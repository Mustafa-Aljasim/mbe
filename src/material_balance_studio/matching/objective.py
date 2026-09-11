"""A fixed-length pressure residual vector; MBE closure is a separate validity gate."""
from dataclasses import dataclass
import numpy as np
from material_balance_studio.solver.simulation import simulate
from material_balance_studio.domain.models import SolverSettings
from .parameters import candidate_tank,validate_parameters
from .weighting import residual_scales


@dataclass(frozen=True)
class Evaluation:
    parameters: tuple
    objective: float
    valid: bool
    failure: str | None
    maximum_mbe_relative: float | None


class PressureObjective:
    def __init__(self,tank,history,specs,observations,mode="Unweighted",default_sigma=None):
        self.tank,self.history,self.specs=tank,tuple(history),tuple(specs)
        validate_parameters(tank,self.specs)
        self.active=tuple(s for s in specs if s.active)
        if len({o.date for o in observations})!=len(observations):
            raise ValueError("Duplicate pressure observation dates.")
        actual={h.date:h.observed_pressure for h in history}
        for o in observations:
            if o.date not in actual or actual[o.date] is None or actual[o.date]!=o.pressure:
                raise ValueError("Observations must match exact existing history dates and pressures.")
        self.observations=tuple(sorted((o for o in observations if o.include_in_match),key=lambda o:o.date))
        if len(self.observations)<=len(self.active):
            raise ValueError("FAIL: included observations must outnumber active parameters.")
        self.scales=residual_scales(self.observations,mode,default_sigma)
        self.mode,self.default_sigma=mode,default_sigma
        low,high=tank.pvt_model.pressure_bounds
        possible=np.array([max(abs(low-o.pressure),abs(high-o.pressure)) for o in self.observations])/self.scales
        self.penalty=float(max(1.,np.max(possible))*1000.)
        self.records=[]
        self.last_result=None
        self.last_tank=None

    def __call__(self,vector):
        self.last_result,self.last_tank=None,None
        values=()
        failure=None
        maximum=None
        try:
            if len(vector)!=len(self.active):
                raise ValueError("Parameter vector length mismatch.")
            values=tuple(s.decode(float(v)) for s,v in zip(self.active,vector))
            candidate=candidate_tank(self.tank,self.specs,values)
            result=simulate(candidate,self.history)  # fresh initial aquifer state on EVERY call
            if not result.converged:
                raise ValueError(f"Pressure solver failed at {result.failed_timestep.date}: {result.failed_timestep.solution.diagnostics.message}")
            settings=SolverSettings()
            for state in result.states:
                b=state.balance
                if not np.isfinite(b.relative_residual) or not np.isfinite(b.absolute_residual) or b.relative_residual>settings.relative_residual_tolerance or b.absolute_residual>settings.absolute_residual_tolerance:
                    raise ValueError("Forward MBE closure violates the existing solver tolerances.")
            maximum=max((s.balance.relative_residual for s in result.states),default=0.)
            by_date={s.date:s.pressure for s in result.states}
            residual=np.array([by_date[o.date]-o.pressure for o in self.observations])/self.scales
            if not np.all(np.isfinite(residual)):
                raise ValueError("Nonfinite history-match residual.")
            self.last_result,self.last_tank=result,candidate
        except (ValueError,ArithmeticError,KeyError) as exc:
            failure=str(exc)
            residual=np.full(len(self.observations),self.penalty)
        self.records.append(Evaluation(values,float(np.dot(residual,residual)),failure is None,failure,maximum))
        return residual
