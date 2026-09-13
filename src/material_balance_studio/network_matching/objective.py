"""Fixed-length network pressure objective; every trial starts a fresh Phase 5A run."""
from dataclasses import dataclass
import numpy as np
from material_balance_studio.tank_network import simulate_network,NetworkSettings
from material_balance_studio.matching.weighting import residual_scales
from material_balance_studio.network_history import prepare_history
from .parameters import validate_specs,candidate_network
from .diagnostics import closure_summary,data_qc


@dataclass(frozen=True)
class NetworkEvaluation:
    parameters: tuple
    objective: float
    valid: bool
    failure_type: str | None
    failure: str | None
    affected: tuple
    closure: tuple


class NetworkObjective:
    def __init__(self,network,specs,observations=None,plan=None,mode="Unweighted",default_sigma=None,settings=None):
        self.base=network
        self.specs=tuple(specs)
        validate_specs(network,self.specs)
        self.active=tuple(s for s in self.specs if s.active)
        self.prepared=prepare_history(network,plan,observations)
        self.observations=tuple(sorted((o for o in self.prepared.observations if o.include_in_match),key=lambda o:(o.date,o.tank)))
        self.warnings=data_qc(network,specs,self.observations)
        self.mode,self.default_sigma=mode,default_sigma
        self.scales=residual_scales(self.observations,mode,default_sigma)
        self.settings=settings or NetworkSettings()
        standard=NetworkSettings()
        for name in ("absolute_tolerance","relative_tolerance","transfer_absolute_tolerance","transfer_relative_tolerance"):
            if getattr(self.settings,name)>getattr(standard,name):
                raise ValueError("Matching cannot loosen Phase 5A closure tolerances.")
        if self.settings.normalization_floor!=standard.normalization_floor:
            raise ValueError("Matching retains the Phase 5A normalization floor.")
        bounds={t.name:t.pressure_bounds for t in network.tanks}
        worst=[max(abs(p-o.pressure) for p in bounds[o.tank])/scale for o,scale in zip(self.observations,self.scales)]
        self.penalty=1000.*max(1.,max(worst))
        self.records=[]
        self.last_result=None
        self.last_network=None

    def __call__(self,vector):
        self.last_result=self.last_network=None
        values=()
        failure=None
        kind=None
        affected=()
        closure=()
        try:
            if len(vector)!=len(self.active):
                raise ValueError("Network parameter vector length mismatch.")
            values=tuple(s.decode(float(x)) for s,x in zip(self.active,vector))
            candidate=candidate_network(self.prepared.network,self.specs,values)
            result=simulate_network(candidate,self.settings)
            if not result.converged:
                diag=result.failed_timestep.diagnostics
                affected=tuple(diag.dominant_residuals)+tuple(diag.dominant_connections)
                message=diag.message.lower()
                kind="stiffness/work-limit failure" if any(w in message for w in ("stiffness","substep limit","time resolution")) else "aquifer failure" if "aquifer" in message else "PVT failure" if "pvt" in message else "coupled solver failure"
                failed=result.failed_timestep.candidate
                if failed and any(min(abs(t.pressure-lo),abs(t.pressure-hi))<1. for t,(lo,hi) in zip(failed.tanks,diag.pressure_bounds)):
                    kind="pressure-bound failure"
                raise ValueError(diag.message)
            qc=closure_summary(result)
            closure=tuple(qc.items())
            if not qc["valid"]:
                kind="closure failure"
                affected=tuple((t.name,t.residual,t.relative_residual) for s in result.internal_states for t in s.tanks if abs(t.residual)>1e-5 or t.relative_residual>1e-8)
                raise ValueError("Network candidate violates mandatory per-tank MBE or transfer-conservation acceptance.")
            calculated={(s.time.date(),t.name):t.pressure for s in (result.initial_state,*result.states) for t in s.tanks}
            residual=np.array([calculated[o.date,o.tank]-o.pressure for o in self.observations])/self.scales
            if not np.all(np.isfinite(residual)):
                raise ValueError("Nonfinite network pressure residual.")
            self.last_result,self.last_network=result,candidate
        except (ValueError,ArithmeticError,KeyError) as exc:
            failure=str(exc)
            message=failure.lower()
            kind=kind or ("stiffness/work-limit failure" if "stiffness" in message or "work-limit" in message else "aquifer failure" if "aquifer" in message else "PVT failure" if "pvt" in message else "coupled solver failure" if "coupled solver" in message else type(exc).__name__)
            residual=np.full(len(self.observations),self.penalty)
        self.records.append(NetworkEvaluation(values,float(residual@residual),failure is None,kind,failure,affected,closure))
        return residual
