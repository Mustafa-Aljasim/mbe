"""Bounded simultaneous pressures with conservative trapezoidal communication."""
from dataclasses import replace
from datetime import datetime,time,timedelta
import numpy as np
from scipy.optimize import least_squares
from material_balance_studio.mbe.oil_balance import evaluate_balance
from .network import NetworkSettings
from .history import timeline
from .residuals import initial_state,evaluate_candidate,accepted
from .state import NetworkDiagnostics,NetworkFailure,NetworkResult


def storage_slopes(network,previous):
    """Local reservoir d(Fnet-N Et)/dP only; aquifer storage is not counted twice."""
    slopes=[]
    for node,state in zip(network.tanks,previous.tanks):
        lo,hi=node.pressure_bounds
        h=min(max(1.,abs(state.pressure)*1e-5),.1*(hi-lo))
        a,b=max(lo,state.pressure-h),min(hi,state.pressure+h)
        fa=evaluate_balance(a,node.reservoir,state.cumulative).residual
        fb=evaluate_balance(b,node.reservoir,state.cumulative).residual
        slopes.append((fb-fa)/(b-a))
    return np.array(slopes)


def communication_step(network,previous,settings):
    slopes=storage_slopes(network,previous)
    degrees={t.name:0. for t in network.tanks}
    for e in network.connections:
        degrees[e.from_tank]+=e.effective_transmissibility
        degrees[e.to_tank]+=e.effective_transmissibility
    rates=[]
    for t,c in zip(network.tanks,slopes):
        if degrees[t.name]>0 and c<=0:
            raise ValueError(f"{t.name}: nonpositive local storage slope; communication timestep cannot be certified.")
        rates.append(degrees[t.name]/c if degrees[t.name]>0 else 0.)
    bound=2*max(rates,default=0.)
    limit=settings.communication_step_limit/bound if bound>0 else float("inf")
    if settings.max_step_days is not None:
        limit=min(limit,settings.max_step_days*86400)
    return limit,slopes


def solve_step(network,previous,elapsed_seconds,settings,slopes=None):
    bounds=np.array([t.pressure_bounds for t in network.tanks],dtype=float)
    guess=np.array([s.pressure for s in previous.tanks])
    if slopes is None:
        slopes=storage_slopes(network,previous)
    scales=np.maximum.reduce([np.abs(slopes)*1e6,np.array([abs(s.components.withdrawal.net) for s in previous.tanks]),np.ones(len(guess))])
    evaluations=0
    failures=[]
    last=None
    def objective(x):
        nonlocal evaluations,last
        evaluations+=1
        last=None
        try:
            pressure=np.clip(np.asarray(x)*1e6,bounds[:,0],bounds[:,1])
            last=evaluate_candidate(network,previous,pressure,elapsed_seconds,settings)
            vector=np.array([s.residual for s in last.tanks])/scales
            if not np.all(np.isfinite(vector)):
                raise ValueError("Nonfinite coupled residual.")
            return vector
        except (ValueError,ArithmeticError) as exc:
            failures.append(str(exc))
            last=None
            return np.full(len(guess),1e12)
    x=guess/1e6
    residual=objective(x)
    nfev=0
    message="Previous accepted pressures already satisfy this timestep."
    if last is None or not accepted(last,settings):
        fit=least_squares(objective,x,bounds=(bounds[:,0]/1e6,bounds[:,1]/1e6),method="trf",jac="3-point",
                          x_scale="jac",ftol=1e-14,xtol=1e-14,gtol=1e-14,max_nfev=settings.max_nfev)
        x=fit.x
        nfev=fit.nfev
        message=str(fit.message)
        residual=objective(x)
        # Residual-based Newton polishing within the same coupled equations.
        # This never accepts a least-squares minimum unless every balance closes.
        for _ in range(4):
            if last is None or accepted(last,settings) or not fit.success:
                break
            columns=[]
            for j in range(len(x)):
                h=max(1e-6,abs(x[j])*1e-6)
                a,b=x.copy(),x.copy()
                a[j]=max(bounds[j,0]/1e6,x[j]-h)
                b[j]=min(bounds[j,1]/1e6,x[j]+h)
                ra,rb=objective(a),objective(b)
                columns.append((rb-ra)/(b[j]-a[j]))
            try:
                correction=np.linalg.solve(np.column_stack(columns),residual)
            except np.linalg.LinAlgError:
                break
            trial=np.clip(x-correction,bounds[:,0]/1e6,bounds[:,1]/1e6)
            candidate_residual=objective(trial)
            if np.linalg.norm(candidate_residual)>=np.linalg.norm(residual):
                break
            x,residual=trial,candidate_residual
        residual=objective(x)
    valid=last is not None and accepted(last,settings)
    if not valid:
        message="Coupled pressure solve failed individual/network acceptance. "+message
        if failures:
            message+=" Last candidate failure: "+failures[-1]
    dominant=tuple(sorted(((s.name,s.residual,s.relative_residual) for s in last.tanks),key=lambda x:abs(x[1]),reverse=True)) if last else ()
    connections=tuple(sorted(((e.from_tank,e.to_tank,e.incremental_transfer) for e in last.connections),key=lambda x:abs(x[2]),reverse=True)) if last else ()
    diagnostics=NetworkDiagnostics(valid,evaluations,nfev,tuple(float(v) for v in guess),tuple(map(tuple,bounds)),
        float(np.linalg.norm([s.residual for s in last.tanks])) if last else None,
        max(s.relative_residual for s in last.tanks) if last else None,message,len(failures),dominant,connections)
    if last is not None:
        last=replace(last,diagnostics=diagnostics)
    if not valid:
        return NetworkFailure(datetime.combine(network.initial_date,time())+timedelta(seconds=elapsed_seconds),diagnostics,last)
    return last


def simulate_network(network,settings=None):
    settings=settings or NetworkSettings()
    start=initial_state(network,settings)
    previous=start
    outputs=[]
    internal=[]
    warnings=["Lumped reservoir-volume communication only: transferred phases/composition are not tracked; each tank keeps its assigned PVT and initial fluid inventory model.",
              "Cumulative surface streams interpolate linearly between supplied dates, start at zero, and hold after the last record. Observed pressures are never interpolated or imposed."]
    events=timeline(network)
    if not events:
        raise ValueError("Supply at least one history/event date after the common initial date, even for closed equilibration.")
    refined=False
    for day in events:
        end=(day-network.initial_date).days*86400.
        steps=0
        while previous.elapsed_seconds<end:
            if steps>=settings.max_substeps_per_event:
                diag=NetworkDiagnostics(False,0,0,tuple(s.pressure for s in previous.tanks),tuple(t.pressure_bounds for t in network.tanks),None,None,
                    "Extreme communication stiffness: substep limit reached. Review transmissibility, storage and event spacing; no uncoupled fallback.")
                return NetworkResult(network,settings,start,tuple(outputs),tuple(internal),NetworkFailure(previous.time,diag),tuple(warnings))
            try:
                limit,slopes=communication_step(network,previous,settings)
                remaining=end-previous.elapsed_seconds
                dt=min(remaining,limit)
                if not np.isfinite(dt) or dt<=0 or previous.elapsed_seconds+dt==previous.elapsed_seconds:
                    raise ValueError("Communication timestep is below numerical time resolution.")
                target=end if dt>=remaining else previous.elapsed_seconds+dt
                refined=refined or dt<remaining
                result=solve_step(network,previous,target,settings,slopes)
            except (ValueError,ArithmeticError) as exc:
                diag=NetworkDiagnostics(False,0,0,tuple(s.pressure for s in previous.tanks),tuple(t.pressure_bounds for t in network.tanks),None,None,str(exc))
                result=NetworkFailure(previous.time,diag)
            if isinstance(result,NetworkFailure):
                return NetworkResult(network,settings,start,tuple(outputs),tuple(internal),result,tuple(warnings))
            # The complete candidate network is committed atomically only here.
            previous=result
            internal.append(result)
            steps+=1
        outputs.append(previous)
    if refined:
        warnings.append(f"Communication/user step control used {len(internal)} accepted substeps for {len(outputs)} event dates. Inspect internal timestep diagnostics; refinement is deterministic, not an error estimate.")
    return NetworkResult(network,settings,start,tuple(outputs),tuple(internal),warnings=tuple(warnings))
