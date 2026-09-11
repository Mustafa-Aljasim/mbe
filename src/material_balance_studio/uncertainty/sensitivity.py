"""Bound-aware physical finite differences of the exact Phase 4B residual vector."""
import numpy as np
from .context import context,evaluate,objective_label
from .result import IdentifiabilityResult,scenario_identity,matrix
from .jacobian import matrix_diagnostics
from .confidence import local_confidence
from .identifiability import classify


def analyze(scenario,statistical_assumptions=False,relative_step=1e-4):
    if not 1e-6<=relative_step<=1e-2:
        raise ValueError("Relative finite-difference step must be between 1e-6 and 1e-2.")
    _,specs,obj=context(scenario)
    values=np.array([s.initial for s in specs])
    base=evaluate(obj,values)
    if base is None:
        raise ValueError("Matched point no longer evaluates successfully.")
    columns,schemes,errors,warnings=[],[],[],[]
    def derivative(j,h):
        s,x=specs[j],values[j]
        # Central if both sides fit; second-order one-sided near a bound.
        if x-h>=s.lower and x+h<=s.upper:
            offsets=(-h,h)
            scheme="central"
        elif x+2*h<=s.upper:
            offsets=(h,2*h)
            scheme="forward (second order)"
        elif x-2*h>=s.lower:
            offsets=(-h,-2*h)
            scheme="backward (second order)"
        else:
            return None,"no feasible stencil"
        results=[]
        for offset in offsets:
            trial=values.copy()
            trial[j]+=offset
            results.append(evaluate(obj,trial))
        if any(r is None for r in results):
            return None,"failed forward stencil"
        if scheme=="central":
            return (results[1]-results[0])/(2*h),scheme
        sign=1 if offsets[0]>0 else -1
        return sign*(-3*base+4*results[0]-results[1])/(2*h),scheme
    for j,s in enumerate(specs):
        h=min(relative_step*max(abs(values[j]),s.scale),.1*(s.upper-s.lower))
        coarse,fine,scheme=None,None,"unavailable"
        for _ in range(4):
            coarse,scheme=derivative(j,h)
            fine,_=derivative(j,h/2)
            if coarse is not None and fine is not None:
                break
            h/=2
        if coarse is None or fine is None:
            raise ValueError(f"Cannot form a valid sensitivity stencil for {s.name}; failed evaluations: {sum(not e.valid for e in obj.records)}.")
        error=float(np.linalg.norm(fine-coarse)/max(np.linalg.norm(fine),1e-20))
        if error>.05:
            warnings.append(f"{s.name}: derivative changes >5% on halving step; local intervals withheld.")
        columns.append(fine)
        schemes.append((s.name,scheme,float(h/2)))
        errors.append(error)
    jac=np.column_stack(columns)
    raw=jac*obj.scales[:,None]
    scales=np.array([s.scale for s in specs])
    dimensionless=raw*scales[None,:]/scenario.base_tank.initial_pressure
    scaled_jac=jac*scales[None,:]
    diag=matrix_diagnostics(scaled_jac)
    if diag["rank"]<len(specs):
        warnings.append("Rank-deficient scaled Jacobian: parameter compensation prevents independent local estimates.")
    elif diag["condition"]>1e4:
        warnings.append("High scaled condition number (>1e4): small pressure changes may cause large parameter changes.")
    if np.max(diag["norms"])/max(np.min(diag["norms"]),1e-30)>1e6:
        warnings.append("Extreme scaled column-magnitude imbalance (>1e6); inspect low-sensitivity parameters.")
    bounds=[min(v-s.lower,s.upper-v)/(s.upper-s.lower)<1e-4 for s,v in zip(specs,values)]
    conf=local_confidence(diag["geometry"],specs,values,float(base@base),len(base),scenario.weighting_mode!="Unweighted",
                          statistical_assumptions,bounds,stable=max(errors)<=.05 and (diag["condition"] is None or diag["condition"]<=1e4))
    warnings.extend(conf["warnings"])
    pressures=[o.pressure for o in obj.observations]
    statuses=classify(specs,values,dimensionless,diag,len(base),np.ptp(pressures)/scenario.base_tank.initial_pressure)
    if max(errors)>.05:
        statuses=tuple((name,"POORLY CONSTRAINED" if status in ("WELL CONSTRAINED","MODERATELY CONSTRAINED") else status,
                        note+" Finite-difference instability limits local inference.") for name,status,note in statuses)
    return IdentifiabilityResult(scenario_identity(scenario),tuple(s.name for s in specs),tuple(str(o.date) for o in obj.observations),
        scenario.weighting_mode,objective_label(scenario.weighting_mode),matrix(raw),matrix(dimensionless),matrix(jac),matrix(scaled_jac),
        tuple(schemes),tuple(errors),tuple(float(v) for v in diag["singular"]),diag["rank"],diag["condition"],tuple(float(v) for v in diag["weak"]),
        matrix(conf["covariance"]) if conf["covariance"] is not None else None,conf["se"],conf["intervals"],
        matrix(diag["correlation"]) if diag["correlation"] is not None else None,matrix(diag["coupling"]),statuses,tuple(warnings),
        sum(not e.valid for e in obj.records),len(obj.records),statistical_assumptions)
