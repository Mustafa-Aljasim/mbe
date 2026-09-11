"""Conditional local intervals and grid-resolved profile ranges, never posteriors."""
import numpy as np
from scipy.stats import norm,t,chi2,f


def local_confidence(geometry,specs,values,sse,nobs,weighted,assumptions,bound_flags,stable=True):
    empty=dict(covariance=None,se=tuple(None for _ in specs),intervals=tuple(None for _ in specs),warnings=[])
    if not assumptions:
        empty["warnings"]=["Statistical assumptions not accepted: confidence intervals withheld."]
        return empty
    if geometry is None or not stable:
        empty["warnings"]=["Rank deficiency, excessive conditioning or unstable derivatives: covariance/intervals withheld; inspect profiles."]
        return empty
    dof=nobs-len(specs)
    if dof<=0 or (not weighted and sse<=1e-12):
        empty["warnings"]=["Residual variance cannot be estimated reliably (no degrees of freedom or effectively noiseless fit). Use a supplied plausible threshold."]
        return empty
    variance=1. if weighted else sse/dof
    scales=np.array([s.scale for s in specs])
    covariance=geometry*np.outer(scales,scales)*variance
    se=np.sqrt(np.maximum(0,np.diag(covariance)))
    multiplier=norm.ppf(.975) if weighted else t.ppf(.975,dof)
    intervals=[]
    warnings=[]
    for spec,value,error,bound in zip(specs,values,se,bound_flags):
        low,high=value-multiplier*error,value+multiplier*error
        if bound or low<spec.lower or high>spec.upper:
            intervals.append(None)
            warnings.append(f"{spec.name}: symmetric interval reaches/exceeds a bound; withheld, profile the valid range.")
        else:
            intervals.append((float(low),float(high)))
    return dict(covariance=covariance,se=tuple(float(v) for v in se),intervals=tuple(intervals),warnings=warnings)


def profile_threshold(scenario,assumptions=False,plausible_delta=None):
    if plausible_delta is not None:
        if not np.isfinite(plausible_delta) or plausible_delta<=0:
            raise ValueError("Plausible objective increment must be positive and finite.")
        return scenario.final_objective+plausible_delta,"User-specified objective increment; no confidence probability"
    if not assumptions:
        return None,"No statistical threshold: assumptions were not accepted"
    if scenario.weighting_mode=="Pressure uncertainty weighted":
        return scenario.final_objective+float(chi2.ppf(.95,1)),"Approximate 95% profile likelihood: known independent Gaussian sigma; Δχ²(1)"
    dof=sum(o.include_in_match for o in scenario.observations)-sum(s.active for s in scenario.specifications)
    if dof<=0 or scenario.final_objective<=1e-12:
        return None,"Unknown variance not estimable; supply a deterministic plausible increment"
    return scenario.final_objective*(1+float(f.ppf(.95,1,dof))/dof),"Approximate 95% extra-sum-of-squares F(1,n−p), estimated common variance"


def grid_ranges(points,parameter,threshold):
    """Separate grid-resolved regions; never bridge failed or unresolved nodes."""
    if threshold is None:
        return ()
    ranges=[]
    start=end=None
    for point in points:
        x=dict(point.parameters)[parameter]
        inside=point.valid and point.optimizer_status in ("FIXED OTHERS","CONVERGED") and point.objective<=threshold
        if inside:
            if start is None:
                start=x
            end=x
        elif start is not None:
            ranges.append((start,end))
            start=end=None
    if start is not None:
        ranges.append((start,end))
    return tuple(ranges)
