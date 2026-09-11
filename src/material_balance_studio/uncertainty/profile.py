"""Conditional slices or deterministic nuisance-parameter re-optimization."""
from material_balance_studio.matching import history_match
from material_balance_studio.matching.parameters import candidate_tank
from .context import context,grid
from .objective_surface import fixed_point
from .result import GridPoint,ProfileResult,scenario_identity
from .confidence import profile_threshold,grid_ranges


def parameter_profile(scenario,parameter,count=7,reoptimize=True,max_nfev=40,
                      statistical_assumptions=False,plausible_delta=None):
    history,specs,obj=context(scenario)
    names=[s.name for s in specs]
    if parameter not in names:
        raise ValueError("Select an active match parameter.")
    if not isinstance(max_nfev,int) or not 1<=max_nfev<=1000:
        raise ValueError("max_nfev must be an integer from 1 to 1000.")
    index=names.index(parameter)
    values=[s.initial for s in specs]
    remaining=tuple(s for s in specs if s.name!=parameter)
    points=[]
    failed=0
    for value in grid(specs[index],values[index],count):
        trial=values.copy()
        trial[index]=value
        if not reoptimize or not remaining:
            point=fixed_point(obj,trial)
            failed+=not point.valid
        else:
            try:
                base=candidate_tank(scenario.base_tank,specs,trial)
                fit=history_match(base,history,remaining,scenario.observations,scenario.weighting_mode,
                                  scenario.default_sigma,max_nfev=max_nfev)
            except (ValueError,ArithmeticError) as exc:
                points.append(GridPoint(tuple((s.name,float(v)) for s,v in zip(specs,trial)),
                    None,None,False,"FAILED",(),str(exc),0))
                failed+=1
                continue
            fitted=dict(fit.fitted_parameters)
            fitted[parameter]=value
            valid=fit.final_simulation is not None and fit.final_simulation.converged
            point=GridPoint(tuple((s.name,fitted[s.name]) for s in specs),
                fit.final_objective if valid else None,fit.final_metrics["rmse"] if valid else None,
                valid,"CONVERGED" if fit.converged else "NOT CONVERGED",fit.bound_flags,
                None if fit.converged else fit.termination_message,fit.function_evaluations)
            failed+=sum(not e.valid for e in fit.evaluations)
        points.append(point)
    threshold,method=profile_threshold(scenario,statistical_assumptions,plausible_delta)
    warnings=["Ranges resolve accepted grid nodes only; refine resolution before interpreting endpoints. No interpolation across failed or unresolved nodes."]
    if statistical_assumptions and plausible_delta is None:
        from .sensitivity import analyze
        local=analyze(scenario)
        if local.rank<len(specs) or any(lo or hi for _,lo,hi in scenario.bound_flags) or max(local.derivative_errors)>.05 or (local.condition is not None and local.condition>1e4):
            threshold=None
            method="Statistical threshold withheld: rank, conditioning, derivative or boundary regularity fails; use an explicit plausible increment"
        failed+=local.failed_evaluations
    if not reoptimize and remaining:
        method="Conditional fixed-others slice; "+method
        warnings.append("Fixed-others ranges ignore parameter compensation and can understate uncertainty.")
    ranges=grid_ranges(points,parameter,threshold)
    brackets=[]
    if threshold is not None:
        for a,b in zip(points,points[1:]):
            if a.valid and b.valid and a.optimizer_status in ("FIXED OTHERS","CONVERGED") and b.optimizer_status in ("FIXED OTHERS","CONVERGED"):
                if (a.objective<=threshold)!=(b.objective<=threshold):
                    brackets.append((dict(a.parameters)[parameter],dict(b.parameters)[parameter]))
    if any(b-a<=1e-12*(specs[index].upper-specs[index].lower) for a,b in ranges):
        warnings.append("A plausible region contains only one sampled node: its width is unresolved, not zero. Inspect threshold-crossing brackets or increase resolution.")
    if any(min(abs(a-specs[index].lower),abs(b-specs[index].upper))<=1e-4*(specs[index].upper-specs[index].lower) for a,b in ranges):
        warnings.append("Plausible region reaches a user bound; uncertainty may extend beyond the permitted range and is not determined by this constrained profile.")
    if any(not p.valid or p.optimizer_status=="NOT CONVERGED" for p in points):
        warnings.append("Failed or non-converged nodes do not define confidence limits; their regions remain unresolved.")
    if any(p.valid and p.objective<scenario.final_objective-max(1e-8,abs(scenario.final_objective)*1e-6) for p in points):
        warnings.append("Profile improves on the accepted match: inspect the original optimum before interpreting ranges.")
        threshold=None
        ranges=()
        brackets=[]
        method="Withheld: profile found a better optimum than the accepted match"
    return ProfileResult(scenario_identity(scenario),parameter,"Re-optimized" if reoptimize else "Fixed others",
                         tuple(points),threshold,method,ranges,tuple(warnings),failed,tuple(brackets))
