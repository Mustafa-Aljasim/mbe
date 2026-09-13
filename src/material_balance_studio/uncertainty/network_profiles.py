"""Real network objective grids and bounded nuisance-parameter profiles."""
from types import SimpleNamespace
from dataclasses import replace
from .network_context import context, fixed_point, fit_point, NetworkPoint
from .context import grid, objective_label
from .result import SurfaceResult, ProfileResult
from .confidence import profile_threshold, grid_ranges
from .network_sensitivity import analyze_network
from material_balance_studio.network_matching import history_match_network
from material_balance_studio.network_matching.parameters import candidate_network


def network_surface(scenario, parameter_x, parameter_y, count=5):
    specs,obj=context(scenario)
    names=[s.name for s in specs]
    if parameter_x==parameter_y or any(n not in names for n in (parameter_x,parameter_y)):
        raise ValueError('Select two distinct active network parameters.')
    ix,iy=names.index(parameter_x),names.index(parameter_y)
    values=[s.initial for s in specs]
    xs,ys=grid(specs[ix],values[ix],count),grid(specs[iy],values[iy],count)
    points=[]
    for y in ys:
        for x in xs:
            trial=values.copy()
            trial[ix],trial[iy]=x,y
            points.append(fixed_point(obj,trial))
    initial={s.name:s.initial for s in scenario.specifications}
    return SurfaceResult(scenario.scenario_id,(parameter_x,parameter_y),xs,ys,tuple(points),objective_label(scenario.mode),
        (initial[parameter_x],initial[parameter_y]),(values[ix],values[iy]),
        ('Conditional fixed-allocation surface: other parameters fixed at the match; valleys suggest compensation. Failed nodes remain gaps.',),
        sum(not p.valid for p in points))


def network_profile(scenario, parameter, count=5, reoptimize=True, max_nfev=40,
                    statistical_assumptions=False, plausible_delta=None):
    specs,obj=context(scenario)
    names=[s.name for s in specs]
    if parameter not in names:
        raise ValueError('Select an active network parameter.')
    if type(max_nfev) is not int or not 1<=max_nfev<=1000:
        raise ValueError('max_nfev must be an integer from 1 to 1000.')
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
                network=candidate_network(scenario.base_network,specs,trial)
                fit=history_match_network(network,remaining,scenario.observations,scenario.history_plan,
                    scenario.mode,scenario.default_sigma,scenario.settings,max_nfev=max_nfev,
                    scenario_name=f'Profile {parameter}={value}')
                fitted=dict(fit.fitted_parameters)
                fitted[parameter]=value
                point=fit_point(fit,tuple((s.name,fitted[s.name]) for s in specs))
                fixed=specs[index]
                point=replace(point,bound_flags=point.bound_flags+((parameter,(value-fixed.lower)/(fixed.upper-fixed.lower)<1e-4,(fixed.upper-value)/(fixed.upper-fixed.lower)<1e-4),))
                failed+=fit.failed_candidates
            except (ValueError,ArithmeticError) as exc:
                point=NetworkPoint(tuple(zip(names,trial)),None,None,False,'FAILED',(),str(exc),0,(),(),(),type(exc).__name__)
                failed+=1
        points.append(point)
    proxy=SimpleNamespace(final_objective=scenario.final_objective,weighting_mode=scenario.mode,
                          observations=scenario.observations,specifications=scenario.specifications)
    threshold,method=profile_threshold(proxy,statistical_assumptions,plausible_delta)
    warnings=['Fixed allocation: grid ranges are distinct from allocation-assumption sensitivity. Failed/non-converged nodes remain unresolved; no interpolation.']
    if statistical_assumptions and plausible_delta is None:
        info=analyze_network(scenario)
        local=info.local
        if local.rank<len(specs) or any(lo or hi for _,lo,hi in scenario.bound_flags) or max(local.derivative_errors)>.05 or (local.condition is not None and local.condition>1e4) or any(c[1]!='DIRECTLY INFORMED' for c in info.coverage):
            threshold=None
            method='Statistical range withheld: rank, bounds, conditioning or derivative stability fails.'
        failed+=local.failed_evaluations
    if not reoptimize and remaining:
        method='Conditional fixed-others slice; '+method
        warnings.append('Fixed-others ranges ignore compensation and may understate uncertainty.')
    ranges=grid_ranges(points,parameter,threshold)
    brackets=[]
    if threshold is not None:
        for a,b in zip(points,points[1:]):
            if a.valid and b.valid and all(p.optimizer_status in ('FIXED OTHERS','CONVERGED') for p in (a,b)) and (a.objective<=threshold)!=(b.objective<=threshold):
                brackets.append((dict(a.parameters)[parameter],dict(b.parameters)[parameter]))
    if any(a==b for a,b in ranges):
        warnings.append('Single-node plausible width is unresolved, not zero; refine the grid.')
    if any(a==specs[index].lower or b==specs[index].upper for a,b in ranges):
        warnings.append('Plausible region reaches a user bound; its external extent is unknown.')
    if any(not p.valid or p.optimizer_status=='NOT CONVERGED' for p in points):
        warnings.append('Failed or non-converged regions cannot define plausible endpoints.')
    if any(p.valid and p.objective<scenario.final_objective-max(1e-8,abs(scenario.final_objective)*1e-6) for p in points):
        threshold,ranges,brackets=None,(),[]
        method='Withheld: profile found a better optimum; revisit the accepted match.'
        warnings.append(method)
    return ProfileResult(scenario.scenario_id,parameter,'Re-optimized' if reoptimize else 'Fixed others',tuple(points),
                         threshold,method,ranges,tuple(warnings),failed,tuple(brackets))
