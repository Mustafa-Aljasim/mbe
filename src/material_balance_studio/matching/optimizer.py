"""Deterministic bounded trust-region least squares; no base-model writes."""
from math import sqrt
import numpy as np
from scipy.optimize import least_squares
from material_balance_studio.aquifer import NoAquifer
from material_balance_studio.diagnostics.diagnosis import diagnose
from material_balance_studio.diagnostics.pressure_qc import PressureObservation
from material_balance_studio.diagnostics.aquifer_comparison import ComparisonRun,engineering_qc
from .objective import PressureObjective
from .result import HistoryMatchResult,pvt_fingerprint


def metrics(errors):
    if not errors:
        return dict(rmse=None,mae=None,bias=None,maximum=None)
    return dict(rmse=sqrt(sum(e*e for e in errors)/len(errors)),mae=sum(abs(e) for e in errors)/len(errors),
                bias=sum(errors)/len(errors),maximum=max(abs(e) for e in errors))


def history_match(tank,history,specs,observations,mode="Unweighted",default_sigma=None,max_nfev=100):
    if not isinstance(max_nfev,int) or not 1<=max_nfev<=1000:
        raise ValueError("max_nfev must be an integer from 1 to 1000.")
    objective=PressureObjective(tank,history,specs,observations,mode,default_sigma)
    active=objective.active
    x0=[s.encode(s.initial) for s in active]
    bounds=([s.encode(s.lower) for s in active],[s.encode(s.upper) for s in active])
    initial=objective(x0)
    initial_result=objective.last_result
    fitted=least_squares(objective,x0,bounds=bounds,method="trf",jac="3-point",x_scale=1.,
                         ftol=1e-10,xtol=1e-10,gtol=1e-10,max_nfev=max_nfev)
    final=objective(fitted.x)
    final_result,fitted_tank=objective.last_result,objective.last_tank
    valid=final_result is not None
    converged=bool(fitted.success and valid)
    warnings=[]
    nobs=len(objective.observations)
    if nobs<max(5,2*len(active)+1):
        warnings.append("Weakly overdetermined match; few observations relative to active parameters.")
    pressures=[o.pressure for o in objective.observations]
    if (tank.initial_pressure-min(pressures))/tank.initial_pressure<.01 or np.ptp(pressures)<.01*tank.initial_pressure:
        warnings.append("Pressure depletion/spread below 1% of Pi; parameter recovery may be poorly constrained.")
    if any(o.qc_flag or (o.source and o.source!="Average reservoir pressure") for o in objective.observations):
        warnings.append("Included observations carry quality/source metadata; inspect representativeness.")
    values=tuple((s.name,s.decode(float(v))) for s,v in zip(active,fitted.x))
    flags=[]
    for s,(_,v) in zip(active,values):
        lower=(v-s.lower)/(s.upper-s.lower)<1e-4
        upper=(s.upper-v)/(s.upper-s.lower)<1e-4
        flags.append((s.name,bool(lower),bool(upper)))
        if lower or upper:
            warnings.append(f"{s.name} constrained by {'lower' if lower else 'upper'} bound; not well determined.")
    if not fitted.success:
        warnings.append("Optimizer did not converge: "+fitted.message)
    if not valid:
        warnings.append("Final forward candidate invalid: "+str(objective.records[-1].failure))
    failures=sum(not r.valid for r in objective.records)
    if failures:
        warnings.append(f"{failures} candidate evaluations failed; inspect recorded failure reasons.")
    by_initial={s.date:s.pressure for s in initial_result.states} if initial_result else {}
    by_final={s.date:s.pressure for s in final_result.states} if final_result else {}
    included={o.date for o in objective.observations}
    scales=dict(zip((o.date for o in objective.observations),objective.scales))
    series=[]
    for o in observations:
        before=by_initial.get(o.date)
        after=by_final.get(o.date)
        sigma=o.sigma if o.sigma is not None and np.isfinite(o.sigma) and o.sigma>0 else default_sigma
        series.append(dict(date=str(o.date),observed=o.pressure,initial=before,matched=after,included=o.date in included,
                           initial_residual=before-o.pressure if before is not None else None,
                           residual=after-o.pressure if after is not None else None,
                           normalized=(after-o.pressure)/sigma if after is not None and sigma is not None and sigma>0 else None,
                           sigma=sigma,source=o.source,qc_flag=o.qc_flag,note=o.note))
    initial_metrics=metrics([r["initial_residual"] for r in series if r["included"] and r["initial_residual"] is not None])
    final_metrics=metrics([r["residual"] for r in series if r["included"] and r["residual"] is not None])
    errors=np.array([r["residual"] for r in series if r["included"] and r["residual"] is not None])
    if len(errors) and final_metrics["rmse"]>1.:
        if abs(final_metrics["bias"])>.5*final_metrics["rmse"]:
            warnings.append("Systematic pressure residual bias exceeds half the RMSE.")
        if final_metrics["maximum"]>2*final_metrics["rmse"]:
            warnings.append("One dominant residual / large maximum error despite lower RMSE.")
        if len(errors)>=4 and np.mean(errors[:len(errors)//2])*np.mean(errors[len(errors)//2:])<0:
            warnings.append("Early/late residual means have opposite signs; inspect systematic temporal trend.")
    if any(abs(r["normalized"])>3 for r in series if r["included"] and r["normalized"] is not None):
        warnings.append("Included normalized pressure residual exceeds 3 sigma.")
    metadata=tuple(PressureObservation(str(o.date),o.source,o.sigma,o.qc_flag,o.note) for o in observations)
    diagnostics=diagnose(final_result,fitted_tank,history,metadata) if valid else None
    ho=diagnostics["ho"]["fit"]["slope"] if diagnostics else None
    difference=100*(ho-fitted_tank.oil_in_place)/fitted_tank.oil_in_place if ho is not None and valid else None
    if difference is not None and abs(difference)>20:
        warnings.append("H–O and matched OOIP differ by more than 20%; review assumptions.")
    if diagnostics and diagnostics["ho"]["fit"]["status"]!="PASS":
        warnings.append("H–O cross-check carries diagnostic cautions; inspect its screening notes.")
    if valid:
        _,engineering_notes=engineering_qc(ComparisonRun("Matched scenario",fitted_tank.aquifer or NoAquifer(),final_result))
        warnings.extend(engineering_notes)
    maximum=max((r.maximum_mbe_relative for r in objective.records if r.valid),default=0.)
    return HistoryMatchResult(tank,fitted_tank,tank.aquifer.key if tank.aquifer else "none",pvt_fingerprint(tank.pvt_model),
        tuple(specs),tuple(observations),mode,default_sigma,values,int(fitted.status),converged,str(fitted.message),
        len(objective.records),int(fitted.nfev),float(np.dot(initial,initial)),float(np.dot(final,final)),initial_metrics,final_metrics,
        initial_result,final_result,tuple(series),tuple(objective.records),maximum,tuple(flags),
        "FAIL" if not converged else "CAUTION" if warnings else "PASS",tuple(warnings),diagnostics,ho,difference)
