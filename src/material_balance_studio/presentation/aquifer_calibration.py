"""Aquifer comparison orchestration around the validated single-tank matcher."""
from dataclasses import dataclass,replace
import numpy as np
from material_balance_studio.diagnostics.aquifer_comparison import ComparisonRun,compare_aquifers,engineering_qc
from material_balance_studio.matching import history_match,observations_from_history
from material_balance_studio.matching.parameters import validate_parameters
from material_balance_studio.matching.weighting import residual_scales
from material_balance_studio.matching.optimizer import metrics


@dataclass(frozen=True)
class CalibratedComparison:
    run: ComparisonRun
    fit: object | None
    pressure_metrics: tuple
    objective: float | None
    qc: str
    warnings: tuple
    mode: str = 'History-Match Each Aquifer Before Comparison'


def calibrate_aquifers(tank,history,models,specifications,observations=None,mode='Unweighted',default_sigma=None,max_nfev=100):
    history=tuple(history)
    observations=observations_from_history(history) if observations is None else tuple(observations)
    included=tuple(sorted((o for o in observations if o.include_in_match),key=lambda o:o.date))
    # Exact-date and sigma rules are the same boundary rules as PressureObjective.
    actual={h.date:h.observed_pressure for h in history}
    if len({o.date for o in observations})!=len(observations) or any(o.date not in actual or actual[o.date] is None or actual[o.date]!=o.pressure for o in observations):
        raise ValueError('Observations must match exact existing history dates and pressures, without duplicates.')
    if not included:
        raise ValueError('Include at least one observed pressure for aquifer comparison matching.')
    scales=residual_scales(included,mode,default_sigma)
    output=[]
    for name,model in models.items():
        candidate=replace(tank,aquifer=model)
        specs=tuple(specifications.get(name,()))
        try:
            if any(not s.name.startswith('aquifer.') for s in specs):
                raise ValueError('Aquifer comparison may adjust aquifer parameters only; N and m remain fixed.')
            if specs: validate_parameters(candidate,specs)
            if model.key=='none' and specs:
                raise ValueError('No Aquifer has no aquifer parameters to optimize.')
            if specs:
                fit=history_match(candidate,history,specs,observations,mode,default_sigma,max_nfev)
                fitted_model=fit.fitted_tank.aquifer if fit.fitted_tank is not None else model
                run=ComparisonRun(name,fitted_model,fit.final_simulation,None if fit.converged else '; '.join(fit.warnings))
                status,notes=engineering_qc(run)
                warnings=tuple(dict.fromkeys((*notes,*fit.warnings)))
                status='FAIL' if not fit.converged or status=='FAIL' else 'CAUTION' if warnings else 'PASS'
                output.append(CalibratedComparison(run,fit,tuple(fit.final_metrics.items()),fit.final_objective if fit.converged else None,status,warnings))
            else:
                run=compare_aquifers(tank,history,{name:model})[0]
                status,notes=engineering_qc(run)
                bydate={s.date:s.pressure for s in run.result.states} if run.result else {}
                errors=[bydate[o.date]-o.pressure for o in included if o.date in bydate]
                valid=run.result is not None and run.result.converged and len(errors)==len(included) and status!='FAIL'
                objective=float(np.sum((np.array(errors)/scales)**2)) if valid else None
                if model.key!='none': notes=(*notes,'No aquifer parameters selected: retained as a fixed supplied-parameter run.')
                output.append(CalibratedComparison(run,None,tuple(metrics(errors).items()),objective,status,tuple(notes)))
        except (ValueError,ArithmeticError) as exc:
            output.append(CalibratedComparison(ComparisonRun(name,model,None,str(exc)),None,tuple(metrics([]).items()),None,'FAIL',(str(exc),)))
    return tuple(output)
