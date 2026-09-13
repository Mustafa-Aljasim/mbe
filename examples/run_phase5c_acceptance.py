"""Reproducible deterministic Phase 5C engineering evidence."""
from pathlib import Path
from dataclasses import asdict
import json,runpy
from material_balance_studio.uncertainty.network_sensitivity import analyze_network
from material_balance_studio.uncertainty.network_profiles import network_surface,network_profile
from material_balance_studio.uncertainty.allocation_sensitivity import shift_allocation,allocation_study

ROOT=Path(__file__).resolve().parents[1]
B=runpy.run_path(str(ROOT/'examples/network_information_benchmarks.py'))


def evidence():
    output={}
    for kind in ('well','NT_confounding','aquifer_confounding','unobserved','low_dp','redundant'):
        fit=B['truth_match'](kind)
        information=analyze_network(fit,True)
        local=information.local
        data=dict(rmse_pa=dict(fit.final_network_metrics)['rmse'],parameters=dict(fit.fitted_parameters),
            rank=local.rank,singular_values=local.singular_values,condition=local.condition,
            correlation=local.correlation,statuses=local.statuses,coverage=information.coverage,
            intervals=local.intervals,weak_combination=local.weak_combination,
            derivative_errors=local.derivative_errors,by_observed_tank=information.by_tank,
            connection_observability=information.connections,warnings=local.warnings)
        if kind in ('well','NT_confounding','aquifer_confounding'):
            names=local.parameters
            surface=network_surface(fit,names[0],'connection:A:B:T',3)
            fixed=network_profile(fit,'connection:A:B:T',3,False,plausible_delta=3.84)
            profile=network_profile(fit,'connection:A:B:T',3,True,plausible_delta=3.84)
            data.update(surface=asdict(surface),fixed_profile=asdict(fixed),reoptimized_profile=asdict(profile))
            if kind=='well':
                data['refined_profile']=asdict(network_profile(fit,'connection:A:B:T',21,True,plausible_delta=3.84))
        output[kind]=data
        print(kind,local.rank,local.correlation,flush=True)
    for kind in ('allocation','allocation_aquifer'):
        fit=B['truth_match'](kind)
        fractions=(.5,.55,.6,.65,.7) if kind=='allocation' else (.5,.6,.7)
        plans=tuple((f'Oil A={f:.2f}',shift_allocation(fit.history_plan,'np','B','A',f-.6)) for f in fractions)
        study=allocation_study(fit,plans)
        output[kind]=dict(robustness=study.robustness,cases=[dict(name=c.name,drift=c.drift,
            rmse_pa=dict(c.match.final_network_metrics).get('rmse') if c.match else None,
            per_tank_metrics=c.match.final_tank_metrics if c.match else (),
            converged=c.match.converged if c.match else False,bounds=c.match.bound_flags if c.match else (),
            drives=c.drives,transfers=c.transfers,robustness=c.robustness,failure=c.failure) for c in study.cases])
        print(kind,output[kind],flush=True)
    return output


if __name__=='__main__':
    record=evidence()
    (ROOT/'docs/phase_5c_numerical_results.json').write_text(json.dumps(record,indent=2,default=str,allow_nan=False)+'\n',encoding='utf-8')
