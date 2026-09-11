"""Independent synthetic match acceptance and forward regression evidence."""
from dataclasses import replace
from pathlib import Path
import json,runpy
from material_balance_studio.matching import history_match,observations_from_history

ROOT=Path(__file__).resolve().parents[1]
helpers=runpy.run_path(str(ROOT/"examples/history_match_benchmarks.py"))
benchmark_case,specification=helpers["benchmark_case"],helpers["specification"]


def evidence():
    previous=json.loads((ROOT/"docs/phase_4a_numerical_results.json").read_text())
    actual=json.loads(json.dumps(runpy.run_path(str(ROOT/"examples/run_phase4a_acceptance.py"))["evidence"]()))
    assert previous==actual
    results={}
    for name,aq,multi in (("OOIP",False,False),("Fetkovich_J",True,False),("OOIP_and_J",True,True)):
        true,history=benchmark_case(aquifer=aq)
        base=replace(true,oil_in_place=8e5 if not aq or multi else 1e6,
                     aquifer=replace(true.aquifer,productivity_index=4e-10) if aq else None)
        specs=[]
        if not aq or multi:
            specs.append(specification(base,"oil_in_place",5e5,1.5e6))
        if aq:
            specs.append(specification(base,"aquifer.productivity_index",1e-10,3e-9,"log"))
        match=history_match(base,history,specs,observations_from_history(history))
        assert match.converged
        results[name]=match
    true,history=benchmark_case()
    noisy=tuple(replace(h,observed_pressure=h.observed_pressure+(2e6 if i==7 else 0)) for i,h in enumerate(history))
    obs=tuple(replace(o,sigma=5e6 if i==7 else 5e4) for i,o in enumerate(observations_from_history(noisy)))
    base=replace(true,oil_in_place=8e5)
    specs=[specification(base,"oil_in_place",5e5,1.5e6)]
    for mode in ("Unweighted","Pressure uncertainty weighted"):
        results[mode]=history_match(base,noisy,specs,obs,mode)
    results["restrictive_bound"]=history_match(base,history,[specification(base,"oil_in_place",5e5,9e5)],observations_from_history(history))
    repeat=history_match(base,history,specs,observations_from_history(history))
    assert repeat.fitted_parameters==results["OOIP"].fitted_parameters
    from material_balance_studio.units.display import to_display,from_display
    field_specs=[replace(s,initial=from_display(to_display(s.initial,s.unit,"FIELD"),s.unit,"FIELD"),lower=from_display(to_display(s.lower,s.unit,"FIELD"),s.unit,"FIELD"),upper=from_display(to_display(s.upper,s.unit,"FIELD"),s.unit,"FIELD")) for s in specs]
    field=history_match(base,history,field_specs,observations_from_history(history))
    return dict(previous_forward_outputs_identical=True,protected_hashes=actual["protected_sources"],
                repeat_parameters_identical=True,field_si_ooip_difference=field.fitted_tank.oil_in_place-repeat.fitted_tank.oil_in_place,
                maximum_valid_candidate_mbe_relative=max(r.maximum_mbe_relative for r in results.values()),
                matches={name:dict(parameters=dict(r.fitted_parameters),converged=r.converged,initial_metrics_pa=r.initial_metrics,
                    final_metrics_pa=r.final_metrics,initial_objective=r.initial_objective,final_objective=r.final_objective,
                    bounds=r.bound_flags,forward_evaluations=r.function_evaluations,ho_ooip=r.ho_ooip,
                    maximum_mbe_relative=r.maximum_mbe_relative) for name,r in results.items()})


if __name__=="__main__":
    output=evidence()
    (ROOT/"docs/phase_4b_numerical_results.json").write_text(json.dumps(output,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in output.items() if k!="protected_hashes"},indent=2))
