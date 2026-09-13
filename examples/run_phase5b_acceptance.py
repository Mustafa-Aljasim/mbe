"""Recover independent network parameters and verify protected numerical outputs."""
from dataclasses import replace
from pathlib import Path
import hashlib,json,runpy
from material_balance_studio.network_history import AllocationRow
from material_balance_studio.network_matching import history_match_network

ROOT=Path(__file__).resolve().parents[1]
helper=runpy.run_path(str(ROOT/"examples/network_match_benchmarks.py"))
benchmark,specification,SETTINGS=helper["benchmark"],helper["specification"],helper["SETTINGS"]


def evidence():
    previous=json.loads((ROOT/"docs/phase_5a_numerical_results.json").read_text())
    current=json.loads(json.dumps(runpy.run_path(str(ROOT/"examples/run_phase5a_acceptance.py"))["evidence"]()))
    assert previous==current,"Protected Phase 5A or earlier numerical outputs changed"
    hashes=json.loads((ROOT/"docs/phase_5b_protected_baseline.json").read_text())
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==v for p,v in hashes.items())
    matches={}
    for case,kwargs,params in (("T",{},(("connection:A:B:T",3e-10),)),
        ("N_A_N_B_T",{},(("tank:A:N",8e5),("tank:B:N",1.2e6),("connection:A:B:T",3e-10))),
        ("aquifer_J_T",dict(aquifer=True),(("tank:B:aquifer:productivity_index",4e-10),("connection:A:B:T",3e-10))),
        ("m_A",dict(gas_cap=.2),(("tank:A:m",.05),))):
        network,obs,plan=benchmark(**kwargs)
        specs=[specification(network,n,initial=v,transformation="log" if n.startswith("connection") or "aquifer" in n else "linear") for n,v in params]
        matches[case]=history_match_network(network,specs,obs,plan,settings=SETTINGS,scenario_name=case)
    network,obs,plan=benchmark(allocation=True)
    specs=[specification(network,"connection:A:B:T",initial=3e-10,transformation="log")]
    matches["correct_allocation"]=history_match_network(network,specs,obs,plan,settings=SETTINGS)
    wrong=replace(plan,schedules=(AllocationRow(network.initial_date,"wp",(("A",.5),("B",.5))),*[s for s in plan.schedules if s.stream!="wp"]))
    matches["wrong_allocation"]=history_match_network(network,specs,obs,wrong,settings=SETTINGS)
    network,obs,plan=benchmark()
    noisy=tuple(replace(o,pressure=o.pressure+(2e6 if i==len(obs)-1 else 0.),sigma=5e6 if i==len(obs)-1 else 2e4) for i,o in enumerate(obs))
    specs=[specification(network,"connection:A:B:T",initial=3e-10,transformation="log")]
    for mode in ("Unweighted","Pressure uncertainty weighted"):
        matches[mode]=history_match_network(network,specs,noisy,plan,mode,settings=SETTINGS)
    matches["restrictive_T"]=history_match_network(network,[specification(network,"connection:A:B:T",upper=5e-10,initial=3e-10)],obs,plan,settings=SETTINGS)
    network,obs,plan=benchmark(transmissibility=0.)
    matches["zero_T"]=history_match_network(network,[specification(network,"connection:A:B:T",lower=0.,initial=3e-10)],obs,plan,settings=SETTINGS)
    for name,result in matches.items():
        assert result.converged,name+" did not converge"
    return dict(previous_phase5a_and_all_prior_outputs_identical=True,protected_sources=hashes,
        maximum_individual_tank_relative_residual_during_valid_matching=max(r.maximum_valid_relative_residual for r in matches.values()),
        maximum_individual_tank_absolute_residual_during_valid_matching=max(r.maximum_valid_absolute_residual for r in matches.values()),
        maximum_transfer_conservation_error=max(r.maximum_transfer_conservation_error for r in matches.values()),
        minimum_closure_margin=min(r.minimum_closure_margin for r in matches.values()),
        matches={name:dict(parameters=dict(r.fitted_parameters),converged=r.converged,initial_objective=r.initial_objective,final_objective=r.final_objective,
            initial_network_metrics_pa=dict(r.initial_network_metrics),final_network_metrics_pa=dict(r.final_network_metrics),
            initial_tank_metrics_pa={n:dict(m) for n,m in r.initial_tank_metrics},final_tank_metrics_pa={n:dict(m) for n,m in r.final_tank_metrics},
            bounds=r.bound_flags,failed_candidates=r.failed_candidates,evaluations=r.function_evaluations,closure_margin=r.minimum_closure_margin,
            history_fingerprint=r.history_fingerprint,warnings=r.warnings,transfer_rates=[s.connections[0].rate for s in r.final_simulation.states]) for name,r in matches.items()})


if __name__=="__main__":
    output=evidence()
    (ROOT/"docs/phase_5b_numerical_results.json").write_text(json.dumps(output,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in output.items() if k not in ("matches","protected_sources")},indent=2))
    print(json.dumps({n:dict(parameters=r["parameters"],rmse=r["final_network_metrics_pa"]["rmse"]) for n,r in output["matches"].items()},indent=2))
