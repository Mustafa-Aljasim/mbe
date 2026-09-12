"""Forward-network numerical evidence plus exact protection of all prior phases."""
from dataclasses import asdict,replace
from pathlib import Path
from math import exp
import hashlib,json,runpy
from material_balance_studio.tank_network import simulate_network,NetworkSettings
from material_balance_studio.tank_network.diagnostics import maxima,tank_terms
from material_balance_studio.solver.simulation import simulate

ROOT=Path(__file__).resolve().parents[1]
network_case=runpy.run_path(str(ROOT/"examples/network_benchmarks.py"))["network_case"]


def evidence():
    expected=json.loads((ROOT/"docs/phase_4c_numerical_results.json").read_text())
    actual=json.loads(json.dumps(runpy.run_path(str(ROOT/"examples/run_phase4c_acceptance.py"))["evidence"]()))
    assert actual==expected,"Previous numerical outputs changed"
    hashes=json.loads((ROOT/"docs/phase_5a_protected_baseline.json").read_text())
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==v for p,v in hashes.items())
    cases={kind:simulate_network(network_case(kind)) for kind in ("closed","production","aquifer","injection","reversal")}
    cases["zero_T"]=simulate_network(network_case("production",0.))
    cases["large_T"]=simulate_network(network_case("closed",1e-8))
    fine=simulate_network(network_case(),NetworkSettings(max_step_days=.25))
    cases["closed_refined"]=fine
    for result in cases.values():
        assert result.converged
    continuous=(.004*30e6+.008*20e6)/.012+(.008/.012)*10e6*exp(-1e-9*(1/.004+1/.008)*100*86400)
    independent=simulate(cases["production"].network.tanks[0].reservoir,cases["production"].network.tanks[0].history)
    zero_difference=max(abs(s.tanks[j].pressure-r.pressure) for j,node in enumerate(cases["zero_T"].network.tanks)
        for s,r in zip(cases["zero_T"].states,simulate(node.reservoir,node.history).states))
    summaries={}
    for name,result in cases.items():
        summaries[name]=dict(qc=maxima(result),final_tanks={t.name:tank_terms(t) for t in result.states[-1].tanks},
            final_connections=[asdict(c) for c in result.states[-1].connections],warnings=result.warnings,
            pressure_history=[dict(time=str(s.time),pressures_pa={t.name:t.pressure for t in s.tanks},
                rates_m3_per_s=[e.rate for e in s.connections],cumulative_transfer_m3=[e.cumulative_transfer for e in s.connections]) for s in result.states])
    return dict(previous_phase4c_and_all_prior_outputs_identical=True,protected_sources=hashes,
        maximum_individual_tank_absolute_residual=max(maxima(r)["maximum_tank_absolute_residual"] for r in cases.values()),
        maximum_individual_tank_relative_residual=max(maxima(r)["maximum_tank_relative_residual"] for r in cases.values()),
        maximum_network_transfer_conservation_error=max(maxima(r)["maximum_transfer_conservation_error"] for r in cases.values()),
        zero_T_maximum_pressure_difference_pa=zero_difference,
        closed_continuous_analytical_final_A_pa=continuous,
        closed_refined_error_pa=abs(fine.states[-1].tanks[0].pressure-continuous),
        isolated_produced_tank_final_pressure_pa=independent.states[-1].pressure,
        benchmarks=summaries)


if __name__=="__main__":
    data=evidence()
    (ROOT/"docs/phase_5a_numerical_results.json").write_text(json.dumps(data,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in data.items() if k not in ("protected_sources","benchmarks")},indent=2))
