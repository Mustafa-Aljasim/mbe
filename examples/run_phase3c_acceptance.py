"""Reproducible Phase 3C benchmarks, regression hashes and comparison evidence."""
from dataclasses import asdict, replace
from datetime import date
from math import pi, sqrt
from pathlib import Path
import hashlib
import json
import runpy

from material_balance_studio.aquifer import ModifiedVanEverdingenHurstAquifer, VanEverdingenHurstAquifer
from material_balance_studio.aquifer.registry import MODELS
from material_balance_studio.presentation.aquifer import FIELDS
from material_balance_studio.diagnostics.aquifer_comparison import compare_aquifers, pressure_metrics, engineering_qc
from material_balance_studio.domain.models import ReservoirTank, HistoryRecord, CumulativeVolumes
from material_balance_studio.io.history import pvt_from_frame, read_table
from material_balance_studio.pvt.table_model import TablePVTModel

ROOT = Path(__file__).resolve().parents[1]
DAY = 86400.


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,allow_nan=False).encode()).hexdigest()


def evidence():
    baseline = json.loads((ROOT/"docs/phase_3b_numerical_results.json").read_text())
    current = json.loads(json.dumps(runpy.run_path(str(ROOT/"examples/run_aquifer_acceptance.py"))["evidence"]()))
    regression = {name: {"before_sha256":digest(rows), "after_sha256":digest(current["support_comparison"][name]),
                         "identical": rows == current["support_comparison"][name]}
                  for name,rows in baseline["support_comparison"].items()}
    assert all(row["identical"] for row in regression.values())
    aq = ModifiedVanEverdingenHurstAquifer(1000,10,10,.2,.2*.001*1e-9*1000**2/DAY,.001,1e-9)
    state, segments, previous_td, previous_we, benchmark = aq.initial_state(30e6), [], 0., 0., []
    for td, pressure in ((.001,29e6),(.003,27e6),(.006,26e6)):
        segments.append((previous_td,td,state.previous_reservoir_pressure-pressure))
        expected = 2*pi*.2*1e-9*10*1000**2*sum(dp/(end-start)*4/(3*sqrt(pi))*((td-start)**1.5-(td-end)**1.5)
                                                  for start,end,dp in segments)
        step = aq.compute_step(state,state.previous_reservoir_pressure,pressure,(td-previous_td)*DAY)
        benchmark.append(dict(tD=td, expected_cumulative=expected, actual_cumulative=step.cumulative_influx,
                              expected_increment=expected-previous_we, actual_increment=step.incremental_influx,
                              relative_error=abs(step.cumulative_influx-expected)/expected))
        state, previous_td, previous_we = step.updated_state, td, expected
    refinement = {}
    for count in (1,100):
        values = {}
        for cls in (VanEverdingenHurstAquifer,ModifiedVanEverdingenHurstAquifer):
            model = cls(**asdict(aq))
            state = model.initial_state(30e6)
            for i in range(count):
                state = model.compute_step(state,state.previous_reservoir_pressure,30e6-1e6*(i+1)/count,DAY/count).updated_state
            values[model.key] = state.cumulative_influx
        values["relative_difference"] = abs(values["van_everdingen_hurst"]-values["modified_van_everdingen_hurst"])/values["modified_van_everdingen_hurst"]
        refinement[str(count)] = values
    tank = ReservoirTank(date(2020,1,1),30e6,1e6,.2,5e-10,4e-10,0,TablePVTModel(pvt_from_frame(read_table(ROOT/"examples/pvt.csv"))))
    history = tuple(HistoryRecord(date(2020,m,1),CumulativeVolumes(np=n,gp=n*80),p)
                    for m,n,p in ((2,10000,29e6),(3,20000,28e6),(4,30000,27e6)))
    from dataclasses import fields
    models = {name: cls(**{f.name:FIELDS[f.name][2] for f in fields(cls)}) for name,cls in MODELS.items()}
    runs = compare_aquifers(tank,history,models)
    assert all(r.result and r.result.converged for r in runs)
    maximum = max(s.balance.relative_residual for r in runs for s in r.result.states)
    return dict(regression=regression, benchmark=benchmark, refinement=refinement,
                maximum_relative_residual=max(maximum,current["maximum_relative_residual"]),
                comparison={r.name:dict(metrics=pressure_metrics(r.result),qc=engineering_qc(r),
                    final_pressure=r.result.states[-1].pressure,final_we=r.result.states[-1].balance.aquifer_support) for r in runs})


if __name__ == "__main__":
    output = evidence()
    target = ROOT/"docs/phase_3c_numerical_results.json"
    target.write_text(json.dumps(output,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(output,indent=2))
