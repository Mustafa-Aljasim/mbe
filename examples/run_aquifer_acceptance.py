"""Reproduce Phase 3A benchmark evidence; canonical pressure Pa, volume m³, time s."""
from dataclasses import asdict, replace
from datetime import date
from pathlib import Path
from math import log
import hashlib
import json
from material_balance_studio.aquifer import NoAquifer, PotAquifer, SchilthuisAquifer, FetkovichAquifer
from material_balance_studio.domain.models import ReservoirTank, PVTProperties, PVTTable, HistoryRecord, CumulativeVolumes
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.solver.simulation import simulate
from material_balance_studio.diagnostics.balance_closure import inspect_state

DAY=86400.


def evidence():
    output = {"benchmarks": {}}
    fixtures = [
        ("none",NoAquifer(),[(28e6,DAY,0,0,None),(26e6,DAY,0,0,None)]),
        ("pot",PotAquifer(.0025),[(28e6,10*DAY,5000,5000,28e6),(25e6,10*DAY,7500,12500,25e6),(27e6,10*DAY,-5000,7500,27e6)]),
        ("schilthuis",SchilthuisAquifer(100/(DAY*1e6)),[(28e6,10*DAY,1000,1000,30e6),(26e6,20*DAY,6000,7000,30e6),(27e6,5*DAY,1750,8750,30e6)]),
        ("fetkovich",FetkovichAquifer(1e6,1e-9,.001*log(2)/DAY),[(26e6,DAY,1000,1000,29e6),(26e6,DAY,1500,2500,27.5e6),(26e6,DAY,750,3250,26.75e6),(26e6,DAY,375,3625,26.375e6)]),
    ]
    for name,model,steps in fixtures:
        state=model.initial_state(30e6)
        rows=[]
        for pressure,dt,delta,cumulative,pa in steps:
            result=model.compute_step(state,state.previous_reservoir_pressure,pressure,dt)
            rows.append({"expected_increment":delta,"expected_cumulative":cumulative,"expected_aquifer_pressure":pa,
                         "actual":asdict(result),"increment_error":abs(result.incremental_influx-delta),
                         "cumulative_error":abs(result.cumulative_influx-cumulative),
                         "pressure_error":abs(result.updated_state.aquifer_pressure-pa) if pa is not None else 0})
            state=result.updated_state
        output["benchmarks"][name]=rows
    table=TablePVTModel(PVTTable(tuple(PVTProperties(p,1.2+4e-9*(30e6-p),80,bg,1)
        for p,bg in ((10e6,.012),(20e6,.006),(30e6,.004),(40e6,.003)))))
    tank=ReservoirTank(date(2020,1,1),30e6,1e6,.2,0,0,0,table)
    records=[HistoryRecord(date(2020,m,1),CumulativeVolumes(np=n,gp=80*n)) for m,n in ((2,10000),(3,20000),(4,30000))]
    models={"none":NoAquifer(),"pot_weak":PotAquifer(.001),"pot_strong":PotAquifer(.01),
            "schilthuis_weak":SchilthuisAquifer(1e-10),"schilthuis_strong":SchilthuisAquifer(1e-9),
            "fetkovich_weak":FetkovichAquifer(1e6,1e-9,1e-10),"fetkovich_strong":FetkovichAquifer(1e7,1e-9,1e-9)}
    output["support_comparison"]={}
    for name,model in models.items():
        result=simulate(replace(tank,aquifer=model),records)
        assert result.converged
        output["support_comparison"][name]=[inspect_state(s) for s in result.states]
    output["maximum_relative_residual"]=max(s["relative_residual"] for rows in output["support_comparison"].values() for s in rows)
    root=Path(__file__).resolve().parents[1]
    baseline=json.loads((root/"docs/phase_3a_baseline.json").read_text())
    output["physics_hashes"]={path:{"sha256":hashlib.sha256((root/path).read_bytes()).hexdigest(),
        "unchanged":hashlib.sha256((root/path).read_bytes()).hexdigest()==digest} for path,digest in baseline["hashes"].items()}
    return output


if __name__=="__main__":
    path=Path(__file__).resolve().parents[1]/"docs/phase_3a_numerical_results.json"
    result=evidence()
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(path)
    print("Maximum relative residual:",result["maximum_relative_residual"])
    for name,rows in result["support_comparison"].items():
        print(name,"final pressure MPa",rows[-1]["pressure_pa"]/1e6,"We m3",rows[-1]["aquifer_support_m3"])
