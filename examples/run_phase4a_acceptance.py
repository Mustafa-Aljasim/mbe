"""Reproduce independent Phase 4A analytical cases and protected-source hashes."""
from dataclasses import asdict
from datetime import date,timedelta
from pathlib import Path
from math import sqrt
import hashlib,json,runpy
from material_balance_studio.domain.models import PVTProperties,PVTTable,ReservoirTank,HistoryRecord,CumulativeVolumes
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.aquifer import PotAquifer
from material_balance_studio.solver.simulation import simulate
from material_balance_studio.diagnostics.diagnosis import diagnose

ROOT=Path(__file__).resolve().parents[1]


def evidence():
    hashes=json.loads((ROOT/"docs/phase_4a_forward_baseline.json").read_text())
    protected={p:dict(before=h,after=hashlib.sha256((ROOT/p).read_bytes()).hexdigest()) for p,h in hashes.items()}
    assert all(v["before"]==v["after"] for v in protected.values())
    prior=json.loads((ROOT/"docs/phase_3c_numerical_results.json").read_text())
    current=json.loads(json.dumps(runpy.run_path(str(ROOT/"examples/run_phase3c_acceptance.py"))["evidence"]()))
    assert prior==current
    cases={}
    maximum=prior["maximum_relative_residual"]
    for name,capacity,m,inj in (("depletion",0,0,0),("aquifer",.005,0,0),("gas_cap",0,.5,0),("injection",0,0,20000),("mixed",.005,.5,5000)):
        pvt=TablePVTModel(PVTTable(tuple(PVTProperties(p,1.2+4e-9*(30e6-p),80,.004+1e-10*(30e6-p),1) for p in (10e6,15e6,20e6,25e6,30e6,35e6))))
        tank=ReservoirTank(date(2020,1,1),30e6,1e6,.2,0,0,m,pvt,aquifer=PotAquifer(capacity) if capacity else None)
        history=[]
        for i,p in enumerate((28e6,26e6,24e6,22e6,20e6),1):
            dp=30e6-p
            # Expected OOIP is input, not derived using the diagnostics under test.
            oil=(1e6*(4e-9*dp+m*1.2*1e-10*dp/.004)+capacity*dp+inj*i)/(1.2+4e-9*dp)
            history.append(HistoryRecord(tank.initial_date+timedelta(days=30*i),CumulativeVolumes(np=oil,gp=80*oil,winj=inj*i),p))
        result=simulate(tank,history)
        assert result.converged
        data=diagnose(result,tank,history)
        maximum=max(maximum,max(s.balance.relative_residual for s in result.states))
        cases[name]=dict(expected_ooip=1e6,ho_slope=data["ho"]["fit"]["slope"],ho_intercept=data["ho"]["fit"]["intercept"],
                         ho_r2=data["ho"]["fit"]["r2"],max_drive_closure_error=max(abs(r["drive"]["closure_error"]) for r in data["rows"]),
                         final_drive=data["rows"][-1]["drive"]["values"],final_vrr=data["rows"][-1]["vrr"],
                         maximum_pressure_error=max(abs(s.pressure_error) for s in result.states))
    return dict(protected_sources=protected,phase3c_outputs_identical=True,cases=cases,
                maximum_forward_relative_residual=maximum,
                pressure_error_benchmark=dict(errors_pa=[1e6,-2e6,3e6],rmse=sqrt(14/3)*1e6,mae=2e6,bias=2e6/3,maximum=3e6))


if __name__=="__main__":
    output=evidence()
    (ROOT/"docs/phase_4a_numerical_results.json").write_text(json.dumps(output,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in output.items() if k!="protected_sources"},indent=2))
