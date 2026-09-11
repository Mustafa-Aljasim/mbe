"""Independent analytical diagnostics and forward-physics regression protection."""
from dataclasses import replace, asdict
from datetime import date,timedelta
from math import sqrt
from pathlib import Path
from types import SimpleNamespace
import hashlib,json,runpy
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest

from material_balance_studio.domain.models import PVTProperties,PVTTable,ReservoirTank,HistoryRecord,CumulativeVolumes,WithdrawalTerms
from material_balance_studio.mbe.withdrawal import withdrawal_terms
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.solver.simulation import simulate
from material_balance_studio.aquifer import PotAquifer
from material_balance_studio.diagnostics.voidage import voidage
from material_balance_studio.diagnostics.vrr import vrr
from material_balance_studio.diagnostics.drive_indices import drive_indices
from material_balance_studio.diagnostics.pressure_qc import history_qc,diagnostic_history_from_frame,PressureObservation
from material_balance_studio.diagnostics.pressure_match import pressure_match
from material_balance_studio.diagnostics.havlena_odeh import observed_components,havlena_odeh,linear_diagnostic
from material_balance_studio.diagnostics.campbell_cole import campbell
from material_balance_studio.diagnostics.diagnosis import diagnose
from material_balance_studio.presentation.diagnosis import diagnosis_frame,diagnostic_inspector
from material_balance_studio.units.display import to_display,column_label
from material_balance_studio.units.conversions import from_si

ROOT=Path(__file__).resolve().parents[1]


def synthetic(capacity=0.,gas_cap=0.,injection=0.):
    rows=tuple(PVTProperties(p,1.2+4e-9*(30e6-p),80,.004+1e-10*(30e6-p),1)
               for p in (10e6,15e6,20e6,25e6,30e6,35e6))
    tank=ReservoirTank(date(2020,1,1),30e6,1e6,.2,0,0,gas_cap,TablePVTModel(PVTTable(rows)),aquifer=PotAquifer(capacity) if capacity else None)
    records=[]
    for i,p in enumerate((28e6,26e6,24e6,22e6,20e6),1):
        dp=30e6-p
        # Independent algebra: Fprod = N(Eo+mEg) + C dp + water injection.
        eo=4e-9*dp
        eg=1.2*(1e-10*dp)/.004
        winj=injection*i
        np=(1e6*(eo+gas_cap*eg)+capacity*dp+winj)/(1.2+eo)
        records.append(HistoryRecord(tank.initial_date+timedelta(days=30*i),CumulativeVolumes(np=np,gp=80*np,winj=winj),p))
    return tank,tuple(records)


@pytest.mark.parametrize("winj,ginj,expected",[(0,0,0),(12,0,1),(0,600,1),(6,300,1),(3,150,.5)])
def test_vrr_independent_water_gas_combined(winj,ginj,expected):
    pvt=PVTProperties(20e6,1.2,80,.01,1,bwinj=1,bginj=.02)
    increments=CumulativeVolumes(np=10,gp=800,winj=winj,ginj=ginj)
    cumulative=CumulativeVolumes(np=20,gp=1600,winj=winj*3,ginj=ginj*3)
    state=SimpleNamespace(pvt=pvt,increments=increments,cumulative=cumulative,balance=SimpleNamespace(withdrawal=withdrawal_terms(cumulative,pvt)))
    values=vrr(state)
    assert values["interval"]["injected"]==pytest.approx(expected)
    assert values["cumulative"]["injected"]==pytest.approx(1.5*expected)
    assert values["interval"]["water_injected"]==pytest.approx(winj/12)
    assert values["interval"]["gas_injected"]==pytest.approx(ginj*.02/12)


def test_voidage_components_and_zero_production():
    pvt=PVTProperties(20e6,1.2,80,.01,1,bwinj=1.1,bginj=.02)
    c=CumulativeVolumes(np=10,gp=1000,wp=3,winj=5,ginj=100)
    state=SimpleNamespace(pvt=pvt,increments=c,cumulative=c,balance=SimpleNamespace(withdrawal=withdrawal_terms(c,pvt)))
    d=voidage(state)["cumulative"]
    assert d==pytest.approx(dict(oil=12,free_gas=2,water=3,produced=17,water_injected=5.5,gas_injected=2,injected=7.5,net=9.5))
    state.increments=CumulativeVolumes(winj=10)
    assert vrr(state)["interval"]["injected"] is None
    assert vrr(state)["warnings"]


@pytest.mark.parametrize("supports",[(100,0,0,0,0),(20,80,0,0,0),(10,0,0,90,0),(5,0,0,0,95),(20,10,5,40,25)])
def test_drive_known_supports_close_without_normalization(supports):
    oil,gas,rock,aq,inj=supports
    b=SimpleNamespace(oil_expansion_support=oil,gas_cap_expansion_support=gas,rock_water_expansion_support=rock,aquifer_support=aq,withdrawal=WithdrawalTerms(100,inj*.6,inj*.4))
    d=drive_indices(SimpleNamespace(balance=b))
    assert list(d["values"].values())==pytest.approx([v/100 for v in supports])
    assert d["total"]==pytest.approx(1) and d["stackable"]
    b.withdrawal=WithdrawalTerms(200,inj*.6,inj*.4)
    assert drive_indices(SimpleNamespace(balance=b))["total"]==pytest.approx(.5)


@pytest.mark.parametrize("capacity,gas_cap,injection",[(0,0,0),(.005,0,0),(0,.5,0),(.005,.5,5000),(0,0,20000)])
def test_observed_ho_known_ooip_and_drive_forward_closure(capacity,gas_cap,injection):
    tank,history=synthetic(capacity,gas_cap,injection)
    components=observed_components(tank,history)
    ho=havlena_odeh(components,tank)
    assert ho["fit"]["slope"]==pytest.approx(1e6,rel=1e-10)
    assert ho["fit"]["intercept"]==pytest.approx(0,abs=1e-8)
    assert ho["fit"]["r2"]==pytest.approx(1)
    result=simulate(tank,history)
    assert result.converged
    assert [s.pressure for s in result.states]==pytest.approx([r.observed_pressure for r in history],abs=1e-4)
    for state in result.states:
        assert drive_indices(state)["total"]==pytest.approx(1,abs=1e-10)
    if capacity:
        assert components[-1]["We"]==pytest.approx(capacity*10e6)
    assert all(p["adjusted"]==pytest.approx(1e6) for p in campbell(components)["points"])


def test_gas_cap_transformation_independent_axes():
    tank,_=synthetic(gas_cap=.5)
    # Deliberately non-collinear Eg/Eo values to identify N intercept and Nm slope.
    rows=[dict(date=str(i),pressure=30e6-i*1e6,Eo=.01*i,Efw=.001*i,Eg=.003*i*i,
               Et=.011*i+.5*.003*i*i,adjusted=1e6*(.011*i+.5*.003*i*i)) for i in range(1,6)]
    fit=havlena_odeh(rows,tank,"gas_cap")["fit"]
    assert fit["intercept"]==pytest.approx(1e6)
    assert fit["slope"]==pytest.approx(5e5)


def test_ho_curvature_nonphysical_and_small_spread():
    assert linear_diagnostic([1,2],[2,4])["slope"] is None
    assert linear_diagnostic([1,1,1],[2,3,4])["slope"] is None
    curved=linear_diagnostic([1,2,3,4,5],[1,4,9,16,25])
    assert curved["status"]=="CAUTION"
    assert any("curvature" in n for n in curved["notes"])
    assert any("Nonphysical" in n for n in linear_diagnostic([1,2,3],[-1,-2,-3])["notes"])


def test_missing_observation_does_not_reuse_calculated_pressure_for_ho():
    tank,history=synthetic(capacity=.005)
    history=(history[0],replace(history[1],observed_pressure=None),*history[2:])
    rows=observed_components(tank,history)
    assert rows[1]["F"] is None and rows[2]["We"] is None
    assert rows[2]["adjusted"] is None and rows[2]["F"] is not None
    assert "missing" in rows[2]["note"]


def test_pressure_match_exact_metrics_sigma_and_missing():
    states=tuple(SimpleNamespace(date=date(2020,1,i),pressure_error=e,pressure=10e6+e,observed_pressure=10e6) for i,e in enumerate((1e6,-2e6,3e6),1))
    result=SimpleNamespace(states=states)
    metadata=(PressureObservation("2020-01-01",sigma=.5e6),PressureObservation("2020-01-02",sigma=0))
    match=pressure_match(result,metadata)
    assert match["metrics"]==pytest.approx(dict(count=3,rmse=sqrt(14/3)*1e6,mae=2e6,bias=2e6/3,maximum=3e6))
    assert [r["normalized"] for r in match["rows"]]==[2,None,None]


def test_raw_history_qc_preserves_flagged_rows():
    raw=pd.DataFrame(dict(date=["2020-01-01","2020-01-01","2019-01-01","2022-01-01"],np=[10,9,20,30],gp=[0]*4,wp=[0]*4,winj=[2,1,3,4],observed_pressure=[30e6,-1,None,40e6]))
    before=raw.copy(deep=True)
    flags=history_qc(raw,30e6,date(2020,1,1))
    assert {"Chronology","Duplicate date","Pressure","Pressure above Pi","Survey gap","Abrupt change","np","winj"}<=set(flags.Rule)
    assert "FAIL" in set(flags.Status)
    pd.testing.assert_frame_equal(raw,before)


@pytest.mark.parametrize("units",["SI","FIELD"])
def test_metadata_and_unit_equivalence(units):
    tank,history=synthetic(injection=3000)
    raw=pd.DataFrame([dict(date=str(r.date),np=from_si(r.cumulative.np,"liquid_volume",units),gp=from_si(r.cumulative.gp,"gas_volume",units),wp=0,winj=from_si(r.cumulative.winj,"liquid_volume",units),observed_pressure=from_si(r.observed_pressure,"pressure",units),pressure_sigma=from_si(1e5,"pressure",units),pressure_source="PBU-derived pressure",pressure_quality="CAUTION",pressure_note="Check representativeness") for r in history])
    parsed,metadata=diagnostic_history_from_frame(raw,units)
    assert metadata[0].sigma==pytest.approx(1e5)
    assert metadata[0].note=="Check representativeness"
    result=simulate(tank,parsed)
    data=diagnose(result,tank,parsed,metadata)
    assert data["ho"]["fit"]["slope"]==pytest.approx(1e6)
    original=simulate(tank,history)
    assert [vrr(s)["cumulative"]["injected"] for s in result.states]==pytest.approx([vrr(s)["cumulative"]["injected"] for s in original.states])
    frame=diagnosis_frame(data,units)
    assert frame.iloc[-1][column_label("cumulative produced","reservoir_volume",units)]==pytest.approx(to_display(result.states[-1].balance.withdrawal.production,"reservoir_volume",units))


def test_protected_sources_and_previous_aquifer_outputs_unchanged():
    hashes=json.loads((ROOT/"docs/phase_4a_forward_baseline.json").read_text())
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in hashes.items())
    previous=json.loads((ROOT/"docs/phase_3c_numerical_results.json").read_text())
    current=json.loads(json.dumps(runpy.run_path(str(ROOT/"examples/run_phase3c_acceptance.py"))["evidence"]()))
    assert current==previous


def test_reservoir_diagnosis_ui_saved_units_and_state():
    app=AppTest.from_file(str(ROOT/"app.py"),default_timeout=60).run()
    app.button[0].click().run()
    before=asdict(app.session_state["simulation"][0])
    app.button(key="build_diagnosis").click().run()
    assert not app.exception and not app.error
    assert any("Reservoir Diagnosis"==t.label for t in app.tabs)
    assert len(app.get("plotly_chart"))>=12
    app.selectbox(key="display_units").select("FIELD").run()
    assert not app.exception and not app.error
    assert asdict(app.session_state["simulation"][0])==before
    assert any("cumulative produced (rb)" in d.value for d in app.dataframe)
