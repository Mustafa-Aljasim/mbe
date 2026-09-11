"""Independent convolution benchmarks and comparison acceptance."""
from dataclasses import asdict, replace, fields
from datetime import date, timedelta
from math import pi, sqrt
from pathlib import Path
import json
import runpy

import pytest
from scipy.integrate import quad
from streamlit.testing.v1 import AppTest

from material_balance_studio.aquifer import ModifiedVanEverdingenHurstAquifer, VanEverdingenHurstAquifer, PotAquifer, NoAquifer
from material_balance_studio.aquifer.modified_veh import integrated_response
from material_balance_studio.aquifer.veh_response import veh_infinite_water_influx_dimensionless as response, VEH_INFINITE_TABLE
from material_balance_studio.diagnostics.aquifer_comparison import compare_aquifers, pressure_metrics, engineering_qc, support_fraction, sensitivity_runs, ComparisonRun
from material_balance_studio.presentation.aquifer_comparison import comparison_frame, comparison_figures
from material_balance_studio.presentation.tables import equation_inspector_frame
from material_balance_studio.domain.models import ReservoirTank, HistoryRecord, CumulativeVolumes
from material_balance_studio.io.history import pvt_from_frame, read_table
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.solver.simulation import simulate
from material_balance_studio.units.display import to_display, from_display, column_label
from material_balance_studio.presentation.aquifer import FIELDS
from material_balance_studio.aquifer.registry import MODELS

ROOT = Path(__file__).resolve().parents[1]
DAY = 86400.


def model(cls=ModifiedVanEverdingenHurstAquifer):
    return cls(1000, 10, 10, .2, .2*.001*1e-9*1000**2/DAY, .001, 1e-9)


def fixture():
    tank = ReservoirTank(date(2020,1,1), 30e6, 1e6, .2, 5e-10, 4e-10, 0,
                         TablePVTModel(pvt_from_frame(read_table(ROOT/"examples/pvt.csv"))))
    history = tuple(HistoryRecord(date(2020,m,1), CumulativeVolumes(np=n,gp=n*80), p)
                    for m,n,p in ((2,10000,29e6),(3,20000,28e6),(4,30000,27e6)))
    return tank, history


def test_modified_reference_multistep_early_time():
    # C2.9 convolution with classical R(u)=2 sqrt(u/pi), all ages < .01.
    # Independent closed-form integral: 4/(3 sqrt(pi))*(b^1.5-a^1.5).
    aq = model()
    state = aq.initial_state(30e6)
    segments, previous_td, previous_we = [], 0., 0.
    for td, pressure in ((.001,29e6),(.003,27e6),(.006,26e6)):
        segments.append((previous_td,td,state.previous_reservoir_pressure-pressure))
        expected = 2*pi*.2*1e-9*10*1000**2*sum(
            dp/(end-start)*4/(3*sqrt(pi))*((td-start)**1.5-(td-end)**1.5)
            for start,end,dp in segments)
        step = aq.compute_step(state,state.previous_reservoir_pressure,pressure,(td-previous_td)*DAY)
        assert step.cumulative_influx == pytest.approx(expected, rel=2e-13)
        assert step.incremental_influx == pytest.approx(expected-previous_we, rel=2e-13)
        assert dict(step.updated_state.model_variables)["diag_tD"] == pytest.approx(td)
        state, previous_td, previous_we = step.updated_state, td, expected


@pytest.mark.parametrize("bounds", [(0,.003),(.004,.03),(.01,1),(1,100),(90,200),(9000,12000),(20,20.000001)])
def test_integrated_shared_response_against_independent_quadrature(bounds):
    low, high = bounds
    points = [t for t,_ in VEH_INFINITE_TABLE if low<t<high]
    expected = quad(response, low, high, points=points, epsabs=1e-9, epsrel=1e-10, limit=200)[0]
    assert integrated_response(low,high) == pytest.approx(expected, rel=2e-8)


def test_modified_candidates_repeatable_and_rejected_history_not_committed():
    aq = model()
    state = aq.compute_step(aq.initial_state(30e6),30e6,29e6,DAY).updated_state
    snapshot = asdict(state)
    first = aq.compute_step(state,29e6,27e6,DAY)
    aq.compute_step(state,29e6,35e6,DAY)
    assert first == aq.compute_step(state,29e6,27e6,DAY)
    assert asdict(state) == snapshot
    assert dict(first.updated_state.model_variables)["segment_count"] == 2
    with pytest.raises(ValueError, match="changed"):
        replace(aq,permeability=aq.permeability*2).compute_step(state,29e6,27e6,DAY)
    with pytest.raises(ValueError, match="changed"):
        model(VanEverdingenHurstAquifer).compute_step(state,29e6,27e6,DAY)


def test_original_modified_refinement_and_large_step_difference():
    values = []
    for count in (1,100):
        totals = []
        for cls in (VanEverdingenHurstAquifer, ModifiedVanEverdingenHurstAquifer):
            aq = model(cls)
            state = aq.initial_state(30e6)
            for i in range(count):
                state = aq.compute_step(state,state.previous_reservoir_pressure,30e6-1e6*(i+1)/count,DAY/count).updated_state
            totals.append(state.cumulative_influx)
        values.append(abs(totals[0]-totals[1])/totals[1])
    assert values[1] < .02
    assert values[0] > .2


def test_field_si_equivalence_and_inspector():
    aq = model()
    field = {k: from_display(to_display(v,FIELDS[k][1],"FIELD"),FIELDS[k][1],"FIELD") for k,v in asdict(aq).items()}
    tank, history = fixture()
    result = simulate(replace(tank,aquifer=aq),history)
    converted = simulate(replace(tank,aquifer=ModifiedVanEverdingenHurstAquifer(**field)),history)
    assert result.converged and converted.converged
    for a,b in zip(result.states,converted.states):
        assert a.pressure == pytest.approx(b.pressure,rel=1e-12)
        assert a.balance.aquifer_support == pytest.approx(b.balance.aquifer_support,rel=1e-12)
    inspector = equation_inspector_frame(result.states[0],"FIELD")
    assert "Linear between accepted endpoints" in set(inspector.Value)
    assert "psia/day" in set(inspector.Unit)


def test_comparison_alone_equivalence_metrics_and_display():
    tank, history = fixture()
    aq = model()
    original = asdict(aq)
    runs = compare_aquifers(tank,history,{"Modified":aq,"None":NoAquifer()})
    assert runs[0].result == simulate(replace(tank,aquifer=aq),history)
    assert asdict(aq) == original and tank.aquifer is None
    errors = [s.pressure-s.observed_pressure for s in runs[0].result.states]
    metrics = pressure_metrics(runs[0].result)
    assert metrics["rmse"] == pytest.approx(sqrt(sum(e*e for e in errors)/3))
    assert metrics["mae"] == pytest.approx(sum(abs(e) for e in errors)/3)
    assert metrics["bias"] == pytest.approx(sum(errors)/3)
    assert metrics["maximum"] == max(abs(e) for e in errors)
    frame = comparison_frame(runs,"FIELD")
    assert frame.iloc[0][column_label("Pressure RMSE","pressure","FIELD")] == pytest.approx(to_display(metrics["rmse"],"pressure","FIELD"))
    assert frame.iloc[0]["QC status"] == "CAUTION"
    assert len(comparison_figures(runs,"SI")) == 7
    assert pressure_metrics(None)["rmse"] is None
    assert support_fraction(runs[0].result.initial_state) is None


def test_qc_fail_and_pass_and_missing_observations():
    assert engineering_qc(ComparisonRun("bad",NoAquifer(),None,"bad geometry"))[0] == "FAIL"
    tank, history = fixture()
    result = simulate(replace(tank,aquifer=NoAquifer()),history)
    exact = tuple(replace(s,observed_pressure=s.pressure) for s in result.states)
    assert engineering_qc(ComparisonRun("none",NoAquifer(),replace(result,states=exact)))[0] == "PASS"
    missing = tuple(replace(s,observed_pressure=None) for s in result.states)
    assert engineering_qc(ComparisonRun("none",NoAquifer(),replace(result,states=missing)))[0] == "CAUTION"


def test_sensitivity_preserves_baseline_and_rejects_inactive_parameter():
    tank, history = fixture()
    aq = model()
    before = asdict(aq)
    runs = sensitivity_runs(tank,history,aq,"permeability",[aq.permeability*.5,aq.permeability])
    assert len(runs) == 2 and asdict(aq) == before
    assert runs[1].result == simulate(replace(tank,aquifer=aq),history)
    with pytest.raises(ValueError):
        sensitivity_runs(tank,history,aq,"radius_ratio",[5,10])


def test_all_previous_phase3b_fixed_outputs_unchanged():
    expected = json.loads((ROOT/"docs/phase_3b_numerical_results.json").read_text())
    actual = json.loads(json.dumps(runpy.run_path(str(ROOT/"examples/run_aquifer_acceptance.py"))["evidence"]()))
    assert actual["support_comparison"] == expected["support_comparison"]
    assert actual["benchmarks"] == expected["benchmarks"]


def test_comparison_ui_run_units_and_modified_switch():
    app = AppTest.from_file(str(ROOT/"app.py"),default_timeout=60).run()
    app.radio(key="aquifer_type").set_value("Modified Van Everdingen-Hurst").run()
    app.button[0].click().run()
    assert not app.exception and not app.error
    assert app.session_state["simulation"][0].converged
    app.multiselect(key="comparison_models").set_value(list(MODELS)).run()
    app.button(key="run_comparison").click().run()
    assert not app.exception and not app.error
    saved = app.session_state["aquifer_comparison"]
    assert len(saved[0]) == 7
    app.selectbox(key="display_units").select("FIELD").run()
    assert not app.exception
    assert app.session_state["aquifer_comparison"] == saved
    assert app.number_input(key="cmp_Pot_capacity").value == pytest.approx(to_display(.005,"aquifer_capacity","FIELD"))
    app.radio(key="aquifer_type").set_value("None").run()
    app.button[0].click().run()
    assert app.session_state["simulation"][0].initial_state.aquifer_state.model_key == "none"


@pytest.mark.parametrize("name", list(MODELS))
def test_every_model_identical_alone_and_in_comparison(name):
    tank, history = fixture()
    cls = MODELS[name]
    aq = cls(**{f.name:FIELDS[f.name][2] for f in fields(cls)})
    assert compare_aquifers(tank,history,{name:aq})[0].result == simulate(replace(tank,aquifer=aq),history)


def test_modified_zero_decline_rebound_and_solver_failure_state():
    from material_balance_studio.domain.models import SolverSettings
    aq = model()
    initial = aq.initial_state(30e6)
    step = aq.compute_step(initial,30e6,30e6,DAY)
    assert step.cumulative_influx == 0
    rebound = aq.compute_step(initial,30e6,32e6,DAY)
    assert rebound.cumulative_influx < 0 and rebound.warnings
    tank, history = fixture()
    result = simulate(replace(tank,aquifer=aq),history,SolverSettings(max_iterations=1))
    assert not result.converged
    assert result.initial_state.aquifer_state == initial
    assert not result.states
    assert engineering_qc(ComparisonRun("modified",aq,result))[0] == "FAIL"
