from dataclasses import asdict, replace
from datetime import date
from io import StringIO
from pathlib import Path

import pandas as pd
import pytest

from material_balance_studio.diagnostics.qc import balance_qc_status
from material_balance_studio.domain.models import CumulativeVolumes, HistoryRecord, ReservoirTank
from material_balance_studio.io.history import history_from_frame, pvt_from_frame, read_table
from material_balance_studio.io.reservoir import tank_from_inputs
from material_balance_studio.presentation.charts import closure_figure, pressure_figure
from material_balance_studio.presentation.tables import (
    RELATIVE_LABEL, advanced_residuals_frame, engineering_results_frame, equation_inspector_frame,
    equation_sides_frame, history_display_frame, pvt_display_frame, reservoir_display_frame,
)
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.solver.simulation import simulate
from material_balance_studio.units.conversions import from_si
from material_balance_studio.units.display import display_unit, from_display, to_display


@pytest.fixture
def example_case():
    examples = Path(__file__).resolve().parents[1] / "examples"
    pvt = pvt_from_frame(read_table(examples / "pvt.csv"))
    history = history_from_frame(read_table(examples / "history.csv"))
    tank = ReservoirTank(date(2020, 1, 1), 30e6, 1e6, .2, 5e-10, 4e-10, 0., TablePVTModel(pvt))
    return tank, pvt, history, simulate(tank, history)


@pytest.mark.parametrize("quantity,value", [
    ("pressure", 30e6), ("oil_volume", 1e6), ("water_volume", 10000),
    ("gas_volume", 80e6), ("bo", 1.2), ("bw", 1.01), ("bg", .004),
    ("rs", 80.), ("compressibility", 5e-10), ("reservoir_volume", 12345.),
    ("bwinj", 1.02), ("bginj", .005), ("expansion", .04),
])
@pytest.mark.parametrize("units", ["SI", "FIELD"])
def test_engineering_display_round_trips(quantity, value, units):
    shown = to_display(value, quantity, units)
    assert from_display(shown, quantity, units) == pytest.approx(value, rel=1e-14)


def test_pressure_reference_and_distinct_si_file_and_display_units():
    assert to_display(30e6, "pressure", "FIELD") == pytest.approx(4351.1321319, abs=1e-6)
    assert to_display(30e6, "pressure", "SI") == 30
    assert from_display(30, "pressure", "SI") == 30e6
    assert from_si(30e6, "pressure", "SI") == 30e6  # Existing SI file contract is unchanged.
    assert to_display(5e-10, "compressibility", "SI") == pytest.approx(.0005)


def test_pvt_field_display_regression(example_case):
    _, table, _, _ = example_case
    original = asdict(table)
    frame = pvt_display_frame(table, "FIELD")
    at_pi = frame.iloc[4]
    assert at_pi["Pressure (psia)"] == pytest.approx(4351.1321319, abs=1e-6)
    # Independent known SI/field factors: rs / 0.17810760668, bg / 5.61458333333.
    assert at_pi["Rs (scf/STB)"] == pytest.approx(449.1666666667)
    assert at_pi["Bg (rb/scf)"] == pytest.approx(.000712430426716)
    assert at_pi["Bo (rb/STB)"] == 1.2
    assert at_pi["Bw (rb/STB)"] == 1.0
    assert pvt_display_frame(table, "SI").iloc[4]["Pressure (MPa)"] == 30
    assert asdict(table) == original


def test_missing_observations_and_water_gas_history_labels():
    history = [HistoryRecord(date(2020, 1, 1), CumulativeVolumes(np=1, gp=1, wp=1, winj=1, ginj=1))]
    frame = history_display_frame(history, "FIELD")
    assert frame.iloc[0]["Cumulative Oil Production (STB)"] == pytest.approx(6.2898107704)
    assert frame.iloc[0]["Cumulative Water Production (bbl)"] == pytest.approx(6.2898107704)
    assert frame.iloc[0]["Cumulative Gas Production (scf)"] == pytest.approx(35.3146667215)
    assert pd.isna(frame.iloc[0]["Observed Pressure (psia)"])


@pytest.mark.parametrize("units,pressure_unit,volume_unit", [("SI", "MPa", "reservoir m³"), ("FIELD", "psia", "rb")])
def test_engineering_result_export_inspector_and_plot_agree(example_case, units, pressure_unit, volume_unit):
    tank, _, _, result = example_case
    frame = engineering_results_frame(result, units)
    assert f"Calculated Pressure ({pressure_unit})" in frame
    assert f"Production Withdrawal ({volume_unit})" in frame
    assert not any("_pa" in label or "_m3" in label for label in frame)
    assert set(frame["Balance QC"]) == {"PASS"}
    exported = pd.read_csv(StringIO(frame.to_csv(index=False)))
    pd.testing.assert_frame_equal(exported, frame, check_exact=False)
    fig = pressure_figure(result, units)
    assert fig.layout.yaxis.title.text == f"Average reservoir pressure ({pressure_unit})"
    assert list(fig.data[0].y) == frame[f"Calculated Pressure ({pressure_unit})"].tolist()
    assert list(fig.data[1].y) == frame[f"Observed Pressure ({pressure_unit})"].tolist()
    inspector = equation_inspector_frame(result.states[-1], units).set_index("Quantity")
    assert inspector.loc["Calculated Pressure", "Unit"] == pressure_unit
    assert inspector.loc["N × Eo", "Unit"] == volume_unit
    assert inspector.loc["Balance Error", "Unit"] == volume_unit
    assert inspector.loc["Solver Convergence", "Value"] == "Converged"
    assert {"Eo", "mEg", "Efw", "N × Eo", "N × mEg", "N × Efw", "Expansion Support"} <= set(inspector.index)
    sides = equation_sides_frame(result.states[-1], units)
    assert sides.iloc[0, 1] == pytest.approx(sides.iloc[1, 1], abs=1e-6)
    assert sides.iloc[2, 1] == to_display(result.states[-1].balance.residual, "reservoir_volume", units)
    assert f"Absolute Balance Residual ({volume_unit})" in advanced_residuals_frame(result, units)
    setup = reservoir_display_frame(tank, units).set_index("Quantity")
    assert setup.loc["Initial Pressure", "Unit"] == pressure_unit


def test_all_display_views_leave_canonical_state_unchanged(example_case):
    tank, pvt, history, result = example_case
    before = asdict(result)
    for units in ["FIELD", "SI", "FIELD", "SI"]:
        pvt_display_frame(pvt, units)
        history_display_frame(history, units)
        engineering_results_frame(result, units)
        reservoir_display_frame(tank, units)
        pressure_figure(result, units)
        for state in result.states:
            equation_inspector_frame(state, units)
            equation_sides_frame(state, units)
    assert asdict(result) == before


@pytest.mark.parametrize("error,status", [(0., "PASS"), (1e-14, "PASS"), (1e-8, "PASS"),
                                         (1.0001e-8, "WARNING"), (1e-5, "WARNING"), (1.0001e-5, "FAIL"),
                                         (float("nan"), "FAIL"), (float("inf"), "FAIL"), (None, "FAIL")])
def test_qc_thresholds(error, status):
    assert balance_qc_status(error) == status


def test_failed_solver_never_receives_pass():
    assert balance_qc_status(0., converged=False) == "FAIL"


def test_closure_plot_uses_tolerance_scale_without_altering_residuals(example_case):
    result = example_case[-1]
    chart = closure_figure(result)
    assert chart.layout.yaxis.range[0] == 0
    assert chart.layout.yaxis.range[1] >= 1e-8
    assert list(chart.data[0].y) == [state.balance.relative_residual for state in result.states]
    assert chart.layout.yaxis.title.text == RELATIVE_LABEL
    assert any(shape.y0 == 1e-8 for shape in chart.layout.shapes)
    assert chart.data[0].y[0] == 0  # No log-scale epsilon substitution.


def test_multistep_field_and_si_inputs_produce_same_canonical_results(example_case):
    si_tank, table, history, si_result = example_case
    field_pvt = pd.DataFrame({key: from_si(getattr(row, key), key, "FIELD")
                             for key in ("pressure", "bo", "rs", "bg", "bw")} for row in table.rows)
    field_history = pd.DataFrame({
        "date": row.date.isoformat(),
        **{name: from_si(getattr(row.cumulative, name), "gas_volume" if name in ("gp", "ginj") else "liquid_volume", "FIELD")
           for name in ("np", "gp", "wp", "winj", "ginj")},
        "observed_pressure": None if row.observed_pressure is None else from_si(row.observed_pressure, "pressure", "FIELD"),
    } for row in history)
    field_tank = tank_from_inputs(
        initial_date=si_tank.initial_date, initial_pressure=from_si(si_tank.initial_pressure, "pressure", "FIELD"),
        oil_in_place=from_si(si_tank.oil_in_place, "liquid_volume", "FIELD"), swc=si_tank.swc,
        cf=from_si(si_tank.cf, "compressibility", "FIELD"), cw=from_si(si_tank.cw, "compressibility", "FIELD"),
        m=si_tank.m, units="FIELD", pvt_model=TablePVTModel(pvt_from_frame(field_pvt, "FIELD")))
    converted_history = history_from_frame(field_history, "FIELD")
    field_result = simulate(field_tank, converted_history)
    assert si_result.converged and field_result.converged
    for si_state, field_state in zip(si_result.states, field_result.states, strict=True):
        assert field_state.pressure == pytest.approx(si_state.pressure, rel=0, abs=.01)
        assert asdict(field_state.cumulative) == pytest.approx(asdict(si_state.cumulative), rel=1e-12)
        assert asdict(field_state.increments) == pytest.approx(asdict(si_state.increments), rel=1e-12)
        assert field_state.pvt.bo == pytest.approx(si_state.pvt.bo, rel=1e-12)
        assert field_state.pvt.rs == pytest.approx(si_state.pvt.rs, rel=1e-12)
        assert field_state.pvt.bg == pytest.approx(si_state.pvt.bg, rel=1e-12)
        assert field_state.balance.total_expansion_support == pytest.approx(si_state.balance.total_expansion_support, abs=1e-6)
        assert field_state.balance.withdrawal.net == pytest.approx(si_state.balance.withdrawal.net, abs=1e-6)
        assert field_state.balance.residual == pytest.approx(si_state.balance.residual, abs=1e-6)
        assert field_state.balance.relative_residual <= 1e-8
