"""Expected terms are fixed hand-calculated fixture values, not engine outputs."""

import json
from datetime import date
from pathlib import Path

import pandas as pd
import pytest

from material_balance_studio.diagnostics.balance_closure import inspect_state
from material_balance_studio.domain.models import ReservoirTank
from material_balance_studio.io.history import history_from_frame, pvt_from_frame
from material_balance_studio.mbe.oil_balance import evaluate_balance
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.solver.simulation import simulate


@pytest.fixture
def hand_case():
    fixture = json.loads((Path(__file__).parent / "fixtures" / "volumetric_hand_case.json").read_text())
    inputs = fixture["tank"].copy()
    inputs["initial_date"] = date.fromisoformat(inputs["initial_date"])
    tank = ReservoirTank(**inputs, pvt_model=TablePVTModel(pvt_from_frame(pd.DataFrame(fixture["pvt"]))))
    return tank, history_from_frame(pd.DataFrame(fixture["history"])), fixture["expected"]


def test_hand_calculated_volumetric_closure(hand_case):
    tank, history, expected = hand_case
    balance = evaluate_balance(expected["pressure_pa"], tank, history[-1].cumulative)
    assert balance.withdrawal.production == pytest.approx(5222.0, abs=1e-10)
    assert balance.withdrawal.water_injection == pytest.approx(200.0)
    assert balance.withdrawal.gas_injection == pytest.approx(22.0)
    assert balance.withdrawal.net == pytest.approx(5000.0)
    assert balance.total_expansion_support == pytest.approx(5000.0)
    assert balance.absolute_residual <= 1e-9


def test_hand_benchmark_recovers_pressure_and_every_term(hand_case):
    tank, history, expected = hand_case
    result = simulate(tank, history)
    assert result.converged
    actual = inspect_state(result.states[-1])
    assert actual["pressure_pa"] == pytest.approx(expected["pressure_pa"], rel=0, abs=.01)
    for key, value in expected.items():
        if key != "pressure_pa":
            assert actual[key] == pytest.approx(value, rel=1e-12, abs=1e-9), key
    assert actual["relative_residual"] <= 1e-12
