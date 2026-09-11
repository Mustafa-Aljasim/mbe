from dataclasses import replace
from datetime import date

import pytest

from conftest import analytic_oil_production
from material_balance_studio.diagnostics.balance_closure import inspect_timestep, results_frame
from material_balance_studio.domain.models import CumulativeVolumes, HistoryRecord
from material_balance_studio.solver.simulation import simulate


def benchmark_history():
    records = []
    for month, pressure in [(2, 28e6), (3, 25e6), (4, 20e6)]:
        np = analytic_oil_production(pressure)
        records.append(HistoryRecord(date(2020, month, 1), CumulativeVolumes(np=np, gp=80 * np)))
    return records


def test_d_chronological_state_carry_and_closure(tank):
    result = simulate(tank, benchmark_history())
    assert result.converged
    assert [state.pressure for state in result.states] == pytest.approx([28e6, 25e6, 20e6], abs=.01)
    previous = result.initial_state
    for state in result.states:
        assert state.previous_pressure == previous.pressure
        assert state.solver.initial_guess == previous.pressure
        assert state.step_index == previous.step_index + 1
        assert state.elapsed_days == (state.date - previous.date).days
        assert state.cumulative.np == pytest.approx(previous.cumulative.np + state.increments.np)
        assert state.balance.relative_residual <= 1e-8
        assert state.balance.absolute_residual <= 1e-5
        previous = state


def test_no_production_multiple_timesteps(tank):
    result = simulate(tank, [HistoryRecord(date(2020, month, 1)) for month in (1, 2, 3)])
    assert result.converged
    assert all(state.pressure == 30e6 and state.balance.residual == 0 for state in result.states)


def test_observations_are_optional_and_do_not_drive_solution(tank):
    history = benchmark_history()
    predicted = simulate(tank, history)
    observed = simulate(tank, [replace(row, observed_pressure=15e6) for row in history])
    assert [s.pressure for s in observed.states] == [s.pressure for s in predicted.states]
    assert predicted.states[0].pressure_error is None
    assert observed.states[0].pressure_error == pytest.approx(13e6)


def test_inspector_exposes_every_term_and_increment(tank):
    result = simulate(tank, benchmark_history())
    row = inspect_timestep(result, date(2020, 4, 1))
    assert row["eo"] == pytest.approx(.04)
    assert row["meg"] == row["efw"] == 0
    assert row["total_expansion_support_m3"] == pytest.approx(40000)
    assert row["net_withdrawal_m3"] == pytest.approx(40000)
    assert row["converged"] and row["increment_np_m3"] > 0
    assert len(results_frame(result)) == 3
    with pytest.raises(KeyError):
        inspect_timestep(result, date(2030, 1, 1))


def test_failure_stops_and_preserves_successful_states(tank):
    history = benchmark_history()[:1] + [
        HistoryRecord(date(2020, 3, 1), CumulativeVolumes(np=900000, gp=72e6)),
        HistoryRecord(date(2020, 4, 1), CumulativeVolumes(np=950000, gp=76e6)),
    ]
    result = simulate(tank, history)
    assert not result.converged
    assert len(result.states) == 1
    assert result.failed_timestep.date == date(2020, 3, 1)
    assert not result.failed_timestep.solution.diagnostics.converged


def test_injection_fvf_assumption_is_reported(tank):
    result = simulate(tank, [HistoryRecord(date(2020, 2, 1), CumulativeVolumes(winj=1000, ginj=10000))])
    assert result.converged
    assert len(result.warnings) == 2
    assert any("exceeds initial pressure" in message for message in result.states[0].warnings)


def test_splitting_history_does_not_change_cumulative_equation(tank):
    history = benchmark_history()
    multi = simulate(tank, history)
    single = simulate(tank, history[-1:])
    assert multi.states[-1].pressure == pytest.approx(single.states[-1].pressure, abs=.01)
