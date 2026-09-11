"""Chronological simulation with an explicit initial state and stop-on-failure policy."""

from collections.abc import Sequence

from material_balance_studio.domain.models import (
    CumulativeVolumes, FailedTimestep, HistoryRecord, ReservoirTank,
    SimulationResult, SimulationState, SolverDiagnostics, SolverSettings,
)
from material_balance_studio.domain.validation import validate_history
from material_balance_studio.mbe.oil_balance import evaluate_balance
from .timestep import advance_timestep
from material_balance_studio.aquifer.none import NoAquifer


def simulate(tank: ReservoirTank, history: Sequence[HistoryRecord],
             settings: SolverSettings | None = None) -> SimulationResult:
    """Predict pressure; observations are QC only and never overwrite reservoir state."""
    settings = settings or SolverSettings()
    validate_history(history, tank.initial_date)
    zero = CumulativeVolumes()
    initial_pvt = tank.pvt_model.properties_at_pressure(tank.initial_pressure)
    initial_state = SimulationState(
        step_index=0, date=tank.initial_date, pressure=tank.initial_pressure,
        previous_pressure=None, elapsed_days=0, cumulative=zero, increments=zero,
        pvt=initial_pvt,
        balance=evaluate_balance(tank.initial_pressure, tank, zero, settings.normalization_floor),
        solver=SolverDiagnostics(True, 0, 1, tank.initial_pressure,
                                 tank.pvt_model.pressure_bounds,
                                 (tank.initial_pressure, tank.initial_pressure), 0.0, "Initial equilibrium."),
        aquifer_state=(tank.aquifer or NoAquifer()).initial_state(tank.initial_pressure),
    )
    warnings = []
    if any(row.cumulative.winj > 0 for row in history) and initial_pvt.bwinj is None:
        warnings.append("Bwinj omitted: water injection uses resident Bw at candidate tank pressure.")
    if any(row.cumulative.ginj > 0 for row in history) and initial_pvt.bginj is None:
        warnings.append("Bginj omitted: gas injection uses resident Bg at candidate tank pressure.")
    previous = initial_state
    states: list[SimulationState] = []
    for record in history:
        step = advance_timestep(tank, previous, record, settings)
        if isinstance(step, FailedTimestep):
            return SimulationResult(initial_state, tuple(states), step, tuple(warnings))
        states.append(step)
        previous = step
    return SimulationResult(initial_state, tuple(states), warnings=tuple(warnings))
