"""A timestep explicitly consumes the previous state and cumulative increments."""

from material_balance_studio.domain.models import (
    CUMULATIVE_FIELDS, CumulativeVolumes, FailedTimestep, HistoryRecord,
    ReservoirTank, SimulationState, SolverSettings,
)
from material_balance_studio.domain.validation import EngineeringValidationError
from .pressure_solver import solve_pressure
from material_balance_studio.aquifer.base import AquiferContext
from material_balance_studio.aquifer.none import NoAquifer
from dataclasses import replace


def cumulative_increments(previous: CumulativeVolumes, current: CumulativeVolumes) -> CumulativeVolumes:
    """Reject decreasing totals, including injection; do not clip small decreases."""
    changes = {}
    for name in CUMULATIVE_FIELDS:
        change = getattr(current, name) - getattr(previous, name)
        if change < 0:
            raise EngineeringValidationError(f"Cumulative {name} decreases from {getattr(previous, name)} "
                                             f"to {getattr(current, name)}.")
        changes[name] = change
    return CumulativeVolumes(**changes)


def advance_timestep(tank: ReservoirTank, previous: SimulationState, record: HistoryRecord,
                     settings: SolverSettings) -> SimulationState | FailedTimestep:
    """Update component totals from the prior state, solve, then commit a new state."""
    if record.date < previous.date or (record.date == previous.date and previous.step_index > 0):
        raise EngineeringValidationError("Timestep dates must advance chronologically.")
    increments = cumulative_increments(previous.cumulative, record.cumulative)
    if record.date == previous.date and any(getattr(increments, name) > 0 for name in CUMULATIVE_FIELDS):
        raise EngineeringValidationError("Nonzero production/injection cannot occur at initial_date.")
    # Reconstruct totals from state + increments. Current-pressure cumulative
    # voidage is recomputed, rather than summing historical reservoir volumes.
    cumulative = CumulativeVolumes(**{
        name: getattr(previous.cumulative, name) + getattr(increments, name)
        for name in CUMULATIVE_FIELDS
    })
    model = tank.aquifer or NoAquifer()
    if previous.aquifer_state is None and previous.step_index > 0 and model.key != "none":
        raise EngineeringValidationError("Previous aquifer state is missing; restart the simulation instead of resetting aquifer memory.")
    prior_aquifer = previous.aquifer_state or model.initial_state(tank.initial_pressure)
    dt = (record.date-previous.date).days*86400.
    # An explicitly recorded zero-volume initial-date row is not a time interval.
    # Preserve the legacy row while leaving aquifer time and state unchanged.
    context = AquiferContext(model, prior_aquifer, previous.pressure, dt) if dt > 0 else None
    solution = solve_pressure(tank if dt > 0 else replace(tank, aquifer=None), cumulative,
                              previous.pressure, settings, aquifer_context=context)
    if not solution.diagnostics.converged:
        return FailedTimestep(record.date, increments, solution)
    assert solution.pressure is not None and solution.balance is not None
    pvt = tank.pvt_model.properties_at_pressure(solution.pressure)
    warnings = list(solution.diagnostics.warnings)
    if solution.aquifer_step is not None:
        warnings.extend(solution.aquifer_step.warnings)
        if abs(solution.balance.aquifer_support) > 2*max(abs(solution.balance.withdrawal.net), settings.normalization_floor):
            warnings.append("Aquifer support exceeds twice the net withdrawal magnitude; review aquifer parameters, injection and pressure recovery.")
    if cumulative.gp < cumulative.np * pvt.rs:
        warnings.append("Cumulative produced GOR is below current Rs. The signed solution-gas "
                        "correction was retained; review gas history and PVT.")
    if cumulative.np > tank.oil_in_place:
        warnings.append("Cumulative oil production exceeds initial oil in place; review the case.")
    if solution.pressure > tank.initial_pressure:
        warnings.append("Injection-supported pressure exceeds initial pressure; confirm PVT coverage "
                        "and validity of the constant-compressibility approximation.")
    return SimulationState(
        step_index=previous.step_index + 1, date=record.date, pressure=solution.pressure,
        previous_pressure=previous.pressure, elapsed_days=(record.date - previous.date).days,
        cumulative=cumulative, increments=increments, pvt=pvt, balance=solution.balance,
        solver=solution.diagnostics, observed_pressure=record.observed_pressure,
        warnings=tuple(warnings),
        aquifer_state=solution.aquifer_step.updated_state if solution.aquifer_step else prior_aquifer,
        aquifer_step=solution.aquifer_step,
    )
