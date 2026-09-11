"""Equation inspector: values come directly from saved engine terms."""

from dataclasses import asdict
from datetime import date
from typing import Any

import pandas as pd

from material_balance_studio.domain.models import CUMULATIVE_FIELDS, SimulationResult, SimulationState
from material_balance_studio.aquifer.transient import diagnostics_from_variables, history_entries


def inspect_state(state: SimulationState) -> dict[str, Any]:
    """Flat auditable values; suffixes identify pressure/volume units."""
    balance = state.balance
    row = {
        "date": state.date.isoformat(), "step_index": state.step_index,
        "pressure_pa": state.pressure, "observed_pressure_pa": state.observed_pressure,
        "pressure_error_pa": state.pressure_error, "previous_pressure_pa": state.previous_pressure,
        "elapsed_days": state.elapsed_days,
        "production_withdrawal_m3": balance.withdrawal.production,
        "water_injection_support_m3": balance.withdrawal.water_injection,
        "gas_injection_support_m3": balance.withdrawal.gas_injection,
        "net_withdrawal_m3": balance.withdrawal.net,
        "eo": balance.expansion.eo, "eg": balance.expansion.eg,
        "meg": balance.expansion.meg, "efw": balance.expansion.efw, "et": balance.expansion.et,
        "oil_expansion_support_m3": balance.oil_expansion_support,
        "gas_cap_expansion_support_m3": balance.gas_cap_expansion_support,
        "rock_water_expansion_support_m3": balance.rock_water_expansion_support,
        "total_expansion_support_m3": balance.total_expansion_support,
        "aquifer_support_m3": balance.aquifer_support,
        "total_support_m3": balance.total_support,
        "residual_m3": balance.residual, "absolute_residual_m3": balance.absolute_residual,
        "relative_residual": balance.relative_residual,
        "converged": state.solver.converged, "iterations": state.solver.iterations,
        "function_calls": state.solver.function_calls,
        "initial_guess_pa": state.solver.initial_guess,
        "pressure_bounds_pa": state.solver.pressure_bounds,
        "bracket_pa": state.solver.bracket, "solver_message": state.solver.message,
        "warnings": " | ".join(state.warnings),
    }
    row.update({f"cumulative_{name}_m3": getattr(state.cumulative, name) for name in CUMULATIVE_FIELDS})
    row.update({f"increment_{name}_m3": getattr(state.increments, name) for name in CUMULATIVE_FIELDS})
    row.update({f"pvt_{key}": value for key, value in asdict(state.pvt).items()})
    aq, step = state.aquifer_state, state.aquifer_step
    row.update({
        "aquifer_model": aq.model_key if aq else "none",
        "incremental_aquifer_influx_m3": step.incremental_influx if step else 0.,
        "cumulative_aquifer_influx_m3": aq.cumulative_influx if aq else 0.,
        "aquifer_average_rate_m3_s": step.average_influx_rate if step else 0.,
        "aquifer_endpoint_rate_m3_s": step.endpoint_influx_rate if step else None,
        "previous_aquifer_pressure_pa": step.previous_aquifer_pressure if step else None,
        "aquifer_pressure_pa": aq.aquifer_pressure if aq else None,
        "aquifer_initial_pressure_pa": aq.initial_pressure if aq else None,
        "aquifer_driving_difference_pa": step.driving_pressure_difference if step else None,
        "aquifer_elapsed_seconds": aq.elapsed_time if aq else 0.,
        "aquifer_model_parameters": dict(aq.model_parameters) if aq else {},
        "aquifer_model_variables": dict(aq.model_variables) if aq else {},
    })
    diagnostics = diagnostics_from_variables(aq.model_variables) if aq else {}
    row.update({
        "aquifer_dimensionless_time": diagnostics.get("tD"),
        "aquifer_dimensionless_pressure": diagnostics.get("dimensionless_pressure"),
        "aquifer_dimensionless_pressure_derivative": diagnostics.get("dimensionless_pressure_derivative"),
        "aquifer_recurrence_term": diagnostics.get("recurrence_term"),
        "aquifer_response_value": diagnostics.get("response_value"),
        "aquifer_active_pressure_steps": diagnostics.get("active_pressure_steps"),
        "aquifer_current_step_contribution_m3": diagnostics.get("current_step_contribution"),
        "aquifer_historical_contribution_m3": diagnostics.get("historical_contribution"),
        "aquifer_pressure_step_count": len(history_entries(aq.model_variables)) if aq else 0,
    })
    return row


def inspect_timestep(result: SimulationResult, selected_date: date) -> dict[str, Any]:
    """Inspect a successfully solved history date; failed attempts remain separate."""
    for state in result.states:
        if state.date == selected_date:
            return inspect_state(state)
    raise KeyError(f"No successfully solved timestep on {selected_date}.")


def results_frame(result: SimulationResult) -> pd.DataFrame:
    """Export history states only (the initial state is separately accessible)."""
    return pd.DataFrame(inspect_state(state) for state in result.states)
