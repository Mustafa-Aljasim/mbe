"""Readable engineering datasets and inspector views built from canonical objects."""

from collections.abc import Sequence
from dataclasses import asdict

import pandas as pd

from material_balance_studio.diagnostics.balance_closure import inspect_state, results_frame
from material_balance_studio.diagnostics.qc import balance_qc_status
from material_balance_studio.domain.models import HistoryRecord, PVTTable, ReservoirTank, SimulationResult, SimulationState
from material_balance_studio.units.conversions import UnitSystem
from material_balance_studio.units.display import column_label, display_unit, format_value, to_display

# Shared schemas ensure plots, tables and CSV exports use identical conversions.
PVT_COLUMNS = {
    "pressure": ("Pressure", "pressure"), "bo": ("Bo", "bo"),
    "rs": ("Rs", "rs"), "bg": ("Bg", "bg"), "bw": ("Bw", "bw"),
    "bwinj": ("Bwinj", "bwinj"), "bginj": ("Bginj", "bginj"),
    "oil_viscosity": ("Oil Viscosity", "viscosity"),
    "gas_viscosity": ("Gas Viscosity", "viscosity"), "z": ("Z Factor", "dimensionless"),
}
HISTORY_COLUMNS = {
    "np": ("Cumulative Oil Production", "oil_volume"),
    "gp": ("Cumulative Gas Production", "gas_volume"),
    "wp": ("Cumulative Water Production", "water_volume"),
    "winj": ("Cumulative Water Injection", "water_volume"),
    "ginj": ("Cumulative Gas Injection", "gas_volume"),
    "observed_pressure": ("Observed Pressure", "pressure"),
}
RESULT_COLUMNS = {
    "pressure_pa": ("Calculated Pressure", "pressure"),
    "observed_pressure_pa": ("Observed Pressure", "pressure"),
    "pressure_error_pa": ("Pressure Error", "pressure"),
    "production_withdrawal_m3": ("Production Withdrawal", "reservoir_volume"),
    "water_injection_support_m3": ("Water Injection Support", "reservoir_volume"),
    "gas_injection_support_m3": ("Gas Injection Support", "reservoir_volume"),
    "net_withdrawal_m3": ("Net Withdrawal", "reservoir_volume"),
    "total_expansion_support_m3": ("Expansion Support", "reservoir_volume"),
    "incremental_aquifer_influx_m3": ("Incremental Aquifer Influx", "reservoir_volume"),
    "cumulative_aquifer_influx_m3": ("Cumulative Aquifer Influx We", "reservoir_volume"),
    "aquifer_average_rate_m3_s": ("Average Aquifer Influx Rate", "aquifer_rate"),
    "aquifer_pressure_pa": ("Aquifer Pressure", "pressure"),
    "total_support_m3": ("Total Expansion and Aquifer Support", "reservoir_volume"),
    "residual_m3": ("Balance Error", "reservoir_volume"),
}
EXPANSION_COLUMNS = {
    "eo": ("Eo", "expansion"), "meg": ("mEg", "expansion"),
    "efw": ("Efw", "expansion"),
    "oil_expansion_support_m3": ("N × Eo", "reservoir_volume"),
    "gas_cap_expansion_support_m3": ("N × mEg", "reservoir_volume"),
    "rock_water_expansion_support_m3": ("N × Efw", "reservoir_volume"),
}
RELATIVE_LABEL = "Absolute Relative Balance Residual (dimensionless)"


def _converted_frame(frame: pd.DataFrame, schema: dict[str, tuple[str, str]],
                     units: UnitSystem | str) -> pd.DataFrame:
    converted = pd.DataFrame(index=frame.index)
    for key, (label, quantity) in schema.items():
        if key in frame:
            converted[column_label(label, quantity, units)] = frame[key].map(
                lambda value: to_display(value, quantity, units)).astype(float)
    return converted


def pvt_display_frame(table: PVTTable, units: UnitSystem | str) -> pd.DataFrame:
    """Convert measured canonical rows, dropping wholly absent optional properties."""
    canonical = pd.DataFrame(asdict(row) for row in table.rows).dropna(axis=1, how="all")
    return _converted_frame(canonical, PVT_COLUMNS, units)


def history_display_frame(history: Sequence[HistoryRecord], units: UnitSystem | str) -> pd.DataFrame:
    canonical = pd.DataFrame({**asdict(row.cumulative), "observed_pressure": row.observed_pressure}
                             for row in history)
    frame = _converted_frame(canonical, HISTORY_COLUMNS, units)
    frame.insert(0, "Date", [row.date.isoformat() for row in history])
    return frame


def engineering_results_frame(result: SimulationResult, units: UnitSystem | str) -> pd.DataFrame:
    """Primary results/export view. Raw variable names stay in advanced diagnostics."""
    raw = results_frame(result)
    frame = _converted_frame(raw, RESULT_COLUMNS, units)
    frame.insert(0, "Date", [row.date.isoformat() for row in result.states])
    frame[RELATIVE_LABEL] = [row.balance.relative_residual for row in result.states]
    frame["Balance QC"] = [balance_qc_status(row.balance.relative_residual,
                                           converged=row.solver.converged) for row in result.states]
    return frame


def equation_inspector_frame(state: SimulationState, units: UnitSystem | str) -> pd.DataFrame:
    """Engineering values in selected units; raw solver JSON is a separate API."""
    raw = inspect_state(state)
    rows = [{"Quantity": "Date", "Value": state.date.isoformat(), "Unit": ""}]
    for schema in (RESULT_COLUMNS, EXPANSION_COLUMNS):
        for key, (label, quantity) in schema.items():
            rows.append({"Quantity": label, "Value": format_value(to_display(raw[key], quantity, units)),
                         "Unit": display_unit(quantity, units)})
    rows.extend([
        {"Quantity": "Aquifer Model", "Value": raw["aquifer_model"], "Unit": ""},
        {"Quantity": "Previous Aquifer Pressure", "Value": format_value(to_display(raw["previous_aquifer_pressure_pa"], "pressure", units)), "Unit": display_unit("pressure", units)},
        {"Quantity": "Aquifer minus Average Reservoir Pressure", "Value": format_value(to_display(raw["aquifer_driving_difference_pa"], "pressure", units)), "Unit": display_unit("pressure", units)},
        {"Quantity": "Endpoint Aquifer Influx Rate", "Value": format_value(to_display(raw["aquifer_endpoint_rate_m3_s"], "aquifer_rate", units)), "Unit": display_unit("aquifer_rate", units)},
        {"Quantity": "Aquifer Dimensionless Time tD", "Value": format_value(raw["aquifer_dimensionless_time"]), "Unit": "dimensionless"},
        {"Quantity": "Aquifer Dimensionless Pressure PD", "Value": format_value(raw["aquifer_dimensionless_pressure"]), "Unit": "dimensionless"},
        {"Quantity": "Aquifer Dimensionless Pressure Derivative", "Value": format_value(raw["aquifer_dimensionless_pressure_derivative"]), "Unit": "dimensionless"},
        {"Quantity": "Aquifer Recurrence Term", "Value": format_value(raw["aquifer_recurrence_term"]), "Unit": "reservoir m³"},
        {"Quantity": "VEH Response WeD", "Value": format_value(raw["aquifer_response_value"]), "Unit": "dimensionless"},
        {"Quantity": "VEH Active Pressure Steps", "Value": format_value(raw["aquifer_active_pressure_steps"]), "Unit": "count"},
        {"Quantity": "VEH Current Pressure-Step Contribution", "Value": format_value(to_display(raw["aquifer_current_step_contribution_m3"], "reservoir_volume", units)), "Unit": display_unit("reservoir_volume", units)},
        {"Quantity": "VEH Historical Superposition Contribution", "Value": format_value(to_display(raw["aquifer_historical_contribution_m3"], "reservoir_volume", units)), "Unit": display_unit("reservoir_volume", units)},
        {"Quantity": "Absolute Relative Balance Residual", "Value": format_value(state.balance.relative_residual), "Unit": "dimensionless"},
        {"Quantity": "Balance QC", "Value": balance_qc_status(state.balance.relative_residual, converged=state.solver.converged), "Unit": ""},
        {"Quantity": "Solver Convergence", "Value": "Converged" if state.solver.converged else "Not converged", "Unit": ""},
        {"Quantity": "Solver Iterations", "Value": str(state.solver.iterations), "Unit": "count"},
    ])
    if raw["aquifer_model"] == "modified_van_everdingen_hurst" and state.aquifer_step is not None:
        variables = dict(state.aquifer_state.model_variables)
        rows.append({"Quantity": "Pressure treatment", "Value": "Linear between accepted endpoints", "Unit": ""})
        for key, label in (("start_pressure", "Reservoir pressure at start"),
                           ("end_pressure", "Reservoir pressure at end"),
                           ("average_pressure", "Timestep-average reservoir pressure")):
            rows.append({"Quantity": label, "Value": format_value(to_display(variables["diag_"+key], "pressure", units)),
                         "Unit": display_unit("pressure", units)})
        rows.append({"Quantity": "Reservoir pressure slope", "Value": format_value(to_display(variables["diag_pressure_slope"]*86400, "pressure", units)),
                     "Unit": display_unit("pressure", units)+"/day"})
        rows.append({"Quantity": "Integrated current dimensionless response", "Value": format_value(variables["diag_integrated_current_response"]), "Unit": "dimensionless"})
    return pd.DataFrame(rows)


def equation_sides_frame(state: SimulationState, units: UnitSystem | str) -> pd.DataFrame:
    """Saved left/right sides and signed difference, without reimplementing physics."""
    return pd.DataFrame({
        "Equation": ["LEFT SIDE · Net Withdrawal", "RIGHT SIDE · N × Et + We", "DIFFERENCE · Residual"],
        column_label("Value", "reservoir_volume", units): [to_display(value, "reservoir_volume", units)
            for value in (state.balance.withdrawal.net, state.balance.total_support, state.balance.residual)],
    })


def reservoir_display_frame(tank: ReservoirTank, units: UnitSystem | str) -> pd.DataFrame:
    """A saved setup snapshot uses the active display basis, not raw Pa."""
    rows = [{"Quantity": "Initial Date", "Value": str(tank.initial_date), "Unit": ""}]
    for label, value, quantity in (
        ("Initial Pressure", tank.initial_pressure, "pressure"),
        ("Initial Oil in Place", tank.oil_in_place, "oil_volume"),
        ("Connate Water Saturation Swc", tank.swc, "dimensionless"),
        ("Rock Compressibility", tank.cf, "compressibility"),
        ("Water Compressibility", tank.cw, "compressibility"),
        ("Gas-cap Ratio m", tank.m, "dimensionless"),
    ):
        rows.append({"Quantity": label, "Value": format_value(to_display(value, quantity, units)),
                     "Unit": display_unit(quantity, units)})
    return pd.DataFrame(rows)


def advanced_residuals_frame(result: SimulationResult, units: UnitSystem | str) -> pd.DataFrame:
    return pd.DataFrame({
        "Date": [row.date.isoformat() for row in result.states],
        column_label("Absolute Balance Residual", "reservoir_volume", units): [
            to_display(row.balance.absolute_residual, "reservoir_volume", units) for row in result.states],
        RELATIVE_LABEL: [row.balance.relative_residual for row in result.states],
    })
