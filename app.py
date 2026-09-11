"""Phase 1.1 testing interface; engineering physics remain in the SI engine."""

from dataclasses import asdict, replace
from datetime import date
from pathlib import Path

import pandas as pd
import streamlit as st

from material_balance_studio.diagnostics.balance_closure import inspect_state, results_frame
from material_balance_studio.diagnostics.qc import balance_qc_status
from material_balance_studio.domain.models import FailedTimestep, HistoryRecord, PVTTable, SolverSettings
from material_balance_studio.domain.validation import EngineeringValidationError
from material_balance_studio.io.history import history_from_frame, pvt_from_frame, read_table
from material_balance_studio.io.reservoir import tank_from_inputs
from material_balance_studio.presentation.charts import closure_figure, pressure_figure
from material_balance_studio.presentation.tables import (
    advanced_residuals_frame, engineering_results_frame, equation_inspector_frame,
    equation_sides_frame, history_display_frame, pvt_display_frame, reservoir_display_frame,
)
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.presentation.pvt_workflow import correlation_workflow, refresh_pvt_units
from material_balance_studio.presentation.aquifer import aquifer_inputs, refresh_aquifer_units, aquifer_figures, aquifer_setup_frame
from material_balance_studio.presentation.aquifer_comparison import comparison_workflow, refresh_comparison_units
from material_balance_studio.presentation.diagnosis import diagnosis_workflow, diagnostic_inspector
from material_balance_studio.presentation.history_matching import history_matching_workflow, refresh_matching_units
from material_balance_studio.presentation.uncertainty import uncertainty_workflow
from material_balance_studio.diagnostics.pressure_qc import diagnostic_history_from_frame, history_qc
from material_balance_studio.solver.simulation import simulate
from material_balance_studio.units.display import column_label, display_unit, format_value, from_display, to_display
from material_balance_studio.units.conversions import UnitSystem

EXAMPLES = Path(__file__).resolve().parent / "examples"
SETUP_FIELDS = {
    "setup_pressure": ("pressure", 30e6), "setup_oil": ("oil_volume", 1e6),
    "setup_cf": ("compressibility", 5e-10), "setup_cw": ("compressibility", 4e-10),
}


def initialize_setup() -> None:
    """Keep an SI setup draft, independent of saved simulation states."""
    st.session_state.setdefault("setup_si", {key: default for key, (_, default) in SETUP_FIELDS.items()})
    units = st.session_state.get("display_units", "SI")
    for key, (quantity, _) in SETUP_FIELDS.items():
        st.session_state.setdefault(key, to_display(st.session_state["setup_si"][key], quantity, units))


def change_display_units() -> None:
    """Switch widget presentation from the SI draft; edited values are preserved."""
    units = st.session_state["display_units"]
    for key, (quantity, _) in SETUP_FIELDS.items():
        st.session_state[key] = to_display(st.session_state["setup_si"][key], quantity, units)
    refresh_pvt_units()
    refresh_aquifer_units()
    refresh_comparison_units()
    refresh_matching_units()


def save_setup_value(key: str) -> None:
    quantity = SETUP_FIELDS[key][0]
    st.session_state["setup_si"][key] = from_display(st.session_state[key], quantity, st.session_state["display_units"])


def setup_input(label: str, key: str, units: UnitSystem, number_format: str = "%.8g") -> float:
    return st.number_input(column_label(label, SETUP_FIELDS[key][0], units), key=key,
                           format=number_format, on_change=save_setup_value, args=(key,))


def engineering_table(frame: pd.DataFrame) -> None:
    """Keep small Bg and residual values legible without rounding exported data."""
    columns = {name: st.column_config.NumberColumn(format="%.8g")
               for name in frame.select_dtypes(include="number").columns}
    st.dataframe(frame, hide_index=True, width="stretch", column_config=columns)


def input_table(kind: str, units: UnitSystem, initial_date=None) -> tuple[PVTTable | tuple[HistoryRecord, ...] | None, str]:
    """Parse original file units once, then present canonical data in project units."""
    is_pvt = kind == "PVT"
    stem = kind.lower()
    upload = st.file_uploader(f"Upload {kind} CSV / Excel", type=["csv", "xlsx"], key=f"{stem}_upload")
    file_units = st.selectbox(f"Uploaded {kind} file unit system", ["SI", "FIELD"],
                              key=f"{stem}_file_units", disabled=upload is None)
    effective_units = file_units if upload is not None else "SI"
    if upload is None:
        st.caption(f"Bundled synthetic example: source SI (Pa); preview and engineering download: {units.value}. "
                   "The uploaded-file selector is inactive until a file is supplied.")
        source = EXAMPLES / f"{stem}.csv"
    else:
        st.caption(f"Uploaded file: {upload.name} · source units: {effective_units} · displayed in: {units.value}.")
        upload.seek(0)
        source = upload
    if not is_pvt:
        st.session_state.pop("history_qc", None)
        st.session_state.pop("pressure_metadata", None)
    with st.expander(f"{kind} import format / original SI template"):
        st.caption("The uploaded-file unit system describes file contents. It does not control project display units.")
        if is_pvt:
            st.write("Required headers: pressure, bo, rs, bg, bw. Optional: bwinj, bginj, oil_viscosity, gas_viscosity, z.")
            st.caption("SI files: pressure in Pa; FVFs and Rs in m³/m³; viscosity in Pa·s. "
                       "FIELD: psia; Bo/Bw/Bwinj in rb/STB; Bg/Bginj in rb/scf; Rs in scf/STB; viscosity in cP. "
                       "Z is dimensionless. Pressure rows must increase. No extrapolation.")
        else:
            st.write("Headers: date, np, gp, wp; optional winj, ginj, observed_pressure. Dates: YYYY-MM-DD, sorted and unique.")
            st.write("Optional pressure metadata: pressure_source, pressure_sigma, pressure_quality, pressure_note. Sigma uses the uploaded pressure unit (Pa or psi). Sources: Average reservoir pressure, Static well pressure, PBU-derived pressure, RFT/MDT pressure, Estimated pressure, Other.")
            st.caption("SI files: surface m³ and Pa. FIELD: oil in STB, water in bbl, gas in scf, pressure in psia. "
                       "Cumulative volumes start at initial date. Blank observed pressure is allowed.")
        st.caption("Engineering downloads have readable unit-bearing headers. Use the original template for imports; "
                   "SI import pressure is Pa, whereas SI display pressure is MPa.")
        st.download_button(f"Download original {kind} import template (SI, Pa)",
                           (EXAMPLES / f"{stem}.csv").read_bytes(), f"{stem}_import_SI_Pa.csv", "text/csv")
    try:
        raw = read_table(source, upload.name if upload is not None else None)
        if is_pvt:
            canonical = pvt_from_frame(raw, effective_units)
        else:
            st.session_state["history_qc"] = history_qc(raw,st.session_state["setup_si"]["setup_pressure"],initial_date,effective_units)
            canonical, metadata = diagnostic_history_from_frame(raw,effective_units)
            st.session_state["pressure_metadata"] = metadata
        shown = pvt_display_frame(canonical, units) if is_pvt else history_display_frame(canonical, units)
        engineering_table(shown)
        if not is_pvt:
            with st.expander("History and pressure quality screening"):
                engineering_table(st.session_state["history_qc"])
                st.caption("Pressure metadata is retained with the history snapshot. Suspicious points are never deleted automatically.")
        st.download_button(f"Download {kind} engineering data (CSV, {units.value})", shown.to_csv(index=False),
                           f"{stem}_engineering_{units.value}.csv", "text/csv")
        return canonical, effective_units
    except (EngineeringValidationError, ValueError) as exc:
        st.error(str(exc))
        if not is_pvt and "history_qc" in st.session_state:
            engineering_table(st.session_state["history_qc"])
        return None, effective_units


def show_failed_step(failed: FailedTimestep, units: UnitSystem) -> None:
    solution = failed.solution
    diagnostics = solution.diagnostics
    st.error(f"FAIL · Simulation stopped at {failed.date}. Pressure solve did not satisfy the convergence requirements.")
    low, high = diagnostics.pressure_bounds
    st.write(f"Attempted pressure bounds: {format_value(to_display(low, 'pressure', units))} to "
             f"{format_value(to_display(high, 'pressure', units))} {display_unit('pressure', units)}.")
    if solution.balance is not None:
        st.write(f"Last absolute relative balance residual: {solution.balance.relative_residual:.3e}.")
    with st.expander("Advanced / Failed Solver Attempt (canonical SI)"):
        st.write(diagnostics.message)
        st.json(asdict(solution))


def show_results(units: UnitSystem) -> None:
    saved = st.session_state.get("simulation")
    if saved is None:
        st.info("Run the simulation to inspect pressure and balance closure.")
        return
    result, snapshot = saved
    st.caption(f"Last-run results · displayed in {units.value}. Changing display units converts the saved view. "
               "Run again after editing engineering inputs.")
    with st.expander("Reservoir inputs used for this run"):
        if "tank" in snapshot:
            engineering_table(reservoir_display_frame(snapshot["tank"], units))
            if snapshot["tank"].aquifer is not None:
                engineering_table(aquifer_setup_frame(snapshot["tank"].aquifer, units))
        else:
            st.caption("This saved run predates Phase 1.1. Run again to rebuild its setup summary.")
        st.caption(f"PVT source units: {snapshot['pvt_units']} · History source units: {snapshot['history_units']}.")
    for warning in result.warnings:
        st.warning(warning)
    if result.converged:
        st.success(f"Converged: {len(result.states)} history timesteps.")
    else:
        show_failed_step(result.failed_timestep, units)
    if not result.states:
        st.info("No history steps converged; the initial state remains available in advanced output.")
    else:
        frame = engineering_results_frame(result, units)
        maximum_error = max(state.balance.relative_residual for state in result.states)
        status = balance_qc_status(maximum_error, converged=result.converged)
        st.metric("Material-balance QC", status)
        st.metric("Maximum absolute relative balance residual", f"{maximum_error:.3e}")
        st.caption("PASS ≤ 1e-8 · WARNING > 1e-8 to 1e-5 · FAIL > 1e-5 or solver failure. "
                   "A failed timestep is reported separately and never treated as a successful state.")
        st.plotly_chart(pressure_figure(result, units), width="stretch", key="pressure_plot")
        st.plotly_chart(closure_figure(result), width="stretch", key="residual_plot")
        for index, figure in enumerate(aquifer_figures(result, units)):
            st.plotly_chart(figure, width="stretch", key=f"aquifer_plot_{index}")
        st.caption("Relative-residual chart includes the PASS tolerance and starts at zero. "
                   "Hover values retain scientific precision; floating-point noise is not magnified.")
        engineering_table(frame)
        st.download_button(f"Download engineering results (CSV, {units.value})", frame.to_csv(index=False),
                           f"simulation_engineering_{units.value}.csv", "text/csv")
        selected = st.selectbox("Equation inspector: date", [state.date for state in result.states])
        state = next(state for state in result.states if state.date == selected)
        st.latex(r"F_{prod}-W_{inj}B_{winj}-G_{inj}B_{ginj}=N(E_o+mE_g+E_{fw})+W_e")
        st.caption(f"Expansion factors: {display_unit('expansion', units)} of original oil. "
                   f"Withdrawal, injection support and expansion support: {display_unit('reservoir_volume', units)}.")
        for warning in state.warnings:
            st.warning(warning)
        engineering_table(equation_sides_frame(state, units))
        inspector = equation_inspector_frame(state, units)
        engineering_table(inspector)
        with st.expander("VRR / Drive-support decomposition"):
            engineering_table(diagnostic_inspector(state, units))
        st.download_button(f"Download equation inspector (CSV, {units.value})", inspector.to_csv(index=False),
                           f"equation_inspector_{selected}_{units.value}.csv", "text/csv")
        with st.expander("Advanced Diagnostics / Absolute Residuals"):
            engineering_table(advanced_residuals_frame(result, units))
        with st.expander("Advanced / Raw Solver State (canonical SI)"):
            st.json(inspect_state(state))
    with st.expander("Advanced / Raw Solver Output (canonical SI)"):
        st.caption("Raw pressure is Pa; volumes are canonical m³. This debug view does not use project display units.")
        raw = results_frame(result)
        st.dataframe(raw, hide_index=True, width="stretch")
        st.json(inspect_state(result.initial_state))
        st.download_button("Download raw solver output (CSV, canonical SI)", raw.to_csv(index=False),
                           "raw_solver_output_SI_Pa_m3.csv", "text/csv")


def main() -> None:
    st.set_page_config(page_title="Material Balance Studio", layout="wide")
    st.title("Material Balance Studio")
    st.caption("Phase 4C · Single black-oil tank · Matching and parameter identifiability")
    initialize_setup()
    units = UnitSystem(st.selectbox("Project / display unit system", ["SI", "FIELD"],
                                   key="display_units", on_change=change_display_units))
    st.caption(f"Setup, previews, results, charts and engineering downloads use {units.value} "
               f"({display_unit('pressure', units)} for pressure). Uploaded-file units are set separately.")
    setup, pvt_tab, history_tab, run_tab, results_tab, comparison_tab, diagnosis_tab, matching_tab, uncertainty_tab = st.tabs(
        ["Reservoir Setup", "PVT Data", "History", "Run Simulation", "Results / Equation Inspector", "Aquifer Comparison", "Reservoir Diagnosis", "History Matching", "Identifiability & Uncertainty"])
    with setup:
        initial_date = st.date_input("Initial date", date(2020, 1, 1))
        setup_input("Initial pressure", "setup_pressure", units)
        setup_input("Initial oil in place", "setup_oil", units)
        swc = st.number_input("Connate water saturation Swc (dimensionless)", value=0.2, format="%.4f")
        setup_input("Rock compressibility", "setup_cf", units, "%.6e")
        setup_input("Water compressibility", "setup_cw", units, "%.6e")
        m = st.number_input("Initial gas-cap / oil reservoir-volume ratio m (dimensionless)", value=0.0, format="%.4f",key="setup_m")
        st.caption("Initial reference is zero cumulative production/injection. Switching display units preserves edited setup values.")
        aquifer = aquifer_inputs(units, st.session_state["setup_si"]["setup_pressure"])
    with pvt_tab:
        pvt_source = st.radio("PVT Source", ["Tabulated PVT", "Correlation PVT", "Matched Correlation PVT"], horizontal=True, key="pvt_source")
        if pvt_source == "Tabulated PVT":
            pvt_table, pvt_units = input_table("PVT", units)
            pvt_model = TablePVTModel(pvt_table) if pvt_table is not None else None
        else:
            pvt_model = correlation_workflow(units, st.session_state["setup_si"]["setup_pressure"], pvt_source)
            pvt_units = f"{pvt_source} (canonical SI)"
    with history_tab:
        history, history_units = input_table("History", units, initial_date)
    with run_tab:
        settings = SolverSettings()
        st.write(f"Solve chronological pressure states; results are displayed in {units.value}.")
        st.caption(f"Solver acceptance requires BOTH |residual| ≤ "
                   f"{to_display(settings.absolute_residual_tolerance, 'reservoir_volume', units):.3e} "
                   f"{display_unit('reservoir_volume', units)} AND relative residual ≤ {settings.relative_residual_tolerance:.0e}. "
                   f"Normalization floor = {format_value(to_display(settings.normalization_floor, 'reservoir_volume', units))} "
                   f"{display_unit('reservoir_volume', units)}. WARNING classification does not relax solver acceptance.")
        if st.button("Run simulation", type="primary", disabled=pvt_model is None or history is None or aquifer is None):
            try:
                values = st.session_state["setup_si"]
                tank = tank_from_inputs(initial_date=initial_date, initial_pressure=values["setup_pressure"],
                                        oil_in_place=values["setup_oil"], swc=swc, cf=values["setup_cf"],
                                        cw=values["setup_cw"], m=m, pvt_model=pvt_model)
                tank = replace(tank, aquifer=aquifer)
                result = simulate(tank, history, settings)
                st.session_state["simulation"] = (result, {
                    "tank": tank, "pvt_units": pvt_units, "history_units": history_units,
                    "history": tuple(history), "history_qc": st.session_state["history_qc"].copy(),
                    "pressure_metadata": st.session_state["pressure_metadata"],
                })
                st.success("Run finished. Open Results / Equation Inspector.")
            except (EngineeringValidationError, ValueError) as exc:
                st.session_state.pop("simulation", None)
                if "outside PVT range" in str(exc):
                    st.error("Initial pressure is outside the measured PVT range shown on the PVT tab.")
                    with st.expander("Advanced / Input Validation (canonical SI)"):
                        st.write(str(exc))
                else:
                    st.error(str(exc))
    with results_tab:
        show_results(units)
    with comparison_tab:
        def comparison_tank():
            values = st.session_state["setup_si"]
            return tank_from_inputs(initial_date=initial_date, initial_pressure=values["setup_pressure"],
                                    oil_in_place=values["setup_oil"], swc=swc, cf=values["setup_cf"],
                                    cw=values["setup_cw"], m=m, pvt_model=pvt_model)
        comparison_workflow(comparison_tank, history, units, pvt_model is not None and history is not None)
    with diagnosis_tab:
        diagnosis_workflow(units)
    with matching_tab:
        history_matching_workflow(units)
    with uncertainty_tab:
        uncertainty_workflow(units)


if __name__ == "__main__":
    main()
