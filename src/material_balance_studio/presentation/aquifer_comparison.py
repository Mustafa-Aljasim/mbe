"""Unit-aware comparison tables, charts and isolated configuration controls."""
from dataclasses import fields
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from material_balance_studio.aquifer.registry import MODELS
from material_balance_studio.diagnostics.aquifer_comparison import (
    compare_aquifers, pressure_metrics, engineering_qc, support_fraction,
    sensitivity_runs, SENSITIVITY_PARAMETERS,
)
from material_balance_studio.units.display import to_display, from_display, column_label
from .aquifer import FIELDS, aquifer_setup_frame
from .tables import equation_inspector_frame


def comparison_frame(runs, units):
    rows = []
    for run in runs:
        metrics = pressure_metrics(run.result)
        status, notes = engineering_qc(run)
        states = run.result.states if run.result else ()
        final = states[-1] if states else None
        row = {"Model": run.name, "Observed points": metrics["count"],
               "Complete history": bool(run.result and run.result.converged), "QC status": status,
               "Engineering notes": " | ".join(notes) or "Screening checks passed; geological review remains necessary."}
        for key, label in (("rmse", "Pressure RMSE"), ("mae", "Pressure MAE"),
                           ("bias", "Mean pressure bias"), ("maximum", "Maximum absolute pressure error")):
            row[column_label(label, "pressure", units)] = to_display(metrics[key], "pressure", units)
        row[column_label("Final cumulative We", "reservoir_volume", units)] = to_display(final.balance.aquifer_support if final else None, "reservoir_volume", units)
        # Only finite lumped aquifers calculate a depleted average aquifer pressure.
        pa = final.aquifer_state.aquifer_pressure if final and run.model.key in ("pot", "fetkovich") else None
        row[column_label("Final aquifer pressure", "pressure", units)] = to_display(pa, "pressure", units)
        row["Maximum relative MBE residual"] = max((s.balance.relative_residual for s in states), default=None)
        row["Final aquifer support fraction"] = support_fraction(final) if final else None
        rows.append(row)
    return pd.DataFrame(rows)


SERIES = {
    "Calculated pressure": ("pressure", lambda s: s.pressure),
    "Cumulative We": ("reservoir_volume", lambda s: s.balance.aquifer_support),
    "Incremental We": ("reservoir_volume", lambda s: s.aquifer_step.incremental_influx if s.aquifer_step else 0.),
    "Average influx rate": ("aquifer_rate", lambda s: s.aquifer_step.average_influx_rate if s.aquifer_step else None),
    "Pressure error": ("pressure", lambda s: s.pressure_error),
    "Relative MBE residual": ("dimensionless", lambda s: s.balance.relative_residual),
    "Aquifer support fraction": ("dimensionless", support_fraction),
}


def comparison_figures(runs, units):
    figures = []
    for label, (quantity, getter) in SERIES.items():
        fig = go.Figure()
        for run in runs:
            if run.result:
                fig.add_scatter(x=[s.date for s in run.result.states],
                                y=[to_display(getter(s), quantity, units) for s in run.result.states],
                                name=run.name, mode="lines+markers")
        if label == "Calculated pressure":
            # Use the longest available common-history view even if one model failed early.
            result = max((r.result for r in runs if r.result), key=lambda r: len(r.states), default=None)
            if result:
                fig.add_scatter(x=[s.date for s in result.states], y=[to_display(s.observed_pressure, "pressure", units) for s in result.states],
                                name="Observed pressure", mode="markers", marker=dict(color="black", symbol="x", size=10))
        fig.update_layout(title=label, xaxis_title="Date", yaxis_title=column_label(label, quantity, units))
        figures.append(fig)
    return figures


def refresh_comparison_units():
    for name, draft in st.session_state.get("comparison_drafts", {}).items():
        for parameter, value in draft.items():
            st.session_state[f"cmp_{name}_{parameter}"] = to_display(value, FIELDS[parameter][1], st.session_state["display_units"])


def save_comparison_value(name, parameter):
    st.session_state["comparison_drafts"][name][parameter] = from_display(
        st.session_state[f"cmp_{name}_{parameter}"], FIELDS[parameter][1], st.session_state["display_units"])


def comparison_workflow(tank_factory, history, units, enabled):
    st.subheader("Aquifer Comparison")
    st.caption("All runs use the current reservoir, PVT and production/injection history. Observed pressure is a common diagnostic target; parameters are supplied by you.")
    selected = st.multiselect("Models to compare", list(MODELS), default=[], key="comparison_models")
    st.session_state.setdefault("comparison_drafts", {})
    models, invalid = {}, False
    for name in selected:
        with st.expander(name+" parameters"):
            if "Everdingen" in name or name == "Carter-Tracy":
                st.caption("Infinite-acting radial. Radius ratio has no effect on this response. Original VEH uses pressure steps; Modified VEH uses linear pressure segments.")
            draft = st.session_state["comparison_drafts"].setdefault(name, {f.name: FIELDS[f.name][2] for f in fields(MODELS[name])})
            for parameter, value in draft.items():
                label, quantity, _, help_text = FIELDS[parameter]
                key = f"cmp_{name}_{parameter}"
                st.session_state.setdefault(key, to_display(value, quantity, units))
                st.number_input(column_label(label, quantity, units), key=key, format="%.8g", help=help_text,
                                on_change=save_comparison_value, args=(name, parameter))
            try:
                models[name] = MODELS[name](**draft)
            except ValueError as exc:
                st.error(f"FAIL · {name}: {exc}")
                invalid = True
    if st.button("Run aquifer comparison", key="run_comparison", disabled=not enabled or not models or invalid):
        try:
            tank = tank_factory()
            runs = compare_aquifers(tank, history, models)
            st.session_state["aquifer_comparison"] = (runs, tank, tuple(history))
        except ValueError as exc:
            st.session_state.pop("aquifer_comparison", None)
            st.error(str(exc))
    saved = st.session_state.get("aquifer_comparison")
    if not saved:
        return
    runs, tank, saved_history = saved
    st.caption("Saved comparison: rerun after editing inputs. Unit switching converts the saved view.")
    st.write("NUMERICAL FIT")
    eligible = [r for r in runs if r.result and r.result.converged and pressure_metrics(r.result)["count"]]
    if eligible:
        best = min(eligible, key=lambda r: pressure_metrics(r.result)["rmse"])
        st.info(f"Best numerical pressure fit: {best.name}. Lowest RMSE does not establish physical correctness.")
    else:
        st.info("No complete run with observed pressures is available for numerical ranking.")
    summary = comparison_frame(runs, units)
    st.dataframe(summary, hide_index=True, width="stretch")
    st.download_button("Download comparison summary", summary.to_csv(index=False), f"aquifer_comparison_{units}.csv")
    st.write("ENGINEERING ACCEPTABILITY")
    st.caption("PASS means screening checks passed, not a verified geological model. CAUTION requires review. FAIL indicates invalid or incomplete results. Support fraction = We / net withdrawal; unavailable below 1 reservoir m³.")
    for run in runs:
        status, notes = engineering_qc(run)
        with st.expander(f"{run.name} · {status}", expanded=status != "PASS"):
            st.dataframe(aquifer_setup_frame(run.model, units), hide_index=True)
            for note in notes:
                (st.error if status == "FAIL" else st.warning)(note)
    for i, fig in enumerate(comparison_figures(runs, units)):
        st.plotly_chart(fig, key=f"comparison_plot_{i}", width="stretch")
    dates = sorted({s.date for r in runs if r.result for s in r.result.states})
    if dates:
        selected_date = st.selectbox("Comparison inspector: date", dates, key="comparison_date")
        frames = []
        for run in runs:
            state = next((s for s in run.result.states if s.date == selected_date), None) if run.result else None
            if state:
                frame = equation_inspector_frame(state, units)
                frame.insert(0, "Model", run.name)
                frames.append(frame)
            else:
                st.warning(f"{run.name}: no converged state at {selected_date}.")
        if frames:
            combined = pd.concat(frames, ignore_index=True)
            quantities = ["Calculated Pressure", "Observed Pressure", "Pressure Error", "Pressure treatment",
                          "Reservoir pressure at start", "Reservoir pressure at end", "Timestep-average reservoir pressure",
                          "Reservoir pressure slope", "Aquifer Dimensionless Time tD", "VEH Response WeD",
                          "Incremental Aquifer Influx", "Cumulative Aquifer Influx We", "Average Aquifer Influx Rate",
                          "Absolute Relative Balance Residual"]
            compact = combined[combined.Quantity.isin(quantities)].pivot(index=["Quantity", "Unit"], columns="Model", values="Value").reset_index()
            st.dataframe(compact, hide_index=True, width="stretch")
            with st.expander("Full comparison equation terms"):
                st.dataframe(combined, hide_index=True, width="stretch")
    available = {r.name: r.model for r in runs if r.model.key in SENSITIVITY_PARAMETERS}
    if available:
        with st.expander("One-parameter sensitivity preview"):
            name = st.selectbox("Sensitivity model", list(available), key="sensitivity_model")
            model = available[name]
            parameter = st.selectbox("Sensitivity parameter", SENSITIVITY_PARAMETERS[model.key], key="sensitivity_parameter")
            quantity = FIELDS[parameter][1]
            # Unit-specific key prevents stale values being reinterpreted after conversion.
            raw = st.text_input(column_label("Comma-separated values", quantity, units),
                                value=", ".join(f"{to_display(getattr(model, parameter)*f, quantity, units):.8g}" for f in (.5, 1)),
                                key=f"sensitivity_values_{name}_{parameter}_{units}")
            if st.button("Run sensitivity preview", key="run_sensitivity"):
                try:
                    values = [from_display(float(v.strip()), quantity, units) for v in raw.split(",")]
                    preview = sensitivity_runs(tank, saved_history, model, parameter, values)
                    st.session_state["sensitivity_preview"] = preview
                except ValueError as exc:
                    st.session_state.pop("sensitivity_preview", None)
                    st.error(str(exc))
            if st.session_state.get("sensitivity_preview"):
                preview = st.session_state["sensitivity_preview"]
                st.caption("Saved sensitivity runs use the comparison snapshot. No optimizer or automatic parameter selection is used.")
                st.plotly_chart(comparison_figures(preview, units)[0], key="sensitivity_plot", width="stretch")
                st.dataframe(comparison_frame(preview, units), hide_index=True)
