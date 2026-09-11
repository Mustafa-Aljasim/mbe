"""Unit-aware PVT comparisons and engineering reporting; no equations here."""
import pandas as pd
import plotly.graph_objects as go
from material_balance_studio.pvt.laboratory import QUANTITIES
from material_balance_studio.pvt.correlations.registry import REGISTRY
from material_balance_studio.pvt.correlations.gas_pseudocritical import PSEUDO_CRITICAL_METHODS
from material_balance_studio.pvt.qc.physical_checks import pressure_grid, property_qc_status
from material_balance_studio.units.display import to_display, column_label, display_unit, format_value

PROPERTY_LABELS = {"pb": "Bubble point", "bo": "Bo", "rs": "Rs", "oil_viscosity": "Oil viscosity",
                   "z": "Gas z-factor", "bg": "Bg", "bw": "Bw", "co": "Undersaturated oil compressibility"}


def comparison_figure(raw, matched, lab, property, units):
    quantity = QUANTITIES[property]
    grid = sorted(set(pressure_grid(raw, lab)) | (set(pressure_grid(matched, lab)) if matched is not None else set()))
    fig = go.Figure()
    for name, model, color, dash in (("Raw correlation", raw, "#2563eb", "solid"),
                                      ("Matched correlation", matched, "#16a34a", "dash")):
        if model is None:
            continue
        fig.add_scatter(x=[to_display(float(p), "pressure", units) for p in grid],
                        y=[to_display(model.property_at_pressure(property, float(p)), quantity, units) for p in grid],
                        name=name, mode="lines", line=dict(color=color, dash=dash, width=2.5),
                        hovertemplate="Pressure: %{x:.6g}<br>Value: %{y:.6g}<extra>%{fullData.name}</extra>")
    points = lab.measurements(property)
    if points:
        fig.add_scatter(x=[to_display(p, "pressure", units) for p, _ in points],
                        y=[to_display(v, quantity, units) for _, v in points],
                        name="Laboratory", mode="markers",
                        marker=dict(color="#f59e0b", size=10, symbol="diamond", line=dict(color="#475569", width=1)))
    fig.add_vline(x=to_display(raw.pb, "pressure", units), line_dash="dot", line_color="#94a3b8",
                  annotation_text="Pb (raw = matched)" if matched is not None and matched.pb == raw.pb else "Pb (raw)")
    if matched is not None and matched.pb != raw.pb:
        fig.add_vline(x=to_display(matched.pb, "pressure", units), line_dash="dash", line_color="#94a3b8", annotation_text="Pb (matched)")
    fig.update_layout(title=f"{PROPERTY_LABELS[property]} vs Pressure", xaxis_title=column_label("Pressure", "pressure", units),
                      yaxis_title=column_label(PROPERTY_LABELS[property], quantity, units),
                      legend=dict(orientation="h", y=1.16, x=0), height=430,
                      margin=dict(l=70, r=30, t=100, b=65), hovermode="x unified")
    return fig


def screening_frame(results, selected, fits, units):
    rows = []
    for result in results:
        metric = result.metrics
        quantity = QUANTITIES[result.property]
        rows.append({"Correlation": REGISTRY[result.correlation].name,
                     "Applicable?": result.applicability.status,
                     column_label("RMSE", quantity, units): to_display(metric.rmse, quantity, units) if metric else None,
                     column_label("MAE", quantity, units): to_display(metric.mae, quantity, units) if metric else None,
                     "MAPE (%)": metric.mape if metric else None,
                     column_label("Bias", quantity, units): to_display(metric.bias, quantity, units) if metric else None,
                     "Maximum absolute error (%)": metric.max_ape if metric else None,
                     "Points": metric.n if metric else 0,
                     "Percentage points": metric.percentage_n if metric else 0,
                     "Match available?": result.correlation in fits,
                     "Selected?": result.correlation == selected,
                     "Notes": "; ".join(result.applicability.reasons + result.errors)})
    return pd.DataFrame(rows)


def fit_frame(records, units):
    return pd.DataFrame([{"Property": r.property, "Correlation": REGISTRY[r.correlation].name,
                          "A (display units)": to_display(r.a, QUANTITIES[r.property], units), "B": r.b, "n": r.n,
                          "Property unit": column_label(r.property, QUANTITIES[r.property], units),
                          "Pre RMSE": to_display(r.before.rmse, QUANTITIES[r.property], units),
                          "Post RMSE": to_display(r.after.rmse, QUANTITIES[r.property], units),
                          "Pre bias": to_display(r.before.bias, QUANTITIES[r.property], units),
                          "Post bias": to_display(r.after.bias, QUANTITIES[r.property], units),
                          "Maximum full-range change": to_display(r.maximum_full_range_deviation, QUANTITIES[r.property], units),
                          "Maximum full-range change (%)": r.maximum_full_range_deviation_percent,
                          column_label("Pressure at maximum change", "pressure", units): to_display(r.maximum_deviation_pressure, "pressure", units),
                          "QC before": r.qc_before, "QC after": r.qc_after,
                          "Formulation": r.formulation} for r in records])


def selected_model_frame(model, qc):
    rows = []
    matched = {t.property for t in model.transforms}
    for property, key in vars(model.selection).items():
        treatment = "Matched" if property in matched else "Raw"
        if property == "pb" and model.use_measured_pb and model.fluid.pb is not None:
            treatment = "Measured anchor"
        if property in ("rs", "bo", "oil_viscosity") and property not in matched:
            treatment = "Raw, saturation anchored"
        if property == "co" and "bo" in matched:
            treatment = "Derived from matched Bo"
        if property == "z" and "bg" in matched:
            treatment = "Scaled by matched Bg"
        rows.append({"Property": PROPERTY_LABELS[property], "Method": REGISTRY[key].name,
                     "Treatment": treatment, "QC": property_qc_status(qc, property)})
    rows.insert(4, {"Property": "Gas pseudo-critical", "Method": PSEUDO_CRITICAL_METHODS[model.pseudo_critical_method].name,
                    "Treatment": "Pseudo-critical properties", "QC": property_qc_status(qc, "pseudo_critical")})
    return pd.DataFrame(rows)


def property_detail_frame(model, property, record, units):
    """Readable vertical details avoid a wide, horizontally clipped fit summary."""
    quantity = QUANTITIES[property]
    unit = display_unit(quantity, units)
    shown = lambda value: f"{format_value(to_display(value, quantity, units))} {unit}" if value is not None else "Unavailable"
    transform = next((t for t in model.transforms if t.property == property), None)
    rows = [{"Item": "Matched", "Value": "Yes" if transform else "No"}]
    if property == "pb":
        try:
            calculated = shown(model.predicted_pb)
        except (ValueError, ArithmeticError, TypeError) as exc:
            calculated = f"Unavailable: {exc}"
        rows += [{"Item": "Measured Pb", "Value": shown(model.fluid.pb)},
                 {"Item": "Calculated Pb (untuned)", "Value": calculated},
                 {"Item": "Active Pb", "Value": shown(model.pb)}]
    if transform:
        rows += [{"Item": "A", "Value": shown(transform.a)}, {"Item": "B", "Value": format_value(transform.b)}]
        if record is not None:
            rows += [{"Item": "RMSE before", "Value": shown(record.before.rmse)},
                     {"Item": "RMSE after", "Value": shown(record.after.rmse)},
                     {"Item": "Maximum full-range change", "Value": shown(record.maximum_full_range_deviation)},
                     {"Item": "QC at matching (before → after)", "Value": f"{record.qc_before} → {record.qc_after}"},
                     {"Item": "Laboratory points", "Value": str(record.n)}]
    return pd.DataFrame(rows)
