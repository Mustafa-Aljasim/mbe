"""Plotly charts with shared engineering labels and tolerance-scaled closure QC."""

import plotly.graph_objects as go

from material_balance_studio.diagnostics.qc import PASS_RELATIVE_TOLERANCE, WARNING_RELATIVE_TOLERANCE
from material_balance_studio.domain.models import SimulationResult
from material_balance_studio.units.conversions import UnitSystem
from material_balance_studio.units.display import column_label, display_unit
from .tables import RELATIVE_LABEL, engineering_results_frame


def pressure_figure(result: SimulationResult, units: UnitSystem | str) -> go.Figure:
    frame = engineering_results_frame(result, units)
    figure = go.Figure()
    for label, name, mode in (("Calculated Pressure", "Calculated", "lines+markers"),
                              ("Observed Pressure", "Observed", "markers")):
        figure.add_scatter(x=frame["Date"].tolist(), y=frame[column_label(label, "pressure", units)].tolist(),
                           mode=mode, name=name)
    figure.update_layout(title="Average reservoir pressure", xaxis_title="Date",
                         yaxis_title=f"Average reservoir pressure ({display_unit('pressure', units)})")
    return figure


def closure_figure(result: SimulationResult) -> go.Figure:
    """Linear zero-based scale includes the PASS limit, so roundoff stays near zero.

    Values are neither clipped nor rounded. WARNING/FAIL thresholds and hover
    scientific notation retain meaningful differences without magnifying noise.
    """
    errors = [state.balance.relative_residual for state in result.states]
    figure = go.Figure(go.Scatter(
        x=[state.date.isoformat() for state in result.states], y=errors,
        mode="lines+markers", name="Absolute relative residual",
        hovertemplate="%{x}<br>Relative residual: %{y:.3e}<extra></extra>",
    ))
    figure.add_hline(y=PASS_RELATIVE_TOLERANCE, line_dash="dash", line_color="green",
                     annotation_text="PASS limit · 1e-8")
    upper = max(PASS_RELATIVE_TOLERANCE, max(errors, default=0.0)) * 1.2
    if upper >= WARNING_RELATIVE_TOLERANCE:
        figure.add_hline(y=WARNING_RELATIVE_TOLERANCE, line_dash="dash", line_color="orange",
                         annotation_text="WARNING limit · 1e-5")
    figure.update_layout(title="Material-balance closure QC", xaxis_title="Date",
                         yaxis_title=RELATIVE_LABEL,
                         yaxis={"range": [0, upper], "tickformat": ".1e"})
    return figure
