"""Aquifer configuration widgets and plots; committed state stays in the engine."""
from dataclasses import asdict
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from material_balance_studio.aquifer.registry import MODELS
from material_balance_studio.units.display import from_display, to_display, column_label, display_unit, format_value


FIELDS = {
    "capacity": ("Aquifer capacity Caq", "aquifer_capacity", .005,
                 "Connected initial aquifer water volume × total compressibility. Reservoir volume released per unit pressure decline."),
    "productivity_index": ("Aquifer productivity index J", "aquifer_productivity", 1e-9,
                           "Reservoir water-flow rate per aquifer-to-reservoir pressure difference. Schilthuis retains fixed reference pressure; Fetkovich depletes."),
    "initial_water_volume": ("Connected initial aquifer water volume Wi", "reservoir_volume", 5e6,
                             "Water volume at initial aquifer conditions, including only the connected aquifer fraction. This is a reservoir volume."),
    "total_compressibility": ("Aquifer total compressibility ct", "compressibility", 1e-9,
                             "Sum of aquifer water and pore-volume compressibilities; Caq = Wi × ct."),
    "inner_radius": ("Reservoir/aquifer inner radius", "length", 1800.,
                     "Reference radius at the reservoir-aquifer contact used in tD and B."),
    "radius_ratio": ("Aquifer outer/inner radius ratio", "dimensionless", 10.,
                     "Outer aquifer radius divided by the reservoir/aquifer contact radius. Phase 3B uses infinite-acting response functions but validates geometry."),
    "thickness": ("Aquifer thickness", "length", 20.,
                  "Net aquifer thickness communicating with the reservoir."),
    "porosity": ("Aquifer porosity", "dimensionless", .2,
                 "Connected aquifer pore-volume fraction."),
    "permeability": ("Aquifer permeability", "permeability", 1e-13,
                     "Aquifer permeability used in the radial dimensionless time."),
    "water_viscosity": ("Aquifer water viscosity", "viscosity", .001,
                        "Water viscosity in the aquifer."),
    "encroachment_angle": ("Encroachment angle", "angle", 360.,
                           "Full-circle aquifer is 360 degrees; partial aquifers scale the radial aquifer constant."),
}

TRANSIENT_NAMES = ("inner_radius", "radius_ratio", "thickness", "porosity", "permeability",
                   "water_viscosity", "total_compressibility", "encroachment_angle")


def refresh_aquifer_units():
    draft = st.session_state.get("aquifer_draft", {})
    for name, (_, quantity, _, _) in FIELDS.items():
        if name in draft:
            st.session_state["aq_"+name] = to_display(draft[name], quantity, st.session_state["display_units"])


def save_aquifer_value(name):
    st.session_state["aquifer_draft"][name] = from_display(st.session_state["aq_"+name], FIELDS[name][1], st.session_state["display_units"])


def aquifer_inputs(units, initial_pressure):
    st.subheader("Aquifer")
    selected = st.radio("Aquifer type", list(MODELS), horizontal=True, key="aquifer_type")
    st.session_state.setdefault("aquifer_draft", {k:v[2] for k,v in FIELDS.items()})
    names = {"None": (), "Pot": ("capacity",), "Schilthuis": ("productivity_index",),
             "Fetkovich": ("initial_water_volume", "total_compressibility", "productivity_index"),
             "Carter-Tracy": TRANSIENT_NAMES, "Van Everdingen-Hurst": TRANSIENT_NAMES,
             "Modified Van Everdingen-Hurst": TRANSIENT_NAMES}[selected]
    if names:
        st.caption(f"Initial aquifer pressure equals initial reservoir pressure: {format_value(to_display(initial_pressure, 'pressure', units))} {display_unit('pressure', units)}. Each run starts a fresh aquifer state.")
        st.caption("Signed flow is retained: pressure recovery can return water to the aquifer. Negative influx and decreasing cumulative influx are flagged for review.")
        if selected in ("Carter-Tracy", "Van Everdingen-Hurst", "Modified Van Everdingen-Hurst"):
            st.caption("Infinite-acting radial response. Radius ratio is geometry metadata only; finite boundaries are not modelled.")
        if selected == "Van Everdingen-Hurst":
            st.caption("Original Van Everdingen-Hurst — infinite-acting radial; pressure-step superposition.")
        if selected == "Modified Van Everdingen-Hurst":
            st.caption("Modified Van Everdingen-Hurst — linear-pressure-history formulation, infinite-acting radial; PETEX Appendix C C2.9.")
    for name in names:
        label, quantity, _, help_text = FIELDS[name]
        st.session_state.setdefault("aq_"+name, to_display(st.session_state["aquifer_draft"][name], quantity, units))
        st.number_input(column_label(label, quantity, units), key="aq_"+name, format="%.8g",
                        help=help_text, on_change=save_aquifer_value, args=(name,))
    try:
        return MODELS[selected](**{name:st.session_state["aquifer_draft"][name] for name in names})
    except ValueError as exc:
        st.error(str(exc))
        return None


def aquifer_setup_frame(model, units):
    rows = [{"Quantity": "Aquifer model", "Value": model.key, "Unit": ""}]
    for name, value in asdict(model).items():
        label, quantity, _, _ = FIELDS[name]
        rows.append({"Quantity": label, "Value": format_value(to_display(value, quantity, units)), "Unit": display_unit(quantity, units)})
    return pd.DataFrame(rows)


def aquifer_figures(result, units):
    initial = result.initial_state.aquifer_state
    if initial is None or initial.model_key == "none":
        return []
    states = (result.initial_state, *result.states)
    series = [("Cumulative Aquifer Influx We", "reservoir_volume", [s.balance.aquifer_support for s in states])]
    if initial.model_key in ("schilthuis", "fetkovich", "carter_tracy", "van_everdingen_hurst", "modified_van_everdingen_hurst"):
        series.append(("Average Aquifer Influx Rate", "aquifer_rate", [s.aquifer_step.average_influx_rate if s.aquifer_step else None for s in states]))
    if initial.model_key in ("fetkovich", "carter_tracy", "van_everdingen_hurst"):
        series.append(("Aquifer Pressure", "pressure", [s.aquifer_state.aquifer_pressure for s in states]))
    if initial.model_key in ("carter_tracy", "van_everdingen_hurst", "modified_van_everdingen_hurst"):
        series.append(("Aquifer Dimensionless Time tD", "dimensionless",
                       [dict(s.aquifer_state.model_variables).get("diag_tD") if s.aquifer_step else None for s in states]))
    figures = []
    for label, quantity, values in series:
        figure = go.Figure(go.Scatter(x=[s.date.isoformat() for s in states],
            y=[to_display(v, quantity, units) for v in values], mode="lines+markers", name=label))
        figure.update_layout(title=label+" vs Time", xaxis_title="Date", yaxis_title=column_label(label, quantity, units))
        figures.append(figure)
    return figures
