"""Streamlit PVT workflow. Canonical drafts survive project unit changes."""
from dataclasses import asdict, replace
from hashlib import sha256
from pathlib import Path
import pandas as pd
import streamlit as st
from material_balance_studio.io.history import read_table
from material_balance_studio.pvt.fluid import BlackOilFluid
from material_balance_studio.pvt.laboratory import LaboratoryData, lab_from_frame, lab_to_frame, PROPERTIES, QUANTITIES
from material_balance_studio.pvt.correlations.registry import candidates, REGISTRY
from material_balance_studio.pvt.correlations.gas_pseudocritical import PSEUDO_CRITICAL_METHODS
from material_balance_studio.pvt.models.correlation_model import CorrelationPVTModel, PropertySelection
from material_balance_studio.pvt.matching.ranking import screen_property, recommendation
from material_balance_studio.pvt.matching.regression import fit_property
from material_balance_studio.pvt.matching.metrics import error_metrics
from material_balance_studio.pvt.qc.physical_checks import physical_qc, failures
from material_balance_studio.units.display import to_display, from_display, column_label
from material_balance_studio.units.temperature import temperature_to_kelvin, temperature_from_kelvin
from .pvt import comparison_figure, screening_frame, fit_frame, selected_model_frame, property_detail_frame, PROPERTY_LABELS

DIMENSIONAL = {"temperature": "temperature", "pb": "pressure", "rsb": "rs", "low": "pressure", "high": "pressure"}
DEFAULTS = {"temperature": 363.15, "pb": 18e6, "rsb": 90., "low": 2e6, "high": 35e6}


def _shown(value, quantity, units):
    return temperature_from_kelvin(value, units) if quantity == "temperature" else to_display(value, quantity, units)


def refresh_pvt_units():
    if "pvt_fluid_si" not in st.session_state:
        return
    for name, quantity in DIMENSIONAL.items():
        st.session_state[f"fluid_{name}"] = _shown(st.session_state["pvt_fluid_si"][name], quantity, st.session_state["display_units"])


def _save(name):
    quantity, units = DIMENSIONAL[name], st.session_state["display_units"]
    value = st.session_state[f"fluid_{name}"]
    try:
        st.session_state["pvt_fluid_si"][name] = (temperature_to_kelvin(value, units) if quantity == "temperature"
                                                  else from_display(value, quantity, units))
        st.session_state.pop(f"fluid_error_{name}", None)
    except ValueError as exc:
        st.session_state[f"fluid_error_{name}"] = str(exc)


def _number(name, label, units):
    quantity = DIMENSIONAL[name]
    st.session_state.setdefault(f"fluid_{name}", _shown(st.session_state["pvt_fluid_si"][name], quantity, units))
    label = f"{label} ({'°C' if str(units) in ('SI', 'UnitSystem.SI') else '°F'})" if quantity == "temperature" else column_label(label, quantity, units)
    st.number_input(label, key=f"fluid_{name}", format="%.8g", on_change=_save, args=(name,))
    if st.session_state.get(f"fluid_error_{name}"):
        raise ValueError(st.session_state[f"fluid_error_{name}"])
    return st.session_state["pvt_fluid_si"][name]


def _lab_editor(units):
    st.caption("Sparse measurements: leave unavailable properties blank. Pressure is mandatory for a populated row. Duplicate pressures are rejected.")
    upload = st.file_uploader("Upload laboratory CSV / Excel", type=["csv", "xlsx"], key="lab_upload")
    source_units = st.selectbox("Laboratory file units", ["SI", "FIELD"], key="lab_file_units", disabled=upload is None)
    st.caption("File headers: pressure, bo, rs, oil_viscosity, z, bg, bw. SI files use Pa and Pa·s; FIELD files use psia, cP, scf/STB and rb/scf. Editor pressure uses project display units.")
    st.session_state.setdefault("pvt_lab", LaboratoryData())
    signature = None if upload is None else sha256(upload.getvalue() + source_units.encode()).hexdigest()
    if signature != st.session_state.get("lab_upload_signature"):
        if upload is not None:
            upload.seek(0)
            st.session_state["pvt_lab"] = lab_from_frame(read_table(upload, upload.name), source_units)
        st.session_state["lab_upload_signature"] = signature
        st.session_state.pop("lab_editor_context", None)
    if st.button("Load bundled sparse lab example", key="load_lab_example"):
        source = Path(__file__).resolve().parents[3] / "examples" / "correlation_lab_SI.csv"
        st.session_state["pvt_lab"] = lab_from_frame(pd.read_csv(source))
        st.session_state.pop("lab_editor_context", None)
    context = (str(units), signature)
    if st.session_state.get("lab_editor_context") != context:
        # Keep editor input fixed until its context changes; repeated added rows
        # must not be reapplied to the editor's already-edited output.
        frame = lab_to_frame(st.session_state["pvt_lab"], units)
        if str(units) in ("SI", "UnitSystem.SI"):
            frame["pressure"] = frame["pressure"] / 1e6
        st.session_state["lab_editor_base"] = frame
        st.session_state["lab_editor_canonical"] = st.session_state["pvt_lab"]
        st.session_state["lab_editor_revision"] = st.session_state.get("lab_editor_revision", 0) + 1
        st.session_state["lab_editor_context"] = context
    edited = st.data_editor(st.session_state["lab_editor_base"], num_rows="dynamic", hide_index=True,
        key=f"lab_editor_{st.session_state['lab_editor_revision']}",
        column_config={k: st.column_config.NumberColumn(column_label(k, QUANTITIES[k], units), format="%.8g")
                       for k in ("pressure", *PROPERTIES)}, width="stretch")
    imported = edited.copy()
    if str(units) in ("SI", "UnitSystem.SI"):
        imported["pressure"] *= 1e6
    lab = (st.session_state["lab_editor_canonical"] if edited.equals(st.session_state["lab_editor_base"])
           else lab_from_frame(imported, units))
    st.session_state["pvt_lab"] = lab
    st.download_button("Download lab data in import format", lab_to_frame(lab, units).to_csv(index=False),
                       f"laboratory_{units.value}.csv", "text/csv")
    st.caption("Download uses file units (SI pressure Pa), with blank measurements preserved.")
    return lab


def correlation_workflow(units, initial_pressure, source):
    overview, fluid_tab, lab_tab, screening_tab, match_tab, selected_tab, qc_tab = st.tabs(
        ["PVT Overview", "Fluid Inputs", "Laboratory Data", "Correlation Screening", "Correlation Matching", "Selected PVT Model", "PVT QC"])
    with overview:
        st.write("Define fluid → add available measurements → screen each property → select correlations → optionally match → review QC.")
        st.caption("Defaults are synthetic. Supply at least measured Pb or Rsb; initial pressure is not assumed to be Pb. Predictions are continuous inside the requested range.")
    try:
        with fluid_tab:
            st.session_state.setdefault("pvt_fluid_si", dict(DEFAULTS))
            api = st.number_input("Oil API gravity", value=35., key="fluid_api")
            gg = st.number_input("Gas specific gravity (air = 1)", value=.75, key="fluid_gg")
            temperature = _number("temperature", "Reservoir temperature", units)
            st.write(f"Initial reservoir pressure from setup: {to_display(initial_pressure, 'pressure', units):.8g}")
            has_pb = st.checkbox("Measured Pb available", value=True, key="has_pb")
            pb = _number("pb", "Measured Pb", units) if has_pb else None
            has_rsb = st.checkbox("Measured Rsb available", value=True, key="has_rsb")
            rsb = _number("rsb", "Measured Rsb", units) if has_rsb else None
            oil_sg = st.number_input("Measured oil SG (0 = unavailable)", value=0., key="fluid_oil_sg")
            water_sg = st.number_input("Water specific gravity", value=1., key="fluid_water_sg")
            salinity = st.number_input("Water salinity (mass fraction, 0–0.35)", value=0., key="fluid_salinity", format="%.6g")
            gas_basis = st.radio("Gas composition basis", ["sweet", "sour"], horizontal=True, key="fluid_gas_basis")
            low, high = _number("low", "Minimum project pressure", units), _number("high", "Maximum project pressure", units)
            fluid = BlackOilFluid(api, gg, temperature, initial_pressure, pb, rsb, oil_sg or None, water_sg, salinity, gas_basis)
            st.caption("Vasquez–Beggs uses gas gravity on its 100-psig separator basis. McCain Bw does not explicitly use the entered water SG or salinity.")
        with lab_tab:
            lab = _lab_editor(units)
        with screening_tab:
            use_pb = st.checkbox("Use measured Pb as the final bubble-point anchor", value=True, key="use_measured_pb", disabled=pb is None)
            pseudo_method = st.selectbox("Gas pseudo-critical property method", list(PSEUDO_CRITICAL_METHODS),
                                        format_func=lambda key: PSEUDO_CRITICAL_METHODS[key].name, key="select_pseudo_critical")
            st.caption("The pseudo-critical method supplies Ppc and Tpc. The independently selected gas z-factor correlation uses Pr and Tr.")
            selections = {}
            for name in vars(PropertySelection()):
                options = [c.key for c in candidates(name)]
                selections[name] = st.selectbox(f"Selected {name} correlation", options, format_func=lambda key: REGISTRY[key].name, key=f"select_{name}")
            raw = CorrelationPVTModel(fluid, (low, high), PropertySelection(**selections), use_measured_pb=use_pb, pseudo_critical_method=pseudo_method)
            signature = repr((fluid, (low, high), raw.selection, lab, use_pb, pseudo_method, "phase2.1"))
            if signature != st.session_state.get("pvt_fit_signature"):
                had_fits = bool(st.session_state.get("pvt_fits"))
                st.session_state["pvt_fits"] = {}
                st.session_state["pvt_fit_signature"] = signature
                st.session_state["pvt_screening_results"] = {}
                if had_fits:
                    st.info("Fluid, range, lab or selection changed. Previous matches were cleared; refit the updated model.")
            fits = st.session_state["pvt_fits"]
            property = st.selectbox("Property to screen / match", ["pb", *PROPERTIES], index=1, key="screen_property")
            results = screen_property(raw, lab, property)
            st.session_state["pvt_screening_results"][property] = results
            recommended = recommendation(results)
            st.write(f"Recommended: {REGISTRY[recommended].name}" if recommended else "No numerical recommendation: independent measurements are unavailable or candidates are inapplicable.")
            st.caption("Default ranking: VALID first, CAUTION next, then lowest RMSE within each tier. Click column headings to sort. Recommendations do not change the selections above. Zero measurements are omitted only from percentage metrics.")
            table = screening_frame(results, selections[property], [r.correlation for r in fits.values()], units)
            st.dataframe(table, hide_index=True, width="stretch")
            st.download_button("Download screening results", table.to_csv(index=False), f"screening_{property}_{units.value}.csv", "text/csv")
            import json
            st.download_button("Download full screening predictions (JSON, canonical SI)",
                               json.dumps([asdict(r) for r in results], indent=2), f"screening_{property}_SI.json", "application/json")
            with st.expander("Selected formulation and reference"):
                entry = REGISTRY[selections[property]]
                st.write(entry.formulation)
                st.write(entry.notes)
                st.markdown(f"[Correlation source]({entry.reference}) · Exact equations: docs/pvt_phase_2.md")
        with match_tab:
            st.caption("Affine: A + B × raw. Rs uses an Rsb-preserving slope; Bg/Pb use scale only. Bg matching also scales z. All accepted fits must pass physical QC over the project range.")
            pending = st.session_state.pop("pvt_activate_pending", None)
            if pending:
                st.session_state[f"activate_match_{pending}"] = True
            for name in fits:
                st.session_state.setdefault(f"activate_match_{name}", True)
            active_names = [name for name in fits if st.checkbox(f"Use matched {name}", key=f"activate_match_{name}")]
            def assembled():
                from material_balance_studio.pvt.models.correlation_model import Transform
                return replace(raw, transforms=tuple(Transform(name, fits[name].a, fits[name].b) for name in active_names if name in fits))
            comparison = assembled()
            if st.button(f"Match selected {property}", key="fit_selected"):
                # A changed upstream curve invalidates downstream fits and their metrics.
                dependents = {"pb": set(fits), "rs": {"bo", "oil_viscosity"}, "z": {"bg"}, "bg": {"z"}}.get(property, set())
                candidate_base = replace(comparison, transforms=tuple(t for t in comparison.transforms if t.property not in dependents))
                attempt = fit_property(candidate_base, lab, property)
                if attempt.accepted:
                    for name in dependents:
                        fits.pop(name, None)
                    fits[property] = attempt.record
                    st.session_state["pvt_activate_pending"] = property
                    st.session_state["pvt_last_match"] = f"Accepted {property}: RMSE {attempt.record.before.rmse:.6g} → {attempt.record.after.rmse:.6g} (canonical units)."
                    st.rerun()
                else:
                    st.error("Match rejected: " + " · ".join(attempt.reasons))
                    if attempt.record is not None:
                        st.dataframe(fit_frame([attempt.record], units), hide_index=True)
            if fits:
                st.success("Accepted matches are retained separately from the raw model.")
                st.dataframe(fit_frame(list(fits.values()), units), hide_index=True, width="stretch")
                st.caption("Fit records describe each regression at acceptance. Current combined model errors are shown in QC; disabling an upstream match can change them.")
            else:
                st.info("No accepted matches. Add measurements and match a selected property.")
        comparison = assembled()
        active = comparison if source == "Matched Correlation PVT" else raw
        raw_qc = physical_qc(raw, lab)
        qc = physical_qc(active, lab)
        with selected_tab:
            st.write(f"Active source: {source}")
            st.caption(f"{column_label('Active Pb', 'pressure', units)}: {to_display(active.pb, 'pressure', units):.8g} · "
                       f"{column_label('Rsb', 'rs', units)}: {to_display(active.rsb, 'rs', units):.8g}")
            st.dataframe(selected_model_frame(active, qc), hide_index=True, width="stretch")
            for name in selections:
                with st.expander(f"{PROPERTY_LABELS[name]} · {REGISTRY[selections[name]].name}", expanded=name == "pb"):
                    st.dataframe(property_detail_frame(active, name, fits.get(name), units), hide_index=True, width="stretch")
            st.caption("Plots compare the raw curve with enabled matches. The active source above determines what the MBE engine uses.")
            for name in PROPERTIES:
                try:
                    st.plotly_chart(comparison_figure(raw, comparison if comparison.transforms else None, lab, name, units),
                                    width="stretch", key=f"pvt_plot_{name}")
                except (ValueError, ArithmeticError, TypeError) as exc:
                    st.error(f"{name} curve unavailable: {exc}")
            export = {"source": source, "canonical_fluid": asdict(fluid), "pressure_bounds_Pa": [low, high],
                      "selections": selections, "use_measured_pb": use_pb, "pseudo_critical_method": pseudo_method,
                      "property_summary": selected_model_frame(active, qc).to_dict(orient="records"),
                      "active_transforms": [asdict(t) for t in active.transforms],
                      "fit_records": [asdict(r) for r in fits.values()]}
            import json
            st.download_button("Download selected model and match records (JSON, canonical SI)", json.dumps(export, indent=2), "pvt_model.json", "application/json")
        with qc_tab:
            st.write("Active model QC")
            st.dataframe(pd.DataFrame([asdict(item) for item in qc]), hide_index=True, width="stretch")
            with st.expander("Raw model QC before matching"):
                st.dataframe(pd.DataFrame([asdict(item) for item in raw_qc]), hide_index=True, width="stretch")
            errors = []
            for name in PROPERTIES:
                data = lab.measurements(name)
                if data:
                    try:
                        metric = error_metrics([v for _, v in data], [active.property_at_pressure(name, p) for p, _ in data])
                        errors.append({"Property": name, "MAPE (%)": metric.mape, "Maximum error (%)": metric.max_ape,
                                       "RMSE (display units)": to_display(metric.rmse, QUANTITIES[name], units),
                                       "n": metric.n})
                    except (ValueError, ArithmeticError, TypeError):
                        pass  # already reported as a QC failure
            if errors:
                st.dataframe(pd.DataFrame(errors), hide_index=True)
            st.caption("QC checks 1,001 evenly spaced pressures, lab pressures, 0.99/0.999/1/1.001/1.01 Pb and Pb ±1e-6. Maximum full-range change means |matched − raw|, not error against unknown laboratory values. QC does not establish accuracy outside measured or published ranges.")
        if failures(qc):
            st.error("PVT QC has failures. Resolve them before running simulation.")
            return None
        if source == "Matched Correlation PVT" and not comparison.transforms:
            st.warning("Matched mode requires an accepted, enabled match. Fit a measured property or choose Correlation PVT.")
            return None
        st.session_state["pvt_active_model"] = active
        return active
    except (ValueError, ArithmeticError, TypeError) as exc:
        st.error(f"PVT input: {exc}")
        return None
