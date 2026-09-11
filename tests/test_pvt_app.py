from dataclasses import asdict
from pathlib import Path
import json
import pytest
from streamlit.testing.v1 import AppTest
from material_balance_studio.pvt.models.correlation_model import CorrelationPVTModel
from material_balance_studio.pvt.laboratory import LaboratoryData, LabPoint

ROOT=Path(__file__).resolve().parents[1]


def correlation_app():
    app=AppTest.from_file(str(ROOT/'app.py'),default_timeout=40).run()
    app.radio(key="pvt_source").set_value("Correlation PVT").run()
    assert not app.exception and not app.error
    return app


def run_button(app):
    return next(b for b in app.button if b.label=="Run simulation")


def fit_example(app):
    app.button(key="load_lab_example").click().run()
    app.selectbox(key="screen_property").select("bo").run()
    app.button(key="fit_selected").click().run()
    assert not app.exception and not app.error
    record=app.session_state["pvt_fits"]["bo"]
    assert record.after.rmse<record.before.rmse


def test_raw_correlation_ui_runs_mbe_and_six_property_plots():
    app=correlation_app()
    assert len(app.get("plotly_chart"))==6
    run_button(app).click().run()
    assert not app.exception and not app.error
    result,snapshot=app.session_state["simulation"]
    assert result.converged and len(result.states)==6
    assert isinstance(snapshot["tank"].pvt_model,CorrelationPVTModel)
    assert not snapshot["tank"].pvt_model.transforms


def test_matched_ui_retains_raw_lab_and_matched_and_runs_mbe():
    app=correlation_app()
    fit_example(app)
    app.radio(key="pvt_source").set_value("Matched Correlation PVT").run()
    assert not run_button(app).disabled
    bo_plot=json.loads(app.get("plotly_chart")[0].proto.spec)
    assert {t["name"] for t in bo_plot["data"]}=={"Raw correlation","Matched correlation","Laboratory"}
    assert bo_plot["layout"]["shapes"][0]["x0"]==18
    run_button(app).click().run()
    assert not app.exception and not app.error
    result,snapshot=app.session_state["simulation"]
    assert result.converged and len(snapshot["tank"].pvt_model.transforms)==1
    app.checkbox(key="activate_match_bo").uncheck().run()
    assert run_button(app).disabled


def test_pvt_unit_switch_preserves_lab_fluid_matches_and_saved_results():
    app=correlation_app()
    app.number_input(key="fluid_temperature").set_value(95.).run()
    fit_example(app)
    app.radio(key="pvt_source").set_value("Matched Correlation PVT").run()
    run_button(app).click().run()
    fluid=dict(app.session_state["pvt_fluid_si"])
    lab=app.session_state["pvt_lab"]
    fits=dict(app.session_state["pvt_fits"])
    saved=asdict(app.session_state["simulation"][0])
    app.selectbox(key="display_units").select("FIELD").run()
    assert not app.exception and not app.error
    assert app.number_input(key="fluid_temperature").value==pytest.approx(203)
    assert dict(app.session_state["pvt_fluid_si"])==fluid
    assert app.session_state["pvt_lab"]==lab
    assert app.session_state["pvt_fits"]==fits
    assert asdict(app.session_state["simulation"][0])==saved
    app.selectbox(key="display_units").select("SI").run()
    assert app.number_input(key="fluid_temperature").value==95
    assert app.session_state["pvt_lab"]==lab
    assert app.session_state["pvt_fits"]==fits


def test_manual_override_clears_stale_matches():
    app=correlation_app()
    fit_example(app)
    app.selectbox(key="select_bo").select("glaso_bo").run()
    assert not app.exception and not app.error
    assert app.session_state["pvt_active_model"].selection.bo=="glaso_bo"
    assert not app.session_state["pvt_fits"]


def test_refit_existing_property_and_raw_mode_selection():
    app=correlation_app()
    fit_example(app)
    app.checkbox(key="activate_match_bo").uncheck().run()
    app.button(key="fit_selected").click().run()
    assert not app.exception and not app.error
    assert app.checkbox(key="activate_match_bo").value
    assert not app.session_state["pvt_active_model"].transforms  # source is Raw


def test_missing_anchor_and_sour_gas_fail_gracefully():
    app=correlation_app()
    app.checkbox(key="has_pb").uncheck().run()
    app.checkbox(key="has_rsb").uncheck().run()
    assert not app.exception and run_button(app).disabled
    assert any("Provide measured" in e.value for e in app.error)
    app.checkbox(key="has_pb").check().run()
    app.checkbox(key="has_rsb").check().run()
    app.radio(key="fluid_gas_basis").set_value("sour").run()
    assert not app.exception and run_button(app).disabled


def test_invalid_temperature_callback_is_recoverable():
    app=correlation_app()
    app.number_input(key="fluid_temperature").set_value(-300.).run()
    assert not app.exception and run_button(app).disabled
    app.number_input(key="fluid_temperature").set_value(90.).run()
    assert not app.exception and not app.error and not run_button(app).disabled


def test_manual_editor_changes_are_canonical_and_survive_reruns():
    app=correlation_app()
    revision=app.session_state["lab_editor_revision"]
    key=f"lab_editor_{revision}"
    # AppTest has no direct dataframe-cell API; exercise Streamlit's editor delta.
    app.session_state[key]={"edited_rows":{},"added_rows":[{"pressure":10.,"bo":1.2}],"deleted_rows":[]}
    app.run()
    assert not app.exception and not app.error
    assert app.session_state["pvt_lab"]==LaboratoryData((LabPoint(10e6,bo=1.2),))
    # AppTest does not retain a browser-side data-editor delta automatically.
    app.session_state[key]={"edited_rows":{},"added_rows":[{"pressure":10.,"bo":1.2}],"deleted_rows":[]}
    app.run()
    assert len(app.session_state["pvt_lab"].points)==1
    app.selectbox(key="display_units").select("FIELD").run()
    assert app.session_state["pvt_lab"]==LaboratoryData((LabPoint(10e6,bo=1.2),))


def test_matched_mode_without_fit_does_not_silently_run_raw():
    app=correlation_app()
    app.radio(key="pvt_source").set_value("Matched Correlation PVT").run()
    assert not app.exception and run_button(app).disabled
    assert any("requires an accepted" in w.value for w in app.warning)
