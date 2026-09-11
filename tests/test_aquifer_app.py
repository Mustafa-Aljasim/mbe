from pathlib import Path
from dataclasses import asdict
import pytest
from streamlit.testing.v1 import AppTest
from material_balance_studio.units.display import to_display

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("model,key,fields,charts", [
    ("None","none",set(),2), ("Pot","pot",{"aq_capacity"},3),
    ("Schilthuis","schilthuis",{"aq_productivity_index"},4),
    ("Fetkovich","fetkovich",{"aq_initial_water_volume","aq_total_compressibility","aq_productivity_index"},5),
    ("Carter-Tracy","carter_tracy",{"aq_inner_radius","aq_radius_ratio","aq_thickness","aq_porosity","aq_permeability",
                     "aq_water_viscosity","aq_total_compressibility","aq_encroachment_angle"},6),
    ("Van Everdingen-Hurst","van_everdingen_hurst",{"aq_inner_radius","aq_radius_ratio","aq_thickness","aq_porosity","aq_permeability",
                             "aq_water_viscosity","aq_total_compressibility","aq_encroachment_angle"},6),
])
def test_aquifer_ui_model_inputs_run_and_charts(model,key,fields,charts):
    app = AppTest.from_file(str(ROOT/"app.py"),default_timeout=30).run()
    app.radio(key="aquifer_type").set_value(model).run()
    assert not app.exception
    assert {w.key for w in app.number_input if w.key and w.key.startswith("aq_")} == fields
    app.button[0].click().run()
    assert not app.exception and not app.error
    result = app.session_state["simulation"][0]
    assert result.converged and result.initial_state.aquifer_state.model_key == key
    assert len(app.get("plotly_chart")) == charts
    if model != "None":
        assert result.states[-1].balance.aquifer_support > 0
    assert any("+W_e" in latex.value for latex in app.latex)


def test_aquifer_unit_switch_preserves_inputs_and_saved_state():
    app = AppTest.from_file(str(ROOT/"app.py"),default_timeout=30).run()
    app.radio(key="aquifer_type").set_value("Fetkovich").run()
    app.number_input(key="aq_initial_water_volume").set_value(8e6).run()
    app.number_input(key="aq_total_compressibility").set_value(.0012).run()
    app.number_input(key="aq_productivity_index").set_value(100).run()
    app.button[0].click().run()
    assert not app.exception and not app.error
    snapshot = asdict(app.session_state["simulation"][0])
    draft = dict(app.session_state["aquifer_draft"])
    app.selectbox(key="display_units").select("FIELD").run()
    assert not app.exception
    assert app.number_input(key="aq_initial_water_volume").value == pytest.approx(to_display(8e6,"reservoir_volume","FIELD"))
    assert app.number_input(key="aq_productivity_index").value == pytest.approx(to_display(draft["productivity_index"],"aquifer_productivity","FIELD"))
    assert "rb/(day·psi)" in app.number_input(key="aq_productivity_index").label
    assert asdict(app.session_state["simulation"][0]) == snapshot
    app.selectbox(key="display_units").select("SI").run()
    assert app.number_input(key="aq_productivity_index").value == pytest.approx(100)
    assert dict(app.session_state["aquifer_draft"]) == draft


def test_ui_model_switch_starts_fresh_state():
    app = AppTest.from_file(str(ROOT/"app.py"),default_timeout=30).run()
    app.radio(key="aquifer_type").set_value("Fetkovich").run()
    app.button[0].click().run()
    assert app.session_state["simulation"][0].states[-1].aquifer_state.cumulative_influx > 0
    app.radio(key="aquifer_type").set_value("None").run()
    app.button[0].click().run()
    result = app.session_state["simulation"][0]
    assert result.initial_state.aquifer_state.cumulative_influx == 0
    assert all(s.aquifer_state.cumulative_influx == 0 for s in result.states)
    app.radio(key="aquifer_type").set_value("Pot").run()
    app.number_input(key="aq_capacity").set_value(7500).run()
    app.radio(key="aquifer_type").set_value("None").run()
    app.radio(key="aquifer_type").set_value("Pot").run()
    assert app.number_input(key="aq_capacity").value == 7500


def test_invalid_aquifer_blocks_run_with_engineering_message():
    app = AppTest.from_file(str(ROOT/"app.py"),default_timeout=30).run()
    app.radio(key="aquifer_type").set_value("Fetkovich").run()
    app.number_input(key="aq_total_compressibility").set_value(-.001).run()
    assert not app.exception and app.button[0].disabled
    assert any("compressibility" in item.value for item in app.error)
