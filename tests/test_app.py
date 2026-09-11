from pathlib import Path
from dataclasses import asdict
import json

import pytest

from streamlit.testing.v1 import AppTest

ROOT = Path(__file__).resolve().parents[1]


def test_streamlit_example_runs_and_displays_results():
    app = AppTest.from_file(str(ROOT / "app.py"), default_timeout=30).run()
    assert not app.exception
    app.button[0].click().run()
    assert not app.exception
    assert not app.error
    assert any("Converged: 6 history timesteps" in item.value for item in app.success)
    assert len(app.get("plotly_chart")) == 2
    assert len(app.session_state["simulation"][0].states) == 6
    app.selectbox[-1].select(app.session_state["simulation"][0].states[-1].date).run()
    assert not app.exception


def test_streamlit_invalid_tank_shows_error():
    app = AppTest.from_file(str(ROOT / "app.py"), default_timeout=30).run()
    app.number_input[2].set_value(1.0)
    app.button[0].click().run()
    assert not app.exception
    assert any("swc" in item.value for item in app.error)


def test_project_units_convert_example_previews_results_and_plots():
    app = AppTest.from_file(str(ROOT / "app.py"), default_timeout=30).run()
    assert app.number_input[0].label == "Initial pressure (MPa)"
    assert app.number_input[0].value == 30
    assert app.selectbox[1].disabled and app.selectbox[2].disabled
    assert app.dataframe[0].value.iloc[4]["Pressure (MPa)"] == 30
    app.button[0].click().run()
    canonical_before = asdict(app.session_state["simulation"][0])
    app.selectbox(key="display_units").select("FIELD").run()
    assert not app.exception and not app.error
    assert app.number_input[0].label == "Initial pressure (psia)"
    assert app.number_input[0].value == pytest.approx(4351.1321319)
    assert app.dataframe[0].value.iloc[4]["Pressure (psia)"] == pytest.approx(4351.1321319)
    assert app.dataframe[0].value.iloc[4]["Rs (scf/STB)"] == pytest.approx(449.1666666667)
    assert "Cumulative Water Production (bbl)" in app.dataframe[1].value
    results = next(item.value for item in app.dataframe if "Calculated Pressure (psia)" in item.value)
    assert results.iloc[0]["Calculated Pressure (psia)"] == pytest.approx(4351.1321319)
    chart = json.loads(app.get("plotly_chart")[0].proto.spec)
    assert chart["layout"]["yaxis"]["title"]["text"] == "Average reservoir pressure (psia)"
    assert chart["data"][0]["y"][0] == pytest.approx(4351.1321319)
    assert chart["data"][1]["y"][0] == pytest.approx(4351.1321319)
    assert asdict(app.session_state["simulation"][0]) == canonical_before
    assert any(metric.value == "PASS" for metric in app.metric)
    assert all("Raw" in expander.label for expander in app.expander if len(expander.json))
    app.selectbox(key="display_units").select("SI").run()
    assert app.number_input[0].value == 30
    assert asdict(app.session_state["simulation"][0]) == canonical_before


def test_switching_units_preserves_edited_reservoir_inputs():
    app = AppTest.from_file(str(ROOT / "app.py"), default_timeout=30).run()
    app.number_input(key="setup_pressure").set_value(28).run()
    app.number_input(key="setup_oil").set_value(2000000).run()
    app.number_input(key="setup_cf").set_value(.0007).run()
    draft_before = dict(app.session_state["setup_si"])
    app.selectbox(key="display_units").select("FIELD").run()
    assert app.number_input(key="setup_pressure").value == pytest.approx(4061.05665645)
    assert app.number_input(key="setup_oil").value == pytest.approx(12579621.5409)
    app.selectbox(key="display_units").select("SI").run()
    assert app.number_input(key="setup_pressure").value == 28
    assert app.number_input(key="setup_oil").value == 2000000
    assert app.number_input(key="setup_cf").value == pytest.approx(.0007)
    assert dict(app.session_state["setup_si"]) == draft_before
    app.button[0].click().run()
    assert not app.exception and not app.error
    assert app.session_state["simulation"][0].initial_state.pressure == 28e6
