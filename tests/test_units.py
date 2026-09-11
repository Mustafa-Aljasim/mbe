from datetime import date

import pandas as pd
import pytest

from material_balance_studio.domain.models import CumulativeVolumes
from material_balance_studio.io.history import history_from_frame, pvt_from_frame
from material_balance_studio.io.reservoir import tank_from_inputs
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.solver.pressure_solver import solve_pressure
from material_balance_studio.units.conversions import from_si, to_si


@pytest.mark.parametrize("quantity", ["pressure", "liquid_volume", "gas_volume", "compressibility",
                                      "bo", "bw", "bg", "rs", "viscosity"])
def test_field_si_roundtrip(quantity):
    assert from_si(to_si(123.4, quantity, "FIELD"), quantity, "FIELD") == pytest.approx(123.4)


def test_field_gas_basis_is_explicit():
    assert to_si(1, "liquid_volume", "FIELD") == pytest.approx(.158987294928)
    assert to_si(1, "gas_volume", "FIELD") == pytest.approx(.028316846592)
    assert to_si(.001, "bg", "FIELD") == pytest.approx(.00561458333333)
    assert to_si(1000, "rs", "FIELD") == pytest.approx(178.107606679)


def test_complete_field_case_matches_si_case(tank):
    fields = [dict(pressure=from_si(r.pressure, "pressure", "FIELD"), bo=r.bo,
                   rs=from_si(r.rs, "rs", "FIELD"), bg=from_si(r.bg, "bg", "FIELD"), bw=r.bw)
              for r in tank.pvt_model.table.rows]
    converted = tank_from_inputs(
        initial_date=date(2020, 1, 1), initial_pressure=from_si(30e6, "pressure", "FIELD"),
        oil_in_place=from_si(1e6, "liquid_volume", "FIELD"), swc=.2, cf=0, cw=0, m=0,
        pvt_model=TablePVTModel(pvt_from_frame(pd.DataFrame(fields), "FIELD")), units="FIELD")
    volumes = CumulativeVolumes(np=20000, gp=1600000, wp=100, winj=2000, ginj=10000)
    row = {"date": "2020-02-01", **{name: from_si(getattr(volumes, name),
           "gas_volume" if name in ("gp", "ginj") else "liquid_volume", "FIELD")
           for name in ("np", "gp", "wp", "winj", "ginj")}}
    history = history_from_frame(pd.DataFrame([row]), "FIELD")
    si_result = solve_pressure(tank, volumes, 30e6)
    field_result = solve_pressure(converted, history[0].cumulative, 30e6)
    assert si_result.diagnostics.converged and field_result.diagnostics.converged
    assert si_result.pressure == pytest.approx(field_result.pressure, abs=.01)
