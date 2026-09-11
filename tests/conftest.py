from datetime import date

import pytest

from material_balance_studio.domain.models import PVTProperties, PVTTable, ReservoirTank
from material_balance_studio.pvt.table_model import TablePVTModel


@pytest.fixture
def pvt():
    # Analytic Bo = 1.2 + 4e-9 (Pi - p); no Rs change. Synthetic, not field data.
    return TablePVTModel(PVTTable(tuple(
        PVTProperties(float(p), 1.2 + 4e-9 * (30e6 - p), 80.0, bg, 1.0)
        for p, bg in [(10e6, .012), (20e6, .006), (30e6, .004), (40e6, .003)]
    )))


@pytest.fixture
def tank(pvt):
    return ReservoirTank(date(2020, 1, 1), 30e6, 1e6, .2, 0.0, 0.0, 0.0, pvt)


def analytic_oil_production(pressure: float, oil_in_place: float = 1e6) -> float:
    """Independent closed form for the constant-Rs, incompressible-rock benchmark."""
    bo = 1.2 + 4e-9 * (30e6 - pressure)
    return oil_in_place * (bo - 1.2) / bo
