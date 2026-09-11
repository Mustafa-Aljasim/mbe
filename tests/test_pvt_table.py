from dataclasses import replace

import pytest

from material_balance_studio.domain.models import PVTProperties, PVTTable
from material_balance_studio.domain.validation import EngineeringValidationError, PVTOutOfRangeError
from material_balance_studio.pvt.table_model import TablePVTModel


def test_exact_nodes_and_linear_interpolation(pvt):
    for row in pvt.table.rows:
        assert pvt.properties_at_pressure(row.pressure) == row
    midpoint = pvt.properties_at_pressure(25e6)
    assert midpoint.bo == pytest.approx(1.22)
    assert midpoint.bg == pytest.approx(.005)
    assert midpoint.rs == 80
    assert midpoint.bw == 1


@pytest.mark.parametrize("pressure", [9e6, 41e6])
def test_no_silent_extrapolation(pvt, pressure):
    with pytest.raises(PVTOutOfRangeError, match="Extrapolation is disabled"):
        pvt.properties_at_pressure(pressure)


@pytest.mark.parametrize("field,value", [("bo", 0), ("bg", -1), ("bw", 0), ("rs", -1),
                                         ("pressure", float("nan")), ("z", 0), ("bwinj", -1)])
def test_bad_properties_rejected(pvt, field, value):
    with pytest.raises(EngineeringValidationError):
        replace(pvt.table.rows[0], **{field: value})


def test_duplicate_unsorted_short_and_incomplete_tables(pvt):
    row, other = pvt.table.rows[:2]
    for rows in [(row,), (row, row), (other, row), (row, replace(other, z=.9))]:
        with pytest.raises(EngineeringValidationError):
            PVTTable(rows)


def test_optional_columns_interpolate():
    rows = (PVTProperties(10e6, 1.2, 50, .01, 1, 1.01, .009, .002, .00002, .8),
            PVTProperties(20e6, 1.3, 70, .005, 1.02, 1.03, .005, .004, .00004, 1))
    actual = TablePVTModel(PVTTable(rows)).properties_at_pressure(15e6)
    assert actual.bwinj == pytest.approx(1.02)
    assert actual.bginj == pytest.approx(.007)
    assert actual.oil_viscosity == pytest.approx(.003)
    assert actual.gas_viscosity == pytest.approx(.00003)
    assert actual.z == pytest.approx(.9)
