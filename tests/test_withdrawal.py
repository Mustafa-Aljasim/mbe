from dataclasses import replace

import pytest

from material_balance_studio.domain.models import CumulativeVolumes, PVTProperties
from material_balance_studio.mbe.withdrawal import produced_withdrawal, withdrawal_terms


def test_hand_calculated_withdrawal_and_separate_injection():
    pvt = PVTProperties(20e6, 1.2, 50, .005, 1.02, bwinj=1.01, bginj=.006)
    actual = withdrawal_terms(CumulativeVolumes(100, 8000, 20, 10, 1000), pvt)
    assert actual.production == pytest.approx(120 + 15 + 20.4)
    assert actual.water_injection == pytest.approx(10.1)
    assert actual.gas_injection == pytest.approx(6)
    assert actual.net == pytest.approx(139.3)


def test_zero_oil_gas_only_and_water_only_are_safe(pvt):
    props = pvt.properties_at_pressure(20e6)
    assert produced_withdrawal(CumulativeVolumes(), props) == 0
    assert produced_withdrawal(CumulativeVolumes(gp=1000), props) == pytest.approx(6)
    assert produced_withdrawal(CumulativeVolumes(wp=100), props) == 100
    assert withdrawal_terms(CumulativeVolumes(winj=10, ginj=1000), props).net == -16


def test_solution_gas_correction_is_not_clipped(pvt):
    props = pvt.properties_at_pressure(20e6)
    actual = produced_withdrawal(CumulativeVolumes(np=100, gp=7000), props)
    assert actual == pytest.approx(124 - 6)
