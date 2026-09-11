from dataclasses import replace

import pytest

from material_balance_studio.domain.models import CumulativeVolumes, PVTProperties
from material_balance_studio.mbe.expansion import expansion_terms
from material_balance_studio.mbe.oil_balance import balance_residual, evaluate_balance


def test_expansion_terms_against_hand_values(tank):
    tank = replace(tank, m=.5, cf=5e-10, cw=4e-10)
    initial = PVTProperties(30e6, 1.2, 80, .004, 1)
    current = PVTProperties(20e6, 1.24, 60, .006, 1)
    terms = expansion_terms(current, initial, tank)
    assert terms.eo == pytest.approx(.16)
    assert terms.eg == pytest.approx(.6)
    assert terms.meg == pytest.approx(.3)
    # 1.2 * 1.5 * (5e-10 + .2*4e-10) / .8 * 10e6
    assert terms.efw == pytest.approx(.01305)
    assert terms.et == pytest.approx(.47305)


def test_zero_expansion_at_initial_pressure(tank):
    initial = tank.pvt_model.properties_at_pressure(tank.initial_pressure)
    terms = expansion_terms(initial, initial, replace(tank, m=2.0))
    assert (terms.eo, terms.eg, terms.meg, terms.efw, terms.et) == (0, 0, 0, 0, 0)


def test_independent_residual_and_normalization(tank):
    volumes = CumulativeVolumes(np=10000, gp=800000, winj=1000)
    balance = evaluate_balance(20e6, tank, volumes)
    assert balance.total_expansion_support == pytest.approx(40000)
    assert balance.withdrawal.net == pytest.approx(11400)
    assert balance_residual(20e6, tank, volumes) == pytest.approx(-28600)
    assert balance.absolute_residual == pytest.approx(28600)
    assert balance.relative_residual == pytest.approx(28600 / 11400)
    assert evaluate_balance(30e6, tank, CumulativeVolumes()).relative_residual == 0


def test_support_above_pi_has_correct_sign(tank):
    tank = replace(tank, m=.2, cf=5e-10, cw=4e-10)
    terms = evaluate_balance(35e6, tank, CumulativeVolumes())
    assert terms.expansion.eo < 0
    assert terms.expansion.meg < 0
    assert terms.expansion.efw < 0
