from dataclasses import replace

import pytest

from conftest import analytic_oil_production
from material_balance_studio.domain.models import CumulativeVolumes, PVTProperties, PVTTable, SolverSettings
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.solver.pressure_solver import solve_pressure


def test_a_no_production(tank):
    result = solve_pressure(tank, CumulativeVolumes(), tank.initial_pressure)
    assert result.diagnostics.converged
    assert result.pressure == tank.initial_pressure
    assert result.balance.residual == 0


@pytest.mark.parametrize("target", [10e6, 17.345e6, 20e6, 27e6, 30e6])
def test_b_recover_known_depletion_pressure(tank, target):
    np = analytic_oil_production(target)
    result = solve_pressure(tank, CumulativeVolumes(np=np, gp=80 * np), tank.initial_pressure)
    assert result.diagnostics.converged, result.diagnostics.message
    assert result.pressure == pytest.approx(target, abs=.01)
    assert result.balance.relative_residual <= 1e-8
    assert result.diagnostics.pressure_bounds == (10e6, 30e6)


def test_c_water_injection_reduces_decline(tank):
    np = analytic_oil_production(20e6)
    dry = solve_pressure(tank, CumulativeVolumes(np=np, gp=80 * np), 30e6)
    wet = solve_pressure(tank, CumulativeVolumes(np=np, gp=80 * np, winj=10000), 30e6)
    expected = 30e6 - (np * 1.2 - 10000) / (4e-9 * (1e6 - np))
    assert dry.diagnostics.converged and wet.diagnostics.converged
    assert 20e6 < wet.pressure < 30e6
    assert wet.pressure == pytest.approx(expected, abs=.01)
    assert wet.balance.withdrawal.water_injection == 10000


def test_injection_can_raise_pressure_above_initial(tank):
    result = solve_pressure(tank, CumulativeVolumes(winj=20000), 30e6)
    assert result.diagnostics.converged
    assert result.pressure == pytest.approx(35e6, abs=.01)
    assert result.balance.withdrawal.net < 0
    assert result.diagnostics.pressure_bounds == (10e6, 40e6)


def test_gas_injection_support(tank):
    result = solve_pressure(tank, CumulativeVolumes(ginj=1e6), 30e6)
    assert result.diagnostics.converged
    assert result.pressure > 30e6
    assert result.balance.withdrawal.gas_injection > 0
    assert result.balance.withdrawal.water_injection == 0


def test_unbracketable_pressure_is_controlled_failure(tank):
    result = solve_pressure(tank, CumulativeVolumes(np=900000, gp=72e6), 30e6)
    assert not result.diagnostics.converged
    assert "No pressure root" in result.diagnostics.message
    assert "PVT extrapolation is disabled" in result.diagnostics.message
    assert result.diagnostics.final_residual is not None
    assert result.diagnostics.bracket is None


def test_iteration_limit_and_closure_failure_are_reported(tank):
    model = TablePVTModel(PVTTable((PVTProperties(10e6, 1.25, 20, .012, 1),
                                   PVTProperties(30e6, 1.2, 80, .004, 1))))
    tank = replace(tank, pvt_model=model)
    volumes = CumulativeVolumes(np=150000, gp=12e6)
    failed = solve_pressure(tank, volumes, 30e6, SolverSettings(max_iterations=1))
    assert not failed.diagnostics.converged
    assert failed.diagnostics.iterations == 1
    assert "iteration limit" in failed.diagnostics.message
    passed = solve_pressure(tank, volumes, 30e6)
    assert passed.diagnostics.converged
    strict = solve_pressure(tank, volumes, 30e6, SolverSettings(
        absolute_residual_tolerance=1e-30, relative_residual_tolerance=1e-30))
    assert not strict.diagnostics.converged
    assert "Large material-balance residual" in strict.diagnostics.message


def test_empty_user_bounds_fail_without_crash(tank):
    result = solve_pressure(tank, CumulativeVolumes(np=1), 30e6,
                            SolverSettings(minimum_pressure=35e6))
    assert not result.diagnostics.converged
    assert result.pressure is None
    assert "No admissible" in result.diagnostics.message


def test_combined_gas_cap_compressibility_and_injection_benchmark(tank):
    model = TablePVTModel(PVTTable((
        PVTProperties(10e6, 1.26, 40, .012, 1.02, 1.03, .013),
        PVTProperties(20e6, 1.24, 60, .006, 1.01, 1.02, .007),
        PVTProperties(30e6, 1.20, 80, .004, 1.00, 1.01, .005),
    )))
    tank = replace(tank, pvt_model=model, m=.5, cf=5e-10, cw=4e-10)
    # At 20 MPa: N*Et = 160000 + 300000 + 13050 = 473050 m³.
    # Fnet = 1.48*Np + 1010 - 2040 - 70; solve by hand for Np.
    np = 474150 / 1.48
    result = solve_pressure(tank, CumulativeVolumes(np, 100 * np, 1000, 2000, 10000), 30e6)
    assert result.diagnostics.converged, result.diagnostics.message
    assert result.pressure == pytest.approx(20e6, abs=.01)
    assert result.balance.oil_expansion_support == pytest.approx(160000)
    assert result.balance.gas_cap_expansion_support == pytest.approx(300000)
    assert result.balance.rock_water_expansion_support == pytest.approx(13050)
    assert result.balance.withdrawal.water_injection == pytest.approx(2040)
    assert result.balance.withdrawal.gas_injection == pytest.approx(70)
