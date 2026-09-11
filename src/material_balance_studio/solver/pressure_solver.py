"""Brent pressure solve constrained to measured PVT coverage."""

import numpy as np
from scipy.optimize import brentq

from material_balance_studio.domain.models import (
    BalanceTerms, CumulativeVolumes, PressureSolution, ReservoirTank,
    SolverDiagnostics, SolverSettings,
)
from material_balance_studio.domain.validation import EngineeringValidationError, finite_value
from material_balance_studio.mbe.oil_balance import evaluate_balance
from material_balance_studio.aquifer.base import AquiferContext


def closure_passes(balance: BalanceTerms, settings: SolverSettings) -> bool:
    """Pressure convergence alone is insufficient; both volume QC limits must pass."""
    return (balance.absolute_residual <= settings.absolute_residual_tolerance
            and balance.relative_residual <= settings.relative_residual_tolerance)


def solve_pressure(tank: ReservoirTank, cumulative: CumulativeVolumes,
                   previous_pressure: float, settings: SolverSettings | None = None,
                   *, aquifer_context: AquiferContext | None = None) -> PressureSolution:
    """Find a root near previous pressure, returning controlled failure diagnostics.

    Scan each PVT interval before selecting a sign-changing bracket. Injection
    permits a search above Pi, only as far as the table/user maximum allows.
    """
    settings = settings or SolverSettings()
    finite_value("previous_pressure", previous_pressure, positive=True)
    pvt_low, pvt_high = tank.pvt_model.pressure_bounds
    low = max(pvt_low, settings.minimum_pressure)
    has_injection = cumulative.winj > 0 or cumulative.ginj > 0
    high = min(pvt_high, pvt_high if has_injection else tank.initial_pressure)
    if settings.maximum_pressure is not None:
        high = min(high, settings.maximum_pressure)
    guess = min(max(previous_pressure, low), high) if low <= high else previous_pressure
    calls = 0
    warnings: list[str] = []
    cache: dict[float, BalanceTerms] = {}
    aquifer_steps = {}

    def balance_at(pressure: float) -> BalanceTerms:
        nonlocal calls
        if pressure not in cache:
            calls += 1
            step = aquifer_context.evaluate(pressure) if aquifer_context is not None else None
            aquifer_steps[pressure] = step
            cache[pressure] = evaluate_balance(pressure, tank, cumulative, settings.normalization_floor,
                                               aquifer_influx=step.cumulative_influx if step else 0.)
        return cache[pressure]

    def result(pressure: float | None, iterations: int, bracket: tuple[float, float] | None,
               message: str, *, force_failure: bool = False) -> PressureSolution:
        balance = None if pressure is None else balance_at(pressure)
        converged = balance is not None and closure_passes(balance, settings) and not force_failure
        if balance is not None and not closure_passes(balance, settings):
            message += (f" Large material-balance residual: |R|={balance.absolute_residual:.6g} m³, "
                        f"relative={balance.relative_residual:.6g}; limits are "
                        f"{settings.absolute_residual_tolerance:g} m³ and "
                        f"{settings.relative_residual_tolerance:g}.")
        return PressureSolution(pressure, balance, SolverDiagnostics(
            converged=converged, iterations=iterations, function_calls=calls,
            initial_guess=guess, pressure_bounds=(low, high), bracket=bracket,
            final_residual=None if balance is None else balance.residual,
            message=message, warnings=tuple(warnings),
        ), aquifer_steps.get(pressure) if converged else None)

    if low > high:
        return result(None, 0, None, "No admissible pressure interval. Check PVT and solver bounds.")
    try:
        if aquifer_context is not None:
            if aquifer_context.previous_reservoir_pressure != previous_pressure:
                raise EngineeringValidationError("Aquifer context disagrees with solver previous pressure.")
            if tank.aquifer is not None and tank.aquifer != aquifer_context.model:
                raise EngineeringValidationError("Aquifer context model differs from tank configuration.")
        if tank.aquifer is not None and tank.aquifer.key != "none" and aquifer_context is None:
            raise EngineeringValidationError("An active aquifer requires prior state and a positive timestep; use simulate or supply aquifer_context.")
        if closure_passes(balance_at(guess), settings):
            return result(guess, 0, (guess, guess), "Previous-pressure guess satisfies balance closure.")
        knots = sorted({low, high, guess, *(p for p in tank.pvt_model.pressure_nodes if low <= p <= high)})
        samples = set(knots)
        for left, right in zip(knots, knots[1:]):
            samples.update(float(p) for p in np.linspace(left, right, settings.bracket_subdivisions + 1))
        pressures = sorted(samples)
        residuals = [balance_at(p).residual for p in pressures]
        brackets = [(p, p) for p, residual in zip(pressures, residuals) if residual == 0.0]
        brackets.extend((left, right) for left, right, fl, fr in
                        zip(pressures, pressures[1:], residuals, residuals[1:])
                        if (fl < 0 < fr) or (fr < 0 < fl))
        if not brackets:
            best = min(pressures, key=lambda p: abs(balance_at(p).residual))
            if closure_passes(balance_at(best), settings):
                return result(best, 0, (best, best), "Sampled pressure satisfies both closure tolerances.")
            return result(best, 0, None,
                          f"No pressure root could be bracketed in [{low:g}, {high:g}] Pa. "
                          f"Residuals at bounds: {residuals[0]:.6g}, {residuals[-1]:.6g} m³. "
                          "Check cumulative history, N, PVT coverage and injection. "
                          "PVT extrapolation is disabled.", force_failure=True)
        if len(brackets) > 1:
            warnings.append("Multiple root brackets found; selected the bracket nearest previous pressure. "
                            "Review PVT consistency and the physical pressure branch.")
        bracket = min(brackets, key=lambda pair: min(abs(pair[0] - guess), abs(pair[1] - guess)))
        if bracket[0] == bracket[1]:
            return result(bracket[0], 0, bracket, "Exact root at a sampled pressure.")
        root, info = brentq(lambda p: balance_at(p).residual, *bracket,
                            xtol=settings.pressure_tolerance, rtol=1e-14,
                            maxiter=settings.max_iterations, full_output=True, disp=False)
        return result(float(root), info.iterations, bracket,
                      "Brent solver converged." if info.converged else "Brent solver reached its iteration limit.",
                      force_failure=not info.converged)
    except (EngineeringValidationError, ValueError, RuntimeError, OverflowError) as exc:
        # A failed step never becomes a successful reservoir state.
        return PressureSolution(None, None, SolverDiagnostics(
            False, 0, calls, guess, (low, high), None, None,
            f"Pressure solve failed: {exc}", tuple(warnings),
        ))
