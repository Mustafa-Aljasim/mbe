"""Carter-Tracy infinite radial constant-terminal-rate pressure response."""
from __future__ import annotations

from math import log, sqrt

from material_balance_studio.domain.validation import finite_value


def carter_tracy_dimensionless_pressure(tD: float) -> float:
    finite_value("Carter-Tracy dimensionless time", tD, positive=True)
    if tD <= 100.0:
        root = sqrt(tD)
        return (370.528*root + 137.582*tD + 5.69549*tD*root) / (
            328.834 + 265.488*root + 45.2157*tD + tD*root)
    return 0.5*(log(tD) + 0.80907)


def carter_tracy_dimensionless_pressure_derivative(tD: float) -> float:
    finite_value("Carter-Tracy dimensionless time", tD, positive=True)
    if tD > 100.0:
        return 0.5/tD
    step = max(1e-7, tD*1e-5)
    left = max(tD-step, tD*0.5)
    right = tD+step
    return (carter_tracy_dimensionless_pressure(right)-carter_tracy_dimensionless_pressure(left))/(right-left)
