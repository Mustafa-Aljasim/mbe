"""Original VEH infinite-acting radial constant-terminal-pressure response."""
from __future__ import annotations

from bisect import bisect_left
from math import exp, log, sqrt, pi

from material_balance_studio.domain.validation import EngineeringValidationError, finite_value

# Dimensionless cumulative water influx WeD for an infinite radial aquifer.
# Values are the classical van Everdingen-Hurst table as reproduced in Ahmed-style
# reservoir engineering texts. Interpolation is in log(tD)-log(WeD) space.
VEH_INFINITE_TABLE: tuple[tuple[float, float], ...] = (
    (0.0, 0.0), (0.01, 0.112), (0.05, 0.276), (0.10, 0.404),
    (0.15, 0.520), (0.20, 0.617), (0.30, 0.796), (0.40, 0.949),
    (0.50, 1.020), (0.60, 1.169), (0.70, 1.275), (0.80, 1.395),
    (1.0, 1.569), (2.0, 2.447), (3.0, 3.202), (4.0, 3.893),
    (5.0, 4.539), (6.0, 5.154), (7.0, 5.765), (8.0, 6.390),
    (9.0, 6.965), (10.0, 7.537), (15.0, 10.235), (20.0, 12.747),
    (25.0, 15.182), (30.0, 17.482), (40.0, 21.928), (50.0, 26.166),
    (60.0, 30.249), (70.0, 34.224), (80.0, 37.966), (90.0, 41.843),
    (100.0, 43.025), (150.0, 56.098), (200.0, 71.691), (300.0, 102.39),
    (400.0, 132.72), (500.0, 162.62), (600.0, 192.30), (700.0, 221.82),
    (800.0, 251.15), (900.0, 280.32), (1000.0, 309.37), (1500.0, 453.37),
    (2000.0, 596.28), (5000.0, 1253.1), (10000.0, 2388.0),
)


def _edwardson_infinite(tD: float) -> float:
    root = sqrt(tD)
    if tD <= 200:
        return (1.12838*root + 1.19328*tD + 0.269872*tD*root + 0.00855294*tD*tD) / (
            1 + 0.616599*root + 0.0413008*tD)
    return (2.02566*tD - 4.29881)/log(tD)


def veh_infinite_water_influx_dimensionless(tD: float) -> float:
    finite_value("VEH dimensionless time", tD, nonnegative=True)
    if tD == 0:
        return 0.0
    if tD < 0.01:
        return 2.0*sqrt(tD/pi)
    times = [row[0] for row in VEH_INFINITE_TABLE]
    index = bisect_left(times, tD)
    if index < len(times) and times[index] == tD:
        return VEH_INFINITE_TABLE[index][1]
    if index == len(times):
        return _edwardson_infinite(tD)
    if index == 0:
        raise EngineeringValidationError("VEH interpolation failed at the table origin.")
    t0, w0 = VEH_INFINITE_TABLE[index-1]
    t1, w1 = VEH_INFINITE_TABLE[index]
    if w0 <= 0 or w1 <= 0:
        return _edwardson_infinite(tD)
    fraction = (log(tD)-log(t0))/(log(t1)-log(t0))
    return exp_interp(w0, w1, fraction)


def exp_interp(y0: float, y1: float, fraction: float) -> float:
    return exp(log(y0) + fraction*(log(y1)-log(y0)))
