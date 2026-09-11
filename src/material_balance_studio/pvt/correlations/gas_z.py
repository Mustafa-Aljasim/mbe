"""Dranchuk–Abou-Kassem (1975) z from reduced pressure and temperature.

Original: Dranchuk and Abou-Kassem (1975), DOI 10.2118/75-03-03.
The equation and coefficients are also reproduced in Sutton (1985), SPE 14265,
equations 12–13. GasCompressibility-py supplies a secondary reference example.
Only supercritical reduced temperature is supported; no sour-gas corrections.
"""

from math import exp

from scipy.optimize import brentq

from material_balance_studio.domain.validation import EngineeringValidationError
# Compatibility import for Phase 2 callers; implementation/selection are separate.
from .gas_pseudocritical import sutton_pseudocritical


def dak_z(pr: float, tr: float) -> float:
    """Solve the DAK density equation. Reference: Pr=3.1995, Tr=1.5006 → .7730934971."""
    if pr <= 0 or tr <= 1:
        raise EngineeringValidationError("DAK requires Pr > 0 and Tr > 1; subcritical gas is not supported.")
    a = .3265 - 1.07 / tr - .5339 / tr ** 3 + .01569 / tr ** 4 - .05165 / tr ** 5
    b = .5475 - .7361 / tr + .1844 / tr ** 2
    c = -.1056 * (-.7361 / tr + .1844 / tr ** 2)

    def eos(rho: float) -> float:
        z = 1 + a * rho + b * rho ** 2 + c * rho ** 5
        z += .6134 * (1 + .7210 * rho ** 2) * rho ** 2 / tr ** 3 * exp(-.7210 * rho ** 2)
        return rho * z - .27 * pr / tr

    rho = brentq(eos, 1e-12, 20.0, xtol=1e-13)
    return .27 * pr / (rho * tr)
