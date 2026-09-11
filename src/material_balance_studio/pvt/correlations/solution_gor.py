"""Saturated Rs; FIELD equations wrapped by the registry into canonical SI.

Standing: doi:10.2118/947275-G, 18.2/1.4 formulation, as in rNodal oil.Rs.standing.
Vasquez–Beggs: doi:10.2118/6719-PA; Restec Al-Marhoun comparative study, slide 24.
Glasø: doi:10.2118/8016-PA; inverse quadratic as in rNodal Rs.glaso.
Exact equations, gas-gravity basis and links: docs/pvt_phase_2.md.
"""

from math import exp, log10, sqrt

from material_balance_studio.domain.validation import EngineeringValidationError


def standing_rs(p: float, t: float, api: float, gas_gravity: float) -> float:
    """psia, °F, API, SG → scf/STB, saturated branch only."""
    return gas_gravity * ((p / 18.2 + 1.4) * 10 ** (.0125 * api - .00091 * t)) ** (1 / .83)


def vb_coefficients(api: float) -> tuple[float, float, float]:
    return (.0362, 1.0937, 25.7240) if api <= 30 else (.0178, 1.1870, 23.9310)


def vasquez_beggs_rs(p: float, t: float, api: float, gas_gravity: float) -> float:
    """SG is the correlation's gas gravity at the 100-psig separator basis."""
    a, b, c = vb_coefficients(api)
    return a * gas_gravity * p ** b * exp(c * api / (t + 460))


def glaso_rs(p: float, t: float, api: float, gas_gravity: float) -> float:
    """Lower quadratic root on the monotonic saturated branch; no acid-gas correction."""
    discriminant = 1.7447 ** 2 + 4 * .30218 * (1.7669 - log10(p))
    if t <= 0 or discriminant < 0:
        raise EngineeringValidationError("Glaso Rs is outside its real-valued domain.")
    log_x = (1.7447 - sqrt(discriminant)) / (2 * .30218)
    return gas_gravity * (10 ** log_x * api ** .989 / t ** .172) ** (1 / .816)
