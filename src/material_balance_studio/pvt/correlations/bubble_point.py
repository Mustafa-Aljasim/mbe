"""Published bubble pressures in psia from Rsb in scf/STB and temperature °F.

Standing (1947), Vasquez–Beggs (1980), Glasø (1980); see docs/pvt_phase_2.md.
Vasquez–Beggs is inverted from its published saturated Rs equation without
rounding its inverse coefficients.
"""

from math import exp, log10

from .solution_gor import vb_coefficients


def standing_pb(rsb: float, t: float, api: float, gas_gravity: float) -> float:
    return 18.2 * ((rsb / gas_gravity) ** .83 * 10 ** (.00091 * t - .0125 * api) - 1.4)


def vasquez_beggs_pb(rsb: float, t: float, api: float, gas_gravity: float) -> float:
    a, b, c = vb_coefficients(api)
    return (rsb / (a * gas_gravity * exp(c * api / (t + 460)))) ** (1 / b)


def glaso_pb(rsb: float, t: float, api: float, gas_gravity: float) -> float:
    x = log10((rsb / gas_gravity) ** .816 * t ** .172 / api ** .989)
    return 10 ** (1.7669 + 1.7447 * x - .30218 * x * x)
