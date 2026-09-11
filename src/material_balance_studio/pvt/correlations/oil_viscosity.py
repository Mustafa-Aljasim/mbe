"""Beggs–Robinson (1975) dead/live oil with Vasquez–Beggs above Pb.

Pengtools documented SG variant, SPE-5434-PA and SPE-6719-PA. °F, scf/STB → cP.
The SI registry boundary converts cP to Pa.s. docs/pvt_phase_2.md gives sources.
"""

from math import exp


def beggs_robinson_dead(t: float, oil_gravity: float) -> float:
    return 10 ** (t ** -1.163 * exp(13.108 - 6.591 / oil_gravity)) - 1


def beggs_robinson_live(rs: float, t: float, oil_gravity: float) -> float:
    return 10.715 * (rs + 100) ** -.515 * beggs_robinson_dead(t, oil_gravity) ** (5.44 * (rs + 150) ** -.338)


def undersaturated_viscosity(mu_b: float, p: float, pb: float) -> float:
    """Vasquez–Beggs multiplier; p and Pb psia, viscosity units unchanged."""
    exponent = 2.6 * p ** 1.187 * exp(-11.513 - 8.98e-5 * p)
    return mu_b * (p / pb) ** exponent
