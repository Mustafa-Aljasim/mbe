"""Saturated Bo correlations (rb/STB = m³/m³).

Standing uses the documented 0.972/1.47e-4/1.175 formulation, not the alternate
0.9759/0.00012/1.2 refit. Glasø uses 0.968*T. References and exact variants in
docs/pvt_phase_2.md. Inputs: Rs scf/STB, temperature °F, gravities dimensionless.
"""

from math import log10

from material_balance_studio.domain.validation import EngineeringValidationError
from .base import CorrelationInput


def standing_bo(rs: float, t: float, oil_gravity: float, gas_gravity: float) -> float:
    return .972 + 1.47e-4 * (rs * (gas_gravity / oil_gravity) ** .5 + 1.25 * t) ** 1.175


def vasquez_beggs_bo(rs: float, t: float, api: float, gas_gravity: float) -> float:
    a, b, c = (4.677e-4, 1.751e-5, -1.811e-8) if api <= 30 else (4.670e-4, 1.100e-5, 1.337e-9)
    return 1 + a * rs + (t - 60) * api / gas_gravity * (b + c * rs)


def glaso_bo(rs: float, t: float, oil_gravity: float, gas_gravity: float) -> float:
    x = log10(rs * (gas_gravity / oil_gravity) ** .526 + .968 * t)
    return 1 + 10 ** (-6.58511 + 2.91329 * x - .27683 * x * x)


def vasquez_beggs_compressibility(inputs: CorrelationInput) -> float:
    """Published SI reformulation: Rs m³/m³, T K, p MPa → co 1/Pa.

    Pengtools / Afanasyev et al. (2004): numerator/(1e5*p_MPa) gives 1/MPa.
    The model integrates co=K/p exactly, giving Bo=Bob*(p/Pb)^(-K).
    """
    if inputs.rs is None:
        raise EngineeringValidationError("Oil compressibility requires Rsb.")
    f = inputs.fluid
    numerator = 28.1 * inputs.rs + 30.6 * f.temperature - 1180 * f.gas_gravity + 1784 / f.sg_oil - 10910
    return numerator / (1e5 * inputs.pressure)
