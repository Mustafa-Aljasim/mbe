"""Gas pseudo-critical methods, independent of the z-factor registry.

Sutton, R.P. (1985), Compressibility Factors for High-Molecular-Weight
Reservoir Gases, SPE-14265-MS, doi:10.2118/14265-MS.
Published equations use psia and °R; method results here use Pa and K.
"""
from collections.abc import Callable
from dataclasses import dataclass
from material_balance_studio.domain.validation import finite_value, EngineeringValidationError
from material_balance_studio.units.conversions import to_si


@dataclass(frozen=True)
class PseudoCriticalProperties:
    pressure: float  # Pa
    temperature: float  # K

    def __post_init__(self):
        finite_value("Pseudo-critical pressure", self.pressure, positive=True)
        finite_value("Pseudo-critical temperature", self.temperature, positive=True)


@dataclass(frozen=True)
class PseudoCriticalMethod:
    key: str
    name: str
    evaluate: Callable[[float], PseudoCriticalProperties]
    reference: str
    gas_gravity_range: tuple[float, float]
    notes: str = "Sweet gas; no composition correction."


def sutton_pseudocritical(gas_gravity: float) -> tuple[float, float]:
    """Published field equation: gas SG (air=1) → (Ppc psia, Tpc °R)."""
    # Sutton (1985), Eq. 14: 131.0, not the 131.07 found in some libraries.
    return (756.8 - 131.0 * gas_gravity - 3.6 * gas_gravity ** 2,
            169.2 + 349.5 * gas_gravity - 74 * gas_gravity ** 2)


def _sutton(gas_gravity):
    pc, tc = sutton_pseudocritical(gas_gravity)
    return PseudoCriticalProperties(to_si(pc, "pressure", "FIELD"), tc / 1.8)


PSEUDO_CRITICAL_METHODS = {
    "sutton": PseudoCriticalMethod("sutton", "Sutton", _sutton, "https://doi.org/10.2118/14265-MS", (.57, 1.68)),
}


def register_pseudocritical(method: PseudoCriticalMethod):
    if method.key in PSEUDO_CRITICAL_METHODS:
        raise EngineeringValidationError(f"Duplicate pseudo-critical method: {method.key}")
    PSEUDO_CRITICAL_METHODS[method.key] = method
