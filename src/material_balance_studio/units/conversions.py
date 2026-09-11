"""SI/field conversion factors, including the distinct field gas-volume basis.

Field convention: psia, STB liquids, scf gas, rb/STB Bo/Bw,
rb/scf Bg/Bginj, scf/STB Rs, 1/psi compressibility, cP viscosity.
SI: Pa absolute, m³ at corresponding reference conditions, 1/Pa, Pa.s.
"""

from enum import StrEnum

from material_balance_studio.domain.validation import EngineeringValidationError

PA_PER_PSI = 6894.757293168
M3_PER_STB = 0.158987294928
M3_PER_SCF = 0.028316846592
M_PER_FT = 0.3048
M2_PER_MD = 9.869233e-16


class UnitSystem(StrEnum):
    SI = "SI"
    FIELD = "FIELD"


_FIELD_TO_SI = {
    "pressure": PA_PER_PSI,
    "liquid_volume": M3_PER_STB,
    "gas_volume": M3_PER_SCF,
    "compressibility": 1.0 / PA_PER_PSI,
    "bo": 1.0,
    "bw": 1.0,
    "bwinj": 1.0,
    "bg": M3_PER_STB / M3_PER_SCF,
    "bginj": M3_PER_STB / M3_PER_SCF,
    "rs": M3_PER_SCF / M3_PER_STB,
    "viscosity": 0.001,
    "length": M_PER_FT,
    "permeability": M2_PER_MD,
    "dimensionless": 1.0,
    "angle": 1.0,
}


def conversion_factor(quantity: str, units: UnitSystem | str) -> float:
    """Return a factor to SI, rejecting unknown quantity/unit names."""
    if quantity not in _FIELD_TO_SI:
        raise EngineeringValidationError(f"Unknown unit quantity: {quantity}.")
    try:
        system = UnitSystem(units)
    except ValueError as exc:
        raise EngineeringValidationError(f"Unknown unit system: {units}.") from exc
    return 1.0 if system == UnitSystem.SI else _FIELD_TO_SI[quantity]


def to_si(value: float, quantity: str, units: UnitSystem | str) -> float:
    return value * conversion_factor(quantity, units)


def from_si(value: float, quantity: str, units: UnitSystem | str) -> float:
    return value / conversion_factor(quantity, units)
