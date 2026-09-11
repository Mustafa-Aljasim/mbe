"""Absolute temperature conversions. Canonical temperature is kelvin."""

from material_balance_studio.domain.validation import finite_value
from .conversions import UnitSystem


def temperature_to_kelvin(value: float, units: UnitSystem | str) -> float:
    """Engineering SI input is °C; FIELD is °F (not °R)."""
    finite_value("temperature", value)
    kelvin = value + 273.15 if UnitSystem(units) == UnitSystem.SI else (value - 32) * 5 / 9 + 273.15
    finite_value("absolute temperature", kelvin, positive=True)
    return kelvin


def temperature_from_kelvin(value: float, units: UnitSystem | str) -> float:
    finite_value("absolute temperature", value, positive=True)
    return value - 273.15 if UnitSystem(units) == UnitSystem.SI else (value - 273.15) * 9 / 5 + 32
