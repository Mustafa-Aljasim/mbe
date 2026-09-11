"""Engineering presentation units, separate from the canonical file-import contract.

SI uploads still contain Pa and 1/Pa. SI presentation uses MPa and 1/MPa.
These helpers create converted values; they never change domain objects.
"""

from .conversions import UnitSystem, conversion_factor

# Name: (canonical/file conversion quantity, SI display label, FIELD label).
_QUANTITIES = {
    "pressure": ("pressure", "MPa", "psia"),
    "oil_volume": ("liquid_volume", "stock-tank m³", "STB"),
    "water_volume": ("liquid_volume", "surface m³", "bbl"),
    "gas_volume": ("gas_volume", "standard m³", "scf"),
    "reservoir_volume": ("liquid_volume", "reservoir m³", "rb"),
    "bo": ("bo", "m³/m³", "rb/STB"),
    "bw": ("bw", "m³/m³", "rb/STB"),
    "bwinj": ("bwinj", "m³/m³", "rb/STB"),
    "bg": ("bg", "m³/m³", "rb/scf"),
    "bginj": ("bginj", "m³/m³", "rb/scf"),
    "rs": ("rs", "m³/m³", "scf/STB"),
    "expansion": ("bo", "m³/m³", "rb/STB"),
    "compressibility": ("compressibility", "1/MPa", "1/psi"),
    "viscosity": ("viscosity", "Pa·s", "cP"),
    "length": ("length", "m", "ft"),
    "permeability": ("permeability", "m²", "mD"),
    "angle": ("angle", "degrees", "degrees"),
    "dimensionless": ("dimensionless", "dimensionless", "dimensionless"),
    "aquifer_capacity": ("dimensionless", "reservoir m³/MPa", "rb/psi"),
    "aquifer_productivity": ("dimensionless", "reservoir m³/(day·MPa)", "rb/(day·psi)"),
    "aquifer_rate": ("dimensionless", "reservoir m³/day", "rb/day"),
}


def display_unit(quantity: str, units: UnitSystem | str) -> str:
    """Label matching the returned display values, including surface/reservoir basis."""
    _, si_label, field_label = _QUANTITIES[quantity]
    return si_label if UnitSystem(units) == UnitSystem.SI else field_label


def _display_factor(quantity: str, units: UnitSystem | str) -> float:
    if quantity in ("aquifer_capacity", "aquifer_productivity", "aquifer_rate"):
        pressure_factor = 1e6 if UnitSystem(units) == UnitSystem.SI else conversion_factor("pressure", "FIELD")
        volume_factor = 1. if UnitSystem(units) == UnitSystem.SI else 1/conversion_factor("liquid_volume", "FIELD")
        return volume_factor * (1 if quantity == "aquifer_rate" else pressure_factor) * (1 if quantity == "aquifer_capacity" else 86400)
    canonical_quantity = _QUANTITIES[quantity][0]
    if UnitSystem(units) == UnitSystem.FIELD:
        return 1.0 / conversion_factor(canonical_quantity, UnitSystem.FIELD)
    return {"pressure": 1e-6, "compressibility": 1e6}.get(quantity, 1.0)


def to_display(value: float | None, quantity: str, units: UnitSystem | str) -> float | None:
    """Convert canonical SI to presentation, preserving missing observations."""
    factor = _display_factor(quantity, units)
    return None if value is None else value * factor


def from_display(value: float, quantity: str, units: UnitSystem | str) -> float:
    """Convert an engineering setup value back to the canonical SI boundary."""
    return value / _display_factor(quantity, units)


def column_label(label: str, quantity: str, units: UnitSystem | str) -> str:
    return f"{label} ({display_unit(quantity, units)})"


def format_value(value: object) -> str:
    """Readable inspector values without rounding the underlying numeric datasets."""
    if value is None:
        return "—"
    if isinstance(value, float):
        if value != 0 and (abs(value) < 1e-3 or abs(value) >= 1e9):
            return f"{value:.6e}"
        return f"{value:,.8g}"
    return str(value)
