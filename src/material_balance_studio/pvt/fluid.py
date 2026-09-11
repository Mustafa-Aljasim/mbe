"""Canonical fluid definition; saturation information may initially be incomplete."""

from dataclasses import dataclass

from material_balance_studio.domain.validation import EngineeringValidationError, finite_value
from material_balance_studio.units.conversions import UnitSystem, to_si
from material_balance_studio.units.temperature import temperature_to_kelvin


@dataclass(frozen=True)
class BlackOilFluid:
    api: float
    gas_gravity: float
    temperature: float  # K
    initial_pressure: float  # Pa absolute
    pb: float | None = None  # Pa absolute, measured
    rsb: float | None = None  # standard gas m³ / stock-tank oil m³, measured
    oil_gravity: float | None = None
    water_gravity: float = 1.0
    salinity: float = 0.0  # mass fraction; stored, not used by the McCain Bw approximation
    gas_basis: str = "sweet"  # no acid-gas correction in this release

    def __post_init__(self) -> None:
        for key in ("api", "gas_gravity", "temperature", "initial_pressure", "water_gravity"):
            finite_value(key, getattr(self, key), positive=True)
        for key in ("pb", "rsb", "oil_gravity"):
            if getattr(self, key) is not None:
                finite_value(key, getattr(self, key), positive=True)
        finite_value("salinity mass fraction", self.salinity, nonnegative=True)
        if self.salinity > 0.35:
            raise EngineeringValidationError("Salinity must be a mass fraction between 0 and 0.35, not ppm.")
        if self.oil_gravity is not None and abs(self.oil_gravity - self.sg_oil) > .005:
            raise EngineeringValidationError("Oil specific gravity conflicts with API gravity (tolerance 0.005 SG).")
        if self.gas_basis not in ("sweet", "sour"):
            raise EngineeringValidationError("Gas basis must be sweet or sour.")

    @property
    def sg_oil(self) -> float:
        return 141.5 / (self.api + 131.5)


def fluid_from_inputs(*, api: float, gas_gravity: float, temperature: float,
                      initial_pressure: float, units: UnitSystem | str = UnitSystem.SI,
                      pb: float | None = None, rsb: float | None = None,
                      **optional) -> BlackOilFluid:
    """File/API boundary: SI pressures Pa and temperature °C; FIELD psia and °F."""
    return BlackOilFluid(api, gas_gravity, temperature_to_kelvin(temperature, units),
                         to_si(initial_pressure, "pressure", units),
                         None if pb is None else to_si(pb, "pressure", units),
                         None if rsb is None else to_si(rsb, "rs", units), **optional)
