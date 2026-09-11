"""Registry contract: every equation has property, reference, units and bounds."""

from collections.abc import Callable
from dataclasses import dataclass

from material_balance_studio.domain.validation import EngineeringValidationError
from material_balance_studio.units.conversions import from_si
from material_balance_studio.units.temperature import temperature_from_kelvin
from ..fluid import BlackOilFluid


@dataclass(frozen=True)
class CorrelationInput:
    fluid: BlackOilFluid
    pressure: float
    rs: float | None = None
    z: float | None = None
    pseudo_critical_method: str = "sutton"

    @property
    def reduced_state(self) -> tuple[float, float]:
        """Resolve the independently selected pseudo-critical method before z."""
        from .gas_pseudocritical import PSEUDO_CRITICAL_METHODS
        pc = PSEUDO_CRITICAL_METHODS[self.pseudo_critical_method].evaluate(self.fluid.gas_gravity)
        return self.pressure / pc.pressure, self.fluid.temperature / pc.temperature

    @property
    def p(self) -> float:
        return from_si(self.pressure, "pressure", "FIELD")

    @property
    def t(self) -> float:
        return temperature_from_kelvin(self.fluid.temperature, "FIELD")

    @property
    def r(self) -> float:
        if self.rs is None:
            raise EngineeringValidationError("This correlation requires Rs/Rsb.")
        return from_si(self.rs, "rs", "FIELD")


@dataclass(frozen=True)
class ValidityRange:
    key: str
    minimum: float
    maximum: float
    unit: str


@dataclass(frozen=True)
class Correlation:
    key: str
    property: str
    name: str
    evaluate: Callable[[CorrelationInput], float]
    reference: str
    formulation: str
    ranges: tuple[ValidityRange, ...] = ()
    notes: str = ""
