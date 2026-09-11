"""Continuous black-oil branches and independent property selections; no MBE code."""
from dataclasses import dataclass, replace
from functools import cached_property, lru_cache
from math import isfinite
from material_balance_studio.domain.models import PVTProperties
from material_balance_studio.domain.validation import EngineeringValidationError, PVTOutOfRangeError
from material_balance_studio.units.conversions import from_si
from ..fluid import BlackOilFluid
from ..correlations.base import CorrelationInput
from ..correlations.registry import REGISTRY
from ..correlations.gas_pseudocritical import PSEUDO_CRITICAL_METHODS
from ..correlations.oil_viscosity import undersaturated_viscosity


@dataclass(frozen=True)
class PropertySelection:
    pb: str = "standing_pb"
    rs: str = "standing_rs"
    bo: str = "standing_bo"
    oil_viscosity: str = "beggs_robinson_mu"
    z: str = "sutton_dak_z"
    bg: str = "real_gas_bg"
    bw: str = "mccain_bw"
    co: str = "vasquez_beggs_co"


@dataclass(frozen=True)
class Transform:
    property: str
    a: float
    b: float

    def __post_init__(self):
        if self.property not in ("bo", "rs", "oil_viscosity", "z", "bg", "bw", "pb"):
            raise EngineeringValidationError("Unsupported matching property.")
        if not isfinite(self.a) or not isfinite(self.b) or self.b <= 0:
            raise EngineeringValidationError("Matching requires finite A and positive B to preserve trends.")


@dataclass(frozen=True)
class CorrelationPVTModel:
    fluid: BlackOilFluid
    pressure_bounds: tuple[float, float]
    selection: PropertySelection = PropertySelection()
    transforms: tuple[Transform, ...] = ()
    use_measured_pb: bool = True
    pseudo_critical_method: str = "sutton"

    def __post_init__(self):
        if self.pseudo_critical_method not in PSEUDO_CRITICAL_METHODS:
            raise EngineeringValidationError(f"Unknown gas pseudo-critical method: {self.pseudo_critical_method}")
        low, high = self.pressure_bounds
        if not all(isfinite(v) for v in (low, high)) or not 0 < low < high:
            raise EngineeringValidationError("PVT range requires 0 < minimum < maximum pressure.")
        if not low <= self.fluid.initial_pressure <= high:
            raise EngineeringValidationError("Initial pressure is outside PVT range.")
        for property, key in vars(self.selection).items():
            if key not in REGISTRY or REGISTRY[key].property != property:
                raise EngineeringValidationError(f"Invalid {property} selection: {key}")
        names = [t.property for t in self.transforms]
        if len(set(names)) != len(names) or {"bg", "z"} <= set(names):
            raise EngineeringValidationError("Use one transform per property; matching Bg and z simultaneously is ambiguous.")
        for t in self.transforms:
            if t.property in ("bg", "pb") and t.a != 0:
                raise EngineeringValidationError(f"{t.property} matching is scale-only (A=0).")
            if t.property == "pb" and self.use_measured_pb and self.fluid.pb is not None:
                raise EngineeringValidationError("Measured Pb is a fixed anchor; select predicted Pb before fitting it.")
        if self.fluid.rsb is None and (self.fluid.pb is None or not self.use_measured_pb):
            raise EngineeringValidationError("Provide measured Pb or Rsb; API, gas gravity, temperature and Pi cannot determine both.")
        if self.pb <= 0 or not isfinite(self.pb) or self.rsb <= 0 or not isfinite(self.rsb):
            raise EngineeringValidationError("Computed Pb and Rsb must be finite and positive.")
        for t in self.transforms:
            if t.property == "rs" and abs(t.a - self.rsb * (1-t.b)) > 1e-9 * max(1, self.rsb):
                raise EngineeringValidationError("Rs matching must preserve the Rsb anchor: A=Rsb(1−B).")

    def _evaluate(self, property, p, rs=None, z=None):
        if self.fluid.gas_basis != "sweet" and property in ("z", "bg"):
            raise EngineeringValidationError("Sour gas is not applicable: composition correction is not implemented.")
        return REGISTRY[getattr(self.selection, property)].evaluate(CorrelationInput(self.fluid, p, rs, z, self.pseudo_critical_method))

    @cached_property
    def pseudo_critical_properties(self):
        return PSEUDO_CRITICAL_METHODS[self.pseudo_critical_method].evaluate(self.fluid.gas_gravity)

    def _transform(self, property, value):
        t = next((t for t in self.transforms if t.property == property), None)
        return value if t is None else t.a + t.b * value

    @cached_property
    def rsb(self):
        return self.fluid.rsb if self.fluid.rsb is not None else self._evaluate("rs", self.fluid.pb)

    @cached_property
    def pb(self):
        if self.use_measured_pb and self.fluid.pb is not None:
            return self.fluid.pb
        return self._transform("pb", self._evaluate("pb", self.fluid.initial_pressure, self.rsb))

    @property
    def predicted_pb(self):
        """Independent untuned Pb estimate, unavailable without measured Rsb."""
        if self.fluid.rsb is None:
            return None
        return self._evaluate("pb", self.fluid.initial_pressure, self.fluid.rsb)

    @property
    def pressure_nodes(self):
        low, high = self.pressure_bounds
        return tuple(sorted({low, high, *([self.pb] if low < self.pb < high else [])}))

    def with_selection(self, property, key):
        return replace(self, selection=replace(self.selection, **{property: key}), transforms=())

    def with_transform(self, transform):
        return replace(self, transforms=tuple(t for t in self.transforms if t.property != transform.property) + (transform,))

    @lru_cache(maxsize=8192)
    def property_at_pressure(self, property: str, pressure: float) -> float:
        if not isfinite(pressure) or not self.pressure_bounds[0] <= pressure <= self.pressure_bounds[1]:
            raise PVTOutOfRangeError(f"Pressure {pressure:g} Pa is outside PVT range {self.pressure_bounds}.")
        p, pb = pressure, self.pb
        if property == "pb":
            return pb
        if property == "rs":
            raw = self.rsb if p >= pb else self.rsb * self._evaluate("rs", p) / self._evaluate("rs", pb)
            value = self._transform("rs", raw)
        elif property in ("bo", "oil_viscosity"):
            # At/below Pb: live-oil equation. Above: continuous branch anchored at Pb.
            rs = self.property_at_pressure("rs", p) if p < pb else self.rsb
            value = self._evaluate(property, min(p, pb), rs)
            if p > pb:
                if property == "bo":
                    k = self._evaluate("co", pb, self.rsb) * pb
                    if k <= 0:
                        raise EngineeringValidationError("Nonpositive undersaturated oil compressibility.")
                    value *= (p / pb) ** -k
                else:
                    value = undersaturated_viscosity(value, from_si(p, "pressure", "FIELD"), from_si(pb, "pressure", "FIELD"))
            value = self._transform(property, value)
        elif property == "z":
            value = self._transform("z", self._evaluate("z", p))
            # Bg scale matching scales z too, preserving Bg = Z*T*Psc/(P*Tsc).
            value = self._transform("bg", value)
        elif property == "bg":
            value = self._evaluate("bg", p, z=self.property_at_pressure("z", p))
        elif property == "bw":
            value = self._transform("bw", self._evaluate("bw", p))
        elif property == "co":
            if p < pb:
                raise EngineeringValidationError("Undersaturated oil compressibility is not applicable below Pb.")
            # Above Pb, report −d ln(Bo_matched)/dp consistently with affine Bo.
            value = self._evaluate("co", p, self.rsb)
            transform = next((t for t in self.transforms if t.property == "bo"), None)
            if transform is not None:
                bo = self.property_at_pressure("bo", p)
                value *= (bo - transform.a) / bo
        else:
            raise EngineeringValidationError(f"Unknown PVT property: {property}")
        if not isfinite(value) or value < 0 or (property != "rs" and value == 0):
            raise EngineeringValidationError(f"Nonphysical {property} at {p:g} Pa.")
        return value

    def properties_at_pressure(self, pressure: float) -> PVTProperties:
        return PVTProperties(pressure=pressure, **{k: self.property_at_pressure(k, pressure)
                             for k in ("bo", "rs", "bg", "bw", "oil_viscosity", "z")})
