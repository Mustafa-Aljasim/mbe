"""Shared radial transient-aquifer calculations in canonical SI units."""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, pi

from material_balance_studio.domain.validation import EngineeringValidationError, finite_value


@dataclass(frozen=True)
class RadialAquiferGeometry:
    inner_radius: float  # reservoir/aquifer contact radius, m
    radius_ratio: float  # outer aquifer radius / inner radius
    thickness: float  # m
    porosity: float
    permeability: float  # m2
    water_viscosity: float  # Pa.s
    total_compressibility: float  # 1/Pa
    encroachment_angle: float = 360.0  # degrees

    def __post_init__(self):
        for name in ("inner_radius", "radius_ratio", "thickness", "permeability",
                     "water_viscosity", "total_compressibility"):
            finite_value("aquifer "+name, getattr(self, name), positive=True)
        finite_value("aquifer porosity", self.porosity, positive=True)
        finite_value("aquifer encroachment angle", self.encroachment_angle, positive=True)
        if self.porosity >= 1:
            raise EngineeringValidationError("Aquifer porosity must be less than 1.")
        if self.radius_ratio <= 1:
            raise EngineeringValidationError("Aquifer radius ratio must be greater than 1.")
        if self.encroachment_angle > 360:
            raise EngineeringValidationError("Aquifer encroachment angle cannot exceed 360 degrees.")

    @property
    def outer_radius(self):
        return self.inner_radius*self.radius_ratio

    @property
    def encroachment_radians(self):
        return 2*pi*self.encroachment_angle/360.0

    @property
    def aquifer_constant(self):
        """Reservoir volume per pressure drop, m3/Pa, for radial edge-water influx."""
        return self.encroachment_radians*self.porosity*self.total_compressibility*self.thickness*self.inner_radius**2

    def dimensionless_time(self, elapsed_seconds: float) -> float:
        finite_value("aquifer elapsed time", elapsed_seconds, nonnegative=True)
        return self.permeability*elapsed_seconds/(self.porosity*self.water_viscosity*
                                                  self.total_compressibility*self.inner_radius**2)


def transient_warnings(*, tD: float, response_name: str) -> tuple[str, ...]:
    if not isfinite(tD) or tD < 0:
        raise EngineeringValidationError("Invalid aquifer dimensionless time.")
    warnings = []
    if tD > 10000:
        warnings.append(f"{response_name} evaluated beyond the tabulated benchmark range; large-time approximation was used.")
    if 0 < tD < 1e-8:
        warnings.append(f"{response_name} evaluated at very small dimensionless time; early-time approximation was used.")
    return tuple(warnings)


def history_entries(variables: tuple[tuple[str, float], ...]) -> tuple[tuple[float, float], ...]:
    pairs = []
    index = 0
    values = dict(variables)
    while f"step_time_{index}" in values:
        pairs.append((values[f"step_time_{index}"], values[f"step_dp_{index}"]))
        index += 1
    return tuple(pairs)


def with_history_diagnostics(history: tuple[tuple[float, float], ...], **diagnostics: float) -> tuple[tuple[str, float], ...]:
    rows = []
    for index, (time_seconds, delta_pressure) in enumerate(history):
        rows.append((f"step_time_{index}", float(time_seconds)))
        rows.append((f"step_dp_{index}", float(delta_pressure)))
    rows.extend((f"diag_{key}", float(value)) for key, value in diagnostics.items() if value is not None)
    return tuple(rows)


def diagnostics_from_variables(variables: tuple[tuple[str, float], ...]) -> dict[str, float]:
    return {name[5:]: value for name, value in variables if name.startswith("diag_")}
