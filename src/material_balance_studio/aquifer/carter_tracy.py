"""Carter-Tracy transient radial water influx recurrence."""
from dataclasses import dataclass
from typing import ClassVar

from material_balance_studio.domain.validation import EngineeringValidationError

from .base import AquiferBase
from .carter_tracy_response import (carter_tracy_dimensionless_pressure,
                                    carter_tracy_dimensionless_pressure_derivative)
from .transient import RadialAquiferGeometry, transient_warnings, with_history_diagnostics


@dataclass(frozen=True)
class CarterTracyAquifer(AquiferBase):
    inner_radius: float
    radius_ratio: float
    thickness: float
    porosity: float
    permeability: float
    water_viscosity: float
    total_compressibility: float
    encroachment_angle: float = 360.0
    key: ClassVar[str] = "carter_tracy"

    def __post_init__(self):
        self.geometry

    @property
    def geometry(self):
        return RadialAquiferGeometry(self.inner_radius, self.radius_ratio, self.thickness,
                                     self.porosity, self.permeability, self.water_viscosity,
                                     self.total_compressibility, self.encroachment_angle)

    def compute_step(self, previous, previous_reservoir_pressure, candidate_reservoir_pressure, dt):
        average = self.validate_step(previous, previous_reservoir_pressure, candidate_reservoir_pressure, dt)
        geometry = self.geometry
        previous_tD = geometry.dimensionless_time(previous.elapsed_time)
        current_tD = geometry.dimensionless_time(previous.elapsed_time+dt)
        if current_tD <= previous_tD:
            raise EngineeringValidationError("Carter-Tracy dimensionless time must advance.")
        pressure_drop = previous.initial_pressure-candidate_reservoir_pressure
        pD = carter_tracy_dimensionless_pressure(current_tD)
        pD_prime = carter_tracy_dimensionless_pressure_derivative(current_tD)
        denominator = pD - previous_tD*pD_prime
        if denominator <= 0:
            raise EngineeringValidationError("Carter-Tracy recurrence denominator is nonpositive.")
        term = (geometry.aquifer_constant*pressure_drop - previous.cumulative_influx*pD_prime)/denominator
        cumulative = previous.cumulative_influx + (current_tD-previous_tD)*term
        variables = with_history_diagnostics((), tD=current_tD, previous_tD=previous_tD,
                                             dimensionless_pressure=pD,
                                             dimensionless_pressure_derivative=pD_prime,
                                             recurrence_term=term,
                                             pressure_drop=pressure_drop,
                                             aquifer_constant=geometry.aquifer_constant)
        result = self.step_result(previous, candidate_reservoir_pressure, dt, cumulative,
                                  previous.initial_pressure, average,
                                  driving_difference=previous.initial_pressure-average,
                                  model_variables=variables)
        warnings = transient_warnings(tD=current_tD, response_name="Carter-Tracy pressure response")
        return result.__class__(result.incremental_influx, result.cumulative_influx, result.updated_state,
                                result.average_influx_rate, result.endpoint_influx_rate,
                                result.previous_aquifer_pressure, result.average_reservoir_pressure,
                                result.driving_pressure_difference, tuple((*result.warnings, *warnings)))
