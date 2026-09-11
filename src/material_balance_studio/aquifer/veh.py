"""Original van Everdingen-Hurst infinite radial water influx by pressure-step superposition."""
from dataclasses import dataclass
from typing import ClassVar

from .base import AquiferBase
from .transient import RadialAquiferGeometry, history_entries, transient_warnings, with_history_diagnostics
from .veh_response import veh_infinite_water_influx_dimensionless


@dataclass(frozen=True)
class VanEverdingenHurstAquifer(AquiferBase):
    inner_radius: float
    radius_ratio: float
    thickness: float
    porosity: float
    permeability: float
    water_viscosity: float
    total_compressibility: float
    encroachment_angle: float = 360.0
    key: ClassVar[str] = "van_everdingen_hurst"

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
        elapsed = previous.elapsed_time + dt
        history = history_entries(previous.model_variables)
        delta_pressure = previous_reservoir_pressure - candidate_reservoir_pressure
        if delta_pressure:
            history = (*history, (previous.elapsed_time, delta_pressure))
        contributions = []
        cumulative = 0.0
        warnings = []
        for index, (start_time, dp) in enumerate(history):
            age = elapsed - start_time
            if age < 0:
                raise ValueError("Accepted VEH pressure history is not chronological.")
            tD = geometry.dimensionless_time(age)
            weD = veh_infinite_water_influx_dimensionless(tD)
            contribution = geometry.aquifer_constant*dp*weD
            contributions.append((index, age, tD, contribution))
            cumulative += contribution
            warnings.extend(transient_warnings(tD=tD, response_name="VEH infinite response"))
        total_tD = geometry.dimensionless_time(elapsed)
        current = contributions[-1] if contributions else (-1, 0.0, 0.0, 0.0)
        historical = cumulative-current[3] if contributions else 0.0
        variables = with_history_diagnostics(history, tD=total_tD,
                                             active_pressure_steps=len(history),
                                             current_step_contribution=current[3],
                                             historical_contribution=historical,
                                             response_value=veh_infinite_water_influx_dimensionless(total_tD) if total_tD else 0.0)
        result = self.step_result(previous, candidate_reservoir_pressure, dt, cumulative, previous.initial_pressure,
                                  average, driving_difference=previous.initial_pressure-average,
                                  model_variables=variables)
        return result.__class__(result.incremental_influx, result.cumulative_influx, result.updated_state,
                                result.average_influx_rate, result.endpoint_influx_rate,
                                result.previous_aquifer_pressure, result.average_reservoir_pressure,
                                result.driving_pressure_difference, tuple((*result.warnings, *warnings)))
