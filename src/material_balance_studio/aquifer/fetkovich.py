"""Fetkovich (1971), SPE 2603, Eqs 6–9: finite capacity, interval-average boundary pressure."""
from dataclasses import dataclass
from math import expm1, isclose
from typing import ClassVar
from material_balance_studio.domain.validation import EngineeringValidationError, finite_value
from .base import AquiferBase


@dataclass(frozen=True)
class FetkovichAquifer(AquiferBase):
    initial_water_volume: float  # initial reservoir m³, already connected fraction
    total_compressibility: float  # water + aquifer pore compressibility, 1/Pa
    productivity_index: float  # reservoir m³/(Pa s)
    key: ClassVar[str] = "fetkovich"

    def __post_init__(self):
        for name in ("initial_water_volume", "total_compressibility", "productivity_index", "capacity"):
            finite_value("aquifer "+name, getattr(self, name), positive=True)

    @property
    def capacity(self):
        return self.initial_water_volume*self.total_compressibility

    def compute_step(self, previous, previous_reservoir_pressure, candidate_reservoir_pressure, dt):
        average = self.validate_step(previous, previous_reservoir_pressure, candidate_reservoir_pressure, dt)
        prior_pa = previous.initial_pressure-previous.cumulative_influx/self.capacity
        if previous.aquifer_pressure is None or not isclose(prior_pa, previous.aquifer_pressure, rel_tol=1e-12, abs_tol=1e-6):
            raise EngineeringValidationError("Aquifer pressure and cumulative influx are inconsistent with finite capacity.")
        difference = prior_pa-average
        fraction = -expm1(-self.productivity_index*dt/self.capacity)
        cumulative = previous.cumulative_influx+self.capacity*difference*fraction
        current_pa = previous.initial_pressure-cumulative/self.capacity
        return self.step_result(previous, candidate_reservoir_pressure, dt, cumulative, current_pa,
                                average, self.productivity_index*(current_pa-candidate_reservoir_pressure), difference)
