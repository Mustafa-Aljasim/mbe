"""Schilthuis steady influx; constant reference pressure and trapezoidal time integral."""
from dataclasses import dataclass
from typing import ClassVar
from material_balance_studio.domain.validation import finite_value
from .base import AquiferBase


@dataclass(frozen=True)
class SchilthuisAquifer(AquiferBase):
    productivity_index: float  # reservoir m³/(Pa s)
    key: ClassVar[str] = "schilthuis"

    def __post_init__(self):
        finite_value("aquifer productivity index", self.productivity_index, positive=True)

    def compute_step(self, previous, previous_reservoir_pressure, candidate_reservoir_pressure, dt):
        average = self.validate_step(previous, previous_reservoir_pressure, candidate_reservoir_pressure, dt)
        difference = previous.initial_pressure-average
        cumulative = previous.cumulative_influx+self.productivity_index*difference*dt
        return self.step_result(previous, candidate_reservoir_pressure, dt, cumulative, previous.initial_pressure,
                                average, self.productivity_index*(previous.initial_pressure-candidate_reservoir_pressure), difference)
