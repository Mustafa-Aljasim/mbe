"""Equilibrated finite pot: We = Caq (Pi - p); Caq = Wi ct [m³/Pa]."""
from dataclasses import dataclass
from typing import ClassVar
from material_balance_studio.domain.validation import finite_value
from .base import AquiferBase


@dataclass(frozen=True)
class PotAquifer(AquiferBase):
    capacity: float
    key: ClassVar[str] = "pot"

    def __post_init__(self):
        finite_value("aquifer capacity (m³/Pa)", self.capacity, positive=True)

    def compute_step(self, previous, previous_reservoir_pressure, candidate_reservoir_pressure, dt):
        average = self.validate_step(previous, previous_reservoir_pressure, candidate_reservoir_pressure, dt)
        cumulative = self.capacity*(previous.initial_pressure-candidate_reservoir_pressure)
        return self.step_result(previous, candidate_reservoir_pressure, dt, cumulative,
                                candidate_reservoir_pressure, average)
