from dataclasses import dataclass
from typing import ClassVar
from .base import AquiferBase


@dataclass(frozen=True)
class NoAquifer(AquiferBase):
    key: ClassVar[str] = "none"

    def compute_step(self, previous, previous_reservoir_pressure, candidate_reservoir_pressure, dt):
        average = self.validate_step(previous, previous_reservoir_pressure, candidate_reservoir_pressure, dt)
        if previous.cumulative_influx != 0:
            raise ValueError("NoAquifer requires zero previous influx.")
        return self.step_result(previous, candidate_reservoir_pressure, dt, 0., None, average)
