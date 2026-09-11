"""Immutable candidate/committed records: Pa absolute, reservoir m³, seconds."""
from dataclasses import dataclass
from material_balance_studio.domain.validation import finite_value


@dataclass(frozen=True)
class AquiferState:
    model_key: str
    cumulative_influx: float
    aquifer_pressure: float | None
    previous_reservoir_pressure: float
    elapsed_time: float
    initial_pressure: float
    model_parameters: tuple[tuple[str, float], ...] = ()
    model_variables: tuple[tuple[str, float], ...] = ()

    def __post_init__(self):
        finite_value("cumulative aquifer influx", self.cumulative_influx)
        finite_value("aquifer elapsed time", self.elapsed_time, nonnegative=True)
        for name in ("previous_reservoir_pressure", "initial_pressure"):
            finite_value(name, getattr(self, name), positive=True)
        if self.aquifer_pressure is not None:
            finite_value("aquifer pressure", self.aquifer_pressure, positive=True)
        for field in ("model_parameters", "model_variables"):
            entries = tuple((str(k), float(v)) for k, v in getattr(self, field))
            for key, value in entries:
                finite_value(key, value)
            object.__setattr__(self, field, entries)


@dataclass(frozen=True)
class AquiferStepResult:
    incremental_influx: float
    cumulative_influx: float
    updated_state: AquiferState
    average_influx_rate: float  # delta We / dt, not the endpoint Darcy rate
    endpoint_influx_rate: float | None
    previous_aquifer_pressure: float | None
    average_reservoir_pressure: float
    driving_pressure_difference: float | None
    warnings: tuple[str, ...] = ()
