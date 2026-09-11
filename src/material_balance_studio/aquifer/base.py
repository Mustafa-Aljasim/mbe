"""Common pure aquifer interface and candidate bookkeeping."""
from dataclasses import asdict, dataclass
from math import isclose
from typing import Protocol
from material_balance_studio.domain.validation import EngineeringValidationError, finite_value
from .state import AquiferState, AquiferStepResult


class AquiferModel(Protocol):
    key: str
    def initial_state(self, reservoir_pressure: float) -> AquiferState: ...
    def compute_step(self, previous: AquiferState, previous_reservoir_pressure: float,
                     candidate_reservoir_pressure: float, dt: float) -> AquiferStepResult: ...


class AquiferBase:
    """Initial aquifer/reservoir equilibrium; signed reversible flow, never clipped."""
    @property
    def parameters(self):
        return tuple(sorted(asdict(self).items()))

    def initial_state(self, reservoir_pressure):
        if reservoir_pressure < 1.:
            raise EngineeringValidationError("Initial aquifer pressure must be at least 1 Pa absolute.")
        return AquiferState(self.key, 0., None if self.key == "none" else reservoir_pressure,
                            reservoir_pressure, 0., reservoir_pressure, self.parameters)

    def validate_step(self, previous, p_previous, p, dt):
        finite_value("aquifer timestep (seconds)", dt, positive=True)
        finite_value("candidate reservoir pressure", p, positive=True)
        finite_value("previous reservoir pressure", p_previous, positive=True)
        if previous.model_key != self.key or previous.model_parameters != self.parameters:
            raise EngineeringValidationError("Aquifer model or parameters changed: start a new simulation; state cannot be reused.")
        if not isclose(previous.previous_reservoir_pressure, p_previous, rel_tol=0, abs_tol=1e-8):
            raise EngineeringValidationError("Previous reservoir pressure disagrees with aquifer state.")
        return (p_previous+p)/2

    def step_result(self, previous, p, dt, cumulative, aquifer_pressure, average_pressure,
                    endpoint_rate=None, driving_difference=None, model_variables=None):
        finite_value("aquifer cumulative influx", cumulative)
        delta = cumulative-previous.cumulative_influx
        finite_value("incremental aquifer influx", delta)
        finite_value("aquifer average rate", delta/dt)
        if endpoint_rate is not None:
            finite_value("aquifer endpoint rate", endpoint_rate)
        warnings = []
        if delta < 0:
            warnings.append("Negative aquifer influx (efflux): cumulative We decreases under the signed reversible-flow policy; verify pressure rebound.")
        if cumulative < 0:
            warnings.append("Net water transfer is into the aquifer (negative cumulative We); verify injection history and flow direction.")
        if endpoint_rate is not None and endpoint_rate < 0 and delta >= 0:
            warnings.append("Endpoint aquifer rate indicates efflux although interval-average influx is nonnegative; review pressure recovery and timestep size.")
        if aquifer_pressure is not None and aquifer_pressure < 1.:
            raise EngineeringValidationError("Aquifer pressure is below the permitted 1 Pa absolute lower bound.")
        updated = AquiferState(self.key, cumulative, aquifer_pressure, p,
                               previous.elapsed_time+dt, previous.initial_pressure,
                               self.parameters, previous.model_variables if model_variables is None else model_variables)
        return AquiferStepResult(delta, cumulative, updated, delta/dt, endpoint_rate,
                                 previous.aquifer_pressure, average_pressure, driving_difference, tuple(warnings))


@dataclass(frozen=True)
class AquiferContext:
    model: AquiferModel
    previous_state: AquiferState
    previous_reservoir_pressure: float
    dt: float

    def evaluate(self, candidate_pressure):
        return self.model.compute_step(self.previous_state, self.previous_reservoir_pressure,
                                       candidate_pressure, self.dt)
