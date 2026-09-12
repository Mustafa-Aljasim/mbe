"""Candidate and accepted records are immutable; transfer state is edge-oriented."""
from dataclasses import dataclass
from datetime import datetime
from material_balance_studio.domain.models import CumulativeVolumes,BalanceTerms,PVTProperties
from material_balance_studio.aquifer.state import AquiferState,AquiferStepResult


@dataclass(frozen=True)
class ConnectionState:
    from_tank: str
    to_tank: str
    transmissibility: float
    enabled: bool
    pressure_difference: float
    rate: float
    incremental_transfer: float
    cumulative_transfer: float


@dataclass(frozen=True)
class TankState:
    name: str
    pressure: float
    previous_pressure: float
    cumulative: CumulativeVolumes
    increments: CumulativeVolumes
    pvt: PVTProperties
    components: BalanceTerms
    aquifer_state: AquiferState
    aquifer_step: AquiferStepResult | None
    intertank_support: float
    incremental_support: float
    residual: float
    relative_residual: float
    observed_pressure: float | None
    connection_contributions: tuple
    warnings: tuple = ()


@dataclass(frozen=True)
class NetworkDiagnostics:
    converged: bool
    function_evaluations: int
    optimizer_nfev: int
    initial_guess: tuple
    pressure_bounds: tuple
    residual_norm: float | None
    maximum_tank_relative_residual: float | None
    message: str
    failed_evaluations: int = 0
    dominant_residuals: tuple = ()
    dominant_connections: tuple = ()


@dataclass(frozen=True)
class NetworkState:
    time: datetime
    elapsed_seconds: float
    dt: float
    tanks: tuple[TankState,...]
    connections: tuple[ConnectionState,...]
    signed_network_residual: float
    network_relative_residual: float
    transfer_error: float
    incremental_transfer_error: float
    diagnostics: NetworkDiagnostics | None = None


@dataclass(frozen=True)
class NetworkFailure:
    time: datetime
    diagnostics: NetworkDiagnostics
    candidate: NetworkState | None = None


@dataclass(frozen=True)
class NetworkResult:
    network: object
    settings: object
    initial_state: NetworkState
    states: tuple[NetworkState,...]
    internal_states: tuple[NetworkState,...]
    failed_timestep: NetworkFailure | None = None
    warnings: tuple = ()

    @property
    def converged(self):
        return self.failed_timestep is None
