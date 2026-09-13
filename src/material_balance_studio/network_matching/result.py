"""Complete scenario snapshots; applying parameters is a separate explicit action."""
from dataclasses import dataclass


@dataclass(frozen=True)
class NetworkMatchResult:
    scenario_name: str
    scenario_id: str
    base_network: object
    fitted_network: object | None
    history_plan: object
    history_fingerprint: str
    observations: tuple
    specifications: tuple
    settings: object
    mode: str
    default_sigma: float | None
    fitted_parameters: tuple
    converged: bool
    optimizer_status: int
    termination_message: str
    function_evaluations: int
    optimizer_nfev: int
    initial_objective: float
    final_objective: float
    initial_simulation: object | None
    final_simulation: object | None
    initial_pressure_rows: tuple
    final_pressure_rows: tuple
    initial_tank_metrics: tuple
    final_tank_metrics: tuple
    initial_network_metrics: tuple
    final_network_metrics: tuple
    bound_flags: tuple
    evaluations: tuple
    maximum_valid_relative_residual: float | None
    maximum_valid_absolute_residual: float | None
    maximum_transfer_conservation_error: float | None
    minimum_closure_margin: float | None
    failed_candidates: int
    qc_status: str
    warnings: tuple
    allocation_audit: tuple
