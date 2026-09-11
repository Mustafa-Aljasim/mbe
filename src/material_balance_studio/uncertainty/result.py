"""Immutable, serializable analysis records in canonical physical units."""
from dataclasses import dataclass,asdict
import hashlib,json


def scenario_identity(scenario):
    payload=dict(base=asdict(scenario.base_tank),parameters=scenario.fitted_parameters,
                 specifications=[asdict(s) for s in scenario.specifications],observations=[asdict(o) for o in scenario.observations],
                 pvt=scenario.pvt_fingerprint,mode=scenario.weighting_mode,default_sigma=scenario.default_sigma,
                 history=[dict(date=str(s.date),cumulative=asdict(s.cumulative),observed_pressure=s.observed_pressure)
                          for s in scenario.final_simulation.states] if scenario.final_simulation is not None else None)
    return hashlib.sha256(json.dumps(payload,sort_keys=True,default=str).encode()).hexdigest()


def matrix(array):
    return tuple(tuple(float(v) for v in row) for row in array)


@dataclass(frozen=True)
class IdentifiabilityResult:
    scenario_identity: str
    parameters: tuple
    dates: tuple
    match_mode: str
    objective_definition: str
    sensitivity: tuple
    scaled_sensitivity: tuple
    residual_jacobian: tuple
    scaled_jacobian: tuple
    schemes: tuple
    derivative_errors: tuple
    singular_values: tuple
    rank: int
    condition: float | None
    weak_combination: tuple
    covariance: tuple | None
    standard_errors: tuple
    intervals: tuple
    correlation: tuple | None
    column_coupling: tuple
    statuses: tuple
    warnings: tuple
    failed_evaluations: int
    evaluations: int
    statistical_assumptions: bool
    profiles: tuple = ()
    surfaces: tuple = ()


@dataclass(frozen=True)
class GridPoint:
    parameters: tuple
    objective: float | None
    rmse: float | None
    valid: bool
    optimizer_status: str
    bound_flags: tuple
    failure: str | None
    evaluations: int


@dataclass(frozen=True)
class ProfileResult:
    scenario_identity: str
    parameter: str
    mode: str
    points: tuple
    threshold: float | None
    threshold_method: str
    plausible_ranges: tuple
    warnings: tuple
    failed_evaluations: int
    threshold_brackets: tuple = ()


@dataclass(frozen=True)
class SurfaceResult:
    scenario_identity: str
    parameters: tuple
    x: tuple
    y: tuple
    points: tuple
    objective_definition: str
    initial: tuple
    matched: tuple
    warnings: tuple
    failed_evaluations: int
