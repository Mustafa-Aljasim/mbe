"""Immutable scenario envelope; applying it is a separate user action."""
from dataclasses import dataclass,asdict,is_dataclass
import hashlib,json


def pvt_fingerprint(model):
    bounds=model.pressure_bounds
    payload={"type":type(model).__module__+"."+type(model).__qualname__,
             "configuration":asdict(model) if is_dataclass(model) else None,
             "samples":[asdict(model.properties_at_pressure(bounds[0]+(bounds[1]-bounds[0])*i/20)) for i in range(21)]}
    return hashlib.sha256(json.dumps(payload,sort_keys=True,default=str).encode()).hexdigest()


@dataclass(frozen=True)
class HistoryMatchResult:
    base_tank: object
    fitted_tank: object | None
    aquifer_model: str
    pvt_fingerprint: str
    specifications: tuple
    observations: tuple
    weighting_mode: str
    default_sigma: float | None
    fitted_parameters: tuple
    optimizer_status: int
    converged: bool
    termination_message: str
    function_evaluations: int
    optimizer_nfev: int
    initial_objective: float
    final_objective: float
    initial_metrics: dict
    final_metrics: dict
    initial_simulation: object | None
    final_simulation: object | None
    pressure_series: tuple
    evaluations: tuple
    maximum_mbe_relative: float
    bound_flags: tuple
    qc_status: str
    warnings: tuple
    diagnostics: dict | None
    ho_ooip: float | None
    ho_difference_percent: float | None
