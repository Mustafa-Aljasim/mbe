"""Independent immutable runs, diagnostic pressure errors and engineering review."""
from dataclasses import dataclass, replace
from math import sqrt, isfinite

from material_balance_studio.solver.simulation import simulate


@dataclass(frozen=True)
class ComparisonRun:
    name: str
    model: object
    result: object | None
    error: str | None = None


def compare_aquifers(tank, history, models, settings=None):
    """Observations are shared QC targets, never forced into predicted pressure."""
    history = tuple(history)
    runs = []
    for name, model in models.items():
        try:
            runs.append(ComparisonRun(name, model, simulate(replace(tank, aquifer=model), history, settings)))
        except ValueError as exc:
            runs.append(ComparisonRun(name, model, None, str(exc)))
    return tuple(runs)


def pressure_metrics(result):
    errors = [s.pressure_error for s in result.states if s.pressure_error is not None] if result else []
    if not errors:
        return dict(count=0, rmse=None, mae=None, bias=None, maximum=None)
    return dict(count=len(errors), rmse=sqrt(sum(e*e for e in errors)/len(errors)),
                mae=sum(abs(e) for e in errors)/len(errors), bias=sum(errors)/len(errors),
                maximum=max(abs(e) for e in errors))


def support_fraction(state):
    """We / net production voidage; undefined at nonpositive or <1 m3 basis."""
    basis = state.balance.withdrawal.net
    return state.balance.aquifer_support/basis if basis >= 1.0 else None


def engineering_qc(run):
    """Screening thresholds are review triggers, not universal reservoir limits."""
    if run.result is None:
        return "FAIL", (run.error or "Run unavailable.",)
    result, model = run.result, run.model
    failures, notes = [], list(result.warnings)
    if not result.converged:
        failures.append(f"Pressure solve failed at {result.failed_timestep.date}; partial results are not ranked.")
    for s in result.states:
        if not isfinite(s.balance.relative_residual) or s.balance.relative_residual > 1e-8:
            failures.append("Material-balance closure exceeds 1e-8.")
        if s.aquifer_state.aquifer_pressure is not None and s.aquifer_state.aquifer_pressure < 1:
            failures.append("Aquifer pressure is invalid.")
        if s.aquifer_step and s.pressure <= s.previous_pressure and s.aquifer_step.incremental_influx < -1e-8:
            notes.append("Decreasing We during a depletion step; inspect prior rebound and signed flow.")
        fraction = support_fraction(s)
        if fraction is None:
            notes.append("Support fraction undefined: net withdrawal is below 1 reservoir m³.")
        elif fraction > 1.05 or fraction < 0:
            notes.append("Aquifer support lies outside 0–105% of net withdrawal; review injection/rebound and parameters.")
        notes.extend(s.warnings)
        if s.aquifer_step:
            notes.extend(s.aquifer_step.warnings)
        td = dict(s.aquifer_state.model_variables).get("diag_tD", 0)
        if 0 < td < .01 or td > 10000:
            notes.append("Transient response outside tabulated range; approximation used.")
    if hasattr(model, "geometry"):
        notes.append("Infinite-acting radial assumption requires geological justification; radius ratio does not impose a boundary.")
        if not 1e-16 <= model.permeability <= 1e-11:
            notes.append("Permeability outside screening range 0.1–10000 mD (approximately).")
    if hasattr(model, "initial_water_volume") and model.initial_water_volume <= 0:
        failures.append("Invalid connected water volume.")
    metrics = pressure_metrics(result)
    if metrics["count"] == 0:
        notes.append("No observed pressure: numerical fit is unavailable.")
    elif metrics["rmse"] > .05*result.initial_state.pressure:
        notes.append("Pressure RMSE exceeds 5% of initial pressure; review supplied parameters and model assumptions.")
    notes = tuple(dict.fromkeys((*failures, *notes)))
    return ("FAIL" if failures else "CAUTION" if notes else "PASS"), notes


SENSITIVITY_PARAMETERS = {
    "fetkovich": ("productivity_index", "initial_water_volume"),
    "carter_tracy": ("permeability", "encroachment_angle"),
    "van_everdingen_hurst": ("permeability", "encroachment_angle"),
    "modified_van_everdingen_hurst": ("permeability", "encroachment_angle"),
}


def sensitivity_runs(tank, history, model, parameter, values):
    if parameter not in SENSITIVITY_PARAMETERS.get(model.key, ()):
        raise ValueError("This parameter is not supported for sensitivity preview.")
    values = tuple(values)
    if not 2 <= len(values) <= 8 or len(set(values)) != len(values):
        raise ValueError("Supply 2–8 distinct parameter values.")
    models = {f"{parameter}={value:.8g} SI": replace(model, **{parameter: value}) for value in values}
    return compare_aquifers(tank, history, models)
