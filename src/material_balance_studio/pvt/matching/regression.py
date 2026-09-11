"""Least squares with branch-preserving constraints and full-range acceptance QC."""
from dataclasses import dataclass, replace
import numpy as np
from ..models.correlation_model import Transform
from ..qc.physical_checks import physical_qc, failures, pressure_grid, qc_status
from .metrics import Metrics, error_metrics
from .ranking import observations


@dataclass(frozen=True)
class FitRecord:
    property: str
    correlation: str
    a: float
    b: float
    n: int
    formulation: str
    before: Metrics
    after: Metrics
    maximum_full_range_deviation: float | None = None  # |matched − baseline|, canonical property units
    maximum_full_range_deviation_percent: float | None = None
    maximum_deviation_pressure: float | None = None  # Pa
    qc_before: str = "NOT_CHECKED"
    qc_after: str = "NOT_CHECKED"


@dataclass(frozen=True)
class MatchResult:
    accepted: bool
    record: FitRecord | None
    model: object | None
    reasons: tuple[str, ...] = ()


def fit_property(model, lab, property) -> MatchResult:
    """Fits keep other property transforms but replace this property's transform.

    Rs: anchored slope; Bg/Pb: scale only; others: affine. z and Bg cannot both
    be fitted. The input model remains immutable. A failed attempt never returns
    an active model. Scalar measured Pb is used only in predicted-Pb mode.
    """
    record = None
    try:
        baseline = replace(model, transforms=tuple(t for t in model.transforms if t.property != property))
        before_qc = physical_qc(baseline, lab)
        data = observations(baseline, lab, property)
        if not data:
            raise ValueError(f"No measured {property} data.")
        if property == "pb" and baseline.use_measured_pb and baseline.fluid.pb is not None:
            raise ValueError("Measured Pb is already an anchor; select predicted Pb to match it.")
        x = np.array([baseline.property_at_pressure(property, p) for p, _ in data])
        y = np.array([v for _, v in data])
        if property == "rs":
            dx = x-baseline.rsb
            if np.dot(dx, dx) <= 1e-24:
                raise ValueError("Rs matching needs a point below Pb; plateau-only data cannot identify a slope.")
            b = float(np.dot(dx, y-baseline.rsb)/np.dot(dx, dx))
            a = baseline.rsb*(1-b)
            formulation = "Rsb + B (Rs_raw − Rsb); A=Rsb(1−B)"
        elif property in ("bg", "pb"):
            a, b = 0.0, float(np.dot(x, y)/np.dot(x, x))
            formulation = "B × raw (A=0); Bg scale also scales z" if property == "bg" else "B × predicted Pb (A=0)"
        else:
            if len(x) < 2 or np.ptp(x) <= max(1e-15, np.max(np.abs(x))*1e-10):
                raise ValueError("Affine matching requires at least two distinct predicted values.")
            # Centered regression avoids conditioning problems for small Bg/viscosity.
            dx, dy = x-x.mean(), y-y.mean()
            b = float(np.dot(dx, dy)/np.dot(dx, dx))
            a = float(y.mean()-b*x.mean())
            formulation = "A + B × raw"
        # Preserve a reviewable candidate report even when it fails full-range
        # physics. Evaluate the affine arithmetic before guarded model creation.
        grid = pressure_grid(baseline, lab)
        grid_raw = np.array([baseline.property_at_pressure(property, float(p)) for p in grid])
        deviation = np.abs(a + b*grid_raw - grid_raw)
        index = int(np.argmax(deviation))
        relative = None if np.any(grid_raw == 0) else float(np.max(100*deviation/np.abs(grid_raw)))
        record = FitRecord(property, getattr(model.selection, property), a, b, len(x), formulation,
                           error_metrics(y, x), error_metrics(y, a+b*x), float(deviation[index]), relative,
                           float(grid[index]), qc_status(before_qc), "FAIL")
        transformed = baseline.with_transform(Transform(property, a, b))
        after_qc = physical_qc(transformed, lab)
        record = replace(record, qc_after=qc_status(after_qc))
        problems = failures(after_qc)
        if problems:
            return MatchResult(False, record, None, problems)
        if record.after.rmse > record.before.rmse + max(1e-14, record.before.rmse*1e-9):
            return MatchResult(False, record, None, ("Matching increased RMSE.",))
        return MatchResult(True, record, transformed)
    except (ValueError, ArithmeticError, TypeError) as exc:
        return MatchResult(False, record, None, (str(exc),))
