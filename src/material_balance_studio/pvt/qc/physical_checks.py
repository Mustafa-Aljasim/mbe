"""Deterministic full-range sampling, including both sides of the branch boundary."""
from dataclasses import dataclass
import numpy as np
from ..correlations.registry import REGISTRY
from ..correlations.gas_pseudocritical import PSEUDO_CRITICAL_METHODS
from ..laboratory import LaboratoryData, PROPERTIES
from .validity import applicability


@dataclass(frozen=True)
class QCItem:
    status: str
    check: str
    detail: str


def pressure_grid(model, lab=LaboratoryData()):
    low, high = model.pressure_bounds
    special = [model.pb * factor for factor in (.99, .999, 1-1e-6, 1, 1+1e-6, 1.001, 1.01)]
    return np.array(sorted({*np.linspace(low, high, 1001),
                            *(p for p in special if low <= p <= high),
                            *(p.pressure for p in lab.points if low <= p.pressure <= high)}))


def physical_qc(model, lab=LaboratoryData()) -> tuple[QCItem, ...]:
    grid = pressure_grid(model, lab)
    rows = [QCItem("PASS", "Pressure-range coverage", "Continuous evaluation within the requested range; requests outside it are rejected.")]
    if any(not model.pressure_bounds[0] <= p.pressure <= model.pressure_bounds[1] for p in lab.points):
        rows.append(QCItem("FAIL", "Laboratory coverage", "Laboratory pressures lie outside the requested project range."))
    values = {}
    for name in PROPERTIES:
        try:
            values[name] = np.array([model.property_at_pressure(name, float(p)) for p in grid])
            rows.append(QCItem("PASS", f"{name} positivity", f"Finite, {'nonnegative' if name == 'rs' else 'positive'} at {len(grid)} check pressures."))
        except (ValueError, ArithmeticError, TypeError) as exc:
            rows.append(QCItem("FAIL", f"{name} evaluation", str(exc)))
    for name, lo, hi in (("bo", .7, 4), ("bw", .8, 1.5), ("z", .1, 2), ("rs", 0, 1000)):
        if name in values and (np.min(values[name]) < lo or np.max(values[name]) > hi):
            rows.append(QCItem("CAUTION", f"Unusual {name} magnitude", f"Curve leaves the review envelope {lo:g}–{hi:g} in canonical units; verify fluid and extrapolation. This is a heuristic, not a published validity bound."))
    try:
        published_rsb = model._evaluate("rs", model.pb)
        scale = model.rsb / published_rsb if published_rsb > 0 else float("inf")
        rows.append(QCItem("CAUTION" if abs(scale-1) > .1 else "PASS", "Saturation anchors",
                           f"Untuned Rs curve is normalized to the selected Pb/Rsb anchors by factor {scale:.6g}. Factors differing from 1 by >10% need review."))
    except (ValueError, ArithmeticError, TypeError) as exc:
        rows.append(QCItem("FAIL", "Saturation anchors", str(exc)))
    below, above = grid <= model.pb, grid >= model.pb
    if np.any(above):
        try:
            co = np.array([model.property_at_pressure("co", float(p)) for p in grid[above]])
            good = np.all(np.isfinite(co)) and np.all(co > 0)
            rows.append(QCItem("PASS" if good else "FAIL", "co positivity", "Undersaturated co checked at/above Pb; not applicable below Pb."))
        except (ValueError, ArithmeticError, TypeError) as exc:
            rows.append(QCItem("FAIL", "co evaluation", str(exc)))
    else:
        rows.append(QCItem("CAUTION", "co coverage", "No undersaturated pressures in the project range."))
    for name, mask, increasing in (("rs", below, True), ("bo", below, True), ("bo", above, False),
                                    ("oil_viscosity", below, False), ("oil_viscosity", above, True),
                                    ("bg", np.ones(len(grid), dtype=bool), False),
                                    ("bw", np.ones(len(grid), dtype=bool), False)):
        if name not in values:
            continue
        data = values[name][mask]
        tolerance = max(1e-12, np.max(np.abs(data), initial=0) * 1e-9)
        good = np.all(np.diff(data) * (1 if increasing else -1) >= -tolerance)
        rows.append(QCItem("PASS" if good else "FAIL", f"{name} trend",
                           f"Expected {'increasing' if increasing else 'decreasing'} with pressure on the checked branch; {len(data)} samples."))
    if "rs" in values:
        good = np.all(values["rs"] <= model.rsb * (1+1e-9)) and np.allclose(values["rs"][above], model.rsb, rtol=1e-9)
        rows.append(QCItem("PASS" if good else "FAIL", "Rsb plateau", "Rs must not exceed Rsb; above Pb it must equal Rsb."))
    if model.pressure_bounds[0] < model.pb < model.pressure_bounds[1]:
        for name in ("bo", "rs", "oil_viscosity"):
            try:
                v = [model.property_at_pressure(name, min(model.pressure_bounds[1], max(model.pressure_bounds[0], model.pb*(1+e))))
                     for e in (-1e-6, 0, 1e-6)]
                jump = (max(v)-min(v))/max(abs(v[1]), 1e-12)
                rows.append(QCItem("PASS" if jump <= 1e-4 else "FAIL", f"Pb continuity: {name}", f"Relative change across Pb ±1e-6: {jump:.3g}; tolerance 1e-4."))
            except (ValueError, ArithmeticError, TypeError) as exc:
                rows.append(QCItem("FAIL", f"Pb continuity: {name}", str(exc)))
    else:
        rows.append(QCItem("CAUTION", "Pb continuity", "Pb lies at/outside the range; both branches cannot be checked within coverage."))
    for name, key in vars(model.selection).items():
        result = applicability(REGISTRY[key], model)
        rows.append(QCItem({"VALID": "PASS", "CAUTION": "CAUTION", "NOT_APPLICABLE": "FAIL"}[result.status],
                           f"Applicability: {name}", "; ".join(result.reasons) or REGISTRY[key].name))
    method = PSEUDO_CRITICAL_METHODS[model.pseudo_critical_method]
    try:
        model.pseudo_critical_properties
        valid = method.gas_gravity_range[0] <= model.fluid.gas_gravity <= method.gas_gravity_range[1]
        status = "FAIL" if model.fluid.gas_basis != "sweet" else "PASS" if valid else "CAUTION"
        rows.append(QCItem(status, "Gas pseudo-critical method", f"{method.name}; gas SG range {method.gas_gravity_range}; {method.notes}"))
    except ValueError as exc:
        rows.append(QCItem("FAIL", "Gas pseudo-critical method", str(exc)))
    for warning in lab.warnings:
        rows.append(QCItem("CAUTION", "Laboratory support", warning))
    for name in PROPERTIES:
        measured = lab.measurements(name)
        if measured:
            lo, hi = min(p for p, _ in measured), max(p for p, _ in measured)
            if model.pressure_bounds[0] < lo or model.pressure_bounds[1] > hi:
                rows.append(QCItem("CAUTION", f"Extrapolation: {name}", "Requested range extends beyond laboratory support; QC samples are not validation of extrapolated accuracy."))
    rows.append(QCItem("CAUTION", "Water composition", "McCain Bw does not explicitly use entered salinity or water gravity."))
    return tuple(rows)


def failures(items):
    return tuple(f"{item.check}: {item.detail}" for item in items if item.status == "FAIL")


def qc_status(items):
    statuses = {item.status for item in items}
    return "FAIL" if "FAIL" in statuses else "CAUTION" if "CAUTION" in statuses else "PASS"


def property_qc_status(items, property):
    """Summarize direct checks plus dependencies; preserve all detailed QC rows."""
    dependencies = {"bo": {"rs", "co", "pb"}, "oil_viscosity": {"rs", "pb"},
                    "bg": {"z"}, "rs": {"pb"}}
    names = {property, *dependencies.get(property, set())}
    relevant = []
    for item in items:
        if item.check in ("Laboratory coverage", "Pressure-range coverage"):
            relevant.append(item)
        elif any(item.check.startswith(name + " ") or item.check.endswith(": " + name)
                 or item.check == f"Unusual {name} magnitude"
                 or (item.check == "Laboratory support" and item.detail.startswith(name + ":")) for name in names):
            relevant.append(item)
        elif item.check in ("Rsb plateau", "Saturation anchors") and names & {"rs", "pb"}:
            relevant.append(item)
        elif item.check == "Water composition" and property == "bw":
            relevant.append(item)
        elif item.check == "Gas pseudo-critical method" and property in ("z", "bg", "pseudo_critical"):
            relevant.append(item)
        elif item.check == "Pb continuity" and names & {"bo", "rs", "oil_viscosity", "pb"}:
            relevant.append(item)
    return qc_status(relevant)
