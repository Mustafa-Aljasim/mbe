from dataclasses import dataclass, replace
from ..correlations.registry import candidates
from ..correlations.base import CorrelationInput
from ..qc.validity import applicability, Applicability
from ..qc.physical_checks import physical_qc, failures
from .metrics import error_metrics, Metrics


@dataclass(frozen=True)
class ScreeningResult:
    property: str
    correlation: str
    applicability: Applicability
    pressures: tuple[float, ...]
    measured: tuple[float, ...]
    predicted: tuple[float, ...]
    metrics: Metrics | None
    errors: tuple[str, ...] = ()


def observations(model, lab, property):
    if property == "pb":
        return () if model.fluid.pb is None else ((model.fluid.initial_pressure, model.fluid.pb),)
    return lab.measurements(property)


def screen_property(model, lab, property) -> tuple[ScreeningResult, ...]:
    """Keep every candidate; rank applicability tier, then SI RMSE within property.

    Pb screening uses published predictions against scalar measured Pb, not the
    anchored final model. All other predictions include continuous branch logic.
    Unmeasured properties have no numerical recommendation.
    """
    data = observations(model, lab, property)
    pressures, measured = tuple(p for p, _ in data), tuple(v for _, v in data)
    results = []
    for candidate in candidates(property):
        predicted, metric, errors = (), None, ()
        try:
            alternative = model.with_selection(property, candidate.key)
            if property == "pb" and model.fluid.rsb is not None:
                alternative = replace(alternative, use_measured_pb=False)
            validity = applicability(candidate, alternative)  # requested range, not only lab span
            if validity.status != "NOT_APPLICABLE" and data:
                if property == "pb":
                    if model.fluid.rsb is None:
                        raise ValueError("Pb screening needs measured Rsb independent of the measured Pb target.")
                    predicted = (candidate.evaluate(CorrelationInput(model.fluid, pressures[0], model.fluid.rsb)),)
                    if predicted[0] <= 0:
                        raise ValueError("Predicted Pb is nonpositive.")
                else:
                    predicted = tuple(alternative.property_at_pressure(property, p) for p in pressures)
                metric = error_metrics(measured, predicted)
                if property != "pb":
                    relevant = tuple(item for item in physical_qc(alternative, lab)
                        if item.check.startswith(property + " ") or item.check == f"Pb continuity: {property}"
                        or (property == "rs" and item.check == "Rsb plateau"))
                    problems = failures(relevant)
                    if problems:
                        validity = Applicability("NOT_APPLICABLE", validity.reasons + problems)
        except (ValueError, ArithmeticError, TypeError) as exc:
            validity = Applicability("NOT_APPLICABLE", (str(exc),))
            errors = (str(exc),)
        results.append(ScreeningResult(property, candidate.key, validity, pressures, measured, predicted, metric, errors))
    order = {"VALID": 0, "CAUTION": 1, "NOT_APPLICABLE": 2}
    return tuple(sorted(results, key=lambda r: (order[r.applicability.status], r.metrics.rmse if r.metrics else float("inf"), r.correlation)))


def recommendation(results):
    return next((r.correlation for r in results if r.metrics is not None and r.applicability.status != "NOT_APPLICABLE"), None)
