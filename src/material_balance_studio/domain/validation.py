"""Engineering validation; invalid inputs are rejected, never repaired silently."""

from math import isfinite
from collections.abc import Sequence
from datetime import date
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .models import HistoryRecord


class EngineeringValidationError(ValueError):
    """An input violates the supported engineering model."""


class PVTOutOfRangeError(EngineeringValidationError):
    """Pressure is outside measured PVT coverage; extrapolation is disabled."""


def finite_value(name: str, value: float, *, positive: bool = False,
                 nonnegative: bool = False) -> None:
    """Require a finite engineering scalar with optional sign constraints."""
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not isfinite(value):
        raise EngineeringValidationError(f"{name} must be a finite number.")
    if positive and value <= 0:
        raise EngineeringValidationError(f"{name} must be > 0; received {value}.")
    if nonnegative and value < 0:
        raise EngineeringValidationError(f"{name} must be >= 0; received {value}.")


def validate_history(history: Sequence["HistoryRecord"], initial_date: date | None = None) -> None:
    """Require sorted unique dates and nondecreasing totals from the initial reference."""
    from .models import CUMULATIVE_FIELDS, CumulativeVolumes

    if not history:
        raise EngineeringValidationError("History requires at least one record.")
    previous = CumulativeVolumes()
    previous_date = None
    for record in history:
        if previous_date is not None and record.date <= previous_date:
            raise EngineeringValidationError("History must be sorted chronologically with unique dates.")
        if initial_date is not None and record.date < initial_date:
            raise EngineeringValidationError("History cannot precede reservoir initial_date.")
        for name in CUMULATIVE_FIELDS:
            if getattr(record.cumulative, name) < getattr(previous, name):
                raise EngineeringValidationError(f"Cumulative {name} decreases on {record.date}.")
            if record.date == initial_date and getattr(record.cumulative, name) != 0:
                raise EngineeringValidationError("Cumulative production/injection must be zero at initial_date.")
        previous, previous_date = record.cumulative, record.date
