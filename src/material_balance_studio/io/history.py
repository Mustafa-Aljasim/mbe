"""CSV/XLSX inputs with strict schemas, coherent unit conversion and chronology."""

from datetime import date
from pathlib import Path
from typing import BinaryIO

import pandas as pd

from material_balance_studio.domain.models import (
    CUMULATIVE_FIELDS, CumulativeVolumes, HistoryRecord, PVTProperties, PVTTable,
)
from material_balance_studio.domain.validation import EngineeringValidationError, validate_history
from material_balance_studio.units.conversions import UnitSystem, to_si


def read_table(source: str | Path | BinaryIO, filename: str | None = None) -> pd.DataFrame:
    """Read CSV or the first sheet of an XLSX workbook; headers are exact lowercase names."""
    suffix = Path(filename or (str(source) if isinstance(source, (str, Path)) else source.name)).suffix.lower()
    try:
        if suffix == ".csv":
            return pd.read_csv(source)
        if suffix == ".xlsx":
            return pd.read_excel(source, engine="openpyxl")
    except (ValueError, OSError, ImportError) as exc:
        raise EngineeringValidationError(f"Could not read input table: {exc}") from exc
    raise EngineeringValidationError("Supported input formats are .csv and .xlsx (first sheet).")


def _check_columns(frame: pd.DataFrame, required: set[str], optional: set[str]) -> None:
    if frame.empty:
        raise EngineeringValidationError("Input table is empty.")
    if frame.columns.duplicated().any():
        raise EngineeringValidationError("Duplicate column names are not allowed.")
    missing = required - set(frame.columns)
    unknown = set(frame.columns) - required - optional
    if missing or unknown:
        raise EngineeringValidationError(f"Invalid columns. Missing: {sorted(missing)}; "
                                         f"unrecognized: {sorted(unknown)}. Use documented lowercase headers.")


def _number(value: object, column: str, row_number: int) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise EngineeringValidationError(f"Row {row_number}, {column}: expected a number.") from exc


def _date(value: object, row_number: int) -> date:
    if isinstance(value, str):
        try:
            return date.fromisoformat(value)
        except ValueError as exc:
            raise EngineeringValidationError(f"Row {row_number}: date must be YYYY-MM-DD.") from exc
    if isinstance(value, (date, pd.Timestamp)) and not pd.isna(value):
        if isinstance(value, pd.Timestamp) and value != value.normalize():
            raise EngineeringValidationError(f"Row {row_number}: use calendar dates without times.")
        return value if type(value) is date else value.date()
    raise EngineeringValidationError(f"Row {row_number}: missing/invalid calendar date.")


def history_from_frame(frame: pd.DataFrame, units: UnitSystem | str = UnitSystem.SI) -> tuple[HistoryRecord, ...]:
    """Missing injection columns mean zero injection; blank observations mean unobserved."""
    _check_columns(frame, {"date", "np", "gp", "wp"}, {"winj", "ginj", "observed_pressure"})
    records = []
    for row_number, row in enumerate(frame.to_dict("records"), start=2):
        values = {
            name: to_si(_number(row.get(name, 0.0), name, row_number),
                        "gas_volume" if name in ("gp", "ginj") else "liquid_volume", units)
            for name in CUMULATIVE_FIELDS
        }
        observed = row.get("observed_pressure")
        pressure = None if pd.isna(observed) else to_si(
            _number(observed, "observed_pressure", row_number), "pressure", units)
        records.append(HistoryRecord(_date(row["date"], row_number), CumulativeVolumes(**values), pressure))
    validate_history(records)
    return tuple(records)


def pvt_from_frame(frame: pd.DataFrame, units: UnitSystem | str = UnitSystem.SI) -> PVTTable:
    """Parse mandatory PVT and optional injection/viscosity/z columns into SI."""
    required = {"pressure", "bo", "rs", "bg", "bw"}
    optional = {"bwinj", "bginj", "oil_viscosity", "gas_viscosity", "z"}
    _check_columns(frame, required, optional)
    rows = []
    for row_number, row in enumerate(frame.to_dict("records"), start=2):
        values = {}
        for name, value in row.items():
            quantity = "viscosity" if name.endswith("viscosity") else ("dimensionless" if name == "z" else name)
            values[name] = to_si(_number(value, name, row_number), quantity, units)
        rows.append(PVTProperties(**values))
    return PVTTable(tuple(rows))
