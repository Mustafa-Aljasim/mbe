"""Sparse laboratory data boundary; independent of the strict legacy PVT table."""
from dataclasses import dataclass
from math import isfinite
import pandas as pd
from material_balance_studio.domain.validation import EngineeringValidationError
from material_balance_studio.units.conversions import to_si, from_si, UnitSystem

PROPERTIES = ("bo", "rs", "oil_viscosity", "z", "bg", "bw")
QUANTITIES = {**{k: k for k in ("pressure", "bo", "rs", "bg", "bw")},
              "oil_viscosity": "viscosity", "z": "dimensionless", "pb": "pressure", "co": "compressibility"}


@dataclass(frozen=True)
class LabPoint:
    pressure: float
    bo: float | None = None
    rs: float | None = None
    oil_viscosity: float | None = None
    z: float | None = None
    bg: float | None = None
    bw: float | None = None

    def __post_init__(self):
        if not isfinite(self.pressure) or self.pressure <= 0:
            raise EngineeringValidationError("Laboratory pressure must be finite and positive.")
        for name in PROPERTIES:
            value = getattr(self, name)
            if value is not None and (not isfinite(value) or value < 0 or (name != "rs" and value == 0)):
                raise EngineeringValidationError(f"Laboratory {name} must be finite and {'nonnegative' if name == 'rs' else 'positive'}.")


@dataclass(frozen=True)
class LaboratoryData:
    points: tuple[LabPoint, ...] = ()

    def __post_init__(self):
        pressures = [p.pressure for p in self.points]
        if len(set(pressures)) != len(pressures):
            raise EngineeringValidationError("Duplicate laboratory pressures: combine measurements into one row.")

    def measurements(self, property: str) -> tuple[tuple[float, float], ...]:
        return tuple((p.pressure, getattr(p, property)) for p in self.points if getattr(p, property) is not None)

    @property
    def warnings(self) -> tuple[str, ...]:
        result = []
        for name in PROPERTIES:
            data = self.measurements(name)
            if len(data) < 3:
                result.append(f"{name}: {len(data)} laboratory points; sparse or absent support.")
        for p in self.points:
            for name, low, high in (("bo", .7, 4), ("bw", .8, 1.5), ("z", .1, 2), ("bg", 1e-5, 2), ("rs", 0, 1000)):
                v = getattr(p, name)
                if v is not None and not low <= v <= high:
                    result.append(f"Suspicious {name}={v:g} (canonical SI) at {p.pressure:g} Pa; verify units/fluid.")
        return tuple(result)


def lab_from_frame(frame: pd.DataFrame, units: UnitSystem | str = "SI") -> LaboratoryData:
    UnitSystem(units)
    if frame.columns.duplicated().any():
        raise EngineeringValidationError("Duplicate laboratory column headers.")
    if "pressure" not in frame or set(frame.columns) - {"pressure", *PROPERTIES}:
        raise EngineeringValidationError("Use pressure and optional bo, rs, oil_viscosity, z, bg, bw headers; set file units explicitly.")
    points = []
    for row_number, (_, row) in enumerate(frame.iterrows(), 1):
        if row.isna().all():
            continue  # the editor's unused blank row is not a measurement
        if pd.isna(row["pressure"]):
            raise EngineeringValidationError(f"Laboratory row {row_number}: pressure is missing.")
        values = {}
        for name in frame.columns:
            if pd.isna(row[name]):
                continue
            try:
                values[name] = to_si(float(row[name]), QUANTITIES[name], units)
            except (ValueError, TypeError) as exc:
                raise EngineeringValidationError(f"Laboratory row {row_number}: invalid {name}.") from exc
        points.append(LabPoint(**values))
    return LaboratoryData(tuple(sorted(points, key=lambda p: p.pressure)))


def lab_to_frame(lab: LaboratoryData, units: UnitSystem | str = "SI") -> pd.DataFrame:
    return pd.DataFrame([{k: None if getattr(p, k) is None else from_si(getattr(p, k), QUANTITIES[k], units)
                          for k in ("pressure", *PROPERTIES)} for p in lab.points], columns=("pressure", *PROPERTIES))
