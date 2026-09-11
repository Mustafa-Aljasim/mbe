"""Piecewise-linear interpolation preserves positive table data without overshoot."""

from dataclasses import dataclass, fields

import numpy as np

from material_balance_studio.domain.models import PVTProperties, PVTTable
from material_balance_studio.domain.validation import PVTOutOfRangeError, finite_value


@dataclass(frozen=True)
class TablePVTModel:
    table: PVTTable

    @property
    def pressure_nodes(self) -> tuple[float, ...]:
        return tuple(row.pressure for row in self.table.rows)

    @property
    def pressure_bounds(self) -> tuple[float, float]:
        return self.pressure_nodes[0], self.pressure_nodes[-1]

    def properties_at_pressure(self, pressure: float) -> PVTProperties:
        """Interpolate all supplied columns; out-of-range requests always raise."""
        finite_value("PVT pressure", pressure, positive=True)
        low, high = self.pressure_bounds
        if not low <= pressure <= high:
            raise PVTOutOfRangeError(
                f"Pressure {pressure:g} Pa is outside PVT range [{low:g}, {high:g}] Pa. "
                "Extrapolation is disabled; supply measured data covering this pressure."
            )
        values: dict[str, float | None] = {"pressure": pressure}
        for column in fields(PVTProperties):
            name = column.name
            if name == "pressure":
                continue
            samples = [getattr(row, name) for row in self.table.rows]
            values[name] = None if samples[0] is None else float(
                np.interp(pressure, self.pressure_nodes, samples)
            )
        return PVTProperties(**values)
