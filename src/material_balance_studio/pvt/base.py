"""Model abstraction; material balance is independent of property generation."""

from typing import Protocol

from material_balance_studio.domain.models import PVTProperties


class PVTModel(Protocol):
    @property
    def pressure_bounds(self) -> tuple[float, float]:
        """Supported pressure range in Pa absolute."""
        ...

    @property
    def pressure_nodes(self) -> tuple[float, ...]:
        """Knots for bracketing; smooth future models may return only their bounds."""
        ...

    def properties_at_pressure(self, pressure: float) -> PVTProperties:
        """Return validated properties at the requested pressure."""
        ...
