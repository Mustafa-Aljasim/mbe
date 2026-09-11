"""Immutable domain records. Pressure: Pa absolute; volumes: m³; time: dates."""

from dataclasses import dataclass, field
from datetime import date
from typing import TYPE_CHECKING

from .validation import EngineeringValidationError, finite_value

if TYPE_CHECKING:
    from material_balance_studio.pvt.base import PVTModel
    from material_balance_studio.aquifer.base import AquiferModel
    from material_balance_studio.aquifer.state import AquiferState, AquiferStepResult

CUMULATIVE_FIELDS = ("np", "gp", "wp", "winj", "ginj")


@dataclass(frozen=True)
class PVTProperties:
    """FVFs in reservoir m³/surface m³; Rs in standard gas m³/stock-tank oil m³.

    Optional injection FVFs default explicitly to resident Bw/Bg in withdrawal.
    Optional viscosities are in Pa.s; z is dimensionless.
    """

    pressure: float
    bo: float
    rs: float
    bg: float
    bw: float
    bwinj: float | None = None
    bginj: float | None = None
    oil_viscosity: float | None = None
    gas_viscosity: float | None = None
    z: float | None = None

    def __post_init__(self) -> None:
        for name in ("pressure", "bo", "bg", "bw"):
            finite_value(name, getattr(self, name), positive=True)
        finite_value("rs", self.rs, nonnegative=True)
        for name in ("bwinj", "bginj", "oil_viscosity", "gas_viscosity", "z"):
            value = getattr(self, name)
            if value is not None:
                finite_value(name, value, positive=True)


@dataclass(frozen=True)
class PVTTable:
    """Strictly increasing measured pressure rows, with consistent optional columns."""

    rows: tuple[PVTProperties, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "rows", tuple(self.rows))
        if len(self.rows) < 2:
            raise EngineeringValidationError("PVT table requires at least two pressure rows.")
        if any(b.pressure <= a.pressure for a, b in zip(self.rows, self.rows[1:])):
            raise EngineeringValidationError("PVT pressure rows must be strictly increasing and unique.")
        for name in ("bwinj", "bginj", "oil_viscosity", "gas_viscosity", "z"):
            present = [getattr(row, name) is not None for row in self.rows]
            if any(present) and not all(present):
                raise EngineeringValidationError(f"PVT column {name} must be complete or omitted.")


@dataclass(frozen=True)
class ReservoirTank:
    """Single equilibrated tank; m = initial gas-cap volume / initial oil volume."""

    initial_date: date
    initial_pressure: float
    oil_in_place: float
    swc: float
    cf: float
    cw: float
    m: float
    pvt_model: "PVTModel"
    name: str = "Single tank"
    aquifer: "AquiferModel | None" = None  # None uses the explicit NoAquifer implementation

    def __post_init__(self) -> None:
        if type(self.initial_date) is not date:
            raise EngineeringValidationError("initial_date must be a calendar date.")
        for name in ("initial_pressure", "oil_in_place"):
            finite_value(name, getattr(self, name), positive=True)
        for name in ("cf", "cw", "m", "swc"):
            finite_value(name, getattr(self, name), nonnegative=True)
        if self.swc >= 1:
            raise EngineeringValidationError("swc must satisfy 0 <= Swc < 1.")
        low, high = self.pvt_model.pressure_bounds
        if not low <= self.initial_pressure <= high:
            raise EngineeringValidationError(
                f"Initial pressure {self.initial_pressure:g} Pa is outside PVT range "
                f"[{low:g}, {high:g}] Pa. Extend measured PVT coverage."
            )


@dataclass(frozen=True)
class CumulativeVolumes:
    """Surface production/injection totals (also used for nonnegative step increments)."""

    np: float = 0.0
    gp: float = 0.0
    wp: float = 0.0
    winj: float = 0.0
    ginj: float = 0.0

    def __post_init__(self) -> None:
        for name in CUMULATIVE_FIELDS:
            finite_value(name, getattr(self, name), nonnegative=True)


@dataclass(frozen=True)
class HistoryRecord:
    """Cumulative volumes since tank initial_date; observed pressure is optional."""

    date: date
    cumulative: CumulativeVolumes = field(default_factory=CumulativeVolumes)
    observed_pressure: float | None = None

    def __post_init__(self) -> None:
        if type(self.date) is not date:
            raise EngineeringValidationError("History date must be a calendar date.")
        if self.observed_pressure is not None:
            finite_value("observed_pressure", self.observed_pressure, positive=True)


@dataclass(frozen=True)
class ExpansionTerms:
    """Reservoir m³ per stock-tank m³ of initial oil; Efw includes (1 + m)."""

    eo: float
    eg: float
    meg: float
    efw: float

    @property
    def et(self) -> float:
        return self.eo + self.meg + self.efw


@dataclass(frozen=True)
class WithdrawalTerms:
    """Cumulative reservoir m³ evaluated at current tank pressure."""

    production: float
    water_injection: float
    gas_injection: float

    @property
    def net(self) -> float:
        return self.production - self.water_injection - self.gas_injection


@dataclass(frozen=True)
class BalanceTerms:
    """All cumulative material-balance terms and signed closure in reservoir m³."""

    withdrawal: WithdrawalTerms
    expansion: ExpansionTerms
    oil_expansion_support: float
    gas_cap_expansion_support: float
    rock_water_expansion_support: float
    total_expansion_support: float
    residual: float
    absolute_residual: float
    relative_residual: float
    aquifer_support: float = 0.0

    @property
    def total_support(self) -> float:
        return self.total_expansion_support + self.aquifer_support


@dataclass(frozen=True)
class SolverSettings:
    """Both absolute and relative closure limits must pass after root finding."""

    pressure_tolerance: float = 1e-5  # Pa
    absolute_residual_tolerance: float = 1e-5  # reservoir m³
    relative_residual_tolerance: float = 1e-8
    normalization_floor: float = 1.0  # reservoir m³; stabilizes zero-net-voidage QC
    max_iterations: int = 100
    minimum_pressure: float = 1.0  # Pa absolute
    maximum_pressure: float | None = None
    bracket_subdivisions: int = 8

    def __post_init__(self) -> None:
        for name in ("pressure_tolerance", "absolute_residual_tolerance",
                     "relative_residual_tolerance", "normalization_floor", "minimum_pressure"):
            finite_value(name, getattr(self, name), positive=True)
        if self.maximum_pressure is not None:
            finite_value("maximum_pressure", self.maximum_pressure, positive=True)
            if self.maximum_pressure <= self.minimum_pressure:
                raise EngineeringValidationError("Maximum pressure must exceed minimum pressure.")
        for name in ("max_iterations", "bracket_subdivisions"):
            if type(getattr(self, name)) is not int or getattr(self, name) < 1:
                raise EngineeringValidationError(f"{name} must be a positive integer.")


@dataclass(frozen=True)
class SolverDiagnostics:
    converged: bool
    iterations: int
    function_calls: int
    initial_guess: float
    pressure_bounds: tuple[float, float]
    bracket: tuple[float, float] | None
    final_residual: float | None
    message: str
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class PressureSolution:
    pressure: float | None
    balance: BalanceTerms | None
    diagnostics: SolverDiagnostics
    aquifer_step: "AquiferStepResult | None" = None


@dataclass(frozen=True)
class SimulationState:
    step_index: int
    date: date
    pressure: float
    previous_pressure: float | None
    elapsed_days: int
    cumulative: CumulativeVolumes
    increments: CumulativeVolumes
    pvt: PVTProperties
    balance: BalanceTerms
    solver: SolverDiagnostics
    observed_pressure: float | None = None
    warnings: tuple[str, ...] = ()
    aquifer_state: "AquiferState | None" = None
    aquifer_step: "AquiferStepResult | None" = None

    @property
    def pressure_error(self) -> float | None:
        return None if self.observed_pressure is None else self.pressure - self.observed_pressure


@dataclass(frozen=True)
class FailedTimestep:
    date: date
    increments: CumulativeVolumes
    solution: PressureSolution


@dataclass(frozen=True)
class SimulationResult:
    initial_state: SimulationState
    states: tuple[SimulationState, ...]
    failed_timestep: FailedTimestep | None = None
    warnings: tuple[str, ...] = ()

    @property
    def converged(self) -> bool:
        return self.failed_timestep is None
