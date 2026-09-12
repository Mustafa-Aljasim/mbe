"""Immutable tank nodes and oriented physical links, canonical Pa, m³ and seconds."""
from dataclasses import dataclass
import re
from material_balance_studio.domain.models import ReservoirTank,HistoryRecord
from material_balance_studio.domain.validation import validate_history,finite_value


@dataclass(frozen=True)
class NetworkTank:
    name: str
    reservoir: ReservoirTank
    history: tuple[HistoryRecord,...] = ()
    minimum_pressure: float | None = None
    maximum_pressure: float | None = None

    def __post_init__(self):
        if not isinstance(self.name,str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9 _-]{0,39}",self.name) or self.name!=self.name.strip():
            raise ValueError("Tank name must start with a letter and contain at most 40 letters, digits, spaces, underscores or hyphens.")
        object.__setattr__(self,"history",tuple(self.history))
        if self.history:
            validate_history(self.history,self.reservoir.initial_date)
        for value in (self.minimum_pressure,self.maximum_pressure):
            if value is not None:
                finite_value("Tank pressure bound",value,positive=True)
        lo,hi=self.pressure_bounds
        if not lo<hi or not lo<=self.reservoir.initial_pressure<=hi:
            raise ValueError(f"{self.name}: pressure bounds must overlap PVT coverage and contain initial pressure.")

    @property
    def pressure_bounds(self):
        lo,hi=self.reservoir.pvt_model.pressure_bounds
        return max(1.,lo,self.minimum_pressure or 1.),min(hi,self.maximum_pressure or hi)


@dataclass(frozen=True)
class Connection:
    from_tank: str
    to_tank: str
    transmissibility: float
    enabled: bool = True

    def __post_init__(self):
        finite_value("Transmissibility (reservoir m³/(Pa s))",self.transmissibility,nonnegative=True)
        if self.from_tank==self.to_tank:
            raise ValueError("A connection cannot connect a tank to itself.")
        if type(self.enabled) is not bool:
            raise ValueError("Connection enabled flag must be boolean.")

    @property
    def key(self):
        return (self.from_tank,self.to_tank)

    @property
    def effective_transmissibility(self):
        return self.transmissibility if self.enabled else 0.


@dataclass(frozen=True)
class TankNetwork:
    tanks: tuple[NetworkTank,...]
    connections: tuple[Connection,...] = ()

    def __post_init__(self):
        object.__setattr__(self,"tanks",tuple(self.tanks))
        object.__setattr__(self,"connections",tuple(self.connections))
        if not self.tanks:
            raise ValueError("Network requires at least one tank.")
        names=[t.name for t in self.tanks]
        if len(set(names))!=len(names):
            raise ValueError("Tank names must be unique.")
        if len({t.reservoir.initial_date for t in self.tanks})!=1:
            raise ValueError("Phase 5A requires a common initial reference date for every tank.")
        seen=set()
        for edge in self.connections:
            if edge.from_tank not in names or edge.to_tank not in names:
                raise ValueError("Connection references an unknown tank.")
            pair=frozenset(edge.key)
            if pair in seen:
                raise ValueError("Duplicate physical connection, including reversed orientation.")
            seen.add(pair)

    @property
    def initial_date(self):
        return self.tanks[0].reservoir.initial_date


@dataclass(frozen=True)
class NetworkSettings:
    absolute_tolerance: float = 1e-5
    relative_tolerance: float = 1e-8
    normalization_floor: float = 1.
    transfer_absolute_tolerance: float = 1e-8
    transfer_relative_tolerance: float = 1e-12
    max_nfev: int = 150
    max_step_days: float | None = None
    communication_step_limit: float = 0.5
    max_substeps_per_event: int = 4096

    def __post_init__(self):
        for name in ("absolute_tolerance","relative_tolerance","normalization_floor","transfer_absolute_tolerance","transfer_relative_tolerance","communication_step_limit"):
            finite_value(name,getattr(self,name),positive=True)
        if self.communication_step_limit>1:
            raise ValueError("Communication step limit must be <=1 for monotone linear-storage equilibration.")
        if self.max_step_days is not None:
            finite_value("Maximum step days",self.max_step_days,positive=True)
        for name in ("max_nfev","max_substeps_per_event"):
            if type(getattr(self,name)) is not int or getattr(self,name)<1:
                raise ValueError(f"{name} must be a positive integer.")
