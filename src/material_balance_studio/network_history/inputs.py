"""Canonical immutable input records; schedules are engineer supplied, never fitted."""
from dataclasses import dataclass
from datetime import date
from math import fsum
from material_balance_studio.domain.models import CUMULATIVE_FIELDS
from material_balance_studio.domain.validation import finite_value


@dataclass(frozen=True)
class NetworkObservation:
    tank: str
    date: date
    pressure: float
    source: str = ""
    sigma: float | None = None
    qc_flag: str = ""
    include_in_match: bool = True
    note: str = ""

    def __post_init__(self):
        if type(self.date) is not date or type(self.include_in_match) is not bool:
            raise ValueError("Observation requires a calendar date and boolean include flag.")
        finite_value("Observed pressure",self.pressure,positive=True)
        if self.sigma is not None:
            finite_value("Pressure sigma",self.sigma,positive=True)


@dataclass(frozen=True)
class AllocationRow:
    date: date
    stream: str
    fractions: tuple

    def __post_init__(self):
        object.__setattr__(self,"fractions",tuple(tuple(pair) for pair in self.fractions))
        if type(self.date) is not date or self.stream not in CUMULATIVE_FIELDS:
            raise ValueError("Allocation requires a calendar effective date and valid stream.")
        if not self.fractions or len({n for n,_ in self.fractions})!=len(self.fractions):
            raise ValueError("Allocation requires unique tank fractions at each stream/date.")
        for _,v in self.fractions:
            finite_value("Allocation fraction",v,nonnegative=True)
            if v>1:
                raise ValueError("Allocation fraction cannot exceed one.")
        if abs(fsum(v for _,v in self.fractions)-1.)>1e-10:
            raise ValueError("Allocation fractions must sum to one within 1e-10; no renormalization is applied.")


@dataclass(frozen=True)
class HistoryPlan:
    allocated_streams: tuple = ()
    field_history: tuple = ()
    schedules: tuple[AllocationRow,...] = ()
    replace_direct: bool = False

    def __post_init__(self):
        for key in ("allocated_streams","field_history","schedules"):
            object.__setattr__(self,key,tuple(getattr(self,key)))
        if len(set(self.allocated_streams))!=len(self.allocated_streams) or not set(self.allocated_streams)<=set(CUMULATIVE_FIELDS):
            raise ValueError("Select unique valid allocated streams.")
        if type(self.replace_direct) is not bool:
            raise ValueError("Direct-stream replacement requires an explicit boolean choice.")


def observations_from_network(network):
    return tuple(NetworkObservation(t.name,h.date,h.observed_pressure) for t in network.tanks for h in t.history if h.observed_pressure is not None)


def validate_observations(network,observations):
    names={t.name for t in network.tanks}
    if len({(o.tank,o.date) for o in observations})!=len(observations):
        raise ValueError("Duplicate tank/date pressure observation.")
    for o in observations:
        if o.tank not in names or o.date<network.initial_date:
            raise ValueError("Observation references an unknown tank or precedes network initialization.")
