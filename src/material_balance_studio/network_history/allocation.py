"""Allocate interval increments, splitting at piecewise-constant effective dates."""
from dataclasses import dataclass,replace,asdict
from math import fsum,isclose
import hashlib,json
from material_balance_studio.domain.models import HistoryRecord,CumulativeVolumes,CUMULATIVE_FIELDS
from material_balance_studio.tank_network.history import volumes_at
from .inputs import HistoryPlan,validate_observations,observations_from_network


@dataclass(frozen=True)
class AllocationAudit:
    start: object
    end: object
    tank: str
    stream: str
    field_increment: float
    fraction: float
    allocated_increment: float
    tank_cumulative: float
    interval_conservation_error: float
    cumulative_conservation_error: float


@dataclass(frozen=True)
class PreparedHistory:
    network: object
    plan: HistoryPlan
    observations: tuple
    audit: tuple
    fingerprint: str
    warnings: tuple


def prepare_history(network,plan=None,observations=None):
    plan=plan or HistoryPlan()
    observations=observations_from_network(network) if observations is None else tuple(observations)
    validate_observations(network,observations)
    names={t.name for t in network.tanks}
    streams=set(plan.allocated_streams)
    if len({(s.stream,s.date) for s in plan.schedules})!=len(plan.schedules):
        raise ValueError("Duplicate allocation stream/effective date.")
    if {s.stream for s in plan.schedules}!=streams:
        raise ValueError("Schedules must cover exactly the selected field-allocated streams.")
    for s in plan.schedules:
        if {n for n,_ in s.fractions}!=names:
            raise ValueError("Each allocation row must name every network tank exactly once, including zero shares.")
        if s.date<network.initial_date:
            raise ValueError("Allocation effective dates cannot precede network initialization.")
    schedules={s:sorted([r for r in plan.schedules if r.stream==s],key=lambda r:r.date) for s in streams}
    if any(rows[0].date!=network.initial_date for rows in schedules.values()):
        raise ValueError("Every allocated stream requires a schedule effective on the initial date.")
    if streams and not plan.field_history:
        raise ValueError("Allocated streams require field-total history.")
    if plan.field_history and not streams:
        raise ValueError("Field history requires explicitly selected allocated streams; direct mode uses tank histories.")
    if not plan.replace_direct and any(getattr(h.cumulative,s)>0 for n in network.tanks for h in n.history for s in streams):
        raise ValueError("Direct and allocated data conflict for a selected stream. Explicitly enable replacement; streams are never added together.")
    field=replace(network.tanks[0],history=plan.field_history) if streams else None
    dates=sorted({network.initial_date,*[h.date for n in network.tanks for h in n.history],*[h.date for h in plan.field_history],*[s.date for s in plan.schedules],*[o.date for o in observations]})
    obs={(o.tank,o.date):o.pressure for o in observations}
    histories={n.name:[] for n in network.tanks}
    cumulative={(n,s):0. for n in names for s in streams}
    audit=[]
    for i,day in enumerate(dates):
        elapsed=(day-network.initial_date).days*86400.
        if i and streams:
            prior=dates[i-1]
            before=volumes_at(field,(prior-network.initial_date).days*86400.)[0]
            after=volumes_at(field,elapsed)[0]
            for stream in sorted(streams):
                active=next(s for s in reversed(schedules[stream]) if s.date<=prior)
                delta=getattr(after,stream)-getattr(before,stream)
                allocated={n:delta*f for n,f in active.fractions}
                for n,v in allocated.items():
                    cumulative[n,stream]+=v
                error=fsum(allocated.values())-delta
                cum_error=fsum(cumulative[n,stream] for n in names)-getattr(after,stream)
                if not isclose(error,0.,abs_tol=1e-8+1e-10*abs(delta)) or not isclose(cum_error,0.,abs_tol=1e-8+1e-10*abs(getattr(after,stream))):
                    raise ValueError("Allocation conservation failed; no history was committed.")
                audit.extend(AllocationAudit(prior,day,n,stream,delta,f,allocated[n],cumulative[n,stream],error,cum_error) for n,f in active.fractions)
        for node in network.tanks:
            direct,_=volumes_at(node,elapsed)
            totals=CumulativeVolumes(**{s:cumulative[node.name,s] if s in streams else getattr(direct,s) for s in CUMULATIVE_FIELDS})
            histories[node.name].append(HistoryRecord(day,totals,obs.get((node.name,day))))
    prepared=replace(network,tanks=tuple(replace(n,history=tuple(histories[n.name])) for n in network.tanks))
    payload=dict(network=asdict(network),plan=asdict(plan),observations=[asdict(o) for o in observations])
    fingerprint=hashlib.sha256(json.dumps(payload,sort_keys=True,default=str).encode()).hexdigest()
    warnings=[]
    if any(h.observed_pressure is not None for h in plan.field_history):
        warnings.append("Field pressure observations are not allocated; supply tank-specific pressure metadata.")
    if streams:
        warnings.append("Field increments use constant rates between cumulative records. Fractions are piecewise constant from their effective dates; no interpolation or normalization of fractions.")
    return PreparedHistory(prepared,plan,observations,tuple(audit),fingerprint,tuple(warnings))
