"""Pressure metrics, physical acceptance margins and signed network drive terms."""
from dataclasses import dataclass
from math import fsum
from types import SimpleNamespace
from material_balance_studio.matching.optimizer import metrics
from material_balance_studio.matching.weighting import residual_scales
from material_balance_studio.diagnostics.drive_indices import drive_indices
from material_balance_studio.tank_network.residuals import accepted
from material_balance_studio.tank_network import NetworkSettings
from .parameters import registry


@dataclass(frozen=True)
class PressureRow:
    tank: str
    date: object
    observed: float
    calculated: float | None
    residual: float | None
    sigma: float | None
    objective_residual: float | None
    included: bool
    source: str
    qc_flag: str
    note: str
    effective_sigma: float | None = None


def pressure_diagnostics(result,observations,mode="Unweighted",default_sigma=None):
    states={} if result is None else {(s.time.date(),t.name):t.pressure for s in (result.initial_state,*result.states) for t in s.tanks}
    rows=[]
    for o in sorted(observations,key=lambda o:(o.date,o.tank)):
        calculated=states.get((o.date,o.tank))
        residual=None if calculated is None else calculated-o.pressure
        scale=residual_scales([o],mode,default_sigma)[0] if o.include_in_match else None
        rows.append(PressureRow(o.tank,o.date,o.pressure,calculated,residual,o.sigma,residual/scale if residual is not None and scale is not None else None,o.include_in_match,o.source,o.qc_flag,o.note,
            float(scale) if scale is not None and mode=="Pressure uncertainty weighted" else None))
    grouped=tuple((n,tuple(metrics([r.residual for r in rows if r.tank==n and r.included and r.residual is not None]).items())) for n in sorted({o.tank for o in observations}|({t.name for t in result.network.tanks} if result else set())))
    aggregate=tuple(metrics([r.residual for r in rows if r.included and r.residual is not None]).items())
    return tuple(rows),grouped,aggregate


def closure_summary(result):
    settings=NetworkSettings()
    states=result.internal_states or (result.initial_state,)
    maximum=max(t.relative_residual for s in states for t in s.tanks)
    absolute=max(abs(t.residual) for s in states for t in s.tanks)
    conservation=max(max(abs(s.transfer_error),abs(s.incremental_transfer_error)) for s in states)
    valid=result.converged and all(accepted(s,settings) for s in states)
    margin=settings.relative_tolerance-maximum
    return dict(valid=valid,maximum_relative=maximum,maximum_absolute=absolute,transfer_error=conservation,closure_margin=margin,
        status="FAIL" if not valid else "NEAR TOLERANCE" if maximum>=.8*settings.relative_tolerance else "PASS")


def data_qc(network,specs,observations):
    active=[s for s in specs if s.active]
    included=[o for o in observations if o.include_in_match]
    if len(included)<=len(active):
        raise ValueError("FAIL: included network pressure observations must outnumber active parameters.")
    warnings=[]
    if len(included)<max(5,2*len(active)+1):
        warnings.append("Weak network overdetermination: few observations relative to active parameters.")
    targets=registry(network)
    for node in network.tanks:
        count=sum(o.tank==node.name for o in included)
        nparams=sum(targets[s.name].edge is None and targets[s.name].owner==node.name for s in active)
        if count==0 and nparams:
            warnings.append(f"STRONG CAUTION: {node.name} has {nparams} active tank parameters but no included pressure observations; communication does not establish direct constraint.")
        elif nparams and count<=nparams:
            warnings.append(f"{node.name}: limited direct observations relative to its active parameters.")
    for s in active:
        t=targets[s.name]
        if t.edge and not any(o.tank in t.edge for o in included):
            warnings.append(f"STRONG CAUTION: {s.name} has no observations at either endpoint.")
    if any(o.qc_flag or (o.source and o.source!="Average reservoir pressure") for o in included):
        warnings.append("Included pressures carry source/quality metadata; check reservoir representativeness.")
    return tuple(warnings)


def drive_terms(state):
    base=drive_indices(SimpleNamespace(balance=state.components))
    values=dict(base["values"])
    denominator=state.components.withdrawal.production
    values["TDI"]=state.intertank_support/denominator if denominator>1e-12 else None
    total=fsum(values.values()) if all(v is not None for v in values.values()) else None
    return dict(values=values,total=total,closure_error=None if total is None else total-1,
                supports={**base["supports"],"TDI":state.intertank_support},denominator=denominator)


def field_drive(state):
    supports={k:fsum(drive_terms(t)["supports"][k] for t in state.tanks) for k in ("DDI","SDI","CDI","WDI","IDI","TDI")}
    denominator=fsum(t.components.withdrawal.production for t in state.tanks)
    values={k:v/denominator if denominator>1e-12 else None for k,v in supports.items()}
    return dict(supports=supports,denominator=denominator,values=values)


def allocation_crosscheck(plan,rows,network):
    warnings=[]
    for stream in plan.allocated_streams:
        schedules=sorted([s for s in plan.schedules if s.stream==stream],key=lambda s:s.date)
        for old,new in zip(schedules,schedules[1:]):
            for name,fraction in new.fractions:
                if abs(fraction-dict(old.fractions)[name])<.2:
                    continue
                series=sorted([r for r in rows if r.tank==name],key=lambda r:r.date)
                pi=next(t.reservoir.initial_pressure for t in network.tanks if t.name==name)
                for a,b in zip(series,series[1:]):
                    if abs((b.date-new.date).days)<=30 and (abs(b.observed-a.observed)>=.05*pi or (a.residual is not None and b.residual is not None and abs(b.residual-a.residual)>=1e6)):
                        warnings.append(f"{name}: major {stream} allocation change on {new.date} near pressure/residual change at {b.date}; diagnostic only, allocation unchanged.")
    return tuple(dict.fromkeys(warnings))
