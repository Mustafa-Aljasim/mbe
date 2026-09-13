"""Deterministic engineer-selected allocation cases; fractions are never fitted."""
from dataclasses import dataclass, replace
from math import isfinite
from material_balance_studio.network_history import prepare_history
from material_balance_studio.network_matching import history_match_network
from material_balance_studio.network_matching.diagnostics import drive_terms
from .network_context import context


def shift_allocation(plan, stream, donor, receiver, delta, effective_dates=None):
    """Transfer a fraction from donor to receiver, at selected existing effective rows."""
    if stream not in plan.allocated_streams or donor==receiver or not isfinite(delta):
        raise ValueError('Select an allocated stream, distinct tanks and a finite fraction shift.')
    available={r.date for r in plan.schedules if r.stream==stream}
    selected=available if effective_dates is None else set(effective_dates)
    if not selected or not selected<=available:
        raise ValueError('Select existing effective dates for this stream.')
    rows=[]
    for row in plan.schedules:
        if row.stream==stream and row.date in selected:
            fractions=dict(row.fractions)
            if donor not in fractions or receiver not in fractions:
                raise ValueError('Unknown allocation donor/receiver.')
            fractions[donor]-=delta
            fractions[receiver]+=delta
            # AllocationRow rejects impossible shifts; no clipping or renormalization.
            row=replace(row,fractions=tuple(fractions.items()))
        rows.append(row)
    return replace(plan,schedules=tuple(rows))


@dataclass(frozen=True)
class AllocationCase:
    name: str
    plan: object
    match: object | None
    drift: tuple
    drives: tuple
    transfers: tuple
    failure: str | None
    robustness: str


@dataclass(frozen=True)
class AllocationStudy:
    scenario_id: str
    cases: tuple
    robustness: str
    warnings: tuple


def allocation_study(scenario, cases, max_nfev=40):
    specs,_=context(scenario)
    cases=tuple(cases)
    if not cases or len(cases)>21 or len({name for name,_ in cases})!=len(cases):
        raise ValueError('Supply 1–21 uniquely named allocation cases.')
    if type(max_nfev) is not int or not 1<=max_nfev<=1000:
        raise ValueError('max_nfev must be an integer from 1 to 1000.')
    base=dict(scenario.fitted_parameters)
    results=[]
    for name,plan in cases:
        try:
            # A study varies allocation only: field totals and selected stream mode stay fixed.
            if plan.field_history!=scenario.history_plan.field_history or plan.allocated_streams!=scenario.history_plan.allocated_streams or plan.replace_direct!=scenario.history_plan.replace_direct:
                raise ValueError('Allocation studies preserve field history and direct/allocated stream mode.')
            prepare_history(scenario.base_network,plan,scenario.observations)
            fit=history_match_network(scenario.base_network,specs,scenario.observations,plan,
                scenario.mode,scenario.default_sigma,scenario.settings,max_nfev=max_nfev,scenario_name=name)
            if not fit.converged:
                results.append(AllocationCase(name,plan,fit,(),(),(),fit.termination_message,'UNRESOLVED'))
                continue
            fitted=dict(fit.fitted_parameters)
            drift=[]
            magnitudes=[]
            for s in specs:
                change=fitted[s.name]-base[s.name]
                meaningful=abs(base[s.name])>1e-8*s.scale
                percent=100*change/abs(base[s.name]) if meaningful else None
                fraction=abs(change)/abs(base[s.name]) if meaningful else abs(change)/(s.upper-s.lower)
                magnitudes.append(fraction)
                drift.append((s.name,base[s.name],fitted[s.name],change,percent,fraction,
                              'relative to base' if meaningful else 'fraction of bound span; base near zero'))
            maximum=max(magnitudes)
            status='ROBUST' if maximum<=.10 else 'MODERATELY SENSITIVE' if maximum<=.25 else 'HIGHLY SENSITIVE'
            final=fit.final_simulation.states[-1] if fit.final_simulation.states else fit.final_simulation.initial_state
            drives=tuple((t.name,tuple(drive_terms(t)['values'].items())) for t in final.tanks)
            transfers=tuple((e.from_tank,e.to_tank,e.cumulative_transfer) for e in final.connections)
            results.append(AllocationCase(name,plan,fit,tuple(drift),drives,transfers,None,status))
        except (ValueError,ArithmeticError) as exc:
            results.append(AllocationCase(name,plan,None,(),(),(),str(exc),'UNRESOLVED'))
    order={'ROBUST':0,'MODERATELY SENSITIVE':1,'HIGHLY SENSITIVE':2,'UNRESOLVED':3}
    status=max((r.robustness for r in results),key=order.get)
    return AllocationStudy(scenario.scenario_id,tuple(results),status,
        ('Deterministic ranges over selected allocation assumptions, not statistical confidence intervals.',
         'ROBUST ≤10%, MODERATELY SENSITIVE ≤25%, HIGHLY SENSITIVE >25% maximum absolute parameter drift. Near-zero bases use bound-span fractions.',
         'Failed cases make the study UNRESOLVED. Fractions are supplied inputs; no automatic allocation selection or fitting.'))
