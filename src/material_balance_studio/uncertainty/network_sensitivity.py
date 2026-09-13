"""Bound-aware network derivatives, using Phase 4C SVD/confidence/classification."""
from dataclasses import dataclass
import numpy as np
from .network_context import context
from .context import evaluate, objective_label
from .result import IdentifiabilityResult, matrix
from .jacobian import matrix_diagnostics
from .confidence import local_confidence
from .identifiability import classify
from material_balance_studio.network_matching.parameters import registry


@dataclass(frozen=True)
class NetworkInformation:
    local: IdentifiabilityResult
    observation_tanks: tuple
    by_tank: tuple
    coverage: tuple
    connections: tuple
    evaluations: tuple


def analyze_network(scenario, statistical_assumptions=False, relative_step=1e-4):
    if not 1e-6 <= relative_step <= 1e-2:
        raise ValueError('Relative step must lie between 1e-6 and 1e-2.')
    specs,obj = context(scenario)
    values = np.array([s.initial for s in specs])
    base = evaluate(obj,values)
    if base is None:
        raise ValueError('The matched network cannot be evaluated.')
    columns,schemes,errors,warnings = [],[],[],[]
    def derivative(j,h):
        s,x = specs[j],values[j]
        if x-h>=s.lower and x+h<=s.upper:
            offsets,coefficients,label = (-h,h),(-1,1),'central'
        elif x+2*h<=s.upper:
            offsets,coefficients,label = (h,2*h),(4,-1),'forward second order'
        elif x-2*h>=s.lower:
            offsets,coefficients,label = (-h,-2*h),(-4,1),'backward second order'
        else:
            return None,'unavailable'
        rows=[]
        for offset in offsets:
            trial=values.copy()
            trial[j]+=offset
            rows.append(evaluate(obj,trial))
        if any(r is None for r in rows):
            return None,'failed stencil'
        numerator=sum(c*r for c,r in zip(coefficients,rows))
        if label.startswith('forward'):
            numerator-=3*base
        elif label.startswith('backward'):
            numerator+=3*base
        return numerator/(2*h),label
    for j,s in enumerate(specs):
        h=min(relative_step*max(abs(values[j]),s.scale),.1*(s.upper-s.lower))
        for attempt in range(4):
            coarse,scheme=derivative(j,h)
            fine,_=derivative(j,h/2)
            if coarse is not None and fine is not None:
                break
            h/=2
        if coarse is None or fine is None:
            raise ValueError(f'Invalid sensitivity stencil for {s.name}: '+str(obj.records[-1].failure))
        errors.append(float(np.linalg.norm(fine-coarse)/max(np.linalg.norm(fine),1e-20)))
        columns.append(fine)
        schemes.append((s.name,scheme,float(h/2)))
    jac=np.column_stack(columns)
    raw=jac*obj.scales[:,None]
    scales=np.array([s.scale for s in specs])
    pi={t.name:t.reservoir.initial_pressure for t in scenario.base_network.tanks}
    dimensionless=raw*scales/np.array([pi[o.tank] for o in obj.observations])[:,None]
    scaled=jac*scales
    diag=matrix_diagnostics(scaled)
    bounds=[min(v-s.lower,s.upper-v)/(s.upper-s.lower)<1e-4 for s,v in zip(specs,values)]
    stable=max(errors)<=.05 and (diag['condition'] is None or diag['condition']<=1e4)
    confidence=local_confidence(diag['geometry'],specs,values,float(base@base),len(base),
        scenario.mode!='Unweighted',statistical_assumptions,bounds,stable=stable)
    warnings.extend(confidence['warnings'])
    if diag['rank']<len(specs):
        warnings.append('NON-IDENTIFIABLE: rank-deficient scaled network Jacobian; intervals withheld.')
    if not stable:
        warnings.append('Unstable derivatives (>5%) or condition >1e4: local confidence withheld.')
    # Use within-tank pressure span, avoiding artificial information from different initial pressures.
    spread=max((np.ptp([o.pressure for o in obj.observations if o.tank==n])/p
                for n,p in pi.items() if any(o.tank==n for o in obj.observations)),default=0.)
    statuses=classify(specs,values,dimensionless,diag,len(base),spread)
    targets=registry(scenario.base_network)
    coverage=[]
    norms=np.sqrt(np.mean(dimensionless**2,axis=0))
    observed={o.tank for o in obj.observations}
    for j,s in enumerate(specs):
        target=targets[s.name]
        direct=target.owner in observed if target.edge is None else all(n in observed for n in target.edge)
        label='UNINFORMED' if norms[j]<1e-8 else 'WEAKLY INFORMED' if norms[j]<1e-5 else 'DIRECTLY INFORMED' if direct else 'INDIRECTLY INFORMED'
        coverage.append((s.name,label,float(norms[j])))
    for name in pi:
        if name not in observed:
            warnings.append(f'No direct pressure observations exist for {name}. Its parameters are constrained only through network response.')
    # An indirect-only local interval can be misleading even if numerically full rank.
    intervals=list(confidence['intervals'])
    updated=[]
    for j,(name,status,note) in enumerate(statuses):
        if coverage[j][1]!='DIRECTLY INFORMED':
            intervals[j]=None
            note+=' Direct coverage is missing or sensitivity is weak; marginal interval withheld.'
            if status in ('WELL CONSTRAINED','MODERATELY CONSTRAINED'):
                status='POORLY CONSTRAINED'
        if not stable and status in ('WELL CONSTRAINED','MODERATELY CONSTRAINED'):
            status='POORLY CONSTRAINED'
        updated.append((name,status,note))
    tank_rows=tuple(o.tank for o in obj.observations)
    by_tank=tuple((name,tuple(float(v) for v in np.sqrt(np.mean(dimensionless[np.array(tank_rows)==name]**2,axis=0)))) for name in sorted(observed))
    connections=connection_observability(scenario)
    warnings.extend(w for c in connections for w in dict(c)['warnings'])
    local=IdentifiabilityResult(scenario.scenario_id,tuple(s.name for s in specs),tuple(str(o.date) for o in obj.observations),
        scenario.mode,objective_label(scenario.mode),matrix(raw),matrix(dimensionless),matrix(jac),matrix(scaled),
        tuple(schemes),tuple(errors),tuple(float(x) for x in diag['singular']),diag['rank'],diag['condition'],
        tuple(float(x) for x in diag['weak']),matrix(confidence['covariance']) if confidence['covariance'] is not None else None,
        confidence['se'],tuple(intervals),matrix(diag['correlation']) if diag['correlation'] is not None else None,
        matrix(diag['coupling']),tuple(updated),tuple(warnings),sum(not e.valid for e in obj.records),len(obj.records),statistical_assumptions)
    return NetworkInformation(local,tank_rows,by_tank,tuple(coverage),connections,tuple(obj.records))


def connection_observability(scenario):
    observed={o.tank for o in scenario.observations if o.include_in_match}
    last=max(o.date for o in scenario.observations if o.include_in_match)
    states=[s for s in (scenario.final_simulation.initial_state,*scenario.final_simulation.states) if s.time.date()<=last]
    rows=[]
    for edge in scenario.base_network.connections:
        series=[]
        for state in states:
            pressures={t.name:t.pressure for t in state.tanks}
            transfer=next(e for e in state.connections if (e.from_tank,e.to_tank)==edge.key)
            series.append((str(state.time.date()),pressures[edge.from_tank]-pressures[edge.to_tank],transfer.rate))
        pi=max(t.reservoir.initial_pressure for t in scenario.base_network.tanks if t.name in edge.key)
        ratio=max(abs(x[1]) for x in series)/pi
        warnings=[]
        if ratio<1e-3:
            warnings.append(f'{edge.key}: poorly observable connection; maximum |ΔP|/Pi <0.001.')
        if not all(n in observed for n in edge.key):
            warnings.append(f'{edge.key}: pressure coverage missing at an endpoint; T information is indirect.')
        if (last-scenario.base_network.initial_date).days<30:
            warnings.append(f'{edge.key}: observed history shorter than 30 days; communication information may be limited.')
        rows.append(tuple(dict(connection=edge.key,maximum_relative_differential=ratio,series=tuple(series),warnings=tuple(warnings)).items()))
    return tuple(rows)


def coupling_pairs(scenario):
    specs,_=context(scenario)
    targets=registry(scenario.base_network)
    names=[s.name for s in specs]
    pairs=[]
    for i,a in enumerate(names):
        for b in names[i+1:]:
            ta,tb=targets[a],targets[b]
            if ta.edge or tb.edge:
                edge,other=(ta,tb) if ta.edge else (tb,ta)
                if other.owner in edge.edge and (other.field=='oil_in_place' or 'aquifer' in a+b):
                    pairs.append((a,b))
            elif (ta.field==tb.field=='oil_in_place') or (ta.owner==tb.owner and {ta.field,tb.field}=={'oil_in_place','m'}):
                pairs.append((a,b))
    return tuple(pairs)
