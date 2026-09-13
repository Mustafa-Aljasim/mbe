"""Network-qualified identities reuse Phase 4B physical parameter rules."""
from dataclasses import dataclass,replace
from material_balance_studio.matching.parameters import MatchParameter,parameter_registry,parameter_value,validate_parameters,candidate_tank


@dataclass(frozen=True)
class Target:
    owner: str
    field: str
    description: str
    unit: str
    value: float
    edge: tuple | None = None
    enabled: bool = True


def registry(network):
    result={}
    for node in network.tanks:
        for field,(description,unit) in parameter_registry(node.reservoir).items():
            suffix="N" if field=="oil_in_place" else field.replace(".",":")
            result[f"tank:{node.name}:{suffix}"]=Target(node.name,field,description,unit,parameter_value(node.reservoir,field))
    for edge in network.connections:
        result[f"connection:{edge.from_tank}:{edge.to_tank}:T"]=Target(edge.from_tank+" → "+edge.to_tank,"transmissibility","Connection transmissibility","aquifer_productivity",edge.transmissibility,edge.key,edge.enabled)
    return result


def validate_specs(network,specs):
    targets=registry(network)
    if not any(s.active for s in specs):
        raise ValueError("Select at least one active network parameter.")
    if len({s.name for s in specs})!=len(specs):
        raise ValueError("Duplicate network parameter identity.")
    for spec in specs:
        if spec.name not in targets or spec.unit!=targets[spec.name].unit:
            raise ValueError(f"Unapproved network parameter/unit: {spec.name}")
        if not spec.active:
            continue
        t=targets[spec.name]
        if t.edge:
            if spec.lower<0 or not t.enabled:
                raise ValueError("Active transmissibility must be nonnegative and belong to an enabled connection.")
        else:
            tank=next(n.reservoir for n in network.tanks if n.name==t.owner)
            validate_parameters(tank,[replace(spec,name=t.field)])


def candidate_network(network,specs,physical):
    from math import isfinite
    active=[s for s in specs if s.active]
    if len(active)!=len(physical):
        raise ValueError("Network parameter vector length mismatch.")
    targets=registry(network)
    tank_specs={n.name:[] for n in network.tanks}
    tank_values={n.name:[] for n in network.tanks}
    edges={}
    for s,v in zip(active,physical):
        if not isfinite(v) or not s.lower<=v<=s.upper:
            raise ValueError(f"{s.name}: candidate outside explicit bounds.")
        t=targets[s.name]
        if t.edge:
            edges[t.edge]=float(v)
        else:
            tank_specs[t.owner].append(replace(s,name=t.field))
            tank_values[t.owner].append(v)
    nodes=tuple(replace(n,reservoir=candidate_tank(n.reservoir,tank_specs[n.name],tank_values[n.name])) if tank_specs[n.name] else n for n in network.tanks)
    return replace(network,tanks=nodes,connections=tuple(replace(e,transmissibility=edges[e.key]) if e.key in edges else e for e in network.connections))


def default_spec(name,target):
    value=target.value
    zero_allowed=target.edge is not None or target.field=="m"
    scale=max(value,1e-9 if target.edge else .1) if zero_allowed else value
    lower=0. if zero_allowed else value*.2
    upper=max(value*5,1e-8 if target.edge else 1.) if zero_allowed else value*5
    if target.field=="aquifer.porosity": upper=min(upper,.99)
    if target.field=="aquifer.encroachment_angle": upper=min(upper,360.)
    return MatchParameter(name,target.description,value,lower,upper,target.unit,scale)
