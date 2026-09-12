"""Each edge is calculated once; equal/opposite incidence enforces conservation."""
from math import fsum
from material_balance_studio.domain.validation import finite_value
from .state import ConnectionState


def transfers(network,previous,pressures,dt):
    old={s.name:s.pressure for s in previous.tanks}
    prior={(s.from_tank,s.to_tank):s.cumulative_transfer for s in previous.connections}
    edges=[]
    contributions={t.name:[] for t in network.tanks}
    increments={t.name:[] for t in network.tanks}
    for edge in network.connections:
        a,b=edge.key
        dp=pressures[a]-pressures[b]
        rate=edge.effective_transmissibility*dp
        volume=edge.effective_transmissibility*dt*((old[a]-old[b])+dp)/2
        total=prior[edge.key]+volume
        for label,value in (("rate",rate),("incremental volume",volume),("cumulative volume",total)):
            finite_value(f"{a} to {b} transfer {label}",value)
        edges.append(ConnectionState(a,b,edge.transmissibility,edge.enabled,dp,rate,volume,total))
        contributions[a].append((b,-total))
        contributions[b].append((a,total))
        increments[a].append(-volume)
        increments[b].append(volume)
    return tuple(edges),{n:tuple(c) for n,c in contributions.items()},{n:fsum(c) for n,c in increments.items()}
