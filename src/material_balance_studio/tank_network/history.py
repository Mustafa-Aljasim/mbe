"""Union timeline with explicit piecewise-linear cumulative stream alignment."""
from bisect import bisect_left
from material_balance_studio.domain.models import CumulativeVolumes,CUMULATIVE_FIELDS,HistoryRecord


def timeline(network):
    return tuple(sorted({h.date for tank in network.tanks for h in tank.history if h.date>network.initial_date}))


def volumes_at(tank,elapsed_seconds):
    records=tank.history
    if not records:
        return CumulativeVolumes(),None
    times=[(r.date-tank.reservoir.initial_date).days*86400. for r in records]
    index=bisect_left(times,elapsed_seconds)
    if index<len(times) and times[index]==elapsed_seconds:
        return records[index].cumulative,records[index].observed_pressure
    if index==len(records):
        return records[-1].cumulative,None
    before=records[index-1].cumulative if index else CumulativeVolumes()
    t0=times[index-1] if index else 0.
    after=records[index].cumulative
    fraction=(elapsed_seconds-t0)/(times[index]-t0)
    # Convex interpolation of monotonically increasing surface totals.
    return CumulativeVolumes(**{n:getattr(before,n)+(getattr(after,n)-getattr(before,n))*fraction for n in CUMULATIVE_FIELDS}),None


def allocate_field_history(history,shares,streams):
    """Explicit opt-in per stream; unselected streams/observations are not allocated."""
    from math import fsum
    from material_balance_studio.domain.validation import validate_history,finite_value
    validate_history(history)
    if not shares or len(set(streams))!=len(streams) or not set(streams)<=set(CUMULATIVE_FIELDS) or not streams:
        raise ValueError("Choose unique valid cumulative streams and nonempty tank shares.")
    for value in shares.values():
        finite_value("Allocation share",value,nonnegative=True)
    if abs(fsum(shares.values())-1.)>1e-10:
        raise ValueError("Fixed allocation shares must sum to 1 within 1e-10.")
    # No silent normalization. Shares retain their submitted values.
    return {name:tuple(HistoryRecord(h.date,CumulativeVolumes(**{s:getattr(h.cumulative,s)*share for s in streams})) for h in history)
            for name,share in shares.items()}
