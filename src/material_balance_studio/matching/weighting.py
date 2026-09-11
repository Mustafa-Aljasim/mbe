"""Exact-date observation selection and explicit sigma policy."""
from dataclasses import dataclass
from math import isfinite
import numpy as np


@dataclass(frozen=True)
class MatchObservation:
    date: object
    pressure: float
    source: str = ""
    sigma: float | None = None
    qc_flag: str = ""
    include_in_match: bool = True
    note: str = ""


def observations_from_history(history,metadata=()):
    meta={m.date:m for m in metadata}
    rows=[]
    for h in history:
        if h.observed_pressure is None:
            continue
        m=meta.get(str(h.date))
        rows.append(MatchObservation(h.date,h.observed_pressure,m.source if m else "",m.sigma if m else None,
                                     m.quality if m else "",True,m.note if m else ""))
    return tuple(rows)


def residual_scales(observations,mode,default_sigma=None):
    if mode not in ("Unweighted","Pressure uncertainty weighted"):
        raise ValueError("Unknown weighting mode.")
    if default_sigma is not None and (not isfinite(default_sigma) or default_sigma<=0):
        raise ValueError("Explicit default sigma must be positive and finite.")
    if mode=="Unweighted":
        return np.full(len(observations),1e6)  # Pa: fixed scale independent of display units
    values=[]
    for o in observations:
        sigma=o.sigma if o.sigma is not None and isfinite(o.sigma) and o.sigma>0 else default_sigma
        if sigma is None:
            raise ValueError(f"Missing/invalid sigma at {o.date}; exclude this observation or explicitly supply default sigma.")
        if sigma<1e-12:
            raise ValueError("Sigma below 1e-12 Pa is numerically unsupported.")
        values.append(sigma)
    return np.asarray(values)
