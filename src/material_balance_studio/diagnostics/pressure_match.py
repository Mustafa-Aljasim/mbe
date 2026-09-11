"""Unweighted diagnostic metrics; optional sigma normalizes residuals only."""
from .aquifer_comparison import pressure_metrics


def pressure_match(result, observations=()):
    metadata = {m.date:m for m in observations}
    rows = []
    for s in result.states:
        m = metadata.get(str(s.date))
        sigma = m.sigma if m else None
        rows.append(dict(date=str(s.date), calculated=s.pressure, observed=s.observed_pressure,
                         residual=s.pressure_error, sigma=sigma,
                         normalized=s.pressure_error/sigma if s.pressure_error is not None and sigma is not None and sigma>0 else None,
                         source=m.source if m else "", quality=m.quality if m else "", note=m.note if m else ""))
    return dict(metrics=pressure_metrics(result), rows=rows)
