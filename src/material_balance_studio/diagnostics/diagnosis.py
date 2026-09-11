"""Compose read-only reservoir diagnostics; no fitting or forward-state mutation."""
from .voidage import voidage
from .vrr import vrr
from .drive_indices import drive_indices
from .pressure_match import pressure_match
from .havlena_odeh import observed_components,havlena_odeh
from .campbell_cole import campbell


def diagnose(result,tank,history,observations=()):
    components=observed_components(tank,history)
    return dict(rows=[dict(date=str(s.date),voidage=voidage(s),vrr=vrr(s),drive=drive_indices(s)) for s in result.states],
                components=components,ho=havlena_odeh(components,tank),
                ho_gas=havlena_odeh(components,tank,"gas_cap") if tank.m>0 else None,
                campbell=campbell(components),pressure=pressure_match(result,observations))
