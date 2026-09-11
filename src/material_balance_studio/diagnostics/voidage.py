"""Decompose validated withdrawal; never introduce a second gas/withdrawal formula."""
from material_balance_studio.mbe.withdrawal import withdrawal_terms


def decomposition(volumes, pvt, terms):
    oil, water = volumes.np*pvt.bo, volumes.wp*pvt.bw
    return dict(oil=oil, free_gas=terms.production-oil-water, water=water,
                produced=terms.production, water_injected=terms.water_injection,
                gas_injected=terms.gas_injection,
                injected=terms.water_injection+terms.gas_injection, net=terms.net)


def voidage(state):
    """Interval surface increments use endpoint PVT, not differences of repriced totals."""
    return {"cumulative": decomposition(state.cumulative,state.pvt,state.balance.withdrawal),
            "interval": decomposition(state.increments,state.pvt,withdrawal_terms(state.increments,state.pvt))}
