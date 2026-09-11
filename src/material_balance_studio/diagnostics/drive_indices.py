"""Additive gross-production-basis drive indices from accepted balance terms."""

def drive_indices(state):
    b = state.balance
    supports = dict(DDI=b.oil_expansion_support, SDI=b.gas_cap_expansion_support,
                    CDI=b.rock_water_expansion_support, WDI=b.aquifer_support,
                    IDI=b.withdrawal.water_injection+b.withdrawal.gas_injection)
    denominator = b.withdrawal.production
    values = {k:v/denominator if denominator > 1e-12 else None for k,v in supports.items()}
    total = sum(values.values()) if all(v is not None for v in values.values()) else None
    return dict(values=values, supports=supports, denominator=denominator, total=total,
                closure_error=total-1 if total is not None else None,
                stackable=total is not None and abs(total-1)<1e-8 and all(v>=0 for v in values.values()))
