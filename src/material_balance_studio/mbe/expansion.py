"""Expansion per stock-tank volume of original oil; no embedded field factors."""

from material_balance_studio.domain.models import ExpansionTerms, PVTProperties, ReservoirTank


def oil_solution_gas_expansion(current: PVTProperties, initial: PVTProperties) -> float:
    """Eo = (Bo - Boi) + (Rsi - Rs) Bg."""
    return current.bo - initial.bo + (initial.rs - current.rs) * current.bg


def gas_cap_expansion(current: PVTProperties, initial: PVTProperties) -> float:
    """Eg = Boi (Bg / Bgi - 1); multiply by m when assembling Et."""
    return initial.bo * (current.bg / initial.bg - 1.0)


def rock_connate_water_expansion(pressure: float, tank: ReservoirTank,
                                initial: PVTProperties) -> float:
    """Efw includes the gas-cap pore volume through (1 + m).

    First-order constant-compressibility approximation; Swc applies to the tank.
    """
    return (initial.bo * (1.0 + tank.m) * (tank.cf + tank.cw * tank.swc)
            / (1.0 - tank.swc) * (tank.initial_pressure - pressure))


def expansion_terms(current: PVTProperties, initial: PVTProperties,
                    tank: ReservoirTank) -> ExpansionTerms:
    eg = gas_cap_expansion(current, initial)
    return ExpansionTerms(
        eo=oil_solution_gas_expansion(current, initial), eg=eg, meg=tank.m * eg,
        efw=rock_connate_water_expansion(current.pressure, tank, initial),
    )
