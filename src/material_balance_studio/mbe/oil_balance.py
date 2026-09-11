"""Independent balance API shared by the solver, inspector, and tests."""

from material_balance_studio.domain.models import BalanceTerms, CumulativeVolumes, ReservoirTank
from material_balance_studio.domain.validation import finite_value
from .expansion import expansion_terms
from .withdrawal import withdrawal_terms


def evaluate_balance(pressure: float, tank: ReservoirTank, cumulative: CumulativeVolumes,
                     normalization_floor: float = 1.0, *, aquifer_influx: float = 0.0) -> BalanceTerms:
    """Evaluate cumulative Fnet - N Et - We, with volumes in reservoir m³."""
    finite_value("normalization_floor", normalization_floor, positive=True)
    finite_value("cumulative aquifer influx", aquifer_influx)
    initial = tank.pvt_model.properties_at_pressure(tank.initial_pressure)
    current = tank.pvt_model.properties_at_pressure(pressure)
    withdrawal = withdrawal_terms(cumulative, current)
    expansion = expansion_terms(current, initial, tank)
    total_support = tank.oil_in_place * expansion.et
    residual = withdrawal.net - total_support - aquifer_influx
    absolute = abs(residual)
    terms = BalanceTerms(
        withdrawal=withdrawal, expansion=expansion,
        oil_expansion_support=tank.oil_in_place * expansion.eo,
        gas_cap_expansion_support=tank.oil_in_place * expansion.meg,
        rock_water_expansion_support=tank.oil_in_place * expansion.efw,
        total_expansion_support=total_support, residual=residual,
        absolute_residual=absolute,
        relative_residual=absolute / max(abs(withdrawal.net), normalization_floor),
        aquifer_support=aquifer_influx,
    )
    finite_value("material-balance residual", terms.residual)
    finite_value("relative material-balance residual", terms.relative_residual)
    return terms


def balance_residual(pressure: float, tank: ReservoirTank,
                     cumulative: CumulativeVolumes) -> float:
    """Signed independent root function: Residual(p) = Fnet(p) - N Et(p)."""
    return evaluate_balance(pressure, tank, cumulative).residual
