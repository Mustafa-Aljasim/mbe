"""Reservoir setup adapter: convert boundary values before constructing the domain."""

from datetime import date

from material_balance_studio.domain.models import ReservoirTank
from material_balance_studio.pvt.base import PVTModel
from material_balance_studio.units.conversions import UnitSystem, to_si


def tank_from_inputs(*, initial_date: date, initial_pressure: float, oil_in_place: float,
                     swc: float, cf: float, cw: float, m: float, pvt_model: PVTModel,
                     units: UnitSystem | str = UnitSystem.SI) -> ReservoirTank:
    """Convert pressure, stock-tank oil volume and compressibilities to internal SI."""
    return ReservoirTank(
        initial_date=initial_date, initial_pressure=to_si(initial_pressure, "pressure", units),
        oil_in_place=to_si(oil_in_place, "liquid_volume", units), swc=swc,
        cf=to_si(cf, "compressibility", units), cw=to_si(cw, "compressibility", units),
        m=m, pvt_model=pvt_model,
    )
