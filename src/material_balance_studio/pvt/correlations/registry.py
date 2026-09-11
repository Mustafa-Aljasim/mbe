"""Property registry. Public evaluators accept/return canonical SI (T in K).

The primitive equations use psia, °F, scf/STB, cP. See pvt_phase_2.md for
equations, original publications, reference calculations and validity sources.
"""
from material_balance_studio.units.conversions import to_si
from .base import Correlation, ValidityRange as R
from .bubble_point import standing_pb, vasquez_beggs_pb, glaso_pb
from .solution_gor import standing_rs, vasquez_beggs_rs, glaso_rs
from .oil_fvf import standing_bo, vasquez_beggs_bo, glaso_bo, vasquez_beggs_compressibility
from .oil_viscosity import beggs_robinson_live
from .gas_z import dak_z
from .gas_fvf import real_gas_fvf
from .water import mccain_bw

REGISTRY: dict[str, Correlation] = {}


def register(correlation: Correlation) -> None:
    if correlation.key in REGISTRY:
        raise ValueError(f"Duplicate correlation: {correlation.key}")
    REGISTRY[correlation.key] = correlation


def candidates(property: str) -> tuple[Correlation, ...]:
    return tuple(c for c in REGISTRY.values() if c.property == property)


def _ranges(api, gg, t, rs, pb):
    return tuple(R(key, *bounds, unit) for key, bounds, unit in (
        ("api", api, "°API"), ("gas_gravity", gg, "air=1"), ("t", t, "°F"),
        ("rsb", rs, "scf/STB"), ("pb", pb, "psia")))


FAMILIES = (
    ("standing", "Standing", standing_pb, standing_rs,
     lambda i: standing_bo(i.r, i.t, i.fluid.sg_oil, i.fluid.gas_gravity),
     "https://doi.org/10.2118/947275-G",
     _ranges((16.5, 63.8), (.59, .95), (100, 258), (20, 1425), (130, 7000))),
    ("vasquez_beggs", "Vasquez–Beggs", vasquez_beggs_pb, vasquez_beggs_rs,
     lambda i: vasquez_beggs_bo(i.r, i.t, i.fluid.api, i.fluid.gas_gravity),
     "https://doi.org/10.2118/6719-PA",
     _ranges((15.3, 59.5), (.511, 1.351), (70, 295), (0, 2199), (15, 6055))),
    ("glaso", "Glasø", glaso_pb, glaso_rs,
     lambda i: glaso_bo(i.r, i.t, i.fluid.sg_oil, i.fluid.gas_gravity),
     "https://doi.org/10.2118/8016-PA",
     _ranges((22.3, 48.1), (.65, 1.273), (80, 280), (90, 2637), (165, 7142))),
)
for key, name, pb_fn, rs_fn, bo_fn, source, ranges in FAMILIES:
    notes = ("Gas gravity must use the published 100-psig separator basis; no separator correction is applied."
             if key == "vasquez_beggs" else "Sweet black oil; no composition correction.")
    register(Correlation(f"{key}_pb", "pb", name,
        lambda i, fn=pb_fn: to_si(fn(i.r, i.t, i.fluid.api, i.fluid.gas_gravity), "pressure", "FIELD"),
        source, f"{key} saturated bubble pressure; FIELD primitive → Pa", ranges, notes))
    register(Correlation(f"{key}_rs", "rs", name,
        lambda i, fn=rs_fn: to_si(fn(i.p, i.t, i.fluid.api, i.fluid.gas_gravity), "rs", "FIELD"),
        source, f"{key} saturated solution GOR; FIELD primitive → m³/m³", ranges, notes))
    register(Correlation(f"{key}_bo", "bo", name, bo_fn, source,
        f"{key} saturated oil FVF; FIELD primitive → m³/m³", ranges, notes))

register(Correlation("beggs_robinson_mu", "oil_viscosity", "Beggs–Robinson / Vasquez–Beggs",
    lambda i: to_si(beggs_robinson_live(i.r, i.t, i.fluid.sg_oil), "viscosity", "FIELD"),
    "https://doi.org/10.2118/5434-PA",
    "SG variant dead/live oil, VB undersaturated multiplier; cP → Pa·s",
    (R("t", 70, 295, "°F"), R("rsb", 20, 2070, "scf/STB"), R("sg_oil", .75, .96, "water=1"),
     R("p", 0, 5250, "psia"))))


def _z(i):
    return dak_z(*i.reduced_state)


# Historical key retained for Phase 2 configurations; it no longer identifies
# a pseudo-critical method. The model's pseudo_critical_method selects that.
register(Correlation("sutton_dak_z", "z", "Dranchuk–Abou-Kassem", _z,
    "https://doi.org/10.2118/75-03-03",
    "DAK implicit density EOS using independently supplied Pr and Tr; dimensionless",
    (R("tr", 1, 3, "dimensionless"), R("pr", .2, 30, "dimensionless")),
    "Sweet gas only; Tr must exceed 1. No sour-gas correction."))
register(Correlation("real_gas_bg", "bg", "Real-gas law (selected z)",
    lambda i: real_gas_fvf(i.pressure, i.fluid.temperature, i.z if i.z is not None else _z(i)),
    "https://manual.whitson.com/modules/bot/", "Bg=Z T Psc/(P Tsc); Zsc=1, 60°F and 101325 Pa; m³/m³",
    notes="Inherits selected z applicability and fixed standard conditions."))
register(Correlation("mccain_bw", "bw", "McCain", lambda i: mccain_bw(i.p, i.t),
    "https://doi.org/10.2118/18571-PA",
    "Bw=(1+ΔVp)(1+ΔVT); psia/°F polynomial → m³/m³",
    (R("p", 0, 5000, "psia"), R("t", 32, 260, "°F")),
    "Salinity and water SG are stored but not explicit terms; 32°F lower bound is an implementation guard."))
register(Correlation("vasquez_beggs_co", "co", "Vasquez–Beggs (SI reformulation)", vasquez_beggs_compressibility,
    "https://doi.org/10.2118/6719-PA",
    "co=[28.1Rs+30.6TK−1180SGg+1784/SGo−10910]/(1e5 pPa); 1/Pa",
    FAMILIES[1][-1], "Used only above Pb, integrated exactly with constant Rsb and T."))
