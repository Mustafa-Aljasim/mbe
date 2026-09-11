from dataclasses import dataclass
from material_balance_studio.units.conversions import from_si
from ..correlations.base import CorrelationInput
from ..correlations.gas_pseudocritical import PSEUDO_CRITICAL_METHODS


@dataclass(frozen=True)
class Applicability:
    status: str
    reasons: tuple[str, ...] = ()


def applicability(correlation, model, pressures=None) -> Applicability:
    """Published population limits are cautions, mathematical/composition domains fail."""
    f = model.fluid
    ps = tuple(pressures) if pressures is not None else model.pressure_bounds
    if f.gas_basis != "sweet":
        return Applicability("NOT_APPLICABLE", ("No sour-fluid composition corrections are implemented.",))
    context = CorrelationInput(f, model.pb, model.rsb)
    try:
        pseudo = model.pseudo_critical_properties
        pc = from_si(pseudo.pressure, "pressure", "FIELD")
        tr = f.temperature / pseudo.temperature
    except ValueError as exc:
        if correlation.property in ("z", "bg"):
            return Applicability("NOT_APPLICABLE", (str(exc),))
        pc, tr = 1, 0
    if correlation.property in ("z", "bg") and (pc <= 0 or tr <= 1):
        return Applicability("NOT_APPLICABLE", ("DAK/Sutton requires positive pseudo-criticals and Tr > 1.",))
    if context.t <= 0:
        return Applicability("NOT_APPLICABLE", ("This implementation requires temperature > 0°F.",))
    values = {"api": (f.api,), "gas_gravity": (f.gas_gravity,), "t": (context.t,),
              "rsb": (context.r,), "sg_oil": (f.sg_oil,), "pb": (context.p,),
              "p": tuple(from_si(p, "pressure", "FIELD") for p in ps),
              "tr": (tr,), "pr": tuple(from_si(p, "pressure", "FIELD") / pc for p in ps) if pc > 0 else (0,)}
    reasons = []
    if correlation.property == "z":
        method = PSEUDO_CRITICAL_METHODS[model.pseudo_critical_method]
        lo, hi = method.gas_gravity_range
        if not lo <= f.gas_gravity <= hi:
            reasons.append(f"{method.name} pseudo-critical method gas gravity outside {lo:g}–{hi:g}.")
    for bound in correlation.ranges:
        if any(not bound.minimum <= v <= bound.maximum for v in values[bound.key]):
            reasons.append(f"{bound.key} outside {bound.minimum:g}–{bound.maximum:g} {bound.unit}.")
    if correlation.key.startswith("vasquez_beggs"):
        reasons.append("Verify gas gravity on the 100-psig separator basis; no correction is applied.")
    if correlation.property == "bg":
        from ..correlations.registry import REGISTRY
        inherited = applicability(REGISTRY[model.selection.z], model, ps)
        if inherited.status == "NOT_APPLICABLE":
            return inherited
        reasons.extend(inherited.reasons)
    return Applicability("CAUTION" if reasons else "VALID", tuple(reasons))
