"""Physical parameter registry and immutable candidate construction."""
from dataclasses import dataclass,replace
from math import log,exp,isfinite


def parameter_registry(tank):
    registry={"oil_in_place":("OOIP N","oil_volume"),"m":("Initial gas-cap ratio m","dimensionless")}
    key=tank.aquifer.key if tank.aquifer else "none"
    applicable={"none":(),"pot":("capacity",),"schilthuis":("productivity_index",),
                "fetkovich":("productivity_index","initial_water_volume")}
    transient=("permeability","inner_radius","thickness","porosity","encroachment_angle")
    names=applicable.get(key,transient if key in ("carter_tracy","van_everdingen_hurst","modified_van_everdingen_hurst") else ())
    definitions={"capacity":("Aquifer capacity","aquifer_capacity"),"productivity_index":("Aquifer productivity J","aquifer_productivity"),
                 "initial_water_volume":("Connected aquifer water volume","reservoir_volume"),"permeability":("Aquifer permeability","permeability"),
                 "inner_radius":("Reservoir/aquifer inner radius","length"),"thickness":("Aquifer thickness","length"),
                 "porosity":("Aquifer porosity","dimensionless"),"encroachment_angle":("Encroachment angle","angle")}
    registry.update({"aquifer."+name:definitions[name] for name in names})
    return registry


def parameter_value(tank,name):
    return getattr(tank.aquifer,name.split(".")[1]) if name.startswith("aquifer.") else getattr(tank,name)


@dataclass(frozen=True)
class MatchParameter:
    name: str
    description: str
    initial: float
    lower: float
    upper: float
    unit: str
    scale: float
    transformation: str = "linear"
    active: bool = True

    def __post_init__(self):
        if not all(isfinite(v) for v in (self.initial,self.lower,self.upper,self.scale)):
            raise ValueError("Parameter values, scales and bounds must be finite.")
        if self.scale<=0 or not self.lower<=self.initial<=self.upper or self.lower>=self.upper:
            raise ValueError("Require positive scale and lower ≤ initial ≤ upper with distinct bounds.")
        if self.transformation not in ("linear","log"):
            raise ValueError("Transformation must be linear or log.")
        if self.transformation=="log" and self.lower<=0:
            raise ValueError("Log parameters require strictly positive bounds.")

    def encode(self,value):
        return log(value/self.scale) if self.transformation=="log" else value/self.scale

    def decode(self,value):
        if value==self.encode(self.lower):
            return self.lower
        if value==self.encode(self.upper):
            return self.upper
        return self.scale*exp(value) if self.transformation=="log" else self.scale*value


def validate_parameters(tank,specs):
    registry=parameter_registry(tank)
    if len({s.name for s in specs})!=len(specs):
        raise ValueError("Duplicate match parameter.")
    for s in specs:
        if s.name not in registry or s.unit!=registry[s.name][1]:
            raise ValueError(f"Parameter or unit not approved for this model: {s.name}.")
        if not s.active:
            continue
        if (s.name=="m" and s.lower<0) or (s.name!="m" and s.lower<=0):
            raise ValueError("Bounds must remain physical (positive, or nonnegative for m).")
        if s.name=="aquifer.porosity" and s.upper>=1:
            raise ValueError("Porosity upper bound must be below one.")
        if s.name=="aquifer.encroachment_angle" and s.upper>360:
            raise ValueError("Encroachment upper bound cannot exceed 360 degrees.")
    if not any(s.active for s in specs):
        raise ValueError("Select at least one active parameter.")


def candidate_tank(tank,specs,physical):
    reservoir,aquifer={},{}
    active=[s for s in specs if s.active]
    if len(active)!=len(physical):
        raise ValueError("Parameter vector length mismatch.")
    for spec,value in zip(active,physical):
        if not isfinite(value) or not spec.lower<=value<=spec.upper:
            raise ValueError(f"{spec.name} is outside explicit bounds.")
        if spec.name.startswith("aquifer."):
            aquifer[spec.name.split(".")[1]]=float(value)
        else:
            reservoir[spec.name]=float(value)
    if aquifer:
        reservoir["aquifer"]=replace(tank.aquifer,**aquifer)
    return replace(tank,**reservoir)
