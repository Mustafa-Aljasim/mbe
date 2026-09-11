"""Physical grids evaluated by the unchanged history-match objective."""
import numpy as np
from .context import context,evaluate,grid,objective_label
from .result import GridPoint,SurfaceResult,scenario_identity


def fixed_point(obj,values):
    residual=evaluate(obj,values)
    flags=tuple((s.name,(v-s.lower)/(s.upper-s.lower)<1e-4,(s.upper-v)/(s.upper-s.lower)<1e-4)
                for s,v in zip(obj.active,values))
    return GridPoint(tuple((s.name,float(v)) for s,v in zip(obj.active,values)),
        None if residual is None else float(residual@residual),
        None if residual is None else float(np.sqrt(np.mean((residual*obj.scales)**2))),
        residual is not None,"FIXED OTHERS",flags,obj.records[-1].failure,1)


def objective_surface(scenario,parameter_x,parameter_y,count=7):
    _,specs,obj=context(scenario)
    names=[s.name for s in specs]
    if parameter_x==parameter_y or parameter_x not in names or parameter_y not in names:
        raise ValueError("Select two different active match parameters.")
    ix,iy=names.index(parameter_x),names.index(parameter_y)
    values=[s.initial for s in specs]
    xs,ys=grid(specs[ix],values[ix],count),grid(specs[iy],values[iy],count)
    points=[]
    for y in ys:
        for x in xs:
            candidate=values.copy()
            candidate[ix],candidate[iy]=x,y
            points.append(fixed_point(obj,candidate))
    failed=sum(not p.valid for p in points)
    warnings=["Conditional surface: all other parameters held at matched values; a coarse grid can miss a narrow minimum."]
    if failed:
        warnings.append(f"{failed} failed forward nodes are invalid regions, not objective minima.")
    initial={s.name:s.initial for s in scenario.specifications}
    return SurfaceResult(scenario_identity(scenario),(parameter_x,parameter_y),xs,ys,tuple(points),
        objective_label(scenario.weighting_mode),(initial[parameter_x],initial[parameter_y]),
        (values[ix],values[iy]),tuple(warnings),failed)
