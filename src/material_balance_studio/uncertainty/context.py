"""Reuse the exact Phase 4B objective and immutable match scenario snapshot."""
from dataclasses import replace
import numpy as np
from material_balance_studio.domain.models import HistoryRecord
from material_balance_studio.matching.objective import PressureObjective


def context(scenario):
    if not scenario.converged or scenario.final_simulation is None:
        raise ValueError("Select a converged match with valid forward results.")
    history=tuple(HistoryRecord(s.date,s.cumulative,s.observed_pressure) for s in scenario.final_simulation.states)
    fitted=dict(scenario.fitted_parameters)
    specs=tuple(replace(s,initial=fitted[s.name]) for s in scenario.specifications if s.active)
    objective=PressureObjective(scenario.base_tank,history,specs,scenario.observations,scenario.weighting_mode,scenario.default_sigma)
    return history,specs,objective


def evaluate(objective,physical):
    vector=[s.encode(float(v)) for s,v in zip(objective.active,physical)]
    residual=objective(vector)
    return residual if objective.last_result is not None else None


def grid(spec,matched,count):
    if not isinstance(count,int) or not 3<=count<=21:
        raise ValueError("Grid resolution must be an integer between 3 and 21.")
    values=np.geomspace(spec.lower,spec.upper,count) if spec.transformation=="log" else np.linspace(spec.lower,spec.upper,count)
    # Include the actual matched coordinate, not just the nearest coarse node.
    separation=1e-12*(spec.upper-spec.lower)
    return tuple(sorted(set([float(v) for v in values if abs(v-matched)>separation or v in (spec.lower,spec.upper)]+[float(matched)])))


def objective_label(mode):
    return "Σ[(Pcalc−Pobs)/sigma]²" if mode=="Pressure uncertainty weighted" else "Σ[(Pcalc−Pobs)/1 MPa]²"
