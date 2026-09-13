"""Network adapters for the shared Phase 4C diagnostic mathematics."""
from dataclasses import dataclass, replace
from material_balance_studio.network_matching.objective import NetworkObjective
from material_balance_studio.network_matching.diagnostics import pressure_diagnostics, closure_summary
from .context import evaluate
from .result import GridPoint


def context(scenario):
    if not scenario.converged or scenario.final_simulation is None:
        raise ValueError("Select a converged network match with valid forward results.")
    fitted = dict(scenario.fitted_parameters)
    specs = tuple(replace(s, initial=fitted[s.name]) for s in scenario.specifications if s.active)
    obj = NetworkObjective(scenario.base_network, specs, scenario.observations,
                           scenario.history_plan, scenario.mode, scenario.default_sigma, scenario.settings)
    return specs, obj


@dataclass(frozen=True)
class NetworkPoint(GridPoint):
    tank_metrics: tuple
    closure: tuple
    affected: tuple
    failure_type: str | None


def fixed_point(obj, values):
    residual = evaluate(obj, values)
    record = obj.records[-1]
    valid = residual is not None
    _, tanks, metrics = pressure_diagnostics(obj.last_result, obj.prepared.observations,
                                            obj.mode, obj.default_sigma)
    flags = tuple((s.name, (v-s.lower)/(s.upper-s.lower)<1e-4,
                   (s.upper-v)/(s.upper-s.lower)<1e-4) for s,v in zip(obj.active,values))
    return NetworkPoint(tuple((s.name,float(v)) for s,v in zip(obj.active,values)),
        float(residual@residual) if valid else None, dict(metrics).get('rmse') if valid else None,
        valid, 'FIXED OTHERS', flags, record.failure, 1, tanks,
        record.closure, record.affected, record.failure_type)


def fit_point(fit, parameters):
    valid = fit.final_simulation is not None and fit.final_simulation.converged
    record = fit.evaluations[-1]
    return NetworkPoint(tuple(parameters), fit.final_objective if valid else None,
        dict(fit.final_network_metrics).get('rmse') if valid else None,
        valid, 'CONVERGED' if fit.converged else 'NOT CONVERGED', fit.bound_flags,
        None if fit.converged else fit.termination_message, fit.function_evaluations,
        fit.final_tank_metrics, tuple(closure_summary(fit.final_simulation).items()) if valid else record.closure,
        record.affected, record.failure_type)
