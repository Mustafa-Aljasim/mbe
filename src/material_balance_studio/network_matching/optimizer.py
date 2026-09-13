"""Bounded multi-tank matching through the unchanged Phase 5A forward engine."""
import hashlib,json
from dataclasses import asdict
from scipy.optimize import least_squares
from .objective import NetworkObjective
from .diagnostics import pressure_diagnostics,allocation_crosscheck
from .parameters import registry,candidate_network
from .result import NetworkMatchResult


def history_match_network(network,specs,observations=None,plan=None,mode="Unweighted",default_sigma=None,
                          settings=None,max_nfev=60,scenario_name="Network match"):
    if type(max_nfev) is not int or not 1<=max_nfev<=1000:
        raise ValueError("max_nfev must be an integer from 1 to 1000.")
    objective=NetworkObjective(network,specs,observations,plan,mode,default_sigma,settings)
    active=objective.active
    initial=objective([s.encode(s.initial) for s in active])
    initial_sim=objective.last_result
    initial_rows,initial_tank,initial_network=pressure_diagnostics(initial_sim,objective.prepared.observations,mode,default_sigma)
    fit=least_squares(objective,[s.encode(s.initial) for s in active],bounds=([s.encode(s.lower) for s in active],[s.encode(s.upper) for s in active]),
                      method="trf",jac="3-point",x_scale=1.,ftol=1e-10,xtol=1e-10,gtol=1e-10,max_nfev=max_nfev)
    final=objective(fit.x)
    final_sim,fitted_network=objective.last_result,objective.last_network
    converged=bool(fit.success and final_sim is not None)
    values=tuple((s.name,s.decode(float(x))) for s,x in zip(active,fit.x))
    flags=tuple((s.name,(v-s.lower)/(s.upper-s.lower)<1e-4,(s.upper-v)/(s.upper-s.lower)<1e-4) for s,(_,v) in zip(active,values))
    final_rows,final_tank,final_network=pressure_diagnostics(final_sim,objective.prepared.observations,mode,default_sigma)
    warnings=list(objective.warnings)+list(objective.prepared.warnings)
    warnings.append("A good network pressure fit does not establish unique parameters; multi-tank identifiability is deferred to Phase 5C.")
    for name,lo,hi in flags:
        if lo or hi:
            warnings.append(f"{name}: {'lower' if lo else 'upper'} bound limited; weak communication/parameter constraint may remain.")
    failed=sum(not e.valid for e in objective.records)
    if failed:
        warnings.append(f"{failed} failed network candidates received fixed finite penalties; inspect failure audit.")
    if not converged:
        warnings.append("Network match not converged: "+str(fit.message)+(": "+str(objective.records[-1].failure) if final_sim is None else ""))
    for n,items in final_tank:
        old=dict(dict(initial_tank).get(n,()))
        new=dict(items)
        if old.get("rmse") is not None and new.get("rmse") is not None and new["rmse"]>old["rmse"]+1.:
            warnings.append(f"{n}: pressure RMSE worsened despite network optimization; review tank-level trade-off.")
    warnings.extend(allocation_crosscheck(objective.prepared.plan,final_rows,network))
    closures=[dict(e.closure) for e in objective.records if e.valid]
    maxrel=max((c["maximum_relative"] for c in closures),default=None)
    maxabs=max((c["maximum_absolute"] for c in closures),default=None)
    maxtransfer=max((c["transfer_error"] for c in closures),default=None)
    margin=min((c["closure_margin"] for c in closures),default=None)
    if maxrel is not None and maxrel>=8e-9:
        warnings.append("NEAR TOLERANCE: at least one valid candidate used ≥80% of the 1e-8 relative closure limit; tolerance was not relaxed.")
    identity=hashlib.sha256(json.dumps([scenario_name,objective.prepared.fingerprint,[asdict(s) for s in objective.specs],asdict(objective.settings),values,mode,default_sigma],sort_keys=True,default=str).encode()).hexdigest()
    return NetworkMatchResult(scenario_name,identity,network,fitted_network,objective.prepared.plan,objective.prepared.fingerprint,
        objective.prepared.observations,objective.specs,objective.settings,mode,default_sigma,values,converged,int(fit.status),str(fit.message),len(objective.records),fit.nfev,
        float(initial@initial),float(final@final),initial_sim,final_sim,initial_rows,final_rows,initial_tank,final_tank,initial_network,final_network,
        flags,tuple(objective.records),maxrel,maxabs,maxtransfer,margin,failed,"FAIL" if not converged else "CAUTION" if len(warnings)>1 else "PASS",tuple(warnings),objective.prepared.audit)


def apply_matched_parameters(current,result):
    if not result.converged:
        raise ValueError("Cannot apply an unsuccessful network match.")
    def topology(n):
        return (tuple((t.name,t.reservoir.aquifer.key if t.reservoir.aquifer else "none") for t in n.tanks),tuple((e.key,e.enabled) for e in n.connections))
    if topology(current)!=topology(result.base_network):
        raise ValueError("Current network topology/aquifer selection differs from the scenario; cannot apply parameters.")
    return candidate_network(current,result.specifications,[dict(result.fitted_parameters)[s.name] for s in result.specifications if s.active])
