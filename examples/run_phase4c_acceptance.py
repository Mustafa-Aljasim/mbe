"""Reproducible engineering evidence; verifies all prior numerical outputs."""
from dataclasses import asdict
from pathlib import Path
import hashlib,json,runpy
import numpy as np
from material_balance_studio.uncertainty import analyze,parameter_profile,objective_surface

ROOT=Path(__file__).resolve().parents[1]
synthetic_match=runpy.run_path(str(ROOT/"examples/identifiability_benchmarks.py"))["synthetic_match"]


def evidence():
    previous=json.loads((ROOT/"docs/phase_4b_numerical_results.json").read_text())
    actual=json.loads(json.dumps(runpy.run_path(str(ROOT/"examples/run_phase4b_acceptance.py"))["evidence"]()))
    assert actual==previous,"Protected Phase 4B numerical results changed"
    hashes=json.loads((ROOT/"docs/phase_4c_matching_baseline.json").read_text())
    assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==v for p,v in hashes.items())
    records={}
    fits={kind:synthetic_match(kind)[2] for kind in ("single","exact","near","zero")}
    fits["bound"]=synthetic_match(upper=9e5)[2]
    for kind,fit in fits.items():
        r=analyze(fit,True)
        records[kind]=dict(fitted=dict(fit.fitted_parameters),rmse_pa=fit.final_metrics["rmse"],rank=r.rank,
            singular_values=r.singular_values,condition=r.condition,correlation=r.correlation,statuses=r.statuses,
            standard_errors=r.standard_errors,local_intervals=r.intervals,schemes=r.schemes,warnings=r.warnings,
            derivative_relative_changes=r.derivative_errors)
    for kind in ("single","exact","near"):
        fit=fits[kind]
        before=asdict(fit)
        profile=parameter_profile(fit,"oil_in_place",7,True,statistical_assumptions=True,plausible_delta=1. if kind=="exact" else None)
        records[kind]["profile"]=asdict(profile)
        assert asdict(fit)==before
    surface=objective_surface(fits["exact"],"oil_in_place","aquifer.capacity",5)
    records["exact"]["surface"]=asdict(surface)
    records["exact"]["distinct_valley_nodes"]=sum(p.valid and p.objective<1e-8 for p in surface.points)
    assert analyze(fits["near"])==analyze(fits["near"])
    return dict(protected_forward_and_phase4b_outputs_identical=True,protected_matching_hashes=hashes,
        previous_maximum_valid_candidate_mbe_relative=actual["maximum_valid_candidate_mbe_relative"],
        scenario_immutability=True,repeated_sensitivity_identical=True,benchmarks=records)


if __name__=="__main__":
    result=evidence()
    (ROOT/"docs/phase_4c_numerical_results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({k:dict(rank=v["rank"],correlation=v["correlation"],rmse_pa=v["rmse_pa"],intervals=v["local_intervals"],
        plausible_ranges=v.get("profile",{}).get("plausible_ranges")) for k,v in result["benchmarks"].items()},indent=2))
