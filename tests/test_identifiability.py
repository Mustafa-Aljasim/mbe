"""Independent analytical acceptance of local and deterministic uncertainty."""
from dataclasses import asdict,replace
from pathlib import Path
import hashlib,json,runpy
import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest
from material_balance_studio.uncertainty import analyze,parameter_profile,objective_surface
from material_balance_studio.uncertainty.jacobian import matrix_diagnostics
from material_balance_studio.uncertainty.confidence import grid_ranges,profile_threshold
from material_balance_studio.uncertainty.result import GridPoint,scenario_identity
from material_balance_studio.matching import history_match,observations_from_history
from material_balance_studio.presentation.uncertainty import sensitivity_frame,point_frame
from material_balance_studio.units.display import to_display,from_display

ROOT=Path(__file__).resolve().parents[1]
synthetic_match=runpy.run_path(str(ROOT/"examples/identifiability_benchmarks.py"))["synthetic_match"]


@pytest.fixture(scope="module")
def single():
    return synthetic_match()[2]


@pytest.fixture(scope="module")
def exact():
    return synthetic_match("exact")[2]


@pytest.fixture(scope="module")
def near():
    return synthetic_match("near")[2]


def test_single_analytical_sensitivity_jacobian_confidence(single):
    r=analyze(single,True)
    n=dict(single.fitted_parameters)["oil_in_place"]
    npv=np.array([s.cumulative.np for s in single.final_simulation.states])
    # Invert Np*(Boi+c*dp)=N*c*dp; differentiate P=Pi-Boi*Np/[c*(N-Np)].
    oracle=1.2*npv/(4e-9*(n-npv)**2)
    assert np.asarray(r.sensitivity)[:,0]==pytest.approx(oracle,rel=2e-6)
    assert np.asarray(r.residual_jacobian)[:,0]==pytest.approx(oracle/1e5,rel=2e-6)
    assert np.asarray(r.scaled_sensitivity)[:,0]==pytest.approx(oracle*1e6/30e6,rel=2e-6)
    assert r.rank==1 and r.condition==1 and r.statuses[0][1]=="WELL CONSTRAINED"
    assert r.intervals[0][0]<1e6<r.intervals[0][1]
    expected_se=1/np.sqrt(np.sum((oracle/1e5)**2))
    assert r.standard_errors[0]==pytest.approx(expected_se,rel=2e-6)
    profile=parameter_profile(single,"oil_in_place",9,True,statistical_assumptions=True)
    minimum=min(profile.points,key=lambda p:p.objective)
    assert dict(minimum.parameters)["oil_in_place"]==pytest.approx(1e6)
    assert any(a-1e-3<=1e6<=b+1e-3 for a,b in profile.plausible_ranges)
    assert len(profile.threshold_brackets)==2
    assert any("width is unresolved" in w for w in profile.warnings)
    assert profile.points[0].objective>minimum.objective+100


def test_analytical_correlation_sign_magnitude():
    d=matrix_diagnostics([[1,1],[1,2],[1,3]])
    # (J'J)^-1 = [[14,-6],[-6,3]] / 6.
    assert d["geometry"]==pytest.approx(np.array([[14,-6],[-6,3]])/6)
    assert d["correlation"][0,1]==pytest.approx(-6/np.sqrt(42))
    assert np.diag(d["correlation"])==pytest.approx([1,1])


@pytest.mark.parametrize("bad",[[[1,float('nan')]],[],[[float('inf')]]])
def test_nonfinite_or_empty_jacobian_rejected(bad):
    with pytest.raises(ValueError,match="finite nonempty"):
        matrix_diagnostics(bad)


def test_exact_rank_deficiency_withholds_marginal_uncertainty(exact):
    r=analyze(exact,True)
    assert exact.final_metrics["rmse"]<1e-5 and r.rank==1
    assert r.correlation is None and r.covariance is None and r.intervals==(None,None)
    assert all(s[1]=="NON-IDENTIFIABLE" for s in r.statuses)
    assert abs(r.column_coupling[0][1])>1-1e-10
    profile=parameter_profile(exact,"oil_in_place",3,True,statistical_assumptions=True)
    assert profile.threshold is None and "regularity" in profile.threshold_method


def test_near_confounded_correlation_and_wider_profile(near,single):
    r=analyze(near,True)
    assert near.final_metrics["rmse"]<1e-5 and r.rank==2
    assert r.correlation[0][1]<-.99
    assert all(s[1]=="POORLY CONSTRAINED" for s in r.statuses)
    a=parameter_profile(near,"oil_in_place",7,True,statistical_assumptions=True)
    b=parameter_profile(single,"oil_in_place",7,True,statistical_assumptions=True)
    assert max(y-x for x,y in a.plausible_ranges)>max(y-x for x,y in b.plausible_ranges)
    assert len(a.plausible_ranges)==1


def test_exact_reoptimization_matches_independent_compensation_and_is_pure(exact):
    before=asdict(exact)
    prof=parameter_profile(exact,"oil_in_place",5,True,plausible_delta=1.)
    fixed=parameter_profile(exact,"oil_in_place",5,False,plausible_delta=1.)
    for p in prof.points:
        n,c=dict(p.parameters).values()
        assert p.valid and p.optimizer_status=="CONVERGED"
        assert c==pytest.approx(.005+4e-9*(1e6-n),abs=1e-9)
        assert p.objective<1e-10 and p.rmse<1.
    assert max(p.objective for p in fixed.points)>100
    assert prof.plausible_ranges==((5e5,1.5e6),)
    assert asdict(exact)==before


def test_real_surface_analytical_valley_and_invalid_nodes(exact):
    surface=objective_surface(exact,"oil_in_place","aquifer.capacity",5)
    npv=np.array([s.cumulative.np for s in exact.final_simulation.states])
    obs=np.array([s.observed_pressure for s in exact.final_simulation.states])
    for p in surface.points:
        n,c=dict(p.parameters).values()
        # Pressure follows Np*Boi=[(N-Np)*c_o+Caq]*dp.
        pressure=30e6-1.2*npv/(4e-9*(n-npv)+c)
        expected=np.sum(((pressure-obs)/1e5)**2)
        assert p.valid and p.objective==pytest.approx(expected,abs=1e-7,rel=1e-8)
    valley=[p for p in surface.points if p.objective<1e-8]
    assert len(valley)>=3  # diagonal N/C combinations, not just one minimum
    assert len({dict(p.parameters)["oil_in_place"] for p in valley})>=3
    broad=replace(exact,specifications=tuple(replace(s,lower=100.) if s.name=="oil_in_place" else replace(s,lower=1e-8) for s in exact.specifications))
    bad=objective_surface(broad,"oil_in_place","aquifer.capacity",3)
    assert bad.failed_evaluations>0
    assert all(p.objective is None and p.rmse is None and p.failure for p in bad.points if not p.valid)


def test_zero_sensitivity_no_false_tight_interval():
    fit=synthetic_match("zero")[2]
    r=analyze(fit,True)
    assert max(abs(x[1]) for x in r.sensitivity)<1e-10
    assert r.rank==1 and r.intervals[1] is None
    assert r.statuses[1][1]=="NON-IDENTIFIABLE"


def test_bound_onesided_stencil_and_no_outside_interval():
    fit=synthetic_match(upper=9e5)[2]
    r=analyze(fit,True)
    assert r.statuses[0][1]=="BOUND LIMITED"
    assert r.intervals[0] is None and "backward" in r.schemes[0][1]
    p=parameter_profile(fit,"oil_in_place",5,True,plausible_delta=1.)
    assert all(5e5<=dict(q.parameters)["oil_in_place"]<=9e5 for q in p.points)
    assert any("beyond" in w for w in p.warnings)


def test_no_statistical_assumptions_or_zero_variance_withholds(single):
    assert analyze(single).intervals==(None,)
    noiseless=synthetic_match(weighted=False)[2]
    assert analyze(noiseless,True).intervals==(None,)
    assert profile_threshold(noiseless,True)[0] is None
    assert "no confidence probability" in profile_threshold(noiseless,False,2.)[1]


def test_profile_does_not_bridge_failed_unresolved_nodes():
    def p(x,valid=True,status="CONVERGED"):
        return GridPoint((("N",x),),1. if valid else None,0.,valid,status,(),None,1)
    assert grid_ranges([p(0),p(1,False),p(2),p(3,status="NOT CONVERGED"),p(4)],"N",2.)==((0,0),(2,2),(4,4))


def test_weighting_changes_information_consistently(single):
    history=tuple(replace(s,observed_pressure=s.observed_pressure+(2e6 if i==7 else 0)) for i,s in enumerate(runpy.run_path(str(ROOT/"examples/history_match_benchmarks.py"))["benchmark_case"]()[1]))
    obs=tuple(replace(o,sigma=5e6 if i==7 else 5e4) for i,o in enumerate(observations_from_history(history)))
    a=history_match(single.base_tank,history,single.specifications,obs)
    b=history_match(single.base_tank,history,single.specifications,obs,"Pressure uncertainty weighted")
    ra,rb=analyze(a,True),analyze(b,True)
    for fit,r,scales in ((a,ra,np.full(8,1e6)),(b,rb,np.array([o.sigma for o in obs]))):
        assert np.asarray(r.sensitivity)[:,0]==pytest.approx(np.asarray(r.residual_jacobian)[:,0]*scales)
        assert r.covariance is not None and r.intervals[0] is not None
    assert not np.allclose(ra.sensitivity,rb.sensitivity)
    assert ra.standard_errors[0]>rb.standard_errors[0]*10


def test_field_si_complete_equivalence(near):
    from material_balance_studio.diagnostics.pressure_qc import diagnostic_history_from_frame
    from material_balance_studio.units.conversions import from_si
    raw=pd.DataFrame([dict(date=str(h.date),np=from_si(h.cumulative.np,"liquid_volume","FIELD"),gp=from_si(h.cumulative.gp,"gas_volume","FIELD"),wp=0,
        observed_pressure=from_si(h.observed_pressure,"pressure","FIELD"),pressure_sigma=from_si(1e5,"pressure","FIELD")) for h in near.final_simulation.states])
    history,meta=diagnostic_history_from_frame(raw,"FIELD")
    specs=tuple(replace(s,initial=from_display(to_display(s.initial,s.unit,"FIELD"),s.unit,"FIELD"),
        lower=from_display(to_display(s.lower,s.unit,"FIELD"),s.unit,"FIELD"),upper=from_display(to_display(s.upper,s.unit,"FIELD"),s.unit,"FIELD")) for s in near.specifications)
    field=history_match(near.base_tank,history,specs,observations_from_history(history,meta),near.weighting_mode)
    a,b=analyze(near,True),analyze(field,True)
    assert np.asarray(a.scaled_sensitivity)==pytest.approx(np.asarray(b.scaled_sensitivity),rel=1e-6)
    assert np.asarray(a.correlation)==pytest.approx(np.asarray(b.correlation),rel=1e-6)
    assert [s[1] for s in a.statuses]==[s[1] for s in b.statuses]
    for units in ("SI","FIELD"):
        displayed=sensitivity_frame(a,near,units)
        assert displayed.iloc[:,0].to_numpy()/to_display(1.,"pressure",units)*to_display(1.,"oil_volume",units)==pytest.approx(np.asarray(a.sensitivity)[:,0])
    pa,pb=parameter_profile(near,"oil_in_place",3,False),parameter_profile(field,"oil_in_place",3,False)
    sa,sb=objective_surface(near,"oil_in_place","aquifer.capacity",3),objective_surface(field,"oil_in_place","aquifer.capacity",3)
    assert [p.objective for p in pa.points]==pytest.approx([p.objective for p in pb.points],abs=1e-7)
    assert [p.objective for p in sa.points]==pytest.approx([p.objective for p in sb.points],abs=1e-7)
    assert point_frame(sa.points,near,"FIELD").iloc[0,0]==pytest.approx(to_display(dict(sa.points[0].parameters)["oil_in_place"],"oil_volume","FIELD"))


def test_determinism_and_no_mutation(exact):
    before=asdict(exact)
    assert analyze(exact)==analyze(exact)
    assert parameter_profile(exact,"oil_in_place",3,True)==parameter_profile(exact,"oil_in_place",3,True)
    assert objective_surface(exact,"oil_in_place","aquifer.capacity",3)==objective_surface(exact,"oil_in_place","aquifer.capacity",3)
    assert asdict(exact)==before


@pytest.mark.parametrize("count",[0,2,22,3.5])
def test_safe_grid_resolution(single,count):
    with pytest.raises(ValueError,match="resolution"):
        parameter_profile(single,"oil_in_place",count)


def test_scenario_identity_tracks_weighting_bounds(single):
    assert scenario_identity(single)!=scenario_identity(replace(single,weighting_mode="Unweighted"))
    assert scenario_identity(single)!=scenario_identity(replace(single,specifications=(replace(single.specifications[0],upper=2e6),)))
    state=single.final_simulation.states[0]
    other=replace(single,final_simulation=replace(single.final_simulation,states=(replace(state,cumulative=replace(state.cumulative,np=state.cumulative.np+1)),*single.final_simulation.states[1:])))
    assert scenario_identity(single)!=scenario_identity(other)


def test_protected_matching_source_hashes():
    for path,digest in json.loads((ROOT/"docs/phase_4c_matching_baseline.json").read_text()).items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest


def test_log_transform_uses_physical_sensitivity_and_geometric_grid(single):
    history=synthetic_match()[1]
    fit=history_match(single.base_tank,history,[replace(single.specifications[0],transformation="log")],single.observations,single.weighting_mode)
    assert np.asarray(analyze(fit).sensitivity)==pytest.approx(np.asarray(analyze(single).sensitivity),rel=1e-6)
    profile=parameter_profile(fit,"oil_in_place",3,False)
    values=[dict(p.parameters)["oil_in_place"] for p in profile.points]
    assert any(v==pytest.approx(np.sqrt(5e5*1.5e6),rel=1e-12) for v in values)
    assert all(5e5<=v<=1.5e6 for v in values)


def test_excluded_observations_absent_from_jacobian(single):
    observations=tuple(replace(o,include_in_match=i!=3) for i,o in enumerate(single.observations))
    fit=history_match(single.base_tank,synthetic_match()[1],single.specifications,observations,single.weighting_mode)
    result=analyze(fit)
    assert len(result.dates)==7
    assert str(observations[3].date) not in result.dates


def test_profile_failed_and_budget_exhausted_nodes_are_unresolved(exact,monkeypatch):
    import material_balance_studio.uncertainty.profile as module
    limited=parameter_profile(exact,"oil_in_place",3,True,max_nfev=1,plausible_delta=1.)
    assert any(p.optimizer_status=="NOT CONVERGED" for p in limited.points)
    assert all(not (a<=dict(p.parameters)["oil_in_place"]<=b) for p in limited.points if p.optimizer_status=="NOT CONVERGED" for a,b in limited.plausible_ranges)
    def failure(*args,**kwargs):
        raise ValueError("Independent injected forward/optimizer construction failure")
    monkeypatch.setattr(module,"history_match",failure)
    failed=parameter_profile(exact,"oil_in_place",3,True,plausible_delta=1.)
    assert failed.failed_evaluations==len(failed.points)
    assert not failed.plausible_ranges
    assert all(not p.valid and p.objective is None for p in failed.points)


def test_finite_difference_failure_not_misread_as_information(single,monkeypatch):
    import material_balance_studio.uncertainty.sensitivity as module
    original=module.evaluate
    calls=[]
    def unavailable(obj,values):
        calls.append(tuple(values))
        return original(obj,values) if len(calls)==1 else None
    monkeypatch.setattr(module,"evaluate",unavailable)
    with pytest.raises(ValueError,match="valid sensitivity stencil"):
        analyze(single)


def test_profile_evidence_downgrades_classification(single):
    from material_balance_studio.uncertainty.identifiability import classify
    r=analyze(single)
    p=parameter_profile(single,"oil_in_place",3,False,plausible_delta=1e9)
    status=classify(single.specifications,[dict(single.fitted_parameters)["oil_in_place"]],r.scaled_sensitivity,
                    matrix_diagnostics(r.scaled_jacobian),8,.2,(p,))
    assert status[0][1]=="POORLY CONSTRAINED"


def test_ui_existing_scenario_diagnostics_and_unit_switch(near):
    app=AppTest.from_file(str(ROOT/"app.py"),default_timeout=120).run()
    app.session_state["history_match_scenarios"]=[near]
    app.run()
    original=app.session_state["setup_si"].copy()
    assert not app.exception
    assert app.session_state["identifiability_results"]=={}
    app.button(key="unc_run").click().run()
    assert not app.exception and not app.error
    identity=scenario_identity(near)
    assert app.session_state["identifiability_results"][identity].rank==2
    app.slider(key="unc_surface_count_"+identity[:16]).set_value(3).run()
    app.button(key="unc_surface_run").click().run()
    assert not app.exception and not app.error
    app.slider(key="unc_profile_count_"+identity[:16]).set_value(3).run()
    app.button(key="unc_profile_run").click().run()
    assert not app.exception and not app.error
    saved=app.session_state["identifiability_results"][identity]
    assert len(saved.profiles)==len(saved.surfaces)==1
    app.selectbox(key="display_units").select("FIELD").run()
    assert not app.exception and not app.error
    assert app.session_state["identifiability_results"][identity]==saved
    assert app.session_state["setup_si"]==original
    assert app.session_state["history_match_scenarios"]==[near]
