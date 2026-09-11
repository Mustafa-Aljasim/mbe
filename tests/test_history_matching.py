"""Independent bounded-match recovery, failure, purity and UI acceptance."""
from dataclasses import replace,asdict
from pathlib import Path
import runpy,json
import numpy as np
import pytest
from scipy.optimize import minimize_scalar
from streamlit.testing.v1 import AppTest
from material_balance_studio.matching import history_match,observations_from_history,MatchParameter
from material_balance_studio.matching.parameters import parameter_registry,parameter_value
from material_balance_studio.matching.objective import PressureObjective
from material_balance_studio.presentation.history_matching import fitted_parameter_frame
from material_balance_studio.units.display import to_display,from_display

ROOT=Path(__file__).resolve().parents[1]
helpers=runpy.run_path(str(ROOT/"examples/history_match_benchmarks.py"))
benchmark_case,specification=helpers["benchmark_case"],helpers["specification"]


@pytest.mark.parametrize("transform",["linear","log"])
def test_independent_ooip_recovery(transform):
    true,history=benchmark_case()
    base=replace(true,oil_in_place=7e5)
    snapshot=asdict(base)
    result=history_match(base,history,[specification(base,"oil_in_place",5e5,1.5e6,transform)],observations_from_history(history))
    assert result.converged
    assert dict(result.fitted_parameters)["oil_in_place"]==pytest.approx(1e6,rel=2e-6)
    assert result.final_metrics["rmse"]<1
    assert result.maximum_mbe_relative<1e-8
    assert asdict(base)==snapshot and base.oil_in_place==7e5
    assert result.ho_ooip==pytest.approx(1e6)


@pytest.mark.parametrize("multi",[False,True])
def test_independent_fetkovich_j_and_n_recovery(multi):
    true,history=benchmark_case(aquifer=True)
    base=replace(true,oil_in_place=8e5 if multi else 1e6,aquifer=replace(true.aquifer,productivity_index=4e-10))
    specs=[specification(base,"aquifer.productivity_index",1e-10,3e-9,"log")]
    if multi:
        specs.insert(0,specification(base,"oil_in_place",5e5,1.5e6))
    result=history_match(base,history,specs,observations_from_history(history))
    assert result.converged
    assert result.fitted_tank.aquifer.productivity_index==pytest.approx(1e-9,rel=2e-4)
    assert result.fitted_tank.oil_in_place==pytest.approx(1e6,rel=2e-4)
    assert result.final_metrics["rmse"]<1


def test_m_recovery_and_inactive_parameter_fixed():
    true,history=benchmark_case(gas_cap=.4)
    base=replace(true,m=.1)
    specs=[specification(base,"m",0,1),replace(specification(base,"oil_in_place",5e5,2e6),active=False,initial=8e5)]
    result=history_match(base,history,specs,observations_from_history(history))
    assert result.converged and result.fitted_tank.m==pytest.approx(.4,rel=2e-5)
    assert result.fitted_tank.oil_in_place==1e6
    assert list(dict(result.fitted_parameters))==["m"]


def test_weighted_noisy_observation_matches_independent_analytical_oracle():
    true,history=benchmark_case()
    noisy=tuple(replace(h,observed_pressure=h.observed_pressure+(2e6 if i==7 else 0)) for i,h in enumerate(history))
    obs=tuple(replace(o,sigma=5e6 if i==7 else 5e4) for i,o in enumerate(observations_from_history(noisy)))
    base=replace(true,oil_in_place=8e5)
    spec=[specification(base,"oil_in_place",5e5,1.5e6)]
    unweighted=history_match(base,noisy,spec,obs)
    weighted=history_match(base,noisy,spec,obs,"Pressure uncertainty weighted")
    def objective(n,weighted_mode):
        # Independent inversion of Np*(Boi+c*dp)=N*c*dp, no simulator.
        return sum(((30e6-1.2*h.cumulative.np/(4e-9*(n-h.cumulative.np))-o.pressure)/(o.sigma if weighted_mode else 1e6))**2 for h,o in zip(noisy,obs))
    for match,mode in ((unweighted,False),(weighted,True)):
        oracle=minimize_scalar(lambda n:objective(n,mode),bounds=(5e5,1.5e6),method="bounded",options={"xatol":.001})
        assert match.fitted_tank.oil_in_place==pytest.approx(oracle.x,rel=1e-5)
    assert abs(weighted.fitted_tank.oil_in_place-1e6)<abs(unweighted.fitted_tank.oil_in_place-1e6)*.01


def test_restrictive_bound_flag_and_every_trial_in_bounds():
    true,history=benchmark_case()
    base=replace(true,oil_in_place=7e5)
    result=history_match(base,history,[specification(base,"oil_in_place",5e5,9e5)],observations_from_history(history))
    assert result.converged
    assert result.bound_flags==(("oil_in_place",False,True),)
    assert all(5e5<=e.parameters[0]<=9e5 for e in result.evaluations)
    assert result.specifications[0].upper==9e5
    assert result.qc_status=="CAUTION"


def test_failed_candidate_fixed_penalty_and_deterministic_recovery():
    true,history=benchmark_case(aquifer=True)
    base=replace(true,oil_in_place=8e5)
    # Match N with no aquifer so low N cannot sustain withdrawal inside PVT bounds.
    base=replace(base,aquifer=None)
    spec=specification(base,"oil_in_place",100,2e6)
    objective=PressureObjective(base,history,[spec],observations_from_history(history))
    failed=objective([spec.encode(100)])
    assert len(failed)==len(history) and np.all(np.isfinite(failed))
    assert objective.records[-1].failure
    first=objective([spec.encode(1e6)])
    objective([spec.encode(100)])
    assert np.array_equal(first,objective([spec.encode(1e6)]))
    assert objective.last_result is not None


def test_weighting_missing_sigma_requires_explicit_policy_and_selection():
    tank,history=benchmark_case()
    specs=[specification(tank,"oil_in_place",5e5,1.5e6)]
    obs=observations_from_history(history)
    with pytest.raises(ValueError,match="sigma"):
        PressureObjective(tank,history,specs,obs,"Pressure uncertainty weighted")
    obj=PressureObjective(tank,history,specs,obs,"Pressure uncertainty weighted",1e5)
    assert len(obj([1]))==8
    partial=tuple(replace(o,sigma=1e5 if i else None,include_in_match=i!=0) for i,o in enumerate(obs))
    obj=PressureObjective(tank,history,specs,partial,"Pressure uncertainty weighted")
    assert len(obj([1]))==7


def test_insufficient_data_alignment_and_fixed_registry():
    tank,history=benchmark_case(aquifer=True)
    obs=observations_from_history(history)
    spec=specification(tank,"oil_in_place",5e5,2e6)
    with pytest.raises(ValueError,match="outnumber"):
        PressureObjective(tank,history,[spec],[replace(o,include_in_match=i==0) for i,o in enumerate(obs)])
    with pytest.raises(ValueError,match="exact"):
        PressureObjective(tank,history,[spec],[replace(obs[0],pressure=12e6),*obs[1:]])
    for fixed in ("initial_pressure","cf","cw","aquifer.total_compressibility","aquifer.radius_ratio"):
        assert fixed not in parameter_registry(tank)
    with pytest.raises(ValueError,match="approved"):
        PressureObjective(tank,history,[replace(spec,name="cf")],obs)


def test_unit_invariance_repeatability_and_result_display():
    true,history=benchmark_case()
    base=replace(true,oil_in_place=8e5)
    si=specification(base,"oil_in_place",5e5,1.5e6)
    field=replace(si,initial=from_display(to_display(si.initial,"oil_volume","FIELD"),"oil_volume","FIELD"),lower=from_display(to_display(si.lower,"oil_volume","FIELD"),"oil_volume","FIELD"),upper=from_display(to_display(si.upper,"oil_volume","FIELD"),"oil_volume","FIELD"))
    obs=observations_from_history(history)
    a=history_match(base,history,[si],obs)
    b=history_match(base,history,[field],obs)
    c=history_match(base,history,[si],obs)
    assert a.fitted_tank.oil_in_place==pytest.approx(b.fitted_tank.oil_in_place,rel=1e-10)
    assert a.fitted_parameters==c.fitted_parameters
    assert a.pvt_fingerprint==c.pvt_fingerprint
    assert fitted_parameter_frame(a,"FIELD").iloc[0].Fitted==pytest.approx(to_display(a.fitted_tank.oil_in_place,"oil_volume","FIELD"))


def test_ui_scenario_no_automatic_apply_and_explicit_apply():
    app=AppTest.from_file(str(ROOT/"app.py"),default_timeout=60).run()
    app.button[0].click().run()
    original=app.session_state["setup_si"]["setup_oil"]
    app.checkbox(key="hm_configure").check().run()
    app.multiselect(key="hm_active").set_value(["oil_in_place"]).run()
    lower_key=next(w.key for w in app.number_input if w.key and w.key.startswith("hm_lower_"))
    app.number_input(key=lower_key).set_value(300000.).run()
    app.button(key="run_history_match").click().run()
    assert not app.exception and not app.error
    results=app.session_state["history_match_scenarios"]
    assert len(results)==1 and results[0].converged
    assert app.session_state["setup_si"]["setup_oil"]==original
    app.selectbox(key="display_units").select("FIELD").run()
    assert not app.exception
    assert app.number_input(key=lower_key).value==pytest.approx(to_display(300000.,"oil_volume","FIELD"))
    assert app.session_state["history_match_scenarios"][0].fitted_parameters==results[0].fitted_parameters
    app.button(key="apply_history_match").click().run()
    assert not app.exception
    assert app.session_state["setup_si"]["setup_oil"]==results[0].fitted_tank.oil_in_place


@pytest.mark.parametrize("name",["None","Pot","Schilthuis","Fetkovich","Carter-Tracy","Van Everdingen-Hurst","Modified Van Everdingen-Hurst"])
def test_all_aquifers_objective_interface_and_applicable_registry(name):
    from dataclasses import fields
    from material_balance_studio.aquifer.registry import MODELS
    from material_balance_studio.presentation.aquifer import FIELDS
    tank,history=benchmark_case()
    cls=MODELS[name]
    aq=cls(**{f.name:FIELDS[f.name][2] for f in fields(cls)})
    tank=replace(tank,aquifer=aq)
    registry=parameter_registry(tank)
    assert "aquifer.radius_ratio" not in registry
    spec=specification(tank,"oil_in_place",5e5,2e6)
    objective=PressureObjective(tank,history,[spec],observations_from_history(history))
    result=objective([spec.encode(1e6)])
    assert objective.last_result is not None and len(result)==8
    if name=="Fetkovich":
        assert "aquifer.productivity_index" in registry and "aquifer.initial_water_volume" in registry
    if name in ("Carter-Tracy","Van Everdingen-Hurst","Modified Van Everdingen-Hurst"):
        assert "aquifer.permeability" in registry and "aquifer.encroachment_angle" in registry


def test_bad_mbe_closure_is_not_accepted_even_with_exact_pressures(monkeypatch):
    from material_balance_studio.solver.simulation import simulate
    import material_balance_studio.matching.objective as module
    tank,history=benchmark_case()
    result=simulate(tank,history)
    bad=replace(result,states=(replace(result.states[0],balance=replace(result.states[0].balance,absolute_residual=1.,relative_residual=.1)),*result.states[1:]))
    monkeypatch.setattr(module,"simulate",lambda *args:bad)
    objective=PressureObjective(tank,history,[specification(tank,"oil_in_place",5e5,2e6)],observations_from_history(history))
    residual=objective([1.])
    assert objective.last_result is None
    assert "MBE closure" in objective.records[-1].failure
    assert np.all(residual==objective.penalty)


def test_all_failed_match_cannot_be_converged():
    tank,history=benchmark_case()
    tank=replace(tank,oil_in_place=200.)
    result=history_match(tank,history,[specification(tank,"oil_in_place",100,300)],observations_from_history(history))
    assert not result.converged and result.qc_status=="FAIL"
    assert result.fitted_tank is None and result.final_metrics["rmse"] is None


def test_field_history_and_uncertainty_input_fit_equivalence():
    import pandas as pd
    from material_balance_studio.units.conversions import from_si
    from material_balance_studio.diagnostics.pressure_qc import diagnostic_history_from_frame
    tank,history=benchmark_case()
    raw=pd.DataFrame([dict(date=str(h.date),np=from_si(h.cumulative.np,"liquid_volume","FIELD"),gp=from_si(h.cumulative.gp,"gas_volume","FIELD"),wp=0,observed_pressure=from_si(h.observed_pressure,"pressure","FIELD"),pressure_sigma=from_si(1e5,"pressure","FIELD")) for h in history])
    field,meta=diagnostic_history_from_frame(raw,"FIELD")
    tank=replace(tank,oil_in_place=8e5)
    spec=[specification(tank,"oil_in_place",5e5,1.5e6)]
    a=history_match(tank,field,spec,observations_from_history(field,meta),"Pressure uncertainty weighted")
    b=history_match(tank,history,spec,tuple(replace(o,sigma=1e5) for o in observations_from_history(history)),"Pressure uncertainty weighted")
    assert a.fitted_tank.oil_in_place==pytest.approx(b.fitted_tank.oil_in_place,rel=1e-9)
