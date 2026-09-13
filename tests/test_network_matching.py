"""Independent allocation arithmetic and analytical network matching acceptance."""
from dataclasses import asdict,replace
from datetime import timedelta
from pathlib import Path
from math import fsum
import hashlib,json,runpy
import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest
from material_balance_studio.domain.models import CumulativeVolumes,HistoryRecord,CUMULATIVE_FIELDS
from material_balance_studio.network_history import NetworkObservation,AllocationRow,HistoryPlan,prepare_history
from material_balance_studio.network_matching import history_match_network,apply_matched_parameters,registry,MatchParameter
from material_balance_studio.network_matching.objective import NetworkObjective
from material_balance_studio.network_matching.parameters import candidate_network,validate_specs
from material_balance_studio.network_matching.diagnostics import pressure_diagnostics,data_qc,drive_terms,field_drive,closure_summary,allocation_crosscheck
from material_balance_studio.tank_network import simulate_network,NetworkSettings
from material_balance_studio.units.display import from_display,to_display
from material_balance_studio.presentation.network_matching import audit_frame,parameter_frame,parse_observations,parse_schedule,observation_input

ROOT=Path(__file__).resolve().parents[1]
helpers=runpy.run_path(str(ROOT/"examples/network_match_benchmarks.py"))
benchmark,specification,SETTINGS=helpers["benchmark"],helpers["specification"],helpers["SETTINGS"]


@pytest.fixture(scope="module")
def t_match():
    network,obs,plan=benchmark()
    return history_match_network(network,[specification(network,"connection:A:B:T",initial=3e-10,transformation="log")],obs,plan,settings=SETTINGS)


@pytest.fixture(scope="module")
def multi_match():
    network,obs,plan=benchmark()
    specs=[specification(network,name,initial=v,transformation="log" if name.startswith("connection") else "linear") for name,v in (("tank:A:N",8e5),("tank:B:N",1.2e6),("connection:A:B:T",3e-10))]
    return history_match_network(network,specs,obs,plan,settings=SETTINGS)


@pytest.fixture(scope="module")
def aquifer_match():
    network,obs,plan=benchmark(aquifer=True)
    specs=[specification(network,name,initial=v,transformation="log") for name,v in (("tank:B:aquifer:productivity_index",4e-10),("connection:A:B:T",3e-10))]
    return history_match_network(network,specs,obs,plan,settings=SETTINGS)


def allocation_example():
    network,_,_=benchmark()
    initial=network.initial_date
    field=(HistoryRecord(initial+timedelta(days=30),CumulativeVolumes(np=300.,gp=30000.,wp=600.,winj=900.,ginj=1200.)),)
    start=(.6,.2,.3,0.,.9)
    after=(.25,.8,.7,1.,.1)
    schedules=tuple(AllocationRow(day,s,(("A",a),("B",1-a))) for day,shares in ((initial,start),(initial+timedelta(days=10),after)) for s,a in zip(CUMULATIVE_FIELDS,shares))
    plan=HistoryPlan(CUMULATIVE_FIELDS,field,schedules,True)
    return network,plan,start,after


def test_effective_date_allocation_all_streams_independent_conservation():
    network,plan,start,after=allocation_example()
    observation=NetworkObservation("A",network.initial_date+timedelta(days=15),29e6)
    result=prepare_history(network,plan,[observation])
    for stream,a,b in zip(CUMULATIVE_FIELDS,start,after):
        total=getattr(plan.field_history[-1].cumulative,stream)
        expected=total/3*a+total*2/3*b
        history=result.network.tanks[0].history
        row=next(h for h in history if h.date==network.initial_date+timedelta(days=30))
        assert getattr(row.cumulative,stream)==pytest.approx(expected,abs=1e-10)
        for day in sorted({r.end for r in result.audit}):
            records=[r for r in result.audit if r.end==day and r.stream==stream]
            assert fsum(r.allocated_increment for r in records)==pytest.approx(records[0].field_increment,abs=1e-10)
            assert max(abs(r.cumulative_conservation_error) for r in records)<1e-9
    # Change at day 10 affects (10,15], not (0,10]. Fractions never ramp.
    assert next(r for r in result.audit if r.stream=="np" and r.tank=="A" and r.end==network.initial_date+timedelta(days=10)).fraction==.6
    assert next(r for r in result.audit if r.stream=="np" and r.tank=="A" and r.end==observation.date).fraction==.25
    assert sum(h.observed_pressure is not None for n in result.network.tanks for h in n.history)==1
    assert audit_frame(result.audit,"FIELD").iloc[0].Unit


@pytest.mark.parametrize("fractions",[(("A",.2),("B",.2)),(("A",.7),("B",.7)),(("A",-.1),("B",1.1)),(("A",1.2),("B",-.2)),(("A",.5),("A",.5))])
def test_invalid_allocation_fractions_rejected(fractions):
    network,_,_=benchmark()
    with pytest.raises(ValueError): AllocationRow(network.initial_date,"np",fractions)


def test_allocation_unknown_duplicate_missing_and_direct_conflict():
    network,plan,_,_=allocation_example()
    with pytest.raises(ValueError,match="conflict"):
        prepare_history(network,replace(plan,replace_direct=False))
    with pytest.raises(ValueError,match="Duplicate"):
        prepare_history(network,replace(plan,schedules=(*plan.schedules,plan.schedules[0])))
    with pytest.raises(ValueError,match="every network tank"):
        prepare_history(network,replace(plan,schedules=(replace(plan.schedules[0],fractions=(("A",.6),("Z",.4))),*plan.schedules[1:])))
    with pytest.raises(ValueError,match="initial date"):
        prepare_history(network,HistoryPlan(("np",),plan.field_history,(AllocationRow(network.initial_date+timedelta(days=1),"np",(("A",.6),("B",.4))),)))


def test_direct_stream_not_replaced_by_other_stream_allocation():
    network,plan,_,_=allocation_example()
    oil=replace(plan,allocated_streams=("np",),schedules=tuple(s for s in plan.schedules if s.stream=="np"),replace_direct=False)
    prepared=prepare_history(network,oil)
    for node,old in zip(prepared.network.tanks,network.tanks):
        for row in old.history:
            new=next(h for h in node.history if h.date==row.date)
            assert new.cumulative.wp==row.cumulative.wp and new.cumulative.winj==row.cumulative.winj


def test_observation_order_metadata_exclusion_and_timeline():
    network,obs,plan=benchmark()
    obs=(replace(obs[1],include_in_match=False,source="RFT",qc_flag="review",note="depth corrected"),obs[0],
         NetworkObservation("B",network.initial_date+timedelta(days=15),29e6,sigma=2e5),*obs[2:])
    objective=NetworkObjective(network,[specification(network,"connection:A:B:T")],obs,plan,settings=SETTINGS)
    residual=objective([.8])
    assert len(residual)==len(obs)-1
    assert [(o.date,o.tank) for o in objective.observations]==sorted((o.date,o.tank) for o in obs if o.include_in_match)
    assert network.initial_date+timedelta(days=15) in [h.date for h in objective.prepared.network.tanks[0].history]
    rows,_,_=pressure_diagnostics(objective.last_result,obs)
    excluded=next(r for r in rows if not r.included)
    assert excluded.objective_residual is None and excluded.note=="depth corrected"
    assert excluded.source=="RFT" and excluded.qc_flag=="review"


def test_known_t_recovery_and_reversal(t_match):
    r=t_match
    assert r.converged and dict(r.fitted_parameters)["connection:A:B:T"]==pytest.approx(8e-10,rel=2e-4)
    assert dict(r.final_network_metrics)["rmse"]<1.
    rates=[s.connections[0].rate for s in r.final_simulation.states]
    assert min(rates)<0<max(rates)
    assert r.maximum_valid_relative_residual<=1e-8 and r.maximum_transfer_conservation_error<=1e-8
    assert r.minimum_closure_margin>=0


def test_known_n_n_t_recovery(multi_match):
    assert multi_match.converged
    fitted=dict(multi_match.fitted_parameters)
    for name,true in (("tank:A:N",1e6),("tank:B:N",1.5e6),("connection:A:B:T",8e-10)):
        assert fitted[name]==pytest.approx(true,rel=5e-4)
    assert dict(multi_match.final_network_metrics)["rmse"]<1.
    assert set(dict(multi_match.final_tank_metrics))=={"A","B"}


def test_known_aquifer_and_t_recovery(aquifer_match):
    assert aquifer_match.converged
    fitted=dict(aquifer_match.fitted_parameters)
    assert fitted["tank:B:aquifer:productivity_index"]==pytest.approx(1e-9,rel=5e-4)
    assert fitted["connection:A:B:T"]==pytest.approx(8e-10,rel=5e-4)


@pytest.fixture(scope="module")
def m_match():
    network,obs,plan=benchmark(gas_cap=.2)
    return history_match_network(network,[specification(network,"tank:A:m",initial=.05)],obs,plan,settings=SETTINGS)


def test_gas_cap_parameter_recovery(m_match):
    assert m_match.converged and dict(m_match.fitted_parameters)["tank:A:m"]==pytest.approx(.2,rel=5e-4)


def test_wrong_allocation_changes_fit_without_fitting_fractions():
    network,obs,plan=benchmark(allocation=True)
    specs=[specification(network,"connection:A:B:T",initial=3e-10,transformation="log")]
    before=asdict(plan)
    correct=history_match_network(network,specs,obs,plan,settings=SETTINGS)
    wrong=replace(plan,schedules=(AllocationRow(network.initial_date,"wp",(("A",.5),("B",.5))),*[s for s in plan.schedules if s.stream!="wp"]))
    incorrect=history_match_network(network,specs,obs,wrong,settings=SETTINGS)
    assert correct.converged and incorrect.converged
    assert dict(correct.final_network_metrics)["rmse"]<1
    assert dict(incorrect.final_network_metrics)["rmse"]>1e4
    assert asdict(plan)==before
    assert all("allocation" not in n for n,_ in correct.fitted_parameters)


def test_sufficiency_global_and_unobserved_tank_caution():
    network,obs,_=benchmark()
    spec=specification(network,"tank:B:N")
    one_tank=tuple(o for o in obs if o.tank=="A")
    warnings=data_qc(network,[spec],one_tank)
    assert any("STRONG CAUTION" in w and "B" in w for w in warnings)
    with pytest.raises(ValueError,match="outnumber"):
        data_qc(network,[spec],one_tank[:1])


def test_weighted_objective_and_low_quality_survey():
    network,obs,plan=benchmark()
    noisy=tuple(replace(o,pressure=o.pressure+(2e6 if i==len(obs)-1 else 0.),sigma=5e6 if i==len(obs)-1 else 2e4) for i,o in enumerate(obs))
    specs=[specification(network,"connection:A:B:T",initial=3e-10,transformation="log")]
    unweighted=history_match_network(network,specs,noisy,plan,settings=SETTINGS)
    weighted=history_match_network(network,specs,noisy,plan,"Pressure uncertainty weighted",settings=SETTINGS)
    assert unweighted.converged and weighted.converged
    assert abs(dict(weighted.fitted_parameters)[specs[0].name]-8e-10)<abs(dict(unweighted.fitted_parameters)[specs[0].name]-8e-10)*.1
    assert weighted.final_objective==pytest.approx(sum(r.objective_residual**2 for r in weighted.final_pressure_rows if r.included))
    assert all(r.objective_residual==pytest.approx(r.residual/r.sigma) for r in weighted.final_pressure_rows)


def test_sigma_policy_missing_default_and_exclusion():
    network,obs,_=benchmark()
    spec=specification(network,"connection:A:B:T")
    mixed=(replace(obs[0],sigma=None),*obs[1:])
    with pytest.raises(ValueError,match="sigma"):
        NetworkObjective(network,[spec],mixed,mode="Pressure uncertainty weighted")
    objective=NetworkObjective(network,[spec],mixed,mode="Pressure uncertainty weighted",default_sigma=2e5)
    assert objective.scales[0]==2e5
    excluded=NetworkObjective(network,[spec],(replace(mixed[0],include_in_match=False),*mixed[1:]),mode="Pressure uncertainty weighted")
    assert len(excluded.scales)==len(mixed)-1


def test_bounds_and_exact_zero_t_boundary():
    network,obs,plan=benchmark()
    spec=specification(network,"connection:A:B:T",upper=5e-10,initial=3e-10)
    result=history_match_network(network,[spec],obs,plan,settings=SETTINGS)
    assert result.converged and result.bound_flags==((spec.name,False,True),)
    assert all(spec.lower<=e.parameters[0]<=spec.upper for e in result.evaluations)
    network,obs,plan=benchmark(transmissibility=0.)
    spec=specification(network,"connection:A:B:T",lower=0.,initial=3e-10)
    result=history_match_network(network,[spec],obs,plan,settings=SETTINGS)
    assert result.converged and dict(result.fitted_parameters)[spec.name]<2e-13
    assert result.bound_flags[0][1]
    with pytest.raises(ValueError,match="strictly positive"):
        replace(spec,transformation="log")


def test_registry_bounds_and_allocation_not_eligible():
    network,obs,_=benchmark(aquifer=True,gas_cap=.2)
    names=("tank:A:N","tank:A:m","tank:B:aquifer:productivity_index","connection:A:B:T")
    specs=[specification(network,n) for n in names]
    candidate=candidate_network(network,specs,[s.lower for s in specs])
    for s in specs:
        assert registry(candidate)[s.name].value==s.lower
    with pytest.raises(ValueError,match="bounds"):
        candidate_network(network,specs,[specs[0].lower-1,*[s.initial for s in specs[1:]]])
    with pytest.raises(ValueError,match="Unapproved"):
        validate_specs(network,[replace(specs[0],name="allocation:A:wp")])


def test_objective_determinism_no_candidate_state_leak():
    network,obs,plan=benchmark(aquifer=True)
    before=asdict(network)
    spec=specification(network,"connection:A:B:T")
    objective=NetworkObjective(network,[spec],obs,plan,settings=SETTINGS)
    first=objective([.8]).copy()
    saved=asdict(objective.last_result)
    objective([.4])
    again=objective([.8])
    assert np.array_equal(first,again) and asdict(objective.last_result)==saved
    assert asdict(network)==before


def test_failed_pressure_bound_candidate_fixed_penalty():
    network,obs,plan=benchmark()
    network=replace(network,tanks=tuple(replace(t,minimum_pressure=29.99e6) for t in network.tanks))
    spec=specification(network,"connection:A:B:T")
    obj=NetworkObjective(network,[spec],obs,plan,settings=SETTINGS)
    before=asdict(network)
    a,b=obj([.8]),obj([.8])
    assert np.array_equal(a,b) and np.all(a==obj.penalty)
    assert obj.last_result is None and obj.records[-1].failure_type=="pressure-bound failure"
    assert obj.records[-1].affected and asdict(network)==before


@pytest.mark.parametrize("message",["Aquifer candidate failure","Coupled solver failed","Stiffness work-limit failure"])
def test_exception_candidates_never_crash_or_commit(message,monkeypatch):
    import material_balance_studio.network_matching.objective as module
    network,obs,plan=benchmark()
    def fail(*a,**k): raise ValueError(message)
    monkeypatch.setattr(module,"simulate_network",fail)
    obj=NetworkObjective(network,[specification(network,"connection:A:B:T")],obs,plan)
    residual=obj([.8])
    assert np.all(np.isfinite(residual)) and len(residual)==len(obs)
    assert obj.records[-1].failure==message and not obj.records[-1].valid and obj.last_result is None


def test_closure_failure_and_margin_are_not_pressure_fit(monkeypatch):
    import material_balance_studio.network_matching.objective as module
    network,obs,plan=benchmark()
    prepared=prepare_history(network,plan,obs)
    real=simulate_network(prepared.network,SETTINGS)
    first=real.internal_states[0]
    bad=replace(real,internal_states=(replace(first,tanks=(replace(first.tanks[0],residual=1.,relative_residual=.1),first.tanks[1])),*real.internal_states[1:]))
    monkeypatch.setattr(module,"simulate_network",lambda *a,**k:bad)
    obj=NetworkObjective(network,[specification(network,"connection:A:B:T")],obs,plan)
    assert np.all(obj([.8])==obj.penalty)
    assert obj.records[-1].failure_type=="closure failure"
    near=replace(real,internal_states=(replace(first,tanks=(replace(first.tanks[0],relative_residual=9e-9),first.tanks[1])),))
    assert closure_summary(near)["status"]=="NEAR TOLERANCE"
    assert closure_summary(near)["closure_margin"]==pytest.approx(1e-9)
    with pytest.raises(ValueError,match="loosen"):
        NetworkObjective(network,[specification(network,"connection:A:B:T")],obs,settings=NetworkSettings(relative_tolerance=1e-7))


def test_all_failed_match_never_converged(monkeypatch):
    import material_balance_studio.network_matching.objective as module
    network,obs,plan=benchmark()
    def fail(*a,**k): raise ValueError("PVT failure")
    monkeypatch.setattr(module,"simulate_network",fail)
    result=history_match_network(network,[specification(network,"connection:A:B:T")],obs,plan,max_nfev=2)
    assert not result.converged and result.qc_status=="FAIL" and result.fitted_network is None
    assert result.maximum_valid_relative_residual is None


def test_signed_drives_close_and_field_transfer_not_created(t_match):
    positive=negative=False
    for state in t_match.final_simulation.states:
        for tank in state.tanks:
            d=drive_terms(tank)
            if d["total"] is not None:
                assert d["values"]["TDI"]==tank.intertank_support/tank.components.withdrawal.production
                assert d["total"]==pytest.approx(1.,abs=1e-8)
                positive|=d["values"]["TDI"]>0
                negative|=d["values"]["TDI"]<0
        assert field_drive(state)["supports"]["TDI"]==pytest.approx(0.,abs=1e-8)
    assert positive and negative
    initial=t_match.final_simulation.initial_state.tanks[0]
    assert drive_terms(initial)["values"]["TDI"] is None


def test_scenario_apply_only_parameters_no_history_or_allocation_overwrite(t_match):
    current=t_match.base_network
    before=asdict(current)
    changed=apply_matched_parameters(current,t_match)
    assert asdict(current)==before
    assert tuple(t.history for t in changed.tanks)==tuple(t.history for t in current.tanks)
    assert changed.connections[0].transmissibility==dict(t_match.fitted_parameters)["connection:A:B:T"]
    with pytest.raises(ValueError,match="topology"):
        apply_matched_parameters(replace(current,connections=()),t_match)


def test_metadata_and_schedule_excel_dates_parse():
    frame=pd.DataFrame([dict(tank="A",date=pd.Timestamp("2020-01-02"),pressure=30e6,sigma=1e5,source="PBU",include_in_match=False,note="retained")])
    result=parse_observations(frame,"SI")
    assert result[0].note=="retained" and not result[0].include_in_match
    schedule=parse_schedule(pd.DataFrame([dict(date=pd.Timestamp("2020-01-01"),stream="np",tank=n,fraction=.5) for n in ("A","B")]))
    assert schedule[0].fractions==(("A",.5),("B",.5))


def test_allocation_pressure_crosscheck_is_diagnostic_only():
    from material_balance_studio.network_matching.diagnostics import PressureRow
    network,_,plan=benchmark(allocation=True)
    day=network.initial_date+timedelta(days=40)
    rows=(PressureRow("A",day-timedelta(days=10),30e6,30e6,0.,None,0.,True,"","",""),PressureRow("A",day+timedelta(days=1),27e6,30e6,3e6,None,3.,True,"","",""))
    before=asdict(plan)
    assert allocation_crosscheck(plan,rows,network)
    assert asdict(plan)==before


@pytest.mark.parametrize("fixture_name",["multi_match","aquifer_match","m_match"])
def test_field_si_network_matching_equivalence(fixture_name,request):
    from material_balance_studio.io.history import history_from_frame,pvt_from_frame
    from material_balance_studio.pvt.table_model import TablePVTModel
    from material_balance_studio.units.conversions import from_si
    from material_balance_studio.presentation.aquifer import FIELDS
    original=request.getfixturevalue(fixture_name)
    nodes=[]
    for node in original.base_network.tanks:
        tank=node.reservoir
        pvt=pvt_from_frame(pd.DataFrame([{n:from_si(getattr(p,n),n,"FIELD") for n in ("pressure","bo","rs","bg","bw")} for p in tank.pvt_model.table.rows]),"FIELD")
        values={n:from_display(to_display(getattr(tank,n),q,"FIELD"),q,"FIELD") for n,q in (("oil_in_place","oil_volume"),("initial_pressure","pressure"),("m","dimensionless"))}
        history=history_from_frame(pd.DataFrame([dict(date=str(h.date),**{n:from_si(getattr(h.cumulative,n),"gas_volume" if n in ("gp","ginj") else "liquid_volume","FIELD") for n in CUMULATIVE_FIELDS}) for h in node.history]),"FIELD")
        aquifer=tank.aquifer
        if aquifer is not None:
            aquifer=replace(aquifer,**{k:from_display(to_display(v,FIELDS[k][1],"FIELD"),FIELDS[k][1],"FIELD") for k,v in asdict(aquifer).items()})
        nodes.append(replace(node,reservoir=replace(tank,**values,pvt_model=TablePVTModel(pvt),aquifer=aquifer),history=history))
    network=replace(original.base_network,tanks=tuple(nodes),connections=tuple(replace(e,transmissibility=from_display(to_display(e.transmissibility,"aquifer_productivity","FIELD"),"aquifer_productivity","FIELD")) for e in original.base_network.connections))
    specs=tuple(replace(s,**{k:from_display(to_display(getattr(s,k),s.unit,"FIELD"),s.unit,"FIELD") for k in ("initial","lower","upper","scale")}) for s in original.specifications)
    observations=parse_observations(observation_input(original.observations,"FIELD"),"FIELD")
    result=history_match_network(network,specs,observations,original.history_plan,original.mode,original.default_sigma,original.settings)
    assert result.converged
    assert [v for _,v in result.fitted_parameters]==pytest.approx([v for _,v in original.fitted_parameters],rel=2e-5)
    for a,b in zip(result.final_simulation.states,original.final_simulation.states):
        assert [t.pressure for t in a.tanks]==pytest.approx([t.pressure for t in b.tanks],abs=.1)
        assert a.connections[0].cumulative_transfer==pytest.approx(b.connections[0].cumulative_transfer,abs=.01)
    assert result.final_objective==pytest.approx(original.final_objective,abs=1e-8)
    frame=parameter_frame(result,"FIELD")
    assert frame[frame.Parameter==specs[0].name].iloc[0].Fitted==pytest.approx(to_display(dict(result.fitted_parameters)[specs[0].name],specs[0].unit,"FIELD"))


def test_protected_single_and_phase5a_hashes():
    for name,digest in json.loads((ROOT/"docs/phase_5b_protected_baseline.json").read_text()).items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest


def test_ui_network_match_scenario_and_explicit_application(t_match):
    app=AppTest.from_file(str(ROOT/"app.py"),default_timeout=120).run()
    app.session_state["nt_network"]=t_match.base_network
    app.session_state["nh_observations"]=t_match.observations
    app.session_state["network_match_scenarios"]=[t_match]
    original=asdict(t_match.base_network)
    app.checkbox(key="nh_enable").check().run()
    assert not app.exception and not app.error
    assert asdict(app.session_state["nt_network"])==original
    app.selectbox(key="display_units").select("FIELD").run()
    assert not app.exception and not app.error
    app.button(key="nh_apply").click().run()
    assert not app.exception and not app.error
    assert app.session_state["nt_network"].connections[0].transmissibility==dict(t_match.fitted_parameters)["connection:A:B:T"]
    assert app.session_state["network_match_scenarios"]==[t_match]


def test_ui_runs_matching_using_saved_physical_bounds(t_match):
    app=AppTest.from_file(str(ROOT/"app.py"),default_timeout=120).run()
    app.session_state["nt_network"]=t_match.base_network
    app.session_state["nh_observations"]=t_match.observations
    app.session_state["nh_specs"]={s.name:replace(s,initial=dict(t_match.fitted_parameters)[s.name]) for s in t_match.specifications}
    app.checkbox(key="nh_enable").check().run()
    app.multiselect(key="nh_active").set_value(["connection:A:B:T"]).run()
    app.number_input(key="nh_step").set_value(2.).run()
    assert not app.button(key="nh_run").disabled
    app.button(key="nh_run").click().run()
    assert not app.exception and not app.error
    assert len(app.session_state["network_match_scenarios"])==1
    assert app.session_state["network_match_scenarios"][0].converged


def test_ui_explicit_default_sigma_survives_unit_switch():
    network,observations,_=benchmark()
    app=AppTest.from_file(str(ROOT/"app.py"),default_timeout=120).run()
    app.session_state["nt_network"]=network
    app.session_state["nh_observations"]=observations
    app.checkbox(key="nh_enable").check().run()
    app.checkbox(key="nh_default_enabled").check().run()
    key=next(w.key for w in app.number_input if w.key and w.key.startswith("nh_") and w.key.endswith("sigma"))
    app.number_input(key=key).set_value(.25).run()
    app.selectbox(key="display_units").select("FIELD").run()
    assert not app.exception and not app.error
    key=next(w.key for w in app.number_input if w.key and w.key.startswith("nh_") and w.key.endswith("sigma"))
    assert app.number_input(key=key).value==pytest.approx(to_display(2.5e5,"pressure","FIELD"))
    assert app.session_state["nh_default_sigma_si"]==pytest.approx(2.5e5)


def test_ui_allocation_save_preserves_prepared_histories_and_base_network():
    network,observations,plan=benchmark(allocation=True)
    before=asdict(network)
    expected=prepare_history(network,plan,observations)
    app=AppTest.from_file(str(ROOT/"app.py"),default_timeout=120).run()
    app.session_state["nt_network"]=network
    app.session_state["nh_observations"]=observations
    app.session_state["nh_plan"]=plan
    app.checkbox(key="nh_enable").check().run()
    assert not app.exception and not app.error
    app.button(key="nh_save_inputs").click().run()
    assert not app.exception and not app.error
    saved=prepare_history(network,app.session_state["nh_plan"],app.session_state["nh_observations"])
    assert saved.network==expected.network and saved.audit==expected.audit
    assert asdict(app.session_state["nt_network"])==before
