"""Independent storage, transfer and conservation oracles for coupled pressures."""
from dataclasses import asdict,replace,fields
from datetime import timedelta
from pathlib import Path
from math import exp,fsum
import hashlib,json,runpy
import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest
from material_balance_studio.domain.models import HistoryRecord,CumulativeVolumes
from material_balance_studio.tank_network import NetworkTank,Connection,TankNetwork,NetworkSettings,simulate_network
from material_balance_studio.tank_network.residuals import initial_state,evaluate_candidate,accepted
from material_balance_studio.tank_network.history import volumes_at,allocate_field_history
from material_balance_studio.tank_network.diagnostics import maxima,tank_terms
from material_balance_studio.solver.simulation import simulate
from material_balance_studio.units.display import from_display,to_display
from material_balance_studio.units.conversions import from_si
from material_balance_studio.io.history import history_from_frame,pvt_from_frame
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.presentation.tank_network import tank_frame,connection_frame

ROOT=Path(__file__).resolve().parents[1]
helpers=runpy.run_path(str(ROOT/"examples/network_benchmarks.py"))
storage_tank,network_case=helpers["storage_tank"],helpers["network_case"]


@pytest.mark.parametrize("enabled,transmissibility",[(True,0.),(False,1e-9)])
def test_zero_transmissibility_independent_single_tank_equivalence(enabled,transmissibility):
    network=network_case("production",transmissibility)
    network=replace(network,connections=(replace(network.connections[0],enabled=enabled),))
    result=simulate_network(network)
    assert result.converged
    for j,node in enumerate(network.tanks):
        independent=simulate(node.reservoir,node.history)
        assert independent.converged
        for state,old in zip(result.states,independent.states):
            assert state.tanks[j].pressure==pytest.approx(old.pressure,abs=1e-4)
            assert state.tanks[j].components.aquifer_support==pytest.approx(old.balance.aquifer_support,abs=1e-6)
            assert state.tanks[j].intertank_support==0


def test_two_identical_isolated_tanks_reproduce_same_single_tank_solution():
    base=network_case("production",0.)
    a=base.tanks[0]
    network=replace(base,tanks=(a,replace(a,name="B")))
    result=simulate_network(network)
    independent=simulate(a.reservoir,a.history)
    assert result.converged
    for state,expected in zip(result.states,independent.states):
        assert state.tanks[0].pressure==state.tanks[1].pressure
        assert state.tanks[0].pressure==pytest.approx(expected.pressure,abs=1e-4)
        assert state.connections[0].cumulative_transfer==0


def test_empty_history_tank_still_communicates_on_other_tank_timeline():
    base=network_case("production")
    network=replace(base,tanks=(base.tanks[0],replace(base.tanks[1],history=())))
    result=simulate_network(network)
    reference=simulate_network(base)
    assert result.converged and result.states==reference.states
    assert result.states[-1].tanks[1].pressure<network.tanks[1].reservoir.initial_pressure


def test_closed_two_tank_independent_trapezoidal_storage_recurrence():
    network=network_case()
    result=simulate_network(network)
    assert result.converged
    ca,cb=.004,.008
    mean=(ca*30e6+cb*20e6)/(ca+cb)
    difference=10e6
    lam=1e-9*(1/ca+1/cb)
    prior=(30e6,20e6)
    for state in result.internal_states:
        assert state.diagnostics.initial_guess==prior
        difference*= (1-lam*state.dt/2)/(1+lam*state.dt/2)
        a,b=state.tanks
        assert a.pressure==pytest.approx(mean+cb/(ca+cb)*difference,abs=2e-5)
        assert b.pressure==pytest.approx(mean-ca/(ca+cb)*difference,abs=2e-5)
        assert ca*(a.pressure-30e6)+cb*(b.pressure-20e6)==pytest.approx(0.,abs=1e-6)
        assert a.pressure>b.pressure
        assert state.transfer_error==state.incremental_transfer_error==0
        assert abs(a.residual)<1e-8 and abs(b.residual)<1e-8
        prior=(a.pressure,b.pressure)
    assert result.states[-1].tanks[0].pressure<30e6 and result.states[-1].tanks[1].pressure>20e6


def test_refinement_converges_to_independent_continuous_equilibration():
    network=network_case()
    coarse=simulate_network(network)
    fine=simulate_network(network,NetworkSettings(max_step_days=.25))
    exact=(.004*30e6+.008*20e6)/.012+(.008/.012)*10e6*exp(-1e-9*(1/.004+1/.008)*100*86400)
    error_coarse=abs(coarse.states[-1].tanks[0].pressure-exact)
    error_fine=abs(fine.states[-1].tanks[0].pressure-exact)
    assert fine.converged and error_fine<10 and error_fine<error_coarse/100


@pytest.mark.parametrize("a,b,sign",[(30e6,20e6,1),(20e6,30e6,-1),(25e6,25e6,0)])
def test_connection_direction_zero_difference_and_signed_volume(a,b,sign):
    base=network_case()
    network=replace(base,tanks=tuple(replace(n,reservoir=replace(n.reservoir,initial_pressure=p)) for n,p in zip(base.tanks,(a,b))))
    settings=NetworkSettings()
    previous=initial_state(network,settings)
    candidate=evaluate_candidate(network,previous,[a,b],86400.,settings)
    edge=candidate.connections[0]
    assert np.sign(edge.rate)==sign and np.sign(edge.incremental_transfer)==sign
    assert edge.incremental_transfer==pytest.approx(1e-9*86400*(a-b))
    assert candidate.tanks[0].intertank_support==-edge.cumulative_transfer
    assert candidate.tanks[1].intertank_support==edge.cumulative_transfer
    if sign==0:
        assert simulate_network(network).states[-1].connections[0].cumulative_transfer==0


def test_production_support_and_independent_forced_storage_equations():
    network=network_case("production")
    result=simulate_network(network)
    isolated=simulate(network.tanks[0].reservoir,network.tanks[0].history)
    assert result.converged
    a,b=result.states[-1].tanks
    assert isolated.states[-1].pressure<a.pressure<30e6
    assert b.pressure<30e6 and a.intertank_support>0 and b.intertank_support<0
    # Independent two-by-two trapezoidal storage solve at each accepted dt.
    old=np.array([30e6,30e6])
    c=np.diag([.004,.008])
    lap=np.array([[1.,-1.],[-1.,1.]])*1e-9
    for state in result.internal_states:
        delta=np.array([s.increments.wp for s in state.tanks])
        expected=np.linalg.solve(c+state.dt*lap/2,(c-state.dt*lap/2)@old-delta)
        assert [s.pressure for s in state.tanks]==pytest.approx(expected,abs=2e-5)
        old=expected


def test_aquifer_and_injection_support_remain_separate():
    result=simulate_network(network_case("aquifer"))
    assert result.converged
    a,b=result.states[-1].tanks
    assert a.components.aquifer_support==0 and b.components.aquifer_support>0
    assert a.intertank_support>0 and b.intertank_support<0
    assert b.components.withdrawal.water_injection==0
    injection=simulate_network(network_case("injection"))
    assert injection.converged
    a,b=injection.states[-1].tanks
    assert a.components.withdrawal.water_injection==0 and b.components.withdrawal.water_injection>0
    assert a.intertank_support>0 and b.intertank_support<0
    assert b.pressure>30e6  # No artificial cap at initial pressure.
    assert b.components.aquifer_support==0


def test_flow_reversal_uses_one_oriented_cumulative_state():
    result=simulate_network(network_case("reversal"))
    assert result.converged
    edges=[s.connections[0] for s in result.internal_states]
    assert any(e.rate>0 for e in edges) and any(e.rate<0 for e in edges)
    assert all(len(s.connections)==1 for s in result.internal_states)
    total=0.
    for e in edges:
        total+=e.incremental_transfer
        assert e.cumulative_transfer==total
    assert total<0  # Later injection reverses cumulative transfer as well.


def test_arbitrary_sparse_n_tank_network_and_pairwise_conservation():
    template=storage_tank()
    history=(HistoryRecord(template.initial_date+timedelta(days=10)),)
    nodes=tuple(NetworkTank(n,replace(template,initial_pressure=p),history) for n,p in (("A",30e6),("B",25e6),("C",20e6),("D",28e6)))
    network=TankNetwork(nodes,(Connection("A","B",1e-9),Connection("B","C",1e-9)))
    result=simulate_network(network)
    assert result.converged
    state=result.states[-1]
    assert state.tanks[3].pressure==28e6 and state.tanks[3].intertank_support==0
    assert len(state.connections)==2
    for edge in state.connections:
        aa=next(t for t in state.tanks if t.name==edge.from_tank)
        bb=next(t for t in state.tanks if t.name==edge.to_tank)
        assert dict(aa.connection_contributions)[edge.to_tank]+dict(bb.connection_contributions)[edge.from_tank]==0
    assert abs(fsum(t.intertank_support for t in state.tanks))<1e-8


def test_large_transmissibility_strengthens_equilibration_without_oscillation():
    low=simulate_network(network_case(transmissibility=1e-10))
    high=simulate_network(network_case(transmissibility=1e-8))
    assert low.converged and high.converged
    assert abs(high.states[-1].connections[0].pressure_difference)<abs(low.states[-1].connections[0].pressure_difference)*.01
    assert all(s.connections[0].pressure_difference>=-1e-4 for s in high.internal_states)
    assert len(high.internal_states)>len(low.internal_states)


def test_extreme_stiffness_fails_explicitly_at_work_limit():
    result=simulate_network(network_case(transmissibility=1e-5),NetworkSettings(max_substeps_per_event=3))
    assert not result.converged and "stiffness" in result.failed_timestep.diagnostics.message
    assert len(result.internal_states)==3 and len(result.states)==0


def test_candidate_repeat_determinism_and_aquifer_immutability():
    network=network_case("aquifer")
    settings=NetworkSettings()
    previous=initial_state(network,settings)
    before,model_before=asdict(previous),asdict(network)
    a=evaluate_candidate(network,previous,[29e6,29.5e6],10*86400.,settings)
    evaluate_candidate(network,previous,[28e6,29.8e6],10*86400.,settings)
    b=evaluate_candidate(network,previous,[29e6,29.5e6],10*86400.,settings)
    assert a==b and asdict(previous)==before and asdict(network)==model_before
    r1,r2=simulate_network(network),simulate_network(network)
    assert r1==r2 and asdict(network)==model_before


@pytest.mark.parametrize("name",["None","Pot","Schilthuis","Fetkovich","Carter-Tracy","Van Everdingen-Hurst","Modified Van Everdingen-Hurst"])
def test_each_aquifer_independent_assignment_and_zero_t_equivalence(name):
    from material_balance_studio.aquifer.registry import MODELS
    from material_balance_studio.presentation.aquifer import FIELDS
    base=network_case("production",0.)
    cls=MODELS[name]
    aq=cls(**{f.name:FIELDS[f.name][2] for f in fields(cls)})
    # Both produce so the independently assigned B aquifer is exercised.
    b=replace(base.tanks[1],reservoir=replace(base.tanks[1].reservoir,aquifer=aq),history=base.tanks[0].history)
    network=replace(base,tanks=(base.tanks[0],b))
    result=simulate_network(network)
    independent=simulate(b.reservoir,b.history)
    assert result.converged and independent.converged
    for state,old in zip(result.states,independent.states):
        assert state.tanks[0].components.aquifer_support==0
        assert state.tanks[1].pressure==pytest.approx(old.pressure,abs=1e-3)
        assert state.tanks[1].components.aquifer_support==pytest.approx(old.balance.aquifer_support,abs=1e-5)


def test_different_pvt_assigned_per_tank():
    base=network_case("production",0.)
    b=replace(base.tanks[1],reservoir=storage_tank(n=1e6,c=8e-9),history=base.tanks[0].history)
    result=simulate_network(replace(base,tanks=(base.tanks[0],b)))
    assert result.converged
    a,b=result.states[-1].tanks
    assert a.pressure==pytest.approx(30e6-10000/.004,abs=1e-4)
    assert b.pressure==pytest.approx(30e6-10000/.008,abs=1e-4)
    assert a.pvt.bo==pytest.approx(1.21) and b.pvt.bo==pytest.approx(1.21)


def test_union_timeline_interpolation_hold_and_observations_only_at_exact_dates():
    base=network_case("production",0.)
    a=replace(base.tanks[0],history=(HistoryRecord(base.initial_date+timedelta(days=20),CumulativeVolumes(wp=2000),29e6),))
    b=replace(base.tanks[1],history=(HistoryRecord(base.initial_date+timedelta(days=10)),HistoryRecord(base.initial_date+timedelta(days=30))))
    network=replace(base,tanks=(a,b))
    result=simulate_network(network)
    assert result.converged and len(result.states)==3
    assert [s.tanks[0].cumulative.wp for s in result.states]==[1000,2000,2000]
    assert [s.tanks[0].observed_pressure for s in result.states]==[None,29e6,None]
    assert result.states[1].tanks[0].pressure!=29e6
    assert result.states[2].tanks[0].increments.wp==0


def test_individual_balance_gate_cannot_be_hidden_by_field_cancellation():
    network=network_case()
    settings=NetworkSettings()
    start=initial_state(network,settings)
    bad=replace(start,tanks=(replace(start.tanks[0],residual=1.,relative_residual=1.),replace(start.tanks[1],residual=-1.,relative_residual=1.)))
    assert bad.signed_network_residual==0 and not accepted(bad,settings)
    assert not accepted(replace(start,transfer_error=.001),settings)
    assert not accepted(replace(start,incremental_transfer_error=.001),settings)


def test_pressure_bound_failure_preserves_previous_accepted_network():
    base=network_case("production")
    a=replace(base.tanks[0],minimum_pressure=29.99e6)
    result=simulate_network(replace(base,tanks=(a,base.tanks[1])))
    assert not result.converged and not result.states and not result.internal_states
    assert result.initial_state.tanks[0].pressure==30e6
    diag=result.failed_timestep.diagnostics
    assert diag.dominant_residuals and diag.dominant_connections
    assert "failed" in diag.message.lower()


@pytest.mark.parametrize("kind",["self","duplicate","unknown","negative","duplicate_names","date","bounds","name"])
def test_network_validation(kind):
    base=network_case()
    with pytest.raises(ValueError):
        if kind=="self": Connection("A","A",1e-9)
        elif kind=="negative": Connection("A","B",-1.)
        elif kind=="duplicate": replace(base,connections=base.connections+(Connection("B","A",1e-9),))
        elif kind=="unknown": replace(base,connections=(Connection("A","Z",1e-9),))
        elif kind=="duplicate_names": replace(base,tanks=(base.tanks[0],base.tanks[0]))
        elif kind=="date": replace(base,tanks=(base.tanks[0],replace(base.tanks[1],reservoir=replace(base.tanks[1].reservoir,initial_date=base.initial_date-timedelta(days=1)))))
        elif kind=="bounds": replace(base.tanks[0],minimum_pressure=31e6)
        elif kind=="name": replace(base.tanks[0],name="bad/name")


def test_fixed_allocation_explicit_streams_conserve_totals_and_reject_bad_shares():
    initial=storage_tank().initial_date
    history=(HistoryRecord(initial+timedelta(days=10),CumulativeVolumes(np=100,gp=8000,wp=20,winj=50,ginj=600),25e6),)
    result=allocate_field_history(history,{"A":.6,"B":.4},("np","wp","winj"))
    for stream in ("np","wp","winj"):
        assert fsum(getattr(h[0].cumulative,stream) for h in result.values())==getattr(history[0].cumulative,stream)
    assert result["A"][0].cumulative.np==60 and result["B"][0].cumulative.np==40
    assert result["A"][0].cumulative.gp==0 and result["A"][0].observed_pressure is None
    with pytest.raises(ValueError,match="sum to 1"):
        allocate_field_history(history,{"A":.6,"B":.5},("np",))


def test_field_si_full_network_equivalence():
    base=network_case("aquifer")
    nodes=[]
    for node in base.tanks:
        t=node.reservoir
        props={n:from_display(to_display(getattr(t,n),q,"FIELD"),q,"FIELD") for n,q in (("initial_pressure","pressure"),("oil_in_place","oil_volume"),("cf","compressibility"),("cw","compressibility"))}
        rows=pd.DataFrame([{n:from_si(getattr(p,n),n,"FIELD") for n in ("pressure","bo","rs","bg","bw")} for p in t.pvt_model.table.rows])
        pvt=TablePVTModel(pvt_from_frame(rows,"FIELD"))
        history=history_from_frame(pd.DataFrame([dict(date=str(h.date),**{n:from_si(getattr(h.cumulative,n),"gas_volume" if n in ("gp","ginj") else "liquid_volume","FIELD") for n in ("np","gp","wp","winj","ginj")}) for h in node.history]),"FIELD")
        aq=t.aquifer
        if aq:
            aq=replace(aq,initial_water_volume=from_display(to_display(aq.initial_water_volume,"reservoir_volume","FIELD"),"reservoir_volume","FIELD"),total_compressibility=from_display(to_display(aq.total_compressibility,"compressibility","FIELD"),"compressibility","FIELD"),productivity_index=from_display(to_display(aq.productivity_index,"aquifer_productivity","FIELD"),"aquifer_productivity","FIELD"))
        nodes.append(replace(node,reservoir=replace(t,**props,pvt_model=pvt,aquifer=aq),history=history))
    edges=tuple(replace(e,transmissibility=from_display(to_display(e.transmissibility,"aquifer_productivity","FIELD"),"aquifer_productivity","FIELD")) for e in base.connections)
    a,b=simulate_network(base),simulate_network(TankNetwork(tuple(nodes),edges))
    assert a.converged and b.converged
    for sa,sb in zip(a.states,b.states):
        for ta,tb in zip(sa.tanks,sb.tanks):
            assert ta.pressure==pytest.approx(tb.pressure,abs=1e-4)
            assert ta.components.aquifer_support==pytest.approx(tb.components.aquifer_support,abs=1e-6)
            assert abs(ta.residual-tb.residual)<1e-6
        assert sa.connections[0].rate==pytest.approx(sb.connections[0].rate,abs=1e-10)
        assert sa.connections[0].cumulative_transfer==pytest.approx(sb.connections[0].cumulative_transfer,abs=1e-6)
    assert tank_frame(a,"FIELD").iloc[0,2]==pytest.approx(to_display(a.states[0].tanks[0].pressure,"pressure","FIELD"))
    assert connection_frame(a,"FIELD").iloc[0,4]==pytest.approx(to_display(base.connections[0].transmissibility,"aquifer_productivity","FIELD"))


def test_protected_physics_matching_and_uncertainty_hashes():
    for path,digest in json.loads((ROOT/"docs/phase_5a_protected_baseline.json").read_text()).items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest


def test_ui_network_run_inspector_units_and_single_tank_preservation():
    app=AppTest.from_file(str(ROOT/"app.py"),default_timeout=90).run()
    original=dict(app.session_state["setup_si"])
    app.checkbox(key="nt_configure").check().run()
    app.button(key="nt_demo").click().run()
    assert not app.exception and not app.error
    app.button(key="nt_run").click().run()
    assert not app.exception and not app.error
    result=app.session_state["nt_result"]
    assert result.converged
    app.selectbox(key="nt_inspect_tank").select("B").run()
    app.selectbox(key="display_units").select("FIELD").run()
    assert not app.exception and not app.error
    assert app.session_state["nt_result"]==result
    assert app.session_state["setup_si"]==original
    app.text_input(key="nt_new_name").set_value("C").run()
    app.button(key="nt_add").click().run()
    assert not app.exception and len(app.session_state["nt_network"].tanks)==3


@pytest.mark.parametrize("units",["SI","FIELD"])
def test_ui_applying_displayed_connections_preserves_canonical_transmissibility(units):
    app=AppTest.from_file(str(ROOT/"app.py"),default_timeout=90).run()
    app.checkbox(key="nt_configure").check().run()
    app.button(key="nt_demo").click().run()
    app.selectbox(key="display_units").select(units).run()
    app.button(key="nt_edges_apply").click().run()
    assert not app.exception and not app.error
    assert app.session_state["nt_network"].connections[0].transmissibility==pytest.approx(1e-9,rel=1e-12)
    app.button(key="nt_run").click().run()
    assert not app.exception and not app.error
    assert app.session_state["nt_result"].converged
    assert app.session_state["nt_result"].states[-1].tanks[0].pressure==pytest.approx(23586264.055936027,abs=1e-3)


@pytest.mark.parametrize("name",["Carter-Tracy","Van Everdingen-Hurst","Modified Van Everdingen-Hurst"])
def test_transient_aquifer_memory_with_active_connection(name):
    from material_balance_studio.aquifer.registry import MODELS
    from material_balance_studio.presentation.aquifer import FIELDS
    cls=MODELS[name]
    aq=cls(**{f.name:FIELDS[f.name][2] for f in fields(cls)})
    network=network_case("production")
    network=replace(network,tanks=(network.tanks[0],replace(network.tanks[1],reservoir=replace(network.tanks[1].reservoir,aquifer=aq))))
    result=simulate_network(network)
    assert result.converged
    for state in result.internal_states:
        a,b=state.tanks
        assert a.components.aquifer_support==0
        assert b.components.aquifer_support>0 and a.intertank_support>0
        assert b.aquifer_state.elapsed_time==pytest.approx(state.elapsed_seconds)
        assert b.aquifer_state.previous_reservoir_pressure==b.pressure
    assert maxima(result)["maximum_tank_relative_residual"]<=1e-8
