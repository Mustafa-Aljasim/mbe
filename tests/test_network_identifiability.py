"""Independent Phase 5C information, allocation and immutable-workflow acceptance."""
from dataclasses import asdict,replace
from pathlib import Path
from datetime import timedelta
import hashlib,json,runpy
import numpy as np
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest
from material_balance_studio.uncertainty.network_sensitivity import analyze_network,coupling_pairs
from material_balance_studio.uncertainty.network_profiles import network_surface,network_profile
from material_balance_studio.uncertainty.allocation_sensitivity import shift_allocation,allocation_study
from material_balance_studio.uncertainty.network_context import context,fixed_point
from material_balance_studio.uncertainty.context import evaluate
from material_balance_studio.uncertainty.jacobian import matrix_diagnostics
from material_balance_studio.network_history import AllocationRow,HistoryPlan,prepare_history
from material_balance_studio.network_matching import history_match_network
from material_balance_studio.presentation.network_matching import observation_input,parse_observations
from material_balance_studio.presentation.uncertainty import sensitivity_frame
from material_balance_studio.units.display import from_display,to_display

ROOT=Path(__file__).resolve().parents[1]
B=runpy.run_path(str(ROOT/'examples/network_information_benchmarks.py'))


@pytest.fixture(scope='module')
def well(): return B['truth_match']('well')

@pytest.fixture(scope='module')
def information(well): return analyze_network(well,True)

@pytest.fixture(scope='module')
def nt(): return B['truth_match']('NT_confounding')

@pytest.fixture(scope='module')
def aq(): return B['truth_match']('aquifer_confounding')

@pytest.fixture(scope='module')
def allocation(): return B['truth_match']('allocation')

@pytest.fixture(scope='module')
def allocation_results(allocation):
    cases=tuple((f'A={f:.2f}',shift_allocation(allocation.history_plan,'np','B','A',f-.6)) for f in (.5,.55,.6,.65,.7))
    return allocation_study(allocation,cases)


def test_well_identified_network(well,information):
    a=information.local
    assert well.converged and well.final_objective<1e-10
    assert a.rank==3 and a.condition<20
    assert max(abs(a.correlation[i][j]) for i in range(3) for j in range(i))<.7
    assert all(s=='WELL CONSTRAINED' for _,s,_ in a.statuses)
    for (name,value),interval in zip(well.fitted_parameters,a.intervals):
        assert interval is not None and interval[0]<value<interval[1]
    assert len(information.observation_tanks)==16 and len(information.by_tank)==2
    assert len(coupling_pairs(well))==3


def test_jacobian_independent_storage_oracle(well,information):
    # Independently perturb the known linear two-tank recurrence in the benchmark.
    # T changes its graph Laplacian; this is independent of NetworkObjective differentiation.
    helper=B['base']
    lower=helper['benchmark'](transmissibility=8e-10-1e-13)[1]
    upper=helper['benchmark'](transmissibility=8e-10+1e-13)[1]
    truth=np.array([(b.pressure-a.pressure)/2e-13 for a,b in zip(lower,upper)])
    actual=np.array(information.local.sensitivity)[:,2]
    assert actual==pytest.approx(truth,rel=2e-5,abs=5e7)
    sigma=np.array([o.sigma for o in well.observations])
    assert np.array(information.local.residual_jacobian)[:,2]==pytest.approx(truth/sigma,rel=2e-5,abs=1000)


def test_unweighted_objective_scaling(well,information):
    f=history_match_network(well.base_network,well.specifications,well.observations,well.history_plan,
        settings=well.settings)
    a=analyze_network(f).local
    assert np.array(a.sensitivity)==pytest.approx(np.array(information.local.sensitivity),rel=2e-4,abs=1e5)
    assert np.array(a.residual_jacobian)==pytest.approx(np.array(a.sensitivity)/1e6)
    assert a.intervals==(None,None,None)


def test_nt_confounding(nt):
    a=analyze_network(nt,True)
    assert nt.final_objective<1e-10
    assert a.local.correlation[0][1]>.98
    assert all(s=='POORLY CONSTRAINED' for _,s,_ in a.local.statuses)
    assert a.local.intervals[1] is None


def test_aquifer_t_confounding(aq):
    a=analyze_network(aq,True)
    assert aq.final_objective<1e-10
    assert a.local.correlation[0][1]>.94
    assert all(s=='POORLY CONSTRAINED' for _,s,_ in a.local.statuses)
    assert coupling_pairs(aq)==((a.local.parameters[0],a.local.parameters[1]),)


@pytest.mark.parametrize('fixture',['nt','aq'])
def test_real_objective_valley_along_weak_combination(fixture,request):
    scenario=request.getfixturevalue(fixture)
    local=analyze_network(scenario).local
    weak=np.array(local.weak_combination)
    strong=np.array([-weak[1],weak[0]])
    specs,obj=context(scenario)
    values=np.array([s.initial for s in specs])
    scales=np.array([s.scale for s in specs])
    along=fixed_point(obj,values+.02*scales*weak)
    across=fixed_point(obj,values+.02*scales*strong)
    assert along.valid and across.valid
    assert along.objective<.1*across.objective


def test_unobserved_tank_is_indirect():
    f=B['truth_match']('unobserved')
    a=analyze_network(f,True)
    assert all(c[1]=='INDIRECTLY INFORMED' for c in a.coverage)
    assert a.local.intervals==(None,None)
    assert any('No direct pressure observations exist for B' in w for w in a.local.warnings)


@pytest.mark.parametrize('kind,rank',[('low_dp',0),('redundant',1)])
def test_zero_information_and_rank_deficiency(kind,rank):
    f=B['truth_match'](kind)
    a=analyze_network(f,True)
    assert a.local.rank==rank
    assert a.local.covariance is None
    assert all(i is None for i in a.local.intervals)
    assert all(s=='NON-IDENTIFIABLE' for _,s,_ in a.local.statuses)
    assert any('poorly observable' in w for w in a.local.warnings)


def test_independent_correlation_sign():
    # J'J=[[2,1],[1,1]], inverse=[[1,-1],[-1,2]], corr=-1/sqrt(2).
    d=matrix_diagnostics([[1.,0.],[1.,1.],[0.,0.]])
    assert d['correlation'][0,1]==pytest.approx(-1/np.sqrt(2))
    assert d['geometry']==pytest.approx(np.array([[1.,-1.],[-1.,2.]]))
    deficient=matrix_diagnostics([[1.,2.],[2.,4.],[3.,6.]])
    assert deficient['rank']==1 and deficient['geometry'] is None


def test_surface_real_objective_and_distinct_minimum(well):
    surface=network_surface(well,'tank:A:N','connection:A:B:T',3)
    specs,obj=context(well)
    valid=[p for p in surface.points if p.valid]
    minimum=min(valid,key=lambda p:p.objective)
    assert dict(minimum.parameters)['tank:A:N']==dict(well.fitted_parameters)['tank:A:N']
    assert minimum.objective<1e-10
    assert sum(p.objective<1. for p in valid)==1
    for point in surface.points[:3]:
        residual=evaluate(obj,[dict(point.parameters)[s.name] for s in specs])
        assert (residual is not None)==point.valid
        if point.valid:
            assert float(residual@residual)==pytest.approx(point.objective)
            assert dict(point.closure)['valid']


@pytest.fixture(scope='module')
def nt_profiles(nt):
    return (network_profile(nt,'connection:A:B:T',3,False,plausible_delta=3.84),
            network_profile(nt,'connection:A:B:T',3,True,plausible_delta=3.84))


def test_profile_compensation_and_ranges(nt,nt_profiles):
    fixed,profile=nt_profiles
    assert profile.plausible_ranges
    true=dict(nt.fitted_parameters)['connection:A:B:T']
    assert any(a<=true<=b for a,b in profile.plausible_ranges)
    improved=[]
    for a,b in zip(fixed.points,profile.points):
        if a.valid and b.valid and b.optimizer_status=='CONVERGED':
            assert b.objective<=a.objective+1e-6
            improved.append(a.objective-b.objective)
            assert dict(b.closure)['valid'] and b.tank_metrics
            assert len(b.bound_flags)==2
    assert max(improved)>10
    assert len({round(dict(p.parameters)['tank:A:N']) for p in profile.points if p.valid})>1


def test_reoptimized_profile_repeatability(nt,nt_profiles):
    assert network_profile(nt,'connection:A:B:T',3,True,plausible_delta=3.84)==nt_profiles[1]


def test_aquifer_profile_compensation(aq):
    fixed=network_profile(aq,'connection:A:B:T',3,False,plausible_delta=3.84)
    profile=network_profile(aq,'connection:A:B:T',3,True,plausible_delta=3.84)
    assert any(a.valid and b.valid and b.optimizer_status=='CONVERGED' and b.objective<a.objective-1 for a,b in zip(fixed.points,profile.points))
    values=[dict(p.parameters)['tank:B:aquifer:productivity_index'] for p in profile.points if p.valid]
    assert max(values)-min(values)>1e-10


def test_well_profile_width(well):
    p=network_profile(well,'connection:A:B:T',3,True,plausible_delta=3.84)
    true=dict(well.fitted_parameters)['connection:A:B:T']
    assert any(a<=true<=b for a,b in p.plausible_ranges)
    assert p.threshold_brackets
    assert max(b-a for a,b in p.plausible_ranges)<1e-9
    assert any('unresolved' in w for w in p.warnings)  # coarse node is not a zero-width interval


@pytest.mark.parametrize('stream',['np','gp','wp','winj','ginj'])
def test_general_n_conserving_stream_local_shifts(stream):
    from datetime import date
    d=date(2020,1,1)
    rows=tuple(AllocationRow(day,s,(('A',.6),('B',.3),('C',.1))) for s in ('np','gp','wp','winj','ginj') for day in (d,d+timedelta(days=30)))
    plan=HistoryPlan(('np','gp','wp','winj','ginj'),(),rows)
    changed=shift_allocation(plan,stream,'A','C',.05,[d+timedelta(days=30)])
    for old,new in zip(plan.schedules,changed.schedules):
        assert sum(dict(new.fractions).values())==pytest.approx(1.,abs=1e-15)
        if old.stream==stream and old.date!=d:
            assert dict(new.fractions)==pytest.approx({'A':.55,'B':.3,'C':.15})
        else:
            assert old==new


@pytest.mark.parametrize('donor,receiver,delta',[('B','A',.9),('A','B',-.9),('A','A',.1),('Z','B',.1),('A','B',float('nan'))])
def test_invalid_shifts_rejected(allocation,donor,receiver,delta):
    with pytest.raises(ValueError):
        shift_allocation(allocation.history_plan,'np',donor,receiver,delta)


def test_allocation_drift_and_conservation(allocation,allocation_results):
    assert all(c.match and c.match.converged for c in allocation_results.cases)
    middle=allocation_results.cases[2]
    assert dict(middle.match.final_network_metrics)['rmse']<1e-3
    assert all(abs(d[4])<1e-5 for d in middle.drift)
    assert allocation_results.robustness in ('MODERATELY SENSITIVE','HIGHLY SENSITIVE')
    for case in allocation_results.cases:
        prepared=prepare_history(allocation.base_network,case.plan,allocation.observations)
        assert max(abs(r.interval_conservation_error) for r in prepared.audit)<1e-10
        assert max(abs(r.cumulative_conservation_error) for r in prepared.audit)<1e-9
        assert case.transfers and case.drives
        assert all(dict(case.match.fitted_parameters)[s.name]>=s.lower for s in case.match.specifications)
    assert all(dict(c.match.final_network_metrics)['rmse']>dict(middle.match.final_network_metrics)['rmse'] for c in allocation_results.cases if c is not middle)
    for name in ('tank:A:N','tank:B:N','connection:A:B:T'):
        assert max(abs(next(d[4] for d in c.drift if d[0]==name)) for c in allocation_results.cases)>1.


def test_water_allocation_error_independent(allocation):
    plan=shift_allocation(allocation.history_plan,'wp','B','A',.15)
    assert [r for r in plan.schedules if r.stream=='np']==[r for r in allocation.history_plan.schedules if r.stream=='np']
    study=allocation_study(allocation,[('water error',plan)])
    assert study.cases[0].match.converged
    assert any(abs(d[4])>1 for d in study.cases[0].drift)


def test_failed_allocation_case_stays_unresolved(allocation):
    invalid=replace(allocation.history_plan,schedules=())
    result=allocation_study(allocation,[('invalid',invalid)])
    assert result.robustness=='UNRESOLVED'
    assert result.cases[0].match is None and result.cases[0].failure


def test_failed_surface_and_profiles_keep_holes(well,monkeypatch):
    import material_balance_studio.network_matching.objective as module
    monkeypatch.setattr(module,'simulate_network',lambda *a: (_ for _ in ()).throw(ValueError('aquifer failure in B')))
    surface=network_surface(well,'tank:A:N','connection:A:B:T',3)
    assert all(not p.valid and p.objective is None and p.failure_type=='aquifer failure' for p in surface.points)
    profile=network_profile(well,'connection:A:B:T',3,False,plausible_delta=1.)
    assert profile.plausible_ranges==()
    assert all(not p.valid for p in profile.points)


def test_repeated_diagnostics_immutable(well,information):
    before=asdict(well)
    assert analyze_network(well,True)==information
    a=network_profile(well,'connection:A:B:T',3,False)
    b=network_profile(well,'connection:A:B:T',3,False)
    assert a==b
    assert network_surface(well,'tank:A:N','connection:A:B:T',3)==network_surface(well,'tank:A:N','connection:A:B:T',3)
    assert asdict(well)==before


def test_allocation_repeatability_and_no_state_leak(allocation,allocation_results):
    before=asdict(allocation)
    original=allocation_results.cases[1]
    repeat=allocation_study(allocation,[(original.name,original.plan)])
    assert repeat.cases[0]==original
    assert asdict(allocation)==before


def field_roundtrip(scenario):
    from material_balance_studio.io.history import history_from_frame,pvt_from_frame
    from material_balance_studio.pvt.table_model import TablePVTModel
    from material_balance_studio.units.conversions import from_si
    from material_balance_studio.domain.models import CUMULATIVE_FIELDS
    from material_balance_studio.presentation.aquifer import FIELDS
    def histories(records):
        return history_from_frame(pd.DataFrame([dict(date=str(h.date),**{n:from_si(getattr(h.cumulative,n),'gas_volume' if n in ('gp','ginj') else 'liquid_volume','FIELD') for n in CUMULATIVE_FIELDS}) for h in records]),'FIELD') if records else ()
    nodes=[]
    for node in scenario.base_network.tanks:
        tank=node.reservoir
        pvt=pvt_from_frame(pd.DataFrame([{n:from_si(getattr(p,n),n,'FIELD') for n in ('pressure','bo','rs','bg','bw')} for p in tank.pvt_model.table.rows]),'FIELD')
        aquifer=tank.aquifer
        if aquifer is not None:
            aquifer=replace(aquifer,**{k:from_display(to_display(v,FIELDS[k][1],'FIELD'),FIELDS[k][1],'FIELD') for k,v in asdict(aquifer).items()})
        reservoir=replace(tank,oil_in_place=from_display(to_display(tank.oil_in_place,'oil_volume','FIELD'),'oil_volume','FIELD'),
                          initial_pressure=from_display(to_display(tank.initial_pressure,'pressure','FIELD'),'pressure','FIELD'),pvt_model=TablePVTModel(pvt),aquifer=aquifer)
        nodes.append(replace(node,reservoir=reservoir,history=histories(node.history)))
    network=replace(scenario.base_network,tanks=tuple(nodes),connections=tuple(replace(e,transmissibility=from_display(to_display(e.transmissibility,'aquifer_productivity','FIELD'),'aquifer_productivity','FIELD')) for e in scenario.base_network.connections))
    plan=replace(scenario.history_plan,field_history=histories(scenario.history_plan.field_history))
    specs=tuple(replace(s,**{k:from_display(to_display(getattr(s,k),s.unit,'FIELD'),s.unit,'FIELD') for k in ('initial','lower','upper','scale')}) for s in scenario.specifications)
    obs=parse_observations(observation_input(scenario.observations,'FIELD'),'FIELD')
    return history_match_network(network,specs,obs,plan,scenario.mode,scenario.default_sigma,scenario.settings)


def test_field_si_information_and_grids(well,information):
    field=field_roundtrip(well)
    info=analyze_network(field,True)
    assert np.array(info.local.scaled_sensitivity)==pytest.approx(np.array(information.local.scaled_sensitivity),rel=3e-4,abs=1e-8)
    assert np.array(info.local.correlation)==pytest.approx(np.array(information.local.correlation),abs=1e-4)
    assert [s for _,s,_ in info.local.statuses]==[s for _,s,_ in information.local.statuses]
    si=sensitivity_frame(information.local,well,'SI')
    fi=sensitivity_frame(info.local,field,'FIELD')
    for j,s in enumerate(field.specifications):
        factor=to_display(1.,'pressure','FIELD')/to_display(1.,s.unit,'FIELD')/(to_display(1.,'pressure','SI')/to_display(1.,s.unit,'SI'))
        assert fi.iloc[:,j].to_numpy()==pytest.approx(si.iloc[:,j].to_numpy()*factor,rel=3e-4,abs=1e-6)
    for run,args in [(network_profile,('connection:A:B:T',3,False)),(network_surface,('tank:A:N','connection:A:B:T',3))]:
        a,b=run(well,*args),run(field,*args)
        assert [p.valid for p in a.points]==[p.valid for p in b.points]
        assert [p.objective for p in a.points if p.valid]==pytest.approx([p.objective for p in b.points if p.valid],rel=1e-5,abs=1e-6)


def test_field_si_allocation_scenario(allocation,allocation_results):
    field=field_roundtrip(allocation)
    case=allocation_results.cases[1]
    result=allocation_study(field,[(case.name,replace(case.plan,field_history=field.history_plan.field_history))])
    assert [v for _,v in result.cases[0].match.fitted_parameters]==pytest.approx([v for _,v in case.match.fitted_parameters],rel=2e-4)
    assert result.cases[0].robustness==case.robustness


def test_protected_571_tests_and_sources():
    for name,digest in json.loads((ROOT/'docs/phase_5c_protected_baseline.json').read_text()).items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest


def test_ui_information_and_unit_switch(well,information):
    app=AppTest.from_file(str(ROOT/'app.py'),default_timeout=120).run()
    app.session_state['network_match_scenarios']=[well]
    app.session_state['network_identifiability']={well.scenario_id:{'information':information}}
    app.run()
    assert not app.exception and not app.error
    app.selectbox(key='display_units').select('FIELD').run()
    assert not app.exception and not app.error
    assert app.session_state['network_match_scenarios']==[well]


def test_ui_allocation_preparation_is_separate(allocation):
    app=AppTest.from_file(str(ROOT/'app.py'),default_timeout=120).run()
    app.session_state['network_match_scenarios']=[allocation]
    app.run()
    app.button(key='nu_prepare').click().run()
    assert not app.exception and not app.error
    saved=app.session_state['network_identifiability'][allocation.scenario_id]
    assert len(saved['allocation_plans'])==5
    assert app.session_state['network_match_scenarios']==[allocation]
