"""Post-5C workflow acceptance, with unchanged engineering engines."""
from dataclasses import replace, asdict
from pathlib import Path
import hashlib
import json
import runpy
import pandas as pd
import pytest
from streamlit.testing.v1 import AppTest
from material_balance_studio.presentation.engineering_inputs import (
    template_frame, normalize_file_frame, sparse_lab_input, laboratory_file_frame,
    informational_composition)
from material_balance_studio.presentation.history_plots import history_figure, history_values
from material_balance_studio.presentation.pressure_oil import pressure_oil_figure, without_aquifer_reference
from material_balance_studio.presentation.aquifer_calibration import calibrate_aquifers
from material_balance_studio.io.history import pvt_from_frame, history_from_frame
from material_balance_studio.aquifer import NoAquifer, PotAquifer, CarterTracyAquifer
from material_balance_studio.solver.simulation import simulate
from material_balance_studio.matching import observations_from_history

ROOT = Path(__file__).resolve().parents[1]
B = runpy.run_path(str(ROOT/'examples/history_match_benchmarks.py'))


@pytest.mark.parametrize('kind,parser', [('PVT',pvt_from_frame),('History',history_from_frame)])
def test_engineering_templates_canonical_equivalence(kind,parser):
    si,field = template_frame(kind,'SI'),template_frame(kind,'FIELD')
    assert any('[MPa]' in c for c in si) and any('[psia]' in c for c in field)
    if kind=='History':
        assert {'np [STB]','gp [scf]','wp [bbl]','winj [bbl]','ginj [scf]'} <= set(field)
    a,b = parser(normalize_file_frame(si,'SI'),'SI'),parser(normalize_file_frame(field,'FIELD'),'FIELD')
    def numbers(value):
        if isinstance(value,dict): return [n for v in value.values() for n in numbers(v)]
        if isinstance(value,(tuple,list)): return [n for v in value for n in numbers(v)]
        return [value] if isinstance(value,(int,float)) else []
    assert numbers(asdict(a) if kind=='PVT' else [asdict(h) for h in a]) == pytest.approx(numbers(asdict(b) if kind=='PVT' else [asdict(h) for h in b]))


@pytest.mark.parametrize('units',['SI','FIELD'])
def test_sparse_lab_roundtrip_and_strict_final(units):
    raw=pd.DataFrame(dict(pressure=[1e7,2e7],bo=[1.2,1.3],rs=[40,80],gas_viscosity=[None,1e-5]))
    original=sparse_lab_input(raw)
    result=sparse_lab_input(laboratory_file_frame(original.lab,original.gas_viscosity,units),units)
    assert [p.bo for p in result.lab.points]==[1.2,1.3]
    assert all(p.z is None and p.oil_viscosity is None for p in result.lab.points)
    assert result.gas_viscosity[0]==pytest.approx((2e7,1e-5))
    with pytest.raises(ValueError): pvt_from_frame(raw.drop(columns='gas_viscosity'))


def test_tagged_unit_conflict_and_legacy_pressure():
    raw=pd.DataFrame({'pressure':[2e7]})
    assert normalize_file_frame(raw,'SI').equals(raw)
    with pytest.raises(ValueError): normalize_file_frame(pd.DataFrame({'pressure [psia]':[2000]}),'SI')
    assert informational_composition(1,2,3)['co2_mole_percent']==2
    with pytest.raises(ValueError): informational_composition(50,40,20)


@pytest.mark.parametrize('units',['SI','FIELD'])
@pytest.mark.parametrize('x,y', [('Date','Cumulative Oil Production'),('Date','Cumulative Water Injection'),('Date','Observed Reservoir Pressure'),('Cumulative Oil Production','Observed Reservoir Pressure'),('Cumulative Water Production','Cumulative Oil Production')])
def test_history_pairs(x,y,units):
    _,history=B['benchmark_case']()
    history=tuple(replace(h,observed_pressure=None) if i%2 else h for i,h in enumerate(history))
    before=history
    fig,paired=history_figure(history,x,y,units)
    assert len(paired)==(4 if y=='Observed Reservoir Pressure' else 8)
    assert len(fig.data[0].x)==len(paired)
    if y=='Observed Reservoir Pressure':
        assert ('MPa' if units=='SI' else 'psia') in fig.layout.yaxis.title.text
        assert fig.data[0].mode=='markers'
    assert history==before
    assert history_figure(history,x,y,units,'Line')[0].data[0].mode=='lines+markers'
    assert history_figure(history,x,y,units,'Scatter')[0].data[0].mode=='markers'
    assert list(history_values(history).rp)==pytest.approx([80]*8)


@pytest.mark.parametrize('model,name',[(NoAquifer(),'None'),(B['benchmark_case'](True)[0].aquifer,'Fetkovich'),(CarterTracyAquifer(1000,10,20,.2,1e-13,.001,1e-9),'Carter-Tracy')])
def test_pressure_oil_reference(model,name):
    tank,history=B['benchmark_case']()
    tank=replace(tank,aquifer=model)
    active=simulate(tank,history)
    snapshot=asdict(active)
    reference_tank,reference=without_aquifer_reference(tank,history)
    assert replace(reference_tank,aquifer=model)==tank
    assert all(s.balance.aquifer_support==0 for s in reference.states)
    fig,_=pressure_oil_figure(active,tank,history,'FIELD',reference)
    assert [t.name for t in fig.data]==['Observed Pressure','Calculated — Without Aquifer']+([] if name=='None' else ['Calculated — With '+name])
    assert asdict(active)==snapshot and tank.aquifer==model


@pytest.mark.parametrize('family',['Pot','Fetkovich'])
def test_independent_aquifer_calibration(family):
    tank,history=B['benchmark_case'](family=='Fetkovich')
    if family=='Pot':
        truth=.001
        tank=replace(tank,aquifer=PotAquifer(truth))
        records=[]
        for h in history:
            dp=tank.initial_pressure-h.observed_pressure
            oil=(tank.oil_in_place*4e-9*dp+truth*dp)/(1.2+4e-9*dp)
            records.append(replace(h,cumulative=replace(h.cumulative,np=oil,gp=80*oil)))
        history=tuple(records)
        parameter='aquifer.capacity'
        initial=PotAquifer(truth*.5)
    else:
        truth=1e-9
        parameter='aquifer.productivity_index'
        initial=replace(tank.aquifer,productivity_index=truth*.5)
    candidate=replace(tank,aquifer=initial)
    spec=B['specification'](candidate,parameter,truth*.1,truth*2)
    observations=observations_from_history(history)
    observations=tuple(replace(o,include_in_match=i!=0) for i,o in enumerate(observations))
    cases=calibrate_aquifers(tank,history,{family:initial,'None':NoAquifer()},{family:[spec]},observations,'Pressure uncertainty weighted',1e5)
    fit=cases[0].fit
    assert fit.converged and dict(fit.fitted_parameters)[parameter]==pytest.approx(truth,rel=2e-4)
    assert dict(cases[0].pressure_metrics)['rmse']<10
    assert fit.fitted_tank.oil_in_place==tank.oil_in_place and fit.fitted_tank.m==tank.m
    assert fit.fitted_tank.pvt_model==tank.pvt_model and tank.aquifer!=initial
    assert cases[1].fit is None and cases[1].run.model.key=='none'
    assert max(s.balance.relative_residual for s in fit.final_simulation.states)<1e-6
    limited=replace(spec,upper=truth*.7)
    bound=calibrate_aquifers(tank,history,{family:initial},{family:[limited]})[0]
    assert bound.fit.converged and bound.fit.bound_flags and bound.qc=='CAUTION'
    rejected=calibrate_aquifers(tank,history,{family:initial},{family:[B['specification'](tank,'oil_in_place')]})[0]
    assert rejected.qc=='FAIL' and 'N and m remain fixed' in rejected.warnings[0]


def test_post5c_protected_sources():
    baseline=json.loads((ROOT/'docs/post5c_baseline/source_hashes.json').read_text())
    intentional=json.loads((ROOT/'docs/post5c_intentional_changes.json').read_text())
    for name,digest in baseline.items():
        actual=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        if name in intentional:
            assert name=='app.py' or name.startswith('src/material_balance_studio/presentation/')
            assert intentional[name]['before']==digest
            assert intentional[name]['after']==actual
        else: assert actual==digest,name


def test_comparison_controls_and_unit_switch():
    app=AppTest.from_file(str(ROOT/'app.py'),default_timeout=60).run()
    assert app.radio(key='comparison_mode').value=='Compare with Supplied Parameters'
    app.multiselect(key='comparison_models').set_value(['None','Pot']).run()
    app.radio(key='comparison_mode').set_value('History-Match Each Aquifer Before Comparison').run()
    app.multiselect(key='cmp_fit_active_Pot').set_value(['aquifer.capacity']).run()
    assert not app.exception
    widget=next(w for w in app.number_input if w.key and w.key.startswith('cmp_fit_Pot_') and w.key.endswith('_Lower'))
    widget.set_value(widget.value*1.1).run()
    drafts=json.loads(json.dumps(app.session_state['cmp_fit_drafts']))
    app.selectbox(key='display_units').select('FIELD').run()
    assert not app.exception and app.session_state['cmp_fit_drafts']==drafts
    app.selectbox(key='display_units').select('SI').run()
    assert not app.exception and app.session_state['cmp_fit_drafts']==drafts
    app.button(key='cmp_fit_run').click().run()
    assert not app.exception
    assert 'calibrated_aquifer_comparison' in app.session_state


def test_new_upload_default_and_persistent_source_units(monkeypatch):
    from types import SimpleNamespace
    import streamlit
    from material_balance_studio.presentation.engineering_inputs import default_upload_units
    state={}
    monkeypatch.setattr(streamlit,'session_state',state)
    upload=SimpleNamespace(name='history.csv',getvalue=lambda:b'first source')
    default_upload_units(upload,'file_units','FIELD')
    assert state['file_units']=='FIELD'
    state['file_units']='SI'
    default_upload_units(upload,'file_units','FIELD')
    assert state['file_units']=='SI'
    default_upload_units(upload,'file_units','SI')
    assert state['file_units']=='SI'
    upload.getvalue=lambda:b'new source'
    default_upload_units(upload,'file_units','FIELD')
    assert state['file_units']=='FIELD'


def test_comparison_observation_policy():
    tank,history=B['benchmark_case']()
    obs=observations_from_history(history)
    with pytest.raises(ValueError,match='sigma'):
        calibrate_aquifers(tank,history,{'None':NoAquifer()},{},obs,'Pressure uncertainty weighted')
    with pytest.raises(ValueError,match='Include'):
        calibrate_aquifers(tank,history,{'None':NoAquifer()},{},[replace(o,include_in_match=False) for o in obs])
    with pytest.raises(ValueError,match='exact'):
        calibrate_aquifers(tank,history,{'None':NoAquifer()},{},[replace(obs[0],pressure=1e6)])



