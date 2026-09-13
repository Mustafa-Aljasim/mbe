"""Temporary independently calibrated aquifer scenarios; no Apply operation."""
from dataclasses import replace,asdict
import json
import pandas as pd
import streamlit as st
from material_balance_studio.matching import MatchParameter,observations_from_history
from material_balance_studio.matching.parameters import parameter_registry,parameter_value
from material_balance_studio.units.display import to_display,from_display,display_unit
from .aquifer_calibration import calibrate_aquifers
from .history_matching import fitted_parameter_frame


def refresh_calibration_units():
    units=st.session_state['display_units']
    for prefix,draft in st.session_state.get('cmp_fit_drafts',{}).items():
        for field in ('Initial','Lower','Upper'):
            st.session_state[prefix+'_'+field]=to_display(draft[field],draft['quantity'],units)
    if 'cmp_fit_sigma_si' in st.session_state:
        st.session_state['cmp_fit_sigma']=to_display(st.session_state['cmp_fit_sigma_si'],'pressure',units)


def save_bound(prefix,field):
    draft=st.session_state['cmp_fit_drafts'][prefix]
    draft[field]=from_display(st.session_state[prefix+'_'+field],draft['quantity'],st.session_state['display_units'])


def save_sigma():
    st.session_state['cmp_fit_sigma_si']=from_display(st.session_state['cmp_fit_sigma'],'pressure',st.session_state['display_units'])


def calibrated_workflow(tank_factory,history,models,units,enabled):
    st.caption('Independently adjusts the selected aquifer parameters within their engineering bounds against the same observed pressure history before comparing the models. N, m, PVT, production/injection and compressibilities remain fixed.')
    if not enabled:
        st.info('Provide valid reservoir/PVT/history inputs and select comparison models.')
        return
    tank=tank_factory()
    specifications={}
    invalid=False
    for name,model in models.items():
        candidate=replace(tank,aquifer=model)
        registry={n:v for n,v in parameter_registry(candidate).items() if n.startswith('aquifer.')}
        with st.expander(name+' matching controls',expanded=True):
            if not registry:
                st.caption('No Aquifer: fixed no-aquifer baseline; no parameters to optimize.')
                specifications[name]=()
                continue
            key='cmp_fit_active_'+name
            selected=st.multiselect(name+' aquifer parameters allowed to change',list(registry),key=key)
            specs=[]
            for parameter in selected:
                label,quantity=registry[parameter]
                initial=parameter_value(candidate,parameter)
                prefix=f'cmp_fit_{name}_{parameter}_{initial}'
                st.write(label+' · '+parameter)
                values=[]
                upper=initial*4
                if parameter.endswith('porosity'): upper=min(.99,upper)
                if parameter.endswith('encroachment_angle'): upper=min(360.,upper)
                draft=st.session_state.setdefault('cmp_fit_drafts',{}).setdefault(prefix,dict(Initial=initial,Lower=initial*.25,Upper=upper,quantity=quantity))
                for field in ('Initial','Lower','Upper'):
                    st.session_state.setdefault(prefix+'_'+field,to_display(draft[field],quantity,units))
                    st.number_input(f'{field} [{display_unit(quantity,units)}]',format='%.8g',key=prefix+'_'+field,on_change=save_bound,args=(prefix,field))
                    values.append(draft[field])
                transform=st.selectbox('Transformation',['linear','log'],key=prefix+'_transform')
                try:
                    specs.append(MatchParameter(parameter,label,*values,quantity,initial,transform))
                except ValueError as exc:
                    invalid=True
                    st.error(str(exc))
            specifications[name]=tuple(specs)
    obs=observations_from_history(history,st.session_state.get('pressure_metadata',()))
    selected=[]
    with st.expander('Common pressure observations and quality',expanded=True):
        for o in obs:
            include=st.checkbox(str(o.date),value=True,key='cmp_fit_include_'+str(o.date))
            selected.append(replace(o,include_in_match=include))
        st.dataframe(pd.DataFrame([dict(Date=o.date,Pressure=to_display(o.pressure,'pressure',units),Sigma=to_display(o.sigma,'pressure',units),Source=o.source,QC=o.qc_flag,Included=o.include_in_match,Note=o.note) for o in selected]),hide_index=True)
        st.caption('Pressure and sigma: '+display_unit('pressure',units)+'. All models use these same inclusion and weighting choices.')
    mode=st.radio('Comparison residual weighting',['Unweighted','Pressure uncertainty weighted'],key='cmp_fit_weighting')
    default=None
    if mode!='Unweighted' and st.checkbox('Supply default sigma for missing uncertainties',key='cmp_fit_sigma_enabled'):
        st.session_state.setdefault('cmp_fit_sigma_si',1e5)
        st.session_state.setdefault('cmp_fit_sigma',to_display(st.session_state['cmp_fit_sigma_si'],'pressure',units))
        st.number_input('Default sigma ['+display_unit('pressure',units)+']',min_value=1e-12,key='cmp_fit_sigma',on_change=save_sigma)
        default=st.session_state['cmp_fit_sigma_si']
    limit=st.number_input('Per-model optimizer evaluation limit',min_value=1,max_value=1000,value=100,key='cmp_fit_limit')
    st.caption('The model with the lowest pressure error is not automatically the physically correct aquifer. No selected aquifer parameters means a fixed comparison run.')
    if st.button('History-match aquifers and compare',key='cmp_fit_run',disabled=invalid or not selected):
        try:
            cases=calibrate_aquifers(tank,history,models,specifications,tuple(selected),mode,default,int(limit))
            st.session_state['calibrated_aquifer_comparison']=(cases,tank,tuple(history),tuple(selected),mode,default)
        except (ValueError,ArithmeticError) as exc:
            st.error(str(exc))
    saved=st.session_state.get('calibrated_aquifer_comparison')
    if not saved: return
    cases,base,saved_history,observations,weighting,sigma=saved
    st.caption(f'Saved temporary comparison · {weighting}. Rerun after editing inputs. Reservoir Setup and accepted single/network scenarios are unchanged.')
    eligible=[c for c in cases if c.objective is not None and dict(c.pressure_metrics)['rmse'] is not None]
    if eligible:
        best=min(eligible,key=lambda c:dict(c.pressure_metrics)['rmse'])
        st.info('Best numerical pressure fit: '+best.run.name+'. Engineering QC is reported separately.')
    rows=[]
    for c in cases:
        states=c.run.result.states if c.run.result else ()
        rows.append(dict(Model=c.run.name,Comparison_Mode=c.mode,**{n:to_display(v,'pressure',units) for n,v in c.pressure_metrics},
            Pressure_unit=display_unit('pressure',units),Final_We=to_display(states[-1].balance.aquifer_support,'reservoir_volume',units) if states else None,
            We_unit=display_unit('reservoir_volume',units),Maximum_MBE_relative=max((s.balance.relative_residual for s in states),default=None),
            Bound_Status=str(c.fit.bound_flags) if c.fit else 'Fixed',Engineering_QC=c.qc,Notes=' | '.join(c.warnings)))
    summary=pd.DataFrame(rows)
    st.dataframe(summary,hide_index=True)
    st.download_button('Download calibrated comparison summary',summary.to_csv(index=False),'calibrated_aquifer_comparison.csv')
    for c in cases:
        with st.expander(c.run.name+' · '+c.qc,expanded=c.qc!='PASS'):
            if c.fit: st.dataframe(fitted_parameter_frame(c.fit,units),hide_index=True)
            for warning in c.warnings:
                (st.error if c.qc=='FAIL' else st.warning)(warning)
    from .aquifer_comparison import comparison_figures
    for i,figure in enumerate(comparison_figures([c.run for c in cases],units)):
        st.plotly_chart(figure,key=f'cmp_calibrated_plot_{i}')
    st.caption('Observation markers retain measured data. Summary metrics and objectives use only included observations; missing sigma follows the existing matching policy.')
    st.download_button('Download calibrated comparison record',json.dumps([asdict(c) for c in cases],default=str,indent=2),'calibrated_aquifer_comparison.json')
