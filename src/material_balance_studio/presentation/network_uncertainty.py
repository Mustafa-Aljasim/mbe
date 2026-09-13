"""Fixed-allocation information and separate deterministic allocation comparisons."""
from dataclasses import asdict, replace
import json
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from material_balance_studio.uncertainty.network_sensitivity import analyze_network,coupling_pairs
from material_balance_studio.uncertainty.network_profiles import network_profile,network_surface
from material_balance_studio.uncertainty.allocation_sensitivity import shift_allocation,allocation_study
from material_balance_studio.units.display import to_display,display_unit
from .uncertainty import sensitivity_frame,point_frame
from .network_matching import parse_schedule,parameter_frame


def network_identifiability_workflow(units):
    st.subheader('Network Identifiability')
    st.info('A good pressure fit does not prove unique N, m, T or aquifer strength. Fixed-allocation parameter information and allocation-assumption sensitivity are separate diagnostics.')
    scenarios=st.session_state.get('network_match_scenarios',[])
    if not scenarios:
        st.caption('Create a converged scenario in Network History & Matching first.')
        return
    i=st.selectbox('Network diagnostic scenario',range(len(scenarios)),format_func=lambda i:f'{i+1}: {scenarios[i].scenario_name}',key='nu_scenario')
    scenario=scenarios[i]
    if not scenario.converged:
        st.warning('Select a converged scenario.')
        return
    key=scenario.scenario_id
    saved=st.session_state.setdefault('network_identifiability',{}).setdefault(key,{})
    specs={s.name:s for s in scenario.specifications if s.active}
    names=list(specs)
    assumptions=st.checkbox('Accept independent Gaussian error and local-linearity assumptions for approximate confidence',key='nu_assumptions_'+key)
    st.caption('Weighted: sigma is known. Unweighted: common variance is estimated. These assumptions are not established by the fit.')
    if st.button('Run network sensitivity and identifiability',key='nu_run'):
        try:
            with st.spinner('Evaluating real network sensitivity stencils…'):
                saved['information']=analyze_network(scenario,assumptions)
        except (ValueError,ArithmeticError) as exc:
            st.error(str(exc))
    info=saved.get('information')
    st.metric('Network match RMSE',f"{to_display(dict(scenario.final_network_metrics)['rmse'],'pressure',units):.5g} {display_unit('pressure',units)}")
    if info:
        result=info.local
        st.caption(f'Fixed allocation • {result.objective_definition} • recorded statistical assumptions: {result.statistical_assumptions}')
        st.dataframe(pd.DataFrame(result.statuses,columns=['Parameter','Identifiability','Evidence']),hide_index=True)
        for warning in result.warnings:
            st.warning(warning)
        with st.expander('Sensitivity and observation information',expanded=True):
            raw=sensitivity_frame(result,scenario,units)
            raw.index=[f'{t} · {d}' for t,d in zip(info.observation_tanks,result.dates)]
            st.dataframe(raw)
            st.caption('Physical dP/dparameter in displayed units; dimensionless sensitivity = dP/dparameter × physical parameter scale / observed tank Pi.')
            st.plotly_chart(go.Figure(go.Heatmap(z=np.abs(result.scaled_sensitivity),x=result.parameters,y=list(raw.index),colorscale='Viridis')),key='nu_information')
            st.caption('Objective information below retains each survey’s sigma weighting (or the fixed 1 MPa unweighted normalization).')
            st.dataframe(pd.DataFrame(np.abs(result.scaled_jacobian),columns=result.parameters,index=list(raw.index)))
            st.dataframe(pd.DataFrame(dict(info.by_tank),index=result.parameters).T)
            st.dataframe(pd.DataFrame(info.coverage,columns=['Parameter','Observation coverage','RMS dimensionless sensitivity']),hide_index=True)
            st.caption('Coverage labels describe direct/indirect information; rank and sensitivity determine whether this information separates parameters.')
        with st.expander('Parameter correlation',expanded=True):
            if result.correlation is None:
                st.warning('Rank deficient: covariance correlation is withheld. Scaled Jacobian column coupling is shown as a different diagnostic.')
                corr=result.column_coupling
            else:
                corr=result.correlation
            st.plotly_chart(go.Figure(go.Heatmap(z=corr,x=names,y=names,zmin=-1,zmax=1,colorscale='RdBu')),key='nu_corr')
            rows=[dict(Parameter_A=a,Parameter_B=b,Correlation=corr[j][k],QC='STRONG COUPLING' if abs(corr[j][k])>.9 else 'MODERATE COUPLING' if abs(corr[j][k])>=.7 else 'LOW COUPLING') for j,a in enumerate(names) for k,b in enumerate(names) if k>j]
            st.dataframe(pd.DataFrame(sorted(rows,key=lambda r:-abs(r['Correlation']))),hide_index=True)
            st.caption('Thresholds: |correlation| ≥0.7 moderate, >0.9 strong. Compensation is not causation.')
        with st.expander('Approximate confidence'):
            st.write('LOCAL / LINEARIZED APPROXIMATION — fixed allocation only. Missing intervals are withheld, not zero uncertainty.')
            st.dataframe(pd.DataFrame([dict(Parameter=n,Unit=display_unit(specs[n].unit,units),
                Lower=None if interval is None else to_display(interval[0],specs[n].unit,units),
                Upper=None if interval is None else to_display(interval[1],specs[n].unit,units)) for n,interval in zip(names,result.intervals)]),hide_index=True)
        with st.expander('SVD and advanced diagnostics'):
            st.write(dict(rank=result.rank,active_parameters=len(names),condition=result.condition,singular_values=result.singular_values))
            st.dataframe(pd.DataFrame(sorted(zip(names,result.weak_combination),key=lambda p:-abs(p[1])),columns=['Parameter','Weak singular-vector coefficient']),hide_index=True)
            st.caption('Dominant terms describe a weak scaled parameter combination permitted by the local pressure response; no causal interpretation.')
            st.dataframe(pd.DataFrame(result.schemes,columns=['Parameter','Stencil','Physical step']),hide_index=True)
            st.write('Derivative step-halving errors',result.derivative_errors)
            st.dataframe(pd.DataFrame(result.residual_jacobian,columns=result.parameters,index=[f'{t} · {d}' for t,d in zip(info.observation_tanks,result.dates)]))
            st.caption('Residual Jacobian dr/dx in canonical physical parameter units. Full validity and failure evaluations are in the export.')
        with st.expander('Connection observability'):
            for entry in info.connections:
                row=dict(entry)
                st.write(str(row['connection']))
                dates,dp,q=zip(*row['series'])
                fig=go.Figure(go.Scatter(x=dates,y=[to_display(v,'pressure',units) for v in dp],name='ΔP'))
                fig.add_trace(go.Scatter(x=dates,y=q,name='Transfer rate [m³/s]',yaxis='y2'))
                fig.update_layout(yaxis_title=f"ΔP [{display_unit('pressure',units)}]",yaxis2=dict(title='q [m³/s]',overlaying='y',side='right'))
                st.plotly_chart(fig,key='nu_edge_'+str(row['connection']))
                for warning in row['warnings']:
                    st.warning(warning)
    with st.expander('Objective surfaces'):
        if len(names)<2:
            st.caption('Select a match with at least two active parameters.')
        else:
            presets=coupling_pairs(scenario)
            chosen=st.selectbox('Coupling preset',['Custom',*presets],format_func=str,key='nu_preset_'+key)
            x=st.selectbox('Parameter X',names,key='nu_x_'+key)
            y=st.selectbox('Parameter Y',[n for n in names if n!=x],key='nu_y_'+key)
            if chosen!='Custom':
                x,y=chosen
            count=st.selectbox('Surface grid', [3,5,9],format_func=lambda n:{3:'Coarse (3)',5:'Medium (5)',9:'Fine (9)'}[n],key='nu_grid_'+key)
            st.caption(f'Up to {(count+1)**2} forward runs, including the matched point; other parameters fixed.')
            if st.button('Run network objective surface',key='nu_surface_run'):
                try:
                    saved['surface']=network_surface(scenario,x,y,count)
                except (ValueError,ArithmeticError) as exc:
                    st.error(str(exc))
            surface=saved.get('surface')
            if surface:
                x,y=surface.parameters
                xs=[to_display(v,specs[x].unit,units) for v in surface.x]
                ys=[to_display(v,specs[y].unit,units) for v in surface.y]
                z=np.array([p.objective if p.valid else np.nan for p in surface.points]).reshape(len(ys),len(xs))
                fig=go.Figure(go.Heatmap(x=xs,y=ys,z=z,connectgaps=False,colorbar=dict(title='Objective')))
                for label,point in [('Initial',surface.initial),('Matched',surface.matched)]:
                    fig.add_trace(go.Scatter(x=[to_display(point[0],specs[x].unit,units)],y=[to_display(point[1],specs[y].unit,units)],mode='markers',name=label))
                bad=[p for p in surface.points if not p.valid]
                if bad:
                    fig.add_trace(go.Scatter(x=[to_display(dict(p.parameters)[x],specs[x].unit,units) for p in bad],y=[to_display(dict(p.parameters)[y],specs[y].unit,units) for p in bad],mode='markers',marker_symbol='x',name='Failed'))
                fig.update_layout(xaxis_title=f'{x} [{display_unit(specs[x].unit,units)}]',yaxis_title=f'{y} [{display_unit(specs[y].unit,units)}]')
                st.plotly_chart(fig,key='nu_surface')
                st.caption('Grid edges are the original user bounds. Failed regions are unresolved; no interpolated validity.')
                st.dataframe(point_frame(surface.points,scenario,units),hide_index=True)
    with st.expander('Parameter profiles'):
        parameter=st.selectbox('Profile parameter',names,key='nu_profile_'+key)
        count=st.selectbox('Profile grid',[3,5,9],key='nu_profile_count_'+key)
        reopt=st.checkbox('Re-optimize other parameters within original bounds',value=True,key='nu_reopt_'+key)
        budget=st.number_input('Profile optimizer evaluation limit',min_value=1,max_value=1000,value=40,key='nu_budget_'+key)
        use_delta=st.checkbox('Use a deterministic plausible objective increment',key='nu_delta_enable_'+key)
        delta=st.number_input('Plausible objective increment',min_value=1e-12,value=.01,format='%.8g',key='nu_delta_'+key)
        st.caption(f'Up to {count+1} profile points. Re-optimization cost varies; each three-point Jacobian adds about {2*max(0,len(names)-1)} forward calls per optimizer step.')
        if st.button('Run network parameter profile',key='nu_profile_run'):
            try:
                saved['profile']=network_profile(scenario,parameter,count,reopt,int(budget),assumptions,delta if use_delta else None)
            except (ValueError,ArithmeticError) as exc:
                st.error(str(exc))
        profile=saved.get('profile')
        if profile:
            name=profile.parameter
            fig=go.Figure(go.Scatter(x=[to_display(dict(p.parameters)[name],specs[name].unit,units) for p in profile.points],
                y=[p.objective if p.valid and p.optimizer_status in ('CONVERGED','FIXED OTHERS') else None for p in profile.points],mode='lines+markers',connectgaps=False))
            if profile.threshold is not None:
                fig.add_hline(y=profile.threshold)
            fig.update_layout(xaxis_title=f'{name} [{display_unit(specs[name].unit,units)}]',yaxis_title='Objective')
            st.plotly_chart(fig,key='nu_profile_plot')
            st.write(profile.threshold_method)
            st.write('Grid-resolved plausible ranges',[(to_display(a,specs[name].unit,units),to_display(b,specs[name].unit,units)) for a,b in profile.plausible_ranges])
            st.dataframe(point_frame(profile.points,scenario,units),hide_index=True)
            st.dataframe(pd.DataFrame([dict(Parameters_SI=str(p.parameters),Tank_metrics_Pa=str(p.tank_metrics),Closure=str(p.closure),Affected=str(p.affected),Failure=p.failure) for p in profile.points]),hide_index=True)
            for warning in profile.warnings:
                st.warning(warning)
    with st.expander('Allocation assumption sensitivity',expanded=True):
        st.info('This compares separate fits under fixed, engineer-selected schedules. It does not fit or infer allocation fractions, and these ranges are not confidence intervals.')
        if not scenario.history_plan.allocated_streams:
            st.caption('This match uses direct histories. Create a field-allocation match to compare allocation schedules.')
        else:
            stream=st.selectbox('Allocation stream',scenario.history_plan.allocated_streams,key='nu_stream_'+key)
            tanks=[t.name for t in scenario.base_network.tanks]
            donor=st.selectbox('Transfer fraction from',tanks,key='nu_donor_'+key)
            receiver=st.selectbox('Transfer fraction to',[n for n in tanks if n!=donor],key='nu_receiver_'+key) if len(tanks)>1 else None
            dates=sorted(r.date for r in scenario.history_plan.schedules if r.stream==stream)
            selected=st.multiselect('Effective intervals (empty means entire history)',dates,key='nu_dates_'+key)
            shifts=st.text_input('Fraction shifts, comma separated',value='-0.10,-0.05,0,0.05,0.10',key='nu_shifts_'+key)
            case_plans=saved.setdefault('allocation_plans',{})
            if st.button('Prepare allocation shifts',key='nu_prepare') and receiver:
                try:
                    values=[float(v.strip()) for v in shifts.split(',')]
                    if len(values)>21:
                        raise ValueError('At most 21 cases.')
                    prepared={f'{stream}: {donor} → {receiver}, δ={v:g}':shift_allocation(scenario.history_plan,stream,donor,receiver,v,selected or None) for v in values}
                    case_plans.clear()
                    case_plans.update(prepared)
                except ValueError as exc:
                    st.error(str(exc))
            st.caption('Positive shift subtracts from donor and adds to receiver. Selected dates identify intervals through the next effective date.')
            upload=st.file_uploader('Custom complete schedule CSV (date, stream, tank, fraction)',type=['csv'],key='nu_upload_'+key)
            custom_name=st.text_input('Custom allocation case name',value='Custom allocation',key='nu_custom_'+key)
            if st.button('Add custom allocation case',key='nu_custom_add') and upload is not None:
                try:
                    plan=replace(scenario.history_plan,schedules=parse_schedule(pd.read_csv(upload)))
                    from material_balance_studio.network_history import prepare_history
                    prepare_history(scenario.base_network,plan,scenario.observations)
                    case_plans[custom_name]=plan
                except ValueError as exc:
                    st.error(str(exc))
            st.write('Prepared cases',list(case_plans))
            if case_plans:
                st.dataframe(pd.DataFrame([dict(Case=n,Date=r.date,Stream=r.stream,Tank=t,Fraction=f) for n,p in case_plans.items() for r in p.schedules for t,f in r.fractions]),hide_index=True)
            limit=st.number_input('Allocation-case optimizer evaluation limit',min_value=1,max_value=1000,value=40,key='nu_alloc_budget_'+key)
            st.caption(f'{len(case_plans)} independent bounded fits; forward calls depend on convergence and active parameter count.')
            if st.button('Run allocation comparison',key='nu_allocation_run',disabled=not case_plans):
                try:
                    saved['allocation']=allocation_study(scenario,tuple(case_plans.items()),int(limit))
                except (ValueError,ArithmeticError) as exc:
                    st.error(str(exc))
            study=saved.get('allocation')
            if study:
                st.write('Allocation robustness:',study.robustness)
                for warning in study.warnings:
                    st.caption(warning)
                rows=[]
                for case in study.cases:
                    for name,base,value,change,percent,magnitude,basis in case.drift:
                        rows.append(dict(Case=case.name,Parameter=name,Base=to_display(base,specs[name].unit,units),Fitted=to_display(value,specs[name].unit,units),Change=to_display(change,specs[name].unit,units),Percent=percent,Unit=display_unit(specs[name].unit,units),Basis=basis,Robustness=case.robustness))
                st.dataframe(pd.DataFrame(rows),hide_index=True)
                for name in names:
                    subset=[r for r in rows if r['Parameter']==name]
                    st.plotly_chart(go.Figure(go.Scatter(x=[r['Case'] for r in subset],y=[r['Fitted'] for r in subset],mode='lines+markers')).update_layout(title=name,yaxis_title=display_unit(specs[name].unit,units)),key='nu_drift_'+name)
                summaries=[]
                for c in study.cases:
                    rmse=dict(c.match.final_network_metrics).get('rmse') if c.match else None
                    summaries.append(dict(Case=c.name,Network_RMSE=None if rmse is None else to_display(rmse,'pressure',units),
                        Pressure_unit=display_unit('pressure',units),Tank_metrics_Pa=str(c.match.final_tank_metrics) if c.match else '',
                        Bounds=str(c.match.bound_flags) if c.match else '',Drive_indices=str(c.drives),Cumulative_transfers_m3=str(c.transfers),Failure=c.failure))
                st.dataframe(pd.DataFrame(summaries),hide_index=True)
                chosen_case=st.selectbox('Inspect allocation case',range(len(study.cases)),format_func=lambda i:study.cases[i].name,key='nu_case_'+key)
                case=study.cases[chosen_case]
                if case.match:
                    st.dataframe(parameter_frame(case.match,units),hide_index=True)
                    st.dataframe(pd.DataFrame([dict(Tank=n,**{metric:None if value is None else to_display(value,'pressure',units) for metric,value in metrics}) for n,metrics in case.match.final_tank_metrics]),hide_index=True)
                    st.caption(f"Per-tank pressure metrics in {display_unit('pressure',units)}. Parameter table includes inactive properties; allocation was held fixed throughout this fit.")
                if case.failure:
                    st.warning(case.failure)
    if saved:
        payload={k:asdict(v) for k,v in saved.items() if k!='allocation_plans'}
        st.download_button('Export network diagnostic record',json.dumps(payload,default=str,indent=2),file_name='network_identifiability.json',key='nu_export')
