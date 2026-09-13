"""Historical allocation and matching workspace, separate from forward setup."""
from dataclasses import asdict,replace
from datetime import date
import json
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from material_balance_studio.domain.models import CUMULATIVE_FIELDS
from material_balance_studio.io.history import read_table,history_from_frame,_date
from material_balance_studio.units.conversions import from_si,to_si
from material_balance_studio.units.display import from_display,to_display,display_unit
from material_balance_studio.network_history import NetworkObservation,AllocationRow,HistoryPlan,observations_from_network,prepare_history
from material_balance_studio.network_matching import registry,default_spec,MatchParameter,history_match_network,apply_matched_parameters
from material_balance_studio.network_matching.parameters import validate_specs
from material_balance_studio.network_matching.diagnostics import pressure_diagnostics,data_qc,closure_summary,drive_terms,field_drive,allocation_crosscheck
from material_balance_studio.tank_network import simulate_network,NetworkSettings
from .tank_network import _history_input,tank_frame,connection_frame


def observation_input(observations,units):
    return pd.DataFrame([dict(tank=o.tank,date=str(o.date),pressure=from_si(o.pressure,"pressure",units),source=o.source,
        sigma=None if o.sigma is None else from_si(o.sigma,"pressure",units),qc_flag=o.qc_flag,include_in_match=o.include_in_match,note=o.note) for o in observations],
        columns=["tank","date","pressure","source","sigma","qc_flag","include_in_match","note"])


def parse_observations(frame,units):
    required={"tank","date","pressure"}
    optional={"source","sigma","qc_flag","include_in_match","note"}
    if not required<=set(frame.columns) or set(frame.columns)-required-optional:
        raise ValueError("Pressure table needs tank, date, pressure and optional source, sigma, qc_flag, include_in_match, note.")
    result=[]
    for r in frame.to_dict("records"):
        def text(k): return "" if pd.isna(r.get(k)) else str(r[k])
        included=r.get("include_in_match",True)
        if pd.isna(included): included=True
        if not isinstance(included,bool): raise ValueError("include_in_match must be boolean, not text.")
        result.append(NetworkObservation(str(r["tank"]),_date(r["date"],2),to_si(float(r["pressure"]),"pressure",units),text("source"),
            None if pd.isna(r.get("sigma")) else to_si(float(r["sigma"]),"pressure",units),text("qc_flag"),included,text("note")))
    return tuple(result)


def parse_schedule(frame):
    if set(frame.columns)!={"date","stream","tank","fraction"}:
        raise ValueError("Allocation columns must be date, stream, tank, fraction.")
    if frame.duplicated(["date","stream","tank"]).any():
        raise ValueError("Duplicate tank/date/stream allocation row.")
    rows=[]
    for (day,stream),group in frame.groupby(["date","stream"],sort=True,dropna=False):
        rows.append(AllocationRow(_date(day,2),str(stream),tuple((str(r["tank"]),float(r["fraction"])) for r in group.to_dict("records"))))
    return tuple(rows)


def audit_frame(audit,units):
    rows=[]
    for a in audit:
        q="gas_volume" if a.stream in ("gp","ginj") else "oil_volume" if a.stream=="np" else "water_volume"
        row=dict(Start=a.start,End=a.end,Tank=a.tank,Stream=a.stream,Fraction=a.fraction,Unit=display_unit(q,units))
        row.update({k:to_display(getattr(a,k),q,units) for k in ("field_increment","allocated_increment","tank_cumulative","interval_conservation_error","cumulative_conservation_error")})
        rows.append(row)
    return pd.DataFrame(rows)


def parameter_frame(result,units,current=None):
    targets=registry(result.base_network)
    fitted=dict(result.fitted_parameters)
    bounds={n:(lo,hi) for n,lo,hi in result.bound_flags}
    rows=[]
    specifications={s.name:s for s in result.specifications}
    for name,t in targets.items():
        s=specifications.get(name)
        value=fitted.get(name,t.value)
        initial=s.initial if s is not None and s.active else t.value
        row=dict(Parameter=name,Owner=t.owner,Kind="Connection" if t.edge else "Tank",Unit=display_unit(t.unit,units),
            Initial=to_display(initial,t.unit,units),Lower=to_display(s.lower,t.unit,units) if s else None,Upper=to_display(s.upper,t.unit,units) if s else None,
            Fitted=to_display(value,t.unit,units),Change_percent=None if initial==0 else 100*(value-initial)/initial,
            At_bound=any(bounds.get(name,(False,False))),Active=s.active if s else False)
        if current is not None and name in current:
            row["Current base"]=to_display(current[name].value,t.unit,units)
            row["Will change"]=bool(s and s.active and current[name].value!=value)
        rows.append(row)
    return pd.DataFrame(rows)


def comparison_frame(result,units):
    initial={(r.tank,r.date):r for r in result.initial_pressure_rows}
    rows=[]
    for r in result.final_pressure_rows:
        a=initial[r.tank,r.date]
        rows.append(dict(Tank=r.tank,Date=r.date,Included=r.included,Source=r.source,QC=r.qc_flag,Note=r.note,
            Pobs=to_display(r.observed,"pressure",units),Pcalc_initial=to_display(a.calculated,"pressure",units),Pcalc_matched=to_display(r.calculated,"pressure",units),
            Residual_initial=to_display(a.residual,"pressure",units),Residual_matched=to_display(r.residual,"pressure",units),Sigma=to_display(r.sigma,"pressure",units),Effective_sigma=to_display(r.effective_sigma,"pressure",units),
            Objective_residual_initial=a.objective_residual,Objective_residual_matched=r.objective_residual))
    return pd.DataFrame(rows)


def _metrics_table(before,after,units):
    old,new=dict(before),dict(after)
    return pd.DataFrame([{"Tank":name,**{f"{metric} {label}":to_display(dict(source.get(name,())).get(metric),"pressure",units)
        for metric in ("rmse","mae","bias","maximum") for label,source in (("initial",old),("matched",new))}} for name in sorted(set(old)|set(new))])


def _remember_sigma(widget_key,units):
    st.session_state["nh_default_sigma_si"]=from_display(st.session_state[widget_key],"pressure",units)


def network_history_workflow(units):
    st.subheader("Network History & Matching")
    if not st.checkbox("Configure network historical analysis",key="nh_enable"):
        return
    network=st.session_state.get("nt_network")
    if network is None:
        st.info("Create a network in Multi-Tank Forward first.")
        return
    plan=st.session_state.get("nh_plan",HistoryPlan())
    observations=st.session_state.get("nh_observations",observations_from_network(network))
    rev=st.session_state.get("nh_revision",0)
    prefix=f"nh_{rev}_{units}"
    input_units=st.selectbox("History/pressure input table units (SI pressure is Pa)",["SI","FIELD"],key="nh_input_units")
    key=prefix+input_units
    with st.expander("Direct histories or field allocation",expanded=True):
        st.caption("Direct tank history remains the default. Selected field streams replace their direct counterparts only with the explicit replacement choice. Fractions apply from the effective date until the next row.")
        streams=st.multiselect("Field-allocated streams",list(CUMULATIVE_FIELDS),default=list(plan.allocated_streams),key=key+"streams")
        replacement=st.checkbox("Explicitly replace direct values for selected allocated streams",value=plan.replace_direct,key=key+"replace")
        field_frame=st.data_editor(_history_input(plan.field_history,input_units),num_rows="dynamic",key=key+"field",hide_index=True)
        field_upload=st.file_uploader("Field-total history CSV/XLSX (optional, used on save)",type=["csv","xlsx"],key=key+"field_upload")
        schedules=pd.DataFrame([dict(date=str(s.date),stream=s.stream,tank=n,fraction=v) for s in plan.schedules for n,v in s.fractions],columns=["date","stream","tank","fraction"])
        edited=st.data_editor(schedules,num_rows="dynamic",key=key+"schedule",hide_index=True,column_config={"stream":st.column_config.SelectboxColumn(options=list(CUMULATIVE_FIELDS)),"tank":st.column_config.SelectboxColumn(options=[t.name for t in network.tanks])})
        allocation_upload=st.file_uploader("Allocation CSV/XLSX: date, stream, tank, fraction",type=["csv","xlsx"],key=key+"allocation_upload")
        st.caption("Each stream/date needs every tank, including zero fractions, summing to 1 within 1e-10. Start each selected stream on the network initial date. No normalization or fraction fitting.")
    with st.expander("Tank pressure observations and metadata",expanded=True):
        obs_frame=st.data_editor(observation_input(observations,input_units),num_rows="dynamic",key=key+"observations",hide_index=True)
        obs_upload=st.file_uploader("Tank pressure observations CSV/XLSX",type=["csv","xlsx"],key=key+"obs_upload")
        st.caption("Dates may differ by tank. sigma is observational standard deviation. Missing include flags default to included; all dates enter the timeline, including excluded surveys.")
    unsaved_inputs=(tuple(streams)!=plan.allocated_streams or replacement!=plan.replace_direct or not field_frame.equals(_history_input(plan.field_history,input_units)) or not edited.equals(schedules) or not obs_frame.equals(observation_input(observations,input_units)) or any(u is not None for u in (field_upload,allocation_upload,obs_upload)))
    if unsaved_inputs:
        st.info("Save the edited historical inputs before evaluating or matching the network.")
    if st.button("Save histories, allocation and observations",key="nh_save_inputs"):
        try:
            for upload in (field_upload,allocation_upload,obs_upload):
                if upload is not None: upload.seek(0)
            proposed=HistoryPlan(tuple(streams),history_from_frame(read_table(field_upload) if field_upload else field_frame,input_units) if streams else (),
                parse_schedule(read_table(allocation_upload) if allocation_upload else edited) if streams else (),replacement)
            obs=parse_observations(read_table(obs_upload) if obs_upload else obs_frame,input_units)
            prepare_history(network,proposed,obs)
            st.session_state["nh_plan"],st.session_state["nh_observations"]=proposed,obs
            st.session_state["nh_revision"]=rev+1
            st.rerun()
        except (ValueError,TypeError) as exc:
            st.error(str(exc))
    try:
        prepared=prepare_history(network,plan,observations)
    except ValueError as exc:
        st.error(str(exc))
        return
    if prepared.audit:
        with st.expander("Allocation audit and effective-date plots",expanded=True):
            audit=audit_frame(prepared.audit,units)
            st.dataframe(audit,hide_index=True)
            st.download_button("Download allocation audit",audit.to_csv(index=False),"allocation_audit.csv","text/csv")
            for stream in plan.allocated_streams:
                rows=sorted([s for s in plan.schedules if s.stream==stream],key=lambda s:s.date)
                end=max(h.date for n in prepared.network.tanks for h in n.history)
                fig=go.Figure()
                for node in network.tanks:
                    fig.add_trace(go.Scatter(x=[s.date for s in rows]+[end],y=[dict(s.fractions)[node.name] for s in rows]+[dict(rows[-1].fractions)[node.name]],name=node.name,line_shape="hv"))
                fig.update_layout(title=stream+" effective allocation",yaxis_title="Fraction",yaxis_range=[0,1])
                st.plotly_chart(fig,key="nh_allocation_"+stream)
    mode=st.selectbox("Network pressure objective",["Unweighted","Pressure uncertainty weighted"],key="nh_mode")
    use_default=st.checkbox("Explicitly supply default sigma for missing uncertainties",key="nh_default_enabled")
    sigma=from_display(st.number_input(f"Default sigma [{display_unit('pressure',units)}]",min_value=1e-8,value=float(to_display(st.session_state.get("nh_default_sigma_si",1e5),"pressure",units)),key=prefix+"sigma",on_change=_remember_sigma,args=(prefix+"sigma",units)),"pressure",units) if use_default else None
    step=st.number_input("Network matching maximum internal step (days; 0 uses Phase 5A control)",min_value=0.,value=0.,key="nh_step")
    settings=NetworkSettings(max_step_days=step or None)
    st.caption("Objective: Σ[(Pcalc−Pobs)/1 MPa]² for unweighted mode; Σ[(Pcalc−Pobs)/sigma]² for weighted mode. Pressure RMSE is reported separately. Closure tolerances are unchanged.")
    if st.button("Evaluate initial network pressure match",key="nh_evaluate",disabled=unsaved_inputs):
        try:
            with st.spinner("Evaluating the unchanged network simulator…"):
                sim=simulate_network(prepared.network,settings)
                rows,tank_metrics,aggregate=pressure_diagnostics(sim,observations,mode,sigma)
                st.session_state["nh_initial"]=(prepared.fingerprint,mode,sigma,settings,sim,rows,tank_metrics,aggregate)
        except ValueError as exc:
            st.error(str(exc))
    prior=st.session_state.get("nh_initial")
    if prior and prior[:4]==(prepared.fingerprint,mode,sigma,settings):
        with st.expander("Initial per-tank pressure metrics",expanded=True):
            st.dataframe(_metrics_table(prior[6],prior[6],units),hide_index=True)
            st.write("Network metrics",{k:to_display(v,"pressure",units) for k,v in prior[7]})
            st.write("Closure",closure_summary(prior[4]))
            for warning in allocation_crosscheck(plan,prior[5],network): st.warning(warning)
    targets=registry(network)
    active=st.multiselect("Active network match parameters",[n for n,t in targets.items() if t.enabled],key="nh_active")
    stored=st.session_state.get("nh_specs",{})
    specs=[stored.get(n,default_spec(n,targets[n])) for n in active]
    frame=pd.DataFrame([dict(Parameter=s.name,Owner=targets[s.name].owner,Unit=display_unit(s.unit,units),Initial=to_display(s.initial,s.unit,units),Lower=to_display(s.lower,s.unit,units),Upper=to_display(s.upper,s.unit,units),Scale=to_display(s.scale,s.unit,units),Transformation=s.transformation) for s in specs],columns=["Parameter","Owner","Unit","Initial","Lower","Upper","Scale","Transformation"])
    edited_specs=st.data_editor(frame,hide_index=True,key=prefix+str(tuple(active))+"specs",disabled=["Parameter","Owner","Unit"],column_config={"Transformation":st.column_config.SelectboxColumn(options=["linear","log"],required=True)})
    st.caption("Linear T can include exactly zero. Log T requires an explicitly positive lower bound. Bounds never expand automatically. Run uses saved bounds; save edited bounds before running or switching units.")
    if st.button("Save network matching bounds",key="nh_save_bounds"):
        try:
            built=[]
            for r in edited_specs.to_dict("records"):
                target=targets[r["Parameter"]]
                vals=[from_display(float(r[k]),target.unit,units) for k in ("Initial","Lower","Upper","Scale")]
                built.append(MatchParameter(r["Parameter"],target.description,*vals[:3],target.unit,vals[3],r["Transformation"]))
            validate_specs(network,built)
            st.session_state["nh_specs"]={**stored,**{s.name:s for s in built}}
            st.session_state["nh_revision"]=rev+1
            st.rerun()
        except (ValueError,TypeError) as exc:
            st.error(str(exc))
    if active:
        try:
            for warning in data_qc(network,specs,observations): st.warning(warning)
        except ValueError as exc:
            st.error(str(exc))
    budget=st.slider("Network optimizer maximum evaluations",10,200,60,key="nh_budget")
    name=st.text_input("Network match scenario name",value="Network match",key="nh_name")
    if st.button("Run bounded network history match",key="nh_run",disabled=not active or unsaved_inputs or not edited_specs.equals(frame)):
        try:
            with st.spinner("Matching network parameters using fresh coupled simulations…"):
                result=history_match_network(network,specs,observations,plan,mode,sigma,settings,budget,name)
            st.session_state.setdefault("network_match_scenarios",[]).append(result)
        except ValueError as exc:
            st.error(str(exc))
    scenarios=st.session_state.get("network_match_scenarios",[])
    if scenarios:
        index=st.selectbox("Saved network match scenario",range(len(scenarios)),index=len(scenarios)-1,format_func=lambda i:f"{i+1} · {scenarios[i].scenario_name}",key="nh_scenario")
        show_match(scenarios[index],network,units)


def show_match(result,current,units):
    st.subheader(result.scenario_name)
    st.write(f"{result.qc_status} · optimizer status {result.optimizer_status} · {result.function_evaluations} candidate evaluations · {result.failed_candidates} failed candidates")
    st.caption(f"{result.mode}; history/allocation fingerprint {result.history_fingerprint[:16]}; result uses its saved input snapshot.")
    for warning in result.warnings: st.warning(warning)
    st.write("Objective before / after",result.initial_objective,result.final_objective)
    st.dataframe(_metrics_table(result.initial_tank_metrics,result.final_tank_metrics,units),hide_index=True)
    st.write("Network pressure metrics",{label:{k:to_display(v,"pressure",units) for k,v in items} for label,items in (("initial",result.initial_network_metrics),("matched",result.final_network_metrics))})
    parameters=parameter_frame(result,units,registry(current))
    st.write("Tank and aquifer parameters")
    st.dataframe(parameters[parameters.Kind=="Tank"],hide_index=True)
    st.write("Connection transmissibilities")
    st.dataframe(parameters[parameters.Kind=="Connection"],hide_index=True)
    st.caption("Explicit application changes only the active parameters shown above. Current histories, allocation schedules and pressure metadata remain unchanged.")
    if st.button("Apply matched network parameters",key="nh_apply",disabled=not result.converged):
        try:
            st.session_state["nt_network"]=apply_matched_parameters(current,result)
            st.session_state["nt_revision"]=st.session_state.get("nt_revision",0)+1
            st.rerun()
        except ValueError as exc: st.error(str(exc))
    comparison=comparison_frame(result,units)
    st.caption(f"Pressure values, errors and sigma columns use {display_unit('pressure',units)}. Effective sigma includes an explicitly configured default when applicable; objective residuals retain the scenario's matching definition.")
    names=[t.name for t in result.base_network.tanks]
    name=st.selectbox("Network pressure plot tank",names,key="nh_plot_tank")
    subset=comparison[comparison.Tank==name]
    fig=go.Figure()
    for col,label in (("Pobs","Observed"),("Pcalc_initial","Initial"),("Pcalc_matched","Matched")):
        fig.add_trace(go.Scatter(x=subset.Date,y=subset[col],mode="markers" if col=="Pobs" else "lines+markers",name=label))
    fig.update_layout(title=name+" pressure match",yaxis_title=display_unit("pressure",units))
    st.plotly_chart(fig,key="nh_pressure")
    fig=go.Figure()
    for column in ("Residual_initial","Residual_matched"):
        fig.add_trace(go.Scatter(x=subset.Date,y=subset[column],name=column,mode="lines+markers"))
    fig.update_layout(title=name+" pressure residuals",yaxis_title=display_unit("pressure",units))
    st.plotly_chart(fig,key="nh_residuals")
    st.dataframe(comparison,hide_index=True)
    with st.expander("History-match observation inspector"):
        if not comparison.empty:
            row=st.selectbox("Observed pressure row",range(len(comparison)),format_func=lambda i:f"{comparison.iloc[i].Tank} · {comparison.iloc[i].Date}",key="nh_obs_inspect")
            st.dataframe(comparison.iloc[[row]],hide_index=True)
    if result.final_simulation is not None:
        final=result.final_simulation
        with st.expander("Matched transfer histories"):
            st.dataframe(connection_frame(final,units),hide_index=True)
            for field,q in (("pressure_difference","pressure"),("rate","aquifer_rate"),("cumulative_transfer","reservoir_volume")):
                fig=go.Figure()
                for j,edge in enumerate(final.network.connections):
                    for label,sim in (("initial",result.initial_simulation),("matched",final)):
                        if sim:
                            fig.add_trace(go.Scatter(x=[s.time for s in sim.states],y=[to_display(getattr(s.connections[j],field),q,units) for s in sim.states],name=f"{edge.from_tank} → {edge.to_tank} {label}"))
                fig.update_layout(title=field,yaxis_title=display_unit(q,units))
                st.plotly_chart(fig,key="nh_transfer_"+field)
        with st.expander("Signed drive indices and equation inspector",expanded=True):
            j=names.index(name)
            fig=go.Figure()
            for index in ("DDI","SDI","CDI","WDI","IDI","TDI"):
                fig.add_trace(go.Scatter(x=[s.time for s in final.states],y=[drive_terms(s.tanks[j])["values"][index] for s in final.states],name=index))
            fig.update_layout(title=name+" signed drive contributions",yaxis_title="Contribution / produced voidage")
            st.plotly_chart(fig,key="nh_drives")
            st.caption("TDI=X/Fprod is positive for received support and negative for export. Undefined at zero production. No normalization or 100% stacking. Field net transfer is internal redistribution, not added field support.")
            if final.states:
                ix=st.selectbox("Matched network equation date",range(len(final.states)),format_func=lambda i:str(final.states[i].time),key="nh_equation_date")
                state=final.states[ix]
                st.dataframe(tank_frame(replace(final,states=(state,)),units),hide_index=True)
                st.dataframe(pd.DataFrame([dict(Tank=t.name,**drive_terms(t)["values"],Drive_sum=drive_terms(t)["total"],Closure_error=drive_terms(t)["closure_error"]) for t in state.tanks]),hide_index=True)
                st.write("Field drive values",field_drive(state)["values"])
                st.write("Field net transfer",to_display(field_drive(state)["supports"]["TDI"],"reservoir_volume",units),display_unit("reservoir_volume",units))
                st.dataframe(pd.DataFrame([{"Tank":t.name,"Neighbor":n,"Signed contribution":to_display(v,"reservoir_volume",units)} for t in state.tanks for n,v in t.connection_contributions]),hide_index=True)
    with st.expander("Closure margin, failure audit and scenario export"):
        st.write(dict(maximum_valid_relative_residual=result.maximum_valid_relative_residual,maximum_valid_absolute_residual=result.maximum_valid_absolute_residual,
            minimum_closure_margin=result.minimum_closure_margin,maximum_transfer_error=result.maximum_transfer_conservation_error))
        st.dataframe(pd.DataFrame([dict(Objective=e.objective,Valid=e.valid,Failure_type=e.failure_type,Failure=e.failure,Affected=str(e.affected),Parameters=str(e.parameters)) for e in result.evaluations]),hide_index=True)
        st.download_button("Download network match scenario",json.dumps(asdict(result),default=str,indent=2,allow_nan=False),"network_match.json","application/json")
