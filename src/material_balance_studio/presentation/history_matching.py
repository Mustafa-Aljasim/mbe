"""Scenario configuration and explicit acceptance; optimization stays in matching."""
from dataclasses import asdict,replace
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from material_balance_studio.matching import MatchParameter,observations_from_history,history_match
from material_balance_studio.matching.parameters import parameter_registry,parameter_value
from material_balance_studio.units.display import to_display,from_display,column_label,display_unit
from material_balance_studio.aquifer.registry import MODELS
from .diagnosis import diagnosis_frame


def refresh_matching_units():
    for suffix,draft in st.session_state.get("hm_parameter_drafts",{}).items():
        for field in ("initial","lower","upper"):
            st.session_state[f"hm_{field}_{suffix}"]=to_display(draft[field],draft["quantity"],st.session_state["display_units"])
    if "hm_default_sigma_si" in st.session_state:
        st.session_state["hm_default_sigma"]=to_display(st.session_state["hm_default_sigma_si"],"pressure",st.session_state["display_units"])


def save_matching_value(suffix,field):
    draft=st.session_state["hm_parameter_drafts"][suffix]
    draft[field]=from_display(st.session_state[f"hm_{field}_{suffix}"],draft["quantity"],st.session_state["display_units"])


def save_matching_sigma():
    st.session_state["hm_default_sigma_si"]=from_display(st.session_state["hm_default_sigma"],"pressure",st.session_state["display_units"])


def fitted_parameter_frame(result,units,current_values=None):
    fitted=dict(result.fitted_parameters)
    flags={n:(lo,hi) for n,lo,hi in result.bound_flags}
    rows=[]
    for s in result.specifications:
        value=fitted.get(s.name,parameter_value(result.base_tank,s.name))
        row={"Parameter":s.name,"Active":s.active,"Unit":display_unit(s.unit,units),
             "Initial":to_display(s.initial,s.unit,units),"Lower":to_display(s.lower,s.unit,units),
             "Upper":to_display(s.upper,s.unit,units),"Fitted":to_display(value,s.unit,units),
             "Change %":100*(value-s.initial)/s.initial if s.initial else None,
             "At lower bound":flags.get(s.name,(False,False))[0],"At upper bound":flags.get(s.name,(False,False))[1]}
        if current_values is not None:
            row["Current base value"]=to_display(current_values[s.name],s.unit,units)
        rows.append(row)
    return pd.DataFrame(rows)


def apply_matched_to_setup(result):
    """Called only by the explicit apply button, before Streamlit widgets rerender."""
    if not result.converged or result.fitted_tank is None:
        raise ValueError("Cannot apply an unsuccessful match.")
    selected=st.session_state.get("aquifer_type","None")
    if MODELS[selected].key!=result.aquifer_model:
        raise ValueError("Select the scenario's aquifer model before applying its parameters.")
    units=st.session_state["display_units"]
    quantities=parameter_registry(result.base_tank)
    for name,value in result.fitted_parameters:
        if name=="oil_in_place":
            st.session_state["setup_si"]["setup_oil"]=value
            st.session_state["setup_oil"]=to_display(value,"oil_volume",units)
        elif name=="m":
            st.session_state["setup_m"]=value
        else:
            parameter=name.split(".")[1]
            st.session_state["aquifer_draft"][parameter]=value
            st.session_state["aq_"+parameter]=to_display(value,quantities[name][1],units)
    st.session_state["hm_applied"]=True


def history_matching_workflow(units):
    st.subheader("Bounded Reservoir History Matching")
    saved=st.session_state.get("simulation")
    if not saved or "history" not in saved[1]:
        st.info("Run the base simulation first to save the reservoir/PVT/history snapshot used for matching.")
        return
    _,snapshot=saved
    tank,history=snapshot["tank"],snapshot["history"]
    st.caption("Matching uses the last simulation's reservoir, aquifer, PVT and history snapshot. Rerun the base simulation after editing inputs. Only selected parameters can change. PVT, initial pressure, compressibilities and histories remain fixed.")
    if "FAIL" in set(snapshot["history_qc"].Status):
        st.error("History QC contains FAIL flags; resolve these before matching.")
        return
    # Keep normal simulation controls uncluttered and avoid running fits on rerender.
    if not st.checkbox("Configure history match",key="hm_configure"):
        return
    registry=parameter_registry(tank)
    if any(name not in registry for name in st.session_state.get("hm_active",[])):
        st.session_state["hm_active"]=[name for name in st.session_state["hm_active"] if name in registry]
    active=st.multiselect("Parameters allowed to change",list(registry),default=[],key="hm_active")
    specs=[]
    for name in active:
        description,quantity=registry[name]
        initial=parameter_value(tank,name)
        lower=initial*.25 if initial>0 else 0.
        upper=initial*4 if initial>0 else 1.
        if name=="aquifer.encroachment_angle":
            upper=min(360.,upper)
        if name=="aquifer.porosity":
            upper=min(.99,upper)
        with st.expander(description+" — physical bounds",expanded=True):
            suffix=f"{tank.aquifer.key if tank.aquifer else 'none'}_{name}_{initial}"
            draft=st.session_state.setdefault("hm_parameter_drafts",{}).setdefault(suffix,dict(initial=initial,lower=lower,upper=upper,quantity=quantity))
            for field in ("initial","lower","upper"):
                st.session_state.setdefault(f"hm_{field}_{suffix}",to_display(draft[field],quantity,units))
            start=st.number_input(column_label("Initial value",quantity,units),format="%.8g",key="hm_initial_"+suffix,on_change=save_matching_value,args=(suffix,"initial"))
            lo=st.number_input(column_label("Lower bound",quantity,units),format="%.8g",key="hm_lower_"+suffix,on_change=save_matching_value,args=(suffix,"lower"))
            hi=st.number_input(column_label("Upper bound",quantity,units),format="%.8g",key="hm_upper_"+suffix,on_change=save_matching_value,args=(suffix,"upper"))
            transform=st.selectbox("Parameter transformation",["linear","log"],key="hm_transform_"+name)
            try:
                specs.append(MatchParameter(name,description,from_display(start,quantity,units),from_display(lo,quantity,units),from_display(hi,quantity,units),quantity,max(abs(initial),.1) if name=="m" else initial,transform))
            except ValueError as exc:
                st.error(str(exc))
    observations=observations_from_history(history,snapshot["pressure_metadata"])
    selected=[]
    with st.expander("Pressure observation inclusion and quality",expanded=True):
        st.caption("All measured observations are retained. Choose exclusions deliberately. Survey values, sigma, source and notes come from the uploaded history.")
        for o in observations:
            include=st.checkbox(str(o.date),value=True,key=f"hm_include_{o.date}")
            selected.append(replace(o,include_in_match=include))
        table=pd.DataFrame([dict(Date=str(o.date),Pressure=to_display(o.pressure,"pressure",units),Sigma=to_display(o.sigma,"pressure",units),Source=o.source,QC=o.qc_flag,Included=o.include_in_match,Note=o.note) for o in selected])
        st.caption("Pressure and sigma: "+display_unit("pressure",units))
        st.dataframe(table,hide_index=True,width="stretch")
        st.dataframe(snapshot["history_qc"],hide_index=True,width="stretch")
    mode=st.radio("Pressure residual weighting",["Unweighted","Pressure uncertainty weighted"],key="hm_weighting")
    default=None
    if mode=="Pressure uncertainty weighted" and st.checkbox("Use an explicit default sigma for missing/invalid uncertainties",key="hm_use_default_sigma"):
        st.session_state.setdefault("hm_default_sigma_si",1e5)
        st.session_state.setdefault("hm_default_sigma",to_display(st.session_state["hm_default_sigma_si"],"pressure",units))
        default=from_display(st.number_input(column_label("Default sigma","pressure",units),key="hm_default_sigma",on_change=save_matching_sigma),"pressure",units)
    st.caption("Unweighted objective = Σ[(Pcalc−Pobs)/1 MPa]². Weighted objective = Σ[(Pcalc−Pobs)/sigma]². MBE residual is a separate strict forward-simulation validity check.")
    evaluations=st.number_input("Maximum optimizer function evaluations",min_value=1,max_value=1000,value=100,step=1,key="hm_max_nfev")
    count=sum(o.include_in_match for o in selected)
    if count<=len(specs):
        st.error("FAIL: observations must outnumber active parameters.")
    elif count<max(5,2*len(specs)+1):
        st.warning("Weakly overdetermined: limited pressure data for the selected parameters.")
    if st.button("Run bounded history match",key="run_history_match",disabled=not active or len(specs)!=len(active) or count<=len(specs)):
        try:
            with st.spinner("Running bounded forward simulations…"):
                result=history_match(tank,history,specs,selected,mode,default,int(evaluations))
            scenarios=st.session_state.setdefault("history_match_scenarios",[])
            scenarios.append(result)
            st.session_state["hm_scenario_index"]=len(scenarios)-1
        except ValueError as exc:
            st.error(str(exc))
    scenarios=st.session_state.get("history_match_scenarios",[])
    if not scenarios:
        return
    st.write("In-session match scenarios")
    st.dataframe(pd.DataFrame([dict(Scenario=i+1,Aquifer=r.aquifer_model,Parameters=", ".join(n for n,_ in r.fitted_parameters),Weighting=r.weighting_mode,RMSE=to_display(r.final_metrics["rmse"],"pressure",units),Objective=r.final_objective,QC=r.qc_status) for i,r in enumerate(scenarios)]),hide_index=True)
    st.caption("RMSE units: "+display_unit("pressure",units)+". Objectives with different weighting/observations are not directly comparable.")
    index=st.selectbox("Match scenario",range(len(scenarios)),format_func=lambda i:f"Scenario {i+1}",key="hm_scenario_index")
    result=scenarios[index]
    st.write(f"History match: {'CONVERGED' if result.converged else 'NOT CONVERGED'} · QC {result.qc_status}")
    st.write(f"Observations used: {sum(o.include_in_match for o in result.observations)} / {len(result.observations)}. Aquifer: {result.aquifer_model}.")
    st.write(result.termination_message)
    st.write(dict(Initial_objective=result.initial_objective,Final_objective=result.final_objective,
                  Forward_evaluations=result.function_evaluations,Optimizer_nfev=result.optimizer_nfev))
    st.metric("Maximum forward MBE relative residual (valid candidates)",f"{result.maximum_mbe_relative:.3e}")
    st.caption("A low pressure mismatch never relaxes MBE closure. Fitted parameters are not established as unique.")
    for warning in result.warnings:
        st.warning(warning)
    st.dataframe(fitted_parameter_frame(result,units),hide_index=True,width="stretch")
    st.dataframe(pd.DataFrame([{"Stage":"Initial",**{k:to_display(v,"pressure",units) for k,v in result.initial_metrics.items()}},{"Stage":"Matched",**{k:to_display(v,"pressure",units) for k,v in result.final_metrics.items()}}]),hide_index=True)
    st.caption("Pressure metric units: "+display_unit("pressure",units))
    for residual_plot in (False,True):
        fig=go.Figure()
        if not residual_plot:
            for label,simulation in (("Initial",result.initial_simulation),("Matched",result.final_simulation)):
                if simulation:
                    fig.add_scatter(x=[s.date for s in simulation.states],y=[to_display(s.pressure,"pressure",units) for s in simulation.states],name=label,mode="lines")
            for included in (True,False):
                rows=[r for r in result.pressure_series if r["included"]==included]
                fig.add_scatter(x=[r["date"] for r in rows],y=[to_display(r["observed"],"pressure",units) for r in rows],name="Observed included" if included else "Observed excluded",mode="markers",marker_symbol="circle" if included else "x")
        else:
            for key in ("initial_residual","residual"):
                fig.add_scatter(x=[r["date"] for r in result.pressure_series],y=[to_display(r[key],"pressure",units) for r in result.pressure_series],name=key,mode="lines+markers")
            fig.add_hline(y=0,line_dash="dash")
        fig.update_layout(xaxis_title="Date",yaxis_title=column_label("Pressure residual" if residual_plot else "Pressure","pressure",units))
        st.plotly_chart(fig,key=f"hm_pressure_{residual_plot}")
    with st.expander("Evaluation audit and pressure residuals"):
        st.caption("Every objective call is recorded, including finite differences and failed trials. These are not accepted optimizer iterations and need not improve monotonically.")
        audit=pd.DataFrame(asdict(e) for e in result.evaluations)
        st.dataframe(audit,hide_index=True)
        residuals=pd.DataFrame(result.pressure_series)
        for key in ("observed","initial","matched","initial_residual","residual","sigma"):
            residuals[key]=residuals[key].map(lambda v:to_display(v,"pressure",units))
            residuals.rename(columns={key:column_label(key,"pressure",units)},inplace=True)
        st.dataframe(residuals,hide_index=True)
        st.download_button("Download match residuals",residuals.to_csv(index=False),f"match_residuals_{index+1}.csv")
        st.write("PVT fingerprint: "+result.pvt_fingerprint)
    if result.diagnostics:
        with st.expander("Matched VRR, drive indices and H–O cross-check",expanded=True):
            st.dataframe(diagnosis_frame(result.diagnostics,units),hide_index=True,width="stretch")
            st.write({"Matched OOIP":to_display(result.fitted_tank.oil_in_place,"oil_volume",units),"H–O diagnostic OOIP":to_display(result.ho_ooip,"oil_volume",units),"Difference %":result.ho_difference_percent})
            st.caption("OOIP units: "+display_unit("oil_volume",units)+". H–O uses all available survey points and the fitted aquifer; excluded match observations remain available for this independent diagnostic review.")
            for note in result.diagnostics["ho"]["fit"]["notes"]:
                st.warning(note)
    st.write("Apply scenario — review proposed base-model changes")
    current={s.name:st.session_state["setup_si"]["setup_oil"] if s.name=="oil_in_place" else st.session_state.get("setup_m",0.) if s.name=="m" else st.session_state["aquifer_draft"].get(s.name.split(".")[1]) for s in result.specifications}
    st.dataframe(fitted_parameter_frame(result,units,current),hide_index=True,width="stretch")
    same_model=MODELS[st.session_state.get("aquifer_type","None")].key==result.aquifer_model
    if not same_model:
        st.warning("Select this scenario's aquifer type in Reservoir Setup before applying.")
    st.button("Apply matched parameters to reservoir model",key="apply_history_match",on_click=apply_matched_to_setup,args=(result,),disabled=not result.converged or not same_model)
    if st.session_state.get("hm_applied"):
        st.info("Matched parameters were explicitly applied to setup. Run simulation again to create new base results.")
