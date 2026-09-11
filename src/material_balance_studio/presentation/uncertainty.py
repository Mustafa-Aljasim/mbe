"""Scenario diagnostics; calculations and stored results remain in canonical SI."""
from dataclasses import asdict,replace
import json
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from material_balance_studio.uncertainty import analyze,parameter_profile,objective_surface
from material_balance_studio.uncertainty.result import scenario_identity
from material_balance_studio.uncertainty.jacobian import matrix_diagnostics
from material_balance_studio.uncertainty.identifiability import classify
from material_balance_studio.units.display import to_display,display_unit


def sensitivity_frame(result,scenario,units):
    specs={s.name:s for s in scenario.specifications}
    return pd.DataFrame({f"{name} [{display_unit('pressure',units)}/{display_unit(specs[name].unit,units)}]":
        np.asarray(result.sensitivity)[:,j]*to_display(1.,"pressure",units)/to_display(1.,specs[name].unit,units)
        for j,name in enumerate(result.parameters)},index=result.dates)


def point_frame(points,scenario,units):
    specs={s.name:s for s in scenario.specifications}
    rows=[]
    for p in points:
        row={f"{n} [{display_unit(specs[n].unit,units)}]":to_display(v,specs[n].unit,units) for n,v in p.parameters}
        row.update(Objective=p.objective,RMSE=to_display(p.rmse,"pressure",units) if p.rmse is not None else None,
                   Valid=p.valid,Status=p.optimizer_status,Bounds=str(p.bound_flags),Failure=p.failure,Evaluations=p.evaluations)
        rows.append(row)
    return pd.DataFrame(rows)


def uncertainty_workflow(units):
    st.subheader("Identifiability & Uncertainty")
    st.info("A low pressure RMSE does not imply unique reservoir parameters. These diagnostics assess information in the selected history; they do not change the reservoir or accepted match.")
    scenarios=st.session_state.get("history_match_scenarios",[])
    if not scenarios:
        st.caption("Create a scenario in History Matching first.")
        return
    index=st.selectbox("Diagnostic match scenario",range(len(scenarios)),format_func=lambda i:f"Scenario {i+1} · {scenarios[i].aquifer_model} · {scenarios[i].weighting_mode}",key="unc_scenario")
    scenario=scenarios[index]
    if not scenario.converged:
        st.warning("Select a converged scenario with a valid final simulation.")
        return
    identity=scenario_identity(scenario)
    key=identity[:16]
    storage=st.session_state.setdefault("identifiability_results",{})
    assumptions=st.checkbox("Accept statistical assumptions for approximate confidence",key="unc_assumptions_"+key,
        help="Independent Gaussian errors; known observational sigma for weighted fits, otherwise common variance estimated from residuals; adequate observations and local linearity. These assumptions are not established automatically.")
    if st.button("Run local identifiability",key="unc_run"):
        try:
            with st.spinner("Evaluating bounded pressure sensitivities…"):
                storage[identity]=analyze(scenario,statistical_assumptions=assumptions)
        except ValueError as exc:
            st.error(str(exc))
    result=storage.get(identity)
    if result is None:
        return
    st.caption(f"Stored calculation: {result.match_mode}; objective {result.objective_definition}; statistical assumptions {'accepted' if result.statistical_assumptions else 'not accepted'}. Rerun local diagnostics to change this setting.")
    st.metric("History-match pressure RMSE",f"{to_display(scenario.final_metrics['rmse'],'pressure',units):.5g} {display_unit('pressure',units)}")
    specs=[s for s in scenario.specifications if s.active]
    fitted=dict(scenario.fitted_parameters)
    st.dataframe(pd.DataFrame([dict(Parameter=n,Estimate=to_display(fitted[n],next(s.unit for s in specs if s.name==n),units),
        Unit=display_unit(next(s.unit for s in specs if s.name==n),units),Status=status,Evidence=note) for n,status,note in result.statuses]),hide_index=True)
    for warning in result.warnings:
        st.warning(warning)
    with st.expander("Sensitivity and informative observation dates",expanded=True):
        scaled=np.asarray(result.scaled_sensitivity)
        ranking=pd.DataFrame({"Parameter":result.parameters,"RMS dimensionless sensitivity":np.sqrt(np.mean(scaled**2,axis=0))}).sort_values("RMS dimensionless sensitivity",ascending=False)
        st.dataframe(ranking,hide_index=True)
        fig=go.Figure(go.Bar(x=ranking["Parameter"],y=ranking["RMS dimensionless sensitivity"]))
        fig.update_layout(yaxis_title="RMS (dP/dx × parameter scale / Pi)")
        st.plotly_chart(fig,key="unc_rank")
        selected=st.selectbox("Sensitivity parameter",result.parameters,key="unc_sens_"+key)
        raw=sensitivity_frame(result,scenario,units)
        j=result.parameters.index(selected)
        fig=go.Figure(go.Scatter(x=result.dates,y=raw.iloc[:,j],mode="lines+markers"))
        fig.update_layout(xaxis_title="Included pressure observation date",yaxis_title=raw.columns[j])
        st.plotly_chart(fig,key="unc_time")
        st.caption("Near-zero sensitivity means that survey carries little local information about this parameter. This does not predict the benefit of future surveys.")
        st.plotly_chart(go.Figure(go.Heatmap(x=result.parameters,y=result.dates,z=scaled,colorscale="RdBu",zmid=0)),key="unc_sens_heat")
        st.dataframe(raw)
    with st.expander("Parameter correlation",expanded=True):
        if result.correlation is None:
            st.warning("Parameter correlation is unavailable: the scaled Jacobian is rank deficient. A pseudo-inverse would conceal unconstrained directions.")
        else:
            corr=np.asarray(result.correlation)
            st.plotly_chart(go.Figure(go.Heatmap(z=corr,x=result.parameters,y=result.parameters,zmin=-1,zmax=1,colorscale="RdBu")),key="unc_corr")
            st.dataframe(pd.DataFrame([{"Parameter A":a,"Parameter B":b,"Correlation":corr[i,j],"QC":
                "STRONG COUPLING" if abs(corr[i,j])>.9 else "CAUTION" if abs(corr[i,j])>=.7 else "Generally acceptable"}
                for i,a in enumerate(result.parameters) for j,b in enumerate(result.parameters) if j>i]),hide_index=True)
            st.caption("Strong absolute correlation means parameter combinations can compensate in the local match. It is not physical causation. Thresholds 0.7 and 0.9 are engineering diagnostics, not universal statistical laws.")
    with st.expander("Approximate confidence"):
        st.write("LOCAL / LINEARIZED APPROXIMATION — approximate marginal 95% intervals only when the recorded assumptions and numerical checks permit them.")
        st.dataframe(pd.DataFrame([{"Parameter":s.name,"Unit":display_unit(s.unit,units),
            "Standard error":None if result.standard_errors[j] is None else to_display(result.standard_errors[j],s.unit,units),
            "Lower":None if result.intervals[j] is None else to_display(result.intervals[j][0],s.unit,units),
            "Upper":None if result.intervals[j] is None else to_display(result.intervals[j][1],s.unit,units)} for j,s in enumerate(specs)]),hide_index=True)
        st.caption("Missing values mean withheld, not zero uncertainty. Bound-limited intervals are not clipped into artificial certainty. Use profiles to inspect compensation and the permitted direction.")
    with st.expander("Parameter profiles"):
        parameter=st.selectbox("Profile parameter",result.parameters,key="unc_profile_parameter_"+key)
        reopt=st.checkbox("Re-optimize other active parameters",value=True,key="unc_reopt_"+key)
        count=st.slider("Profile grid points",3,21,7,key="unc_profile_count_"+key)
        budget=st.slider("Optimizer evaluations per profile node",10,100,40,key="unc_budget_"+key)
        custom=st.checkbox("Use a deterministic plausible objective increment",key="unc_custom_"+key)
        delta=st.number_input("Plausible objective increment (same objective units)",min_value=1e-8,value=1.,format="%.6g",key="unc_delta_"+key) if custom else None
        st.caption(f"At most {count+1} grid nodes (including the matched coordinate). Approximately up to {(count+1)*(budget*(1+2*(len(specs)-1))+2) if reopt and len(specs)>1 else count+1} forward evaluations. Serial, fresh simulations.")
        if st.button("Run parameter profile",key="unc_profile_run"):
            try:
                with st.spinner("Evaluating profile…"):
                    profile=parameter_profile(scenario,parameter,count,reopt,budget,result.statistical_assumptions,delta)
                profiles=tuple(p for p in result.profiles if not (p.parameter==parameter and p.mode==profile.mode))+(profile,)
                diag=matrix_diagnostics(result.scaled_jacobian)
                pressures=[o.pressure for o in scenario.observations if o.include_in_match]
                statuses=classify(specs,[fitted[s.name] for s in specs],result.scaled_sensitivity,diag,len(pressures),np.ptp(pressures)/scenario.base_tank.initial_pressure,profiles)
                if max(result.derivative_errors)>.05:
                    statuses=tuple((n,"POORLY CONSTRAINED" if status in ("WELL CONSTRAINED","MODERATELY CONSTRAINED") else status,
                        note+" Finite-difference instability limits local inference.") for n,status,note in statuses)
                result=replace(result,profiles=profiles,statuses=statuses)
                storage[identity]=result
                st.rerun()
            except ValueError as exc:
                st.error(str(exc))
        for i,profile in enumerate(result.profiles):
            spec=next(s for s in specs if s.name==profile.parameter)
            st.write(f"{profile.parameter} — {profile.mode}")
            fig=go.Figure(go.Scatter(x=[to_display(dict(p.parameters)[spec.name],spec.unit,units) for p in profile.points],
                y=[p.objective if p.valid and p.optimizer_status!='NOT CONVERGED' else None for p in profile.points],mode="lines+markers",connectgaps=False))
            if profile.threshold is not None:
                fig.add_hline(y=profile.threshold,line_dash="dash",annotation_text="Selected threshold")
            fig.update_layout(xaxis_title=f"{spec.name} [{display_unit(spec.unit,units)}]",yaxis_title=result.objective_definition)
            st.plotly_chart(fig,key=f"unc_profile_plot_{i}")
            st.caption(profile.threshold_method)
            st.write("Grid-resolved plausible regions:",[(to_display(a,spec.unit,units),to_display(b,spec.unit,units)) for a,b in profile.plausible_ranges])
            st.write("Threshold-crossing brackets (endpoints unresolved within each bracket):",[(to_display(a,spec.unit,units),to_display(b,spec.unit,units)) for a,b in profile.threshold_brackets])
            for warning in profile.warnings:
                st.warning(warning)
            st.dataframe(point_frame(profile.points,scenario,units),hide_index=True)
    with st.expander("Objective surfaces · N / aquifer ambiguity"):
        if len(specs)<2:
            st.caption("A surface requires two active matched parameters.")
        else:
            aquifers=[s.name for s in specs if s.name.startswith("aquifer.")]
            classic=st.checkbox("Classic OOIP N versus aquifer",value=bool(aquifers and "oil_in_place" in result.parameters),disabled=not(aquifers and "oil_in_place" in result.parameters),key="unc_classic_"+key)
            x="oil_in_place" if classic else st.selectbox("Surface parameter 1",result.parameters,key="unc_x_"+key)
            y=st.selectbox("Aquifer strength parameter" if classic else "Surface parameter 2",aquifers if classic else [n for n in result.parameters if n!=x],key="unc_y_"+key)
            count=st.slider("Surface grid points per axis",3,21,7,key="unc_surface_count_"+key)
            st.caption(f"At most {(count+1)**2} forward evaluations. Other parameters stay at matched values. Bounds span the complete axes.")
            if st.button("Run objective surface",key="unc_surface_run"):
                try:
                    with st.spinner("Evaluating real forward objective grid…"):
                        surface=objective_surface(scenario,x,y,count)
                    result=replace(result,surfaces=tuple(s for s in result.surfaces if s.parameters!=(x,y))+(surface,))
                    storage[identity]=result
                except ValueError as exc:
                    st.error(str(exc))
            for i,surface in enumerate(result.surfaces):
                sx,sy=[next(s for s in specs if s.name==name) for name in surface.parameters]
                xs=[to_display(v,sx.unit,units) for v in surface.x]
                ys=[to_display(v,sy.unit,units) for v in surface.y]
                z=np.array([p.objective if p.valid else np.nan for p in surface.points]).reshape(len(ys),len(xs))
                fig=go.Figure(go.Heatmap(x=xs,y=ys,z=z,colorbar_title="Objective",colorscale="Viridis",hoverongaps=False))
                for label,values,symbol in (("Initial",surface.initial,"x"),("Matched",surface.matched,"star")):
                    fig.add_trace(go.Scatter(x=[to_display(values[0],sx.unit,units)],y=[to_display(values[1],sy.unit,units)],mode="markers",name=label,marker=dict(symbol=symbol,size=13,color="red")))
                invalid=[p for p in surface.points if not p.valid]
                if invalid:
                    fig.add_trace(go.Scatter(x=[to_display(dict(p.parameters)[sx.name],sx.unit,units) for p in invalid],y=[to_display(dict(p.parameters)[sy.name],sy.unit,units) for p in invalid],mode="markers",name="Failed forward node",marker=dict(symbol="x",color="gray")))
                fig.update_layout(title=surface.objective_definition,xaxis_title=f"{sx.name} [{display_unit(sx.unit,units)}]",yaxis_title=f"{sy.name} [{display_unit(sy.unit,units)}]")
                st.plotly_chart(fig,key=f"unc_surface_plot_{i}")
                for warning in surface.warnings:
                    st.caption(warning)
                st.dataframe(point_frame(surface.points,scenario,units),hide_index=True)
    with st.expander("Advanced numerical diagnostics"):
        st.write(f"Effective rank {result.rank}/{len(specs)}; condition indicator: {result.condition if result.condition is not None else 'rank deficient (unbounded)'}. SVD uses the parameter-scaled objective Jacobian.")
        st.write("Singular values",result.singular_values)
        st.dataframe(pd.DataFrame({"Parameter":result.parameters,"Weak-direction coefficient (scaled coordinates)":result.weak_combination,"Derivative relative change on halving step":result.derivative_errors}),hide_index=True)
        st.write("Finite difference stencils and physical steps",result.schemes)
        st.caption("Weak-vector sign is arbitrary. Opposing coefficients can indicate compensation; the vector is not a causal relationship.")
        st.write("Residual Jacobian, per canonical physical parameter")
        st.dataframe(pd.DataFrame(result.residual_jacobian,index=result.dates,columns=result.parameters))
        if result.covariance is not None:
            st.write("Local covariance in canonical physical parameter units")
            st.dataframe(pd.DataFrame(result.covariance,index=result.parameters,columns=result.parameters))
        st.write("Sensitivity column cosine coupling (not parameter covariance correlation)")
        st.dataframe(pd.DataFrame(result.column_coupling,index=result.parameters,columns=result.parameters))
        st.write(f"Local forward evaluations: {result.evaluations}; failed: {result.failed_evaluations}. Profile and surface counts are retained in their records.")
        st.download_button("Download diagnostic record (canonical SI)",json.dumps(asdict(result),default=str,indent=2),"identifiability.json","application/json")
