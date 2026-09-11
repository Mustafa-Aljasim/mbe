"""Reservoir Diagnosis workspace, engineering conversions and auditable plots."""
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st
from material_balance_studio.diagnostics.diagnosis import diagnose
from material_balance_studio.diagnostics.voidage import voidage
from material_balance_studio.diagnostics.vrr import vrr
from material_balance_studio.diagnostics.drive_indices import drive_indices
from material_balance_studio.units.display import to_display, column_label, display_unit
from .aquifer_comparison import comparison_frame


def diagnostic_inspector(state,units):
    volumes,ratios,drive=voidage(state),vrr(state),drive_indices(state)
    rows=[]
    def add(label,value,quantity):
        rows.append(dict(Quantity=label,Value=to_display(value,quantity,units),Unit=display_unit(quantity,units)))
    for basis,values in volumes.items():
        for name,value in values.items():
            add(basis+" "+name.replace("_"," "),value,"reservoir_volume")
    for name,value in drive["supports"].items():
        add(name+" support",value,"reservoir_volume")
        add(name,value=drive["values"][name],quantity="dimensionless")
    add("Drive sum",drive["total"],"dimensionless")
    add("Drive closure error (sum−1)",drive["closure_error"],"dimensionless")
    add("VRR instantaneous",ratios["interval"]["injected"],"dimensionless")
    add("VRR cumulative",ratios["cumulative"]["injected"],"dimensionless")
    add("Pressure residual",state.pressure_error,"pressure")
    return pd.DataFrame(rows)


def diagnosis_frame(data,units):
    rows=[]
    for r in data["rows"]:
        row={"Date":r["date"]}
        for basis,values in r["voidage"].items():
            for name,value in values.items():
                row[column_label(basis+" "+name.replace("_"," "),"reservoir_volume",units)]=to_display(value,"reservoir_volume",units)
        for basis in ("interval","cumulative"):
            for name,value in r["vrr"][basis].items():
                row[f"VRR {basis} {name}"]=value
        row.update(r["drive"]["values"])
        row["Drive sum"]=r["drive"]["total"]
        row["Drive closure error"]=r["drive"]["closure_error"]
        rows.append(row)
    return pd.DataFrame(rows)


def time_figure(title,dates,series,quantity,units,reference=None):
    fig=go.Figure()
    for name,values in series.items():
        fig.add_scatter(x=dates,y=[to_display(v,quantity,units) for v in values],name=name,mode="lines+markers",connectgaps=False)
    fig.update_layout(title=title,xaxis_title="Date",yaxis_title=column_label(title,quantity,units))
    if reference is not None:
        fig.add_hline(y=reference,line_dash="dash",annotation_text=f"Reference = {reference}")
    return fig


def diagnosis_workflow(units):
    st.subheader("Reservoir Diagnosis")
    saved=st.session_state.get("simulation")
    if saved is None or "history" not in saved[1]:
        st.info("Run a simulation to create a reservoir/history snapshot for diagnosis.")
        return
    result,snapshot=saved
    qc=snapshot["history_qc"]
    st.write("History QC")
    st.dataframe(qc,hide_index=True,width="stretch")
    st.caption("Screening retains all rows. FAIL blocks derived diagnostics; CAUTION calls for engineering review.")
    if "FAIL" in set(qc.Status):
        st.error("Resolve history FAIL flags before building diagnostics.")
        return
    if not result.converged:
        st.warning("Forward simulation is incomplete. Diagnostics describe accepted states only; pressure matching covers a partial history.")
    if st.button("Build reservoir diagnosis",key="build_diagnosis"):
        st.session_state["reservoir_diagnosis"]=(diagnose(result,snapshot["tank"],snapshot["history"],snapshot["pressure_metadata"]),result,snapshot)
    stored=st.session_state.get("reservoir_diagnosis")
    if not stored:
        return
    data,diagnosed_result,diagnosed_snapshot=stored
    if diagnosed_result is not result:
        st.info("Simulation changed. Build reservoir diagnosis again to review the new snapshot.")
        return
    st.caption("Saved simulation snapshot. H–O/Campbell use observed pressure; VRR and drive indices use accepted calculated-pressure terms. No parameter optimization is performed.")
    match=data["pressure"]
    metrics=match["metrics"]
    st.write("Engineering summary")
    st.write(f"History: {'CAUTION' if 'CAUTION' in set(qc.Status) else 'PASS'} · {metrics['count']} pressure observations. Aquifer: {result.initial_state.aquifer_state.model_key}.")
    st.metric(column_label("Pressure RMSE","pressure",units),"Unavailable" if metrics["rmse"] is None else f"{to_display(metrics['rmse'],'pressure',units):.6g}")
    if data["rows"]:
        last=data["rows"][-1]
        st.write(f"Latest cumulative VRR: {last['vrr']['cumulative']['injected']}. Current model support fractions:",last["drive"]["values"])
    st.write(f"H–O screening: {data['ho']['fit']['status']}.")
    for note in data["ho"]["fit"]["notes"]:
        st.warning(note)
    frame=diagnosis_frame(data,units)
    with st.expander("Pressure / Production / Injection and voidage"):
        st.caption("Oil and water components are explicit; free-gas contribution is the residual of the validated production withdrawal. Signed gas is retained. Interval voidage uses surface increments at endpoint PVT, not differences of cumulative volumes valued at different pressures.")
        st.dataframe(frame,hide_index=True,width="stretch")
        st.download_button("Download reservoir diagnostics",frame.to_csv(index=False),f"reservoir_diagnosis_{units}.csv")
        dates=[r["date"] for r in data["rows"]]
        st.plotly_chart(time_figure("Cumulative reservoir volume",dates,{k:[r["voidage"]["cumulative"][k] for r in data["rows"]] for k in ("produced","water_injected","gas_injected","net")},"reservoir_volume",units),key="diagnosis_voidage")
    with st.expander("VRR",expanded=True):
        st.caption("VRR = injected reservoir volume / produced reservoir voidage. VRR = 1 does not guarantee constant pressure: aquifer support, expansion, compressibility and reservoir communication also matter.")
        for basis,title in (("interval","Instantaneous VRR"),("cumulative","Cumulative VRR")):
            st.plotly_chart(time_figure(title,dates,{k:[r["vrr"][basis][k] for r in data["rows"]] for k in ("water_injected","gas_injected","injected")},"dimensionless",units,1),key="diagnosis_vrr_"+basis)
        for warning in dict.fromkeys(w for r in data["rows"] for w in r["vrr"]["warnings"]):
            st.warning(warning)
        combined=make_subplots(specs=[[{"secondary_y":True}]])
        combined.add_scatter(x=dates,y=[to_display(s.pressure,"pressure",units) for s in result.states],name="Calculated pressure",secondary_y=False)
        combined.add_scatter(x=dates,y=[to_display(s.observed_pressure,"pressure",units) for s in result.states],name="Observed pressure",mode="markers",secondary_y=False)
        combined.add_scatter(x=dates,y=[r["vrr"]["cumulative"]["injected"] for r in data["rows"]],name="Cumulative VRR",secondary_y=True)
        combined.add_hline(y=1,secondary_y=True,line_dash="dash")
        combined.update_yaxes(title_text=column_label("Pressure","pressure",units),secondary_y=False)
        combined.update_yaxes(title_text="VRR (dimensionless)",secondary_y=True)
        st.plotly_chart(combined,key="diagnosis_pressure_vrr")
    with st.expander("Havlena–Odeh",expanded=True):
        st.caption("Observed-pressure points only. Injection is subtracted: Fnet = Fproduction − water injection − gas injection. We is replayed from the supplied aquifer against complete observed pressure history. The current gas-cap ratio and compressibilities remain assumptions; high R² is not physical validation and OOIP is never overwritten.")
        for index,ho in enumerate(h for h in (data["ho"],data["ho_gas"]) if h):
            st.code(ho["equation"],language=None)
            fit=ho["fit"]
            st.write(f"{fit['status']} · Usable points: {fit['count']} · R²: {fit['r2']}")
            st.write({"Slope": to_display(fit["slope"],"oil_volume",units),
                      "Intercept": to_display(fit["intercept"],ho["y_quantity"],units),
                      "OOIP inferred from":ho["inferred_from"]})
            st.caption(f"Slope units: {display_unit('oil_volume',units)}; intercept units: {display_unit(ho['y_quantity'],units)}.")
            for note in fit["notes"]:
                st.warning(note)
            points=ho["points"]
            fig=go.Figure(go.Scatter(x=[to_display(p["x"],ho["x_quantity"],units) for p in points],y=[to_display(p["y"],ho["y_quantity"],units) for p in points],text=[p["date"] for p in points],mode="markers",name="Observed-pressure transform"))
            if fit["slope"] is not None:
                xs=sorted(p["x"] for p in points if p["x"] is not None and p["y"] is not None)
                fig.add_scatter(x=[to_display(x,ho["x_quantity"],units) for x in xs],y=[to_display(fit["slope"]*x+fit["intercept"],ho["y_quantity"],units) for x in xs],mode="lines",name="Diagnostic trend")
            fig.update_layout(xaxis_title=column_label(ho["x_label"],ho["x_quantity"],units),yaxis_title=column_label(ho["y_label"],ho["y_quantity"],units))
            st.plotly_chart(fig,key=f"diagnosis_ho_{index}")
        components=pd.DataFrame(data["components"])
        for col in ("pressure","F","Eo","Eg","Efw","Et","We","injection","adjusted"):
            quantity="pressure" if col=="pressure" else "expansion" if col in ("Eo","Eg","Efw","Et") else "reservoir_volume"
            components[col]=components[col].map(lambda v:to_display(v,quantity,units))
            components.rename(columns={col:column_label(col,quantity,units)},inplace=True)
        st.dataframe(components,hide_index=True,width="stretch")
        for row in data["components"]:
            if row["note"]:
                st.warning(row["date"]+": "+row["note"])
    with st.expander("Campbell / Cole"):
        st.caption("Oil-reservoir Campbell plot. Dry-gas Cole is not applicable to this black-oil tank. X = net cumulative withdrawal; Y = Fnet/Et (apparent oil volume), with aquifer-adjusted overlay. Points require Et > 1e-12. Trends alone do not identify the drive mechanism.")
        st.code(data["campbell"]["equation"],language=None)
        st.write(data["campbell"]["cue"])
        points=data["campbell"]["points"]
        fig=go.Figure()
        for key,label in (("unadjusted","Fnet/Et"),("adjusted","(Fnet-We)/Et")):
            fig.add_scatter(x=[to_display(p["x"],"reservoir_volume",units) for p in points],y=[to_display(p[key],"oil_volume",units) for p in points],text=[p["date"] for p in points],mode="lines+markers",name=label,connectgaps=False)
        fig.update_layout(xaxis_title=column_label("Fnet","reservoir_volume",units),yaxis_title=column_label("Apparent oil volume","oil_volume",units))
        st.plotly_chart(fig,key="diagnosis_campbell")
    with st.expander("Drive Indices",expanded=True):
        st.caption("Common denominator Fproduction: DDI=N·Eo/Fproduction; SDI=N·mEg/Fproduction; CDI=N·Efw/Fproduction; WDI=We/Fproduction; IDI=(water+gas injection support)/Fproduction. Sum−1 = −MBE residual/Fproduction. No renormalization or clipping.")
        series={k:[r["drive"]["values"][k] for r in data["rows"]] for k in ("DDI","SDI","CDI","WDI","IDI")}
        st.plotly_chart(time_figure("Drive/support indices",dates,series,"dimensionless",units,1),key="diagnosis_drive")
        active=[r for r in data["rows"] if r["drive"]["total"] is not None]
        if active and all(r["drive"]["stackable"] for r in active):
            fig=go.Figure()
            for key in series:
                fig.add_scatter(x=[r["date"] for r in active],y=[100*r["drive"]["values"][key] for r in active],stackgroup="support",name=key)
            fig.update_layout(yaxis_title="Support (%) — unrenormalized",xaxis_title="Date")
            st.plotly_chart(fig,key="diagnosis_drive_stack")
        else:
            st.warning("Stacked view withheld: support fractions are undefined, signed, or do not close within 1e-8.")
    with st.expander("Pressure Match",expanded=True):
        st.caption("Missing observations remain missing. Sigma normalizes residuals for review only; regression and metrics are unweighted.")
        st.write({k:to_display(v,"pressure",units) for k,v in metrics.items() if k!="count"})
        pressures=pd.DataFrame(match["rows"])
        for col in ("calculated","observed","residual","sigma"):
            pressures[col]=pressures[col].map(lambda v:to_display(v,"pressure",units))
            pressures.rename(columns={col:column_label(col,"pressure",units)},inplace=True)
        st.dataframe(pressures,hide_index=True,width="stretch")
        for label,keys,quantity in (("Pressure",("calculated","observed"),"pressure"),("Pressure residual",("residual",),"pressure"),("Normalized pressure residual",("normalized",),"dimensionless")):
            st.plotly_chart(time_figure(label,dates,{k:[r[k] for r in match["rows"]] for k in keys},quantity,units,0 if "residual" in label else None),key="diagnosis_match_"+label)
    if result.states:
        selected=st.selectbox("Diagnosis equation inspector: date",[s.date for s in result.states],key="diagnosis_date")
        st.dataframe(diagnostic_inspector(next(s for s in result.states if s.date==selected),units),hide_index=True,width="stretch")
    with st.expander("Aquifer Comparison"):
        st.caption("Use the Aquifer Comparison tab to configure/run models. This view reuses its saved results; its snapshot may differ from the current diagnosis.")
        comparison=st.session_state.get("aquifer_comparison")
        if comparison:
            st.dataframe(comparison_frame(comparison[0],units),hide_index=True,width="stretch")
