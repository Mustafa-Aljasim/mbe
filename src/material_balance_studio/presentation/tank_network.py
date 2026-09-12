"""Separate multi-tank forward workspace; canonical snapshots never alter single tank."""
from dataclasses import asdict,replace,fields
from datetime import date,timedelta
import json
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from material_balance_studio.domain.models import ReservoirTank,PVTProperties,PVTTable,HistoryRecord,CumulativeVolumes,CUMULATIVE_FIELDS
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.io.history import read_table,history_from_frame,pvt_from_frame
from material_balance_studio.aquifer.registry import MODELS
from material_balance_studio.tank_network import NetworkTank,Connection,TankNetwork,NetworkSettings,simulate_network
from material_balance_studio.tank_network.diagnostics import tank_terms,maxima
from material_balance_studio.units.display import from_display,to_display,display_unit
from material_balance_studio.units.conversions import from_si
from .aquifer import FIELDS
from .tables import history_display_frame


def example_network():
    pvt=TablePVTModel(PVTTable(tuple(PVTProperties(p,1.2+4e-9*(30e6-p),80.,.004,1.) for p in (1e6,10e6,20e6,30e6,40e6))))
    a=ReservoirTank(date(2020,1,1),30e6,1e6,.2,0.,0.,0.,pvt)
    history=tuple(HistoryRecord(a.initial_date+timedelta(days=d)) for d in (10,20,40,70,100))
    return TankNetwork((NetworkTank("A",a,history),NetworkTank("B",replace(a,initial_pressure=20e6,oil_in_place=2e6),history)),(Connection("A","B",1e-9),))


def tank_frame(result,units,internal=False):
    rows=[]
    for state in (result.internal_states if internal else result.states):
        for tank in state.tanks:
            row={"Time":state.time,"Tank":tank.name}
            for name,value in tank_terms(tank).items():
                quantity="pressure" if name=="pressure" else "dimensionless" if name in ("relative_residual","transfer_support_ratio") else "reservoir_volume"
                row[f"{name} [{display_unit(quantity,units)}]"]=to_display(value,quantity,units)
            rows.append(row)
    return pd.DataFrame(rows)


def connection_frame(result,units,internal=False):
    rows=[]
    for state in (result.internal_states if internal else result.states):
        for edge in state.connections:
            rows.append({"Time":state.time,"From":edge.from_tank,"To":edge.to_tank,"Enabled":edge.enabled,
                f"T [{display_unit('aquifer_productivity',units)}]":to_display(edge.transmissibility,"aquifer_productivity",units),
                f"ΔP [{display_unit('pressure',units)}]":to_display(edge.pressure_difference,"pressure",units),
                f"q [{display_unit('aquifer_rate',units)}]":to_display(edge.rate,"aquifer_rate",units),
                f"ΔV last substep [{display_unit('reservoir_volume',units)}]":to_display(edge.incremental_transfer,"reservoir_volume",units),
                f"Cumulative signed transfer [{display_unit('reservoir_volume',units)}]":to_display(edge.cumulative_transfer,"reservoir_volume",units)})
    return pd.DataFrame(rows)


def _save(nodes,edges):
    network=TankNetwork(tuple(nodes),tuple(edges))
    st.session_state["nt_network"]=network
    st.session_state["nt_revision"]=st.session_state.get("nt_revision",0)+1


def _history_input(history,units):
    return pd.DataFrame([dict(date=str(h.date),**{n:from_si(getattr(h.cumulative,n),"gas_volume" if n in ("gp","ginj") else "liquid_volume",units) for n in CUMULATIVE_FIELDS},
        observed_pressure=None if h.observed_pressure is None else from_si(h.observed_pressure,"pressure",units)) for h in history],
        columns=["date",*CUMULATIVE_FIELDS,"observed_pressure"])


def multi_tank_workflow(units):
    st.subheader("Multi-Tank Forward Model")
    st.caption("Separate compartment pressures and explicit connections. Forward modeling only; the single-tank project and accepted matches remain independent.")
    if not st.checkbox("Configure multi-tank network",key="nt_configure"):
        return
    st.info("Communication transfers lumped reservoir volume, not tracked fluid phases/composition. Each tank retains its assigned PVT. Histories interpolate cumulative streams linearly onto a common timeline and hold after their last date; observed pressures are not imposed.")
    if st.button("Load two-tank equilibration example",key="nt_demo"):
        demo=example_network()
        _save(demo.tanks,demo.connections)
        st.rerun()
    network=st.session_state.get("nt_network")
    nodes=list(network.tanks) if network else []
    edges=list(network.connections) if network else []
    name=st.text_input("New tank name",value="Tank C" if nodes else "Tank A",key="nt_new_name")
    if st.button("Add tank",key="nt_add"):
        try:
            snapshot=st.session_state.get("simulation")
            tank=snapshot[1]["tank"] if snapshot else example_network().tanks[0].reservoir
            if nodes:
                tank=replace(tank,initial_date=nodes[0].reservoir.initial_date)
            _save(nodes+[NetworkTank(name,tank)],edges)
            st.rerun()
        except ValueError as exc:
            st.error(str(exc))
    st.caption("New tanks copy the saved single-tank reservoir/PVT/aquifer when available, otherwise a simple example PVT. New histories start empty. Tanks must share an initial reference date.")
    if not network:
        return
    rev=st.session_state.get("nt_revision",0)
    names=[t.name for t in nodes]
    selected=st.selectbox("Configure tank",names,key=f"nt_selected_{rev}")
    index=names.index(selected)
    node=nodes[index]
    prefix=f"nt_{rev}_{selected}_{units}"
    if st.button("Remove selected tank and its connections",key="nt_remove"):
        remaining=[n for n in nodes if n.name!=selected]
        if remaining:
            _save(remaining,[e for e in edges if selected not in e.key])
        else:
            st.session_state.pop("nt_network",None)
        st.rerun()
    with st.expander("Tank properties and pressure bounds",expanded=True):
        st.caption("Apply edits before switching display units. Pressure bounds intersect this tank's PVT coverage; pressures above Pi are allowed.")
        quantities={"initial_pressure":"pressure","oil_in_place":"oil_volume","m":"dimensionless","swc":"dimensionless","cf":"compressibility","cw":"compressibility"}
        with st.form(prefix+"properties"):
            values={n:from_display(st.number_input(f"{n} [{display_unit(q,units)}]",value=float(to_display(getattr(node.reservoir,n),q,units)),format="%.10g",key=prefix+n),q,units) for n,q in quantities.items()}
            low,high=node.pressure_bounds
            low=from_display(st.number_input(f"Minimum pressure [{display_unit('pressure',units)}]",value=float(to_display(low,"pressure",units)),format="%.10g",key=prefix+"low"),"pressure",units)
            high=from_display(st.number_input(f"Maximum pressure [{display_unit('pressure',units)}]",value=float(to_display(high,"pressure",units)),format="%.10g",key=prefix+"high"),"pressure",units)
            newdate=st.date_input("Common initial date (applied to all tanks)",value=network.initial_date,key=prefix+"date")
            if st.form_submit_button("Apply tank properties"):
                try:
                    updated=[replace(n,reservoir=replace(n.reservoir,initial_date=newdate)) for n in nodes]
                    updated[index]=replace(node,reservoir=replace(node.reservoir,initial_date=newdate,**values),minimum_pressure=low,maximum_pressure=high)
                    _save(updated,edges)
                    st.rerun()
                except ValueError as exc:
                    st.error(str(exc))
    with st.expander("Assign PVT"):
        choices=["Keep assigned PVT"]+[f"Share PVT from {n}" for n in names if n!=selected]
        snapshot=st.session_state.get("simulation")
        if snapshot:
            choices.append("Use saved single-tank PVT")
        source=st.selectbox("PVT reference",choices,key=prefix+"pvt_ref")
        if st.button("Assign selected PVT reference",key=prefix+"assign_pvt"):
            try:
                model=node.reservoir.pvt_model if source==choices[0] else snapshot[1]["tank"].pvt_model if source=="Use saved single-tank PVT" else nodes[names.index(source.removeprefix("Share PVT from "))].reservoir.pvt_model
                updated=nodes.copy()
                updated[index]=replace(node,reservoir=replace(node.reservoir,pvt_model=model),minimum_pressure=None,maximum_pressure=None)
                _save(updated,edges)
                st.rerun()
            except ValueError as exc:
                st.error(str(exc))
        upload=st.file_uploader("Independent PVT CSV/XLSX",type=["csv","xlsx"],key=prefix+"pvt_file")
        file_units=st.selectbox("PVT file units (SI uses Pa)",["SI","FIELD"],key=prefix+"pvt_units")
        if st.button("Load tank PVT",disabled=upload is None,key=prefix+"pvt_load"):
            try:
                upload.seek(0)
                model=TablePVTModel(pvt_from_frame(read_table(upload),file_units))
                updated=nodes.copy()
                updated[index]=replace(node,reservoir=replace(node.reservoir,pvt_model=model),minimum_pressure=None,maximum_pressure=None)
                _save(updated,edges)
                st.rerun()
            except ValueError as exc:
                st.error(str(exc))
        st.caption(f"Assigned {type(node.reservoir.pvt_model).__name__}; PVT range {to_display(node.reservoir.pvt_model.pressure_bounds[0],'pressure',units):g}–{to_display(node.reservoir.pvt_model.pressure_bounds[1],'pressure',units):g} {display_unit('pressure',units)}. Assigning PVT resets custom pressure bounds to its coverage.")
    with st.expander("Assign aquifer"):
        model_name=next(n for n,c in MODELS.items() if (node.reservoir.aquifer.key if node.reservoir.aquifer else "none")==c.key)
        aq_name=st.selectbox("Tank aquifer",list(MODELS),index=list(MODELS).index(model_name),key=prefix+"aq_name")
        cls=MODELS[aq_name]
        with st.form(prefix+aq_name+"aquifer"):
            aq_values={}
            for f in fields(cls):
                label,q,default,help_text=FIELDS[f.name]
                value=getattr(node.reservoir.aquifer,f.name,default) if aq_name==model_name else default
                aq_values[f.name]=from_display(st.number_input(f"{label} [{display_unit(q,units)}]",value=float(to_display(value,q,units)),format="%.10g",help=help_text,key=prefix+aq_name+f.name),q,units)
            if st.form_submit_button("Apply tank aquifer"):
                try:
                    updated=nodes.copy()
                    updated[index]=replace(node,reservoir=replace(node.reservoir,aquifer=cls(**aq_values)))
                    _save(updated,edges)
                    st.rerun()
                except ValueError as exc:
                    st.error(str(exc))
    with st.expander("Tank production / injection / observed-pressure history"):
        file_units=st.selectbox("History input units (SI pressure is Pa)",["SI","FIELD"],key=prefix+"history_units")
        st.caption("Required columns: date, np, gp, wp. Optional: winj, ginj, observed_pressure. Add zero-volume dated rows for closed equilibration or additional events.")
        upload=st.file_uploader("Tank history CSV/XLSX",type=["csv","xlsx"],key=prefix+"history_file")
        if st.button("Load tank history",disabled=upload is None,key=prefix+"history_load"):
            try:
                upload.seek(0)
                history=history_from_frame(read_table(upload),file_units)
                updated=nodes.copy()
                updated[index]=replace(node,history=history)
                _save(updated,edges)
                st.rerun()
            except ValueError as exc:
                st.error(str(exc))
        edited=st.data_editor(_history_input(node.history,file_units),num_rows="dynamic",key=prefix+file_units+"history_edit",hide_index=True)
        if st.button("Apply edited tank history",key=prefix+"history_apply"):
            try:
                updated=nodes.copy()
                updated[index]=replace(node,history=history_from_frame(edited,file_units))
                _save(updated,edges)
                st.rerun()
            except ValueError as exc:
                st.error(str(exc))
    with st.expander("Connections",expanded=True):
        st.caption("Each physical link appears once. Positive rate and signed volume follow From → To; the same edge can reverse. Disabled links transfer nothing.")
        rows=pd.DataFrame([dict(From=e.from_tank,To=e.to_tank,T=to_display(e.transmissibility,"aquifer_productivity",units),Enabled=e.enabled) for e in edges],columns=["From","To","T","Enabled"])
        edited=st.data_editor(rows,num_rows="dynamic",hide_index=True,key=prefix+"edges",column_config={
            "From":st.column_config.SelectboxColumn(options=names,required=True),"To":st.column_config.SelectboxColumn(options=names,required=True),
            "T":st.column_config.NumberColumn(f"T [{display_unit('aquifer_productivity',units)}]",min_value=0.,required=True),"Enabled":st.column_config.CheckboxColumn(default=True,required=True)})
        if st.button("Apply connections",key="nt_edges_apply"):
            try:
                proposed=[Connection(str(r["From"]),str(r["To"]),from_display(float(r["T"]),"aquifer_productivity",units),r["Enabled"]) for r in edited.to_dict("records")]
                _save(nodes,proposed)
                st.rerun()
            except (ValueError,TypeError) as exc:
                st.error(str(exc))
    st.dataframe(pd.DataFrame([dict(Tank=n.name,Initial_date=n.reservoir.initial_date,PVT=type(n.reservoir.pvt_model).__name__,Aquifer=n.reservoir.aquifer.key if n.reservoir.aquifer else "none",History_rows=len(n.history),Neighbors=", ".join(e.to_tank if e.from_tank==n.name else e.from_tank for e in edges if e.enabled and n.name in e.key)) for n in nodes]),hide_index=True)
    max_days=st.number_input("Maximum internal step (days; 0 uses communication control only)",min_value=0.,value=0.,key="nt_max_days")
    st.caption("Communication stiffness may refine steps further. Every tank advances on every substep, including tanks with no production. Large/stiff networks can require many coupled evaluations.")
    if st.button("Run multi-tank forward simulation",key="nt_run",type="primary"):
        try:
            with st.spinner("Solving simultaneous compartment pressures…"):
                st.session_state["nt_result"]=simulate_network(network,NetworkSettings(max_step_days=max_days or None))
        except ValueError as exc:
            st.error(str(exc))
    result=st.session_state.get("nt_result")
    if result is not None:
        show_network_results(result,units)


def show_network_results(result,units):
    st.subheader("Multi-Tank Results")
    st.caption("Results belong to the saved run snapshot. Later configuration edits require another run.")
    if result.converged:
        st.success("PASS: all tanks close individually and paired transfer is conserved at every accepted substep.")
    else:
        st.error(f"FAIL at {result.failed_timestep.time}: {result.failed_timestep.diagnostics.message}")
        st.write(asdict(result.failed_timestep.diagnostics))
        if result.failed_timestep.candidate:
            st.write("Failed candidate tank residuals (canonical SI)",[(s.name,s.residual,s.relative_residual) for s in result.failed_timestep.candidate.tanks])
    for warning in result.warnings:
        st.warning(warning)
    stats=maxima(result)
    st.write("Accepted-state QC",stats)
    states=(result.initial_state,*result.internal_states)
    fig=go.Figure()
    for j,node in enumerate(result.network.tanks):
        fig.add_trace(go.Scatter(x=[s.time for s in states],y=[to_display(s.tanks[j].pressure,"pressure",units) for s in states],mode="lines",name=node.name))
        obs=[h for h in node.history if h.observed_pressure is not None]
        if obs:
            fig.add_trace(go.Scatter(x=[h.date for h in obs],y=[to_display(h.observed_pressure,"pressure",units) for h in obs],mode="markers",name=node.name+" observed"))
    fig.update_layout(title="Tank pressures",yaxis_title=display_unit("pressure",units))
    st.plotly_chart(fig,key="nt_pressure_plot")
    with st.expander("Tank balances / aquifer / inter-tank support",expanded=True):
        st.dataframe(tank_frame(result,units),hide_index=True)
        for j,node in enumerate(result.network.tanks):
            fig=go.Figure()
            for label,values in (("Aquifer We",[s.tanks[j].components.aquifer_support for s in states]),("Inter-tank X",[s.tanks[j].intertank_support for s in states]),("Direct injection",[s.tanks[j].components.withdrawal.water_injection+s.tanks[j].components.withdrawal.gas_injection for s in states])):
                fig.add_trace(go.Scatter(x=[s.time for s in states],y=[to_display(v,"reservoir_volume",units) for v in values],name=label))
            fig.update_layout(title=node.name+" · separate supports",yaxis_title=display_unit("reservoir_volume",units))
            st.plotly_chart(fig,key="nt_support_"+node.name)
        st.caption("X/produced voidage is a signed transfer-support diagnostic; undefined for zero production. Phase 4A drive indices are unchanged.")
    with st.expander("Inter-tank transfer / connectivity"):
        st.dataframe(connection_frame(result,units),hide_index=True)
        for label,field,quantity in (("Pressure difference","pressure_difference","pressure"),("Endpoint transfer rate","rate","aquifer_rate"),("Cumulative signed transfer","cumulative_transfer","reservoir_volume")):
            fig=go.Figure()
            for j,edge in enumerate(result.network.connections):
                fig.add_trace(go.Scatter(x=[s.time for s in states],y=[to_display(getattr(s.connections[j],field),quantity,units) for s in states],name=edge.from_tank+" → "+edge.to_tank))
            fig.update_layout(title=label,yaxis_title=display_unit(quantity,units))
            st.plotly_chart(fig,key="nt_connection_"+field)
    with st.expander("Production / injection histories"):
        for node in result.network.tanks:
            st.write(node.name)
            if node.history:
                st.dataframe(history_display_frame(node.history,units),hide_index=True)
    with st.expander("Equation inspector and network inspector",expanded=True):
        if result.internal_states:
            selected=st.selectbox("Network timestep",range(len(result.internal_states)),format_func=lambda i:str(result.internal_states[i].time),key="nt_inspect_time")
            state=result.internal_states[selected]
            name=st.selectbox("Inspect tank",[t.name for t in state.tanks],key="nt_inspect_tank")
            tank=next(t for t in state.tanks if t.name==name)
            st.write("Fnet = N Eo + N mEg + N Efw + We + X. Positive X is net support received from adjacent tanks.")
            st.dataframe(tank_frame(replace(result,states=(state,)),units).query("Tank == @name"),hide_index=True)
            st.dataframe(pd.DataFrame([{"Neighbor":n,f"Signed cumulative contribution into {name} [{display_unit('reservoir_volume',units)}]":to_display(v,"reservoir_volume",units)} for n,v in tank.connection_contributions]),hide_index=True)
            for warning in tank.warnings:
                st.warning(warning)
            st.dataframe(connection_frame(replace(result,states=(state,)),units),hide_index=True)
            st.write(f"Sum tank cumulative transfer = {to_display(state.transfer_error,'reservoir_volume',units):.6e} {display_unit('reservoir_volume',units)}; expected ≈ 0.")
            st.write(asdict(state.diagnostics))
    with st.expander("Advanced closure / internal timestep records"):
        fig=go.Figure()
        for j,node in enumerate(result.network.tanks):
            fig.add_trace(go.Scatter(x=[s.time for s in states],y=[s.tanks[j].relative_residual for s in states],name=node.name))
        fig.add_hline(y=result.settings.relative_tolerance,line_dash="dash")
        fig.update_layout(title="Individual tank relative MBE closure",yaxis=dict(rangemode="tozero",title="Absolute relative residual"))
        st.plotly_chart(fig,key="nt_closure")
        fig=go.Figure(go.Scatter(x=[s.time for s in states],y=[to_display(abs(s.transfer_error),"reservoir_volume",units) for s in states],name="Cumulative transfer error"))
        fig.add_hline(y=to_display(result.settings.transfer_absolute_tolerance,"reservoir_volume",units),line_dash="dash")
        fig.update_layout(title="Network transfer-conservation error",yaxis=dict(rangemode="tozero",title=display_unit("reservoir_volume",units)))
        st.plotly_chart(fig,key="nt_conservation")
        st.dataframe(tank_frame(result,units,True),hide_index=True)
        st.dataframe(connection_frame(result,units,True),hide_index=True)
        st.download_button("Download network run (canonical SI)",json.dumps(asdict(result),default=str,indent=2,allow_nan=False),"network_run.json","application/json")
