"""Independent 2x2 storage/Fetkovich recurrence, not the network solver or optimizer."""
from dataclasses import replace
from datetime import date,timedelta
from math import exp
import numpy as np
from material_balance_studio.domain.models import ReservoirTank,PVTProperties,PVTTable,HistoryRecord,CumulativeVolumes
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.aquifer import FetkovichAquifer
from material_balance_studio.tank_network import NetworkTank,Connection,TankNetwork,NetworkSettings
from material_balance_studio.network_history import NetworkObservation,AllocationRow,HistoryPlan
from material_balance_studio.network_matching import registry,MatchParameter

SETTINGS=NetworkSettings(max_step_days=2.)


def benchmark(aquifer=False,transmissibility=8e-10,gas_cap=0.,allocation=False):
    initial=date(2020,1,1)
    n_a,n_b,j=1e6,1.5e6,1e-9
    c,cg,bgi,boi=4e-9,1e-10,.004,1.2
    pvt=TablePVTModel(PVTTable(tuple(PVTProperties(p,boi+c*(30e6-p),80.,bgi+cg*(30e6-p),1.) for p in (1e6,10e6,20e6,30e6,40e6))))
    a=ReservoirTank(initial,30e6,n_a,.2,0.,0.,gas_cap,pvt)
    b=replace(a,oil_in_place=n_b,m=0.,aquifer=FetkovichAquifer(5e6,1e-9,j) if aquifer else None)
    capacities=np.diag([n_a*(c+gas_cap*boi*cg/bgi),n_b*c])
    lap=transmissibility*np.array([[1.,-1.],[-1.,1.]])
    history={"A":[],"B":[]}
    observations=[]
    field=[]
    pressure=np.array([30e6,30e6])
    pa=30e6
    we=0.
    totals=np.zeros((2,2))  # Wp, Winj
    days=(10,20,30,40,50,60,80,100)
    for end in range(2,101,2):
        start=end-2
        if allocation:
            fraction=.8 if start<40 else .3
            rates=np.array([[100*fraction,0.],[100*(1-fraction),150. if start>=50 else 0.]])
        else:
            rates=np.array([[40. if start<30 else 120.,0.],[250. if start<30 else 0.,150. if start>=50 else 0.]])
        delta=rates*2
        totals+=delta
        dt=2*86400.
        beta=.005*(1-exp(-j*dt/.005)) if aquifer else 0.
        matrix=capacities+dt*lap/2+np.diag([0.,beta/2])
        rhs=(capacities-dt*lap/2)@pressure-(delta[:,0]-delta[:,1])+np.array([0.,beta*(pa-pressure[1]/2)])
        new=np.linalg.solve(matrix,rhs)
        if aquifer:
            we+=beta*(pa-(pressure[1]+new[1])/2)
            pa=30e6-we/.005
        pressure=new
        if end in days:
            day=initial+timedelta(days=end)
            for i,name in enumerate(("A","B")):
                history[name].append(HistoryRecord(day,CumulativeVolumes(wp=totals[i,0],winj=totals[i,1])))
                observations.append(NetworkObservation(name,day,float(pressure[i]),"Average reservoir pressure",5e4))
            field.append(HistoryRecord(day,CumulativeVolumes(wp=float(totals[:,0].sum()),winj=float(totals[:,1].sum()))))
    true=TankNetwork((NetworkTank("A",a,tuple(history["A"])),NetworkTank("B",b,tuple(history["B"]))),(Connection("A","B",transmissibility),))
    plan=HistoryPlan()
    if allocation:
        plan=HistoryPlan(("wp","winj"),tuple(field),(AllocationRow(initial,"wp",(("A",.8),("B",.2))),
            AllocationRow(initial+timedelta(days=40),"wp",(("A",.3),("B",.7))),AllocationRow(initial,"winj",(("A",0.),("B",1.)))),True)
    return true,tuple(observations),plan


def specification(network,name,lower=None,upper=None,transformation="linear",initial=None):
    target=registry(network)[name]
    value=target.value if initial is None else initial
    if target.edge:
        lo,hi,scale=1e-11,2e-9,1e-9
    elif target.field=="m":
        lo,hi,scale=0.,.8,.2
    elif target.field.startswith("aquifer"):
        lo,hi,scale=2e-10,2e-9,1e-9
    else:
        lo,hi,scale=6e5,2e6,1e6
    return MatchParameter(name,target.description,value,lo if lower is None else lower,hi if upper is None else upper,target.unit,scale,transformation)
