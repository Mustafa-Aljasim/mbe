"""Independent oil/water storage recurrence for Phase 5C allocation acceptance."""
from dataclasses import replace
from datetime import timedelta
from math import exp
import runpy
from pathlib import Path
import numpy as np
from material_balance_studio.domain.models import PVTProperties,PVTTable,HistoryRecord,CumulativeVolumes
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.network_history import NetworkObservation,AllocationRow,HistoryPlan
from material_balance_studio.network_matching import history_match_network

base=runpy.run_path(str(Path(__file__).with_name('network_match_benchmarks.py')))
specification,SETTINGS=base['specification'],base['SETTINGS']


def oil_allocation_benchmark(aquifer=False):
    network,_,_=base['benchmark'](aquifer=aquifer)
    pvt=TablePVTModel(PVTTable(tuple(PVTProperties(p,1.2+4e-9*(30e6-p),0.,.004,1.) for p in (1e6,10e6,20e6,30e6,40e6))))
    network=replace(network,tanks=tuple(replace(t,reservoir=replace(t.reservoir,pvt_model=pvt)) for t in network.tanks))
    pressure=np.array([30e6,30e6])
    capacities=np.array([.004,.006])
    lap=8e-10*np.array([[1.,-1.],[-1.,1.]])
    totals=np.zeros((2,3))  # Np,Wp,Winj
    histories={'A':[],'B':[]}
    field=[]
    obs=[]
    pa=30e6
    we=0.
    for end in range(2,101,2):
        start=end-2
        oil=80. if start<40 else 160.
        water=30. if start<60 else 100.
        rate=np.array([[.6*oil,.2*water,0.],[.4*oil,.8*water,150. if start>=50 else 0.]])
        delta=rate*2
        old_storage=capacities-4e-9*totals[:,0]
        totals+=delta
        new_storage=capacities-4e-9*totals[:,0]
        dt=2*86400.
        beta=.005*(1-exp(-1e-9*dt/.005)) if aquifer else 0.
        matrix=np.diag(new_storage)+dt*lap/2+np.diag([0.,beta/2])
        rhs=(np.diag(old_storage)-dt*lap/2)@pressure+(new_storage-old_storage)*30e6-(1.2*delta[:,0]+delta[:,1]-delta[:,2])+np.array([0.,beta*(pa-pressure[1]/2)])
        new=np.linalg.solve(matrix,rhs)
        if aquifer:
            we+=beta*(pa-(pressure[1]+new[1])/2)
            pa=30e6-we/.005
        pressure=new
        day=network.initial_date+timedelta(days=end)
        for i,name in enumerate(('A','B')):
            histories[name].append(HistoryRecord(day,CumulativeVolumes(np=totals[i,0],wp=totals[i,1],winj=totals[i,2])))
            if end in (10,20,30,40,50,60,80,100):
                obs.append(NetworkObservation(name,day,float(pressure[i]),'Average reservoir pressure',2e4))
        field.append(HistoryRecord(day,CumulativeVolumes(np=float(totals[:,0].sum()),wp=float(totals[:,1].sum()),winj=float(totals[:,2].sum()))))
    network=replace(network,tanks=tuple(replace(t,history=tuple(histories[t.name])) for t in network.tanks))
    plan=HistoryPlan(('np','wp','winj'),tuple(field),tuple(AllocationRow(network.initial_date,s,(('A',f),('B',1-f))) for s,f in (('np',.6),('wp',.2),('winj',0.))),True)
    return network,tuple(obs),plan


def truth_match(kind='well'):
    if kind.startswith('allocation'):
        network,obs,plan=oil_allocation_benchmark(aquifer=kind.endswith('aquifer'))
        params=[('tank:A:N',1e6),('tank:B:N',1.5e6),('connection:A:B:T',8e-10)]
        if kind.endswith('aquifer'):
            params.append(('tank:B:aquifer:productivity_index',1e-9))
    elif kind in ('low_dp','redundant'):
        network,_,plan=base['benchmark']()
        node=network.tanks[0]
        reservoir=replace(node.reservoir,m=.2)
        history=tuple(HistoryRecord(network.initial_date+timedelta(days=d),CumulativeVolumes(wp=100.*d)) for d in (10,20,30,40,50,60,80,100))
        network=replace(network,tanks=tuple(replace(t,reservoir=reservoir,history=history) for t in network.tanks))
        capacity=1e6*(4e-9+.2*1.2*1e-10/.004)
        obs=tuple(NetworkObservation(t.name,h.date,30e6-h.cumulative.wp/capacity,'Average reservoir pressure',2e4) for h in history for t in network.tanks)
        params=[('connection:A:B:T',8e-10)] if kind=='low_dp' else [('tank:A:N',1e6),('tank:A:m',.2)]
    else:
        network,obs,plan=base['benchmark'](aquifer=kind=='aquifer_confounding')
        params=[('tank:A:N',1e6),('tank:B:N',1.5e6),('connection:A:B:T',8e-10)]
        if kind=='NT_confounding':
            obs=tuple(o for o in obs if o.tank=='A' and (o.date-network.initial_date).days<=40)
            params=[('tank:A:N',1e6),('connection:A:B:T',8e-10)]
        elif kind=='aquifer_confounding':
            obs=tuple(o for o in obs if o.tank=='A' and (o.date-network.initial_date).days<=40)
            params=[('tank:B:aquifer:productivity_index',1e-9),('connection:A:B:T',8e-10)]
        elif kind=='unobserved':
            obs=tuple(o for o in obs if o.tank=='A')
            params=[('tank:B:N',1.5e6),('connection:A:B:T',8e-10)]
    specs=[specification(network,name,initial=value) for name,value in params]
    return history_match_network(network,specs,obs,plan,'Pressure uncertainty weighted',settings=SETTINGS,scenario_name=kind)
