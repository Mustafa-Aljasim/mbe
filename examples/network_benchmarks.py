"""Known linear storage capacities embedded in existing black-oil component calls."""
from dataclasses import replace
from datetime import date,timedelta
from material_balance_studio.domain.models import ReservoirTank,PVTProperties,PVTTable,CumulativeVolumes,HistoryRecord
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.aquifer import FetkovichAquifer
from material_balance_studio.tank_network import NetworkTank,Connection,TankNetwork


def storage_tank(pressure=30e6,n=1e6,c=4e-9):
    pvt=TablePVTModel(PVTTable(tuple(PVTProperties(p,1.2+c*(30e6-p),80.,.004,1.) for p in (1e6,10e6,20e6,30e6,40e6))))
    return ReservoirTank(date(2020,1,1),pressure,n,.2,0.,0.,0.,pvt)


def network_case(kind="closed",transmissibility=1e-9):
    a=storage_tank()
    b=storage_tank(20e6,2e6) if kind=="closed" else storage_tank(n=2e6)
    if kind=="aquifer":
        b=replace(b,aquifer=FetkovichAquifer(5e6,1e-9,1e-9))
    ha,hb=[],[]
    for day in (10,20,40,70,100):
        produced=day*100. if kind in ("production","aquifer","injection","reversal") else 0.
        # Water production has pressure-independent Bw=1 so analytic linear storage
        # benchmarks remain possible without recomputing oil-phase inventory.
        injected=max(0.,day-20)*400. if kind in ("injection","reversal") else 0.
        if kind=="reversal":
            # Produce B early (A -> B), then inject B (B -> A).
            ca=CumulativeVolumes()
            cb=CumulativeVolumes(wp=min(day,20)*100.,winj=injected)
        else:
            ca=CumulativeVolumes(wp=produced)
            cb=CumulativeVolumes(winj=injected)
        ha.append(HistoryRecord(a.initial_date+timedelta(days=day),ca))
        hb.append(HistoryRecord(a.initial_date+timedelta(days=day),cb))
    return TankNetwork((NetworkTank("A",a,tuple(ha)),NetworkTank("B",b,tuple(hb))),
                       (Connection("A","B",transmissibility),))
