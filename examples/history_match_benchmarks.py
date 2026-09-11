"""Analytical synthetic histories, independent of the forward simulator/optimizer."""
from dataclasses import replace
from datetime import date,timedelta
from math import exp
from material_balance_studio.domain.models import PVTProperties,PVTTable,ReservoirTank,HistoryRecord,CumulativeVolumes
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.aquifer import FetkovichAquifer
from material_balance_studio.matching.parameters import MatchParameter,parameter_registry,parameter_value


def benchmark_case(aquifer=False,gas_cap=0.):
    n,j,capacity=1e6,1e-9,.005
    table=TablePVTModel(PVTTable(tuple(PVTProperties(p,1.2+4e-9*(30e6-p),80,.004+1e-10*(30e6-p),1) for p in (5e6,10e6,20e6,30e6,40e6))))
    tank=ReservoirTank(date(2020,1,1),30e6,n,.2,0,0,gas_cap,table,aquifer=FetkovichAquifer(5e6,1e-9,j) if aquifer else None)
    records=[]
    we,pa,previous,day_before=0.,30e6,30e6,0
    for day,p in zip((10,30,60,110,180,270,390,540),(29.7e6,29.1e6,28.3e6,27.6e6,26.9e6,26.1e6,25.4e6,24.5e6)):
        if aquifer:
            we+=capacity*(pa-(previous+p)/2)*(1-exp(-j*(day-day_before)*86400/capacity))
            pa=30e6-we/capacity
        dp=30e6-p
        expansion=4e-9*dp+gas_cap*1.2*1e-10*dp/.004
        np=(n*expansion+we)/(1.2+4e-9*dp)
        records.append(HistoryRecord(tank.initial_date+timedelta(days=day),CumulativeVolumes(np=np,gp=80*np),p))
        previous,day_before=p,day
    return tank,tuple(records)


def specification(tank,name,lower=None,upper=None,transformation="linear"):
    description,unit=parameter_registry(tank)[name]
    value=parameter_value(tank,name)
    return MatchParameter(name,description,value,lower if lower is not None else value*.2,
                          upper if upper is not None else value*5,unit,max(value,.1) if name=="m" else value,transformation)
