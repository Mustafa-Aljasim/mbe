"""Synthetic pressure histories derived analytically, without calling simulate."""
from dataclasses import replace
from datetime import date,timedelta
from material_balance_studio.domain.models import PVTProperties,PVTTable,ReservoirTank,HistoryRecord,CumulativeVolumes
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.aquifer import PotAquifer
from material_balance_studio.matching import history_match,observations_from_history,MatchParameter


def synthetic_match(kind="single",weighted=True,upper=1.5e6):
    # Eo=c*dp for exact confounding: F=(N*c+Caq)*dp. Curvature breaks the
    # exact degeneracy slightly, retaining a deliberately narrow trade-off valley.
    curvature=.03 if kind=="near" else 0.
    pressures=(5e6,10e6,20e6,24.5e6,25.4e6,26.1e6,26.9e6,27.6e6,28.3e6,29.1e6,29.7e6,30e6,40e6)
    def bo(p):
        dp=30e6-p
        return 1.2+4e-9*dp+curvature*(dp/1e7)**2
    pvt=TablePVTModel(PVTTable(tuple(PVTProperties(p,bo(p),80,.004,1.) for p in pressures)))
    coupled=kind in ("exact","near")
    true=ReservoirTank(date(2020,1,1),30e6,1e6,.2,0.,0.,0.,pvt,aquifer=PotAquifer(.005) if coupled else None)
    history=[]
    for day,p in zip((10,30,60,110,180,270,390,540),reversed(pressures[3:11])):
        dp=30e6-p
        np=(1e6*(bo(p)-1.2)+(.005*dp if coupled else 0))/bo(p)
        history.append(HistoryRecord(true.initial_date+timedelta(days=day),CumulativeVolumes(np=np,gp=80*np),p))
    base=replace(true,oil_in_place=8e5)
    specs=[MatchParameter("oil_in_place","OOIP",8e5,5e5,upper,"oil_volume",1e6)]
    if coupled:
        base=replace(base,aquifer=PotAquifer(.0058))
        specs.append(MatchParameter("aquifer.capacity","Aquifer capacity",.0058,.001,.009,"aquifer_capacity",.005))
    if kind=="zero":
        specs.append(MatchParameter("m","Gas-cap ratio",.3,0.,1.,"dimensionless",1.))
        base=replace(base,m=.3)
    observations=tuple(replace(o,sigma=1e5) for o in observations_from_history(history))
    fit=history_match(base,history,specs,observations,"Pressure uncertainty weighted" if weighted else "Unweighted")
    return true,tuple(history),fit
