"""Transparent interface-equivalence fixture, not a new empirical correlation.

Temporarily supplies identical table functions through the correlation registry.
The real CorrelationPVTModel dispatch and unchanged MBE solver execute both runs.
Gas-cap and rock/water terms are deliberately nonzero. Run in a separate process.
"""
from dataclasses import replace
from datetime import date
import json
import numpy as np
from material_balance_studio.domain.models import PVTProperties, PVTTable, ReservoirTank, CumulativeVolumes, HistoryRecord
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.pvt.fluid import fluid_from_inputs
from material_balance_studio.pvt.models.correlation_model import CorrelationPVTModel, PropertySelection
from material_balance_studio.pvt.correlations.base import Correlation
from material_balance_studio.pvt.correlations.registry import REGISTRY, register
from material_balance_studio.solver.simulation import simulate


def comparison():
    table=TablePVTModel(PVTTable(tuple(PVTProperties(p,bo,80,bg,bw) for p,bo,bg,bw in (
        (10e6,1.28,.012,1.008),(20e6,1.24,.006,1.004),(30e6,1.2,.004,1.),(40e6,1.16,.003,.996)))))
    selection=vars(PropertySelection()).copy()
    temporary=[]
    try:
        for name in ("bo","rs","bg","bw"):
            key="acceptance_reference_"+name
            register(Correlation(key,name,"Identical reference function",
                lambda i,n=name:getattr(table.properties_at_pressure(i.pressure),n),
                "Engineering interface acceptance fixture", "Exactly the same continuous tabular function"))
            temporary.append(key)
            selection[name]=key
        fluid=fluid_from_inputs(api=35,gas_gravity=.75,temperature=90,initial_pressure=30e6,pb=40e6,rsb=80)
        model=CorrelationPVTModel(fluid,(10e6,40e6),PropertySelection(**selection))
        property_error=max(abs(getattr(table.properties_at_pressure(float(p)),name)-model.property_at_pressure(name,float(p)))
                           for p in np.linspace(10e6,40e6,1001) for name in ("bo","rs","bg","bw"))
        tank=ReservoirTank(date(2020,1,1),30e6,1e6,.2,5e-10,4e-10,.15,table)
        history=tuple(HistoryRecord(date(2020,month,1),CumulativeVolumes(np=n,gp=100*n,wp=.1*n,winj=.05*n,ginj=n))
                      for month,n in ((2,5000),(3,20000),(4,50000)))
        a,b=simulate(tank,history),simulate(replace(tank,pvt_model=model),history)
        if not a.converged or not b.converged:
            raise RuntimeError("Acceptance comparison did not converge.")
        def terms(state):
            return {"pressure_Pa":state.pressure,"Fnet_m3":state.balance.withdrawal.net,
                    "Eo":state.balance.expansion.eo,"mEg":state.balance.expansion.meg,
                    "Efw":state.balance.expansion.efw,"residual_m3":state.balance.residual}
        rows=[]
        for left,right in zip(a.states,b.states):
            l,r=terms(left),terms(right)
            rows.append({"date":str(left.date),"table":l,"correlation":r,
                         "absolute_difference":{name:abs(l[name]-r[name]) for name in l}})
        return {"maximum_property_difference":property_error,"both_converged":True,"steps":rows,
                "tolerances":{"pressure_Pa":1e-6,"Fnet_m3":1e-7,"Eo":1e-12,"mEg":1e-12,"Efw":1e-12,"residual_m3":1e-7},
                "scope":"Synthetic identical-functions interface check; no empirical accuracy claim."}
    finally:
        for key in temporary:
            REGISTRY.pop(key)


if __name__=="__main__":
    print(json.dumps(comparison(),indent=2))
