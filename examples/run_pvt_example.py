"""Reproducible Phase 2 screening/matching/QC and unchanged MBE integration."""
from dataclasses import asdict, replace
from datetime import date
from pathlib import Path
import json
import pandas as pd
from material_balance_studio.domain.models import ReservoirTank
from material_balance_studio.io.history import history_from_frame
from material_balance_studio.pvt.fluid import fluid_from_inputs
from material_balance_studio.pvt.laboratory import lab_from_frame
from material_balance_studio.pvt.models.correlation_model import CorrelationPVTModel
from material_balance_studio.pvt.matching.ranking import screen_property
from material_balance_studio.pvt.matching.regression import fit_property
from material_balance_studio.pvt.qc.physical_checks import physical_qc, failures
from material_balance_studio.solver.simulation import simulate


def main():
    folder=Path(__file__).resolve().parent
    fluid=fluid_from_inputs(api=35,gas_gravity=.75,temperature=90,initial_pressure=30e6,pb=18e6,rsb=90)
    raw=CorrelationPVTModel(fluid,(2e6,35e6))
    lab=lab_from_frame(pd.read_csv(folder/'correlation_lab_SI.csv'))
    for candidate in screen_property(raw,lab,'bo'):
        print(candidate.correlation,candidate.applicability.status,candidate.metrics)
    fit=fit_property(raw,lab,'bo')
    if not fit.accepted:
        raise RuntimeError(fit.reasons)
    print(json.dumps(asdict(fit.record),indent=2))
    tank=ReservoirTank(date(2020,1,1),30e6,1e6,.2,5e-10,4e-10,0,raw)
    history=history_from_frame(pd.read_csv(folder/'history.csv'))
    for name,model in (('raw',raw),('matched',fit.model)):
        problems=failures(physical_qc(model,lab))
        if problems:
            raise RuntimeError(problems)
        result=simulate(replace(tank,pvt_model=model),history)
        print(name,'QC failures:',problems,'MBE converged:',result.converged,
              'last pressure MPa:',result.states[-1].pressure/1e6)
        if not result.converged:
            raise RuntimeError(result.failed_timestep)


if __name__=='__main__':
    main()
