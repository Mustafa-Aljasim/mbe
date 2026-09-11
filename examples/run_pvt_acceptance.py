"""Write reproducible Phase 2.1 numerical acceptance evidence (canonical SI)."""
from collections import Counter
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path

from material_balance_studio.pvt.fluid import fluid_from_inputs
from material_balance_studio.pvt.models.correlation_model import CorrelationPVTModel, Transform
from material_balance_studio.pvt.laboratory import LaboratoryData, LabPoint
from material_balance_studio.pvt.matching.regression import fit_property
from material_balance_studio.pvt.qc.physical_checks import physical_qc, pressure_grid
from pvt_acceptance_case import comparison


def acceptance():
    root = Path(__file__).resolve().parents[1]
    raw = CorrelationPVTModel(fluid_from_inputs(api=35, gas_gravity=.75, temperature=90,
        initial_pressure=30e6, pb=18e6, rsb=90), (2e6, 35e6))
    lab = LaboratoryData(tuple(LabPoint(p, bo=.04+1.02*raw.property_at_pressure("bo", p))
                              for p in (4e6, 8e6, 12e6, 18e6, 24e6, 30e6)))
    fit = fit_property(raw, lab, "bo")
    assert fit.accepted
    qc = {}
    for label, model in (("raw", raw), ("matched", fit.model)):
        items = physical_qc(model, lab)
        qc[label] = {"pressure_count": len(pressure_grid(model, lab)),
                     "counts": dict(Counter(item.status for item in items)),
                     "checks": [asdict(item) for item in items]}
    probes = {}
    for family in ("standing", "vasquez_beggs", "glaso"):
        base = replace(raw, selection=replace(raw.selection, rs=family+"_rs", bo=family+"_bo"))
        for matched in (False, True):
            model = replace(base, transforms=(Transform("rs", .1*base.rsb, .9),
                Transform("bo", .04, 1.02), Transform("oil_viscosity", .0001, 1.1))) if matched else base
            rows = []
            for factor in (.99, .999, 1, 1.001, 1.01):
                p = factor*model.pb
                row = {"p_over_pb": factor, "pressure_Pa": p}
                row.update({name: model.property_at_pressure(name, p) for name in ("rs", "bo", "oil_viscosity")})
                row["co"] = model.property_at_pressure("co", p) if p >= model.pb else None
                rows.append(row)
            probes[family+("_matched" if matched else "_raw")] = rows
    bad_lab = LaboratoryData(tuple(LabPoint(p, bo=-12+10*raw.property_at_pressure("bo", p))
                                  for p in (16e6, 18e6, 20e6)))
    rejected = fit_property(raw, bad_lab, "bo")
    assert not rejected.accepted
    original = json.loads((root/"docs"/"phase_2_core_hashes.json").read_text())
    core = []
    for entry in original["files"]:
        path = Path(entry["Path"])
        digest = hashlib.sha256(path.read_bytes()).hexdigest().upper()
        core.append({"path": str(path.relative_to(root)), "sha256": digest, "unchanged": digest == entry["Hash"]})
    assert all(item["unchanged"] for item in core)
    return {"fixture": "Synthetic 35 API, gas SG .75, 90 C, Pi 30 MPa, Pb 18 MPa, Rsb 90 m3/m3; range 2–35 MPa",
            "units": "Canonical SI; co=null means not applicable below Pb",
            "qc": qc, "bubble_point_probes": probes, "accepted_bo_fit": asdict(fit.record),
            "rejected_bo_fit": {"record": asdict(rejected.record), "reasons": rejected.reasons,
                                "candidate_bo_at_2_MPa": -12+10*raw.property_at_pressure("bo", 2e6)},
            "mbe_equivalence": comparison(), "core_files": core, "all_core_files_unchanged": True}


if __name__ == "__main__":
    output = Path(__file__).resolve().parents[1]/"docs"/"phase_2_1_numerical_results.json"
    output.write_text(json.dumps(acceptance(), indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(output)
