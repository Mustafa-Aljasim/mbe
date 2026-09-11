"""Phase 2.1 acceptance; the existing suite retains one corrected Sutton reference."""
from dataclasses import replace
import numpy as np
import pandas as pd
import pytest
from material_balance_studio.pvt.fluid import fluid_from_inputs
from material_balance_studio.pvt.models.correlation_model import CorrelationPVTModel, Transform
from material_balance_studio.pvt.laboratory import LaboratoryData, LabPoint
from material_balance_studio.pvt.qc.physical_checks import pressure_grid, physical_qc, failures, property_qc_status
from material_balance_studio.pvt.matching.regression import fit_property
from material_balance_studio.pvt.correlations.registry import REGISTRY
from material_balance_studio.pvt.correlations.gas_pseudocritical import PSEUDO_CRITICAL_METHODS, PseudoCriticalMethod, PseudoCriticalProperties
from material_balance_studio.pvt.correlations.oil_fvf import standing_bo
from material_balance_studio.presentation.pvt import selected_model_frame, property_detail_frame, comparison_figure, fit_frame
from material_balance_studio.domain.validation import EngineeringValidationError


@pytest.fixture
def model():
    return CorrelationPVTModel(fluid_from_inputs(api=35,gas_gravity=.75,temperature=90,
        initial_pressure=30e6,pb=18e6,rsb=90),(2e6,35e6))


@pytest.mark.parametrize("family",["standing","glaso","vasquez_beggs"])
@pytest.mark.parametrize("matched",[False,True])
def test_dedicated_five_pressure_pb_acceptance(model,family,matched):
    model=replace(model,selection=replace(model.selection,rs=f"{family}_rs",bo=f"{family}_bo"))
    if matched:
        model=replace(model,transforms=(Transform("rs",.1*model.rsb,.9),Transform("bo",.04,1.02),Transform("oil_viscosity",.0001,1.1)))
    factors=np.array([.99,.999,1,1.001,1.01])
    ps=factors*model.pb
    rs=np.array([model.property_at_pressure("rs",float(p)) for p in ps])
    bo=np.array([model.property_at_pressure("bo",float(p)) for p in ps])
    mu=np.array([model.property_at_pressure("oil_viscosity",float(p)) for p in ps])
    assert np.all(np.diff(rs[:3])>0)
    assert np.all(rs[2:]==pytest.approx(model.rsb))
    assert bo[0]<bo[1]<bo[2]>bo[3]>bo[4]
    assert mu[0]>mu[1]>mu[2]<mu[3]<mu[4]
    # Finite changes shrink towards the boundary; test both sides independently.
    for v in (rs,bo,mu):
        assert abs(v[1]-v[2]) <= .12*abs(v[0]-v[2])+1e-12
        assert abs(v[3]-v[2]) <= .12*abs(v[4]-v[2])+1e-12
        assert max(abs(v[1]/v[2]-1),abs(v[3]/v[2]-1))<.01
    co=[model.property_at_pressure("co",float(p)) for p in ps[2:]]
    assert all(c>0 for c in co) and abs(co[1]/co[0]-1)<.002
    assert not failures(physical_qc(model))


@pytest.mark.parametrize("factor",[.99,.999])
def test_undersaturated_co_not_reported_on_saturated_branch(model,factor):
    with pytest.raises(EngineeringValidationError,match="not applicable below Pb"):
        model.property_at_pressure("co",factor*model.pb)


@pytest.mark.parametrize("matched",[False,True])
def test_co_is_derivative_of_actual_undersaturated_bo(model,matched):
    if matched:
        model=model.with_transform(Transform("bo",.04,1.02))
    for p in (model.pb,1.001*model.pb,1.01*model.pb,30e6):
        h=1.
        # One-sided at Pb; symmetric above. co below Pb is deliberately undefined.
        left=p if p==model.pb else p-h
        right=p+h
        finite_difference=-(np.log(model.property_at_pressure("bo",right))-np.log(model.property_at_pressure("bo",left)))/(right-left)
        assert model.property_at_pressure("co",p)==pytest.approx(finite_difference,rel=2e-6)


def test_pseudocritical_method_is_independent_and_affects_gas_only(model,monkeypatch):
    reference=model.pseudo_critical_properties
    method=PseudoCriticalMethod("test_pc","Test pseudo-critical method",lambda sg: PseudoCriticalProperties(reference.pressure*1.1,reference.temperature),
                                "test-only",(.5,1.5))
    monkeypatch.setitem(PSEUDO_CRITICAL_METHODS,method.key,method)
    other=replace(model,pseudo_critical_method=method.key)
    assert other.selection.z==model.selection.z
    assert other.property_at_pressure("z",20e6)!=model.property_at_pressure("z",20e6)
    assert other.property_at_pressure("bg",20e6)!=model.property_at_pressure("bg",20e6)
    for name in ("bo","rs","oil_viscosity","bw"):
        assert other.property_at_pressure(name,20e6)==model.property_at_pressure(name,20e6)
    assert REGISTRY[model.selection.z].name=="Dranchuk–Abou-Kassem"


def test_z_method_can_change_without_changing_pseudocritical(model,monkeypatch):
    entry=replace(REGISTRY[model.selection.z],key="test_z",name="Test z",evaluate=lambda inputs:1.)
    monkeypatch.setitem(REGISTRY,entry.key,entry)
    other=model.with_selection("z",entry.key)
    assert other.pseudo_critical_method==model.pseudo_critical_method
    assert other.pseudo_critical_properties==model.pseudo_critical_properties
    assert other.property_at_pressure("z",20e6)==1


def test_dense_grid_covers_requested_range_and_pb_probe_pressures(model):
    grid=pressure_grid(model)
    assert len(grid)>=1001 and grid[0]==2e6 and grid[-1]==35e6
    assert np.max(np.diff(grid)) <= (35e6-2e6)/1000+1e-8
    assert all(factor*model.pb in grid for factor in (.99,.999,1,1.001,1.01))


def test_rejected_match_reports_better_lab_fit_and_bad_extrapolation(model):
    ps=(16e6,18e6,20e6)
    lab=LaboratoryData(tuple(LabPoint(p,bo=-12+10*model.property_at_pressure("bo",p)) for p in ps))
    attempt=fit_property(model,lab,"bo")
    assert not attempt.accepted and attempt.model is None
    record=attempt.record
    assert record.a==pytest.approx(-12) and record.b==pytest.approx(10)
    assert record.after.rmse < 1e-10 < record.before.rmse
    assert record.qc_before in ("PASS","CAUTION") and record.qc_after=="FAIL"
    assert record.maximum_full_range_deviation>0
    assert any("Nonphysical bo" in reason for reason in attempt.reasons)
    assert -12+10*model.property_at_pressure("bo",2e6)<0


def test_accepted_match_reports_full_range_change_at_unmeasured_pressure(model):
    lab=LaboratoryData(tuple(LabPoint(p,bo=.04+1.02*model.property_at_pressure("bo",p)) for p in (4e6,8e6,12e6,24e6,30e6)))
    attempt=fit_property(model,lab,"bo")
    assert attempt.accepted
    record=attempt.record
    assert record.maximum_deviation_pressure==model.pb  # absent from lab points
    assert record.maximum_full_range_deviation==pytest.approx(.04+.02*model.property_at_pressure("bo",model.pb))
    assert record.qc_after in ("PASS","CAUTION")
    for p in (2e6,6e6,17e6,19e6,35e6):
        assert attempt.model.property_at_pressure("bo",p)>0


def test_negative_slope_attempt_retains_reviewable_coefficients(model):
    lab=LaboratoryData(tuple(LabPoint(p,bo=3-model.property_at_pressure("bo",p)) for p in (4e6,10e6,18e6)))
    result=fit_property(model,lab,"bo")
    assert not result.accepted and result.record.b==pytest.approx(-1) and result.record.qc_after=="FAIL"


def test_selected_summary_separates_gas_methods_and_shows_property_qc(model):
    qc=physical_qc(model)
    summary=selected_model_frame(model,qc).set_index("Property")
    assert summary.loc["Gas pseudo-critical","Method"]=="Sutton"
    assert summary.loc["Gas z-factor","Method"]=="Dranchuk–Abou-Kassem"
    assert set(summary.QC)<={"PASS","CAUTION","FAIL"}
    assert summary.loc["Bubble point","Treatment"]=="Measured anchor"


@pytest.mark.parametrize("units",["SI","FIELD"])
def test_matched_summary_reports_pb_fit_coefficients_errors_and_units(model,units):
    lab=LaboratoryData(tuple(LabPoint(p,bo=.04+1.02*model.property_at_pressure("bo",p)) for p in (4e6,18e6,30e6)))
    fit=fit_property(model,lab,"bo")
    details=property_detail_frame(fit.model,"bo",fit.record,units).set_index("Item")
    assert details.loc["Matched","Value"]=="Yes"
    assert {"A","B","RMSE before","RMSE after","Maximum full-range change"}<=set(details.index)
    pb=property_detail_frame(model,"pb",None,units).set_index("Item")
    assert {"Measured Pb","Calculated Pb (untuned)","Active Pb"}<=set(pb.index)
    assert ("MPa" if units=="SI" else "psia") in pb.loc["Active Pb","Value"]
    table=fit_frame([fit.record],units)
    assert "Maximum full-range change" in table and "QC after" in table


def test_property_summary_carries_failures_from_dependencies(model):
    invalid=model.with_transform(Transform("bo",-1.15,1.))
    qc=physical_qc(invalid)
    assert property_qc_status(qc,"bo")=="FAIL"
    assert property_qc_status(qc,"bw")!="FAIL"


@pytest.mark.parametrize("property",["bo","rs","oil_viscosity","z","bg","bw"])
@pytest.mark.parametrize("units",["SI","FIELD"])
def test_engineering_plot_lab_raw_match_and_units(model,property,units):
    lab=LaboratoryData(tuple(LabPoint(p,**{property:model.property_at_pressure(property,p)}) for p in (4e6,18e6,30e6)))
    figure=comparison_figure(model,model,lab,property,units)
    assert [trace.name for trace in figure.data]==["Raw correlation","Matched correlation","Laboratory"]
    assert "vs Pressure" in figure.layout.title.text
    assert ("MPa" if units=="SI" else "psia") in figure.layout.xaxis.title.text
    assert len(figure.data[0].x)>=1001 and len(figure.data[2].x)==3
    assert figure.layout.shapes[0].type=="line"


def test_published_standing_bo_worked_example():
    # T.A. Blasingame, Properties of Reservoir Fluids appendix, A-50,
    # equations A-37/A-38: Rs 569.54 scf/STB, gas SG .786, oil SG .8217,
    # T 220°F -> Bo 1.3687 rb/STB (rounded published worked solution).
    assert standing_bo(569.54,220,.8217,.786)==pytest.approx(1.3687,abs=5e-5)


def test_mbe_equivalence_with_nonzero_gas_cap_and_rock_water_terms():
    import runpy
    from pathlib import Path
    example=runpy.run_path(str(Path(__file__).resolve().parents[1]/"examples"/"pvt_acceptance_case.py"))
    result=example["comparison"]()
    assert result["maximum_property_difference"]==0 and result["both_converged"]
    for step in result["steps"]:
        assert step["table"]["mEg"]>0 and step["table"]["Efw"]>0 and step["table"]["Eo"]>0
        for name,tolerance in result["tolerances"].items():
            assert step["absolute_difference"][name]<=tolerance
    assert not any(key.startswith("acceptance_reference_") for key in REGISTRY)
