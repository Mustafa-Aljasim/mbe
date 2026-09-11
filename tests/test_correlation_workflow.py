from dataclasses import replace, asdict
from datetime import date
from io import BytesIO
import numpy as np
import pandas as pd
import pytest
from material_balance_studio.domain.models import CumulativeVolumes, HistoryRecord, ReservoirTank
from material_balance_studio.domain.validation import EngineeringValidationError, PVTOutOfRangeError
from material_balance_studio.io.history import read_table, pvt_from_frame
from material_balance_studio.pvt.fluid import BlackOilFluid, fluid_from_inputs
from material_balance_studio.pvt.laboratory import LaboratoryData, LabPoint, lab_from_frame, lab_to_frame, PROPERTIES
from material_balance_studio.pvt.models.correlation_model import CorrelationPVTModel, PropertySelection, Transform
from material_balance_studio.pvt.correlations.registry import REGISTRY
from material_balance_studio.pvt.correlations.base import Correlation
from material_balance_studio.pvt.matching.metrics import error_metrics
from material_balance_studio.pvt.matching.ranking import screen_property, recommendation
from material_balance_studio.pvt.matching.regression import fit_property
from material_balance_studio.pvt.qc.physical_checks import physical_qc, failures
from material_balance_studio.pvt.qc.validity import applicability
from material_balance_studio.pvt.correlations.gas_fvf import real_gas_fvf
from material_balance_studio.solver.simulation import simulate
from material_balance_studio.units.conversions import from_si
from material_balance_studio.units.temperature import temperature_from_kelvin


@pytest.fixture
def model():
    return CorrelationPVTModel(fluid_from_inputs(api=35, gas_gravity=.75, temperature=90, initial_pressure=30e6,
                                                pb=18e6, rsb=90), (2e6,35e6))


def synthetic(model, property, a=0, b=1, pressures=(4e6,8e6,12e6,18e6,24e6,30e6)):
    return LaboratoryData(tuple(LabPoint(p, **{property: a+b*model.property_at_pressure(property,p)}) for p in pressures))


def test_field_and_si_fluid_models_are_equivalent(model):
    f = model.fluid
    field = fluid_from_inputs(api=f.api, gas_gravity=f.gas_gravity, temperature=temperature_from_kelvin(f.temperature,"FIELD"),
        initial_pressure=from_si(f.initial_pressure,"pressure","FIELD"), pb=from_si(f.pb,"pressure","FIELD"),
        rsb=from_si(f.rsb,"rs","FIELD"), units="FIELD")
    other = replace(model, fluid=field)
    for p in (3e6, 18e6, 30e6):
        for name in PROPERTIES:
            assert other.property_at_pressure(name,p) == pytest.approx(model.property_at_pressure(name,p), rel=1e-12)


@pytest.mark.parametrize("family", ["standing", "vasquez_beggs", "glaso"])
def test_pb_branches_and_continuity(model, family):
    m = replace(model, selection=replace(model.selection, pb=f"{family}_pb", rs=f"{family}_rs", bo=f"{family}_bo"))
    p = m.pb
    below, at, above = [m.properties_at_pressure(v) for v in (p*.9, p, p*1.1)]
    assert below.rs < at.rs == above.rs == m.rsb
    assert below.bo < at.bo > above.bo
    assert below.oil_viscosity > at.oil_viscosity < above.oil_viscosity
    for name in ("rs", "bo", "oil_viscosity"):
        assert m.property_at_pressure(name,p*(1-1e-8)) == pytest.approx(m.property_at_pressure(name,p*(1+1e-8)),rel=1e-6)
    assert not failures(physical_qc(m))


@pytest.mark.parametrize("pb,rsb", [(None,90), (18e6,None)])
def test_one_measured_anchor_is_sufficient(model, pb, rsb):
    m = replace(model, fluid=replace(model.fluid,pb=pb,rsb=rsb))
    assert m.pb > 0 and m.rsb > 0
    assert m.property_at_pressure("rs", 30e6) == m.rsb


def test_missing_both_anchors_rejected(model):
    with pytest.raises(EngineeringValidationError, match="Provide measured"):
        replace(model, fluid=replace(model.fluid,pb=None,rsb=None))


@pytest.mark.parametrize("p", [0., -1., 1e6, 36e6, float("nan"), float("inf")])
def test_model_never_extrapolates(model,p):
    with pytest.raises(PVTOutOfRangeError):
        model.properties_at_pressure(p)


@pytest.mark.parametrize("property", ["bo","rs"])
def test_known_correlation_ranks_first_and_keeps_all_results(model, property):
    lab = synthetic(model, property)
    results = screen_property(model,lab,property)
    assert len(results) == 3
    assert recommendation(results) == getattr(model.selection,property)
    assert results[0].metrics.rmse < 1e-12
    for result in results:
        assert result.pressures == tuple(p.pressure for p in lab.points)
        assert len(result.predicted) == 6 and result.metrics.n == 6
    assert model.selection == PropertySelection()


def test_ranking_prefers_validity_tier_over_perfect_extrapolated_fit(model):
    vb = model.with_selection("bo", "vasquez_beggs_bo")
    results = screen_property(model,synthetic(vb,"bo"),"bo")
    winner_by_error = min(results,key=lambda r:r.metrics.rmse)
    assert winner_by_error.correlation == "vasquez_beggs_bo"
    assert winner_by_error.applicability.status == "CAUTION"
    assert results[0].applicability.status == "VALID"


def test_pb_screening_compares_unanchored_predictions(model):
    results = screen_property(model,LaboratoryData(),"pb")
    assert len({r.predicted for r in results}) == 3
    assert any(r.metrics.rmse > 0 for r in results)
    no_rsb = replace(model,fluid=replace(model.fluid,rsb=None))
    assert recommendation(screen_property(no_rsb,LaboratoryData(),"pb")) is None


def test_no_measurements_no_recommendation(model):
    assert recommendation(screen_property(model,LaboratoryData(),"bo")) is None


@pytest.mark.parametrize("property,a,b", [("bo",.04,1.02),("bw",.02,.98),("oil_viscosity",.0001,1.05),("z",.01,1.01),("bg",0,1.03)])
def test_matching_improves_sparse_property_and_preserves_raw(model,property,a,b):
    before = asdict(model.properties_at_pressure(10e6))
    lab = synthetic(model,property,a,b)
    attempt = fit_property(model,lab,property)
    assert attempt.accepted, attempt.reasons
    assert attempt.record.a == pytest.approx(a,abs=1e-12)
    assert attempt.record.b == pytest.approx(b,rel=1e-10)
    assert attempt.record.after.rmse < 1e-12
    assert attempt.record.before.rmse > attempt.record.after.rmse
    assert not failures(physical_qc(attempt.model,lab))
    assert asdict(model.properties_at_pressure(10e6)) == before
    if property in ("bg","z"):
        for p in (2e6,18e6,35e6):
            v=attempt.model.properties_at_pressure(p)
            assert v.bg == pytest.approx(real_gas_fvf(p,model.fluid.temperature,v.z))


def test_rs_matching_retains_rsb_plateau(model):
    lab=synthetic(model,"rs",model.rsb*.1,.9)
    result=fit_property(model,lab,"rs")
    assert result.accepted, result.reasons
    assert result.model.property_at_pressure("rs",30e6) == pytest.approx(model.rsb)
    assert result.record.after.rmse < 1e-12


def test_pb_scalar_scale_matching(model):
    m=replace(model,use_measured_pb=False)
    result=fit_property(m,LaboratoryData(),"pb")
    assert result.accepted, result.reasons
    assert result.model.pb == pytest.approx(model.fluid.pb)
    assert not fit_property(model,LaboratoryData(),"pb").accepted


def test_refitting_property_does_not_compound_transforms(model):
    lab=synthetic(model,"bo",.04,1.02)
    one=fit_property(model,lab,"bo")
    two=fit_property(one.model,lab,"bo")
    assert two.accepted and two.record.a == one.record.a and two.record.b == one.record.b


def test_negative_slope_match_is_rejected(model):
    result=fit_property(model,synthetic(model,"bo",3.,-1.),"bo")
    assert not result.accepted and result.model is None
    assert "positive B" in " ".join(result.reasons)


def test_positive_lab_but_negative_extrapolation_rejected(model):
    # All lab values positive near Pb, while this regression goes negative at 2 MPa.
    lab=synthetic(model,"bo",-12.,10.,pressures=(16e6,18e6,20e6))
    result=fit_property(model,lab,"bo")
    assert not result.accepted and result.model is None
    assert "bo" in " ".join(result.reasons)


def test_rs_all_above_pb_cannot_identify_slope(model):
    result=fit_property(model,synthetic(model,"rs",pressures=(20e6,25e6,30e6)),"rs")
    assert not result.accepted and "below Pb" in " ".join(result.reasons)


def test_affine_single_point_rejected(model):
    assert not fit_property(model,synthetic(model,"bo",pressures=(18e6,)),"bo").accepted


def test_invalid_transform_constraints(model):
    with pytest.raises(EngineeringValidationError, match="preserve"):
        model.with_transform(Transform("rs",0,2))
    with pytest.raises(EngineeringValidationError,match="scale-only"):
        model.with_transform(Transform("bg",1,1))
    with pytest.raises(EngineeringValidationError,match="simultaneously"):
        replace(model,transforms=(Transform("bg",0,1),Transform("z",0,1)))


def test_metrics_zero_measurements_are_excluded_only_from_percentages():
    m=error_metrics([0,2,4],[1,3,2])
    assert m.rmse == pytest.approx(2**.5) and m.mae == pytest.approx(4/3)
    assert m.mape == 50 and m.max_ape == 50 and m.bias == 0 and m.percentage_n == 2
    assert error_metrics([0],[1]).mape is None


@pytest.mark.parametrize("frame,match", [
    ({"pressure":[1e6,1e6],"bo":[1.1,1.2]},"Duplicate"),
    ({"pressure":[None],"bo":[1.1]},"missing"),
    ({"pressure":[0],"bo":[1.1]},"positive"),
    ({"pressure":[1e6],"bo":[-1]},"positive"),
    ({"pressure":[1e6],"z":[float('inf')]},"finite"),
    ({"pressure":[1e6],"rs":[-1]},"nonnegative"),
    ({"pressure_MPa":[1],"bo":[1.1]},"headers"),
    ({"pressure":[1e6],"bo":["bad"]},"invalid"),
])
def test_bad_lab_data_is_rejected(frame,match):
    with pytest.raises(EngineeringValidationError,match=match):
        lab_from_frame(pd.DataFrame(frame))


def test_sparse_csv_excel_roundtrip_and_units():
    frame=pd.DataFrame({"pressure":[1e6,2e6,3e6],"bo":[1.1,None,1.2],"rs":[None,0,None]})
    lab=lab_from_frame(frame)
    assert len(lab.measurements("bo"))==2 and lab.points[1].bo is None and lab.points[1].rs==0
    assert all(p.bg is None for p in lab.points)
    field=lab_from_frame(lab_to_frame(lab,"FIELD"),"FIELD")
    assert field.points[0].pressure==pytest.approx(1e6)
    for extension in ("csv","xlsx"):
        buffer=BytesIO()
        if extension=="csv":
            buffer.write(frame.to_csv(index=False).encode())
        else:
            frame.to_excel(buffer,index=False)
        buffer.seek(0)
        assert lab_from_frame(read_table(buffer,"lab."+extension))==lab
    with pytest.raises(EngineeringValidationError):
        pvt_from_frame(frame)  # existing full-table contract stays strict


def test_blank_editor_rows_ignored_and_suspicious_data_flagged():
    lab=lab_from_frame(pd.DataFrame({"pressure":[1e6,None],"bo":[8,None]}))
    assert len(lab.points)==1 and any("Suspicious" in w for w in lab.warnings)


def test_qc_reports_lab_extrapolation_and_outside_range(model):
    qc=physical_qc(model,synthetic(model,"bo"))
    assert any(row.check=="Extrapolation: bo" and row.status=="CAUTION" for row in qc)
    assert failures(physical_qc(model,LaboratoryData((LabPoint(40e6,bo=1.2),))))


def test_applicability_caution_and_not_applicable(model):
    hot=replace(model,fluid=replace(model.fluid,temperature=450))
    assert applicability(REGISTRY["standing_bo"],hot).status=="CAUTION"
    sour=replace(model,fluid=replace(model.fluid,gas_basis="sour"))
    assert applicability(REGISTRY["sutton_dak_z"],sour).status=="NOT_APPLICABLE"
    assert failures(physical_qc(sour))


def test_correlation_model_runs_real_mbe_depletion_and_injection(model):
    tank=ReservoirTank(date(2020,1,1),30e6,1e6,.2,5e-10,4e-10,.1,model)
    history=[HistoryRecord(date(2020,month,1),CumulativeVolumes(np=n,gp=n*90,wp=n*.1,winj=n*.05,ginj=n*2))
             for month,n in ((2,1000),(3,10000),(4,40000))]
    result=simulate(tank,history)
    assert result.converged
    assert all(s.balance.relative_residual<=1e-8 and s.balance.absolute_residual<=1e-5 for s in result.states)
    assert result.states[-1].pressure < result.states[0].pressure < 30e6


def test_identical_property_functions_give_identical_mbe_results(model,pvt,tank,monkeypatch):
    # Inject reference equations into the registry, not into the MBE or model's
    # pressure dispatch. The real CorrelationPVTModel evaluates them continuously.
    selection=vars(model.selection).copy()
    for name in ("rs","bo","bg","bw"):
        key="test_table_"+name
        entry=Correlation(key,name,"Reference table equation",
            lambda i,n=name:getattr(pvt.properties_at_pressure(i.pressure),n),"test reference","identical tabular function")
        monkeypatch.setitem(REGISTRY,key,entry)
        selection[name]=key
    corr=replace(model,fluid=replace(model.fluid,pb=40e6,rsb=80),pressure_bounds=(10e6,40e6),selection=PropertySelection(**selection))
    for p in np.linspace(10e6,40e6,31):
        one,two=pvt.properties_at_pressure(float(p)),corr.properties_at_pressure(float(p))
        assert (one.bo,one.rs,one.bg,one.bw)==(two.bo,two.rs,two.bg,two.bw)
    history=[HistoryRecord(date(2020,month,1),CumulativeVolumes(np=n,gp=80*n,wp=100,winj=50,ginj=1000))
             for month,n in ((2,5000),(3,20000),(4,50000))]
    original=simulate(tank,history)
    alternate=simulate(replace(tank,pvt_model=corr),history)
    assert original.converged and alternate.converged
    for a,b in zip(original.states,alternate.states):
        assert a.pressure==pytest.approx(b.pressure,abs=1e-6)
        for key,value in asdict(a.balance).items():
            if isinstance(value,(int,float)):
                assert value==pytest.approx(getattr(b.balance,key),abs=1e-7)
