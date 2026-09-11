"""Independent Phase 3A benchmarks, chronological state and solver safety."""
from dataclasses import asdict, FrozenInstanceError, replace
from datetime import date
from math import log
from pathlib import Path
import hashlib
import json
import pytest

from material_balance_studio.aquifer import (NoAquifer, PotAquifer, SchilthuisAquifer, FetkovichAquifer,
                                             CarterTracyAquifer, VanEverdingenHurstAquifer)
from material_balance_studio.aquifer.base import AquiferContext
from material_balance_studio.aquifer.veh_response import veh_infinite_water_influx_dimensionless
from material_balance_studio.aquifer.registry import MODELS
from material_balance_studio.domain.models import CumulativeVolumes, HistoryRecord, ReservoirTank, SolverSettings
from material_balance_studio.io.history import read_table, pvt_from_frame, history_from_frame
from material_balance_studio.mbe.oil_balance import evaluate_balance
from material_balance_studio.pvt.table_model import TablePVTModel
from material_balance_studio.solver.pressure_solver import solve_pressure
from material_balance_studio.solver.simulation import simulate
from material_balance_studio.units.display import from_display, to_display
from material_balance_studio.diagnostics.balance_closure import inspect_state
from material_balance_studio.presentation.tables import equation_inspector_frame, equation_sides_frame

ROOT = Path(__file__).resolve().parents[1]
DAY = 86400.
MODELS_UNDER_TEST = [NoAquifer(), PotAquifer(.0025), SchilthuisAquifer(100/(DAY*1e6)),
                     FetkovichAquifer(1e6, 1e-9, .001*log(2)/DAY),
                     CarterTracyAquifer(1000,10,10,.2,.2*.001*1e-9*1000**2/DAY,.001,1e-9),
                     VanEverdingenHurstAquifer(1000,10,10,.2,.2*.001*1e-9*1000**2/DAY,.001,1e-9)]


def history():
    return [HistoryRecord(date(2020, m, 1), CumulativeVolumes(np=n, gp=80*n))
            for m, n in ((2, 10000), (3, 20000), (4, 30000))]


def test_exact_no_aquifer_pre_phase_3_baseline():
    baseline = json.loads((ROOT/"docs/phase_3a_baseline.json").read_text())
    tank = ReservoirTank(date(2020,1,1),30e6,1e6,.2,5e-10,4e-10,0,
        TablePVTModel(pvt_from_frame(read_table(ROOT/"examples/pvt.csv"))), aquifer=NoAquifer())
    result = simulate(tank, history_from_frame(read_table(ROOT/"examples/history.csv")))
    actual = json.loads(json.dumps(asdict(result), default=str))
    def compare(original, current):
        if isinstance(original, dict):
            for key, value in original.items():
                compare(value, current[key])
        elif isinstance(original, list):
            assert len(original) == len(current)
            for a,b in zip(original,current):
                compare(a,b)
        else:
            assert original == current  # exact; no numerical tolerance
    compare(baseline["result"], actual)
    assert all(s.balance.aquifer_support == 0 for s in result.states)


def test_pvt_and_original_expansion_withdrawal_files_unchanged():
    baseline = json.loads((ROOT/"docs/phase_3a_baseline.json").read_text())
    for path, expected in baseline["hashes"].items():
        if not path.endswith("mbe/oil_balance.py"):
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == expected


def test_pot_hand_benchmark_and_reversal():
    model = PotAquifer(.0025)  # 2500 reservoir m³/MPa
    state = model.initial_state(30e6)
    for pressure, expected, increment in ((28e6,5000,5000),(25e6,12500,7500),(27e6,7500,-5000)):
        result = model.compute_step(state, state.previous_reservoir_pressure, pressure, 10*DAY)
        assert result.cumulative_influx == expected
        assert result.incremental_influx == increment
        assert result.updated_state.aquifer_pressure == pressure
        assert bool(result.warnings) == (increment < 0)
        state = result.updated_state


def test_schilthuis_hand_trapezoidal_benchmark():
    model = SchilthuisAquifer(100/(DAY*1e6))  # 100 reservoir m³/(day MPa)
    state = model.initial_state(30e6)
    for pressure, days, increment, cumulative, endpoint in (
            (28e6,10,1000,1000,200),(26e6,20,6000,7000,400),(27e6,5,1750,8750,300)):
        result = model.compute_step(state, state.previous_reservoir_pressure, pressure, days*DAY)
        assert result.incremental_influx == pytest.approx(increment)
        assert result.cumulative_influx == pytest.approx(cumulative)
        assert result.average_influx_rate*DAY == pytest.approx(increment/days)
        assert result.endpoint_influx_rate*DAY == pytest.approx(endpoint)
        assert result.updated_state.aquifer_pressure == 30e6
        state = result.updated_state
    assert state.elapsed_time == 35*DAY


def test_fetkovich_independent_half_life_multistep_benchmark():
    # C=1000 m³/MPa, J=C ln(2)/day: exact exponential fraction=1/2 per day.
    # First boundary average=(30+26)/2=28 MPa; then boundary holds at26 MPa.
    model = FetkovichAquifer(1e6,1e-9,.001*log(2)/DAY)
    state = model.initial_state(30e6)
    increments = []
    for increment, cumulative, pressure in ((1000,1000,29e6),(1500,2500,27.5e6),
                                            (750,3250,26.75e6),(375,3625,26.375e6)):
        result = model.compute_step(state,state.previous_reservoir_pressure,26e6,DAY)
        assert result.incremental_influx == pytest.approx(increment, abs=1e-10)
        assert result.cumulative_influx == pytest.approx(cumulative, abs=1e-10)
        assert result.updated_state.aquifer_pressure == pytest.approx(pressure, abs=1e-8)
        increments.append(result.incremental_influx)
        state = result.updated_state
    assert increments[1] > increments[2] > increments[3]
    assert state.elapsed_time == 4*DAY


@pytest.mark.parametrize("model", MODELS_UNDER_TEST)
def test_candidate_state_is_immutable_order_independent_and_repeatable(tank, model):
    state = model.initial_state(30e6)
    before = asdict(state)
    context = AquiferContext(model,state,30e6,DAY)
    first = context.evaluate(25e6)
    context.evaluate(20e6)
    assert context.evaluate(25e6) == first
    assert asdict(state) == before
    with pytest.raises(FrozenInstanceError):
        state.cumulative_influx = 10
    volumes = CumulativeVolumes(np=10000,gp=800000)
    def residual(p):
        return evaluate_balance(p,tank,volumes,aquifer_influx=context.evaluate(p).cumulative_influx).residual
    a = residual(25e6)
    residual(20e6)
    assert residual(25e6) == a and asdict(state) == before


@pytest.mark.parametrize("model", MODELS_UNDER_TEST)
def test_chronological_commit_and_cumulative_consistency(tank, model):
    result = simulate(replace(tank,aquifer=model),history())
    assert result.converged
    previous = result.initial_state
    for current in result.states:
        step = current.aquifer_step
        assert step.updated_state == current.aquifer_state
        assert step.cumulative_influx == pytest.approx(previous.aquifer_state.cumulative_influx+step.incremental_influx)
        assert current.aquifer_state.previous_reservoir_pressure == current.pressure
        assert current.aquifer_state.elapsed_time == (current.date-tank.initial_date).days*DAY
        assert current.balance.aquifer_support == step.cumulative_influx
        assert current.balance.residual == current.balance.withdrawal.net-current.balance.total_expansion_support-step.cumulative_influx
        assert current.balance.relative_residual <= 1e-8
        previous = current


@pytest.mark.parametrize("model", MODELS_UNDER_TEST)
def test_initial_date_zero_row_does_not_advance_aquifer(tank, model):
    result = simulate(replace(tank,aquifer=model),[HistoryRecord(tank.initial_date),*history()])
    assert result.converged
    assert result.states[0].aquifer_state == result.initial_state.aquifer_state
    assert result.states[0].aquifer_step is None


@pytest.mark.parametrize("weak,strong", [
    (PotAquifer(.001),PotAquifer(.01)),
    (SchilthuisAquifer(1e-10),SchilthuisAquifer(1e-9)),
    (FetkovichAquifer(1e6,1e-9,1e-10),FetkovichAquifer(1e7,1e-9,1e-9)),
])
def test_stronger_aquifer_support_reduces_pressure_decline(tank, weak, strong):
    runs = [simulate(replace(tank,aquifer=m),history()) for m in (NoAquifer(),weak,strong)]
    assert all(run.converged for run in runs)
    for none,w,s in zip(*(run.states for run in runs)):
        assert none.pressure < w.pressure < s.pressure < tank.initial_pressure


def test_failed_step_does_not_commit_aquifer(tank):
    model = FetkovichAquifer(1e6,1e-9,1e-10)
    active = replace(tank,aquifer=model)
    records = history()[:1]+[HistoryRecord(date(2020,3,1),CumulativeVolumes(np=950000,gp=76e6)),history()[-1]]
    # Keep cumulative history valid after the failing row; the final row must never run.
    records[-1] = HistoryRecord(date(2020,4,1),CumulativeVolumes(np=960000,gp=76.8e6))
    result = simulate(active,records)
    prefix = simulate(active,records[:1])
    assert not result.converged and len(result.states) == 1
    assert result.states == prefix.states
    assert result.failed_timestep.solution.aquifer_step is None


def test_brent_trials_do_not_depend_on_sampling_density(tank):
    model = FetkovichAquifer(1e6,1e-9,1e-9)
    state = model.initial_state(30e6)
    context = AquiferContext(model,state,30e6,31*DAY)
    solutions = [solve_pressure(replace(tank,aquifer=model),history()[0].cumulative,30e6,
                    SolverSettings(bracket_subdivisions=n),aquifer_context=context) for n in (1,8,31)]
    assert all(s.diagnostics.converged for s in solutions)
    assert [s.pressure for s in solutions] == pytest.approx([solutions[0].pressure]*3,abs=1e-5)
    assert state.cumulative_influx == 0 and state.elapsed_time == 0


@pytest.mark.parametrize("bad", [0,-1,float("nan"),float("inf")])
@pytest.mark.parametrize("factory", [lambda x:PotAquifer(x),lambda x:SchilthuisAquifer(x),
    lambda x:FetkovichAquifer(x,1e-9,1e-9),lambda x:FetkovichAquifer(1e6,x,1e-9),lambda x:FetkovichAquifer(1e6,1e-9,x)])
def test_invalid_aquifer_parameters_rejected(factory,bad):
    with pytest.raises(ValueError):
        factory(bad)


@pytest.mark.parametrize("kwargs", [
    {"inner_radius":1000,"radius_ratio":1,"thickness":10,"porosity":.2,"permeability":1e-13,"water_viscosity":.001,"total_compressibility":1e-9},
    {"inner_radius":1000,"radius_ratio":10,"thickness":10,"porosity":1,"permeability":1e-13,"water_viscosity":.001,"total_compressibility":1e-9},
    {"inner_radius":1000,"radius_ratio":10,"thickness":10,"porosity":.2,"permeability":0,"water_viscosity":.001,"total_compressibility":1e-9},
    {"inner_radius":1000,"radius_ratio":10,"thickness":10,"porosity":.2,"permeability":1e-13,"water_viscosity":0,"total_compressibility":1e-9},
    {"inner_radius":1000,"radius_ratio":10,"thickness":10,"porosity":.2,"permeability":1e-13,"water_viscosity":.001,"total_compressibility":0},
    {"inner_radius":1000,"radius_ratio":10,"thickness":10,"porosity":.2,"permeability":1e-13,"water_viscosity":.001,"total_compressibility":1e-9,"encroachment_angle":361},
])
@pytest.mark.parametrize("factory", [CarterTracyAquifer, VanEverdingenHurstAquifer])
def test_invalid_transient_aquifer_geometry_rejected(factory,kwargs):
    with pytest.raises(ValueError):
        factory(**kwargs)


@pytest.mark.parametrize("dt", [0,-1,float("nan"),float("inf")])
@pytest.mark.parametrize("model", MODELS_UNDER_TEST)
def test_nonpositive_or_nonfinite_time_rejected(model,dt):
    with pytest.raises(ValueError):
        model.compute_step(model.initial_state(30e6),30e6,25e6,dt)


def test_invalid_state_and_model_switch_are_not_reused():
    model = FetkovichAquifer(1e6,1e-9,1e-9)
    state = model.initial_state(30e6)
    with pytest.raises(ValueError,match="inconsistent"):
        model.compute_step(replace(state,cumulative_influx=1000),30e6,25e6,DAY)
    for other in (PotAquifer(.001),replace(model,productivity_index=2e-9)):
        with pytest.raises(ValueError,match="changed"):
            other.compute_step(state,30e6,25e6,DAY)
    with pytest.raises(ValueError,match="disagrees"):
        model.compute_step(state,29e6,25e6,DAY)
    with pytest.raises(ValueError):
        replace(state,aquifer_pressure=-1)
    with pytest.raises(ValueError,match="lower bound"):
        PotAquifer(.001).compute_step(PotAquifer(.001).initial_state(30e6),30e6,.5,DAY)


@pytest.mark.parametrize("model", MODELS_UNDER_TEST[1:])
def test_injection_efflux_is_signed_and_warned(tank,model):
    result = simulate(replace(tank,aquifer=model),[HistoryRecord(date(2020,2,1),CumulativeVolumes(winj=1000))])
    assert result.converged and result.states[0].pressure > tank.initial_pressure
    assert result.states[0].aquifer_step.incremental_influx < 0
    assert any("efflux" in w for w in result.states[0].warnings)


@pytest.mark.parametrize("units", ["SI","FIELD"])
def test_equation_inspector_separates_all_aquifer_terms(tank,units):
    result = simulate(replace(tank,aquifer=FetkovichAquifer(1e6,1e-9,1e-9)),history())
    state = result.states[-1]
    raw = inspect_state(state)
    assert raw["cumulative_aquifer_influx_m3"] == raw["aquifer_support_m3"]
    assert raw["previous_aquifer_pressure_pa"] > raw["aquifer_pressure_pa"]
    assert raw["total_support_m3"] == raw["total_expansion_support_m3"]+raw["aquifer_support_m3"]
    frame = equation_inspector_frame(state,units).set_index("Quantity")
    assert {"Incremental Aquifer Influx","Cumulative Aquifer Influx We","Previous Aquifer Pressure","Aquifer minus Average Reservoir Pressure"} <= set(frame.index)
    sides = equation_sides_frame(state,units)
    assert sides.iloc[1,1] == to_display(state.balance.total_support,"reservoir_volume",units)
    assert sides.iloc[1,0].endswith("+ We")


@pytest.mark.parametrize("quantity,si,field", [("aquifer_capacity",1e6,6894.757293168/.158987294928),
    ("aquifer_productivity",86400e6,86400*6894.757293168/.158987294928),
    ("aquifer_rate",86400,86400/.158987294928)])
def test_engineering_units_have_independent_scales(quantity,si,field):
    assert to_display(1,quantity,"SI") == pytest.approx(si)
    assert to_display(1,quantity,"FIELD") == pytest.approx(field)
    assert from_display(field,quantity,"FIELD") == pytest.approx(1)


@pytest.mark.parametrize("label", ["Pot","Schilthuis","Fetkovich"])
def test_field_si_aquifer_inputs_produce_equivalent_simulations(tank,label):
    kwargs = {"Pot":{"capacity":.003},"Schilthuis":{"productivity_index":1e-9},
              "Fetkovich":{"initial_water_volume":3e6,"total_compressibility":1e-9,"productivity_index":1e-9}}[label]
    quantities = {"capacity":"aquifer_capacity","productivity_index":"aquifer_productivity", "initial_water_volume":"reservoir_volume","total_compressibility":"compressibility"}
    field_kwargs = {k:from_display(to_display(v,quantities[k],"FIELD"),quantities[k],"FIELD") for k,v in kwargs.items()}
    a,b = [simulate(replace(tank,aquifer=MODELS[label](**kw)),history()) for kw in (kwargs,field_kwargs)]
    assert a.converged and b.converged
    assert [s.pressure for s in a.states] == pytest.approx([s.pressure for s in b.states],abs=1e-5)
    assert [s.balance.aquifer_support for s in a.states] == pytest.approx([s.balance.aquifer_support for s in b.states],abs=1e-7)


@pytest.mark.parametrize("label", ["Carter-Tracy","Van Everdingen-Hurst"])
def test_field_si_transient_aquifer_inputs_produce_equivalent_simulations(tank,label):
    kwargs = {"inner_radius":1800,"radius_ratio":10,"thickness":20,"porosity":.2,
              "permeability":1e-13,"water_viscosity":.001,"total_compressibility":1e-9,
              "encroachment_angle":180}
    quantities = {"inner_radius":"length","radius_ratio":"dimensionless","thickness":"length",
                  "porosity":"dimensionless","permeability":"permeability","water_viscosity":"viscosity",
                  "total_compressibility":"compressibility","encroachment_angle":"angle"}
    field_kwargs = {k:from_display(to_display(v,quantities[k],"FIELD"),quantities[k],"FIELD") for k,v in kwargs.items()}
    a,b = [simulate(replace(tank,aquifer=MODELS[label](**kw)),history()) for kw in (kwargs,field_kwargs)]
    assert a.converged and b.converged
    assert [s.pressure for s in a.states] == pytest.approx([s.pressure for s in b.states],abs=1e-5)
    assert [s.balance.aquifer_support for s in a.states] == pytest.approx([s.balance.aquifer_support for s in b.states],abs=1e-7)


@pytest.mark.parametrize("factory", [CarterTracyAquifer, VanEverdingenHurstAquifer])
def test_transient_physical_sensitivities(factory):
    base = {"inner_radius":1000,"radius_ratio":10,"thickness":10,"porosity":.2,
            "permeability":1e-13,"water_viscosity":.001,"total_compressibility":1e-9,
            "encroachment_angle":360}
    def influx(**updates):
        model = factory(**{**base,**updates})
        return model.compute_step(model.initial_state(30e6),30e6,28e6,30*DAY).cumulative_influx
    assert influx(permeability=2e-13) > influx(permeability=5e-14)
    assert influx(encroachment_angle=180) == pytest.approx(influx()*0.5)
    assert influx(water_viscosity=.002) < influx(water_viscosity=.001)
    assert influx(total_compressibility=2e-9) > 0


def test_active_aquifer_direct_solve_requires_context(tank):
    result = solve_pressure(replace(tank,aquifer=PotAquifer(.001)),history()[0].cumulative,30e6)
    assert not result.diagnostics.converged and "prior state" in result.diagnostics.message


def test_fetkovich_tiny_timestep_uses_stable_exponential():
    model = FetkovichAquifer(1e6,1e-9,1e-30)
    result = model.compute_step(model.initial_state(30e6),30e6,20e6,1)
    assert result.incremental_influx == pytest.approx(5e-24,rel=1e-12,abs=0)


def test_phase_3b_models_retained_with_phase_3c_registration():
    # Phase 3C explicitly supersedes the Phase 3B exclusion of Modified VEH.
    assert list(MODELS) == ["None","Pot","Schilthuis","Fetkovich","Carter-Tracy","Van Everdingen-Hurst",
                            "Modified Van Everdingen-Hurst"]


def test_veh_reference_response_nodes_and_small_time():
    for tD, expected in ((0,0),(.01,.112),(.1,.404),(1,1.569),(10,7.537),(100,43.025),(1000,309.37)):
        assert veh_infinite_water_influx_dimensionless(tD) == pytest.approx(expected, rel=2e-5, abs=1e-12)
    assert veh_infinite_water_influx_dimensionless(1e-6) == pytest.approx(2*(1e-6/3.141592653589793)**.5)


def test_veh_pressure_step_superposition_independent_benchmark():
    k = .2*.001*1e-9*1000**2/DAY  # tD advances by one per day.
    model = VanEverdingenHurstAquifer(1000,10,10,.2,k,.001,1e-9)
    state = model.initial_state(30e6)
    first = model.compute_step(state,30e6,28e6,DAY)
    assert first.cumulative_influx == pytest.approx(39433.27098785909)
    second = model.compute_step(first.updated_state,28e6,27e6,DAY)
    assert second.cumulative_influx == pytest.approx(81216.45328060335)
    variables = dict(second.updated_state.model_variables)
    assert variables["diag_active_pressure_steps"] == 2
    assert variables["diag_current_step_contribution"] == pytest.approx(19716.635493929546)
    assert variables["diag_historical_contribution"] == pytest.approx(61499.8177866738)


def test_carter_tracy_first_step_independent_recurrence_benchmark():
    k = .2*.001*1e-9*1000**2/DAY
    model = CarterTracyAquifer(1000,10,10,.2,k,.001,1e-9)
    result = model.compute_step(model.initial_state(30e6),30e6,28e6,DAY)
    assert result.cumulative_influx == pytest.approx(31331.833883943957)
    diagnostics = dict(result.updated_state.model_variables)
    assert diagnostics["diag_tD"] == pytest.approx(1)
    assert diagnostics["diag_dimensionless_pressure"] == pytest.approx(.8021471491841932)
    assert diagnostics["diag_dimensionless_pressure_derivative"] > 0


@pytest.mark.parametrize("model", MODELS_UNDER_TEST[1:])
@pytest.mark.parametrize("matched", [False,True])
def test_aquifer_works_with_raw_and_matched_correlation_pvt(tank,model,matched):
    from material_balance_studio.pvt.fluid import fluid_from_inputs
    from material_balance_studio.pvt.models.correlation_model import CorrelationPVTModel, Transform
    fluid = fluid_from_inputs(api=35,gas_gravity=.75,temperature=90,initial_pressure=30e6,pb=18e6,rsb=90)
    pvt = CorrelationPVTModel(fluid,(2e6,35e6))
    if matched:
        pvt = pvt.with_transform(Transform("bo",.04,1.02))
    records = [replace(r,cumulative=replace(r.cumulative,gp=90*r.cumulative.np)) for r in history()]
    result = simulate(replace(tank,aquifer=model,pvt_model=pvt),records)
    assert result.converged
    assert all(s.balance.relative_residual<=1e-8 for s in result.states)
    assert result.states[-1].balance.aquifer_support>0


def test_mismatched_solver_context_is_controlled_failure(tank):
    model = PotAquifer(.001)
    context = AquiferContext(model,model.initial_state(30e6),30e6,DAY)
    a = solve_pressure(replace(tank,aquifer=model),history()[0].cumulative,29e6,aquifer_context=context)
    b = solve_pressure(replace(tank,aquifer=PotAquifer(.002)),history()[0].cumulative,30e6,aquifer_context=context)
    assert not a.diagnostics.converged and not b.diagnostics.converged
    assert a.aquifer_step is None and b.aquifer_step is None
