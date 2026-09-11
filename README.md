# Material Balance Studio — Phase 3A

An auditable calculation foundation for a single equilibrated black-oil tank with
tabulated or correlation PVT, cumulative production/injection, and optional stateful aquifer support. The Streamlit
interface is a functional test harness. The engine runs without Streamlit.
Phase 1.1 stabilizes units and engineering diagnostics; the Phase 1 equations,
PVT interpolation and reservoir expansion/withdrawal physics are preserved. Phase 3A
extends domain state and the chronological solver to include aquifer influx explicitly.

Phase 2 adds independent property correlations, sparse laboratory imports/editing,
screening, transparent matching, and PVT QC. See the [PVT engineering specification](docs/pvt_phase_2.md)
for equations, sources, validity ranges, matching safeguards and the full file inventory.
The PVT section now offers **Tabulated PVT**, **Correlation PVT**, and
**Matched Correlation PVT**. The default remains the original tabular example.

Available families: Standing, Vasquez–Beggs and Glasø for Pb/Rs/Bo;
Beggs–Robinson/Vasquez–Beggs oil viscosity; Sutton pseudo-critical properties
and independently selected Dranchuk–Abou-Kassem z; real-gas Bg;
McCain Bw; Vasquez–Beggs undersaturated oil compressibility. Temperature input
is °C in SI or °F in FIELD, stored in K. At least measured Pb or Rsb is needed.
Source-file pressures remain Pa/psia; the lab editor uses MPa/psia.

Phase 2.1 adds dense before/after QC, full-range matching-deviation reports,
per-property model summaries, and browser visual acceptance. See the
[original-correlation reference audit](docs/pvt_reference_audit.md) and
[Phase 2.1 acceptance report](docs/phase_2_1_acceptance.md). Sutton's Ppc coefficient
is corrected from 131.07 to the original paper's 131.0. Undersaturated co is
reported only at/above Pb and is consistent with matched Bo. That acceptance pass made no MBE core changes.

Phase 3A adds **None, Pot, Schilthuis and finite Fetkovich** aquifers. Select an
aquifer under Reservoir Setup; inputs and results use FIELD/SI engineering units.
Immutable aquifer state advances only after convergence. The equation inspector
shows `Fnet = N Et + We`, with separate water influx, rate and pressure diagnostics.
Signed pressure-reversal efflux is retained with warnings. See the
[aquifer equations and references](docs/aquifer_phase_3a.md) and
[Phase 3A acceptance report](docs/phase_3a_acceptance.md). PVT physics is unchanged.

## Run (Python 3.12)

From this repository in PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m streamlit run app.py
```

The UI initially loads `examples/pvt.csv` and `examples/history.csv`; click **Run
simulation**, then open **Results / Equation Inspector**. Defaults describe a
synthetic tank with Pi = 30 MPa, N = 1,000,000 stock-tank m³, Swc = 0.2,
cf = 5e-10 /Pa, cw = 4e-10 /Pa, m = 0. Approximate synthetic pressure observations
illustrate comparison; they are not fitted or used to set calculated pressure.

```powershell
# Full tests, including Streamlit AppTest and Excel import:
.\.venv\Scripts\python.exe -m pytest -q -p no:cacheprovider
# Independent command-line engineering example:
.\.venv\Scripts\python.exe examples/run_example.py
# Phase 2 screening, matching and unchanged MBE integration:
.\.venv\Scripts\python.exe examples/run_pvt_example.py
# Aquifer hand benchmarks and pressure-support comparisons:
.\.venv\Scripts\python.exe examples/run_aquifer_acceptance.py
```

On macOS/Linux use `.venv/bin/python` in place of the Windows executable path.
Dependencies are declared in `pyproject.toml`; `requirements.txt` installs the
package and its test extra. Validated dataclasses provide immutable domain
records; Pydantic is not needed for this phase.

## Architecture and file inventory

The conventional `src/material_balance_studio/` package avoids ambiguous top-level
imports such as Python's built-in `io`. There are no engineering equations in
`app.py`.

```text
app.py                                      Streamlit input, plots, inspector
pyproject.toml / requirements.txt            Package, dependencies, pytest setup
src/material_balance_studio/
  domain/models.py                          Inputs, terms, states, solver/result records
  domain/validation.py                      Scalar and chronological history validation
  units/conversions.py                      Explicit SI/field boundary conversions
  units/display.py                          SI-to-engineering display conversion and labels
  pvt/base.py                               PVTModel protocol
  pvt/table_model.py                        Strict piecewise-linear table interpolation
  mbe/withdrawal.py                         Produced voidage and separate injection support
  mbe/expansion.py                          Eo, Eg, Efw, Et
  mbe/oil_balance.py                        Independent balance and residual APIs
  solver/pressure_solver.py                 Bounded Brent root solve and closure checks
  solver/timestep.py                        Prior-state transition and cumulative increments
  solver/simulation.py                      Initial state, time march, stop on failure
  diagnostics/balance_closure.py            Equation inspector and result DataFrame
  diagnostics/qc.py                         PASS / WARNING / FAIL presentation categories
  presentation/tables.py                    Readable converted engineering datasets and inspector
  presentation/charts.py                    Unit-aware pressure plot and relative-residual QC
  io/history.py                            CSV/XLSX schemas and unit conversion
  io/reservoir.py                          Reservoir input conversion adapter
examples/pvt.csv / history.csv               Small synthetic demonstration (SI)
examples/run_example.py                      Headless example
tests/conftest.py                            Independent analytic benchmark fixture
tests/test_pvt_table.py                      Nodes, interpolation, coverage, invalid data
tests/test_withdrawal.py                     Voidage, injection, zero-oil edge cases
tests/test_expansion.py                      Hand-calculated expansion and residual checks
tests/test_pressure_solver.py                Analytic recovery, injection, failure diagnostics
tests/test_simulation.py                     State carry, closure, observations, inspector
tests/test_history.py                        Validation, CSV/XLSX, reservoir safeguards
tests/test_units.py                          Gas basis and end-to-end SI/field equivalence
tests/test_app.py                            Bundled UI run and invalid-input handling
tests/test_presentation.py                   Display regression, state invariance and SI/FIELD equivalence
tests/test_hand_benchmark.py                 Independent fixed-value engineering benchmark
tests/fixtures/volumetric_hand_case.json      Hand-specified inputs and expected terms
docs/phase_1_1_validation.md                  Hand calculations and stabilization validation record
```

Each package directory also includes `__init__.py`; `.gitignore` excludes generated
environments, bytecode, test caches and build artifacts.

## Units and input formats

Internally: **Pa absolute**, surface component volumes in **m³ at their specified
surface reference conditions**, reservoir volumes in **reservoir m³**,
compressibility in **1/Pa**, viscosity in **Pa·s**, and calendar dates.
N/Np/Wp/Winj are stock-tank liquid m³; Gp/Ginj are standard gas m³. All surface
gas data and PVT must share the same reference temperature/pressure. Unit
conversion does not reconcile different standard conditions.

| Quantity | Canonical SI / SI file input | FIELD file input |
|---|---|---|
| Pressure | Pa absolute | psia |
| N, Np, Wp, Winj | stock-tank m³ | STB |
| Gp, Ginj | standard gas m³ | scf (**not** Mscf) |
| Bo, Bw, Bwinj | reservoir m³ / liquid surface m³ | rb/STB |
| Bg, Bginj | reservoir m³ / standard gas m³ | rb/scf (**not** ft³/scf) |
| Rs | standard gas m³ / stock-tank oil m³ | scf/STB |
| cf, cw | 1/Pa | 1/psi |
| Optional viscosities | Pa·s | cP |
| Swc, m, z | dimensionless | dimensionless |

**Project / display unit system** controls setup widgets, all engineering
previews, pressure plots, results, the equation inspector and engineering CSV
downloads. SI display pressure is **MPa** and compressibility is **1/MPa**;
FIELD display pressure is **psia**, oil is **STB**, water is **bbl**, gas is
**scf**, and reservoir withdrawal/support volumes are **rb**. SI surface and
reservoir volumes are labeled by their reference basis. Swc and m are
dimensionless. Every engineering table has readable, unit-bearing columns.

**Uploaded PVT/History file unit system** describes the values in each uploaded
file and is independent of display units. The canonical file contracts above
are unchanged: an SI file uses **Pa**, not MPa. Import first converts to
canonical SI; display then converts from that SI data to the project units.
Bundled examples are always parsed as SI and converted for display. File-unit
selectors are disabled while using bundled examples. Changing the display
system preserves edited reservoir values and does not recompute or modify a
saved simulation result.

Engineering CSV downloads contain converted values and readable unit-bearing
headers, including MPa in SI mode. They are reports, not input templates. For
imports, use the original SI/Pa templates provided in the import-format
expanders, or construct files with the exact import schemas below. The raw SI
result DataFrame/CSV is still available only in **Advanced / Raw Solver Output**.

**PVT:** CSV or first sheet of `.xlsx`, at least two rows, exact lowercase headers:

```csv
pressure,bo,rs,bg,bw
20000000,1.24,80,0.006,1.004
30000000,1.20,80,0.004,1.000
```

Optional columns: `bwinj,bginj,oil_viscosity,gas_viscosity,z`. Supplied optional
columns must be complete; omit unused columns. Pressure rows must already be
strictly increasing and unique. Piecewise-linear interpolation avoids cubic
overshoot. Out-of-range calls raise `PVTOutOfRangeError`; the solver never
extrapolates. Initial pressure outside coverage is rejected before simulation.
Unsorted tables, missing values, nonfinite values, nonpositive FVFs/viscosities/z,
negative Rs and unknown columns are rejected. Legacy `.xls` is not supported;
save it as `.xlsx` first.

**History:** exact lowercase headers `date,np,gp,wp`, optional
`winj,ginj,observed_pressure`. Missing injection columns mean zero injection;
blank cells in a supplied cumulative column are errors. Observations may be blank
or the entire column omitted. CSV dates are `YYYY-MM-DD`; Excel date cells work.
Dates must already be sorted and unique; input is never silently reordered.
Negative/decreasing/nonfinite cumulative production **or injection** is rejected.

The reservoir initial date defines **zero** cumulative production/injection.
A history row on that date must contain zeros. The first row may be later and
contain nonzero totals accumulated since the initial date. No history may
precede the initial date. The engine derives increments from these totals.

## Material-balance formulation

For a black-oil tank, with initial properties evaluated from the same PVT model
at Pi, the implemented equations are:

```text
Fprod = Np Bo + (Gp - Np Rs) Bg + Wp Bw
Iwater = Winj Bwinj
Igas = Ginj Bginj
Fnet = Fprod - Iwater - Igas

Eo = (Bo - Boi) + (Rsi - Rs) Bg
Eg = Boi (Bg/Bgi - 1)
Efw = Boi (1 + m) (cf + cw Swc)/(1 - Swc) (Pi - p)
Et = Eo + m Eg + Efw

ExpansionSupport = N Et
Residual(p) = Fnet(p) - N Et(p)
AbsoluteResidual = abs(Residual)
RelativeResidual = abs(Residual) / max(abs(Fnet), normalization_floor)
```

The expansion convention follows the black-oil relation in [Thomas Blasingame,
Texas A&M PETE 324, material balance lecture, slide 6](https://blasingame.engr.tamu.edu/z_zCourse_Archive/P324_07A/P324_07A_l_Lessons/070216_Lec_Material_Balance.pdf).
In this implementation **Efw already includes (1 + m)**; it must not be
multiplied by that factor again. m is initial gas-cap reservoir volume divided
by initial oil reservoir volume. No field conversion factors enter these equations.

The division-free production expression is algebraically equivalent to
`Np [Bo + (Rp - Rs) Bg] + Wp Bw`, and also works for zero oil production.
The signed `Gp - Np Rs` correction is retained. A negative correction produces
a review warning at a solved timestep rather than being silently clipped.

Injection FVFs are evaluated at **candidate tank pressure**, not wellhead
conditions. Optional Bwinj/Bginj tables allow separate FVFs; when omitted,
resident Bw/Bg are used and that assumption is recorded as a result warning.

## Stateful time march and solver

`simulate(tank, history, settings)` constructs State 0 at Pi with zero cumulative
volumes. `advance_timestep` consumes the previous state, derives increments,
updates totals, uses previous pressure as the initial guess, solves the current
balance, and commits the new immutable state. Observations are used only for
reported `calculated - observed` pressure error.

The no-aquifer equation remains a **cumulative equilibrium balance referenced to
initial PVT**. Every step re-evaluates all cumulative voidage at its candidate
pressure. Summing past reservoir-volume withdrawals at past pressures would
change this equation. Therefore a final cumulative history endpoint is invariant
to adding intermediate dates in this Phase 1 model, while the execution and
state transition architecture are chronological. Elapsed days and increments
are stored for future transient models; no rate-dependent model is claimed here.

Bounds are the intersection of PVT coverage, a positive pressure floor and any
user maximum. Without cumulative injection the upper bound is Pi. Injection
permits pressure above Pi within measured coverage; the residual must still
support that root. The solver samples PVT intervals, including the previous
pressure, selects a sign-changing bracket nearest previous pressure, then uses
`scipy.optimize.brentq`. Multiple detected brackets generate a warning. Sampling
does not guarantee detection of all tangential/closely spaced roots in arbitrary
nonphysical tables; engineers must inspect questionable PVT/pressure branches.

Default settings:

| Setting | Default |
|---|---|
| Pressure absolute tolerance | 1e-5 Pa |
| Brent pressure relative tolerance | 1e-14 |
| Absolute balance tolerance | 1e-5 reservoir m³ |
| Relative balance tolerance | 1e-8 |
| Relative normalization floor | 1 reservoir m³ |
| Maximum Brent iterations | 100 |
| Minimum pressure | 1 Pa (also limited by table coverage) |
| Bracket subdivisions per knot interval | 8 |

**Both** absolute and normalized relative closure tests must pass; pressure
convergence alone is not acceptance. The normalization floor prevents division
by zero when production and injection cancel. It is configurable for the scale
of the modeled reservoir. Residuals retain a sign as well as absolute magnitude.

No bracket, iteration exhaustion or failed closure returns non-converged solver
diagnostics. Simulation stops at that date, retaining previous successful states
and a separate `failed_timestep`; failed attempts are not passed to later steps.
The diagnostic includes attempted bounds, selected bracket, iteration/evaluation
counts and final residual. Invalid inputs raise `EngineeringValidationError`,
which the UI displays. No engineering inputs are silently corrected.

## Inspector API

The primary UI inspector uses engineering tables for the saved left side,
right side, signed difference and all expansion/withdrawal contributions.
**Advanced / Raw Solver State** retains the original canonical JSON API.

The primary closure plot shows the **absolute relative balance residual**,
with a zero-based linear axis that includes the 1e-8 PASS limit. This keeps tiny
floating-point residuals near zero rather than amplifying them. Data are not
clipped, rounded or replaced with a plotting epsilon. Exact small values remain
available in scientific notation in hovers, tables and advanced diagnostics.

| Presentation QC | Absolute relative balance residual |
|---|---|
| PASS | ≤ 1e-8 |
| WARNING | > 1e-8 and ≤ 1e-5 |
| FAIL | > 1e-5, unavailable/nonfinite error, or solver failure |

Presentation QC does not relax numerical acceptance: the solver still requires
both its original absolute and relative tolerances. Failed attempts remain
separate from successful history states. Absolute residuals are available in
**Advanced Diagnostics / Absolute Residuals**, in the selected engineering units.

```python
from material_balance_studio.diagnostics.balance_closure import inspect_timestep
from material_balance_studio.mbe.oil_balance import balance_residual
from material_balance_studio.solver.simulation import simulate

result = simulate(tank, history)
if result.converged:
    detail = inspect_timestep(result, result.states[-1].date)
    residual = balance_residual(result.states[-1].pressure, tank,
                                result.states[-1].cumulative)
```

The inspector exposes Fprod, each injection support, Eo, Eg, mEg, Efw, Et,
each N-scaled expansion support, N Et, Fnet, signed/absolute/relative residual,
calculated/observed pressure and pressure error, PVT values, cumulative totals,
increments, previous pressure, elapsed days and solver diagnostics. Result tables
contain successful history rows; `initial_state` is separately accessible.

## Validation and engineering assumptions

See [Phase 1.1 validation and the explicit hand benchmark](docs/phase_1_1_validation.md)
for the fixed-input 20 MPa case, independently calculated terms, display
regressions and canonical-SI invariance checks.

Tests include the six required acceptance cases: A zero production, B analytic
volumetric depletion, C water injection support, D closure at every timestep,
E decreasing history rejection, and F exact-node/midpoint PVT interpolation.
Additional checks cover a hand-calculated combined gas-cap/rock-water/injection
case, gas-only and water-only flows, pressure above Pi with injection, solver
failure retention, optional observations, SI/field equivalence, CSV/XLSX imports,
the stateful inspector and the running Streamlit interface. Benchmark expected
pressures/terms use independent analytical or hand-calculated values.

Assumptions: uniform average pressure and fixed temperature; supplied black-oil
PVT represents the reservoir; fixed N and initial m; constant cf/cw with a
first-order rock/connate-water expansion approximation; initial Swc applies to
the tank; surface quantities use consistent standard conditions; injected gas
provides explicit voidage support but does not alter the supplied fluid PVT.
This is not a compositional or miscible-injection simulator. Mathematical closure
does not independently establish physical validity of measured inputs.

See [Phase 2 validation](docs/phase_2_validation.md) for retained/new tests and
the original core-file hash comparison. The synthetic Bo matching example
recovers A=.04, B=1.02 and reduces RMSE from .0644754 to numerical zero.

**Intentionally deferred:** Carter–Tracy, VEH and modified VEH aquifers, aquifer
parameter matching, reservoir-parameter history matching,
multi-tank systems, transmissibility, forecasting, DCA, uncertainty analysis
and PDF reporting.
