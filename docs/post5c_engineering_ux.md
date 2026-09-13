# Post-Phase-5C engineering workflow revision

This revision changes presentation, input adapters and orchestration. It adds no
Phase 6 forecasting, correlations, aquifer equations, network physics, allocation
optimizer or uncertainty engine.

## Regression evidence

Before application edits, the entire post-5C suite passed: **608 tests in
995.76 seconds**. A fresh latest Phase 5C acceptance run, including its 21-point
re-optimized profile, exactly matched the saved numerical results. Evidence and
pre-edit source/test SHA-256 hashes are preserved in `post5c_baseline/`.

Final acceptance: **632 tests passed in 783.93 seconds** (608 existing + 24 new).
After the compact-browser layout refinement, **46 relevant UI/workflow tests
passed in 40.69 seconds**. The fresh latest **Phase 5C and Phase 5B outputs match
their baselines exactly**; recursive Phase 5A and Phase 4C checks also passed.
Hash verification confirms **140 protected files unchanged**, with exactly the
three intentional application/presentation differences listed in the manifest.
The local app health endpoint returns `ok`; the browser default forward run
converges for all six timesteps with material-balance QC PASS.

The final suite and numerical acceptance are recorded in
`post5c_final_pytest.log` and `post5c_final_acceptance.log`. The repeatable runner
is `examples/run_post5c_acceptance.py`. It compares complete JSON structures,
without numerical tolerances, with the pre-edit 5C and 5B snapshots. The 5B runner
recursively reruns and verifies 5A and 4C acceptance, including earlier evidence.

`post5c_intentional_changes.json` records before/after hashes for the three
intentionally revised pre-existing application/presentation files. The existing
5C hash manifest renews only its two affected presentation entries. Its original
version remains in `post5c_baseline/`. The new regression test verifies every
pre-existing source/test against the fresh post-5C snapshot, allowing only these
explicit presentation changes. No pre-existing test or numerical-engine source
is edited.

## PVT and file handling

All three PVT sources remain available. The PVT workspace now follows Fluid
Inputs → Laboratory Data → Correlations → Match Data / Regression → Final PVT
Model → PVT Tables → PVT QC. In tabulated mode, fluid inputs are metadata; the
complete table remains the authoritative forward model.

Rsb is the scalar solution GOR at bubble point. Rs(P) is the pressure-dependent
solution ratio in the PVT model. Rp = Gp/Np is calculated only from history;
it stays missing when Np is zero. Separate field names and models prevent the
history ratio from becoming a fluid input.

API, gas gravity, temperature, salinity and Rsb are exposed. H2S, CO2 and N2 mole
percentages are stored separately and exported as informational composition.
They do not enter the correlation model or matching signature. Their total may
not exceed 100%. Existing applicability/QC rules remain authoritative.

Laboratory input accepts pressure with any subset of Bo, Rs, oil viscosity, z,
Bg, Bw and gas viscosity. Missing measurements remain missing. Gas viscosity is
retained and exported as informational lab data because no validated regression
for it exists. Existing lab validation still rejects invalid measurements and
duplicate pressures. Final tabulated PVT still requires pressure, Bo, Rs, Bg
and Bw, with all existing positivity and range checks.

Primary templates follow project units, including injection columns. SI uses
explicit MPa headers and FIELD uses psia, STB, bbl and scf. A tagged header must
agree with the selected file units. Legacy plain SI pressure remains Pa. New
uploads default to project units; explicit overrides persist for the same file.
Changing display units does not reinterpret an imported source. Computation
continues in canonical SI.

## History and Results

Enable **Explore history plots** in History to choose X, Y and Auto/Line/Scatter.
Available axes cover date, cumulative production/injection, measured pressure,
and history Rp. Date defaults to lines with markers where appropriate; numeric
X defaults to scatter. Only valid paired rows are shown. Pressure observations
are never filled, interpolated or replaced by calculated pressure. Axis/hover
units and large-volume scaling follow the project display system.

Results retains the temporal pressure diagnostic and adds a dedicated second
pressure-versus-Np panel. No Aquifer shows observed pressure and the without-
aquifer calculation only. An active model adds its actual name, for example
Calculated — With Fetkovich or Calculated — With Carter-Tracy.

The temporary reference is a frozen tank copy with only the aquifer replaced by
NoAquifer. The same N, m, PVT, history, injection, compressibilities and initial
conditions enter the original simulator. Observed points use Np on the same
history date. A failed reference is reported and only its converged states are
shown. Active simulation and accepted single/network scenarios are untouched.

## Aquifer Comparison

**Compare with Supplied Parameters** remains the default and calls the same
validated comparison simulator with exactly the entered models.

**History-Match Each Aquifer Before Comparison** calls the existing bounded
single-tank matcher separately for each formulation. Only explicitly selected
`aquifer.*` parameters from that model's registry may vary. Shared N, m, PVT,
production/injection and compressibilities remain fixed. No Aquifer is a fixed
baseline. A real model with no selected parameters also remains fixed.

The common observed-pressure list supports inclusion/exclusion and the existing
unweighted or pressure-uncertainty weighting policy. Missing sigma requires an
explicit default or exclusion. Initial values, bounds, transformations, fitted
values and bound flags are visible. Bounds and sigma retain canonical values
when display units change.

Results include RMSE, MAE, bias, maximum error, final We, MBE closure, bound
status, engineering QC and warnings. Bound-limited fits remain visible. The
best numerical pressure fit is explicitly distinguished from a physically
correct aquifer. CSV and JSON exports preserve the temporary comparison.
There is no automatic Apply operation or rerun of identifiability analyses.

## Files in this revision

- `app.py`: neutral header, templates, PVT organization, History and Results wiring.
- `presentation/pvt_workflow.py`: fluid/lab workflow, informational fields and tables.
- `presentation/aquifer_comparison.py`: default mode and optional calibrated workflow.
- `presentation/engineering_inputs.py`: explicit-unit templates and sparse lab adapter.
- `presentation/history_plots.py`: paired history plots and scaling.
- `presentation/pressure_oil.py`: immutable reference and pressure/Np panel.
- `presentation/aquifer_calibration.py`: existing-engine comparison orchestration.
- `presentation/aquifer_calibration_workflow.py`: bounds, observations and comparison UI.
- `tests/test_engineering_ux.py`: 24 added acceptance cases.
- `examples/run_post5c_acceptance.py`: deterministic regression verification.
- README, this report, baseline snapshots, manifests and final evidence logs.

Presentation paths above are relative to `src/material_balance_studio/`.
