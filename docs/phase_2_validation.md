# Phase 2 validation record — 10 September 2026

Phase 2 implementation and verification completed in `E:\Code_Folder\MBE_tool`.
Python 3.12.10; existing `.venv` and installed project dependencies.

## Automated acceptance

Command:

```powershell
.\.venv\Scripts\python.exe -m pytest -q -p no:cacheprovider
```

Result: **228 passed in 12.00 seconds**. The cache plugin is disabled because
the existing Windows pytest-cache directory has permission restrictions.

| Test group | Tests passed |
|---|---:|
| Original Phase 1/1.1 tests, retained unchanged | 136 |
| New correlation/reference and metadata tests | 37 |
| New model, lab, ranking, matching, QC and MBE integration tests | 46 |
| New Streamlit PVT workflow regression tests | 9 |
| **Total new tests** | **92** |
| **Total tests** | **228** |

Coverage includes fixed reference values for every registered correlation;
published rounded viscosity/compressibility and DAK examples; unit conversions;
Rs/Bo/viscosity behavior below, at and above Pb; each supported anchor configuration;
ranking by applicability before error; exact lab-pressure predictions; all fit
treatments; failed regression safeguards; sparse and invalid CSV/XLSX inputs;
curve and coverage QC; MBE depletion/injection; identical-properties regression;
and tabular/raw/matched UI operation.

UI regression checks confirm that manual overrides clear stale matches, unit
changes preserve canonical fluid/lab/fit/saved-result values, matched mode needs
an enabled accepted fit, raw mode keeps its untuned model, sparse editor changes
are interpreted in displayed units, invalid temperature/anchors are recoverable,
and sour gas blocks the final model. No visual browser screenshot audit is claimed;
UI verification uses Streamlit AppTest and Plotly figure specifications.

## Headless example

Command:

```powershell
.\.venv\Scripts\python.exe examples/run_pvt_example.py
```

Synthetic fluid: API 35, gas SG .75, 90°C, Pi 30 MPa, measured Pb 18 MPa,
Rsb 90 m³/m³. Range 2–35 MPa. Six sparse Bo observations use
`Bo_lab=.04+1.02 Bo_raw`.

| Quantity | Observed result |
|---|---:|
| Fitted A | .03999999999999937 m³/m³ |
| Fitted B | 1.0200000000000002 |
| Laboratory points | 6 |
| Before RMSE | .06447544321627238 m³/m³ |
| After RMSE | 4.1540741810552243e-16 m³/m³ |
| Before bias | −.06445772161037733 m³/m³ |
| After bias | −4.070817756958907e-16 m³/m³ |
| Before MAPE | 5.0182910222282855% |
| Raw/matched physical QC failures | 0 / 0 |
| Raw/matched MBE runs | Both converged |
| Raw final pressure | 16.860877871149906 MPa |
| Matched final pressure | 16.725425122873613 MPa |

These differences in MBE pressure arise solely from the different PVT curves.
The separate identical-properties regression supplies the same Bo/Rs/Bg/Bw
through `CorrelationPVTModel` and `TablePVTModel`; pressures agree within 1e-6 Pa
and balance terms within numerical tolerance using the unchanged solver.

The example also demonstrates the ranking policy: Standing and Glasø are VALID;
Vasquez–Beggs has a lower Bo error than Standing here but is CAUTION until its
separator gas-gravity basis is verified. It does not silently become the
recommended model on error alone; the engineer can select it manually.

## Core preservation

SHA-256 comparisons against the pre-implementation values found **zero changed
files among all 13 original core files**: every Python file under `mbe/`,
`solver/`, `domain/`, and `pvt/base.py` / `pvt/table_model.py`, including package
initializers. See [machine-readable hash record](phase_2_core_hashes.json).
Consequently the Phase 1 withdrawal, expansion, residual, time marching,
pressure solver, domain validation and tabular interpolation are byte-for-byte
unchanged. Existing tests were not edited.

All added/changed files, correlations grouped by property, sources, exact
equations, applicability limits, lab workflow, ranking methodology, matching
formulations and QC rules are listed in the [PVT engineering specification](pvt_phase_2.md).

## Practical limits

This is a sweet black-oil correlation workflow with optional lab calibration.
It does not implement composition/acid-gas or separator-condition corrections,
gas viscosity, explicit salinity effects, or new injection-fluid physics.
Correlation ranges and sampled physical QC are not proof of fluid accuracy.
The matching example is synthetic. Aquifers, reservoir-parameter history
matching, multiple tanks, forecasting, DCA and uncertainty remain deferred.
