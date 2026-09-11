# Phase 3B acceptance report

Command evidence:

```powershell
.\.venv\Scripts\python.exe -m pytest -q -p no:cacheprovider
.\.venv\Scripts\python.exe examples\run_aquifer_acceptance.py
```

Result: **390 tests passed**. The 349 pre-Phase-3B tests remain passing, with 41
new Phase 3B tests covering transient response functions, benchmarks,
immutability, UI integration, unit equivalence and validation.

## Implemented models

Phase 3B adds `CarterTracyAquifer` and `VanEverdingenHurstAquifer` under
`src/material_balance_studio/aquifer/`. Both use the existing Phase 3A
`AquiferModel` interface and are evaluated by the solver through `AquiferContext`.
No model-specific logic was added to the main time-marching solver.

Original VEH is implemented as an infinite-acting radial constant-terminal-pressure
solution with pressure-step superposition. Finite VEH is not claimed as supported.
Carter-Tracy is implemented as the recursive transient aquifer formulation using
dimensionless pressure and its derivative.

Modified VEH was not added.

## Benchmarks

The generated artifact is
`docs/phase_3b_numerical_results.json`.

| Benchmark | Expected | Actual |
| --- | ---: | ---: |
| VEH day-1 cumulative We, m3 | 39433.27098785909 | 39433.27098785909 |
| VEH day-2 cumulative We, m3 | 81216.45328060335 | 81216.45328060335 |
| Carter-Tracy day-1 cumulative We, m3 | 31331.833883943957 | 31331.833883943957 |

The shared synthetic material-balance case gives:

| Model | Final pressure, MPa | Final We, m3 |
| --- | ---: | ---: |
| Carter-Tracy | 29.696609794 | 34822.845999 |
| Infinite VEH | 29.716605367 | 34900.428822 |

The trends are consistent and the difference is documented as formulation
difference: Carter-Tracy uses a recurrence based on constant-terminal-rate
pressure response, while VEH uses pressure-step superposition of the
constant-terminal-pressure cumulative influx response.

Maximum relative material-balance residual in the acceptance comparison:
`1.636292862885668e-14`.

## Validation

New tests verify:

- VEH response nodes and small-time behavior.
- VEH multi-pressure-step superposition with current and historical
  contribution diagnostics.
- Carter-Tracy recurrence terms and dimensionless pressure diagnostics.
- Candidate evaluations remain immutable and repeatable.
- Accepted state is committed only after convergence.
- FIELD/SI transient-aquifer inputs produce equivalent simulations.
- Invalid transient geometry and rock/fluid inputs are rejected.
- Permeability, viscosity and encroachment-angle sensitivities are physically
  consistent.
- UI selectors expose only relevant Carter-Tracy and VEH inputs and produce the
  expected transient plots.

PVT source-file integrity tests remain passing. No PVT correlation, expansion,
withdrawal, injection or material-balance residual definitions were changed.

## Changed files

Key additions:

- `src/material_balance_studio/aquifer/transient.py`
- `src/material_balance_studio/aquifer/veh_response.py`
- `src/material_balance_studio/aquifer/carter_tracy_response.py`
- `src/material_balance_studio/aquifer/veh.py`
- `src/material_balance_studio/aquifer/carter_tracy.py`
- `docs/aquifer_phase_3b.md`
- `docs/phase_3b_acceptance.md`
- `docs/phase_3b_numerical_results.json`

Key integrations:

- Aquifer registry and package exports.
- Presentation-unit conversions for length, permeability and angle.
- Streamlit aquifer setup and plots.
- Equation inspector transient diagnostics.
- Phase 3B tests in `tests/test_aquifer.py` and `tests/test_aquifer_app.py`.
