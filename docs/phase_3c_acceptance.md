# Phase 3C engineering acceptance

Implemented PETEX User Guide Appendix C C2.9 linear-pressure-history convolution,
using the original VEH dimensionless response. Geometry remains infinite-acting
radial; no finite boundary response or Hurst Modified Steady-State is claimed.
See [equations, references and benchmark derivation](aquifer_phase_3c.md).

The independent early-time, three-segment reference-equation benchmark has maximum
relative cumulative-influx error 1.464e-16. It validates pressure slopes,
dimensionless time, incremental and cumulative influx. This is an analytic
reference-equation test, not a commercial software equivalence certification.

For a 1 MPa linear decline over tD=1, one original pressure step gives
19716.635493929545 m³ and the modified segment gives 12522.52104974676 m³
(57.4494% difference relative to modified). With 100 steps the original gives
12618.232523658779 m³ and modified gives 12522.521049746765 m³ (0.7643%).
Original VEH applies the full drop at each step start; the modified model spreads
it through the interval. The comparison demonstrates discretization effects,
not arbitrary superiority.

All nine Phase 3B fixed support cases covering every previous model match their
saved outputs exactly after JSON normalization. Before/after SHA-256 hashes are
recorded in [numerical evidence](phase_3c_numerical_results.json). Previous benchmark
records also match exactly. Maximum relative MBE residual across the acceptance
cases is **1.636292862885668e-14**.

Validation completed: **414 passed in 102.33 seconds**, comprising the 390
existing cases and 24 new Phase 3C cases. One obsolete registry assertion
that explicitly prohibited Modified VEH was updated because Phase 3C requires its
registration; its assertions for all six previous models remain. Added coverage
includes analytic benchmark, independent quadrature, repeatable/rejected candidates,
rebound, solver failure, model switching, FIELD/SI equivalence, all models alone
versus comparison, metrics/dataframes/charts, QC, immutable sensitivity, exact
regression outputs, and a Streamlit comparison/unit-switching workflow.

Run validation:

```
.\.venv\Scripts\python.exe -m pytest -q -p no:cacheprovider
.\.venv\Scripts\python.exe examples/run_phase3c_acceptance.py
```

Comparison includes all requested models plus NoAquifer, independent parameter
sets, seven aligned output charts, observed-pressure errors, summary CSV, same-date
inspector, visible QC and optional sensitivity. Numerical fit is explicitly
separate from engineering acceptability. Automatic aquifer/reservoir history
matching and all other deferred optimization/forecasting features are not included.
