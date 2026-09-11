# Phase 4A acceptance

**439 tests passed in 42.34 seconds: all 414 existing cases plus 25 new cases.**
No previous test was changed. UI coverage builds Reservoir Diagnosis, displays
the diagnostic plots/tables, switches to FIELD units, and verifies unchanged
saved simulation state.

Implemented modules: pressure/history QC and observation metadata, voidage,
VRR, Havlena–Odeh, oil Campbell, drive indices, pressure match and diagnosis
composition. The dedicated Reservoir Diagnosis tab includes row flags, voidage
export, instantaneous/cumulative/pressure-overlay VRR plots, observed-pressure
H–O trends and regression QC, Campbell overlays, individual/stacked drive indices,
pressure residuals/normalized residuals, a same-date inspector, and reuse of saved
aquifer comparison results. The existing equation inspector includes VRR and
drive decomposition.

See [the engineering specification](diagnostics_phase_4a.md) for equations,
pressure bases, references, thresholds, unit conventions and limitations.

## Independent deterministic evidence

* VRR tests: 10 stock-tank m³ oil at Bo=1.2 and Gp=Np Rs produce 12 reservoir m³.
  Injecting 12 m³ water at Bw=1 gives VRR=1; injecting 600 standard m³ gas at
  Bginj=.02 also gives VRR=1; 6 m³ water plus 300 m³ gas gives VRR=1. Half these
  combined volumes give .5. Independent cumulative totals check a different
  ratio (1.5× the interval ratio). Zero production produces no numeric ratio.
* Component benchmark: produced oil=12, free gas=2, water=3, total=17 m³;
  injected water=5.5, gas=2, total=7.5 m³; net=9.5 m³.
* Five synthetic pressure/volume histories constructed algebraically for known
  N=1,000,000 m³ cover depletion, aquifer, gas cap, injection and mixed support.
  H–O slopes range from 999999.9999999987 to 1000000.0 m³ with R²=1 and
  negligible intercepts. Aquifer-adjusted Campbell values equal the known N.
* A separate non-collinear gas-cap transform recovers intercept N=1,000,000 and
  slope Nm=500,000, without using the tested diagnostic to generate expected data.
* Drive tests cover pure depletion, gas support, aquifer, injection and five-term
  mixed support (20%,10%,5%,40%,25%). A deliberate nonclosing example remains
  at sum=.5, demonstrating that no arbitrary normalization is applied.
* Forward synthetic drive sums close within **2.665e-15**. In the aquifer case
  final WDI=.555556; gas-cap case SDI=.789474; injection case IDI=.714286.
* Pressure errors [+1,−2,+3] MPa give RMSE=2.160246899 MPa, MAE=2 MPa,
  bias=.666666667 MPa, maximum=3 MPa. First sigma=.5 MPa gives normalized
  residual=2; zero/missing sigma leaves normalization unavailable.
* FIELD/SI tests convert volume, gas, pressure and sigma input data, compare VRR
  and recovered H–O slope, and verify display values. QC tests retain suspicious
  input rows. Missing observations never borrow calculated pressure for H–O.

## Regression protection

The maximum forward-model relative MBE residual across retained and new acceptance
cases is **1.636292862885668e-14**. All protected source hashes (domain, MBE,
PVT, aquifer and solver) match the pre-Phase-4A manifest. The entire Phase 3C
acceptance output remains exactly identical after JSON normalization, including
all prior aquifer output hashes and modified VEH benchmarks.

[Numerical evidence and hashes](phase_4a_numerical_results.json) are reproducible:

```
.\.venv\Scripts\python.exe -m pytest -q -p no:cacheprovider
.\.venv\Scripts\python.exe examples/run_phase4a_acceptance.py
```

No automatic history matching, OOIP/gas-cap/aquifer parameter updates, weighted
optimization, forecasts or uncertainty analysis were added. The H–O line is a
diagnostic trend only. The supported Campbell transformation is explicitly for
black oil; a dry-gas Cole plot is not mislabeled as an implemented oil diagnostic.
