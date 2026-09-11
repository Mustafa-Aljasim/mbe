# Phase 3A acceptance report

Completed 11 September 2026. **349 tests pass: all 262 existing tests retained
unchanged, plus 87 new tests.** Full-suite command:

```powershell
.\.venv\Scripts\python.exe -m pytest -q -p no:cacheprovider
# 349 passed in 51.67s
.\.venv\Scripts\python.exe examples/run_aquifer_acceptance.py
```

## Implementation and preserved behavior

Created a common aquifer protocol, immutable configuration/state/result records,
and NoAquifer, Pot, Schilthuis and finite Fetkovich implementations under
`src/material_balance_studio/aquifer/`. Validation is colocated with the common
base and records. Model-specific equations remain outside the MBE module.

ReservoirTank now holds the aquifer configuration. SimulationState retains the
committed aquifer state and interval diagnostics. State contains cumulative We,
aquifer pressure, previous reservoir pressure, elapsed seconds, initial pressure,
model identity/parameters and an immutable model-variable extension tuple.

The solver receives an AquiferContext for the prior committed state and timestep.
Candidate pressure evaluations are pure; only a successful pressure solution
returns a proposed aquifer state for commitment. Failed trials never advance
aquifer time or influx. Models/parameters cannot reuse incompatible state.

`evaluate_balance` accepts a separate reservoir-volume We and computes
`Fnet − N Et − We`. Total expansion retains its original meaning; total support
adds We explicitly. NoAquifer reproduces every pre-existing field in the captured
six-step baseline **exactly**, including pressure, balance terms, PVT, solver
iterations/calls and warnings. No existing test expectations were changed.

All **25 PVT Python files** and the existing expansion/withdrawal implementation
files match pre-change SHA-256 hashes. The MBE balance adapter, domain records
and solver were extended as required; their files are not claimed to be unchanged.
PVT equations, correlations, matching, interpolation and validation are unchanged.

## Equations and references

The [engineering specification](aquifer_phase_3a.md) documents publications,
parameter definitions, units, initial-equilibrium assumption, flow reversal,
numerical treatment and limitations. Exact implemented forms:

- None: delta We=We=0.
- Pot: `We=Caq(Pi−p)`, with `Caq=Wi ct`.
- Schilthuis: `delta We=J[Pi−(p_previous+p)/2]dt`, added to prior We.
- Fetkovich: `delta We=Caq[pa_previous−(p_previous+p)/2][1−exp(−Jdt/Caq)]`;
  `pa=Pi−We/Caq` after adding the increment.

References are Schilthuis (1936), AIME 118, 33–52; Fetkovich (1971), SPE 2603,
JPT 23(7), 814–828; and the cited NTNU/USGS compressible-storage treatments for
Pot. Original-paper indexed equations and primary engineering references were
used; source-access limits are stated in the specification. No open-source code
served as the sole engineering authority.

## Independent benchmark outcomes

| Case | Expected hand result | Outcome |
|---|---|---|
| NoAquifer | Exact captured Phase 2.1 six-step result, zero We | Exact match |
| Pot, Caq=2500 m³/MPa; p=28,25,27 MPa from Pi=30 | We=5000,12500,7500 m³ | Exact values; reversal retained and warned |
| Schilthuis, J=100 m³/(day·MPa); intervals 10,20,5 days | delta We=1000,6000,1750; cumulative=1000,7000,8750 m³ | Matches within floating-point tolerance |
| Fetkovich, Caq=1000 m³/MPa, J=Caq ln2/day | delta We=1000,1500,750,375; cumulative=1000,2500,3250,3625 m³ | Zero recorded numerical error in benchmark artifact |
| Same Fetkovich case | pa=29,27.5,26.75,26.375 MPa | Zero recorded numerical error |

Fetkovich's first interval changes boundary pressure from 30 to 26 MPa; subsequent
intervals hold 26 MPa. After that transition, support halves each day as finite
aquifer pressure depletes. Expected values are constants derived by hand using
an exponential fraction of exactly 1/2, not calls to production functions.

## Qualitative pressure-support comparison

Synthetic table case: Pi=30 MPa, N=1,000,000 m³, constant Rs=80,
Bo=1.2+4e−9(30e6−p), m=cf=cw=0. Cumulative oil production is 10000, 20000,
30000 m³ on 1 February, 1 March and 1 April 2020; Gp=80 Np. All runs converge.

| Model | Weak final pressure MPa | Strong final pressure MPa | Weak final We m³ | Strong final We m³ |
|---|---:|---:|---:|---:|
| NoAquifer | 20.721649485 | — | 0 | — |
| Pot | 22.622950820 | 27.406340058 | 7377.049180 | 25936.599424 |
| Schilthuis | 21.591479547 | 25.967659864 | 3374.940640 | 20354.520271 |
| Fetkovich | 21.402050408 | 25.127255418 | 2639.955585 | 17093.751021 |

Weak/strong Pot Caq=.001/.01 m³/Pa; Schilthuis J=1e−10/1e−9 m³/(Pa·s).
Fetkovich uses Wi=1e6/1e7 m³, ct=1e−9 /Pa, J=1e−10/1e−9 m³/(Pa·s).
At every history step, pressure satisfies **none < weak < strong < Pi** for each
family. Maximum relative MBE residual over these **21 solved steps** is
**6.452864909078926e−15**, below the 1e−8 acceptance limit. This maximum describes
the documented comparison fixtures, not an estimate for arbitrary user cases.

## Presentation and QC

Reservoir Setup shows only relevant inputs for the selected model, with physical
help text and explicit initial-equilibrium/reversible-flow assumptions. Inputs
and outputs use MPa/psia, reservoir m³/rb, per-day rates, compressibility and PI
units. SI/FIELD switching preserves edited inputs and saved canonical results.

Results include cumulative We versus time, average influx rate for Schilthuis and
Fetkovich, and aquifer pressure for Fetkovich. The inspector separately reports
production, water/gas injection, expansion, delta We, We, total right-hand support,
residual, previous/current aquifer pressure and driving pressure difference.
Raw state remains in Advanced expanders.

Automated UI tests cover all four selectors, relevant-field visibility, runs,
chart counts, invalid-input blocking, edited-value preservation, unit conversion,
and model switching with fresh state. The running app's default Fetkovich example
was also inspected in the browser, including all three aquifer plots in SI and
FIELD: readable reservoir m³/rb, reservoir m³/day/rb/day and MPa/psia axes.
Six history rows converged, maximum relative
residual 5.651e−15; cumulative influx and declining finite pressure are visible.
The final SI inspector displayed Fnet and N Et+We as 60648.293 m³, with signed
residual 1.4551915e−11 m³.

Negative/reversed flow is not clipped: it is retained and warned. Invalid
capacity, J, volume, compressibility, pressure, time, state consistency or context
produce validation/control failures. Extreme support relative to net withdrawal
is flagged. The legacy CSV round-trip tests detected an all-missing pressure
column dtype issue; the numeric presentation conversion was corrected and all
existing tests now pass.

## New test coverage and artifacts

The 87 new tests cover hand benchmarks; exact no-aquifer regression; physics
hashes; candidate immutability/determinism; chronological cumulative consistency;
initial-date reference rows; weak/strong support ordering; failed-step retention;
Brent sampling-density independence; invalid parameters and timesteps; model and
context mismatch; signed injection efflux; equation-inspector terms; independent
unit scales; FIELD/SI simulation equivalence; tiny-step numerical stability;
raw/matched correlation PVT with all three aquifers; and UI integration.

- [Pre-change baseline and hashes](phase_3a_baseline.json)
- [Reproducible numerical evidence](phase_3a_numerical_results.json)
- [Benchmark runner](../examples/run_aquifer_acceptance.py)
- [Engineering tests](../tests/test_aquifer.py)
- [UI tests](../tests/test_aquifer_app.py)

No Carter–Tracy, VEH, modified VEH or other Phase 3B models were added. No aquifer
history matching, multi-tank flow, forecasting, optimization, DCA, uncertainty
analysis or Phase 4 drive-index system was implemented.
