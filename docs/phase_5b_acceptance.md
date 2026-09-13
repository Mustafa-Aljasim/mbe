# Phase 5B acceptance

**571 tests passed in 341.67 seconds:** all 533 existing tests retained unchanged,
plus 38 new Phase 5B tests. The final combined run completed on 2026-09-12;
its output is retained in [the test log](phase_5b_pytest.log).

Phase 5B adds effective-date, stream-specific allocation and bounded network
history matching around the unchanged Phase 5A simulator. The allocation data
model, parameter identities, optimizer/scaling, pressure objective, closure policy,
scenario application and signed drive diagnostics are documented in
[the engineering method](network_history_matching_phase_5b.md).

## Independent recovery benchmark

The acceptance fixture generates observations with an independent two-by-two
linear compressible-storage recurrence, not with the production network solver or
history matcher. Both compartments have prescribed piecewise-constant water
production/injection rates, with changing pressure gradients and flow reversal.
Eight survey dates provide sixteen pressure observations for the joint fits.

Storage C_A=N_A(c_o+m_A Boi c_g/Bgi), C_B=N_B c_o, and a two-node graph Laplacian
L give the discrete oracle

```
(C + dt L/2) Pnew = (C - dt L/2) Pold - ΔF
```

For the Fetkovich case, independently calculate
beta=Caq(1−exp(−J dt/Caq)), add beta/2 to the B diagonal, add
beta(Paq_old−P_B_old/2) to the B right-hand side, then update cumulative aquifer
influx and P aq from the solved pressure. The independent oracle and forward
matching runs use a two-day maximum step, kept below the communication stiffness
limit throughout the recovery bounds. This makes the discretization controlled
without introducing an alternate optimizer forward model.

True parameters are N_A=1,000,000 m³, N_B=1,500,000 m³, T_AB=8e-10
reservoir m³/(Pa·s), J_B=1e-9 reservoir m³/(Pa·s), and m_A=0.2 in its dedicated
gas-cap case. Perturbed initial values are used for recovery.

| Recovery case | Required verification |
|---|---|
| T only | Relative recovery error ≤2e-4 and pressure RMSE <1 Pa |
| N_A, N_B and T | Each relative recovery error ≤5e-4 and pressure RMSE <1 Pa |
| Fetkovich J_B and T | Each relative recovery error ≤5e-4 |
| Gas-cap m_A | Relative recovery error ≤5e-4 |
| Restrictive T upper bound | Fit reaches upper bound and every evaluated T stays inside bounds |
| Isolated true T=0 | Linear parameterization approaches lower boundary, flags weak/bound-limited communication |

The recorded T-only fit gives **8.000000000000317e-10**. Joint recovery gives
N_A=**1000000.0000000849**, N_B=**1500000.0000000326**, and
T=**8.000000000000073e-10**. The aquifer/connection fit returns
J_B=**1.0000000000005448e-9** and T=**8.000000000000131e-10**. Gas-cap recovery
gives m_A=**0.19999999999994972**. These noiseless recovery fits have pressure RMSE
below 1e-6 Pa. The isolation fit approaches T=**4.338402570792328e-15**, with a
lower-bound caution; it is not silently rounded or clamped to zero.

The T-only fitted simulation includes both positive and negative rates on the same
connection. No duplicate direction-specific link or clipped signed transfer is
introduced during matching.

## Allocation and weighting benchmarks

Independent arithmetic checks allocate all five cumulative streams, with different
fractions by stream and a change inside a field-history interval. They verify each
interval and cumulative sum, including an extra pressure survey between allocation
dates. The old fraction applies through the interval ending at the change date;
the new fraction applies afterwards. Oil allocation is not reused for water or
injection. Negative, above-one, incomplete, duplicate and non-unit-sum schedules
are rejected, as is implicit mixing with direct history.

The wrong-allocation case uses a true water-production share changing from 80/20
to 30/70, with injection assigned to B. Correct allocation recovers the known
transmissibility with pressure RMSE below 1 Pa. Replacing the schedule with a wrong
fixed 50/50 split produces a materially worse fit (RMSE >10,000 Pa) and biased T.
The optimizer never changes any allocation fraction.

Numerically, correct allocation gives T=**7.999999999207774e-10** and RMSE about
**1.24e-5 Pa**. Wrong fixed allocation gives T=**1.0106820269782596e-9** and
RMSE **90,093.50 Pa**.

The weighted benchmark adds +2 MPa to one B pressure observation, assigns it
sigma=5 MPa, and assigns sigma=0.02 MPa to other surveys. Weighted T recovery is
closer to truth than the unweighted fit by more than a factor of ten in absolute
parameter error. Every weighted objective entry is checked against its own
pressure error divided by sigma, and the stored objective equals their squared
sum. Missing sigma requires an explicit default or individual exclusion.

Unweighted noisy fitting returns T=**1.184417119598998e-10**; weighted fitting
returns T=**7.999795453858116e-10**. The weighted fit's raw pressure RMSE is higher
(about 0.500 versus 0.452 MPa), because it deliberately gives the poor-quality
outlier little weight. Recovery and the selected weighted objective, not a claim
of universally lower unweighted RMSE, establish this benchmark's result.

## Safety and diagnostics

Per-tank and aggregate RMSE, MAE, bias and maximum pressure errors are retained.
Metadata, excluded surveys, asynchronous dates and effective default sigma remain
auditable. Global underdetermination fails; active B parameters with only A
observations receive a strong caution. This is data QC, not a claim of quantified
multi-tank identifiability.

Tests verify fixed-length finite penalties for pressure-bound, solver, aquifer,
stiffness and closure failures, and no candidate state commit. All-failed fits
cannot converge or apply. Repeated objective evaluations reproduce residuals and
complete aquifer/connection state while leaving the base network unchanged.
Closure-margin tests distinguish NEAR TOLERANCE without loosening acceptance.

Across the reproducible valid matching candidates:

- Maximum individual tank relative residual: **1.6473222785862163e-10**.
- Maximum individual tank absolute residual: **3.5183594926380124e-10 m³**.
- Maximum network transfer-conservation error: **0.0 m³**.
- Minimum relative closure margin: **9.835267772141379e-9**.

These are maxima for the reported fit benchmarks. The earlier Phase 5A cases
still reproduce their original values, including the near-tolerance 8.42e-9 case;
the new figures do not imply changed equations or relaxed acceptance criteria.

Signed TDI=X/Fprod is independently checked for both import and export, for
undefined zero-production cases, and for signed drive closure. Field transfer
contribution remains zero within conservation tolerance. The existing drive
equations and their normalization are unchanged.

FIELD/SI tests convert PVT, tank properties, histories, pressure/sigma metadata,
aquifer parameters, T, and physical parameter bounds/scales, then compare fitted
N, m, aquifer J and T, predicted pressures, transfer histories and objective values.
UI tests cover scenario inspection/application, actual matching with saved bounds,
unit-switch preservation of explicit default sigma, and saving allocation without
altering the base network. The application table identifies current versus fitted
values and exactly which active parameters will change.

## Regression and reproduction

The Phase 5B acceptance runner recomputes the complete Phase 5A record and its
recursive earlier acceptance cases, comparing outputs exactly. A SHA-256 manifest
protects the original domain, PVT, aquifers, MBE, single-tank solver/matching/
uncertainty, and the complete Phase 5A network implementation.

```powershell
.\.venv\Scripts\python.exe -m pytest -q -p no:cacheprovider
.\.venv\Scripts\python.exe examples/run_phase5b_acceptance.py
```

The [numerical record](phase_5b_numerical_results.json) stores all recovery cases,
per-tank/network metrics, objective values, bounds, closure margins, evaluation
counts, failure counts and protected hashes.

Base networks are not overwritten automatically. Applying parameters is explicit;
current histories, allocation and observations remain unchanged. No allocation
fractions were automatically fitted. No multi-tank uncertainty/correlation/surface,
Bayesian/Monte Carlo, forecasting or phase-transport feature was implemented.
