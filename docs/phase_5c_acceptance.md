# Phase 5C acceptance

**608 tests passed**: all **571 existing tests retained unchanged**, plus
**37 new Phase 5C tests**. The final combined run completed in 2581.17 seconds;
the complete output is retained in [the test log](phase_5c_pytest.log).

The [engineering method](network_identifiability_phase_5c.md) documents network
residual sensitivity, scaling, SVD/rank, covariance, profiles, surfaces, observation
coverage and deterministic allocation robustness. All diagnostics call the actual
Phase 5B objective and Phase 5A simulator. Shared Phase 4C mathematical utilities
are reused without changing their original implementation.

## Independent engineering benchmarks

The direct-history cases reuse the independent two-by-two compressible-storage
recurrence established in Phase 5B. They do not generate observations through the
production matcher. The aquifer case independently advances finite Fetkovich
storage and pressure. True N_A=1,000,000 m³, N_B=1,500,000 m³,
T=8e-10 m³/(Pa·s), and aquifer J=1e-9 m³/(Pa·s).

| Case | Information result |
|---|---|
| Both tanks observed with meaningful changing gradients | Full rank 3/3; maximum absolute correlation about 0.561; well constrained locally |
| Early Tank A surveys, N_A and T active | Correlation +0.985; poorly constrained despite near-zero fit error |
| Early Tank A surveys, B aquifer J and T active | Correlation +0.951; aquifer/communication compensation, poor constraint |
| Tank B unobserved, N_B and T active | Indirect information explicitly reported; local intervals withheld |
| Equal tank pressures throughout history | T sensitivity zero, rank 0/1; connection poorly observable |
| Redundant N_A and m_A storage terms | Rank 1/2; NON-IDENTIFIABLE; covariance intervals withheld |

Physical T derivatives are compared against central differences of the independent
storage oracle. A separate analytic Jacobian has JᵀJ=[[2,1],[1,1]], whose inverse
is [[1,-1],[-1,2]] and covariance correlation is −1/√2. This independently tests
the sign/magnitude convention rather than using the production function to supply
expected correlation.

Conditional surfaces have a distinct minimum at truth in the informative case.
Real forward evaluations along weak versus transverse scaled directions test the
elongated valleys in both confounded cases. Re-optimized profiles explicitly record
compensating N or J values and lower objectives than fixed-others slices. Their
plausible grids include the true T. Coarse single-node accepted regions are reported
as unresolved widths, with crossing brackets; they are never called zero uncertainty.

## Oil allocation benchmark and drift

A new independent recurrence uses true **60/40 oil allocation**, **20/80 water
production allocation**, and injection into B. Field rates change through history.
Oil uses Rs=0 and Bo=1.2+c(Pi−P), with c=4e-9 Pa⁻¹. The oracle accounts for
pressure-dependent oil withdrawal with storage D_i=C_i−c Np_i at each step:

```
(diag(Dnew) + dt L/2) Pnew
  = (diag(Dold) − dt L/2) Pold
    + (Dnew−Dold) Pi − Δ(1.2 Np + Wp − Winj)
```

Finite aquifer terms are independently added in the companion benchmark.
The oracle uses two-day steps and retains cumulative field samples at each step.
No surrogate replaces the production simulator during diagnostic evaluation.

Engineer-selected oil cases 50/50, 55/45, 60/40, 65/35 and 70/30 are each matched
separately for N_A, N_B and T. Companion aquifer cases quantify J drift as well.
The 60/40 case retains the known values and smallest pressure error. The other
cases quantify parameter drift and fit degradation without selecting or fitting
allocation. A separate test leaves oil unchanged and shifts water allocation to
verify independent stream behavior.

Recorded oil-allocation results (N in m³, T in m³/(Pa·s)):

| Oil A/B | Fitted N_A | Fitted N_B | Fitted T | Network RMSE |
|---|---:|---:|---:|---:|
| 50/50 | 852,527 | 1,760,624 | 6.72883e-10 | 9.709 kPa |
| 55/45 | 928,520 | 1,626,269 | 7.30885e-10 | 5.038 kPa |
| 60/40 | 1,000,000 | 1,500,000 | 8.00000e-10 | 1.78e-8 Pa |
| 65/35 | 1,066,283 | 1,382,922 | 8.81127e-10 | 5.357 kPa |
| 70/30 | 1,126,918 | 1,275,758 | 9.74304e-10 | 10.978 kPa |

T drifts from −15.89% to +21.79%, N_A from −14.75% to +12.69%, and N_B from
−14.95% to +17.37%. These selected cases classify as **MODERATELY SENSITIVE**
under the documented maximum-drift thresholds. No parameter bound is active in
these five fits. The result quantifies dependence on the allocation assumption;
it is not a statistical interval or automatic identification of true allocation.

With B aquifer J additionally active, the 50/50 case fits J=1.61381e-9
(**+61.38%**) and T=6.76123e-10, with RMSE **7.906 kPa**. The 70/30 case fits
J=3.81746e-10 (**−61.83%**) and T=9.08612e-10, with RMSE **3.987 kPa**.
The true 60/40 case recovers J=1e-9 and T=8e-10 with RMSE 1.87e-8 Pa.
Thus the aquifer-active cases are **HIGHLY SENSITIVE**: relatively small pressure
errors can coexist with substantially allocation-dependent water-support estimates.

At objective increment 3.84, the coarse re-optimized N–T profile accepts the
entire permitted T grid **[1e-11, 2e-9]**; both endpoints have objective below 1.24.
The aquifer–T profile accepts grid nodes **[8e-10, 1.005e-9]**, with crossing
brackets on each side. In contrast, the well-observed coarse profile rejects all
non-optimal coarse nodes. Its refined profile is retained in the numerical record
to resolve its narrower accepted range. All these grids hold allocation fixed.

The 21-sample refined well-observed profile (with the matched point inserted)
resolves accepted nodes from **7.065e-10 to 9.055e-10 m³/(Pa·s)**, containing
the true T=8e-10. Threshold crossings are bracketed by **[6.070e-10, 7.065e-10]**
and **[9.055e-10, 1.005e-9]**. These are grid-resolved plausible nodes and crossing
brackets at the specified objective increment, not interpolated exact limits.

Each scenario conserves field increments/cumulative totals through the existing
allocation audit. General-N tests transfer a fraction between two tanks while
holding a third unchanged, for all five streams and a selected effective interval.
Negative/out-of-range fractions, unknown tanks, duplicate identity and non-finite
shifts cannot become valid cases. Failed allocation cases remain UNRESOLVED.

## Numerical and workflow protection

- The full saved Phase 5B numerical record and all recursively earlier acceptance
  outputs were recomputed and compared **exactly unchanged**. See
  [the comparison record](phase_5c_prior_output_check.txt).
- A SHA-256 manifest protects all 571 existing tests and original source modules,
  including forward, matching, allocation and Phase 4C mathematical infrastructure.
- FIELD/SI checks round-trip PVT, reservoir/aquifer/connection properties,
  histories, field totals, observation pressure/sigma and physical bounds/scales.
  They compare physical/scaled sensitivity, correlation, statuses, profile/surface
  objectives and allocation-case fitted parameters.
- Repeated sensitivities, fixed and re-optimized profiles, surfaces and allocation
  refits preserve the original immutable scenarios and reproduce results.
- UI tests exercise information display, FIELD switching and preparation of five
  allocation cases without altering the accepted match.
- A final AppTest smoke check also rendered completed allocation results, the
  all-parameter inspector (including inactive properties), and FIELD display,
  verifying that the accepted scenario remained unchanged. The local application
  health endpoint returned `ok`.
- Invalid surface/profile nodes retain failure reasons and validity; they cannot
  bridge plausible regions. Rank deficiency, unstable derivatives, bounds and
  missing coverage withhold unsupported local intervals.

Reproduce the checks with:

```powershell
.\.venv\Scripts\python.exe -m pytest -q -p no:cacheprovider
.\.venv\Scripts\python.exe examples/run_phase5c_acceptance.py
```

The [numerical record](phase_5c_numerical_results.json) retains ranks, singular values,
correlations, weak combinations, coverage, profiles/surfaces and allocation-case
parameters, drift, metrics, drive contributions and transfers. The
[combined test log](phase_5c_pytest.log) records regression validation.

There is no automatic Apply operation in Phase 5C. Base networks, accepted matches
and allocation plans remain unchanged. Allocation fractions were **not** automatically
fitted. No Monte Carlo, Bayesian, ensemble or forecasting method was added.
