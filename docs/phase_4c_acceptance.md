# Phase 4C acceptance

**489 tests passed in 98.21 seconds: all 461 existing tests retained unchanged,
plus 28 new tests.** The full suite includes Streamlit UI interaction checks.

Phase 4C adds an independent `uncertainty` package and an **Identifiability &
Uncertainty** workspace for existing converged match scenarios. The forward
physics and Phase 4B objective/optimizer source files are unchanged. All earlier
tests remain unchanged. Full methods, thresholds, statistical assumptions and
limitations are documented in [the engineering method](identifiability_phase_4c.md).

## Numerical and engineering benchmarks

Eight prescribed pressure observations are converted independently into cumulative
production using F=N Eo+We, with constant Rs and no free-gas production. No forward
simulation or history-match implementation generates the synthetic observations.
The exact Pot case has Eo=4e-9(Pi−P), We=C(Pi−P), hence
**F=(4e-9 N+C)(Pi−P)**. A small quadratic Bo contribution makes the near-confounded
case full rank while retaining strong compensation. These synthetic cases are
acceptance problems, not proposed reservoir-fluid correlations.

| Benchmark | Result |
|---|---|
| Single OOIP, true N=1,000,000 m³ | Fitted N=1,000,000.000000176 m³; rank 1/1; WELL CONSTRAINED |
| Known sigma=0.1 MPa, single OOIP | Local approximate marginal 95% interval 979,140.03–1,020,859.97 m³, containing truth |
| Near-confounded OOIP/Pot capacity | Essentially zero RMSE; rank 2/2; correlation **−0.9976700682**; both POORLY CONSTRAINED |
| Near-confounded re-optimized N profile | All sampled N from 500,000 to 1,500,000 m³ accepted by the approximate threshold; reaches bounds |
| Exactly confounded OOIP/Pot capacity | Rank 1/2; covariance and marginal intervals withheld; NON-IDENTIFIABLE |
| Exact compensation | C=0.005+4e-9(1,000,000−N); re-optimized profiles recover this independent algebraic relationship |
| Exact conditional surface | Multiple distinct N/C nodes form the zero-objective trade-off valley; each valid node is checked against an independent analytical pressure formula |
| Zero sensitivity gas-cap ratio | Constant Bg makes Eg zero; m sensitivity is zero, rank is deficient, no false tight interval |
| Restrictive N≤900,000 m³ | Fitted N≈900,000; BOUND LIMITED; backward finite differences; symmetric interval withheld |

The exactly confounded fit can retain N=800,000 m³ and C=0.0058 m³/Pa with zero
pressure RMSE, even though the generating values are N=1,000,000 and C=0.005.
This demonstrates why numerical match quality cannot establish uniqueness.

The single-OOIP profile has a distinct minimum near one million m³. At a coarse
seven-point full-bound grid, only the matched coordinate is accepted. Its profile
width is explicitly **unresolved**, with threshold-crossing brackets; this is not
a zero-width confidence interval. The local interval above supplies the separate
linearized estimate. The broad confounded profiles are therefore clearly
distinguishable without pretending coarse grid endpoints are exact confidence
limits. User-bound truncation is identified in both profile and local output.

## Independent verification

Tests differentiate the independent closed-form pressure expression
P=Pi−Boi*Np/[c*(N−Np)] to check dP/dN and recover the known-sigma covariance directly.
The controlled Jacobian [[1,1],[1,2],[1,3]] gives correlation
−6/sqrt(42), independently of the production correlation calculation.

Profiles re-optimize compensating parameters using the unchanged Phase 4B
optimizer. Fixed-others and re-optimized results differ as expected. The surface
oracle inverts P=Pi−Boi*Np/[c*(N−Np)+C] and checks every valid objective node.
Failure injection and low optimizer budgets verify that failed/non-converged
points never become plausible minima. Stencil failures stop sensitivity analysis.

FIELD/SI tests convert pressure, production, observational sigma and physical
parameter bounds, then compare sensitivity, local correlation, statuses, profile
objectives and real surface objectives. Raw derivative display uses the correct
pressure/parameter unit ratio. The noisy-survey test checks each weighted and
unweighted Jacobian against its own residual scaling, and confirms different
covariance/intervals and different physical sensitivities at their fitted solutions.

Repeated calculations agree and leave complete scenario snapshots unchanged.
The UI test opens an existing scenario without rematching, runs diagnostics,
profiles and surfaces, switches units, and checks that the base and accepted
match remain unchanged. Log transformations, excluded observations, safe grid
resolution, invalid Jacobians and scenario identity changes are also covered.

## Protected regression and reproduction

`run_phase4c_acceptance.py` executes the earlier acceptance runners and compares
the entire Phase 4B numerical record exactly, including the earlier Phase 4A/3C
results and protected physics hashes. It also checks SHA-256 hashes of every
Phase 4B matching source. These checks pass. The earlier maximum valid matching
candidate relative MBE residual remains **1.0743406159705994e-13**; the original
forward acceptance tolerance is unchanged.

The [machine-readable numerical record](phase_4c_numerical_results.json) contains
fitted values, SVD/rank, confidence results, profile nodes/statuses, surface nodes,
warnings and protected source hashes. Reproduce with:

```powershell
.\.venv\Scripts\python.exe -m pytest -q -p no:cacheprovider
.\.venv\Scripts\python.exe examples/run_phase4c_acceptance.py
```

No Bayesian inference, MCMC, Monte Carlo uncertainty, surrogate model, multi-tank
matching, forecast uncertainty or forecasting was implemented. Correlation and
confidence remain local approximations; grids and local re-optimization do not
prove global uniqueness.
