# Phase 4C — deterministic identifiability

A good pressure match does not establish unique reservoir parameters. This layer
operates on a converged, immutable Phase 4B scenario and its included pressure
observations. Opening the workspace does not rematch the reservoir. Each button
runs explicit temporary calculations; only the diagnostic record is saved in the
session. The base reservoir and accepted scenario are never changed.

## Sensitivity and objective definitions

The existing `PressureObjective` is the only residual builder. Its unchanged
unweighted objective uses **r=(Pcalc−Pobs)/1 MPa**, not pressure numerically in Pa.
This constant normalization preserves the original Phase 4B optimum. Weighted
residuals are **r=(Pcalc−Pobs)/sigma**, including the scenario's original explicit
default-sigma policy. The objective is Q=sum(r²); pressure RMSE is reported
separately in the selected engineering units. Every grid in a scenario uses the
same objective; no switching or mixing weighting.

Physical residual derivatives J=dr/dx come from fresh forward evaluations at the
fitted physical parameters. Multiplying row i by its original residual scale
recovers S=dP/dx. Dimensionless pressure sensitivity is
**Sscaled[i,j]=S[i,j] × parameter_scale[j] / Pi**. Parameter scales are the explicit
positive Phase 4B scales, even for log-transformed optimization. Derivatives and
reported intervals are physical, not derivatives with respect to log parameters.

Default perturbation h=min(1e-4 × max(abs(x),scale), 0.1 × bound span). Central
differences are used when feasible, otherwise second-order forward/backward
differences. Calculations at h and h/2 check step stability; the finer result is
retained. Failed stencils retry up to four smaller steps, never outside bounds.
An unavailable stencil stops the analysis with a named failure, rather than
turning a penalty into a derivative. A relative derivative change above 5%
withholds confidence estimates. API relative steps may range from 1e-6 to 1e-2.

Sensitivity ranking uses RMS dimensionless sensitivity across included survey
dates. Raw derivatives remain visible. Display conversion uses
Sdisplay=Sphysical × pressure_conversion_factor / parameter_conversion_factor.
The canonical records and conclusions do not change with FIELD/SI display units.
Time plots describe existing information only, not future survey optimization.

## Rank, correlation and engineering QC

SVD is performed on **A=J diag(parameter_scale)**. A singular value is retained
when greater than max(1e-10, 1e-6 × largest singular value). Full-rank condition
number is largest/smallest; rank-deficient condition is explicitly unavailable
(unbounded), never a misleading finite number. The weakest right singular vector
is exposed in scaled parameter coordinates; its sign is arbitrary.

For effective full rank, G=V diag(1/s²) Vᵀ, and correlation is
Gij/sqrt(Gii Gjj). This is local parameter covariance geometry, not cosine
similarity of sensitivity columns. Column cosine coupling is available separately
in Advanced Diagnostics. No pseudo-inverse marginal covariance is reported when
rank deficient because it could assign false precision to unconstrained directions.

Classification follows these ordered engineering rules:

| Evidence | Status |
|---|---|
| Within 1e-4 of the allowed span from a bound | BOUND LIMITED |
| Effective rank deficiency or RMS scaled sensitivity <1e-8 | NON-IDENTIFIABLE |
| Maximum absolute parameter correlation >0.9, RMS sensitivity <1e-5, or condition >1e4 | POORLY CONSTRAINED |
| Correlation ≥0.7, fewer than max(5,2p+1) observations, or pressure spread <1% Pi | MODERATELY CONSTRAINED |
| Remaining cases | WELL CONSTRAINED, locally only |

A plausible profile region spanning more than half the permitted parameter range
downgrades an otherwise stronger status to POORLY CONSTRAINED. Disconnected
grid-resolved regions receive a caution. Overall rank deficiency conservatively
flags every active parameter; this is not a proof that each individual coordinate
is unestimable. Removing an inactive/zero-information parameter and rematching
may isolate the estimable coordinates. Diagnostic thresholds are documented
engineering heuristics, not universal statistical laws or evidence of causation.
Extreme column-scale imbalance above 1e6 is also flagged.

## Profiles and surfaces

Two-dimensional surfaces evaluate the real unchanged forward objective at each
pair of physical parameter values. All other parameters remain at fitted values.
The UI offers an OOIP versus applicable aquifer-parameter shortcut, while any
two active parameters may be selected. Initial and matched coordinates, complete
bound axes, failed-node markers, objective and pressure RMSE are shown.

Profiles have two distinct modes: conditional fixed-others slices, and profiles
that re-optimize all remaining active parameters with the unchanged Phase 4B
trust-region optimizer. Each node starts from the accepted fitted values for the
remaining parameters, with their original bounds/scales/transforms. No warm-state
or aquifer-state carryover is used. Compensating parameter estimates, bounds,
convergence, objective, RMSE and failure information are retained at every node.
Non-converged valid nodes are visible in tables but are not treated as confidence
limits. Failed forward nodes have null objective/RMSE, never an attractive penalty.

Resolution is 3–21 nodes per axis, plus the exact fitted coordinate when absent.
Log-transformed parameters receive geometric grids; others linear grids. The UI
shows expected calculation cost. All evaluations are serial. Coarse grids can
miss narrow minima or plausible regions; reported ranges are **accepted sampled
nodes**, not interpolated confidence endpoints. Failures and unresolved nodes
break regions. Boundary-touching regions are truncated by user bounds and may
extend farther in reality. These calculations do not prove global uniqueness.
Threshold-crossing brackets show where an endpoint could lie between evaluated
nodes. A region with one accepted node has unresolved width, not zero uncertainty.
Near-duplicate grid coordinates are merged into the exact fitted coordinate.

## Approximate confidence and plausible ranges

Statistical assumptions default to **not accepted**. An engineer may explicitly
accept independent Gaussian observational errors, suitable variance information,
adequate observations, and local linearity. This action records an assumption;
the software does not establish it from the data. The distinction between known
absolute sigma and variance estimated from residuals follows the
[SciPy covariance documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html).

For known weighted sigma, physical covariance is D G D and local marginal 95%
intervals are x ± 1.959964 SE. For unweighted fits with estimated common variance,
covariance is D G D × Q/(n−p), with t(.975,n−p) intervals. D contains parameter
scales. All are **LOCAL / LINEARIZED APPROXIMATIONS**, not exact nonlinear bounds.
Rank deficiency, unstable derivatives, condition >1e4, insufficient degrees of
freedom, and effectively noiseless unknown-variance fits with Q≤1e-12 withhold
covariance/intervals as appropriate. An interval crossing a bound is withheld,
not clipped; bound-limited estimates never receive invalid symmetric intervals.

For a regular known-sigma re-optimized profile, the approximate marginal 95%
threshold is Qmin + chi-square(.95,1). This is the profile likelihood-ratio
construction described in
[MathWorks parameter confidence documentation](https://ch.mathworks.com/help/simbio/ref/sbioparameterci.html).
For estimated common variance, the approximate extra-sum-of-squares criterion is
(Qprofile−Qmin)/(Qmin/(n−p)) ≤ F(.95,1,n−p). This is a linear-model F construction
used here as an explicitly approximate nonlinear criterion; it is not an exact
nonlinear confidence statement. Both require the accepted statistical assumptions.
Rank/conditioning, derivative or bound regularity failures disable statistical
profile thresholds. A profile finding a better optimum than the accepted match
also withholds ranges until that mismatch is investigated.

An explicit positive **plausible objective increment** is available without
assigning any confidence probability: accept nodes with Q≤Qmin+increment.
This is useful for noiseless, rank-deficient and bound-limited examples. Fixed-others
ranges are labeled conditional and can understate uncertainty through omitted
compensation. Profile optimization is local and can miss alternate minima.

The result retains scenario identity, objective mode, active parameters, matrices,
rank, singular values, stencils, step checks, uncertainty, warnings, profiles,
surfaces and failure counts; download it as canonical-SI JSON. Scenario identity
includes reservoir/PVT, complete cumulative production/injection history,
observations, weighting, fitted parameters and bounds.

No Bayesian posterior, MCMC, Monte Carlo, proxy, forecast uncertainty, global
optimization or future-survey design is implemented in this phase.
