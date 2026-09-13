# Phase 5C — network information and allocation sensitivity

Phase 5C adds network adapters under `uncertainty` and a **Network Identifiability**
workspace. It leaves the Phase 5A forward engine, Phase 5B matching objective,
and all existing source/test files in those layers unchanged. The Phase 4C SVD,
covariance geometry, correlation, local-confidence, classification, parameter-grid,
and profile-threshold/range utilities are reused. A diagnostic does not apply
parameters or alter the base network, accepted scenario, or allocation plan.

## Three distinct questions

1. **Fit quality:** per-tank and network pressure errors in a chosen match.
2. **Parameter identifiability at fixed allocation:** which parameter combinations
   the recorded pressure observations can distinguish locally or along profiles.
3. **Allocation assumption sensitivity:** how separate bounded matches change
   when an engineer specifies alternative allocation schedules.

A small RMSE answers only the first question. The deterministic allocation ranges
are not confidence intervals and the application never selects an inferred true
allocation. No Monte Carlo, Bayesian inference, ensemble, forecast or transport
calculation is added.

## Residual derivatives and parameter scaling

The network adapter reconstructs an immutable Phase 5B `NetworkObjective` from a
converged scenario, with fitted active parameter values, the original physical
bounds/scales, exact history/allocation plan, included observations, weighting,
sigma and network settings. Every evaluation starts the existing simulator fresh.
Observations retain date-then-tank ordering and unique parameter identities such
as `tank:A:N`, `tank:B:aquifer:productivity_index`, and `connection:A:B:T`.

Physical-coordinate differences calculate `J_ij = dr_i/dx_j`. The objective is
exactly Phase 5B: unweighted residuals are `(Pcalc-Pobs)/1 MPa`; weighted residuals
are `(Pcalc-Pobs)/sigma_i`. Physical pressure sensitivity is `dP/dx = J sigma_i`
in weighted mode or `J × 1 MPa` unweighted. The scaled residual Jacobian is
`J diag(parameter_scale)`. Dimensionless engineering sensitivity is
`(dP/dx) × parameter_scale / Pi`, using the **observation tank's** initial pressure
for each row. Raw sensitivities convert to displayed pressure/parameter units;
dimensionless sensitivities, rank, correlation and objective are unit invariant.

The initial stencil step is the smaller of `1e-4 × max(|x|,scale)` and 10% of the
allowed span. The engineer/API can select a relative step from 1e-6 through 1e-2.
Central differences are used inside bounds; second-order one-sided differences
are used near a bound. The step is halved for an independent stability comparison.
Failed stencils shrink deterministically up to four attempts; failure is reported,
never replaced by derivatives of the optimizer's failure penalty. More than 5%
derivative change with step halving withholds local confidence. Physical T=0
continues to work with linear bounds; there is no hidden logarithmic floor.

## Observation information and coverage

Rows of the information heatmap identify both tank and survey date; columns
identify active network parameters. Cell magnitude is absolute dimensionless
sensitivity. Per-observed-tank RMS magnitudes expose cross-tank responses.

Coverage uses RMS dimensionless sensitivity: below 1e-8 is **UNINFORMED**;
below 1e-5 is **WEAKLY INFORMED**. Above these thresholds, a tank parameter with
observations in its own tank is **DIRECTLY INFORMED**; otherwise it is
**INDIRECTLY INFORMED**. A connection requires observations at both endpoints for
the direct label. Missing tank pressure coverage is explicitly reported.
These are engineering screening labels, not proofs of identifiability. Indirect
or weak coverage conservatively withholds marginal local intervals and downgrades
an otherwise favorable parameter classification.

## SVD, correlation and confidence

The shared Phase 4C matrix utility computes SVD of the scaled residual Jacobian.
Singular values above `max(1e-10, 1e-6 × largest singular value)` determine rank.
Rank below the active parameter count is **NON-IDENTIFIABLE**. No pseudo-inverse
covariance or marginal confidence is displayed for a rank-deficient network.
The weakest right singular vector reports compensation coefficients in scaled
parameter coordinates. Its overall sign is arbitrary; its components are not
causal effects.

For full rank, covariance geometry is `V diag(1/s²) Vᵀ`. Correlation comes from
this covariance, not directly from Jacobian column cosines. When covariance is
unavailable, column coupling is explicitly labeled as a different diagnostic.
Absolute correlation ≥0.7 is **MODERATE COUPLING**, >0.9 **STRONG COUPLING**.
Convenience pairs cover adjacent N–T, adjacent aquifer–T, N–N and same-tank m–N;
any two active parameters remain available.

The original Phase 4C classification combines sensitivity, rank, coupling,
conditioning, observation count, bounds and within-tank pressure span. Conditions
above 1e4 or unstable derivatives withhold local confidence. Weighted confidence
requires accepted independent Gaussian errors with known sigma; unweighted
confidence estimates common residual variance with positive degrees of freedom.
Effectively noiseless unweighted residuals cannot estimate variance. Intervals
reaching a user bound are withheld rather than clipped. All intervals are
**LOCAL / LINEARIZED APPROXIMATIONS** under fixed allocation, not global uniqueness.

## Connection observability

For every edge, the observed-history timeline reports predicted pressure difference
and signed instantaneous transfer rate. `max |DeltaP|/Pi < 0.001` flags poor
excitation. Missing endpoint observations and history shorter than 30 days add
cautions. These screening thresholds complement the numerical T sensitivity and
aquifer/T correlation; they are not automatic proof that a connection is absent.

## Objective surfaces and profiles

Surfaces evaluate the actual Phase 5B objective on bounded physical grids with
other parameters held at their matched values. Linear/log grid spacing follows
the parameter specification, with the actual optimum inserted. Initial and matched
points, original bounds and failed nodes are shown. A coarse grid may miss a narrow
minimum. No surrogate or interpolation over failed regions is used.

Profiles either hold all other parameters fixed or re-optimize nuisance parameters
through the existing Phase 5B matcher, inside their original bounds. Each node
records all physical parameter values, bounds, objective, per-tank and network
RMSE, convergence, closure margin, conservation and failure details. Every node
starts from the same matched nuisance values, independent of traversal history.
Non-converged nuisance fits remain unresolved even when their forward run is valid.

Shared Phase 4C thresholds support a user-selected positive objective increment
(no probability attached), or conditional approximate likelihood/F thresholds
when statistical and numerical regularity permit. Statistical ranges are withheld
for rank deficiency, active bounds, instability, poor conditioning or missing direct
coverage. Disconnected grid-resolved regions remain separate. A single accepted
node has **unresolved width**, not zero uncertainty. Bound-reaching ranges have
unknown extent outside the allowed domain. A profile finding a better optimum than
the accepted match withholds ranges until that match is reconsidered.

UI grids offer 3/5/9 samples plus the inserted optimum; the API permits 3–21.
Expected surface forward counts and profile-point counts are shown. Re-optimized
profiles expose an optimizer budget; actual forward counts include finite-difference
Jacobian calls and final verification, so they can exceed optimizer `nfev`.

## Allocation perturbation and robustness

The operation transfers fraction `delta` **from donor to receiver**, leaving every
other tank unchanged. Negative delta reverses the transfer. It applies to one
selected stream and either all its effective rows or selected existing rows, whose
intervals continue through the next effective date. It works for two or general N
tanks. Impossible fractions are rejected, never clipped or renormalized. Each of
oil, gas, water production, water injection and gas injection can be investigated
independently. The UI also accepts complete engineer-defined schedule CSVs.

Every case goes through Phase 5B history preparation and conservation validation.
Field totals, direct/allocated mode and observations stay fixed. Only approved
reservoir/aquifer/connection parameters enter each independent match; allocation
fractions never enter its parameter vector. Up to 21 named cases can be compared.

For each active parameter, results retain base/fitted value, signed absolute change,
and percent change relative to the absolute base value. When `|base| <= 1e-8 × scale`,
percentage change is undefined; absolute drift divided by the original bound span
provides the documented alternative screening measure. The largest absolute drift
across parameters and selected cases gives:

- **ROBUST:** ≤10%.
- **MODERATELY SENSITIVE:** >10% through 25%.
- **HIGHLY SENSITIVE:** >25%.
- **UNRESOLVED:** any failed/non-converged allocation case.

This classification depends on the engineer's chosen plausible cases. It is not
a probability statement. Tables and plots show allocation versus fitted N, m, T
and aquifer parameters, per-tank/network fit degradation, bound warnings, final
signed drive indices and cumulative transfers. Failed cases remain visible with
failure records; they are never silently excluded to make robustness look better.

All records are immutable snapshots and can be exported as JSON. There is no Phase
5C Apply action, and the accepted Phase 5B match and allocation plan remain intact.
