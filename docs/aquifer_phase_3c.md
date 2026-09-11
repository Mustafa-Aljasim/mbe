# Phase 3C formulation and scope

## Formulation decision (before implementation)

Modified Van Everdingen-Hurst (linear-pressure-history formulation) follows
Petroleum Experts, **MBAL Reservoir Engineering Toolkit User Guide**, Appendix C,
section C2.9, equations 1.6a–1.6d. The publisher-authored guide is available in a
[public mirror](https://studylib.net/doc/25780905/mbal-complete).
Section C2.6 also identifies the related Vogt–Wang convolution treatment;
this implementation uses the explicit C2.9 piecewise-linear convolution, not
Hurst Modified Steady-State, and does not claim binary equivalence to MBAL.

Let a = k/(phi mu ct ri²), tD = a t, B = theta phi ct h ri² (theta in radians),
and R(u) be the existing dimensionless cumulative VEH step response.
On accepted segment [tj-1,tj], pressure is
p(t) = pj-1 + (pj-pj-1)(t-tj-1)/(tj-tj-1).
With positive pressure drop dj = pj-1-pj, cumulative influx at T is

    We(T) = B sum_j dj/[a(tj-tj-1)] * integral[a(T-tj), a(T-tj-1)] R(u) du.

Incremental influx is We(T)-We(previous T). The original implementation instead
uses B sum_j dj R(a(T-tj-1)), placing each full drop at the interval start.
Both approach the same convolution as steps shrink for a smooth pressure history;
neither is automatically physically preferable when actual within-step pressure
is unknown. Rate is not constant even for a constant-pressure VEH step.

The modified integral uses exact power-law antiderivatives of the shared log-log
table interpolation, the analytic small-time integral, and adaptive quadrature
only beyond the table. No duplicate response table is introduced.

Supported geometry is homogeneous, single-phase, infinite-acting radial edge
water with constant properties and partial encroachment scaling. Neither VEH
model implements finite boundaries. Radius ratio is retained as geometry metadata
and has **no influence** on infinite-acting influx; it is not offered as an active
sensitivity input. The inherited response table/approximations limit accuracy;
integration accuracy is not independent validation of the full response table.
An infinite aquifer has no calculated lumped average pressure: comparison reports
that quantity as unavailable, rather than treating reference pressure as depletion.

Each accepted state contains immutable (start, end, pressure drop) segments.
Candidate evaluation appends only a temporary segment, evaluates the complete
convolution, and returns a proposed state through the existing solver interface.
No solver, withdrawal, expansion, injection, PVT, or existing aquifer formula changes
are required.

## Comparison methodology and diagnostic metrics

Each model runs the unmodified chronological pressure solver with an immutable
replacement aquifer configuration on the same reservoir, PVT and cumulative
production/injection records. Observations do not overwrite predicted pressure.
Errors use calculated minus observed pressure at matching recorded dates; missing
observations are omitted, without interpolation. RMSE = sqrt(mean(e²)),
MAE = mean(abs(e)), bias = mean(e), maximum error = max(abs(e)). Observation counts
are shown. Runs that fail are visible with partial results and are excluded from
best numerical fit selection. Complete runs therefore share the observation basis.

Comparison controls maintain separate SI parameter drafts per model. Saved results
retain their model, reservoir and history snapshots; changing widgets requires a
new run, while changing display units only converts presentation. Tables and plots
use MPa or psia, reservoir m³ or rb, and the shared rate conversion. The inspector
compares the same calendar date across models.

The provisional support fraction is We / (Fprod - Winj Bwinj - Ginj Bginj).
It is unavailable when net withdrawal is less than 1 reservoir m³, including
negative or zero net withdrawal. Signed We is preserved; no clipping to [0,1].
All original withdrawal, injection and expansion terms remain available for future
drive indices. This is not the final WDI framework.

QC uses FAIL for failed pressure solves or relative closure above 1e-8, and
CAUTION for response approximation, signed-flow warnings, depletion with falling
We, missing observations, RMSE above 5% of initial pressure, or support fractions
outside [0,1.05]. Radial models always request confirmation of the infinite-acting
geological assumption; permeability outside approximately 0.1–10000 mD is flagged.
Constructor validation rejects invalid geometry/parameters before simulation.
Thresholds are explicitly screening triggers, not universal geological bounds.
PASS means those automated checks passed; geological size, connectivity and
property justification remain the engineer's responsibility. A low RMSE can result
from compensating, unjustified parameters, so numerical fit never establishes
physical correctness.

Sensitivity runs use 2–8 distinct user-supplied values and immutable model copies.
Fetkovich supports productivity and connected volume; CT and both VEH variants
support permeability and encroachment angle. No ranking, optimizer, search,
history matching, forecast or uncertainty analysis is added by this preview.

## Independent benchmark derivation

The three-segment benchmark applies C2.9 to the classical early-time response
R(u)=2 sqrt(u/pi), yielding the hand antiderivative
I(l,h)=4/(3 sqrt(pi)) (h^(3/2)-l^(3/2)). This is a constructed reference-equation
benchmark, not a claim to reproduce a published field-data table or MBAL output.
The underlying constant-terminal-pressure solution is from van Everdingen and
Hurst (1949), [original paper record](https://www.onetunnel.org/documents/reservoir-engineering-application-of-the-laplace-transformation-to-flow-problems-in-reservoirs).

Inputs: ri=1000 m, h=10 m, phi=.2, ct=1e-9 /Pa, mu=.001 Pa·s,
full encroachment, k=.2*.001*1e-9*1000²/86400 m², hence a=1/day and
B=.01256637061435917 m³/Pa. Start at 30 MPa and use endpoints
(tD=.001,29 MPa), (.003,27 MPa), (.006,26 MPa). All response ages are below .01.
Sum B*dj/(endD-startD)*I(TD-endD,TD-startD) independently for every segment.

| tD | Expected cumulative We (m³) | Expected incremental We (m³) |
|---|---:|---:|
| .001 | 298.9328648745563 | 298.9328648745563 |
| .003 | 1553.3007300445597 | 1254.3678651700034 |
| .006 | 3357.8640977163855 | 1804.5633676718257 |

The test calculates these expectations solely from the closed-form equation and
input constants. A separate adaptive-quadrature cross-check covers integration
through the existing table and extrapolated tail; that check verifies integration,
not the accuracy of the inherited response data themselves.
