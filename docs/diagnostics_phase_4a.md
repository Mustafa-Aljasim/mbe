# Phase 4A reservoir diagnostics

## Architecture and pressure basis

The new diagnostic modules consume accepted simulation states, history and the
existing balance API. They do not alter the domain records, PVT, withdrawal,
expansion, aquifer equations, pressure solver or time marching. Pressure metadata
is an immutable sidecar keyed by calendar date. The diagnostic import adapter
removes only the four documented metadata columns before calling the original
strict history parser. Original file formats remain valid.

VRR and drive indices describe the accepted calculated-pressure simulation.
Havlena–Odeh and Campbell instead evaluate the existing `evaluate_balance` API at
**observed** pressure. Using solved pressures would largely reconstruct the
assumed model and give misleading straight-line agreement. No observed pressure
is replaced with calculated pressure. Aquifer influx for these observed-pressure
transforms is replayed independently through the same immutable aquifer interface.
After a missing pressure, aquifer-adjusted transforms are unavailable for the rest
of that observed history; no silent interpolation or model-pressure substitution
is performed. Unadjusted Campbell terms remain available at later valid pressures.
NoAquifer requires no replay history. Out-of-PVT-range observations produce visible
FAIL rows; they do not extend the PVT table or modify the forward run.

## Voidage and VRR

Let Fp be validated produced reservoir voidage, Iw the validated water injection
support, Ig the validated gas injection support, and I=Iw+Ig. Define Fnet=Fp−I.
Totals come directly from `WithdrawalTerms`. Oil and water components are Np Bo
and Wp Bw. The signed free-gas component is obtained as Fp minus these two
components, preserving the validated solution-gas accounting. There is no second
implementation of the withdrawal equation. Injection FVFs retain the engine's
resident-fluid fallback when separate injection properties are absent.

Instantaneous (interval-average) VRR = I_interval/Fp_interval. Surface production
and injection increments are evaluated with the **endpoint accepted PVT** by the
existing withdrawal function. It is not a derivative of differently repriced
cumulative totals. Cumulative VRR = I_cumulative/Fp_cumulative, with all cumulative
surface quantities evaluated at the current accepted PVT, consistent with MBE.
Water and gas contributions use that same denominator separately and sum to VRR.
If produced voidage is ≤1e-12 reservoir m³, ratios are unavailable with a warning;
neither infinity nor a fabricated zero is returned. These definitions mean
cumulative reservoir voidage is not the sum of past interval reservoir volumes
when PVT changes. VRR=1 alone does not guarantee pressure maintenance.

## Drive/support indices

The engine equation Fnet=N(Eo+mEg+Efw)+We rearranges to
Fp=N Eo+N mEg+N Efw+We+Iw+Ig. Therefore every index uses **Fp**, not net withdrawal:

| Index | Definition | Contribution |
|---|---|---|
| DDI | N Eo/Fp | Oil and originally dissolved-gas expansion |
| SDI | N mEg/Fp | Initial gas-cap expansion |
| CDI | N Efw/Fp | Rock and connate-water expansion |
| WDI | We/Fp | External aquifer support |
| IDI | (Iw+Ig)/Fp | Operational injection support |

The sum minus one equals −residual/Fp. Values are not normalized or clipped to
force closure. Indices are unavailable for Fp≤1e-12 m³. The stacked percentage
chart uses 100 times the actual indices only when all plotted values are
nonnegative and sums close within 1e-8; undefined dates are omitted. Signed or
nonclosing histories retain individual curves and a visible caution. Efw includes
the existing (1+m) pore-volume factor. IDI equals cumulative total VRR under this
common basis; aquifer support remains distinct.

## Havlena–Odeh variants

The straight-line methodology originates with Havlena and Odeh (1963),
[SPE-559-PA](https://doi.org/10.2118/559-PA). The following transformations are
algebraic rearrangements of this application's validated component balance:

* No gas cap: X=Eo+Efw, Y=Fnet−We. Expected slope N, intercept zero.
* Supplied gas cap: X=Et=Eo+mEg+Efw, Y=Fnet−We. Expected slope N.
* Additional gas-cap view: X=Eg/(Eo+Efw), Y=(Fnet−We)/(Eo+Efw).
  Expected intercept N and slope Nm. Efw is conditioned on the supplied m, so this
  is not an independent automatic gas-cap estimate.

NoAquifer sets We=0. Active aquifers use observed-pressure replay as described
above; inferred volumes remain conditional on supplied aquifer parameters.
Only positive expansion denominators >1e-12 are plotted. Initial equilibrium and
nonpositive expansion points remain visible in the component table but have no
transformed coordinates. No user-specified OOIP or gas-cap ratio is overwritten.

Trend lines are unweighted mathematical regressions with an unrestricted
intercept. Three valid points are required. Centering/scaling the X axis improves
numerical conditioning. R² is unavailable for a constant Y series. QC warns about
<1% pressure depletion; axis spread ≤1e-12 times max(1,|X|); nonpositive inferred
OOIP; intercept >10% of max |Y| in the total-expansion view; R²<0.9; and >20% OOIP
change on omission of the first or last observation. With ≥4 points, material
curvature is flagged when linear RMS error exceeds 5% of Y spread and a quadratic
reduces squared error by >50%. These are screening thresholds, not proof of
geological validity. A high R² can reflect compensating assumptions.

## Campbell / Cole

For oil, Pletcher (2002), *Improvements to Reservoir Material-Balance Methods*,
equations 19–21, describes F/Et versus F; see the
[author's paper hosted by Texas A&M](https://blasingame.engr.tamu.edu/z_zCourse_Archive/P324_03A/Lecture_Refs_%28pdf%29/P324_Mod2_01_Pletcher_%28SPE_75354%29_%28Add%29.pdf).
Here injection is explicitly subtracted: X=Fnet, Y=Fnet/Et=N+We/Et.
An additional aquifer-adjusted curve uses (Fnet−We)/Et, expected to equal N under
the supplied model. Y is a stock-tank oil volume; X is reservoir volume. Et
includes the supplied gas cap and rock/water expansion. This is the oil Campbell
counterpart, not the dry-gas Cole formulation, which is inapplicable to this tank.

Endpoint changes within ±5% are described as approximately horizontal; larger
changes receive upward or downward cues. The wording identifies endpoint behavior,
not a claim that the entire series is monotonic. Scatter, PVT, pressure quality,
gas-cap assumptions and aquifer parameters must also be reviewed. No reservoir
drive mechanism is classified automatically from curve shape.

## Pressure QC and matching

The raw-history screen precedes strict parsing, preserving invalid rows for
review. FAIL flags cover invalid/nonchronological/duplicate dates, negative or
nonfinite cumulative volumes, decreasing cumulative production/injection,
production/injection without positive elapsed time, and nonpositive/invalid
absolute pressures. CAUTION covers missing pressure, fewer than five valid
observations, pressure above Pi, depletion <1%, survey gaps >365 days, and
survey changes >20% Pi. Thresholds are documented engineering review triggers.
FAIL history blocks building diagnosis. Existing strict parser/solver validation
continues independently; no points are deleted or repaired automatically.

Optional upload columns are `pressure_source`, `pressure_sigma`, `pressure_quality`,
and `pressure_note`. Sigma is in the file pressure unit (Pa in SI, psi in FIELD).
Supported source labels are Average reservoir pressure, Static well pressure,
PBU-derived pressure, RFT/MDT pressure, Estimated pressure and Other. Non-average
sources prompt representativeness review; unrecognized labels are retained and
flagged. Nonpositive/invalid sigma is flagged and excluded from normalization.
Metadata is retained with the saved run and displayed with pressure residuals.

Residual e=Pcalc−Pobs; RMSE=sqrt(mean(e²)); MAE=mean(|e|); bias=mean(e);
maximum=max(|e|). Missing observations are excluded, not converted to zero.
Where sigma>0, normalized residual=e/sigma. Metrics and trend lines are unweighted;
sigma is never used for optimization. Saved diagnosis, comparisons and simulation
snapshots are identified separately; editing inputs requires rerunning them.

No reservoir/aquifer parameter fitting, multi-parameter optimization, uncertainty
analysis, forecasting or DCA is implemented in Phase 4A.
