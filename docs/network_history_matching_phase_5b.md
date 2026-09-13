# Phase 5B — allocation and bounded network history matching

Three new layers surround the unchanged Phase 5A engine: `network_history` validates
and prepares historical inputs; `network_matching` constructs bounded candidates,
calls `simulate_network`, and records scenarios; the separate presentation module
provides allocation, pressure, parameter and equation inspectors. The original
single-tank and network equations, aquifers, PVT, pressure solver and timestep
refinement are not changed.

## Allocation and pressure input model

`HistoryPlan` contains explicitly selected field-allocated streams, immutable
field-total history, allocation rows, and the explicit choice to replace conflicting
direct stream values. Direct histories remain the default. Unselected streams
remain direct; selected streams are replaced, never added to direct values.
Conflicting nonzero direct histories require explicit replacement authorization.

Each `AllocationRow` contains an effective calendar date, one stream (`np`, `gp`,
`wp`, `winj`, `ginj`), and all tank/fraction pairs, including zeros. Every stream
has its own schedule beginning at the common initial date. Rows remain active
until the next effective date; fractions are never interpolated. Fractions must
be finite, between zero and one, and sum to one within 1e-10. Unknown/missing tanks,
duplicate dates/streams/tank pairs, and invalid sums are rejected. No renormalization
or automatic allocation fitting is performed.

The common timeline includes direct and field history, allocation changes, and
every pressure-observation date, including excluded surveys. The existing
Phase 5A cumulative-history interpolation supplies field totals at split dates.
For interval (a,b], allocate its field increment using the schedule active at a;
a change effective at b applies to the following interval. Each allocated
increment is added to that tank's prior allocated cumulative total. Applying a new
fraction to the entire historical cumulative volume would be incorrect and is
not done. This assumes constant field surface rates between cumulative samples,
with held totals after the final record, as in Phase 5A.

The audit retains start/end dates, tank, stream, field increment, active fraction,
allocated increment, resulting tank cumulative, and both interval and cumulative
conservation errors. Conservation uses absolute allowance 1e-8 surface m³ plus
1e-10 times the applicable field volume. Display/export converts each stream
according to liquid or standard-gas units. The long-form schedule input columns
are `date, stream, tank, fraction`. Field histories retain the established
`date,np,gp,wp` plus optional `winj,ginj,observed_pressure` schema; field pressure
is explicitly not allocated to compartments.

`NetworkObservation` preserves tank, date, pressure, source, sigma, QC flag,
include flag and note. Different tanks can have different survey dates. Positive
pressure and supplied sigma, valid tank references, and unique tank/date pairs
are required. Additional observation dates interpolate cumulative surface streams,
not pressure. Embedded tank pressure observations are imported by default;
an explicit observation table is authoritative. Excluded observations remain in
the result/inspection timeline but not the objective.

CSV/XLSX table inputs use SI **Pa**, not the UI's MPa display convention; FIELD
uses the established psia/STB/scf contract. Widgets and result tables display
selected engineering units. Input changes and parameter-bound edits must be
saved before evaluation/matching, avoiding stale-input runs.

## Parameter registry and optimizer

Stable identities are qualified by tank or stored edge orientation:

```
tank:A:N
tank:A:m
tank:B:aquifer:productivity_index
connection:A:B:T
```

Tank identifiers cannot contain colons under existing Phase 5A validation, so
these identities remain unambiguous even when tank names contain hyphens/spaces.
The registry reuses Phase 4B's approved per-tank parameters: N and m; Pot capacity;
Schilthuis J; Fetkovich J and initial connected water volume; validated transient
aquifer permeability, inner radius, thickness, porosity and encroachment angle.
It does not invent adjustable inactive radius-ratio or compressibility parameters.
Enabled connection transmissibility is additionally eligible. Disabled connections
cannot have active T matching because their T has no effect. No parameters are
automatically activated in the UI. Allocation fractions are absent from the registry.

Every active parameter has initial/lower/upper physical values, positive scale
and linear or log transformation, reusing Phase 4B's `MatchParameter`. Linear
coordinates are x/scale, log coordinates log(x/scale). T is nonnegative and linear
bounds may contain **exact zero**. Log bounds must be explicitly strictly positive;
no hidden positive floor substitutes for zero. Physical units are always displayed.
The approved per-tank domain restrictions and user bounds are enforced without
expansion for every candidate.

SciPy trust-region reflective bounded least squares uses a three-point numerical
Jacobian, parameter-scaled coordinates, `x_scale=1`, and ftol/xtol/gtol=1e-10.
The configured evaluation limit defaults to 60. This is the existing bounded
least-squares approach, with a network-aware objective and immutable candidate
constructor, not a surrogate network solver. A viable initial candidate remains
important: a flat region of failure penalties may prevent local optimization
from finding a valid region.

## Exact residual definition and data QC

Included observations are ordered by date, then tank name. Unweighted residuals
retain Phase 4B's fixed normalization: **(Pcalc−Pobs)/1 MPa**; weighted residuals
are **(Pcalc−Pobs)/sigma**. Q is the sum of squared objective residuals, while
RMSE/MAE/bias/maximum error remain physical pressure metrics. Weighted mode requires
valid sigma for every included survey or an explicitly supplied positive default;
there is no silent mixing of weighted and unweighted residuals. Both per-tank
and aggregate metrics are stored before optimization and at the final candidate.

Nobs≤Nactive is a pre-fit failure. Fewer than max(5,2Nactive+1) observations generates
weak-overdetermination caution. A tank with active parameters but no included
pressures receives a strong caution; connectivity does not establish direct
constraint. Low tank-level observation counts and connections with no observed
endpoints are also flagged. These are data-sufficiency cautions, not multi-tank
identifiability or correlation calculations. Source/QC metadata remain visible,
and a tank whose RMSE worsens by more than 1 Pa is flagged even if the overall
objective improves.

## Candidate validity, closure margins and reproducibility

The prepared histories, allocations and observations are fixed during matching.
Each decoded candidate replaces only its selected reservoir/aquifer/connection
parameters, then calls the existing `simulate_network` from a fresh initial state.
No pressure, cumulative transfer or aquifer state carries over between evaluations.

A candidate must converge and pass every original per-tank and network acceptance
check at every internal timestep. Matching rejects settings that loosen the
Phase 5A 1e-5 m³ absolute, 1e-8 relative, or transfer-conservation limits, or change
the normalization floor. Pressure RMSE cannot compensate for MBE failure.
Candidate failures return a deterministic, finite, fixed-length penalty vector:
1000 times the largest possible normalized residual within the fixed pressure
bounds, with a minimum penalty of 1000. Failed vectors are recorded separately
with parameters, failure category, message and available affected tank/edge details.
An all-failed optimizer termination is not convergence and cannot be applied.

For valid candidates, margin=1e-8−maximum individual relative tank residual.
The record reports PASS below 80% of tolerance, NEAR TOLERANCE at or above 80%,
and FAIL if any mandatory check fails. It stores the minimum margin across all
valid matching candidates, maximum individual absolute/relative residual and
maximum transfer error. No tolerance is relaxed. Users may explicitly select a
smaller maximum internal step before matching; the unchanged Phase 5A stiffness
refinement remains active. No candidate-specific hidden retry policy changes
the objective discretization.

The structured result stores the original and fitted network, topology, complete
history plan, pressure metadata, input fingerprint, specs/bounds, physical fitted
values, optimizer and evaluation status, before/after pressure records and metrics,
failure audit, closure statistics, warnings and allocation audit. Fingerprints
include reservoir/PVT/history, allocation and observation snapshots. Scenario IDs
also include settings, specifications, objective mode and fitted values.

Matching appends in-session scenarios. Applying a scenario is a separate explicit
action after a table shows current/initial/fitted values and bounds by tank or
connection. Only active parameters are applied to the current compatible topology;
current histories, observations and allocation schedules are preserved. Editing
the project does not rewrite an existing scenario.

## Drive and historical diagnostics

Existing Phase 4A support terms are reused and supplemented by **TDI=X/Fprod**.
Positive TDI means received support; negative TDI means export. DDI, SDI, CDI,
WDI, IDI and TDI remain signed and unnormalized. Their sum minus one is
−R/Fprod, where the denominator is suitable (>1e-12 reservoir m³); otherwise
indices are undefined. Individual curves are used, not a 100% stacked chart.
The field sums the physical contributions before division. Field TDI is therefore
zero within conservation tolerance, because transfer redistributes existing
support rather than creating network volume.

Per-tank before/after pressure and residual plots, per-tank and network metrics,
physical parameter tables, observed-pressure inspector, connection histories,
drive curves, and a date/tank equation inspector expose the match construction.
The inspector retains Fprod, injections, Fnet, each expansion support, We, X,
residuals, drive terms and individual neighbor contributions.

An allocation/pressure cross-check flags a fraction change of at least 0.20 near
(within 30 days of the subsequent survey) an observed pressure change of at least
5% Pi or a residual change of at least 1 MPa. These are documented engineering
screening thresholds. The flag does not infer causation or alter allocation;
well/completion-event metadata are not yet included.

No allocation optimization, multi-tank uncertainty, correlation matrices, objective
surfaces, Bayesian/Monte Carlo methods, forecasts or phase/component transport
were introduced. A good fit is not a claim of unique network parameters.
