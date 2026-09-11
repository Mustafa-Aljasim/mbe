# Phase 4B bounded reservoir history matching

The matching package calls the unchanged `simulate` function for every candidate,
starting with a fresh initial aquifer state. It never duplicates MBE, PVT or aquifer
equations. Configuration, objective construction, weighting, optimizer and scenario
records live outside Streamlit. The UI configures, runs, reviews and explicitly
applies a scenario.

## Parameters and numerical method

The primary optimizer is SciPy `least_squares`, trust-region reflective (`trf`),
linear loss, three-point numerical differences and explicit finite bounds.
[SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html)
describes its vector residual/bounded least-squares interface. The local validated
runtime is SciPy 1.18.1. Solver termination tolerances are 1e-10 for ftol, xtol and
gtol; forward solver tolerances are unchanged. Maximum optimizer evaluations are
user-controlled (1–1000). The separately reported forward evaluation count also
includes numerical differences and initial/final verification.

Every MatchParameter records name, description, initial, lower, upper, engineering
quantity, scale, linear/log transformation and active state. Internal coordinates
are x/scale or log(x/scale); physical values and bounds are always displayed.
Log scaling requires positive bounds. Inactive parameters retain their base-model
values. Candidate construction uses immutable dataclass replacements. Bounds are
never expanded and invalid physical values are never silently clipped.

| Model | Approved active parameters |
|---|---|
| All | `oil_in_place` (N), `m` |
| NoAquifer | No aquifer parameters |
| Pot | `aquifer.capacity` |
| Schilthuis | `aquifer.productivity_index` |
| Fetkovich | `aquifer.productivity_index`, `aquifer.initial_water_volume` |
| Carter-Tracy, original VEH, modified VEH | `aquifer.permeability`, `aquifer.inner_radius`, `aquifer.thickness`, `aquifer.porosity`, `aquifer.encroachment_angle` |

Names have the existing engine meanings. Radius ratio is deliberately absent:
the validated transient responses are infinite-acting and radius ratio has no
effect on their influx. Initial pressure, all compressibilities, PVT calibration,
dates, observed pressure and production/injection are fixed. Positive physical
bounds are enforced, except m may include zero; porosity must stay below one and
encroachment angle cannot exceed 360 degrees.

## Observations and objectives

Only real nonmissing pressure observations enter the selectable list. Every
selection retains date, pressure, source, sigma, QC flag, include state and note.
Exact observation dates already belong to the simulation timeline. The backend
rejects altered pressures, duplicate dates and dates absent from the original
history; it does not interpolate or invent observations. Exclusions affect only
the objective, not the stored history or diagnostic availability.

Unweighted residual: r=(Pcalc−Pobs)/1e6 Pa. The fixed 1 MPa numerical scale is
independent of display units. Weighted residual: r=(Pcalc−Pobs)/sigma, dimensionless.
Weighted mode requires valid sigma for every included observation, or the user's
explicit default sigma for missing/invalid entries. No arbitrary default is
silently applied. Sigma below 1e-12 Pa is numerically unsupported.

The optimizer receives the fixed-length residual vector. Reported objective is
sum(r²); SciPy's internal cost is half this value. Weighted objective is a
chi-square-style diagnostic, not a posterior probability. RMSE, MAE, mean bias
and maximum absolute pressure error use included observations in physical pressure
units. Excluded points remain plotted separately. Normalized residuals retain
sigma information even in unweighted review when available.

## Failure handling and closure

MBE residual Fnet−N Et−We is distinct from pressure mismatch Pcalc−Pobs.
Every valid candidate must converge for the entire history and pass the existing
absolute (1e-5 reservoir m³) and relative (1e-8) closure tolerances at every step.
A low pressure error cannot compensate for bad closure.

Failed candidates return a constant finite penalty vector of the same observation
length. Each element is 1000 times max(1, the maximum possible scaled pressure
error over the PVT pressure bounds). Parameters, objective, validity, failure
reason and maximum valid MBE residual are recorded for every call. This trace
contains finite-difference probes and rejected trials, not accepted iterations;
it is not expected to improve monotonically. A failed final candidate is marked
FAIL even if the optimizer stops on its flat penalty. This deterministic policy
does not guarantee recovery from an entirely invalid initial/bound region.

## QC and scenario acceptance

Nobs ≤ Nactive prevents matching. Fewer than max(5,2*Nactive+1) observations,
pressure depletion/spread below 1% Pi and included source/quality concerns cause
cautions. Bound distance below 1e-4 of the physical parameter range flags the
specific lower/upper bound. Residual review flags bias above half RMSE, maximum
error above twice RMSE, opposing early/late mean signs, and included normalized
residuals exceeding three sigma. These are review triggers; no points are removed.

Matched diagnostics reuse Phase 4A equations, and aquifer engineering screening
reuses Phase 3C QC. The H–O comparison uses all available observations and the
matched aquifer, including observations excluded from the match objective; its
basis is labeled. More than 20% difference from matched OOIP prompts review.
H–O is not forced to agree. Scenarios do not establish unique parameter values.

Each in-session scenario retains parameter specifications, observations, weighting,
initial/final objectives and pressure metrics, calculated histories, residuals,
evaluation trace, closure, QC, recomputed diagnostics and a PVT configuration/
sample fingerprint. Different observation sets or weighting modes do not have
directly comparable objectives. Display-unit changes preserve stored scenarios
and edited parameter bounds.

Running a match does not overwrite the base reservoir. The explicit **Apply
matched parameters to reservoir model** action shows current and fitted physical
values first, then updates only selected setup parameters. It requires the same
aquifer type and a converged valid result. Simulation results remain their original
snapshot until the engineer reruns the forward simulation.

No stochastic/global optimization, uncertainty distributions, correlation matrices,
automated identifiability analysis, multi-tank matching, forecasting or DCA is added.
