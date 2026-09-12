# Phase 5A — multi-tank forward engineering method

Each compartment is a `NetworkTank` referencing an unchanged `ReservoirTank`,
its own immutable cumulative history, and optional pressure bounds. A `TankNetwork`
contains an ordered tuple of tanks and an explicit tuple of `Connection` edges.
Any sparse topology is represented directly; isolated nodes do not exchange fluid.
The engine is general N, using a dense coupled solver. Acceptance exercises two
and four tanks; very large networks have not been performance-qualified.

All tanks share an initial reference date in Phase 5A, but may have different
initial pressures, OOIP, m, Swc, cf/cw, PVT and aquifers. Different initialization
dates are rejected rather than given ambiguous pre-initialization communication.
PVT objects may be shared by reference or independently assigned; the existing
PVT interface is unchanged. All seven validated aquifer choices use their existing
pure `initial_state` and `compute_step` interfaces independently for each tank.

## Transfer equation, physical scope and conservation

For an edge oriented from i to j, positive q and cumulative V mean i → j:

```
q_ij = T_ij (P_i - P_j)
ΔV_ij = T_ij Δt [(P_i_old - P_j_old) + (P_i_new - P_j_new)] / 2
V_ij_new = V_ij_old + ΔV_ij
X_i += -V_ij_new
X_j += +V_ij_new
```

T is reservoir m³/(Pa·s), pressure is absolute Pa, and time is seconds. SI display
uses reservoir m³/(day·MPa); FIELD uses rb/(day·psi). Endpoint rate is distinguished
from interval-average rate and incremental volume. When pressure order reverses,
q changes sign on the same physical edge; an interval's integrated volume may
have a different sign from its endpoint rate. Zero T or a disabled edge gives
zero rate and zero volume. Self-connections, duplicate links in either orientation,
unknown tank names and negative/nonfinite T are rejected.

Each edge is computed once, then applied with exactly opposing signs. Net tank
support X is cumulative reservoir-volume influx, which can be negative for an
exporting compartment. This is the requested **lumped reservoir-volume exchange
model**, not a phase/component transport model. It does not update fluid composition,
phase inventories, saturation or PVT due to transferred fluid. With different
fluids/PVT sets it conserves the specified communication-volume ledger, not an
independently resolved oil/gas/water component mass balance. Those interpretations
must not be conflated; phase transport is outside Phase 5A.

## Coupled residual and solver

The new layer calls the unchanged `evaluate_balance` separately for each tank,
including that tank's candidate cumulative aquifer influx. It then subtracts X:

```
R_i = Fnet_i - N_i Et_i - We_i - X_i
Fnet_i = Fprod_i - water_injection_support_i - gas_injection_support_i
```

The original balance object is retained as component terms, including its
uncoupled residual; the network's explicit `TankState.residual` is the corrected
coupled residual used for acceptance. The single-tank balance is not modified or
misrepresented as closed when X is nonzero. Aquifer, injection and transfer are
separate terms in tables and the equation inspector.

All pressures are solved simultaneously with SciPy `least_squares(method="trf",
jac="3-point")`, using bounded trust-region reflective least squares. It supports
physical bounds and numerical Jacobians without replacing the earlier scalar
solver; see [SciPy's method documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html).
Pressure coordinates are P/1e6 and the fixed per-tank solve scaling is
max(abs(local_storage_slope)×1e6, abs(previous_Fnet), 1 m³). Jacobian column scaling
uses `x_scale="jac"`. The initial guess is each tank's previous accepted pressure.
Default optimizer limits are 150 evaluations and ftol/xtol/gtol=1e-14.

A successful optimizer termination alone does not establish conservation.
Up to four bounded Newton polishing steps can improve physical residual closure
after successful optimizer termination, using the same simultaneous equations.
Final acceptance always applies the unscaled engineering checks below. An initial
guess already satisfying all checks is accepted without redundant optimization.
Invalid candidates receive a high finite residual and retain failure reasons;
no invalid candidate can be committed. No uncoupled fallback is available.

Per-tank bounds intersect its supplied bounds with its PVT coverage and the
positive 1 Pa lower limit. Initial pressure must lie inside them. There is no
automatic upper cap at Pi: injection and communication can increase pressure,
with a validity caution if the constant-compressibility model is being extended.

## Timestep integration and stiffness

Trapezoidal communication can oscillate across equilibrium if a timestep is too
large for a stiff connection. Phase 5A therefore subdivides each event interval
deterministically. At each accepted state, a bound-aware finite difference of the
existing reservoir balance estimates positive local storage slope
C_i=d(Fnet_i−N_i Et_i)/dP. Direct aquifer storage is excluded from this estimate.
The communication rate bound is

```
lambda_bound = 2 max_i [sum_j(T_ij) / C_i]
dt <= communication_step_limit / lambda_bound
```

The default communication step limit is 0.5 and may not exceed 1. This is a
conservative spectral/Gershgorin bound for the linear positive-storage network;
it keeps the linear trapezoidal relaxation factors positive. In nonlinear
black-oil cases it is a **local stiffness heuristic**, not an error estimator or
universal accuracy proof. Smaller user maximum-step durations can be imposed.
Every accepted internal step is retained for inspection and conservation checks.
The numerical benchmark verifies both the exact discrete recurrence and convergence
toward an independent continuous exponential solution as the maximum step shrinks.

A connected tank with nonpositive storage slope is rejected with an explanatory
failure. A default limit of 4096 internal steps per event prevents uncontrolled
work on extreme stiffness; failure preserves earlier accepted states. No pressures
are forced equal, no transfer is clipped, and no negative numerical oscillations
are hidden. Large systems, highly nonlinear storage or long aquifer histories
require refinement studies; the UI reports actual substep counts.

## Histories and immutable state

The timeline is the union of all tank history/event dates after the common initial
date. A zero-volume dated record is a valid event for closed equilibration.
Cumulative Np/Gp/Wp/Winj/Ginj are linearly interpolated between supplied dates,
starting from zero at the reference date and held after the final record. This
explicitly assumes constant surface rates between cumulative samples. Exact input
dates retain exact submitted totals. Observed pressure is retained only at its
own date and never imposed on the model or interpolated. A tank with no change in
production still advances through aquifer response and communication.

Initial-date rows must have zero cumulative streams and their pressure observations
are retained on the initial state. Empty tank histories mean zero direct streams;
at least one network event date is required. All original monotonic-history
validation remains in use. For asynchronous histories, the relevant independent
zero-T comparison is against single-tank runs on the same aligned timeline:
adding aquifer timesteps can change numerical integration relative to a coarser
original history. The tested zero-T regression uses identical supplied dates.

`allocate_field_history(history, shares, streams)` is an optional programmatic
fixed-allocation helper. It requires explicit selected streams and shares summing
to 1 within 1e-10, with no silent normalization. Unselected streams and field
pressure observations are not allocated. Output consists of new per-tank histories,
not modifications to the base project. The UI uses direct per-tank histories;
time-varying allocation is not implemented.

Frozen dataclasses and tuples hold reservoir references, histories, tank states,
aquifer states and edge ledgers. Every solver evaluation reconstructs candidate
aquifer steps and transfer totals from the previous accepted network state. It
does not increment stored state. Only a fully accepted pressure vector commits
all tank and connection states atomically. Failed candidates remain separate;
prior accepted internal steps survive failure of a later event.

## Numerical QC and inspectors

Every accepted internal step must pass:

| Check | Default acceptance |
|---|---|
| Each tank absolute coupled residual | ≤1e-5 reservoir m³ |
| Each tank abs(R)/max(abs(Fnet),1 m³) | ≤1e-8 |
| Signed network residual | abs(sum R) ≤ N_tanks × 1e-5 m³ |
| Network relative closure | abs(sum R)/max(sum abs(Fnet),1 m³) ≤1e-8 |
| Cumulative and incremental sum of tank transfers | absolute error ≤1e-8 m³ **and** relative error ≤1e-12 |

Relative transfer error uses max(sum absolute tank transfer magnitudes,1 m³).
Closure of the field sum cannot conceal a tank imbalance. Endpoint rate,
last-internal-step volume, cumulative transfer, neighboring contributions, direct
injection, We, X and all expansion terms are accessible separately. X/Fprod is
reported as a signed transfer-support diagnostic and is undefined for zero
produced voidage; Phase 4A drive-index equations are unchanged.

Solver diagnostics retain initial pressure guesses, individual bounds, objective
evaluation counts, optimizer nfev, physical residual norm, maximum individual
relative residual and failure reasons. Failed solves list dominant tank residuals
and incident edge increments. Counts describe coupled objective evaluations,
including Jacobian calls; storage-slope estimates are additional component calls.

The separate Multi-Tank Forward workspace supports tank add/remove/configuration,
PVT sharing/import, independent aquifers, direct CSV/XLSX or edited histories,
explicit connection editing and complete results/inspectors. Applying edits is
explicit. Saved runs contain the network snapshot and remain separate from later
edits and the single-tank workflow. Connection tables and neighbor lists provide
the engineering connectivity schematic. Canonical-SI JSON exports retain all
internal states; FIELD/SI display conversion does not mutate them.

No multi-tank matching, transmissibility fitting, uncertainty, forecasting,
time-varying allocation, spatial saturation or relative permeability is implemented.
