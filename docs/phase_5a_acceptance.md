# Phase 5A acceptance

**533 tests passed in 125.75 seconds: all 489 existing tests retained unchanged,
plus 44 new network tests.** The suite includes independent analytical benchmarks,
candidate-state isolation and Streamlit UI interactions.

The dedicated multi-tank forward layer composes the existing black-oil balance,
PVT and aquifer components. It does not replace the single-tank solver. See
[the engineering formulation](tank_network_phase_5a.md) for integration, scaling,
state handling, history alignment, supported scope and limitations.

## Architecture and coupled calculations

The general-N architecture uses immutable tank nodes and explicit oriented edges;
acceptance exercises two tanks and a four-tank sparse network with an isolated
node. Every tank owns independent reservoir properties/history and references
its own or shared PVT. All seven aquifer types are verified against independent
single-tank runs at zero T; transient aquifers also have active-connection memory
tests. Candidate evaluations cannot mutate accepted aquifer or transfer state.

Each edge uses q=T(Pfrom−Pto), positive in its stored orientation, and
ΔV=T dt (ΔPold+ΔPnew)/2. Volume is applied with opposite signs at the two nodes.
The simultaneous bounded SciPy trust-region solver enforces
Fnet−N Et−We−X=0 for every tank. Deterministic internal refinement controls local
communication stiffness; all internal steps and failures remain inspectable.

## Independent engineering results

The closed benchmark embeds constant storage capacities C_A=0.004 and
C_B=0.008 reservoir m³/Pa in the existing black-oil expansion calculation.
Initial pressures are 30 and 20 MPa, T=1e-9 m³/(Pa·s), with no direct production,
injection or aquifer. The independent equilibrium pressure is 23.333333 MPa.
The exact discrete oracle updates pressure difference by
(1−lambda dt/2)/(1+lambda dt/2), lambda=T(1/C_A+1/C_B), and retains the
storage-weighted mean pressure. Every internal step agrees within 2e-5 Pa.

At 100 days the default-refinement pressures are 23.586264 and 23.206868 MPa.
With a 0.25-day maximum step, Tank A reaches 23.594421342 MPa versus the independent
continuous exponential solution 23.594425967 MPa: **4.6252 Pa difference**.
The refinement test shows improvement by more than a factor of 100.

| Case at 100 days | Result |
|---|---|
| Produced A, isolated | P_A=27.500000 MPa |
| Produced A, communicating B | P_A=28.671780 MPa; P_B=29.414110 MPa; B supplies 4,687.118536 reservoir m³ to A |
| Fetkovich aquifer on B only | Direct We_B=1,103.790498 m³; We_A=0; net B→A support=4,891.628401 m³ |
| Water injection into B | Direct injection remains in B; P_A=30.384690 MPa, P_B=32.557655 MPa; B→A support=11,538.758935 m³ |
| Flow reversal | Early A→B reverses after injection into B; final signed A→B cumulative volume=−6,220.991231 m³ |
| Large T=1e-8 | Both pressures approach 23.333333 MPa without forcing equality; residual pressure difference is numerical noise below 1e-6 Pa |
| Zero T | Maximum pressure difference from original independent single-tank runs is 3.3528e-8 Pa |

The produced-tank benchmark also checks each internal timestep against an
independent two-by-two linear storage solve including the known withdrawal.
Pairwise transfer conservation is checked explicitly, as is the sum over a sparse
four-tank network. Disabled links, equal-pressure zero flow, one edge reversing
direction, zero-history tank participation, separate PVT and aquifer assignment,
and exact previous-pressure initial guesses are covered.

Across the reproducible benchmark collection:

- Maximum individual tank absolute MBE residual: **8.418282959610224e-9 m³**.
- Maximum individual tank relative MBE residual: **8.418282959610224e-9**.
- Maximum network cumulative/incremental transfer-conservation error: **0.0 m³**.

The original per-tank limits remain 1e-5 m³ absolute and 1e-8 relative. Every tank
must pass both; offsetting tank residuals cannot pass by canceling in the field
sum. Failure tests verify pressure-bound failure, explicit extreme-stiffness work
limits and no invalid state commit or uncoupled fallback.

## Units, UI and protection

FIELD/SI equivalence tests round-trip reservoir properties, per-tank PVT/history,
aquifer parameters and transmissibility, then compare pressure, edge rate/volume,
aquifer influx and MBE residuals. Separate UI tests apply the displayed SI and
FIELD connection table and verify canonical T remains 1e-9 m³/(Pa·s). The UI run
test exercises tank addition, inspectors, unit switching and preservation of the
single-tank setup and stored network result.

The network workspace supports direct per-tank histories. An optional
programmatic fixed-share allocation helper is tested with 60/40 shares and
explicit streams; it conserves selected totals and rejects invalid shares. It
does not automatically assign field pressure observations or unselected streams.

The acceptance runner recomputes the entire Phase 4C record, recursively including
Phase 4B and all previous numerical acceptance results, and compares them exactly.
All match unchanged. SHA-256 hashes protect the original domain, PVT, aquifer,
MBE, solver, matching and uncertainty source files. No earlier tests were modified.

Reproduce:

```powershell
.\.venv\Scripts\python.exe -m pytest -q -p no:cacheprovider
.\.venv\Scripts\python.exe examples/run_phase5a_acceptance.py
```

[Numerical evidence](phase_5a_numerical_results.json) includes benchmark histories,
final component balances, edge states, QC maxima and protected source hashes.

The model exchanges lumped reservoir volume; it does not simulate phase/component
transport or changing PVT from mixing. No multi-tank history matching,
transmissibility fitting, uncertainty, time-varying allocation or forecasting was
implemented.
