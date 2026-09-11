# Phase 3B - transient aquifer engineering specification

Phase 3B adds two radial transient aquifer models through the same immutable
interface introduced in Phase 3A:

- Carter-Tracy
- Original van Everdingen-Hurst (VEH), infinite-acting radial aquifer

Modified VEH, finite VEH boundary behavior, aquifer matching, forecasting and
multi-tank communication are not implemented in this phase.

## Shared radial transient layer

Both models use `RadialAquiferGeometry` in `src/material_balance_studio/aquifer/transient.py`.
All user inputs are converted to canonical SI before physics is evaluated.

The dimensionless time is

```text
tD = k t / (phi mu_w ct r_i^2)
```

where `r_i` is the reservoir/aquifer contact radius. The radial aquifer
constant used for edge-water influx is

```text
B = theta_rad phi ct h r_i^2
```

with `theta_rad = 2 pi theta_deg / 360`. `B` has units of reservoir m3/Pa.
Partial encroachment therefore scales the effective aquifer support by the
published radial angle factor. Geometry validation rejects radius ratio <= 1,
invalid porosity, nonpositive permeability, viscosity, compressibility,
thickness and angles outside `(0, 360]`.

## Original VEH

`VanEverdingenHurstAquifer` implements the original constant-terminal-pressure
superposition form for an infinite-acting radial aquifer:

```text
We(T) = B sum_j Delta p_j WeD(tD(T - t_j))
```

Each accepted step stores a pressure-step pair `(t_j, Delta p_j)` in immutable
`AquiferState.model_variables`. A root-solver candidate appends a temporary
current pressure step, evaluates cumulative influx, and returns a proposed
state. The accepted pressure history is committed only when the pressure solve
converges.

The infinite radial dimensionless cumulative influx response `WeD(tD)` lives in
`veh_response.py`. It uses classical VEH infinite-aquifer table nodes with
log-time/log-response interpolation, an early-time analytical approximation for
`tD < 0.01`, and the Edwardson infinite-aquifer large-time approximation beyond
the tabulated range. The implementation documents any large-time extrapolation
as an engineering warning.

The validated Phase 3B benchmark uses `tD = 1/day` and `B = 0.012566370614359175
m3/Pa`. For pressure steps 30 -> 28 MPa over day 1 and 28 -> 27 MPa over day 2:

| Quantity | Expected | Actual |
| --- | ---: | ---: |
| Day-1 cumulative We, m3 | 39433.27098785909 | 39433.27098785909 |
| Day-2 cumulative We, m3 | 81216.45328060335 | 81216.45328060335 |
| Day-2 current-step contribution, m3 | 19716.635493929546 | 19716.635493929546 |
| Day-2 historical contribution, m3 | 61499.8177866738 | 61499.8177866738 |

## Carter-Tracy

`CarterTracyAquifer` implements the published recursive transient-aquifer form:

```text
We_n = We_(n-1) + (tD_n - tD_(n-1))
       [B Delta p_n - We_(n-1) pD'(tD_n)]
       / [pD(tD_n) - tD_(n-1) pD'(tD_n)]
```

The pressure drop is `Delta p_n = p_i - p_n`. The model stores `tD`, `pD`,
`pD'`, the recurrence term, pressure drop and aquifer constant as diagnostics in
the proposed state. Candidate evaluations never mutate the accepted recurrence
state.

The pressure function uses the Carter-Tracy/Edwardson infinite radial
piecewise approximation:

```text
pD = (370.528 sqrt(tD) + 137.582 tD + 5.69549 tD sqrt(tD))
     / (328.834 + 265.488 sqrt(tD) + 45.2157 tD + tD sqrt(tD)), tD <= 100

pD = 0.5 [ln(tD) + 0.80907], tD > 100
```

For the same one-day, `tD = 1` benchmark and a 2 MPa pressure drop, the expected
first-step cumulative influx is `31331.833883943957 m3`, matching the model
output.

## Solver and diagnostics

The material balance remains

```text
Fnet(p) = N Et(p) + We(candidate)
```

No PVT, withdrawal, injection, expansion or residual definitions are changed.
The pressure solver still evaluates aquifer support through `AquiferContext`.
Repeated residual calls with the same candidate pressure return identical values.

The engineering inspector now includes transient diagnostics: `tD`, Carter-Tracy
`pD`, `pD'` and recurrence term, and VEH active pressure-step count plus current
and historical superposition contributions. Raw immutable variables remain in
advanced diagnostics.

## References

- A. F. van Everdingen and W. Hurst, "The Application of the Laplace
  Transformation to Flow Problems in Reservoirs", AIME, 1949.
- R. D. Carter and G. W. Tracy, "An Improved Method for Calculating Water
  Influx", Petroleum Transactions AIME, 1960.
- Edwardson-style infinite radial approximations as reproduced in reservoir
  engineering water-influx calculation examples and used only with documented
  range warnings.
