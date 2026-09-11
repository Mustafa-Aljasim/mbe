# Phase 3A — stateful aquifer engineering specification

The supported models are **None, Pot, Schilthuis and Fetkovich**. All use one
interface and the existing chronological pressure solver. PVT functions, reservoir
expansion formulas and production/injection withdrawal formulas are unchanged.
The cumulative balance is now explicitly:

```text
Fnet(p) = N Et(p) + We(p, previous aquifer state, dt)
Residual(p) = Fnet(p) − N Et(p) − We(...)
```

We is **reservoir water volume**, not surface water volume. It is added directly
to reservoir expansion support, with no further multiplication by Bw. Positive
We denotes influx into the oil tank; negative increments denote efflux. The
existing normalization remains `|Residual| / max(|Fnet|, normalization_floor)`.
Both absolute and relative convergence tolerances must pass.

## Architecture and state ownership

`src/material_balance_studio/aquifer/` contains `base.py`, `state.py`, `none.py`,
`pot.py`, `schilthuis.py`, `fetkovich.py`, and `registry.py`. Validation is shared
in `base.py` and the immutable records, with parameter checks in each model,
instead of a separate validation module.

`AquiferModel` defines `initial_state(initial_reservoir_pressure)` and
`compute_step(previous_state, previous_reservoir_pressure, candidate_pressure,
dt_seconds)`. Every implementation returns `AquiferStepResult`. ReservoirTank
holds a model configuration; omission uses the explicit NoAquifer model.

AquiferState carries model key, cumulative We, aquifer pressure where applicable,
previous reservoir pressure, elapsed seconds, initial aquifer pressure, immutable
parameter pairs and an immutable model-variable tuple for later extensions.
AquiferStepResult carries delta We, We, new state, interval-average rate,
endpoint rate where applicable, prior aquifer pressure, average boundary pressure,
driving pressure difference, and engineering warnings. Configurations and states
are frozen dataclasses. State never resides in loose Streamlit variables.

The current phase assumes **initial aquifer/reservoir pressure equilibrium**.
Initial aquifer pressure is the reservoir Pi, shown in the UI and retained in
the state. Different initial pressures or datum offsets are not modeled.

The timestep builds an immutable AquiferContext from the previous committed
state. Each candidate is evaluated independently, and the local solver cache
stores its balance and proposed aquifer step. Repeated/reordered trial pressures
produce identical values. A PressureSolution exposes a proposed aquifer step
only when pressure convergence and both closure tolerances pass. The timestep
then commits it with the reservoir state. Failed steps retain the previous
successful states and stop the run. Changing models or parameters requires a
fresh simulation and incompatible state is rejected.

A zero-volume history row exactly at the initial date is preserved for backward
compatibility; it represents the initial reference, not a zero-length aquifer
step. Its aquifer state is unchanged. `compute_step` itself rejects dt≤0.

## Equations, assumptions and references

| Model | Governing formulation | Reference and assumptions |
|---|---|---|
| NoAquifer | delta We=0; We=0 | Explicit no-support boundary condition; no empirical correlation |
| Pot | `Caq=Wi ct`; `We_n=Caq(Pi−p_n)`; `delta We=We_n−We_(n−1)` | Equilibrium compressible storage; pressure throughout the connected aquifer equals tank pressure; no flow resistance or transient diffusion. Pot storage is illustrated in Stanko's NTNU material-balance treatment [1] and the USGS report [2]. |
| Schilthuis | `q(t)=J(Pi−p(t))`; `pavg=(p_(n−1)+p_n)/2`; `delta We=J(Pi−pavg)dt`; `We_n=We_(n−1)+delta We` | R. J. Schilthuis (1936), *Active Oil and Reservoir Energy* [3]. Steady Darcy inflow, constant J and externally maintained reference pressure; no finite-capacity depletion. Trapezoidal integration is exact for an assumed linear boundary-pressure segment. |
| Fetkovich | `Caq=Wi ct`; `pa_(n−1)=Pi−We_(n−1)/Caq`; `pavg=(p_(n−1)+p_n)/2`; `delta We=Caq(pa_(n−1)−pavg)[1−exp(−J dt/Caq)]`; `We_n=We_(n−1)+delta We`; `pa_n=Pi−We_n/Caq` | M. J. Fetkovich (1971), *A Simplified Approach to Water Influx Calculations—Finite Aquifer Systems*, SPE 2603, Eqs. 6–9 [4]. Finite compressible water storage and constant aquifer PI; interval-average boundary pressure held constant during each exponential update. |

Fetkovich's published notation has `Wei=Wi ct Pi`, so `Wei/Pi=Caq` and
`(qwi)max/Wei=J/Caq`. The implementation evaluates the exponential fraction
with `−expm1(−J dt/Caq)` to preserve small-timestep precision. It is not a
constant-pressure-source Schilthuis approximation: aquifer pressure depends on
all previous cumulative influx. This interval method is not an exact solution
for a continuously varying boundary pressure; shorten history intervals around
rapid changes. A large J and long step equilibrate with the **interval average**,
not necessarily the endpoint reservoir pressure.

The reported average rate is `delta We/dt`. Schilthuis endpoint rate is
`J(Pi−p_n)`; Fetkovich endpoint rate is `J(pa_n−p_n)`. They can differ from the
interval-average rate, especially during pressure recovery. Pot has no
instantaneous Darcy-rate model; only the finite-interval average is reported.

All capacities, J, water volumes and total compressibilities must be finite
and positive. For finite models, Wi is the **connected initial reservoir water
volume**, including any geometry/contact fraction already applied. ct includes
aquifer water and pore-volume compressibilities. Constant parameters and a
single lumped aquifer pressure are assumptions, not calibrated field properties.
The reservoir's cf/cw terms remain separate; aquifer storage must not be counted
again as resident reservoir pore volume.

### Reference access and authority

1. M. Stanko, NTNU, TPG4230 (2021), *Problem 2-v2*, analytical material balance
   with pot aquifer, p. 5, Eqs. 1–3. The aquifer contribution is proportional to
   connected pore volume and water-plus-pore compressibility.
   [Author's teaching reference](https://www.ipt.ntnu.no/~stanko/files/Courses/TPG4230/2021/Class_files/20210204/Problem_2_v2.pdf).
   Pot is a storage approximation derived from compressibility, not a newly
   invented empirical correlation or a claim of unique original authorship.
2. USGS Water-Resources Investigations 76-61 (1976), aquifer influx
   representations including pot and Schilthuis.
   [Government engineering report](https://pubs.usgs.gov/wri/1976/0061/report.pdf).
3. R. J. Schilthuis (1936), *Active Oil and Reservoir Energy*, Transactions AIME
   118, pp. 33–52, DOI 10.2118/936033-G.
   [Publication](https://doi.org/10.2118/936033-G),
   [AIME publisher record](https://www.onemine.org/documents/papers-estimation-of-petroleum-reserves-active-oil-and-reservoir-energy-with-discussion-).
4. M. J. Fetkovich (1971), *A Simplified Approach to Water Influx Calculations—Finite
   Aquifer Systems*, JPT 23(7), pp. 814–828, DOI 10.2118/2603-PA.
   [Original publication](https://doi.org/10.2118/2603-PA),
   [Original-paper copy at NTNU](https://www.ipt.ntnu.no/~curtis/courses/Reservoir-Recovery/2018-TPG4150/Handouts/Material-Balance/Optional/1971-SPE2603-Fetkovich-Simplified-Water-Influx-MB.pdf).

Publisher metadata and publicly indexed original-paper equation text were
consulted. Some direct PDF/publisher fetches were unavailable; no claim is made
that all original full texts were retrieved. The equations were checked through
the primary engineering references above and independent hand cases below.
No open-source library is used as the engineering authority.

## Engineering units

| Parameter/result | Canonical engine | SI display | FIELD display |
|---|---|---|---|
| Pressure | Pa absolute | MPa | psia |
| Wi, We, delta We | reservoir m³ | reservoir m³ | rb |
| Caq | reservoir m³/Pa | reservoir m³/MPa | rb/psi |
| J | reservoir m³/(s·Pa) | reservoir m³/(day·MPa) | rb/(day·psi) |
| ct | 1/Pa | 1/MPa | 1/psi |
| Influx rate | reservoir m³/s | reservoir m³/day | rb/day |
| Time | seconds | dates / days | dates / days |

Conversions use 86400 s/day, 6894.757293168 Pa/psi and .158987294928 m³/rb.
UI unit changes convert an input draft and saved result views; they never modify
the engine's committed aquifer state. File-PVT conversion rules are unchanged.

## Flow reversal, validity and QC

The model uses an explicit **signed reversible-flow policy**. Pressure recovery
can give negative incremental We or cumulative We. Values are retained, with
warnings for efflux, decreasing We or net transfer into the aquifer. A negative
endpoint rate with a nonnegative interval average also raises a warning. There
is no hidden `max(0, ...)` clipping. This policy represents bidirectional hydraulic
communication; a one-way/check-valve aquifer is not implemented.

Nonfinite quantities, nonpositive parameters/dt, inconsistent finite-aquifer
pressure versus We, mismatched model/state, and aquifer pressure below 1 Pa
absolute are rejected. Aquifer pressure <0 cannot become a state. PVT coverage
still bounds the pressure solve; invalid candidate evaluation produces controlled
failure rather than silently skipping or extrapolating the domain.

Accepted support with `|We| > 2 max(|Fnet|, normalization_floor)` raises a review
warning. This is a heuristic for cancellation/extreme support, not a physical
cap. Non-convergence is reported separately from committed history. Negative
influx remains visible even if balance closure passes. Finite-aquifer pressure
and the pressure difference driving each interval are available for inspection.

The inspector retains production withdrawal, water/gas injection support,
Fnet, Eo/mEg/Efw and their N-scaled support, total expansion, incremental We,
cumulative We, total right-hand support, signed/relative residual, rates and
aquifer pressures. No normalized Phase 4 drive index is calculated. The separate
We and production/net-withdrawal terms make future normalization explicit.

## Independent numerical benchmarks

All fixtures begin at Pi=30 MPa. Expected values are hand-derived constants,
not results obtained by calling the tested implementation.

**Pot:** Caq=2500 m³/MPa. At pressures 28, 25 and 27 MPa, cumulative We is
5000, 12500 and 7500 m³; increments are 5000, 7500 and −5000 m³. The final
negative increment is a deliberate pressure-recovery check.

**Schilthuis:** J=100 m³/(day·MPa), boundary pressures 30→28→26→27 MPa,
durations 10, 20 and 5 days. Average pressure deficits are 1, 3 and 3.5 MPa.
Increments are `100×1×10=1000`, `100×3×20=6000`, and
`100×3.5×5=1750 m³`; cumulative We is 1000, 7000 and 8750 m³.
Endpoint rates are 200, 400 and 300 m³/day.

**Fetkovich:** Wi=1,000,000 m³, ct=.001 /MPa, Caq=1000 m³/MPa,
J=1000 ln(2) m³/(day·MPa), one-day steps. Thus the exponential fraction is
exactly 1/2 each day. Boundary pressure goes from 30 to 26 MPa in the first
interval, then stays at 26 MPa. First pavg=28 MPa; later pavg=26 MPa.

| Day | Delta We m³ | Cumulative We m³ | Aquifer pressure MPa |
|---:|---:|---:|---:|
| 1 | 1000 | 1000 | 29 |
| 2 | 1500 | 2500 | 27.5 |
| 3 | 750 | 3250 | 26.75 |
| 4 | 375 | 3625 | 26.375 |

Support halves after the initial boundary transition. For example, day 3 influx
is `1000(27.5−26)/2=750 m³`. The initial two steps intentionally have different
interval-average boundary pressures; support need not decrease across that
transition. Tests require errors ≤1e−10 m³ and ≤1e−8 Pa for this case.

**NoAquifer:** a complete six-step example result was captured before changing
the solver. Every pre-existing saved field is compared exactly, including
pressures, residuals, PVT, withdrawal/expansion, solver iterations/calls and
warnings. The new zero-aquifer fields do not alter the baseline.

Run `python examples/run_aquifer_acceptance.py` to regenerate all benchmark
and pressure-support outputs in `docs/phase_3a_numerical_results.json`.
See [the acceptance report](phase_3a_acceptance.md) for executed results.

## Limits and deferred work

The boundary pressure is the lumped reservoir pressure at a shared datum.
No spatial pressure solution, relative-permeability water invasion, salinity,
changing aquifer geometry, nonlinear compressibility or geometry-derived J is
modeled. Early-time transient diffusion is not represented by these models.
This phase adds no Carter–Tracy, VEH, modified VEH, aquifer matching, multi-tank
communication, optimization, forecasting, DCA or uncertainty engine.
