# Phase 2 PVT engineering specification and validation

Phase 2 adds a property registry, sparse laboratory data, screening, constrained
matching, continuous black-oil branches, and PVT QC. The existing MBE equations,
domain records, chronological solver, PVT protocol and table interpolation are
unchanged. All fluid calculations live under `src/material_balance_studio/pvt/`.
The interface remains `properties_at_pressure(p)` with canonical SI values.

## Available models and workflow

Choose **PVT Source** in the existing PVT section:

1. **Tabulated PVT** uses the original complete measured table and linear
   interpolation. Its import contract and range restrictions are unchanged.
2. **Correlation PVT** uses independent selected property correlations and
   measured saturation anchors, with no regression transforms.
3. **Matched Correlation PVT** applies enabled, accepted property matches.
   At least one enabled fit is required; otherwise simulation is disabled.

The PVT subsections are Overview, Fluid Inputs, Laboratory Data, Correlation
Screening, Correlation Matching, Selected PVT Model and PVT QC. Recommendations
never change selections automatically. The registry supplies the selector
choices. Property plots show laboratory points, raw predictions, accepted
matched predictions and Pb markers. JSON exports preserve the chosen model,
canonical fluid, range, coefficients and fit records. Screening CSVs report
engineering-unit errors; screening JSONs retain every candidate's exact lab
pressures, measurements, predictions, metrics and applicability reasons.

## Inputs and units

`BlackOilFluid` stores API, gas gravity (air=1), temperature K, initial pressure
Pa absolute, optional measured Pb (Pa) and Rsb (standard gas m³/stock-tank oil m³),
optional oil SG, water SG and salinity mass fraction. Oil SG is calculated as
`141.5/(API+131.5)`; an optional measured SG is a consistency check with tolerance
0.005. Water SG and salinity are retained but are not explicit arguments in the
implemented McCain approximation. No separator or sour-gas composition correction
is implemented. Sour fluids receive NOT_APPLICABLE and cannot pass final QC.

At least one saturation anchor, Pb or Rsb, is needed to build the final model.
API, gas gravity, temperature and Pi alone cannot determine both. Pi is never
silently treated as Pb. Both anchors may be supplied; their inconsistency with
the selected untuned correlation is shown through the normalization-factor QC.

| Quantity | Canonical | SI file/API input | FIELD file/API input | UI SI |
|---|---|---|---|---|
| Pressure, Pb | Pa absolute | Pa absolute | psia | MPa |
| Temperature | K | °C | °F | °C |
| Rs/Rsb | m³/m³ | m³/m³ | scf/STB | m³/m³ |
| Bo, Bw | m³/m³ | m³/m³ | rb/STB | m³/m³ |
| Bg | reservoir m³/standard m³ | m³/m³ | rb/scf | m³/m³ |
| Oil viscosity | Pa·s | Pa·s | cP | Pa·s |
| Compressibility | 1/Pa | 1/Pa | 1/psi | 1/MPa |
| z, gravities, salinity | dimensionless | dimensionless | dimensionless | dimensionless |

`fluid_from_inputs` is the file/API conversion boundary; direct construction of
`BlackOilFluid` requires canonical units. Temperature conversions extend the
units package. Existing pressure, GOR, FVF and viscosity conversions are reused.
Project-unit changes preserve the canonical fluid, lab points, matches and saved
simulation. Uploaded file units are independent of project display units.

Sparse CSV/XLSX headers are `pressure` plus any of `bo,rs,oil_viscosity,z,bg,bw`.
Missing cells remain `None`; zero is allowed only for Rs. Pressure is required
for every populated row. Completely blank editor rows are ignored. Duplicate
pressures/headers, unknown columns, nonnumeric/nonfinite values and invalid signs
are rejected. Rows are explicitly sorted at this sparse-data boundary. The old
full-table importer still requires complete, already sorted rows.

The editor displays pressure in project units. Lab CSV downloads use import
units (Pa for SI), with exact lowercase headers and missing cells blank.
Suspicious Bo/Bw/z/Bg/Rs magnitudes generate review warnings. Units cannot be
inferred reliably from plausible numbers; explicit source-unit selection remains
mandatory. Salinity uses mass fraction, not ppm or percent.

## Registry and sources

Fourteen property entries and an independent pseudo-critical method registry are implemented. The viscosity entry contains separate
dead, saturated and undersaturated equations. The z entry evaluates the DAK
equation of state; pseudo-critical properties are selected independently.

| Property | Registry entries / formulations | References |
|---|---|---|
| Pb | Standing, Vasquez–Beggs, Glasø | [S1], [S2], [S3] |
| Saturated Rs | Standing, Vasquez–Beggs, Glasø | [S1], [S2], [S3] |
| Saturated Bo | Standing, Vasquez–Beggs, Glasø | [S1], [S2] |
| Oil viscosity | Beggs–Robinson dead/live; Vasquez–Beggs above Pb | [S4] |
| z | Dranchuk–Abou-Kassem; separate Sutton pseudo-criticals | [S5], [S6] |
| Bg | Real-gas law using selected z | [S7] and gas equation of state |
| Bw | McCain temperature/pressure approximation | [S8] |
| Undersaturated co | Vasquez–Beggs SI reformulation | [S9] |

Sources consulted, including publicly inspectable equation implementations:

- **[S1]** [rNodal author source: oil correlations](https://rdrr.io/github/f0nzie/rNodal/src/R/oil_correlations.R).
  Use the named Standing and Glasø equations; unrelated placeholders in that
  source are not used. Original publications: Standing (1947),
  DOI `10.2118/947275-G`; Glasø (1980), DOI `10.2118/8016-PA`.
- **[S2]** [Al-Marhoun, Black Oils Correlations Comparative Study, RESTEC](https://www.restec.com/wp-content/uploads/2024/01/PETE-205-Black_oils_Correlations_Comparative_Study.pdf),
  slides 24 and 31 for Vasquez–Beggs Rs and Bo. Original Vasquez–Beggs (1980),
  DOI `10.2118/6719-PA`.
- **[S3]** [Energies 14(9), 2653](https://www.mdpi.com/1996-1073/14/9/2653),
  Standing Pb formulation and ranges;
  [published comparative research, Table 5](https://iasj.rdd.edu.iq/journals/uploads/2025/02/05/596b02c184b5c5e4917164ac3e2bcfd3.pdf)
  for the Standing, Vasquez–Beggs and Glasø population ranges.
- **[S4]** [Pengtools: Beggs–Robinson viscosity](https://wiki.pengtools.com/index.php?title=Beggs_and_Robinson_correlation),
  SG formulation, example and applicability. Original Beggs–Robinson (1975),
  SPE-5434-PA; undersaturated branch from SPE-6719-PA.
- **[S5]** [GasCompressibility-py author documentation: DAK theory and example](https://aegis4048.github.io/GasCompressibility-py/theories.html).
  Original Dranchuk–Abou-Kassem (1975), DOI `10.2118/75-03-03`.
- **[S6]** [GasCompressibility-py Sutton source](https://aegis4048.github.io/GasCompressibility-py/_modules/sutton.html).
  The implementation uses the original paper's **positive** `349.5 SG` term in Tpc;
  the library theory page has a conflicting sign. Phase 2.1 also corrects its source-code Ppc coefficient 131.07 to the original Eq. 14 value 131.0. Original Sutton (1985), SPE-14265-MS.
- **[S7]** [Whitson black-oil table documentation](https://manual.whitson.com/modules/bot/)
  for the 60°F, 1-atm standard-condition convention. All gas totals and Rs/Bg
  must use consistent standard conditions; units alone do not reconcile them.
- **[S8]** [Pengtools: water formation volume factor](https://wiki.pengtools.com/index.php?title=Water_formation_volume_factor),
  McCain approximation, example and range; McCain, *Properties of Petroleum Fluids*.
- **[S9]** [Pengtools: Vasquez–Beggs compressibility](https://wiki.pengtools.com/index.php?title=Vasquez_and_Beggs_Oil_Compressibility_correlation),
  including the Afanasyev et al. (2004) SI reformulation and worked example.

## Exact equations and variants

In the field equations below, `p` and `Pb` are psia, `t` is °F, `R` is scf/STB,
`g` is gas SG, `o` is oil SG, and `API` is oil gravity. All logarithms marked
`log10` are base ten; `exp` is the natural exponential. Saturated Bo is rb/STB.
The registry converts the published field outputs to canonical units.

### Standing and Glasø

```text
Standing:
Pb = 18.2 [(Rsb/g)^0.83 × 10^(0.00091t − 0.0125API) − 1.4]
Rs = g [(p/18.2 + 1.4) × 10^(0.0125API − 0.00091t)]^(1/0.83)
Bo = 0.972 + 1.47e-4 [R sqrt(g/o) + 1.25t]^1.175

Glasø:
x = log10[(Rsb/g)^0.816 t^0.172 / API^0.989]
Pb = 10^(1.7669 + 1.7447x − 0.30218x²)
d = 1.7447² + 4(0.30218)(1.7669 − log10 p)
x = [1.7447 − sqrt(d)] / (2 × 0.30218)
Rs = g [10^x API^0.989 / t^0.172]^(1/0.816)
u = log10[R(g/o)^0.526 + 0.968t]
Bo = 1 + 10^(−6.58511 + 2.91329u − 0.27683u²)
```

Standing uses the 18.2/1.4 Pb/Rs variant and the 0.972/1.47e-4/1.175 Bo
variant in [S1], not the alternative 0.9759/0.00012/1.2 Bo refit. Glasø uses
the lower quadratic root on the monotonic branch; a negative discriminant or
nonpositive temperature is rejected. The coefficient is **0.968**, not 0.986.

### Vasquez–Beggs

```text
Rs = a g p^b exp[c API/(t+460)]
Pb = [Rsb / (a g exp[c API/(t+460)])]^(1/b)
Bo = 1 + d R + (t−60)(API/g)(e + f R)

API ≤ 30: a=.0362, b=1.0937, c=25.7240
          d=4.677e-4, e=1.751e-5, f=−1.811e-8
API > 30: a=.0178, b=1.1870, c=23.9310
          d=4.670e-4, e=1.100e-5, f=1.337e-9
```

Pb is the exact inverse of this Rs equation; rounded inverse coefficients are
not substituted. The published gas-gravity basis is a 100-psig separator.
Its verification is always a visible CAUTION because this version does not
convert arbitrary separator conditions. [S2]

### Oil viscosity and compressibility

```text
Beggs–Robinson SG variant, cP:
x = t^(−1.163) exp(13.108 − 6.591/o)
mu_dead = 10^x − 1
mu_sat = 10.715(R+100)^(−.515) mu_dead^[5.44(R+150)^(−.338)]

Above Pb, Vasquez–Beggs:
m(p) = 2.6 p^1.187 exp(−11.513 − 8.98e-5 p)
mu(p) = mu(Pb) (p/Pb)^m(p)

co SI reformulation, using Rsb in m³/m³, T in K, pressure in Pa:
co(p) = [28.1Rsb + 30.6T − 1180g + 1784/o − 10910] / (1e5 p)
```

The numerator of co must be positive. The canonical result is 1/Pa. For fixed
temperature/Rsb the model integrates `co=K/p` exactly, with `K=co(Pb)Pb`.
Thus above Pb, `Bo(p)=Bob(p/Pb)^(-K)`. Viscosity is converted from cP to Pa·s.
Sources: [S4], [S9].

### Gas and water

```text
Sutton: Ppc = 756.8 − 131.0g − 3.6g²  [psia]
        Tpc = 169.2 + 349.5g − 74g²    [°R]
Pr = p/Ppc; Tr = T_K × 1.8/Tpc; rho = 0.27 Pr/(z Tr)

DAK:
A = .3265 − 1.07/Tr − .5339/Tr³ + .01569/Tr⁴ − .05165/Tr⁵
B = .5475 − .7361/Tr + .1844/Tr²
C = −.1056(−.7361/Tr + .1844/Tr²)
z = 1 + A rho + B rho² + C rho⁵
    + [.6134(1+.7210rho²)rho²/Tr³] exp(−.7210rho²)

Bg = z T_K Psc/(p_Pa Tsc Zsc)
Psc=101325 Pa; Tsc=288.7055555556 K (60°F); Zsc=1

McCain Bw (p psia, t °F):
dVp = −1.95301e-9 p t −1.72834e-13 p²t −3.58922e-7 p −2.25341e-10 p²
dVT = −.010001 +1.33391e-4 t +5.50654e-7 t²
Bw = (1+dVp)(1+dVT)
```

DAK solves `rho z(rho) − .27Pr/Tr = 0` using Brent, density bracket
`[1e-12,20]`, tolerance `1e-13`. Only `Pr>0, Tr>1` and sweet gas are supported;
no two-phase root selection is attempted. Bg is reservoir/standard m³ and
inherits the selected z applicability. Sources: [S5]–[S8].

## Validity and saturation branches

Population ranges are guidance rather than guarantees for an individual fluid:

| Family | API | Gas SG | Temperature °F | Rsb scf/STB | Pb psia |
|---|---:|---:|---:|---:|---:|
| Standing | 16.5–63.8 | .59–.95 | 100–258 | 20–1425 | 130–7000 |
| Vasquez–Beggs | 15.3–59.5 | .511–1.351 | 70–295 | 0–2199 | 15–6055 |
| Glasø | 22.3–48.1 | .65–1.273 | 80–280 | 90–2637 | 165–7142 |

These family ranges [S3] are applied as conservative screening guidance to
Pb/Rs/Bo and the VB co branch; they are not claimed to be separately validated
accuracy limits for every derived property. Beggs–Robinson checks 70–295°F,
Rs 20–2070 scf/STB, oil SG .75–.96 and p up to 5250 psia [S4]. Sutton checks
SG .57–1.68; DAK checks Tr 1–3, Pr .2–30 [S5]. McCain checks up to 5000 psia
and 260°F [S8]; 32°F is an additional implementation review boundary.

**VALID** means these published checks pass; **CAUTION** means range/basis review
is required; **NOT_APPLICABLE** means a required input or mathematical/composition
domain is unsupported, or that property's sampled physical checks fail. Cautions
remain selectable. Final physical failures disable the simulation.

Measured Pb is used by default; the engineer can explicitly choose predicted
Pb. Measured Rsb is retained if supplied. With only Pb, Rsb is calculated from
the selected Rs correlation at Pb. With only Rsb, Pb is calculated from the
selected Pb correlation. The final untreated Rs curve is:

```text
p < Pb: Rs(p) = Rsb × Rs_published(p)/Rs_published(Pb)
p ≥ Pb: Rs(p) = Rsb
```

This anchor normalization is explicit even in **raw** mode: raw means no lab
regression, not disregard of measured saturation anchors. It permits independent
Pb/Rs family selection without a jump. The unnormalized published primitives
remain available through the registry. Pb screening always compares the
published predicted Pb against measured Pb and needs an independent measured
Rsb; it does not return the measured anchor for every candidate.

Bo and live viscosity below Pb use the selected Rs curve. Above Pb, Bob and
mu(Pb) anchor the compressibility/viscosity branches described above. Positive
slope transforms preserve continuity; Rs transforms must preserve Rsb.

## Ranking, matching and QC

For each measured property, all registered candidates are evaluated at exactly
its measured pressures. Results retain complete arrays and failure reasons.
Errors use `prediction − measurement`: RMSE, MAE, bias, MAPE and maximum absolute
percentage error. Zero measured values remain in RMSE/MAE/bias; percentage metrics
exclude them and expose `percentage_n`. Empty/zero-only percentage support yields
`None`, not infinity. Unmeasured properties receive no numerical recommendation.

Ranking is lexicographic: VALID, then CAUTION, then NOT_APPLICABLE; within each
tier, increasing RMSE in canonical units, then stable registry key. Errors from
different properties are never compared as one score. An exact but out-of-range
fit can rank below an in-range alternative. The UI supports sorting and manual
selection; selecting a correlation never fits it automatically.

Least squares treatments:

| Property | Matched formulation | Minimum information |
|---|---|---|
| Bo, Bw, oil viscosity, z | `A + B raw`, centered affine regression | ≥2 distinct predicted values |
| Rs | `Rsb + B(raw−Rsb)`; `A=Rsb(1−B)` | ≥1 point below Pb |
| Bg | `B raw`, A=0; applies the same scale to z | ≥1 point |
| Predicted Pb | `B raw`, A=0 | Measured Pb and independent Rsb |

All B values must be positive and coefficients finite. Bg and z cannot both
have independent transforms: this would double-count the same gas state.
Scaling Bg scales z so the real-gas identity remains exact. A fixed measured Pb
cannot be fitted again. Matching Rs to plateau-only data is unidentifiable and
rejected. Sparse properties can be fitted without measurements of the others.

Every attempt returns acceptance/rejection and reasons; accepted records contain
property, correlation, A, B, n, formulation, pre/post RMSE and pre/post bias
(plus other metrics). Models are immutable and the raw model remains available.
Regression coefficients are expressed against canonical raw values. The report
converts A and dimensional errors to engineering units; B is dimensionless.

The complete transformed model is checked over 1,001 uniform range samples, all
in-range lab pressures, Pb and Pb±1e-6 relative pressure. The active UI model is
checked again after enabling/disabling fits. Failures include nonfinite/negative
properties, zero Bo/Bg/Bw/viscosity/z, nonmonotonic expected branches, Rs exceeding
Rsb or changing above Pb, relative Pb jumps >1e-4, or lab pressures outside the
requested range. Expected trends: Rs and Bo increase below Pb; Bo decreases
above Pb; oil viscosity decreases below and increases above Pb; Bg and Bw
decrease with pressure. z is deliberately not constrained to be monotonic.

Cautions report sparse/missing lab data, lab extrapolation, correlation ranges,
water-composition limitations, heuristic magnitude envelopes and an Rs anchor
normalization differing from one by >10%. QC sampling is not proof of accuracy
between samples or outside laboratory coverage. A Pb at/outside project bounds
gets a caution because both sides cannot be checked within coverage.

Changing fluid, bounds, laboratory data, saturation treatment or property
selections clears prior UI fits. Changing an upstream match clears downstream
fits: Pb clears all, Rs clears Bo/viscosity, z and Bg exclude each other. Historical
fit metrics describe the regression at acceptance; current combined-model lab
errors are recalculated in QC, including after disabling an upstream fit.

## Reference calculations and tests

The following fixed reference values were calculated independently with 50-digit
decimal arithmetic from the cited equations. They are stored as constants in
tests and are not generated using the production implementation. They verify
equation transcription/units, not independent fluid accuracy. Common inputs:
API=35, gas SG=.75, T=180°F, p=2500 psia, R/Rsb=500 scf/STB,
oil SG=141.5/(35+131.5).

| Reference calculation | Expected field result |
|---|---:|
| Standing Pb | 2113.543928855771 psia |
| Vasquez–Beggs Pb | 2366.736834072329 psia |
| Glasø Pb | 2451.815344319128 psia |
| Standing Rs | 610.7613681303263 scf/STB |
| Vasquez–Beggs Rs | 533.5913496932905 scf/STB |
| Glasø Rs | 511.6207937670833 scf/STB |
| Standing Bo | 1.292949195140943 rb/STB |
| Vasquez–Beggs Bo | 1.2988436000 rb/STB |
| Glasø Bo | 1.258464495419306 rb/STB |
| Beggs–Robinson saturated viscosity | .6398767305448119 cP |
| McCain Bw | 1.028363968194184 rb/STB |

Additional external/reference cases:

- DAK published example [S5]: Pr=3.1995, Tr=1.5006 gives
  z=.7730934971021096. Sutton SG=.7 gives Ppc=663.336 psia, Tpc=377.59°R.
- Real-gas identity at standard conditions: z=1 gives Bg=1 m³/m³;
  ten times standard pressure and z=.8 gives Bg=.08 at the same temperature.
- Compressibility [S9]: Rs=53.24 m³/m³, T=363 K, oil SG=.85, gas SG=.75,
  p=10 MPa gives .002907667529411765 /MPa by independent arithmetic. The source
  reports .002907 /MPa; the rounded-source comparison uses absolute tolerance
  1e-6 /MPa, with a separate tight exact-arithmetic comparison.
- Viscosity [S4]: T=137°F, API=22, Rs=90 gives approximately 17.44 cP dead
  and 8.24 cP saturated in the source. It reports oil SG rounded to .922;
  the exact API-derived SG gives 17.46338 and 8.24722. Tests allow .03 and .02 cP
  respectively; the independent 180°F reference above uses tight tolerance.
- Undersaturated viscosity: mu(Pb)=1.2 cP, Pb=2500, p=4000 psia gives
  1.409550354932235 cP by independent arithmetic.

The test suite covers every registry entry, both VB API regimes, FIELD/SI
equivalence, branches/continuity, ranking, constrained matching and rejection,
sparse CSV/XLSX/editor workflows, MBE depletion/injection, and all three UI modes.
The interface regression injects identical continuous property functions through
the registry into the real `CorrelationPVTModel`; matching Bo/Rs/Bg/Bw yields
the same MBE pressures to 1e-6 Pa and balance terms to numerical tolerance as the
table model. No solver or equation stubs are used.

## Reproducible matching example

Use Correlation PVT with the synthetic default fluid (API 35, gas SG .75,
90°C, Pi 30 MPa, Pb 18 MPa, Rsb 90 m³/m³; range 2–35 MPa). Load the bundled
sparse lab example, select Bo and click **Match selected bo**. The lab file uses
`Bo_lab=.04+1.02 Bo_raw` at 4, 8, 12, 18, 24 and 30 MPa. Expected fit:

| Quantity | Before | After |
|---|---:|---:|
| RMSE (m³/m³) | .06447544321627238 | approximately 4.2e-16 |
| Bias (m³/m³) | −.06445772161037733 | approximately −4.1e-16 |
| MAPE | 5.018291% | numerical zero |

Coefficients: A=.04 m³/m³, B=1.02, n=6. This is a synthetic regression
recovery check, not validation against a real reservoir fluid. Switch to Matched
Correlation PVT to use the fit in the unchanged MBE engine.

Headless equivalent: `python examples/run_pvt_example.py`. It prints screening,
the match record, QC failures and raw/matched MBE convergence.

## File inventory and integration boundary

Created under `src/material_balance_studio/`:

```text
units/temperature.py
pvt/fluid.py
pvt/laboratory.py
pvt/correlations/{base,bubble_point,solution_gor,oil_fvf,oil_viscosity,
                  gas_z,gas_fvf,water,registry}.py
pvt/models/correlation_model.py
pvt/matching/{metrics,ranking,regression}.py
pvt/qc/{validity,physical_checks}.py
presentation/{pvt,pvt_workflow}.py
```

The new subpackages also contain `__init__.py`. Added examples:
`correlation_lab_SI.csv`, `run_pvt_example.py`. Added tests:
`test_correlation_references.py`, `test_correlation_workflow.py`, `test_pvt_app.py`.
Changed existing files: `app.py` and `README.md` only. Added this document and
`phase_2_validation.md` / `phase_2_core_hashes.json` as validation records.

No MBE/solver/domain/core-PVT file is changed. A raw or matched model supplies
Bo, Rs, Bg, Bw, oil viscosity and z. Gas viscosity is unimplemented (`None`).
Separate injection FVFs are not calculated; the existing resident-fluid FVF
fallback and its existing simulation warnings apply unchanged. Aqueous influx,
reservoir history matching, multiple tanks, forecasting, DCA and uncertainty
remain outside Phase 2.


Phase 2.1 supersedes the original reference hierarchy: see [the original-publication audit](pvt_reference_audit.md). The z registry key `sutton_dak_z` is retained solely for configuration compatibility; the displayed z method is Dranchuk–Abou-Kassem, and `pseudo_critical_method` independently selects Sutton. For matched Bo, reported undersaturated co is `co_raw × (Bo_matched−A)/Bo_matched`; below Pb this undersaturated property is unavailable. Dedicated boundary probes also include 0.99, 0.999, 1, 1.001 and 1.01 Pb. Match records now include maximum absolute full-range change, pressure at maximum change, percentage change where defined, and before/after QC status.
