# Phase 2.1 engineering acceptance

Completed 11 September 2026. **262 tests passed**: the existing 228 plus
34 dedicated acceptance cases. All 136 Phase 1/1.1 tests remain intact.
One existing Phase 2 reference expectation was corrected after comparison with
Sutton's original publication. No aquifers or new empirical correlation families
were added.

## Bubble-point acceptance

Dedicated tests cover Standing, Vasquez–Beggs and Glasø Rs/Bo combinations,
each raw and matched, at 0.99, 0.999, 1, 1.001 and 1.01 Pb. All six combinations
pass: Rs increases to Rsb and stays at Rsb above Pb; Bo peaks at Pb; oil viscosity
reaches its minimum there and increases above it. Changes shrink towards Pb.
Additional ±1e−6 Pb probes pass the relative-change threshold of 1e−4.

Representative raw Standing fixture: API 35, gas SG .75, 90°C, Pi=30 MPa,
Pb=18 MPa, Rsb=90 m³/m³; requested pressure range 2–35 MPa.

| p/Pb | Rs m³/m³ | Bo m³/m³ | Oil viscosity Pa·s | Undersaturated co 1/Pa |
|---:|---:|---:|---:|---:|
| .99 | 88.92723537 | 1.302105730 | .000586498358 | Not applicable |
| .999 | 89.89262523 | 1.304883149 | .000583148902 | Not applicable |
| 1 | 90 | 1.305192279 | .000582779479 | 2.191991303e−9 |
| 1.001 | 90 | 1.305140809 | .000582915800 | 2.189801502e−9 |
| 1.01 | 90 | 1.304679963 | .000584149650 | 2.170288419e−9 |

The co tests use the undersaturated one-sided derivative at Pb and centered
differences above Pb. Both raw and affine-matched Bo agree with `−d ln Bo/dp`
to relative tolerance 2e−6. Calling undersaturated co below Pb raises an explicit
not-applicable error.

## Full-range QC and matching

Each QC pass uses 1,001 uniformly spaced pressures plus in-range lab pressures
and the dedicated Pb probes. This example has **1,013 unique pressures**.
Raw and accepted matched models each produce **26 PASS, 10 CAUTION, 0 FAIL**.
All six plotted properties satisfy their positivity requirements; co is positive
where applicable. Oil branch trends, gas/water FVF trends, the Rsb plateau and
Pb continuity pass. z is not incorrectly constrained to be monotonic.

The ten cautions are retained: saturation-anchor normalization .813945, the
upper end of the water-correlation pressure range, unverified VB separator
gravity basis, five unmeasured properties, Bo extrapolation beyond lab coverage,
and water composition not explicitly modeled. This is an accepted synthetic
engineering fixture with visible cautions, not independently validated fluid data.

| Bo fit diagnostic | Accepted fixture | Rejected extrapolation fixture |
|---|---:|---:|
| A m³/m³ | .04 | −12 |
| B | 1.02 | 10 |
| Lab pressures MPa | 4, 8, 12, 18, 24, 30 | 16, 18, 20 |
| RMSE before m³/m³ | .0644754432163 | .393710912179 |
| RMSE after m³/m³ | 4.1541e−16 | 0 |
| Maximum full-range absolute change m³/m³ | .0661038455895 | 2.261919383337 |
| Pressure at maximum absolute change MPa | 18 | 2 |
| Maximum full-range relative change % | 5.696827067 | 209.048119967 |
| QC before → after | CAUTION → CAUTION | CAUTION → FAIL |
| Result | Accepted | Rejected: Bo(2 MPa)=−1.179910426 |

Maximum change means `max |candidate−raw|` over the project grid. It does not
estimate error against unknown measurements. Relative and absolute maxima can
occur at different pressures. Rejected candidates retain coefficients, metrics,
range-change diagnostics and rejection reasons but cannot become active models.
An additional accepted-fit test omits Pb from laboratory data and still finds
the largest absolute Bo change at Pb.

## Reference audit and corrections

The [reference audit](pvt_reference_audit.md) maps all 14 property entries, the
independent pseudo-critical method and viscosity branches to authors,
publications, exact equation variants, input/output units, applicability and
implementation assumptions. It distinguishes full-text evidence from engineering
cross-checks; open-source implementations are not the sole authority.

- Sutton original SPE 14265 Eq. 14 uses **131.0**, replacing the secondary
  implementation's 131.07. At gas SG=.7, Ppc is now 663.336 psia. The single
  incorrect reference expectation and its constructed reduced-state input were
  updated. DAK's equation and coefficients are unchanged.
- Undersaturated co is now unavailable below Pb and is derivative-consistent
  with affine-matched Bo: `co_raw (Bo_matched−A)/Bo_matched`.
- No FIELD/SI conversion equations changed.
- Sutton pseudo-critical selection and Dranchuk–Abou-Kassem z selection are
  independent. Tests demonstrate that either method can be replaced independently.

## MBE interface equivalence

The executable [comparison fixture](../examples/pvt_acceptance_case.py) supplies
identical continuous Bo/Rs/Bg/Bw functions through TablePVTModel and the actual
CorrelationPVTModel registry dispatch. Temporary fixture entries are removed
afterwards. This tests the abstraction and existing solver, not empirical
correlation accuracy or the physical validity of the synthetic gas identity.

N=1,000,000 m³, initial pressure=30 MPa, m=.15, Swc=.2, cf=5e−10 /Pa,
cw=4e−10 /Pa. Three cumulative history steps include oil/gas/water production
and water/gas injection. The fixture deliberately produces nonzero Eo, mEg and
Efw. All four property functions agree exactly over 1,001 comparison pressures;
both simulations converge.

| Quantity | Maximum absolute difference across three steps | Tolerance |
|---|---:|---:|
| Pressure Pa | 7.4506e−9 | 1e−6 |
| Fnet m³ | 0 | 1e−7 |
| Eo | 0 | 1e−12 |
| mEg | 7.9797e−17 | 1e−12 |
| Efw | 7.5894e−18 | 1e−12 |
| Balance residual m³ | 9.0950e−11 | 1e−7 |

Final table-model values: pressure=25.128402953819 MPa,
Fnet=68204.79444505674 m³, Eo=.019486388184723813,
mEg=.04384437341562875, Efw=.004874032844704065,
residual=1.1641532182693481e−10 m³. Both complete step-by-step sets are retained
in [numerical evidence](phase_2_1_numerical_results.json).

## Running-app visual QC

Inspected the running Streamlit app at `http://127.0.0.1:8501/` in the in-app
browser, 1280×720 viewport, dark theme. Loaded the bundled sparse Bo lab example,
performed its fit, activated Matched Correlation PVT and opened Selected PVT
Model. Visually inspected each of the six complete plots in **both SI and FIELD**.

| Plot | SI ordinate | FIELD ordinate | Observed behavior |
|---|---|---|---|
| Bo | m³/m³ | rb/STB | Six laboratory diamonds on matched curve, distinct raw curve, continuous peak at Pb |
| Rs | m³/m³ | scf/STB | Increasing saturated branch and constant Rsb plateau |
| Oil viscosity | Pa·s | cP | Continuous minimum at Pb, rising undersaturated branch |
| Gas z-factor | Dimensionless | Dimensionless | Smooth minimum; no artificial monotonic constraint |
| Bg | m³/m³ | rb/scf | Smoothly decreasing curve with consistent field conversion |
| Bw | m³/m³ | rb/STB | Smoothly decreasing curve |

Pressure axes display MPa or psia. Titles, legends, axes and Pb indicators are
readable. Raw and matched curves coincide for properties without transforms;
the bundled lab has Bo measurements only. The 12 automated plot cases also
verify lab/raw/matched traces for each property in both unit systems.

The initial dark-theme check exposed low-contrast lab markers and Pb lines.
Laboratory points now use amber diamonds with an outline, and Pb lines use
visible gray. Coincident raw/matched Pb is labeled explicitly. All twelve
rendered plots were reviewed after that correction.

The selected-model table displays separate Sutton and DAK rows, treatments and
per-property PASS/CAUTION/FAIL. Bubble-point details show measured 18 MPa,
untuned calculated 15.145224 MPa and active 18 MPa. The matching screen shows
the recovered A/B, before/after errors and full-range change; vertical property
details make the main fit values accessible without scrolling a wide table.

## Tests and protected core

The 34 new cases cover six raw/matched family boundary cases, co applicability
and derivative consistency, independent gas methods, dense-grid coverage,
accepted/rejected extrapolation, negative slopes, per-property QC propagation,
model summary and Pb metrics in both units, 12 plot combinations, a published
Standing Bo worked example, and nonzero-term MBE equivalence.

Validation commands, run from the repository root:

```powershell
.\.venv\Scripts\python.exe -m pytest -q -p no:cacheprovider
.\.venv\Scripts\python.exe examples/run_pvt_acceptance.py
```

The complete suite passed **262/262 in 28.10 seconds**. After the visual marker
contrast adjustment, the affected acceptance and app suites passed **43/43 in
21.05 seconds**. SHA-256 comparison against the captured Phase 2 baseline confirms
**all 13 MBE, solver, domain, base-PVT and table-PVT files remain byte-for-byte
unchanged**. Hashes and per-file outcomes are in the numerical evidence JSON.

Changed areas: pseudo-critical method registry; correlation input/registry and
model; validity and dense physical QC; matching diagnostics; PVT presentation
and workflow; app phase caption; reference tests; README and engineering docs.
Added `tests/test_pvt_acceptance.py`, `examples/pvt_acceptance_case.py`,
`examples/run_pvt_acceptance.py`, the reference audit and this acceptance record.
