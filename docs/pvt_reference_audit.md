# PVT correlation reference audit — Phase 2.1

Audited 11 September 2026. This table covers all 14 property-registry entries,
the independent Sutton pseudo-critical method, and the viscosity branches.
The [equation specification](pvt_phase_2.md#exact-equations-and-variants) records
every coefficient and canonical-unit conversion. Source attribution does not
constitute validation against independent laboratory data.

## Notation and unit boundary

In the field equations, `p` and `Pb` are **psia**, `t` is **°F**, `R` and `Rsb`
are **scf/STB**, `g` is gas specific gravity relative to air, and `o` is stock-tank
oil specific gravity. `API` is degrees API. `log` below means log10.
The model stores pressure in Pa, temperature in K, Rs in standard m³/stock-tank
m³, FVF in the corresponding m³/m³, viscosity in Pa·s and co in 1/Pa.
Absolute pressure is used throughout. The SI co equation is explicitly identified.

## Correlation reference table

| Method / property | Original author and publication | Equation implemented | Primitive input → output units | Applicability and implementation assumptions |
|---|---|---|---|---|
| Standing / Pb | M. B. Standing, 1947 [1] | `18.2[(Rsb/g)^.83 10^(.00091t−.0125API)−1.4]` | Field notation → psia | Standing population below. Uses the 18.2/1.4 analytical variant. |
| Standing / Rs | Standing, 1947 [1] | `g[(p/18.2+1.4)10^(.0125API−.00091t)]^(1/.83)` | Field notation → scf/STB | Inverse of selected Pb equation. Saturated primitive; the model separately anchors and caps Rs at Rsb. |
| Standing / Bo | Standing, 1947 [1]; engineering worked example [9] | `.972+.000147[R√(g/o)+1.25t]^1.175` | Field notation → rb/STB | Specific analytical variant; not the alternative .9759/.00012/1.2 refit. Published worked example independently tested. |
| Vasquez–Beggs / Rs | M. Vasquez and H. D. Beggs, 1980 [2] | `a g p^b exp[c API/(t+460)]`; two API coefficient regimes in equation specification | Field notation → scf/STB | VB population below. Gas gravity is on a 100-psig separator basis; no separator correction is implemented, so basis verification remains CAUTION. |
| Vasquez–Beggs / Pb | Vasquez and Beggs, 1980 [2] | `[Rsb/(a g exp[c API/(t+460)])]^(1/b)` | Field notation → psia | Exact inverse of implemented Rs, rather than separately rounded inverse coefficients. |
| Vasquez–Beggs / Bo | Vasquez and Beggs, 1980 [2] | `1+dR+(t−60)(API/g)(e+fR)`; API≤30 and API>30 coefficients | Field notation → rb/STB | Same separator basis; coefficients cross-checked with Al-Marhoun's engineering course [10]. |
| Glasø / Pb | Ø. Glasø, 1980 [3] | `10^(1.7669+1.7447x−.30218x²)`, `x=log[(Rsb/g)^.816 t^.172/API^.989]` | Field notation → psia | Glasø population below; positive temperature, valid quadratic branch required. |
| Glasø / Rs | Glasø, 1980 [3] | Lower quadratic root of Pb equation, then `g[10^x API^.989/t^.172]^(1/.816)` | Field notation → scf/STB | Saturated inverse, monotonic root; negative discriminant rejected. |
| Glasø / Bo | Glasø, 1980 [3] | `1+10^(−6.58511+2.91329u−.27683u²)`, `u=log[R(g/o)^.526+.968t]` | Field notation → rb/STB | Coefficient **.968**, not transposed .986; equation cross-check [11]. |
| Beggs–Robinson / dead and saturated oil viscosity | H. D. Beggs and J. R. Robinson, 1975 [4] | `x=t^(−1.163)exp(13.108−6.591/o)`; `μd=10^x−1`; `μs=10.715(R+100)^−.515 μd^[5.44(R+150)^−.338]` | Field notation → cP, converted ×.001 to Pa·s | Uses documented SG reformulation [12], whose rounded coefficients need not be numerically identical to an API-form variant. Guidance: 70–295°F, R 20–2070 scf/STB, o .75–.96, p≤5250 psia. |
| Vasquez–Beggs / undersaturated viscosity branch | Vasquez and Beggs, 1980 [2] | `μ=μ(Pb)(p/Pb)^m`; `m=2.6p^1.187 exp(−11.513−8.98e−5p)` | p/Pb in psia; μ anchor in cP → cP | Used only above Pb, continuous at Pb. Family applicability guidance, not separate validation of every branch. |
| Vasquez–Beggs / undersaturated co | Vasquez and Beggs, 1980 [2]; SI reformulation Afanasyev et al., 2004 [13] | `[28.1Rsb+30.6T−1180g+1784/o−10910]/(1e5 p)` | **Rsb m³/m³, T K, p Pa** → 1/Pa | Positive numerator required. Only at/above Pb. The raw Bo branch integrates co=K/p exactly; matched Bo co is its actual logarithmic derivative. |
| Sutton / pseudo-critical pressure and temperature | R. P. Sutton, 1985 [5], Eqs. 14–15 | `Ppc=756.8−131.0g−3.6g²`; `Tpc=169.2+349.5g−74g²` | Gas SG → psia, °R; wrapper converts to Pa, K | Population g=.571–1.679 (screen rounded .57–1.68). Sweet-gas use; no acid-gas correction. Independent method registry. |
| Dranchuk–Abou-Kassem / z | P. M. Dranchuk and J. H. Abou-Kassem, 1975 [6]; reproduced in [5], Eqs. 12–13 | Eleven-coefficient reduced-density EOS; `ρ=.27Pr/(zTr)`; full equation in specification | Dimensionless Pr, Tr → dimensionless z | Published principal range Tr=1–3, Pr=.2–30. Implementation requires Tr>1, Pr>0; out-of-range use receives caution. No subcritical root selection. Brent density solve, bracket [1e−12,20], tolerance 1e−13. |
| Real-gas equation / Bg | Thermodynamic identity; standard-condition convention [14] | `z T Psc/(p Tsc Zsc)` | T K, p Pa → reservoir m³/standard m³ | Psc=101325 Pa, Tsc=60°F, Zsc=1. Inherits selected z validity; gas volumes and Rs must share these standard conditions. Not an empirical named correlation. |
| McCain / Bw | W. D. McCain Jr., 1990 textbook and 1991 review [7] | `(1+ΔVp)(1+ΔVT)`; pressure/temperature polynomials in specification and [15] | p psia, t °F → rb/STB | Guidance p≤5000 psia, t≤260°F. The 32°F lower review boundary is implementation policy. No explicit salinity or water-SG dependence; those inputs are retained and the limitation is visible. |

## Applicability populations

These are reported development-population ranges, collated by comparative
engineering research [8]. They are screening guidance, not an independent
accuracy envelope for each property or the anchor-normalized model.

| Family | API | Gas SG | Temperature °F | Rsb scf/STB | Pb psia |
|---|---:|---:|---:|---:|---:|
| Standing | 16.5–63.8 | .59–.95 | 100–258 | 20–1425 | 130–7000 |
| Vasquez–Beggs | 15.3–59.5 | .511–1.351 | 70–295 | 0–2199 | 15–6055 |
| Glasø | 22.3–48.1 | .65–1.273 | 80–280 | 90–2637 | 165–7142 |

Sutton's original data span 200–12500 psia and 100–360°F, with gas gravity
.571–1.679. These describe that paper's samples, not a substitute for the
separate reduced-state validity screen of the selected z equation.

## Audit findings and corrections

1. **Sutton Ppc coefficient corrected from 131.07 to 131.0.** Original SPE 14265,
   printed page 3, Eq. 14, states 131.0. At gas SG=.7, Ppc changes from 663.287
   to **663.336 psia**. The existing test encoded the secondary implementation's
   coefficient and was corrected with this evidence. The +349.5g Tpc term was
   already correct. DAK coefficients and equation did not change.
2. **Undersaturated co applicability and matched derivative corrected.** It is
   unavailable below Pb. For `Bo_matched=A+B Bo_raw`, report
   `co_matched=co_raw (Bo_matched−A)/Bo_matched` at/above Pb. This is obtained by
   differentiating the actual matched Bo curve; it is not a new empirical co
   correlation. Bo values and Phase 1 MBE equations remain unchanged.
3. Sutton and DAK now have independent selectors/registries. Historical z key
   `sutton_dak_z` is retained for compatibility, but has no combined-method
   meaning in the UI or computation.
4. No FIELD/SI conversion equation was changed. Pressure is absolute;
   temperature affine conversions and viscosity/Bg/co scales were rechecked
   through reference and presentation tests.

The review distinguishes original-publication identification from full-text
inspection. Sutton's original paper was inspected and also supplies the DAK
equation and coefficients. Some original oil articles are behind publisher
access restrictions; their bibliographic identities are recorded below and
implemented forms are cross-checked against named engineering teaching,
research or calculator references. This is **not** a claim that every original
paper's full text was independently re-derived. Open-source rNodal and
GasCompressibility-py were secondary transcription checks only. Their results
do not supersede the original paper when a discrepancy is found.

## Publications and engineering cross-checks

1. Standing, M. B. (1947), *A Pressure-Volume-Temperature Correlation for Mixtures
   of California Oils and Gases*, Drilling and Production Practice, API,
   pp. 275–287. [Publication identifier](https://doi.org/10.2118/947275-G).
2. Vasquez, M.; Beggs, H. D. (1980), *Correlations for Fluid Physical Property
   Prediction*, JPT 32(6), 968–970. [Original publication](https://doi.org/10.2118/6719-PA).
3. Glasø, Ø. (1980), *Generalized Pressure-Volume-Temperature Correlations*,
   JPT 32(5), 785–795. [Original publication](https://doi.org/10.2118/8016-PA).
4. Beggs, H. D.; Robinson, J. R. (1975), *Estimating the Viscosity of Crude Oil
   Systems*, JPT 27(9), 1140–1141. [Original publication](https://doi.org/10.2118/5434-PA).
5. Sutton, R. P. (1985), *Compressibility Factors for High-Molecular-Weight
   Reservoir Gases*, SPE 14265, ATCE Las Vegas, 22–25 September.
   [Original-paper scan](https://yacsgas.wordpress.com/wp-content/uploads/2013/03/00014265.pdf),
   [DOI](https://doi.org/10.2118/14265-MS). Equations 12–15, printed page 3.
6. Dranchuk, P. M.; Abou-Kassem, J. H. (1975), *Calculation of Z Factors for
   Natural Gases Using Equations of State*, JCPT 14(3), 34–36.
   [Original publication](https://doi.org/10.2118/75-03-03).
7. McCain, W. D. Jr. (1990), *The Properties of Petroleum Fluids*, 2nd ed.,
   PennWell; and (1991), *Reservoir-Fluid Property Correlations—State of the
   Art*, SPERE 6(2), 266–272. [Review](https://doi.org/10.2118/18571-PA),
   [author-paper copy hosted by Texas A&M](https://blasingame.engr.tamu.edu/z_zCourse_Archive/P324_03A/Lecture_Refs_%28pdf%29/P324_Mod1_02_McCain_%28SPE_18571%29.pdf).
8. [Comparative engineering research, Table 5](https://iasj.rdd.edu.iq/journals/uploads/2025/02/05/596b02c184b5c5e4917164ac3e2bcfd3.pdf).
9. T. A. Blasingame, *Properties of Reservoir Fluids*, appendix, Eqs. A-37/A-38,
   worked solution A-50. [Texas A&M teaching reference](https://blasingame.engr.tamu.edu/z_zCourse_Archive/P613_23A/P613_Reference/%5B1988%5D_PVT_Appendix_%28Blasingame%29_%28pdf%29.pdf).
   Rs=569.54 scf/STB, t=220°F, o=.8217, g=.786 gives Bo=1.3687 rb/STB.
10. M. Al-Marhoun, *Black Oils Correlations Comparative Study*, slides 24 and 31.
    [Author's engineering course](https://www.restec.com/wp-content/uploads/2024/01/PETE-205-Black_oils_Correlations_Comparative_Study.pdf).
11. *New Computer Program to Calculate Pressure Drop in Pipes*, Sudan University
    of Science and Technology, chapter 3, Eqs. 3.63–3.65.
    [University repository](https://repository.sustech.edu/jspui/bitstream/123456789/15250/1/New%20Computer%20Program%20to%20Calculate%20Pressure%20Drop%20in%20Pipes.pdf).
12. [Pengtools engineering documentation: Beggs–Robinson](https://wiki.pengtools.com/index.php?title=Beggs_and_Robinson_correlation),
    SG formulation and worked example. This is a secondary engineering reference.
13. Afanasyev, V.; Moskvin, I.; Wolcott, K.; McCain, W. D. (2004),
    *Practical PVT Calculations for Black Oils*, YUKOS; SI formulation reproduced
    in [Pengtools compressibility documentation](https://wiki.pengtools.com/index.php?title=Vasquez_and_Beggs_Oil_Compressibility_correlation).
14. [Whitson black-oil documentation](https://manual.whitson.com/modules/bot/),
    standard-condition convention.
15. [Pengtools water FVF documentation](https://wiki.pengtools.com/index.php?title=Water_formation_volume_factor),
    McCain polynomial, example and engineering range.
