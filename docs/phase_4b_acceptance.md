# Phase 4B acceptance

**461 tests passed in 56.87 seconds: all 439 existing tests plus 22 new tests.**
No existing tests or protected forward equations were changed.

Bounded history matching uses the existing forward simulator and SciPy trust-region
reflective least squares. See [method, registry and limitations](history_matching_phase_4b.md).

## Independent recovery benchmarks

Eight pressure surveys and nonuniform intervals are prescribed independently of
the forward simulator. For volumetric oil with constant Rs, Bo=Boi+c(Pi−P),
surface production is defined algebraically by Np=N*c*(Pi−P)/Bo. This supplies
known OOIP without generating observations from the optimizer or simulator.
The Fetkovich benchmark independently calculates accepted influx using
deltaWe=C*(Pa−(Pprevious+Pcurrent)/2)*(1−exp(−J*dt/C)), then sets
Np=(N*c*(Pi−P)+We)/Bo. Capacity and compressibility remain fixed.

| Case | True value | Fitted value |
|---|---:|---:|
| OOIP alone | 1,000,000 m³ | approximately 1,000,000 m³ |
| Fetkovich J alone | 1e-9 m³/(Pa·s) | 9.999999999994715e-10 |
| Joint OOIP | 1,000,000 m³ | 999999.9999999952 |
| Joint J | 1e-9 m³/(Pa·s) | 1.0000000000000129e-9 |

Tests require OOIP-only recovery within 2e-6 relative tolerance, J/joint recovery
within 2e-4, and final noiseless pressure RMSE below 1 Pa. Gas-cap selection and
recovery are tested independently; inactive OOIP remains fixed.

One deliberately noisy survey has a +2 MPa error and sigma=5 MPa; the other seven
have sigma=.05 MPa. Unweighted fitting gives N=1,145,470.029984 m³; weighted fitting
gives N=1,000,019.681747 m³. Both are checked against a separate scalar minimization
of the analytical pressure formula P=Pi−Boi*Np/[c*(N−Np)]. Weighted fitting reduces
the outlier's influence, rather than minimizing the unweighted RMSE.

A restrictive upper bound of 900,000 m³ returns 899999.9999999998 m³ and an upper
bound CAUTION. All recorded candidates remain inside the original physical bounds.
FIELD/SI equivalence covers transformed bounds, converted pressure/volume/sigma
history inputs, displayed results and edited UI bounds. Repeated identical runs
return identical fitted parameter tuples.

Failure tests verify fixed finite penalties and captured failure reasons. A bad
MBE closure candidate is rejected even with matching pressures. A wholly invalid
final candidate cannot be labeled converged. Every supported aquifer passes the
objective-interface check. Observation exclusion, missing-sigma policy, exact-date
alignment, data sufficiency, fixed-parameter enforcement and explicit UI application
are covered.

## Regression and reproducibility

The protected source manifest and entire Phase 4A numerical acceptance output
remain unchanged. The maximum relative MBE residual across valid matching
benchmark candidates is approximately **1.075e-13**, still below the unchanged
1e-8 forward tolerance. This larger candidate set differs from the prior fixed
acceptance cases, whose outputs remain identical.

[Numerical evidence and hashes](phase_4b_numerical_results.json) can be reproduced:

```
.\.venv\Scripts\python.exe -m pytest -q -p no:cacheprovider
.\.venv\Scripts\python.exe examples/run_phase4b_acceptance.py
```

Base reservoirs are not overwritten automatically. Multiple in-session results
remain independent scenarios, with explicit parameter application. No uncertainty
or automated identifiability analysis was implemented.
