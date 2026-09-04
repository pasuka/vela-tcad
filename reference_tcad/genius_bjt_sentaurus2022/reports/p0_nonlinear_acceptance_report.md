# Genius NPN BJT P0 nonlinear acceptance validation

## Outcome

P0 passes. The M1 accepted-state workflow now requires all of the following:

- the Newton solve converges and all 5611 primary-state rows are finite;
- every physically qualified electron and hole continuity row satisfies the
  normalized local residual threshold `eps_row=1e-3`;
- global electron and hole contact-flux/source closure satisfies `1e-6`;
- three-terminal KCL is no larger than `1e-10 A/um`;
- a post-convergence, uncapped Newton probe has no active-carrier update larger
  than `1e-6 V`, using `delta(phi_n-psi)` and `delta(phi_p-psi)` rather than the
  quasi-Fermi updates alone.

The complete VBE=0.70 V, VCE=0.00..3.00 V chain passes 31/31 points. Generated
evidence is under `build-release/reference_tcad/genius_bjt_sentaurus2022/
m1_p0_acceptance/`; its schema-2 `manifest.json` SHA-256 is
`ed7cd0f39fd346f06e2e5e449efef083311df06028a1b08f95aa6b5d6d01b641`.

## Physical row qualification

The previous local-row check qualified rows only from their recombination or
generation source. P0 adds three independent physical supports:

1. carrier density at or above `1e16 m^-3`;
2. local absolute SG flux sum at or above `1e-6` of the carrier-wide maximum;
3. local recombination/generation source at or above `1e-6` of the
   carrier-wide maximum, while retaining the existing local source-to-row
   ratio and absolute floor.

Dirichlet contact rows and gauge rows without a continuity control volume are
explicitly marked inactive. They remain covered by the ordinary nonlinear
residual, but are not misinterpreted as local flux-balance equations.

This distinction was calibrated at VCE=0.8 V. Before the global source
qualification was added, 144 low-hole-density rows were selected solely by
sources of about `1e-18` in normalized units and blocked the 0.7-to-0.8 V
continuation. With the complete qualification, 7984 rows were enforced, 3238
rows remained diagnostic-only, the maximum qualified ratio was
`1.6675e-7`, and the same point converged normally.

## Full-chain maxima

| Quantity | Maximum | Bias |
|---|---:|---:|
| Qualified carrier-row residual ratio | 1.8059893e-5 | 1.6 V |
| Global electron continuity-closure ratio | 8.1714468e-8 | 2.2 V |
| Global hole continuity-closure ratio | 3.5596562e-8 | 1.6 V |
| Absolute three-terminal KCL | 1.7660790e-15 A/um | 2.2 V |
| Post-accept active-carrier raw update | 3.5195528e-10 V | 2.0 V |

All values pass their registered thresholds. At 3 V, the unconstrained
full-domain raw update is `0.02574996 V`, but the physically relevant update is
only `2.4342e-12 V`. This confirms why a global quasi-Fermi maximum must remain
a diagnostic rather than a hard gate in negligible-density regions.

The most extreme full-domain diagnostic occurs at 2.9 V: its uncapped raw
quasi-Fermi update is `1.3549474e5 V`, while the active-carrier update is only
`3.0215e-13 V`. This is retained in the manifest as evidence of severe
ill-conditioning in a negligible-density tail. P0 does not claim to remove
that tail; it prevents it from either hiding a physically supported continuity
error or falsely rejecting an otherwise closed physical state.

## Verification

- `test_newton_solver`: 100 test cases, 1310 assertions passed.
- carrier-row and global-closure focused tests passed.
- `newton_solve_from_state` and `example_runner_newton`: 3/3 CTest cases passed.
- Genius BJT reference/P0 Python regressions: 23 tests passed.
- The production `m1_spatial_vce3.json` strict configuration completed its
  one-point run successfully.
- Full repository CTest regression: 707/707 tests passed.

Generated states, VTK files, probes, and logs are intentionally not committed.
The manifest records the runner, mesh, doping, material, state, parent-state,
and probe SHA-256 values.
