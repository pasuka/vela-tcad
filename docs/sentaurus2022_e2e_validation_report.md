# Sentaurus 2022 workflow validation report

> Historical snapshot: PN2D BV status was superseded by the M5 acceptance record
> (`docs/validation/sentaurus_gmsh_m5_acceptance_2026-08-26.md`). New main sync and
> reference-script coverage are tracked in the 2026-09-05 reference SDE plan.


Validation date: 2026-08-22

Machine-readable closeout: `docs/validation/sentaurus2022_remaining_tasks_20260822.json`.

## Build and regression status

- Windows MSYS2 UCRT64 Release targets were rebuilt and used for all new
  full-mesh simulations.
- Python regression coverage includes frontend provenance, ABA/BAB geometry,
  contact edge splitting, ordered doping replacement, directed Gaussian
  profiles, refinement controls, execution-stage ordering, IV/BV mode mapping,
  contract loading, and machine-readable failure reports.

## Contract executions

| Contract case | Runner | Result | Evidence |
| --- | --- | --- | --- |
| `pn2d_sentaurus2018_coarse7x3_iv` | Debug (earlier baseline) | accepted | 21 nodes, 24 cells, 16 points, last bias 0.3 V, maximum terminal-current imbalance 0.001216 |
| `pn2d_sentaurus2018_coarse7x3_bv` | Debug (earlier baseline) | accepted | reached the contracted -0.2 V stop with `bv_reverse` output mode |
| `pn2d_sentaurus2018_iv` | Release, explicit qualified profile | accepted | 1767 nodes, 3360 cells, 16 accepted points, last bias 0.3 V; maximum relative terminal imbalance 0.003882 and maximum mixed-tolerance ratio 0.01334 |
| `pn2d_sentaurus2018_bv` | Release (pre-strict baseline) | failed acceptance after reaching target | 1767 nodes, 3360 cells, 59 accepted solver points, last bias -20 V with `bv_reverse`; maximum absolute terminal mismatch 6.474e-4 A/um and mixed-tolerance ratio 20.0 |
| `pn2d_sentaurus2018_bv` | Release (strict convergence follow-up) | failed solve, correctly rejected | 7 accepted points, last stable bias -2.326689e-6 V; the next point failed global-continuity closure after adaptive reduction reached the minimum step |
| `pn2d_sentaurus2018_bv` | Release, explicit qualified profile | solved; failed acceptance | 1767 nodes, 3360 cells, 406 accepted solver points, last bias -20 V with `bv_reverse`; maximum relative terminal imbalance 0.11214 and mixed-tolerance ratio 1.3632 |

The contract now selects `pn2d_sentaurus2018_iv_qualified` or
`pn2d_sentaurus2018_bv_qualified` explicitly. The generic importer is no longer
the owner of this case-specific policy. Both profiles restore zero Gummel
iterations, no fallback, frozen mobility-field Jacobian terms for BV, and no
implicit quasi-Fermi/contact-majority policy. The selected profile is recorded
in both the deck and the run manifest.

The full PN2D IV low-current false rejection is closed. Acceptance retains the
raw relative diagnostic but applies a checked-in mixed absolute/relative gate.
The Release result also passes the frozen Sentaurus curve gate over 0.2--0.3 V:
four aligned points, matching trend, 0.001853 decades maximum magnitude delta
and 0.4276% maximum relative error.

## BV continuity follow-up

`global_continuity_closure=report` is now observationally pure: it no longer
changes the Newton line-search merit function. In `enforce` mode the closure
gate is evaluated only after the ordinary Newton convergence criterion. A C++
regression test verifies that `off` and `report` produce identical iteration
counts and current on the same problem. The fallback path also rejects an
unconverged Gummel state rather than exposing it as a converged curve point.

With the restored PN2D BV profile the Release solver again reaches -20 V, but
the result remains correctly `failed` at acceptance. The worst terminal-current
relative mismatch is 0.11214 (mixed-gate ratio 1.3632). At that point the
internal electron closure is about 1.87e-5 while the hole closure is 0.11361,
which matches the terminal imbalance. Tightening Newton `reltol` to 1e-10 does
not change the states, and changing the closure mode to `enforce` merely drives
the adaptive step toward its minimum. The remaining defect is therefore in the
assembled hole-continuity/contact-flux closure, not in curve post-processing or
the continuation settings.

Against the frozen Sentaurus BV curve, the high-field -19.9--5 V interval has
150 aligned points and a maximum delta of 0.8702 decades; the trend is still
mixed rather than monotone. At -20 V Vela reports approximately
-1.147e-16 A/um versus Sentaurus -9.105e-16 A/um. BV is solved but neither the
conservation gate nor the full Sentaurus-curve gate is accepted.

## Sentaurus 2022 official SDE oracle

The configured licensed VM ran Sentaurus Structure Editor T-2022.03-SP2. The
frozen `reference_tcad/pn2d_sentaurus2022_oracle` fixture contains a
variable-expanded `.log.eval`, official `_bnd.tdr`, `_msh.cmd`, `_msh.tdr`,
stdout/stderr, version output, expected Device IR and SHA-256 manifest. The
official mesh has 1943 points, 3680 elements and 3 regions. The VM evidence also
records T-2022.03's deprecation warning for the source deck's legacy
`sdegeo:define-2d-contact` calls.

The evaluated frontend now records the official `_msh.tdr` path and hash. The
device orchestration path imports and requalifies that TDR when present;
`offline-static` output is explicitly marked
`offline-static-approximation`. A frozen-oracle import check produced 1943
nodes, 3680 cells and 17 nodes on each terminal contact.
