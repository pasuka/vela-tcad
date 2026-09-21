# Templates/LDMOS Stage-4 Vg=4 V, Vd=10 mV Newton/Jacobian audit

Date: 2026-09-03

## Outcome

The first nonzero-drain D5 transfer is no longer blocked.  The audit first
found and fixed a configuration-contract defect: with Fermi-Dirac carrier statistics,
`mobility.jacobian_field_derivatives=false` still recomputed field-dependent
mobility inside the numerical transport Jacobian.  After making the switch
actually freeze mobility, the exact 0 -> 10 mV transfer converged in four
Newton iterations with no predictor and no IALMob.

This is a solver-linearization result, not a mobility calibration.  The HFS
mobility used by the residual, the AverageBox carrier couples, material data,
contact reconstruction and all physical parameters are unchanged.

## Frozen failed-state comparison

The original live-mobility Jacobian run moved away after its best third
iteration:

| State | Poisson block | electron continuity block |
| --- | ---: | ---: |
| iteration 0 | 66.5411 | 14497.6 |
| best iteration 3 | 2.20930 | 3229.08 |
| final iteration 40 | 18.1294 | 20034.7 |

The final electron residual is almost entirely transport flux.  Its leading
nodes are 687, 682, 683, 675, 702, 713 and 695 in the heavily doped drain
contact neighbourhood.  SRH/Auger terms are negligible there.  The configured
continuity row scaling is inactive at these rows because the source magnitude
is far below both the absolute source floor and the configured fraction of the
flux sum; the row weights are therefore one.

## Jacobian evidence

A 60-node, frozen one-ring potential JVP audit was run on iteration 0, best
iteration 3 and final iteration 40.  It confirms that the contact-cell
ElectricField fallback depends on triangle nodes beyond the transport edge
endpoints while the current analytic stencil omits some of those columns.  The
maximum electron-block JVP relative error grows from zero at the equilibrium
initial state to 7.224e-4 at the best state and 2.282e-3 at the final state.
This is a real follow-up defect, but its measured global magnitude is too small
to explain the failed transfer by itself.

The more important finding was that the intended frozen-mobility control did
not work: live and frozen configurations initially produced bit-identical
steps.  The Fermi-Dirac branch selected a numerical edge-flux derivative and
continued to recompute mobility regardless of
`jacobian_field_derivatives=false`.  The implementation now passes the current
edge mobility into that derivative when field derivatives are disabled.  A
dedicated Catch2 regression test verifies that residual/finite-difference truth
is unchanged while the analytic Jacobian responds to the switch.

With the repaired control, freezing mobility-field derivatives changes the
Newton direction decisively:

| Frozen state | raw-step ratio, frozen/live | trial electron residual ratio, frozen/live |
| --- | ---: | ---: |
| best iteration 3 | 0.20735 | 0.16545 |
| final iteration 40 | 8.510e5 | 0.25335 |

The enormous uncapped final-state step in both formulations identifies a
remaining near-null carrier mode.  Nevertheless, after the configured QF cap
and globalisation, the lagged-mobility direction reduces the electron residual
instead of amplifying it.

## Reclose gate

The exact D5 0 -> 10 mV sweep with lagged mobility derivatives passed:

- Vd=10 mV converged in 4 Newton iterations.
- Final Poisson block: 7.693e-10.
- Final electron continuity block: 9.271e-13.
- Drain current: 9.254406e-7 A/um.
- Source/drain/substrate terminal balance closes at floating-point level.

The Vg=4 V full curve has therefore been started under the same contract.  It
must still complete and pass the exact-point Id-Vd/KCL gates before Stage 4 is
declared qualified.  Vg=8 V remains gated on that result.  IALMob and predictor
remain disabled.

## 2026-09-04 follow-up: complete HFS stencil and QF-coordinate floor

The omitted contact-cell HFS columns have now been implemented for every
non-endpoint potential DOF in the cell-gradient stencil.  Replaying the same
60 one-ring potential directions reduced the maximum full-Jacobian JVP error
at the initial, best and final frozen states to `8.073e-12`, `6.195e-11` and
`2.522e-10`, respectively.  The previous `7.224e-4`--`2.282e-3` discrepancy
is therefore closed without changing the residual or HFS driving-force
contract.

The full curve then exposed a second, independent representation floor near
`Vd=0.1625 V`.  At the leading electron rows (4084, 4096, 4097, 4118, 4245 and
4246), the requested `phin` corrections were only `0.84e-17`--`2.15e-17 V`,
or about 2.6--6.2 ULP of the stored absolute QF value.  Applying the rounded
step flipped the leading residual signs and increased the electron block from
`1.0994e-11` to `1.4008e-11`.  Extended-precision summation and exact replay
of every incident SG edge did not change the residual, excluding nodal
accumulation order as the cause.

An opt-in `quasi_fermi_recenter_on_initial_state` mode now places each
warm-start's QF reference on the supplied physical state.  It changes only the
internal coordinate representation; reconstructed absolute `phin/phip`,
carrier densities and the residual equations are unchanged.  With this mode:

- the previously rejected `0.1625 V` point converged under the unchanged
  `1e-11` electron ceiling;
- all 25 explicit points from `0.1625` through `0.2 V` passed with the complete
  live HFS Jacobian;
- the qualified lagged-HFS and live-HFS solutions at `0.2 V` differ by only
  `1.066e-17 A/um` (`5.924e-13` relative) in drain current;
- an exact no-predictor production path has passed through `0.8 V`.

The production D5 generator now freezes the HFS field derivative as the
qualified quasi-Newton strategy, recenters QF coordinates at every warm start,
and freezes `min_step=1e-3 V`, `max_retries=12`.  The residual continues to use
live HFS mobility.  The complete live stencil remains covered by the device
JVP audit and Catch2 regression.

The 31-point gate is not yet complete.  The exact no-predictor sweep currently
requires roughly 1.5--3.125 mV accepted internal steps in the startup region
(64 accepted substeps were needed for `0.4 -> 0.6 V`).  A new default-off
`write_state_every_accepted_step_prefix` option now persists these internal
states.  A real-device `0.6 -> 0.65 V` test wrote all 17 accepted states,
including the 15 non-requested intermediate points, and completed normally.
The ordinary `write_state_file` is also refreshed after each internal accepted
step, providing a bounded-size rolling restart checkpoint for production.
The remaining issue is therefore runtime rather than lost progress or the
former 10 mV/0.1625 V convergence blockers; no partial curve is scored.

A `growth_factor=2` control from `0.65 V` repeatedly proposed 5--10 mV steps
and fell back to 2.5 mV, making normalized progress slower.  The fixed
`growth_factor=1` control then completed all 49 accepted states from `0.68` to
`0.8 V`.  Production therefore keeps growth disabled; this negative A/B is not
used to weaken convergence gates or enable a predictor.

## Remaining work

1. Use the accepted-step checkpoints to finish and score the Vg=4 V exact-point
   curve; run Vg=8 V only if it passes.
2. Audit why the strict no-predictor path needs millivolt-scale steps above
   `0.2 V`; do not weaken the frozen block residual ceilings.
3. Audit the near-null carrier mode exposed by the uncapped final-state solve.
4. Only after D5 curve qualification resume IALMob development.

Machine-readable results are in
`reference_tcad/templates_ldmos_sentaurus2022/stage4_vd10mv_newton_audit_summary.json`.
Large probe outputs remain in ignored `reference_staging` directories.
