# Templates/LDMOS Stage-4 Vg=4 V, Vd=10 mV Newton/Jacobian audit

Date: 2026-09-03

## Outcome

The first nonzero-drain D5 transfer is no longer blocked.  The audit found and
fixed a configuration-contract defect: with Fermi-Dirac carrier statistics,
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

## Remaining work

1. Finish and score the Vg=4 V exact-point curve, then run Vg=8 V if it passes.
2. Add the missing non-endpoint potential stencil for contact-cell HFS mobility
   and qualify it independently; do not silently change the current residual.
3. Audit the near-null carrier mode exposed by the uncapped final-state solve.
4. Only after D5 curve qualification resume IALMob development.

Machine-readable results are in
`reference_tcad/templates_ldmos_sentaurus2022/stage4_vd10mv_newton_audit_summary.json`.
Large probe outputs remain in ignored `reference_staging` directories.
