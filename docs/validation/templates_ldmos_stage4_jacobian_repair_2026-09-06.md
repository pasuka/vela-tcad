# Templates/LDMOS Stage-4 mobility Jacobian repair

Date: 2026-09-06. The targeted repair and local controls are complete. Vg=4 V
D5 continuation is running; this report does not qualify the full 31-point
curve. The previous [numerical controls](templates_ldmos_stage4_numerics_validation_2026-09-06.md)
remain a separate, unchanged record of the pre-repair behavior.

## Result and implementation

The fine-perturbation error plateau at node 5569 came from differentiating
the entire live HFS transport flux with a voltage secant that crossed the
narrow mobility transition. The edge-projected bulk QF path now differentiates
the mobility multiplier analytically and retains the existing finite
difference of the frozen-mobility SG flux. At the archived failed state,
the electron JVP relative error at h=1e-9 V falls from 1.52565e-4 to
2.80279e-10. This fixes the identified plateau; it is not a claim that every
Jacobian block is exact.

For each adjacent transport cell, write
`mu_k = mu0_k * (1 + r_k^beta)^(-1/beta)`, where
`r_k = mu0_k * |Delta_phi_qf| / (h * vsat)` with the configured field-unit
conversion. The code uses the same cached low-field mobilities and arithmetic
cell average as the residual. If `F = mu * F0`, the added feedback is

```
dF_mobility / d(Delta_phi_qf)
  = -(F / Delta_phi_qf) * mean(mu_k * r_k^beta / (1+r_k^beta)) / mu
```

The two endpoint QF columns receive opposite signs. At zero QF drop, the
feedback has zero limit. The existing surface, cell-vector and contact-field
fallback stencils remain in use outside this bulk edge-projected path.
`constant_field` now participates in field differentiation even when contact
fallback is disabled; its omission from the selector was a separate defect.

Implementation: [CoupledDDAssembler.cpp](../../src/equation/CoupledDDAssembler.cpp).
Contract documentation: [mobility configuration](../config_schema.md).
The configured physical mobility, nonlinear residual, convergence ceilings,
production template defaults and low-bias derivative policy were not changed.
Diagnostic and script changes already present in the worktree were preserved.

## Reproduction context

- Worktree: `D:\code-repo\vela-tcad\.worktrees\templates-ldmos-phase-a`.
- Branch: `codex/templates-ldmos-phase-a`; HEAD
  `345fda30d541003f4aa3a08341dca59418f06ded`. Changes are uncommitted.
- Exact imported mesh: 10241 nodes, 30022 edges, 19782 triangles; 16237
  qualified external AverageBox transport couples. Barycentric node volumes,
  material-local Poisson charge and `legacy_node_local` contacts are retained.
- 300 K, Fermi statistics, OldSlotboom, SRH/Auger, `constant_field` HFS,
  edge-projected QF drive and the existing contact electric-field fallback.
  Predictor, IALMob, avalanche, quantum correction and heating remain off.
- Units: input geometry in um, `unit_scaling`; reported terminal current is
  `current_total_A_per_um`. Raw block residuals are solver quantities and must
  not be interpreted as terminal currents or compared directly to SDevice RHS.
- Block gates: psi <=5e-8, electron <=1e-11, hole <=3e-10, enforced.
  Newton budget 160. Physical parameters, damping/update-limit settings and
  the residual filter are identical to the prior controls.
- Windows UCRT64 Release preset. HDF5 and SuiteSparse capabilities are detected;
  actual device controls explicitly select Eigen SparseLU, COLAMD ordering and
  L2 row/column equilibration.
- Repaired runner SHA-256:
  `89ed3a95e828e4190adf460dfa0881fd90a8749bef203bb4ed5adc02ffa757fa`.
  The previous runner was preserved as
  `build-release/vela_example_runner_pre_chainrule_20260906.exe`, SHA-256
  `04a3c2ac2282819abef4d4bbe383e54c005e349dab54caa2228aebf1a17e5076`.

Generated evidence is in
`reference_staging/templates_ldmos_jacobian_repair_20260906/`.
The [analysis script](../../reference_staging/templates_ldmos_jacobian_repair_20260906/analyze_repair.py)
checks unchanged D5 gates, matched checkpoint hashes, frozen-path equivalence
and the exact remaining reference-voltage grid. Its aggregate is
[repair_summary.json](../../reference_staging/templates_ldmos_jacobian_repair_20260906/repair_summary.json).

## Hotspot derivative validation

The same archived parent and failed states, 12-node neighborhood, three
potential directions and perturbation scan were reused: 204 JVP evaluations.
Lagged Jacobians use a matched frozen-mobility residual perturbation. The
table uses `absolute_error / finite_difference_norm`; the diagnostic CSV's
floor-normalized `relative_error` column is not substituted for this metric.

| Node 5569 electron-QF direction, failed state | Old live error | Repaired live error |
| --- | ---: | ---: |
| h=1e-6 V | 2.066e-13 | 1.52588e-4 |
| h=1e-7 V | 1.43923e-4 | 8.65821e-6 |
| h=1e-8 V | 1.52540e-4 | 2.50010e-8 |
| h=1e-9 V | 1.52565e-4 | 2.80279e-10 |

The old h=1e-6 agreement arose because the diagnostic reused the production
secant scale. With an independent analytic mobility derivative, coarse
perturbations correctly reveal nonlinear finite-difference truncation error;
the fine-step plateau disappears. The matched frozen cases remain identical.
The small minority-hole directions still show large relative errors in tiny
absolute terms, and the unchanged psi/contact/surface/vector stencils have
not received a new general accuracy qualification.

## Low- and high-bias activation controls

Each pair is one direct fixed-target solve from an identical checkpoint,
without a preliminary reclose, predictor or adaptive intermediate points.
The manifest records the physical parent-to-target delta; the single-point
runner's startup `attempted_step_V=0` is bookkeeping, not a zero voltage change.

| Parent -> target Vd, Vg=4 V | Frozen mobility derivatives | Repaired live derivatives |
| --- | ---: | ---: |
| 0 -> 0.0025 V, accepted zero-bias state | pass, 3 updates | pass, 13 updates |
| 0.2 -> 0.2025 V | pass, 8 updates | pass, 12 updates |
| 1.33333333333333 -> 1.33583333333333 V | pass, 26 updates | pass, 12 updates |
| 3.369166666666683 -> 3.371666666666683 V | pass, 85 updates | pass, 13 updates |
| 3.851666666666673 -> 3.854166666666673 V | fail after 153 updates, next line search rejected | pass, 13 updates |

The corresponding live/frozen Id relative differences for the four passing
pairs are 2.4071e-10, 3.4925e-12, 2.4980e-14 and 2.4014e-13. At the floor
transfer the live electron residual is 9.30655e-12 and normalized terminal
KCL is 1.83674e-13. The old live path needed 14 updates there; most of the
85-to-13 improvement comes from retaining mobility feedback, not from this
additional accuracy repair alone.

The frozen high-bias iteration CSVs and all emitted final/rejected state CSVs
are byte-identical to the pre-repair controls. Thus the high-bias frozen
baseline and its failure remain reproducible.

These samples support retaining frozen derivatives at startup/low bias and
using repaired live derivatives for the current high-bias continuation. They
do not establish an optimal global switching voltage, so no automatic switch
or all-bias default change was introduced.

Two initial 0->10 mV seed controls and a second pair using the accepted 0 V
checkpoint both failed. The former used a gate-prebias seed, not a qualified
D5 parent. They remain archived and are not conflated with the controlled
2.5 mV comparisons above. A 10 mV direct transfer also omits the intermediate
steps of the historical successful 0->10 mV sweep.

## Floating-point update resolution

This was tested independently with derivatives frozen, using the archived
failed final state at 3.8541666666666727 V and solving again at that same bias.

| QF reference arrangement | Result | Final electron residual |
| --- | --- | ---: |
| contact-basin references, initial-state recenter off | failed after 29 updates | 8.15742e-11 |
| initial-state recenter on | passed after 1 update | 4.10550e-13 |

Both read the same checkpoint, SHA-256
`16519194daa3624d56148e4d14be34264dea111bd28fe02f120b2d12a52acea4`.
Repartitioning reference plus increment changes finite-precision coordinates;
this is not an assertion that both internal initial residuals or coordinate
vectors are identical. The saved split representation is preserved during
loading and repartitioning.

In the recentered case, node 5569 starts with electron residual 5.08793e-12.
Its raw 1.93183e-17 V correction is actually applied, and the trial row
residual drops to 1.14081e-15. The accepted solution has Id
1.8299877286054715e-4 A/um and normalized KCL 1.32253e-15; its Id differs
from the live-Jacobian accepted floor solution by only 1.82077e-13 relatively.

This establishes a strict same-bias restart/recenter remedy for this archived
floor, separate from mobility derivatives. It uses the existing
`quasi_fermi_recenter_on_initial_state` implementation. No arbitrary minimum
update, rounded-up correction, residual exemption, automatic mid-Newton
reference change or new solver heuristic was added. It is not a universal
proof that recentering resolves every subsequent stall. Production continuation
still starts from accepted checkpoints, rather than promoting a failed state.

## D5 continuation

The accepted 3.975416666666673 V production checkpoint was reclosed and
continued with live derivatives, 2.5 mV initial steps, growth factor 1,
minimum retry step 1 mV and the unchanged gates. Vd=4 V was accepted with
20 updates in the final internal transfer; the whole short segment used
135 updates across 11 accepted solves, including its startup reclose, and
had zero rejected attempts.

At the exact Vd=4 V point:

- psi/electron/hole residuals: 9.80882e-10 / 7.99897e-12 / 3.65150e-25.
- Id=1.857642151059825e-4 A/um, versus Sentaurus 1.68591810156894e-4 A/um:
  +10.18579%. This is one exact point, not a full-curve median/P95 verdict.
- KCL divided by the largest terminal magnitude: 3.96566e-13
  (3.96566e-11 percent).
- Accepted checkpoint: [point_bias_4p000000.csv](../../reference_staging/templates_ldmos_jacobian_repair_20260906/near4_fine/point_bias_4p000000.csv),
  SHA-256 `215c03232802296b39b4e4560707b3c4208dc50c9b8dfceda42b34209de54904`.

An isolated 24.5833 mV direct transfer to 4 V failed. From the accepted 4 V
checkpoint, an isolated 10 mV transfer also failed, while the 2.5 mV sweep
passed 4.01 V with 53 updates across five solves including startup. The
large-step equivalence controls therefore did not pass and no larger
production step was adopted.

The ongoing [continuation config](../../reference_staging/templates_ldmos_jacobian_repair_20260906/vg4_continue_4p01_to_40/reclose.json)
starts from the accepted 4.01 V checkpoint and retains every remaining exact
reference point from 5.33333333333333 V to 40 V. Its directory contains
accepted-step checkpoints and live `reclose.log`; attempt/iteration CSVs may
be buffered during execution. `summary.json` is written on runner completion.

At 13:43 China time on 2026-09-06, process 7356 was running, with a complete
accepted checkpoint at 4.1725 V. This is a timestamped observation, not a
permanent process-status claim. Inspect the directory for a newer checkpoint
or completion summary before resuming; do not launch a duplicate scan.
The read-only `inspect_progress.py` helper and `progress_snapshot.json` in
the evidence directory distinguish complete checkpoints from buffered CSV
rows; verify the operating-system process as well as the files.

The full 31-point Vg=4 score, zero-bias KCL issue, Vg=8 curve and two-gate ratio
remain outstanding. Preserve the zero-bias KCL ratio in scoring; no exemption
was approved. Vg=8 remains gated on Vg=4 qualification, and IALMob remains gated
on D5 qualification. On a failed transfer, use the last accepted checkpoint
and the rejected-state/iteration evidence; do not score an interpolated or
unconverged point or loosen the residual gates.

## Verification

- Matching UCRT64 Release configuration and affected targets built successfully.
- New narrow-transition JVP test: both carriers, both drop signs and two
  independent perturbation sizes; 16 assertions pass at relative tolerance
  2e-6. Its deliberately low test saturation velocities expose the secant
  defect without changing device physics.
- Newton solver: 111 cases / 1460 assertions passed.
- DC sweep: 98 cases / 3479 assertions passed.
- Material-local Poisson charge: 5 cases / 14 assertions passed.
- Templates/LDMOS Phase 2/3 script regressions: 27 tests passed.
- 204 real-device JVP controls, five low/high-bias A/B pairs, independent
  precision restart controls and two fine-step continuation segments completed.
  Deliberately larger-step failures are retained as negative evidence.
- Frozen trace/state equivalence and unchanged D5 gate/grid assertions passed.
  Local report links and `git diff --check` were checked. No changes committed.

Replay helpers and summaries are generated, ignored artifacts; their local
links require this workspace's fixtures. Source, test and schema changes are
reviewable in the worktree diff.
