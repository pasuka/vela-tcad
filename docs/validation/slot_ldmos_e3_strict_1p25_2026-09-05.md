# SLOT-LDMOS E3 error budget and strict 1.25 V checkpoint

Date: 2026-09-05

## Outcome

The error budget is complete for the frozen 1.25 V failure and its final
restoration trial. The effective correction combines the existing
`quasi_fermi_reference: "contact_basin"` coordinates with the new opt-in
`global_continuity_closure.restoration_physical_rows_only: true`. A strict
1.25 V checkpoint passes the unchanged handoff acceptance contract after two
real Newton iterations. Avalanche, IALMob and external resistance remain off.
This is a low-voltage checkpoint, not a Vela BVDS result.

This report supersedes the pending E3 work and no-checkpoint status in the
[earlier E4/E1/E2 audit](slot_ldmos_closure_restoration_audit_2026-09-05.md).
The local main merge remains `dbe6d3d`, incorporating main tip `8e99e57`.
The implementation and this report are worktree changes beyond that commit.

## Frozen experiment and error budget

All paths below are relative to the worktree. Evidence is under
`build/slot_ldmos_e3_20260905/`. The input is the rejected 1.25 V final state
from `slot_ldmos_noavalanche_1p25_resume02`, SHA-256
`13bba7c57646d25b57eeb3f0fb188c56a950168acc1aa83bbd3dc801001f5edc`.
It has 13881 nodes. Mesh, doping, material files and protected configuration
are hashed in `manifest.json` in the evidence directory.

The E3 diagnostic first freezes the assembled residual and compares ordinary,
pairwise, Kahan, Neumaier, long-double, increasing/decreasing magnitude and
deterministically permuted sums against a Boost quad reference. Separately,
it changes the actual edge scatter order (forward, reverse, fixed permutation)
and uses ordinary or compensated per-node accumulation. Every assembly is
repeated to check bitwise determinism. These are diagnostics, not changes to
the Newton residual or its thresholds.

| Measurement | Electron | Hole |
|---|---:|---:|
| Initial physical closure mismatch | 1.8465419608e-12 | 3.4286899885e-12 |
| Initial fixed-vector maximum reduction error vs quad | 1.0909079921e-23 | 5.5463816397e-24 |
| Alpha=1 fixed-vector maximum reduction error vs quad | 8.1341291288e-23 | 3.5542545350e-23 |
| Initial actual assembly global-sum span | 1.8789278367e-24 | 1.6410326313e-24 |
| Alpha=1 actual assembly global-sum span | 1.8787258883e-24 | 1.6407274007e-24 |
| Alpha=1 auxiliary gauge free-row sum | 8.2953872210e-14 | 3.4058967228e-13 |
| Inactive free gauge rows | 1471 | 1471 |

Closure and residual sums above use the solver's existing internal units,
not terminal A/um. No physical unit reinterpretation is made. The maximum
per-node full-residual change is 4.8260480056e-15; this includes Poisson rows
and must not be compared directly to a carrier global sum. Repeated identical
assemblies have zero bitwise changes. Contact/source mismatch is unchanged
across all six assembly variants at each frozen state. The quad reference
quantifies reduction of the already assembled double values; it is not a
quad-precision reevaluation of the device physics.

The large contact/source versus free-row discrepancy is explained by inactive
gauge equations. The old projection can reduce the sum of physical and gauge
residuals by moving unconstrained auxiliary carrier coordinates, while leaving
physical imbalance. This is not a roundoff floor. The measured reduction and
assembly variations give no reason to loosen the absolute closure limit 1e-12.

## Correction and controlled A/B

The final physical-row projection excludes inactive gauge rows from the two
constraint sums and excludes their unknowns from the minimum-norm QF update.
Dirichlet rows/unknowns remain excluded. Its default is false; the main Newton
equations, row replacement, diagonal floor and outer regularization are unchanged.
The successful mode retains the stabilized final sensitivity. Removing its
stabilization is not required for this checkpoint.

Contact-basin coordinates keep the dominant electron update representable.
The prior audit showed a nominal electron step smaller than half an ULP of a
stored biased coordinate. With the combined correction, electron nominal and
actually applied L2 norms both measure approximately 1.0334e-14. Remaining
lost updates are much smaller and do not prevent closure. No QF component is
capped in the successful final trial.

The following fixed-state tests use zero Newton iterations to isolate the
corrector. They are experiments, not strict checkpoints. For failing modes the
table reports the alpha=1 trial, not the unchanged failed-state IV record.

| Reference coordinates | Physical rows only | Electron mismatch | Hole mismatch | Final trial |
|---|---|---:|---:|---|
| None | No | 1.84654e-12 | 3.40590e-13 | Fails electron gate |
| Contact basin | No | 8.29541e-14 | 3.40590e-13 | Passes |
| None | Yes | 1.84654e-12 | 2.35925e-19 | Fails electron gate |
| Contact basin | Yes | 2.48575e-19 | 2.35925e-19 | Passes |

Contact-basin coordinates alone also pass a separate four-Newton-iteration
test (`basin_only_strict`). The combined mode is retained because it addresses
both the lost update and the incorrectly included auxiliary constraints.
Neither change alters physics or increases a tolerance.

## Strict checkpoint and restart

Final run: `build/slot_ldmos_e3_20260905/final_strict/`.

| Quantity | Measured value |
|---|---:|
| Drain bias | 1.25 V |
| Real Newton iterations | 2 |
| Convergence reason | `final_global_continuity_restoration` |
| Drain terminal current | 6.5889954852964324e-16 A/um |
| Normalized scalar residual | 1.2946469730360184e-09 |
| Electron physical closure mismatch | 2.4601233527902148e-19 |
| Hole physical closure mismatch | 2.3598226144186885e-19 |
| Electron/hole closure limit in this regime | 1e-12 |
| Contact/source versus free-sum gap, electron/hole | 4.728059e-25 / 5.512098e-25 |

Accepted state: `final_strict/final_state.csv` within the evidence root.
SHA-256: `8de098c8b1416210a0ef05e1e4269104a8f555e48d01ef99e9a308d96818b10a`.

Strict means the six requirements in section 10 of the handoff: real Newton
iterations, the existing carrier-row gate, the unchanged two-regime global
closure gate, a saved restartable accepted state, an explicit convergence
reason, and unchanged physical settings. The carrier-row gate has **zero
qualified rows** at its existing source threshold, so this result is not a
claim that every local row passes an independently enforced absolute test.
The existing IV diagnostic also reports 1186 QF-bound warnings; these are
reported by the existing validation path and are not a rejecting gate in this
contract. The configured per-step QF limits remain 0.1/0.05 V.

Restart testing exposed an additional precision issue: computing
`old_reference + increment - new_reference` loses part of a tiny increment
even with long-double intermediates. The initial restart passed closure but
changed terminal current by 3.33e-5 relatively. Reference repacking now computes
`(old_reference - new_reference) + increment` in both solver warm starts and
diagnostic packing. Matching references cancel exactly before the increment is
added. Both electron and hole regression tests retain 1e-22 V increments at
nonzero reference biases. The final checkpoint and restart below use this fix;
the earlier `combined_strict`/`combined_restart` artifacts are exploratory only.

The final same-bias reload (`final_restart`) passes
`initial_stall_residual_floor` with zero additional Newton iterations. This is
a reload verification of the two-iteration checkpoint, not a second independent
strict checkpoint. Terminal current is exactly unchanged at the reported double
precision; all saved electron and hole densities are also exactly unchanged.
Electron/hole closure mismatches are
`2.4601233689460862e-19` /
`2.3598225982628171e-19`.
Maximum electron/hole increment changes are 2.096e-31 / 5.275e-32 V, consistent
with conversion to and from scaled coordinates. The redundant reconstructed
physical electron QF column changes by at most 8.882e-16 V; the stored reference
and increment fields are the precision-preserving restart coordinates.
`restart_validation.json` records the complete column comparison.
Restart state SHA-256: `d4871cbdb9df24cbf583d09db7d735fe914dfc5c946d520014c75097add6007c`.

## Validation and reproducibility

- UCRT64 Debug build: `test_newton_solver` and `vela_example_runner`.
- All Newton tests pass: 1610 assertions in 114 cases, including assembly-order
  auditing, gauge exclusion, default behavior and referenced restart precision.
- Test log: `build/restart_precision_tests.log`.
- Budget summary: `build/slot_ldmos_e3_20260905/noise_budget.json`.
- Fixed-state evidence: `baseline`, `basin_only`, `physical_rows_only`, `combined`.
- Final checkpoint and same-bias reload evidence: `final_strict`, `final_restart`.
- Each case contains its rendered `config.json`, `iv.csv`, audit and runner logs;
  accepted cases include saved state CSVs with reference and increment columns.
- The tracked [restoration contract](contracts/slot_ldmos_physical_closure_restoration_v1.json)
  records required options, gates, frozen input hashes and final evidence hashes.

From a PowerShell session in this worktree, prepend
`D:\msys64\ucrt64\bin;D:\msys64\usr\bin` to PATH, build the two targets, then run:

```powershell
.\build\vela_example_runner.exe --config build/slot_ldmos_e3_20260905/final_strict/config.json
.\build\vela_example_runner.exe --config build/slot_ldmos_e3_20260905/final_restart/config.json
```

Audit JSONL appends, so a new replay should redirect outputs into a fresh case
directory to keep evidence runs distinct. Generated outputs remain ignored in
`build/` and are not committed.

## Next gate

Continue from the saved strict 1.25 V state to 1.5 V and then 2 V under the
same restoration contract. Each new voltage needs real Newton iterations and
the same acceptance evidence. Only after those gates should avalanche-off
voltage continuation, IIC, self-consistent avalanche and external-resistor BVDS
resume. IALMob comparison remains later work.
