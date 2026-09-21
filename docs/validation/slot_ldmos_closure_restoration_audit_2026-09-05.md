# SLOT-LDMOS 1.25 V closure-restoration audit

Date: 2026-09-05

Follow-up: [E3 budget and strict 1.25 V checkpoint](slot_ldmos_e3_strict_1p25_2026-09-05.md)
completes the next steps below. It identifies auxiliary gauge rows as the
contact/source versus free-sum discrepancy, validates an effective correction,
and obtains a restart-verified strict 1.25 V checkpoint. The no-checkpoint
statements below describe this earlier E4/E1/E2 experiment only.

## Outcome

The deterministic E4/E1/E2 replay does **not** support stabilized-Jacobian
contamination as the dominant failure mechanism at the saved 1.25 V state.
The nominal minimum-norm correction solves both linear constraints, but adding
it to the stored double state loses most of the effective electron correction.
The reassembled free-row sum follows the physical Jacobian prediction of the
**actually represented** update. Contact/source closure and the free-row sum
also separate at this accuracy. This is evidence for a representability audit
and the two-level E3 assembly/reduction audit, not grounds to relax a tolerance.

No strict 1.25 V checkpoint or Vela BVDS result was obtained. All tested final
restoration candidates fail the unchanged electron global-closure gate.

## Code and input baseline

- Worktree: `.worktrees/slot-ldmos-bvds-next`.
- Branch: `codex/slot-ldmos-bvds-next`.
- Local main tip merged: `8e99e57`.
- Merge commit: `dbe6d3d` (no conflicts).
- Original handoff tip: `77202ad`.
- Frozen input: `build/slot_ldmos_noavalanche_1p25_resume02/outputs/stages/01_delaunay_noavalanche_adaptive/attempt_001_v_1p25/rejected_states/attempt_1_bias_1p250000_final.csv`.
- Input SHA-256: `13bba7c57646d25b57eeb3f0fb188c56a950168acc1aa83bbd3dc801001f5edc`.
- Original rendered config SHA-256: `98bdae0f838ef9434eb2688190b03682278303dc4f7d0a9e0df74c103aba2a04`.

The Delaunay/mixed-Voronoi geometry, Masetti field mobility, recombination,
Old Slotboom, QF limits (0.1/0.05 V), diagonal floor (scale 10), outer row
regularization (1e-8), and all convergence thresholds are retained. Avalanche,
IALMob and external resistance remain disabled. The absolute closure limit is
1e-12; the relative limit is 1% above the 1e-11 source/current transition.

Two fixed-state configurations use zero Newton iterations solely to isolate
final restoration. They differ only in restoration Jacobian mode and output
paths. These zero-iteration experiments are not physical checkpoints. A
separate four-iteration replay starts from the original resume02 input and
reproduces the handoff final scalar residual, 1.2948077136185354e-9.
The persisted four-step final state is byte-identical to the frozen input above
(the SHA-256 hashes match), confirming trajectory parity after the local-main
merge and with the final diagnostic implementation.

## E4: replacement rows

Both carriers replace the equation at node 0, located at
(7.092727668385104, 0) um. With N=13881, the electron/hole equation indices are
13881 and 27762. The selected electron density is below intrinsic density.

| Accepted Newton iteration | Electron physical residual | Electron rank | Hole physical residual | Hole rank |
|---|---:|---:|---:|---:|
| Initial | -4.4150166470e-13 | 308 | 3.5076607373e-23 | 9216 |
| 1 | -4.4066647876e-13 | 311 | 9.9931148080e-23 | 9144 |
| 2 | -4.4016848934e-13 | 313 | 3.7876456637e-23 | 9376 |
| 3 | -4.4010159442e-13 | 313 | 5.4799632597e-23 | 9325 |
| 4 | -4.4006788588e-13 | 313 | 5.9275505216e-23 | 9319 |

Ranks compare absolute assembled residuals within each carrier block. Neither
replacement row is in the top ten or is node 6248, 5849 or 8063. This does not
prove row replacement harmless for all states, but it supplies no evidence for
redesigning the linear system around these reported hotspots.

## E1/E2: final-restoration A/B

The initial physical electron/hole mismatches are 1.846541960790157e-12 and
3.4286899885251264e-12, both in absolute mode. The frozen-mode ideal alpha lower
bound is approximately 0.7083. Only alpha=1 in the existing halving sequence is
ideally feasible for both carriers; the other twelve trials are labelled.

| Quantity | Stabilized sensitivity | Physical sensitivity |
|---|---:|---:|
| Electron row norm | 180.9537540567 | 180.9537307533 |
| Hole row norm | 121.6899494941 | 121.6899418589 |
| Normalized two-row condition number | approximately 1 | approximately 1 |
| Projected correction L2 | 2.9966604936e-14 | 2.9966607046e-14 |
| QF capped DOFs | 0 | 0 |
| Alpha=1 scalar residual | 1.2948098501e-9 | 1.2948098501e-9 |
| Alpha=1 physical electron mismatch | 1.8465419608e-12 | 1.8465419608e-12 |
| Alpha=1 physical hole mismatch | 3.4058988818e-13 | 3.4058971239e-13 |
| Alpha=1 rejection | global_closure | global_closure |

The correction-norm difference is about 7.04e-8 relative. Both normalized row
singular values are approximately 1, sin(theta)=1, and the multipliers are of
order 1e-14. Neither conditioning nor cap activation explains this replay.
All thirteen candidates pass the device floor. Contrary to the proposed
failure chain in the original handoff, alpha=1 does not break the device gate.
The configured carrier-row gate passes with zero qualified rows at these
candidates; this does not establish accuracy of every local continuity row.
The final-only corrector does not invoke Poisson recorrection.

### Nominal versus represented update

For stabilized sensitivity at alpha=1:

- Nominal electron step L2: 1.0204496559993189e-14.
- Actually represented electron step L2: 2.1628695523766688e-15.
- Nominal and represented hole step L2: 2.8175621756594644e-14.
- Total addition-rounding error L2: 9.972649865623315e-15.
- Nominal predicted electron free-row sum: approximately zero.
- Prediction using the represented step: -1.76358808775e-12.
- Physical-Jacobian prediction using the represented step: -1.76358808858e-12.
- Reassembled electron free-row sum: -1.76358808858e-12.

Only about 4.49% of the ideal electron closure change survives state addition.
At node 10237, the electron solver coordinate is 48.35215883979197 and the
nominal delta is -2.8704875976464014e-15. The directional ULP is
7.105427357601002e-15. This delta is below half an ULP and is rounded away.
The trace records eight largest lost updates, their node/carrier identities,
state values, deltas and directional ULPs. Counts of lost updates include
arbitrarily small cross-carrier components, so the sensitivity loss is the
more useful measure of impact.

This is new evidence about final-correction representability. It does not
invalidate the earlier contact-bias lifting or restart-coordinate audits.

### Contact/source closure versus assembled free sums

At alpha=1, the physical-Jacobian prediction and the reassembled free-row sum
agree closely for the represented update. The authoritative contact/source
closure is still a different numerical observation:

| Carrier | Reassembled free-row sum | Contact/source mismatch | Sum of these quantities |
|---|---:|---:|---:|
| Electron | -1.7635880886e-12 | 1.8465419608e-12 | 8.2953872210e-14 |
| Hole | -2.1589481116e-19 | 3.4058988818e-13 | 3.4058967228e-13 |

The audit's `free_sum_inferred_mismatch` is the initial contact/source mismatch
minus the predicted change in the physical-Jacobian free-row sum. It is not
an independently differentiated contact/source functional. This distinction
prevents interpreting the observed gap as Jacobian contamination.

These observations do not measure the complete assembly noise floor. E3 must
still separate fixed-residual reduction error, per-node assembly/scatter error,
and repeatability before any change to the 1e-12 contract is considered.

## Implementation and validation

- Opt-in `global_continuity_closure.audit_jsonl_file` appends replacement-row,
  trigger and trial records, including failed candidates.
- `restoration_jacobian: physical` excludes both forms of artificial carrier
  stabilization only from final restoration; the main Newton solve is unchanged.
- The default remains `stabilized`. The Gram solve, QF cap order, alpha schedule,
  scalar/device gate, carrier gate and physical-closure gate are unchanged.
- The assembler accepts an explicit option to omit its diagonal floor for the
  physical Jacobian. A test compares it with a floor-disabled assembler and
  verifies that subsequent stabilized assembly and residuals remain unchanged.
- Catch2 covers exact projection, near dependence, degeneracy, post-cap
  predictions, the ideal feasible-alpha bound, configuration validation,
  unchanged Newton trajectories with auditing, replacement rows and rejected
  final nonlinear candidates.
- Final Debug build: successful with MSYS2 UCRT64.
- Final `test_newton_solver.exe`: **1571 assertions in 111 test cases passed**.

Local artifacts are under `build/slot_ldmos_closure_audit_20260905/final/` in
`stabilized`, `physical`, and `newton_replay`. Each contains `config.json`,
`closure_audit.jsonl`, `iv.csv`, iteration/failure diagnostics and `run.log`.
`summary.json` at the final artifact root consolidates the verified records.
The parent `manifest.json` seals the input and runner hashes. Initial
exploratory records outside `final/` are superseded, including their unconverted
coordinate labels. Generated outputs remain ignored and are not committed.

Rebuild and run from this worktree in PowerShell:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --preset windows-ucrt64-debug
cmake --build build --target test_newton_solver vela_example_runner --parallel 2
Push-Location build
.\test_newton_solver.exe
Pop-Location
.\build\vela_example_runner.exe --config build/slot_ldmos_closure_audit_20260905/final/stabilized/config.json
.\build\vela_example_runner.exe --config build/slot_ldmos_closure_audit_20260905/final/physical/config.json
.\build\vela_example_runner.exe --config build/slot_ldmos_closure_audit_20260905/final/newton_replay/config.json
```

Runner exit code 1 is expected for these deliberately unresolved closure
experiments. A failed IV row's zero terminal-current columns remain unavailable
results, not physical zero-current measurements. Use fresh output paths when
repeating, since JSONL diagnostics append.

## Revised next task

1. Complete E3 at this exact saved state: deterministic alternate reductions,
   reversed/permuted edge scatter, compensated per-node assembly, repeated
   assembly, and a separate contact/source versus free-row accounting check.
2. Audit representable final QF updates in the drain basin. Compare the predicted
   effect of nominal and representable corrections, including ULPs and closure
   sensitivity. Consider a precision-preserving reference/increment scheme only
   after establishing its end-to-end assembly and restart semantics.
3. Keep the physical Jacobian mode as an explicit diagnostic A/B option; this
   replay does not justify promoting it as a numerical fix.
4. Revisit cap-aware projection, row conditioning or row replacement only if new
   evidence activates those hypotheses. Preserve all existing physics and gates.
5. Resume strict 1.25, 1.5 and 2 V continuation only after an actual qualified
   correction is demonstrated. Avalanche, IALMob and high-voltage work remain
   downstream of those gates.
