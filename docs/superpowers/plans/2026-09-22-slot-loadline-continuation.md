# Reproducible SLOT external-load-line continuation

User authorization: execute stage 2 of the approved sequential task list in
`docs/validation/slot_ldmos_next_tasks_2026-09-22.md`. Stage 1 direct-bias
drivers have completed 20/37 V with original gates. This is a bounded
automation of the existing successful projected-SG 37.02/37.05 V procedure,
not a new physical model or a relaxation of any final qualification rule.

## Design decision

Add a focused Python driver around existing direct device helpers. Keep the
C++ circuit/physics implementations unchanged. Compared with extending the
one-direction direct-bias CLI, a separate scalar controller makes inner and
outer voltage explicit and can recover in either inner-voltage direction.
Repairing the C++ coupled globalization remains a separate investigation.

The research wrapper's single-attempt budget abandoned an already-queued
midpoint after the 37.1 V failure. Formal recovery must preserve failed evidence,
reduce the inner step toward its qualified seed, and never promote a failed
candidate. Distinguish qualified device evaluations from accepted external
checkpoints. Only the latter satisfy the load line and circuit confirmation.

Use signed `current_total_A_per_um`, R=1e12 ohm*um and residual tolerance
1e-6 V for this campaign. Root proposals use bounded secants; recognize a
bracket only from opposite signs. Fall back to bounded exploration/midpoints
for singular slopes, out-of-bracket steps and numerical failures. Do not infer
monotonicity, uniqueness, or a root from same-sign samples. Bound attempts,
retries, minimum representable steps, and wall time per solver process.

An imported seed is requalified by a direct full solver and all-row audit at
its known inner bias before campaign acceptance. This avoids trusting old
`qualified` booleans or placeholder circuit residuals. Reuse identical-bias
qualified states; each new inner bias requires at least one cumulative real
coupled Newton update on its own candidate lineage. Poisson and carrier-only
updates are not coupled updates. All original full solver/global closure and
local gates apply; validate active node IDs, not just row counts.

Each external acceptance must pass recomputed signed load-line residual and
confirmation under the same physical deck with external circuit enabled.
Byte-identical confirmation states inherit the bound direct evidence. A
changed state must undergo independent full requalification and load-line
recomputation; it cannot silently inherit the old state evidence.

Persist atomic checkpoint JSON and immutable per-attempt evidence. Freeze the
runner, original physical/numerical profile, mesh/material/doping/IALMob input
hashes, circuit contract, state/evidence hashes, pending targets and counters.
Resume rejects identity mismatches and ignores interrupted/unqualified output.
Runtime/DLL/invalid-input failures stop; numerical failures allow bounded step
reduction. Exit success only when every requested external target is accepted.

## Implementation scope and ownership

One implementer owns new `scripts/run_slot_ldmos_loadline_continuation.py`,
its regression test(s), and a focused CMake test registration. Reuse
`render_attempt`, `render_poisson_prewarm`, `execute_config`,
`audit_local_candidate`, `correct_local_candidate`, and guarded Poly restore.
Do not change existing helper semantics or broad-refactor the old driver.
Resolve original relative inputs including nested IALMob geometry before
relocation; keep per-attempt output paths separate. Root owns this plan,
validation documentation, real runs, and integration/commits.

## Ordered tasks and checks

- [x] Add meaningful failing tests for scalar solving/recovery, signed units,
  original gates, state lineage, evidence binding and resume integrity.
- [x] Implement pure controller and checkpoint layer, then real device adapter
  and CLI. Keep science-specific row expectations explicit in the frozen
  campaign profile; count all active rows, including weak rows.
- [x] Test same-sign/nonfinite/near-zero secant denominators, decreasing inner
  corrections, failed-step midpoint recovery, outer retry queues and budgets.
- [x] Test missing/duplicate/omitted rows, zero-new-bias updates, placeholder
  circuit proofs, failed final audit, wrong state hash/bias, absolute closure
  criterion, runtime errors, interrupted resume and changed input identities.
- [x] Run new tests plus existing continuation/Poly regression suites. Review
  spec compliance first, then code quality; fix findings before acceptance.
- [x] Reproduce the prior external checkpoints, then qualify 37.075/37.1 V
  using the frozen Release runner and original final gates. Preserve failed
  trials. Continue requested higher targets only from qualified checkpoints.
- [x] Record actual runtime/backend and all inner/outer/current/closure
  evidence; commit reviewed implementation locally.

No full C++ build is required solely for this Python orchestration change.
Reconfigure the existing Release preset only as needed for test registration.
All simulation outputs remain ignored under `build/slot_resume_20260922/`.
Original AvalDensGradQF identification and full BVDS remain subsequent stages.
