# Templates/LDMOS Stage-4 D5 sequential continuation

Date: 2026-09-05. Status: **completed through exact Vd=2.66666666666667 V;
80-iteration blocker at 3.369166666666683 V resolved and continuation toward
Vd=4 V restarted with a qualified 160-iteration budget; full curve not qualified**.

This follows section 8 of `templates_ldmos_stage4_handoff_2026-09-05.md`.
The worktree remains `D:/code-repo/vela-tcad/.worktrees/templates-ldmos-phase-a`,
branch `codex/templates-ldmos-phase-a`, with solver HEAD
`345fda30d541003f4aa3a08341dca59418f06ded`.

## Iteration-budget blocker resolved

The 2.66666666666667 -> 4 V interval stopped at approximately 17:11 on
2026-09-05 with 344 accepted attempts and two rejected attempts. It retained
the final accepted state at 3.3679166666666829 V. A 2.5 mV proposal at
3.2141666666666584 V first exhausted 80 Newton iterations; its 1.25 mV fallback
passed. At 3.3691666666666831 V the 1.25 mV proposal also exhausted 80
iterations. The next halving would be 0.625 mV, below the frozen 1 mV minimum,
so step control stopped rather than using the remaining retry allowance.

This was **slow convergence cut off by the iteration budget**, not a proven
floating-point residual floor. The failed electron block decreased from
5.43787e-11 at iteration 70 to 1.03568e-11 at iteration 80, with the leading
electron row at node 3132. The update and leading residual continued to
contract. Poisson and hole residuals met their frozen ceilings.

The recovered checkpoint contains 10,241 unique finite node rows, nonnegative
carrier densities, and matches the rolling state byte for byte. Its SHA-256
is `f7bb2ec10c99e11f07932f42451815724be6509c79afb16bf129cd5f1282e17b`.
The saved drain QF differs from the attempt's exact applied bias by one ULP;
restart commands retain the exact bias recorded in the attempt metadata.

Controlled replay evidence is under
`reference_staging/templates_ldmos_stage4_blocker_3p369167_20260905/`:

| Replay | Result |
| --- | --- |
| Fresh QF recenter, 80 iterations | Fails; final electron block 1.055218834688413e-11 |
| Same replay, budget 160 | Passes at iteration **81**; electron block 8.928110286173847e-12 |
| Trace-prefix comparison | All 83 recorded rows through the 80-iteration run are identical |
| Previously passing 2.6666667 -> 2.6691667 V control, budgets 80/160 | Both converge in 54 iterations; curve, terminal table, iteration trace and state files are byte-identical |

The successful 3.369166666666683 V replay has Poisson block
8.87975385836046e-10, hole block 3.2111993710217456e-25, drain current
1.7265845733373104e-4 A/um, KCL/max-terminal-current 3.711395103693046e-13,
and zero QF violations. Its checkpoint SHA-256 is
`d748971a29980be73d52f260c0c233921ddad8e205d275483316f3a076c72c3a`.
This internal diagnostic point is not added to the 31-point reference score.

The short diagnostic interval needed an explicit `sweep.max_step=0.0025`:
otherwise the runner defaulted max_step to its 1.25 mV reporting span and
rejected the 2.5 mV initial-step setting before solving. The failed preflight
is preserved in `recenter80/`; the actual control is `recenter80_v2/`.
The explicit maximum was used in both diagnostic arms. It is not added to
the resumed full-reference-interval production config.

The reclose CLI now exposes `--max-iter` with its original default of 80.
The resumed sequence explicitly uses **160**. This is the only qualified
numerical-setting amendment: no physical parameter, step-control setting,
mobility Jacobian, QF policy or residual acceptance ceiling was relaxed.

The new live job is under
`reference_staging/templates_ldmos_stage4_sequence_budget160_20260905/`.
Its `status.json` replaces the stopped earlier controller as the current
progress source. It starts from the accepted 3.369166666666683 V replay,
targets exactly 4 V, and retains the previous 0/1.333333/2.666667 V reference
points with provenance. A subsequent 2.5 mV production step has been accepted.
Vg=8 remains gated on Vg=4 qualification, and IALMob is disabled.

## KCL scoring correction and unresolved equilibrium issue

The scorer now uses `abs(sum(I))/max(abs(I))`, matching the frozen validation
plan and handoff. Previously it used `sum(abs(I))` in the denominator, which
could understate a nonzero-bias imbalance by about twofold. No KCL threshold
was changed. It rejects incomplete, duplicate-contact or non-finite terminal
data and reports the zero-bias absolute imbalance separately from the
nonzero-bias maximum. The 27 phase23 regressions, including the added KCL
normalization and malformed-input cases, pass.

For the saved startup zero point, the absolute imbalance is
5.294309596804935e-26 A/um and maximum terminal magnitude is
3.392302652261677e-26 A/um. The corrected raw relative ratio is 156.068315%.
The nonzero startup points' maximum ratio is only 6.390891e-10 percent.
These values are recorded in `kcl_scoring_audit.json` in the blocker evidence
directory. The near-equilibrium ratio remains ill-conditioned; the job does
not invent an absolute resolution floor or silently exempt the zero point
from the gate. The raw gate therefore still prevents a final pass with this
zero-point input, even if future nonzero points pass. This scoring-policy
issue must remain explicit and is distinct from the resolved iteration-budget
blocker. No full-curve or D5 qualification is claimed.

## Second nonzero reference point completed

The recovery interval 1.6483333333333232 -> 2.66666666666667 V completed on
2026-09-05 at approximately 14:34 local time with return code 0, **409 accepted
attempts, 0 rejected attempts, and 409 accepted-step checkpoints**.

| Endpoint metric | Value |
| --- | ---: |
| Exact Vd | 2.66666666666667 V |
| Vela Id | 1.5378988888696833e-4 A/um |
| Sentaurus Id | 1.4179418763901000e-4 A/um |
| Signed current error | +8.4599386% |
| Absolute log current error | 0.03526935 decade |
| Newton iterations | 50 |
| Final Poisson block | 8.70552657587233e-10 |
| Final electron block | 8.540644918018925e-12 |
| Final hole block | 3.4761723414310205e-25 |
| Absolute terminal KCL | 1.3135084054246668e-18 A/um |
| KCL / maximum absolute terminal current | 8.540928242623689e-15 |
| KCL / sum of absolute terminal currents | 4.270464121311863e-15 |
| QF bound violations | 0 |

All three residual blocks satisfy the frozen ceilings. These are pointwise
results; no median/P95 or full-curve accuracy verdict is asserted.

The exact checkpoint is
`reference_staging/templates_ldmos_stage4_sequence_20260905/vg4_interval_02/point_bias_2p666667.csv`,
SHA-256 `bb6556482df88c9b13edfcf7f4fce8ca8eabeffa7af3bcf6a653d95c0149ba21`.
Its 10,241 unique node rows contain only finite values. The same directory's
`completed_interval_audit.json` records the raw endpoint metrics and hashes
for configuration, curve, terminal balance, Newton attempts and iterations.

The next interval, 2.66666666666667 -> 4.0 V, started automatically in
`vg4_interval_03/`. Its same-bias initial reclose converged in one Newton
iteration, and a subsequent accepted internal checkpoint has been saved.
Vg remains 4 V; predictor and IALMob remain disabled. Current live progress
is in the sequence's `status.json` rather than this dated snapshot.

## Interrupted interval and verified recovery

The first 1.33333333333333 -> 2.66666666666667 V run stopped without a
completion summary. No solver process remained when inspected at approximately
12:30 local time. The last checkpoint was written at 12:25:02. The available
attempt records contain no rejected transfer, and there is no rejected-state
bundle. The interruption cause is unknown; it is not evidence of a numerical
convergence failure. Buffered curve/terminal outputs were not finalized, so this
interval contributes no new exact reference point to scoring.

The verified recovery file is:

`reference_staging/templates_ldmos_stage4_vg4_1p333333_to_2p666667_20260905/accepted_step_bias_1p648333.csv`

- SHA-256: `671fd0100e307e7bde4370ad896341c0257db41780a20e7f992f4052d733c104`.
- All 10,241 node rows are present, with the original node IDs in order.
- All numerical values are finite, and carrier densities are nonnegative.
- The file equals the rolling `state.csv` byte for byte.
- Every drain-contact electron QF value equals **1.6483333333333232 V**.
  This also agrees with 126 successive 2.5 mV increments from the starting
  bias. The rounded checkpoint filename was not used as the restart voltage.

The original base-config and 1.333333 V checkpoint hashes were verified against
the handoff. UCRT64 Release configuration/build and the handoff regression set
passed before the first launch: 5/14 material-local cases/assertions, 98/3479
DC-sweep cases/assertions, 110/1438 Newton cases/assertions, and 25 Python
phase23 tests. No solver source or frozen threshold was changed.

## Original 80-iteration sequential execution

The local execution bundle is under the ignored directory:

`reference_staging/templates_ldmos_stage4_sequence_20260905/`

`run_sequence.py` runs one interval at a time through the existing
`run_templates_ldmos_stage4_vd10mv_frozen_mobility_reclose.py` runner. Its first
interval resumes the verified state at 1.6483333333333232 V and targets exactly
2.66666666666667 V. Subsequent intervals use adjacent points from the actual
31-point Sentaurus reference CSV, with a new output directory per interval.

The frozen contract is unchanged: Vg=4 V initially, no predictor, no IALMob,
live HFS residual with lagged field Jacobian, QF recentering, 2.5 mV initial
internal step, 1 mV minimum step, growth factor 1, 12 retries, and original
block residual ceilings. Each interval retains rolling, accepted-step and
exact-point states and the complete Newton/terminal diagnostics.

The launcher checks solver HEAD and executable hash before every interval.
The executable SHA-256 is
`0569876d846792777390a3d134986bf6df8314af85acfaf87a5f0166d12abb36`.
An exclusive `sequence.lock` prevents accidental duplicate launches into the
same evidence directory. Do not rerun the launcher against this directory.

After all Vg=4 exact points exist, the job assembles one curve and terminal
table, then uses the scoring functions and frozen limits from
`scripts/analyze_templates_ldmos_stage4_d5.py`. The zero point comes from the
completed lagged/recenter startup run and 1.33333333333333 V from the completed
post-main run; point provenance is retained. Internal bridge points are not
scored. Missing exact points, unconverged endpoints, residual ceiling violations
and QF bound violations fail closed.

Vg=8 starts only if all Vg=4 **final** single-gate limits pass. It uses the
previously sealed Vg=8, Vd=0 prebias checkpoint and the same D5 contract. Once
both curves exist, the original two-gate analyzer evaluates the current ratio
and full engineering/final gates. IALMob development is not launched by this
job. Any runner failure stops the sequence and records diagnostic paths in
`failure.json`; it does not trigger calibration, threshold relaxation, or
automatic alternate physics.

## Inspecting the original sequence evidence

- `status.json`: current interval, process IDs, most recent checkpoint,
  completed intervals, and overall qualification status; refreshed every 15 s
  while an interval is running. Check timestamp and live processes before
  interpreting a saved `running` status as current.
- `interrupted_run_audit.json`: recovery provenance.
- `vg4_accepted_exact_points.json` / `vg8_accepted_exact_points.json`: current,
  Newton count, residual blocks, KCL and source for each collected reference
  point.
- `vg4_interval_02/`: recovery interval to 2.66666666666667 V.
- `vg4_qualification.json`: generated only after the full Vg=4 reference grid.
- `two_gate_analysis/`: generated only after both curves have completed.
- `controller.stderr.log` and `failure.json`: failure investigation.

The local sequence preflight passed, including reference lattice and recovered
state integrity checks. Synthetic scoring checks confirmed that an exact
reference curve passes and that a missing endpoint, excessive endpoint error,
electron residual ceiling violation, and QF violation are rejected. The
recovery run has produced new accepted checkpoints beyond 1.648333 V.
These startup checks do not qualify the unfinished physical curve.

## Historical scoring audit item found during preparation

The existing zero-drain startup row has source current
`3.3923026522616769e-26`, drain current `1.9020069447711415e-26`, gate current
zero, and substrate current `-2.2788304958386476e-36` A/um. The analyzer's
`abs(sum(I))/sum(abs(I))` normalization therefore produces nearly **100%**
at this numerically tiny equilibrium point. As written, the existing analyzer
includes this row in the KCL maximum, so this input will prevent its KCL gate
from passing even if all future nonzero points pass. No equilibrium-point
exemption or denominator floor has been introduced by this continuation job.
This is a scoring/resolution audit item, not a rejected solver transfer.

Also retain the definition distinction: the handoff's quoted 1.333333 V KCL
value (`5.117277e-14`) uses the maximum absolute terminal current as denominator,
whereas the existing analyzer uses the sum of absolute terminal currents
(`2.558639e-14` at that point). Both are small here; the definitions must not
be silently interchanged near a qualification limit. This job preserves the
requested analyzer behavior and will not promote D5 on an incomplete or
failed result. These issues remain to be resolved explicitly before claiming
final Stage-4 qualification.

## Runtime-log preservation follow-up

The reclose wrapper used to overwrite `reclose.log` with captured stdout after
the runner exited. The runner itself writes its detailed runtime log to that
same path, while stdout can contain only a short summary. This explains why
some completed historical intervals have very short `reclose.log` files.

The wrapper now saves captured stdout to `runner_stdout.log` and preserves any
existing `reclose.log`, using stdout as a fallback only when no runtime log was
created. A focused smoke check verified both branches. This changes logging
only; solver executable, input physics and numerical settings are unchanged.
The interval-02 wrapper was loaded before this edit, so its available
detailed runtime log is additionally retained as
`vg4_interval_02/reclose_runtime_snapshot.log`. That file is a live-run snapshot,
not proof of a complete interval log. Subsequent intervals load the corrected
wrapper.

`audit_completed_interval.py` in the local execution bundle independently
checks completed interval outputs and records exact reference-point current,
Newton count, residual ceilings, both KCL denominators, checkpoint integrity,
and artifact hashes. Replaying the completed 1.0 -> 1.33333333333333 V interval
reproduced the handoff's current error and KCL exactly. It refuses an interval
without a successful completion summary and does not score the unfinished
31-point curve.
