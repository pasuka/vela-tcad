# SLOT-LDMOS BVDS next-chat handoff

Date: 2026-09-05  
Repository: `D:\code-repo\vela-tcad`  
Worktree: `D:\code-repo\vela-tcad\.worktrees\slot-ldmos-bvds-next`  
Branch: `codex/slot-ldmos-bvds-next`  
Implementation baseline before this document: `d5a3f8a`  
Primary solver change: `e19d388`  

## Current qualification reference (2026-09-07)

Read [SLOT-LDMOS BVDS validation](slot_ldmos_bvds_validation.md) and the
[IALMob implementation and controls](slot_ldmos_ialmob_alignment_2026-09-07.md)
before the historical notes below. The old v1 contract selected zero local
carrier rows; its strict global-closure checkpoints do not qualify spatial
continuity under the current all-active profile. Masetti 0--2 V replay and
3/5/10 V checkpoints now pass the stronger local and global gates.

Original-coefficient IALMob has separate qualified no-avalanche 2/3/5/10 V
controls, and self-consistent avalanche 2/3/5 V activation controls. These do
not establish original Sentaurus spatial parity or BVDS. Original
`AvalDensGradQF` also differs from the local SG/GSS-Laux source-current path;
the current report records that distinction and the independent oracle deck.

Numerical work now uses Release. The reported missing `libcholmod.dll` was
an omitted UCRT64 runtime search path and is resolved; use the repository
`scripts/with_ucrt64_path.cmd` launcher. Runtime loader failures are recorded
separately and immediately stop voltage bisection.

## Latest follow-up: strict 1.25 V obtained (2026-09-05)

Local main tip `8e99e57` was merged as `dbe6d3d`. The
[E4/E1/E2 audit](slot_ldmos_closure_restoration_audit_2026-09-05.md) and the
[completed E3 budget and checkpoint report](slot_ldmos_e3_strict_1p25_2026-09-05.md)
supersede the failure status and hypothesis ordering in the historical sections
below. Fixed-vector reduction errors are at most about 8.2e-23 in the audited
baseline trials; actual assembly-order spans are about 1.9e-24. The larger
contact/source versus free-sum discrepancy came from inactive auxiliary gauge
rows, while the dominant electron correction was lost on addition to biased
absolute QF coordinates.

The retained correction combines `quasi_fermi_reference: "contact_basin"` with
`global_continuity_closure.restoration_physical_rows_only: true`. It preserves
the stabilized sensitivity, main Newton equations, physical models, QF limits
and closure tolerances. Reference repacking also subtracts references before
adding saved increments, eliminating the measured restart current drift.

The final strict 1.25 V run executes two real Newton iterations and converges
with `final_global_continuity_restoration`. Drain current is
`6.5889954852964324e-16 A/um`; electron/hole closure mismatches are
`2.4601233527902148e-19` and `2.3598226144186885e-19`. The same-bias reload
retains exactly the reported current and densities. Acceptance uses the existing
carrier-row qualification threshold (zero qualified rows here); QF-bound
diagnostics remain warnings, as documented in the report.

Restart from
`build/slot_ldmos_e3_20260905/final_strict/final_state.csv`, SHA-256
`8de098c8b1416210a0ef05e1e4269104a8f555e48d01ef99e9a308d96818b10a`.
The [frozen restoration contract](contracts/slot_ldmos_physical_closure_restoration_v1.json)
records the options, unchanged gates and evidence hashes. All Newton tests pass
(1610 assertions, 114 cases). Next obtain strict 1.5 V and 2 V checkpoints with
this contract before resuming the high-voltage sequence. Vela BVDS is still
unavailable. Code and reports after the merge remain uncommitted worktree changes.

## 1. Historical handoff outcome before follow-up

The Sentaurus reference is complete and reproducible. Vela has strict,
avalanche-off checkpoints through 1.0 V, but it does not yet have a credible
SLOT-LDMOS BVDS curve. The next task is not avalanche or IALMob development.
It is a deterministic audit of the 1.25 V final global-continuity restoration.

The latest 1.25 V same-bias restart reduces the scalar normalized residual to
`1.2948077136e-9`, but the carrier block norms remain approximately
`1.47536e-8` and `1.02587e-8`. The electron/hole global closure mismatches are
`1.84654e-12` and `3.42869e-12`, above the low-current absolute limit
`1e-12`. No 1.25 V point has been accepted as a physical checkpoint.

Current physics for this debug path:

| Item | Setting |
|---|---|
| Mesh | Delaunay assembled fallback plus mixed-Voronoi volumes |
| Mobility | `masetti_field` |
| IALMob / Lombardi | Off |
| Impact ionization | Off |
| External resistor | Off; direct drain-voltage continuation |
| Recombination | SRH(DopingDep) plus Auger |
| Band-gap narrowing | Old Slotboom |
| High-field drive | Quasi-Fermi gradient, transport-cell vector |
| QF update | `componentwise_poisson_recorrect` |
| QF limits | 0.1 V majority, 0.05 V minority |
| Carrier diagonal floor | Enabled, scale 10, minority ratio 1 |
| Newton globalization | Merit line search |
| Closure | Final gate; absolute below `1e-11`, relative above it |
| Absolute/relative limits | `1e-12` and 1% |
| Final restoration | `final_linear_qf_projection` |
| Stall floor | `1e-8` |

## 2. Sentaurus oracle

The original project is:

`D:\gpt-tmp\slot-LDMOS\V3E\_NLDMOS\_40V\_Copy`

Sentaurus release: T-2022.03-SP2.

| Stage | Outer drain V | Inner drain V | Drain current A/um |
|---|---:|---:|---:|
| Equilibrium | 0 | 0 | `-5.9567e-24` |
| 1 Tohm-um, 1 V calibration | 1 | `0.999409604` | `5.90255e-16` |
| Direct 1 V | 1 | 1 | `5.90523e-16` |
| Avalanche off | 60 | `59.9934120` | `6.58796e-15` |
| IIC postprocessing | 60 | `59.9934120` | `6.58796e-15` |
| Self-consistent avalanche | 60 | `38.5048827` | `2.14951e-11` |
| Final external resistor | `109076.963` | `38.5215816` | `1.09038e-7` |

Interpolation at `|Id|=1e-7 A/um` gives the authoritative internal BVDS:

`BVDS = 38.5209012 V`.

The approximately 100 kV outer voltage is the voltage drop required by the
`1e12 ohm*um` series resistor. It is not the device BVDS. The 2-D unit
calibration establishes that Vela must use `resistance_ohm_um=1e12` without a
second width factor.

Reference document:

`docs/validation/slot_ldmos_bvds_sentaurus_reference_2026-08-20.md`

## 3. Valid Vela results

| Checkpoint | Id A/um | Electron mismatch | Hole mismatch | Status |
|---|---:|---:|---:|---|
| Strict 0.5 V | `3.0535948881e-16` | `5.29563e-13` | `3.58819e-13` | Accepted by final restoration |
| Strict 0.625 V | `3.7207300754e-16` | `3.20996e-13` | `5.11552e-13` | 8 Newton steps plus final restoration |
| Strict 0.8125 V | `4.67737102146e-16` | `9.62699e-13` | `7.41985e-13` | Accepted intermediate point |
| Strict 0.90625 V | `5.63087057643e-16` | `3.96311e-14` | `4.25603e-13` | Accepted intermediate point |
| Strict 1.0 V | `5.1772757291e-16` | `8.98421e-13` | `7.52880e-13` | Accepted |

The strict 1 V current differs from the direct-voltage Sentaurus current by
`-0.0571355223 decade`, inside the 0.3-decade pre-breakdown gate.

The old 0.5 V provisional state with current `2.9187277463e-16 A/um` is not a
strict checkpoint. It was accepted with zero Newton iterations through
`initial_stall_residual_floor`. The strict 0.5 V point supersedes it.

The old legacy-profile high-voltage Vela curves are invalidated and must not be
used to quote Vela BVDS.

## 4. Work already completed

1. Generated the final SProcess TDR on the Sentaurus VM.
2. Imported topology, coordinates, node doping, regions, materials and five
   contacts into Vela.
3. Resolved SLOT and substrate contact semantics.
4. Calibrated 2-D current and resistance units.
5. Transcribed Old Slotboom, SRH(DopingDep), Auger and Van Overstraeten data.
6. Audited obtuse versus non-Delaunay elements. The mesh contains many obtuse
   triangles but Sentaurus reports only two non-Delaunay cells. Vela now has a
   Delaunay/mixed-Voronoi path with local fallback.
7. Implemented the augmented load-line Newton infrastructure.
8. Added checkpoint persistence and explicit last-accepted failure state
   handling.
9. Added the two-regime global-continuity closure criterion.
10. Corrected ordinary globalization so `report/enforce` closure is a final
    gate rather than an intermediate merit constraint.
11. Froze closure normalization at bias-point scope instead of recomputing a
    moving denominator on every iteration.
12. Added `final_linear_qf_projection` and validated strict checkpoints through
    1 V.

The implementation checkpoint is described in:

- `docs/validation/slot_ldmos_final_only_globalization_progress_2026-08-30.md`
- `docs/validation/contracts/slot_ldmos_globalization_baseline_v1.json`

## 5. Latest failed runs

| Attempt | Scalar residual | Psi block | Electron block | Hole block | Electron mismatch | Hole mismatch |
|---|---:|---:|---:|---:|---:|---:|
| 1.0 to 2.0 V, 10 steps | `2.2946e-1` | `2.4592e-1` | `8.44183e2` | `1.572e-2` | `-3.888e3` | `3.347e-8` |
| 1.0 to 1.5 V, 10 steps | `1.1346e-3` | `1.2913e-3` | `2.80068` | `2.579e-8` | `12.6239` | `4.447e-12` |
| 1.0 to 1.25 V, 10 steps | `1.2387e-8` | `1.4736e-8` | `1.7614e-8` | `7.912e-9` | `1.358e-12` | `4.452e-12` |
| 1.25 V resume 1 | `4.0325e-9` | `4.7971e-9` | `1.4811e-8` | `9.487e-9` | `1.203e-12` | `4.313e-12` |
| 1.25 V resume 2 | `1.2948e-9` | `1.5403e-9` | `1.4754e-8` | `1.0259e-8` | `1.847e-12` | `3.429e-12` |

Latest directory:

`build/slot_ldmos_noavalanche_1p25_resume02/outputs/stages/01_delaunay_noavalanche_adaptive/attempt_001_v_1p25/`

Important local files:

- `iv.csv`
- `iterations.csv`
- `local_updates.csv`
- `iv_newton_failure_diagnostics.json`
- `rejected_states/attempt_1_bias_1p250000_final.csv`

These are ignored local validation outputs and are intentionally not committed.
In a failed `iv.csv` row, terminal current fields are zero because no accepted
terminal state was reported; that is not a physical zero-current result.

## 6. Code-level findings that guide the next audit

### F1. The halving trial schedule is incompatible with the linear closure gate

Final restoration starts with
`feasibilityRestorationFactor=1` and halves alpha for every retry. For a
linear equality correction, the remaining mismatch is approximately
`(1-alpha) * g0`.

At the latest 1.25 V state:

- Electron requires
  `alpha >= 1 - 1e-12 / 1.84654e-12 = 0.4585`.
- Hole requires
  `alpha >= 1 - 1e-12 / 3.42869e-12 = 0.7083`.

Therefore the current sequence `1, 0.5, 0.25, ...` has only one trial that can
simultaneously satisfy both carrier gates in the ideal linear model: alpha 1.
The remaining trials do not constitute a meaningful nonlinear search for this
state.

Relevant code begins near `src/solver/NewtonSolver.cpp:7323`.

### F2. Restoration sensitivity is constructed from the stabilized Jacobian

The final restoration reassembles the same Jacobian used by Newton, applies
carrier row regularization, and then constructs the two global rows by summing
free carrier rows by column. With `carrier_diagonal_floor scale=10`, weak
deep-depletion rows receive artificial diagonal sensitivity. This can make the
minimum-norm QF correction concentrate on degrees of freedom with little
physical terminal-current observability.

This gives a plausible complete failure chain:

1. Stabilized global sensitivity predicts an effective closure correction.
2. The physical nonlinear residual does not change by the predicted amount.
3. The correction is larger than necessary or concentrated in depleted QFs.
4. Componentwise QF clipping changes the correction again.
5. Alpha 1 breaks the device-floor gate, while the halved trials cannot meet
   the hole closure gate even in the linear model.

Relevant code:

- `src/solver/NewtonSolver.cpp:608`: carrier row regularization.
- `src/solver/NewtonSolver.cpp:1341`: global row construction.
- `src/solver/NewtonSolver.cpp:7336`: final-restoration Jacobian path.

### F3. The absolute gate is already close at 0.8125 and 1 V

The strict electron closure mismatch is `9.62699e-13` at 0.8125 V and
`8.98421e-13` at 1 V, compared with an absolute limit of `1e-12`.

This is evidence that the available margin is small, but it does not by itself
prove that the gate is impossible. The sum of already assembled double
residuals and the error introduced while assembling those residuals are two
different floors and must be measured separately before changing the contract.

### F4. `replace_free_row_sum` replaces one physical row per carrier

The stabilized linear solve replaces a selected free electron row and a
selected free hole row with the corresponding global row-sum equations. The
ordinary merit does not include closure, so subsequent line-search damping may
apply only a small fraction of the linearized closure correction. The replaced
physical row is also not directly solved in that linear system.

The selected replacement row IDs should be compared with the observed hot
residual nodes, including electron nodes 6248 and 5849 and hole node 8063.

Relevant code begins near `src/solver/NewtonSolver.cpp:1118` and its use near
`src/solver/NewtonSolver.cpp:6749`.

### F5. The current Gram degeneracy threshold is too weak as a diagnostic

After normalization, the 2x2 determinant corresponds to a squared sine of the
angle between the two rows. A threshold close to machine epsilon rejects only
almost exactly parallel rows. A poorly conditioned pair can proceed with large
multipliers without an actionable diagnostic. Record singular values,
condition number and multiplier magnitudes; do not infer good conditioning
from the current determinant guard.

## 7. Hypothesis ordering

| Rank | Hypothesis | Current status |
|---:|---|---|
| 1 | Stabilized Jacobian contaminates physical closure sensitivity | Directly supported by code path; needs A/B replay |
| 2 | Device/block stopping floors and `1e-12` global cancellation gate have inconsistent attainable accuracy | Strong evidence; assembly floor not yet measured |
| 3 | Alpha-halving schedule has no feasible trial after alpha 1 | Demonstrated for the latest e/h mismatches |
| 4 | Unweighted minimum-norm projection over weakly observable QFs is ill-conditioned | Plausible; Gram/SVD diagnostics absent |
| 5 | Project-then-cap breaks closure | Code order confirms it; may be secondary to an oversized polluted correction |
| 6 | Replaced physical row causes or anchors a carrier-block hotspot | Cheap to test; not yet confirmed |
| 7 | Absolute/relative mode transition | Low priority; both carriers currently use the absolute mode |
| 8 | Adaptive-driver strategy | Real workflow issue but not the root of final closure failure |

Do not reopen the already closed investigations into contact/QF bias lifting,
inner versus outer voltage, warm-start reference coordinates, predictor
coordinates or the obtuse-triangle policy unless new evidence contradicts the
existing audits.

## 8. Immediate development plan

### E4 first: replacement-row audit

This has the lowest implementation cost. Record for both carriers:

- Replaced equation index.
- Carrier and node ID.
- Coordinates and density regime.
- Physical residual before and after every accepted iteration.
- Whether the row appears in the top residual list.

If a replacement row is node 6248, 5849, 8063 or a persistent neighboring
hotspot, retain the evidence. Do not immediately redesign the solver; complete
E1 and E2 first.

### E1: trial-level final-restoration diagnostics

For every restoration trigger and candidate alpha, record:

- Initial physical electron/hole mismatch and criterion mode.
- Predicted mismatch from `g + A*delta`.
- Reassembled nonlinear physical mismatch.
- Scalar residual and psi/electron/hole block norms.
- Carrier-row qualification.
- Raw, projected, capped and applied correction L2/Linf norms.
- Number of capped QF DOFs, maximum over-limit value, carrier and node IDs.
- Gram entries, row norms, normalized row angle or `sin(theta)`, singular
  values, condition number and multiplier values.
- Poisson-recorrect before/after values.
- Explicit rejection-reason enum.

Minimum tests:

1. Independent rows, no cap: exact linear restoration.
2. Nearly dependent rows: conditioning is reported.
3. Active cap: the implementation does not claim the original constraint is
   still satisfied.
4. Regularized and physical Jacobian modes remain distinguishable.
5. Alpha values below the ideal closure-feasible lower bound are labelled and
   skipped or reported as infeasible.

### E2: physical-Jacobian A/B

Keep the main Newton trajectory unchanged. Add an audit/restoration-only option
that constructs the closure sensitivity from an unregularized physical
Jacobian or independent physical JVP. Compare, from the same saved state:

- Predicted versus actual closure change.
- Correction norms.
- QF cap count.
- Candidate device/block residuals.
- Final qualification.

This experiment must not disable the diagonal floor for the main Newton solve.
Only the restoration sensitivity changes.

### E3: two-level noise-floor measurement

Do not estimate the complete floor solely by reordering the final global sum.
Separate:

1. Global reduction floor: fixed node residuals summed in ordinary double,
   pairwise, Kahan/Neumaier, long double, ascending/descending magnitude and
   deterministic random order.
2. Assembly floor: reassemble the same state with reversed/permuted edge
   scatter order and compare ordinary versus compensated per-node
   accumulation.
3. Repeatability: assemble the same state more than once and check bitwise and
   numerical differences.

Only after measuring both floors may the `1e-12` absolute contract be changed.
If changed, derive the new limit from measured floor plus a documented safety
factor and requalify every 0.5 to 1.25 V checkpoint consistently.

## 9. Expected minimal fix order

Apply one change per A/B experiment:

1. Use the unregularized physical Jacobian/JVP for closure sensitivity while
   retaining stabilized Newton.
2. Replace the current determinant guard with a two-row SVD/QR diagnostic and
   robust solve.
3. Weight QF corrections by physical closure observability and optionally
   restrict them to a drain/current-path basin.
4. Replace project-then-cap with a cap-aware active set: solve, freeze capped
   DOFs at their bounds, and resolve the two-row system on the remaining DOFs.
5. Replace blind halving with a closure-feasible interval. For the latest
   state, search only in approximately `[0.7083, 1]`, with a numerical margin.
6. If replacement-row auditing proves that physical rows are being stranded,
   consider a bordered/KKT constraint solve that retains all physical rows.

PTC is not the first choice at 1.25 V: the ordinary Newton path is still
reducing its scalar residual. The present blocker is final qualification and
restoration geometry, not a failure to find a descending nonlinear path.

## 10. Acceptance gates before high-voltage work

Obtain strict 1.25, 1.5 and 2 V checkpoints. Each must:

- Execute real Newton iterations.
- Pass carrier-row convergence.
- Pass the same documented two-regime closure rule.
- Save a restartable accepted final state.
- Report a non-ambiguous convergence reason.
- Preserve all physical settings used for the strict 1 V point.

For `max_iterations` failures, allow another same-bias chunk while residual
reduction is material and line search remains productive. Return to voltage
bisection only after a documented plateau, repeated severe damping or a
nonfinite/physical failure. A restartable failed state must never be labelled
as a converged checkpoint.

## 11. Work after 2 V

1. Continue avalanche-off points at 5, 10, 20, 40 and 60 V.
2. At 60 V compare with the Sentaurus current `6.58796e-15 A/um`; require an
   error no larger than 0.3 decade.
3. Perform frozen-state/IIC avalanche-source mapping.
4. Use avalanche coupling homotopy.
5. Enable the already developed augmented series-resistor/load-line path and
   cross the folded branch.
6. Reach `|Id| >= 1e-7 A/um` and interpolate internal BVDS versus log current.
7. Only then perform strict IALMob off/on A/B. Keep every setting bytewise
   identical except `masetti_field` versus `masetti_field_lombardi`.

IALMob must not be enabled during the immediate E1-E4 debug work.

## 12. Files to read first in the next chat

1. This document.
2. `docs/validation/slot_ldmos_final_only_globalization_progress_2026-08-30.md`
3. `docs/validation/contracts/slot_ldmos_globalization_baseline_v1.json`
4. `docs/validation/slot_ldmos_bvds_sentaurus_reference_2026-08-20.md`
5. `src/solver/NewtonSolver.cpp`
6. `include/vela/solver/NewtonSolver.h`
7. `tests/test_newton_solver.cpp`
8. `scripts/run_slot_ldmos_noavalanche_continuation.py`
9. Latest local `iv.csv`, `iterations.csv`, `local_updates.csv` and
   `iv_newton_failure_diagnostics.json` listed above.

Other useful code locations:

- `src/solver/NewtonSolver.cpp:608`: carrier-row regularization.
- `src/solver/NewtonSolver.cpp:1118`: replacement-row selection.
- `src/solver/NewtonSolver.cpp:1341`: two global sensitivity rows and Gram
  solve.
- `src/solver/NewtonSolver.cpp:6312`: closure normalization and merit setup.
- `src/solver/NewtonSolver.cpp:6749`: free-row-sum replacement.
- `src/solver/NewtonSolver.cpp:7323`: final-only restoration.

## 13. Verification state

At the implementation baseline `d5a3f8a`, the targeted Newton test executable
passes:

`All tests passed (1409 assertions in 103 test cases)`

No generated simulation output should be added to Git. Keep new diagnostics in
`build/` unless a small reviewed contract or Markdown summary is intentionally
promoted.

## 14. Suggested opening prompt for the next chat

> Work in `D:\code-repo\vela-tcad\.worktrees\slot-ldmos-bvds-next` on branch
> `codex/slot-ldmos-bvds-next`. Read
> `docs/validation/slot_ldmos_bvds_next_chat_handoff_2026-09-05.md` completely,
> then inspect the referenced implementation and latest 1.25 V outputs. Begin
> with E4 and E1: add replacement-row and trial-level final-restoration
> diagnostics plus Catch2 coverage. Do not change physics, closure tolerances,
> QF limits, carrier diagonal floor, avalanche, external resistance or IALMob
> during the diagnostic replay. Rebuild with MSYS2 UCRT64, run the targeted
> Newton tests, and replay the exact 1.25 V last-accepted state. Report whether
> predicted versus actual closure changes confirm stabilized-Jacobian
> contamination before implementing a numerical fix.
