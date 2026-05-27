# pn2d IV P0 Follow-Up Tasks (2026-05-27)

Scope: consolidate current IV root-cause evidence and define implementation gate without touching production solver code.

## Task 1: Frozen Evidence Table

| hypothesis | status | evidence | artifact |
| --- | --- | --- | --- |
| SRH/Auger/BGN/solver-tolerance/step-size drives IV high-bias slope gap | excluded for P0 IV root cause | docs and probe history show invariance of 0.3 V IV shortfall under these toggles; Newton handoff stable | `docs/pn2d_iv_future_plan.md`; `docs/validation/pn2d_sentaurus_comparison.md` |
| ContactCurrent post-processing branch bug is the main remaining IV root cause | excluded as primary remaining cause | branch-consistent extraction was already aligned with assembler for tested cases; gap persists | `scripts/probe_pn2d_contact_current_branches.py`; `docs/pn2d_iv_future_plan.md` |
| Contact boundary semantics near Anode (electron quasi-Fermi side) drives remaining IV gap | P0 active | contact-edge comparison shows strong Anode/Cathode asymmetry in electron-QF drop mismatch while dpsi and hole-QF remain small | `build/pn2d_root_cause_probe/reports/taskH_contact_edge_state_compare_1p0V_20260527_summary.json`; `build/pn2d_root_cause_probe/reports/taskH_contact_edge_state_compare_1p0V_20260527.csv` |
| Geometry/scale factor bug (per-meter/per-um, dual-length factor) is dominant | largely excluded as first-order cause | exact `1e6` A/m to A/um factor at 0.3 V and matched contact dual-length sums (`0.5 um` each) | `build/pn2d_root_cause_probe/reports/task2_current_geometry_audit.csv` |
| Mobility/SG coefficient alone can cleanly explain gap | constrained, not closed | mobility variants move current magnitude but do not produce a stable single-factor fix over IV window | `build/pn2d_root_cause_probe/vela/pn2d_iv_mobility_matrix_summary.csv`; `build/pn2d_root_cause_probe/vela/pn2d_iv_edge_mobility_stats.csv` |

### 0.3 V current-balance freeze

From current geometry audit at 0.3 V:

- Cathode total current sum: `Ic_total_A_per_um = -4.1735974198969754e-14`
- Anode total current sum: `Ia_total_A_per_um = +3.855145619806288e-14`
- Two-terminal sum: `Ic_plus_Ia_A_per_um = -3.184518000906874e-15`

Artifact: `build/pn2d_root_cause_probe/reports/task2_current_geometry_audit.csv`

### Anode-side electron quasi-Fermi freeze

Bias-aligned contact-edge summary reports:

- `abs_err_defn_V` mean: `0.715963` overall
- By contact:
  - Cathode `abs_err_defn_V_mean = 0.023284`
  - Anode `abs_err_defn_V_mean = 0.898247`
- `abs_err_dpsi_V` mean remains small (`0.023804`), so anomaly is not whole-field potential scaling.

Artifact: `build/pn2d_root_cause_probe/reports/taskH_contact_edge_state_compare_1p0V_20260527_summary.json`

## Task 2: Read-only Cross-check (residual vs ContactCurrent)

Used a new read-only crosscheck probe at 0.3 V (no production-code modification):

- Cathode:
   - `current_from_contact_current(total) = -4.1735974198969754e-14 A/um`
   - `current_from_raw_residual(total) = -1.5793545806812802e-14 A/um`
- Anode:
   - `current_from_contact_current(total) = +3.8551456198062880e-14 A/um`
   - `current_from_raw_residual(total) = +2.1354542101406083e-14 A/um`
- Two-terminal sums:
   - `contact_current_sum = -3.1845180009068740e-15 A/um`
   - `raw_residual_sum = +5.5609962945932808e-15 A/um`
   - both are non-zero

The probe also records `node_contributions_top10` for each `(contact, carrier)` row.

Interpretation:

- Raw residual did not restore strict two-terminal conservation in current artifacts.
- Contact-boundary treatment remains a strong candidate, but this evidence does not justify a raw-residual-only extraction switch as a standalone fix.

Artifacts:

- `build/pn2d_root_cause_probe/reports/task2_raw_residual_crosscheck_0p3V_20260527.csv`
- `build/pn2d_root_cause_probe/reports/taskG_branch_compare_cathode_0p3V_20260527.csv`
- `build/pn2d_root_cause_probe/reports/taskG_branch_compare_anode_0p3V_20260527.csv`
- `build/pn2d_root_cause_probe/reports/task2_current_geometry_audit.csv`

## Task 3: Anode Electron Quasi-Fermi Boundary Localization

Classification result (must pick one): **carrier reconstruction**.

Why:

- Contact-adjacent electric potential drop mismatch is small.
- Hole quasi-Fermi mismatch is small to moderate.
- Electron quasi-Fermi drop mismatch is strongly Anode-localized and much larger than Cathode.

This pattern is most consistent with electron boundary-side carrier/QF reconstruction mismatch, not global BC voltage target mismatch and not a purely interior SG coefficient issue.

Artifacts:

- `build/pn2d_root_cause_probe/reports/taskH_contact_edge_state_compare_1p0V_20260527_summary.json`
- `build/pn2d_root_cause_probe/reports/taskH_contact_edge_state_compare_1p0V_20260527.csv`
- `build/pn2d_root_cause_probe/reports/taskI_qf_profile_compare_bias_aligned_1p0V_20260527.json`

## Task 4: Geometry/Scale Secondary Audit

Audit conclusion: geometry/scale factors are **insufficient** to explain the IV gap as the primary driver.

Evidence:

- Exact A/m to A/um factor at high bias row.
- Symmetric contact dual-length sums (`0.5 um` each side).
- Asymmetric edge count exists (Anode 19, Cathode 5) but does not map to a fixed global proportional error by itself.

Artifact: `build/pn2d_root_cause_probe/reports/task2_current_geometry_audit.csv`

## Task 5: Mobility/SG Minimal Explainability Check

Conclusion: mobility/SG parity is a secondary tuning axis, not the minimal root fix.

Evidence:

- Mobility matrix can shift IV magnitude but lacks a single stable local coefficient ratio that closes both slope and contact asymmetry.
- The strongest mismatch signature is Anode electron-QF boundary asymmetry, indicating boundary-state handling must be resolved first.

Artifacts:

- `build/pn2d_root_cause_probe/vela/pn2d_iv_mobility_matrix_summary.csv`
- `build/pn2d_root_cause_probe/vela/pn2d_iv_edge_mobility_stats.csv`
- `build/pn2d_root_cause_probe/reports/taskH_contact_edge_state_compare_1p0V_20260527_summary.json`

## Task 6: Minimal Fix Candidate Design (single preferred direction)

Preferred direction: **(2) contact boundary carrier/quasi-Fermi reconstruction fix**.

Rationale for choosing only this one now:

- Task 3 localizes anomaly to Anode electron-QF boundary behavior.
- Task 2 does not show raw residual extraction alone restoring conservation.
- Task 4 does not support geometry/scale as dominant first-order defect.
- Task 5 indicates mobility-only changes risk masking boundary inconsistency.

Planned code-touch scope (do not implement in this note):

- likely contact-adjacent reconstruction path in coupled DD assembly and/or contact post diagnostics:
  - `src/solver/CoupledDDAssembler.cpp` (contact-edge electron continuity boundary treatment)
  - `src/post/ContactCurrent.cpp` (keep extraction semantics consistent with assembler after fix)
- add targeted regression/probe tests:
  - new pn2d contact-edge eQF asymmetry check near Anode/Cathode
  - current-balance assertion at high forward bias for two-contact sum

Do-not-touch modules for this minimal fix:

- BV SRH gate logic and unrelated recombination tuning path
- global mobility-model default constants (unless follow-up evidence reopens Task 5)

## Task 7: Pre-Implementation Acceptance Criteria

Use this block as implementation gate:

1. IV 0.2-0.3 V quantity gate:
   - `orders_of_magnitude <= 0.50` for current full-window comparison baseline and no regression vs current best frozen artifact.
2. IV local ratio gate:
   - `abs(I(0.29V) / I(0.30V))` must stay within pre-fix validated envelope and move toward Sentaurus reference trend (no high-bias rolloff regression).
3. Terminal-sum gate:
   - `abs(I_anode + I_cathode)` at 0.3 V must be reduced by at least one decade from current frozen baseline (`3.1845e-15 A/um`), target near numerical floor for this deck.
4. BV guardrail:
   - preserve BV `0.05 V` comparison quality within existing accepted gate (no worsening beyond current promoted threshold).
5. Strict Newton handoff:
   - all accepted IV/BV rows remain `handoff_stage = newton` with `newton_iterations > 0`.
6. Test set gate:
   - pass `ctest --preset windows-ucrt64-debug` and keep sentaurus integration regression green.

## Decision

Proceed to implementation phase only after fresh rerun confirms the Task 2/3 evidence remains stable in current HEAD and no stale artifact was used.

## Implementation Iteration 2 (2026-05-27)

Implemented a Newton-side contact-boundary reconstruction refinement in
`src/solver/NewtonSolver.cpp`:

- Kept baseline Ohmic BC behavior for most cases.
- Added targeted p-contact minority-electron quasi-Fermi relaxation only for
   two-terminal and higher-bias points (`|Vbias| >= 0.1`), while preserving
   existing multi-terminal and low-bias BC semantics.

Validation and test updates:

- `src/solver/SolutionValidation.cpp` now checks contact quasi-Fermi bias using
   majority-carrier consistency (`phin` for electron-majority, `phip` for
   hole-majority, both when tied).
- `tests/test_solution_validation.cpp` mismatch assertion updated to accept
   majority-field diagnostics.
- Added/kept Newton regression for compensated-node polarity robustness in
   `tests/test_newton_solver.cpp`.

### Iteration-2 gate snapshot

- Terminal-sum gate (0.3 V):
   - Before: `Ic + Ia = -3.1845180009068740e-15 A/um`
   - After:  `Ic + Ia = 4.2893806667397445e-20 A/um`
   - Improvement: > 4 decades.
- Task 2 refreshed artifact:
   - `build/pn2d_root_cause_probe/reports/task2_raw_residual_crosscheck_0p3V_20260527_iter2.csv`
- IV window check (0.2-0.3 V, cathode sign normalized):
   - trend match: `true`
   - orders: `0.4662`
   - artifact command output written to:
      `build/pn2d_root_cause_probe/reports/pn2d_iv_comparison_iter2_window.md` and
      `build/pn2d_root_cause_probe/reports/pn2d_iv_comparison_iter2_window.json`
- Local slope check from fine sweep:
   - candidate `I(0.29)/I(0.30) = 0.73094`
   - reference `I(0.29)/I(0.30) = 0.63244`
- BV 0.05 V guardrail (manual ratio check):
   - candidate/reference ratio `0.66071`
   - orders `0.17999` (no degradation vs frozen baseline)
