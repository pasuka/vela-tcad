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

## Two-Week Execution Plan (M1-M3)

### M1: Freeze regression entry and explicit strategy boundary

- Add pn2d high-bias regression gate in Sentaurus integration coverage with:
   - IV `0.2-0.3 V` trend match.
   - IV window `orders_of_magnitude <= 0.50`.
   - Local slope gate at `I(0.29)/I(0.30)` relative to the frozen baseline delta
      versus Sentaurus (`|candidate-ref| <= |0.7309-0.6324|`).
   - Two-terminal sum gate at `0.3 V`:
      `abs(I_anode + I_cathode) <= 1e-18 A/um`.
- Expose an explicit Newton solver config key:
   `contact_boundary_reconstruction`, defaulting to the currently accepted
   behavior (`dominant_signed_contact_mean`).
- Document that this key controls only contact-boundary quasi-Fermi/carrier
   reconstruction and is not a mobility-tuning parameter.

Acceptance:
- New/updated tests pass.
- pn2d frozen metrics do not regress.
- Full preset test gate remains green (`ctest --preset windows-ucrt64-debug`).

### M2: High-bias local slope refinement

- Restrict tuning to anode-adjacent minority-electron quasi-Fermi/carrier
   reconstruction (A/B only).
- Candidate matrix:
   - current strategy;
   - more conservative relaxation;
   - stronger relaxation;
   - bias-threshold micro-adjustments.
- For each candidate, record:
   - `I(0.29)/I(0.30)`;
   - IV window orders;
   - terminal sum at `0.3 V`;
   - BV `0.05 V` orders;
   - strict Newton handoff status.

Target:
- Reduce local slope absolute gap to Sentaurus from about `0.0985` to
   `<= 0.075`.
- Keep IV window orders `<= 0.4662` (target `< 0.45`).
- Keep terminal sum `<= 1e-18 A/um`.
- No BV guardrail regression.

Rollback:
- Keep iter2 as baseline;
- reject candidates that improve one-point slope but degrade window/BV gates.

### M2 execution result (2026-05-27)

Implemented M2 as explicit Newton knobs plus automated candidate scan:

- New Newton config knobs (default behavior preserved):
   - `contact_boundary_minority_electron_relaxation`
   - `contact_boundary_minority_electron_relaxation_bias_threshold_V`
   - `contact_boundary_minority_electron_relaxation_two_terminal_only`
   - `contact_boundary_minority_electron_relaxation_contact_side`
- Added scan automation:
   - `scripts/scan_pn2d_contact_relax_candidates.py`
   - Output root: `build/pn2d_contact_relax_scan`
   - Summary artifacts:
      - `build/pn2d_contact_relax_scan/pn2d_contact_relax_summary.csv`
      - `build/pn2d_contact_relax_scan/pn2d_contact_relax_summary.json`

Candidate matrix evaluated (strategy x contact-side):

- `baseline` (iter2 current behavior)
- `dominant_p_only`
- `dominant_n_only`
- `dominant_both`
- `legacy_p_only`
- `legacy_n_only`
- `legacy_both`

Observed metric snapshot (from summary CSV):

- Local slope ratio `I(0.29)/I(0.30)`:
   - baseline / `*_p_only` / `*_both`: `0.7309387` (`delta ~= 0.098503`)
   - `dominant_n_only` / `legacy_n_only`: `0.7025914` (`delta ~= 0.070156`)
- IV window orders (`0.2-0.3 V`):
   - baseline / `*_p_only` / `*_both`: `0.4662111`
   - `dominant_n_only` / `legacy_n_only`: `0.4945596`
- Terminal sum at 0.3 V:
   - all candidates remain near numerical floor (`8.5e-20` to `1.7e-19 A/um`)
- BV 0.05 V orders:
   - all candidates `0.1109109`
- Strict Newton handoff:
   - all candidates `true` for IV/BV/fine sweeps

M2 decision:

- `n_contact_only` branches improve local slope gap (`0.0985 -> 0.0702`) but
  regress IV window orders (`0.4662 -> 0.4946`), violating the current gate.
- Reconstruction mode (`dominant` vs `legacy`) does not materially change this
  outcome in the tested matrix.
- Keep iter2 promoted baseline (`dominant_p_only` equivalent behavior) and mark
  `n_contact_only` as a promising but currently unqualified direction.

### M2 round2 result: n-only threshold refinement (2026-05-27)

A focused threshold micro-sweep was executed for the `n_contact_only` branch
using:

- `scripts/scan_pn2d_n_only_thresholds.py`
- Summary artifacts:
   - `build/pn2d_contact_relax_scan/pn2d_contact_relax_round2_n_only_summary.csv`
   - `build/pn2d_contact_relax_scan/pn2d_contact_relax_round2_n_only_summary.json`

Round2 candidate set:

- `baseline`
- `n_only_th0p08`
- `n_only_th0p1`
- `n_only_th0p12`
- `n_only_th0p15`
- `n_only_th0p2`

Observed outcome (stable across all tested thresholds):

- All `n_only_th*` points are numerically identical in this deck:
   - `I(0.29)/I(0.30) = 0.7025914`
   - `delta_to_reference = 0.0701557`
   - `iv_window_orders(0.2-0.3V) = 0.4945596`
   - `terminal_sum_abs_A_per_um_at_0p3 = 1.7293e-19`
   - `bv_orders_at_0p05 = 0.1109109`
   - `strict_newton_handoff_all = true`
- Relative to baseline:
   - local slope improves (`delta: 0.0985 -> 0.0702`),
   - but IV-window orders remain worse (`0.4662 -> 0.4946`).

Round2 decision:

- The threshold axis (`0.08-0.20 V`) does not provide additional separation for
  the current `n_contact_only` branch.
- No round2 threshold candidate is promotable under the existing IV-window gate.
- Keep iter2 promoted baseline unchanged; treat `n_contact_only` as a known
  trade-off branch pending a new non-threshold mechanism.

### M2 round3 result: n-only edge-threshold refinement (2026-05-27)

To close the threshold-axis uncertainty near the comparison window, the same
scan script was rerun with high thresholds:

- Command:
    - `python scripts/scan_pn2d_n_only_thresholds.py --thresholds 0.24,0.26,0.28,0.29,0.295 --summary-prefix pn2d_contact_relax_round3_n_only_edge_summary --candidate-dirname candidates_round3`
- Summary artifacts:
    - `build/pn2d_contact_relax_scan/pn2d_contact_relax_round3_n_only_edge_summary.csv`
    - `build/pn2d_contact_relax_scan/pn2d_contact_relax_round3_n_only_edge_summary.json`

Round3 candidate set:

- `baseline`
- `n_only_th0p24`
- `n_only_th0p26`
- `n_only_th0p28`
- `n_only_th0p29`
- `n_only_th0p295`

Observed outcome (again identical for all `n_only_th*` points):

- `I(0.29)/I(0.30) = 0.7025914`
- `delta_to_reference = 0.0701557`
- `iv_window_orders(0.2-0.3V) = 0.4945596`
- `terminal_sum_abs_A_per_um_at_0p3 = 1.7293e-19`
- `bv_orders_at_0p05 = 0.1109109`
- `strict_newton_handoff_all = true`

Round3 decision:

- Extending thresholds up to `0.295 V` still does not separate candidate
   behavior in this deck.
- The threshold axis is now considered exhausted for `n_contact_only` in M2;
   next progress must come from a genuinely non-threshold mechanism.

### M2 round4 result: n-only strength refinement (2026-05-27)

After adding a continuous minority-relaxation strength knob, a small strength
matrix was executed with a fixed `0.1 V` activation threshold:

- Command:
   - `python scripts/scan_pn2d_n_only_thresholds.py --strengths 0.0,0.25,0.5,0.75,1.0 --strength-threshold 0.1 --summary-prefix pn2d_contact_relax_round4_n_only_strength_summary --candidate-dirname candidates_round4`
- Summary artifacts:
   - `build/pn2d_contact_relax_scan/pn2d_contact_relax_round4_n_only_strength_summary.csv`
   - `build/pn2d_contact_relax_scan/pn2d_contact_relax_round4_n_only_strength_summary.json`

Round4 candidate set:

- `baseline`
- `n_only_str0p0`
- `n_only_str0p25`
- `n_only_str0p5`
- `n_only_str0p75`
- `n_only_str1p0`

Observed outcome:

- All `n_only_str*` points are again numerically identical in this deck:
   - `I(0.29)/I(0.30) = 0.7025914`
   - `delta_to_reference = 0.0701557`
   - `iv_window_orders(0.2-0.3V) = 0.4945596`
   - `terminal_sum_abs_A_per_um_at_0p3 = 1.7293e-19`
   - `bv_orders_at_0p05 = 0.1109109`
   - `strict_newton_handoff_all = true`
- Baseline remains unchanged.

Round4 decision:

- The new continuous strength axis does not create separation either.
- For this pn2d deck, both threshold and strength sweeps now point to the
  same fixed `n_contact_only` response.
- The remaining path forward must be a different non-threshold mechanism, not
  another scalar tweak of the same relaxation policy.

### M3: Generalization verification and documentation convergence

Validation matrix:
- pn2d current reference import;
- simplified pn_diode deck;
- at least one nmos2d/pmos2d two-terminal or mixed-terminal smoke;
- existing BV/SRH gates;
- regression finite-output path.

Deliverables:
- Update `docs/validation/pn2d_sentaurus_comparison.md` with fresh pn2d IV/BV
   import metrics and probe paths.
- Keep this follow-up plan file synchronized with accepted gate values and
   unresolved issues.

### M3 execution result (2026-05-27)

M3 generalization validation was executed on current HEAD with:

- Full preset gate:
   - `ctest --preset windows-ucrt64-debug`
   - Result: `274/274` passed, `0` failed.
- Engineering regression matrix artifact:
   - `build/regression_output/regression_summary.json`

Focused matrix evidence (all `passed=true`, all rows converged):

| case | rows | converged_rows | final_current_total |
| --- | ---: | ---: | ---: |
| `pn_diode_iv` | 3 | 3 | `3.481904227298678e-03` |
| `pn_diode_bv` | 3 | 3 | `-5.78531939292131e-11` |
| `nmos2d_dd_iv` | 3 | 3 | `5.009164630732429` |
| `pmos2d_dd_iv` | 3 | 3 | `-1.7810363137092906` |
| `nmos2d_mos_dd_iv` | 3 | 3 | `1.62055767856109e-06` |
| `pmos2d_mos_dd_iv` | 3 | 3 | `-5.77996479489852e-07` |
| `nmos2d_mos_dd_bv` | 3 | 3 | `1.72038242915746e-07` |
| `pmos2d_mos_dd_bv` | 3 | 3 | `-6.3310014841849e-08` |

M3 status decision:

- pn2d + simplified pn/PMOS/NMOS + mixed-material smoke matrix is green under
   the current promoted iter2 baseline.
- M3 acceptance is marked complete for this cycle; unresolved work remains the
   pn2d high-bias local slope improvement branch (new M2 direction, not
   threshold-only tuning).
- Preserve reproducibility commands and artifact paths.

Open items to keep visible:
- external source alignment for Sentaurus default SRH parameters;
- Avalanche/OkutoCrowell parity;
- broader contact model generalization beyond pn2d.
