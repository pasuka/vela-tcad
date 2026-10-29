# Vela Documentation

This directory is organized around current project behavior: build it, choose
a configuration schema, import reference fixtures, and inspect validation
evidence.

## Current References

- [SLOT actual SRH/Auger branches and local chain checks, 2026-10-27](validation/slot_ldmos_srh_auger_branch_chain_2026-10-27.md): live source operands and partials, fixed-branch directional checks, high-precision replay and separate native-consumer limits; production acceptance unchanged.

- [SLOT native 80-bit QF/DOS chain checks, 2026-10-26](validation/slot_ldmos_native80_qf_dos_chain_2026-10-26.md): typed cache binding, exact eta replay, frozen-state directional checks and asymmetric weak recombination responses; full native Newton chain remains separate.

- [SLOT P1 acceptance checklist, 2026-10-21](validation/slot_ldmos_p1_acceptance_checklist_2026-10-21.md): fixed shadow-diagnostics scope, verified evidence, implementation and coverage gaps, closure conditions, and separate P2/native-alignment follow-ups; no change to acceptance rules.
- [SLOT P1 closeout batch, 2026-10-21](validation/slot_ldmos_p1_closeout_batch_2026-10-21.md): independent 24-hour observer build, update-chain and rejected-trial capture, paired validation, and explicitly limited historical substitution.

- [SLOT native low-bias convergence replay, 2026-10-21](validation/slot_ldmos_native_convergence_replay_2026-10-21.md): Sentaurus iteration capture, SLOT update-formula replay, P1 shadow limits, current-definition checks and an isolated IdVd doping-mismatch attribution; production criteria unchanged.

- [SLOT P1 shadow convergence diagnostics, 2026-10-20](validation/slot_ldmos_p1_shadow_convergence_2026-10-20.md): default-off observer, single-thread invariance pairing, exact archived-state replay and explicit P1 coverage limits; legacy acceptance unchanged.

- [SLOT weak-row classification, update metrics and bounded repeats, 2026-10-20](validation/slot_ldmos_weak_rows_repeat_adjoint_2026-10-20.md): saved-state edge classification, actual-update action, three U metrics, regional conservation, independent adjoint checks and an ongoing four-hour B repeat campaign.
- [SLOT scale binding, bounded span and Fermi-state audit, 2026-10-20](validation/slot_ldmos_scale_span_fermi_audit_2026-10-20.md): matched BVDS archive/build constants, diagnostic controls, same-parent span probes, IdVd population audit, and corrected weak-row physical-current units.

- [SLOT BVDS three-arm overnight attempt, 2026-10-19](validation/slot_ldmos_bvds_three_arm_overnight_2026-10-19.md): two fresh zero-bias starts with different high-voltage step caps and one qualified 170 mV restart; common twelve-hour deadline, full qualification and explicit resistor confirmation.

- [SLOT BVDS 170 mV tail and load-line confirmation, 2026-10-19](validation/slot_ldmos_bvds_170mv_tail_loadline_2026-10-19.md): one extra qualified update per path, full readbacks and production resistor seed confirmation; prior current path screen clarified as report-only for this case.

- [SLOT BVDS structure-preserving adjacent points, 2026-10-19](validation/slot_ldmos_bvds_poly_structure_adjacent_2026-10-19.md): full qualification and zero-update readback at 160/170 mV, unchanged homogeneous electron block, and a stopped current path-reproduction screen; no long scan.

- [SLOT BVDS homogeneous polysilicon structure-preserving solve, 2026-10-19](validation/slot_ldmos_bvds_poly_structure_preserving_2026-10-19.md): exact eligibility guards, original full-system defect checks, and paired bounded updates with repeated full qualification and archive readback; no new bias scan.

- [SLOT BVDS first-direction capture and thread replay, 2026-10-19](validation/slot_ldmos_bvds_first_direction_thread_replay_2026-10-19.md): actual PARDISO CSR/RHS and affine lift identity, one/four-thread fixed-system repeats, and homogeneous polysilicon direction drift; no new bias continuation.

- [SLOT BVDS production-cache identity and step-failure association, 2026-10-19](validation/slot_ldmos_bvds_cache_identity_step_association_2026-10-19.md): same-run population/SG generation evidence, source derivative identity checks, and an independent read-only audit of repeated small steps near 0.15576 V.

- [SLOT BVDS gap/tangent probes, 2026-10-19](validation/slot_ldmos_bvds_gap_tangent_preparation_2026-10-19.md): saved-operand symmetric decomposition and bounded VM finite perturbations, with explicit state/readback and branch-evidence limits; no Newton updates or new bias scan.

- [SLOT BVDS same-state operand replay, 2026-10-19](validation/slot_ldmos_bvds_same_state_operand_probe_2026-10-19.md): frozen 25/50 mV native SG/source captures, independent arithmetic replay, population/drive attribution and verified raw contact-reaction ledgers; no Newton updates.

- [SLOT BVDS qualified 50 mV twelve-hour continuation, 2026-10-19](validation/slot_ldmos_bvds_50mv_resume12_2026-10-19.md): independent qualified-parent continuation, staged adaptive device steps, explicit production resistor confirmation and fixed global deadline.

- [SLOT BVDS 5 mV cold-parent short path, 2026-10-19](validation/slot_ldmos_bvds_5mv_cold_path_2026-10-19.md): one-hour independent zero-bias requalification, adaptive device-voltage continuation with 25/50 mV checkpoints, unchanged full qualification and separate supply-voltage ledger.

- [SLOT BVDS overnight bounded directions, 2026-10-18](validation/slot_ldmos_bvds_overnight_directions_2026-10-18.md): four independent fixed-device, augmented-circuit, direction/coordinate and operand-arithmetic groups; fixed 12-hour window and unchanged full qualification.

- [SLOT BVDS charge-constant compatibility and local residuals, 2026-10-18](validation/slot_ldmos_bvds_charge_compat_local_residuals_2026-10-18.md): isolated q/epsilon profile with fixed Vt, 25/50 mV population/QF contracts, exact residual ledgers and remaining transport-sensitive weak rows.

- [SLOT BVDS 25/50 mV same-state matrix comparison, 2026-10-18](validation/slot_ldmos_bvds_25_50mv_matrix_comparison_2026-10-18.md): matched actual device biases, population/QF state contracts, physical Jacobian/RHS blocks and augmented contact equations.

- [SLOT BVDS cross-block direction and finite-step effects, 2026-10-18](validation/slot_ldmos_bvds_crossblock_direction_effect_2026-10-18.md): latest frozen-state population/QF contracts, source-branch attribution limits, fixed-RHS derivative controls and full-physics residual/port evaluations.

- [SLOT BVDS weak-floor report policy and short paths, 2026-10-18](validation/slot_ldmos_bvds_floor_report_short_path_2026-10-18.md): explicit user-authorized weak-row reporting policy; non-floor rows, closure and boundary protections remain enforced.

- [SLOT BVDS ordered contract and short validation, 2026-10-18](validation/slot_ldmos_bvds_ordered_contract_validation_2026-10-18.md): fixed evidence contract, one-direction archived trials, constant-material SG coordinates, contact-voltage low parts and original-qualification gates.

- [SLOT BVDS boundary projection and weak-row response, 2026-10-18](validation/slot_ldmos_bvds_boundary_projection_weak_response_2026-10-18.md): default-off synchronized-boundary diagnostic, strict paired-direction stop, original qualification, and frozen PolySilicon common-mode Jacobian row-sum audit.

- [SLOT BVDS 75 mV path and bounded nonzero-supply diagnostic, 2026-10-18](validation/slot_ldmos_bvds_75mv_path_bordered_2026-10-18.md): strict path stop, poly weak-row and actual low-state audit, opt-in bordered interface regressions and contact-voltage rounding limitation.

- [SLOT BVDS independent high-precision directions and fresh cold qualification, 2026-10-18](validation/slot_ldmos_bvds_high_precision_direction_cold_2026-10-18.md): 113/160-bit frozen-system references, closure-sum rounding audit and bounded new-boundary cold qualification.

- [SLOT BVDS pivot controls and split boundary targets, 2026-10-18](validation/slot_ldmos_bvds_pivot_split_boundary_2026-10-18.md): frozen-system pivot/scaling controls, target-row defect audit, split Dirichlet high/low targets, voltage-column regression and strict repeated-read qualification.

- [SLOT BVDS frozen-system refinement and boundary synchronization, 2026-10-18](validation/slot_ldmos_bvds_fixed_refinement_boundary_sync_2026-10-18.md): unsuccessful reused-factor corrections, old/new-bias qualification controls, first-read Dirichlet low-tail clearing and idempotent repeated read.

- [SLOT BVDS regularization strength paired updates, 2026-10-18](validation/slot_ldmos_bvds_regularization_strength_2026-10-18.md): frozen-parent 1e-8/1e-10/0 controls, strict qualification failures, tail linear defects and zero-update drain low-state changes.

- [SLOT BVDS equilibrium reference and weak-row regularization, 2026-10-17](validation/slot_ldmos_bvds_equilibrium_weak_rows_2026-10-17.md): same-potential common-QF controls, 220-digit SG replay, explicit contact-constraint ledger correction, and saved-direction explanation of residual retention on all 3,878 weak rows.

- [SLOT BVDS screening off/on relinearization, 2026-10-17](validation/slot_ldmos_bvds_ial_screening_relinearization_2026-10-17.md): same-parent three-update paired control, 39,150-dimensional device/circuit solves, unchanged strict qualification failures, exact low-state and boundary audits, and verified evidence package.

- [SLOT BVDS opt-in screening compatibility, 2026-10-17](validation/slot_ldmos_bvds_ial_screening_compatibility_2026-10-17.md): default-preserving compatibility selection, moving-cutoff thermal chain rule, one-sided branch validation, actual source derivatives and explicitly reproduced baseline test failure.

- [SLOT BVDS IALMob screening and derivative-cache identity, 2026-10-17](validation/slot_ldmos_bvds_ial_screening_derivative_identity_2026-10-17.md): directly observed screening cutoff/G, independent one-step reconstruction, scattering/HFS/source replay, and normalized own-node source derivative identities.

- [SLOT BVDS low-field/HFS and carrier-source decomposition, 2026-10-17](validation/slot_ldmos_bvds_hfs_source_detail_2026-10-17.md): native vertex low-field operands, same-contract HFS reconstruction, actual carrier-resolved source caches, and resolved branch-controlled source derivatives.

- [SLOT BVDS state contract, SG/HFS and avalanche pairing, 2026-10-17](validation/slot_ldmos_bvds_sg_hfs_state_contract_2026-10-17.md): same-system cache phases, actual SG and mobility operands, explicit AvalDensGradQF semantics, source scatter and branch-aware derivative diagnostics.

- [SLOT BVDS population generation and QF cache audit, 2026-10-17](validation/slot_ldmos_bvds_population_qf_cache_2026-10-17.md): same-population versus generation-input-QF controls, runtime band/DOS/ni operands, cache conversion stages, and finite-drive source value/derivative comparisons.

- [SLOT BVDS direct source and same-input SG arithmetic, 2026-10-17](validation/slot_ldmos_bvds_direct_source_same_input_2026-10-17.md): actual SRH/Auger numerator, denominator and scatter, independently matched support/lifetime, primary-variable equilibrium consistency, and 18 frozen PolySilicon edge/carrier arithmetic controls.

- [SLOT BVDS actual SG operand back-substitution, 2026-10-17](validation/slot_ldmos_bvds_sg_operand_backsubstitution_2026-10-17.md): runtime mobility operands, 23 active-edge stencils, actual-RHS row reconstruction, source/scatter remainder limits, and native default versus DirectCurrent contact-current definitions.

- [SLOT BVDS SG diagnostic export repair, 2026-10-17](validation/slot_ldmos_bvds_sg_diagnostic_export_fix_2026-10-17.md): carrier-specific production prefactors, independent 200-digit arithmetic reference, exact split-input sidecar, zero-support availability, and unchanged fixed-state RHS/current evidence.

- [SLOT BVDS split ledger and native SG operands, 2026-10-17](validation/slot_ldmos_bvds_split_ledger_native_sg_operands_2026-10-17.md): actual RHS/contact-reaction ledger, runtime SG operands and ordered accumulation, independent decimal arithmetic, and the observer/production population-prefactor boundary.

- [SLOT BVDS zero-bias RHS and bounded follow-up, 2026-10-17](validation/slot_ldmos_bvds_zero_rhs_followup_2026-10-17.md): constitutive/contact-reaction separation, six descending coupled updates, unchanged row/closure qualification failures, explicit boundary audit and separately passing equilibrium references.

- [SLOT BVDS zero-bias bordered iteration, 2026-10-17](validation/slot_ldmos_bvds_zero_bordered_iteration_2026-10-17.md): real 39,150-dimensional PARDISO device/circuit system, affine carrier constraints, staged Poisson/Jacobian comparison and an unqualified single-update diagnostic.

- [SLOT BVDS weak flux, finite tail and equilibrium reference, 2026-10-16](validation/slot_ldmos_bvds_weak_flux_tail_equilibrium_2026-10-16.md): independent electron/hole stable-flux checks, source balance, rejected-tail writeback audit and a separately qualified zero-bias equilibrium reference.

- [SLOT BVDS near-zero Auger relinearization and manual audit, 2026-10-16](validation/slot_ldmos_bvds_auger_relinearization_manual_2026-10-16.md): residual-preserving three-arm, three-update diagnostics, unchanged strict qualification, archived low-part statistics and verified native Auger formula/branch documentation.

- [SLOT BVDS Auger equilibrium arithmetic, 2026-10-16](validation/slot_ldmos_bvds_auger_equilibrium_arithmetic_2026-10-16.md): direct-library zero-point properties, 100-digit arithmetic references, native source/derivative transition brackets and explicit limits on inferring the native internal zero convention.

- [SLOT BVDS hole-source branches and cross-carrier support, 2026-10-16](validation/slot_ldmos_bvds_compatibility_hole_source_branches_2026-10-16.md): paired hole-population directions, distinct Auger and avalanche branches, actual element support and local absolute derivative scales.

- [SLOT BVDS thermal/ni compatibility, 2026-10-16](validation/slot_ldmos_bvds_thermal_ni_compatibility_2026-10-16.md): isolated version-bound four-arm experiment, fV-scale constitutive alignment, source/matrix comparison and explicit thermal-baseline regression limitation.

- [SLOT BVDS native constitutive and source finite differences, 2026-10-16](validation/slot_ldmos_bvds_native_constitutive_source_fd_2026-10-16.md): actual DOS/BGN outputs, independent band-reference mapping reconstruction, validated native state perturbations and source-only ablation.

- [SLOT BVDS constitutive and contact-chain audit, 2026-10-16](validation/slot_ldmos_bvds_constitutive_contact_chain_audit_2026-10-16.md): corrected eliminated-contact chain transformation, thermal-scale mapping trend, fixed-population finite differences and current VM task stops.

- [SLOT BVDS same-state assembly and actual PARDISO capture, 2026-10-16](validation/slot_ldmos_bvds_same_state_pardiso_capture_2026-10-16.md): matched native populations, QF-to-density Jacobian differences, native-QF control, and same-process CSR/RHS/returned-direction verification at the actual solver interface.

- [SLOT BVDS split-state 24-hour continuation and native single iteration, 2026-10-16](validation/slot_ldmos_bvds_C_resume24_native_iteration_2026-10-16.md): independent strict-policy hot continuation, preserved original deadline after a packed-predictor failure, and compensated-doping native binary64 input/direction/writeback capture.

- [SLOT BVDS C-arm tail and load-line validation, 2026-10-16](validation/slot_ldmos_bvds_C_tail_loadline_validation_2026-10-16.md): same-target diagnostic block corrections, original full qualification, actual external-resistor control and decoded-state confirmation.

- [SLOT BVDS C-arm split-state adjacent validation, 2026-10-16](validation/slot_ldmos_bvds_C_split_adjacent_2026-10-16.md): same-parent nearby targets and sampled path reproduction under strict qualification, exact residual reload, and checkpoint QF-frame restoration.

- [SLOT BVDS C-arm production low-state validation, 2026-10-16](validation/slot_ldmos_bvds_C_production_low_state_2026-10-16.md): opt-in IALMob/SRH/Auger/avalanche integration, strict qualification and exact reload, explicit low-zeroing control, and a same-target warm-start reproduction.

- [SLOT BVDS C-arm compensated Poisson validation, 2026-10-16](validation/slot_ldmos_bvds_C_compensated_poisson_2026-10-16.md): same-direction hi/lo versus ordinary coordinates, nonlinear Poisson reduction, loss on projection into ordinary DD, and local contact-guard evidence without changing qualification.

- [SLOT BVDS C-arm paired and frozen-state validation, 2026-10-16](validation/slot_ldmos_bvds_C_paired_poisson_contact_validation_2026-10-16.md): matched hot/secant failures, direct packed-writeback quantization evidence, and source-edge current/continuity sensitivity with unchanged qualification.

- [SLOT BVDS C-arm contact-guard stop analysis, 2026-10-16](validation/slot_ldmos_bvds_C_guard_failure_2026-10-16.md): eight bounded failures at 29.233 V, independent source-edge QF reconstruction, Poisson plateau, and unchanged qualification rules.

- [SLOT BVDS disk-stop resumption and 24-hour extension, 2026-10-16](validation/slot_ldmos_bvds_disk_resume_extension_2026-10-16.md): six independently supervised warm continuations, source-specific frozen solvers and acceptance rules, preserved pending spans, and original deadlines extended by 24 hours.

- [SLOT compensated BVDS hot-start 1 V ceiling, 2026-10-15](validation/slot_ldmos_compensated_bvds_hot1_cold_2026-10-15.md): independent 24-hour zero-bias cold start, same B strategy with only the span ceiling increased.

- [SLOT compensated BVDS parallel cold starts, 2026-10-15](validation/slot_ldmos_compensated_bvds_parallel_cold_2026-10-15.md): four independently supervised 24-hour cold-start arms comparing hot/secant initialization and step ceilings under unchanged strict qualification.

- [SLOT compensated BVDS adjacent validation, 2026-10-15](validation/slot_ldmos_compensated_bvds_adjacent_2026-10-15.md): strict cold-to-100 mV qualification with one bounded step reduction, four accepted points and unchanged physics.

- [SLOT compensated cold recovery, 2026-10-15](validation/slot_ldmos_compensated_cold_recovery_2026-10-15.md): bounded recovery obtains a strictly qualified zero state; paired continuation fails, and frozen probes quantify Poisson writeback loss without changing acceptance rules.
- [SLOT compensated cold failure analysis, 2026-10-15](validation/slot_ldmos_compensated_cold_failure_analysis_2026-10-15.md): coupled cold initialization, residual plateau, weak rows, gate contact edge and the unused recovery path; no restart or tolerance change.
- [SLOT compensated-doping BVDS cold campaign, 2026-10-15](validation/slot_ldmos_compensated_bvds_cold_campaign_2026-10-15.md): independent 24-hour VM cold-start run with original active compensation, current discretization, strict BVDS qualification and detached supervision.
- [SLOT original active-doping node audit, 2026-10-15](validation/slot_ldmos_active_doping_node_audit_2026-10-15.md): original As/P/B active components match the old compensated Vela input; the recent net-split input differs despite matching net doping.
- [SLOT BVDS total-doping cold initialization, 2026-10-15](validation/slot_ldmos_bvds_total_doping_cold_2026-10-15.md): independent native control and compensated-doping zero states, preserved input fields, BGN agreement and scoped mobility comparison to the old Vela model.
- [SLOT BVDS tdx doping comparison, 2026-10-15](validation/slot_ldmos_bvds_tdx_doping_comparison_2026-10-15.md): GRD/DAT conversion independently confirms input versus cold-state donor/acceptor differences, with source-neighborhood coordinates and conversion limitations.
- [SLOT BVDS September/October history comparison, 2026-10-15](validation/slot_ldmos_bvds_history_comparison_2026-10-15.md): verified historical high-voltage states, changed doping/discretization/reference frame, and conditional Poisson/contact stall evidence.
- [SLOT IdVd 1.10 V final analysis and plan, 2026-10-15](validation/slot_ldmos_idvd_1p10_analysis_plan_2026-10-15.md): completed research chain, early-stop versus local recovery, remaining material-specific errors and scoped next steps.
- [SLOT overnight and daily work summary, 2026-10-14 to 15](validation/slot_ldmos_work_summary_2026-10-14_night_to_15.md): source-backed daily-report material separating completed diagnostics, numerical limitations and ongoing tasks.
- [SLOT BVDS Poisson writeback audit, 2026-10-15](validation/slot_ldmos_bvds_poisson_writeback_audit_2026-10-15.md): frozen linear directions, sub-ULP packed updates, residual readback and independent contact guard evidence.
- [SLOT BVDS local double/binary128 solve, 2026-10-15](validation/slot_ldmos_bvds_precision_local_solve_2026-10-15.md): packed PolySilicon restore consistency repair; paired local solves retain Poisson/contact limitations despite carrier closure passing.
- [SLOT BVDS extended-precision prefactor repair, 2026-10-15](validation/slot_ldmos_bvds_precision_prefactor_2026-10-15.md): corrected material-reference populations, unchanged double baseline, gauge/ULP limitations and frozen contact/load-line sensitivity.
- [SLOT IdVd path reproduction and BVDS Poisson/contact guard diagnostics, 2026-10-15](validation/slot_ldmos_path_reproduction_poisson_guard_2026-10-15.md): 1.030 V path comparison, bounded continuation status, matched hot/secant failures, ULP sensitivity and contact-edge threshold evidence.
- [SLOT IdVd adjacent results and BVDS stall analysis, 2026-10-15](validation/slot_ldmos_idvd_adjacent_bvds_stall_analysis_2026-10-15.md): recovered dilute populations, remaining polysilicon rows, and the contact-edge QF veto at the BVDS Poisson residual plateau.
- [SLOT BVDS overnight analysis and 24-hour continuation, 2026-10-15](validation/slot_ldmos_bvds_overnight_analysis_resume_2026-10-15.md): budget versus numerical stops, Poisson/contact-QF stall evidence, and three independently resumed strict-policy arms.
- [SLOT coordinated-source initial-state and adjacent-point validation, 2026-10-15](validation/slot_ldmos_source_coordination_adjacent_2026-10-15.md): bounded research campaign, previous-job status and contact-consistent warm-start checks; see report for completion status.
- [SLOT transport/source linearization coordination, 2026-10-15](validation/slot_ldmos_source_transport_coordination_2026-10-15.md): same-seed four-arm diagnostics, released-row balance and unchanged residual checks.
- [SLOT dilute mapping and temporary density hold, 2026-10-15](validation/slot_ldmos_inverse_hold_parallel_validation_2026-10-15.md):
  extreme-density arithmetic repair, four bounded same-seed diagnostic arms,
  and original versus actual update-direction/source-derivative attribution.

- [SLOT BVDS aligned cold start and three matched targets, 2026-10-14](validation/slot_ldmos_bvds_aligned_cold_three_points_2026-10-14.md):
  independent AverageBox and effective-doping configuration, strict fresh-cold
  qualification, and bounded nonzero-bias recovery with four-thread PARDISO.

- [SLOT BVDS zero-bias same-state comparison, 2026-10-14](validation/slot_ldmos_bvds_zero_state_comparison_2026-10-14.md):
  matched populations expose compensated-doping/BGN and contact differences;
  physical matrix comparison separates native drain auxiliary equations.

- [SLOT PhuMob hotspot and secant validation, 2026-10-14](validation/slot_ldmos_phumob_hotspot_secant_2026-10-14.md):
  exact invariant caching removes repeated G-minimum searches; fixed-state matrix,
  direction and field checks precede bounded two-point drain continuation.

- [SLOT overnight cost audit and native references, 2026-10-14](validation/slot_ldmos_overnight_cost_native_reference_2026-10-14.md):
  first-update residual guards limit IdVd spans; repeated bounded recovery suppresses
  BVDS growth; independent native reference jobs export matrices and physical fields.

- [SLOT BVDS guarded linear refinement, 2026-10-13](validation/slot_ldmos_bvds_guarded_refinement_2026-10-13.md):
  opt-in ordinary Newton refinement preserves old BVDS physics and qualification;
  fixed-matrix correction passes, while the paired failed target remains limited by nonlinear closure.

- [SLOT VM refined adjacent points, 2026-10-13](validation/slot_ldmos_vm_refined_adjacent_validation_2026-10-13.md):
  three first-attempt diagnostic continuations to 42 mV with retained-factor refinement,
  unchanged linear guards and exact saved-state residual/current readback.

- [SLOT VM refinement and shared diagnostics, 2026-10-13](validation/slot_ldmos_vm_refinement_shared_diagnostics_2026-10-13.md):
  unchanged linear guards pass with retained-factor refinement; bounded 6 to
  13.008 mV validation and bitwise paired shared-evaluation checks, with measured costs.

- [SLOT VM linear guard and diagnostic cost audit, 2026-10-13](validation/slot_ldmos_vm_linear_guard_cost_audit_2026-10-13.md):
  six rejected steps share one interface hole row; frozen-system residual
  correction passes unchanged linear guards, with dual-system timing limits.

- [SLOT Linux VM campaign and PARDISO threading, 2026-10-13](validation/slot_ldmos_linux_vm_long_campaign_2026-10-13.md):
  four-thread same-parent replay passed; detached twelve-hour IdVd and independent
  old-policy BVDS campaign launched, with remote deadline supervision.

- [SLOT Linux VM short validation, 2026-10-13](validation/slot_ldmos_linux_vm_short_validation_2026-10-13.md):
  GCC16 and Linux PARDISO verified at three adjacent drain biases with the
  diagnostic transport linearization; final readback completed in a separate
  zero-update audit. Sentaurus execution awaits its license configuration.

- [Remote compute VM](office_vm_workflow.md): private SSH configuration,
  Micromamba GCC 16 builds, isolated Sentaurus environment, file transfer and
  detached job submission/status through `scripts/Invoke-VelaVm.ps1`.

- [SLOT native derivative-switch audit, 2026-10-13](validation/slot_ldmos_native_derivative_switch_2026-10-13.md):
  same-state Sentaurus mobility/avalanche derivative controls preserve the solve RHS;
  111 negative hole diagonals and omitted terms match the Vela diagnostic closely,
  while wider source coupling and historical scan context remain distinct.

- [SLOT mobility-derivative adjacent paths, 2026-10-13](validation/slot_ldmos_mobility_adjacent_path_2026-10-13.md):
  bounded local/T470p continuation from 6 mV, different-path field and contact
  reaction checks, and field changes after the first RHS-or-update stopping point;
  diagnostic linearization only, with original residuals and guards retained.

- [SLOT projection and mobility-derivative recovery, 2026-10-13](validation/slot_ldmos_projection_mobility_recovery_2026-10-13.md):
  local/T470p bounded diagnostics distinguish repeated positivity projection,
  an unsuccessful activity-set release rule, and transport derivative feedback;
  original residuals and acceptance thresholds remain unchanged.

- [SLOT retained-potential three-step trial, 2026-10-13](validation/slot_ldmos_keep_psi_three_steps_2026-10-13.md):
  retain 126 contact-potential constraints in the 39,149-dimensional PARDISO
  system; identical initial assembly/RHS, three accepted adaptive drain points,
  and remaining low-density hole-QF path sensitivity.

- [SLOT IdVd three-day handoff, 2026-10-11 to 13](validation/slot_ldmos_idvd_handoff_2026-10-11_to_13.md):
  evidence-backed configuration and acceptance history, linear replay versus
  nonlinear field differences, current adaptive-step failures, and bounded
  diagnostic questions for an independent model review.

- [SLOT native-style steps toward 1 V, 2026-10-13](validation/slot_ldmos_native_steps_to_1v_2026-10-13.md):
  bounded local 6 mV start, iteration-aware growth and halving with Fermi QF
  extrapolation; six accepted points to 13.002 mV, numerical guards and budget stops.

- [SLOT RHS-or-update stopping policy, 2026-10-13](validation/slot_ldmos_rhs_or_update_policy_2026-10-13.md):
  current SLOT opt-in entry reports closure/boundary thresholds, retains numerical
  and saved-state guards, and verifies RHS-only and update-only acceptance.

- [SLOT single-process population loop, 2026-10-13](validation/slot_ldmos_population_memory_2026-10-13.md):
  opt-in in-memory updates and linear guards; fixed-parent 5/12.1 mV replay,
  unchanged fields and acceptance, bounded timing and output comparison.

- [SLOT acceptance and diagnostic efficiency, 2026-10-12](validation/slot_ldmos_population_acceptance_efficiency_2026-10-12.md):
  update-or-RHS acceptance with unchanged safeguards; paired frozen-state checks
  of compact matrix output and removal of repeated port/residual diagnostics.

- [SLOT population projection pair, 2026-10-12](validation/slot_ldmos_population_projection_pair_2026-10-12.md):
  fixed PARDISO four-update comparison of 0.001/0.01 positivity projection;
  0.001 qualifies, while 0.01 improves some hole populations but retains a Poisson residual gap.

- [SLOT same-direction population update discrimination, 2026-10-12](validation/slot_ldmos_same_direction_update_discrimination_2026-10-12.md):
  native and current rules diverge at first-step hole positivity projection;
  a single-factor 0.001/0.01 replay preserves the other fields but does not restore closure.

- [SLOT native system through production linear interfaces, 2026-10-12](validation/slot_ldmos_native_linear_interface_replay_2026-10-12.md):
  PARDISO and UMFPACK reproduce native float64 directions; single-thread results
  match direct library calls bitwise, with independent residual and input audits.

- [SLOT 5 mV refinement and contact-hole ledger, 2026-10-12](validation/slot_ldmos_5mv_linear_refinement_hole_ledger_2026-10-12.md):
  explicit residual correction repairs the frozen UMFPACK linear error;
  contact source accounting revises the source-hole comparison from 26.1% to 5.16 ppm.
  Low-density population projection and field differences remain unresolved.

- [SLOT exact 5 mV revalidation versus native float64, 2026-10-12](validation/slot_ldmos_5mv_native64_revalidation_2026-10-12.md):
  PARDISO passes scoped acceptance in four updates; UMFPACK is rejected by the
  third linear guard. Drain current differs by 7.59 ppm, with remaining hole-field differences.

- [SLOT 6 mV final states versus Sentaurus, 2026-10-12](validation/slot_ldmos_6mv_native_final_comparison_2026-10-12.md):
  exact-bias native comparison finds 7.59 ppm drain-current difference but
  substantial low-density hole-QF and source hole-current differences.

- [SLOT near-zero linear guard and backend pair, 2026-10-12](validation/slot_ldmos_mixed_linear_guard_pair_2026-10-12.md):
  opt-in mixed residual guard retains block and constraint checks; same-parent
  6 mV UMFPACK/PARDISO relinearization, with unchanged nonlinear acceptance.

- [SLOT T470p population update and strict elimination, 2026-10-12](validation/slot_ldmos_t470p_population_elimination_2026-10-12.md):
  paired UMFPACK/PARDISO pass scoped 10 mV acceptance; adjacent steps expose
  a finite-population inverse-mapping guard and a separate linear-error outlier.

- [SLOT population-update linear backend comparison, 2026-10-11](validation/slot_ldmos_population_linear_backends_2026-10-11.md):
  SparseLU and UMFPACK repeat four stabilization combinations against PARDISO;
  scoped acceptance and current agree while low-density hole-QF trajectories differ.

- [SLOT regularization and closure row replacement comparison, 2026-10-11](validation/slot_ldmos_stabilization_ablation_2026-10-11.md):
  four-arm population finite-step comparison with independent stabilization
  controls, unchanged physical equations, and preserved closure acceptance.

- [SLOT population finite-step verification, 2026-10-11](validation/slot_ldmos_population_finite_step_verification_2026-10-11.md):
  bounded relinearization reaches the scoped acceptance at 10 mV; a separate
  trajectory reproduces drain current but retains low-density hole-QF differences.

- [SLOT component population update and frozen Charon-inspired kernels, 2026-10-11](validation/slot_ldmos_component_population_charon_diagnostics_2026-10-11.md):
  three-step same-parent update comparison, stagnation and backward-error audit,
  and separate low-density model sensitivities; no new qualified continuation point.

- [SLOT local Charon and DEVSIM Newton review, 2026-10-11](validation/slot_ldmos_charon_devsim_newton_review_2026-10-11.md):
  source-verified population positivity, potential damping, scaling and stagnation
  checks; separates numerical updates from low-density driving-field model changes.

- [SLOT Bank-Rose population-domain recovery, 2026-10-11](validation/slot_ldmos_bank_rose_population_recovery_2026-10-11.md):
  bounded rejection of typed invalid-population trials preserves the parent and
  physics; the paired uncapped continuation audit records qualification separately.

- [SLOT uncapped Bank-Rose comparison, 2026-10-11](validation/slot_ldmos_row_report_bank_rose_2026-10-11.md):
  replacing backtracking with Bank-Rose leaves the population-domain exception unresolved;
  a larger initial K permits six tiny updates but does not establish global closure.

- [SLOT row-report continuation without QF caps, 2026-10-11](validation/slot_ldmos_row_report_nocap_2026-10-11.md):
  10 mV and the 7.5 mV retry exit during the first trial population evaluation;
  the initial 10 mV systems match the capped baseline, with no new accepted updates.

- [SLOT row-report neighboring voltages and path reproduction, 2026-10-11](validation/slot_ldmos_row_report_neighbors_2026-10-11.md):
  three new-policy points reach 20 mV with closure enforced; a different voltage path
  reproduces current and potential, while low-density carrier QF fields remain path dependent.

- [SLOT row-report stopping policy and 5 mV verification, 2026-10-11](validation/slot_ldmos_row_report_policy_2026-10-11.md):
  the user-selected case-local policy stops at update nine with original closure enforced;
  strict row diagnostics and local minority-field differences remain separately reported.

- [SLOT PARDISO/SparseLU three-update comparison, 2026-10-11](validation/slot_ldmos_mkl_relinearization_2026-10-11.md):
  an opt-in MKL backend reduces some weak-row violations while unchanged QF updates
  leave the dominant hole imbalance unqualified; initial matrices and RHS are identical.

- [SLOT oneMKL PARDISO native-system replay, 2026-10-11](validation/slot_ldmos_mkl_native_replay_2026-10-11.md):
  raw float64 native systems reproduce the saved direction closely; common population
  finite steps separate linear accuracy from remaining nonlinear qualification.

- [SLOT native-style OR stopping diagnostic, 2026-10-11](validation/slot_ldmos_native_gate_trial_2026-10-11.md):
  both eliminated backends stop at update nine on the relaxed RHS branch;
  unchanged trajectories and original zero-update audits distinguish diagnostic success from row qualification.

- [SLOT native/Vela convergence criteria comparison, 2026-10-11](validation/slot_ldmos_convergence_criteria_comparison_2026-10-11.md):
  native update-error stopping and shadow population-update metrics distinguish
  global convergence indicators from unresolved Vela carrier-row qualification.

- [SLOT eliminated-system QF cap on/off comparison, 2026-10-11](validation/slot_ldmos_eliminated_cap_toggle_2026-10-11.md):
  both uncapped precisions exit on a first-trial SG population-domain exception;
  no accepted updates, distinguishing failed trial evaluation from exhausted backtracking.

- [SLOT eliminated-system 5 mV Newton trials, 2026-10-11](validation/slot_ldmos_eliminated_newton_trial_2026-10-11.md):
  opt-in double/long-double elimination backends pass solver tests; bounded
  nonlinear trials satisfy closure but both fail original carrier-row qualification.

- [SLOT fixed-system constraint elimination comparison, 2026-10-11](validation/slot_ldmos_constraint_elimination_precision_2026-10-11.md):
  exact affine elimination agrees with the high-precision reference; double fixes
  the 10280 sign but increases other direction errors, so elimination alone is not a recovery.

- [SLOT high-precision fixed-system direction reference, 2026-10-11](validation/slot_ldmos_high_precision_direction_reference_2026-10-11.md):
  113/160-bit references and independent decimal residual checks quantify weak
  direction errors; double has the wrong sign at 10280, without implying nonlinear recovery.

- [SLOT in-memory matrix and direction capture, 2026-10-11](validation/slot_ldmos_in_memory_direction_capture_2026-10-11.md):
  unchanged production trajectories and exact matrix replay locate the 10280 sign
  split in the finite-precision base solve; both nonlinear paths remain unqualified.

- [SLOT same-state precision direction layers, 2026-10-11](validation/slot_ldmos_float80_direction_layers_2026-10-11.md):
  regularization, Woodbury, QF caps and Poisson recorrection isolated; three-neighbor
  swaps reproduce the 10280 row difference, while original qualification still fails.

- [SLOT mixed-precision three-update comparison, 2026-10-11](validation/slot_ldmos_float80_relinearization_2026-10-11.md):
  identical-parent production updates reduce some carrier-row violations with
  long double, while dominant residuals and closure remain nearly unchanged.

- [SLOT float80 trial and optional mixed-precision backend, 2026-10-11](validation/slot_ldmos_float80_trial_2026-10-11.md):
  checked long-double precision, global type migration blockers, identical-matrix
  residual comparison, and zero-update SLOT probes; production defaults unchanged.

- [SLOT Bank–Rose implementation and bounded comparison, 2026-10-10](validation/slot_ldmos_bank_rose_2026-10-10.md):
  opt-in Algorithm Global damping, persistent K and NotDamped activation;
  numerical tests and a same-parent 5 mV comparison retain original qualification.

- [SLOT undamped-window and OR stopping comparison, 2026-10-10](validation/slot_ldmos_native_style_stopping_2026-10-10.md):
  the diagnostic stops on RHS after nine updates but fails original hole-row qualification;
  distinguishes raw-direction update error, applied updates and unchanged production gates.

- [SLOT native 64-bit and 80-bit precision comparison, 2026-10-10](validation/slot_ldmos_native_precision_2026-10-10.md):
  matched 5 mV native solves both stop after four updates; current agrees closely,
  while residual and near-zero quantities retain precision sensitivity.

- [SLOT linear backends and Newton algorithm review, 2026-10-10](validation/slot_ldmos_linear_backends_newton_review_2026-10-10.md):
  six backend controls share identical matrices and RHS, with finite-step checks;
  version-bound Sentaurus and PETSc references separate linear accuracy from globalization.

- [SLOT boundary elimination, stabilization and finite steps, 2026-10-10](validation/slot_ldmos_boundary_stabilization_finite_steps_2026-10-10.md):
  compares both matrix origins in a common free QF space, audits stable closure-row
  solves, and separates finite population accuracy from actual row and conservation quality.

- [SLOT cross-block flux and source decomposition, 2026-10-10](validation/slot_ldmos_crossblock_flux_source_decomposition_2026-10-10.md):
  corrects omitted contact-column coordinate terms and separates smooth SG transport
  from near-zero avalanche-current directional sensitivity using matched-state probes.

- [SLOT zero and exact 5 mV initial matrix comparison, 2026-10-10](validation/slot_ldmos_initial_zero_5mv_matrix_comparison_2026-10-10.md):
  original matrix statistics include a subsequently corrected contact-column mapping
  error; see the erratum and cross-block decomposition above. RHS statistics are unaffected.

- [SLOT matched-system population update comparison, 2026-10-10](validation/slot_ldmos_pure_population_update_comparison_2026-10-10.md):
  three-step controls reuse the original stabilized direction, reduce hole-row failures
  without full qualification, and distinguish population updates from Poisson recorrection.

- [SLOT contact-preserving relative population inverse, 2026-10-10](validation/slot_ldmos_contact_preserving_relative_qf_2026-10-10.md):
  exact zero-change and contact preservation pass, with manufactured increment recovery;
  native population matching still fails row qualification and exposes sub-ULP QF sensitivity.

- [SLOT population inverse mapping and flux representation, 2026-10-10](validation/slot_ldmos_population_inverse_flux_representation_2026-10-10.md):
  frozen controls distinguish exact QF re-anchoring, density round-trip perturbations,
  boundary remapping and stable SG arithmetic; all qualification gates remain unchanged.

- [SLOT regional closure, relinearization and constitutive profile, 2026-10-10](validation/slot_ldmos_four_stage_relinearization_2026-10-10.md):
  three-step diagnostic replays separate transient electron closure from persistent row failures;
  bounded population directions and an independently built constants profile remain unqualified,
  with a frozen population-matching control exposing charge/transport sensitivity.

- [SLOT constitutive constants and population update coordination, 2026-10-10](validation/slot_ldmos_mapping_population_coordination_2026-10-10.md):
  native-final-state formulas explain remaining mapping and Poisson constant terms;
  bounded population and cross-block regularization trials improve some rows but
  fail full qualification, with explicit linear-backend sensitivity controls.

- [SLOT native final-state qualification and direction layers, 2026-10-10](validation/slot_ldmos_native_qualification_direction_layers_2026-10-10.md):
  zero-update qualification distinguishes residual mapping from stopping criteria;
  identical physical matrices expose regularization and clipping effects, with
  an independently reconstructed first production update and unchanged gates.

- [SLOT effective-band transport implementation and bounded validation, 2026-10-10](validation/slot_ldmos_effective_band_transport_validation_2026-10-10.md):
  opt-in Fermi edge transport, complete small-mesh derivative tests and SLOT column checks;
  frozen Si electron residual agreement improves, but both nonlinear controls stop
  on carrier-row line-search rejection with no qualified new state.

- [SLOT Fermi SG edge and constitutive mapping, 2026-10-10](validation/slot_ldmos_fermi_sg_edge_constitutive_mapping_2026-10-10.md):
  frozen edge substitutions isolate a Fermi transport-weight difference; multi-node
  thermal-voltage and base-ni checks explain spatial QF offsets without changing production defaults.

- [SLOT full residual, SG transport and contact mapping, 2026-10-10](validation/slot_ldmos_full_residual_transport_mapping_2026-10-10.md):
  frozen population/QF controls separate constitutive driving-field differences from
  transport assembly; contact sets and source/transport units are audited without new accepted states.

- [SLOT source-kernel comparison with joint stopping rules, 2026-10-10](validation/slot_ldmos_source_trend_stop_revalidation_2026-10-10.md):
  opt-in row-quality/raw-direction/line-search monitoring and longer bounded controls;
  both arms stop on original line-search failure before time limits, with no new qualified point.

- [SLOT SG source-kernel alignment, 2026-10-10](validation/slot_ldmos_sg_source_alignment_2026-10-10.md):
  opt-in effective-band SG current reconstruction with partial-layer current/force
  projection, live population derivatives, source-support checks and bounded
  self-consistent controls. Frozen source agreement is separate from device qualification.

- [SLOT native source support and residual audit, 2026-10-09](validation/slot_ldmos_native_source_residual_audit_2026-10-09.md):
  same-state source splitting verifies runtime source volumes against actual RHS;
  SG current reconstruction and interface projection differ from Vela's density-gradient source.
  Fresh native assembly passes two fixed-population psi finite differences; earlier iteration
  matrices are not interchangeable with the independently reassembled physical derivative.

- [SLOT PMI source and cross-block comparison, 2026-10-09](validation/slot_ldmos_pmi_source_crossblock_comparison_2026-10-09.md):
  matched-population diagnostics reduce QF mapping differences to sub-microvolt
  levels, but hotspot source and psi cross-block differences remain. Two Vela
  derivative columns pass finite-difference checks; native source support remains unresolved.

- [SLOT native integral BGN PMI comparison, 2026-10-09](validation/slot_ldmos_native_integral_bgn_pmi_2026-10-09.md):
  GCC 11 PMI passes scalar/derivative tests and two zero-update native observations.
  Explicit BW plus normalized-integral statistics matches all exported Si/Poly
  nodes; no new self-consistent device state or full Jacobian equivalence is claimed.

- [SLOT gate continuation and PolySilicon BGN audit, 2026-10-09](validation/slot_ldmos_gate_continuation_poly_bgn_2026-10-09.md):
  four new gate points and independent reproduction pass to 0.160859375 V.
  Native normalized quadrature aligns PolySilicon population relations; a
  NoFermi diagnostic isolates the remaining Si statistical BGN correction.

- [SLOT aligned gate path and native integral route, 2026-10-09](validation/slot_ldmos_aligned_gate_integral_route_2026-10-09.md):
  three small gate points and a separate direct cold-to-target reproduction pass
  the original gates. Native RelativeToDOSMass quadrature aligns the Si population
  integral, but BGN statistical correction and PolySilicon remain unmatched.

- [SLOT AverageBox cold state and constitutive alignment, 2026-10-09](validation/slot_ldmos_average_box_cold_constitutive_alignment_2026-10-09.md):
  a separate configuration uses the actual AverageBox Poisson path and passes
  fresh cold-state qualification. The 8.214 ppm ratio is explained by physical
  constants; native DOS matches, while Fermi approximation and BGN mapping
  differences remain. No gate/drain continuation or stabilization ablation.

- [SLOT constitutive, Poisson and source audit, 2026-10-09](validation/slot_ldmos_constitutive_poisson_source_audit_2026-10-09.md):
  the signed per-cell Poisson override bypasses AverageBox; a frozen control
  aligns potential stencils, while charge scaling and Fermi/QF mapping still
  differ. Two hotspot cross-block derivatives pass local finite differences.
  Native intermediate fields are usable despite final-output timeout; no Vela
  continuation point is accepted.

- [SLOT same-state residual and Jacobian layers, 2026-10-09](validation/slot_ldmos_same_state_matrix_layers_2026-10-09.md):
  native free populations are matched before a full Fermi chain transformation.
  Contact constitutive mismatches and hotspot cross-block differences remain;
  weak-row regularization is quantified separately. No recovery step is accepted.

- [SLOT connected population trials with neighbor guards, 2026-10-09](validation/slot_ldmos_connected_population_guard_2026-10-09.md):
  twelve same-state frozen trials compare restricted and locally recomputed
  population directions. No trial passes true neighbor guards; new electron
  failures are already visible in linear predictions, motivating explicit
  protected-row constraints rather than soft least-squares weights alone.

- [SLOT population updates with Fermi inverse mapping, 2026-10-09](validation/slot_ldmos_population_fermi_update_2026-10-09.md):
  same-state bounded population trials improve hotspot recovery and Poisson
  residuals but create thousands of new hole-row failures. Mapping and frozen
  population checks pass; no trial is accepted for continuation.

- [SLOT second-step population recovery and QF clipping, 2026-10-09](validation/slot_ldmos_second_step_recovery_2026-10-09.md):
  exact first-step replay and a bounded second update distinguish a positive
  recovery direction clipped at node 10849 from continued depletion and
  cross-block coupling at 10237. Full qualification remains unmet; comparisons
  across different states and variable coordinates do not establish a Jacobian bug.

- [SLOT native Newton matrices and iteration comparison, 2026-10-09](validation/slot_ldmos_native_newton_matrix_2026-10-09.md):
  Sentaurus `-d` exports five full sparse systems; NewtonPlot captures a
  four-update 6 mV recovery with density positivity handling. Matrix/vector
  reconstruction and observer controls pass; cross-solver variable transforms
  and Vela final qualification remain distinct, unfinished checks.

- [SLOT bounded joint least-squares directions, 2026-10-09](validation/slot_ldmos_bounded_joint_direction_2026-10-09.md):
  equilibrated constrained subproblems yield millivolt-scale coupled updates and
  reduce all three residual blocks. Original failed rows and closure remain
  unqualified: hole source reduction has not established balancing transport.
  No recovery step or continuation state is approved.

- [SLOT regularization and cross-block direction comparison, 2026-10-09](validation/slot_ldmos_regularization_block_direction_2026-10-09.md):
  four regularization strengths and conditional hole-block directions do not
  provide an admissible 6 mV first step. Removing regularization yields enormous
  directions; clipping and block damping retain thousands of failed rows.
  All trials are diagnostic, with unchanged physics and final qualification.

- [SLOT 6 mV direction, clipping and merit decomposition, 2026-10-09](validation/slot_ldmos_first_step_direction_merit_2026-10-09.md):
  zero-update replay reproduces the first accepted trial and separates carrier
  clipping from Poisson recorrection. Global merit masks hole-row deterioration;
  source growth and regularized cross-block directions require coordinated repair.
  No production physics, final threshold or qualified-point count changes.

- [SLOT zero-state population/source and 6 mV first-step audit, 2026-10-09](validation/slot_ldmos_zero_state_source_first_step_2026-10-09.md):
  frozen native fields show up to 1.35% population remapping, while density-only
  substitution changes avalanche integration by 0.278%. Two bounded first-step
  arms retain thousands of failed rows; global merit decreases while hole
  residuals increase. No qualification or production physics change is implied.

- [SLOT zero-drain native field comparison, 2026-10-09](validation/slot_ldmos_zero_drain_native_comparison_2026-10-09.md):
  exact Vg=5 V / Vd=0 native TDR and Vela state share node coordinates,
  connectivity and materials. Channel-surface potential differs by at most
  1.8 µV; localized silicon differences reach 0.346 mV and 2.71% in density.
  Near-zero currents and source outputs require separate interpretation;
  this comparison does not resolve the nonzero-drain bootstrap failure.

- [SLOT 6 mV drain bootstrap short validation, 2026-10-09](validation/slot_ldmos_drain_6mv_short_validation_2026-10-09.md):
  the latest 5 V gate / zero-drain state revalidates, but 6, 3 and 1.5 mV
  warm starts all reach their bounded soft deadlines with thousands of failed
  carrier rows, mainly holes. No new drain state is qualified; neighboring
  extrapolation and reproduction remain conditional on establishing a first point.

- [SLOT gate warning-only continuation, 2026-10-09](validation/slot_ldmos_gate_warning_resume_2026-10-09.md):
  the predictor now logs changes above 0.1 V without rejecting them; the adaptive
  policy retains its 0.25 V maximum and halves failed spans. The separately
  authorized three-hour gate-only continuation starts from a revalidated
  0.496289 V state. See the report and live campaign evidence for completion.

- [SLOT gate 0.25 V predictor comparison, 2026-10-09](validation/slot_ldmos_gate_025_comparison_2026-10-09.md):
  one 0.25 V step and a 0.1/0.1/0.05 V path both pass the original full gates
  from Vg=0.246289 V to 0.496289 V at Vd=0. Maximum potential difference is
  8.22e-15 V; elapsed costs are 157 versus 391 seconds. This local diagnostic
  did not itself change the then-production 0.1 V predictor cap or establish
  full-ramp stability; the subsequent warning-only change is documented above.

- [SLOT IdVd three-hour continuation and diagnosis, 2026-10-08](validation/slot_ldmos_idvd_three_hour_resume_2026-10-08.md):
  original gates passed at 2.4 V and reproduced through a different fine voltage
  path; SG current agreement is 4.1e-10 relative. Frozen observations distinguish
  failed electron-row inflow from source terms. The budgeted run stopped at
  3.339357 V with 81 new main-path qualified points; three reproduction points
  are counted separately. Final record and hash checks pass; no 60 V claim.

- [SLOT IdVd small-step secant path, 2026-10-08](validation/slot_ldmos_idvd_small_step_path_2026-10-08.md):
  four new qualified points through 1.194279 V without fallback; endpoint
  reproduction from an earlier qualified seed pair passes the original gates.
  Total solver wall time is 816 s; original physics and thresholds retained.

- [SLOT IdVd electron weak rows and bounded continuation, 2026-10-08](validation/slot_ldmos_idvd_electron_weak_row_2026-10-08.md):
  two new qualified drain voltages through 0.921167 V with same-target reproduction;
  frozen HFS response dominates one old failed row, while connected restoration
  is rejected for neighbor damage. Original physics and qualification retained.

- [SLOT IdVd AverageBox implementation and bounded recomputation, 2026-10-07](validation/slot_ldmos_idvd_average_box_2026-10-07.md):
  opt-in unweighted 2-D face intersections reproduce native zero edges and
  element measures; 49 C++ cases and 19 Python checks pass. Three bounded
  controls remain unqualified with carrier-row and closure failures.
- [SLOT IdVd native box geometry and update decomposition, 2026-10-07](validation/slot_ldmos_idvd_update_decomposition_geometry_2026-10-07.md):
  native coefficients and measures observed with zero updates; first-step row
  expansion reproduced; factor-one block filter is not a nonworsening guard.
  Two bounded arms remain unqualified with original final gates preserved.
- [SLOT IdVd and SimpleMOS alignment controls, 2026-10-07](validation/slot_ldmos_idvd_simplemos_alignment_2026-10-07.md):
  non-Delaunay edges block direct geometry transfer; three native ablations
  complete, while six bounded local physics/numerics arms remain unqualified.
  Frozen edge/source/port checks preserve original gates and ledger identity.
- [SLOT IdVd existing failed-row ratio guards, 2026-10-07](validation/slot_ldmos_idvd_failed_row_ratio_guard_2026-10-07.md):
  real guards detect nine worsening failed rows; five constrained inner solves
  remain unverified despite first-order descent directions. No accepted update
  or qualified point; original physics and final gates unchanged.
- [SLOT IdVd nonlinear SG neighbor guards, 2026-10-07](validation/slot_ldmos_idvd_nonlinear_sg_guard_2026-10-07.md):
  protects 4034/4036 in same-parent controls and accepts two diagnostic updates;
  still 147 failed rows and worsening 4046 ratio from faster scale reduction.
  Original gates retained; no qualified point or scan restart.
- [SLOT IdVd SG/mobility curvature split, 2026-10-07](validation/slot_ldmos_idvd_sg_mobility_curvature_2026-10-07.md):
  fixed-mobility SG factors explain 98.77%/98.15% of local flux curvature at
  4034/4036; mobility and interaction terms separately reconciled to full flux.
  Frozen probes only; no accepted update or qualified point.
- [SLOT IdVd fixed-direction derivative check, 2026-10-07](validation/slot_ldmos_idvd_fixed_direction_derivative_2026-10-07.md):
  central differences on 4034/4036 converge quadratically to Jd within
  2.57/3.03 ppm; flux curvature explains the finite-step prediction failure.
  Frozen diagnostics only, no accepted Newton update or qualified point.
- [SLOT IdVd edge subdomains and prediction margins, 2026-10-07](validation/slot_ldmos_idvd_edge_patch_margin_2026-10-07.md):
  three connected 63-DOF patches and reoptimized smaller steps remain rejected;
  isolates nonlinear prediction errors at 4034/4036 and changing row scale at
  4046. No accepted update, qualified point or production restart.
- [SLOT IdVd carrier and conservation guards, 2026-10-07](validation/slot_ldmos_idvd_carrier_conservation_guard_2026-10-07.md):
  block-balanced directions and real closure checks accept three diagnostic
  updates across two chains; scaled trial directions with nonlinear damage
  are rejected. Still 147 failed rows and no qualified point.
- [SLOT IdVd guarded coupled carrier/potential repair, 2026-10-07](validation/slot_ldmos_idvd_guarded_coupled_repair_2026-10-07.md):
  matched 63-DOF controls prevent new failed rows; four diagnostic updates
  reduce Poisson residual 42-fold but leave 147 failed rows and worsen hole
  closure. No qualified point or production restart.
- [SLOT IdVd density mapping and hotspot initialization, 2026-10-07](validation/slot_ldmos_idvd_density_mapping_hotspot_seed_2026-10-07.md):
  independent Fermi quadrature localizes the mapping discrepancy; 55-node
  population-consistent seeds increase transport but disturb Poisson and
  neighbor rows. No common fixed-state loading scale, no qualified point.
- [SLOT IdVd contact reaction, population representation and native HFS audit, 2026-10-07](validation/slot_ldmos_idvd_contact_reaction_population_hfs_2026-10-07.md):
  zero-update native observer restored and identity verified; contact source
  reactions and SG fluxes separated, saved/reconstructed populations matched.
  No Newton update or qualified point; initialization remains conditional.
- [SLOT IdVd source ledger and HFS formation audit, 2026-10-07](validation/slot_ldmos_idvd_source_ledger_hfs_audit_2026-10-07.md):
  independent edge/source conservation closes internally; native signed-Si
  source integration exposes a contact-accounting question. Native zero-update
  observer times out at startup; conditional loading remains unstarted.
- [SLOT IdVd coupled restoration, source homotopy and native RefDens controls, 2026-10-07](validation/slot_ldmos_idvd_restoration_homotopy_refdens_2026-10-07.md):
  connected 63-variable restoration remains unqualified; verified homotopy
  directions reveal a hotspot population barrier. Two native smoothing arms
  complete and the explicit-zero arm times out; no production model change.
- [SLOT IdVd low-density mechanism and native-control review, 2026-10-07](validation/slot_ldmos_idvd_low_density_mechanism_review_2026-10-07.md):
  symbolic residual-barrier analysis, primary literature/code references, and
  148 re-extracted native node records; version-specific RefDens controls proposed.
  Read-only investigation, no new solve or qualified point.
- [SLOT IdVd local nonlinear trust-region control, 2026-10-07](validation/slot_ldmos_idvd_local_trust_region_2026-10-07.md):
  fresh residuals and Jacobians enable three/two guarded diagnostic steps;
  thirty/twenty-nine failed rows remain. One closure-probe timeout, no qualified
  point or production restart; all original gates unchanged.
- [SLOT IdVd directly constrained local directions, 2026-10-07](validation/slot_ldmos_idvd_constrained_qf_direction_2026-10-07.md):
  feasible active-set directions avoid uniform clipping, but nonlinear row
  checks still fail; directional differences support finite-step curvature.
  Twelve frozen probes, no Newton restart or qualified point.
- [SLOT IdVd bounded local QF repair, 2026-10-07](validation/slot_ldmos_idvd_local_qf_repair_2026-10-07.md):
  six frozen local trials with 30/50 variables; neighboring-row guards avoid
  new failures but do not resolve the original thirty. No Newton restart or
  qualified point; unchanged physics and acceptance gates.
- [SLOT IdVd QF mobility-feedback audit, 2026-10-07](validation/slot_ldmos_idvd_qf_feedback_audit_2026-10-07.md):
  stabilized matrix/direction and finite-difference evidence, frozen single-row
  scan, and a residual-preserving diagnostic linearization; the bounded solve
  remains unqualified, with original gates and ledgers unchanged.
- [SLOT IdVd LLM technical handoff, 2026-10-06](validation/slot_ldmos_idvd_llm_handoff_2026-10-06.md):
  self-contained context, qualification boundaries, source mapping controls,
  failure trajectories and node/material/angle evidence for external analysis.
- [SLOT IdVd contact-field fallback bounded control, 2026-10-06](validation/slot_ldmos_idvd_contact_field_fallback_2026-10-06.md):
  205 contact-adjacent cells audited; source-side holes substantially restored,
  but three bounded solves produce no qualified state. Thirty local rows remain;
  original gates and ledgers unchanged, no long continuation or float80 claim.
- [SLOT IdVd spatial, model-switch and precision controls, 2026-10-06](validation/slot_ldmos_idvd_spatial_ablation_precision_2026-10-06.md):
  four matched drain points through 2.3683165 V; native avalanche/HFS/64-bit
  controls complete. Frozen probes prioritize contact avalanche driving-force
  mapping; no changed-model qualified continuation or Vela float80 claim.
- [SLOT IdVd 60 V secant campaign, 2026-10-06](validation/slot_ldmos_idvd_60v_campaign_2026-10-06.md):
  user-authorized nine-hour local run from Vd=0.310540 V, with twenty-minute
  monitoring; consult live summary and exit evidence for progress, not this launch record.
- [SLOT IdVd opt-in secant CLI integration, 2026-10-06](validation/slot_ldmos_idvd_secant_integration_2026-10-06.md):
  explicit drain predictor and validated history import; two new drain points
  qualify through 0.310540 V in 326 s with unchanged gates. Saved-evidence audit
  passes; no long scan or numerical failure-recovery claim.
- [SLOT IdVd secant neighbors and same-target reproduction, 2026-10-06](validation/slot_ldmos_idvd_secant_neighbors_2026-10-06.md):
  three adjacent drain points qualify through 0.210913 V; alternate secant
  reproduction qualifies, while anchor warm-start reproduction fails the original
  row and closure gates; diagnostic only, no production ledger update.
- [SLOT IdVd warm-start/secant bounded comparison, 2026-10-06](validation/slot_ldmos_idvd_secant_comparison_2026-10-06.md):
  secant directly qualifies Vd=0.1134313 V in 182 s; warm-start exhausts its
  1200 s arm budget after qualifying only the intermediate point. Original
  gates remain unchanged; this is not a paired-endpoint speedup measurement.

- [SLOT IdVd bounded resume, 2026-10-06](validation/slot_ldmos_idvd_resume_2026-10-06.md):
  resumes after the user pause from qualified Vg=2 V / Vd=0 V, with explicit
  parent revalidation, separate imported/new counts and the remaining budget.
  Consult the new campaign summary for current progress.

- [SLOT IdVd gate-initialization integration and bounded full-path run, 2026-10-05](validation/slot_ldmos_idvd_gate_path_integration_2026-10-05.md):
  opt-in Fermi Poly/Poisson initialization is integrated with original strict
  acceptance and retry isolation. New cold and first gate point pass; the
  campaign was paused by the user at Vg=2 V; see the October 6 resume above.

- [SLOT IdVd Fermi Poly initialization validation, 2026-10-05](validation/slot_ldmos_idvd_fermi_poly_initialization_2026-10-05.md):
  rigid Poly bias prediction passes both direct-Newton and Poisson-then-Newton
  qualification at Vg=0.1 V, followed by qualified 0.2/0.3 V points at Vd=0.
  Original gates are unchanged; the full gate/drain path and native agreement
  remain unvalidated. The predictor is explicit and not yet in the default driver.

- [SLOT IdVd combination implementation and cold validation, 2026-10-05](validation/slot_ldmos_idvd_combination_validation_2026-10-05.md):
  explicit Fermi/PhuMob/Lombardi/HFS plus Si-only avalanche and Poly constant
  mobility pass targeted regressions and cold-equilibrium gates. Four bounded
  nonzero-gate attempts fail the original gates; no qualified IdVd current curve
  is available. This element-box profile does not reuse the external AverageBox
  transport geometry, whose compatibility guard remains in force.

- [SLOT IdVd native baseline and compatibility preflight, 2026-10-05](validation/slot_ldmos_idvd_preflight_2026-10-05.md):
  native main-deck Vg=5 V / Vd=0--60 V completes with 93 drain records; Vela
  preflight identified Fermi/PhuMob/avalanche and material-scope compatibility
  blocks. See the newer implementation report above for the subsequent work;
  this preflight produced zero qualified Vela IdVd points.

- [SLOT IDL independent materials and cold-path validation, 2026-10-05](validation/slot_ldmos_idl_cold_material_validation_2026-10-05.md):
  a separate Bennett-Wilson material bundle completes cold equilibrium, Vd=0.1 V
  and Vg=0--5 V; 89 accepted states, unchanged gates, 0.01111% maximum current
  difference at 11 native matched-gate samples. Native geometry remains required.

- [SLOT IDL PolySilicon material attribution, 2026-10-05](validation/slot_ldmos_idl_poly_material_attribution_2026-10-05.md):
  coherent intrinsic density and energy-reference mapping explains the previous
  threshold-region current excess; local interventions precede the cold path above.

- [SLOT IDL Bennett-Wilson continuation, 2026-10-05](validation/slot_ldmos_idl_bennett_continuation_2026-10-05.md):
  Vd=0.1 V and the 0--5 V gate path complete under the original gates; six matched-bias
  samples expose 28--29% subthreshold current differences despite 0.135% agreement at 5 V.

- [SLOT IDL Bennett-Wilson model mapping, 2026-10-05](validation/slot_ldmos_idl_bennett_wilson_mapping_2026-10-05.md):
  explicit model, parameter/unit mapping and related regressions pass; cold zero-bias
  integration passes on the SLOT mesh. Nonzero-drain and gate results are in the continuation report above.

- [SLOT IDL reference and initialization diagnostics, 2026-10-05](validation/slot_ldmos_idl_initial_validation_2026-10-05.md):
  native 0--5 V gate reference completed; runtime Bennett/Wilson BGN differs from the
  provisional Vela none configuration. No matched-physics IDL pass is claimed.

- [SLOT bounded breakdown continuation, 2026-10-04](validation/slot_ldmos_s04_breakdown_goal_2026-10-04.md):
  38 accepted points reach 6.27688e-8 A/um under the original gates; the next point exhausts
  bounded retries. Evidence audit passes, but the 1e-7 current goal remains unmet.

- [SLOT damped adjacent points and same-supply validation, 2026-10-03](validation/slot_ldmos_s03_damped_adjacent_same_supply_2026-10-03.md):
  full recovery followed by half-step polishing passes two target points and one independent
  sample; the third target stops at the original local-patch cap, with all costs retained.

- [SLOT carrier-block damping comparison, 2026-10-03](validation/slot_ldmos_s03_carrier_damping_comparison_2026-10-03.md):
  half steps suppress the alternating correction at both qualified inputs; eight cases
  pass original gates, with actual-supply confirmation kept distinct from same-supply reproduction.

- [SLOT symmetric precision revalidation, 2026-10-03](validation/slot_ldmos_s03_symmetric_precision_revalidation_2026-10-03.md):
  two original-gate repair cycles reveal alternating current differences; own-bias controls
  reproduce the effect, tighter solves fail, and no production points are promoted.

- [SLOT same-parent workflow and binary comparison, 2026-10-03](validation/slot_ldmos_s03_same_parent_workflow_comparison_2026-10-03.md):
  both binaries reproduce the old recovery result; current recovery also passes and
  closes the fixed target, with all attempt costs and remaining path sensitivity reported.

- [SLOT native-settings and row-threshold review, 2026-10-03](validation/slot_ldmos_s03_native_settings_row_review_2026-10-03.md):
  original deck/log and archived successful branch checked; six physical inputs match,
  original row gate previously crossed the criterion, workflow regression prioritized.

- [SLOT recovered-state target correction and independent confirmation, 2026-10-03](validation/slot_ldmos_s03_recovered_target_reproduction_2026-10-03.md):
  both three-map arms pass original device and fixed-supply gates at 1666.4974 V;
  recovery costs retained, 108 regressions pass, diagnostic only and no ledger promotion.

- [SLOT saved Poisson-failure state recovery, 2026-10-03](validation/slot_ldmos_s03_poisson_failure_recovery_2026-10-03.md):
  one carrier correction and two Poisson updates pass original device and actual-supply
  gates; target supply remains 0.642 V away, diagnostic state only, no ledger promotion.

- [SLOT next continuation stopped at pre-carrier Poisson repair, 2026-10-03](validation/slot_ldmos_s03_port_continuation_stop_2026-10-03.md):
  first 10% target has a converged coupled candidate but Poisson line-search
  failure; zero new points, original gates preserved, last qualified parent retained.

- [SLOT port-corrected three-point continuation and same-supply samples, 2026-10-03](validation/slot_ldmos_s03_port_corrected_continuation_2026-10-03.md):
  three original-gate new pilot points and two independent samples pass through
  1.480e-9 A/um; full cost recorded, recovery retained, no long scan started.

- [SLOT same-supply port correction and independent reproduction, 2026-10-03](validation/slot_ldmos_s03_same_supply_reproduction_2026-10-03.md):
  both recovered seeds pass original fixed-supply and device gates after three
  bounded maps each; five coupled exits require recovery, no new bias advance.

- [SLOT carrier and Poisson recovery of independent-seed failures, 2026-10-03](validation/slot_ldmos_s03_carrier_poisson_recovery_2026-10-03.md):
  both saved failures recover through original device and actual-supply gates;
  one carrier-block correction resolves the rows, while same-inner current
  differences remain and no new bias points are added.

- [SLOT predictive-selection switch and independent-seed validation, 2026-10-03](validation/slot_ldmos_s03_switch_independent_validation_2026-10-03.md):
  same-seed on/off both pass three cases; two different-seed reproductions fail
  row and patch-budget gates, so no new bias advance; original criteria preserved.

- [SLOT predictive neighboring-row selection and three-point validation, 2026-10-03](validation/slot_ldmos_s03_predictive_neighbor_validation_2026-10-03.md):
  default-off first-cycle FD-impact expansion; three original-gate consecutive
  points accepted through 1.112e-9 A/um, unchanged production ledger and gates.

- [SLOT patch-boundary propagation and anticipatory repair, 2026-10-02](validation/slot_ldmos_s02_patch_boundary_2026-10-02.md):
  same-seed six/eight-row initial patches recover the remaining 10% input with
  original gates; boundary transfer quantified, diagnostic states only.

- [SLOT bounded neighboring-row repair, 2026-10-02](validation/slot_ldmos_s02_neighbor_row_repair_2026-10-02.md):
  default-off local QF repair recovers two of three failed inputs with unchanged
  gates and actual-supply confirmation; one remains rejected, no scan restart.

- [SLOT row-tolerance, ULP and local-repair diagnostic, 2026-10-02](validation/slot_ldmos_s02_row_tolerance_diagnostic_2026-10-02.md):
  relaxed tolerances merely accept the unchanged seed; a damped two-node repair
  passes original device gates with unchanged evaluated drain current; PETSc/Wolfram
  parameter comparison, no production promotion or scan restart.

- [SLOT adjacent-point, fixed-supply and goal continuation, 2026-10-02](validation/slot_ldmos_s02_sequence_port_goal_2026-10-02.md):
  three original-gate adjacent points pass; a same-seed repair-path control resolves
  the fixed-supply discrepancy to 0.100 microvolt; bounded goal continuation is launched.

- [SLOT six-row repair and merit controls, 2026-10-02](validation/slot_ldmos_s02_six_row_repair_2026-10-02.md):
  one carrier correction followed by Poisson repair passes original device and
  actual-supply gates in both coordinate controls; target-supply reproduction remains separate.

- [SLOT combined-path replay and actual-supply confirmation, 2026-10-02](validation/slot_ldmos_s02_repeat_loadline_2026-10-02.md):
  deterministic replay and original actual-supply confirmation pass; independent
  parent-seeded control still fails six electron rows, with no ledger promotion.

- [SLOT merit trials and poly/Poisson recovery controls, 2026-10-02](validation/slot_ldmos_s02_merit_f3_2026-10-02.md):
  rejected full trial improves Si rows; bounded combined repair passes original
  fixed-inner device checks, while target load-line confirmation remains open.

- [SLOT ULP, fixed-inner and bounded coupled-depth diagnostics, 2026-10-02](validation/slot_ldmos_s02_ulp_f1_f2_2026-10-02.md):
  controlled QF sensitivity, five unqualified fixed-inner trials, and a one-update
  coupled improvement blocked by merit on further correction; original gates retained.

- [SLOT fold current-target control and same-supply controls, 2026-10-01](validation/slot_ldmos_fold_current_target_2026-10-01.md):
  opt-in fold-safe proposal/history, 360 recorded-state decisions, and bounded
  A/B/C controls; all candidates fail local-row qualification, no new accepted point.

- [SLOT full-state Jacobian tangent predictor, 2026-10-01](validation/slot_ldmos_full_state_tangent_2026-10-01.md):
  default-off predictor, unchanged strict gates, four qualified same-parent
  controls; independent same-bias reclosure exhausts its bounded budget.

- [SLOT short continuation after validation reuse, 2026-09-30](validation/slot_ldmos_reuse_short_continuation_2026-09-30.md):
  two six-step paths and two independent controls pass unchanged gates;
  real bounded recovery succeeds, while full-run speedup remains unproven.

- [SLOT validation reuse and bounded candidate correction, 2026-09-30](validation/slot_ldmos_validation_reuse_candidate_correction_2026-09-30.md):
  unified proof adoption, content-checked process-local reuse, 145 regressions,
  and fixed-candidate controls; extra correction does not improve cost.

- [SLOT candidate-current bias and repeated validation cost, 2026-09-30](validation/slot_ldmos_candidate_bias_validation_cost_2026-09-30.md):
  read-only balance decomposition and state/proof comparison identify incomplete
  candidate closure, a repeated validation, and reusable parsing costs; no new solve.

- [SLOT conditional row repair and 10/12.5 mV short batches, 2026-09-30](validation/slot_ldmos_conditional_row_repair_batch_2026-09-30.md):
  default-off one-pass repair, 137 regressions, 12 qualified branch points and
  two independent controls; row quality improves without a throughput gain.

- [SLOT nodes 3595/3603 correction-cost check, 2026-09-30](validation/slot_ldmos_nodes3595_3603_correction_cost_2026-09-30.md):
  two one-pass carrier repairs restore row-quality margin at unchanged bias and
  gates; measured full diagnostic overhead motivates conditional integration.

- [SLOT point-6 reclosure and step-quality gate, 2026-09-30](validation/slot_ldmos_point6_reclosure_step_gate_2026-09-30.md):
  bounded reclosure removes the repeatability gap; 8/10 mV same-parent trials pass
  final gates, but the 10 mV hole-row quality prevents further step expansion.

- [SLOT matched-feedback same-parent short batch, 2026-09-30](validation/slot_ldmos_matched_feedback_batch_2026-09-30.md):
  six accepted steps per mode on T470p, 11.94% higher local supply advance per
  second, working candidate-quality backoff and a same-bias repeatability caveat.

- [SLOT matched-supply feedback implementation, 2026-09-30](validation/slot_ldmos_matched_feedback_fix_2026-09-30.md):
  opt-in frozen-model evaluation at accepted supply, conservative fallback,
  identity binding, 119 regressions and 18-point replay; no new device solve.

- [SLOT prediction feedback coordinate review, 2026-09-30](validation/slot_ldmos_prediction_feedback_review_2026-09-30.md):
  causal replay of 18 accepted steps separates actual-supply migration from
  same-supply model error; production policy and final qualification unchanged.

- [SLOT finite step expansion from 5 to 10 mV, 2026-09-30](validation/slot_ldmos_finite_step_expansion_2026-09-30.md):
  three qualified same-parent trials, six strict advances and two independent
  controls; port prediction feedback limits sustained step growth.

- [SLOT point-3 bounded same-bias reclosure, 2026-09-30](validation/slot_ldmos_point3_bounded_reclosure_2026-09-30.md):
  ordinary restart performed no update; one carrier correction followed by a full
  Newton update removed the current discrepancy under unchanged final gates.

- [SLOT current-scaled stepping and bounded growth validation, 2026-09-30](validation/slot_ldmos_current_growth_batch_2026-09-30.md):
  four same-parent step trials, 12 qualified advances with quality-triggered
  backoff, two independent controls, and a residual-related repeatability caveat.

- [SLOT daily report and figures, Sept 29–30 morning](validation/slot_ldmos_daily_report_2026-09-30.md):
  management-ready summary, actual device cross-section, 12-point continuity and
  fixed-parent step comparison; qualified voltage is not a BVDS result.

- [SLOT CLI integration and short continuation, 2026-09-30](validation/slot_ldmos_short_cli_batch_2026-09-30.md):
  12 qualified T470p advances, two independent controls, three same-parent step
  controls, 47 regressions and downloaded evidence independently revalidated.

- [SLOT higher-bias candidate conservation and repair cost, 2026-09-29](validation/slot_ldmos_candidate_conservation_high_parents_2026-09-29.md):
  three local parent advances, independent same-bias controls, complete qualification,
  and measured candidate/strict-reclosure costs.

- [SLOT frozen candidate conservation, 2026-09-29](validation/slot_ldmos_candidate_conservation_2026-09-29.md):
  contact reference/increment preservation, independent candidate closure,
  unchanged strict reclosure and bounded local verification.

- [SLOT progress/predictor repair and bounded diagnostics, 2026-09-29](validation/slot_ldmos_progress_predictor_diagnostics_2026-09-29.md):
  opt-in stagnation guard, matched predictor voltage, same-state conservation
  diagnostics, strict local step controls, and remaining candidate-current mismatch.

- [SLOT zero-device-RHS option and three-parent validation, 2026-09-29](validation/slot_ldmos_zero_rhs_multi_parent_2026-09-29.md):
  default-off candidate option, unchanged strict final gates, three additional
  T470p controls, and revalidated 35 V / latest high-voltage restart windows.

- [SLOT candidate forcing comparison on T470p, 2026-09-29](validation/slot_ldmos_candidate_forcing_comparison_2026-09-29.md):
  same-parent strict/inexact controls, unchanged final qualification, successful
  zero-device-RHS candidate at the original tolerance, and independent reclosure.

- [SLOT local 35 V matched-bias review and high-voltage continuation, 2026-09-28](validation/slot_ldmos_home_35v_to_bv_2026-09-28.md):
  matched current/field/source comparison, independently requalified origin-36
  parent window and local continuation toward the original current criterion.

- [SLOT bounded recovery and independently audited continuation through 35 V, 2026-09-28](validation/slot_ldmos_bounded_recovery_2026-09-28.md):
  unchanged qualification gates, actual-supply circuit confirmation, reference-frame
  control and 44 accepted T470p points from 31.3 to 35 V external bias.

- [SLOT hybrid prediction through 30 V on T470p, 2026-09-28](validation/slot_ldmos_hybrid_continuation_30v_2026-09-28.md):
  28 new resistor points and one frame requalification, actual rational-candidate
  use, independent audits and quantified post-frame step-recovery cost.

- [SLOT hybrid port predictor and T470p pilot, 2026-09-28](validation/slot_ldmos_hybrid_port_predictor_2026-09-28.md):
  opt-in limited endpoint baseline, causal rational-candidate selection, 87
  regressions and independently audited 26.77-to-27.2 V forward validation.

- [SLOT port step controller and offline prediction, 2026-09-27](validation/slot_ldmos_port_step_prediction_2026-09-27.md):
  opt-in independent outer, accepted-inner and trial-inner caps, qualified
  three-point slope fitting, cost/error feedback and read-only evidence replay.

- [SLOT native-goal adaptive BV campaign, 2026-09-27](validation/slot_ldmos_native_goal_campaign_2026-09-27.md):
  launched from the qualified 3 V resistor state toward 200000 V external or
  the 1e-7 A/um current criterion, with qualified-parent prediction and adaptive steps.

- [SLOT predicted resistor continuation through 3 V, 2026-09-27](validation/slot_ldmos_predictor_continuation_3v_2026-09-27.md):
  conservative adaptive continuation from the qualified 2 V state, observed
  first-step recovery, unchanged gates and full attempt-cost accounting.

- [SLOT state prediction and adaptive stepping, 2026-09-27](validation/slot_ldmos_predictor_adaptive_validation_2026-09-27.md):
  opt-in qualified-parent prediction and success-based growth, unchanged fallback
  and acceptance gates, with separate fixed-path and adaptive controls through 2 V.

- [SLOT resistor continuation through 2 V, 2026-09-27](validation/slot_ldmos_low_voltage_continuation_2026-09-27.md):
  serial continuation from the recovered 0.8 V state, unchanged device/circuit
  gates, and bounded Windows checkpoint replacement recovery.

- [SLOT candidate failure classification and 0.8 V recovery, 2026-09-27](validation/slot_ldmos_candidate_failure_recovery_2026-09-27.md):
  bound numerical-failure protocol, original-state replay and qualified 0.7–0.8 V
  resistor recovery with unchanged device and circuit gates; full 0-to-BV scan remains open.

- [SLOT full resistor campaign, 2026-09-26](validation/slot_ldmos_full_resistor_campaign_2026-09-26.md):
  in-progress local equilibrium-to-BV campaign, serial execution and checkpoint recovery;
  this work has not yet qualified a complete new-configuration 0-to-BV path.

- [SLOT native transport resistor branch and BVDS comparison, 2026-09-26](validation/slot_ldmos_native_transport_bv_branch_2026-09-26.md):
  qualified 38.2/38.4 V controls, audited resistor branch across the current criterion,
  independent threshold refinement and interpolation comparison with Sentaurus.

- [SLOT native transport 38 V reclosure and continuous path, 2026-09-26](validation/slot_ldmos_native_transport_38v_continuation_2026-09-26.md):
  independent SC/IIC controls, qualified 25-point 35–38 V continuation,
  matching endpoints and external transport resource checks; the subsequent BVDS
  criterion comparison is linked above.

- [SLOT full native transport geometry and qualified 35 V reclosure, 2026-09-26](validation/slot_ldmos_native_transport_35v_2026-09-26.md):
  shared opt-in source guard, assembled/Jacobian checks and independent IIC/SC
  reclosure; the subsequent 38 V and continuous-path evidence is linked above.

- [SLOT near-surface hole source and transport localization, 2026-09-26](validation/slot_ldmos_near_surface_hole_transport_2026-09-26.md):
  same-state source replay, density/field/current factor budget, and frozen
  edge-weight plus mobility-measure control; independent reclosure remains pending.

- [SLOT circumcentric continuation through 38 V, 2026-09-26](validation/slot_ldmos_circumcentric_continuation_38v_2026-09-26.md):
  qualified continuous 35–38 V path, bounded density/QF candidate recovery,
  and remaining local hole-source versus non-Delaunay measure diagnostics.

- [SLOT circumcentric volume localization and 35–38 V controls, 2026-09-26](validation/slot_ldmos_circumcentric_volume_35_38v_2026-09-26.md):
  opt-in charge/recombination measures, independent 35 V gates and high-voltage continuation evidence.

- [SLOT per-cell dielectric operator and independent 35 V reclosure, 2026-09-26](validation/slot_ldmos_per_cell_reclosure_35v_2026-09-26.md):
  opt-in operator, full-device replica acceptance and qualified IIC/SC points;
  source/current improvement with remaining bulk-potential and IIC-source differences.

- [SLOT 35 V IIC nonlinearity and H2 frozen controls, 2026-09-26](validation/slot_ldmos_iic_nonlinearity_h2_35v_2026-09-26.md):
  state-path cancellation checks and charge-volume/permittivity controls;
  promising joint variant, replica-gate limitations, no production policy change.

- [SLOT 35 V adjoint localization and H1 control, 2026-09-26](validation/slot_ldmos_iic_loop_gain_localization_35v_2026-09-26.md):
  verified linear probes, two-segment diagnostic attribution, unsuccessful H1 promotion
  and the remaining IIC linearization discrepancy; no production policy change.

- [SLOT IIC source budget at 35/37/38 V, 2026-09-25](validation/slot_ldmos_iic_source_budget_2026-09-25.md):
  integration-support checks, frozen production-kernel replay, independently
  qualified no-feedback states and a self-consistent source-amplification budget.

- [SLOT 38.2/38.4 V comparison and high-resistance branch, 2026-09-25](validation/slot_ldmos_near_breakdown_branch_2026-09-25.md):
  qualified direct continuation, native field/source comparison and independently
  audited resistor points; the local breakdown-current criterion is not yet crossed.

- [SLOT local-origin continuation through 38 V and resistor start, 2026-09-25](validation/slot_ldmos_local_origin_38v_loadline_2026-09-25.md):
  eight newly qualified voltages, native 37.5/38 V comparisons, and an independently
  audited original-resistor start plus a 0.01 V outer-voltage advance.

- [SLOT local-origin continuation through 37 V, 2026-09-25](validation/slot_ldmos_local_origin_37v_2026-09-25.md):
  57 newly qualified voltages, 63 audited records, native 35/37 V field/source
  comparisons, origin 35 V restart, and tighter-floor checkpoint revalidation.

- [SLOT local-origin continuation through 30 V, 2026-09-25](validation/slot_ldmos_local_origin_30v_2026-09-25.md):
  81 newly qualified voltages, 87 audited records, serial native 25 V reference,
  25/30 V field/source comparisons, origin 15 V restart and residual-floor limits.

- [SLOT local-origin continuation through 20 V, 2026-09-24](validation/slot_ldmos_local_origin_20v_2026-09-24.md):
  48 newly qualified voltages, 52 audited records, guarded Poly restart support,
  native 15/20 V current, field and avalanche-source comparisons, and restart anchors.

- [SLOT local-origin continuation through 10 V, 2026-09-24](validation/slot_ldmos_local_origin_high_voltage_2026-09-24.md):
  36 newly qualified voltages, 40 audited qualification records, serial Release
  continuation in origin 5 V, native 5/10 V comparisons and final restart entry.

- [SLOT reference-frame continuation tools, 2026-09-24](validation/slot_ldmos_reference_frame_tools_2026-09-24.md):
  explicit frame conversion, inactive contact constraints, requalification of both
  predictor anchors, physical-voltage reporting and Release low-bias integration.

- [SLOT 10 V electrical-frame audit, 2026-09-24](validation/slot_ldmos_10v_frame_audit_2026-09-24.md):
  shifted-frame qualification, original-frame revalidation failure, contact-floor
  semantics and inactive contact projection; not continuous cold-origin 10 V or BVDS completion.

- [SLOT local-origin low-voltage completion, 2026-09-24](validation/slot_ldmos_local_origin_low_voltage_2026-09-24.md):
  36 qualified points from local equilibrium to 2.125 V including restart recovery,
  six-point path comparison, serial 92-test tool verification; high-voltage/BVDS remains open.

- [SLOT HDF5 continuation and native low-bias comparison, 2026-09-24](validation/slot_ldmos_hdf5_continuation_2026-09-24.md):
  three-driver state migration, block-probe output fix, Release regression,
  checkpoint confirmation, repeat-chain progress and native 2.025 V evidence.

- [SLOT status and handoff, 2026-09-24](validation/slot_ldmos_status_and_handoff_2026-09-24.md):
  original-deck function and parameter matrix, completed versus unfinished
  validation, current HDF5 anchors and next-stage entry; separates projected-SG
  high-voltage results from the controlled new-source 2.125 V qualification.

- [SLOT local main integration, 2026-09-23](validation/slot_ldmos_main_integration_2026-09-23.md):
  HDF5 restart compatibility, default Jacobian structure reuse, retained SLOT
  physics fixes and bounded Release revalidation; continuation driver migration
  remains separate from high-voltage and BVDS qualification.

- [SLOT mixed-contact guard and escape controls, 2026-09-23](validation/slot_ldmos_contact_guard_and_escape_controls_2026-09-23.md):
  local contact-majority classification fix, T1 update traces, HFS-off and
  single-node residual controls, ULP sensitivity limits and unchanged-gate
  device revalidation. A 0.025 V continuation now qualifies new-source points
  through 2.125 V; the 10 V observed-source checkpoint remains rejected.

- [SLOT comparison analysis for LLM review, 2026-09-23](validation/slot_ldmos_comparison_analysis_for_llm_2026-09-23.md):
  Chinese evidence summary, separate projected-SG and observed-source qualification,
  qualified new-source 2 V, failed 2.125/10 V controls, Jacobian coverage limits,
  units, outstanding hypotheses and concrete questions for independent analysis.

- [SLOT observed element source integration, 2026-09-23](validation/slot_ldmos_observed_source_integration_2026-09-23.md):
  complete fixed-mesh native selector coverage, independent screening-floor
  discrepancy evidence, controlled live source/Jacobian integration and explicit
  qualified 2 V and failed 2.125/10 V coupled controls. Strict original BVDS remains unqualified.

- [SLOT contact avalanche kernel fix, 2026-09-23](validation/slot_ldmos_contact_kernel_fix_2026-09-23.md):
  contact/interface full-gradient correction, analytic regression, full local
  nine-column live-mobility Jacobian checks and native 10/37 V replay.
  Controlled coupled assembly and 2 V qualification are recorded above;
  strict BVDS remains unqualified.

- [SLOT PMI callback follow-up, 2026-09-23](validation/slot_ldmos_pmi_callback_followup_2026-09-23.md):
  complete callback-to-element matching, contact full-gradient current evidence,
  effective current scale, saved-state precision limits and qualified frozen
  AvalDensGradQF-off control. The contact-kernel correction is now verified in
  the contact-kernel report; subsequent controlled assembly is recorded above.

- [SLOT native avalanche identification, 2026-09-23](validation/slot_ldmos_native_avalanche_identification_2026-09-23.md):
  compiler-independent frozen controls, extrema harmonic density, cosine alpha
  smoothing, best-vertex mobility and 10/37 V interior replay with explicit
  geometry-choice uncertainty, native truncated source volumes and a standalone
  nine-column Jacobian kernel. GCC 11 PMI loading and frozen-output parity now
  pass, with direct near-tie vertex-choice evidence. Subsequent fixed-mesh selector
  completion and controlled coupled integration are recorded above.
  Original AvalDensGradQF and complete BVDS remain unqualified.

- [SLOT load-line driver review, 2026-09-23](validation/slot_ldmos_loadline_review_2026-09-23.md):
  guarded independent-device continuation, original all-row/circuit gates,
  reviewed restart identity and budgets, and real 37.075/37.1 V replay evidence.
  Original AvalDensGradQF and complete BVDS remain separate.

- [SLOT completed internal continuation, 2026-09-22](validation/slot_ldmos_continuation_completion_2026-09-22.md):
  both resumed segments complete through 20/37 V, 49 qualified new points,
  and 20/30/35/37 V field/current overlap with independent checkpoints.
  Cold-origin replay and original avalanche/BVDS validation remain separate.

- [SLOT commit handoff and next tasks, 2026-09-22](validation/slot_ldmos_next_tasks_2026-09-22.md):
  committed fixes, a timestamped live-run snapshot, and concrete acceptance
  criteria for continuation, load-line automation, original avalanche physics,
  error budgeting, and final BVDS qualification.

- [SLOT high-voltage continuation, 2026-09-22](validation/slot_ldmos_high_voltage_resume_2026-09-22.md):
  independently qualified projected-SG checkpoints through 37 V, complete
  Jacobian checks, guarded Poly restarts, compensated-QF restart repairs, and
  qualified 37.02/37.05 V external-load-line controls. Original AvalDensGradQF source
  averaging and full BVDS qualification remain unresolved.

- [SLOT energy reference and source volume, 2026-09-21](validation/slot_ldmos_energy_reference_and_source_volume_2026-09-21.md):
  band-derived poly reference, material-local SRH/Auger integration, qualified
  2 V error budget and independent local-equilibrium continuation status.

- [SLOT Auger and poly independent controls, 2026-09-21](validation/slot_ldmos_auger_poly_independent_controls_2026-09-21.md):
  density-enhanced Auger controls at 2/5/10 V, isolated poly transport and
  material-scoped SRH/Auger, qualified 2 V restoration and remaining energy-reference gap.

- [SLOT element IALMob execution, 2026-09-21](validation/slot_ldmos_element_ialmob_execution_2026-09-21.md):
  original isotropic parameter and element replay, strict 2 V checkpoint,
  single-factor error budget and complete mobility/avalanche Jacobian work.
  High-voltage and external-resistor qualification are recorded separately.

- [SimpleMOS raw states/logs archive and local main migration](../reference_tcad/simplemos_sentaurus2022/raw_assets_20261001/README.md):
  complete local H5/TDR coverage, 60,954 assets including structured logs and existing bundles,
  persistent VM storage, restore tool and source-to-main commit mapping.

- [SimpleMOS historical evidence archive and recovery](../reference_tcad/simplemos_sentaurus2022/local_evidence_20260930/README.md):
  1,424 small records in Git; complete 1,946-file archive on the Sentaurus VM,
  with exact hashes and recovery of 522 detailed historical artifacts to their original paths.

- [Latest D0 Release closeout, 2026-09-25](validation/templates_ldmos_release_closeout_2026-09-25.md):
  Codespaces Release 858/858 CTest, 590 exact points, 698 identical state pairs,
  62 native local fields and six checkpoints passed. Three paired rounds reduce
  mean dual-gate wall time by 1.10%; all 6,390 archive members verified remotely,
  with the small verification report retrieved locally at the user's request.

- [Linked defaults and electrothermal cross-point linear reuse, 2026-09-25](validation/templates_ldmos_thermal_linear_reuse_2026-09-25.md):
  unified production defaults, exact-pattern reuse and failure/restart invalidation;
  UCRT64 Debug 858-test regression passed; subsequent full-curve timing is recorded
  in the latest D0 Release closeout above.

- [Shared Halley default integration, 2026-09-25](validation/templates_ldmos_halley_default_2026-09-25.md):
  C++ and isothermal/electrothermal JSON defaults unified, explicit legacy retained.
  Codespaces Release 858/858 CTest and seven full D0/G3/D5/D4 curves (217 points)
  passed. Default D0 states match the qualified explicit Halley run exactly.

- [Screening-root full D0 qualification and repeats, 2026-09-25](validation/templates_ldmos_screening_full_2026-09-25.md):
  20 full curves / 620 exact points, 10 native joint checks and 744 state pairs
  passed. Halley reduced aggregate wall time by 6.70% over three fresh
  interleaved dual-gate rounds (5.68%-7.24% per round); all six gate pairs were
  faster. Newton work is unchanged, with small assembly-count differences.
  Evidence was retrieved; default was legacy in that experiment. The subsequent
  shared-default integration above promotes Halley after full regression.

- [Screening-root candidates on Codespaces, 2026-09-24](validation/templates_ldmos_screening_candidates_2026-09-24.md):
  local main merged; Release 856/856 CTest and eight short D0 curves passed.
  Safeguarded Newton/Halley/TOMS748 retain states and Newton work; one-round
  paired wall reductions are 6.36%/7.18%/6.32%. Legacy was the default then;
  full curves and repeated end-to-end qualification are tracked in the newer
  2026-09-25 report above.

- [Codespaces build environment](codespaces.md): container dependencies and
  PowerShell/Codex commands for remote builds, tests, and result retrieval.

- [T470p electrothermal subcost measurements, 2026-09-24](validation/templates_ldmos_thermal_subcost_2026-09-24.md):
  eight full D0 curves and four joint checks passed with identical states/work;
  screening roots dominate the three measured candidates. Thermal scatter
  requires explicit clock-bias qualification; diagnostic timers default off.

- [T470p electrothermal gprof, 2026-09-24](validation/templates_ldmos_thermal_gprof_2026-09-24.md):
  latest default structure reuse, full D0 Release/gprof controls and original gates;
  assembly/IALMob and factorization dominate. Four-equation linear-object reuse
  across points remains distinct from the completed assembly cache.

- [Default electrothermal assembly structure reuse, 2026-09-24](validation/templates_ldmos_thermal_structure_2026-09-24.md):
  four-equation/thermal structure reuse; C560/T470p full regression, D0 curves,
  three paired repetitions and omitted-option checks passed. Default is on;
  measured aggregate timing improvement is small (0.70%).

- [Default DC assembly structure reuse, 2026-09-23](validation/templates_ldmos_structure_default_2026-09-23.md):
  production option, invalidation boundaries and T470p default-path regression.
- [Restart state storage](state_archive.md): current `vela.state/2` fields, units,
  mesh/source identity, atomic checkpoint writes and strict recovery semantics.
- [HDF5 production state migration plan, 2026-09-23](validation/templates_ldmos_hdf5_production_migration_plan_2026-09-23.md):
  authorized scope, state schema, restart transactions, active-fixture migration,
  electrothermal coverage and acceptance gates. See the execution record for final results.
- [HDF5 migration execution, 2026-09-23](validation/templates_ldmos_hdf5_migration_execution_2026-09-23.md):
  T470p Release validation; 850 CTest entries, G3/D4/D0 and frozen-300 K recovery
  passed. Final three-round D5 overhead is 0.4572%, below 3%; C560 Linux CTest passed.
- [HDF5 Linux CI, 2026-09-23](validation/templates_ldmos_linux_ci_2026-09-23.md):
  isolated Ubuntu 26.04/GCC 15 Debug build and 850/850 CTest passed on C560 with
  SparseLU; actual dependency detection and GitHub Actions differences are recorded.

- [LDMOS worktree cleanup, 2026-09-22](validation/templates_ldmos_worktree_cleanup_2026-09-22.md):
  duplicate artifacts and explicitly authorized historical JSON/CSV outputs removed;
  118.68 GB total logical file bytes released. Production dependencies and latest HDF5
  evidence retained; old per-point output paths may no longer be replayable.

- [TDR and Vela HDF5 schema compatibility, 2026-09-22](validation/templates_ldmos_tdr_hdf5_compatibility_2026-09-22.md):
  actual files share the HDF5 container but require different schemas; importer tests
  and bidirectional schema-rejection checks pass. No native TDR write compatibility claimed.

- [T470p HDF5 full-curve repeated pairing, 2026-09-22](validation/templates_ldmos_hdf5_full_repeats_2026-09-22.md):
  three alternating serial rounds complete: twelve full D5 curves / 372 points, six joint
  gates and all-field repeated audits pass. HDF5 adds 0.8695% total wall time versus VDS1,
  meeting the approved 3% criterion; unchanged Release program, numerical gates and defaults.

- [Public DD state formats on T470p, 2026-09-22](validation/templates_ldmos_public_state_formats_2026-09-22.md):
  optional HDF5/HighFive, FlatBuffers and structured NPY interfaces pass local/remote
  checks and twelve paired short curves / 96 exact points; all-field equality confirmed.
  One-round combined wall costs are +0.297%/+0.304%/+0.672% versus matched VDS1.
  No speed ranking or default promotion; HDF5 remains the preferred public-storage candidate.

- [T470p binary restart state and reduced seed writes, 2026-09-22](validation/templates_ldmos_binary_state_2026-09-22.md):
  eight short curves / 64 exact points pass unchanged gates and identical states
  and trajectories. VDS1 saves 12.11%/10.95% versus CSV in one paired round.
  Separate memory predictor seeds remove 60 files per curve but add 11.78%/12.52%
  wall time including transport and reconstruction audit; not promoted. Defaults unchanged.

- [T470p D5 file reuse and cost-per-volt controls, 2026-09-22](validation/templates_ldmos_d5_file_efficiency_2026-09-22.md):
  eight short curves / 64 exact points pass unchanged gates and cold final audits.
  Bounded CSV/hash reuse saves 3.00%/2.61% in one paired round with identical trajectories.
  Separate efficiency control saves six services per curve but increases combined wall time
  0.87%; cost-regression feedback never triggers. No default promotion or full-curve claim.

- [T470p D5 cost-feedback continuation and parent timing, 2026-09-22](validation/templates_ldmos_d5_cost_continuation_2026-09-22.md):
  dual-gate first-eight-point controls pass original state gates, but wall time grows
  10.57%/10.21% and updates 367/349 to 386/379. Candidate not promoted to full curves.
  Nonoverlapping parent timers locate 26.1/26.3 s outside worker calls in the controls;
  budget-per-step feedback does not minimize total continuation cost. Defaults unchanged.

- [T470p D5 assembly/Fermi reuse candidates, 2026-09-22](validation/templates_ldmos_d5_kernel_candidates_2026-09-22.md):
  immutable pattern/scatter reuse passes dual-gate first-eight-point controls;
  one paired round reduces wall time by 5.66%/5.90% with identical states and trajectories.
  Fermi endpoint reuse reduces evaluations without a clear timing gain. Combined
  62-point full curves pass original gates and exact state/trajectory checks;
  full-curve qualification is not paired timing. Defaults unchanged.

- [T470p D5 completed gprof results, 2026-09-22](validation/templates_ldmos_d5_gprof_results_2026-09-22.md):
  62 exact points pass unchanged curve/joint gates, with identical states, trajectories
  and recorded counters to ordinary Release. Cross-checks Fermi edge assembly, repeated
  pattern/scatter construction and point-service overhead; proposes progress/cost-based
  continuation without v3's fixed voltage split. No production-policy changes.

- [T470p D5 baseline gprof and adaptive continuation, 2026-09-21](validation/templates_ldmos_d5_gprof_t470p_2026-09-21.md):
  historical pause record and isolated Release/frozen DLL profiling contract;
  resumed and completed on September 22 (see results above).
  Separates v3's case-specific voltage split from a proposed
  progress/cost-driven adaptive controller; no production changes.

- [T470p D5 protected 0.8 V full curves, 2026-09-21](validation/templates_ldmos_d5_step08_full_2026-09-21.md):
  four full curves and fresh paired states pass unchanged gates. One paired round
  reduces wall time by 9.70%/29.06%, but Vg4 Newton updates increase 5.34% due to
  costly recovery points. Low-voltage and no-history safeguards are required;
  two earlier failures remain separate evidence. No production-default promotion.

- [T470p D5 targeted step/prediction/row controls, 2026-09-21](validation/templates_ldmos_d5_targeted_controls_2026-09-21.md):
  20 bounded paths and 12 node-trace controls qualify; 0.8 V steps with fixed 0.2 V QF
  cap reduce short-window work. Larger jumps can trigger internal recovery. Diagnoses
  one low-voltage prediction counterexample and corrects legacy trace-density units
  in analysis; no full-curve or production-default promotion.

- [T470p D5 Newton and prediction localization, 2026-09-21](validation/templates_ldmos_d5_newton_localization_2026-09-21.md):
  six repeated UMFPACK traces: 251 drain advances versus native 49; average updates per
  advance comparable. Sixteen fixed-state prediction controls pass original gates;
  one low-voltage prediction counterexample identified. No production default changes.

- [LDMOS performance closeout and production defaults, 2026-09-20](validation/templates_ldmos_performance_closeout_2026-09-20.md):
  cross-point main-Newton analysis reuse enabled by default; UMFPACK preferred,
  SparseLU fallback, then STRUMPACK/MUMPS/SuperLU_MT in the available-backend list.
  Portable D5 inputs and dual-gate validation entry; historical experiments remain explicit controls.

- [T470p 1/2/4-thread pilot, 2026-09-20](validation/templates_ldmos_t470p_threads_2026-09-20.md):
  39 matrix replays and 13 Vg8 first-eight-point D5 curves pass; 104 exact points,
  unchanged gates. STRUMPACK 2/1 takes 167.445 s versus 230.187 s at 1/1;
  UMFPACK/MUMPS do not benefit under tested settings. One curve round, defaults unchanged.

- [T470p identical-binary host pilot, 2026-09-20](validation/templates_ldmos_t470p_host_pilot_2026-09-20.md):
  five single-thread reuse backends, Vg8 first eight D5 points on N150 and T470p;
  all 80 points pass, remote wall times 12.79–34.48% lower in this single round.
  Full-curve repeats remain paused; host timings must not be pooled across machines.

- [Single-thread five-backend repeats, 2026-09-20](validation/templates_ldmos_single_thread_repeats_2026-09-20.md):
  STRUMPACK dual-gate qualification first, then all five reuse backends with solver/BLAS
  1+1, three fresh D5 rounds (30 curves); original gates and historical timings preserved.

- [Full-curve backend timing execution, 2026-09-20](validation/templates_ldmos_full_curve_timing_execution_2026-09-20.md):
  all 12 D5 first-round curves and six dual-gate qualifications pass; stopped for shortlist
  review. Future runs use solver/BLAS 1+1 for every backend; historical STRUMPACK 2+1
  results stay separate. Repeats and D4 pending; no production default change.

- [Full-curve backend timing plan, 2026-09-20](validation/templates_ldmos_full_curve_timing_plan_2026-09-20.md):
  revised future thread contract fixes all five solvers and BLAS to one thread.
  Original schedule retained as provenance; thread changes require a new frozen batch,
  not a resume that mixes old STRUMPACK two-thread timings with single-thread repeats.

- [Four-backend cross-point curve validation, 2026-09-19](validation/templates_ldmos_four_backend_reuse_2026-09-19.md):
  SparseLU/UMFPACK/MUMPS/SuperLU_MT, ordinary worker versus analysis reuse,
  all 16 dual-gate D5 first-eight-point curves pass with identical paired states,
  Newton updates and factorization counts. One paired round, solver/BLAS 1+1;
  default remains off, auxiliary analyses remain separate.

- [Cross-point Newton analysis reuse, 2026-09-19](validation/templates_ldmos_cross_point_analysis_2026-09-19.md):
  default-off sequential linear context; all six D5 first-eight-point controls pass.
  Incremental wall savings over ordinary worker are 8.57%/8.78%; main-only services
  add no analyses after reuse. Auxiliary solvers remain separate; one round only.

- [STRUMPACK ordering comparison, 2026-09-19](validation/templates_ldmos_strumpack_ordering_2026-09-19.md):
  METIS/AMD/MMD/AND over 36 captured systems, three analysis policies and three rounds;
  all 1,296 systems pass. METIS wins with reuse, AMD wins with per-matrix rebuilding.
  Ordering/tree timings and fill recorded; factor structure varies across some repeats.
  No production or curve qualification change.

- [STRUMPACK verified 2+1 threads: D5 first-eight-point controls, 2026-09-19](validation/templates_ldmos_strumpack_verified_threads_2026-09-19.md):
  all four dual-gate curves and runtime thread audits pass; end-to-end wall time is
  2.87%/4.60% higher than paired UMFPACK. Faster factorization is offset by analysis
  and solve costs; one paired round only, no full-curve qualification or promotion.

- [BLAS threading and STRUMPACK compression, 2026-09-19](validation/templates_ldmos_blas_compression_2026-09-19.md):
  isolated Eigen BLAS and compression candidates; runtime thread verification and
  13 configurations over 36 matrices, three rounds, all accuracy gates pass.
  Corrects the previous BLAS=1 assumption for OpenMP OpenBLAS; no production promotion.

- [Sparse backend and METIS execution, 2026-09-18](validation/templates_ldmos_sparse_backend_execution_2026-09-18.md):
  optional adapters and linear-contract tests implemented; historical/current matrix screening passes.
  All six single-thread D5 dual-gate first-eight-point controls pass; no new wall-time winner.
  STRUMPACK 2/4-thread curve controls also pass but are slower than matched UMFPACK controls;
  SuperLU_MT multithreaded replay failures retained. No winner triggered the full-curve repeat stage.
  Screening complete; default unchanged.

- [Sparse backend and METIS comparison plan, 2026-09-18](validation/templates_ldmos_sparse_backend_plan_2026-09-18.md):
  UCRT64 MUMPS/SuperLU_MT/STRUMPACK/METIS inventory and five smoke checks pass;
  staged integration and fixed-matrix/curve plan only, no new production backend qualification.

- [D5 diagnostics-disabled repeats and strategy controls, 2026-09-18](validation/templates_ldmos_umfpack_low_noise_2026-09-18.md):
  all 8 curves and original joint gates pass; four wall-time gains of 10.38%–17.93%.
  Repeated states/counters match; strict idle conditions were unavailable.
  Historical frozen-matrix controls distinguish ordering, scaling and pivot effects.

- [Isothermal full-curve backend repeats, 2026-09-18](validation/templates_ldmos_umfpack_full_repeats_2026-09-18.md):
  all 16 D5/D4 curves and original joint gates pass; repeat states are identical.
  Factor fill/ordering/work/memory quantified; D4 wall-time gains repeat, D5 timings vary.

- [Isothermal SparseLU/UMFPACK controls, 2026-09-18](validation/templates_ldmos_isothermal_umfpack_2026-09-18.md):
  optional shared backend with symbolic/numeric reuse and failure checks;
  fixed-system replay and first-eight-point controls, qualification tracked separately from R11.

- [LDMOS group review, September 12–17, 2026](reports/ldmos_review_2026-09-12_17/report.md):
  meeting report with device/model parameters, equations, solver flow and pseudocode,
  11 source-backed figures, R11 curve/local-field results and paired performance evidence.

- [Architecture](architecture.md): source tree map, solver paths, and supported
  implementation boundaries.
- [Config schema](config_schema.md): implementation-aligned JSON field
  reference for Poisson, DC sweeps, Newton, unit-scaling input mode, contacts,
  boundaries, solver options, and regression blocks.
- [Sentaurus import](sentaurus_import.md): HDF5/TDR import prerequisites,
  `sentaurus_import` CLI usage, and end-to-end conversion workflow.
- [Sentaurus VMware SSH workflow](sentaurus_vm_ssh_workflow.md): Host-only
  VMware networking, legacy CentOS SSH setup, and remote Sentaurus run/copy
  commands for the local Sentaurus 2018 VM.
- [Poisson unit-scaling notes](development_poisson_unit_scaling.md): developer
  notes for the scaled Poisson assembly path used by
  `scaling.mode = "unit_scaling"`.
- [PN2D BV validation](validation/pn2d_bv_validation.md): current qualified
  template policy, validation gates, limitations, and evidence map.
- [SLOT-LDMOS BVDS validation](validation/slot_ldmos_bvds_validation.md):
  local and global qualification profiles, current low-voltage evidence,
  and remaining physical-model differences.
- [SLOT/main integration and template comparison, 2026-09-21](validation/slot_ldmos_main_integration_2026-09-21.md):
  retained SLOT compatibility, imported template capabilities, and separate
  numerical qualification boundaries.

- [SimpleMOS engineering reference](../reference_tcad/simplemos_sentaurus2022/engineering/README.md):
  qualified configuration, run entry, external input requirements, 816-point
  evidence and a portable read-only hash/coverage audit.
- [SimpleMOS branch status](validation/simplemos_branch_status.md): explicit
  production options, diagnostic boundaries, local evidence dependencies, and
  the current n23 Id-Vg comparison.
- [Validation evidence](validation/): dated reports and machine-readable
  contracts for checked-in reference TCAD work.
- [LDMOS performance summary, 2026-09-11](validation/templates_ldmos_stage4_performance_summary_2026-09-11.md):
  implemented optimizations, Release verification, and the validated Vg=8 V
  linked-continuation experiment with its reproduction limits.
- [LDMOS current validation](validation/templates_ldmos_current_status.md):
  current qualification boundaries, configuration profiles and reproducible
  linked D4/D5 entry points.
- [IALMob explicit and generated partials, 2026-09-17](validation/templates_ldmos_symbolic_partials_2026-09-17.md):
  high-field analytic partials, branch-preserving low-field code generation,
  independent AD/finite-difference checks and frozen-seed/curve controls; implementation-stage evidence.
- [IALMob kernel costs and same-VM follow-up, 2026-09-17](validation/templates_ldmos_symbolic_vm_2026-09-17.md):
  Release kernel microbenchmarks and deadline-protected full-curve controls; qualification tracked separately from R10.
- [Explicit high-field full repeated qualification, 2026-09-17](validation/templates_ldmos_symbolic_full_2026-09-17.md):
  recommended explicit R11 profile; 12 full curves, exact repeated trajectories, real restart and 62-point physical audit; R10 rollback.
- [R11 G3/D5/D4 isothermal full sweeps, 2026-09-17](validation/templates_ldmos_r11_isothermal_2026-09-17.md):
  final-source isothermal regression: all 155 exact points and both dual-gate joint scores pass; original physics/gates, full continuation distinct from saved-point reclosure.
- [R9 UCRT64 first8 gprof, 2026-09-17](validation/templates_ldmos_gprof_r9_2026-09-17.md):
  Vg8 Release control and cold/drain-only profiles; exact trajectories, IALMob/assembly hotspots and profiler limits.
- [R9 preparation and kernel follow-up, 2026-09-17](validation/templates_ldmos_preparation_followup_2026-09-17.md):
  immutable preparation reuse, thermal high-field reuse controls and fixed-seed candidate-cost diagnosis.
- [Preparation full VM repeat study, 2026-09-17](validation/templates_ldmos_preparation_full_2026-09-17.md):
  qualified explicit R10 preparation profile; two full dual-gate rounds, native pairing and actual pause/resume.
- [LDMOS work summary, 2026-09-16–17](validation/templates_ldmos_work_summary_2026-09-16_17.md):
  daily-report material covering contact repair, root controls, hotspot profiling and R10 qualification.
- [Safeguarded neutral-root Newton, 2026-09-16](validation/templates_ldmos_neutral_root_newton_2026-09-16.md):
  original-equation scalar root acceleration, finite-budget rounding guard and isolated representative controls.
- [Newton and assembly cost controls, 2026-09-16](validation/templates_ldmos_newton_cost_2026-09-16.md):
  explicit R9 contact recommendation, exact neutral-root caching and failed local QF clipping controls.
- [LDMOS production electrothermal reproduction](validation/templates_ldmos_production_reproduction.md):
  recommended explicit R11 high-field profile, R10 rollback and compatible R7 export; Release UMFPACK execution and qualification limits.
- [Extrapolation and continuation research, 2026-09-15](validation/vela_extrapolation_solver_research_2026-09-15.md):
  primary papers and open-source algorithms, with proposed R7/R8 predictor and Newton experiments.
- [Extrapolation A–E execution, 2026-09-15](validation/templates_ldmos_extrapolation_execution_2026-09-15.md):
  guarded predictor, step-growth, tangent, local-fit and NGMRES controls; complete costs,
  frozen evidence and negative results, with R7 retained as the production default.
- [Newton update execution, 2026-09-15](validation/templates_ldmos_newton_update_execution_2026-09-15.md):
  F0 reproduces all 74 R7 attempts exactly; local density projection and NLEQ_ERR-type controls
  fail promotion, and conservative Jacobian refresh shifts cost to extra updates/assemblies.
  R7 remains the production default; projected pseudo-transient mass-matrix research is recorded.
- [Carrier pseudo-transient validation, 2026-09-15](validation/templates_ldmos_pseudo_transient_2026-09-15.md):
  storage signs, temperature/statistics partials and SG diffusion verified; twelve controls
  fail promotion, with fixed-state coefficient isolation identifying non-descent steady-residual directions.
- [Pseudo-transient acceptance validation, 2026-09-15](validation/templates_ldmos_pseudo_acceptance_2026-09-15.md):
  actual backward-Euler defect trials cross initial rejection, but both time controllers fail
  steady gates within 60 updates at each frozen point; R7 remains unchanged.
- [Poisson-consistent initialization controls, 2026-09-15](validation/templates_ldmos_poisson_initialization_2026-09-15.md):
  both preparations close with QF/temperature held; model-controlled PTC reaches near-steady
  states but fails original carrier-row gates, with all preparation costs retained. No promotion.
- [Near-steady mass, density and QF representation controls, 2026-09-15](validation/templates_ldmos_near_steady_isolation_2026-09-15.md):
  state-only local reference changes plus one QF update close both frozen failures under original gates;
  density mapping adds rejected trials. That isolation stage did not test automatic switching or full-path acceleration.
- [Automatic near-steady switch, 2026-09-15](validation/templates_ldmos_near_switch_2026-09-15.md):
  opt-in switching closes both full representative trajectories after Poisson preparation in 32/31 total
  updates; direct raw seeds still fail. Original R7 needs 15/15, so the candidate is not promoted.
- [Standalone R7 high-voltage QF rebase, 2026-09-15](validation/templates_ldmos_r7_local_rebase_2026-09-15.md):
  four points with two paired repeats on each platform pass original gates; Windows/Linux updates
  remain 44/45 per round, with extra reassembly and no reduction of the Linux five-update floor tail.
  Diagnostic only; production defaults unchanged.
- [Finite-contact row arithmetic audit, 2026-09-16](validation/templates_ldmos_contact_row_roundoff_2026-09-16.md):
  exact frozen R7 prefixes isolate the high-voltage floor tail to neutral-contact potential consistency,
  amplified by finite hole exchange; high-precision sum/SG/recombination arithmetic does not remove it.
  Existing contact projection closes the five failed frozen states under original row/block gates;
  that audit alone did not qualify an automatic Newton candidate.
- [Finite-contact algebraic consistency candidate, 2026-09-16](validation/templates_ldmos_contact_consistency_2026-09-16.md):
  an opt-in terminal repair passes 32 paired original-seed high-voltage runs on Windows/Linux;
  Linux's 14-update tail becomes 9 under unchanged gates. Per-round updates fall 44/45 to 40/40,
  but Linux representative-point wall time reverses between repeats. Subsequent full-curve
  qualification is linked below; production R7 defaults are unchanged.
- [Contact-repair sweep validation, 2026-09-16](validation/templates_ldmos_contact_sweep_2026-09-16.md):
  first8 dual-gate controls and real pause/resume pass with exact R7 states and unchanged
  152/137 drain updates; initialization is explicitly protected. All 62 points and two rounds
  of R7/candidate/native controls pass the original joint gates with exactly repeated trajectories.
  Updates fall 366/321 to 344/302; median wall time falls 3.53%/7.20%, with variable per-round
  benefit and higher CPU cost than native. The candidate remains opt-in.
- [Newton update-mapping and damping plan, 2026-09-15](validation/templates_ldmos_newton_update_plan_2026-09-15.md):
  R7 iteration-phase decomposition, the R8 global positivity-alpha failure mechanism, and the
  proposed F0–F5 stage contracts and confirmed design choices (per-node density projection first);
  implementation and local VM experiments are authorized; see the execution report for completed controls.
- [Linux/UCRT64 environment alignment study](validation/vela_linux_ucrt_environment_alignment_2026-09-14.md):
  measured package versions, Linux dependency dry-runs, source-build gaps and unqualified deployment plans.
- [LDMOS daily report, 2026-09-12](validation/templates_ldmos_daily_report_2026-09-12.md):
  physics and electrothermal qualification summary with reproducible static
  current, temperature-rise and nodal-temperature comparison figures.
- [LDMOS daily report, 2026-09-13](validation/templates_ldmos_daily_report_2026-09-13.md):
  finite contact and Auger alignment, native Poisson geometry qualification,
  and serial performance results with reproducible comparison figures.

Optional feature switches used by this repository:

- `VELA_ENABLE_HDF5` (default ON): enables Sentaurus inventory/export support
  when an HDF5 package is found by CMake.
- `VELA_ENABLE_PYTHON` (default OFF): enables the pybind11 Python module and
  `python_api` CTest target.

See `CMakePresets.json` for the shipped Windows UCRT64 preset combinations.

## External Fixture And Test References

- [Regression README](../tests/regression/README.md): engineering regression
  runner behavior, summary JSON fields, and assertion configuration.
- [Reference TCAD README](../reference_tcad/README.md): neutral CSV export
  format and comparison workflow.

## Evidence And History Policy

Dated files under `validation/` are immutable point-in-time evidence and can
contain decisions that were later superseded. Read the current validation
summary first, then follow its links to the relevant evidence. Design specs and
the small number of execution plans still referenced by retained evidence are
provenance records, not active work queues.

Treat the root `README.md`, this index, `config_schema.md`, checked-in
`reference_tcad/` fixtures, CMake targets, tests, and current source code as
the source of truth for current behavior.
