# Vela Documentation

This directory is organized around current project behavior: build it, choose
a configuration schema, import reference fixtures, and inspect validation
evidence.

## Current References

- [Codespaces build environment](codespaces.md): container dependencies and
  PowerShell/Codex commands for remote builds, tests, and result retrieval.

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
