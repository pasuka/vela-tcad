# Vela Documentation

This directory is organized around current project behavior: build it, choose
a configuration schema, import reference fixtures, and inspect validation
evidence.

## Current References

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
- [Safeguarded neutral-root Newton, 2026-09-16](validation/templates_ldmos_neutral_root_newton_2026-09-16.md):
  original-equation scalar root acceleration, finite-budget rounding guard and isolated representative controls.
- [Newton and assembly cost controls, 2026-09-16](validation/templates_ldmos_newton_cost_2026-09-16.md):
  explicit R9 contact recommendation, exact neutral-root caching and failed local QF clipping controls.
- [LDMOS production electrothermal reproduction](validation/templates_ldmos_production_reproduction.md):
  recommended explicit R9 contact profile and compatible R7 export, Release UMFPACK execution, pause/resume and qualification limits.
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
