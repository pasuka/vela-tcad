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
- [LDMOS production electrothermal reproduction](validation/templates_ldmos_production_reproduction.md):
  frozen R7 input export, Release UMFPACK execution, pause/resume and qualification limits.
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
