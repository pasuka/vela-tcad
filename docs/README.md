# Vela Documentation

This directory is organized around how people use the project: build it, choose
a configuration schema, run examples, validate reference fixtures, and inspect
archived planning notes.

## Current References

- [Architecture](architecture.md): source tree map, solver paths, and supported
  implementation boundaries.
- [Config schema](config_schema.md): implementation-aligned JSON field
  reference for Poisson, DC sweeps, Newton, unit-scaling input mode, contacts,
  boundaries, solver options, and regression blocks.
- [Examples](examples.md): support matrix for every checked-in example deck and
  the regression expectations tied to those decks.
- [Poisson unit-scaling notes](development_poisson_unit_scaling.md): developer
  notes for the scaled Poisson assembly path used by
  `scaling.mode = "unit_scaling"`.
- [Validation notes](validation/): trend-validation summaries for checked-in
  reference TCAD fixtures.

## External Fixture And Test References

- [Regression README](../tests/regression/README.md): engineering regression
  runner behavior, summary JSON fields, and assertion configuration.
- [Reference TCAD README](../reference_tcad/README.md): neutral CSV export
  format and comparison workflow.

## Archive Notes

These documents are historical project planning and handoff artifacts:

- [Agent handoff baseline, 2026-05-13](agent_handoff_baseline_2026-05-13.md)
- [Merged implementation backlog, 2026-05-13](multi_agent_merged_backlog_2026-05-13.md)
- [Weekly development summary, 2026-05-14](weekly_development_summary_2026-05-14.md)

They are kept for traceability. Treat `README.md`, `docs/config_schema.md`,
`docs/examples.md`, CMake targets, tests, and current source code as the source
of truth for current behavior.
