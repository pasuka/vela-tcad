# Regression tests

The engineering example regression suite is driven by `scripts/run_regression.py`
and registered with CTest as the `regression` test. The script copies each
example into a build-local working directory, runs it with `vela_example_runner`,
checks generated outputs, and writes `regression_summary.json`.

The optional Python binding has separate CTest coverage. Configure with
`-DVELA_ENABLE_PYTHON=ON` to build the pybind11 module and register the
`python_api` test, which imports `vela`, loads a mesh, runs Poisson and DC sweep
flows, and verifies generated CSV/VTK outputs.

## DC sweep regression configuration

DC sweep decks can tune checks under `regression.dc_sweep` without changing how
the runner is invoked. Supported stability fields include:

- `expected_rows`: exact number of CSV rows expected from the configured sweep.
- `max_abs_attempted_step`: largest absolute attempted continuation step allowed
  in any row's `step_diagnostics`.
- `max_abs_accepted_step`: largest absolute accepted continuation step allowed in
  any row's `step_diagnostics`.
- `max_retry_count`: largest allowed per-row `retry_count` from
  `step_diagnostics`.
- `min_converged_rows`: minimum number of rows with `converged=1`.
- `allow_nonconverged_final_bv_point`: for BV sweeps only, allows the final row
  to be non-converged so a last-stable-before-nonconvergence breakdown row can
  be reported while earlier rows remain strict.
- `require_monotone_abs_current`: requires `abs(current_total)` to be
  non-decreasing across converged rows. Optional tolerances are
  `current_monotone_abs_tolerance` and `current_monotone_rel_tolerance`.
- `require_monotone_max_field`: for BV sweeps, requires
  `max_electric_field_V_per_m` to be non-decreasing across converged rows.
  Optional tolerances are `max_field_monotone_abs_tolerance` and
  `max_field_monotone_rel_tolerance`.
- `min_max_electric_field_V_per_m` / `max_max_electric_field_V_per_m`: lower and
  upper bounds for BV `max_electric_field_V_per_m` on converged rows.
- `allow_zero_capacitance`, `expected_zero_capacitance_rows`, and
  `min_nonzero_capacitance_rows`: CV-specific capacitance checks.

Regression failures include the example name, CSV row index where applicable,
field name, and actual value so threshold regressions can be traced directly to
the failing deck row.

## `regression_summary.json` DC sweep fields

For every example with the `dc_sweep_regression` check, the summary includes:

- `rows`: total CSV rows written.
- `converged_rows`: number of rows with `converged=1`.
- `max_attempted_step`: maximum absolute attempted step seen in
  `step_diagnostics`.
- `max_accepted_step`: maximum absolute accepted step seen in
  `step_diagnostics`.
- `max_retry_count_seen`: maximum retry count seen in `step_diagnostics`.
- `final_current_total`: final row `current_total`.
- `current_trend_checked`: whether `require_monotone_abs_current` was enabled
  and checked.
- `max_electric_field_seen`, `breakdown_detected`, and
  `max_field_trend_checked`: BV-only electric-field and breakdown diagnostics.
- `nonzero_capacitance_rows`: CV-only count of non-zero differential
  capacitance rows after the initial bias point.

The script also keeps the legacy `max_abs_attempted_step`,
`max_abs_accepted_step`, `max_retry_count`, and BV field-list keys for consumers
that already parse older summaries.
