# Regression Tests

The engineering example regression suite is driven by
`scripts/run_regression.py` and registered with CTest as `regression`. It is the
main end-to-end smoke gate for checked-in example decks.

The reference TCAD CSV tools have separate CTest coverage through
`reference_tcad_regression`, implemented in
`tests/regression/test_reference_tcad_tools.py`.

The optional Python binding has separate CTest coverage. Configure with
`-DVELA_ENABLE_PYTHON=ON` to build the pybind11 module and register
`python_api`.

## Runner Behavior

For each entry in `scripts/run_regression.py::EXAMPLES`, the runner:

1. Copies the source example directory into a build-local work directory.
2. Copies the selected named deck to `simulation.json` for check code that
   expects a canonical config name.
3. Runs `vela_example_runner --config <selected deck>`.
4. Verifies configured expected CSV/VTK outputs.
5. Scans generated outputs for NaN/Inf tokens.
6. Applies convergence, sweep, and device-specific trend checks.
7. Writes `regression_summary.json`.

Default work directory:

```text
build/regression_output/
```

Direct run:

```bash
python scripts/run_regression.py --runner build/vela_example_runner
```

CTest run:

```bash
ctest --test-dir build --output-on-failure -R regression
```

## Covered Example Families

The current regression matrix includes:

- PN diode IV, CV, and BV.
- MOS capacitor Poisson interface continuity.
- Legacy simplified NMOS Poisson.
- All-silicon NMOS/PMOS drift-diffusion IV, Id-Vg, and CV.
- Mixed-material NMOS MOS DD IV, Id-Vg, surface-mobility Id-Vg, CV, and BV.
- Mixed-material PMOS MOS DD IV, Id-Vg, CV, and BV.
- LDMOS Poisson, unit-scaling DD-IV, unit-scaling BV, and field-plate BV trend.
- IGBT Poisson, unit-scaling IV, high-injection IV, charge/CV, BV, and
  impact-ionization BV trend.
- Schottky diode IV prototype.

See [../../docs/examples.md](../../docs/examples.md) for the user-facing
support matrix.

## Boundary And Contact Schema Notes

- `contacts[].type` is parsed by Poisson and DC sweep paths.
- Omitted contact `type` remains backward-compatible and defaults to `ohmic`.
- Explicit `boundaries[]` segments are consumed by the Poisson driver for
  Neumann, insulating, and symmetry boundaries.
- DD sweep regressions keep the historical natural zero-flux behavior on
  non-contact boundaries.
- Schottky contacts are a prototype Gummel/Poisson path and are intentionally
  rejected by Newton sweeps.

## DC Sweep Regression Configuration

DC sweep decks can tune checks under `regression.dc_sweep` without changing how
the runner is invoked.

Supported stability fields include:

- `expected_rows`: exact number of CSV rows expected from the configured sweep.
- `max_abs_attempted_step`: largest absolute attempted continuation step
  allowed in any row's `step_diagnostics`.
- `max_abs_accepted_step`: largest absolute accepted continuation step allowed
  in any row's `step_diagnostics`.
- `max_retry_count`: largest allowed per-row retry count.
- `min_converged_rows`: minimum rows with `converged=1`.
- `allow_nonconverged_final_bv_point`: for BV sweeps only, allows the final row
  to be non-converged so a last-stable-before-nonconvergence breakdown row can
  be reported.
- `require_monotone_abs_current`: requires `abs(current_total)` to be
  non-decreasing across converged rows.
- `current_monotone_abs_tolerance` and `current_monotone_rel_tolerance`:
  tolerances for current trend checks.
- `require_monotone_max_field`: for BV sweeps, requires
  `max_electric_field_V_per_m` to be non-decreasing across converged rows.
- `max_field_monotone_abs_tolerance` and
  `max_field_monotone_rel_tolerance`: tolerances for max-field trend checks.
- `min_max_electric_field_V_per_m` and `max_max_electric_field_V_per_m`:
  bounds for converged BV max-field diagnostics.
- `allow_zero_capacitance`, `expected_zero_capacitance_rows`, and
  `min_nonzero_capacitance_rows`: CV-specific capacitance checks.

Regression failures include the example name, CSV row index where applicable,
field name, and actual value.

## Device-Specific Regression Blocks

The runner also consumes top-level `regression` blocks for specialized smoke
checks:

- `mos`: Id-Vd polarity/monotonicity and optional generated Id-Vg trend checks.
- `surface_mobility`: compares a surface-mobility Id-Vg variant against a
  baseline deck.
- `ldmos_iv`: checks LDMOS low-bias drain-current sign and trend.
- `ldmos_fieldplate_trend`: compares field-plate and baseline LDMOS BV max
  field using a ratio limit.
- `igbt_high_injection`: compares high-injection IV against a baseline and can
  require stored-charge monotonicity.
- `igbt_charge_cv`: checks stored charge, multi-terminal charge columns, and
  nonzero quasi-static capacitance rows.
- `igbt_bv`: compares impact-ionization BV smoke output against a baseline at
  a matched bias.
- `schottky_iv`: validates Schottky IV current sign and monotonic trend.

These checks are regression guards, not calibrated silicon acceptance criteria.

## `regression_summary.json`

For every example with `dc_sweep_regression`, the summary includes:

- `rows`: total CSV rows written.
- `converged_rows`: number of rows with `converged=1`.
- `max_attempted_step`: maximum absolute attempted step seen in
  `step_diagnostics`.
- `max_accepted_step`: maximum absolute accepted step seen in
  `step_diagnostics`.
- `max_retry_count_seen`: maximum retry count seen in `step_diagnostics`.
- `final_current_total`: final row `current_total`.
- `current_trend_checked`: whether monotone current checking ran.
- `max_electric_field_seen`, `breakdown_detected`, and
  `max_field_trend_checked`: BV-only electric-field and breakdown diagnostics.
- `breakdown_criterion`, `breakdown_voltage`, `last_stable_bias`, and
  `failed_bias`: BV-only last-breakdown-row semantics.
- `stored_charge_final`: final-row `stored_charge_C_per_m` or
  `stored_charge_C` when stored-charge diagnostics are present.
- `fieldplate_variant_name`: present on the LDMOS field-plate trend check.
- `nonzero_capacitance_rows`: CV-only count of nonzero differential
  capacitance rows after the initial bias point.

The script also preserves legacy summary keys such as
`max_abs_attempted_step`, `max_abs_accepted_step`, and `max_retry_count` for
consumers that parse older summaries.

## Reference TCAD Tool Coverage

`reference_tcad_regression` verifies:

- neutral CSV exports convert to `unit_scaling` Vela decks;
- checked-in PN, NMOS, PMOS, LDMOS, and IGBT fixture assets are complete;
- comparison JSON/Markdown reports can be generated from CSV curves; and
- validation documents retain the expected "trend and order-of-magnitude" and
  "no calibration claim" language.
