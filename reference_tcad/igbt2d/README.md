# IGBT2D Reference Fixture

This directory contains a neutral reference_tcad validation fixture for a
coarse IGBT-like 2-D device using `unit_scaling`.

The checked-in data is an explicit CSV/text export for engineering trend
validation. This is not a calibrated IGBT model and makes no calibration claim.
In short: no calibration claim.

## Structure

- Regions: `p_collector`, `n_buffer`, `n_drift`, `p_base`, and `n_emitter`.
- Contacts: `collector`, `gate`, and `emitter`.
- Coverage: low-current IV, high-injection IV, stored charge proxy CV,
  off-state BV diagnostics, and impact-ionization BV smoke diagnostics.
- Vela decks: generated or derived with `"scaling": {"mode": "unit_scaling"}`.

## Validation Coverage

- low-current collector IV.
- high-injection current trend.
- stored charge and quasi-static CV trend.
- breakdown diagnostic max-field trend.
- impact-ionization BV smoke diagnostic.

The fixture is an engineering trend validation guard only.
