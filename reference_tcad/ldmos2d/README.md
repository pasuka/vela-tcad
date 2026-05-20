This directory contains a neutral reference_tcad validation fixture for a
coarse LDMOS-like mixed Si/SiO2 device using `unit_scaling`.

The checked-in data is an explicit CSV/text export for engineering trend
validation. This is an engineering trend validation fixture. It is not a
calibrated power-device model and makes no calibration claim.

## Structure

- Regions: `p_body`, `n_source`, `n_drift`, `n_drain`, and `gate_oxide`.
- Contacts: `source`, `body`, `drain`, and `gate`.
- Coverage: low-bias DD-IV, off-state BV/high-field diagnostics, and a
  field-plate variant for max field trend checks.
- Vela decks: generated or derived with `"scaling": {"mode": "unit_scaling"}`.
