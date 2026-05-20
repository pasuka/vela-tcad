# IGBT Unit Scaling Validation

Date: 2026-05-20

This page records the checked-in IGBT reference_tcad chain for Vela
`unit_scaling`. It is engineering trend validation only and makes no
calibration claim. In short: no calibration claim.

## Fixture

- Source export: `reference_tcad/igbt2d/nodes.csv`, `elements.csv`,
  `contacts.csv`, and `doping.csv`.
- Structure: coarse vertical IGBT-like stack with `p_collector`, `n_buffer`,
  `n_drift`, `p_base`, and `n_emitter`.
- Contacts: `collector`, `gate`, and `emitter`.
- Decks: `reference_tcad/igbt2d/vela/simulation_iv.json`,
  `simulation_high_injection_iv.json`, `simulation_charge_cv.json`,
  `simulation_bv.json`, and `simulation_bv_ii.json`.
- Scaling: every validation deck uses
  `"scaling": {"mode": "unit_scaling"}`.

## Scope

- low-current collector IV
- high-injection IV with stored charge proxy
- stored charge and quasi-static CV trend
- breakdown diagnostic max-field trend
- impact-ionization BV smoke diagnostic

## Result

The checked-in reports show monotone high-injection current, non-negative stored
charge trend data, and non-decreasing max field over the sampled reverse-bias
points. These are engineering trend validation fixtures, not calibrated IGBT
or calibrated breakdown diagnostic results. This document makes no calibration
claim. In short: no calibration claim.

The comparison uses explicit CSV/text reference data and checked-in Vela
outputs. It validates signs, trends, and key orders of magnitude for this
fixture.
