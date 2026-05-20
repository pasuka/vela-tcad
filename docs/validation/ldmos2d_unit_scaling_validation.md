# LDMOS Unit Scaling Validation

Date: 2026-05-20

This page records the checked-in LDMOS reference_tcad chain for Vela
`unit_scaling`. It is engineering trend validation only and makes no
calibration claim. In short: no calibration claim.

## Fixture

- Source export: `reference_tcad/ldmos2d/nodes.csv`, `elements.csv`,
  `contacts.csv`, and `doping.csv`.
- Structure: coarse mixed Si/SiO2 LDMOS-like lateral device with `p_body`,
  `n_source`, `n_drift`, `n_drain`, and `gate_oxide`.
- Contacts: `source`, `body`, `drain`, and `gate`.
- Decks: `reference_tcad/ldmos2d/vela/simulation_iv.json`,
  `simulation_bv.json`, and `simulation_bv_fieldplate.json`.
- Scaling: every validation deck uses
  `"scaling": {"mode": "unit_scaling"}`.

## Scope

- low-bias DD-IV finite current trend
- off-state breakdown diagnostic using max field trend
- field-plate variant comparison for relative max field behavior

## Result

The checked-in reports show monotone low-bias current and non-decreasing max
field over the sampled reverse-bias points. The field-plate deck is retained as
an engineering trend validation fixture, not a calibrated LDMOS model. This
document makes no calibration claim.

The comparison uses explicit CSV/text reference data and checked-in Vela
outputs. It validates signs, trends, and key orders of magnitude for this
fixture.
