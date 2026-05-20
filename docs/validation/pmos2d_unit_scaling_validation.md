# PMOS2D Unit Scaling Validation

Date: 2026-05-20

This page records the checked-in PMOS mixed Si/SiO2 reference_tcad validation
chain for Vela `unit_scaling`. It is a trend and order-of-magnitude validation
only; there is no calibration claim.

## Fixture

- Source export: `reference_tcad/pmos2d/nodes.csv`, `elements.csv`,
  `contacts.csv`, and `doping.csv`.
- Structure: `n_body`, `p_source`, `p_drain`, `gate_oxide`, source/drain/body
  ohmic contacts, and a `metal_gate` contact.
- Physics knobs covered: fixed interface charge, interface trap charge in CV,
  and a surface mobility Id-Vg variant.
- all-Si MOS baseline: `examples/pmos2d_dd` is retained as the stable all-Si
  MOS deck family.

## Results

| Check | Result |
| --- | --- |
| Id-Vg | Pass: drain current magnitude increases under negative gate bias from `8.689e-8` to `3.881e-7` A/m. |
| Id-Vd | Pass: drain current uses the expected positive reported sign for the PMOS drain-current convention. |
| multi-terminal CV | Pass: gate, drain, source, and body charge/capacitance columns are present and finite. |
| BV max field | Pass: max field rises from `6.575e4` to `8.575e4` to `1.058e5` V/cm. |
| Reference comparison | Pass: Id-Vd, Id-Vg, CV, and BV report JSON files have matching trends. |

The CV comparison is quasi-static trend validation only. It is not compared to
an AC small-signal matrix.
