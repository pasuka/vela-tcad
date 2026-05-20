# NMOS2D reference fixture

This directory contains a neutral reference_tcad validation fixture for a
mixed Si/SiO2 NMOS deck using `unit_scaling`.

## Structure

- Device: coarse 2D mixed Si/SiO2 NMOS prototype.
- Silicon regions: `p_body`, `n_source`, and `n_drain`.
- Oxide region: `gate_oxide` with material `SiO2`.
- Contacts: `source`, `drain`, `body`, and a `metal_gate` gate contact.
- interface charge: fixed and trapped charge diagnostics are exercised on the
  `p_body` / `gate_oxide` interface.
- surface mobility: Id-Vg surface-mobility deck uses
  `caughey_thomas_field_surface` at the Si/SiO2 interface.
- Vela decks: generated or derived with `"scaling": {"mode": "unit_scaling"}`.
- all-Si MOS baseline: the existing `examples/nmos2d_dd` decks remain the
  stable all-Si MOS reference baseline.

## Validation coverage

- Id-Vd drain sweep at fixed positive gate bias.
- Id-Vg gate sweep with current increasing under positive gate bias.
- multi-terminal quasi-static CV for gate, drain, source, and body charges.
- Off-state drain BV/high-field diagnostic.

The first dataset checks signs, trends, and key orders of magnitude only. It is
not calibrated to a process or external simulator.
