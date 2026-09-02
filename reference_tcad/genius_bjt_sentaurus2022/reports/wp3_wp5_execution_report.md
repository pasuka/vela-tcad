# Genius BJT Vela WP3-WP5 execution report

Date: 2026-09-02

## Outcome

WP3 structure conversion, WP4 Vela model-ladder execution, and WP5
Sentaurus/Vela current comparison are complete. Both Vela collector sweeps
converged through 3 V. The exact 31 requested collector voltages, direct
three-terminal currents, and KCL operational gates pass for M0 and M1.

Numerical parity is reported as characterization-only. WP0-WP2 intentionally
reserved the Vela thresholds, so this first run must not invent a post-hoc
tolerance.

## Common input

- 5611 nodes and 10940 triangles, imported from the accepted Sentaurus SDE
  structure.
- Exact node donor and acceptor concentrations; no profile refit.
- Base, emitter, and collector contact-node sets are preserved exactly.
- Emitter is 0 V, base is 0.70 V, and collector is swept from 0 to 3 V.
- Currents are reported in A/um with positive sign into the device.

## Vela continuation

M0 uses a zero-to-0.70 V base ramp followed by the collector sweep. A direct
M1 base ramp from equilibrium encountered a continuity line-search
non-decrease near VBE=2.022 mV after reducing the voltage increment to about
1.4e-8 V. This is a near-zero-current numerical continuation failure.

The accepted M1 path starts from the converged M0 VBE=0.70 V state, performs a
fixed-bias relaxation with the complete M1 physics, and then starts the
collector sweep from that M1 state. The relaxation converged in three Newton
iterations. No M1 physics term was disabled or weakened.

## Comparison result

| Model | Points | Max relative Vela KCL | Sentaurus Ic @ 3 V (A/um) | Vela Ic @ 3 V (A/um) | Ic ratio | Sentaurus Ib @ 3 V (A/um) | Vela Ib @ 3 V (A/um) | Ib ratio | Sentaurus beta | Vela beta |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| M0 | 31 | 2.086e-9 | 1.843315e-6 | 3.397337e-6 | 1.84306 | 1.603892e-8 | 7.127638e-9 | 0.444397 | 114.928 | 476.643 |
| M1 | 31 | 3.889e-9 | 2.807532e-6 | 2.679006e-6 | 0.954221 | 6.018761e-8 | 2.741779e-7 | 4.55539 | 46.6464 | 9.77105 |

Over VCE=0.5-3.0 V, the M1 median absolute log10 error is 0.0204 decades for
Ic and 0.6603 decades for Ib. Thus M1 collector transport is already close to
the Sentaurus oracle, while base-current/recombination behavior is the primary
remaining discrepancy. M0 is useful as a numerical baseline but does not show
comparable quantitative agreement.

## Reproduction

Run the four Vela stages in order with `build-release/vela_example_runner.exe`:

1. `vela/configs/m0_base_ramp.json`
2. `vela/configs/m0_collector_sweep.json`
3. `vela/configs/m1_model_relaxation.json`
4. `vela/configs/m1_collector_sweep.json`

Then run `scripts/compare_genius_bjt_sentaurus_vela.py` with the fixture root,
the ignored Vela run directory, and `comparison/` as its three arguments. Raw
Vela states, diagnostics, and adaptive-step curves remain under
`build-release/reference_tcad/genius_bjt_sentaurus2022/vela_wp3_wp5`.
