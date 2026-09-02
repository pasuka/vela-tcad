# Genius BJT Vela WP3-WP5 execution report

Date: 2026-09-02

## Outcome

WP3 structure conversion, WP4 Vela model-ladder execution, and WP5
Sentaurus/Vela current comparison are complete. Both Vela collector sweeps
converged through 3 V. The exact 31 requested collector voltages, direct
three-terminal currents, and KCL operational gates pass for M0 and M1.

The M1 numerical-parity gate also passes after aligning the electron
Scharfetter maximum lifetime with the T-2022.03-SP2 Silicon default. M0 remains
an intentionally minimal characterization-only baseline.

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
| M1 | 31 | 2.738e-9 | 2.807532e-6 | 2.790669e-6 | 0.993993 | 6.018761e-8 | 6.332852e-8 | 1.05219 | 46.6464 | 44.0665 |

Over VCE=0.5-3.0 V, the M1 maximum absolute log10 errors are 0.002694 decades
for Ic, 0.022092 decades for Ib, and 0.024735 decades for beta. All are below
the pre-registered 0.05-decade limit. This corresponds to maximum magnitude
ratio deviations of about 0.62%, 5.21%, and 5.54%, respectively. M0 is useful
as a numerical baseline but does not show comparable quantitative agreement.

## M1 SRH parameter closure

The original Vela M1 deck used an electron `tau_max` of `3e-8 s`. A direct
`sdevice -P:Silicon` export from T-2022.03-SP2 gives the Scharfetter pair
`taumax = 1e-5, 3e-6 s` for electrons and holes. The remaining Scharfetter
inputs already matched: `taumin=0`, `Nref=1e16 cm^-3`, `gamma=1`, and
`Etrap=0 eV`. The exported `models.par` diagnostic has SHA-256
`aab018ee48a57521decb4d70fb3be4f67f02a1c6e082c930330e5a838a810b4f`.

The numerical gate was registered before the full rerun in
`contracts/comparison_thresholds.json`: over VCE=0.5-3.0 V, the maximum
absolute log10 errors for Ic, Ib, and beta must each be no greater than 0.05
decades. The comparison command now returns failure when any asserted gate
fails.

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
