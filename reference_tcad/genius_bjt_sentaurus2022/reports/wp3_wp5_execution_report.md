# Genius BJT Vela WP3-WP5 execution report

Date: 2026-09-03

## Outcome

WP3 structure conversion, WP4 Vela model-ladder execution, and WP5
Sentaurus/Vela current comparison are complete. Both Vela collector sweeps
converged through 3 V. The exact 31 requested collector voltages, direct
three-terminal currents, and KCL operational gates pass for M0 and M1.

The M1 numerical-parity gate also passes after aligning the electron
Scharfetter maximum lifetime and Nc/Nv with T-2022.03-SP2 Silicon. M0 remains
an intentionally minimal characterization-only baseline. A newly registered
spatial-state gate also passes for potential, electron density, and hole
density. Transport and recombination gates are now asserted as well; the final
acceptance is false only because the global hole-current-density gate still
fails. Electron current density, SRH, and Auger pass.

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
| M1 | 31 | 9.605e-10 | 2.807532e-6 | 2.785499e-6 | 0.992152 | 6.018761e-8 | 5.989697e-8 | 0.995171 | 46.6464 | 46.5048 |

Over VCE=0.5-3.0 V, the M1 maximum absolute log10 errors are 0.003437 decades
for Ic, 0.002102 decades for Ib, 0.003407 decades for Ie, and 0.001349 decades
for beta. All are below the pre-registered 0.05-decade limit. This corresponds
to maximum magnitude-ratio deviations of about 0.79%, 0.49%, 0.78%, and 0.31%,
respectively. M0 is useful
as a numerical baseline but does not show comparable quantitative agreement.

## Spatial-state acceptance and hole-density localization

The exact common mesh at VBE=0.70 V and VCE=3.00 V is used without
interpolation. Potential is gated on all nodes. Density gates use only nodes
where the Sentaurus reference carrier concentration is at least 1e10 cm^-3;
full-domain metrics remain visible as characterization.

| Field | Selected nodes | RMSE | P95 absolute error | Maximum absolute error | Gate |
|---|---:|---:|---:|---:|---:|
| Potential (V) | 5611 | 0.0010965 | 0.0022881 | 0.0150609 | pass |
| Electron density (decade) | 5196 | 0.0029782 | 0.0065252 | 0.0232410 | pass |
| Hole density (decade) | 2207 | 0.0069382 | 0.0077658 | 0.100667 | pass |

The full-domain hole-density maximum of 2.8708 decades is caused by reference
concentrations of roughly 1-100 cm^-3, primarily outside the emitter window and
in the bulk. At the largest-error node Sentaurus gives 32.03 cm^-3 and Vela
2.3793e4 cm^-3. Both are physically negligible compared with the registered
1e10 cm^-3 relevance floor. The emitter-window P95 error is 0.00289 decade.

## Transport and recombination acceptance

The reproducible 3 V Vela export contains electron and hole current-density
vectors and independent SRH/Auger rates. Reference-relative masks plus
integral, normalized-shape, and direction metrics now form asserted initial
engineering regression gates. Single-node peak locations remain diagnostic
only because a nearly flat maximum makes `argmax` location unstable.

| Quantity | Main diagnostic result |
|---|---|
| Electron current density | cosine 0.897; normalized RMSE 0.564; pass |
| Hole current density | cosine 0.788; normalized RMSE 0.815; fail |
| SRH recombination | integral ratio 0.779; shape TV 0.106; pass |
| Auger recombination | integral ratio 0.988; shape TV 0.00658; pass |

This records an important residual: terminal currents agree within the asserted
limits, while local reconstructed transport/source fields are not yet close
enough to justify a separate numerical-parity claim.

A two-node Auger audit showed that SDevice and Vela carrier densities agreed
within about 0.4% while the former Vela rate was about one half of SDevice.
The new opt-in classical `n*p-ni_eff^2` mode raises the Auger integral ratio
from 0.575 to 0.988 and reduces normalized L1 from 0.425 to 0.0122. The full
31-point curve remains accepted. On the real 3 V state, the source-only
analytic Jacobian agrees with centered finite differences to `5.07e-9`
relative error.

The base-collector current audit covers 516 vertically oriented edges. Vela's
nodal hole-current reconstruction agrees with the SDevice nodal projection in
this window (normalized RMSE 0.0628; cosine 0.99925), but the production SG
edge flux has RMSE 0.7189 and requires a 3.53 fitted scale. Direction is already
consistent. The remaining hole-current failure is therefore assigned to nodal
current recovery/semantics; a conservative dual-face recovery is the next
implementation target.

## M1 SRH parameter closure

The original Vela M1 deck used an electron `tau_max` of `3e-8 s`. A direct
`sdevice -P:Silicon` export from T-2022.03-SP2 gives the Scharfetter pair
`taumax = 1e-5, 3e-6 s` for electrons and holes. The remaining Scharfetter
inputs already matched: `taumin=0`, `Nref=1e16 cm^-3`, `gamma=1`, and
`Etrap=0 eV`. The exported `models.par` diagnostic has SHA-256
`aab018ee48a57521decb4d70fb3be4f67f02a1c6e082c930330e5a838a810b4f`.

The numerical gate was registered before the full rerun in
`contracts/comparison_thresholds.json`: over VCE=0.5-3.0 V, the maximum
absolute log10 errors for Ic, Ib, Ie, and beta must each be no greater than 0.05
decades. The comparison command now returns failure when any asserted gate
fails.

## Reproduction

Run the four Vela stages in order with `build-release/vela_example_runner.exe`:

1. `vela/configs/m0_base_ramp.json`
2. `vela/configs/m0_collector_sweep.json`
3. `vela/configs/m1_model_relaxation.json`
4. `vela/configs/m1_collector_sweep.json`

Generate the unique accepted-state chain and its transport/source products in
this order:

1. `scripts/run_genius_bjt_accepted_state_pipeline.py`
2. `scripts/export_genius_bjt_accepted_transport_sources.py`
3. `scripts/compare_genius_bjt_spatial_fields.py`
4. `scripts/compare_genius_bjt_transport_fields.py`
5. `scripts/compare_genius_bjt_sentaurus_vela.py`

The script CLIs require the fixture, ignored Sentaurus/Vela data, and output
paths shown in their `--help` text. Raw Vela states, diagnostics, VTK exports,
and adaptive-step curves remain under
`build-release/reference_tcad/genius_bjt_sentaurus2022/`, including the
`vela_wp3_wp5` run directory.
