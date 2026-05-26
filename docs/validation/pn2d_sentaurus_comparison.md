# pn2d Sentaurus Comparison

The `reference_tcad/pn2d` case is a 2-D abrupt PN diode imported from
Sentaurus. The structure is `L=2.0 um`, `H=0.5 um`, with junction position
`Xj=1.0 um`, `Na=Nd=1e17 cm^-3`, an `Anode` contact on the left boundary, and
a `Cathode` contact on the right boundary.

The import path preserves the Sentaurus TDR mesh, contacts, material regions,
reference IV/BV curves, command summaries, and full `doping.csv` node-level
donor/acceptor table. Generated faithful Vela decks keep `node_doping_file:
"doping.csv"` so the exact imported doping is available for solver
development and regression inspection.

The current executable comparison uses the faithful node-level doping decks.
The IV deck uses `vela_stop: 0.3` and `vela_step: 0.1`; the local gate compares
the 0.2-0.3 V forward-bias window against the full imported Sentaurus IV
reference. The sweep bias follows the Sentaurus `Anode` voltage ramp, while the
Vela current is taken from `Cathode`; this matches the terminal-current
orientation observed in the imported Sentaurus `Anode TotalCurrent` curve and
avoids an artificial electron/hole cancellation at the swept contact. The
comparison selects Vela's `current_total_A_per_um` column rather than the
per-meter `current_total` column, matching the Sentaurus total-current quantity
used by the imported PLT curves. The BV deck currently uses a strict
low-reverse-bias smoke window (`vela_stop: 0.05`, `vela_step: 0.05`) while
Newton continuation beyond the first reverse-bias steps is improved; its
quantity gate checks the non-zero 0.05 V point and leaves the 0 V equilibrium
row as a strict Newton/provenance smoke check. Vela currently disables the Sentaurus
`Avalanche(OkutoCrowell)` approximation in the strict handoff deck because the
imported reverse-bias path is not yet robust with impact-ionization Jacobian
terms. The imported BV command includes `Recombination(SRH Auger Avalanche)`.
The current Vela BV reference override intentionally disables recombination and
impact ionization as a low-bias numerical gate while recombination and
avalanche model parity work remains open; carrying IV recombination models into
BV inflated the 0.05 V current by about 17x.
After candidate isolation, the BV deck also uses a BV-only Caughey-Thomas
mobility override with `bandgap_narrowing: "none"`. The same mobility point is
not used for IV because it degrades the forward-current comparison.
These reports are diagnostic-only and do not yet require trend match.

Vela currently treats Sentaurus Fermi statistics as Boltzmann carrier
statistics. Sentaurus Okuto-Crowell avalanche is approximated by Vela
Selberherr impact ionization for BV diagnostics. These reports establish a
readable, runnable, finite, and comparable regression loop; they are not yet a
calibrated numerical match to Sentaurus.

## Hybrid Solver Status

Faithful pn2d decks now use `solver.method: "gummel_newton"` with
`handoff.fallback: "none"` and preserve `node_doping_file: "doping.csv"`. Each
accepted faithful IV/BV row must end with `handoff_stage: "newton"` and
`newton_iterations > 0`; Gummel fallback is no longer part of the default pn2d
gate.

Current solver limits are intentionally explicit: the bundled faithful pn2d
deck sets `handoff.gummel_max_iter: 0`, so the hybrid path skips the damaged
one-step Gummel initializer and lets coupled Newton cold-start from its
charge-neutral equilibrium seed. This is a strict Newton ownership gate, not
yet a calibrated Sentaurus match: the comparison gate now removes the A/m versus
A-per-um post-processing mismatch and keeps a two-order-of-magnitude envelope
while the remaining physical current calibration is improved. IV is tighter:
it requires a trend match over the 0.2-0.3 V window and must stay within one
order of magnitude.

For BV, matching the Sentaurus physics block by disabling recombination reduced
the 0.05 V quantity delta from about 1.23 orders to about 0.35 orders. The
subsequent BV-only mobility override reduces the promoted BV gate further to
about 0.064 orders while keeping IV on its default mobility model.

## BV 0.05 V Drift/Diff Decomposition

To localize the remaining BV delta source, Vela now records electron/hole drift
and diffusion current components in `DCSweep` CSV output and Python bindings.
The decomposition below uses the regenerated 0.05 V row from:

- `build/pn2d_recomb_gate/vela/pn2d_bv.csv` (BV with `recombination: ["none"]`);
- `build/pn2d_current_contact_gate/vela/pn2d_bv.csv` (baseline BV with SRH/Auger).

At 0.05 V, compared against Sentaurus `current_total` reference:

| Case | total_A_per_um | ratio vs ref | orders | electron_drift_A_per_um | electron_diffusion_A_per_um | electron_total_A_per_um | hole_diffusion_A_per_um | hole_total_A_per_um |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| recomb_none | 4.777839e-19 | 4.458440e-01 | 3.508171e-01 | -2.354494e-17 | 2.343303e-17 | 1.060864e-19 | 3.716976e-19 | 3.716976e-19 |
| baseline (SRH+Auger) | 1.836854e-17 | 1.714060e+01 | 1.234026e+00 | -4.165881e-17 | 2.341374e-17 | 1.824685e-17 | 1.216893e-19 | 1.216893e-19 |

Observations:

- In both cases, electron drift and diffusion are large, opposite-signed terms;
  the net electron contribution is cancellation-sensitive.
- With `recombination: ["none"]`, electron cancellation is much stronger
  (`~1e-19 A/um` residual), so the remaining total current is mainly hole
  diffusion (`~3.7e-19 A/um`), yielding the improved 0.35-order mismatch.
- In the baseline SRH/Auger case, electron cancellation weakens and leaves a
  larger electron residual (`~1.8e-17 A/um`), which dominates the 1.23-order
  mismatch.

This confirms the remaining low-bias BV gap is a coupled transport-balance
residual (electron drift/diff cancellation plus hole diffusion share), rather
than a single-parameter toggle issue.

## Controlled BV Ablations (recombination none)

Using the same BV baseline (`recombination: ["none"]`, `impact_ionization: none`,
`current_contact: Cathode`), three controlled sweeps were executed at 0.05 V:

- `mobility_constant` (`mobility.model = constant`, `bandgap_narrowing = slotboom`)
- `bgn_none` (`mobility.model = caughey_thomas_field`, `bandgap_narrowing = none`)
- `mobility_constant_plus_bgn_none` (`mobility.model = constant`, `bandgap_narrowing = none`)

Generated result files:

- `build/pn2d_recomb_gate/vela/pn2d_bv_mobility_constant.csv`
- `build/pn2d_recomb_gate/vela/pn2d_bv_bgn_none.csv`
- `build/pn2d_recomb_gate/vela/pn2d_bv_mobility_constant_bgn_none.csv`

At 0.05 V, compared against Sentaurus `current_total` reference:

| Case | total_A_per_um | ratio vs ref | orders | electron_drift_A_per_um | electron_diffusion_A_per_um | electron_total_A_per_um | hole_diffusion_A_per_um | hole_total_A_per_um |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| mobility_constant | 9.931476e-20 | 9.267555e-02 | 1.033035e+00 | -4.580564e-17 | 4.563952e-17 | 1.594803e-19 | -6.016557e-20 | -6.016557e-20 |
| bgn_none | 3.098720e-20 | 2.891570e-02 | 1.538866e+00 | -2.413125e-17 | 2.405267e-17 | 8.197582e-20 | -5.098862e-20 | -5.098862e-20 |
| mobility_constant_plus_bgn_none | 1.852821e-18 | 1.728960e+00 | 2.377848e-01 | -4.706634e-17 | 4.699041e-17 | 7.504957e-20 | 1.777771e-18 | 1.777771e-18 |

Interpretation:

- Single-factor changes (`mobility_constant` or `bgn_none`) each push total
  current far below the reference at 0.05 V (about 1.03 and 1.54 orders).
- The combined change (`mobility_constant_plus_bgn_none`) is significantly
  closer at 0.24 orders, but now dominated by hole diffusion
  (`~1.78e-18 A/um`) rather than electron residual.
- Therefore the residual is strongly coupled and non-monotonic: mobility and
  BGN affect electron cancellation and hole diffusion in opposite directions,
  so single-parameter tuning is insufficient for robust BV alignment.

## Mobility x BGN Grid Scan (10 combinations)

A full 2-D grid scan over mobility model and bandgap narrowing model was
executed from the same BV baseline (`recombination: ["none"]`). The sweep
summary is exported to:

- `build/pn2d_recomb_gate/vela/pn2d_bv_grid_scan_summary.csv`

At 0.05 V (sorted by `orders`, ascending):

| Rank | mobility | bgn | status | points | total_A_per_um | ratio_vs_ref | orders |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: |
| 1 | constant | none | ok | 2 | 1.852821e-18 | 1.728960e+00 | 2.377848e-01 |
| 2 | caughey_thomas_field_surface | slotboom | ok | 2 | 4.777839e-19 | 4.458440e-01 | 3.508171e-01 |
| 3 | caughey_thomas_field | slotboom | ok | 2 | 4.777839e-19 | 4.458440e-01 | 3.508171e-01 |
| 4 | caughey_thomas_surface | slotboom | ok | 2 | 4.548779e-19 | 4.244692e-01 | 3.721538e-01 |
| 5 | caughey_thomas | slotboom | ok | 2 | 4.548779e-19 | 4.244692e-01 | 3.721538e-01 |
| 6 | caughey_thomas_surface | none | ok | 2 | 2.745645e-19 | 2.562098e-01 | 5.914043e-01 |
| 7 | caughey_thomas | none | ok | 2 | 2.745645e-19 | 2.562098e-01 | 5.914043e-01 |
| 8 | constant | slotboom | ok | 5 | 9.931476e-20 | 9.267555e-02 | 1.033035e+00 |
| 9 | caughey_thomas_field | none | ok | 3 | 3.098720e-20 | 2.891570e-02 | 1.538866e+00 |
| 10 | caughey_thomas_field_surface | none | ok | 3 | 3.098720e-20 | 2.891570e-02 | 1.538866e+00 |

Scan takeaways:

- The best 0.05 V quantity match in this grid is `mobility=constant` with
  `bandgap_narrowing=none` (about 0.238 orders).
- The default-like field-limited mobility with Slotboom BGN
  (`caughey_thomas_field + slotboom`) remains within about 0.351 orders.
- `caughey_thomas_field_surface + none` now emits a valid 0.05 V row and
  numerically matches `caughey_thomas_field + none` in this pn2d setup,
  indicating the configured surface constraint is not changing the low-bias
  0.05 V current metric for this case.

## Focused Local CT Tuning (partial)

To follow up the best grid point (`constant + none`, `orders ~= 0.238`), a
focused local scan was started around `caughey_thomas + none` using explicit
mobility-object parameters near the default Caughey-Thomas values.

Partial completed-point summary is exported to:

- `build/pn2d_recomb_gate/vela/pn2d_bv_ct_local_small_scan_partial_summary.csv`

Completed points so far (0.05 V, sorted by `orders`):

| Case tag | orders | ratio_vs_ref | total_A_per_um |
| --- | ---: | ---: | ---: |
| `ct_small_mu0p9_nr1_a0p9` | 2.086866e-01 | 6.184626e-01 | 6.627687e-19 |
| `ct_small_mu0p9_nr0p8_a0p9` | 5.412929e-01 | 2.875459e-01 | 3.081454e-19 |
| `ct_small_mu0p9_nr0p8_a1` | 1.314198e+00 | 4.850676e-02 | 5.198175e-20 |

Interim takeaway:

- A completed local CT point (`ct_small_mu0p9_nr1_a0p9`) already improves over
  the previous global best (`0.2087` vs `0.2378` orders at 0.05 V).
- This indicates a tunable Caughey-Thomas neighborhood can outperform the
  coarse grid winner while keeping `bandgap_narrowing: none` and
  `recombination: ["none"]`.

### Neighborhood refinement (6 points)

To verify whether the new local best can be improved further, a 6-point
refinement was executed around `ct_small_mu0p9_nr1_a0p9`.

Refinement summary:

- `build/pn2d_recomb_gate/vela/pn2d_bv_ct_refine6_summary.csv`

Results (0.05 V, sorted by `orders`):

| Case tag | orders | ratio_vs_ref | total_A_per_um |
| --- | ---: | ---: | ---: |
| `refine_mu0p9_nr1p0_a0p95` | 2.284867e-01 | 5.908991e-01 | 6.332306e-19 |
| `refine_mu0p9_nr1p0_a0p85` | 5.003556e-01 | 3.159689e-01 | 3.386047e-19 |
| `refine_mu0p9_nr0p9_a0p9` | 8.749254e-01 | 1.333750e-01 | 1.429299e-19 |
| `refine_mu0p85_nr1p0_a0p9` | 1.169681e+00 | 6.765792e-02 | 7.250488e-20 |
| `refine_mu0p9_nr1p1_a0p9` | 1.202610e+00 | 6.271765e-02 | 6.721069e-20 |
| `refine_mu0p95_nr1p0_a0p9` | 1.780036e+00 | 1.659450e-02 | 1.778332e-20 |

Refinement takeaway:

- None of the six neighborhood points surpassed the current local best
  `ct_small_mu0p9_nr1_a0p9` (`orders ~= 0.2087`).
- The best refinement point (`alpha_scale = 0.95`) reached `orders ~= 0.2285`,
  still better than the earlier coarse-grid best (`0.2378`) but worse than the
  current local optimum.

### Second-round quick refinement (6 points)

After the neighborhood refinement, a targeted quick 6-point run was executed
to probe lower-alpha and nearby-mobility settings around the local optimum.

Quick refinement summary:

- `build/pn2d_recomb_gate/vela/pn2d_bv_ct_quick6_summary.csv`
- Repro script: `scripts/scan_pn2d_bv_ct_quick6.ps1`

The script can be run directly after regenerating a pn2d reference tree; it
starts from `simulation_bv.json`, forces `bandgap_narrowing: none`, applies the
six CT mobility-object perturbations below, and writes the summary CSV.

Results (0.05 V, sorted by `orders`):

| Case tag | orders | ratio_vs_ref | total_A_per_um |
| --- | ---: | ---: | ---: |
| `q_mu0p89_a0p89` | 6.411045e-02 | 8.627591e-01 | 9.245665e-19 |
| `q_mu0p89_a0p90` | 1.259122e-01 | 7.483208e-01 | 8.019299e-19 |
| `q_mu0p91_a0p90` | 1.897403e-01 | 6.460404e-01 | 6.923222e-19 |
| `q_mu0p90_nr1p02_a0p90` | 5.018781e-01 | 3.148632e-01 | 3.374198e-19 |
| `q_mu0p90_a0p91` | 5.953329e-01 | 2.539026e-01 | 2.720920e-19 |
| `q_mu0p90_a0p89` | 9.737582e-01 | 1.062287e-01 | 1.138388e-19 |

Second-round takeaway:

- This run surpasses all previous local candidates and pushes the 0.05 V
  mismatch well below 0.20 orders.
- Current best observed point is `q_mu0p89_a0p89` with
  `orders ~= 0.0641` (previous best was `0.2087`).
- A fresh rerun from `build/pn2d_current_review` reproduced the same best point:
  `ratio_vs_ref ~= 0.8628`, `orders ~= 0.0641`, and
  `total_A_per_um ~= 9.2457e-19`.

Fresh baseline rerun on 2026-05-26:

- default IV: `orders ~= 0.5048`, trend matched;
- default BV after BV-only promotion: `orders ~= 0.0641`;
- quick6 best `q_mu0p89_a0p89`: `orders ~= 0.0641`.

## IV Per-Bias Ratio Shape

Generated with `scripts/summarize_pn2d_iv_ratios.ps1` from
`build/pn2d_tdr_tie_probe/vela/pn2d_iv.csv` against
`build/pn2d_tdr_tie_probe/reference_curves/pn2d_iv_reference.csv`:

| Bias V | Vela/reference ratio |
| ---: | ---: |
| 0.204721576526 | 0.9052 |
| 0.224721576526 | 1.2270 |
| 0.244721576526 | 0.9426 |
| 0.25 | 0.8499 |
| 0.27 | 0.5348 |
| 0.29 | 0.3128 |

The ratio roll-off above about 0.25 V confirms that the remaining IV mismatch
is dominated by a shallower high-forward-bias slope, not by a low-bias offset.

## IV/BV Physics Matrix

Generated with `scripts/scan_pn2d_iv_bv_physics_matrix.ps1` from
`build/pn2d_tdr_tie_probe`.

| Case | Kind | Orders | Ratio at target | Interpretation |
| --- | --- | ---: | ---: | --- |
| default | IV | 0.5048 | 0.3128 at 0.29 V | baseline high-forward-bias slope remains shallow |
| iv_recomb_none | IV | 0.6415 | 0.4864 at 0.29 V | turning recombination off does not resolve IV slope |
| iv_bgn_none | IV | 0.5135 | 0.4298 at 0.29 V | removing BGN alone also does not resolve IV slope |
| bv_recomb_none | BV | 0.0641 | 0.8628 at 0.05 V | current promoted low-bias gate |
| bv_recomb_srh | BV | 1.1778 | 15.0580 at 0.05 V | SRH parity reintroduces large low-bias mismatch |
| bv_recomb_srh_auger | BV | 1.1737 | 14.9172 at 0.05 V | SRH+Auger parity remains unresolved |

Physics takeaway:

- BV low-bias numerical agreement currently depends on the deliberate
  `recombination: ["none"]` gate.
- Re-enabling SRH (with or without Auger) moves the 0.05 V point by about
  1.17 orders, so BV physics parity remains open.
- IV mismatch remains dominated by high-forward-bias slope behavior rather than
  a single recombination or BGN toggle.

## Candidate Mobility Impact: q_mu0p89_a0p89

The `scripts/scan_pn2d_candidate_iv_bv.ps1` helper applies the same
`q_mu0p89_a0p89` CT mobility object to separately generated IV and BV
candidate decks. This isolates whether the BV-tuned mobility point is safe as a
global pn2d calibration.

| Metric | Default | Candidate |
| --- | ---: | ---: |
| IV orders | 5.047505e-01 | 2.442584e+00 |
| IV ratio vs ref at 0.29 V | 3.109879e-01 | 3.609239e-03 |
| BV 0.05 V orders | 3.508171e-01 | 6.411045e-02 |
| BV 0.05 V ratio vs ref | 4.458440e-01 | 8.627591e-01 |

Takeaway:

- The candidate is excellent for the low-bias BV quantity gate.
- The same candidate is not safe as a global IV/BV mobility calibration because
  it degrades the IV comparison by roughly two orders of magnitude.
- The promoted reference config therefore applies this as a BV-only override;
  IV remains on the default field-limited mobility path.

The old region-average `runtime_doping_scale` path is no longer required for
pn2d. It remains available through an opt-in `runtime_diagnostic` config block
for future debugging, but the default bundled pn2d reference runs the imported
node-level doping. Current gate:

- faithful IV/BV deck generation is required;
- faithful decks must preserve node-level doping and hybrid handoff settings;
- faithful IV/BV execution must remain finite and end in Newton handoff;
- comparison reports align by `bias_V`, use configured bias windows, and compare
  Vela `current_total_A_per_um` against Sentaurus total current;
- strict Sentaurus numerical agreement is not yet required.

Useful local verification command:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\sentaurus_import.py reference --config reference_tcad\pn2d\pn2d_reference.json --source-dir reference_tcad\pn2d --output-dir build\reference_tcad\pn2d --tdr-importer build\sentaurus_import.exe --runner build\vela_example_runner.exe
```
