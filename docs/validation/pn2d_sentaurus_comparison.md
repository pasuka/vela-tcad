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

## Sentaurus 2018 PN2D Fixture

The Sentaurus 2018 calibration case requested for the next pn2d precision pass
is checked in under `reference_tcad/pn2d_sentaurus2018/source/`. It preserves
the original `D:\pn2d` artifacts, including `.cmd`, `.log`, `.tdr`, `.plt`,
`tdx`-generated `.grd/.dat`, and command backups. Its import config is
`reference_tcad/pn2d_sentaurus2018/pn2d_sentaurus2018_reference.json`.

This fixture has three simulations:

- `0v`: equilibrium state import from `pn2d_0v_des.tdr`, reference curve export
  from `pn2d_0v.plt`, and Sentaurus field export under `sim_fields/0v`.
- `iv`: forward-bias import from `pn2d_iv_des.tdr` and `pn2d_iv.plt`.
- `bv`: reverse-bias import from `pn2d_bv_des.tdr` and `pn2d_bv.plt`.

All generated Vela decks use the current strict coupled-Newton handoff:
`solver.method: "gummel_newton"`, `handoff.gummel_max_iter: 0`, and
`handoff.fallback: "none"`. This is the repository's present way to express a
Newton-only solve of Poisson plus electron and hole continuity while retaining
the existing runner schema. The decks keep `node_doping_file: "doping.csv"`.

The independent input-data gate is `scripts/compare_sentaurus_tdr_tdx.py`.
It compares the TDR-derived neutral export against `tdx` text output for mesh
vertex count, total and bulk element counts, region/contact names, contact edge
counts, and matching `.dat` datasets such as `DopingConcentration`,
`DonorConcentration`, and `AcceptorConcentration` when present. State-field
parity can be run with `scripts/compare_sentaurus_fields.py` for fields such as
`ElectrostaticPotential`, `eDensity`, `hDensity`, quasi-Fermi potentials,
`ElectricField`, SRH recombination, and BV avalanche/impact-ionization exports
where available.

The zero-bias state diagnostic is `scripts/compare_pn2d_0v_state.py`. It
derives a temporary Vela 0V probe deck from the imported
`vela/simulation_0v.json`, forces a single strict-Newton IV point at 0 V, writes
VTK, and compares that state to `pn2d_0v_des.tdr` fields exported under
`sim_fields/0v/fields`. The report includes potential, electron/hole
quasi-Fermi potentials, carrier densities reconstructed from
`psi/phin/phip/ni_eff`, SRH recombination, raw VTK carrier-density cross-checks,
centerline/contact/junction-local statistics, and Anode/Cathode near-zero
terminal-current conservation. The 0 V terminal check follows the Sentaurus
log table convention: Anode and Cathode currents must be equal in magnitude and
opposite in sign, so the report gates both absolute near-zero current and
`abs(I_anode + I_cathode) / max(abs(I_anode), abs(I_cathode))`. The diagnostic
matrix records the current
priority axes: `ni`, OldSlotboom/BGN, Ohmic contact boundary semantics,
quasi-Fermi definitions, carrier formulas, and current units.

The terminal-current root-cause diagnostic is
`scripts/diagnose_pn2d_0v_current_balance.py`. It derives one 0 V probe deck,
enables `sweep.diagnostics.terminal_balance` and multi-contact
`sweep.diagnostics.contact_edge`, then computes Anode and Cathode currents on
the same converged `DDSolution`. Its report is the authoritative source for
0 V terminal-current conservation because it removes the previous ambiguity
from running two independent probes with different `current_contact` values.
The report writes both `electron_minus_hole` and `electron_plus_hole`
conventions, finite-volume flux-link sums, contact coverage, top contributing
links, endpoint `psi/phin/phip/n/p/ni/mu` state, SG continuity fluxes, and a
root-cause classification such as `total_current_sign_convention`,
`contact_current_aggregation`, `contact_edge_coverage`, or
`contact_flux_formula`. With `--require-balanced`, the script becomes a hard
regression gate and returns non-zero unless the same-solution terminal-current
classification is `balanced`.

`balanced` is intentionally only a Vela two-terminal conservation statement. The
same report also parses the Sentaurus `pn2d_0v.plt` final Coupled current table
(falling back to `pn2d_0v.log_des.log`) and writes
`sentaurus_current_reference` plus `sentaurus_current_parity`. That parity block
compares Vela `electron_minus_hole_A_per_um` against Sentaurus `TotalCurrent`,
records the absolute-current ratio and sign relation for each contact, and
flags component-level evidence such as Vela `Anode:hole` or `Cathode:electron`
being numerically zero while Sentaurus reports non-zero components. A report can
therefore be `classification=balanced` while
`sentaurus_current_parity.status=mismatch`; this means current conservation has
been hardened, but absolute Sentaurus current parity remains a separate debug
target.

The follow-up residual-floor probe is
`scripts/probe_pn2d_0v_newton_residual_current.py`. It derives multiple 0 V
probe decks from the same imported fixture, varies Newton `reltol`, `abstol`,
and `max_iter`, runs the same multi-terminal current diagnostics, and ranks
candidates by Sentaurus total-current magnitude only after preferring successful
terminal-balanced candidates. This prevents a numerically unbalanced early-stop
state from being mistaken for Sentaurus parity simply because one terminal
current happens to have a similar absolute magnitude. Current evidence shows
that strict `reltol<=1e-9` candidates remain balanced but at the `1e-33 A/um`
floor, while looser candidates around `3e-9` and above raise residual currents
but fail the two-terminal balance gate.

It also records the Sentaurus/TDR mesh reference
counts separately: `pn2d_0v_des.grd` has `nb_elements=3712`, which includes
`3680` bulk `R.Si` triangle elements plus `16` Cathode and `16` Anode contact
boundary elements. Vela's contact-current diagnostic `flux_link_count` is a
finite-volume post-processing count from contact Dirichlet nodes into the bulk,
so it can differ from the Sentaurus boundary segment count. Sentaurus log lines
516-518 and 547-549 in `pn2d_0v.log_des.log` remain the reference behavior:
Anode and Cathode 0 V conduction currents should be equal magnitude and
opposite sign.

The 0 V deck uses a local tighter Newton tolerance than the IV/BV calibration
decks: `reltol: 1e-10`, `abstol: 1e-24`, and `newton_max_iter: 80`. This keeps
equilibrium terminal currents below the numerical floor where relative balance
is ill-conditioned. The terminal gate therefore accepts either the relative
pair-balance threshold or an absolute conservation threshold (`1e-24 A/um` by
default). This preserves the old failure for visible `~1e-21 A/um` imbalance
while allowing `~1e-33 A/um` roundoff-floor currents from the tightened 0 V
Newton solve.

The field-convention diagnostic is
`scripts/diagnose_pn2d_0v_field_conventions.py`. It tests whether the remaining
0 V state-field mismatch is explained by electrostatic-potential sign, offset,
or affine convention, then ranks carrier-density formulas and `ni/BGN`
candidates. This script is diagnostic-only: it does not change solver physics
or import mappings. Current pn2d 0 V evidence is not a simple sign flip; the
best potential fit is affine but still has a large residual, so the next debug
step is to separate node/field mapping issues from material `ni_eff`/OldSlotboom
differences.

The field-mapping diagnostic is `scripts/diagnose_pn2d_0v_field_mapping.py`.
It compares Vela VTK fields to Sentaurus `sim_fields/0v/fields` by both direct
node-id pairing and nearest-coordinate pairing, while testing whether VTK
coordinates need conversion from meters back to the imported Sentaurus micrometer
mesh. Current pn2d 0 V evidence rules out a node-ordering mismatch: direct and
nearest-coordinate pairings produce the same field errors, and the best
coordinate alignment is `m_to_um`. Carrier-density fields are classified as
Vela VTK `m^-3` against Sentaurus `cm^-3` (`vtk_m3_to_cm3`), but that conversion
only removes the leading `1e6` unit factor; large local density residuals remain
after conversion. This leaves physical field conventions and material/contact
parameters as the next state-parity targets.

The dedicated 0 V electric-field distribution report is
`scripts/compare_pn2d_0v_electric_field.py`. It reads Sentaurus
`ElectricField_region0.csv`, derives Vela field magnitude from the linear
triangle gradient of VTK `Potential`, and writes JSON, Markdown, and per-node
CSV reports. The script tests `V/m`, `V/cm`, and `V/um` candidates and selects
the unit by matching the field-magnitude distribution before reporting
all-node, centerline, junction-near, and contact-local error statistics. Current
pn2d 0 V evidence selects `V_per_cm`: Sentaurus median/max are about
`144 V/cm` and `1.02e5 V/cm`, while Vela-derived median/max are about
`308 V/cm` and `1.05e5 V/cm`. The mean magnitudes are close, but local residuals
remain large near the junction, so this report is a distribution/localization
diagnostic rather than a pass/fail state-parity gate.

The density-decomposition diagnostic is
`scripts/diagnose_pn2d_0v_density_decomposition.py`. It checks whether Vela's
VTK carrier densities are internally consistent with Vela `Potential`,
`ElectronQuasiFermi`, `HoleQuasiFermi`, `ni`, and BGN, then compares those
densities to Sentaurus. Current pn2d 0 V evidence classifies the mismatch as
`vela_state_differs_from_sentaurus`: Vela density export is self-consistent to
about `2.3e-5` relative error after the VTK `m^-3` to `cm^-3` conversion, but
Sentaurus-vs-Vela density residuals remain large. The inferred equilibrium
`sqrt(n*p)` median also differs: Sentaurus is about `1.66e10 cm^-3`, while Vela
is about `1.13e10 cm^-3`, so `ni_eff`/OldSlotboom parity remains a material
parameter target in addition to the electrostatic-potential state mismatch.

The ni/BGN probe is `scripts/diagnose_pn2d_0v_ni_bgn_probe.py`. It creates
temporary 0 V decks under the report directory, writes per-candidate
`materials_file` overrides, toggles `solver.bandgap_narrowing`, runs the same
strict Newton probe, and ranks Sentaurus state parity. Current pn2d 0 V scan
over `ni={1e10,1.45e10,1.6556207295e10} cm^-3` and
`BGN={none,slotboom}` chooses `ni=1.6556207295e10 cm^-3`, `BGN=none` as the
best local 0 V material match. This exactly matches the Sentaurus
`sqrt(n*p)` median, but only reduces the density median log10 error to about
`0.252` and leaves p95 log10 error about `5.44`; potential RMS remains about
`0.138 V`. Therefore `ni_eff` is a confirmed material-parity lever, but the
remaining 0 V mismatch is dominated by electrostatic/contact state parity rather
than a pure intrinsic-density or OldSlotboom toggle.

The QF-driver probe is `scripts/probe_pn2d_0v_qf_drivers.py`. It runs a small
strict-Newton 0 V matrix from the imported `simulation_0v.json` and records
per-variant QF span, terminal current, terminal-balance, and VTK state outputs.
Current evidence selects BGN/effective-ni consistency as the direct driver of
the 0 V quasi-Fermi split: baseline `slotboom` keeps `qf_max_span_V =
0.0043930610509`, disabling recombination still keeps `qf_max_span_V =
0.00439298359789`, and switching `bandgap_narrowing` to `none` collapses the
span to `8.39006e-13 V` while retaining strict Newton handoff
(`gummel_iterations=0`, `handoff_stage=newton`, `newton_iterations=10`). The
`l2_residual` variant does not remove the split (`0.004457594 V`), and the
intentionally over-tight `tight_block_scales` variant fails by
`line_search_non_decrease`; this is recorded as a diagnostic branch failure,
not as a matrix execution failure.

Reproducible import and comparison:

```powershell
python scripts\sentaurus_import.py reference --config reference_tcad\pn2d_sentaurus2018\pn2d_sentaurus2018_reference.json --source-dir reference_tcad\pn2d_sentaurus2018\source --output-dir build\reference_tcad\pn2d_sentaurus2018 --tdr-importer build\sentaurus_import.exe --runner build\vela_example_runner.exe
python scripts\compare_sentaurus_tdr_tdx.py --tdr-export build\reference_tcad\pn2d_sentaurus2018 --tdx-dir reference_tcad\pn2d_sentaurus2018\source --output-dir build\reference_tcad\pn2d_sentaurus2018\reports
python scripts\diagnose_pn2d_0v_current_balance.py --reference-root build\reference_tcad\pn2d_sentaurus2018 --runner build\vela_example_runner.exe --output-dir build\reference_tcad\pn2d_sentaurus2018\reports\0v_current_balance --require-balanced
python scripts\diagnose_pn2d_0v_field_conventions.py --reference-root build\reference_tcad\pn2d_sentaurus2018 --runner build\vela_example_runner.exe --output-dir build\reference_tcad\pn2d_sentaurus2018\reports\0v_field_conventions --ni-cm3 1e10 --ni-cm3 1.45e10
python scripts\diagnose_pn2d_0v_field_mapping.py --reference-root build\reference_tcad\pn2d_sentaurus2018 --runner build\vela_example_runner.exe --output-dir build\reference_tcad\pn2d_sentaurus2018\reports\0v_field_mapping --fields ElectrostaticPotential:Potential --fields eQuasiFermiPotential:ElectronQuasiFermi --fields hQuasiFermiPotential:HoleQuasiFermi --fields eDensity:Electrons --fields hDensity:Holes
python scripts\compare_pn2d_0v_electric_field.py --reference-root build\reference_tcad\pn2d_sentaurus2018 --output-dir build\reference_tcad\pn2d_sentaurus2018\reports\0v_electric_field
python scripts\diagnose_pn2d_0v_density_decomposition.py --reference-root build\reference_tcad\pn2d_sentaurus2018 --runner build\vela_example_runner.exe --output-dir build\reference_tcad\pn2d_sentaurus2018\reports\0v_density_decomposition --ni-cm3 1e10 --ni-cm3 1.45e10
python scripts\diagnose_pn2d_0v_ni_bgn_probe.py --reference-root build\reference_tcad\pn2d_sentaurus2018 --runner build\vela_example_runner.exe --output-dir build\reference_tcad\pn2d_sentaurus2018\reports\0v_ni_bgn_probe --ni-cm3 1e10 --ni-cm3 1.45e10 --ni-cm3 1.6556207295e10 --bgn none --bgn slotboom
python scripts\probe_pn2d_0v_qf_drivers.py --reference-root build\reference_tcad\pn2d_sentaurus2018 --runner build\vela_example_runner.exe --output-dir build\reference_tcad\pn2d_sentaurus2018\reports\0v_qf_drivers
python scripts\compare_pn2d_0v_state.py --reference-root build\reference_tcad\pn2d_sentaurus2018 --runner build\vela_example_runner.exe --output-dir build\reference_tcad\pn2d_sentaurus2018\reports\0v_state
```

Known unsupported or approximate physics remains explicit: Sentaurus avalanche
models are not a calibrated numerical match, Fermi statistics are still treated
with Vela's current carrier-statistics implementation, and BV impact-ionization
comparison is diagnostic until the coupled Jacobian and continuation path are
fully calibrated. The executable BV deck therefore keeps strict coupled Newton
but sets `impact_ionization.model: "none"`; avalanche/impact-ionization fields
remain available through the imported Sentaurus TDR/tdx state exports.

The current executable comparison uses the faithful node-level doping decks.
The IV deck uses `vela_stop: 0.3` and `vela_step: 0.02`; the local gate compares
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
mobility override with `bandgap_narrowing: "none"`. The IV deck now declares a
separate `caughey_thomas_field` mobility override (silicon defaults), matching
the Sentaurus `DopingDep`+`HighFieldSaturation` declaration; this is the
Sentaurus-faithful model and improves the IV forward-current comparison.
The BV-only Caughey-Thomas constant point is still not reused for IV because it
degrades the forward-current comparison relative to the field-dependent point.
These reports are diagnostic-only and do not yet require trend match.

### IV high-bias physical quantities

The June 2026 high-bias IV state comparison found that the previous 1 V
diagnostic was not a same-bias comparison: Sentaurus fields came from the 1.0 V
terminal state, while the Vela probe stopped at `0.8265625 V` and then failed
validation at `0.828125 V` with `contact 'Anode' node 2 phin=0 does not match
bias`. The direct cause was Newton's default p-contact minority-electron
relaxation, which intentionally relaxed the Anode minority `phin` to 0 V at
high forward bias. That behavior is useful as a numerical option, but it is not
appropriate for the Sentaurus-faithful pn2d IV calibration state, where both
contact quasi-Fermi potentials should remain tied to the Ohmic contact voltage.

The pn2d Sentaurus2018 IV reference override now sets
`contact_boundary_minority_electron_relaxation: false`. With that override, the
1 V probe converges through `1.0000000000000002 V`; Anode and Cathode terminal
currents remain balanced, and the high-bias current ratio improves to
Vela/Sentaurus `0.826676` at 1 V. The full `0.20547013066..1.0 V` IV curve
comparison is still not a calibrated pass (`max_relative_error=1.17508`,
`orders_of_magnitude=0.337475`), mostly because the lower forward-bias points
remain under-calibrated, but the previous boundary-condition artifact is
removed.

Same-bias 1.0 V field comparison artifacts live under
`build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/fixed/physical_quantity_compare/`.
The dominant field errors after disabling relaxation are much smaller:
`electron_qf_V` max error drops from the diagnostic 1.0 V discontinuity to
`0.015289 V`; `eDensity` and `hDensity` mean absolute errors drop to about
`1.42e16 cm^-3`; and `ElectricField` p95 error drops to about `333 V/cm`.
Remaining IV mismatch is therefore a current-magnitude calibration problem
(mobility/current-density magnitude, recombination/effective-ni calibration, or
width/unit convention), not a contact-boundary failure.

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
about 0.064 orders.

## IV Mobility Promotion (2026-05-29)

A dedicated IV mobility candidate scan
(`scripts/scan_pn2d_iv_mobility_candidates.py`, summary
`build/pn2d_iv_mobility_scan/pn2d_iv_mobility_summary.csv`) compared the default
(no-mobility) IV deck against Caughey-Thomas variants while leaving the BV deck
untouched. With the imported reference local ratio `I(0.29)/I(0.30) = 0.632436`:

- `baseline` (no mobility): window orders `0.5858` (fails the 0.50 gate).
- `caughey_thomas` (silicon): window orders `0.4459`.
- `caughey_thomas_field` (silicon): window orders `0.4214` (best).
- `caughey_thomas` / `caughey_thomas_field` with BV constants: `0.4364` / `0.4340`.

All Caughey-Thomas variants pass the 0.50 IV window gate, keep the 0.3 V
terminal-current sum near the numerical floor (`<1e-18 A/um`), keep strict
Newton handoff `true`, and leave the BV 0.05 V gate unchanged at `0.1109`.
The winning `caughey_thomas_field` point was promoted into the IV deck via a
per-simulation `vela_solver.mobility` override in
`reference_tcad/pn2d/pn2d_reference.json`. Caughey-Thomas closes the IV
magnitude axis (window orders `-0.15`) but the local slope delta remains
`~0.099` (target `0.075` not yet met).

## Phase 2B: Local-Slope vs Magnitude Tradeoff (2026-05-29)

To test whether combining the promoted `caughey_thomas_field` mobility with the
anode minority-electron quasi-Fermi relaxation knob could also close the local
slope, the contact-relaxation scan
(`scripts/scan_pn2d_contact_relax_candidates.py`, summary
`build/pn2d_phase2b_scan/pn2d_contact_relax_summary.csv`) was rerun on the
promoted base deck. With reference ratio `I(0.29)/I(0.30) = 0.632436`:

- `baseline` (CT_field, no relaxation): local slope delta `0.0987`,
  IV window orders `0.4214`.
- `n_contact_only` relaxation (thresholds `0.08`-`0.20` V, all identical):
  local slope delta `0.0702` (meets the `0.075` target), but IV window orders
  rise to `0.4946`.

This confirms the magnitude and local-slope axes trade against each other with
these two levers: the n-contact-only relaxation meets the local-slope target but
erases the Caughey-Thomas magnitude gain, leaving only `0.0054` of margin below
the `0.50` IV-window gate (regression-fragile). BV (`0.1109`), the 0.3 V
terminal-current sum (`~1e-19 A/um`), and strict Newton handoff are preserved in
all candidates. Decision: keep the `caughey_thomas_field` promotion alone (with
its healthy `0.4214` magnitude margin) and carry the local slope delta `~0.099`
as a documented known gap; the n-contact-only relaxation is not promoted because
the resulting IV-window margin is too small to be regression-stable.

## M2 Strategy x Contact-Side Scan (2026-05-27)

To test a non-threshold M2 direction for the remaining pn2d high-bias slope
gap, a dedicated strategy x contact-side scan was executed with:

- Script: `scripts/scan_pn2d_contact_relax_candidates.py`
- Output root: `build/pn2d_contact_relax_scan`
- Summary files:
  - `build/pn2d_contact_relax_scan/pn2d_contact_relax_summary.csv`
  - `build/pn2d_contact_relax_scan/pn2d_contact_relax_summary.json`

Candidates:

- `baseline` (iter2 behavior)
- `dominant_p_only`, `dominant_n_only`, `dominant_both`
- `legacy_p_only`, `legacy_n_only`, `legacy_both`

Result summary:

- Baseline/`*_p_only`/`*_both` stay at
  `I(0.29)/I(0.30) ~= 0.7309387` (`delta ~= 0.098503`).
- `*_n_only` improves local ratio to `0.7025914`
  (`delta ~= 0.070156`).
- IV window orders:
  - baseline/`*_p_only`/`*_both`: `0.4662111`
  - `*_n_only`: `0.4945596` (worse than promoted baseline gate)
- Terminal current sum at 0.3 V remains near numerical floor
  (`8.5e-20` to `1.7e-19 A/um`).
- BV 0.05 V orders remain `0.1109109`.
- Strict Newton handoff remains `true` for IV/BV/fine sweeps in all cases.

Conclusion:

- `n_contact_only` is a valid improvement direction for local slope but is not
  yet promotable because it regresses IV-window orders beyond the current gate.
- `dominant` vs `legacy` reconstruction did not materially change the result in
  this matrix.
- Promoted baseline remains unchanged (`dominant_p_only` equivalent behavior).

## M2 Round2 N-Only Threshold Refinement (2026-05-27)

A dedicated threshold micro-sweep was run to test whether `n_contact_only`
becomes promotable when only the bias threshold is adjusted.

- Script: `scripts/scan_pn2d_n_only_thresholds.py`
- Output root: `build/pn2d_contact_relax_scan`
- Summary files:
  - `build/pn2d_contact_relax_scan/pn2d_contact_relax_round2_n_only_summary.csv`
  - `build/pn2d_contact_relax_scan/pn2d_contact_relax_round2_n_only_summary.json`

Candidates:

- `baseline`
- `n_only_th0p08`
- `n_only_th0p1`
- `n_only_th0p12`
- `n_only_th0p15`
- `n_only_th0p2`

Result summary:

- All `n_only_th*` thresholds collapse to the same metric point in this deck:
  - `I(0.29)/I(0.30) = 0.7025914` (`delta = 0.0701557`)
  - IV window orders (`0.2-0.3 V`) = `0.4945596`
  - Terminal sum at 0.3 V = `1.7293e-19 A/um`
  - BV orders at 0.05 V = `0.1109109`
  - Strict Newton handoff = `true`
- Baseline remains:
  - `I(0.29)/I(0.30) = 0.7309387` (`delta = 0.0985029`)
  - IV window orders (`0.2-0.3 V`) = `0.4662111`

Conclusion:

- Threshold tuning in `0.08-0.20 V` does not add discrimination for the current
  `n_contact_only` branch.
- The local slope gain persists, but IV-window regression persists as well, so
  no threshold candidate is promotable.
- Promoted baseline remains unchanged.

## M2 Round3 N-Only Edge-Threshold Refinement (2026-05-27)

To test whether thresholds near the IV comparison window change the trade-off,
the same scan was rerun with high thresholds:

- Script: `scripts/scan_pn2d_n_only_thresholds.py`
- Command:
  - `python scripts/scan_pn2d_n_only_thresholds.py --thresholds 0.24,0.26,0.28,0.29,0.295 --summary-prefix pn2d_contact_relax_round3_n_only_edge_summary --candidate-dirname candidates_round3`
- Summary files:
  - `build/pn2d_contact_relax_scan/pn2d_contact_relax_round3_n_only_edge_summary.csv`
  - `build/pn2d_contact_relax_scan/pn2d_contact_relax_round3_n_only_edge_summary.json`

Candidates:

- `baseline`
- `n_only_th0p24`
- `n_only_th0p26`
- `n_only_th0p28`
- `n_only_th0p29`
- `n_only_th0p295`

Result summary:

- All `n_only_th*` points remain identical:
  - `I(0.29)/I(0.30) = 0.7025914` (`delta = 0.0701557`)
  - IV window orders (`0.2-0.3 V`) = `0.4945596`
  - Terminal sum at 0.3 V = `1.7293e-19 A/um`
  - BV orders at 0.05 V = `0.1109109`
  - Strict Newton handoff = `true`
- Baseline remains:
  - `I(0.29)/I(0.30) = 0.7309387` (`delta = 0.0985029`)
  - IV window orders (`0.2-0.3 V`) = `0.4662111`

Conclusion:

- Raising the threshold up to `0.295 V` does not recover IV-window quality.
- The threshold axis is exhausted for this branch; future M2 progress needs a
  non-threshold boundary mechanism.

## M2 Round4 N-Only Strength Refinement (2026-05-27)

To test a continuous non-threshold relaxation axis, the scan was rerun with a
fixed `0.1 V` activation threshold and varying minority-relaxation strength:

- Script: `scripts/scan_pn2d_n_only_thresholds.py`
- Command:
  - `python scripts/scan_pn2d_n_only_thresholds.py --strengths 0.0,0.25,0.5,0.75,1.0 --strength-threshold 0.1 --summary-prefix pn2d_contact_relax_round4_n_only_strength_summary --candidate-dirname candidates_round4`
- Summary files:
  - `build/pn2d_contact_relax_scan/pn2d_contact_relax_round4_n_only_strength_summary.csv`
  - `build/pn2d_contact_relax_scan/pn2d_contact_relax_round4_n_only_strength_summary.json`

Candidates:

- `baseline`
- `n_only_str0p0`
- `n_only_str0p25`
- `n_only_str0p5`
- `n_only_str0p75`
- `n_only_str1p0`

Result summary:

- All `n_only_str*` points remain identical:
  - `I(0.29)/I(0.30) = 0.7025914` (`delta = 0.0701557`)
  - IV window orders (`0.2-0.3 V`) = `0.4945596`
  - Terminal sum at 0.3 V = `1.7293e-19 A/um`
  - BV orders at 0.05 V = `0.1109109`
  - Strict Newton handoff = `true`
- Baseline remains:
  - `I(0.29)/I(0.30) = 0.7309387` (`delta = 0.0985029`)
  - IV window orders (`0.2-0.3 V`) = `0.4662111`

Conclusion:

- The strength axis does not separate candidate behavior in this deck.
- Threshold and strength sweeps now both collapse to the same n-only response.
- Further progress requires a different non-threshold boundary mechanism.

## M3 Cross-Device Generalization Snapshot (2026-05-27)

Generalization checks were executed on the same HEAD after M2 integration:

- Full preset quality gate:
  - `ctest --preset windows-ucrt64-debug`
  - `274/274` passed, `0` failed.
- Regression matrix artifact:
  - `build/regression_output/regression_summary.json`

Representative cross-device rows (all passed with converged sweeps):

| case | rows | converged_rows | final_current_total |
| --- | ---: | ---: | ---: |
| `pn_diode_iv` | 3 | 3 | `3.481904227298678e-03` |
| `pn_diode_bv` | 3 | 3 | `-5.78531939292131e-11` |
| `nmos2d_dd_iv` | 3 | 3 | `5.009164630732429` |
| `pmos2d_dd_iv` | 3 | 3 | `-1.7810363137092906` |
| `nmos2d_mos_dd_iv` | 3 | 3 | `1.62055767856109e-06` |
| `pmos2d_mos_dd_iv` | 3 | 3 | `-5.77996479489852e-07` |
| `nmos2d_mos_dd_bv` | 3 | 3 | `1.72038242915746e-07` |
| `pmos2d_mos_dd_bv` | 3 | 3 | `-6.3310014841849e-08` |

Interpretation:

- The pn2d promoted baseline and M2 parameter plumbing did not introduce
  regressions in the broader NMOS/PMOS and mixed-material smoke matrix.
- Remaining open item is pn2d high-bias local slope closure, which requires a
  new tuning direction beyond threshold-only contact-relaxation scans.

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

## IV Resolution Scan (Task 1)

Generated with `scripts/scan_pn2d_iv_resolution.ps1` from
`build/pn2d_tdr_tie_probe`:

| Case | Step (V) | Accepted rows (`handoff_stage=newton`) | Orders | Max relative error | Ratio at 0.25 V | Ratio at 0.27 V | Ratio at 0.29 V | Gate impact |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| promoted | 0.10 | 4 | 0.5048 | 0.6872 | 0.8499 | 0.5348 | 0.3128 | baseline |
| step0p02 | 0.02 | 16 | 0.4979 | 0.6822 | 0.4385 | 0.3683 | 0.3178 | no material IV gain |
| step0p01 | 0.01 | 31 | 0.5207 | 0.6985 | 0.4194 | 0.3508 | 0.3015 | no material IV gain |

Resolution takeaway:

- Finer IV steps increase accepted strict-Newton rows but do not improve IV
  orders by the `>= 0.15` decision threshold.
- The high-forward-bias slope mismatch remains after interpolation pressure
  reduction, so the IV gate remains unchanged.

## SRH Lifetime Sweep In Existing Physics Matrix (Task 2)

Extended matrix generated with `scripts/scan_pn2d_iv_bv_physics_matrix.ps1`.

| Case | Kind | Orders | Ratio at target | Gate impact |
| --- | --- | ---: | ---: | --- |
| default | IV | 0.5048 | 0.3128 at 0.29 V | current IV baseline |
| iv_srh_tau1e-6 | IV | 0.9417 | 0.1144 at 0.29 V | IV worsened |
| iv_srh_tau1e-8 | IV | 0.8338 | 1.5467 at 0.29 V | IV worsened |
| bv_recomb_none | BV | 0.0641 | 0.8628 at 0.05 V | current promoted BV gate |
| bv_recomb_srh | BV | 1.1778 | 15.0580 at 0.05 V | BV parity fails |
| bv_srh_tau1e-6 | BV | 0.2650 | 1.8406 at 0.05 V | BV improves vs default SRH, still worse than promoted gate |
| bv_srh_tau1e-8 | BV | 2.1721 | 148.6262 at 0.05 V | BV severely worsened |

SRH sweep takeaway:

- `taun=taup=1e-6` reduces the SRH-enabled BV mismatch from about `1.18` orders
  to about `0.26` orders, but it still does not outperform the promoted
  recombination-disabled BV gate (`0.064`).
- The same SRH lifetime change degrades IV materially, so SRH lifetime alone is
  not a viable combined IV/BV promotion candidate.
- Based on this decision check, the next step is exposing Auger coefficients in
  solver JSON (Task 3) for controlled recombination parameter probes.

## Auger JSON Surface (Task 3)

The solver config surface now accepts the Auger coefficient keys in both
Gummel and Newton paths:

- `auger_cn_m6_per_s`
- `auger_cp_m6_per_s`

Coverage status:

- Gummel JSON parsing coverage added in `tests/test_mobility.cpp`.
- Newton JSON parsing coverage added in `tests/test_newton_solver.cpp`.
- Negative-coefficient guard remains enforced by `RecombinationModel`
  validation (`tests/test_recombination.cpp`).

Gate impact:

- No IV/BV gate values are promoted by Task 3 itself; this is a config-surface
  prerequisite for the next recombination matrix probe.

## Auger And Combined Recombination Matrix (Task 4)

The physics matrix now includes:

- SRH+Auger with Auger coefficients scaled down (`0.5x`) and up (`2.0x`)
- Auger-only cases (`recombination: ["auger"]`) for SRH-decoupled behavior

Measured with the same summary columns and target biases as prior matrix runs.

| Case | Kind | Orders | Ratio at target | Gate impact |
| --- | --- | ---: | ---: | --- |
| iv_srh_auger_half | IV | 1.7843 | 0.0164 at 0.29 V | IV strongly worsened |
| iv_srh_auger_double | IV | 0.5031 | 0.3139 at 0.29 V | essentially baseline IV |
| iv_auger_only | IV | 0.8929 | 0.1287 at 0.29 V | IV worsened |
| bv_srh_auger_half | BV | 1.1701 | 14.7956 at 0.05 V | no meaningful BV recovery |
| bv_srh_auger_double | BV | 1.1752 | 14.9681 at 0.05 V | no meaningful BV recovery |
| bv_auger_only | BV | 0.2801 | 0.4753 at 0.05 V | improved vs SRH+Auger, still worse than promoted BV gate |

Task 4 takeaway:

- Within SRH+Auger runs, scaling Auger coefficients over `0.5x` to `2.0x` does
  not materially change the BV mismatch (still about `1.17` orders).
- Auger-only improves BV compared with SRH+Auger, but it degrades IV and still
  underperforms the promoted BV gate (`recombination: ["none"]`, `0.064`
  orders).
- No recombination parameter set in this matrix satisfies both IV and BV
  criteria together, so recombination tuning remains insufficient as a unified
  explanation for the IV slope and BV SRH jump.

## Narrow Recombination Diagnostic (Task 5)

`DCSweep` now supports an opt-in per-bias recombination diagnostic path under
existing `solver.diagnostics` for `newton` and `gummel_newton` sweeps. When
enabled, CSV output appends:

- `recombination_max_abs_rate_m3_per_s`
- `recombination_mean_abs_rate_m3_per_s`
- `carrier_product_max_np_over_ni2`

The diagnostic is disabled by default and verified by
`tests/test_dc_sweep.cpp` to be both opt-in and finite.

### BV 0.05 V Diagnostic Comparison

Using diagnostic-enabled reruns at `0.05 V`:

| Case | total_A_per_um | max \|R\| (m^-3 s^-1) | mean \|R\| (m^-3 s^-1) | max np/ni^2 |
| --- | ---: | ---: | ---: | ---: |
| `bv_recomb_none` | 9.245665e-19 | 0.0 | 0.0 | 13.8270 |
| `bv_recomb_srh` | 1.613673e-17 | 6.304272e+21 | 4.743757e+20 | 6.3793 |
| `bv_srh_tau1e-6` | 1.972474e-18 | 6.304274e+20 | 4.743764e+19 | 5.4077 |
| `bv_recomb_srh_auger` | 1.598583e-17 | 6.304273e+21 | 4.743763e+20 | 5.5690 |
| `bv_auger_only` | 5.623094e-19 | 1.232900e+13 | 8.291935e+11 | 5.4032 |

Task 5 takeaway:

- The BV SRH jump tracks very large SRH recombination magnitude (`~1e21`) at
  the target bias.
- Increasing SRH lifetime from `1e-7` to `1e-6` lowers recombination magnitude
  by about one order and correspondingly reduces the BV current jump.
- Auger-only recombination stays many orders smaller than SRH at this bias,
  and SRH+Auger is nearly identical to SRH-only, so Auger is not the dominant
  source of the BV jump in this window.
- This points to SRH lifetime/model parity as the primary sensitivity axis;
  no local numerical instability signal (NaN/Inf) is observed in the new
  diagnostics.

## Task 6: Fermi Statistics Decision

Degeneracy relevance estimate for pn2d (`Na=Nd=1e17 cm^-3`, Si defaults
`Nc=2.8e19 cm^-3`, `Nv=1.04e19 cm^-3`):

- `n/Nc ~= 3.57e-3`
- `p/Nv ~= 9.62e-3`
- `Ec-Ef ~= 0.1457 eV` (n side, majority estimate)
- `Ef-Ev ~= 0.1201 eV` (p side, majority estimate)

These values remain comfortably in the non-degenerate regime for the baseline
doping scale. Combined with Task 5 diagnostics (`max np/ni^2` in the tested BV
window staying single-digit to low-tens), there is no direct indication that
missing Fermi-Dirac carrier statistics is the dominant cause of the current BV
low-bias mismatch.

Decision for this plan stage:

- IV high-bias slope mismatch remains open after Tasks 1-4.
- BV SRH jump is now localized primarily to SRH recombination sensitivity
  (lifetime/model parity axis), not Auger dominance and not a numerical
  instability signature.
- Therefore, a broad Fermi-statistics implementation is **not** started in this
  plan stage, and a separate Fermi implementation plan is not yet justified by
  current evidence.

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

## HEAD Rerun After Contact-Current QF Alignment

Rerun date: 2026-05-26 on `3e84042`.

The latest `ContactCurrent` change makes uniform-`ni` contact edges use the same
QF SG branch as `CoupledDDAssembler`. A direct cathode probe at IV 0.30 V shows:

| Case | Density total (A/um) | QF total (A/um) | Assembler-style residual total (A/um) |
| --- | ---: | ---: | ---: |
| default Slotboom | 9.786651e-15 | -1.897805e-14 | -1.897805e-14 |
| BGN off | 2.359720e-14 | 4.718129e-16 | 4.718129e-16 |

Takeaway:

- The default IV cathode CSV already matches the assembler-consistent residual,
  so the remaining IV shortfall is not a cathode current-extraction artifact.
- Older BGN-off and promoted BV numbers that relied on density-form extraction
  are stale under current HEAD.
- Sentaurus IV/BV logs report no lifetime file, no model-parameter file, and no
  field-, doping-, or temperature-dependent SRH lifetimes. The SRH parity task
  should target Sentaurus built-in constant lifetime/trap defaults first, not
  Scharfetter doping-dependent lifetime support.
- Equal-lifetime SRH scan:
  - IV target ratio is closest between `taun=taup=3e-8 s` (`0.666x`) and
    `1e-8 s` (`1.704x`).
  - BV target ratio is closest at `taun=taup=3e-6 s` (`0.935x`).
  - No single equal lifetime matches both IV and BV.
- Tested asymmetric pairs also do not match both; all BV target ratios remain
  `3.45x` or worse.

Artifacts:

- `build/pn2d_tdr_tie_probe/vela/pn2d_taugrid_summary.csv`
- `build/pn2d_tdr_tie_probe/vela/pn2d_taupair_summary.csv`

Useful local verification command:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\sentaurus_import.py reference --config reference_tcad\pn2d\pn2d_reference.json --source-dir reference_tcad\pn2d --output-dir build\reference_tcad\pn2d --tdr-importer build\sentaurus_import.exe --runner build\vela_example_runner.exe
```

## 2026-05-27 Root-Cause Execution Update (Tasks 1-8)

Fresh baseline regeneration was completed in `build/pn2d_root_cause_probe`.

Baseline metrics:

- IV comparison orders: `0.49456` (`reports/pn2d_iv_comparison.json`).
- BV comparison orders: `0.17999` (`reports/pn2d_bv_comparison.json`).
- IV and BV accepted rows remained strict Newton handoff (`handoff_stage=newton`, `newton_iterations>0`).

### Excluded/Constrained Hypotheses

- Unit scaling bug (`A/m` to `A/um`) was not supported by data.
  At `0.3 V`, `current_total / current_total_A_per_um = 1e6` exactly.
- Contact geometry weighting mismatch was not supported by first-order perimeter/couple checks.
  Contact boundary `sum_couple_um` was `0.5` for both electrodes.
- However, two-terminal extraction still showed non-zero `Ic + Ia` residual at forward bias
  (see `build/pn2d_root_cause_probe/reports/task2_current_geometry_audit.csv`), so
  contact extraction semantics remain an open IV root-cause candidate.

### Compensated Doping Policy A/B

`dominant_signed_region` remained the better policy for pn2d tie nodes (33 compensated nodes).

- `dominant_signed_region`: IV `0.49456`, BV `0.17999`.
- `reported`: IV `0.77727`, BV `1.39817` (comparison gate fail).

Artifact: `build/pn2d_policy_ab_summary.csv`.

### Field Profile Comparison Boundary

Sentaurus field export from `pn2d_des.tdr` is the final `1.0 V` state. Vela nodewise comparison
used `0.3 V` VTK, so the profile comparison is trend-only and not strict `0.3 V` pointwise proof.

Artifact: `build/pn2d_field_compare/field_profile_error_summary.md`.

### Mobility/Flux Evidence

IV mobility variant runs showed that raising mobility (for example constant mobility) can increase
high-bias current, but worsens combined 0.2-0.3 V match quality; this is not a clean single-factor
global fix.

Artifacts:

- `build/pn2d_root_cause_probe/vela/pn2d_iv_mobility_matrix_summary.csv`
- `build/pn2d_root_cause_probe/vela/pn2d_iv_edge_mobility_stats.csv`

### BV Recombination Matrix and Confirmed Sensitivity

At `0.05 V`, BV current sensitivity was dominated by SRH lifetime settings.
Matrix evidence (`build/pn2d_root_cause_probe/vela/pn2d_bv_recomb_matrix_summary.csv`) showed:

- `recombination=[none]`: ratio `0.6607`.
- `recombination=[srh], taun=taup=1e-6`: ratio `1.2910` (closest among tested SRH points in this run).
- `recombination=[srh], taun=taup=1e-8`: ratio `147.63` (strong over-current).

## Minimal Fix Implemented (Task 7)

Chosen confirmed root cause axis: BV SRH parameter/config parity.

Changes:

- `reference_tcad/pn2d/pn2d_reference.json` (BV `vela_solver`):
  - `recombination: ["none"] -> ["srh"]`
  - add `taun: 1.0e-6`, `taup: 1.0e-6`
- `tests/regression/test_sentaurus_sample_integration.py`:
  - expect BV deck `recombination == ["srh"]`
  - assert `taun` and `taup` are `1e-6`
  - tighten BV comparison gate in test from `< 0.20` to `< 0.15`

Failure-first evidence:

- Before config update, targeted test failed on expected `['srh']` vs actual `['none']`.

Post-fix evidence:

- Targeted regression `sentaurus_sample_integration` passed.
- Fresh import in `build/pn2d_task7_after` produced:
  - IV orders: `0.49456` (unchanged)
  - BV orders: `0.11091` (improved from `0.17999`)

## Final Validation Status (Task 8)

- Build: `cmake --build --preset windows-ucrt64-debug` completed (`ninja: no work to do`).
- Full test preset: `ctest --preset windows-ucrt64-debug --output-on-failure` executed.
  Summary reported `99% tests passed`, `3` failures in total (pn2d change-specific sentaurus integration passed).
- Required narrow sets were exercised and passed in this run output:
  `sentaurus_sample_integration`, `reference_tcad_regression`, `dc_sweep`, `ascii_sources`.

Remaining limitation:

- IV residual root cause is still open; this change only addresses the confirmed BV SRH sensitivity axis.

## 2026-05-27 Follow-up Closure (Residual Full-Suite Failures)

After Task 8, three unrelated full-suite failures remained in this workspace run
(`122`, `262`, `270`). They are now closed with minimal-scope fixes:

- PMOS regression polarity alignment:
  - `examples/pmos2d_dd/simulation_iv.json`: `drain_current_sign` set to `-1.0`
  - `examples/pmos2d_mos_dd/simulation_iv.json`: `drain_current_sign` set to `-1.0`
- BV finite-output checker robustness:
  - `scripts/run_regression.py` now applies CSV-aware finite checks and allows
    `inf` only for BV diagnostic column `current_jump_ratio`.

Verification:

- `ctest --preset windows-ucrt64-debug -I 122,122 --output-on-failure`: pass
- `ctest --preset windows-ucrt64-debug -I 262,262 --output-on-failure`: pass
- `ctest --preset windows-ucrt64-debug -I 270,270 --output-on-failure`: pass
- `ctest --preset windows-ucrt64-debug --output-on-failure`: **100% tests passed (273/273)**

## 2026-05-27 Next-Phase Root-Cause Run (G/H/I)

This follow-up executed the next planned IV residual investigations with currently
available artifacts.

### G. Geometry / Unit-Factor Audit (refresh)

Reused current geometry audit data confirms no first-order A/m to A/um scaling mismatch:

- At `0.3 V`, `scale_A_per_m_over_A_per_um = 1000000.0`.
- Contact dual-length totals stay matched: `sum_couple_um = 0.5` for both Cathode and Anode.

Artifact: `build/pn2d_root_cause_probe/reports/task2_current_geometry_audit.csv`.

Additional branch probe at `0.3 V` on `dc_sweep_0003_0.3V.vtk`:

- Cathode `assembler_residual_total = -1.57935e-14 A/um` (matches the qf branch, not density branch).
- Anode `assembler_residual_total = +2.13545e-14 A/um` (qf and density nearly identical in this case).

Artifacts:

- `build/pn2d_root_cause_probe/reports/taskG_branch_compare_cathode_0p3V_20260527.csv`
- `build/pn2d_root_cause_probe/reports/taskG_branch_compare_anode_0p3V_20260527.csv`

### H. Mobility / Flux-Coefficient Evidence (refresh)

Re-ran contact-side decomposition summary and confirmed persistent electron drift/diffusion cancellation
at forward bias, with growing cancellation ratio at high bias:

- default case at `0.3 V`: `I_e=-2.876e-14`, `I_e_drift=+4.985e-14`, `I_e_diff=-2.109e-14`,
  `|I_e_drift|/|I_e|=1.73`.
- recombination-off case at `0.3 V`: `|I_e_drift|/|I_e|=1.99`.

This supports that the remaining IV gap is not removed by simple recombination toggles and retains
strong sensitivity to transport-state details near contact-adjacent edges.

Artifact: `build/pn2d_root_cause_probe/reports/taskH_contact_decomposition_20260527.txt`.

Bias-aligned (`1.0 V`) contact-edge follow-up was also executed using node-mapped Sentaurus exports
from `pn2d_forward_des.tdr` and Vela `dc_sweep_0001_1V.vtk`:

- Edge report:
  `build/pn2d_root_cause_probe/reports/taskH_contact_edge_state_compare_1p0V_20260527.csv`
- Summary:
  `build/pn2d_root_cause_probe/reports/taskH_contact_edge_state_compare_1p0V_20260527_summary.json`

Key edge-level findings:

- Contact-adjacent electric-potential drop mismatch stays small:
  `abs_err_dpsi_V` mean `0.0238 V` (p95 `0.0587 V`).
- Hole quasi-Fermi edge-drop mismatch is also small:
  `abs_err_defp_V` mean `0.0194 V` (max `0.0749 V`).
- Electron quasi-Fermi edge-drop mismatch is strongly contact-asymmetric:
  overall `abs_err_defn_V` mean `0.716 V`, but by contact:
  - Cathode edges mean `0.0233 V`
  - Anode edges mean `0.898 V`

This reinforces that the residual mismatch is localized near anode-side electron-QF boundary behavior,
not a whole-domain potential/geometry scaling mismatch.

### I. Quasi-Fermi Profile Check (bias-aligned completed)

Per-simulation Sentaurus neutral field export is now enabled in the reference flow, producing
`sim_fields/<simulation>/...` artifacts (including potential and quasi-Fermi fields).

Implementation + verification:

- `scripts/sentaurus_import.py` now exports simulation TDR fields to
  `sim_fields/iv` and `sim_fields/bv` during `reference` import.
- A dedicated probe script was added:
  `scripts/probe_pn2d_quasifermi_profile_compare.py`.
- New import run generated QF/potential fields:
  `build/pn2d_taskI_qf/sim_fields/iv/fields/eQuasiFermiPotential_region{0,1}.csv`,
  `.../hQuasiFermiPotential_region{0,1}.csv`, `.../ElectrostaticPotential_region{0,1}.csv`.

Trend-only comparison versus Vela `0.3 V` VTK (node-mapped):

- Report: `build/pn2d_root_cause_probe/reports/taskI_qf_profile_compare_20260527.json`.
- `ElectrostaticPotential` abs-diff: mean `0.379 V`, p95 `0.656 V`.
- `eQuasiFermiPotential` abs-diff: mean `0.092 V`, p95 `0.144 V`.
- `hQuasiFermiPotential` abs-diff: mean `0.631 V`, p95 `0.674 V`.

Important boundary:

- The imported Sentaurus IV TDR state is the final quasistationary point (typically `1.0 V`),
  while the compared Vela state is `0.3 V`; this run therefore closes Task I as a trend check,
  not a strict bias-aligned parity proof.

Strict bias-aligned follow-up (`1.0 V`) was then executed:

- Sentaurus state source: `reference_tcad/pn2d/pn2d_forward_des.tdr`
  exported to `build/pn2d_taskI_bias_align/forward_fields`.
- Vela aligned run: `build/pn2d_taskI_bias_align/simulation_iv_1p0.json`
  producing `build/pn2d_taskI_bias_align/dc_sweep_0001_1V.vtk`.
- Bias-aligned report:
  `build/pn2d_root_cause_probe/reports/taskI_qf_profile_compare_bias_aligned_1p0V_20260527.json`.

Bias-aligned metrics (node-mapped):

- `ElectrostaticPotential` abs-diff: mean `0.00869 V`, p95 `0.0464 V`.
- `eQuasiFermiPotential` abs-diff: mean `0.0312 V`, p95 `0.0780 V`, max `0.9315 V`.
- `hQuasiFermiPotential` abs-diff: mean `0.0112 V`, p95 `0.0387 V`.

Outlier localization (`eQuasiFermiPotential`):

- Largest mismatches cluster near the high-bias anode-side boundary (for example node `876`
  at approximately `(x,y)=(1.0, 0.0)` um, abs diff `0.9315 V`).

Interpretation:

- Most of the domain shows close potential/QF agreement at matched bias; remaining large discrepancy
  is concentrated in a small boundary-adjacent subset and is consistent with contact/boundary
  treatment differences rather than a whole-domain quasi-Fermi mismatch.

### P1. SRH Default Lifetime Extraction Status

A repository-local evidence pass was completed for SRH default lifetime lookup.

- Evidence report:
  `build/pn2d_root_cause_probe/reports/taskP1_srh_default_evidence_20260527.md`
- Checked local Sentaurus artifacts:
  `reference_tcad/pn2d/pn2d_sdevice.cmd`,
  `reference_tcad/pn2d/pn2d_bv_sdevice.cmd`,
  `reference_tcad/pn2d/pn2d_des.log`,
  `reference_tcad/pn2d/pn2d_bv_des.log`.

Confirmed locally:

- `no Lifetime file`, `no ModelParameters file` in both IV and BV logs.
- SRH is active and configured without field/doping/temperature-dependent lifetimes.
- No explicit numeric built-in SRH lifetime value is printed in available local logs/cmd files.

Status:

- P1 numeric default extraction remains blocked by missing external model-parameter source.
- To close P1, we need either Sentaurus documentation for this model set or a Sentaurus run mode
  that prints resolved SRH numeric lifetime parameters.

#### P1 follow-up: constant-lifetime parity at the BV gate (2026-05-29)

Because the BV deck declares plain `Recombination(SRH ...)` with no
`DopingDependence` keyword, the Sentaurus SRH lifetime is a *constant* equal to
the silicon-default `taumax` (documented `1.0e-5 s`), not a doping-dependent
Scharfetter relation. A direct parity sweep at the validated 0.05 V BV gate
(`scripts/probe_pn2d_bv_srh_lifetime.py`, summary
`build/pn2d_bv_srh_probe/pn2d_bv_srh_summary.csv`) shows:

| `taun=taup` (s) | source | BV orders @ 0.05 V | I(0.05 V) A/µm |
| --- | --- | --- | --- |
| `1.0e-6` | current BV deck | `0.1109` (best) | `1.38e-18` |
| `3.0e-6` | prior tau-grid "best" | `0.4562` | `3.75e-19` |
| `1.0e-5` | Sentaurus silicon default `taumax` | `0.4816` | `-3.54e-19` |

Conclusion: at the validated low-reverse-bias gate the current `1.0e-6 s`
lifetime already gives the best constant-lifetime parity; moving toward the
documented `1.0e-5 s` default *degrades* the gate to `~0.48` orders and even
flips the current sign. So the residual BV mismatch is **not** a constant-SRH
lifetime parity issue at the low-bias gate — it lives in the higher-reverse-bias
avalanche/continuation regime (Newton reverse-bias robustness with the
impact-ionization Jacobian), which is a separate solver-robustness track. The
current BV deck lifetime is therefore left unchanged; no solver or deck change
is warranted from this research pass.

## 2026-05-30 IV Residual Gap: Localization And Development Plan

This section continues the pn2d IV retrospective and turns current evidence
into an execution plan for the remaining forward-bias mismatch.

### Consolidated Evidence (What Is Already Excluded)

- Unit scaling mismatch is excluded (`A/m` vs `A/um` ratio consistently `1e6`).
- Resolution pressure is excluded as a primary cause (finer IV steps did not
  materially improve window orders).
- Single-toggle physics fixes are excluded for IV slope closure in the current
  matrix (`recombination none`, BGN off, and SRH/Auger coefficient-only scans
  do not produce a robust combined IV/BV improvement).
- Whole-domain potential mismatch is excluded at bias-aligned checks (`1.0 V`);
  dominant residuals are boundary-localized.

### Working Root-Cause Hypothesis (IV)

The remaining IV gap is a boundary-localized transport-balance error in the
high-forward-bias window (`>=0.25 V`), dominated by contact-adjacent electron
quasi-Fermi handling and drift/diffusion cancellation sensitivity. Current
evidence points to a mechanism mismatch near the anode-side electron QF
boundary treatment rather than a domain-wide Poisson/DD mismatch.

### Localization Plan (Priority-Ordered)

1. Contact-edge discrete operator parity (highest priority)
   - Build a per-edge parity table between contact current assembly and coupled
     DD assembly at matched bias points (`0.27`, `0.29`, `0.30`, and `1.0 V`).
   - Add edge-level diagnostics for SG coefficients, Bernoulli arguments,
     reconstructed carrier values, and selected branch IDs.
   - Decision gate: identify one or more stable edge-level terms with
     consistent signed bias-correlated drift against Sentaurus trend.

2. Anode-side electron-QF boundary mechanism A/B
   - Isolate anode-side electron QF boundary behavior with strict minimal A/B
     switches (one mechanism change per run).
   - Keep cathode extraction and mobility fixed to current promoted IV baseline
     to avoid confounding variables.
   - Decision gate: recover local slope ratio `I(0.29)/I(0.30)` toward target
     without regressing window orders.

3. Cancellation robustness audit under controlled perturbations
   - Run small deterministic perturbations on boundary-state reconstruction and
     verify monotonicity/conditioning of electron drift and diffusion residuals.
   - Decision gate: remove non-physical high sensitivity where small
     reconstruction perturbations cause outsized net-current changes.

4. Limited cross-device sanity replay
   - Replay only a narrow NMOS/PMOS smoke subset after each candidate promote,
     then run full preset only on merge candidate.
   - Decision gate: zero new regressions in targeted matrix before full-suite.

### Development Plan (Three Iterations)

Iteration A: Instrumentation And Repro Stabilization (1-2 days)

- Add an opt-in IV boundary diagnostic block that exports per-contact-edge:
  drift/diff split, Bernoulli arguments, edge normal/current sign, and selected
  carrier reconstruction branch.
- Add deterministic repro scripts for the four anchor biases and produce one
  compact summary artifact for each run.
- Add/extend tests to assert diagnostic schema stability and finite outputs.

Exit criteria:

- Repeated runs produce stable edge-level signatures (noise below decision
  threshold), and suspected mismatch edges are consistently identified.

Iteration B: Mechanism A/B And Candidate Narrowing (2-4 days)

- Implement two to three minimal anode electron-QF boundary variants behind
  explicit config flags (default behavior unchanged).
- Compare candidates on:
  - IV window orders (`0.2-0.3 V`)
  - local slope ratio delta (`I(0.29)/I(0.30)`)
  - strict Newton handoff ownership
  - BV low-bias gate non-regression (`0.05 V` check)

Exit criteria:

- At least one candidate improves local slope meaningfully while preserving IV
  window margin and BV gate.

Iteration C: Promotion, Hardening, And Full Validation (1-2 days)

- Promote the best candidate into the pn2d reference deck and keep all other
  provisional toggles off.
- Add regression expectations for the newly improved IV metric and preserve
  existing BV guardrails.
- Run full preset tests and regression suite; archive before/after summaries.

Exit criteria:

- Promoted candidate passes strict Newton ownership checks, maintains BV gate,
  and improves IV local-slope behavior with measurable margin.

### Metrics And Promotion Gates

Mandatory gates for candidate promotion:

- IV window orders (`0.2-0.3 V`) must not regress beyond current promoted
  baseline margin.
- IV local slope delta (via `I(0.29)/I(0.30)`) must improve relative to current
  promoted baseline.
- BV `0.05 V` order gate must remain non-regressed.
- Accepted rows must remain strict Newton handoff (`handoff_stage=newton`,
  `newton_iterations>0`).
- Contact current sum near numerical floor at target forward bias remains
  bounded.

### Immediate Next Actions

1. ~~Implement the new opt-in IV contact-edge diagnostic payload.~~ **Done (2026-05-30)**
2. Generate the four-bias parity bundle and lock a baseline signature.
3. Start anode electron-QF boundary A/B with minimal isolated switches.
4. Promote only candidates that satisfy both slope and margin gates.

## Iteration A Completion (2026-05-30)

The opt-in IV contact-edge diagnostic payload was implemented and verified.

### Changes

- `include/vela/post/ContactCurrent.h`: Added `ContactCurrentEdgeDiagnostic` and
  `ContactCurrentDetailedResult` structs; added `computeDetailed()` method to the
  `ContactCurrent` class. Existing `compute()` is a thin wrapper over it, so all
  existing callers remain unaffected.
- `src/post/ContactCurrent.cpp`: Refactored the per-edge loop to emit a
  `ContactCurrentEdgeDiagnostic` per contact-adjacent edge. Fields exported:
  edge geometry (`edge_id`, `node0/1`, `length`, `couple`), orientation
  (`outward_sign`), Bernoulli arguments (`u`, `B+`, `B-`), branch selection
  (`electron_branch: quasi_fermi | density`, `hole_branch`), and all six current
  components (drift, diffusion, total) per carrier.
- `include/vela/simulation/DCSweep.h`: Added `ContactEdgeDiagnosticsConfig` and
  `SweepDiagnosticsConfig` structs; added `diagnostics` field to `DCSweepConfig`.
- `src/simulation/DCSweep.cpp`: Added parsing of `sweep.diagnostics.contact_edge`
  JSON block; if `enabled: true`, a separate CSV (`_contact_edges.csv` by default
  or a user-specified `csv_file`) is written with one row per contact-adjacent edge
  per converged bias point. The per-micron column is emitted when `unit_scaling` is
  active. Default behavior (diagnostics disabled) is unchanged.
- `tests/test_dc_sweep.cpp`: Added `[contact_edge]` test case asserting that the
  diagnostic file is absent when not configured, the CSV schema is correct when
  enabled, and any data rows contain finite values and valid branch labels.

### Validation

- `ctest --preset windows-ucrt64-debug --output-on-failure`: **275/275 passed**.
- No regression in IV/BV gate values, existing CSV schema, or handoff provenance
  columns.

### Next Step (Iteration B)

Generate the four anchor-bias parity bundle (0.27, 0.29, 0.30, 1.0 V) using the
new contact-edge diagnostic and identify specific anode-side edges with persistent
drift/diffusion imbalance relative to Sentaurus. Use this signature to design the
minimal anode electron-QF boundary A/B variants.

---

## Iteration B Completion (2026-05-30)

Contact-edge diagnostics were used to characterize the pn2d IV mismatch and a
config-only tau-lifetime sweep identified the primary cause of the ideality-factor
divergence.

### Diagnostic Setup

Additional scripts created (not committed; generated artefacts live in
`build/pn2d_probe/`):

- `scripts/probe_anode_analysis.py`: Cathode vs anode contact-edge comparison at
  anchor biases (0.27, 0.29, 0.30 V).
- `scripts/ideality_compare.py`: Per-step ideality factor n and orders-gap vs bias.
- `scripts/analyze_iter_b.py`: Comparative IV/ideality table across four tau variants.
- `build/pn2d_probe/vela/make_iter_b_configs.py`: Config generator for B1/B2/B3 variants.

### Bug Fix: DCSweep.cpp — contact-edge write outside VTK block

The contact-edge CSV writer block in `src/simulation/DCSweep.cpp` was
accidentally nested inside `if (converged && sweep.writeVtk)`. Because
`write_vtk: false` in probe configs, the CSV always had only a header row.
Fixed by moving the contact-edge write block before the VTK block as an
independent `if (converged && contactEdgeCsv != nullptr)` check.
Post-fix: 2/2 dc_sweep tests pass; 275/275 total.

### Anode vs Cathode Edge Diagnostics (0.30 V)

| Contact  | Edges | I_total (A/m) | e_drift  | e_diff   | h_drift   | h_diff   | drift_frac |
|----------|-------|---------------|----------|----------|-----------|----------|------------|
| Cathode  | 3     | −1.72e−08     | +2.85e−8 | −2.11e−8 | +2.4e−18  | +9.79e−9 | 0.575      |
| Anode    | 10    | +1.72e−08     | −5e−16   | +5e−11   | −1.61e−7  | +1.44e−7 | ≈0.000     |

**Cathode interpretation**: Dominant electron SG current (drift−diffusion balance)
plus significant minority-hole diffusion (9.79e-9 A/m = 57% of total). This is
correct for a **short diode** regime: with `taup=1e-7 s`, `Lp ≈ 11 μm >> L_n=1 μm`,
so minority holes traverse the n-region nearly intact and recombine at the cathode.

**Anode interpretation**: Hole current (−1.61e−7 drift + 1.44e−7 diffusion) carries
essentially 100% of the anode current. Electron contribution to anode current is
negligible (~5e-11 A/m vs 1.72e-8 A/m total). This is physically correct for a
p-type ohmic contact at forward bias.

**Key finding**: Neither contact shows an abnormal SG branch label or extreme
Bernoulli argument (u≈0 at all contact edges). The per-edge SG decomposition is
numerically healthy; the IV mismatch originates upstream in the recombination model,
not in the contact-edge transport formulas.

### Ideality Factor Analysis

Vela ideality factor n monotonically decreases from ~1.66 at 0.15 V to ~1.24 at
0.30 V. Sentaurus maintains n ≈ 1.00 throughout 0.05–0.70 V. This is the classic
signature of a depletion-region SRH current (n≈2 component) competing with the
diffusion current (n≈1 component). The crossover occurs at ~0.17 V:

- Below 0.17 V: Vela > Sentaurus (SRH-dominated, excess n≈2 component)
- Above 0.17 V: Vela < Sentaurus (diffusion-dominated, but n_eff > 1 because
  residual SRH still inflates the effective denominator)

At 0.17 V, the depletion SRH current equals the diffusion current. Analytical
estimate with tau=1e-7, W≈0.12 μm, ni=9.65e9 cm-3:

- J_depl ≈ q·ni·W/(2·tau) · exp(V/2Vt) ≈ 1.6e-8 A/m at V=0.30 V
- J_diff ≈ q·Dp·ni²/Nd/L_n · exp(V/Vt) ≈ 1.0e-8 A/m at V=0.30 V

Both are comparable, confirming a mixed regime at 0.30 V with effective n≈1.24.

### Iteration B Tau Sweep Results

| Variant                | n(0.30V)  | orders@0.27V | orders@0.30V | I(0.29)/I(0.30) Δ |
|------------------------|-----------|--------------|--------------|-------------------|
| Baseline (tau=1e-7)    | 1.235     | 0.489        | 0.624        | **0.099** > target |
| B1: tau=1e-6 (10×)     | 1.035     | 0.787        | 0.836        | **0.056** ✓ <0.06  |
| B2: tau=1e-5 (100×)    | 1.005     | 0.832        | 0.864        | **0.048** ✓ <0.06  |
| B3: no SRH/Auger       | 1.002     | 0.837        | 0.868        | **0.047** ✓ <0.06  |
| Sentaurus (reference)  | ≈1.004    | —            | —            | 0.632              |

### Conclusions

1. **SRH lifetime is the primary cause of the ideality factor divergence.**
   Setting `taun=taup=1e-6` (10× the default 1e-7) brings n to 1.035 and
   slope delta to 0.056, meeting the `<0.06` target.

2. **Slope delta target is met with tau=1e-6.** The BV simulation already
   uses `taun=taup=1e-6` (verified in `test_sentaurus_sample_integration.py`
   line 112). Updating the IV simulation to the same value is consistent and
   safe with respect to BV non-regression.

3. **A residual ~0.84 orders absolute gap in I₀ remains even with no SRH.**
   The diffusion saturation current I₀ in Vela is ~7× lower than Sentaurus
   across the full forward-bias range (consistent factor from 0.17 to 0.30 V).
   This is not caused by tau and requires further investigation in Iteration C
   (candidates: BGN OldSlotboom vs Slotboom formula, effective ni value,
   contact BC ohmic minority-carrier equilibrium level, or diffusivity model).

4. **Contact-edge SG decomposition is correct.** The bernoulli_u≈0 at all
   contact edges confirms numerical health. The minority-carrier diffusion at
   cathode (57% of total at 0.30 V) is physically expected for the short-diode
   geometry (L_n=1 μm << Lp=11 μm at tau=1e-7).

5. **All 275 tests pass after the DCSweep.cpp brace fix.**

### Next Step (Iteration C)

Investigate the ~7× I₀ discrepancy between Vela and Sentaurus when SRH is
disabled. Candidates:
- Sentaurus OldSlotboom vs Vela "slotboom" BGN formula differences — confirm
  by running with `bandgap_narrowing` disabled in B3-style config and checking
  if Vela ni matches Sentaurus ni (expected value at 300 K for silicon).
- Effective intrinsic carrier density ni: Sentaurus may use ni=1.45e10 vs
  Vela's ni=9.65e9 cm-3 (factor 1.5², = 2.25× in I₀ — partial explanation).
- Contact ohmic BC: verify that minority-carrier equilibrium density at ohmic
  nodes is consistent between Vela and Sentaurus, especially with BGN applied.
- Diffusivity model: confirm Einstein relation Dp=μp·kT/q is used correctly.

---

## Iteration C Completion (2026-05-30)

A 2×2 factorial experiment (ni × BGN, all no-SRH) cleanly decomposes the I₀
diffusion-current gap. **The gap is fully explained by physical parameter
choices, not a solver defect.**

### Experiment Design

All four variants disable SRH/Auger to expose the pure diffusion saturation
current I₀. Generated by `build/pn2d_probe/vela/make_iter_c_configs.py`:

| Variant | ni (cm⁻³) | BGN      |
|---------|-----------|----------|
| c0      | 1.0e10    | slotboom |
| c1      | 1.45e10   | slotboom |
| c2      | 1.0e10    | none     |
| c3      | 1.45e10   | none     |

The `ni=1.45e10` variants use a custom `materials_file` (`si_ni145.json`).

**Important unit-scaling note (cost me one bad run):** under
`scaling.mode = unit_scaling`, the `materials_file` interprets `mun`/`mup` in
**cm²/V·s** (`mobilityToSI` ×1e-4) and `ni`/`Nc`/`Nv` in **cm⁻³**
(`concentrationToSI` ×1e6). The first attempt wrote `mun: 0.135` expecting SI
(m²/V·s); it was read as 0.135 cm²/V·s = 1.35e-5 m²/V·s, collapsing mobility
by 1e4 and dropping the current ~3000×. The corrected file uses
`mun: 1350, mup: 480` cm²/V·s (= 0.135 / 0.048 m²/V·s) and `Nc: 2.8e19,
Nv: 1.04e19` cm⁻³.

### Factorial Results (I₀ at 0.30 V, A/µm)

| Variant                  | I(0.30V)   | gap vs Sentaurus |
|--------------------------|------------|------------------|
| Sentaurus (reference)    | 7.229e−14  | —                |
| c0 (ni=1e10, slotboom)   | 9.804e−15  | 7.37×            |
| c1 (ni=1.45e10, slotboom)| 2.059e−14  | **3.51×**        |
| c2 (ni=1e10, no BGN)     | 7.668e−15  | 9.43×            |
| c3 (ni=1.45e10, no BGN)  | 1.610e−14  | 4.49×            |

### Clean Factor Separation

- **ni effect = 2.100×** (both c1/c0 and c3/c2), exactly matching the
  analytic prediction (1.45/1.0)² = 2.10 for I₀ ∝ ni². This confirms Vela's
  diffusion-current machinery scales correctly with ni.
- **BGN effect = 1.279×** (both c0/c2 and c1/c3): Vela's `slotboom` model
  raises I₀ by 1.28× at 1e17 doping by increasing ni_eff.
- The two effects are orthogonal and multiplicative (2.10 × 1.28 = 2.69),
  reducing the gap from 7.37× (c2 baseline-ni no-BGN would be 9.43×) to 3.51×.

### Conclusions

1. **ni is the dominant contributor.** Adopting the textbook/Sentaurus
   ni=1.45e10 cm⁻³ shrinks the I₀ gap from 7.37× to **3.51× (0.545 orders)**,
   comfortably within the `<1.0 orders` IV gate.

2. **Vela's SG diffusion-current formulation is correct.** The exact 2.10×
   ni-squared scaling rules out a bug in the discretization, Einstein relation,
   or ohmic minority-carrier BC.

3. **Residual 3.51× (0.545 orders) is the BGN-model gap.** Vela's `slotboom`
   gives only 1.28× enhancement, while Sentaurus `OldSlotboom` at 1e17 doping
   plus Fermi-Dirac statistics produces a larger effective ni_eff. This residual
   is within the quantitative gate and does not require a solver change.

4. **No source change made or required.** The current default
   `MaterialDatabase` ni=1.0e16 m⁻³ (1.0e10 cm⁻³) keeps the IV gate passing at
   0.87 orders. Changing the global default would perturb all 275 tests'
   reference values and is out of scope. The recommended path to match
   Sentaurus is a deck-level `materials_file` with ni=1.45e10 cm⁻³, applied
   only to the pn2d reference case.

### Combined Iteration A→C Outcome

| Lever                   | Effect on pn2d IV parity        |
|-------------------------|---------------------------------|
| DCSweep.cpp brace fix   | Enables contact-edge diagnostics |
| tau 1e-7 → 1e-6         | Ideality n 1.24 → 1.04; slope Δ 0.099 → 0.056 (meets <0.06) |
| ni 1.0e10 → 1.45e10     | I₀ gap 7.37× → 3.51× (0.87 → 0.55 orders) |
| Slotboom BGN (existing) | +1.28× I₀ (already enabled)     |

All 275 tests remain green; no solver/source files were modified in
Iteration C (only new analysis scripts and generated `build/` artifacts).

## Stage 0 Completion (2026-05-30): Promote-Lever Gate Quantification

Stage 0 quantifies how each candidate promote lever moves the **official IV
gate metric** (the `iv_report.orders_of_magnitude` value the regression test
asserts `<= 0.50`) before any deck change is made. Four IV variants were run on
the faithful pn2d deck, all keeping official physics
(`mobility=caughey_thomas_field`, `bandgap_narrowing=slotboom`,
`recombination=[srh,auger]`) and differing only in the two levers:

| Variant | tau (s) | ni (cm⁻³) |
|---------|---------|-----------|
| s0 baseline | 1e-7 (default) | 1.0e10 |
| s1 tau1e6 | 1e-6 | 1.0e10 |
| s2 ni145 | 1e-7 | 1.45e10 |
| s3 tau1e6_ni145 | 1e-6 | 1.45e10 |

Tooling: `build/pn2d_probe/vela/make_stage0_configs.py` (fine step 0.01),
`make_stage0_gate_configs.py` (deck step 0.1), and `scripts/analyze_stage0.py`
(imports `compare_reference_curves` so the gate is computed identically to the
regression path: `max |log10(|cand·(-1)| / |ref|)|` aligned at the reference's
own bias points in `[0.2, 0.3]`).

### Fine-Step Gate (vela_step=0.01 — physically faithful)

| Variant | GATE orders | I(0.30) A/µm | gap× | slope `I(0.29)/I(0.30)` | slope Δ vs Sentaurus | handoff |
|---------|------------:|-------------:|-----:|------------------------:|---------------------:|---------|
| s0 baseline | 0.5606 | 1.720e−14 | 4.20 | 0.7311 | 0.0987 | newton ✓ |
| s1 tau1e6 | 0.7995 | 1.054e−14 | 6.86 | 0.6880 | 0.0556 | newton ✓ |
| s2 ni145 | **0.2968** | 3.198e−14 | 2.26 | 0.7219 | 0.0894 | newton ✓ |
| s3 tau1e6_ni145 | 0.4868 | 2.173e−14 | 3.33 | 0.6859 | **0.0534** | newton ✓ |

Sentaurus reference: `I(0.30)=7.229e−14 A/µm`, `slope=0.6324`.

### Authoritative Gate At Deck Step (vela_step=0.1)

| Variant | GATE orders | decision |
|---------|------------:|----------|
| s0g baseline | 0.4212 | PASS |
| s1g tau1e6 | 0.3085 | PASS |
| s2g ni145 | 0.6668 | FAIL (>0.50) |
| s3g tau1e6_ni145 | 0.5715 | FAIL (>0.50) |

### Critical Finding: The Step-0.1 Gate Is An Interpolation Artifact

The two step settings give **opposite rankings**. The cause is mechanical, not
physical. At `vela_step=0.1` the candidate curve has points only at
`0.0/0.1/0.2/0.3`, so the gate (which aligns at reference biases
`0.2047…0.29`) interpolates the candidate **linearly between 0.2 and 0.3**.
The true IV is a steep convex exponential, so the straight chord massively
**overshoots** at the low edge of the window. Per-bias breakdown of the coarse
baseline:

| ref bias | ref I (A/µm) | coarse cand I (A/µm) | order |
|---------:|-------------:|---------------------:|------:|
| 0.2047 | 1.72e−15 | 2.98e−15 | 0.238 |
| 0.2247 | 3.70e−15 | 9.77e−15 | **0.421** ← max |
| 0.2447 | 7.99e−15 | 1.66e−14 | 0.317 |
| 0.2900 | 4.57e−14 | 3.19e−14 | 0.156 |

The gate maximum (`0.4212`) sits at `0.2247 V`, where the linear chord sits
`2.6×` **above** both the true Vela curve and the reference. Raising the current
with `ni=1.45e10` makes this low-edge overshoot worse (`0.6668`), so the coarse
gate *penalizes* a physically-correct current increase. The fine-step gate
removes the artifact and shows the opposite, correct behavior: `ni=1.45e10`
genuinely halves the orders gap (`0.5606 → 0.2968`).

### Stage 0 Decisions

1. **The official step-0.1 IV gate is not a trustworthy promotion metric.** It
   is dominated by linear-interpolation overshoot of the steep forward-bias
   exponential at the low edge of `[0.2, 0.3]`. Any promote evaluation must use
   a fine IV step (`<= 0.02`). Recommended deck change before promotion:
   lower `vela_step` for the pn2d IV simulation (e.g. `0.02`) so the gate
   reflects physics, not chord overshoot.

2. **`ni=1.45e10 cm⁻³` is the only single lever that reduces the true orders
   gap** (fine-step `0.5606 → 0.2968`, gap `4.20× → 2.26×`). It barely changes
   slope (Δ `0.099 → 0.089`).

3. **`tau=1e-6` improves slope but worsens the true orders gap** (fine-step
   `0.5606 → 0.7995`, gap `4.20× → 6.86×`). Increasing the SRH lifetime lowers
   the depletion-region recombination current and so lowers I, moving the
   already-low Vela curve further from the reference magnitude. This confirms
   the Iteration B tension: the slope and magnitude levers pull in opposite
   directions.

4. **`tau=1e-6 + ni=1.45e10` is the best combined candidate**: fine-step gate
   `0.4868` (passes `<0.50`), best slope Δ `0.0534`, strict Newton handoff
   preserved. It satisfies both the slope and the (fine-step) magnitude gate,
   but with thin margin against an artifact-sensitive threshold.

5. **All four variants converge with strict Newton handoff** (`handoff=newton`,
   `newton_iterations>0`) and matching trend, so none introduces a solver
   regression.

### Implication For The Plan

- **Stage 2 promotion is blocked on a deck-step / gate-methodology fix first.**
  Promoting `tau+ni` against the current step-0.1 gate would *fail* the test
  (`0.5715 > 0.50`) for a non-physical reason. The correct sequence is: (a)
  move the pn2d IV deck to a finer step, (b) re-baseline the documented gate,
  (c) then promote `tau+ni` with the fine-step gate as the guard.
- The residual `2.26×` (`0.30` orders) at `ni=1.45e10` remains the BGN-model
  gap identified in Iteration C and is within the `<1.0` order envelope.

No source or deck files were modified in Stage 0 (only new analysis scripts and
generated `build/pn2d_probe/` artifacts). All 275 tests remain green.

## Stage 2 Promotion (2026-05-31): Fine-Step Deck + `ni=1.45e10` Bundle (Path B)

Following the Stage 0 finding that the step-`0.1` IV gate is an
interpolation artifact, Stage 2 promotes the physically-correct `ni=1.45e10`
magnitude lever together with a fine deck step so the regression gate measures
physics rather than chord overshoot. This is **Path B**: bundle the fine-step
deck change with the `ni` promote and tighten the guard.

### Changes Landed

1. **New material override** `reference_tcad/pn2d/pn2d_iv_materials.json` — a
   minimal `materials` array overriding only Si `ni = 1.45e10 cm⁻³`
   (`concentrationToSI` ⇒ `1.45e16 m⁻³`). All other Si fields (mun/mup, Nc/Nv,
   eps_r, bandgap) inherit the built-in defaults, so the override is fully
   isolated to `ni` and scoped to the pn2d IV deck only. The global
   `MaterialDatabase` default (`ni = 1.0e16 m⁻³`) is untouched, so the other
   274 tests are unaffected.

2. **Importer channel** `scripts/sentaurus_import.py` — added a
   `vela_materials_file` simulation field. When present, the import flow copies
   the named file from `--source-dir` into the generated `vela/` directory and
   sets the deck's top-level `materials_file` to its basename (resolved at
   runtime against the deck directory via `cwd=deck_path.parent`). The channel
   is per-simulation, so it applies to IV without touching BV.

3. **Deck** `reference_tcad/pn2d/pn2d_reference.json` IV block —
   `vela_step: 0.1 → 0.02`, added `vela_materials_file: "pn2d_iv_materials.json"`,
   and tightened `comparison.max_orders_of_magnitude: 1.0 → 0.4`. BV block and
   global defaults unchanged.

4. **Regression test** `tests/regression/test_sentaurus_sample_integration.py` —
   tightened the IV gate assertion `<= 0.50 → <= 0.40`, added assertions that
   the IV deck step is `0.02`, that `materials_file == "pn2d_iv_materials.json"`
   and the copied file exists, and that the BV deck has **no** `materials_file`
   (isolation guard). BV assertions are unchanged.

### Validated Outcome (regenerated reference tree)

| metric | before (step 0.1, ni 1e10) | after (step 0.02, ni 1.45e10) |
|---|---|---|
| IV `orders_of_magnitude` | `0.4212` (artifact) | **`0.2731`** |
| IV gate guard | `<= 0.50` | `<= 0.40` (now honest) |
| IV trend_match | yes | yes |
| IV deck step | `0.1` | `0.02` |
| BV `orders_of_magnitude` | `0.1109` | `0.1109` (unchanged) |
| BV `materials_file` | none | none (isolated) |

The honest fine-step gate dropped from the true baseline `0.5606`
(step 0.01, ni 1e10) to `0.2731` at the deck step `0.02` with `ni=1.45e10` —
within the Iteration C BGN-residual envelope (`~0.30` orders, `2.26×`). All
IV/BV rows keep strict Newton handoff (`handoff_stage=newton`,
`newton_iterations>0`).

### Step Selection

Step `0.02` was chosen from the Stage 2 step-sensitivity sweep
(`scripts/analyze_stage2_step.py`): it is the coarsest deck step whose baseline
gate (`0.5386`) is within `0.022` of the true fine value (`0.5606`) while
`ni=1.45e10` cleanly ranks below baseline, at roughly half the solve cost of
step `0.01`.

### Test Status

`ctest -R sentaurus` → 3/3 passed (`sentaurus_sample_integration` 103 s).
Full suite `ctest --preset windows-ucrt64-debug` → **275/275 passed**, confirming
no regression across the other 272 tests (the `ni` override is isolated to the
pn2d IV deck).

## PN2D Sentaurus2018 0V Current Debug (2026-06-12)

The 0V current-related debug pass generated a five-group comparison report under:

```text
build/reference_tcad/pn2d_sentaurus2018/reports/0v_current_related
```

### Current Definition Findings

The baseline terminal signs match Sentaurus, but the absolute current does not:

| contact | Vela total (A/um) | Sentaurus `.plt` total | Sentaurus/Vela |
|---|---:|---:|---:|
| Anode | `-6.5533928359887347e-18` | `-7.17389811693691e-25` | `1.0946845849894114e-07` |
| Cathode | `6.5556001087542772e-18` | `7.17389811693687e-25` | `1.0943160043207828e-07` |

Simple width conversions do not explain the baseline mismatch. Comparing
Sentaurus internal current definitions shows a separate reference-definition
ambiguity:

| contact | `.plt TotalCurrent` | `ContactCurrentFlux` | boundary `TotalCurrentDensity` integral |
|---|---:|---:|---:|
| Anode | `-7.17389811693691e-25` | `-1.45982e-19` | `1.4598200000000005e-15 A/cm-width` |
| Cathode | `7.17389811693687e-25` | `6.65229e-20` | `6.652300000000003e-16 A/cm-width` |

`ContactCurrentFlux` is consistent with the boundary current-density integral
after a `1e4` width conversion, while `.plt TotalCurrent` remains more than
three orders smaller. This means no Vela current-scaling fix should be made from
`.plt` parity alone.

### Contact Driver Findings

Contact quasi-Fermi values are not the source of the baseline current mismatch:

| contact | max eQF contact diff | max hQF contact diff |
|---|---:|---:|
| Anode | `1.83689e-16 V` | `4.59224e-17 V` |
| Cathode | `2.29612e-16 V` | `0.0 V` |

The contact majority densities match Sentaurus, while minority densities differ:

| contact | mean e-density ratio Vela/Sentaurus | mean h-density ratio Vela/Sentaurus |
|---|---:|---:|
| Anode | `0.46664818246822415` | `0.9999999999999999` |
| Cathode | `0.9999999999999999` | `0.46664818246822415` |

Contact-to-interior QF deltas in the Vela edge-current path are tiny
(`~1e-11 V` max), so the contact Dirichlet QF is not the immediate source of the
large baseline edge current.

### BGN Isolation Finding

The `no_bgn` variant is decisive:

| variant | QF max span | Anode Vela current | Cathode Vela current | Vela A/m vs Sentaurus |
|---|---:|---:|---:|---:|
| baseline | `0.0043930610509000005 V` | `-6.553392835988735e-18 A/um` | `6.555600108754277e-18 A/um` | `~9.14e12` |
| no_bgn | `8.39006e-13 V` | `-1.0647573418528946e-30 A/um` | `-8.8957506957498075e-31 A/um` | `~1.2-1.5` |

Turning off BGN collapses both the body QF split and the Vela current in the
physical A/m convention to Sentaurus `.plt` scale. The selected root-cause
hypothesis is therefore BGN/effective-ni consistency across the Newton state,
continuity residual, density reconstruction, and contact-current post-process.

### Debug Decision

Do not apply a current-unit fix yet. The Sentaurus `.plt`, TDR
`ContactCurrentFlux`, and boundary current-density integral are not mutually
consistent after simple conversions, while `no_bgn` strongly implicates the BGN
state path. The next implementation phase should add a focused failing
equilibrium test for BGN-enabled 0V QF flatness, then inspect only the
BGN/effective-ni path in `NewtonSolver`, `CoupledDDAssembler`, and
`ContactCurrent`.

### Follow-up Fix: Variable-ni Quasi-Fermi SG Flux

The follow-up debug implemented a variable-intrinsic-density quasi-Fermi
Scharfetter-Gummel flux. The old coupled Newton path used the balanced
quasi-Fermi flux only when `ni_i == ni_j`; BGN makes `ni_eff` node dependent, so
the code fell back to density SG on BGN edges. That density fallback does not
cancel flat quasi-Fermi levels when `ni_eff` varies, producing the 0V QF split
and artificial terminal current.

The fix adds variable-ni QF flux functions and uses them in both the coupled
continuity residual/Jacobian and `ContactCurrent` post-processing. New unit
coverage verifies:

- variable-ni SG electron/hole flux is zero for flat e/h quasi-Fermi levels;
- `CoupledDDAssembler` BGN continuity residuals vanish for flat e/h
  quasi-Fermi levels;
- BGN coupled residuals intentionally diverge from the older density-form
  `DDAssembler` reference on nonuniform-`ni` edges.

Refreshed 0V diagnostics after the fix:

| metric | before fix | after fix |
|---|---:|---:|
| current-balance status | `diagnostic_fail` | `pass` |
| classification | `contact_boundary_qf_state` | `balanced` |
| baseline QF max span | `0.0043930610509000005 V` | `1.827779e-08 V` |
| Anode total current | `-6.553392835988735e-18 A/um` | `-1.656945621904068e-27 A/um` |
| Cathode total current | `6.555600108754277e-18 A/um` | `3.3547470416890355e-28 A/um` |
| terminal abs pair sum | `2.2072727655425378e-21 A/um` | `1.3214709177351644e-27 A/um` |
| terminal balance gate | relative pass | absolute floor pass |

The current-related baseline report now shows edge currents near numerical
zero: Anode mean absolute edge current `9.746738952376869e-23 A/m`, Cathode
mean absolute edge current `1.9733806127582562e-23 A/m`. The remaining
Sentaurus `.plt` absolute-current mismatch should not be used as a Vela unit
fix trigger because the Sentaurus `.plt`, TDR `ContactCurrentFlux`, and boundary
current-density integral definitions remain mutually inconsistent.

## PN2D BV -20 V Blocked Status

Validation date: 2026-06-18.

The Sentaurus-default BV parity target is now explicit in the generated Vela
deck:

```json
"impact_ionization": {
  "model": "van_overstraeten",
  "driving_force": "quasi_fermi_gradient",
  "generation": "current_density",
  "current_approximation": "density_gradient"
}
```

The committed reference gate remains low-bias only with `vela_stop = -0.05` and
`vela_step = -0.05`. Do not promote the pn2d BV gate to `-20 V` yet.

Fresh execution artifacts are under
`build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution`.
The no-impact branch converged restart points at `-8 V`, `-10 V`, and `-13.2 V`.
The Sentaurus-default SG edge-current branch converged single-point solves at
`-13.2 V` and `-20 V` from the high-bias restart states and wrote C++ SG
edge-source dumps.

The high-bias gate fails. In the focus-edge region:

| branch | bias | median log10 Vela/Sentaurus electron density | focus-edge log10 Vela/Sentaurus generation |
|---|---:|---:|---:|
| no-impact | `-10 V` | `-0.297` | about `-0.303` |
| no-impact | `-13.2 V` | `+2.654` | about `+2.511` on active generation edges |
| SG avalanche | `-13.2 V` | `+2.654` | `+2.511` |
| SG avalanche | `-20 V` | `+1.787` | `+1.382` |

The `-13.2 V` mismatch is already present without impact ionization, so
avalanche source tuning is not the first fix. The selected blocker is
high-bias quasi-Fermi/state anchoring between about `-10 V` and `-13.2 V`.

SG edge-source implementation consistency is good at `-13.2 V`:
`log10(total C++ source / Python source) = +0.00804`, with both paths assigning
about `96.6%` of source to interior-bulk edges. At `-20 V`, the total source is
still within about `0.347` decades, but the C++ diagnostic assigns `54.2%` of
source to contact edges while the Python reconstruction classifies almost all
source as interior bulk. Resolve this contact-edge/source-reporting discrepancy
before any `-20 V` promotion.

Ionization-integral breakdown analysis is a future post-processing diagnostic
only. It must not replace the self-consistent Sentaurus-default drift-diffusion
avalanche comparison for this validation gate.

A lightweight edge-local ionization-integral proxy was added as
`scripts/diagnose_pn2d_bv_ionization_integral.py`. It is not a field-line
integrator; it reports dominant mesh-edge `alpha * dx` values after converting
Vela VTK high-field scalars from `V/cm` scale to `V/m`. On the current
Sentaurus-default SG states:

| bias | max edge-local integral | electron max | hole max |
|---:|---:|---:|---:|
| `-13.2 V` | `0.08193` | `0.08193` | `0.02804` |
| `-20 V` | `1.80364` | `1.80364` | `1.40371` |

This explains why `-20 V` is deep in the avalanche-prone regime, but it does
not pass the Sentaurus-default BV gate because the self-consistent local source
and carrier-density parity checks still fail.

### No-Impact Branch Follow-Up

Additional no-impact probes on 2026-06-18 narrowed the remaining `-13.2 V`
blocker.

The previous explicit `-10 V -> -13.2 V` no-impact restart was re-run as a
small-step continuation from a `-10 V` VTK restart with `step = -0.05 V` and no
explicit `bias_points`. The run reached `-13.2 V` in 65 points, but the
focus-edge median electron-density error remained high:

| probe | `-13.2 V` median log10 Vela/Sentaurus electron density |
|---|---:|
| explicit `-10 -> -13.2 V` restart | `+2.654` |
| small-step restart from `-10 V` | `+2.597` |

The no-impact high-current branch begins near `-12.75 V`, where the terminal
current leaves the `~1e-11 A/m` leakage scale and enters the `~1e-9` to
`~1e-8 A/m` electron-current branch. A Gummel-initialized handoff from the
pre-jump `-12.7 V` state changed the failure mode, shrinking near
`-12.9466659 V`, but did not recover Sentaurus-like carrier densities.

A local material `ni` override to `1.6556153e10 cm^-3` confirms that material
intrinsic density is a real low-bias calibration axis, not the high-bias branch
fix:

| probe | `-10 V` median log10 e-density ratio | `-13.2 V` median log10 e-density ratio |
|---|---:|---:|
| default `ni = 1.0e10 cm^-3` | `-0.297` | `+2.597` |
| local `ni = 1.6556153e10 cm^-3` | `-0.072` | `+2.995` |

Contact quasi-Fermi Dirichlet values were checked directly and match the
electrode biases at both `-10 V` and `-13.2 V`; the contact potential offset is
the known material-`ni` built-in difference of about `13 mV`. The high-bias
error is internal: near the focus edge at `-13.2 V`, Vela `psi` is about
`0.14 V` above Sentaurus while `phin` differs by only about `0.03 V`.

Updated blocker: the remaining first-order task is the no-impact coupled
continuity/Poisson branch selection around `-12.7 V` to `-13.0 V`, especially
the electron-continuity residual/Jacobian balance that admits the high-density
internal solution. Do not promote the `-20 V` BV gate, and do not treat material
`ni`, SRH lifetime, mobility, or avalanche tuning as the next acceptance fix.

Follow-up high-precision residual probes show that this high-density branch is
not a false convergence caused by VTK output truncation or continuity residual
scaling. The diagnostic first confirmed that VTK-roundtripped states can show a
spurious `~0.1` Poisson block residual; converting the final high-precision
`latest_state.csv` at `-13.2 V` instead gives `psi ~= 1.37e-8`,
`phin ~= 7.27e-12`, and `phip ~= 7.27e-12`.

The no-impact branch was then re-run from the same `-10 V` restart to selected
target biases, saving high-precision final states and probing them with the
current C++ Newton residual evaluator:

| bias (V) | electron current (A/m) | `psi` block | `phin` block |
|---:|---:|---:|---:|
| `-12.70` | `-4.475e-11` | `3.28e-09` | `6.33e-12` |
| `-12.75` | `-6.424e-10` | `1.17e-09` | `1.01e-11` |
| `-12.80` | `-3.124e-09` | `1.01e-10` | `1.14e-11` |
| `-12.85` | `-9.001e-09` | `2.24e-08` | `1.99e-11` |
| `-12.90` | `-1.389e-08` | `1.17e-11` | `6.72e-12` |
| `-13.20` | `-1.375e-08` | `1.37e-08` | `7.27e-12` |

At focus nodes `351/986`, `log10(electrons_m^-3)` jumps from about `9.88` at
`-12.70 V` to about `11.14` at `-12.75 V`, then saturates near `12.29` after
`-12.90 V`. Since the accepted high-density states are internally
residual-balanced, the next implementation work should instrument
accepted-step Newton history and coupled Jacobian/continuation behavior around
`-12.70 V -> -12.75 V`, rather than loosening/tightening residual thresholds or
tuning avalanche, SRH, mobility, contact relaxation, or material `ni`.

Accepted-step Newton history diagnostics were then added and run on the same
branch jump. For the direct `-12.70 V -> -12.75 V` no-impact step, Newton took
five full accepted iterations at `-12.75 V`; line-search damping stayed at
`1.0` and the dominant residual was the Poisson block. The first accepted
iteration at `-12.75 V` had `psi block ~= 3.81`, while `phin block ~=
1.95e-10`; the final accepted iteration reached `psi block ~= 1.64e-8` and
`phin block ~= 1.38e-13`.

A smaller `-0.005 V` continuation step from the same `-12.70 V` restart delayed
but did not remove the transition:

| bias | focus `log10(electrons_m^-3)` | note |
|---:|---:|---|
| `-12.75 V` | `~10.71` | lower than the direct `-0.05 V` jump |
| `-12.77 V` | `>= 11` | first threshold crossing |
| `-12.84 V` | `>= 12` | high-density branch established |
| `-13.20 V` | `~12.34` | still high-density at target bias |

Across the `-0.005 V` run, line-search damping remained `1.0`, maximum Newton
iterations per point were `3`, and final residuals were again dominated by the
Poisson block. This points away from line-search failure or electron-continuity
residual imbalance and toward electrostatic Newton branch control:
pseudo-transient continuation, trust-region/max-update policy, or
Sentaurus-like coupled-variable extrapolation should be the next opt-in solver
experiments.

The existing `solver.max_update` cap was then tested as a trust-region proxy.
All variants used the no-impact `-12.70 V` high-precision restart,
`step = -0.005 V`, and ran to `-13.20 V`:

| `max_update` | final focus `log10(electrons_m^-3)` | final electron current (A/m) | max Newton iters/point |
|---:|---:|---:|---:|
| `5.0` | `12.343` | `-1.568e-08` | 3 |
| `1.0` | `12.343` | `-1.568e-08` | 4 |
| `0.5` | `12.343` | `-1.568e-08` | 6 |
| `0.2` | `12.343` | `-1.568e-08` | 8 |
| `0.1` | `12.343` | `-1.568e-08` | 14 |
| `0.05` | `12.343` | `-1.568e-08` | 25 |

Tighter values, `0.02` and `0.01`, fail at the `-12.70 V` restart point with
`max_iterations`; the remaining residual is still Poisson-dominated. Therefore
plain Newton update capping is not enough to recover the Sentaurus-like
low-density branch. The next useful solver work should be a true
continuation-control experiment: pseudo-transient/homotopy, coupled-variable
predictor/extrapolation control, or a block-aware trust region.

An external linear-predictor proxy was tested before adding any production
predictor code. It uses two prior high-precision restart states and solves each
target bias from:

```text
x_pred = x_curr + alpha * (x_curr - x_prev)
```

This proxy is a real branch-control lever. With `-0.005 V` target spacing, the
focus electron density evolves as:

| bias | focus `log10(electrons_m^-3)` | terminal current state |
|---:|---:|---|
| `-12.75 V` | `~10.73` | finite |
| `-12.80 V` | `~11.61` | finite |
| `-12.85 V` | `~12.10` | finite |
| `-12.90 V` | `~12.29` | finite |
| `-12.95 V` | `~11.96` | terminal current columns collapse to zero |
| `-13.20 V` | `~10.88` | terminal current columns remain zero |

So predictor/extrapolation can move the solution away from the previously
stable high-density branch, but the current external proxy lands on a suspicious
low-current branch and is not an acceptance fix. A production predictor should
be opt-in and TDD-covered, with explicit diagnostics for predicted initial
state quality, terminal-current consistency, and no-impact branch parity.

## BV High-Bias Calibration Baseline (2026-06-21)

The materials-aligned reverse-bias deck
(`build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/materials_aligned_to_m13p2/`)
is frozen as the accepted pn2d BV calibration baseline. It is a strict
coupled-Newton (`gummel_newton`) `bv_reverse` sweep that enables the full
Sentaurus-faithful physics block (`van_overstraeten` avalanche with
`quasi_fermi_gradient` driving force, `masetti_field` mobility, `old_slotboom`
BGN, SRH recombination). Reruns are bit-reproducible across local rebuilds: the
`-13.2 V` terminal current reproduces to four digits (`-5.8554e-17` vs
`-5.8558e-17 A/um`).

### Field parity at -13.2 V

Same-bias field comparison
(`scripts/compare_pn2d_bv_multibias_fields.py`, report under
`reports/stepE_bv_rerun_compare/`) confirms the depletion-region state matches
Sentaurus closely on every term except the minority quasi-Fermi level:

- `Potential`: RMS `0.0034 V`, centerline `0.00097 V` (excellent, ~3 mV).
- `ElectricField`: junction-local error `~11%` (the `relative_p95 0.775` metric
  is inflated by near-zero low-field tail nodes and is not representative).
- `eDensity` / `hDensity`: log10 p95 `0.156` / `0.170` decades, the signature of
  the `~7 mV` depletion quasi-Fermi-level offset.
- `SRHRecombination`: log10 p95 `0.0046` decades (effectively identical).
- `AvalancheGeneration`: meaningful comparison is `p99` `1.22e15` vs `3.19e15`
  (the nodal `0.38x` tail-diluted ratio); the reported `12.8` decade log error is
  a near-zero-floor artifact, not a physical discrepancy.

### Full-range IV (0 to -20 V)

The `-20 V` extension deck
(`reports/stepB_bv_minus20/simulation_bv_minus20.json`, generated from the frozen
`-13.2 V` deck by changing only `sweep.stop` and the output paths) converges over
the entire ramp with no breakdown trigger by `-20 V`
(`max_electric_field ~5.6e7 V/m`). The log10-magnitude IV RMS versus the
Sentaurus reference curve (`scripts/compute_bv_iv_rms_minus20.py`) is:

| segment | log10 RMS (decades) | notes |
|---|---:|---|
| `0..-5 V` | `0.0300` | excellent (~7% current) |
| `-5..-10 V` | `0.0376` | excellent (~9%) |
| `-10..-13.2 V` | `0.1161` | calibration band |
| `-13.2..-20 V` | `0.3916` | diverges approaching breakdown |
| full `0..-20 V` | `0.2631` | |

Representative ratios `Vela/Sentaurus`: `-1 V` `1.09`, `-5 V` `1.04`,
`-10 V` `0.81`, `-13.2 V` `0.70`, `-15 V` `0.64`, `-18 V` `0.35`, `-20 V` `0.13`.
The agreement is excellent up to `-10 V` and widens monotonically toward
breakdown. This is the expected multiplication-integral sensitivity: with
`M = 1/(1 - I_ion)`, Sentaurus reaches `I_ion -> 1` (avalanche runaway) at a
lower reverse bias than Vela, so a sub-percent `I_ion` difference is amplified
without bound as `M -> infinity`. Sentaurus current rises ~2x/V at `-20 V`
(approaching its breakdown voltage just beyond the swept range) while Vela's
multiplication grows more gently, giving Vela a slightly higher effective
breakdown voltage.

### BV Acceptance Scope After `avaljac`

The `avaljac` branch demonstrates that the Sentaurus-faithful BV physics block can
be continued to `-20 V` without Newton failure. This is a convergence milestone,
not a final BV parity acceptance. The full-curve shape gate remains open because
the current Vela curve does not reproduce the Sentaurus one-volt growth knee in
the `-18 V..-20 V` region.

Accepted status is limited to:

- SG avalanche source Jacobian completeness for the current production path.
- SRH/Auger local derivative coverage in the coupled residual/Jacobian.
- End-to-end continuation robustness to `-20 V` for the current branch.

Open status remains:

- Full high-bias (`-18 V`, `-19 V`, `-20 V`) real-state replay after the default `avaljac` curve/config/state artifacts are regenerated.
- Curve-shape parity over `-10 V..-20 V`.
- The physical cause of the missing high-bias one-volt current-growth knee.

Reproduce:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
$base = "build-release\reference_tcad\pn2d_sentaurus2018"
build-release\vela_example_runner.exe --config "$base\reports\stepB_bv_minus20\simulation_bv_minus20.json"
python scripts\compute_bv_iv_rms_minus20.py --reference "$base\reference_curves\pn2d_sentaurus2018_bv_reference.csv" --candidate "$base\reports\stepB_bv_minus20\materials_aligned_minus20.csv"
```

### PN2D BV Jacobian Audit Baseline

Current BV production physics uses `van_overstraeten`, `driving_force:
quasi_fermi_gradient`, `generation: current_density`, and
`current_approximation: density_gradient`. In this SG edge-current avalanche
path, `CoupledDDAssembler::assembleJacobian` finite-differences the combined
edge avalanche source with respect to the six endpoint potentials, so the
matrix includes carrier-density, alpha driving-field, and local field-dependent
edge-mobility derivatives for that source discretization.

The non-SG node-local avalanche path remains intentionally approximate: it
includes local carrier-density derivatives but omits driving-field and mobility
derivatives. That path is not the current PN2D BV production path and must not
be used as evidence that the SG BV Jacobian is incomplete.

### PN2D BV Real-State Jacobian Block Replay (2026-06-22)

A real-state replay was run on the available Vela BV restart fields reconstructed
from `bv_newton_residual_states` at `-10 V` and `-13.2 V`. The initial replay
localized a non-avalanche mismatch to the ordinary transport block when
quasi-Fermi-gradient high-field mobility was active. After adding the missing
transport mobility potential sensitivity, the available full-state block audit
reports:

| bias | poisson | transport | srh_auger | sg_avalanche | dirichlet_or_gauge |
|---:|---:|---:|---:|---:|---:|
| `-10 V` | `9.039746e-06` | `1.427413e-05` | `7.623126e-15` | `9.379560e-19` | `2.173359e-10` |
| `-13.2 V` | `1.255905e-05` | `1.932885e-05` | `7.428910e-15` | `2.725849e-11` | `2.083001e-10` |

Artifact: `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_real_state_jacobian_audit/jacobian_blocks_real_state.csv`.

A follow-up high-bias replay converted the existing
`vela/sg_edge_current_vtk/dc_sweep_0006_-20V.vtk` state into restart CSV and ran
the same block probe at `-20 V` under the `sg_edge_current` BV configuration:

| bias | poisson | transport | srh_auger | sg_avalanche | dirichlet_or_gauge |
|---:|---:|---:|---:|---:|---:|
| `-20 V` | `2.487927e-05` | `3.122208e-05` | `7.349736e-15` | `2.544407e-05` | `1.944513e-10` |

Artifacts:

- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_real_state_jacobian_audit_sg_edge_current/states/bv_state_bias_m20p000000.csv`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_real_state_jacobian_audit_sg_edge_current/jacobian_blocks_real_state.csv`

Decision: for the available real BV states, including the VTK-derived `-20 V`
state, the SG avalanche, SRH/Auger, and transport Jacobian blocks are no longer
the leading BV discrepancy. The remaining BV work should stay on curve-shape,
state-path parity, and absolute quasi-Fermi level alignment.

Caveat: after the high-field transport Jacobian fix, rerunning the old
`simulation_bv_minus20_sg_edge_current_probe.json` continuation from `0 V` no
longer reproduces the checked-in `-20 V` CSV path; it stops near
`-2.4414e-5 V` with `max_iterations`. Therefore the checked-in historical curve
and VTK states remain useful for block replay, but current end-to-end
continuation robustness must be re-established before promoting any BV curve.

### PN2D BV Near-0V Early-Stop Globalization Check (2026-06-22)

The near-0V restart blocker was reproduced with a shortened deck derived from
`simulation_bv_minus20_sg_edge_current_probe.json` under:

- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_zero_early_stop_debug/`

With the current deck setting `solver.max_update = 5`, the sweep accepts only
`0 V` and `-2.44140625e-5 V`, then fails at
`-2.4421215057373049e-5 V` with `max_iterations`. Increasing the global cap to
`20` fails earlier, at `-1.220703125e-5 V`. The failed states have positive and
finite carriers, full damping, and residuals dominated by Poisson/hole blocks,
so this is a Newton globalization/stalling issue rather than an immediate
carrier-domain failure.

Disabling the global cap (`solver.max_update = 0`) removes the near-0V early
stop for the short `0 -> -0.05 V` deck, but a longer `0 -> -1 V` replay then
fails near `-0.091002197265625 V` with `nonfinite_residual`. A
`newton_step_probe` from the last stable state shows why: the uncapped Newton
step has non-finite norm and a local `phip` update of order `5.7e221 V` near
`x ~= 1.078 um, y ~= 0.188 um`; the trial hole density reaches
`~1.6e233 m^-3`, and the trial Poisson residual reaches `~1e210`. This localizes
the post-0V failure to an unconstrained quasi-Fermi update around the
junction/drift transition.

Physics isolation:

| variant | result |
|---|---|
| `impact_ionization.model = none` | converges to `-0.05 V` without Newton failure |
| constant mobility | converges to `-0.05 V` without Newton failure |
| `max_update = 0`, no QF-specific limit | passes near-0V, then fails near `-0.091 V` |
| `max_update = 0`, `quasi_fermi_update_limit_V = 0.0259` | converges to `-1 V`, 56 rows, max retry 1 |
| `max_update = 0`, `quasi_fermi_update_limit_V = 0.05` | converges to `-1 V`, 73 rows, max retry 1 |
| `max_update = 0`, `quasi_fermi_update_limit_V = 0.1` | converges to `-1 V`, 21 rows, max retry 0 |

Decision: the immediate near-0V early stop is caused by the global
`max_update` cap interacting with the corrected high-field transport Jacobian.
The safer continuation policy is not an unrestricted Newton step, but a
QF-specific update limit with the global cap disabled. The next BV gate should
run a full or staged `-20 V` sweep with `max_update = 0` and a
`quasi_fermi_update_limit_V` candidate, starting with `0.1 V` because it reached
`-1 V` with the fewest accepted rows and no retries.

### PN2D BV QF-Limit Staged Sweep Follow-Up (2026-06-22)

The recommended `quasi_fermi_update_limit_V` candidates were extended from
`-1 V` toward `-5 V`. This confirmed that the near-0V and `-0.091 V`
nonfinite-update blockers are cleared, but exposed a new continuation blocker
before `-5 V`:

| variant | last stable bias | failed bias | failure | residual at failure |
|---|---:|---:|---|---:|
| `qflim0p1_to5V` | `-2.92471424 V` | `-2.924716955 V` | `line_search_non_decrease` | `6.3757e-9` |
| `qflim0p05_to5V` | `-2.101646658 V` | `-2.101646777 V` | `line_search_non_decrease` | `1.2914e-9` |
| `qflim0p0259_to5V` | `-1.908734980 V` | `-1.908735992 V` | `line_search_non_decrease` | `1.1568e-9` |

All three failures have finite line-search trial residuals and positive finite
carriers; the residual is Poisson-dominated and sits just above the previous
hard-coded stall floor (`1e-9`) for the smaller QF limits. This is no longer the
same failure as the original near-0V `max_iterations` stall or the uncapped
`phip` overflow near `-0.091 V`.

A new solver knob, `solver.stall_residual_floor`, was added so this numerical
floor can be swept without changing code. The default remains `1e-9`, preserving
existing behavior. Diagnostic runs show the knob must be used narrowly:

| variant | result |
|---|---|
| `qflim0p1_abstol1e8_to5V` | advances to `-2.93479955 V`, then fails at residual `1.3444e-8` |
| `qflim0p1_abstol1e7_to5V` | accepts too coarse a path and fails earlier near `-2.85 V` with residual `1.91e-3` |
| `qflim0p1_stall2e8_to5V` | also accepts too coarse a path and fails near `-2.8498 V` with residual `4.95e-4` |
| `qflim0p05_stall1p5e9_to5V` | moves slightly past the strict `qflim0p05` failure, then stalls at `-2.103990259 V` with residual `2.36e-9` |

Decision: do not promote a loose residual floor as the BV fix. The useful part
of this pass is the solver instrumentation: `stall_residual_floor` is now
configurable for controlled diagnostics. The next implementation target should
be continuation/globalization quality around `-2 V..-3 V`, especially why the
raw Newton step remains large (`~38..91` in the smaller-QF-limit failures) when
the accepted residual is already near the Poisson numerical floor. Full `-20 V`
reruns should wait until this staged `-5 V` gate is stable without coarse-path
pollution.


### PN2D BV Newton Continuation Regression Localization (2026-06-23)

The regression boundary is `c04edbf -> 51d8bf1`. With the same imported PN2D BV
deck, `max_update=0`, and `quasi_fermi_update_limit_V=0.1`, `c04edbf` reaches
`-3 V` in 61 accepted points. Its Newton-history tail at `-3 V` has residual
`4.0711e-13` and raw step norm `4.5011`. The same deck on `51d8bf1` fails at
the next continuation point after `-0.3 V` with residual `~7e-9..8e-9` and raw
step norm `93.65..124.30`; the `c04edbf` to `51d8bf1` tail raw-step ratio is
`20.806`. The current branch (`85924e5` plus later changes) fails similarly
after `-0.40625 V`, with the comparison report showing residual `8.2502e-9`,
raw step norm `108.8725`, and tail raw-step ratio `24.188`.

The no-impact control localizes the blocker away from avalanche physics. With
`impact_ionization.model = "none"`, the canonical current-branch deck reaches
only `-0.8165 V` and then stops at `-0.8165000001609326 V` with
`max_iterations`, residual `8.5834e-9`, raw step norm `113.1072`, and positive
finite carriers. The impact-on canonical deck fails at `-0.40625000011175866 V`
with `line_search_non_decrease`, residual `1.9337e-8`, raw step norm
`107.5935`, and positive finite carriers. These controls reject the
`85924e5` contact driving-field fallback as the primary regression source.

A candidate clamp around Poisson recorrection after QF clipping was also
rejected: it passed the synthetic Newton-step test but did not move either
canonical `-3 V` gate. The solver patch was not kept. The next root-cause target
is therefore the raw linear solve and block coupling state before globalization
at the last stable points, especially why near-floor carrier or Poisson
residuals still produce a Newton step of order `100`.

### PN2D BV Carrier-Block Step Probe Follow-Up (2026-06-23)

The next raw-step probe used the last converged restart state and the failing
next contact bias, matching the continuation handoff state. This confirms that
the `O(100)` step is not driven by the Poisson solve or by Poisson recorrection.
For the impact-on gate, the full capped Newton step is `108.668`, while the
Poisson-only block step is only `1.7984e-8` and the carrier-only capped step is
`108.093`. For the no-impact gate, the full capped step is `113.108`, the
Poisson-only block step is `2.5885e-8`, and the carrier-only capped step is
`113.108`.

The carrier rows are already pathological before the QF cap. In the impact-on
probe, raw carrier deltas reach `7.47368e10 V` for `phin` and `666.932 V` for
`phip`; after capping, 545 electron-QF nodes and 172 hole-QF nodes sit exactly
at the `0.1 V` limit. In the no-impact probe, raw carrier deltas reach
`1.90984e217 V` for `phin` and `9.7053e208 V` for `phip`; after capping, 527
`phin` nodes and 328 `phip` nodes sit at the same limit. The weakest carrier-row
diagonal dominance is effectively zero: electron diagonal/row-sum is
`1.64e-13` in the impact-on probe and `6.66e-220` in the no-impact probe, while
the no-impact hole row also reaches `1.92e-211`.

The finite-difference Jacobian block audit does not point to a gross Jacobian
implementation error. Around the same states, analytic-vs-FD relative
mismatches are about `6.0e-7` to `9.3e-7` for `poisson`, `8.6e-6` to `1.1e-5`
for `transport`, and near numerical zero for `dirichlet_or_gauge`; no-impact
`sg_avalanche` is exactly zero as expected. The remaining target is therefore
carrier-block conditioning and row policy near depleted/floor carrier states,
not the avalanche source term and not a Poisson-block correction.

Artifacts were written under `build-release/bv_localization/canonical_probe/`,
including `probe_summary.json`, `jacobian_block_summary.json`, and the per-case
`probes/*.csv` files.

### PN2D BV Frozen High-Field Mobility Jacobian Gate (2026-06-23)

The carrier-block probe narrowed the continuation failure to the high-field
mobility Jacobian path rather than to avalanche, Poisson recorrection, or a
global carrier regularization. A low-field `masetti` control keeps the same
SRH/BGN/impact setup but removes high-field mobility limiting; both impact-on
and no-impact canonical `-3 V` gates then converge. A global
`carrier_regularization_scale` trial can reduce the one-shot carrier step, but
it destabilizes the full sweep at very low bias, so it is rejected as a
production fix.

Vela now exposes `solver.mobility.jacobian_field_derivatives`. The default is
`true`, preserving the analytic/finite-difference transport-Jacobian behavior
added for high-field mobility. The PN2D BV reference deck sets it to `false` for
the `masetti_field`/`quasi_fermi_gradient` mobility object, and also sets
`max_update=0` with `quasi_fermi_update_limit_V=0.1` for this BV solver path.
The residual still uses the Sentaurus-like high-field mobility while the Newton
matrix freezes mobility's field sensitivity. With those settings, the
frozen-Jacobian diagnostic passes the canonical gates:

| case | stop | points | result | last Newton iterations |
|---|---:|---:|---|---:|
| impact-on | `-3 V` | 61 | converged | 3 |
| no-impact | `-3 V` | 61 | converged | 3 |
| impact-on | `-5 V` | 101 | converged | 3 |
| impact-on | `-10 V` | 201 | converged | 3 |
| impact-on | `-20 V` | 401 | converged | 3 |

The `-20 V` terminal row from
`build-release/bv_localization/canonical_probe/frozen_mobility_jacobian/impact_on_to_20V/iv.csv`
has `current_total_A_per_um = -1.168307486e-16` and
`max_electric_field_V_per_cm = 5.607485472e5`. This is a Newton-continuation
milestone only; it does not by itself accept BV current magnitude or knee-shape
parity.


### PN2D BV Frozen-Jacobian Acceptance Refresh (2026-06-23)

Using the frozen high-field mobility Jacobian base config, the visual acceptance
run wrote VTK/state outputs at `0, -0.5, -2, -5, -10, -20 V` and converged all
6 requested points. The terminal `-20 V` row has
`current_total_A_per_um = -1.168116404e-16`,
`max_electric_field_V_per_cm = 5.607485472e5`, and 3 Newton iterations.
Artifacts are under
`build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_acceptance_visual/`.

Current-curve parity is mixed. The `-13.2..-13.0 V` current-window gate now
passes with ratios `0.8034, 0.8040, 0.8044`, and the sampled current errors are
small through `-10 V`: `0.0355` decades at `-2 V`, `0.0160` decades at `-5 V`,
and `0.0924` decades at `-10 V`. At `-20 V`, however, Vela is still low by
`0.8918` decades (`1.1681e-16 A/um` versus Sentaurus `9.1046e-16 A`).

The refreshed knee-shape gate remains diagnostic rather than accepted. In the
`-10..-20 V` window, Sentaurus first reaches a one-volt current-growth ratio
above `1.5` at `-19 V` and above `2.0` at `-20 V`; the frozen-Jacobian Vela
curve reaches neither threshold. The max absolute log-current error in this
window is `0.8917` decades. Field/state parity also remains open: at `-20 V`,
field compare reports electron-density and hole-density `log10_p95` errors of
`0.9207` and `0.9570`, electric-field `relative_p95` error of `0.8685`, and
thresholded avalanche-generation `log10_p95` error of `13.0551`.

This accepts the frozen-Jacobian path as a continuation fix, not as final BV
physics parity. The next physics target is the high-bias avalanche/current
feedback that creates the missing `-19..-20 V` knee.

### PN2D BV High-Bias Feedback Localization (2026-06-23)

The frozen high-field mobility Jacobian path now reaches `-20 V`, so the next
blocker is curve shape rather than continuation reachability. A field/source
summary under
`build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/`
shows that electric-field magnitude is already close at `-20 V` (`max` log10
Vela/Sentaurus ratio about `-8.6e-5`, p95 ratio about `0.005`, and sum-proxy
ratio about `0.011`). The same diagnostic shows avalanche/source generation is
low by about `1.16` decades at `-20 V`, so the missing knee is a feedback/current
problem rather than an alpha(E) or peak-field problem.

The continuity-feedback diagnostic was corrected to accept `--material-ni-m3`,
then rerun with the actual PN2D material value
`1.4638914958767616e16 m^-3`. With that correction, focus-edge quasi-Fermi
fields still match closely at `-20 V` (`5.5318e7` versus `5.5210e7 V/m` for
electrons, `5.6468e7` versus `5.6440e7 V/m` for holes), while edge fluxes and
source density remain low (`log10 Vela/Sentaurus generation = -0.8401` on focus
edge 2886). The focus endpoint effective intrinsic densities match Sentaurus
back-inferred values (`~1.9623e10 cm^-3` and `~1.6556e10 cm^-3`), and using the
Sentaurus `psi/phin/phip` state with Vela's corrected `ni_eff` reconstructs the
Sentaurus carrier densities. Therefore material `ni` and OldSlotboom BGN are not
the remaining high-bias discrepancy.

The remaining localized mismatch is the absolute state: at focus endpoints near
`-20 V`, Vela has `psi-phin` lower by about `47..48 mV` and `phip-psi` lower by
about `56 mV`. That produces electron-density errors of about `-0.79..-0.81`
decades and hole-density errors of about `-0.94..-0.95` decades, matching the
missing flux/source feedback. The next experiment should be a minimal
state-feedback or continuation/state-alignment probe around absolute
quasi-Fermi/carrier-density branch selection. Do not promote alpha, material-ni,
or hidden source-scale changes from this evidence.

### PN2D BV Absolute-State Feedback Probe (2026-06-23)

A minimal post-processing probe now quantifies the high-bias state hypothesis
without changing production physics. The script
`scripts/diagnose_pn2d_bv_absolute_state_feedback.py` reads the corrected
continuity-feedback edge/node CSVs and scales Vela edge source density by the
endpoint carrier-density factors reconstructed from Vela `ni_eff` with the
Sentaurus `psi/phin/phip` absolute state. This is a diagnostic source proxy, not
a replacement flux discretization.

Using
`build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/continuity_feedback_material_ni/`
as input, the probe wrote
`build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/absolute_state_feedback_probe/`.
For active edges at `-20 V`, the median original source gap is `-0.8401`
decades, the state-scaled gap is `+0.04295` decades, and the recovered gap is
`0.8827` decades. On focus edge 2886, the electron state factor is `6.286`, the
hole state factor is `8.848`, the flux-weighted state factor is `7.639`, and the
source proxy moves from `8.0748e21` to `6.1687e22 m^-3 s^-1` versus Sentaurus
`5.5879e22 m^-3 s^-1`. At `-10 V`, the corresponding active-edge recovered gap
is only `0.0610` decades, matching the already-small low-bias source mismatch.

This confirms the prior localization: the missing high-bias knee is explained
at the source-feedback level by the absolute quasi-Fermi/carrier-density branch
offset. The next production-facing investigation should inspect why the coupled
continuation settles on a lower-density absolute state at high reverse bias,
including state initialization, branch selection, and any gauge/contact policy
that can shift `psi-phin` and `phip-psi`. Do not promote this post-processing
source proxy, alpha(E) retuning, material `ni` changes, BGN changes, or hidden
source scaling as a fix.

### PN2D BV Absolute Branch Offset Probe (2026-06-23)

A full-node branch-offset probe now compares Vela and Sentaurus absolute states
by node class, contact membership, doping sign, and impact-active support. The
script `scripts/diagnose_pn2d_bv_absolute_branch_offsets.py` writes node rows and
group summaries under
`build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/absolute_branch_offsets/`.
It reports the absolute offsets `delta_psi`, `delta_phin`, `delta_phip` and the
gauge-invariant carrier exponents `delta(psi-phin)` and `delta(phip-psi)`.

The high-bias offset is not a contact Dirichlet or global electrostatic gauge
error. At `-20 V`, contact nodes remain aligned: contact median
`delta(psi-phin)` is about `-2.45e-5 V`, and the Anode/Cathode contact averages
are only about `+/-4.9e-5 V` in the two carrier exponents. In contrast,
non-contact impact-active nodes show a large internal quasi-Fermi branch offset:
median `delta_psi = 4.96e-6 V`, median `delta_phin = +0.04609 V`, median
`delta_phip = -0.05398 V`, median `delta(psi-phin) = -0.04733 V`, and median
`delta(phip-psi) = -0.05547 V`. The corresponding median density errors are
`-0.7949` decades for electrons and `-0.9318` decades for holes.

The onset is high-bias specific. For impact-active nodes, the median
`delta(psi-phin)` / `delta(phip-psi)` values are about `+0.00015/-0.00031 V` at
`-2 V`, `-0.00088/-0.00078 V` at `-5 V`, `-0.00318/-0.00408 V` at `-10 V`, and
`-0.04733/-0.05547 V` at `-20 V`. Splitting the `-20 V` impact-active support by
doping sign shows the same branch issue on both sides: p-type active nodes have
median density errors of `-0.574/-0.953` decades for electron/hole density, and
n-type active nodes have `-0.892/-0.859` decades.

This narrows the next production-facing experiment: contact boundary values and
electrostatic potential alignment are not the leading cause. The next probe
should target the interior high-field carrier-continuity branch, for example by
replaying the `-20 V` state with controlled quasi-Fermi shifts or adding a
continuation experiment that constrains the high-field active support toward the
Sentaurus absolute carrier-density branch while preserving contact Dirichlet
values. Any production change should be gated by moving the `-10..-20 V` knee
shape, not merely by improving a local source proxy.

### PN2D BV Active-Support QF Shift Replay (2026-06-23)

The controlled branch-alignment replay has now been run on the current frozen
high-field mobility-Jacobian `-20 V` state. The diagnostic uses the current
visual VTK state, regenerates the 99th-percentile avalanche support and SG edge
source decomposition, then replays the carrier-continuity residual proxy with
the median active-state offsets from the branch-offset probe:
`delta(psi-phin) = -0.04732618818171197 V` and
`delta(phip-psi) = -0.05547180240410299 V`. The residual-proxy script now
supports `--qf-shift-scope support_nodes` so the shift can be applied only to
thresholded active-support nodes rather than globally.

The support comparison itself is an important constraint. At the 99th
percentile, Sentaurus and Vela active avalanche nodes have zero overlap:
`20` Sentaurus-only false-negative nodes, `20` Vela-only false-positive nodes,
and peak separation `0.04752 um`. The Vela active-support integral is also much
smaller (`7.873498e16 cm^-3 s^-1`) than the Sentaurus active-support integral
(`1.1225339212493828e18 cm^-3 s^-1`). This means a local support replay is a
cause probe, not a production correction.

For Sentaurus-only active nodes, the baseline Vela-state transport over the
Sentaurus node source is only `0.0955` for electrons and `0.1185` for holes,
with residual/source medians `-0.927` and `-0.904`. A global QF branch shift
moves those ratios to `0.596` and `1.013`, nearly closing the hole-side residual
but leaving electron transport low. A support-only shift moves the electron
ratio closer to parity (`0.844`) but overshoots the hole ratio badly (`14.20`).
For Vela-only false-positive nodes, support-only replay overshoots electron
transport even more (`8.424` versus a Sentaurus-state reference of `0.891`).

This confirms the causal link between the internal high-field QF/carrier-density
branch and the missing avalanche feedback, but rejects hard active-node shifts
as a production path. The next useful experiment should use a smooth
active-region or continuation-level branch control that preserves contact
Dirichlet values and is judged by the curve-level `-10..-20 V` knee gate, not by
local residual-proxy improvement alone.
### PN2D BV Smooth Branch-Control Backscan (2026-06-23)

A continuation-level version of the active-support branch probe was then run via
`scripts/prepare_pn2d_bv_smooth_branch_state.py`. The helper reads the converged
frozen-Jacobian `-20 V` Vela VTK state, applies a Gaussian distance-weighted QF
shift around selected support nodes, preserves all contact nodes at zero weight,
reconstructs carrier densities with the original Vela inferred `ni_eff`, and
writes a DCSweep-compatible `initial_state_file`. This makes the experiment a
real Newton/continuation probe rather than another local source-only replay.

The first run targeted only Sentaurus false-negative active-support nodes with
`decay_length_um = 0.05`, `electron_qf_shift_v = -0.04732618818171197`, and
`hole_qf_shift_v = +0.05547180240410299`. It selected `20` support nodes, kept
`34` contact nodes at zero weight, and assigned nonzero smooth weights to `1008`
interior nodes. The shifted state was used as the initial state for a true
DCSweep backscan over `-20, -19, ..., -10 V`; all `11` points converged.

The curve-level result is negative: the smooth branch-control initial state does
not move the high-bias knee. At `-20 V`, current changes from the frozen visual
baseline `-1.1681164e-16 A/um` to `-1.1690341e-16 A/um`, only `0.00034` decades.
The error versus Sentaurus remains `-0.8914` decades. The `-20 -> -19 V` smooth
current-growth ratio is `1.0007`, while Sentaurus is `2.204`; no Sentaurus-like
`-18..-20 V` knee appears. By `-10 V`, the smooth curve is also unchanged to
within `~2e-6` decades versus the frozen visual baseline.

Residual probes explain why this does not become a new branch. Comparing the
zero-shift and smooth-shift `-20 V` states in `newton_residual_probe`, the global
block residual remains Poisson dominated (`psi = 0.2443387`) while carrier blocks
remain tiny. On Sentaurus false-negative support nodes, median absolute residuals
change from `phin = 1.28e-17`, `phip = 7.61e-18` to only `phin = 3.98e-14`,
`phip = 5.19e-14`; the local Poisson residual is unchanged (`5.80e-4` median).
Thus a smooth active-support QF initialization is accepted by Newton but is
pulled back to the same electrostatic/charge-consistent low-current branch.

The next probe should therefore inspect the Poisson/space-charge consistency of
the desired high-density active-support state, not apply stronger pointwise QF
shifts. A useful next experiment is a mixed-state residual audit that combines
Sentaurus-like active carrier densities with Vela potential and measures the
Poisson residual/charge imbalance by active support, contact distance, and doping
sign before proposing any production branch-control policy.
### PN2D BV Mixed-State Charge Audit (2026-06-23)

The follow-up mixed-state charge audit tests whether the desired high-density
active-support branch is blocked by Poisson/space-charge consistency. The new
script `scripts/diagnose_pn2d_bv_mixed_state_charge_audit.py` keeps the Vela
potential and fixed doping, replaces selected active-support carrier densities
with Sentaurus densities, and integrates the resulting mobile/net charge changes
with the mesh control volumes. It reports per-node rows and group summaries by
support class, doping sign, and contact bucket.

For the `-20 V` frozen-Jacobian state, replacing only Sentaurus false-negative
active-support nodes changes the selected mobile/net charge by only
`1.4116e-23 C/m`, while the same selected nodes carry
`3.9116e-11 C/m` of baseline net charge. The ratio is `3.61e-13`. Replacing both
false-negative and false-positive support nodes still changes only
`3.1196e-23 C/m` against the same `3.9116e-11 C/m` baseline, or `7.98e-13` in
absolute-charge ratio. The false-positive compensated nodes show a larger
relative change only because their baseline net charge is near zero; their
absolute change is still `1.7080e-23 C/m`.

This rules out a meaningful Poisson/space-charge obstruction from the
Sentaurus-like active carrier densities themselves. The smooth branch-control
initial state did not move the curve because the coupled solve returns to the
same low-current carrier-continuity/flux branch, not because the desired active
carrier densities would violate electrostatic charge balance. The next useful
probe should therefore inspect carrier-continuity flux/Jacobian balance around
high-field active edges: whether the residual/Jacobian is insensitive to the
absolute QF density lever, whether the SG current path is damping the shifted
state back to the low-current branch, or whether a coupled predictor must move
QF gradients/current and density together.

That follow-up probe has now been executed. The replay first exposed a stale
SG-edge loader assumption: current `sg_avalanche_edges.csv` files write
`electron_flux_abs`, `hole_flux_abs`, and `edge_area_m2`, while the diagnostic
loader only read the older `*_proxy` names. After making the loader accept both
schemas, the real `-20 V` active-edge mixed-state replay reports nonzero active
support (`40` active x-edges per support class). On false-negative support,
Vela baseline generation is `0.1380x` Sentaurus, particle flux is `0.1296x`, and
a uniform Vela QF branch shift using the measured offsets recovers generation to
`0.9618x` and particle flux to `0.9644x`. This proves the local SG source path is
not insensitive to the absolute QF density lever.

The complementary restart-state relaxation probe then captured why the smooth
branch-control deck still leaves the curve unchanged. A one-point `-20 V`
restart from `smooth_branch_state.csv` converges in two Newton iterations to the
same terminal current as the backscan (`-1.169034088445e-16 A/um`). Relative to
the original frozen visual baseline, false-negative support starts with median
QF shifts `phin=-0.047326 V`, `phip=+0.055472 V` and carrier-density boosts
`n=6.238x`, `p=8.548x`; the converged state retains only `phin=-0.009381 V`,
`phip=+0.014140 V`, with `n=1.438x`, `p=1.728x`. The retained absolute QF shift
is only `0.198x` for electrons and `0.255x` for holes. Therefore the active
source gap is causal, but a local absolute-QF initialization is mostly relaxed
away by the coupled carrier-continuity solve. The next experiment should be a
coupled QF-gradient/current-density branch movement or predictor, not a stronger
local density seed or Poisson charge correction.
### PN2D BV Knee-Shape Acceptance Gate

BV parity is not accepted solely because the `-20 V` sweep converges. The next
acceptance gate requires the Vela knee to move toward the Sentaurus knee
window, approximately `-18 V` to `-19 V`, and requires the `-10 V` to `-20 V`
curve to avoid artificial plateaus or early step transitions near `-11 V` to
`-12 V`.

The default `pn2d_sentaurus2018_bv_minus20_avaljac.csv` candidate is not present
in the current `build-release` workspace, so the default `avaljac` knee gate was
not refreshed. The same gate was refreshed with `--output-json` for two existing
candidate curves:

| curve | first 1 V growth ratio > 1.5 | first 1 V growth ratio > 2.0 | max abs log10 current error |
|---|---:|---:|---:|
| Sentaurus | `-19.0 V` | `-20.0 V` | n/a |
| Vela `pn2d_sentaurus2018_bv_minus20_sg_edge_current.csv` | `-11.0 V` | `-13.0 V` | `2.39217` decades |
| Vela dense release curve | `-13.0 V` | `-13.0 V` | `2.4189` decades |

Artifacts:

- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_knee_shape/pn2d_bv_minus20_sg_edge_current_knee_shape.json`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_knee_shape/pn2d_bv_release_dense_curve_knee_shape.json`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_real_state_jacobian_audit_sg_edge_current/knee_shape_current.json`

This keeps the next BV work focused on curve shape and avalanche feedback
parity, not only Newton convergence or local Jacobian completeness.

### PN2D BV Next Physics Decision

No production physics configuration is promoted from this audit. The current
evidence rejects the unsafe options listed in the BV follow-up plan:

- Do not rewrite the core SG flux-divergence discretization as a BV calibration
  step.
- Do not use `source_geometry_scale` as a hidden calibration factor.
- Do not accept a change that only improves one bias point while preserving the
  wrong high-bias curve shape.
- Do not disable SRH for the Sentaurus-faithful BV comparison.

The supporting diagnostics are:

- `tests/test_impact_ionization.cpp` now covers SG avalanche Jacobian
  consistency with low-density quasi-Fermi-to-electric-field interpolation, so
  the current SG source Jacobian already includes the interpolation sensitivity.
- `stepB_loop_gain_sensitivity_m13p2` reports
  `M_electron_qf_vela_over_sentaurus = 1.00734`; the avalanche multiplication
  integral is close, while the carrier generation ratios close to the
  Boltzmann prediction from `~6..9 mV` quasi-Fermi level offsets.
- `active_edge_flux_factors_m13p2` keeps the active-edge particle-flux ratio at
  about `0.73x`, with electric-field and quasi-Fermi-gradient ratios near
  unity. That points at carrier-density/level alignment, not an SG edge flux
  formula mismatch.
- `edge_direction_source_policy_m13p2` and the focused restart variant both
  select active-edge averaging as the closest diagnostic policy, with median
  ratios around `0.76x`, but this remains source-support semantics evidence,
  not a production correction.

The next implementation target should therefore be a real BV-state Jacobian
block export probe and curve-shape experiments around avalanche feedback and
absolute quasi-Fermi level alignment. Reference configuration changes should
wait until they move the `-10 V..-20 V` knee shape toward Sentaurus rather than
only changing a local source magnitude.

### PN2D BV Jacobian Block Probe

`scripts/diagnose_pn2d_bv_jacobian_block_audit.py` now delegates to
`pn2d_jacobian_block_audit` when the C++ probe executable is present. The probe
constructs a small PN coupled-DD fixture and writes finite analytic-vs-FD norms
for `poisson`, `transport`, `srh_auger`, `sg_avalanche`, and
`dirichlet_or_gauge` blocks. It uses a short strong-field fixture for the SG
avalanche block so Van Overstraeten coefficients are active, and a milder
fixture for SRH/Auger so recombination derivatives are not judged on an
avalanche-stiff finite-difference scale.

This is a real assembler-backed Jacobian block audit, but it is still a fixture
probe rather than a replay of the full `pn2d_sentaurus2018` BV restart state.
The remaining upgrade is to export or reconstruct the large BV state and run
the same block decomposition on that state.

The coupled QF-gradient/current-density predictor experiment was then executed as
a stricter version of the branch-control probe. The new helper
`scripts/prepare_pn2d_bv_coupled_qf_predictor_state.py` starts from the Vela
`-20 V` state, selects active SG edges around requested support classes, blends
both endpoints' `phin/phip` to the Sentaurus endpoint QF pattern, and
reconstructs carriers from Vela inferred-ni while keeping `psi` fixed. The
mixed-state replay now accepts explicit `--state-csv-variant` inputs, so the
predictor initial state can be evaluated directly at the active-edge source
level.

For `false_negative` support only, the predictor touches `40` active edges and
`60` endpoint nodes. Its initial active-edge replay is source-effective:
false-negative generation is `1.0126x` Sentaurus and particle flux is `0.9826x`.
However the `-20 V` single-point restart again converges in two Newton iterations
to `-1.169034089095e-16 A/um`, effectively unchanged from the smooth-branch
restart. Endpoint relaxation shows the QF shift is mostly removed: false-negative
support retains only `0.2216x` electron shift and `0.2462x` hole shift, with
carrier boosts reduced from `7.136x/8.158x` to `1.546x/1.677x`.

For the all-support predictor (`false_negative + false_positive`), the helper
touches `66` active edges and `92` endpoint nodes. Its initial replay restores
both active classes near Sentaurus (`false_negative` generation `1.0126x`,
`false_positive` generation `1.0224x`), but the single-point restart still lands
at `-1.169034109498e-16 A/um`. The all-support endpoint relaxation is the same
pattern: retained QF shift remains only about `0.22..0.26x`. This rules out
"move both edge endpoints to the Sentaurus QF/current-density pattern" as a
sufficient branch-control predictor. The next useful probe is the first Newton
linear solve/update on this predictor state: identify which residual block or
Jacobian coupling drives the 75-80% QF-density rollback.

The first-Newton-step audit of the coupled predictor state was executed with
`scripts/diagnose_pn2d_bv_predictor_first_step_audit.py`. The helper converts the
predictor state to runner probe fields, runs `newton_step_probe` and
`newton_block_step_probe`, and reports rollback of the intended
`psi-phin`/`phip-psi` branch shift by support class.

For the false-negative-only predictor, the first full Newton step rolls back
`0.438x` of the electron branch shift and `0.418x` of the hole branch shift on
false-negative support. The carrier-only block step gives the same medians
(`0.438x` electron, `0.418x` hole), while the Poisson-only step is effectively
zero (`~2e-5` rollback fraction). The block residuals show why the solve still
accepts this direction: the predictor state is Poisson dominated
(`psi = 0.2443387`) but carrier residuals are already tiny (`phin = 1.41e-12`,
`phip = 2.01e-12`); the full step combines the Poisson correction with the
carrier block rollback.

The all-support predictor repeats the same pattern. False-negative support has
full-step rollback `0.438x/0.418x` for electron/hole, false-positive support has
`0.437x/0.425x`, and carrier-only matches full Newton. Poisson-only remains near
zero rollback. Thus the first-step rollback is not caused by Poisson/gauge or
Dirichlet movement; it is encoded in the carrier-continuity Newton block itself.
The remaining follow-up should inspect local carrier-row Jacobian coefficients
and RHS signs on the active endpoints, especially why a source-effective QF
branch is treated as a carrier-continuity residual reduction direction back
toward the low-density branch.

The active-endpoint carrier-row audit has now been executed with
`scripts/diagnose_pn2d_bv_predictor_carrier_row_audit.py`. The helper compares
baseline, coupled predictor initial state, and first Newton trial state through
`newton_carrier_row_probe` plus `newton_carrier_term_probe`, so the rollback can
be read at the row/term level instead of only from state deltas.

For false-negative-only support, the predictor increases median electron flux by
`7.09e-14` and hole flux by `9.68e-14` over baseline. The impact term moves in
the compensating direction, but only by `-8.24e-15`, leaving positive residual
deltas of `6.31e-14` electron and `8.85e-14` hole. The raw carrier-row update
therefore rolls the intended branch back by `0.438x/0.418x`, matching the
previous carrier-block step audit. The first trial state reduces the flux and
impact magnitudes, but still has positive residual deltas and `0.341x/0.336x`
raw rollback on the same nodes.

For all-support prediction, false-negative support has median flux deltas
`2.20e-14` electron and `2.07e-14` hole, impact compensation `-8.24e-15`, and
positive residual deltas `1.32e-14`/`1.27e-14`. False-positive support repeats
the pattern with flux deltas `1.88e-14`/`1.54e-14`, impact compensation
`-8.45e-15`, and residual deltas `1.04e-14`/`6.99e-15`.

Conclusion: the active QF branch is source-effective but also overdrives the SG
continuity flux more than the present impact source feedback cancels. Newton is
not discarding a good branch because of Poisson/gauge constraints; the local
carrier residual itself asks for lower carrier density. The next useful probe is
an impact-source feedback sensitivity, ideally at the same carrier-row/term
level, to determine whether the mismatch is impact sign/scale/Jacobian coupling
or an SG flux balance difference.
The impact-source feedback sensitivity was then run on the same carrier-row
artifacts by extending `scripts/diagnose_pn2d_bv_predictor_carrier_row_audit.py`
with `--impact-scale`. The sensitivity keeps the state and SG flux fixed and
asks how much the impact term would need to be multiplied to close each carrier
row residual, using `term_sum + (scale - 1) * impact`.

Results reinforce the carrier-row diagnosis. For the false-negative-only
predictor, the median required scale is `7.97x` electron and `10.26x` hole, and
the first trial still needs `8.47x`/`10.43x`. The single-support predictor is
therefore far from a balanced continuity row even though it restores the local
source proxy.

For the all-support predictor, the required scale is much closer to plausible
feedback-mismatch territory: predictor false-negative support needs `2.31x`/
`2.37x`, predictor false-positive support needs `2.07x`/`1.71x`, and the first
trial remains around `1.61..2.16x`. This makes the next comparison concrete:
inspect whether Vela's impact source contribution is low by about a factor of
two at the active high-field rows because of source support, current weighting,
unit/area scaling, sign convention, or missing Jacobian coupling.
## Impact Feedback Semantic Audit

The production-facing impact feedback audit has now been executed. The local
reference-code check found no sign reversal in the basic generation semantics:
Genius DDM forms `GII = alpha_n * |Jn| / q + alpha_p * |Jp| / q`, distributes it
with directional current weights, and injects the same generated pair source into
electron and hole continuity rows over the finite-volume truncated partial
volume. Charon computes `alpha_n * |Je| + alpha_p * |Jh|` and subtracts avalanche
from total recombination, i.e. treats it as generation. Vela likewise subtracts
the SG edge-current avalanche source from both carrier residual rows. The local
`devsim` tree did not expose an equivalent built-in impact-ionization assembly to
compare.

The concrete all-support predictor comparison is now:

- Active-edge replay: predictor false-negative generation is `1.01257x`
  Sentaurus and false-positive generation is `1.02244x` Sentaurus on the selected
  active x-edges.
- Source geometry replay, after updating
  `scripts/diagnose_pn2d_bv_source_geometry.py` to accept current C++
  `sg_avalanche_edges.csv` columns (`edge_area_m2`, `electron_flux_abs`,
  `hole_flux_abs`, `source_integral`), reports active endpoint area fraction
  `0.5` for both support classes.
- Multiplying those two facts gives effective active-support feedback of
  `0.506x` on false-negative and `0.511x` on false-positive nodes, matching the
  carrier-row sensitivity that needed roughly `2.31x/2.37x` and `2.07x/1.71x`
  electron/hole impact feedback to close the rows.

Artifacts:

- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/source_geometry_all_support_blend1/`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/active_edge_replay_all_support_blend1/`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/impact_scale_sensitivity_all_support_blend1/`

Conclusion: the all-support predictor's local active-edge avalanche physics is
already Sentaurus-sized; the remaining factor-of-two is explained by effective
finite-volume/support feedback into the carrier rows, not by a raw avalanche
coefficient sign or local QF branch-strength error. Do not introduce a hidden
source multiplier. The next task should inspect whether Vela's endpoint
`0.5 * edge_area` carrier-row injection should be compared to Genius/Sentaurus
truncated partial-volume ownership on active edges, and separately whether the
SG edge-current avalanche Jacobian should include source derivatives rather than
omitting them in the current Newton path.
## Impact Feedback Ownership Policy Probe

A follow-up ownership summary has been executed with the new helper
`scripts/summarize_pn2d_bv_impact_feedback_ownership.py`. The helper joins three
already generated all-support predictor diagnostics:

- active-edge replay generation ratio versus Sentaurus,
- source-geometry active endpoint area fraction,
- carrier-row impact scale needed to close the local residual.

Artifact:

- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/impact_feedback_ownership_all_support_blend1/`

Real `-20 V` all-support results:

| support class | active-edge generation / Sentaurus | active endpoint area fraction | endpoint feedback / Sentaurus | full active-edge feedback / Sentaurus | required e/h impact scale |
|---|---:|---:|---:|---:|---:|
| false_negative | `1.01257x` | `0.5` | `0.50628x` | `1.01257x` | `2.30979x / 2.36797x` |
| false_positive | `1.02244x` | `0.5` | `0.51122x` | `1.02244x` | `2.07195x / 1.71073x` |

This makes the factor-of-two failure mode concrete: the local active-edge source
strength is already approximately Sentaurus-sized, while endpoint-half ownership
leaves the carrier row with only about half of that active feedback. The product
`endpoint_feedback * required_scale` lands near unity (`1.17/1.20` for
false-negative and `1.06/0.87` for false-positive), so the previous `~2x`
impact-scale sensitivity is consistent with source-support ownership rather than
a raw ionization coefficient error.

Current C++ status: the SG edge-current avalanche residual path injects a source,
but the analytic Jacobian path still explicitly omits those nonlocal edge-source
derivatives. That is a separate Newton-coupling issue. The next production probe
should therefore be ordered as:

1. First test a focused source-ownership variant for SG edge-current avalanche,
   comparing Vela endpoint-half ownership against a directional/truncated-volume
   policy analogous to Genius/Sentaurus active-edge ownership. This must remain a
   gated experiment, not a hidden multiplier.
2. Then test source-derivative Jacobian completion for the SG edge-current path
   against finite-difference block probes.
3. Accept neither change without the curve-level `-10..-20 V` knee-shape gate and
   local carrier-row residual gate moving in the same direction.
## SG Edge-Box Source Volume Probe

Executed the focused production probe proposed above by adding an explicit gated
configuration knob, `impact_ionization.source_volume_policy`. The default remains
`edge_half_box`, preserving the previous SG edge-current avalanche source
support. The probe value `edge_box` changes only the SG edge source support from
`0.5 * h * edge.couple` to `1.0 * h * edge.couple`; it is not a hidden avalanche
coefficient multiplier and does not add the still-missing SG source-derivative
Jacobian terms.

Real `-20 V` all-support single-point deck:

- Baseline config: `single_m20_all_support_blend1/simulation_bv_coupled_qf_predictor_all_support_blend1_single_m20.json`
- Probe config: `single_m20_all_support_edge_box/simulation_bv_coupled_qf_predictor_all_support_edge_box_single_m20.json`
- Sentaurus reference current at `-20 V`: `-9.10455666344e-16 A`

| case | converged | Newton iterations | current_total_A_per_um | abs ratio vs baseline | abs ratio vs Sentaurus | decade error vs Sentaurus |
|---|---:|---:|---:|---:|---:|---:|
| endpoint-half baseline | `1` | `2` | `-1.16903410949798e-16` | `1.0000x` | `0.128401x` | `-0.891432` |
| `edge_box` probe | `1` | `2` | `-6.39910455999440e-16` | `5.47384x` | `0.702846x` | `-0.153140` |

The max field is effectively unchanged (`560748.547233267` to
`560748.547233097 V/cm`), so the terminal-current movement is carrier-continuity
feedback from the source ownership policy rather than a field-state movement. The
probe closes about `0.738` decades of the `0.891` decade gap and leaves a
remaining Sentaurus multiplier of `1.42279x`.

Conclusion: source ownership is now confirmed as a production-relevant direction,
but the `edge_box` probe alone is not yet an acceptance change. The next ordered
work is: run the same gated policy through the `-10..-20 V` knee-shape gate and
local carrier-row residual gate, then independently test SG edge-current
avalanche source-derivative Jacobian completion against finite-difference block
probes.

## SG Edge-Box Backscan Knee And Carrier-Row Gate

The `source_volume_policy=edge_box` experiment was extended from a single `-20 V`
point to a matched `-20 -> -10 V` all-support predictor backscan. Two decks were
generated from the same all-support restart family:

- endpoint-half baseline: `all_support_blend1_backscan/simulation_bv_coupled_qf_predictor_all_support_blend1_backscan_m20_to_m10.json`
- edge-box probe: `all_support_edge_box_backscan/simulation_bv_coupled_qf_predictor_all_support_edge_box_backscan_m20_to_m10.json`

Both backscans converged with `11` integer-bias points. The knee-shape gate now
reports:

| curve | first 1 V growth ratio > 1.5 | first 1 V growth ratio > 2.0 | max abs log10 current error |
|---|---:|---:|---:|
| Sentaurus | `-19.0 V` | `-20.0 V` | n/a |
| endpoint-half baseline backscan | none | none | `0.891432` decades |
| `edge_box` backscan | `-16.0 V` | `-19.0 V` | `0.519580` decades |

The probe therefore moves the curve-level knee in the right direction but does
not pass acceptance: the `>2.0` growth point is still one volt early, and the
`>1.5` threshold appears too early at `-16 V`. Pointwise, `edge_box` is close at
some mid-window biases but overshoots around `-19 V` (`3.308x` Sentaurus) while
remaining low at `-20 V` (`0.703x` Sentaurus).

A local carrier-row audit was also run under the `edge_box` policy with the
matched `-20 V` states:

- `carrier_row_audit_all_support_edge_box_policy/`

On the active support classes, the endpoint-half state evaluated with `edge_box`
already has required impact scales below unity (`false_negative` e/h
`0.8898/0.8774`, `false_positive` e/h `0.7784/0.7033`). The converged `edge_box`
state pushes the same rows further source-strong (`false_negative` e/h
`0.6452/0.5778`, `false_positive` e/h `0.6568/0.4142`) and lowers carrier raw
rollback magnitudes. This confirms source ownership is a real lever, but the
full-edge policy is too coarse as a production acceptance change.

Next ordered task: introduce and test an intermediate/truncated ownership factor
rather than binary `0.5` versus `1.0`, or move to the independent SG source-
derivative Jacobian probe if the goal is Newton coupling rather than curve-shape
calibration. Any intermediate policy must be gated by the same backscan knee
summary and carrier-row required-scale summary.
