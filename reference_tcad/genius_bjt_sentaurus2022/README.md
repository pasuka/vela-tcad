# Genius NPN BJT / Sentaurus 2022 reference

This fixture reproduces the two-dimensional NPN BJT shipped with Genius TCAD
(`examples/BJT/step1.inp` and `step2.inp`) as a Sentaurus SDE/SDevice project,
then imports the accepted structure into Vela for the WP3-WP5 common-input
comparison.

## Frozen source contract

- Device: 6 um x 2 um silicon rectangle, 300 K.
- Genius extrusion width: 100 um. Sentaurus and Vela comparisons use current
  per unit width in A/um; the 100 um width is metadata, not a current scale
  silently applied to either solver.
- Top contacts: base x=[1.25,2.00] um and emitter x=[2.75,4.25] um.
- Bottom contact: collector x=[0,6] um.
- Bias: emitter=0 V, base=0.70 V, collector=0..3 V in 0.1 V increments.
- Doping is evaluated with the exact Genius analytical-profile definition in
  device coordinates. Outside each rectangular window the profile decays as
  `exp(-(distance/char)^2)`; inside it is flat. Overlapping profiles add.

The machine-readable form is in `contracts/geometry_contract.json` and the
physics split is in `contracts/physics_contract.json`.

## WP0-WP2 contents

- `source/bjt_sde.cmd`: SDE geometry, contacts, analytical doping, and mesh.
- `source/bjt_m0_des.cmd`: minimal common drift-diffusion model.
- `source/bjt_m1_des.cmd`: realistic common drift-diffusion model.
- `contracts/`: immutable input, physics, and acceptance contracts.
- `reference_curves/`: normalized all-terminal DC curves produced from the
  Sentaurus PLT files.
- `reports/`: structure, convergence, current-sign, and KCL audits.
- `manifests/`: source hashes, tool versions, run identity, and output hashes.

Large/raw Sentaurus products are deliberately kept under
`build-release/reference_tcad/genius_bjt_sentaurus2022/` and are not intended
for source control.

## WP3-WP5 Vela contents

- `vela/input/`: the exact 5611-node/10940-triangle Sentaurus structure and
  node doping converted to Vela input.
- `vela/materials_sentaurus2022.json`: explicit 300 K silicon material values.
- `vela/configs/`: M0 base ramp and collector sweep, plus M1 fixed-bias model
  relaxation, collector sweep, and reproducible 3 V spatial export.
- `comparison/`: exact 31-point all-terminal comparison tables and the
  machine-readable/human-readable summaries.

The M1 path deliberately relaxes the converged M0 state at VBE=0.70 V before
its collector sweep. A direct M1 ramp from zero bias was rejected near 2 mV by
the continuity line search after the adaptive step fell to approximately
1e-8 V. The fixed-bias model ladder reaches the identical requested boundary
condition without changing mesh, doping, contacts, or M1 physics.

## Model ladder

M0 intentionally uses Boltzmann statistics, default low-field mobility, and
plain SRH recombination. M1 enables Fermi statistics, OldSlotboom band-gap
narrowing, doping-dependent mobility, doping-dependent SRH, and Auger. Neither
model enables high-field mobility or impact ionization; that keeps the first
cross-solver comparison focused on equilibrium, injection, transport, and
terminal-current bookkeeping.

The collector sweep is solved only after Poisson, equilibrium coupled DD, and
a base ramp from 0 to 0.70 V. Currents for collector, base, and emitter are
read directly from Sentaurus; base current must never be reconstructed by KCL.

## Validated WP0-WP2 result

The accepted run is `genius_bjt_wp0_wp2_20260902/run_v2` on host `tcad` with
Sentaurus T-2022.03-SP2. It produced a 5611-node, 10940-element mesh and 31
collector-bias points for each model. Structure, voltage-grid, direct-terminal
current, and KCL gates all pass; see `reports/wp0_wp2_validation.json`.

Sentaurus terminal current is stored with the solver sign convention: positive
current flows into the device. Consequently, forward-active collector and base
currents are positive and emitter current is negative. The three directly read
terminal currents sum to zero within the recorded numerical residual.

Run-local raw TDR, PLT, and log files are retained outside source control at
`build-release/reference_tcad/genius_bjt_sentaurus2022/wp0_wp2_20260902/run_v2`.
The remote copy is under
`/root/sentaurus_runs/vela_oracle_2022/genius_bjt_wp0_wp2_20260902/run_v2`.

## Validated WP3-WP5 terminal and spatial result

Both Vela collector sweeps reached 3 V with 74 accepted/internal points. The
comparison tool extracts only the exact 31 requested 0.1 V points and reads
collector, base, and emitter currents directly. Convergence, voltage matching,
and three-terminal KCL pass for both model ladders.

At VCE=3 V, M0 gives an Ic magnitude ratio Vela/Sentaurus of 1.8431 and beta
values of 476.64 versus 114.93; it remains a characterization-only numerical
baseline. The original M1 base-current excess was traced to an electron
Scharfetter `tau_max` of `3e-8 s` instead of the T-2022.03-SP2 Silicon default
`1e-5 s`. After correcting that input and aligning the effective density of
states to the T-2022.03-SP2 Silicon values (`Nc=2.8567e19 cm^-3`,
`Nv=3.1046e19 cm^-3`), the unique accepted-state chain gives Ic and Ib ratios
of 0.992152 and 0.995171 at 3 V, with beta values of 46.505 versus 46.646.

The asserted M1 gate covers VCE=0.5-3.0 V and limits the maximum absolute
log10 error of Ic, Ib, Ie, and beta to 0.05 decades each. Observed maxima are
0.003437, 0.002102, 0.003407, and 0.001349 decades.

The spatial-state gate uses the exact common mesh at VBE=0.70 V and VCE=3.00 V.
Potential is checked on all 5611 nodes. Electron and hole decade errors are
checked where the Sentaurus reference density is at least 1e10 cm^-3. The
potential, electron-density, and hole-density gates all pass. The former
1.6836-decade full-domain hole outlier is confined to nodes with reference
hole concentrations of roughly 1-100 cm^-3; the asserted hole-density maximum
is 0.10079 decade.

The accepted-state VTK path now exports separate SRH, Auger, total electron and
hole current density, and legacy-scale Vela drift/diffusion diagnostics. At 3 V the electron
current-density gate passes. The hole current-density gate fails, while SRH
passes its magnitude, integral, and shape checks but misses the peak-location
tolerance by 0.021875 um. Auger fails magnitude, integral, and peak-location
checks. Consequently the final `overall_pass` is now false and depends on all
terminal, KCL, spatial-state, current-density, SRH, and Auger gates. See
`contracts/comparison_thresholds.json`, `comparison/comparison_summary.md`, and
`reports/transport_source_acceptance_report.md`.

Vela representative exports cover VCE=0, 1, 2, and 3 V. Existing SDevice TDRs
support direct comparisons at 0, 1, and 3 V; the 2 V SDevice state remains
missing because the earlier deck saved 0.1 V instead and the VM currently
rejects the documented SSH key. The updated deck now requests 0, 1, 2, and 3 V
states on its next trusted rerun.

## Comparison figures

Publication-ready mesh, M1 spatial-field, terminal-curve, and parity-error
figures are stored in `figures/`. Their source/output hashes and field-error
metrics are recorded in `figures/figure_manifest.json`; see
`figures/README.md` for plotting conventions and the reproduction command.

## Hole-density solver diagnosis

The historical pre-Nc/Nv-alignment full-domain low-hole outlier is localized in
`hole_density_diagnosis/`. Five exact-mesh state-difference maps show that the
1.6836-decade density error is explained by a -100 mV `hQF-psi` plateau to
within 0.0003 decade RMSE where the Sentaurus hole density is below 1e2 cm^-3.
Carrier-term and Newton-step probes show that a further limited +0.1 V hQF
update removes this plateau, whereas disabling the update limit is unstable.

Strict `carrier_row_convergence=enforce` with `eps_row=1e-4` removes the
original plateau but does not converge globally: after 80 iterations, 148
negligible-density rows still violate the local relative criterion. This is a
diagnostic result, not a replacement baseline. Reproduce the solver probes and
figures with:

```powershell
D:\msys64\ucrt64\bin\python.exe scripts\run_genius_bjt_hole_density_diagnostics.py --phase all
D:\msys64\ucrt64\bin\python.exe scripts\diagnose_genius_bjt_hole_density.py
```
