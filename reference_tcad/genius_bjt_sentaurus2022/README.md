# Genius NPN BJT / Sentaurus 2022 reference

This fixture reproduces the two-dimensional NPN BJT shipped with Genius TCAD
(`examples/BJT/step1.inp` and `step2.inp`) as a Sentaurus SDE/SDevice project.
It is the WP0-WP2 oracle for a later, separately scoped Vela import and
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

