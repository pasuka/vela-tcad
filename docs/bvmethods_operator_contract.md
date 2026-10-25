# BVMethods explicit coupled operator

`solver.coupled_operator: "bvmethods_tri3_qf"` selects an opt-in coupled
operator. The default remains `"legacy"`. This profile does not change the
PN2D BV template bundle. Its numerical qualification is recorded in
[the implementation report](validation/bvmethods_production_profile_2026-10-24.md).

A [full Release sweep](validation/bvmethods_full_profile_2026-10-25.md) reached
the configured stop current using the same core library and a local continuation
driver with a predicted inner-voltage step guard and accepted-state retry.
The default CLI continuation has not acquired that strategy. Native zero-state
operator consistency and breakdown-region step convergence remain separate
qualifications.

## State and statistics

The unknowns are ordinary physical electrostatic and electron/hole
quasi-Fermi potentials, packed as `(psi, fn, fp)/V0`. Only zero-reference
coordinates are qualified; reference-tail and split-DD coordinates are rejected.
Populations are reconstructed at every evaluation from these potentials, the
material DOS, effective intrinsic density and BGN offsets. No imported population,
native SG cache, source array or native Jacobian is a constitutive input.

Select `solver.carrier_statistics.model: "fermi_dirac_joyce_dixon"` explicitly.
The Joyce-Dixon forward, inverse and derivatives form one contract, including
contact neutrality, density reconstruction, FD-BGN and generalized SRH. The
existing `fermi_dirac` defining-integral model remains unchanged. Selecting
Joyce-Dixon in the Newton JSON parser outside this profile is rejected.
Material files and the global CODATA constants remain explicit inputs; this
profile does not substitute Sentaurus thermal-voltage or charge constants.

## Discrete operator

- Poisson and contact/gauge rows retain the existing cell-material assembler.
- Carrier transport uses Tri3 coefficient-derived box weights, total-impurity
  Masetti mobility, distance-gradient Enhanced Lombardi and Canali field
  mobility with live derivatives.
- SG is evaluated with log-density and `expm1` quasi-Fermi differences and the
  matching generalized Einstein factor. Flat quasi-Fermi fields give zero
  equilibrium current algebraically; no small-current cutoff is used.
- AV uses the representative vertex and its two incident SG edges to recover
  the current vector, the existing contact driving-field rule, an unsmoothed
  Euclidean norm and vertex-box source integration. The same pair-generation
  residual enters both continuity equations. Local jets differentiate that
  residual; no captured current direction is supplied.
- SRH uses the same statistics and coefficient-derived silicon control volume.
- The residual, full Jacobian and carrier-term ledger share this runtime.
  Terminal totals and their state derivatives use its source-free transport
  residual directly. Drift/diffusion partitions and legacy edge diagnostics
  are not qualified for this operator. Legacy drift/diffusion CSV columns
  retain default placeholders; use carrier totals and terminal total only.

At exactly zero current the norm has directional derivatives, not a unique
Frechet derivative. The existing zero-norm jet convention is explicit. A finite
difference crossing that point cannot establish the local derivative of an
arbitrarily small nonzero native current. Native zero-state bitwise agreement
is a separate qualification from autonomous nonlinear convergence.

## Qualified configuration boundary

Version 1 requires classical 300 K silicon Tri3 DD, cell-material Poisson
geometry, `masetti_field_lombardi`, `total_impurity`,
`quasi_fermi_gradient` mobility with field derivatives, and SRH or no
recombination. AV is either disabled (for source-free transport evaluation) or
self-consistent Van Overstraeten with GradQF/current-density generation and
`local_ad`. AV cutoffs, reference-density interpolation and source overrides
are rejected. Quantum, density diagonal floors, thermionic contacts, Auger,
BTBT and split DD are not qualified. The DC sweep rejects VTK and legacy
contact-edge/QF-floor diagnostics, whose field/current operators differ.

The profile is available through JSON and C++ configuration. SDevice deck
translation does not automatically enable it. To return to the former physics
path select both `coupled_operator: "legacy"` and the previous statistics model.
Newton and load-line stopping rules and budgets are unchanged.

## Verification entry points

`test_bvmethods_av_contract` contains `[bvmethods_profile]` tests for full Jv,
source pairing, contact rows, ledger identity, equilibrium and domain guards,
and terminal current. `test_joyce_dixon_contract` covers the public statistics
and FD-BGN dispatch.

The uninstalled `scripts/diagnostics/bvmethods_profile_probe.cpp` links against
the current Release `vela_core`. It exports zero and finite-current operators
and drives the production load-line Newton API through multiple outer voltages.
Build it explicitly with
`cmake --build --preset windows-ucrt64-release --target vela_bvmethods_profile_probe`.
It is excluded from the normal build and installation.

The full-sweep controller is retained separately in
`scripts/diagnostics/bvmethods_full_sweep_probe.cpp`. Build its uninstalled target
with `cmake --build --preset windows-ucrt64-release --target vela_bvmethods_full_sweep_probe`,
then run `build-release/vela_bvmethods_full_sweep_probe manifest.json output-directory`.
The trusted local manifest contains `config` (a JSON configuration path),
`previous` and `current` physical checkpoints, and `initial_step_V`.
Each checkpoint contains `outer_V`, `inner_V`, and mesh-order arrays
`psi_V`, `fn_V`, `fp_V`. Use Vela's own accepted states, retaining input hashes
and rechecking the restart anchor; the manifest does not supply populations,
currents, AV sources or a Jacobian.

This frozen case-specific controller uses 1e7 ohm um, a predicted inner-voltage
step cap of 0.02 V, outer step bounds 1e-5--2500 V, growth at most 1.41,
halving after both seed attempts fail, and 40 Newton iterations per attempt.
Its strict equation/load tolerances are 1e-8 and 1e-6 V. It stops at
`abs(Id) >= 1.3e-3 A/um` or outer voltage 100000 V, with a 512-attempt budget.
It records every completed attempt and accepted physical state. These constants
are diagnostic controls, not newly configurable production defaults.
`scripts/run_bvmethods_profile_regression.py` validates the state-only manifest
and rejects an incomplete or failed continuation. Every rejected solve remains
in the evidence; step halving retries from the last accepted state with the same
iteration budget and stopping rules.
