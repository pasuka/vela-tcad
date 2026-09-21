# SLOT-LDMOS BVDS validation

This is the current qualification reference. Read it before the dated SLOT
handoff and execution reports. Evidence is local to
`build/low_voltage_20260906/`; generated simulation states are not checked in.

## Qualification profiles

The historical
[v1 closure profile](contracts/slot_ldmos_physical_closure_restoration_v1.json)
qualified global closure through 2 V. Its carrier source threshold selected
zero local rows. The later all-active audit found large cancelling local
residuals, so those points do not qualify spatial continuity.

The opt-in
[v2 all-active profile](contracts/slot_ldmos_all_active_local_v2.json)
requires full coupled convergence, global closure, and every active electron
and hole continuity row to satisfy `eps_row = 0.001`. It explicitly uses the
existing node-local contact doping boundary and the new state-based QF
reference coordinates. It is not an unchanged-physics continuation of v1.
Global solver defaults and the PN2D BV atomic profile are unaffected.

The continuation driver keeps coupled Newton iterations separate from bounded
carrier-only correction steps. A correction probe only produces a candidate;
the complete original coupled deck and a fresh all-active audit must accept
it before the continuation state advances. A zero-step same-bias reload is
revalidation, not additional Newton work. Failed attempts and the largest
attempted voltage do not establish a qualified voltage range.

## Current evidence

- A continuous qualified state chain covers 0--2 V, with 42 positive-bias
  points: `qualified_0_to_2_chain.json`. It combines the accepted 0--0.5 V
  prefix of `reproducible_0_to_2/` with `dense_0p5_to_2/`. The first run
  exhausted its attempt budget at 0.75 V and has `completed: false`.
- `single_run_0_to_2_replay/` independently replays the full accepted bias
  list in one invocation: completed, 42 positive-bias points, zero rejected
  attempts, all 42 saved states byte-identical to the qualified chain.
  `repeatability.json` and the final summary record this result;
  `invocation.json` records arguments and the driver hash.
- Independent Si intrinsic-density controls at 1.25, 1.5 and 2 V pass both
  local and global gates. Their current differences from the same-Masetti
  Sentaurus reference are -1.965%, -2.063% and -2.748%.
- `intrinsic_aligned_2_to_10/` and `sentaurus_noav_10/` hold the higher-voltage
  no-avalanche controls. The 3 V and 5 V local checkpoints pass all gates,
  with current differences -3.0091% and -4.8941%. The original local run
  stopped at qualified 5.728515625 V; shorter-path recovery `build/h7`
  reached 9.5625 V before the requested Release switch. `build/h8` then
  completed through 10 V, with full coupled, global and all-active local
  qualification. At 10 V, Id is `2.9274254645767804e-15 A/um`, versus
  Sentaurus `2.98614241990607e-15 A/um`; maximum local ratio is
  `0.0008471901590258922`. Its final summary records 27 attempts, including
  18 accepted points. Surface mobility, avalanche and the external resistor
  require separate qualification.
- The independent original-coefficient IALMob controls `build/in2/in3/in5/in10`
  pass the same gates at 2/3/5/10 V after the true-distance normal correction.
  Their current differences from original Sentaurus IALMob are -4.5756%,
  -4.8391%, -6.9082%, and -3.9737%. Separate avalanche activation controls
  `ia2r/ia3/ia5` pass at 2/3/5 V. Continuous IALMob replay, source-current
  alignment and the external load line remain separate work; none of these
  checkpoints establishes BVDS. See the IALMob report below for exact values
  and the original `AvalDensGradQF` versus local SG/GSS-Laux distinction.

See the [local error budget](slot_ldmos_local_error_budget_2026-09-06.md) for
configuration, units, source evidence, numerical changes and tests, and the
[mobility derivative report](slot_ldmos_mobility_execution_2026-09-06.md) for
the endpoint, density and adjacent-cell derivative coverage.

## Physical scope and limits

The mesh has 13,881 nodes and 27,275 triangles. Runs use 300 K, mixed-Voronoi
node volumes, Masetti plus GradQF high-field saturation, SRH with doping
dependence, Auger and OldSlotboom. Newton uses Eigen SparseLU/COLAMD. State
CSV densities are SI `m^-3`; terminal comparison currents are `A/um`.

The base v2 profile retains `ni = 1e10 cm^-3` for Si. The separately named
intrinsic-aligned variant uses `1.0749877484352535e10 cm^-3`. Both retain the
electrostatic-only PolySilicon simplification and metal-gate boundary. The
Sentaurus reference has mobile, doped PolySilicon and a material energy
reference offset, so neither variant is a fully identical physical model.
The partial PolySilicon transport control must retain that limitation.

The original Sentaurus deck uses IALMob. The new local `ialmob` and
`ialmob_field` implementations are separate from Enhanced Lombardi. Model
integration and spatial parity are still under validation; historical
Lombardi runs remain labeled Lombardi. See the
[IALMob implementation and comparison](slot_ldmos_ialmob_alignment_2026-09-07.md).
Matching a model name or obtaining a nearby terminal current does not
establish formula, parameter or spatial-field parity.

Fixed-state volume/source sensitivities, self-consistent material controls,
and solver stopping controls are distinct experiments. Do not add their
percentage changes into one causal budget. Do not infer a verified BVDS or
post-breakdown branch from a no-avalanche sweep.
