# Preserve referenced QF states in circuit and arclength restarts

The authorized high-voltage continuation exposed a confirmed restart defect:
`packArclengthState` packs rounded absolute QFs, and arclength factories build
reference fields without the seed. The external circuit additionally constructs
the device with zero base drain bias even when its initial inner bias is 37 V.
Ordinary Newton diagnostics already preserve split references and increments.

## Design

Give the arclength assembler factories and state conversions an optional common
seed reference state. Each factory resolves its immutable reference fields when
called; callbacks must not retain a pointer to caller-owned seed storage. Reuse
the ordinary compensated state-packing helper. Construct the circuit device
with its actual initial inner contact bias, and pass one immutable seed to the
system, packer, unpacker, and terminal-current functional. Use an immutable
anchor similarly for ordinary pseudo-arclength continuation.

This fixes initial handoff representation. Fixed reference coordinates do not
promise arbitrary sub-ULP resolution after a large continuation-voltage change.
No physical coefficients, residual/Jacobian formula, acceptance tolerances,
linear backend, or default density-update policy change.

## Execution and ownership

- [x] Add a failing numerical regression for a nonzero reference with sub-ULP
  carrier increments; record the old packer's loss before implementation.
- [x] Update NewtonSolver factories/conversions and test seed-aware references,
  residual/current consistency, and legacy API compatibility (reference agent).
- [x] Integrate the shared anchor and actual initial contact bias into both
  DCSweep consumers (root).
- [x] Build Release, run the relevant numerical regressions and full Release
  CTest suite because this changes shared solver interfaces.
- [x] Independently review coordinate consistency and seed lifetime.
- [x] Freeze a new Release runner and repeat the identical external 37 V step
  control; preserve failed old-runner evidence and all qualification gates.
- [x] Record outcomes and commit locally. An initial representation repair is
  not a completed BVDS qualification.

Numerical evidence lives under `build/slot_resume_20260922/`; raw simulation
outputs stay ignored. Existing native Grad and SG resistor references remain
separate, with no change to their decks or reference currents.

## Verified results

The original roundtrip test failed both sub-ULP increment checks. After the
repair, 18 arclength/QF-reference test cases and 165 assertions passed, and the
full Release CTest suite passed 890/890 (140.78 s). The real SLOT 37 V replay
restored the initial current to within 9.22e-15 relative of its qualified seed,
replacing the old +18.7777% drift. The same positive-density circuit step still
failed subsequent line search; this is a completed repair verification, not a
qualified high-voltage continuation point. Evidence is in
`build/slot_resume_20260922/arclength_restart_representation_audit/`.

## Review and Poisson follow-up before the final commit

Independent review found that the new arclength packer could overwrite explicit
absolute-QF edits with stale split fields. Match ordinary Newton's per-node,
per-carrier 32-epsilon consistency check, falling back to the supplied absolute
value when inconsistent. Opt in only the restart consumers so unrelated
diagnostic packing semantics remain unchanged.

The actual 20.1/20.25-to-20.5 V predictor retains split increments exactly, but
`solvePoissonOnly(preserve_initial)` discards 102 free-Si electron increments
while keeping their reference and absolute QF unchanged. Use the same guarded
restart packing after contact projection. `contact_basin` must continue to
reset QFs intentionally; no reference-selection or acceptance policy changes.

- [x] Add/run RED tests for Poisson sub-ULP preservation and stale absolute
  arclength edits, with both scaling modes and a contact-basin reset check.
- [x] Implement the guarded restart helper and Poisson preservation path.
- [x] Repeat targeted and full Release verification on the final source.
- [x] Review the repair delta and replay the actual 20.5 V prewarm.

The 890/890 result above is the earlier arclength implementation's evidence;
the follow-up changes now have their own verification: 23 targeted cases / 221
assertions and full Release CTest 892/892 pass (145.44 s). Independent review
found no remaining actionable issues. The same-input 20.5 V Poisson replay
retains all 102 formerly zeroed increments: 85 are bitwise identical, and 17
differ by one increment ULP from the scaled-coordinate roundtrip (maximum
1.9722e-31 V). This is representation evidence, not coupled qualification.
The fresh frozen runner SHA-256 is
`bc2de24c077584a910511e33ef32002876adbad4b2294825d29af4ba43c6f6a6`.
Evidence: `continuation_reference_audit/{green.log,full_ctest.log,final_manifest.json}`
and `preserved_prewarm_20p5/result.json` under the evidence root above.
