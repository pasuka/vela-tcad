# Auxiliary gate QS large-step control

This control uses the same authorized drain 0 -> 13000 V, gate 0 -> 0.2 -> 0 V,
then drain -> 100000 V path as [the initial candidate](../gate_qs_seed_v1/README.md).
Only the two gate ramps change: InitialStep and MaxStep are 1 instead of 0.05;
local `Extrapolate(MinStep=1)` avoids prediction below a full normalized step.
Failed attempts can still reduce the step down to the original MinStep.
Models, parameters, mesh, fixed resistor, stopping threshold, acquisition
and final drain ramp are identical to the initial candidate.

Source SHA-256:
`53e98eef85df96943e70a58e178304f175bbe66b2de1ea9efe3ae4951b11b455`.
Use the same independent auxiliary-path checker and separate qualification.
Results: [zero-QS diagnosis](../../../../docs/validation/bvmethods_zero_qs_diagnosis_2026-10-06.md).

Executed on 2026-10-06: failed. The gate ramp reached only 0.0966954759244808 V;
the final drain current stopped at 0.00140926509666232 A/um. All 191 accepted
states passed the physical closure checks, but the full path failed. Raw
data are sealed in `bvmethods_zero_qs_controls/v2`; current v4 is unchanged.
