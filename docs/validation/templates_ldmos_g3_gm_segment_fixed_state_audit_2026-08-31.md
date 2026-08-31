# Templates/LDMOS G3 maximum-gm segment fixed-state audit (2026-08-31)

## Question and frozen contract

The 31-point G3 qualification misses only the frozen maximum-gm limit. Both
engines place the maximum on the exact `Vg=1.0--1.16666666666667 V` segment,
so this audit asks whether the `24.3076%` slope difference is produced by the
qualified carrier operator or by its self-consistent state.

Both endpoint replays use the same G3 contract as the 31-point run:

- `Vd=0.1 V`, Fermi statistics, OldSlotboom BGN and SRH/Auger;
- IALMob, predictor and quantum potential disabled;
- interior GradQF plus qualified contact-cell ElectricField HFS fallback;
- `legacy_node_local` source-short reconstruction;
- external AverageBox carrier couples and barycentric node volumes.

For each endpoint, `SSS` applies the Vela carrier operator to the complete
Sentaurus `psi/phin/phip` state. `VVV` applies the identical operator to the
converged Vela state. The remaining six combinations replace one or two state
families without changing coefficients, mobility parameters or acceptance
thresholds.

## Sentaurus endpoint capture

The existing G3 state deck ended at `Vg=1.0 V`. An isolated T-2022.03-SP2 VM
run extended only the gate goal and exact plot list to `1.166666666666667 V`.
The mesh, parameter file and physics block were unchanged. The new endpoint
converged with drain current `2.82791611008207e-6 A/um`; the full run used
`246.81 s` wall time and `339 MB` peak memory.

Provenance checks:

- input `n1_fps.tdr` SHA-256:
  `07af6353c47c03a740583b873bac58e733a377a8bb9fa7e08baa3d30bc2c67f3`;
- extended deck SHA-256:
  `a0dd94a6219602ceba5ab3055d56fde69f3b95b4da94adedbb196a033a438f36`;
- endpoint TDR SHA-256:
  `e7840a6a314d43cff4060a2a76a2289f9667c236f9d6b764645266658df11745`;
- downloaded gzip SHA-256:
  `b6fc018bbd7c37f44df7c4f314d04f0bb6fc3de33198c178fe9ffa70f2a301a8`.

The imported endpoint contains 44 field datasets. No TDR postprocessed
mobility field is injected into Vela; mobility is recomputed by the frozen
qualified operator from each selected state.

## Endpoint operator and state decomposition

| Vg (V) | Sentaurus Id (A/um) | SSS Id (A/um) | SSS/Sentaurus | VVV Id (A/um) | VVV/Sentaurus |
|---:|---:|---:|---:|---:|---:|
| 1.000000 | 1.254155163e-6 | 1.254408686e-6 | 1.000202146 | 1.021797549e-6 | 0.814729772 |
| 1.166667 | 2.827916110e-6 | 2.828487876e-6 | 1.000202146 | 2.213015335e-6 | 0.782560391 |

The fixed-state operator error is `8.7782e-5 dex` and `8.7799e-5 dex` at the
two endpoints. Its ratio is effectively constant, so it cannot generate the
gm mismatch. The endpoint slopes make the separation direct:

| Segment slope | A/(um V) | Ratio to Sentaurus | Relative error |
|---|---:|---:|---:|
| Sentaurus terminal | 9.442565684e-6 | 1.000000000 | 0% |
| Vela operator on Sentaurus states (SSS) | 9.444475139e-6 | 1.000202218 | +0.0202% |
| Self-consistent Vela curve (VVV) | 7.147306711e-6 | 0.756924225 | -24.3076% |

Thus the AverageBox coefficient plus qualified HFS mobility operator closes
the complete maximum-gm segment. The failed gm gate is a state-feedback
problem, not a remaining fixed-state coefficient or mobility amplitude error.

## State-family attribution

| Vg (V) | Full feedback (dex) | psi effect (dex) | phin effect (dex) | phip effect (dex) |
|---:|---:|---:|---:|---:|
| 1.000000 | -0.0890742 | +0.0013985 | -0.0904882 | ~0 |
| 1.166667 | -0.1065699 | +0.0013990 | -0.1079843 | 0 |

The electron-QF contribution becomes `0.0174961 dex` more negative across the
segment, while the electrostatic contribution stays constant to about
`5.3e-7 dex` and the hole-QF contribution remains zero. Consistently, the
centered P95 Vela-minus-Sentaurus `phin` difference grows from `2.068 mV` to
`5.417 mV`; the P95 edge `phin`-drop difference grows from `0.162 mV` to
`0.286 mV`. The corresponding centered P95 `psi` values are `4.269 mV` and
`4.531 mV`.

As an additional HFS check, the absolute-flux-weighted drain-cut electron
mobility is `0.0625667` and `0.0625568 m2/(V s)` for the two SSS endpoints, a
change of only `-0.0159%`. The same measure changes by `-0.0092%` on the VVV
states. The external geometric coefficient file is bias invariant; its
flux-weighted drain-cut `couple/length` changes only because the current
weights move. These observations agree with the complete SSS operator closure
and exclude a hidden endpoint HFS roll-off as the gm cause.

A second read-only A/B matrix holds the two converged Vela states fixed and
switches AverageBox, contact-HFS and all-HFS support independently. Removing
HFS changes the absolute current by about one order of magnitude, as expected,
but all five operator variants retain essentially the same endpoint growth:
`2.165795--2.165975`. Their log-current slope spans only
`2.1565e-4 dex/V`, which is `0.205%` of the qualified-to-Sentaurus slope gap
of `0.104974 dex/V`. AverageBox versus mesh-default support likewise changes
the frozen absolute gm but not its log slope. This separates operator amplitude
importance from the missing bias dependence: HFS is required, yet no tested
coefficient/mobility switch generates the observed slope deficit.

## Decision

1. The requested fixed-state coefficient/mobility replay **passes**: frozen
   operator gm error is `0.0202%`, far below the `20%` curve gate.
2. The `24.3076%` self-consistent gm deficit is assigned to electron-QF state
   feedback. It is not authorized as an engine floor and the known-difference
   ledger remains draft.
3. No IALMob, low-field mobility, flatband or threshold adjustment is
   authorized by this result.
4. The next minimal diagnostic is a two-endpoint electron-continuity/Jacobian
   reclose audit focused on the edges whose `phin` drop changes with bias. A
   full per-region Poisson AverageBox implementation remains deferred because
   the `psi` contribution is nearly bias independent on the failed segment.

## Reproducibility artifacts

- Replay implementation: `scripts/audit_templates_ldmos_g3_state_feedback.py`.
- Frozen Vela-state operator A/B:
  `scripts/audit_templates_ldmos_g3_gm_operator.py`; unit test:
  `tests/regression/test_audit_templates_ldmos_g3_gm_operator.py`.
- Vg=1.0 replay:
  `reference_staging/templates_ldmos_g3_averagebox_state_feedback_vg1_20260831/`.
- Vg=1.166667 replay:
  `reference_staging/templates_ldmos_g3_averagebox_state_feedback_vg1p166667_20260831/`.
- Operator A/B output:
  `reference_staging/templates_ldmos_g3_gm_operator_20260831/`.
- VM capture and imported endpoint:
  `reference_staging/templates_ldmos_g3_gm_state_capture_20260831/`.
- The large TDR, imported fields, hybrid states and edge probes remain ignored
  runtime evidence; the frozen metrics and hashes are recorded in this report
  and the AverageBox diagnostic contract.
