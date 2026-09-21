# Templates/LDMOS Stage-4 recovery validation

Date: 2026-09-07 (Asia/Shanghai).

The 4.1 V density discrepancy is resolved. The previously recovered solution
met the absolute block ceilings, but its depleted carrier rows were not
sufficiently balanced for the separate density consistency requirement.
Enforcing the existing local carrier-row criterion exposed an inactive
`block_filter` path. A narrowly scoped filter correction now permits local
progress while preserving every original block ceiling. One additional Newton
update reduces the maximum electron density difference from 1.13385e-6 to
6.04002e-13; current, KCL and potential consistency also pass.

The proposed full-range recovery policy is **not qualified**. Its isolated
Vg=4 V continuation accepts thirteen 100 mV transfers from 4 to 5.3 V, then
fails to reach the next exact reference, 5.33333333333333 V. Recenter/reclose,
step bisection, uniform Newton scaling and a balanced pair of substeps do not
resolve that interval. The validation stopped with a numerical failure;
no complete 31-point curve or full-range speedup is claimed. The original
production sweep remains stopped at the user's request.

## Configuration and provenance

- Worktree: `D:\code-repo\vela-tcad\.worktrees\templates-ldmos-phase-a`;
  branch `codex/templates-ldmos-phase-a`; HEAD
  `345fda30d541003f4aa3a08341dca59418f06ded`. Existing uncommitted Jacobian,
  diagnostics and scoring changes predate this investigation and are retained.
- MSYS2 UCRT64, CMake `windows-ucrt64-release`, `-O3 -DNDEBUG`;
  actual backend Eigen SparseLU/COLAMD (`VELA_LINEAR_SOLVER=sparselu`),
  L2 row/column equilibration. The experiments did not use a Debug runner.
- Pinned corrected runner:
  `build-release/vela_example_runner_recovery_20260907.exe`, SHA256
  `0ffd9ad4731d7f0b130f575aed267722ad80db654787b0253a702214e6e4db7f`.
  The final rebuilt runner has SHA256
  `37e2f099b08864a3d2be8aef2e4047ba479b3f6ef22a11967e08cfb9f6a733b4`.
  Repeating the 4.1 V correction with that final build gives identical Newton
  iteration CSV bytes and zero physical-state difference from the pinned run.
- Exact imported mesh: 10241 nodes, 30022 edges, 19782 triangles; 16237 external
  AverageBox records, including 3548 zero couples. Barycentric volumes,
  material-local Poisson charge, `legacy_node_local` contact reconstruction.
- 300 K; Fermi statistics, OldSlotboom, SRH/Auger; `constant_field` HFS,
  edge-projected QF driving force and contact electric-field fallback; live
  repaired mobility Jacobian. Predictor and IALMob disabled. No avalanche,
  heating or quantum model was introduced.
- `unit_scaling`; geometry in um, potentials in V, state densities in m^-3,
  terminal comparisons use explicit A/um columns. Residual values are solver
  block norms, not terminal currents.
- Original enforced block ceilings: psi 5e-8, electron 1e-11, hole 3e-10.
  Newton maximum 160 updates, QF update cap 0.1 V, initial-state QF recenter
  enabled. Neither ceilings nor endpoint comparison tolerances were relaxed.
- Strict recovery controls additionally enforce the existing qualified
  carrier-row residual ratio <=1e-8, scale floor 1e-30, source qualifications
  set to zero. This is an experimental stronger local convergence contract;
  no template or global default was changed.

The original accepted 11.899167 V production checkpoint still has SHA256
`4600e6a9a34dcf00eebbc6fd8dca10fe515dc2a7d4dc1ddd791c2c05bfd2bf37`.
Its administrative stop is preserved in the
[production directory](../../reference_staging/templates_ldmos_jacobian_repair_20260906/vg4_continue_4p01_to_40/).
The earlier [performance report](templates_ldmos_stage4_step_performance_2026-09-07.md)
remains unchanged.

## Why the 4.1 V endpoint differed

At node 376, the recovered state has electron residual 5.88001e-19 and
local scale 1.22229e-12, giving relative local imbalance 4.81067e-7.
The global electron norm is only 5.93232e-12, below its 1e-11 ceiling.
Thus a globally accepted solution can still differ appreciably in relative
density in a depleted region: twelve nodes exceed the 1e-6 density comparison
limit, and 207 qualified carrier rows fail the extra 1e-8 local criterion.

The electron row 376 directional derivative agrees with a central difference:
at a 1e-6 V phin perturbation, its relative derivative error is 5.31055e-10.
The probe scans 1e-4 through 1e-9 V. This supports the row's existing mobility
derivative at this state; it is not a claim that every Jacobian block is exact.
Use the row-relative errors in `jvp_rows.csv`; the aggregate probe's differently
normalized `relative_error` column is not a relative error for this weak row.

When the local criterion is enforced, the previous block filter has no active
global block to improve and rejects the useful correction. The raw proposed
phin update at node 376 is +2.93122e-8 V, matching the discrepancy from the
fine-path state. The fix in [NewtonSolver.cpp](../../src/solver/NewtonSolver.cpp)
allows a candidate only when both convergence modes are enforced, the current
global blocks already pass, the candidate preserves all block ceilings, and
the maximum qualified local ratio satisfies the filter's sufficient decrease.
Final acceptance still requires the local criterion. Off/report behavior is
preserved; a row-off 100 mV control reproduces the old Newton trace byte for byte.

After one update, the corrected candidate has block norms
1.13256e-9 / 4.07628e-13 / 2.62046e-24, zero local violations, and maximum
local ratio 3.68329e-13. Against an independently reclosed fine-path reference:

| Quantity | Measured | Preserved comparison limit |
| --- | ---: | ---: |
| Maximum electron density relative difference | 6.04002e-13 | 1e-6 |
| Maximum phin absolute difference | 1.59872e-14 V | 1e-7 V |
| Maximum port-current difference / reference maximum terminal magnitude | 1.44501e-16 | 1e-8 |
| Absolute terminal sum / maximum terminal magnitude | 7.86351e-14 | 1e-8 |

All potentials and both densities pass across all 10241 nodes. Comparison
against the archived, unreclosed fine reference also passes.

## Recovery strategy controls and negative results

At 4.3 V, a rejected 0.2 V transfer can satisfy global ceilings while retaining
local imbalance. Recenter alone is insufficient. A diagnostic undamped update
reduces violating rows from 281 to 39 but increases the worst relative row.
Its raw linear solve has normwise backward error 1.07413e-21, yet maximum
componentwise backward error 0.266901 in row 23975, whose denominator is
1.42029e-20. A small global linear error does not establish accuracy of every
weak row. This is evidence for further numerical investigation, not proof
that the linear backend explains all failures.

The existing Gummel density recovery, followed by strict Newton closure,
does recover this 4.3 V state: two density passes in one cycle, three reported
Newton updates, zero local violations, maximum local ratio 2.84140e-10.
Electron density difference from the reclosed fine reference is 3.82457e-10.
An experimental L2 local-merit variant did not help and was reverted before
the final build. Its isolated executable and negative evidence are retained.

Applying density recovery indiscriminately to large-step states with grossly
failed global blocks instead raises a charge-neutral-potential bracketing
exception. Consequently the full strategy permits density recovery only after
all original global ceilings already pass; it does not accept those bad states.

The full strategy uses the actual 31-point Sentaurus CSV, SHA256
`3c9266b2616e2e779b853cee7bd94a0ee13478015d727c5bc2306ed8e152cb6a`.
There are 27 remaining reference points above the accepted 4 V seed.
Internal 100 mV transfers clip at each reference; failed transfers roll back
and bisect, with minimum retry step 2.5 mV. Each attempt first uses Newton,
then one same-bias reclose, then conditionally the density recovery above.
Every accepted transfer must meet original blocks, local rows and KCL <=1e-8.

The thirteen accepted internal transfers end at 5.299999999999995 V. Each passes
these electrical and local checks, but they do not establish full-state
equivalence against an independent fine path at every intermediate bias.
The next reference and all rollback probes fail:

| Parent / target V | Direct Newton updates | Reclose updates | Result |
| --- | ---: | ---: | --- |
| 5.3 / 5.33333333333333 | 160 | 160 | Rejected |
| 5.3 / 5.316666666666663 | 75 | 40 | Rejected |
| 5.3 / 5.308333333333329 | 139 | 75 | Rejected |
| 5.3 / 5.304166666666662 | 160 | 27 | Rejected; next half-step below retry minimum |

The run finished at 09:57:57 +08:00 with 35 child executions, 996 reported
Newton updates, 868.203125 s child CPU time and 1173.261854 s summed child
wall time. These totals include failed attempts, startup and recovery, and
exclude independent diagnostic controls. This incomplete interval has no
meaningful full-curve speedup. No new exact reference above 4 V was accepted.

Two additional controls also fail: uniform QF trust-region scaling for
5.3 to 5.33333333333333 V reaches 160 updates; splitting 5.2 to 5.33333333333333 V
into two equal substeps fails on the first step (61 updates plus 34 reclose
updates). Their final global blocks fail too, so these are not merely local
ratio threshold failures. Failure labels alone must not be read as proof that
global convergence was achieved.

An earlier draft queue accidentally included the old 4.01 V startup as a
reference point. It was stopped and replaced by the CSV-derived queue above;
its output is excluded from the full-strategy result. Other superseded pilots
remain separately recorded as stopped or negative experiments.

## Zero bias, tests and remaining work

The separate zero-bias KCL problem persists. Reclosing the old zero-bias state
gives raw KCL 130.402%; explicitly zeroing its initial physical and split QF
coordinates before re-solving gives 127.824%. Terminal magnitudes are around
1e-25 A/um, but the current scoring contract has no zero-bias exemption.
Neither case qualifies; no scoring denominator, limit or current was altered.

Final Release validation passed: NewtonSolver 112 cases / 1484 assertions;
DCSweep 98 cases / 3479 assertions; Templates/LDMOS Python regression 27 tests;
registered `pn2d_config_templates` CTest 1 selected / 1 passed. The new Newton
test checks correction of a local imbalance inside already satisfied global
ceilings and preservation of reporting-only behavior. The schema documents
the new conditional filter path. Final-build 4.1 V output matches the pinned
corrected runner. These passing tests do not override the full-strategy failure.

All experiment processes in this worktree had exited at the final inspection.
No other worktree's processes were stopped. Generated controls, checkpoints,
logs and audit output are kept in the ignored
[recovery evidence directory](../../reference_staging/templates_ldmos_recovery_20260907/).
The [audit summary](../../reference_staging/templates_ldmos_recovery_20260907/audit_summary.json)
checks input hashes, original gates, accepted state hashes, KCL, comparisons,
test logs and the preserved production checkpoint. Reproduction helpers there
refuse to overwrite existing case directories.

The next numerical blocker is the short-step path from the accepted 5.3 V
state, including its first Newton updates, componentwise limiting/Poisson
recorrection and rollback state selection. It requires a controlled comparison
from that same checkpoint before a new global recovery policy can be adopted.
Full Vg=4 scoring, zero-bias KCL qualification, Vg=8 and the drain-current
ratio between gate biases remain outstanding. IALMob remains gated on D5
qualification.
