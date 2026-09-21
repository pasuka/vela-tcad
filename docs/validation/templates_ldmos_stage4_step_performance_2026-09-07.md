# Templates/LDMOS Stage-4 voltage-step performance controls

Date: 2026-09-07 (Asia/Shanghai). Experiments completed on 2026-09-06;
artifact verification and reporting completed on 2026-09-07.

The requested local performance validation is complete. At Vg=4 V,
11.9 to 12 V, a 100 mV attempt followed by a strict same-bias solve with
the existing QF recenter option takes 10 Newton updates, versus 493 for
the 2.5 mV reference path. Both repeated trials pass the original residual
ceilings and the additional endpoint current, KCL, potential and density
checks. Combined child CPU time is 11.36–11.67 s versus 378.81 s locally.

This does not qualify a global 100 mV policy. Direct 100 mV stepping fails
the electron residual ceiling. The recovery workflow passes the electrical
checks at 4.1 and 8.1 V too, but the 4.1 V electron density difference is
1.13385e-6, above the predeclared 1e-6 consistency limit. No thresholds or
production solver policies were changed. The original D5 sweep remains stopped.

## Reproduction context and preserved production state

- Worktree: `D:\code-repo\vela-tcad\.worktrees\templates-ldmos-phase-a`.
  Branch `codex/templates-ldmos-phase-a`, HEAD
  `345fda30d541003f4aa3a08341dca59418f06ded`. The earlier Jacobian repairs
  are uncommitted worktree changes; HEAD alone does not identify this binary.
- Pinned executable: `build-release/vela_example_runner_step_performance_20260906.exe`.
  SHA256 `89ed3a95e828e4190adf460dfa0881fd90a8749bef203bb4ed5adc02ffa757fa`.
  Windows MSYS2 UCRT64 Release, `-O3 -DNDEBUG`; Eigen SparseLU/COLAMD selected
  by `VELA_LINEAR_SOLVER=sparselu`, with L2 row/column equilibration.
- Exact imported mesh: 10241 nodes, 30022 edges, 19782 triangles;
  16237 external AverageBox transport records, including 3548 zero couples.
  Barycentric node volumes, material-local Poisson charge, and
  `legacy_node_local` contact reconstruction.
- 300 K, Fermi statistics, OldSlotboom, SRH/Auger, `constant_field` HFS,
  edge-projected QF driving force and existing contact electric-field fallback.
  Live mobility Jacobian and `quasi_fermi_recenter_on_initial_state=true`.
  Predictor and IALMob remain disabled; no avalanche, heating or quantum model.
- `unit_scaling` input mode; geometry in um, terminal comparisons use the
  explicit `current_total_A_per_um` CSV column. State comparisons use physical
  `psi`, `phin`, `phip` in V and `electrons_m3`, `holes_m3` in m^-3.
  Residual norms below are solver block norms, not terminal currents.
- Original enforced ceilings: psi <=5e-8, electron <=1e-11, hole <=3e-10.
  `max_iter=160`, `block_filter` line search and 0.1 V QF update limit retained.

The production runner was stopped at the user's request on 2026-09-06
21:41:58 +08:00. Its final complete accepted checkpoint is 11.899167 V,
SHA256 `4600e6a9a34dcf00eebbc6fd8dca10fe515dc2a7d4dc1ddd791c2c05bfd2bf37`.
The stop and checkpoint are preserved in the
[production evidence directory](../../reference_staging/templates_ldmos_jacobian_repair_20260906/vg4_continue_4p01_to_40/).
The stop-related process return code is an administrative termination, not
a new numerical convergence result. No Vela runner was active at the final check.

## Controls and acceptance criteria

All main controls load one accepted 11.9 V seed, prepared from the preserved
11.899167 V checkpoint. Seed SHA256:
`bea5b12f69fdbe551adfe09602a14cc99c8ee5292924b311f047ff57d70ebbe6`.
Preparation costs 20 updates, 46.1924 s wall and 14.453125 s child CPU time;
it is excluded from every 11.9-to-12 V interval timing below.

Each main control performs the same initial 11.9 V reclose (one update),
then follows a fixed maximum step to 12 V. The startup accepted-state hash
is identical (`dc3291a5e02a3674`). Minimum, initial and maximum steps are equal,
growth is 1, retries are 0, and failure stops the run. Only the last step may
clip to the endpoint. Output requests are `[11.9, 12]` in every main case,
with the same diagnostics and accepted-step checkpoint settings.

The [predeclared plan](../../reference_staging/templates_ldmos_step_performance_20260906/plan.json)
adds these consistency checks without changing the original solver ceilings:

| Endpoint comparison | Limit |
| --- | ---: |
| Maximum port-current difference / maximum reference terminal magnitude | 1e-8 |
| Absolute sum of terminal currents / maximum terminal magnitude (KCL) | 1e-8 |
| Maximum absolute difference of psi, phin, phip | 1e-7 V |
| Maximum carrier density relative difference, denominator max(abs(reference), 1 m^-3) | 1e-6 |

Bias matching allows only 32 binary64 ULPs for endpoint bookkeeping and
retains the raw CSV bias; no interpolation is used. For example, the fine
path ends at 11.99999999999998 V. Density comparisons cover all 10241 nodes.

## Fixed-step results at 11.9 to 12 V

Updates include the common startup update. Rejected line-search evaluations
are timed but are not counted as applied Newton updates. Failure rows are
incomplete intervals and cannot be assigned a speedup.

| Fixed step (mV) | Accepted transfers | Updates | Wall (s) | Child CPU (s) | Result / final electron norm |
| ---: | ---: | ---: | ---: | ---: | --- |
| 2.5 | 40 | 493 | 666.04 | 378.81 | Pass; 6.46409e-12 |
| 3 | 34 | 425 | 375.58 | 258.16 | Pass; 5.28247e-12 |
| 3.75 | 0 | 20 | 20.60 | 13.78 | Reject at 11.90375 V; 1.21675e-11 |
| 5 | 0 | 133 | 167.89 | 110.59 | Reject at 11.905 V; 1.67592e4 |
| 10 | 0 | 38 | 60.10 | 35.44 | Reject at 11.91 V; 2.89904e4 |
| 20 | 0 | 90 | 117.87 | 77.12 | Reject at 11.92 V; 4.58266e4 |
| 50 | 0 | 8 | 13.63 | 8.89 | Reject at 11.95 V; 9.66323e-11 |
| 100 | 0 | 9 | 13.87 | 9.83 | Reject at 12 V; 2.27093e-10 |
| 100 repeat | 0 | 9 | 14.53 | 8.97 | Same rejection and norm |

All rejected cases report `block_absolute_convergence_line_search_rejected`.
Behavior is not monotonic in step size: 5–20 mV enter poor nonlinear
trajectories, while 50–100 mV approach a much smaller residual before rejection.
Thus neither a universal maximum step nor a universal roundoff explanation
can be inferred from these failures.

The 3 mV path passes every endpoint check: normalized current difference
2.67732e-13, KCL 1.85824e-14, maximum potential difference 2.48690e-14 V,
maximum density relative difference 1.58141e-13. It uses 13.8% fewer updates
than the baseline (493/425=1.16), with observed CPU speedup 1.47x.
Only one 3 mV run and one fine baseline were timed.

## Strict same-bias recovery of the 100 mV trial

The rejected 12 V final state is used only as an input checkpoint to another
solve at 12 V. The existing initial-state QF recenter option is enabled, as
it is in all other controls. The state is accepted only after that new solve
passes every original gate. This is an external experimental workflow;
automatic recovery inside the production sweep has not been implemented.

Each trial uses 1 startup + 8 direct-step + 1 recovery update = 10 updates.
The following totals include the unsuccessful direct attempt, its output,
and the subsequent recovery process:

| Trial | Total wall (s) | Total child CPU (s) | CPU speedup vs fine path | Update ratio |
| --- | ---: | ---: | ---: | ---: |
| First | 16.1430 | 11.671875 | 32.46x | 49.3x |
| Repeat | 18.2358 | 11.359375 | 33.35x | 49.3x |

Both direct iteration CSVs, rejected final states, recovered iteration CSVs,
and recovered final states are byte-identical between repeats. At 12 V:

| Quantity | Recovered value |
| --- | ---: |
| psi / electron / hole residual | 1.30087e-9 / 1.93385e-12 / 2.12167e-21 |
| Drain current | 2.244483331269852e-4 A/um |
| Fine-reference drain current | 2.244483331271577e-4 A/um |
| Maximum normalized port-current difference | 7.68536e-13 |
| KCL | 5.23372e-13 |
| Maximum physical potential difference | 1.03588e-10 V |
| Maximum carrier density relative difference | 4.00743e-9 |

All endpoint criteria pass. The direct 100 mV attempt had rejected a line
search after 13 trials; its last global raw update norm alone cannot prove
that every component was below floating-point resolution. The successful
recentered restart establishes a useful local remedy, while the prior
[Jacobian/precision investigation](templates_ldmos_stage4_jacobian_repair_2026-09-06.md)
provides separate evidence about update precision.

Timing is from sequential runs with other host workloads present. Child
CPU time is Windows `GetProcessTimes` user+kernel time; wall time includes
process startup and output. Scheduling, CPU frequency and cache effects
remain uncontrolled, and the fine baseline was not repeated. With common
seed preparation included in both paths, the CPU ratio is 15.05–15.24x,
not 32–33x. Neither ratio predicts a complete 0-to-40 V run or SDevice parity.

## Cross-bias checks and the remaining limitation

The [cross-bias plan](../../reference_staging/templates_ldmos_step_performance_20260906/cross_plan.json)
uses archived accepted 2.5 mV-path endpoints at 4.1 and 8.1 V, each reclosed
at the same bias under the original gates. These are accuracy checks, not
new paired fine-path timing benchmarks. Each direct 100 mV attempt fails;
one recovery update then passes the original residual ceilings.

| Interval (V) | Recovered electron norm | Normalized current difference | KCL | Maximum density relative difference | All extra checks |
| --- | ---: | ---: | ---: | ---: | --- |
| 4 to 4.1 | 4.41195e-13 | 6.41584e-14 | 9.38077e-14 | 1.13385e-6 | Fail density |
| 8 to 8.1 | 9.28590e-13 | 1.56850e-13 | 1.65001e-15 | 2.94314e-10 | Pass |
| 11.9 to 12 | 1.93385e-12 | 7.68536e-13 | 5.23372e-13 | 4.00743e-9 | Pass twice |

At 4.1 V, maximum phin difference is 2.93123e-8 V, within the 1e-7 V
potential limit. The maximum relative electron-density mismatch occurs at
node 376: candidate 1.0847695672230946e10 versus reference
1.0847683372590759e10 m^-3. Its small absolute difference in a depleted
region does not waive the predeclared relative-density criterion.

An additional strict reclose of the already accepted 4.1 V candidate
performs zero updates, returning `initial_block_abstol`; the initial
residual already passes the solver's existing acceptance test. The density
and physical-potential discrepancies therefore remain. This agrees with
the initial-convergence branch in [NewtonSolver.cpp](../../src/solver/NewtonSolver.cpp).
Original electrical acceptance and the stronger endpoint density comparison
are distinct checks; both outcomes are retained.

Before promoting automatic large-step recovery, the next investigation
should resolve the 4.1 V state discrepancy under unchanged production
ceilings, then validate bounded intervals across bias and full D5 reference
points. Current evidence does not justify changing the global step policy.
Vg=8 V, the two-gate current ratio, and IALMob remain behind the full D5 gate.

## Cost attribution, evidence and verification

The fine-path profiler reports 304.22 s in `linear.factorize`, 194.82 s in
`newton.jacobian`, 6.46 s in `linear.solve`, 59.10 s in `dd.residual`, and
58.25 s in `newton.line_search`. These wall-time scopes overlap and must
not be added. `dc.record_point` is 2.18 s but does not isolate all checkpoint
output. Avoid claiming that all I/O is negligible. Reducing repeated Newton
assembly and factorization is a supported performance direction; these
measurements alone do not establish a linear-solver accuracy defect.

The [experiment directory](../../reference_staging/templates_ldmos_step_performance_20260906/)
contains control decks, process status/timing, raw iteration/attempt CSVs,
terminal currents, state CSVs, profiler output, and the execution scripts.
The [comparison](../../reference_staging/templates_ldmos_step_performance_20260906/comparison.json)
includes failed-attempt cost in each recovered total. The
[final audit](../../reference_staging/templates_ldmos_step_performance_20260906/final_audit.json)
checks all 19 completed cases and 19 unique complete input/accepted state
files, configuration/checkpoint hashes, no predictor/retries, fixed-step
bounds, original residual gates, exact repeat identity and cross-bias
comparisons. It records fixture and evidence SHA256 values and preserves
the expected 4.1 V negative result.

Read-only numerical reanalysis from the worktree root:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python reference_staging/templates_ldmos_step_performance_20260906/analyze_controls.py
python reference_staging/templates_ldmos_step_performance_20260906/finalize_evidence.py
```

Both analyses pass their integrity checks. Helper scripts also pass Python
syntax compilation; report links and `git diff --check` were checked.
This task adds the report and isolated experiment/audit artifacts only;
no solver source or production template was changed, and no solver tests
were rerun solely for the report. Generated outputs remain in ignored staging.
