# Templates/LDMOS Stage-4 D5 Vg=4 step-growth A/B

Date: 2026-09-08 (Asia/Shanghai). Branch `codex/templates-ldmos-phase-a`,
source base `959c1c9` plus the uncommitted optional step-growth implementation.

Both independent Vg=4 V drain sweeps completed from the same Vd=0 equilibrium
checkpoint through 40 V. All 31 exact shared reference points pass the existing
single-gate engineering and final limits, including zero-bias KCL. All 31 paired
current and full-state equivalence checks pass. Neither arm requires rollback,
bisection or Gummel density recovery.

This experiment does **not establish a performance improvement**: Newton-aware
growth takes 6.29% more voltage transfers and 1.60% more applied Newton updates.
Its summed child wall time is 0.45% lower, while child CPU time is 2.31% higher.
This is one sequential paired experiment, without repeated timing samples.

This validates the **outer checkpoint adapter** of the growth policy. Each C++
child is a single-target solve with native `step_growth_mode="fixed"`; the
Python controller chooses the next physical voltage. It is not a full LDMOS
validation of the native C++ multi-point `newton_iterations` sweep path.

## Reproducible controls

The immutable [run plan](../../reference_staging/templates_ldmos_d5_vg4_step_growth_20260908/plan.json)
records the script, helper, runner, reference and seed hashes. The
[input inventory](../../reference_staging/templates_ldmos_d5_vg4_step_growth_20260908/input_inventory.json)
pins materials, exact mesh, doping and AverageBox transport couples.
The source diff and source-file hashes are saved beside the plan. The Release
runner SHA256 is `9638411ec0d38434f05afc765d573060adde0d406c6818f3428c76777a066391`.

Both arms load `templates_ldmos_shortstep_20260907/zero_closed_psi/state.csv`,
SHA256 `edc93d5ebe8c6b42eeb1bab215bafdf1669e937af70b8fcefe38923e6c1f541e`,
then independently reclose and advance from zero. Both initial closures take
zero applied Newton updates and produce zero terminal currents. No high-bias
states or previously scored prefix points are spliced into either new curve.
Preparation of this Vg=4 equilibrium seed, including gate ramping, is excluded
from the measured drain sweeps.

The physical model remains D5 without IALMob: Fermi-Dirac statistics,
OldSlotboom narrowing, SRH/Auger and the qualified constant-field HFS mobility,
with live repaired mobility derivatives and QF initial recentering. Predictor,
IALMob and avalanche remain off. The original exact mesh has 10,241 nodes,
30,022 edges and 19,782 triangles, with 16,237 external AverageBox couples
(3,548 zero couples), barycentric node volumes, material-local Poisson charge
and legacy node-local contact reconstruction. Geometry is in um, potentials
in V, state densities in m^-3 and reported terminal currents in A/um.
The build is UCRT64 Release, Eigen SparseLU/COLAMD, with unit scaling and
L2 row/column linear equilibration.

Both arms start at 2.5 mV, cap proposals at 100 mV and use a 160-update budget
per child. Every transfer uses a QF update cap equal to its actual drain
increment. A rejected direct attempt gets one strict same-bias reclose.
Density recovery is eligible only after all original global blocks pass and
local carrier rows still fail. Unqualified transfers roll back to the last
accepted parent and halve the actual step, with a 2.5 mV retry minimum.
Normal clipping to an exact reference point may be smaller than this retry
minimum; it does not change any acceptance threshold.

The fixed arm proposes `min(0.1, previous_proposal * 1.35)`, preserving the
legacy proposal across a clipped reference endpoint. The Newton-aware arm uses

```text
G(N) = 1 + 0.35 * max(0, 1 - max(N - 1, 0) / (0.75 * 160))
next_step = clamp(actual_accepted_step * G(N), 0.0025, 0.1)
```

Here N includes all reported applied updates in the logical transfer's direct
solve and reclose, plus eligible density recovery if used. Ordinary reclose
work is counted; it does not by itself suppress growth. Density recovery
would suppress growth. Children run sequentially, with arm order alternating
at each reference interval. Bulky per-node trace output is disabled equally
for both arms; state, terminal, attempt, iteration and profiler evidence remains.

## Complete sweep costs

| Full Vd=0 to 40 V cost | Fixed growth | Newton-aware growth |
| --- | ---: | ---: |
| Exact reference points | 31 | 31 |
| Accepted physical voltage transfers | 429 | 456 |
| Child executions, including initial closure | 858 | 908 |
| Reported applied Newton updates | 4,761 | 4,837 |
| Profiler Newton iteration bodies | 5,189 | 5,288 |
| Updates per accepted transfer, min / median / max | 5 / 10 / 57 | 3 / 10 / 37 |
| Directly accepted transfers | 1 | 5 |
| Transfers accepted after strict reclose | 428 | 451 |
| Rollbacks / density recoveries | 0 / 0 | 0 / 0 |
| Sum of child wall time | 5,944.764 s / 99.079 min | 5,917.893 s / 98.632 min |
| Sum of child CPU user + kernel time | 4,320.281 s / 72.005 min | 4,420.219 s / 73.670 min |

All direct attempts, including their rejected work, and all recloses are included
in time totals. A final rejected Newton iteration can perform assembly and
linear work without applying an update; hence profiler iteration bodies and
reported applied updates differ. Neither count is the number of output points.

The combined controller ran from 09:48:26 to 13:10:51 +08:00, 12,145.833 s
(202.431 min) including both arms, bookkeeping, comparisons and scoring.
It exited with code 0. Per-arm summed child time excludes controller overhead,
the common seed's historical preparation and final postprocessing. This report
does not infer a SDevice runtime ratio from this paired Vela experiment.

The wall-clock improvement is only 26.871 s across about 99 minutes per arm,
and changes sign at intermediate matched endpoints. CPU work increases by
99.938 s. The results therefore provide no clear speedup evidence.

## Numerical qualification

The original accepted-state limits are unchanged: psi <=5e-8, electron <=1e-11,
hole <=3e-10, zero local carrier-row violations with maximum ratio <=1e-8,
and terminal KCL/max-terminal <=1e-8. All accepted transfers satisfy them.

| Maximum over all accepted nonzero transfers | Fixed | Newton-aware | Limit |
| --- | ---: | ---: | ---: |
| Psi block | 8.09951e-9 | 6.58947e-9 | 5e-8 |
| Electron block | 9.29197e-12 | 9.61479e-12 | 1e-11 |
| Hole block | 2.56595e-19 | 2.64711e-19 | 3e-10 |
| Qualified local carrier-row ratio | 8.87265e-9 | 9.75531e-9 | 1e-8 |
| Normalized terminal KCL | 2.11451e-11 | 1.23894e-11 | 1e-8 |

All saved accepted states contain 10,241 unique finite node records and
nonnegative densities. Across the 31 paired exact points, maximum normalized
terminal-current difference is 2.28322e-12 (limit 1e-8). Maximum psi/phin/phip
differences are 4.97380e-14 / 1.63347e-10 / 8.33023e-11 V (limit 1e-7 V).
Maximum relative electron/hole density differences, with the existing
1 m^-3 denominator floor, are 6.31866e-9 / 3.22228e-9 (limit 1e-6).

The unchanged scorer uses the frozen D5-no-IALMob SDevice Vg=4 curve. Actual
solved biases remain in provenance; the existing <=32 binary64 ULP endpoint
matching convention is the only alignment operation. No interpolation is used.

| Existing single-gate score | Fixed | Newton-aware | Final limit |
| --- | ---: | ---: | ---: |
| Median absolute current error | 2.13030% | 2.13030% | 5% |
| P95 absolute current error | 9.31228% | 9.31228% | 12% |
| Low-voltage resistance error | 4.14501% | 4.14501% | 10% |
| 40 V current error | 2.03117% | 2.03117% | 10% |
| Maximum normalized KCL over 31 points | 1.75295e-10% | 1.91128e-10% | 0.1% |
| Engineering / final Vg=4 verdict | Pass / Pass | Pass / Pass | |

Zero-bias terminal currents and KCL are exactly zero in both new runs. Relative
current statistics use the 30 resolved nonzero points; KCL includes all 31.
At 40 V both Vela currents round to 2.29399922075054e-4 A/um, versus
2.24833171719662e-4 A/um in the reference. The score reproduces the previously
qualified Vg=4 curve, now using two independent complete drain sweeps.
Vg=8 and two-gate ratio qualification are outside this rerun; the outstanding
findings in the [two-gate report](templates_ldmos_stage4_two_gate_validation_2026-09-07.md)
are not cleared by this result.

## Interpretation and evidence

Both arms hit a clipped reference endpoint 30 times. The Newton-aware policy
then regrows from that smaller accepted increment. It has 338 transfers at
100 mV versus 386 for fixed growth, and 118 smaller transfers versus 43.
Its smallest accepted endpoint increment is 0.942299 mV, with no failed-step
subdivision. This records a concrete source of extra work in the current
adapter. The comparison includes both iteration-aware growth and the switch
from legacy proposal memory to the actual clipped increment; it does not
isolate the growth formula alone.

At the 100 mV ceiling, both successful-step policies retain that ceiling.
Neither can exploit larger voltage steps in this experiment. Strict reclose
also remains the usual acceptance path, costing 1,137.470 / 1,178.549 s of
child wall time. Changing growth alone does not remove this work.
The profiler records approximately 1,412 / 1,414 s in Jacobian assembly and
1,801 / 1,768 s in linear factorization. These are inclusive wall-time scopes;
nested profiler stages must not be added as independent costs.

A useful next isolated experiment is to preserve the free-running proposal
across mandatory output-point clipping, while retaining Newton-dependent growth
and all original gates. Larger max-step policies or integrating strict reclose
inside the native solver need separate validation. No physical defaults or
acceptance limits were changed by this rerun.

- [Controller and run plan](../../reference_staging/templates_ldmos_d5_vg4_step_growth_20260908/run_ab.py)
- [Full audit](../../reference_staging/templates_ldmos_d5_vg4_step_growth_20260908/audit_summary.json), covering 1,766 children and 885 accepted transfers
- [Numerical and profiler summary](../../reference_staging/templates_ldmos_d5_vg4_step_growth_20260908/result_summary.json)
- [Matched interval costs](../../reference_staging/templates_ldmos_d5_vg4_step_growth_20260908/interval_costs.csv) and [accepted-step data](../../reference_staging/templates_ldmos_d5_vg4_step_growth_20260908/accepted_transfers.csv)
- [Fixed score](../../reference_staging/templates_ldmos_d5_vg4_step_growth_20260908/fixed/score/summary.json) and [Newton-aware score](../../reference_staging/templates_ldmos_d5_vg4_step_growth_20260908/newton_iterations/score/summary.json)
- [Figure PDF](../../reference_staging/templates_ldmos_d5_vg4_step_growth_20260908/vg4_growth_comparison.pdf)

![D5 Vg=4 complete paired curves and costs](../../reference_staging/templates_ldmos_d5_vg4_step_growth_20260908/vg4_growth_comparison.png)
