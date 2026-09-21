# Templates/LDMOS Stage-4 D5 Vg=8 full drain-sweep rerun

Date: 2026-09-08 (Asia/Shanghai). Branch `codex/templates-ldmos-phase-a`,
source base `959c1c9` plus the unchanged optional step-growth implementation.

The fresh Vg=8 V sweep completed all 31 exact Vd=0-to-40 V reference points.
Every accepted state passes the original global residual, local carrier-row
and terminal KCL criteria. The complete integrity audit passes. There are no
rollbacks; two early transfers require the existing eligible density recovery.
The planned, verified coordinate change at physical 26.666667 V allows the
run to pass the historical near-28 V precision obstruction without interruption.

The **D5 accuracy verdict remains engineering fail and final fail**. Vg=8
current-error median is 17.28090%, P95 18.64330%, and endpoint error 16.20605%.
Combining this new curve with the fresh fixed-growth Vg=4 curve gives a 40 V
two-gate ratio error of 13.89269%, above the unchanged 8% final limit.
Predictor and IALMob remain off.

## Configuration and full-run provenance

The run uses the same UCRT64 Release runner as the
[Vg=4 A/B experiment](templates_ldmos_stage4_vg4_step_growth_ab_2026-09-08.md),
SHA256 `9638411ec0d38434f05afc765d573060adde0d406c6818f3428c76777a066391`.
The runtime backend is Eigen SparseLU/COLAMD, explicitly selected by
`VELA_LINEAR_SOLVER=sparselu`, with unit scaling and L2 row/column equilibration.
Source and preset hashes were verified against that Vg=4 run.

The exact imported mesh has 10,241 nodes, 30,022 edges and 19,782 triangles,
with 16,237 external AverageBox carrier couples, barycentric node volumes,
material-local Poisson charge and legacy node-local contact reconstruction.
The D5 physics remains 300 K Fermi/OldSlotboom, SRH/Auger and constant-field
GradQF HFS with live repaired mobility derivatives and QF initial recentering.
Predictor, IALMob, quantum, avalanche and heating are off. Geometry is in um,
potentials in V, state densities in m^-3 and terminal currents in A/um.

The initial state is the previously closed Vg=8, Vd=0 equilibrium checkpoint,
SHA256 `d1e48ccd2e089bbe4ef3fa8d3d6a25c904b003d104012c76d848484efc906ccf`.
The new initial coupled solve independently accepts it with zero applied
updates and exactly zero terminal currents. Each subsequent scored point
comes from this new accepted-state chain. Historical gate-ramp/seed preparation
is outside the measured drain sweep.

The outer checkpoint driver uses fixed growth: initial step 2.5 mV, multiplier
1.35, maximum 100 mV, with clipping at every exact reference point. Each C++
child is a single-target Newton solve with a 160-update budget. Its physical
QF update cap equals the actual voltage increment. A rejected direct attempt
gets one strict same-bias reclose. Gummel density recovery is eligible only
when all original global block ceilings pass and local carrier rows still
fail. Unqualified transfers would roll back and halve the step, with the
declared 2.5 mV retry minimum. No tolerance, physical model or solver source
was changed during this rerun.

The [immutable run plan](../../reference_staging/templates_ldmos_d5_vg8_rerun_20260908/plan.json)
pins the runner, seed, material/mesh/doping/couple inputs, controller/helper
sources, frozen SDevice references and the fresh Vg=4 scoring inputs.

## Completion and measured cost

The controller ran from 13:43:43.682 to 16:10:21.907 +08:00 and exited with
code 0. All child executions are complete; no active child remains in its ledger.

| Full new drain sweep | Measured |
| --- | ---: |
| Exact reference points | 31 |
| Accepted physical voltage transfers | 429 |
| Child executions | 861 |
| Reported applied Newton updates | 4,985 |
| Updates per physical transfer: min / median / max | 5 / 11 / 135 |
| Rollbacks | 0 |
| Transfers requiring eligible density recovery | 2 |
| Sum of child wall time | 8,569.160 s / 142.819 min |
| Sum of child CPU user + kernel time | 5,278.234 s / 87.971 min |
| Controller elapsed time | 8,798.224 s / 146.637 min |

The 861 children comprise 429 direct attempts, 428 strict recloses, two density
recovery children, one initial closure and one frame-equivalence closure.
All new child costs, including rejected attempts and recovery, are included.
Controller time additionally includes bookkeeping, comparisons and scoring;
the child totals exclude these outer costs and historical seed preparation.
This is an observed full-run cost, not a paired performance comparison against
the old Vg=8 run or Sentaurus.

Density recovery is used at physical 0.84621190996608 V and
1.1462119099660801 V. Their full transfer costs are 18 and 135 applied updates.
Both qualify under the original gates after recovery. The final physical 40 V
raw state SHA256 is `358db151b8046821bd61114ac2574799162a6f45629205c9e3dfcdab1af2a7de`.

## Verified coordinate change

At the accepted physical 26.6666666666667 V checkpoint, the controller subtracts
28 V from physical potentials, QF references and all contact biases, preserving
the existing conventions for inactive oxide carrier placeholders and increments.
One coupled same-bias solve qualifies the translated state in one applied update.
After inverse translation for comparison, the maxima are:

| Equivalence quantity | Observed | Preserved limit |
| --- | ---: | ---: |
| Normalized terminal-current difference | 1.34338e-13 | 1e-8 |
| Psi / electron QF / hole QF absolute difference | each 1.42109e-14 V | 1e-7 V |
| Electron density relative difference | 3.90673e-13 | 1e-6 |
| Hole density relative difference | 2.06336e-13 | 1e-6 |
| KCL ratio before / after | 1.51370e-13 / 3.23583e-14 | 1e-8 |

Only after those checks pass does the new raw state become the parent of the
remaining sweep. Raw source/substrate are -28 V and gate is -20 V; raw drain
then advances to 12 V. Physical Vg stays 8 V and physical Vd is raw drain minus
raw source. All raw states and solver biases retain this frame in provenance.
The physical curve uses directly computed terminal currents. It crosses 28 V
and reaches 40 V with no retry subdivision.

## Integrity and numerical score

The audit checks all 861 child configurations, attempt records, parent hashes
and frame/voltage lineage, plus 431 distinct accepted state files. Every state
has 10,241 unique finite node records and nonnegative carrier densities.

| Maximum over accepted states | Observed | Original limit |
| --- | ---: | ---: |
| Psi block residual | 4.17082e-8 | 5e-8 |
| Electron block residual | 9.86674e-12 | 1e-11 |
| Hole block residual | 7.64989e-20 | 3e-10 |
| Qualified local carrier-row ratio | 8.17210e-9 | 1e-8 |
| Terminal KCL/max-terminal ratio | 6.27902e-12 | 1e-8 |

Scoring preserves the existing <=32 binary64 ULP endpoint matching convention
after exact raw-to-physical coordinate conversion. There is no interpolation.
Current statistics use 30 resolved nonzero points; KCL includes all 31 points.
The zero-bias currents and KCL are exactly zero, without an exemption.

| Vg=8 score | Observed | Engineering limit | Final limit |
| --- | ---: | ---: | ---: |
| Median absolute current error | 17.28090% | 15% | 5% |
| P95 absolute current error | 18.64330% | 25% | 12% |
| Low-voltage resistance error | 3.98476% | 20% | 10% |
| 40 V current error | 16.20605% | 20% | 10% |
| Maximum normalized KCL at 31 reference points | 8.01838e-11% | 1% | 0.1% |
| Vg=8 verdict | Fail | Fail | Fail |

At 40 V, Vela gives 4.869356375653821e-4 A/um and the frozen SDevice reference
gives 4.19027790569749e-4 A/um. The 40 V two-gate ratio error, using the fresh
Vg=4 fixed-growth curve, is 13.89269%: the 15% engineering ratio limit passes,
while the 8% final ratio limit fails. The complete two-gate D5 score fails
both engineering and final acceptance because the current-accuracy criteria
remain unsatisfied.

The maximum nonzero-point current difference from the previous complete Vg=8
curve is 9.19133e-13 relative. This rerun reproduces the existing accuracy gap
while demonstrating uninterrupted physical drain continuation with a verified
frame change. Full-field equivalence was checked at the frame change; the
comparison to the historical full curve is a terminal-current comparison.

## Evidence

- [Completed controller ledger](../../reference_staging/templates_ldmos_d5_vg8_rerun_20260908/fixed/ledger.json)
- [Full integrity audit and numerical summary](../../reference_staging/templates_ldmos_d5_vg8_rerun_20260908/audit_summary.json)
- [Vg=8 Id-Vd CSV](../../reference_staging/templates_ldmos_d5_vg8_rerun_20260908/vg8_score/curve.csv) and [score/provenance](../../reference_staging/templates_ldmos_d5_vg8_rerun_20260908/vg8_score/summary.json)
- [Complete two-gate D5 score](../../reference_staging/templates_ldmos_d5_vg8_rerun_20260908/d5_score/summary.json)
- [Source verification](../../reference_staging/templates_ldmos_d5_vg8_rerun_20260908/source_verification.json)
- [Figure PDF](../../reference_staging/templates_ldmos_d5_vg8_rerun_20260908/vg8_rerun_comparison.pdf)

![Vg=8 full rerun, two-gate ratios and physical continuation steps](../../reference_staging/templates_ldmos_d5_vg8_rerun_20260908/vg8_rerun_comparison.png)
