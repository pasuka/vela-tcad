# Templates/LDMOS Stage-4 matched-step validation

Date: 2026-09-07 (Asia/Shanghai). Completed on source HEAD `959c1c9`.

The 5.3 V short-step blocker is resolved by matching the physical quasi-Fermi
(QF) update cap to the actual drain-voltage increment, followed by strict
same-bias closure. Two repetitions give identical Newton traces and final
states. The resumed Vg=4 curve now covers all 31 exact reference points
through 40 V and passes both engineering and final single-gate limits.
Current error median is 2.13030%, P95 9.31228%, endpoint error 2.03117%;
the existing all-point KCL gate also passes, including physical zero bias.

This updates the outstanding short-step and zero-bias findings in the earlier
[recovery report](templates_ldmos_stage4_recovery_validation_2026-09-07.md).
That dated report remains unchanged. No solver source, model default or
acceptance threshold was changed in this investigation. The result qualifies
this checkpoint-resumed Vg=4 validation workflow; it does not yet qualify the
two-gate D5 comparison. Vg=8 and the two-gate current ratio remain pending,
and IALMob remains disabled.

## Short-step diagnosis

The comparison uses the same accepted 5.299999999999995 V state, SHA256
`9b166fe7c7b2ac70a57a52f0121e70aad5f71c201d2f48827c3ae09458394918`,
and the same Release runner/Jacobian/physics. Only the physical QF update cap
changes (`solver.quasi_fermi_update_limit_V`). The target is 5.33333333333333 V, so the applied voltage increment is
0.033333333333334991 V.

At the first iteration, node 5558's raw electron QF update is 0.0709255465 V.
With the old 0.1 V cap the line search applies half of that, 0.0354627733 V,
overshooting the applied drain increment. Its electron row becomes -2145.49.
Matching the cap to the drain increment applies 0.0333333333 V and leaves
that row at -3.93768e-10. At node 702 the corresponding applied update falls
from 0.0424897928 to 0.0333333333 V and its row from -3301.14 to -1.16621e-10.
Both experiments use the same raw Newton direction; no new Jacobian change
is needed to demonstrate this difference.

| Eight-update prefix | QF cap | Electron block after eight updates |
| --- | ---: | ---: |
| 33.333 mV transfer, previous cap | 0.1 V | 1.10978e4 |
| Same transfer, matched cap | 0.033333333333334991 V | 9.81135e-11 |
| Same transfer, smaller arbitrary cap | 0.01 V | 1.47218e3 |

The full 160-update-budget matched-cap run stops after eight applied updates
at a remaining precision floor. One strict same-bias reclose, using the existing
`quasi_fermi_recenter_on_initial_state=true`, accepts the state: block norms 1.35272e-9 / 5.75717e-13 / 5.57570e-25, zero local
violations, worst local ratio 3.48873e-13. Both transfer and reclose traces
are byte-identical in the repeated experiment.

Against the independently reclosed fine-path reference, the maximum electron
density relative difference is 7.11310e-14, maximum phin difference is
1.99840e-15 V, normalized port-current difference is 1.18707e-13 and KCL
ratio is 1.77438e-13. These pass the earlier 1e-6 density, 1e-7 V potential,
1e-8 normalized current and 1e-8 KCL limits. A matched 4.166667 mV step also
passes the full residual/local/KCL checks in five updates.

A native two-point sweep reproduces the same failed short-step initial
residual and early behavior as the single-target restart. The problem is
therefore not resolved by changing only the wrapper's startup convention.
The first raw linear solve has normwise backward error 2.55555e-16, while a
very weak row has a large componentwise error; this does not prove all weak
linear equations are accurate. The controlled cap change establishes a
globalization remedy for the tested transfer, not a universal Jacobian proof.

Larger matched caps were probed separately with 12-update diagnostic prefixes.
In that tested sequence, a 0.2 V prefix plus strict reclose needs the existing
conditional Gummel density recovery to qualify. A 0.5 V prefix followed by
a 160-update-budget reclose does not qualify. A 1.333333 V diagnostic prefix
also remains far from convergence. These are bounded controls, not exhaustive
tests of larger-step policies; none is used in the complete continuation.

## Applied continuation strategy

The production fine-path sweep remains administratively stopped. A separate
queue loads its last complete accepted 11.899167 V checkpoint, SHA256
`4600e6a9a34dcf00eebbc6fd8dca10fe515dc2a7d4dc1ddd791c2c05bfd2bf37`,
and recloses it with the stronger local criterion. Source bias is read from
the common physical drain-contact QF, not the rounded checkpoint filename.

Each nominal transfer is at most 100 mV and clips at the next exact Sentaurus
reference. The QF cap is its actual target-minus-parent voltage increment,
including floating-point endpoint bookkeeping. Direct failure is followed
by one same-bias strict reclose. Existing Gummel density recovery is eligible
only if every original global block ceiling already passes and local rows
still fail. Otherwise the controller rolls back and halves the step, stopping
below the declared 2.5 mV retry minimum. Normal clipping at a requested exact
point is distinct from retry subdivision.

All accepted states satisfy psi <=5e-8, electron <=1e-11 and hole <=3e-10;
qualified local carrier ratios <=1e-8; terminal KCL/max-terminal <=1e-8.
No acceptance tolerance or physical parameter is changed. The original
160-update budget, live repaired mobility Jacobian and QF recenter option are
retained. Predictor and IALMob remain off. The physical parent/target/cap and hashes are written in
the controller ledger; a child's startup `attempted_step_V=0` is bookkeeping,
not evidence of a physical zero-voltage transfer.

At 12 V the new path passes the independent 2.5 mV-path reference comparison;
maximum electron density difference is 1.50327e-9, below the preserved 1e-6
limit. Old exact states below 12 V are reclosed under the same current
solver and stricter row guard for scoring. They are not interpolated from
internal steps. Consequently the final score describes a checkpoint-resumed
curve, not an uninterrupted fresh 0-to-40 V run.

## Completed continuation and single-gate score

The queue completed at 15:33:11 +08:00 and exited with code 0. It accepted
296 physical transfers from 11.899166666666455 to 40 V, reaching all 22
remaining exact reference points. Of those transfers, one accepts directly
and 295 require a strict same-bias reclose. No rollback/bisection or Gummel
density recovery was needed. The last 33.333333 mV transfer to 40 V takes
six direct updates plus one reclose update and passes all gates.

The [completed ledger](../../reference_staging/templates_ldmos_shortstep_20260907/continuation_from_11p899/ledger.json)
records 592 child executions including initial closure, with 3243 applied
Newton updates. Transfer totals including reclose have minimum 3, median 10,
and maximum 40 updates. Direct attempts have median 9 updates; reclose has
median 1 and maximum 3. Failed direct-attempt costs are included. This is
not evidence that a direct-only 100 mV sweep succeeds.

| Timing for the resumed segment | Measured |
| --- | ---: |
| Controller elapsed wall time (14:18:02 to 15:33:11) | 4508.673 s / 75.145 min |
| Sum of timed child wall times | 4458.745 s / 74.312 min |
| Sum of child CPU user + kernel times | 3300.469 s / 55.008 min |

Summed child time includes initial closure, direct attempts and recloses.
It excludes the independent fine-reference 12 V reclose, prefix preparation,
short-step diagnostics, scoring and historical seed preparation. Controller
elapsed time also includes its bookkeeping and the independent 12 V check.
This resumed interval is not a paired full-sweep performance benchmark;
neither a full-curve speedup nor SDevice runtime parity is inferred. The prior
[local performance controls](templates_ldmos_stage4_step_performance_2026-09-07.md)
remain the bounded paired timing evidence.

The final continuation audit checks all 592 completed children, 296 accepted
transfers, original residual ceilings, live mobility settings, parent hashes,
cap/voltage lineage and 1185 distinct hashed files. Maximum KCL ratio over all
accepted internal transfers is 2.47515e-12, below 1e-8; maximum qualified local
carrier-row ratio is 9.11653e-9, below 1e-8. All original block ceilings pass.
The 40 V state SHA256 is
`e5fe06b9a086de8364932e2250811e390b4aff4e46925e18a750188d30998be4`.

For scoring, nine independently closed prefix points from 0 through
10.666667 V are combined with the 22 new exact points. All 31 state files
have 10241 unique finite node records with nonnegative densities. The score
checks 124 terminal-current rows and preserves the actual raw solver bias in
per-point provenance; only the pre-existing <=32 binary64 ULP endpoint
matching convention is applied to the reference key. There is no interpolation.
The zero-bias state is prepared and qualified as described below.

| Existing Vg=4 metric | Measured | Engineering limit | Final limit |
| --- | ---: | ---: | ---: |
| Median absolute current error | 2.13030% | 15% | 5% |
| P95 absolute current error | 9.31228% | 25% | 12% |
| Low-voltage resistance error | 4.14501% | 20% | 10% |
| 40 V current error | 2.03117% | 20% | 10% |
| Maximum normalized KCL, all 31 points | 1.90040e-10% | 1% | 0.1% |
| Vg=4 verdict | **Pass** | **Pass** | **Pass** |

Maximum current error is 10.18579% at Vd=4 V. At 40 V, Vela gives
2.2939992207505347e-4 A/um and Sentaurus gives 2.24833171719662e-4 A/um.
The scorer uses 30 nonzero points for relative current errors, and all
31 points for KCL. Its historical `low_vd_differential_resistance` field is
the ratio of V/I at the first resolved nonzero point (1.333333 V), not a newly
estimated derivative. The formulas and limits in
[the existing scorer](../../scripts/analyze_templates_ldmos_stage4_d5.py)
are unchanged. Two-gate ratio limits are not evaluated using one curve.

The early shoulder around 4 to 6.7 V retains an approximately 8 to 10% current
difference from Sentaurus. This is within the existing Vg=4 aggregate limits;
it is not removed by the convergence remedy. Full-field equivalence against
an independent fine path was checked at the local control endpoints (including
5.333333 and 12 V), not at every new internal high-bias state.

![Vg=4 comparison at all 31 exact shared points](../../reference_staging/templates_ldmos_shortstep_20260907/vg4_score/vg4_comparison.png)

## Zero-bias equilibrium closure

The earlier zero-bias reclose had converged Poisson potential but tiny
nonconstant QFs and currents around 1e-25 A/um, producing a failed normalized
KCL ratio. With all ohmic contacts at 0 V, equal carrier QFs are the exact
equilibrium condition for the enabled isothermal, unilluminated model.

The new input keeps the already closed Poisson potential and restores both
physical QFs and their reference/increment coordinates to zero. Stored input
currents are never edited; the coupled solver independently evaluates its
residuals, reconstructs densities and computes terminal currents. It accepts
the state at the initial residual check with the original block limits and
the added local guard. All four computed terminal currents are exactly zero.
The unchanged KCL scorer reports zero, with no zero-bias exemption or new
denominator floor. A second solve produces a byte-identical state.

This is restricted to the physical zero-drain equilibrium. From that state,
a 2.5 mV nonzero-drain solve plus strict reclose also succeeds, with
Id=2.3156755854e-7 A/um and KCL ratio 9.96151e-15. The qualification does not
depend on suppressing finite-bias current.

## Reproduction context

Branch `codex/templates-ldmos-phase-a`, source HEAD `959c1c9` at task start.
No solver source change was needed for these controls. Windows UCRT64 Release,
actual Eigen SparseLU/COLAMD, L2 row/column equilibration. Runner SHA256
`37e2f099b08864a3d2be8aef2e4047ba479b3f6ef22a11967e08cfb9f6a733b4`.

Exact imported mesh: 10241 nodes, 30022 edges, 19782 triangles; 16237 external
AverageBox carrier couples, including 3548 zero couples. Barycentric volumes,
material-local Poisson charge, legacy node-local contact reconstruction.
300 K Fermi/OldSlotboom, SRH/Auger, constant-field HFS with live repaired
mobility derivatives, edge QF drive and existing contact-field fallback.
No predictor, IALMob, quantum, avalanche or heating. Geometry is in um;
reported terminal currents use A/um, state densities m^-3 and potentials V.

Generated evidence stays in `reference_staging/templates_ldmos_shortstep_20260907/`.
The input harness refuses to overwrite case directories. `continue_reference.py`
records the physical state lineage. `audit_continuation.py` checks completed
children without treating incomplete intervals as passes. `score_vg4.py` reuses
the existing D5 formulas and limits and refuses to score an unfinished queue.


## Evidence and checks performed

The following local evidence records are generated artifacts and remain in
ignored staging, rather than being committed simulation outputs:

- [Short-step controlled comparison](../../reference_staging/templates_ldmos_shortstep_20260907/shortstep_cause_summary.json),
  [fine-state equivalence](../../reference_staging/templates_ldmos_shortstep_20260907/matched_cap_comparison.json),
  and [12 V comparison](../../reference_staging/templates_ldmos_shortstep_20260907/continuation_from_11p899/comparison_12.json).
- [Zero-bias preparation plan](../../reference_staging/templates_ldmos_shortstep_20260907/zero_equilibrium_plan.json)
  and [nonzero entry check](../../reference_staging/templates_ldmos_shortstep_20260907/zero_entry_validation.json).
- [Continuation plan](../../reference_staging/templates_ldmos_shortstep_20260907/continuation_from_11p899/plan.json)
  and [completed continuation audit](../../reference_staging/templates_ldmos_shortstep_20260907/continuation_from_11p899/integrity_audit.json).
- [31-point score and provenance](../../reference_staging/templates_ldmos_shortstep_20260907/vg4_score/summary.json)
  and [per-point comparison CSV](../../reference_staging/templates_ldmos_shortstep_20260907/vg4_score/point_comparison.csv).
- [Final audit](../../reference_staging/templates_ldmos_shortstep_20260907/final_audit.json):
  passing local controls and prefix closures, input/configuration hashes,
  identical repeated Newton/state files, exact preservation of non-QF
  zero-seed fields, 31 scored checkpoints and 124 terminal records.

The current schema was checked before using historical reports. Python helper
syntax was checked; the complete real-device queue, 31-point score, repeated
short-step controls, zero-equilibrium/nonzero-entry controls, and both final
audits pass. The comparison plot was visually inspected. Documentation links
and `git diff --check` were checked. No C++ source or tracked Python scorer
changed, so no additional solver rebuild or unit-test rerun was performed
solely for this experiment/report. The preceding committed repair's tests are
recorded separately in the earlier recovery report.

From the worktree root, the audit commands can be repeated against the saved
results without running a new device sweep:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
$env:PYTHONUTF8 = '1'
python reference_staging/templates_ldmos_shortstep_20260907/audit_continuation.py
python reference_staging/templates_ldmos_shortstep_20260907/finalize_evidence.py
```

These commands refresh only audit metadata. The original execution and score
helpers deliberately refuse to overwrite existing output directories; a new
experiment must use a separate output directory and preserve its input lineage.
The reusable rule is the recorded physical delta-matched cap plus strict
same-bias closure, not a change to the global C++ or template defaults.

The original production checkpoint and its administrative-stop evidence remain
unchanged. The new continuation has finished and no experiment runner remains
active in this worktree. The next D5 work is Vg=8 and the two-gate current ratio;
a fresh uninterrupted 0-to-40 V strategy qualification and broader performance
claims would require separate evidence. IALMob remains behind complete D5
qualification.
