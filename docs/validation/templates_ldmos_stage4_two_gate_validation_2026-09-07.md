# Templates/LDMOS Stage-4 two-gate D5 validation

Date: 2026-09-07 (Asia/Shanghai). Source HEAD `959c1c9`.

The Vg=8 V drain sweep completed all 31 exact reference points from 0 to 40 V. Combined with the previously qualified Vg=4 curve, the unchanged D5 scorer gives engineering **fail**, final **fail**.

IALMob remains disabled because complete D5 final qualification has not passed.

## Frozen configuration

Windows UCRT64 Release; runner SHA256 `37e2f099b08864a3d2be8aef2e4047ba479b3f6ef22a11967e08cfb9f6a733b4`. Actual Eigen SparseLU/COLAMD selected by `VELA_LINEAR_SOLVER=sparselu`, with L2 row/column equilibration. Exact imported mesh: 10,241 nodes, 19,782 triangles and 30,022 edges. External AverageBox carrier couples, barycentric volumes, material-local Poisson charge and legacy node-local contact reconstruction are retained.

300 K Fermi/OldSlotboom, SRH/Auger and constant-field GradQF HFS with live repaired mobility derivatives. Predictor, IALMob, quantum, avalanche and heating are off. Geometry is in um; terminal current is A/um, state potentials V and densities m^-3.

Original global residual ceilings remain psi <=5e-8, electron <=1e-11, hole <=3e-10. The additional local carrier-row ratio and terminal KCL/max-current limits remain <=1e-8. No physical parameter, acceptance tolerance, C++ source or reference baseline changed.

## Complete two-gate score

| Metric | Vg=4 | Vg=8 | Final limit |
| --- | ---: | ---: | ---: |
| Median current error | 2.130296% | 17.280895% | 5% |
| P95 current error | 9.312277% | 18.643299% | 12% |
| Low-voltage resistance error | 4.145013% | 3.984756% | 10% |
| 40 V current error | 2.031173% | 16.206049% | 10% |

The 40 V two-gate current-ratio error is **13.892691%**, against the unchanged 8% final limit. Maximum normalized KCL over both complete curves is 1.9004e-10% (zero bias included). Relative current metrics use the 30 nonzero points; no interpolation is used to score Vela.

The resistance metric retains the historical first nonzero V/I definition. Gate-ratio acceptance uses the 40 V endpoint. Raw solver voltages and the gauge offset are preserved in per-point provenance. Physical drain voltage is raw drain minus raw source. After that exact coordinate conversion, only the existing <=32 binary64 ULP endpoint matching convention is used.

![Both gate curves and ratios](../../reference_staging/templates_ldmos_vg8_20260907/final_review/two_gate_comparison.png)

## Initialization and path checks

The archived Vg=8, Vd=0 gate-prebias state was reclosed with the current solver. Its first attempt reached the original global block ceilings after five updates but failed the strict local carrier-row guard. Its converged Poisson potential was retained, both carrier QFs and their reference/increment coordinates were restored to the exact zero-drain equilibrium condition, and a new full coupled solve independently accepted the state at iteration zero. Repetition gave a byte-identical state and four exactly zero computed terminal currents. No terminal output was edited.

Each accepted physical transfer uses at most 0.1 V and clips to the next exact reference point. Its QF cap equals its actual target-minus-parent voltage. Failed direct attempts undergo a strict same-bias reclose; density recovery is eligible only after the original global ceilings pass and local rows still fail. Unqualified results are rolled back and subdivided. Failed child costs are included in timing.

At 0.1 V the repeated 100 mV control has identical Newton traces and state. Its comparison with 40 x 2.5 mV transfers passes the preserved field/current criteria: maximum electron density relative difference 2.84486e-14, normalized terminal-current difference 1.272e-12.

At 4 V an independent 2.5 mV subdivision of the final transfer also passes: maximum electron density relative difference 5.66132e-10, normalized terminal-current difference 6.98486e-14. This is a local endpoint control, not an independent complete 0-to-4 V history. It does not explain away the Sentaurus current discrepancy.

## Runtime and verification

The Vg=8 production path's summed child wall time was **5355.807 s (89.263 min)** and child CPU time 4542.922 s. There were 424 accepted physical transfers, 868 child executions, 5172 applied Newton updates and 8 rollbacks.

The original controller interval was 3721.137 s; the separate 28-to-40 V controller interval was 1837.621 s. A 7589.325 s investigation gap separates those intervals. The full path is therefore not reported as one uninterrupted run.

The shifted-frame controller also had an administrative interruption after a qualified 29.533333333333303 V solve: Windows rejected an atomic progress-file replacement. Both the committed and pending ledgers were preserved. The pending accepted state, parent hash, original block norms and KCL were verified before resuming. A bounded retry was added to the generated progress writer; no solver state or physics was changed. The shifted-controller interval includes that interruption, while child wall/CPU sums do not.

Sum includes failed production trials and the chosen frame=28 baseline/28 V transfer; excludes frame=14 controls, other diagnostics, old gate prebias and idle investigation gap. Controller intervals are reported separately. Not a controlled matched performance benchmark.

Python helper syntax checks pass. The eight `SentaurusAblationSummaryTest` regressions pass. The complete real-device queue, exact-grid score, both independent local controls and the final child/configuration/state/lineage audit pass their integrity checks. Numerical score failures remain visible separately from integrity checks. No solver rebuild was performed solely for generated experiments.

## Evidence

- [Original controller ledger, preserving its failure](../../reference_staging/templates_ldmos_vg8_20260907/continuation_r2/ledger.json)
- [Gauge continuation and composed exact points](../../reference_staging/templates_ldmos_vg8_20260907/gauge_continuation/ledger.json)
- [Complete integrity audit](../../reference_staging/templates_ldmos_vg8_20260907/gauge_continuation/integrity_audit.json)
- [Unchanged two-gate scorer result](../../reference_staging/templates_ldmos_vg8_20260907/d5_score/summary.json)
- [Vg=8 per-point provenance](../../reference_staging/templates_ldmos_vg8_20260907/vg8_score/provenance.json)
- [0.1 V path comparison](../../reference_staging/templates_ldmos_vg8_20260907/local_controls/comparison.json)
- [4 V endpoint comparison](../../reference_staging/templates_ldmos_vg8_20260907/midbias_control/comparison.json)
- [Sentaurus D5 replay input manifest](../../reference_staging/templates_ldmos_vg8_20260907/sdevice_bundle/manifest.json)

The Sentaurus Vg=8 Extrapolate on/off replays and the independent 4 V endpoint solve are separate D5-only reference controls; they do not enable IALMob or replace the frozen scoring reference. All three completed and their results follow below. Generated simulation outputs remain ignored.

## High-voltage precision control and coordinate continuation

The original ground-referenced controller stopped at 27.99270833333334 V after eight rollbacks; a further retry would have violated the predeclared 2.5 mV retry minimum. Its final attempted reclose had electron block norm 1.3489701393e-11 against the unchanged 1e-11 ceiling, despite zero violating local carrier rows. The failed state was not accepted.

At the last accepted state, node 739 has psi about 28.567860 V, electron QF about 27.992705 V, and a stored electron QF increment of -1.64e-16 V. Two isolated controls subtract 14 V and 28 V, respectively, from every physical potential and contact bias. All terminal voltage differences remain unchanged. Both same-bias coupled solves qualify and match the original state under the preserved field/current criteria.

Both controls then reach physical 28 V after one strict reclose. Their cross-comparison gives maximum carrier-density relative difference 2.77132e-13, maximum potential difference 1.06581e-14 V, and normalized terminal-current difference 4.87153e-13. All equivalence and KCL criteria pass. This is evidence that finite-precision coordinate representation contributes to the continuation floor; it does not identify the source of the converged Vg=8 current discrepancy against Sentaurus.

The accepted 28 V control in the 28 V shifted frame seeds the remaining physical 28-to-40 V segment. Its raw source/substrate are -28 V, gate -20 V, and drain 0-to-12 V. Raw solver states keep this coordinate system. Inverse-translated states are used only for field comparison; scorer currents come directly from qualified raw solves.

Inactive oxide carrier QFs are zero placeholders. Their physical/reference/increment slots must remain internally consistent when translating a state. Two preliminary input files reset reference slots without resetting increments and were rejected before Newton; the converter was corrected and the retries use new evidence directories. These input errors are separate from convergence failures and are excluded from the production path.

- [Baseline coordinate equivalence](../../reference_staging/templates_ldmos_vg8_20260907/gauge_validated_comparison.json)
- [Qualified 28 V arrivals](../../reference_staging/templates_ldmos_vg8_20260907/gauge_transfer_r2.json)
- [28 V cross-gauge comparison](../../reference_staging/templates_ldmos_vg8_20260907/gauge_28_cross_comparison.json)

## Independent Sentaurus reference controls

Three D5-only controls completed in the authorized Sentaurus T-2022.03-SP2 VM. Uploaded input hashes matched the prepared decks. The downloaded result archive SHA256 is `f25ca3f6ec1d0f70326b9ddaaed96503d1ffa443a46c934245058f5ed9935f50`. All runs use one assembly thread and one solver thread. The original scoring reference is unchanged.

| Control | Forward steps | Forward Newton updates | Forward solver time | Process wall time |
| --- | ---: | ---: | ---: | ---: |
| adaptive_on | 49 | 216 | 97.23 s | 320.40 s |
| adaptive_off | 63 | 656 | 213.34 s | 418.44 s |
| endpoint4_on | 40 | 114 | 40.50 s | 244.18 s |

The Extrapolate-on replay reproduces all 31 frozen reference currents exactly. Turning Extrapolate off changes the nonzero reference currents by at most 2.61857e-08%. The independently solved 4 V endpoint differs from the frozen reference by 1.98175e-11%; Vela remains 12.407682% higher than that independent endpoint.

The measured current discrepancy is therefore not explained by Extrapolate or reference-point sampling. Forward solver timers exclude startup/license/output costs; process wall timers include them. The VM and host ran other validation work during this interval, so these are observed execution costs, not a controlled cross-program performance benchmark.

- [Sentaurus currents, iteration and timing audit](../../reference_staging/templates_ldmos_vg8_20260907/sdevice_comparison.json)

## Independent field-seed limitation

The Release HDF5 importer successfully exported the independent 4 V Sentaurus endpoint. All 10,241 vertex coordinates match the Vela mesh exactly. The existing qualified field mapper was reused: potentials in V, carrier densities converted from cm^-3 to m^-3 by 1e6; no fitted reference shift was applied.

The accepted Sentaurus and Vela fields differ on 5,723 silicon nodes: the P95 absolute electrostatic-potential difference is 0.218192 V and its maximum is 0.251614 V. These field differences do not by themselves identify the responsible operator or model term.

Loading the Sentaurus field state as a Vela initial guess did not qualify in the bounded control: the direct attempt used 160 updates and one strict reclose stopped after 10 more. The direct final block norms were approximately 5.7866 / 19726.97 / 34.1332; the reclose was also far above the original ceilings. The conditional density recovery was therefore ineligible. This failed control is retained and is not included in the curve score. It cannot prove that the two initial states converge to the same Vela solution.

- [Independent seed plan and field differences](../../reference_staging/templates_ldmos_vg8_20260907/sentaurus_seed4/plan.json)
- [Bounded independent-seed result](../../reference_staging/templates_ldmos_vg8_20260907/sentaurus_seed4/comparison.json)

The next investigation should separate initial-guess robustness from the accepted Vg=8 current accuracy gap using the saved 4 V fields and the unchanged classical contract. IALMob remains gated while complete D5 final acceptance fails. The current T-2022.03-SP2 `Siliconc100.par` was downloaded for preparation only; no IALMob simulation or solver model change was made.
