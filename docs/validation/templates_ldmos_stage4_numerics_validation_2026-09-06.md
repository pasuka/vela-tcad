# Templates/LDMOS Stage-4 high-bias numerical controls

Date: 2026-09-06. Status: all four requested numerical controls completed.
The SDevice VM runs finished with exit code 0. This is diagnostic evidence, not
qualification of the 31-point D5 curve or a change to its production settings.

## Result

The high-bias A/B establishes that omitting mobility-field derivatives is a
major cause of the long Newton tail in the two tested transfers. Enabling the
existing live derivative path reduces 85 updates to 13 at 3.371667 V and
converges in 14 updates at the 3.854167 V transfer where the lagged path fails.
No predictor, IALMob, residual-ceiling change, or physical-parameter change was
used. The first pair's drain-current relative difference is 1.34e-13.

The SDevice Extrapolate controls also finish at 40 V with unchanged physics and
stopping criteria. Disabling extrapolation increases Newton updates by 1.65x
on the fixed common sequence and by 2.70x on adaptive sweeps. It is a meaningful
initial-guess advantage, but does not by itself explain the Vela long tails.
The Vela mobility-linearization and state-update-resolution findings above
remain the directly demonstrated local causes.

The live Jacobian still has a measurable finite-difference-step sensitivity
near the drain. These results do not establish a universally accurate Jacobian
or qualify a production switch across the complete voltage range.

## Configuration and provenance

- Worktree: `D:\code-repo\vela-tcad\.worktrees\templates-ldmos-phase-a`.
- Branch: `codex/templates-ldmos-phase-a`.
- HEAD: `345fda30d541003f4aa3a08341dca59418f06ded`; local diagnostic changes are
  uncommitted. Pre-existing script/test/evidence changes were preserved.
- Mesh: exact imported topology, 10241 nodes, 30022 edges, 19782 triangles.
  Carrier transport uses the existing external AverageBox couples, with
  material-local Poisson charge volumes and `legacy_node_local` contacts.
- Physics: 300 K, Fermi statistics, OldSlotboom, SRH/Auger and live residual HFS;
  no IALMob, quantum correction, avalanche or heating.
- Units: `unit_scaling`; external geometry in um. Terminal currents below are
  explicitly taken from `current_total_A_per_um`. Raw solver residual norms are
  not terminal currents and are not directly comparable to SDevice RHS norms.
- Absolute block ceilings: psi 5e-8, electron 1e-11, hole 3e-10, all enforced.
- Newton budget: 160; QF recenter enabled; predictor absent/off.
- Build: Windows UCRT64 Release. HDF5, SuiteSparse SPQR/UMFPACK capabilities
  detected. Actual A/B linear backend was explicitly pinned to Eigen SparseLU,
  with COLAMD and the configured L2 row/column equilibration.
- Phase 1 used the unchanged original runner, SHA-256
  `0569876d846792777390a3d134986bf6df8314af85acfaf87a5f0166d12abb36`.
  It is preserved as `build-release/vela_example_runner_stage4_baseline_345fda3.exe`.

All generated inputs, states, traces and machine-readable summaries are under
`reference_staging/templates_ldmos_numerics_validation_20260906/`. The aggregate
is [validation_summary.json](../../reference_staging/templates_ldmos_numerics_validation_20260906/validation_summary.json).
It distinguishes original, profiling, and later row-localization binary hashes.

## 1. Fixed-transfer mobility-Jacobian A/B

Each pair loads the identical accepted parent checkpoint and performs one
Newton solve directly at the target voltage. There is no preliminary reclose,
intermediate voltage, predictor or retry subdivision. A zero `attempted_step_V`
in this single-point runner's attempt CSV denotes its startup bookkeeping; the
physical parent-to-target delta is 2.5 mV, recorded in the experiment manifest.
Both variants have identical initial-state hashes and initial residuals.

| Parent -> target Vd (V), Vg=4 V | Lagged derivatives | Live derivatives |
| --- | --- | --- |
| 3.369166666666683 -> 3.371666666666683 | pass, 85 updates, 40.16 s host wall | pass, 13 updates, 7.71 s |
| 3.851666666666673 -> 3.854166666666673 | fail on attempt 154, after 153 accepted updates, 71.37 s | pass, 14 updates, 8.02 s |

| Case | Final electron block | Id (A/um) | KCL / max terminal magnitude |
| --- | ---: | ---: | ---: |
| tail, lagged | 9.92593e-12 | 1.727167047567871e-4 | 2.32442e-13 |
| tail, live | 8.12779e-12 | 1.727167047567640e-4 | 9.85774e-14 |
| floor, lagged | 1.00780e-11 | no accepted endpoint | no accepted endpoint |
| floor, live | 7.28009e-12 | 1.829987728605768e-4 | 1.60103e-13 |

The tail pair's relative Id difference is 1.33893e-13. Every passing case also
satisfies the unchanged Poisson and hole ceilings. The failed case reproduces
the previous 3.854167 V failure exactly: 13 unsuccessful line-search trials at
Newton attempt 154, with `block_absolute_convergence_line_search_rejected`.

This is evidence for retaining mobility feedback in the high-bias
linearization. It does not overturn the historical low-bias observation that
the lagged path was more robust from equilibrium.

Evidence: [phase1_metrics.json](../../reference_staging/templates_ldmos_numerics_validation_20260906/phase1_metrics.json)
and [phase1_status.json](../../reference_staging/templates_ldmos_numerics_validation_20260906/phase1_status.json).

## 2. Hotspot JVP and perturbation-size audit

The existing JVP probe always compared against the live residual. A new,
opt-in `freeze_transport_mobility` diagnostic holds each edge's base-state
mobility fixed for both perturbations. It requires lagged mobility derivatives,
edge SG and no coupled avalanche. The ordinary residual and default JVP mode
retain their previous behavior.

The audit covers nodes 3132 and 5569 and their one-ring patch: 12 nodes, with
individual psi, phin and phip directions. Both the accepted parent and archived
failed-final states were checked, using live/live and lagged/frozen comparisons.
At the failed state, both hotspots were also scanned at perturbations from
1e-4 to 1e-9 V. Total: 204 directional evaluations.

Relative errors below are actual `||Jv-FD|| / ||FD||`, not the probe's legacy
normalization using `max(1, ||FD||)`.

| Failed-state direction family, h=1e-7 V | Maximum matched lagged error | Maximum live error |
| --- | ---: | ---: |
| psi, 12 directions | 3.1363e-9 | 3.1363e-9 |
| phin, 12 directions | 8.0192e-10 | 1.5424e-4 |

At node 5569, the live phin error is about 1.44e-4 at h=1e-7 V and tends toward
1.53e-4 at h=1e-8 and 1e-9 V. It is approximately 2e-13 when the diagnostic
step equals the production finite-difference step, h=1e-6 V. Matching that
same secant step is not independent proof of the infinitesimal derivative;
the finer-step plateau remains an open Jacobian-accuracy finding. The matched
lagged phin error at h=1e-7 V is only 4.56e-12 at this node.

At node 3132, live and lagged phin JVP norms are 1.36536e-6 and 1.11155e-5 for
the same 1e-7 V perturbation. The magnitude of the omitted mobility response
there is consistent with a substantial change in the Newton correction. A
local column-norm ratio is not an eigenvalue or a proof of a global rate.

Hole directions have larger relative errors, up to 0.67, but the largest-error
rows are Poisson rows with changes/errors around 1e-20 to 1e-23. These directions
probe extremely small minority-carrier effects through a much larger residual.
They are retained as resolution-limited findings, not reported as a global
Jacobian pass or evidence of the electron-tail cause.

Evidence: [phase2_by_variable.json](../../reference_staging/templates_ldmos_numerics_validation_20260906/phase2_by_variable.json),
the four `phase2/*/jvp.csv` files, and their `sample_rows.csv` files.

## 3. SDevice Extrapolate controls

All four Vg=4 V controls completed to Vd=40 V with exit code 0 and no rejected
voltage steps. The VM was checked after completion and had no remaining
SDevice process. The installed executable is T-2022.03-SP2; logs confirm double
precision, one assembly thread, one solver thread, blocked decomposition and
no incomplete Newton. The installed T-2022.03 user guide confirms that full
mobility derivatives are the default and `-Derivatives` disables them;
extrapolation is off by default. The original D5 deck explicitly enables it.

| Case | Accepted forward steps | Forward Newton updates | Median / step | Range / step | Forward solve time (s) | Process wall time (s) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| fixed_on | 49 | 237 | 5 | 2--7 | 142.91 | 393.32 |
| fixed_off | 49 | 391 | 8 | 5--19 | 171.72 | 406.54 |
| adaptive_on | 49 | 239 | 6 | 2--7 | 70.59 | 272.98 |
| adaptive_off | 63 | 645 | 8 | 5--23 | 181.07 | 382.64 |

Iteration zero is an initial residual evaluation, not a Newton update. The
fixed runs have another 49 same-bias t=0 updates plus one loaded-state reclose
each, giving all-update totals 287/441. Adaptive runs have two such startup
updates each, giving totals 241/647. Every printed forward damping factor is 1.
Forward solve time sums the per-step `Total time` lines. Process wall time
also includes about 194--195 s of startup/license acquisition in each run;
approximate time after license checkout is 198.32/211.54 s for fixed on/off
and 77.98/188.64 s for adaptive on/off. These are single-run observations,
not a repeated benchmark or a hardware-normalized comparison with Vela.

### Fixed common sequence

Both decks use the same 49 target voltages reconstructed from the original
*printed* internal-step log. Those displayed voltages are rounded; this is
not a replacement for the original 31 exact output points. Each target uses
one Quasistationary segment with InitialStep=MinStep=MaxStep=1, so failures
cannot silently introduce intermediate voltages. The on log confirms 48
history reuses across segments and 48 actual extrapolations; off has zero.
Each segment also performs the same-bias solve accounted for above.

Disabling extrapolation raises the forward count from 237 to 391 (1.65x).
The off outlier of 19 updates occurs at 0.062336 -> 0.092432 V; most high-bias
off steps take 8--9 updates. At 3.24136 -> 3.963552 V, the count is 5/8 and
initial RHS is 8.00e8/4.16e12 for on/off. This directly shows the improved
initial guess, without changing the Jacobian or stopping criteria.

At common nonzero biases, maximum relative drain-current difference is
2.1833e-8 and its median is 5.3677e-14. Raw current rows retain both arrival
and subsequent same-bias reclose values. Common-point comparisons take the
first arrival; duplicate rows and their changes are explicitly preserved.

### Adaptive sequence

Both use the original controls: normalized InitialStep=2.5e-4, MinStep=1e-6,
MaxStep=0.05, Increment=1.35 and 31 exact current-output points over 0--40 V.
These correspond to a 10 mV initial step, 40 uV minimum, and 2 V maximum
before clipping to the output lattice. The realized maximum step is about
1.333333 V in both runs. Rounded internal log values must not be interpreted
as higher-precision bias evidence.

Off uses 63 accepted steps versus 49 for on, with no rejected/retried steps.
Newton updates rise from 239 to 645 (2.70x). Much of the extra work is below
1 V: on/off take 12/22 steps and 35/292 updates there. Several off steps near
0.09--0.23 V take 23 updates. The step-size median is about 1.33332/0.395032 V
for on/off, as measured from the printed internal sequence.

The on internal voltage sequence and all per-step iteration counts exactly
match the archived original audit (49 steps, 239 updates); its current PLT
is byte-identical to the original. Across the 31 exact common output points,
off/on maximum nonzero relative drain-current difference is 2.8356e-10, with
median 3.5717e-14. Both give Id(40 V)=2.24833171719662e-4 A/um at the stored
output precision.

### Conservation, input control and evidence

KCL uses all four finite terminal TotalCurrent values and a compensated sum.
Across every nonzero raw row, KCL/max-terminal is at most 6.5763e-15 in the
fixed pair and 4.3694e-15 in the adaptive pair; maximum absolute KCL is below
9.53e-19 A/um. All four zero-bias rows are identical: absolute KCL
3.2840e-34 A/um and relative KCL 8.9868e-16. No zero-bias exemption or
denominator floor was added. These are terminal-current checks, not proof of
identical local residual norms between SDevice and Vela.

All controls load the same saved Vg=4 V, Vd=0 state, mesh and parameters.
Shared input SHA-256 hashes match between original and revised runs. Deck
on/off pairs differ only by local Extrapolate lines. Digits=6, Iterations=25,
NotDamped=100, ErrRef(Electron/Hole)=1e8 and the original D5 physics are retained.
Vela predictor and IALMob remain off, and Vela's absolute residual gates are
unchanged. For SDevice, the accepted solves stop at update error < 1, which
is not interchangeable with Vela's absolute block-residual gates.

After explicit user authorization, the original 40 KiB input bundle was
uploaded and its hash verified on the VM. The first fixed-on run completed,
but MinStep=1 triggered 237 automatic NewtonPlot TDR dumps. Its following
fixed-off run was deliberately stopped; this is not a convergence failure.
All four final controls ran in `quiet/` with AutoCNPMinStepFactor=0 and
AutoNPMinStepFactor=0 applied equally. The controller's status path was also
made absolute. The quiet fixed-on has identical printed numeric Newton fields
and a byte-identical current PLT versus the original fixed-on, while removing
the 237 dumps (process wall time 1049.32 -> 393.32 s). Original evidence remains
under `phase3_initial_results/`; final evidence is under `phase3_results/`.

Remote final directory:
`sentaurus:/root/sentaurus_runs/vela_oracle_2022/templates_ldmos/numerics_validation_20260906/quiet/`.
Original bundle SHA-256:
`d41f1c173dc84567b5fe7fef680914afeee2d63dd69e499012da69b4b235d116`.
Final bundle SHA-256:
`5a7093fea82d1e90cf3477f810bd5aa197412054b10d3d0623cd8adf44ee9e10`.
Retrieved final archive SHA-256:
`ca45068404dc43a2f3bede1e726b3e58216587cd5358d93137d57b151ca8d765`.

Evidence: [phase3_metrics.json](../../reference_staging/templates_ldmos_numerics_validation_20260906/phase3_metrics.json),
[input manifest](../../reference_staging/templates_ldmos_numerics_validation_20260906/phase3_quiet_manifest.json),
[analysis script](../../reference_staging/templates_ldmos_numerics_validation_20260906/analyze_sdevice_ab.py),
and per-case `phase3_results/*/steps.csv` and `currents.csv`.

![SDevice forward Newton counts and adaptive steps](../../reference_staging/templates_ldmos_numerics_validation_20260906/sdevice_extrapolate_comparison.png)

## 4. Linear solve, floating-point resolution and runtime

Opt-in local-update CSVs now include the raw system's normwise and componentwise
linear backward errors, plus the worst component's row and absolute terms.
They are evaluated after undoing row weighting/equilibration and before Newton
step limiting. The ordinary solve and convergence policy are unchanged.

The four complete profiling replays reproduce all numeric trace/curve fields
exactly. Accepted and rejected state files and iteration traces are
byte-identical to the unchanged-binary controls. Only explicit output-file paths
differ in the failed attempt's metadata. Later 12-iteration row-localization
probes also reproduce their corresponding trace prefixes exactly; their
iteration-budget failures are diagnostic, not additional convergence failures.

The maximum normwise backward error across the four full runs is 3.66e-16.
Maximum componentwise relative errors are 0.26 to 0.73, so the global norm
alone is insufficient. Localization finds the maxima in electron rows 80 or
127, with absolute linear residuals from about 1.7e-24 to 8.2e-20. These maxima
do not coincide with the 5569 residual plateau and are much smaller in absolute
value than the blocking electron residual. A conditioning/forward-error audit
has not been performed, so this is not a universal linear-solver accuracy claim.

At node 5569 in the failed lagged run, iterations 120--154 request raw physical
QF updates of about 2.36e-19 V. Their line-search-selected updates are zero,
while the local linear closure error is at most a few e-27. The unscaled
increment is about 0.0025 V; one observed update quantum in the internally
scaled coordinates corresponds to about 3.59e-19 V. Damping by 0.5 or less
therefore rounds this correction away. The node residual remains
5.087928087335752e-12 while the combined electron block sits just above 1e-11.
This is direct evidence of a state-update resolution problem in this tail.
The live path reaches the unchanged block gate before being trapped there.

Timing from the full instrumented tail-pair replay:

| Stage | Lagged, 85 updates | Live, 13 updates |
| --- | ---: | ---: |
| Jacobian assembly | 11.38 s | 2.33 s |
| Linear factorization | 19.29 s | 2.85 s |
| Linear back substitution | 0.462 s | 0.068 s |
| Residual evaluations | 2.71 s | 0.642 s |
| Local-update diagnostics | 0.146 s | 0.023 s |
| Record final point, diagnostics and output | 0.579 s | 0.465 s |

Profiler scopes overlap; residual time must not be added again to a parent
line-search time. Final-point recording includes computation and output rather
than pure disk time. These measurements compare the two local runs; they are
not a hardware-normalized comparison with the VM. Full derivatives cost more
per assembly but avoid most repeated assemblies and factorizations.

## Verification and remaining work

- Release configuration and affected targets built successfully.
- Newton solver: 110 cases, 1444 assertions passed.
- DC sweep: 98 cases, 3479 assertions passed.
- Material-local Poisson charge: 5 cases, 14 assertions passed.
- Final row-diagnostic addition: 24 diagnostic cases, 449 assertions passed,
  and four exact trace-prefix replays passed.
- Four full unchanged-binary/instrumented comparisons passed numeric/state
  equivalence. No tolerance or production mobility setting was changed.

- Four final SDevice controls: exit code 0, full 40 V endpoints, zero rejected
  steps, expected actual extrapolation behavior, input/deck hashes verified.
- Log parser cross-checked against the original 49-step audit; adaptive-on
  reproduces its full per-step counts and byte-identical 31-point current PLT.
- Quiet/original fixed-on trace-field and current-PLT equivalence passed.
- Analysis assertions, Python syntax checks, figure inspection and local
  report-link checks passed. No solver rebuild was needed for this final
  SDevice execution and reporting step.

All four requested validation items are completed. The live-Jacobian step
sensitivity and update-resolution floor remain follow-up engineering findings.
The production Vg=4 31-point sweep has not been resumed or qualified by this diagnostic task;
Vg=8 and IALMob remain gated.
