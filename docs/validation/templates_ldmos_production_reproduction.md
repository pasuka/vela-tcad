# LDMOS production electrothermal reproduction

This is the run guide for the frozen D0 production configurations. R9 contact
consistency is now the recommended explicit profile; R7 remains the compatible
exporter default and rollback baseline. See the R9 section below.
The [current status](templates_ldmos_current_status.md), profile-specific R9
evidence below and [frozen R7 evidence](templates_ldmos_production_electrothermal_2026-09-14.md#r7完成后的暂停记录)
define their separate numerical and performance qualifications. The ordinary `dc_sweep`
entry remains isothermal.

## Inputs and provenance

The versioned profile
[`d0_production_r7.json`](../../reference_tcad/templates_ldmos_sentaurus2022/thermal/d0_production_r7.json)
freezes the two prepared SI inputs, their production decks, mesh and IALMob
geometry by SHA256. These large inputs remain external evidence, not committed
simulation outputs. The profile also records the qualified Linux runner hash
and the complete archive hash.

The existing extracted evidence is at
`reference_staging/templates_ldmos_joint_20260914/linux_environment/copied_vm`.
Alternatively extract the matching `r7_final_evidence.tgz` into a new evidence
directory; preserve `cases_r7/` and `cases/data/`. The exporter checks all six
required files before writing anything. The full archive includes Linux shared
libraries, but these are not needed when exporting inputs for Windows.

[`export_templates_ldmos_production.py`](../../scripts/export_templates_ldmos_production.py)
only rewrites the declared mesh/geometry file paths and deck input/output paths.
It preserves the material model, contacts, native geometry coefficients,
temperature semantics, gates, predictor, iteration budgets and all 31 exact
bias values. Both gates initialize independently from neutral 300 K, followed
by prebias to the requested gate voltage. Exporting does not start a solver.

## Windows Release commands

Run from the current checkout/worktree, using an unused output directory:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --preset windows-ucrt64-release
cmake --build --preset windows-ucrt64-release --parallel 2 --target vela_example_runner
python -X utf8 scripts/export_templates_ldmos_production.py --evidence-root reference_staging/templates_ldmos_joint_20260914/linux_environment/copied_vm --output build-release/reference_tcad/ldmos-r7
build-release/vela_example_runner.exe --config build-release/reference_tcad/ldmos-r7/vg4.json
build-release/vela_example_runner.exe --config build-release/reference_tcad/ldmos-r7/vg8.json
```

Check that configuration actually reports UMFPACK enabled. The exported inputs
explicitly select `electrothermal_linear_solver: "umfpack"`; missing support
fails instead of silently selecting another backend. Debug or a different
backend does not inherit the frozen Release performance qualification.

The generated `manifest.json` records input and dependency hashes, replacements,
and every exported file hash. Dependencies use absolute paths so the runner can
start from another working directory. Regenerate the bundle after moving it.
Existing output directories are never overwritten by the exporter. Each
generated deck writes to its own `results_vg4` or `results_vg8` directory.

## Pause, resume and outputs

Create a `STOP` file inside an active result directory to pause between attempts.
Alternatively set a positive `pause_after_attempts` in the deck. A controlled
pause returns exit code 1 with `ledger.json` status `stopped_at_checkpoint`;
inspect the ledger to distinguish it from failure. For continuation remove
`STOP`, set `resume: true`, and remove/reset `pause_after_attempts`. Input, mesh,
sweep and initialization must still match the checkpoint. Interrupted
initialization requires a fresh result directory.

Each attempt retains full input, state, residuals, current and heat outputs.
`curve.csv` uses V, A/µm and K; `terminal_balance.csv` retains all four contacts.
The ledger includes accepted/rejected attempts, initialization, timings and
exact-point summaries. External geometry and executable identity are frozen by
the evidence/export manifests; the C++ checkpoint itself does not verify every
external file or executable hash. Keep the bundle and runtime dependencies fixed
when resuming a qualified run.

## Qualification and remaining limits

R7 Linux completed both 31-point curves, original electrical/thermal gates,
the approved 2% carrier-density RMS and 30 meV band-edge gates, and real
pause/resume plus production/point-service agreement. Same-VM serial external
wall times were 626.65/516.07 s versus native 610.06/493.71 s (1.027/1.045×).
Those values describe the first pair per gate. The subsequent
[R7 repeatability study](templates_ldmos_r7_stability_newton_2026-09-14.md)
completed three pairs per gate: checked numerical fields were exactly equal
across repeats and all six paired wall ratios were below 1.5. Absolute Vela
wall time still varied substantially (sample CV 26.4%/30.2%); CPU time and
Newton counts remain higher than native. Compare each new candidate with its
own paired native run, and retain the distinction between numerical
repeatability and variable wall time.

The export preserves that configuration, but a changed program, machine or
model requires its own qualification. Shared R6 physics checks include G3 and
D5/D4 reclosure and the 300 K material/transport limit; no new complete cold
device curve is claimed for R7. Original Solve did not activate the optional
quantum, Thermodynamic, Peltier or RecGenHeat extensions; they remain outside
this qualification.

Exporter regressions use synthetic provenance/path fixtures. The real R7 input
was additionally exported and run through frozen Windows Release neutral
initialization and the zero-bias point; this checks the generated file paths
and service entry, not a new Windows full-curve performance qualification.

## Recommended explicit R9 contact-consistency profile

[`d0_production_r9_contact.json`](../../reference_tcad/templates_ldmos_sentaurus2022/thermal/d0_production_r9_contact.json)
freezes the actual qualified contact-repair inputs and decks, without generating
new numerical overrides. It is recommended following the
[62-point joint and repeated timing qualification](templates_ldmos_contact_sweep_2026-09-16.md).
The ordinary exporter default remains R7 for existing callers; select R9 explicitly.
The underlying C++ contact-repair default remains false for other devices/configurations.

The R9 archive contains the full matrix, not another copy of the unchanged
mesh/IALMob dependencies. Extract `full_final_evidence.tgz` beneath
`extrapolation_20260915/contact_sweep_20260916` in the evidence root and add
the hash-checked `cases/data` dependencies from the R7 archive. Both archive
hashes and all six required input/dependency hashes are recorded in the profile.
The existing copied evidence already has this layout:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts/export_templates_ldmos_production.py --profile reference_tcad/templates_ldmos_sentaurus2022/thermal/d0_production_r9_contact.json --evidence-root reference_staging/templates_ldmos_contact_sweep_20260916/copied_vm --output build-release/reference_tcad/ldmos-r9
```

The qualified Linux runner hash is in the profile. A rebuilt Windows binary or
an additional optimization does not inherit that Linux timing qualification.
R9 protects initialization/gate prebias, passes actual pause/resume, and preserves
the original electrical/thermal/2% density/30 meV gates. Two complete numerical
trajectories repeat exactly. Drain updates are 344/302; median end-to-end times
are 520.71/511.32 s, or 1.018/1.063 times paired native medians. CPU is still
1.625/1.763 times native; per-round speedup varies. The recommendation rests on
qualified behavior and repeatable removal of contact closure work, not a promise
of fixed wall-time acceleration. Root caching passed separate exact-trajectory controls but has no demonstrated
stable speedup; local QF limiting failed screening. Neither experiment is
enabled by this profile. See the [cost study](templates_ldmos_newton_cost_2026-09-16.md).

## Optional R8 density-coordinate profile

The explicit
[`d0_production_r8_density_guard.json`](../../reference_tcad/templates_ldmos_sentaurus2022/thermal/d0_production_r8_density_guard.json)
freezes the independently verified R8 candidate. It enables density-coordinate
Newton updates only at drain bias ≤1 V when an accepted-history prediction is
actually used. Initialization and gate prebias use the original QF updates.
The ordinary export command still defaults to R7.

Use the verified `r8_initguard_final_20260915.tgz` archive (SHA256 in the profile),
extracted with its `r8_density_low_20260915/` and `cases/data/` directories:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python -X utf8 scripts/export_templates_ldmos_production.py --profile reference_tcad/templates_ldmos_sentaurus2022/thermal/d0_production_r8_density_guard.json --evidence-root reference_staging/templates_ldmos_joint_20260914/linux_environment/copied_r8_initguard --output build-release/reference_tcad/ldmos-r8-density-guard
```

This only prepares the two decks. Run them with the corresponding Release
runner after checking its actual UMFPACK support; exporting or rebuilding on
another platform does not transfer the frozen Linux qualification.

R8 completed all 62 electrical, thermal and local-field gates. Drain updates
fell from 366/321 to 354/314, and line-search candidates from 496/446 to 448/401.
Its one same-VM timing pair per gate measured Vela 965.50/1050.28 s against
native 920.03/745.56 s (1.049/1.409×), both within 1.5. Vg8 gained an extra
rejected attempt after the altered low-bias path caused earlier step growth.
This does not establish stable wall-time acceleration or transfer R7's
repeatability and actual pause/resume qualification. Details and remaining
predictor work are in the
[R8 results](templates_ldmos_r7_stability_newton_2026-09-14.md#r8完整曲线计时与联合验收结果).
