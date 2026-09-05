# Templates/LDMOS Stage-4 verification handoff

Date: 2026-09-05  
Worktree: `D:/code-repo/vela-tcad/.worktrees/templates-ldmos-phase-a`  
Branch: `codex/templates-ldmos-phase-a`

## 1. Handoff outcome

Phase-A Stage 4 is still **in progress**. The historical Vg=4 V, Vd=10 mV
solver blocker and the later 0.1625 V quasi-Fermi representation floor are
closed. After merging local `main`, the strict no-predictor D5 path was resumed
from the sealed 0.8 V checkpoint and reached the first nonzero common
Sentaurus/Vela reference point at Vd=1.33333333333333 V.

At that point Vela is 4.3243% above Sentaurus and the terminal KCL ratio is
5.12e-14. This individual point is inside the final current and KCL limits,
but the 31-point Vg=4 V curve is not complete and must not be declared
qualified. Vg=8 V and IALMob remain gated on completion of the Vg=4 V D5
curve.

## 2. Repository state

The implementation was committed before integration:

- `832d68e` -- `Stabilize LDMOS drain continuation`.

The latest local `main` was then merged:

- local `main`: `8e99e57` -- `Promote validated TCAD fixtures and remove examples`;
- integration commit: `15a013e` -- `Merge main into LDMOS phase A`;
- `15a013e` parents: `832d68e` and `8e99e57`.

The merge deleted the obsolete example-dependent `test_ldmos_dd` and
`test_mos_mixed_material` targets with `main`. Branch-specific Poisson charge
volume and non-finite-carrier assertions were retained in the focused,
synthetic-mesh `test_material_local_poisson_charge` target.

Large simulation products stay below ignored `reference_staging/` and are not
committed. The tracked handoff state is recorded in
`reference_tcad/templates_ldmos_sentaurus2022/stage4_vd10mv_newton_audit_summary.json`.

## 3. Frozen D5 physical and numerical contract

Do not change these settings while completing the current curve:

| Contract item | Frozen value |
| --- | --- |
| Gate bias | 4 V |
| Carrier equations | Poisson + electron + hole |
| Temperature | isothermal 300 K |
| Carrier statistics/BGN | Fermi + OldSlotboom |
| Recombination | SRH + Auger |
| High-field mobility | live HFS residual, GradQF interior with qualified contact fallback |
| Mobility Jacobian | lagged field, `jacobian_field_derivatives=false` |
| Carrier couples | `templates_ldmos_external_averagebox` |
| Contact reconstruction | `legacy_node_local` |
| Poisson charge volume | material-local barycentric |
| Predictor | disabled |
| IALMob/QP/avalanche/self-heating | disabled |
| Warm-start coordinate | `quasi_fermi_recenter_on_initial_state=true` |
| Internal continuation | initial step 2.5 mV, minimum 1 mV, growth factor 1, 12 retries |
| State persistence | rolling state plus every accepted internal step |

The QF recenter changes the internal coordinate origin only. Absolute
quasi-Fermi potentials, densities and residual equations remain unchanged.
The lagged mobility Jacobian is a qualified quasi-Newton strategy; the
residual continues to use live HFS mobility.

## 4. Post-merge build and regression evidence

The Release tree was reconfigured with UCRT64 Ninja and built successfully.

| Test | Result |
| --- | ---: |
| `test_material_local_poisson_charge` | 5 cases, 14 assertions, pass |
| `test_dc_sweep` | 98 cases, 3479 assertions, pass |
| `test_newton_solver` | 110 cases, 1438 assertions, pass |
| `tests.regression.test_templates_ldmos_phase23` | 25 tests, pass |

The first new test should be retained after future `main` integrations because
it replaces the relevant assertions formerly embedded in deleted example
tests.

## 5. Latest Vg=4 V continuation evidence

Three restartable intervals were run on integration commit `15a013e`:

| Interval | Accepted | Rejected | Endpoint Id (A/um) | Endpoint Newton |
| --- | ---: | ---: | ---: | ---: |
| 0.8 -> 0.9 V | 41 | 1 | 7.166383918583358e-5 | 18 |
| 0.9 -> 1.0 V | 41 | 0 | 7.818799291571467e-5 | 20 |
| 1.0 -> 1.33333333333333 V | 135 | 0 | 9.803614361917374e-5 | 24 |

The only rejected attempt was an isolated 5 mV proposal at 0.805 V. Its
automatic 2.5 mV fallback and every later transfer passed. Therefore 5 mV is
not qualified as the production step in this region.

At the exact common reference point:

| Metric | Value |
| --- | ---: |
| Drain bias | 1.33333333333333 V |
| Vela Id | 9.803614361917374e-5 A/um |
| Sentaurus Id | 9.3972532439272e-5 A/um |
| Vela/Sentaurus | 1.04324254198989 |
| Relative current error | +4.324254% |
| Absolute log error | 0.0183853 decade |
| Final Poisson block | 8.41392e-10 |
| Final electron block | 6.89020e-12 |
| Final hole block | 1.97369e-25 |
| Relative KCL | 5.11728e-14 |
| QF bound violations | 0 |

The single point passes the final nonzero-current and KCL magnitude limits.
It does not establish the curve median, P95, low-Vd differential resistance,
40 V endpoint or two-gate current-ratio gates.

## 6. Evidence locations and hashes

All paths below are relative to the worktree root.

| Evidence | Path | SHA-256 |
| --- | --- | --- |
| D5 base config | `reference_staging/templates_ldmos_stage4_d5_20260903/vela_run_r9_bounded40_no_guard/d5_idvd_vg4.json` | `fc8dcc8079980607b8c16b2dc55acdc53aca177790d996ea205c2104165557d9` |
| Starting 0.8 V state | `reference_staging/templates_ldmos_stage4_vg4_fixed_step_0p68_to_0p8_20260904/point_bias_0p800000.csv` | `88137ed5e45ff3ee3d0282c8732cdd49ee358ffb1b7ce37e806064be2006c6a7` |
| Latest 1.333333 V state | `reference_staging/templates_ldmos_stage4_vg4_post_main_1p0_to_1p333333_20260905/point_bias_1p333333.csv` | `18693c1c50d6fce75cdd5aaaf7dd9720d99533155f0d98564387603e30141bd9` |
| Latest endpoint curve | `reference_staging/templates_ldmos_stage4_vg4_post_main_1p0_to_1p333333_20260905/curve.csv` | `287a3aa85737cd6c484a89d68bccf2b12e066a02ab53c7e31998e3af49538bad` |
| Sentaurus Vg=4 D5 curve | `reference_staging/templates_ldmos_stage4_idvd_ablation_20260903/analysis_r2/normalized/D5-no-IALMob/IdVd_Vg1_n4_des_drain_curve.csv` | `3c9266b2616e2e779b853cee7bd94a0ee13478015d727c5bc2306ed8e152cb6a` |

The latest run directory also contains `terminal_balance.csv`,
`newton_attempts.csv`, `newton_iterations.csv`, `state.csv`, every accepted
checkpoint and `reclose.log`.

## 7. Reproduction and continuation commands

Initialize the required Windows toolchain:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
Set-Location "D:\code-repo\vela-tcad\.worktrees\templates-ldmos-phase-a"
```

Rebuild and repeat the focused regression set:

```powershell
cmake -S . -B build-release -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --parallel 8
build-release\test_material_local_poisson_charge.exe
build-release\test_dc_sweep.exe
build-release\test_newton_solver.exe
python -m unittest tests.regression.test_templates_ldmos_phase23
```

Resume from the latest state to the next exact Sentaurus point. Use a new
output directory because the runner rejects an existing destination:

```powershell
python scripts\run_templates_ldmos_stage4_vd10mv_frozen_mobility_reclose.py `
  --runner build-release\vela_example_runner.exe `
  --base-config reference_staging\templates_ldmos_stage4_d5_20260903\vela_run_r9_bounded40_no_guard\d5_idvd_vg4.json `
  --initial-state reference_staging\templates_ldmos_stage4_vg4_post_main_1p0_to_1p333333_20260905\point_bias_1p333333.csv `
  --gate-bias 4 `
  --mobility-jacobian frozen `
  --recenter-qf `
  --bias-points 1.33333333333333 2.66666666666667 `
  --output-dir reference_staging\templates_ldmos_stage4_vg4_1p333333_to_2p666667_YYYYMMDD `
  --min-step 0.001 `
  --max-retries 12 `
  --initial-step 0.0025 `
  --growth-factor 1
```

At the measured rate, a full 1.333333 V reference interval can take roughly
one to two hours. Do not use a larger production step merely to reduce wall
time. A larger-step experiment, if desired, must use a separate output
directory and be treated as a state/current equivalence A/B, not as an
implicit contract change.

## 8. Required next sequence

1. Continue Vg=4 V D5 through every exact Sentaurus drain-bias point, preserving
   accepted-step checkpoints and recording each point's current, Newton count,
   residual blocks and KCL.
2. Assemble one Vg=4 V 31-point curve on the exact reference lattice. Do not
   score interpolated internal points.
3. Run `scripts/analyze_templates_ldmos_stage4_d5.py` and evaluate median/P95
   nonzero-current error, low-Vd differential resistance, 40 V endpoint and
   KCL against the frozen Stage-4 gates.
4. If and only if Vg=4 passes, repeat the same contract for Vg=8 V and evaluate
   the two-gate current-ratio gate.
5. Only after D5 qualification resume IALMob development. Predictor, parameter
   calibration, QP, avalanche, self-heating and accepted-difference ledger
   promotion remain out of scope.

If a transfer fails, stop curve scoring at the last accepted checkpoint and
classify the failure from `newton_attempts.csv`, `newton_iterations.csv` and
the rejected-state bundle. Do not weaken the frozen block residual ceilings.

## 9. Known open issues

- Runtime remains high because the strict branch needs approximately 2.5 mV
  internal steps above 0.2 V.
- The near-null carrier mode seen in uncapped diagnostic solves remains an
  investigation item, although it has not caused a rejection from 0.9 through
  1.333333 V.
- The complete 31-point Vg=4 accuracy distribution is unknown; one passing
  point cannot establish median or P95.
- Vg=8 V has not been rerun after the `main` merge.
- IALMob is known to be material for the final template, but remains correctly
  disabled until D5 closure.

Canonical background and solver evidence are in:

- `docs/validation/templates_ldmos_stage4_idvd_execution_2026-09-03.md`;
- `docs/validation/templates_ldmos_stage4_vd10mv_newton_audit_2026-09-03.md`;
- `docs/superpowers/plans/2026-08-26-templates-ldmos-sentaurus-vela-validation-plan.md`.
