# Templates/LDMOS G3 WP3 T1 residual-scaling audit

## Decision

The fixed-state replay completed for the exact eight same-contract G3 Id-Vg
states at `Vd=0.1 V`.  T1 cannot lock H1, H2, or H3 by itself: the
pre-registered collinearity guard is exceeded and the density/current/gradient
single-predictor confidence intervals are not pairwise distinguishable.  The
fits below are therefore descriptive only.

The frozen-edge-kernel sensitivity audit also does not support strong
Bernoulli saturation as the explanation for the node-3721 blind direction.
Across the 24 evaluable silicon transport edges in the frozen 40-edge hotspot
set, no edge reaches `|eta| >= 10`.  At `Vg=1.1666667 V`, the largest `|eta|`
is `1.64609`, the smallest downwind/upwind endpoint sensitivity ratio is
`0.96357`, and the node-3721 minimum is `0.981037`.

## Frozen replay contract

- Gate biases: `0, 1/6, 1/3, 1/2, 2/3, 5/6, 1, 7/6 V`; drain bias:
  `0.1 V`.
- Physics: Fermi statistics, OldSlotboom BGN, SRH/Auger, GradQF
  `edge_projection` HFS with the qualified contact ElectricField fallback,
  external AverageBox carrier couples, and `legacy_node_local` contact
  reconstruction.
- IALMob, predictor, QP, and impact ionization remain disabled.  No Id-Vd,
  breakdown, or full-physics state is consumed, and no production default is
  changed.
- Each state is replayed without a nonlinear solve through both
  `newton_carrier_term_probe` and `sg_edge_flux_probe`.
- Residuals are converted to `A/um` with the per-state ratio between SG particle
  line flux and scaled electron flux.  The maximum scale spread is
  `2.22e-16`; the reconstructed-density spread is zero.
- All eight state files contain 10,241 ordered node IDs and have the same node-ID
  SHA-256 as the mesh.

## Reproduction

From the worktree root, with the MSYS2 UCRT64 Python first on `PATH`, run:

```powershell
D:\msys64\ucrt64\bin\python.exe scripts\audit_templates_ldmos_g3_residual_scaling.py `
  --runner build-release\vela_example_runner.exe `
  --baseline-config reference_staging\templates_ldmos_g3_averagebox_state_feedback_vg1p166667_20260831\vg_1p166667\SSS\sg_probe.json `
  --reference-curve reference_staging\templates_ldmos_sentaurus2022\phase01_original_20260826_02\stage1_v4\phase23_t2022_contract_v2\sentaurus_ablation_report\normalized\G3-no-IALMob\IdVg_n2_des_drain_curve.csv `
  --mesh reference_staging\templates_ldmos_sentaurus2022\phase01_original_20260826_02\stage1_v4\vela_exact_topology\mesh.json `
  --state 0=reference_staging\templates_ldmos_sentaurus2022\phase01_original_20260826_02\stage1_v4\qualification\sentaurus_idvg_vg0_vd0p1_state.csv `
  --state 0.166666666666667=reference_staging\templates_ldmos_g3_state_feedback_averagebox_node_local_20260831\vg_0p166667\sentaurus_state.csv `
  --state 0.333333333333333=reference_staging\templates_ldmos_g3_state_feedback_averagebox_node_local_20260831\vg_0p333333\sentaurus_state.csv `
  --state 0.5=reference_staging\templates_ldmos_g3_state_feedback_averagebox_node_local_20260831\vg_0p500000\sentaurus_state.csv `
  --state 0.666666666666667=reference_staging\templates_ldmos_g3_state_feedback_averagebox_node_local_20260831\vg_0p666667\sentaurus_state.csv `
  --state 0.833333333333333=reference_staging\templates_ldmos_g3_contact_hfs_coefficient_ab_20260831\vg_0p833333\sentaurus_state.csv `
  --state 1=reference_staging\templates_ldmos_g3_averagebox_state_feedback_vg1_20260831\vg_1p000000\sentaurus_state.csv `
  --state 1.16666666666667=reference_staging\templates_ldmos_g3_averagebox_state_feedback_vg1p166667_20260831\vg_1p166667\sentaurus_state.csv `
  --output-dir reference_staging\templates_ldmos_g3_wp3_t1_residual_scaling_20260901
```

The script rejects any set other than the exact eight biases and verifies every
input/config contract before launching the probes.

## Residual anchors

| Vg (V) | seven-node L2 (A/um) | interface-band L2 (A/um) |
| ---: | ---: | ---: |
| 0 | 7.16542e-18 | 2.43247e-13 |
| 0.166667 | 1.09765e-17 | 1.90142e-13 |
| 0.333333 | 2.82812e-14 | 9.63003e-13 |
| 0.5 | 9.59348e-15 | 1.73434e-12 |
| 0.666667 | 5.68362e-13 | 6.67655e-11 |
| 0.833333 | 1.74628e-11 | 2.52745e-9 |
| 1 | 1.24146950452e-9 | 3.34151e-8 |
| 1.166667 | 6.97873112028e-9 | 8.50925e-8 |

The last two seven-node values reproduce the frozen-plan anchors.

## Scaling, condition, VIF, and LOOCV

Both responses use the predictors `log(n_interface)`, `log(Id)`, `Vg`, and
`log(|grad phin|)`.  The standardized design condition number is `94.9172`.
The VIFs are `484.78`, `1428.13`, `102.91`, and `533.81`, respectively, all
well above the pre-registered limit of 10.

| response / single predictor | exponent | 95% CI | R2 | LOOCV log-RMSE |
| --- | ---: | ---: | ---: | ---: |
| seven-node / `log(n_interface)` | 2.94581 | [1.44361, 4.44801] | 0.7933 | 6.1401 |
| seven-node / `log(Id)` | 0.976249 | [0.758235, 1.19426] | 0.9524 | 2.0497 |
| seven-node / `log(|grad phin|)` | 1.50029 | [0.854112, 2.14647] | 0.8432 | 3.8465 |
| interface band / `log(n_interface)` | 1.81282 | [0.436568, 3.18906] | 0.6339 | 5.3944 |
| interface band / `log(Id)` | 0.662291 | [0.473749, 0.850833] | 0.9249 | 1.9411 |
| interface band / `log(|grad phin|)` | 1.11076 | [0.934117, 1.28740] | 0.9753 | 1.0531 |

The full standardized model is unstable under leave-one-out validation:
seven-node and interface-band LOOCV log-RMSE values are `19.3394` and
`7.57932`.  High in-sample R2 must not be interpreted causally.

## Hotspot-edge sensitivity scope

The output preserves all 40 unique topology edges incident on the frozen
hotspot nodes for every state (320 rows total); 24 per state are evaluable
silicon transport edges and 16 are inactive/nontransport rows.  Each row gives
both analytic endpoint `phin` sensitivities, Bernoulli eta, upwind/downwind
roles, physical edge flux, and saturation classification.

The derivative is an analytic frozen-edge-kernel diagnostic.  Endpoint Fermi
compressibility is included, while the already assembled HFS mobility and
generalized-Einstein secant factor are held fixed.  It is sufficient to test
Bernoulli endpoint suppression, but it is not the complete production
Jacobian.

## Artifacts

- `reference_staging/templates_ldmos_g3_wp3_t1_residual_scaling_20260901/summary.json`
- `reference_staging/templates_ldmos_g3_wp3_t1_residual_scaling_20260901/residual_scaling_states.csv`
- `reference_staging/templates_ldmos_g3_wp3_t1_residual_scaling_20260901/hotspot_edge_phin_sensitivity.csv`
- `reference_staging/templates_ldmos_g3_wp3_t1_residual_scaling_20260901/state_manifest.json`
- `reference_staging/templates_ldmos_g3_wp3_t1_residual_scaling_20260901/report.md`
- Per-bias copied configs, carrier/SG probe CSVs, and logs under the same output
  directory.

## Verification

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest tests.regression.test_audit_templates_ldmos_g3_residual_scaling -v
D:\msys64\ucrt64\bin\python.exe -m py_compile scripts\audit_templates_ldmos_g3_residual_scaling.py tests\regression\test_audit_templates_ldmos_g3_residual_scaling.py
```

The three focused tests cover exact eight-state contract rejection, analytic
saturated-edge endpoint classification, and downgrade of collinear fits to
descriptive evidence.
