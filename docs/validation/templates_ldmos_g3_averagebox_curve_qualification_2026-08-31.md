# Templates/LDMOS G3 node-local + AverageBox 31-point qualification (2026-08-31)

> Historical result: the remaining maximum-gm failure documented here was
> closed on 2026-09-02 by the independently qualified material-local Poisson
> charge-volume policy. See
> `templates_ldmos_g3_material_local_charge_volume_qualification_2026-09-02.md`.

## Scope and frozen controls

This run executes the exact 31 Sentaurus CurrentPlot gate biases at
`Vd=0.1 V`. It retains the qualified interior GradQF plus contact-cell
ElectricField HFS support and the WP1.5 cached-mobility Jacobian repair. The
only retained contract additions relative to the previous production curve are:

- `contact_boundary_reconstruction=legacy_node_local`, required because the
  source metal physically shorts p+ body pickup and n+ source nodes;
- the explicit, default-off `templates_ldmos_external_averagebox`
  carrier-transport couple profile with 16,237 archived edge coefficients.

Node volumes remain barycentric. Predictor, IALMob, quantum potential, impact
ionization and every acceptance threshold remain disabled or unchanged. No
unscored warm-start points or curve interpolation are used.

## Exact-point result

All 31 points converge. The run uses 177 total Newton iterations, 5--7 per
point. The worst final block residuals are `1.143e-9` psi, `9.851e-12`
electron continuity and `3.378e-17` hole continuity. The worst resolved KCL
relative error is `1.143e-9` at `Vg=0 V`.

| Metric | Result | Frozen limit | Status |
|---|---:|---:|:---:|
| Exact resolved points | 31/31 | 31/31 | pass |
| Median absolute log-current error | 0.069890 dex | 0.10 dex | pass |
| P95 absolute log-current error | 0.099697 dex | 0.20 dex | pass |
| Maximum absolute log-current error | 0.106482 dex | report only | pass |
| Strong-inversion endpoint relative error | 13.327% | 20% | pass |
| Diagnostic 1e-8 A/um Vth error | 0.0387 mV | 100 mV | pass |
| Maximum-gm relative error | 24.308% | 20% | **fail** |
| Worst resolved KCL relative error | 1.143e-9 | 1% | pass |

The P95 curve error improves from `0.432753 dex` to `0.099697 dex`, and the
Vg=0.5 V error remains closed at about `0.0130 dex`. The current ratio is near
unity through `Vg=0.8333 V`, then becomes `0.815` at `1.0 V`, `0.783` at
`1.1667 V`, and gradually recovers to `0.867` at `5 V`.

Both curves place maximum gm on the same exact segment,
`Vg=1.0--1.1667 V` (midpoint `1.0833 V`). The Sentaurus and Vela slopes are
`9.44257e-6` and `7.14731e-6 A/(um V)`, respectively. The gm failure is
therefore not caused by peak-position selection, curve interpolation or an
unresolved bias point.

## Contact-HFS frozen replay and state feedback

The contact-HFS decision was rerun as a true same-contract pair. Each point
uses the same frozen Sentaurus state, node-local contact contract, and external
AverageBox transport couples; the only switch is the contact-node-cell
ElectricField fallback. This removes the previous comparison's mesh-default
coefficient confounder.

| Vg (V) | Contact HFS off / Sentaurus | Contact HFS on / Sentaurus | Error improvement (dex) |
|---:|---:|---:|---:|
| 0.166667 | 24.66599 | 1.74731 | 1.14973 |
| 0.500000 | 11.00098 | 1.00071 | 1.04112 |
| 0.833333 | 10.98829 | 1.00020 | 1.04084 |

The median frozen-state error falls from `1.04143 dex` to `0.000308 dex`.
Contact HFS therefore remains a necessary operator correction; it was not the
cause of the former self-consistent `2.6--2.7x` current platform.

The eight-combination `psi/phin/phip` replay was then repeated at
`Vg=1/6, 1/3, 1/2, 2/3 V` using the new 31-point states and the same corrected
operator:

| Vg (V) | VVV/Sentaurus | SSS/Sentaurus | Feedback (dex) | Largest state-family effect |
|---:|---:|---:|---:|:---:|
| 0.166667 | 0.97364 | 1.74731 | -0.253971 | phin (cancels operator excess) |
| 0.333333 | 1.00860 | 1.00356 | 0.002178 | psi |
| 0.500000 | 1.03047 | 1.00071 | 0.012727 | phin |
| 0.666667 | 1.02125 | 1.00018 | 0.009052 | phin |

The median feedback is now only `0.005615 dex` (`1.0130x`), versus the old
`0.403147 dex` (`2.530x`). At `Vg=1/6 V`, a signed `phin` feedback term cancels
most of the remaining fixed-state operator excess, so reporting only unsigned
"error recovery" would be misleading. The old conclusion that `phin`
self-consistency generated the curve-wide platform is withdrawn for the
corrected contract.

## Decision

1. The node-local source-short and carrier-only AverageBox contracts are
   confirmed to close the former subthreshold amplitude failure. Five of six
   frozen curve gates pass.
2. Stage-3 remains **failed** solely on the frozen 20% maximum-gm gate. The
   gate is not relaxed, and no mobility, flatband or threshold calibration is
   authorized.
3. The follow-up fixed-state coefficient/mobility replay at the exact `1.0 V`
   and `1.1667 V` states closes the operator slope to `+0.0202%`, while the
   self-consistent slope remains `-24.3076%`. The remaining failure is
   electron-QF state feedback, not AverageBox or HFS operator amplitude. See
   `templates_ldmos_g3_gm_segment_fixed_state_audit_2026-08-31.md`.
4. A full per-region Poisson AverageBox implementation remains deferred: the
   present run does not provide evidence that it is the gm root cause.

## Reproducibility artifacts

- Curve runner: `scripts/run_templates_ldmos_g3_averagebox_curve.py`.
- Exact-point scorer: `scripts/analyze_templates_ldmos_g3_idvg.py`.
- Same-contract contact-HFS replay:
  `scripts/run_templates_ldmos_g3_contact_hfs_replay.py`.
- Corrected state-feedback replay:
  `scripts/audit_templates_ldmos_g3_state_feedback.py`.
- Maximum-gm segment replay report:
  `docs/validation/templates_ldmos_g3_gm_segment_fixed_state_audit_2026-08-31.md`.
- Frozen profile contract:
  `reference_tcad/templates_ldmos_sentaurus2022/contracts/diagnostics/templates_ldmos_external_averagebox_profile.json`.
- Ignored run directory:
  `reference_staging/templates_ldmos_g3_averagebox_curve_node_local_20260831/`.
- Run contract: `run_summary.json`; exact score: `qualification.json` and
  `qualification.md`; curve: `g3_averagebox_idvg.csv`; numerical evidence:
  `newton_history.csv`, `newton_iterations.csv`, `terminal_balance.csv` and
  `srh_balance.csv`.
- Ignored replay directories:
  `reference_staging/templates_ldmos_g3_contact_hfs_replay_averagebox_node_local_v2_20260831/`
  and
  `reference_staging/templates_ldmos_g3_state_feedback_averagebox_node_local_20260831/`.
- Ignored maximum-gm endpoint evidence:
  `reference_staging/templates_ldmos_g3_averagebox_state_feedback_vg1_20260831/`,
  `reference_staging/templates_ldmos_g3_averagebox_state_feedback_vg1p166667_20260831/`
  and `reference_staging/templates_ldmos_g3_gm_state_capture_20260831/`.
