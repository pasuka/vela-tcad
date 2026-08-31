# Templates/LDMOS G3 node-local + AverageBox 31-point qualification (2026-08-31)

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

## Decision

1. The node-local source-short and carrier-only AverageBox contracts are
   confirmed to close the former subthreshold amplitude failure. Five of six
   frozen curve gates pass.
2. Stage-3 remains **failed** solely on the frozen 20% maximum-gm gate. The
   gate is not relaxed, and no mobility, flatband or threshold calibration is
   authorized.
3. The known-difference ledger remains draft. Before IALMob, the next minimal
   experiment is a fixed-state coefficient/mobility replay at the exact
   `1.0 V` and `1.1667 V` states to separate AverageBox transport support from
   the remaining HFS/state-feedback slope difference.
4. A full per-region Poisson AverageBox implementation remains deferred: the
   present run does not provide evidence that it is the gm root cause.

## Reproducibility artifacts

- Curve runner: `scripts/run_templates_ldmos_g3_averagebox_curve.py`.
- Exact-point scorer: `scripts/analyze_templates_ldmos_g3_idvg.py`.
- Frozen profile contract:
  `reference_tcad/templates_ldmos_sentaurus2022/contracts/diagnostics/templates_ldmos_external_averagebox_profile.json`.
- Ignored run directory:
  `reference_staging/templates_ldmos_g3_averagebox_curve_node_local_20260831/`.
- Run contract: `run_summary.json`; exact score: `qualification.json` and
  `qualification.md`; curve: `g3_averagebox_idvg.csv`; numerical evidence:
  `newton_history.csv`, `newton_iterations.csv`, `terminal_balance.csv` and
  `srh_balance.csv`.
