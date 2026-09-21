# Templates/LDMOS G3 contact-HFS qualification (2026-08-31)

> 2026-08-31 correction: the fixed-state contact-HFS replay has now been
> repeated with the same node-local + external-AverageBox operator on both
> sides of the HFS switch. At Vg=0.5/0.833333 V, HFS off gives
> 11.00098x/10.98829x and HFS on gives 1.00071x/1.00020x. The HFS conclusion
> remains valid, but the earlier mesh-default baseline must not be mixed with
> AverageBox when assigning effect size. See
> `templates_ldmos_g3_averagebox_curve_qualification_2026-08-31.md`.

## Scope and frozen controls

This qualification tests one physics change only: the G3 electron/hole
high-field mobility drive remains `GradQuasiFermi` in the device interior and
falls back to the Tri3 cell electrostatic-gradient magnitude for transport
cells touching a contact node. IALMob, predictor, impact ionization, and all
acceptance thresholds remain disabled or unchanged. The original 31
Sentaurus CurrentPlot biases are the only curve-scoring points.

## Three-point fixed-state decision gate

The archived Sentaurus states at Vg = 1/6, 1/2, and 5/6 V were replayed through
the production Vela SG operator before and after the Jacobian repair. The
repeated result is unchanged because the repair affects only Newton
linearization, not the residual operator.

| Vg (V) | Baseline Vela/Sentaurus | Contact-HFS Vela/Sentaurus | Error improvement (dex) |
|---:|---:|---:|---:|
| 0.166667 | 24.4976 | 1.79505 | 1.13505 |
| 0.500000 | 11.0772 | 1.04621 | 1.02481 |
| 0.833333 | 11.0647 | 1.04570 | 1.02453 |

The median magnitude error falls from 1.04443 dex to 0.01962 dex. This closes
the fixed-state contact-HFS coefficient/operator decision gate.

## WP1.5 node-702 root cause and repair

The first 400-iteration curve stopped at Vg = 2/3 V with electron continuity
residual 1.57e-10. Psi and hole continuity were already inside their frozen
ceilings, every block-filter step was accepted at damping 1, and the maximum
electron row stayed at interior node 702 next to two drain contact nodes.

The frozen final state was replayed through carrier-row, carrier-term,
edge-mobility, and directional-Jacobian probes. The decisive observations
were:

- node 702 electron residual: -1.52e-10, entirely from cancellation of edge
  fluxes whose absolute sum is 1.38e-5;
- 400/400 full Newton steps accepted, with no line-search rejection;
- the node-702 electron-QF JVP had analytic/finite-difference norms
  13.448/1.234 at a 1e-6 V direction;
- disabling HFS while preserving the state reduced maximum JVP relative error
  to 8.3e-8, excluding contact-basin storage and block-filter logic.

The residual mobility path treats `constant_field` as a high-field model, but
the Jacobian's cached mobility predicate omitted `constant_field`. The
Jacobian therefore differentiated edge fluxes with low-field mobility while
the residual used velocity-saturated mobility. Adding `constant_field` to the
cached high-field predicate reduces the node-702 electron diagonal from
-3.47e5 to -2.59e4 and makes the 1e-6 V QF JVP agree to about 1.1e-13. A
dedicated NewtonSolver regression test freezes this contract. No convergence
threshold, predictor, mobility calibration, or third-vertex experiment is
part of the retained repair.

## Self-consistent 31-point result

The repaired Release run used all 31 exact biases, no predictor, and the
frozen block absolute ceilings: psi 5e-8, electron continuity 1e-11, and hole
continuity 3e-10. All points converge in 186 total Newton iterations (5--8 per
point, median 6). The worst final block residuals are 1.16e-8, 8.82e-12, and
2.41e-11 respectively.

| Metric | Result | L2 limit | Status |
|---|---:|---:|:---:|
| Exact resolved points | 31/31 | 31/31 | pass |
| Median absolute log-current error | 0.007462 dex | 0.10 dex | pass |
| P95 absolute log-current error | 0.432753 dex | 0.20 dex | **fail** |
| Strong-inversion endpoint relative error | 1.513% | 20% | pass |
| Diagnostic 1e-8 A/um Vth error | 43.206 mV | 100 mV | pass |
| Maximum-gm relative error | 6.592% | 20% | pass |
| Worst resolved KCL relative error | 7.96e-6 | 1% | pass |

The maximum log-current error is 0.436285 dex at Vg = 1/3 V. The four points
from 1/6 through 2/3 V retain a nearly multiplicative 2.61--2.73x current
factor. Thus the Jacobian repair closes the numerical trajectory and KCL but
does not close the remaining low-current coefficient/state-feedback
difference.

## Decision

- The contact ElectricField fallback is qualified for fixed-state G3 SG replay.
- The WP1.5 `constant_field` HFS Jacobian repair is qualified; the 31-point
  no-predictor trajectory passes every frozen numerical convergence and KCL
  gate.
- Stage-3 L2 remains **failed** solely on the 0.20 dex P95 curve gate. IALMob
  remains disabled and no threshold shift or mobility calibration is allowed.
- The known-difference ledger remains draft. The remaining 2.6--2.7x
  subthreshold factor must be tested at the coefficient/discretization and
  self-consistent state-feedback levels before it can be classified as an
  engine floor.

## Reproducibility artifacts

- Node-702 audit: `scripts/audit_templates_ldmos_g3_contact_hfs_node702.py`.
- Fixed-state runner: `scripts/run_templates_ldmos_g3_contact_hfs_replay.py`.
- Curve runner: `scripts/run_templates_ldmos_g3_contact_hfs_curve.py`.
- Fixed-state evidence (ignored staging):
  `reference_staging/templates_ldmos_g3_contact_hfs_replay_jacobian_fix_20260831/summary.json`.
- Node-702 evidence (ignored staging):
  `reference_staging/templates_ldmos_g3_contact_hfs_node702_audit_fixed_20260831/summary.json`.
- Final curve, Newton history, KCL, and exact-point score (ignored staging):
  `reference_staging/templates_ldmos_g3_contact_hfs_curve_jacobian_fix_20260831/`.
