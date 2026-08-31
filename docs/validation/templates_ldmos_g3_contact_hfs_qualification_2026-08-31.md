# Templates/LDMOS G3 contact-HFS qualification (2026-08-31)

## Scope and frozen controls

This qualification tests one change only: the G3 electron/hole high-field
mobility drive remains `GradQuasiFermi` in the device interior and falls back
to the Tri3 cell electrostatic-gradient magnitude for transport cells touching
a contact node. IALMob, predictor, impact ionization, and threshold changes are
out of scope and remain disabled. The original 31 Sentaurus CurrentPlot biases
are the only curve-scoring points.

## Three-point fixed-state decision gate

The archived Sentaurus states at Vg = 1/6, 1/2, and 5/6 V were replayed through
the Vela SG operator before any self-consistent curve was attempted. The gate
passed and therefore authorized the more expensive curve run.

| Vg (V) | Baseline Vela/Sentaurus | Contact-HFS Vela/Sentaurus | Error improvement (dex) |
|---:|---:|---:|---:|
| 0.166667 | 24.4976 | 1.79505 | 1.13505 |
| 0.500000 | 11.0772 | 1.04621 | 1.02481 |
| 0.833333 | 11.0647 | 1.04570 | 1.02453 |

The median magnitude error fell from 1.04443 dex to 0.01962 dex. This is a
decisive fixed-state coefficient/operator improvement; it does not by itself
qualify the self-consistent state trajectory.

## Self-consistent 31-point qualification result

The Release run used the original 31 exact biases, no predictor, and the
frozen block absolute convergence ceilings:

- psi: 5e-8;
- electron continuity: 1e-11;
- hole continuity: 3e-10.

The baseline 100-iteration budget failed at Vg = 1/6 V. Increasing only the
iteration budget allowed 1/6 V to pass at iteration 126. A 400-iteration run
then accepted 0, 1/6, 1/3, and 1/2 V, but stopped at 2/3 V:

| Vg (V) | Status | Newton iterations | Vela/Sentaurus current | log error (dex) |
|---:|:---:|---:|---:|---:|
| 0.000000 | pass | 19 | 2.1184 | 0.32601 |
| 0.166667 | pass | 126 | 2.7198 | 0.43453 |
| 0.333333 | pass | 243 | 2.7304 | 0.43623 |
| 0.500000 | pass | 368 | 2.6871 | 0.42928 |
| 0.666667 | fail | 400 | not scored | not scored |

At the failed point, psi (1.00e-9) and hole continuity (4.70e-13) were already
inside their ceilings, while electron continuity remained 1.57e-10, 15.7x
above its ceiling. The maximum row remained node 702 in the contact-neighbour
region. The final linear/state increment was approximately 1.5e-13 while the
electron residual continued to decay only about 2.8 percent per iteration.
Consequently, raising `max_iter` again was rejected: it would mask a contact
branch/Jacobian slow mode rather than close the qualification honestly.

An unscored 1/12 V warm-start experiment and an explicit adjacent-cell
third-vertex mobility-Jacobian experiment did not remove the slow mode. The
third-vertex experiment was therefore not retained as a claimed fix.

## Decision

- The contact ElectricField fallback is **qualified for fixed-state G3 SG
  replay**.
- The self-consistent 31-point G3 curve is **not qualified** and has no P95,
  Vth, gm, or final KCL acceptance result.
- The known-difference ledger remains draft. The four accepted self-consistent
  points already show a roughly 0.43 dex current factor, so the fixed-state
  improvement must not be presented as curve closure.
- The next WP1.5 task is to audit node 702's contact-neighbour electron row,
  including contact-basin branch protection, mobility-field Jacobian support,
  and block-filter envelope decisions. IALMob and predictor remain disabled.

## Reproducibility artifacts

- Fixed-state runner: `scripts/run_templates_ldmos_g3_contact_hfs_replay.py`.
- Curve runner: `scripts/run_templates_ldmos_g3_contact_hfs_curve.py`.
- Fixed-state evidence (ignored staging):
  `reference_staging/templates_ldmos_g3_contact_hfs_replay_release_20260831/summary.json`.
- Final partial curve and diagnostics (ignored staging):
  `reference_staging/templates_ldmos_g3_contact_hfs_curve_release_iter400_20260831/`.
