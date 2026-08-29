# Templates/LDMOS G3 WP1.5 linear-equilibration qualification

Date: 2026-08-29

Status: solver qualification passed; 31-point execution passed; L2 curve
acceptance failed on P95 log error and low-current KCL.

## Scope and invariants

This follow-up closes the registered node-4601 Poisson/Jacobian scaling task.
It retains the frozen G3 physical contract, exact 10,241-node topology,
Sentaurus CurrentPlot bias grid, `5e-8` Poisson ceiling, no predictor, no IALMob,
and no quantum/thermal/avalanche additions.

The two frozen input states are:

- drain-path best state at `Vd=0.00214466666666667 V`;
- exact Sentaurus `Vg=0 V, Vd=0.1 V` seed best state.

## Fixed-state node-4601 audit

The production Jacobian was replayed through the new read-only
`newton_poisson_linear_probe`. The probe compares the existing SparseLU solve
with a mathematically equivalent one-pass two-sided L2 equilibration after the
existing continuity-row scaling.

| Metric | Drain best state | Exact Id-Vg seed |
| --- | ---: | ---: |
| node-4601 one-ring raw condition | 2775.22 | 2770.66 |
| node-4601 one-ring equilibrated condition | 1.79390 | 1.72829 |
| full row-norm spread | 1.87362e22 | 3.15921e21 |
| full column-norm spread | 4.72804e31 | 5.90180e25 |
| raw linear closure norm | 2.86681e-8 | 2.77824e-9 |
| equilibrated linear closure norm | 8.53094e-11 | 1.52715e-11 |
| relative step difference | 1.34315 | 0.0703052 |

The closure reduction and finite-precision step differences show that the
unscaled factorization did not return an accurate enough Newton increment.
This is a linear algebra defect, not a Poisson formula, material, or line-search
threshold difference.

Generated evidence is below the short, ignored staging root
`reference_staging/templates_ldmos_g3_node4601_audit_20260829/`. Its
`frozen_inputs/` copies are immutable probe inputs; the drain input is
reproduced in an isolated `linear_equilibration=off` run so later production
runs cannot overwrite the registered failure state.

## Implementation

- `linear_equilibration.mode = l2_row_column` performs one row-L2 and one
  column-L2 scaling pass, solves the scaled system, and maps the increment back
  to the original Newton coordinates.
- The default is `off`; only generated high-field G3 decks enable it.
- Existing source-aware continuity-row scaling remains first in the chain.
- `newton_poisson_linear_probe` exports per-node Poisson/full row and column
  norms, raw/equilibrated psi increments, local condition estimates, and full
  linear closures without changing `solve()` state.

## Frozen entrance qualification

The final executable produced:

| Route | Result | Key result |
| --- | --- | --- |
| exact Sentaurus seed reclose | pass | 5 Newton iterations; psi 1.97296e-9; electron 7.94147e-11; hole 9.54337e-12 |
| Save/Load same-bias reclose | pass | `initial_block_abstol`, zero Newton iterations |
| exact drain prebias path | fail, retained diagnostic | node-4601 psi reduced from about 5.42e-7 to 5.90256e-10; failure moved to node-702 electron continuity 2.67202e-8 versus 2e-9 ceiling |

The production Id-Vg route is the independently qualified exact Sentaurus seed
followed by two same-bias recloses. The failed drain path is not substituted as
the production seed and remains an open WP1.5 carrier-continuity diagnostic.

## 31-point G3 Id-Vg result

All 31 exact shared bias points from 0 to 5 V converged. Every finite-bias point
used the original 0.1666667 V step, with no retry, predictor, or carrier-row
violation. The run took about 30.5 minutes from the first to final state write,
which consumes most of the current phase-A lower-bound budget and should be
used in the next budget re-freeze.

Against the sealed `G3-no-IALMob` Sentaurus curve:

| L2 metric | Result | Limit | Gate |
| --- | ---: | ---: | --- |
| resolved median absolute log error | 0.00657606 dex | 0.10 dex | pass |
| resolved P95 absolute log error | 0.287219 dex | 0.20 dex | **fail** |
| strong-inversion endpoint relative error | 1.45137% | 20% | pass |
| diagnostic 1e-8 A/um Vth absolute error | 29.6261 mV | 100 mV | pass |
| maximum gm relative error | 6.42497% | 20% | pass |
| worst resolved KCL relative error | 6.16280% | 1% | **fail** |

`Vg=0 V` is unresolved because the drain current is only 5.02 times the KCL
residual. The worst resolved KCL point is `Vg=0.1666667 V`. The five largest
curve errors lie from 0.1666667 to 0.8333333 V and are 0.2466--0.2920 dex,
consistent with the measured 29.6 mV subthreshold horizontal displacement.
The median, endpoint, and gm agreement must not mask these two hard failures.

The machine-readable staging result is `g3_idvg_qualification.json`; the
tracked rendering is `templates_ldmos_g3_idvg_qualification_2026-08-29.md`.

## Decision and next task

WP1.5 node-4601 Poisson scaling is closed and the 31-point execution path is
restored. Phase 3 is not accepted at L2. The next single-factor work is:

1. audit the 0.1667--0.8333 V electrostatic/reference-energy mapping that
   produces the 29.6 mV threshold displacement;
2. audit the `Vg=0.1667 V` source/drain/body current cancellation and terminal
   integration responsible for the 6.16% KCL ratio;
3. keep IALMob and predictor disabled until the classical G3 P95 and KCL gates
   pass or an evidence-backed known difference is approved.
