# SLOT-LDMOS final-only closure globalization progress (2026-08-30)

## Scope

This checkpoint implements the first D0-D3 development batch for the
Delaunay, avalanche-off, `masetti_field` SLOT-LDMOS branch.  It does not claim
a converged 0.5 V replacement checkpoint or a Vela BVDS result.

## Solver semantic correction

`global_continuity_closure.mode = report/enforce` no longer adds the closure
ratio to an ordinary merit or device-block line-search norm.  Closure affects
the nonlinear path only through the explicit `block_closure_filter` mode.
`enforce` still rejects a final state that fails the configured two-regime
absolute/relative closure criterion.

The Newton regression verifies that `off`, `report`, and permissive `enforce`
produce the same merit trajectory through final convergence.

## Diagnostic additions

The opt-in local update CSV now records aggregate current, raw-predicted,
capped-predicted, and selected block norms; raw/capped/applied step norms; the
number and largest voltage magnitude of QF cap changes; and signed free
electron/hole residual sums.  The adaptive continuation script can enable the
trace for repeatable hotspot node ids and override QF limits and the existing
depleted-carrier diagonal floor without editing rendered decks.

## Strict 0.5 V reclosure matrix

All runs started from the provisional 0.5 V state whose original acceptance
reason was `initial_stall_residual_floor` with zero Newton iterations.

| Case | Result | Evidence |
|---|---|---|
| componentwise cap 0.1/0.05 V, no diagonal floor | accepted seven merit steps in 300 s, but did not finish | first raw step norm `4.46e11`; 5755 QF DOFs capped; line search now progresses instead of rejecting immediately |
| QF cap disabled | failed on iteration 1 with `nonfinite_residual` after 13 backtracks | the unconstrained raw direction cannot be used as a nonlinear trial |
| uniform trust region 0.1/0.05 V | accepted only extremely small changes | a few depleted QF DOFs force the global scaling ratio to approximately `1e-11`; device residual is nearly unchanged |
| componentwise cap plus `carrier_diagonal_floor` scale 1 | first three full steps accepted | first raw step norm reduced to `7.83e3`; capped DOFs reduced to 954 |
| componentwise cap plus `carrier_diagonal_floor` scale 10 | completed 16 accepted iterations, then failed only the final closure gate | final device blocks: psi `4.73e-9`, electron `4.01e-11`, hole `1.05e-7`; final hole free-row sum approximately `-3.57e-12` |

The scale-10 run confirms that the previous immediate blocker was acceptance
geometry plus depleted-carrier Jacobian conditioning.  Removing closure from
the merit path is necessary but not sufficient: the original provisional
state is not a strict steady solution, and its raw Newton direction contains
large, numerically inactive depleted-carrier components.

## Checkpoint persistence correction

The scale-10 run used `return_best_feasible_on_failure`.  Because the initial
state passed the absolute closure gate while later lower-device-residual
iterates temporarily did not, `best feasible` selected the initial state.
Consequently the emitted rejected `final` state was byte-identical to the
initial state and the 16 accepted device iterations could not be resumed.

The continuation driver now makes its CLI failure-state policy authoritative
even when the base deck enables `return_best_feasible_on_failure`. It also
copies the solver's `rejected_final_state_file` into the adaptive summary and
labels it as either `last_accepted_iterate` or `best_feasible`. This permits
bounded restart chunks without treating a rejected bias point as converged.

An eight-iteration repeat with the best-feasible rollback disabled reduced the
normalized device residual from `3.51465e-6` to `6.54018e-9` and persisted a
different state hash. A subsequent four-iteration chunk reduced it to
`4.83723e-9`. The hole device block decreased from `1.48377e-7` to
`1.13092e-7`, while the absolute hole closure mismatch changed from
`3.39261e-12` to `3.63401e-12`. This separates device convergence from the
final conservation mode and justifies a final-only corrector.

## Final-only restoration and strict checkpoints

`final_linear_qf_projection` recomputes the Jacobian only after the device
residual has reached `stall_residual_floor`, constructs a minimum-norm
free-QF correction for the two global rows, reapplies QF update bounds, and
backtracks against the full nonlinear residual. It accepts a state only when
the device floor, carrier rows, and enforced global closure all pass.

| Bias transition | Newton result | Final restoration | Id (A/um) | Closure mismatch e/h (A/um-equivalent residual) |
|---|---:|---:|---:|---:|
| provisional 0.5 V -> strict 0.5 V | 13 accepted device steps in restartable chunks plus one final validation step | correction norm `2.98276e-14`, alpha 1 | `3.05359e-16` | `5.29563e-13 / 3.58819e-13` |
| strict 0.5 V -> 0.625 V | 8 accepted steps; initial electron block `584.627` reduced to `9.98001e-10` | correction norm `4.23558e-14`, alpha 1 | `3.72073e-16` | `3.20996e-13 / 5.11552e-13` |

The strict 0.5 V current differs from the provisional value by approximately
`0.0196` decade, below the contract limit of `0.05` decade. The 0.625 V run no
longer fails on its first nonlinear step, directly validating the corrected
globalization semantics on the previously blocked transition.

Generated simulation outputs remain under ignored `build/` directories. The
authoritative strict states for this run are:

1. `build/slot_ldmos_globalization_debug_final_restore01/.../final_state.csv`
2. `build/slot_ldmos_globalization_debug_0p625/.../final_state.csv`

They are local validation artifacts and are not repository deliverables.

For final-only closure experiments, production runs should:

1. disable `return_best_feasible_on_failure` so the last accepted device
   iterate is available for restart; or
2. add a separate last-accepted-iterate checkpoint that does not reinterpret a
   failed point as converged.

## Next gate

Continue the avalanche-off ladder from the strict 0.625 V state to 1 V and
then larger bounded voltage points. Keep the same diagonal-floor,
componentwise-QF, merit, and final-only-restoration settings so each point is a
single-variable continuation experiment. Do not enable avalanche or IALMob
until the 1 V and 60 V pre-breakdown leakage gates have been evaluated.
