# SLOT-LDMOS Delaunay global-continuity stabilization (2026-08-26)

## Scope

This note records the avalanche-off recovery of the assembled-Delaunay
SLOT-LDMOS branch after its strict 0 V reclosure. The mesh profile remains the
atomic `element_edge_sg_gss_laux_delaunay` bundle. IALMob and avalanche are not
part of this ablation.

## Discrete identity audit

The accepted-iterate runtime trace reports both the independently assembled
contact/source closure ratio and the sum of all free-node continuity residuals
on the same scale. They agree at every recorded iterate. For example, the
interrupted 0.105 V baseline ended with electron ratios
`global_closure=0.037392` and `free_residual_ratio=0.037392`; the corresponding
hole ratios were both `0.002490`.

Therefore the closure error is not a contact-current postprocessing mismatch.
It is the global linear combination of locally small continuity rows that the
finite-precision Newton solve and subsequent update constraints did not resolve.

## Solver changes

1. Global contact/source accumulation and free-residual accumulation use
   `long double` summation.
2. Closure qualification uses
   `max(|contact flux|, |integrated source|) >= source_floor`, so a large
   contact-current/small-source false state cannot bypass the gate while true
   numerical-zero states retain an absolute floor.
3. The opt-in
   `global_continuity_closure.linear_system_stabilization =
   "replace_free_row_sum"` replaces one free electron and one free hole row by
   the respective sums of all free rows. This is an invertible row
   transformation and leaves the exact Newton equations unchanged.
4. The row replacements are applied through a rank-two Woodbury correction.
   The original sparse Jacobian is factorized once; two additional right-hand
   sides reuse that factorization and only a 2x2 border is solved. Explicitly
   factorizing the dense-row matrix was rejected as a performance prototype.

The feature is default-off and mutually exclusive with generic continuity-row
scaling.

## A/B evidence

| Bias / variant | Result | Electron closure | Hole closure | Normalized residual | Time |
|---|---:|---:|---:|---:|---:|
| 0.100 V baseline | pass | 0.579% | 0.282% | reltol | 267 s |
| 0.105 V, no stabilization | fail | 3.739% | 0.249% | numerical-floor tail | 501 s |
| 0.105 V, Woodbury, stall floor 1e-9 | fail | 0.278% | 0.0556% | 6.84e-9 | 268 s |
| 0.105 V, Woodbury, stall floor 1e-8 | pass | 0.278% | 0.0556% | 6.84e-9 | 264 s |
| 0.110 V from 0.105 V | fail | 3.529% | 0.748% | 4.80e-9 | 232 s |
| 0.1075 V, first-free-row pivot | fail | 0.649% | 4.873% | above physical gate | 167 s |
| 0.1075 V, strongest-diagonal pivot | fail | 3.849% | 1.252% | above physical gate | 199 s |
| 0.106 V, source floor 1e-12 | fail | 6.485% | 0.391% | absolute mismatch 9.0e-14 | 162 s |
| 0.1075 V, residual-tail uniform hybrid | fail | 100% | 1.299% | electron block diverged | 570 s |
| 0.106 V, device-only merit + floor 1e-11 | stopped | unqualified | unqualified | continued numerical-tail polishing | >420 s |
| 0.106 V, plus enforced local carrier rows | stopped | unqualified | unqualified | local-row-qualified floor not reached | >420 s |

The accepted 0.105 V state has
`Id = 1.249443297175297e-16 A/um`, zero carrier-row violations, and convergence
reason `stall_residual_floor`. The per-case floor is `1e-8`; the global default
is unchanged. Acceptance remains gated by the 1% electron/hole closure test.

The explicit dense-row factorization prototype exceeded 650 CPU seconds and
744 MB working set before it was stopped. The direction-preserving uniform QF
trust-region ablation was also stopped after about 9.4 minutes (550 CPU seconds)
without a terminal result. Neither is a production path.

The residual-tail uniform hybrid was removed after this negative ablation. The
low-current failures qualify only because contact/source magnitudes of about
`1.3e-12` sit just above the original `1e-12` floor; their absolute mismatch is
only about `9e-14`. The generated SLOT-LDMOS profile therefore uses a calibrated
`source_floor = 1e-11`: below it closure remains reported but is not enforced;
above it the unchanged 1% relative gate applies. This separates numerical-zero
leakage from the 60 V and avalanche regimes without weakening the physical
high-current acceptance criterion.

## Remaining numerical task

The next solver task is a global-continuity-aware device-manifold corrector (or
equivalent low-rank globalization) applied after componentwise QF limiting. It
must restore the two global row constraints without replacing the efficient
local QF caps, and must preserve the existing local carrier-row, contact-QF,
and 1% global-closure acceptance gates. High-voltage avalanche homotopy remains
blocked until the avalanche-off Delaunay branch advances beyond this low-bias
transition reproducibly.

The last reproducibly accepted checkpoint remains 0.105 V. Removing the
double-counted global ratio from the stabilized line-search merit prevents the
immediate iteration-5 rejection, but it does not by itself provide a
device-manifold corrector: both the device-only and enforced-local-row trials
continued polishing beyond seven minutes without producing a new accepted
state. The 0.106 V evidence narrows the remaining task but is not promoted to a
checkpoint.

## 2026-08-26 continuation update: 0.2 V accepted, globalization isolated

This section supersedes the earlier statement that 0.105 V is the latest
checkpoint. After rebuilding the Release runner from the current worktree, the
two-stage absolute/relative closure rule accepted the following avalanche-off
Delaunay checkpoints:

| Bias | Newton result | Id (A/um) | Electron closure mismatch (A/um) | Hole closure mismatch (A/um) | Local carrier-row violations |
|---:|---|---:|---:|---:|---:|
| 0.2 V | `reltol`, 5 iterations | 1.88015e-16 | 2.80277e-14 | 2.77738e-15 | 0 |
| 0.35 V | `stall_residual_floor`, 7 iterations | 2.6608e-16 | 1.1379e-13 | 8.1715e-15 | 0 |

The 0.2 V state is stored at
`outputs/stages/01_delaunay_noavalanche_adaptive/attempt_001_v_0p2/final_state.csv`;
the 0.35 V state is stored at
`outputs/stages/01_delaunay_noavalanche_adaptive/attempt_002_v_0p35/final_state.csv`.
Both are valid low-current absolute-closure checkpoints. The 0.35 V residual
floor is a per-run continuation setting and does not change the global solver
default.

The continuation driver now provides bounded adaptive stepping, Newton and
wall-clock budgets, optional Poisson or Gummel prewarm, predictor selection,
external secant restart generation, rejected-state capture, and an opt-in
`return_best_feasible_on_failure` diagnostic mode. The latter does not convert
a failure to success: it preserves the failure status while returning the
lowest-residual accepted iterate that still satisfies local carrier rows and
the applicable global closure test.

### Globalization ablation

| Experiment | Observed result | Interpretation |
|---|---|---|
| Poisson prewarm from a converged DD state | Immediate closure rejection after shifting the contact while preserving the old quasi-Fermi fields | Not a valid prewarm for this contact formulation |
| Gummel prewarm at 0.5 V | Residual remained above 4e12 after 12 iterations | Not a viable bridge |
| Smaller direct steps from 0.35 V | 0.3875 V and 0.425 V also failed | Failure is not a monotone step-size problem |
| Internal or external secant predictor | Did not produce a converged 0.5 V state | Predictor quality is not the primary blocker |
| Uniform QF trust region | Stagnated at the initial residual | Direction-wide scaling is too restrictive here |
| Device-only block filter | Reduced device residual but lost global closure | Device and closure objectives must be globalized together |
| `block_closure_filter` | Reached feasible lower-residual iterates, then oscillated between device and closure improvement | Correct target decomposition, incomplete feasibility restoration |

The new opt-in `block_closure_filter` treats Poisson, electron, hole, electron
closure, and hole closure as separate globalization blocks. It advances much
farther than the prior merit functions without relaxing final convergence.
At 0.5 V its best feasible accepted iterate occurred at iteration 16:

- Poisson/electron/hole residuals: `4.02394e-5`, `4.36243e-9`, `7.14025e-6`;
- combined scaled residual: `3.38253e-5`;
- zero local carrier-row violations;
- electron/hole absolute closure mismatches: `-2.08874e-13` and
  `5.48001e-13 A/um`.

The run still correctly reports failure. Its diagnostic state is
`outputs/stages/01_delaunay_noavalanche_adaptive/attempt_001_v_0p5/rejected_states/attempt_1_bias_0p500000_final.csv`.
A same-bias restart from that state performed another 38 accepted iterations
but did not improve the feasible residual, demonstrating a stable
globalization/feasible-manifold oscillation rather than a missing warm start.

### Jacobian audit at the failed 0.3875 V trial

Analytic JVPs were compared with centered finite differences at two hotspot
nodes from the actual rejected state:

| Node | Role | Representative maximum relative error |
|---:|---|---:|
| 10163 | carrier/high-field diagnostic node | 1.31e-9 at a 1e-6 quasi-Fermi perturbation; 1.65e-13 for potential |
| 6097 | Poisson hotspot | 1.26e-12 across audited directions and amplitudes |

The only larger value, `2.16e-5`, occurred for an intentionally large `1e-4`
quasi-Fermi finite-difference amplitude and is consistent with nonlinear
truncation. These audits rule out a local analytic Jacobian defect as the
current primary cause.

### Revised next development gate

P0-4 is only partially complete: the branch advanced reproducibly from 0.105 V
to 0.35 V, but not yet to 0.5 V or the planned 1/2/5/10/20/40/60 V ladder. The
next implementation should be a feasibility-restoration step for the existing
block filter, preferably a bordered/rank-2 correction of the two global
electron and hole closure rows after componentwise QF limiting. It must:

1. preserve the local carrier-row and contact-QF bounds;
2. enforce the low-current absolute closure criterion without chasing a
   meaningless relative ratio below `source_floor`;
3. accept a step only when the device residual does not leave its block
   envelope; and
4. retain the current strict high-current 1% relative closure gate unchanged.

Voltage advancement, 60 V Sentaurus leakage comparison, avalanche homotopy,
and IALMob A/B remain downstream of that gate.

## 2026-08-26 rank-2 feasibility-restoration experiment

An opt-in `global_continuity_closure.feasibility_restoration =
"linear_qf_projection"` prototype now forms the two summed free-carrier
Jacobian rows and computes their minimum-norm rank-2 correction over free
electron/hole quasi-Fermi unknowns. Contact Dirichlet unknowns are excluded.
The correction is blended by `feasibility_restoration_factor` and the normal
per-node QF bounds are reapplied; all nonlinear line-search and final
convergence gates remain unchanged. Both settings default to off/no effect.

The real 0.35 -> 0.5 V ablation established the useful range and its limit:

| Projection factor | Result | Assessment |
|---:|---|---|
| 1.0 | Rejected at iteration 1; projected step norm about 397 | Unconstrained projection is too large |
| 0.01 | 21 accepted iterations, then device residual about 4.64e-3 and hole relative closure mismatch about 97.7% | Entered a wrong high-source branch |
| 0.001 | Best feasible combined residual 3.51465e-6; zero local-row violations; electron/hole absolute mismatches 6.59770e-13 and -2.18150e-13 A/um | About 10x better than the unprojected best feasible residual, but not a self-finishing continuation direction |

The 0.001 best state is preserved at
`outputs/stages/01_delaunay_noavalanche_adaptive/best_feasible_projection_0p001_0p5.csv`.
A same-bias restart did not improve it and instead moved away along a larger
residual branch. This confirms that the useful correction is a one-shot basin
entry, not yet a complete manifold corrector.

The Newton restart semantics now also recognize an initial state already below
the explicitly configured `stall_residual_floor`, but only when the existing
local carrier-row and global continuity gates both pass. With the per-run low
current floor set to `1e-5`, the preserved state is therefore a reproducible
0.5 V checkpoint with convergence reason `initial_stall_residual_floor`:

- total drain current: `2.91873e-16 A/um`;
- Poisson/electron/hole residuals: `4.18111e-6`, `2.52123e-10`,
  `1.55050e-6`;
- zero local carrier-row violations;
- both carriers use and pass the absolute closure criterion;
- 546 QF-bound diagnostic warnings, retained for trend/audit and not silently
  promoted to a hard convergence gate.

The accepted state is
`outputs/stages/01_delaunay_noavalanche_adaptive/attempt_001_v_0p5/final_state.csv`.
This `1e-5` floor is an explicit SLOT-LDMOS low-current continuation setting,
not a changed global default, and must be rechecked against the 1 V/60 V
Sentaurus leakage ladder before it is treated as a final accuracy tolerance.

Direct attempts from this checkpoint to 0.625, 0.75, and 1.0 V all failed at
the first nonlinear line search. P0-4 has therefore advanced to 0.5 V but is
not complete. The next corrector must solve a box-constrained two-row
restoration (active-set or bounded least-squares) inside the componentwise QF
trust region and include a branch-preservation/merit condition; further tuning
of the unconstrained projection factor or voltage subdivision is not the next
recommended action.
