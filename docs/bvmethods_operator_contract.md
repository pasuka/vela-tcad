# BVMethods explicit coupled operator

The [precision budget](validation/bvmethods_precision_budget_2026-10-28.md)
separates statistical tangents, near-zero current directions and affine boundary
RHS terms. Native `-AvalDerivatives` retains the AV source; its RHS must not be
labeled a same-physics no-AV comparison against Vela `model=none`. The isolated
current-value controls diagnose matrix contributions without qualifying an
autonomous corrected operator.

The [primitive-state and Poisson baseline follow-up](validation/bvmethods_primitive_contract_2026-10-28.md)
adds an un-dispatched experimental normalized band/state interface and corrects
the shared explicit JD inverse coefficients from independent series reversion.
The matched-constant full Poisson residual closes in isolated controls; the
complete zero-bias AV and weak-hole RHS remain unqualified.

The [population SG arithmetic investigation](validation/bvmethods_population_sg_arithmetic_2026-10-28.md)
adds explicit `population_sg_current: band_population` for controlled comparison.
It retains the physical model and derivatives, bypasses QF only in SG flux
arithmetic, and does not establish complete native zero-bias equivalence.

`solver.coupled_operator: "bvmethods_tri3_qf"` selects an opt-in coupled
operator. The default remains `"legacy"`. This profile does not change the
PN2D BV template bundle. Its numerical qualification is recorded in
[the implementation report](validation/bvmethods_production_profile_2026-10-24.md).

A [full Release sweep](validation/bvmethods_full_profile_2026-10-25.md) reached
the configured stop current using the same core library and a local continuation
driver with a predicted inner-voltage step guard and accepted-state retry.
The default CLI continuation remains unchanged. An explicit
`sweep.external_circuit.coupled_continuation` now makes the same predictor guard
and accepted-state retry available to the production sweep. Availability of
this option alone does not qualify an autonomous full CMD sweep. The
[JSON / CMD three-step check](validation/bvmethods_cmd_e2e_2026-10-25.md)
establishes matching configuration, accepted states and Newton attempts for its
explicit BTBT-off control, but both autonomous paths stop at the same low-bias
minimum-step failure. Native zero-state
operator consistency and breakdown-region step convergence remain separate
qualifications.

The uninstalled [adaptive population/derivative cold-start experiment](validation/bvmethods_t470p_adaptive_bv_2026-10-26.md)
accepted 277 points on T470p before interruption under independently confirmed
system memory pressure from a hung log query. Its last inner voltage was
5.487691822 V; it did not reach the breakdown current. Low-bias conservation
reports and missing final iteration diagnostics remain qualification gaps.
The T470p rerun remains paused. A separately authorized
[vela VM comparison](validation/bvmethods_vm_predictor_comparison_2026-10-26.md)
completed its three-hour cold-start enforce/report-only predictor controls.
The [postprocessing analysis](validation/bvmethods_vm_predictor_analysis_2026-10-26.md)
records 240/49 accepted points and time-limit termination of both arms, without
complete BV qualification. Report-only advances farther but repeatedly returns
to bootstrap during AV-derivative restoration; residual, update and restoration
conditions remain distinct. Neither experiment changes production defaults.
The [same-checkpoint restoration control](validation/bvmethods_vm_restoration_verification_2026-10-26.md)
replays the original ten-update prefix and exports identical RHS and full
matrices for the competing choices. Holding full derivatives reaches the same
joint gates in two updates and remains settled through eight, while the current
selection returns to bootstrap. Six selected-column finite differences are
recorded separately from the unresolved, roundoff-sensitive near-null Newton
direction check. This is a local diagnostic result, without a production
restoration-policy or full-BV qualification.
The [automatic derivative-switch control](validation/bvmethods_vm_derivative_switch_2026-10-27.md)
keeps restoration selection separate from the final raw-update RMS. A diagnostic
bootstrap/trial/full state machine uses physical residual model quality and the
original physical gates, confirms restoration, and checks explicit injected
rejection and cooldown. It restores the same checkpoint in two updates, remains
settled through eight, and carries full mode through five additional points with
two updates each. The matched-target baseline and finite-state differences are
recorded separately. This uninstalled diagnostic does not change production
defaults or establish complete cold-start BV convergence.
A separately authorized [automatic cold-start full-sweep task](validation/bvmethods_vm_automatic_full_start_2026-10-27.md)
uses `bvmethods_adaptive_bv_probe` with explicit
`derivative_restoration_mode: "automatic"`. It inherits full mode only from an
accepted point and preserves independent corrector budgets on failed predictors
and accepted-state retries. The linked record is a startup snapshot, with the
first five nonzero points accepted; it is not a complete-BV qualification.
The [evaluation-reuse optimization](validation/bvmethods_evaluation_reuse_2026-10-27.md)
adds explicit diagnostic switches for candidate diagnostics, physical evaluations,
shared terminal assemblers, and boundary-only parameter differences. On T470p
and C560, two Release pairs in the cold first-five-point window preserve field
bits, complete attempt records and sampled operator files, with measured
speedups of 1.77 and 2.26 respectively. The bounded, single-owner physical cache
requires frozen mesh, doping and material inputs, compares full state/reference
bits, separates derivative modes and reapplies boundaries on every call.
All switches default off; stopping criteria, production/CMD defaults and the
currently frozen VM sweep remain unchanged. This does not qualify complete-BV
speedup or alter the outstanding nonlinear-convergence evidence.
The separately frozen [four-switch VM full sweep](validation/bvmethods_vm_four_switch_full_2026-10-27.md)
subsequently completed its cold start with 103 accepted points and exit code 0,
reaching the configured stop current in about 113 minutes. It used global
fraction-to-boundary updates and automatic derivative restoration, without the
later adaptive population-candidate selector or iteration-local reuse. Low-bias
conservation report violations and native-operator consistency remain separate
qualifications; the completed baseline does not qualify the later selector's
complete-BV performance.

## State and statistics

The [adaptive-candidate VM failure analysis](validation/bvmethods_vm_adaptive_failure_2026-10-27.md)
records 24 accepted points before the frozen pre-reduction scan stopped at
outer 1.079671328 V. Positive minority populations did not guarantee a finite
population-coordinate map: electron node 498 reached about 6.98e-301 cm^-3,
and the fixed density-reference column scale overflowed before the next
Jacobian solve. Bounded reuse on/off controls retain identical operators and
states. This identifies the direct exception and its state history, without
qualifying a scaling repair or the physical origin of repeated negative
minority directions. The core reduction qualification remains separate.

The [core prescribed-row reduction](validation/bvmethods_core_constraint_elimination_2026-10-27.md)
is now enabled by default for ordinary `LinearSolver` calls and the normal direct
bordered path. Device QF row IDs survive population-column transformations and
the circuit border; the backend receives the 6765-dimensional free system for
this mesh while state and physical residual arrays retain 7801 entries. Affine
RHS correction and recovery are part of the solve. `VELA_CONSTRAINT_ELIMINATION=0`
provides the unreduced control. This algebraic change does not qualify native
Jacobian support or complete BV convergence, and the frozen VM task uses its
previous source snapshot. Explicit QR/emergency fallback policies remain separate.

The unknowns are ordinary physical electrostatic and electron/hole
quasi-Fermi potentials, packed as `(psi, fn, fp)/V0`. Only zero-reference
coordinates are qualified; reference-tail and split-DD coordinates are rejected.
Populations are reconstructed at every evaluation from these potentials, the
material DOS, effective intrinsic density and BGN offsets. No imported population,
native SG cache, source array or native Jacobian is a constitutive input.

Select `solver.carrier_statistics.model: "fermi_dirac_joyce_dixon"` explicitly.
The Joyce-Dixon forward, inverse and derivatives form one contract, including
contact neutrality, density reconstruction, FD-BGN and generalized SRH. The
existing `fermi_dirac` defining-integral model remains unchanged. Selecting
Joyce-Dixon in the Newton JSON parser outside this profile is rejected.
Material files and the global CODATA constants remain explicit inputs; this
profile does not substitute Sentaurus thermal-voltage or charge constants.

Finite QF updates are not equivalent to finite population Newton updates.
The [population/acceptance audit](validation/bvmethods_population_acceptance_2026-10-26.md)
observes local native population positivity protection even with factor one,
and reconstructs the native RMS error from raw tangents and actual post-state
populations. That diagnostic does not install a population update rule in this
profile or qualify the current absolute residual as a native RHS stopping scale.
Equal-population replay after the first native step still has unresolved
weak-hole operator differences.
The [joint weak-hole audit](validation/bvmethods_weak_hole_joint_2026-10-26.md)
separates full-diagonal RHS amplification from genuine matrix differences.
Native continuous-iteration and fresh-Coupled linear systems are not assumed
interchangeable: phase controls preserve physical primary fields and exported
RHS values while changing continuity entries, concentrated in previously
protected population columns. The internal native mechanism is unqualified;
neither that matrix processing nor population projection is a production rule here.
The [continuous RHS directional audit](validation/bvmethods_continuous_rhs_jv_2026-10-26.md)
verifies selected population columns inside the same native Coupled stage.
Actual RHS responses agree with the fresh full-derivative matrix, while the
observed continuous effective system numerically matches a frozen-mobility,
no-AV-derivative control. Unprotected columns differ too; association with prior
population protection alone does not identify the internal mechanism. Keep
physical residual derivatives and a solver's approximate linear system distinct.
This evidence adds no approximate-linearization default to this profile.
The [linearization/population control](validation/bvmethods_linearization_population_control_2026-10-26.md)
keeps the physical RHS live while independently dropping transport mobility
tangents and AV tangents in an uninstalled probe. Combined with finite population
updates, it reaches the original residual/load gates and reported RMS below one
at two checkpoints; full-derivative population controls do not settle within
40 updates. The successful combined trajectories do not trigger the 1% positivity
fallback. Carrier-row and global-closure reports still contain violations.
This local result installs neither an approximate-Jacobian selector nor a
population update policy in the production solver, and does not qualify full BV.
The [population-Newton trials](validation/bvmethods_population_newton_trials_2026-10-26.md)
solve directly in scaled population columns and independently check physical
population finite differences. Global fraction-to-boundary startup controls
fail linear-accuracy guards on the recorded paths. Staged positive-population
warm starts permit full AV and mobility derivatives, with eight full-step
checks and five additional small-bias diagnostic points. A target-aware
restoration control uses the existing residual/load tolerances and reported
RMS, while the original two-metric baseline comparison is retained separately.
High-precision sums of the stored ledgers reproduce the outstanding electron
closure ratios: changing global summation alone does not remove the remaining
free-equation imbalance. Local state/operator precision and stopping scales
remain unqualified. These trials add no production coordinate, feasibility,
derivative-restoration, or full-CMD default.
The [population-protection comparison](validation/bvmethods_population_protection_comparison_2026-10-27.md)
adds explicit `population_protection: "local_nonpositive"` to the uninstalled
adaptive probe; its default remains `"fraction_to_boundary"`. The local policy
replaces only nonpositive finite population candidates with 1% of their old
population, retaining positive candidates and the existing circuit update cap.
An explicit `protection_checkpoint` mode replays common imported initial states
without qualifying autonomous cold start or full BV. On T470p and C560, both
policies have identical first-step operators: local protection passes the
predictor checkpoint in seven updates, whereas the fraction policy fails;
starting from the prior accepted state reverses the pass/fail outcome. These
results retain the original residual, raw-update and conservation controls,
and establish no universally preferred policy or production default change.
The [candidate-selection and iteration-reuse control](validation/bvmethods_candidate_selection_2026-10-27.md)
adds opt-in `population_protection: "adaptive"` to the diagnostic probe. It
prefers candidates satisfying the original joint gates, otherwise comparing
device and load residuals normalized by their existing stopping tolerances.
The earlier fixed-row L2 candidate control remains available explicitly and its
failed predictor replay is retained. On both hosts the tolerance-based selector
passes the predictor and prior-accepted-state checks, followed by five matched
points. Optional `direction_reuse` and `iteration_linearization_reuse` avoid
duplicate directions, initial-state linearizations and same-state residuals.
Their bounded scopes preserve sampled matrices, iteration records and physical
state bits, including reversed ordinary-Release timing pairs without matrix
output. The fraction policy remains the default and both new reuse switches
default off;
this is not a complete-BV or native-operator qualification.

In the explicitly selected Joyce-Dixon dilute limit, DOS-weighted populations
and derivatives are evaluated in logarithmic scale before the normalized
population becomes subnormal. The matching density inverse and BVMethods SRH
degeneracy law use the same limit. Actual positive, finite populations remain
required; no density floor is introduced. See the [extreme-density recovery
audit](validation/bvmethods_full_recovery_2026-10-25.md).

## Discrete operator

- Poisson and contact/gauge rows retain the existing cell-material assembler.
- Carrier transport uses Tri3 coefficient-derived box weights, total-impurity
  Masetti mobility, distance-gradient Enhanced Lombardi and Canali field
  mobility with live derivatives.
- SG is evaluated with log-density and `expm1` quasi-Fermi differences and the
  matching generalized Einstein factor. Flat quasi-Fermi fields give zero
  equilibrium current algebraically; no small-current cutoff is used.
  Extreme population ratios use a difference of logarithms; large QF
  differences use the equivalent two-Bernoulli flux, avoiding an overflowing
  exponential intermediate while retaining both endpoint derivatives.
- AV uses the representative vertex and its two incident SG edges to recover
  the current vector, the existing contact driving-field rule, an unsmoothed
  Euclidean norm and vertex-box source integration. The same pair-generation
  residual enters both continuity equations. Local jets differentiate that
  residual; no captured current direction is supplied.
- SRH uses the same statistics and coefficient-derived silicon control volume.
- The residual, full Jacobian and carrier-term ledger share this runtime.
  Terminal totals and their state derivatives use its source-free transport
  residual directly. Drift/diffusion partitions and legacy edge diagnostics
  are not qualified for this operator. Legacy drift/diffusion CSV columns
  retain default placeholders; use carrier totals and terminal total only.

At exactly zero current the norm has directional derivatives, not a unique
Frechet derivative. The existing zero-norm jet convention is explicit. A finite
difference crossing that point cannot establish the local derivative of an
arbitrarily small nonzero native current. Native zero-state bitwise agreement
is a separate qualification from autonomous nonlinear convergence.

## Qualified configuration boundary

Version 1 requires classical 300 K silicon Tri3 DD, cell-material Poisson
geometry, `masetti_field_lombardi`, `total_impurity`,
`quasi_fermi_gradient` mobility with field derivatives, and SRH or no
recombination. AV is either disabled (for source-free transport evaluation) or
self-consistent Van Overstraeten with GradQF/current-density generation and
`local_ad`. AV cutoffs, reference-density interpolation and source overrides
are rejected. Quantum, density diagonal floors, thermionic contacts, Auger,
BTBT and split DD are not qualified. The DC sweep accepts explicit basic-state
VTK fields (`ElectrostaticPotential`, `eDensity`, `hDensity`, `SpaceCharge`,
`DopingConcentration`). It rejects implicit/all-field VTK output and other
current/source fields, plus legacy contact-edge/QF-floor diagnostics, whose
field/current operators differ.

The profile is available through JSON and C++ configuration. SDevice deck
translation does not automatically enable it. An explicit sidecar selector
`bvmethods_operator_profile: "bvmethods_tri3_qf"` enables the atomic coordinate,
statistics and geometry bundle for supported 300 K Fermi/OldSlotboom,
Masetti/Enormal/GradQF, SRH and GradQF Van Overstraeten CMD inputs with
`AvalDerivatives` and BTBT disabled. Physical models and coefficients remain
owned by CMD/PAR. Incompatible Physics or unsupported Plot fields are errors.
Conservation reporting/enforcement must be selected explicitly in the sidecar;
the selector does not silently relax those gates. Outer adaptive bounds and
growth/shrink settings are derived from CMD, without turning MaxStep into a
mandatory output grid when no Plot/CurrentPlot times were requested.
For this explicit profile the frontend defaults both
`quasi_fermi_update_limit_V` and `quasi_fermi_update_limit_minority_V` to zero,
including when a preceding warm-start stage supplied an implicit cap. An
explicit sidecar value is retained for controlled comparisons. Raw updates
remain diagnostic outputs; the separate circuit voltage update limit and
the selected acceptance gates still apply. The CMD profile uses `full_step`
for undamped Sentaurus settings, resolving Math and local Coupled controls.
Enabled vendor line search or Bank–Rose damping is explicitly unsupported,
instead of substituting a merit rule. Explicit sidecar modes remain diagnostic
overrides. This does not
change the global C++ defaults or PN2D template policy. The bounded
[same-state QF audit](validation/bvmethods_qf_direction_2026-10-25.md) records
the effect and its multi-step limitations.
For a controlled full-step investigation, explicitly set
`sweep.external_circuit.coupled_line_search_mode: "report_only"`. The merit
backtracking recommendation is recorded but does not scale the applied step;
finite full steps may increase the merit. Device and inner-voltage update
limits are separate controls. Invalid full steps fail without a damped
substitution, and all final acceptance gates remain active. A diagnostic CSV
can be selected through `coupled_update_convergence.diagnostic_csv` with
`mode: "off"` to retain residual-only stopping.
To return to the former physics
path select both `coupled_operator: "legacy"` and the previous statistics model.
Newton and load-line stopping rules and budgets are unchanged.

## Verification entry points

`test_bvmethods_av_contract` contains `[bvmethods_profile]` tests for full Jv,
source pairing, contact rows, ledger identity, equilibrium and domain guards,
and terminal current. `test_joyce_dixon_contract` covers the public statistics
and FD-BGN dispatch.

The uninstalled `scripts/diagnostics/bvmethods_profile_probe.cpp` links against
the current Release `vela_core`. It exports zero and finite-current operators
and drives the production load-line Newton API through multiple outer voltages.
Build it explicitly with
`cmake --build --preset windows-ucrt64-release --target vela_bvmethods_profile_probe`.
It is excluded from the normal build and installation.

The full-sweep controller is retained separately in
`scripts/diagnostics/bvmethods_full_sweep_probe.cpp`. Build its uninstalled target
with `cmake --build --preset windows-ucrt64-release --target vela_bvmethods_full_sweep_probe`,
then run `build-release/vela_bvmethods_full_sweep_probe manifest.json output-directory`.
The trusted local manifest contains `config` (a JSON configuration path),
`previous` and `current` physical checkpoints. Optional `initial_step_V`
overrides the restart's first outer step; otherwise the script supplies it.
Each checkpoint contains `outer_V`, `inner_V`, and mesh-order arrays
`psi_V`, `fn_V`, `fp_V`. Use Vela's own accepted states, retaining input hashes
and rechecking the restart anchor; the manifest does not supply populations,
currents, AV sources or a Jacobian.

The original 0.02 V predictor cap was an empirical stabilization choice for
the recorded BVMethods run, not a material law or a universal device scale.
It has not undergone a breakdown-region step-convergence qualification.
The diagnostic now resolves its cap in this order:

1. Manifest `inner_step_control.max_step_V`: explicit physical inner volts.
2. Manifest `inner_step_control.max_step_fraction` times the absolute span
   between `initial_voltage_V` and `goal_voltage_V` in that same object.
   These endpoints are device voltages. Only without an enabled series
   resistor may missing endpoints come from `sweep.start/stop`.
3. An explicitly present `sweep.external_circuit.max_inner_voltage_step_V`.
   This reuses the script's inner-voltage scale for prediction and also
   preserves its ordinary corrector-update limit. Its parser-default 0.1 V
   is not treated as a script-provided predictor setting.
4. `sweep.max_step` when the input is a direct-voltage sweep without an
   enabled external resistor. This value is already in volts after frontend
   conversion; do not multiply it by the voltage span a second time.
5. The fallback 0.02 V, only when no applicable setting or complete
   device-voltage conversion is available.

For example, either of these manifest additions selects 0.05 V:

```json
{"inner_step_control": {"max_step_V": 0.05}}
```

```json
{"inner_step_control": {
  "max_step_fraction": 0.005,
  "initial_voltage_V": 0.0,
  "goal_voltage_V": 10.0
}}
```

The first explicit form takes precedence within that object. An incomplete or
zero device span permits the remaining resolution rules and eventual fallback;
malformed, nonpositive or nonfinite explicit step sizes raise errors. The
summary records `inner_step_guard_V`, `inner_step_guard_source` and
`inner_step_guard_reason`, including a reason when fallback is needed.
For a resistor sweep, existing `InitialStep/MaxStep` controls belong to the
external voltage. They do not uniquely determine an independent inner cap.

Outer bounds, growth, shrink and retry settings now come from the script's
`coupled_*` overrides, then the corresponding `sweep` fields, with the former
diagnostic constants as missing-field defaults. The predictor limits the next
outer step using the accepted secant's absolute inner/outer voltage slope;
this permits large external steps and reversed inner motion. If the cap
requires a step below the configured outer minimum, it reports an incompatible
step limit rather than silently exceeding the requested cap.

The case-specific diagnostic still uses 1e7 ohm um and 40 Newton iterations.
Its strict equation/load tolerances are 1e-8 and 1e-6 V. It stops at
`abs(Id) >= 1.3e-3 A/um` or outer voltage 100000 V, with a 512-attempt budget.
It records every completed attempt and accepted physical state. These remain
diagnostic controls; the default CMD/CLI path still does not activate the
predictor cap or seed retry. The explicit production option has the same cap
resolution order, using its `coupled_continuation.inner_step_control` object
instead of the local manifest object. Malformed explicit settings are errors.
`scripts/run_bvmethods_profile_regression.py` validates the state-only manifest
and rejects an incomplete or failed continuation. Every rejected solve remains
in the evidence; step halving retries from the last accepted state with the same
iteration budget and stopping rules.
