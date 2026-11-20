# Config Schema Reference

## Bounded nonzero-supply circuit diagnostic

`simulation_type: newton_bounded_supply_load_line_iteration_probe` permits one
bordered device/contact-voltage update, with `abs(supply_voltage_V) <= 0.2` and
`0 < max_abs_device_voltage_V <= 0.2` (default 0.2 V). All other terminals must
be at zero and `resistance_ohm_um` must be positive. It exports the physical and
stabilized systems and reports `Vdevice + R*I - Vsupply`; it never qualifies a
state or starts a sweep. Its uniform direction cap and merit trial are diagnostic
controls, not a replacement for the standard Newton update/qualification policy.
With the default policy the contact voltage remains binary64: a below-ULP voltage correction can round
away while split node coordinates retain it, causing trial boundary residuals.
This interface therefore does not yet establish a production circuit continuation
contract for voltage low parts.
The opt-in `trial_boundary_policy: rounded_voltage_projection` instead projects
all prescribed node high/low coordinates to the boundary targets built from the
actually rounded scalar voltage after each candidate. Its default is
`linear_direction`, preserving the existing trial. Projection requires split
state and changes only this diagnostic finite update, not the frozen matrix,
RHS, raw direction, physical equations, or qualification rules. It is not a
contact-voltage low-part implementation. `trace_weak_rows: true` adds actual
flux/source/residual samples for initially failing rows and the unprojected
linear prediction at every tested alpha. Requested/applied voltage changes,
boundary residual and full high/low coordinate changes are reported separately;
the generic probe completion flag never grants device qualification.
`paired_trial_audit: true` (default false) retains one assembly and one captured
linear direction, tests both `linear_direction` and `rounded_voltage_projection`
at all 13 powers-of-two trial scales, and archives all 26 unaccepted candidates.
Each record includes all-device row qualification, bad-row identities and raw
reaction electron/hole/total terminal currents. It performs no accepted update.
The optional constant-material SG directional audit excludes source-bearing and
IAL/HFS-supported rows; it is not a replacement Jacobian for silicon or a solver
acceptance criterion.

`trial_boundary_policy: split_contact_voltage` is an independent experimental
policy for the bounded probe. It retains the voltage increment in a high/low
pair, constructs boundary targets from both parts, uses the same voltage column,
and evaluates the load line with high/low voltage arithmetic. It requires split
state and unrelaxed Dirichlet targets. In paired-audit mode this policy archives
its own 13 trials instead of the 26-policy comparison. The terminal current is
still a double-valued production result; split voltage does not imply arbitrary
precision terminal current. Archives record `contact_biases_low_V` separately.
This remains a diagnostic circuit interface, not an enabled long-sweep driver.
`constant_relative_coordinates: true` is an additional default-off diagnostic
restricted to a split paired audit with zero carrier regularization. It forms
one common coordinate and relative coordinates for each eligible constant-material
carrier block. Common columns use stable SG directional derivatives; other
columns, RHS, constraints and closure-row semantics are retained. Groups touching
the drain circuit-reaction derivative are rejected. Common and relative increments
are applied separately to the split state so their sum is not prematurely rounded.
`coordinate_backend.mtx` is the actual matrix supplied to the direct backend in
this mode; `backend.mtx` remains the legacy-system comparator. Coordinate-system
defects and separate increment components are exported. This diagnostic is a
different linear representation, not part of the same-matrix update-policy pair,
and still grants no qualification or sweep acceptance.
The original zero-supply probes retain their zero-supply and 10 mV restrictions.
Use a new output directory and an external hard deadline. Qualify any resulting
state separately with the original device gates before using it as a parent.

## Bounded SLOT population memory loop

For explicit three-block `[psi, phin, phip]` comparisons,
`VELA_LINEAR_SOLVER=pardiso_mkl_eliminated_keep_psi` retains the complete
potential block, including Dirichlet equations, and eliminates only exact
carrier diagonal-singleton constraints. The SparseLU counterpart is
`sparselu_eliminated_keep_psi`. These opt-in backends require equal contiguous
nodal blocks; divisibility by three is checked, but the caller must ensure
the semantic layout. They are not appropriate for Poisson-only, interleaved,
or extended multiphysics systems. Existing full-elimination backends and
global defaults are unchanged. Affine carrier constraints still correct the
reduced RHS, and all directions are lifted to the original layout.

`simulation_type: slot_population_memory` is an opt-in, single-target restart
driver, not a change to default Newton. It keeps the loaded mesh, materials and
reference/increment state in one process. `state_file` must contain referenced
Fermi QFs without a packed-low state; carrier regularization and closure-row
replacement must be off. Use the validated SLOT unit-scaled 300 K configuration.
The entry rejects an incompatible temperature, scaling, RHS tolerance, line
search or QF cap instead of silently replacing that policy. Initial boundary
residuals are always measured; they block only under `primary_with_guards`.

`population_memory` requires `output_directory` (a new directory) and
`support_file` (JSON arrays `[psi_node_ids, electron_node_ids, hole_node_ids]`).
The explicit support preserves the existing physical RMS update-error node set;
it must be bound to the mesh/configuration by the calling workflow. Optional
`max_updates` is 1–40 (default 40), `soft_seconds` is in (0,1080] (default 1080).
The soft bound is checked between updates; an external hard deadline is still
required for unattended runs.

The loop applies the established 0.001 concentration positivity projection,
relative Fermi inverse and 0.0259 V logarithmic potential damping.
`population_memory.dilute_log_inverse: true` (default false) selects a stable
log-Fermi-ratio finite inverse when both initial and target reduced populations
are at most 1e-12. It retains the production dilute fugacity correction and
falls back to the original relative inverse outside that domain. This is an
arithmetic selection boundary, not a density floor. The 1e-7 inverse guard,
positivity projection, physics and acceptance policy are unchanged. The option
is independent of the existing small-tangent finite-update arithmetic.
`population_memory.update_diagnostics: true` (default false) enables the
production-entry observer: each completed step writes `update_stages_K.csv`,
`weak_rows_K.csv`, and a `convergence_report` in history. Raw and actual RMS
use the original support, post-update denominator and 1e5 normalization;
actual RMS never participates in acceptance. QF requested finite displacement,
increment-coordinate writeback and physical high-field change are distinct.
References do not change inside this update; a zero physical-high change does
not imply zero increment change. The psi block in this CSV uses the same voltage
columns; density columns are zero/not applicable there. Population deltas and
denominators use internal concentration, with the conversion in JSON.
This option also enables the existing projection diagnostics. Weak rows use the
configured production eligibility gate, separately from the existing all-active
row count. A disabled gate emits an empty ledger with `enabled=false`, not a
claim of physical balance. Incident materials include both sides of interfaces.
External load-line and file-reread checks are marked unevaluated until their
controller performs them; no successful result is inferred from missing data.

For `newton_solve_from_state` and `newton_zero_update_qualification_probe`,
top-level `convergence_diagnostics_directory` optionally creates a new directory
containing `report.json` and `weak_rows.csv`. It separates the solver outcome,
configured row qualification and closure from externally owned checks. It does
not assign the population-loop U definition to coupled Newton. The ledger uses
the actual gate implementation, records absolute current residuals in A/µm and
retains internal flux/source scales and the effective global qualification
thresholds. This output does not change acceptance. Existing directories are
rejected to avoid overwriting evidence.
The standalone weak-row ledger is a final-state re-evaluation of the returned
solution; the solver's `finalCarrierRowConvergence` result remains the original
acceptance decision. Compare them explicitly when auditing restart arithmetic.
Two independent update diagnostics are default off:
`population_memory.stable_finite_update` evaluates the nonlinear remainder of
small, unprojected Fermi density steps separately from their linear QF component,
avoiding cancellation between the electrostatic step and inverse population
step. Large, projected or held steps retain the relative inverse. The existing
population reconstruction guard still applies. This is a local arithmetic
experiment, not evidence of full-device qualification.
`population_memory.consistent_psi_writeback` additionally uses the electrostatic
increment actually stored in binary64 when restoring QF; it requires the first
option. It does **not** retain sub-ULP electrostatic increments or enable the
BVDS split runtime for Fermi statistics. Its effect on weak transport must be
verified separately. Each history entry reports the requested/stored psi
difference, changed QF-step magnitude and number of stable finite updates.
`population_memory.acceptance_policy` selects the SLOT-only acceptance:

- `primary_only` (current default for this entry): RHS <= 1e-10 **or** raw
  physical update RMS < 1 is sufficient for nonlinear convergence. Closure,
  initial/final boundary residuals, boundary preservation and excluded QF
  directions are reports, not vetoes. Existing contact update masks are
  retained; reporting a violation does not change the boundary equations.
- `primary_with_guards`: reproduce the earlier additional closure/boundary
  vetoes. Use this explicitly when replaying older configurations.

Both policies retain numerical finiteness, positive-population,
inverse-mapping and mixed linear-error protections. Row ratios are reported.
Saved-state RHS is diagnostic, while saved-state residual consistency remains
mandatory. `qualified` is relative to the selected policy; inspect
`strict_qualification_satisfied` for the former closure/boundary qualification,
and `rhs_satisfied` / `raw_update_satisfied` for the stopping branch. These are
SLOT-specific rules, not a change to default production Newton. The Python
acceptance helper requires explicit `policy="primary_only"`; its old default
is preserved for callers outside this entry and frozen historical batches.

Normal updates write small progress records, not matrix/state snapshots.
Qualification or abnormal exit saves a complete matrix and state; the matrix
has its own `*_matrix_parent.h5` because it belongs to the state before the
last update. `validation_snapshots: true` additionally saves each update and
matrix for paired tests only. This entry does not implement bias adaptation,
secant prediction or a long scan.
`validation_snapshot_max_updates` (default 40) limits these additional snapshots.
`projection_diagnostics: true` writes per-node raw/used population increments,
projection intersections/streaks, density and QF extrema, original row terms,
and contact currents. The fixed [-0.5,5.5] V SLOT QF box is explicitly a shadow
diagnostic, not a physical bound or acceptance gate. PARDISO pivot perturbation
and refinement counts are retained when available.

`diagnostic_population_active_set: true` is a separate bounded experiment,
default off. Previously projected free carriers with nonpositive original
population direction are held using `dQF=dpsi`; a positive original direction
releases them. A separate modified matrix replaces those carrier rows, while
the original matrix, RHS and direction remain available. Finite recovery uses
the applied (log-damped) psi change to preserve held population. This assumes
the fixed DOS/bandgap of this SLOT configuration. Raw update error continues
to use the **unconstrained** direction, and original residuals are evaluated.
Both linear systems must pass the unchanged linear protections. No constrained
state receives `qualified=true` or becomes a production continuation parent,
even if `original_policy_satisfied` is reported. Stable held sets, small
constrained directions and the shadow QF box do not certify original equations.

The separate `population_memory.diagnostic_freeze_transport_qf_mobility: true`
arm uses the existing diagnostic transport derivative omission for its update
direction. At every step it also assembles/solves the unmodified system, checks
bitwise equality of the two residual/RHS vectors, retains the original direction
for raw U, and checks both linear solves. It cannot be combined with the active
set flag and never returns a qualified production parent. Density, source and
Poisson derivative chains and the nonlinear residual model stay enabled. This
is an opt-in inexact-Newton experiment, not a default physics or Jacobian change.
`population_memory.diagnostic_stop_on_primary: true` may be combined only with
this mobility diagnostic. It stops at the first original R/U success and writes
the diagnostic snapshot with `stop_reason=diagnostic_primary_satisfied`; it still
returns `qualified=false`. An explicit diagnostic continuation controller must
check numerical guards and saved-state consistency before using that snapshot.
The default remains false, retaining bounded tail-update experiments.
`diagnostic_density_hold` is a separate opt-in object for the mobility
diagnostic, with `enabled` (default false), positive `threshold_cm3` (default
1e-6), and `consecutive_nonpositive` (at least 2, default 2). A free carrier
below that threshold is held only after repeated unconstrained proposals
would make its population nonpositive. This is a trigger, not a density floor:
no population is raised to the threshold. The actual diagnostic system is
resolved with `dQF=dpsi` on held rows; a positive unconstrained population
direction releases the constraint. Finite recovery uses the applied damped
psi displacement. Original equations, their residuals and original raw U
remain measured; held-row residuals and activation/release counts are reported.
This option cannot be combined with `diagnostic_population_active_set` and
does not qualify a production parent. A stable held set is not convergence.
`diagnostic_continue_after_primary` (default false, diagnostic only) permits
bounded tail updates after R/U passes. `capture_steps` lists update numbers
at which to retain matrices, pre-update states, actual and unconstrained
directions, and detailed population/row diagnostics. Actual held systems have
different rows/RHSs; compare their original physical residuals separately.
`diagnostic_source_linearization` (integer 0..2, default 0) is available only
with the shared transport-mobility diagnostic pair and PhuMob element avalanche.
Mode 0 retains the complete avalanche Jacobian; mode 1 omits only the direct
QF-field chain of mobility inside the avalanche source (retaining contact
electric-field, density, alpha and SG-current derivatives); mode 2 lags all
avalanche derivatives for the actual direction. All modes retain the original
live source values, physical residual, original matrix/direction and raw U.
Both linear solves and RHS bitwise identity are checked before an update.
These are explicit inexact-Newton diagnostic variants, not physics switches or
production convergence certificates. No production Newton configuration reads
this option. Temporary density hold may be combined with either variant, but
held-row residuals must still be reported and release requires an unconstrained
positive population proposal from that arm's actual diagnostic system.
Optional `linear_row_refinement_sweeps` (integer 0–4, default 0) tries bounded
coordinate residual correction on mixed-error failing nonsingleton rows. It
uses long-double defect accumulation and commits the corrected double direction
only if all original mixed, block and constraint guards pass; otherwise the
original direction is retained and the original guard rejects it. No tolerance
or residual/Jacobian coefficient changes. Both diagnostic directions are checked.
The independent `linear_factor_refinement_steps` option (0–3, default 0)
instead solves correction RHSs from long-double defects with the original
system's retained direct factors. Vendor elimination retains its reduced solver
and recomputes affine correction RHSs. Failed corrections leave the original
direction intact for rejection. This option cannot be combined with row
refinement; it retains the same acceptance guards and checks both systems.
PhuMob caches the G-clamp minimum by exact temperature, carrier mass ratio and
all seven G coefficients in a bounded thread-local cache. This quantity does not
depend on live populations; fluxes, mobility derivatives and residuals are still
reevaluated. `VELA_PHUMOB_G_MINIMUM_CACHE=0` before process startup selects the
original search for diagnostic comparison. It does not select another model.
The SLOT memory loop also saves `linear_quality.profile` and
`evaluation_profile`, including nested assembly and linear-solver stages;
parent/child stage times must not be added together as disjoint costs.

`shared_diagnostics` (default false, mobility diagnostic only) shares assembler
preparation across the original/altered systems, rechecks RHS bitwise equality,
and evaluates all contact functionals from one raw residual. Original-system U
semantics remain unchanged. `compact_diagnostics` (default false) retains scalar
diagnostics each step and detailed CSVs at primary success or bounded end points;
matrix snapshots store one unchanged operator with `matrix_aliases.json` naming
its identical physical/regularized views. Exception snapshots retain the state
and matrix even when no completed-step CSV exists. Per-step `timing` separates
system evaluation, quality/refinement, carriers, closure, ports and output;
`wall_s` includes step output and `elapsed_s` includes final snapshots. The
in-progress JSON can precede its own final write-time measurement; use the
completed result for timing totals. These options are off by default.
Use `result.json` and the process exit status as the completion record;
`progress.json` is an intermediate snapshot. Configuration and update-error
support are copied into the output directory. Input mesh/material/state hashes
and the externally enforced hard deadline remain the caller's responsibility.

## Diagnostic linear system snapshots

`simulation_type: newton_load_line_iteration_probe` is a bounded zero-bias
diagnostic. All configured contacts and `supply_voltage_V` must be zero.
It adds the actual series-resistor equation using `contact` (default `drain`),
`resistance_ohm_um` (default `1e12`) and `bordered_linear_backend` (default
`pardiso_mkl`; also `sparselu` or `umfpack`, with no fallback). It retains every
potential equation and eliminates carrier constraints by affine substitution,
including their dependence on the shared contact voltage. It exports physical,
stabilized, reduced and equilibrated systems to a fresh
`linear_system_output_dir`, plus the elimination map, direction, trial state and
summary. Existing regularization and free-row replacement are retained; the
same row transform is applied to the voltage column. Split state low parts are
preserved. At most one descending coupled update is retained, with uniform
direction-preserving QF damping and up to 12 halvings. This is an independent
diagnostic update policy, not a production sweep or a device qualification.
Original row/closure results remain reports and `qualified` is always false.

`newton_zero_supply_load_line_iteration_probe` provides the same independent
single-update diagnostic while retaining a small nonzero internal contact
voltage from the previous diagnostic state. External supply and all other
terminal biases must remain zero; both initial and trial internal voltages
must lie within 10 mV of zero. It is not a nonzero-supply sweep mode. The caller
must bind the saved state to its actual internal contact voltage.

`newton_zero_update_qualification_probe` runs the original device qualification
without a Newton update or post-loop feasibility restoration. It retains the
ordinary initial boundary projection, so the optional `output_state_file` may
differ from the input at contact nodes; audit that change separately. Device
convergence does not qualify an external circuit: the saved state remains
marked diagnostic and unqualified, and its load-line equation needs a separate
check. Existing residual, row, closure and boundary tolerances are unchanged.

`simulation_type: newton_linear_system_probe` loads `state_file` without
accepting an update. `linear_system_output_dir` must name a new directory.
It exports Matrix Market matrices before and after carrier regularization and
`replace_free_row_sum`, plus residual, solve RHS, raw direction and replacement
row IDs. Matrix coordinates are internal scaled coordinates; `direction_V`
uses the reported potential scale. The first matrix includes boundary/gauge
equations; use `newton_carrier_term_probe` for pre-boundary physical rows.
This probe rejects row scaling, linear equilibration and iterative refinement
rather than silently auditing a different solve. Caps, Poisson recorrection,
feasibility projection and line search are not applied.

`linear_system_output_profile` is `full` by default. Opt-in `solve_only`
omits `physical.mtx` and `regularized.mtx`; it retains `stabilized.mtx`,
`vectors.csv`, replacement rows and all actual assembly/solve operations.
It changes diagnostic output only, not the direction or acceptance criteria.
Other profile names are rejected.

Only this probe accepts top-level
`diagnostic_freeze_transport_qf_mobility: true`. It removes the direct QF-gradient
mobility chain from the element-box PhuMob **transport Jacobian**, preserving
the residual model, density/Poisson chains and avalanche source derivatives.
The dedicated `newton_diagnostic_transport_qf_solve` mode applies the same
diagnostic linearization to a bounded `state_file` restart. It is not a
production solver option or a qualification waiver. The original
`mobility.jacobian_field_derivatives` compatibility guards remain in force.

These diagnostics do not implement `float80` solving or low-density driving
force interpolation. A successful probe does not qualify a device state.
## Physical constants profiles

`VELA_PHYSICAL_CONSTANTS_PROFILE` is a CMake build contract:
`sentaurus_u2022_12` (fresh-build default) or explicit `codata2018` rollback.
Existing explicit CMake cache selections remain unchanged. Compatibility uses older documented constants
for compatibility; it is not a newer CODATA standard. Use separate builds:

```powershell
cmake --preset windows-ucrt64-release-sentaurus2022
cmake --build --preset windows-ucrt64-release-sentaurus2022
./build-release/sentaurus2022/vela_example_runner.exe --physical-constants
```

`vela device --physical-constants-profile sentaurus_u2022_12` and
`python -m vela_sdevice.materials --physical-constants-profile sentaurus_u2022_12`
regenerate Formula-1 `Nc`, `Nv`, and `ni` under the selected constants, including
the existing active OldSlotboom offset. They emit `physical_constants_profile`
in materials, the material report, and generated decks. Point `VELA_BIN` at
the matching build. The CMD frontend defaults to `sentaurus_u2022_12` when its
sidecar selects `bvmethods_operator_profile: "bvmethods_tri3_qf"`. Select a
compatible runner with `VELA_BIN`, or explicitly choose `codata2018` for a
matching CODATA2018 build. Other CMD profiles also default to compatible constants; use explicit CODATA2018 for historical inputs.
The JSON field declares an expected profile; it cannot
switch constants inside a running binary. Explicit mismatches fail before
assembly. Untagged external material files retain CODATA2018 semantics and
are rejected by the compatibility binary. Regenerate them from `.par`;
adding a tag to old density values does not convert those values.

| Constant | CODATA2018 | U-2022.12 compatibility |
| --- | --- | --- |
| q (C) | 1.602176634e-19 | 1.602192e-19 |
| epsilon0 (F/m) | 8.8541878128e-12 | 8.8542e-12 |
| kB (J/K) | 1.380649e-23 | 1.380662e-23 |
| h (J s) | 6.62607015e-34 | 6.626176e-34 |
| m0 (kg) | 9.1093837015e-31 | 9.109534e-31 |
| DOS pi | mathematical pi | 3.141592654 |

DOS pi affects DOS evaluation only. Geometry, trigonometry, and hbar retain
mathematical pi. Relative dielectric parameters (Si 11.7), bandgap laws,
Formula-1 masses, BGN and acceptance gates are unchanged. Built-in Si/SiO2
have separate regenerated no-BGN 300 K snapshots. The independent nonlocal
BTBT source-constant preset remains explicit and is not silently switched;
the BV verification keeps BTBT off. Parameter provenance is recorded in the
[native constants audit](validation/bvmethods_native_parameter_constants_audit_2026-11-03.md).

## SDevice command entry point and generated v2 decks

DF-ISE Interface line regions are preserved as optional `mesh.json.interfaces`
metadata (`name`, `node_ids`, `edge_node_ids`, `adjacent_regions` with region names).
They are not volume regions/materials. See [P2 grid checks](sde_sdevice_checks.md).
Contacts spanning several volume regions additionally retain `adjacent_regions`
as sorted numeric region IDs. Their legacy `region_id` is the smallest adjacent
ID; edge/node sets remain authoritative. See [P5 modeling](sde_bv_modeling.md).

`vela device` reads a preprocessed single-device SDevice command file, imports
its grid and declared `File{Parameter}`, validates them, then emits one
`format_version: 2` JSON deck per stage. The source `.cmd` owns electrodes,
physics, goals and filenames; the optional `.vela.json` sidecar owns only the
documented numerical/diagnostic overrides. `--check` performs those gates and
writes `stage_plan.json` without starting the solver. It is separate from
passing a hand-written JSON deck to `vela_example_runner`.

See the [current workflow and executable example](sentaurus_sdevice_workflow.md#210-可执行中性示例与-sidecar),
[sidecar schema](../schemas/vela.sdevice_sidecar.v1.schema.json), and
[current package acceptance and limits](validation/sdevice_native_package_e2e_2026-10-05.md).
The CLI and exit codes are in workflow section 2; deployment targets are
Windows UCRT64 and Ubuntu 26.04 x86_64.

Other examples using `inputs/...`, `runs/...` or generic filenames below are
configuration patterns whose paths must be supplied by the caller. The linked
workflow section 2.10 provides the concrete, executed source-tree example.

Generated decks use the public v2 unit conventions; conversion to internal SI
belongs to the C++ loader. A `materials_file` is independently schema-tagged
`vela.templates_ldmos.materials.v1`: its native `_cm3`/`_cm2_V_s` keys are not
rewritten as deck keys. The frontend always supplies that file. The C++ Si/SiO2
built-ins match the library's no-BGN 300 K scalar snapshot as of WP13; see
[material alignment and baseline audit](validation/builtin_material_alignment_2026-10-05.md).
Decks without `materials_file`, and omitted fields in partial overrides, use
these defaults. Explicit material fields still take precedence.

For Formula-1 electron DOS masses, an active OldSlotboom constant bandgap
offset enters both gaps in `mt = a*(Eg(0)+dEg0)/(Eg(T)+dEg0)`, as well as the
intrinsic-density exponent. Without that model the offset is zero. The emitted
`bandgap_eV` remains the unshifted Varshni gap; the material report records the
two effective DOS gaps separately. Fresh builds and material emission default to Sentaurus U-2022.12; CODATA2018 remains an explicit rollback. See the
[native parameter audit](validation/bvmethods_material_export_2026-10-23.md).

For SDevice output, the emitter sets `sweep.vtk_format: vtu`, requested
`vtk_fields`, exact `vtk_bias_points`, `vtk_field_data` metadata and
`csv_all_contacts`. Current-target/BreakCriteria termination uses the documented
`breakdown.current_{abs,min,max}_A_per_um` plus `current_contact`. These options
are described under Sweep below. CLI `--threads` takes precedence over
`Math{NumberOfThreads}`; absence of both means one solver thread.

`vtk_field_data` goal status describes the individual stage, and exit 0 may
include an incomplete later stage when `ExitOnFailure` is absent. Verify the
run manifest as well as the exit code. Accepted CSV/PLT and VTU/PVD publication
does not imply reference-qualified fields at every bias; the PN2D BV scope is
documented in workflow section 3.3.

Explicit element-box transport and the restricted PhuMob candidate:

- `mesh_geometry.cell_box_policy: average_intersection` is an opt-in,
  unweighted 2-D Tri3 AverageBox construction. For each edge endpoint it
  intersects the perpendicular bisector with the half-planes of that
  endpoint's mesh neighbors, assigns each face to the incident cell's side
  of the edge, then averages the endpoint contributions. This is not clipping
  each local cotangent or clipping faces to the finite triangle. Non-Delaunay
  faces may vanish; neighboring coefficients and box measures can change.
  It requires `poisson_permittivity_policy: cell_material`, rejects supplied
  cell coefficients, degenerate/nonmanifold/overlapping edge fans and non-Tri3
  cells, and does not implement weighted Voronoi geometry or remeshing.
  Calculations use long double translated coordinates; stored values remain
  double. Use explicit `transport_edge_geometry: element_box` for transport.
  Historical whole-edge coefficients and mesh node volumes remain unchanged.
- `discretization.poisson_charge_volume_policy: material_local_box` and
  `solver.recombination_volume_policy: material_local_box` separately select
  the semiconductor-only sums of `M(T,i)=sum_j(g(T,ij)*length(ij)^2)/4`.
  They require cell-material geometry. Their assembled transport-node volumes
  must be positive and finite. These options align charge and SRH/Auger with
  the selected element box, including non-Delaunay adjustments. They do not
  silently override any existing volume policy or avalanche support option.
  Non-Delaunay box volume sums need not equal geometric triangle areas.
- `mesh_geometry.cell_box_policy`: `legacy_positive` (default) or
  `delaunay_transfer`. The latter requires `poisson_permittivity_policy:
  cell_material` and no supplied `poisson_cell_edge_coefficients`. It computes
  signed Tri3 cotangents in long double and transfers a negative local edge
  contribution to its same-region neighbor, preserving the signed edge sum.
  It rejects globally non-Delaunay, nonmanifold and obtuse region-boundary
  edges, subject only to a 128-double-epsilon roundoff bound. This implements
  the unweighted 2-D Delaunay subset, not the general native non-Delaunay
  intersection algorithm. Mesh node volumes and historical whole-edge
  couplings remain unchanged; explicit Poisson/transport policies consume
  the generated per-cell geometry.
- `solver.region_resolved_interface_assembly.transport_edge_geometry`:
  `inherit` (default) or `element_box`. The latter sums only transport-cell
  local coefficients times edge length. It rejects combination with the old
  `transport_edge_coupling` flag or explicit edge ratios.
- `solver.mobility.edge_averaging`: `legacy` (default), `element_box`, or
  the explicit experimental `element_box_phumob` candidate.
  The new option requires `model: constant` or `masetti`,
  `doping_concentration_basis: total_impurity`, `element_box` transport and
  `cell_material` geometry. It first averages each node's model mobility with
  `M(T,i)=sum_j(g(T,ij)*length(ij)^2)/4`, divided by the actual sum of these
  measures. Edge mobility is `sum_T(g(T,e)*mu_T)/sum_T(g(T,e))` over transport
  cells. Zero transport geometry yields zero mobility. Field-dependent,
  surface and carrier-dependent models are rejected for this averaging mode.

These two carrier controls currently require coupled Newton; Gummel's carrier
assembler rejects them. Coupled residuals, analytic Jacobian, conservative
port extraction and edge mobility probes share the new constitutive function.
For the supported doping-only models, mobility has no state derivatives.
This does not change SRH quadrature or Poisson charge volumes. The existing
`signed_transport` Poisson charge policy remains independently selected.
`simulation_type: element_box_probe` with `output_csv` exports per-cell local
coefficients, node measures in m² and electron/hole mobility in m²/(V s).

`simulation_type: phumob_source_probe` is a read-only, zero-update export for
the explicit `phumob_element_density_grad_qf` avalanche profile. It requires
`state_file` and `output_csv` and reads the saved physical fields and carrier
populations without Newton repacking or contact projection. Each Si cell,
vertex and carrier receives a row with contact adjacency, electric and QF
gradient vectors, selected ionization field, alpha, mobility, density, source
measure, generation rate and integral. Values use the active solver units;
source measure is coordinate area. Split low components are not evaluated.
The output reports `qualified: false`; successful export is not device
convergence or qualification. Unsupported source profiles fail explicitly.

For plain low-field `solver.mobility.model: phumob` with `edge_averaging:
legacy` and Boltzmann statistics, coupled Newton uses explicit carrier chain
derivatives of the PhuMob scalar equations. It differentiates the screening
functions, carrier scattering and Matthiessen combination, using zero `dG/dP`
on the existing G-floor branch. Log-density partials preserve weak populations
without subtracting two complete edge fluxes. These mobility terms are added
to fixed-mobility SG derivatives, including the doping-only BGN potential
partials and both electron/hole cross blocks. This changes the Jacobian only;
the residual, mobility values, G-floor rule and numerical acceptance are unchanged.
The density-exponent clamp boundary/outside range, Fermi statistics and
PhuMob field/surface variants retain their existing finite-difference paths.
The original `element_box` mode still rejects PhuMob. The separate
`element_box_phumob` candidate supports `model: phumob` and the explicit
element-distance Lombardi/HFS variants described below, nodal
`total_impurity`, `element_box` transport geometry and `cell_material`
permittivity geometry. Coupled Newton supports Boltzmann and Fermi statistics;
the Boltzmann density exponent must remain smooth (unclamped). Fermi density
partials include the Fermi-integral response. Avalanche requires the explicit
Si-only `phumob_element_density_grad_qf` profile below. This does not extend
the separate split/multiprecision runtime's Boltzmann-only scope.
It evaluates PhuMob from each vertex's live ND, NA, n and p, averages these
values using the same box measures above, and assembles every adjacent
cell vertex's psi, electron and hole columns conservatively. It uses the
unchanged precise G minimum; native cutoff discrepancies are not fitted away.
Full nodal carrier vectors are mandatory in the constitutive interface;
endpoint-only fallback is rejected. `element_box_probe` additionally requires
`state_file` for this mode. HDF5 states are read with the explicit mesh/reference
context and reject incompatible mesh identities. The cell probe evaluates the
stored double physical fields; it does not restore split low components. Its
JSON reports `state_evaluation: physical_double_fields`, `input_has_split_state`
and `split_low_components_used: false`. Compare with a full-split edge ledger
separately before treating the cell result as a production-state reconstruction.
Coupled residuals, SG edge diagnostics, contact
currents and the cell/edge mobility probes consume the same live populations.
The old endpoint-only `transport_edge_jacobian_probe` rejects this candidate;
use `newton_jvp_probe` to inspect its full cell stencil. Gummel remains
unsupported. This candidate does not change any global model default.

For `simulation_type: sg_edge_flux_probe`, the split Boltzmann runtime exports
carrier-specific population prefactors separately from the legacy common `ni`
columns. The electron/hole SG weights use those production prefactors. Appended
`split_sg_reference_decimal_digits=200` identifies a fresh diagnostic evaluation
of populations and direct Bernoulli products from the full packed high/low state;
production mobility is held fixed. It is an arithmetic reference, not an
independent mobility or device-qualification check. Zero denotes unavailable;
in particular, an edge without transport support has no evaluated SG argument.

Optional `sg_operand_output_json` writes an exact-input sidecar. Optional
`sg_operand_nodes` is an integer node list selecting incident edges (omitted or
empty selects all). This sidecar is supported only by the split runtime and
rejects an empty selection. Primitive JSON numbers represent binary64 inputs:
a reader must recover the binary64 value before exact higher-precision conversion,
rather than treating the printed decimal as the original real input. Packed
high/low coordinates (ordered psi0, psi1, phin0, phin1, phip0, phip1), QF
references (phin0, phin1, phip0, phip1) and the potential scale are kept separately;
mobility and 200-digit reference values are decimal strings. Concentrations,
mobilities and fluxes in this sidecar use internal units, with a concentration
to-SI factor at the top level; potentials and thermal voltage are in volts.
`production_transport_active: false` means no SG kernel was evaluated, rather
than a measured zero argument or a validated zero-current physical solution.
For split-runtime replay, use the availability-tagged `electron_sg_*`,
`hole_sg_*` fields and this sidecar. Legacy generic Fermi/edge-observer columns
that the split exporter does not populate are not authoritative operands;
a default zero in those columns must not be interpreted as a measured value.
These optional exports do not change residuals, Jacobians, state or acceptance.

Poisson permittivity assembly can be selected independently of node volumes
and carrier transport with `mesh_geometry.poisson_permittivity_policy`:

- `legacy_average` (default): retain the arithmetic material permittivity
  average multiplied by the existing whole-edge box geometry.
- `cell_material`: sum `epsilon_cell * (d/l)_cell` over adjacent Tri3 cells.
  Poisson-only, Gummel Poisson, coupled Newton residual/Jacobian, and the
  coupled electrode reaction-charge diagnostic share this coefficient.
  Without additional data, each local geometry contribution uses the existing
  box builder's cotangent/fallback policy. This does not reproduce arbitrary
  external TCAD box modifications automatically.

With `cell_material`, optional `mesh_geometry.poisson_cell_edge_coefficients`
accepts a nonempty array of records in sequential mesh cell-id order:
`{"cell_id":0,"node_ids":[0,1,2],"coefficients":[0.2,0.3,0.4]}`.
Provide exactly one record per cell. The three finite nonnegative,
dimensionless coefficients are `d/l` along local edges `(0,1)`, `(1,2)`,
`(2,0)`. Node order must exactly match that mesh cell. These values supply
geometry, not permittivity; current material parameters still supply epsilon.
Degenerate or non-Tri3 meshes, topology/value/count mismatches, and combination
with the legacy `solver.region_resolved_interface_assembly.poisson_edge_coupling`
flag in coupled Newton are rejected. Data are mesh-specific and require
independent provenance and qualification. Empty/missing data in the C++ options
select computed local geometry; an explicitly empty JSON array is rejected.
Carrier couplings, charge/SRH volumes, units, and convergence gates remain
independent. Existing decks and global defaults are unchanged.

Newton numerical and geometry controls added for explicit validation profiles:

- `solver.state_variables`: `potential_quasi_fermi` (legacy default) or explicitly
  `potential_density` (BVMethods default). Explicitly bind `coupled_operator` and `carrier_statistics`.
  The `bvmethods_tri3_qf` profile retains its 300 K Joyce–Dixon contract.
  Classical `legacy` Boltzmann/Fermi coordinates use the existing SLOT operator,
  including its configured IALMob/PhuMob, BGN, recombination and avalanche models;
  this option selects coordinates, not a replacement physical model.
  The SLOT adapter stores authoritative n/p and derives QF consumer fields.
  It retains Ohmic potential/QF constraints expressed through those populations.
  Quantum, split coordinates, thermionic contacts, BVMethods-only band-population
  SG/bootstrap, QF closure replacement/restoration and binary128 QF Poisson
  arithmetic remain unsupported and are rejected. Population protection is
  selected with `population_update`, not QF `population_step_control`.
  The QF-only `slot_population_memory` and `evaluateStep` diagnostics reject this
  mode; use the production `dc_sweep` and population update reports.
  Integration and actual qualification are recorded in the
  [SLOT integration report](validation/slot_ldmos_main_density_integration_2026-10-29.md).
  Classical Masetti and Masetti-field population configurations use the direct
  PN2D transport columns, while the SLOT adapter retains its QF consumer.
  Packed unknowns are ψ/V0, n and p in internal concentration units. Positive
  semiconductor n/p are authoritative; derived QF fields never overwrite them.
  Population DDSolution retains the original normalized psi separately from QF
  packed coordinates. Optional HDF5 `population_psi_normalized` requires
  `population_potential_schema: vela.population-potential.v1`, positive
  `population_potential_scale_V`, `population_potential_unit: V` and
  `population_potential_origin_V` matching the archive reference origin. It also
  requires the existing population internal-density bundle. Same-scale restores
  reuse normalized psi; changed scales use one physical conversion. Incomplete,
  stale or mixed QF/population coordinate bundles are rejected. Legacy archives
  remain readable through their physical fields. Reference changes require an
  explicit seed translation; changing metadata alone is invalid. The Python
  archive translator rebases population normalized psi from its authoritative
  value in extended precision, rounds once, updates the physical output and
  binds the new origin. It preserves n/p; QF packed translation remains rejected.
  Inactive oxide slots are zero and eliminated as constraints. QF caps,
  references, split coordinates, legacy row scaling and quantum/Gummel recovery
  combinations are rejected by JSON parsing as well as C++ construction. The
  population CMD emitter omits inherited QF caps/reference/trust-region options;
  explicit active QF controls fail preflight. Old merit/stall acceptance, legacy
  linear equilibration/finite-difference Jacobian selection and active legacy
  block-filter searches are not implemented by population Newton and fail instead
  of being silently ignored. Population positivity, physical residuals and
  configured report/enforce gates remain. See the
  [coordinate-policy audit](validation/potential_density_policy_audit_2026-10-29.md).
  `sd-sim --pn2d-population <cmd>` enables this path for classical silicon
  voltage stages with ohmic Anode and Cathode contacts. It retains the CMD
  physics and the BV sidecar. Coupled stages select Unit coordinates and
  relative-update stopping with CMD Digits and ErrRef values. With omitted
  ErrRef, the reference values are 0.0258 V for Poisson and 1e10 cm⁻³ for
  each carrier. The relative error uses the existing population RMS evaluator
  in physical units, with no circuit-voltage term for fixed-voltage solves.
  The RhsMin branch requires a native RHS calibration, which is unavailable
  for PN2D. Carrier-row, global-continuity and configured absolute-block
  thresholds are reported. The emitter uses full Newton steps when the CMD
  damping settings disable both damping methods within the iteration budget.
  The emitter translates normalized InitialStep,
  MinStep and MaxStep to volts and uses Increment with Newton-iteration growth
  and reciprocal Decrement on retries. `sweep.stop_below_minimum_step: true`
  stops an ordinary voltage sweep when a retry would require less than MinStep.
  The CMD `Extrapolate` setting maps to `sweep.continuation.predictor.mode`.
  Linear prediction uses two accepted states. For concentration coordinates,
  it extrapolates positive carrier concentrations in log space.
  PN2D IV and BV sweeps select `sweep.preserve_step_after_output_clip: true` to
  reuse the proposed step after a first-attempt landing at an intermediate
  `vtk_bias_points` voltage shortened the step. A rejected step uses the
  configured reduction and the accepted retry determines subsequent growth.
  Each CMD Goal and requested output bias is landed exactly. Different Newton
  outcomes can produce different adaptive trajectories. Full PN2D IV/BV and
  pseudo-arclength qualification require separate numerical evidence.
  Population-mode coupled load-line secant prediction uses linear ψ/internal
  voltage and logarithmic n/p extrapolation between accepted states:
  `c_pred = exp(log(c_current) + ratio*(log(c_current)-log(c_previous)))`.
  This predicts the next initial state; Newton solves n/C0 and p/C0 by default; physical n/p remain the public state contract.
  Zero inactive pairs remain zero, and unchanged carrier slots retain their
  exact values. Invalid domains, overflow/underflow or a failed predictor guard
  fall back to the last accepted state; no population floor is introduced.
  QF-mode prediction retains its existing behavior.
- `solver.primary_variable_coordinates`: `unit` (BVMethods population default and
  selected by the PN2D CMD mode) or `physical`. Internally psi/VT, n/C0, p/C0
  and inner voltage/VT; output, archive
  populations, positivity and ErrRef remain physical. Existing device row scales
  remain; the scalar row is divided by VT. Unit directions require one original
  defect correction and the original physical precision audit before writeback.
- `solver.compensated_population_state` enables compensated potential and
  carrier concentrations for the legacy, unit-scaled Boltzmann SG solver
  with Ohmic contacts, Masetti mobility, optional SRH, and optional
  element-edge SG/GSS-Laux avalanche. It defaults to `false`. This path
  uses physical residual ceilings of `8e-28 C/m` (Poisson) and
  `5e-22 A/µm` (each continuity equation), with backtracking on the live
  residual. A continuity row can also meet its binary64 roundoff bound
  of 16 machine epsilons times the sum of absolute flux and source terms.
  The correction uses the SG/Poisson Jacobian with mobility
  and generation/recombination coefficients frozen during linearization.
  HDF5 checkpoints store the three low components as one bundle. It is
  supported for voltage-controlled sweeps. Terminal totals use the compensated
  transport fluxes. Spatial fields in VTU retain the usual binary64 output
  representation.
- `solver.population_stopping_rule`: `native_approx` (BVMethods population default) or
  `legacy`. Approximate stop is valid E<1 OR a valid mapped RHS infinity norm
  <1e-5. E uses the physical raw direction and actual protected post-state.
  Mapping is limited to the frozen 300 K, 2600/5020 reference normalization,
  compatible constants and material scales. Outside it only E can accept, with
  an unavailable-mapping diagnostic. This is not native RHSMin equivalence.
- `solver.population_update_convergence`: optional `digits`,
  `potential_reference_V`, `electron_reference_cm3` and `hole_reference_cm3`
  for the relative update test in the fixed-voltage population solve.
- `sweep.external_circuit.physical_acceptance_mode`: `report` or `enforce`.
  Omitted selects report for native_approx, enforce otherwise. Load, carrier,
  continuity and KCL remain recorded; nonfinite states always reject. Linear
  precision gates are independent and continue to enforce their tolerances.
- `solver.population_sg_current`: `quasi_fermi` (default) or experimental
  `band_population`. The latter requires `potential_density` and its qualified
  BVMethods physics. It constructs SG arguments directly from ψ, authoritative
  n/p, effective ni, DOS and the same Joyce–Dixon statistical law, with paired
  Bernoulli evaluation and canonical node ordering. Both transport and AV
  current use that path with complete AD tangents. It does not reconstruct the
  SG carrier populations through QF. GradQF mobility and avalanche driving fields
  continue to use derived QF as required by the selected physical models.
  The ni-based band combination is algebraically equivalent to the previous
  profile; it does not add a separate material electron-affinity/energy-state
  contract. This option is an arithmetic experiment, not native zero-bias
  operator-equivalence qualification or a new global default.
- `solver.population_writeback`: `independent` (default) or opt-in
  `boltzmann_psi_compensated`. Compensation requires the classical `legacy`
  operator, Boltzmann statistics and `potential_density` coordinates. After the
  original population protection, let `e` be actual physical ψ change minus
  the intended `alpha*delta_psi`. Unprojected carrier candidates are evaluated
  as `(n+alpha*delta_n)*exp(e/Vt)` and `(p+alpha*delta_p)*exp(-e/Vt)` with
  extended intermediates, then stored in binary64. This preserves the intended
  finite Boltzmann QF change under fixed statistical prefactors. Contact
  populations use the same map and remain subject to their original QF
  boundary equations. Inactive slots and locally projected carriers retain
  the original policy; positivity, inverse consistency and acceptance checks
  remain. Zero direction is an exact identity. This is an experimental finite
  update map, not a change to the residual/Jacobian or a Fermi update formula.
- `solver.population_sg_bootstrap_linearization`: `unchanged` (default) or
  experimental `floor_adjacent_edges`. Requires BVMethods `potential_density`,
  `quasi_fermi` SG and a positive configured minimum density. Only approximate
  bootstrap matrices freeze population feedback in both the SG prefactor and
  argument on a carrier edge with either endpoint at that density floor
  (16 machine epsilons allow unit-conversion rounding). Potential derivatives,
  direct population terms, SRH, full Jacobian, physical RHS and port derivatives
  remain unchanged. Edge selection is recomputed per state and per carrier;
  full-derivative recovery restores all feedback, and an edge above the floor
  also restores it. This does not change the automatic recovery controller or
  guarantee convergence. It has no mutable cross-stage selection history.
- `solver.mobility.ialmob_transport_feedback_derivatives`: boolean, default
  `true`. Setting `false` is an opt-in quasi-Newton approximation limited to
  IALMob: omit its transport mobility feedback terms from the Jacobian.
  Mobility values, source AD caches, avalanche Jacobian, physical residual,
  contact current and stopping gates are unchanged. Keep the existing
  `jacobian_field_derivatives=true` requirement for the IALMob avalanche
  profile; this option does not bypass it or enable BVMethods-only recovery.
- `solver.population_update`: population-mode-only object. `protection` is
  `local_nonpositive` (default) or `fraction_to_boundary`. The former changes
  only free nodes with c+αδc≤0 to `local_retained_fraction*c` (default 0.01).
  The latter uniformly scales the coupled update using `fraction_tau` (default
  0.99). Both fractions must be finite and strictly between zero and one.
  Every Newton trial repeats the check; invalid predictions fall back to the
  last accepted state. This policy does not change the configured residual,
  load-line, KCL or report/enforce gates. RMS diagnostics retain the raw direction.
- `solver.population_derivative_recovery`: `after_protection` (population default), `full`, or explicit
  `automatic`. `after_protection` starts each target with Full; after an accepted
  protected writeback, only the next linearization omits mobility and AV
  derivatives. It retains the complete physical RHS, SRH derivatives and
  statistical SG edge feedback. An unprotected accepted write restores Full;
  rejected candidates do not alter the event policy. No derivative-selection
  history transfers to a new target. Automatic mode requires population coordinates, a direct bordered
  solve and `full_step` acceptance. Its bootstrap directly assembles transport
  with mobility derivatives frozen and omits AV derivatives, retaining full SRH
  derivatives and the complete physical RHS. A bounded trial/full controller
  uses the model response to the actual applied step; invalid candidates or
  confirmed model failure return to bootstrap with cooldown. Traces distinguish
  bootstrap, trial and full. A failed bootstrap bordered solve probes the full
  Jacobian once at the same state and target, without spending an update or
  changing acceptance gates. Failure of both solves returns to continuation;
  a rejected full candidate does not renew the probe until a step is applied.
  The update CSV marks this event with `bootstrap_linear_retry`.
  This is an approximate Newton policy, not a
  replacement of the physical equations or a guarantee of BV convergence.
  The population direct solve audits the original system at componentwise
  backward error 1e-10 and scaled linear residual 1e-8·max(1,‖b_scaled‖∞), using
  up to at least three correction solves when needed. Legacy QF defaults remain.
  CMD selects these options through its Vela sidecar; no Sentaurus keyword is
  invented. Population HDF5 archives retain SI output fields plus authoritative
  internal n/p coordinates and their explicit SI conversion factor. Same-unit
  restores recover the original carrier bits, while cross-unit/legacy restores
  convert the SI fields. Mesh identity and mixed QF-buffer checks still apply.
  See the [storage contract](state_archive.md).
- `solver.linear_refinement_iterations`: integer 0--10, default 0. After the
  direct Newton solve, recompute the linear defect in 100-digit arithmetic and
  solve each correction with the same row weights, configured L2 row/column
  equilibration and backend. This does not
  alter nonlinear acceptance criteria or enable refinement in other solvers.
- `solver.linear_factor_refinement_steps`: integer 0--3, default 0. Opt-in guarded
  refinement for the ordinary three-block Newton loop, separate from the SLOT
  population option. Audits the actual regularized/closure-replaced equations
  with the SLOT mixed, block and singleton-constraint linear guards; failed
  correction leaves the original direction intact and stops the linear solve.
  Uses long-double defect accumulation and the selected direct backend. For
  Woodbury closure replacement, the first correction factors the explicit
  replaced matrix; subsequent corrections reuse it. This is not reuse of the
  original unreplaced factors. Requires a separate arm from MP100 refinement,
  split state, row scaling and equilibration. Emits `VELA_LINEAR_FACTOR_AUDIT`
  records; does not change nonlinear tolerances, row/closure qualification or
  circuit validation, and does not refine an augmented external-circuit solve.
- `solver.legacy_edge_jacobian_steps`: optional object with positive finite
  `transport` (default `1e-6`) and `avalanche` (default `1e-7`). For the legacy
  edge transport and finite-difference edge avalanche source Jacobian, the
  physical column perturbation is the configured value times
  `max(1 V, abs(column potential))`. Smaller values are diagnostic controls
  requiring a same-direction residual check; they can amplify cancellation.
  They do not change the physical residual or element-box/GSS derivatives.
- `newton_step_probe.directional_audit_scale`: positive finite multiplier of
  the constrained Newton correction, or zero (default) to skip the central
  residual audit. The probe uses production row scaling, L2 equilibration,
  refinement, regularization and exact boundary updates. Its CSV includes
  raw/constrained updates, 100-digit linear defects and per-row backward
  errors, central residual changes and free continuity-row masks. These
  diagnostics do not establish nonlinear or external-load convergence.
- `solver.arithmetic_precision` selects `double` (default) or `extended80`.
  Use `vela_example_runner` for `double` and `vela_example_runner_extended80`
  for `extended80`. The extended runner is built when the compiler provides
  a 64-bit significand and an 80-bit extended exponent range in `long double`.
  Legacy mixed-precision paths still use binary64 state and linear interfaces.
- `solver.require_end_to_end_precision` defaults to `false` in native JSON.
  `sd-sim` sets it when a command file requests `ExtendedPrecision`.
  The bound SLOT LDMOS potential-density BVDS stages are identified by
  `solver.extended80_stage_profile: "slot_population_bvds_v1"`. The runner uses
  wide state, SparseLU, original-system refinement, and exact HDF5 state
  for these stages. The profile covers a zero-bias Poisson/coupled start,
  constant-seed series-resistor voltage continuation, and accepted-state
  reload. The profile is emitted by the SLOT binding. The runner checks the
  stage and its inputs before accepting or reloading a state.
  Use `--arithmetic-capabilities` to inspect the scope and backends.
  Strict requests outside validated scope are rejected before acceptance.
  The profile exports a subset of the Plot fields as Float64 display arrays
  and reports skipped fields. HDF5 remains the authoritative restart state.
  This profile enables coupled equation normalization when any normalization
  reference is supplied. For example, set
  `"scaling": {"mode":"unit_scaling","reference_concentration_cm3":"auto"}`
  or use `solver.normalization.reference_concentration_cm3` in a v2 deck.
  The solve uses `psi / V_T`, `n / C0` and `p / C0`. The reference temperature
  determines `V_T`. The concentration reference `C0` accepts a positive
  concentration in cm^-3 or `"auto"` (the larger of the maximum absolute net
  doping and the intrinsic-density fallback). The length and mobility
  references determine the Poisson and continuity row scales. Output concentrations
  use cm^-3, and the load line uses volts and A/um. The normalized HDF5 state
  stores its concentration scale and resumes with matching normalization references.
  A new normalized run starts from 0 V to produce its own accepted archive.
  The series-resistor stage accepts
  `sweep.external_circuit.coupled_rhs_growth`. For example,
  `{"mode":"report","factor":1e10}` reports RHS growth without affecting the
  outcome. With `"mode":"enforce"`, a strict exceedance returns
  `rhs_growth_limit` and the quasistationary controller retries from the last
  accepted state. Omitting the object disables the check. The factor must be
  finite and at least one. The norm is the maximum absolute augmented residual,
  including the load-line equation, with the initial trial's row scales held
  fixed. It is evaluated after each accepted Newton update and before the update
  convergence check. A zero initial norm has no growth ratio and does not
  trigger the guard. Each enabled trial logs the initial and peak norms, peak
  ratio, and the first update above the limit. This is the Vela norm definition.
- `solver.poisson_residual_precision`: `double` (default) or `binary128`.
  The latter forms physical potential and Boltzmann charge from the packed
  state and retained quasi-Fermi references in software binary128 arithmetic,
  using the same separate electron/hole reference prefactors as production
  (including `intrinsic_reference_potential_V`),
  then assembles the Poisson edge and charge terms before rounding the residual
  to double. State storage, geometry, source volumes and material coefficients
  keep their existing precision. This explicit option currently excludes
  Fermi-Dirac statistics, quantum potentials, thermionic contacts and feedback
  substitutions. Residual probes use the same assembler option as Newton.
- `solver.stable_merit_comparison`: boolean, default `false`. For `l2` or `block`
  residual norms, preserving the configured weights and resolved block scales,
  with `merit` line search and global-continuity merit `off`, compare the sign
  of the residual squared-norm difference without losing weak block changes.
  A compensated estimate falls back to exact binary integer arithmetic, or
  rational arithmetic for weighted blocks, when necessary. Existing zero-residual
  handling and convergence thresholds remain. Read-only residual probes may
  disable this unused line-search option when reporting global closure.
- `solver.split_dd_state`: boolean, default `false`. Experimental, explicit
  production integration of the qualified three-block split-coordinate state.
  Classical warm-start Newton accumulates low components for potential and both
  quasi-Fermi increments. The residual, live PhuMob mobility, SRH, conservative
  edge/row probes and terminal current consume the same state using 100 decimal
  digits before projecting public values to double. Directional residual
  differences are formed before that projection. The existing analytic double
  Jacobian, linear backend, update caps and acceptance thresholds are retained
  and require independent qualification. The PhuMob profile supports 300 K
  `phumob` with `element_box_phumob`, Boltzmann carriers and SRH. The independent
  IALMob profile supports classical 300 K edge SG, configured material-local
  SRH/Auger (including density enhancement and generation clipping), constant
  non-Si material transport and the explicit `ialmob_element_density_grad_qf`
  local-AD source contract. It preserves material electron/hole population
  prefactors separately from the recombination intrinsic density. Local psi,
  electron-QF and hole-QF differences are formed from high/low coordinates before
  projection into the existing binary64 IAL mobility/HFS/source kernels; this
  is not an arbitrary-precision rewrite of those kernels. Poisson, population,
  SG, recombination and row/port accumulation use wide arithmetic. The existing
  analytic Jacobian remains double, evaluated with reconstructed populations;
  local tiny-field derivatives require independent directional checks.
  Both profiles reject quantum, thermionic contacts, feedback substitution,
  BTBT and Gummel recovery. Other avalanche profiles remain unsupported.
  Arclength updates are not supported. Geometry and source-volume policies are
  configured separately; this option does not select a geometry bundle.
- `solver.split_contact_bias_low_V`: optional map of configured contact names to
  finite low voltage parts in volts; default empty. Requires `split_dd_state`.
  Each pair must be normalized (`high + low == high` in binary64). Fixed-bias
  Dirichlet construction uses built-in potential, applied high voltage and low
  voltage separately. Unknown contacts, Schottky boundaries and relaxed minority
  boundary targets are rejected. In a bounded circuit probe nonempty low parts
  require `split_contact_voltage`; its reported low voltage must accompany the
  archived node state during independent qualification. This does not add split
  contact-voltage support to arclength or production circuit continuation.
- `simulation_type: "poisson_arithmetic_model_probe"` with `state_file` and a
  fresh `output_model_json` exports the actual Boltzmann Poisson edge couplings,
  material electron/hole reference populations, charge volumes, packed high
  coordinates and Dirichlet equations. The `vela.poisson-arithmetic-model.v1`
  node tuple is `(psi,n,p,xpsi,xn,xp,potentialScale,eref,href,niElectron,niHole,vt,
  doping,voln,volp,vold,q,area,scale,interfaceRhs)`; edge tuples are `(i,j,g)` and
  Dirichlet tuples `(i,value)` use normalized potential. `zero_low_residual`
  strings retain the binary128 arithmetic result to 40 significant digits;
  `zero_low_residual_double` records the actual production binary64 projection,
  including its Dirichlet subtraction. Do not reconstruct binary64 halfway
  rounding decisions from the finite decimal strings. This read-only export enables
  independent fixed-QF hi/lo Poisson diagnostics with IALMob/Auger/avalanche
  configurations, but does not implement split transport, enable Newton
  acceptance, or bypass the restrictions of `split_dd_state`. Fermi, quantum
  and thermionic configurations are rejected. Density conversion to SI is
  recorded separately as `density_internal_to_m3`.
- `solver.exact_dirichlet_updates`: boolean, default `false`. After the direct
  solve and optional refinement, set the step of each Dirichlet identity row
  to minus that row's residual. Free-row updates are unchanged before the
  existing update caps and line search. Nonzero `carrier_regularization_scale`
  is rejected for this option. It does not add a new convergence criterion.

Warm restarts subtract the old and new quasi-Fermi references before adding
the stored increment. This preserves sub-reference-rounding increments and
does not change the CSV state format.

With `poisson_residual_precision: binary128`, Newton additionally retains its
original normalized coordinates in the optional restart CSV columns
`packed_psi`, `packed_electron_qf_increment`, `packed_hole_qf_increment`, and
`packed_potential_scale_V`. These four columns require the existing QF reference
and increment columns. Loading validates their consistency with physical fields;
Newton and residual probes reuse them when the potential scale and references
match. Changed physical fields invalidate the hint; the writer omits stale hints.
Ordinary restart files and their subnormal-to-zero convention remain supported.
Packed hints are omitted for subnormal coordinates or physical increments.

Split-state CSV checkpoints additionally require `packed_psi_low`,
`packed_electron_qf_increment_low`, `packed_hole_qf_increment_low`,
`split_state_schema` (`vela.split-dd-state.v1`) and `split_mesh_fingerprint`.
These columns must appear together with the packed coordinates and QF reference
fields. Low components retain finite subnormal values. A stale or incomplete
split payload is rejected rather than omitted. Resume requires the same mesh,
potential scale and QF reference frame, with `split_dd_state` enabled. Physical
columns remain rounded projections; the complete packed payload defines the
state. Split-mode contact diagnostic fields named `*LongDoubleReference` project
the same wide evaluator and do not constitute an independent reference.

- `solver.region_resolved_interface_assembly.poisson_charge_node_volume`:
  `inherit` (default) or `signed_transport`. The latter assigns signed
  transport-material AverageBox volume to all three Poisson charge terms
  (electrons, holes, dopants). Dielectric edge coupling and continuity/SRH
  volume are unchanged. Cannot combine with the per-term Poisson volume flags.
- `solver.region_resolved_interface_assembly.transport_edge_coupling_ratios`:
  optional array of finite nonnegative multipliers in mesh edge-id order,
  exactly one per edge. Coupled assembly, terminal integration, and mobility
  probes use the same effective coupling. Cannot combine with
  `transport_edge_coupling`. The array is mesh-specific; an empty array keeps
  the existing policy. There is no global mesh-volume or model default change.

The two explicit geometry controls require the coupled Newton assembler;
passing them directly to the Gummel DD assembler raises an error. Their use
does not by itself qualify arbitrary meshes or signed negative node measures.

This document is the implementation-aligned reference for JSON config files used
by Vela. It describes fields currently parsed by the C++ Poisson, DC sweep, and
single-bias Newton paths.

Use this file as the field-level reference. Use
[reference_tcad/README.md](../reference_tcad/README.md) for validated fixture status and
[architecture.md](architecture.md) for solver path boundaries.

Scope and conventions:
- No `scaling` field keeps the legacy SI input behavior used by existing decks.
- In legacy SI mode, all numeric values are SI unless a field name includes a unit suffix such as `_eV`.
- The optional public input-unit mode is `scaling.mode: "unit_scaling"`; it is a unit interpretation and numeric scaling foundation only, not a calibration feature.
- Relative paths are resolved from the directory of the config JSON file.
- Legacy decks remain supported where noted.
- Prototype features are marked explicitly.

Frozen candidate conservation diagnostic: `simulation_type:
"newton_candidate_conservation_probe"` reads `state_file` and writes a JSON
audit to `output_csv` (the generic runner output-path key). It requires
`solver.global_continuity_closure.mode: "enforce"` and uses that policy unchanged.
The saved quasi-Fermi reference/increment coordinates are retained at contacts;
physical carrier equations are evaluated before Dirichlet replacement. Contact
sets come from the configured electron/hole boundary conditions. Flux, integrated
source and mismatch use assembler continuity units; their ratio is dimensionless.
The process succeeds when the diagnostic executes, even if `satisfied` is false.
`device_qualified` is always false: this probe neither solves nor certifies the
full device. The SLOT campaign opt-in `--candidate-conservation-runner PATH`
records this audit before poly recovery and still requires strict device
reclosure, all active rows, global closure and load-line confirmation afterwards.
Its probe executable and protocol are bound into the campaign identity.

Classical DC assembly reuse: `solver.reuse_jacobian_structure`
(boolean, default `true`) shares an immutable zero-valued sparse pattern, offset
map and scatter tables between sequential Newton point services on the same
`DCSweep` instance. It does not share numerical Jacobian values, residuals or
carrier states. An exact compatibility key checks node numbering, geometry,
transport support/stencils, material statistics, temperature, cell stencil mode
and constrained rows. A prepared-input cache miss or an explicit context clear
also invalidates the cache. `jacobian.structure_cache_hits/misses` and
`jacobian.structure_cache_check` report reuse and compatibility-check cost;
`jacobian.pattern_build_calls` counts actual builds. Set the option to `false`
to clear the shared structure and rebuild independently. The historical
`solver.diagnostic_reuse_jacobian_structure` alias remains valid when the
production option is absent; an explicit production option takes precedence.
This cache is separate
from linear-solver symbolic-analysis reuse and does not cover the independent
four-equation electrothermal service. It requires a persistent DC worker to
reuse across separate point requests.

`solver.diagnostic_fermi_node_cache` (boolean, default `false`) selects a bounded
per-node/per-carrier endpoint-statistics cache within each classical Jacobian
assembly. Neighboring edges can share identical Fermi trial evaluations. Keys
compare the exact potential, local/relative QF, reference origin and material
logarithm; all entries expire at the end of the assembly. Flux formulas and
finite-difference derivatives are unchanged. The counter
`jacobian.fermi_node_cache_assemblies` confirms execution; existing endpoint hit
and miss counters include this cache when enabled. It is independent of structure
reuse and is an experimental option, not the production default.
- Field names with historical SI suffixes are kept for compatibility. In
  `unit_scaling` mode the numeric interpretation is described explicitly below.

## Top-level fields

| Field | Type | Required | Notes |
| --- | --- | --- | --- |
| simulation_type | string | Optional | Common values: `poisson`, `dc_sweep`, `newton`, `newton_solve_from_state`. The selected CLI/tool path determines which fields are consumed. |
| mesh_file | string | Yes | Input mesh JSON path. |
| materials_file | string | Optional | Optional material override file. Supported shapes: top-level array, object with `materials` array, or object map keyed by material name. |
| output_csv | string | Optional | Default CSV output path for DC sweep. Can be overridden by `sweep.csv_file`. |
| output_vtk | string | Poisson: Yes | Poisson VTK output file path. Optional for `newton` and `newton_solve_from_state`. |
| output_vtk_prefix | string | Optional | Default VTK prefix for DC sweep point outputs. Can be overridden by `sweep.vtk_prefix`. |
| current_contact | string | Optional for `newton_solve_from_state` | Explicit configured terminal for `current_total_A_per_um`; invalid or empty names fail. Output includes `current_contact` and `current_contact_selection`. Omission preserves legacy largest absolute **effective boundary potential** selection, which may select a metal gate and is not a drain-current contract. This top-level field is separate from `sweep.current_contact`. |
| output_state_file | string | Optional | Restart-state CSV written by `newton_solve_from_state`. Uses Vela restart format: `node_id,psi,phin,phip,electrons_m3,holes_m3`. |
| state_fields_dir | string | Required for `newton_solve_from_state` and probe-style external-state tools | Directory containing Sentaurus-export-style scalar field CSVs. |
| runtime_log | object | Optional | Runtime log control. Default is enabled and writes `<config stem>.log` next to the config file. |
| scaling | object | Optional | Legacy input unit interpretation. Omit for legacy SI behavior, or set `mode` to `unit_scaling`; see below. Not accepted together with `format_version: 2`. |
| format_version | integer | Optional | Deck format version. The only supported value is `2`, which interprets the whole deck in the TCAD internal units and forbids `scaling`; see below. |
| mesh_geometry | object | Optional | Mesh box-geometry options; see below. |
| doping | array | Yes | Region doping definitions; see below. |
| regions | array | Optional | Region-level fixed charge definitions; see below. |
| interfaces | array | Optional | Interface sheet/fixed/trap charge definitions; see below. |
| contacts | array | Yes | Contact bias and type definitions; see below. |

### Isothermal material energy reference

An optional material field `intrinsic_reference_potential_V` defaults to zero.
At 300 K the coupled/Newton/Gummel transport paths use
`n = ni * exp((psi - phin - delta)/Vt)` and
`p = ni * exp((phip - psi + delta)/Vt)` in the Boltzmann case, with the
corresponding density-of-states prefactors for Fermi statistics. Physical `ni`
in SRH/Auger is unchanged. Ohmic neutrality uses
`psi_contact = bias + equilibrium_potential + delta`; quasi-Fermi potentials
remain physical terminal voltages. The fixed material offset is differentiated
through carrier density, charge, flux, mobility and recombination.

This option represents a material energy reference, not a fitted contact bias.
Derive relative offsets from consistent affinity, band gap, density of states
and a common energy origin; do not infer them from a target current. Nonzero
offsets require 300 K, reject shared semiconductor nodes with incompatible
references, and are not supported with quantum correction or Schottky contacts.
Avalanche with nonzero offsets is restricted to the audited element IALMob
GSS-Laux configuration where offsets occur only in non-Si materials outside
its source support. General heterojunction, electrothermal and nonzero-Si
avalanche-reference extensions require separate work.
### Built-in materials

Without `materials_file`, Si/SiO2 use the generated no-BGN **300 K** snapshot
[`default_materials.json`](../share/vela/materials/default_materials.json).
Si has `ni=10750176656.471613 cm^-3`, `Nc=2.8566791451911786e19 cm^-3`,
`Nv=3.104657033455896e19 cm^-3`, `Eg=1.1241592307692307 eV`, `chi=4.05 eV`,
`mun=1417` and `mup=470.5 cm^2/(V*s)`, with `eps_r=11.7`.
SiO2 has `eps_r=3.9`, `Eg=9 eV`, `chi=0.9 eV` and zero carrier density/mobility.
Legacy SI storage converts concentration by `1e6` and mobility by `1e-4`.

Material overrides inherit omitted fields from these built-ins. Pin every
required scalar explicitly to reproduce an older calibration; see the
[WP13 audit and frozen baselines](validation/builtin_material_alignment_2026-10-05.md).
The existing C++ temperature law holds the gap fixed, scales `Nc/Nv` by
`(T/300)^1.5` and mobility by `(T/300)^-2.2`; equivalence to `.par` emission is
qualified at 300 K only. Alignment does not activate additional models.

### Explicit electrothermal DC entry

The four-equation point input also accepts `diagnostic_electrothermal_cost`
(`off`, default; `heat`; `screening`; or `linear`). This is opt-in wall-clock
instrumentation only. `heat` splits the thermal matrix preparation, derivative
scatter, and finalization specifically during residual-only coupled requests;
`screening` measures complete screening-minimum calls; `linear` separates
solver creation, actual symbolic analysis, and numerical factorization.
Results are in `subcost_diagnostics` with seconds and call counts. These are
nested costs, not independent additions to existing assembly/factorization
totals. Per-call clock overhead must be measured, especially for thermal
scatter and screening; diagnostic times are not promises of removable cost.

The four-equation entry has a separate `reuse_jacobian_structure` option
(default `true`), accepted in the prepared point input or at sweep
deck level (the deck wins). It retains a zero-valued coupled/thermal sparse
pattern and checked scatter positions in the sweep preparation context. An
enabled deck retains that context even when `reuse_static_preparation` is absent
or explicitly `false`. To disable all sweep-context caches, set `reuse_jacobian_structure`,
`reuse_static_preparation`, and `reuse_linear_analysis` to `false`. Linear-analysis
reuse also retains immutable preparation so its exact input identity can be checked.
Exact topology/geometry, IALMob support, and constrained-row identities guard
reuse; ordinary bias and unknown-temperature changes refill values. Contact
row type changes rebuild the coupled pattern. Residuals, all physical values,
temperature partials and numeric factorizations are still evaluated. A new
process starts fresh; explicit disabling drops the shared structural cache.
This entry has its own electrothermal qualification; it does not change the
classical DC implementation. See the [electrothermal structure validation](validation/templates_ldmos_thermal_structure_2026-09-24.md) for correctness and timing limits.

`simulation_type: "electrothermal_dc_sweep"` selects the in-process four-equation
silicon service. This is a separate explicit schema; the legacy DC fields above
do not select or configure it. Its initial scope is the audited LDMOS SI model,
not arbitrary semiconductor materials. Integrated numerical qualification is
tracked in [current LDMOS validation](validation/templates_ldmos_current_status.md).
The [R7 production run guide](validation/templates_ldmos_production_reproduction.md)
provides a versioned, hash-checked input export and complete execution commands.

| Field | Meaning |
|---|---|
| `input_file` | Prepared four-equation SI JSON: mesh, native Poisson/transport/volume geometry, doping, audited material configuration, thermodes and contact boundaries. Relative to this deck. |
| `output_directory` | New output directory, relative to this deck. Existing directories require explicit `resume: true`. |
| `sweep.bias_points_V` | Strictly increasing drain biases beginning at exactly zero. Every requested point must pass the numerical, terminal KCL and heat-balance gates. |
| `sweep.initial_step_V`, `minimum_step_V`, `maximum_step_V` | Defaults 0.1, 0.0001 and 4/3 V. Failed attempts halve the attempted increment. |
| `sweep.max_newton`, `growth_newton` | Defaults 60 and 8. Accepted steps with at most `growth_newton` updates grow by 1.5; more than20 updates halve the next step. |
| `sweep.predictor` | `linear` (default), `none`, or the experimental guarded options below. Linear prediction uses accepted states only, referenced QF increments and an extrapolation ratio at most3. |
| `sweep.predictor_guard` | Experimental `none` (default) or `residual`. Screen a constant-state fallback against the linear prediction under one Jacobian row scale and the original electrical block/row limits. Record all candidate preparation costs; not a convergence criterion. |
| `sweep.predictor_history` | Experimental `legacy` (default) searches older accepted history when needed; `nearest` only uses the nearest accepted predecessor and declines extrapolation when its step ratio exceeds3. |
| `sweep.predictor="tangent_guarded"` | Experimental alternative to default `linear`. Freshly assemble/factor the accepted source-state Jacobian, solve `J dx/dV=-F_V`, then screen the tangent candidate against the original secant/constant seed using common residual scaling and electrical row/block protection. Additional preparation, residual calls and factorization are counted. Failed preparation retains the original seed. |
| `sweep.predictor="local_guarded"` | Experimental anchored local fit, only when the nearest accepted-history step ratio exceeds3. Uses three/four accepted states, normalized bias, QR rank checking, quadratic-to-linear downgrade with weight L1 bound8, and common residual/row/block screening. No physical tolerance changes. |
| `sweep.step_policy` | Experimental `adaptive` (default), `fixed_targets`, or `actual_step`. Fixed-target mode directly attempts each requested bias and saves/stops on failure. Actual step growth uses the successfully executed advance when an exact output target clips a planned step; original startup, growth factor, bounds and failure halving remain. |
| `sweep.density_update_maximum_bias_V` | Optional finite nonnegative drain-bias ceiling. Above it, disables the point input's density-coordinate update. Omitted: no ceiling. Does not change the bias grid or convergence gates. |
| `sweep.density_update_requires_prediction` | Boolean, default `false`. If enabled, density-coordinate updates require an actually used accepted-state predictor; initialization/prebias, zero bias and the first unpredicted drain step use the original QF update. |
| `initialization.mode` | `provided_state` (default) uses the prepared explicit state. `neutral_300K` constructs a neutral initial guess, solves initial Poisson/coupled equilibrium, performs Poisson gate prebias, then restores coupled equations. |
| `initialization.gate_voltage_V` | Required positive final gate voltage for `neutral_300K`; gate boundary work-function/reference offsets are retained from the prepared input. Prebias passes through4 V when the final voltage is higher. |
| `initialization.max_newton` | Default100 for initialization/prebias point solves. |
| `resume` | Reload an accepted checkpoint; reject changed input, mesh, sweep or initialization configuration. |
| `pause_after_attempts` | Optional positive count for a controlled checkpoint pause during drain continuation. A file named `STOP` in the output directory also pauses between attempts. |
| `reuse_static_preparation` | Default false. Retains immutable mesh, edge/volume mapping, doping, heat geometry and IALMob interface geometry across initialization and drain point services within this sweep. Exact serialized configuration and complete mesh/IALMob source bytes determine reuse, including units, crystal axes and contact geometry. Each point creates a new assembler and state caches; changed carrier state or temperature never reuses an old root or transport state. A resumed process creates a fresh context. |

The prepared point input may explicitly set `diagnostic_ngmres_recovery: true`
(default false). This experimental option allows one residual-space recovery
after a stagnation-watch or line-search failure, with at most four prior states
and four trial residual evaluations. Histories from other QF references are
converted and reassembled before mixing; bias/model changes start a new point
and cannot share history. Original voltage/temperature update caps and electrical
row/block guards remain; a candidate must reduce the common scaled norm by more
than1%. This progress safeguard is not a relaxed convergence tolerance. Failed
recovery preserves the original failure state. `ngmres_recovery.updates` is
reported separately from `newton_updates`; total state changes include both.

The source input and mesh are preserved as JSON snapshots. Each attempt records
its full input, solver history and output (referenced state, temperature, residuals,
terminal currents, lattice source and thermode heat). `ledger.json` atomically
records accepted states, rejected attempts, exact points and elapsed wall time.
`curve.csv` uses A/µm for terminal current and K for temperature;
`terminal_balance.csv` records all four terminals. Source and boundary heat are
in W/m (unit device width). Failed attempts never replace the last accepted state.
Incomplete initialization currently requires a new output directory; completed
initialization and drain-continuation checkpoints can be resumed.

The prepared `input_file` can explicitly enable these performance controls (all
default to `false`): `defer_recentered_candidate_jacobian` computes residuals only
for line-search candidates that must be reassembled after QF recentering, and
`skip_equilibrium_poisson_transport` omits zero edge currents during Poisson
prebias only when both QF fields and temperature are globally constant and all
carrier/temperature rows are constrained. Coupled transport and convergence
gates are unchanged. Curve summaries are retained in exact-point ledger entries
to avoid repeatedly parsing all prior full-state files; old checkpoints retain
the full-state fallback.

`reuse_ialmob_local_preparation` reuses field-independent low-field mobility
terms across cell vertices with identical model, doping, carrier density,
interface distance and temperature, within one transport-state preparation.
`residual_ialmob_values_only` omits spatial and temperature mobility derivatives
for residual-only candidates. A subsequent full Jacobian request rebuilds any
incomplete cached state. Both require the IALMob path to have an effect; neither
changes physical parameters, step selection or convergence criteria.

`reuse_ialmob_thermal_high_field` (default false) reuses high-field mobility values
and local partials between spatial and temperature directions of the same element
evaluation. Temperature chain derivatives are still applied separately. It never
shares high-field results between elements or states. `diagnostic_ialmob_kernel_timing`
(default false) times value-only, spatial and temperature passes, in that order
in `performance.ialmob_pass_seconds`; element clocks perturb timing, so use a
separate run for this diagnosis. Pass counts, high-field evaluations/reuses,
screening cache requests/hits and local preparation hits/builds are also recorded.
These point-scoped counters include preparation and final outputs, whereas the
existing solve-loop counters retain their narrower scope.

`solver.mobility.ialmob.screening_method` (isothermal DC) and
`mobility_SI.ialmob.screening_method` (electrothermal input) select the screening
minimum root algorithm. Omission uses safeguarded `halley`, including direct C++
construction. Explicit `legacy` retains the original algorithm; `newton` and
`toms748` remain comparison alternatives. The electrothermal top-level
`diagnostic_ialmob_screening_method` remains an explicit alias; conflicting
formal and alias declarations are rejected. The alternatives solve the same
stationarity equation in a logarithmic variable with an analytic bracket,
finite/sign checks and fallback to the qualified legacy path. No temperature
quantization or cross-temperature root reuse is permitted. The algorithm is part
of both static preparation identity and exact screening-cache keys. Results record
`ialmob_screening_method`, `ialmob_screening_candidate_calls`,
`ialmob_screening_function_evaluations` and `ialmob_screening_fallbacks` in
`performance`. Candidate function counts cover transformed-equation evaluations;
they do not count legacy derivative evaluations and must not be compared as
equivalent work. These four root algorithms do not change the physical model or
acceptance gates. Comparison runs must explicitly select `legacy` as their baseline instead
of relying on an omitted option.

The explicit `u2022_12_single_newton` selection is a separate, experimental
numerical constitutive compatibility policy, **disabled by omission**. It uses
the version-bound, one-step screening cutoff observed in Sentaurus U-2022.12
binary64, rather than solving for a stationary minimum. Runtime value identity
has been checked at 300 K for Si masses 1 and 1.258; the implementation permits
only those masses and 250--500 K. Other temperatures in that interval are
internal consistency checks, not a claim of native validation. Unsupported
inputs fail explicitly; there is no fallback to another cutoff convention.
Both 3-D Coulomb screening and the 2-D accumulation term use the selected G.
The reference differentiated kernel includes the moving-cutoff temperature
chain rule. The generated low-field kernel is not supported with this policy
and rejects that combination. The splice is continuous in value; because the
one-step cutoff is not stationary, its two branch slopes need not agree.
At equality the clamped branch is selected. Check one-sided derivatives rather
than interpreting a central difference across the splice as a smooth Jacobian.
This opt-in choice does not change any solver acceptance gate or an existing
running configuration. Frozen value and derivative evidence is documented in
[the screening compatibility validation](validation/slot_ldmos_bvds_ial_screening_compatibility_2026-10-17.md).

`diagnostic_ialmob_explicit_high_field` and
`diagnostic_ialmob_generated_low_field` (both default false) select independent
experimental mobility kernels in the electrothermal point service. The first
uses explicit high-field partials and a value-only residual path; the second
uses offline SymPy-generated low-field values and seven SI partials. The original
AD kernels remain the reference. Both switches participate in static preparation
identity; generated and reference low-field caches are separate. They do not
change the screening root, branch conventions, physical model or solver gates,
and remain disabled in the R10 rollback profile. The explicitly selected,
qualified R11 profile enables only the high-field switch; generated low-field
remains an experimental candidate. See the
[R11 qualification](validation/templates_ldmos_symbolic_full_2026-09-17.md).

Preparation timing separates identity checking, mesh loading, input mapping,
IALMob geometry and new assembler construction. `point_preparation_seconds` is
inclusive; its components must not be added to it. Assembly seconds are split
into full-Jacobian and residual-only requests; nested IALMob timers overlap them.

The shared isothermal linear solver also accepts the process environment override
`VELA_LINEAR_SOLVER=umfpack` when built with detected UMFPACK support. Without an
override, its preference order is `umfpack`, `sparselu`, `strumpack`, `mumps`,
`superlu_mt`, excluding unavailable optional libraries. SparseLU is always
available, so it is the default when UMFPACK is absent. Direct backends reuse exact sparse patterns and identical
numeric factors; unavailable UMFPACK and failed solves raise errors rather than
silently changing backends. This override is independent of the electrothermal
configuration below; it does not select a different physical model.

`VELA_LINEAR_SOLVER=sparselu_long_double` explicitly selects an experimental
Eigen SparseLU factorisation and solve in `long double`. The assembled matrix
and RHS remain binary64 and the returned direction is rounded back to binary64;
this is mixed precision, not extended-precision residual assembly or state
storage. Platforms where `long double` has no more significand bits than `double`
are rejected. The backend does not cache factors and does not change the default
backend preference or nonlinear acceptance criteria.

On Windows and Linux, `VELA_LINEAR_SOLVER=pardiso_mkl` opts into the oneMKL
PARDISO runtime backend. Set `VELA_MKL_RT_LIBRARY` to an absolute runtime
library path (`mkl_rt.3.dll` on Windows or `libmkl_rt.so` on Linux).
The existing `VELA_MKL_RT_DLL` variable remains a fallback when the new variable
is unset or empty. Dependent runtime libraries must be available to the loader;
on Windows use the process PATH. Linux links the system dynamic-loader library,
not MKL itself. No oneMKL library
is loaded for other backends; this option does not enter the default preference
list. The interface is LP64, real nonsymmetric binary64, with METIS ordering,
scaling, matching and up to two iterative-refinement steps. It reuses exact
patterns and identical factors, and reports errors without backend fallback.
`VELA_LINEAR_THREADS` controls its thread-local MKL setting. Selection applies
to shared auxiliary solves as well as the Newton system; assembly, state
precision, nonlinear updates and acceptance criteria remain unchanged.

`sparselu_eliminated`, `sparselu_eliminated_long_double`, `umfpack_eliminated`
and `pardiso_mkl_eliminated` are explicit
experimental shared-solver backends. They eliminate exact numerical diagonal
singleton rows using their actual diagonal (including regularization), subtract
the constrained-column contribution from the free RHS, and lift the solution.
All compute affine RHS correction in long double, then round it to binary64.
The UMFPACK and PARDISO variants use the corresponding binary64 backend on
the reduced system and require the same package/runtime as that backend.
Only the long-double SparseLU variant changes factorization/solve precision. Explicit stored
zeros do not count as couplings. These backends do not cache factors, change
assembly or acceptance, or become defaults. The selection applies to all shared
linear solves, including auxiliary Poisson and Woodbury RHS solves.
The oneMKL PARDISO adapter is compiled by default on Windows and Linux. Set
`-DVELA_ENABLE_MKL_PARDISO=OFF` to exclude the adapter, mock runtime libraries
and dynamic-loader link dependency. Existing CMake caches retain their setting.
Explicit selection in a disabled build raises an error naming the CMake option.

In an enabled build, `VELA_LINEAR_SOLVER=pardiso_mkl` explicitly selects the
runtime backend. Set `VELA_MKL_RT_LIBRARY` to the absolute path
of the installed MKL runtime; `VELA_MKL_RT_DLL` is a compatibility fallback
when that variable is unset or empty. The runtime's dependent libraries must
also be discoverable (on Windows, include their directories in the process
`PATH`). No MKL headers or link-time libraries are required, and other backends
do not load MKL. Missing libraries, symbols and failed solves raise errors
without switching backends. The default preference order is unchanged.

This backend uses LP64, real nonsymmetric binary64, METIS ordering, scaling,
matching and up to two iterative-refinement steps. `VELA_LINEAR_THREADS`
accepts `1` (default), `2`, or `4`; each call restores the previous thread-local
MKL setting. Its reported thread count is the requested setting, not a measured
worker count. Exact sparse patterns and identical numerical factors are reused.
Selection applies to shared auxiliary solves too; it does not select the
separate electrothermal backend or change assembly, nonlinear updates or
acceptance criteria. This is not extended-precision solving or a guarantee of
Sentaurus-equivalent results. The core `LinearSolver` now eliminates prescribed
constraints by default before invoking the backend and recovers the full direction
afterward. `VELA_CONSTRAINT_ELIMINATION=0` disables this reduction for controls;
`1` (or unset) enables it. Ordinary calls identify numerical diagonal prescribed
rows, ignoring stored zeros without a magnitude cutoff. Coupled load-line and
population-coordinate calls retain the row IDs from the device Jacobian before
adding the circuit border or changing coordinates, so affine contact dependencies
are also eliminated. The selected block must have nonzero diagonal and acyclic
dependencies; unsupported explicit metadata raises an error. The full RHS
correction and affine back-substitution are retained. Original residual/update
acceptance gates are unchanged. Explicit QR and emergency alternate factorization
paths retain their existing policies; changing backend alone does not prove
Sentaurus matrix agreement.
An explicit `VELA_LINEAR_SOLVER=pardiso_mkl` also selects PARDISO for the full
`direct_bordered` load-line matrix and its refinement corrections, including
the direct path of `auto`; a PARDISO error is returned without backend fallback.

For the `bvmethods_tri3_qf` coupled sweep with a directly bordered PARDISO
solve, the original-system defect calculation uses compensated EFT
accumulation. Finalization uses the selected arithmetic precision.
Out-of-range rows use exact rational arithmetic in `double` and widened
arithmetic in `extended80`. A detected discarded-tail cancellation is
recomputed from the exact binary64 inputs using rational arithmetic.
The weighted absolute defect, componentwise backward error and original
correction rules retain their existing tolerances. Floating-point
environments unsupported by EFT use the original MP100 path.
Set `VELA_BVMETHODS_EFT_MP100_FINISH=1` to
compare against EFT with MP100 finalization, or `VELA_BVMETHODS_MP100_DEFECT=1`
to select the original all-MP100 defect path.
The explicitly selected `direct_bordered_qr` policy retains its QR preference.
In an embedding process, MKL must permit LP64 initialization before other MKL
calls; preinitialized NumPy/MKL and Python API compatibility need separate
qualification. Runtime paths belong in local environment configuration.

BVMethods device and terminal transport Jacobians reuse identical long-double
Lombardi powers within each element, including the approximate device
Jacobian. Coupled residual and terminal-current assembly also reuse these
powers in potential-density mode with quasi-Fermi SG.
Set `VELA_BVMETHODS_DISABLE_LOMBARDI_POW_REUSE=1` before launch to compare
against uncached evaluation. The cache is scoped to one element evaluation.

BVMethods potential-density mode with quasi-Fermi SG uses scalar Tri3
transport by default when only the transport value is requested. Set
`VELA_BVMETHODS_TRI3_VALUE_MODE=dual` before launch to use the Dual path for
those calls, or set it to `shadow` to return the Dual result and compare
both value vectors bit for bit with the scalar result. Jacobian requests
and unsupported state or SG modes use the Dual path.

Scalar Tri3 transport computes the Canali mobility value in long double.
Set `VELA_BVMETHODS_CANALI_VALUE_ONLY=0` before launch to use the full Canali
evaluation for scalar transport, or set `VELA_BVMETHODS_CANALI_VALUE_AUDIT=1`
to compute both paths and check the mobility bit pattern at every evaluation.
BVMethods avalanche residual assembly also computes only the Canali mobility
value when it does not assemble avalanche derivatives. Full Canali evaluations
in BVMethods transport and avalanche assembly reuse the common long-double power
for the mobility and its derivatives.
Set `VELA_BVMETHODS_DISABLE_CANALI_AV_VALUE=1` to restore full Canali evaluation
for avalanche residuals, or `VELA_BVMETHODS_DISABLE_CANALI_SHARED_POWER=1` to
restore separate power evaluations. `VELA_BVMETHODS_CANALI_OPTIMIZATION_AUDIT=1`
compares the active optimized evaluations with the original bit patterns and
checks the assembled avalanche source.

BVMethods potential-density mode with quasi-Fermi SG also uses scalar avalanche
residual assembly by default when no avalanche derivatives are requested.
The scalar device transport passes the two carrier mobilities at each
element's representative vertex to avalanche assembly for the same state.
Set `VELA_BVMETHODS_DISABLE_AV_MOBILITY_REUSE=1` before launch to recompute
those mobilities. Set `VELA_BVMETHODS_AV_VALUE_MODE=dual` to restore dual
avalanche residual assembly, or `shadow` to return the dual result and
check scalar source terms, generation rates, zero-current diagnostics and
reused mobilities bit for bit. Derivative requests and cells outside the
scalar arithmetic range use dual assembly. Calls without eligible scalar
device transport compute avalanche mobility independently.

BVMethods potential-density contact and circuit current evaluations reuse the
last source-free transport vector when the packed state matches bit for bit
within the same assembler runtime. Set `VELA_BVMETHODS_TRI3_PORT_CACHE=0`
before launch to disable this reuse. Set it to `shadow` to recompute each cache
hit and check bitwise equality. Jacobian assembly uses its own path.

Optional shared-solver candidates also accept `sparselu_metis`, `umfpack_metis`,
`mumps`, `mumps_metis`, `superlu_mt`, `superlu_mt_metis`, and `strumpack`.
They require their detected CMake packages (`VELA_ENABLE_METIS`,
`VELA_ENABLE_MUMPS`, `VELA_ENABLE_SUPERLU_MT`, `VELA_ENABLE_STRUMPACK`).
An enabled option alone does not prove availability. These candidates do not
change the default and do not extend the separate electrothermal selector.
METIS uses the undirected structural union, not numerical symmetrization.
MUMPS uses its non-MPI OpenMP double-precision library and nonsymmetric mode;
SuperLU_MT uses COLAMD or an external METIS column ordering with internal
postordering; STRUMPACK uses METIS, maximum diagonal-product matching/scaling,
uncompressed DIRECT mode, and disables tiny-pivot replacement and GPU execution.
`VELA_LINEAR_THREADS` accepts `1` (default), `2`, or `4` for those three threaded
adapters. Control BLAS threads separately; the installed OpenMP OpenBLAS ignores
`OPENBLAS_NUM_THREADS`, so that variable alone does not isolate BLAS threads.
The isolated acceleration study uses `openblas_set_num_threads` and records the
runtime query. The shared runner also offers the opt-in experiment
`VELA_ENABLE_OPENBLAS_THREAD_CONTROL=ON` (CMake default OFF) with
`VELA_BLAS_THREADS=1`, `2`, or `4`. It requires OpenBLAS, verifies the reported
thread count after every solve, and records one `linear.openblas_threads`
observation per successful solve, including fully prescribed systems. Constructor
checks do not add observations. Other
requested values are rejected. MUMPS, SuperLU_MT and STRUMPACK additionally record
their OpenMP thread configuration and restrict active nesting to one level.
Without the environment switch the shared runner does not set BLAS threads.
UMFPACK's sparse algorithm remains serial; its dense kernels may use threaded
BLAS. For this installed OpenMP BLAS, also allow the requested OpenMP team size
with `OMP_NUM_THREADS`; the API query alone does not prove every kernel launches
that team. The production Eigen SparseLU path has no BLAS or parallel build
switch enabled by these runtime variables. Performance qualification is
configuration-specific; no thread setting changes a default physical model.
See the [backend screening record](validation/templates_ldmos_sparse_backend_execution_2026-09-18.md)
for tested ranges. SuperLU_MT at 2/4 threads encountered replay accuracy failures
in that historical screening; small-system unit tests do not override them.
The later [T470p thread pilot](validation/templates_ldmos_t470p_threads_2026-09-20.md)
passed its three matrix rounds and bounded D5 curves, a separate environment-specific
qualification rather than a general guarantee for every multithreaded configuration.

The default `DCSweep` and CLI `vela_example_runner --dc-worker` retain the
main Newton linear context across points and sequential JSON-lines DC requests.
`--dc-worker-linear-reuse` remains an alias; `--dc-worker-no-linear-reuse` is the
explicit control that retains prepared inputs only. Each new request invalidates numeric
factors; changed prepared inputs, sparse structure, or failed/invalid requests
invalidate analysis as appropriate. Changing the selected backend creates a new
context. STRUMPACK also retains its value-dependent
matching/scaling until reanalysis. Auxiliary Poisson/recovery solvers keep their
existing lifetimes. No nonlinear state is implicitly reused. Use
`DCSweep(false, false)` to disable both caches explicitly. The portable
[LDMOS reference entry](../reference_tcad/templates_ldmos_sentaurus2022/validation_d5/README.md)
defaults to UMFPACK and reuse. Historical experiment drivers keep explicit control
flags and reject a resume with changed cache policy. See the
[cross-point experiment](validation/templates_ldmos_cross_point_analysis_2026-09-19.md).

The sequential worker loads explicit HDF5 state files. The former experimental
`initial_state_vds_hex` transport is rejected. State reuse and sparse symbolic
analysis reuse are independent policies.

`VELA_LINEAR_CAPTURE_DIR` is an opt-in diagnostic directory, which must already
exist. It records at most four large three-block solver inputs per process as
`VELALU02` files and rejects existing output names. These files contain the matrix,
RHS and solution in solver coordinates, not the unscaled physical Jacobian.
Capture runs must be excluded from performance comparisons. Leave unset normally.

`VELA_LINEAR_FACTOR_STATISTICS=0` disables the optional numeric-factor statistics
and SparseLU structural work estimate in this shared solver. Basic stage timers,
solve/cache counters and numerical checks remain active when profiling is enabled.
Unset, empty, or `1` retains the existing statistics; other values are rejected.
The switch is read at solver construction and does not change factorization controls.
Use it for paired end-to-end timing; gather detailed factor diagnostics separately.

`electrothermal_linear_solver` selects `umfpack`, `sparselu_colamd`, or
`sparselu_amd` for this four-equation service only. Omission prefers `umfpack`
when compiled in, otherwise `sparselu_colamd`. `umfpack`
requires a build with detected SuiteSparse UMFPACK; unavailable or unknown
selections fail explicitly. Unless the explicit adaptive-Jacobian experiment is enabled,
all choices repeat numerical factorization and reuse symbolic analysis only for identical compressed sparse indices when
enabled. Backend changes can alter Newton trajectories through roundoff and
require independent numerical and performance qualification.

`reuse_linear_analysis` defaults to `true`. A sweep retains its direct solver
across accepted points; a standalone point needs an explicit preparation context
to share it. The point's `reuse_sparselu_symbolic` option now defaults to `true`
and applies to all supported backends. Explicit `false` forces analysis on each
factorization and prevents cross-point retention. Deck-level values for these
options and `electrothermal_linear_solver` override the prepared point input.
Reuse checks dimensions and every compressed sparse index, as well as the backend
(including ordering) and immutable preparation identity. Changed coefficients
always receive fresh numerical factors unless the separate adaptive-Jacobian
experiment requests lagging. Failed/throwing solves or rejected sweep gates drop
the linear context. Restart starts with a fresh context; no factors are serialized. The checkpoint
records the effective linear policy and rejects a changed or missing policy,
including on an already completed scan.
`performance.symbolic_analyses` is a per-point delta, including separate tangent
preparation, not a lifetime total. `linear_object_reused` reports object retention,
not a guarantee of a symbolic hit when the matrix pattern changes.


`diagnostic_density_update_iterations` is a nonnegative integer in the point
input, default `0` (disabled). During at most this many initial coupled Newton
iterations, and only while the scaled residual exceeds `1e-6`, it transforms
the existing Jacobian direction to nodal density increments using the full
psi/QF/temperature derivatives. It bounds negative increments to retain positive
densities, inverts the same local Fermi/DOS/BGN model at the candidate psi/T,
and applies the original residual line search. Eight unsuccessful density
trials fall back to the original QF direction. Dirichlet carrier unknowns,
Poisson-only prebias and final convergence gates retain their existing behavior.
The sweep guards above can further restrict its use. Broad-bias use has failed
representative LDMOS controls; do not infer qualification from enabling this flag.

Additional prepared-point experimental controls, all disabled by default:

| Field | Behavior |
| --- | --- |
| `diagnostic_iteration_trace` | Records before/after block and row gates, direction maxima, node temperatures, local density changes, and each ordinary line-search candidate's alpha, merit, acceptance, Jacobian request and inclusive trial time. Does not participate in numerical decisions; `iteration_trace_seconds` reports gate/direction diagnostic overhead. Trial durations include actual trial work and are not additional overhead. |
| `diagnostic_density_projection` | `off` (default), `v1`, or `v2`. Coupled silicon carrier updates use per-node positive density targets with relative floor 0.01 and absolute floor `min(old_density,1e-250 m^-3)`. V1 applies throughout coupled iterations; V2 returns to QF updates when every active raw QF direction is below 0.01 local Vt. Poisson prebias remains unchanged. Unlike the earlier R8 window, these variants deliberately do not require a predictor or low bias. |
| `diagnostic_natural_damping` | NLEQ_ERR-type corrector test using the existing LU and fixed scaling, with extra back-solves counted in `natural_damping`. Preserves the original near-floor behavior and terminal convergence gates; does not implement the full NLEQ_ERR predictor/corrector algorithm. |
| `diagnostic_adaptive_jacobian` | Reuses a factorization for at most two updates following an accepted full step with residual ratio at most 0.1. QF recentering forces refresh; failed stale-Jacobian attempts retry fresh within the existing iteration budget. Candidate residuals are evaluated before preparing another Jacobian. |
| `diagnostic_pseudo_transient` | Adds `d(0,-q*A_cont*n,+q*A_cont*p,0)/dx / tau` to the coupled steady Jacobian, using `recombination_area_m2` (or its documented silicon-area fallback). Constrained rows have no storage; a finite-exchange hole contact remains active. Initial tau is the median positive finite carrier `abs(M_ii/J_ii)` times `diagnostic_pseudo_time_scale` (default 1). SER uses initial steady row scales throughout, with growth factor clamped to [0.5,4], tau bounded to [1e-6,1e12] times its initial value. A failed line search reduces tau tenfold, with at most three retries within the original iteration budget. Poisson prebias ignores this option. Incompatible with natural damping, adaptive Jacobian, and NGMRES. Steady acceptance gates are unchanged; this is not physical transient qualification. |
| `diagnostic_pseudo_direction_audit` | Requires `diagnostic_pseudo_transient` and exactly one coupled diagnostic iteration. Independently compares full, frozen mobility-derivative, frozen bulk-recombination-derivative, and combined directions with/without carrier mass. Reports true steady merit slopes and projected trial residuals. Coefficient values, residuals and contact laws remain unchanged; audited directions never replace the actual update. |
| `diagnostic_pseudo_acceptance` | `steady` (default), `defect_ser`, or `defect_model`. The latter two require coupled pseudo transient and check the actual backward-Euler defect `G(y)=F(y)+(rho(y)-rho(x))/tau` with initial fixed row weights and norm Armijo coefficient 1e-4. The previous accepted state x and tau stay fixed during every backtracking trial. `defect_ser` retains the steady-residual SER time controller; `defect_model` multiplies tau by 2/0.5/1 when `norm(W*G-(1-alpha)*W*F)/(alpha*norm(W*F))` is <=0.25/>0.75/otherwise. The stagnation watch uses the within-step defect reduction when accepted this way. The original near-floor trial rule and all steady terminal gates are retained. Poisson prebias is unchanged. |
| `diagnostic_near_steady_qf_switch` | Default false. Independent candidate requiring referenced coupled `provided_state`, PTC, V1 projection and `defect_model`, without direction audit. Once the original scaled merit is <1e-9 and all three original electrical block gates pass, rebase physical QFs within +/-1 mV to zero references, rebuild the full residual/Jacobian and row scales, clear stagnation history, and permanently retire carrier mass and density mapping for the rest of this point. Records the event and inclusive preparation time. Does not accept a state, reset the attempt budget, relax any terminal gate, or change the production default. |
| `diagnostic_near_steady_qf_rebase` | Default false. Standalone R7 representation experiment, requiring referenced coupled `provided_state` QF Newton, isolated from PTC, density mapping, natural damping, adaptive Jacobian, NGMRES and the PTC switch. If merit <1e-9 and all original electrical blocks pass while the row gate still fails, once per point move physical QFs within +/-1 mV to zero references and fully reassemble. Retains original terminal gates, stagnation history, update limits and iteration budget. Reports changed references and inclusive reassembly time; no production default change. |
| `diagnostic_near_steady_contact_consistency` | Default false. Isolated referenced coupled `provided_state` QF Newton only; cannot combine with density mapping, PTC, natural damping, adaptive Jacobian, NGMRES or either near-steady QF experiment. Once per point, when merit <1e-9 and all original electrical blocks pass while the row gate fails, set finite-hole-contact potentials to the existing neutral target at their current T. QFs, references and T remain fixed. Try a full residual/Jacobian/scaling rebuild only for a nonzero potential correction <=1e-8 V. Accept as a terminal repair only if original merit, all row and electrical block gates pass; otherwise discard the candidate and retain ordinary Newton state, scales, history and budget. Reports the event, correction and inclusive preparation/reassembly time (overlaps assembly profiling). Adds no Newton update/factorization; the repair and its assembly are reported separately. Outer current/heat/state qualification remains required; production default unchanged. |
| `diagnostic_local_qf_limiter` | Default false. Isolated coupled `provided_state` QF experiment; clips each raw carrier direction to the existing voltage update cap, while global alpha is limited by potential and temperature only. Retains real-residual backtracking and original terminal gates. Incompatible with density mapping, PTC, natural damping, adaptive Jacobian, NGMRES and near-steady QF experiments; compatible with contact consistency. The sweep disables it throughout initialization/prebias. Eight-point LDMOS screening failed three points; not a recommended configuration. |
| `reuse_neutral_contact_roots` | Default false. Exact per-node memoization of the neutral-contact potential and its temperature partial, keyed by bias (including zero sign), temperature, donors and acceptors. Any unequal key, including one temperature ULP, recomputes the original root. Scope is one point assembler, with no cross-point or checkpoint cache. Profiling records actual root solves and hits in `neutral_root_counts`; this changes computation reuse only, not equations or tolerances. |
| `diagnostic_neutral_root_newton` | Default false. Uses the existing charge derivative in a bracket-protected scalar Newton solve for each neutral target. Preserves the original bias +/-1.5 V bracket and strict charge comparison. Newton proposals must be finite, have positive derivative, lie strictly inside the bracket and move no more than half its current width; otherwise bisect. A rounded zero step may probe one adjacent float, without repeated neighbor walking. Terminates only at an unsplittable floating-point bracket; after 70 unsuccessful trials or when root spacing is below the original interval width / 2^70, recomputes the original 70-bisection target. Temperature derivatives use the original implicit formula at the returned root. Does not alter the outer Newton algorithm or gates; applies during initialization too. `neutral_root_iteration_counts` records density evaluations, Newton proposals, bisections, neighbor probes and full legacy fallbacks, including all fallback work. Qualification is configuration-specific. |
| `diagnostic_hole_row_audit_nodes` | Optional list of distinct active silicon hole-row IDs, including finite hole contacts; only allowed for zero-update coupled `provided_state` diagnostics. Preserves the raw input by bypassing load-time finite-contact potential projection and QF recentering. Records every oriented hole edge current, frozen endpoint statistics/state/mobility, SRH/Auger rates, source area charge, residual and absolute flux sum. Finite contacts additionally record the outward current, neutral target and temperature derivative, small-delta operands and Fermi derivatives. Independently checks that observation leaves the full residual and hole flux sums identical; this extra observation assembly is outside solve-loop profiling counters. It is a read-only arithmetic audit, not a convergence or model option. |

For an explicit `electrothermal_dc_sweep`, the contact-consistency candidate is
disabled during neutral initialization and every gate-prebias stage, including
coupled reclosure. The requested flag is restored for drain continuation and is
part of the input snapshot checked on resume. Direct point-service calls retain
the isolation guards above.

Projection covers free carrier unknowns only, retaining Dirichlet/contact
constraints and the candidate-temperature statistics inversion. Invalid targets
are rejected; eight unsuccessful mapped trials fall back to QF updates. These
controls have not qualified as production improvements; see the
[Newton update execution report](validation/templates_ldmos_newton_update_execution_2026-09-15.md).
The mass-matrix controls and coefficient-isolation findings are recorded in the
[pseudo-transient validation report](validation/templates_ldmos_pseudo_transient_2026-09-15.md).

`diagnostic_voltage_update_limit_V` is a positive finite point-input value,
default `0.2`. It bounds voltage corrections in the original global damping;
the temperature bound remains 30 K. Legacy density-coordinate carriers instead
use the global density positivity bound; the explicit local projection variants
project each free carrier independently. This is a numerical update limit, separate
from both the outer drain-bias step and the convergence tolerances. Increasing
it has not improved the tested LDMOS representative-point totals.

These solver gates are independent of reference acceptance: native electrical
curves, thermal fields, the user-approved 2% local-density RMS and same-environment
performance still need separate scoring. The legacy `electrothermal_probe` calls
the same point-solve service for compatibility and diagnostics.
| boundaries | array | Optional | Explicit Poisson boundary segments (Neumann/insulating/symmetry); see below. |
| solver | object | Optional | Gummel/Newton settings for DD sweep and Newton solve. |
| sweep | object | Required for dc_sweep | Sweep mode, range, outputs, and diagnostics controls. |

## Simulation type dispatch

`vela_example_runner --config <file>` dispatches by `simulation_type`:

- omitted or `poisson`: Poisson driver.
- `dc_sweep`: adaptive curve sweep driver.
- `newton`: single-bias coupled Newton solve.
- `newton_solve_from_state`: single-bias coupled Newton solve initialized from `state_fields_dir`. The directory must contain `ElectrostaticPotential_region0.csv`, `eQuasiFermiPotential_region0.csv`, and `hQuasiFermiPotential_region0.csv`; each file must include `node_id` and `component0` columns with every mesh node exactly once.
- `newton_jacobian_block_probe`: read `state_file` and write the requested
  Jacobian `blocks` to `output_csv`. `finite_difference_step` sets the
  relative perturbation. The optional `finite_difference_mode` is
  `double_symmetric` by default; the audit-only
  `multiprecision_branch_resolved` mode evaluates the canonical
  element-edge avalanche source with a 50-decimal-digit scalar and requires a
  caller-selected step below every nonzero active-branch margin.

For `dc_sweep`, omitting `solver.method` keeps the default Gummel path. Set
`solver.method` (or legacy alias `solver.type`) to `newton` to use the coupled
Newton sweep path where supported. Use `poisson_only` (alias `aba_poisson`) for
Sentaurus-style ABA: the first point is a coupled equilibrium solve; later
points solve only nonlinear Poisson after reconstructing electron and hole
quasi-Fermi fields with Vela's nearest-majority-contact-basin approximation.

Runtime log CLI override:
- `--log auto`: keep config/default behavior.
- `--log off`: disable runtime log for this run.
- `--log <path>`: force-enable and write to the given path.

Optional runtime log profile CLI override:
- `--log-profile minimal|default|debug`

## format_version

`format_version` declares the deck format. The only supported value is `2`.

```json
{ "format_version": 2 }
```

A `format_version: 2` deck expresses every value in the TCAD internal units:
`um`, `cm^-3`, `cm^-2`, `cm^2/(V s)`, `cm/s`, `V/cm`, `cm^-1`, `cm/V`,
`cm^6/s`, `V`, `K`, `s`, and `eV`. This is the same numeric interpretation as
`scaling.mode: "unit_scaling"`, expressed as a deck-format declaration instead
of a solver option.

Rules:

- `format_version: 2` and a `scaling` block together are rejected. The unit
  mode is implied by the version, and `scaling` no longer has any remaining
  role in a v2 deck.
- Any `format_version` other than `2`, or a non-integer value, is rejected.
- A v2 deck uses the v2 key names. Every key whose name still claims an SI unit
  is renamed, and the old name is rejected rather than aliased; see below.
- Omitting `format_version` keeps the legacy behavior described under
  `scaling`: no `scaling` block means SI input, and
  `scaling.mode: "unit_scaling"` means TCAD input. This transitional form is
  scheduled for removal.

### v2 key names

A v2 deck renames the keys whose names claim an SI unit while their value is
already read in the TCAD internal units. The rename is by unit suffix
(`_m3` to `_cm3`, `_m2` to `_cm2`, `_m2_V_s` to `_cm2_V_s`, `_m_s` to `_cm_s`,
`_m_inv` to `_cm_inv`, `_m_per_s` to `_cm_per_s`, `_m_per_V` to `_cm_per_V`,
`_m6_per_s` to `_cm6_per_s`, `_V_per_m` and `_V_m` to `_V_per_cm`), plus the
unsuffixed names `ni`, `mun`, `mup`, `donors`, `acceptors`, `contact_radius`,
`surface_recombination_velocity` and `depth_m`, which become `ni_cm3`,
`mun_cm2_V_s`, `mup_cm2_V_s`, `donors_cm3`, `acceptors_cm3`,
`contact_radius_um`, `surface_recombination_velocity_cm_per_s` and `depth_um`.
The v1 spelling of any renamed key is rejected in a v2 deck with a message
naming its replacement. A materials file loaded by a v2 deck follows the same
contract.

Two blocks keep their names: `solver.band_to_band`, whose parameters are
genuinely SI, and `solver.normalization`, which is already spelled in v2 units.
Values are unchanged by the rename: a v2 deck is numerically identical to the
same deck with `scaling.mode: "unit_scaling"`.

The equation-normalization references that used to live inside `scaling` move
to `solver.normalization` in a v2 deck:

```json
{
  "format_version": 2,
  "solver": {
    "normalization": {
      "characteristic_length_um": "auto",
      "reference_concentration_cm3": "auto",
      "reference_mobility_cm2_V_s": "auto"
    }
  }
}
```

Each key accepts `"auto"` (the default, derived from the mesh, doping, and
material data) or a positive number. `solver.normalization` requires
`format_version: 2`; supplying it in a legacy deck is rejected. These values
control the internal non-dimensionalization of the equations and are
independent of the input unit system.

## scaling

The `scaling` object selects how numeric input values are interpreted before
internal solver scaling. No `scaling` field means legacy SI input behavior:
lengths use `m`, concentrations use `m^-3`, mobilities use `m^2/(V s)`,
electric fields use `V/m`, interface sheet densities use `m^-2`, voltages use
`V`, temperatures use `K`, and energies use `eV`.

The only supported public mode name is:

```json
"scaling": { "mode": "unit_scaling" }
```

In `unit_scaling` mode, input values use common external TCAD units:

| Quantity | legacy SI input | `unit_scaling` input |
| --- | --- | --- |
| length | m | um |
| concentration | m^-3 | cm^-3 |
| mobility | m^2/(V s) | cm^2/(V s) |
| velocity | m/s | cm/s |
| electric field | V/m | V/cm |
| sheet density | m^-2 | cm^-2 |
| voltage | V | V |
| temperature | K | K |
| energy | eV | eV |

`scaling.mode` does not accept `si` or non-public aliases.
Field names with explicit unit suffixes remain the stable public names; this
mode defines how their numeric values are read at the schema boundary.
The interpretation is applied while reading mesh coordinates, material override
files, doping, region fixed charge, interface sheet/trap charge, mobility
settings, and electric-field-related solver settings. In `unit_scaling`, these
values remain in TCAD internal units (um, cm^-3, cm^-2, cm^2/(V s),
cm/s, V/cm, cm^-1, cm/V) and the assemblers apply named composite factors where
geometry and charge/current units meet.

Poisson driver note:
- In legacy mode (no `scaling`), Poisson uses the historical SI assembly path.
- In `unit_scaling` mode, Poisson uses a scaled (dimensionless) assembly path
  and restores physical potential before producing outputs.
- VTK `potential_V` remains the physical potential in volts in both modes.

CSV output keeps the legacy column names for compatibility. In `unit_scaling`
mode, DC sweep CSV output also appends convenience columns using common
external TCAD display units, such as `current_total_A_per_um`,
`current_electron_A_per_um`, `current_hole_A_per_um`, `charge_C_per_um`,
`capacitance_F_per_um`, and `max_electric_field_V_per_cm`.

DC sweep CSVs include solver provenance columns:
- `solver_method`: selected path, such as `gummel`, `newton`, `poisson_only`,
  `frozen_state`, or `gummel_newton`
- `gummel_iterations`: iterations used by the Gummel stage for this bias point
- `newton_iterations`: iterations used by the coupled Newton stage for this bias point
- `handoff_stage`: final accepted stage or failure stage, such as `newton`, `gummel_failed`, `newton_failed`, or `gummel_fallback`

## mesh_geometry

The optional `mesh_geometry` object controls box-method geometry derived from
the input mesh. Omit it to keep the legacy default.

```json
"mesh_geometry": {
  "node_volume_policy": "mixed_voronoi",
  "fallback_negative_cotangent": true
}
```

Supported `node_volume_policy` values:
- `barycentric` (default): each Tri3 cell contributes one third of its area to
  each adjacent node.
- `mixed_voronoi`: acute and right triangles use the cotangent Voronoi vertex
  area; obtuse triangles assign half the area to the obtuse vertex and one
  quarter to each other vertex. Edge couplings still use the existing
  cotangent box coefficients and negative-cotangent fallback policy.

`fallback_negative_cotangent` defaults to `true`. For an obtuse-cell local
cotangent contribution, `true` substitutes the existing positive barycentric
box contribution. `false` removes that local contribution (zero weight); it
does not retain a signed negative edge coupling. The option is intended for
explicit mesh/discretization qualification.

SDevice mesh qualification and the box builder use the same per-triangle
threshold. A cotangent with magnitude at most `1.7453292519943296e-12`
(an angle within `1e-10` degrees of a right angle) contributes zero to
the opposite edge. Other triangles sharing that edge keep their own
contributions. Mixed-Voronoi shares in a triangle with a rounded
cotangent are rescaled to sum to its area. A cotangent below the negative
threshold still triggers `require_non_obtuse` or the configured fallback
when that guard is disabled.

The node volume policy controls the areas used by Poisson charge,
carrier continuity volume terms, stored/terminal charge, and SG avalanche
source density normalization. It does not change donor/acceptor values used for
Ohmic contact boundary conditions.

The template-private `templates_ldmos_external_averagebox` carrier-transport
profile is an explicit, default-off diagnostic. It replaces only the SG carrier
edge couple; Poisson geometry and source volumes remain unchanged:

```json
"mesh_geometry": {
  "node_volume_policy": "barycentric",
  "carrier_transport_couple_profile":
    "templates_ldmos_external_averagebox",
  "external_averagebox_couples_file": "transport_couples.csv",
  "external_averagebox_expected_edges": 16237
}
```

The CSV header must be exactly `node0,node1,couple_m`; every selected mesh edge
must be listed exactly once, including selected zero-couple records. The record
count must equal `external_averagebox_expected_edges`; unlisted edges retain
their mesh coupling. The profile is qualified
only for the exact Templates/LDMOS classical G3 diagnostic with impact
ionization, surface mobility/IALMob, quantum potential, and predictor disabled.
It requires barycentric node volumes, is not a production/global default, and
does not inherit the PN2D `element_edge_sg_gss_laux` bundle.

The explicit `slot_ldmos_external_averagebox` diagnostic reads the same CSV
contract while preserving the selected mesh node volumes (including
`mixed_voronoi`). This permits a SLOT transport-couple-only comparison and
surface/element IALMob. Quantum potential remains rejected. Avalanche is allowed
only with IALMob, `van_overstraeten`,
`current_approximation: "ialmob_element_density_grad_qf"`, and
`source_jacobian: "local_ad"`; all existing kernel-specific constraints still
apply, including explicit observed source geometry and complete mobility
derivatives. Both `postprocess_only` and `self_consistent` use the same contract
in probes and DC sweeps. The caller must supply SLOT's own native geometry.
For a full native transport comparison, use matching native element edge
coefficients **and** vertex measures in `mobility.ialmob.geometry_file`;
changing only the edge-couple CSV does not align mobility averaging. The IALMob
loader checks coefficient sums against the active carrier couples.
It changes no Poisson couples or source volumes, and makes no curve
qualification claim. The Templates profile retains its original restrictions.

## discretization.poisson_charge_volume_policy

The optional Poisson charge-volume policy is independent of the pure mesh
geometry contract:

```json
"discretization": {
  "poisson_charge_volume_policy": "material_local"
}
```

Supported values are `global` (default) and `material_local`. `global` uses
the mesh node volume. `material_local` sums one-third of each adjacent Tri3
transport-cell area and uses that volume only for the Poisson mobile-carrier
and ionized-dopant charge terms, including their Jacobian entries. It does not
change continuity source volumes, fixed/interface charge integration, BTBT,
SRH/Auger, impact ionization, stored charge, or the carrier SG couple.

The material-local policy is explicit and default-off. Its first qualified
implementation requires `mesh_geometry.node_volume_policy: "barycentric"`;
mixed-Voronoi composition is rejected rather than inferred. A mesh containing
only transport material is returned to the existing node-volume vector so the
single-material path remains bitwise unchanged.

`material_local_mixed_voronoi` is a separate, explicit policy requiring
`mesh_geometry.node_volume_policy: "mixed_voronoi"`. It sums the same local
mixed-Voronoi shares as the mesh builder, restricted to transport cells, for
Poisson mobile and dopant charge and their derivatives. It preserves
continuity sources and edge couples. This policy is useful for isolating
material ownership on SLOT without also changing the geometric volume rule;
it does not import Sentaurus AverageBox measures.

`material_local_circumcentric` is an opt-in diagnostic that sums signed Tri3
circumcentric vertex shares over transport materials. It does not replace obtuse
shares by the mixed-Voronoi half/quarter rule. Every assembled transport-node
measure must be positive and finite; invalid geometry is rejected, not clipped.
The policy computes its own charge measures even on an all-semiconductor mesh;
mesh node volumes, carrier edge couples and avalanche measures remain unchanged.
It is independent of the dielectric coupling choice and is not a claim of exact
Sentaurus AverageBox equivalence on non-Delaunay geometry.

`solver.recombination_volume_policy: "material_local_circumcentric"` selects the
same material-local measure separately for SRH/Auger and their derivatives.
It does not implicitly follow the Poisson volume policy; select each explicitly
when conducting a controlled comparison. Existing defaults are unchanged.

## discretization.poisson_permittivity_coupling

`edge_mean` (default) preserves the existing adjacent-cell arithmetic mean
permittivity times the configured edge couple when the mesh permittivity policy
is not `cell_material`. With `mesh_geometry.poisson_permittivity_policy:
"cell_material"`, that default path instead sums each cell's permittivity times
its configured box edge coefficient, including AverageBox when selected.
For an independent SLOT cold-state example, see
`configs/regression/slot_ldmos_idvd_average_box_aligned_control.json`.
`per_cell` is an opt-in diagnostic
for the coupled Newton assembler and its probes, including its Poisson-block
initialization. It requires nondegenerate Tri3 cells and sums signed
`epsilon_cell * cot(opposite_angle) / 2` contributions. Negative contributions
are retained; the mesh edge fallback is not applied to this dielectric operator.
Transport couples, charge volumes and source terms are unchanged. Pure Poisson
simulation and Gummel/DDAssembler reject the option rather than silently using
a different operator. This option does not claim equivalence to native AverageBox
or change the PN2D template contract.

The same conductance is used in the coupled residual, its full Jacobian and
`newton_poisson_term_probe`. That probe appends `dielectric_operand_sum`, the
sum of `abs(conductance * potential_scale) * (abs(x_i) + abs(x_j))` in the same
row normalization as `dielectric_flux`. Use this operand scale, together with
the absolute charge terms, when checking roundoff in nearly cancelling fluxes;
the net flux alone does not represent the arithmetic operand magnitudes.

## Doping, regions, interfaces

### doping[] entries

Required per entry:
- region: string
- donors: number (m^-3)
- acceptors: number (m^-3)

Optional per entry:
- fixed_charge_m3: number (signed elementary-charge density in m^-3)

Notes:
- Doping uses net doping `Nd - Na` per region.
- `fixed_charge_m3` may be specified either in `doping[]` or `regions[]` for a region, but not both.
- With `scaling.mode: "unit_scaling"`, `donors`, `acceptors`, and
  `fixed_charge_m3` numeric inputs are read and kept internally as `cm^-3`.

### node_doping_file

Optional top-level field:
- node_doping_file: string

`node_doping_file` points to a CSV with columns
`node_id,donors_cm3,acceptors_cm3`. When present, it overrides region-average
`doping[]` entries for drift-diffusion DC sweeps. The file must contain exactly
one row for each mesh node id; missing rows, duplicate rows, invalid ids,
quoted fields, malformed concentrations, and non-finite numeric values are
rejected.

With `scaling.mode: "unit_scaling"`, the donor and acceptor concentrations in
the CSV use the same external concentration convention as deck-level doping:
`cm^-3` and are kept internally as `cm^-3`. In legacy SI mode, the same values are
read through the legacy concentration path.

`node_doping_decimal_input: true` opts an authoritative decimal CSV into
extended-precision parsing in the `extended80` runner. The Poisson charge
assembly subtracts donor and acceptor concentrations before rounding the
resulting charge contribution. The default remains false. This option is
for decimal source data. A CSV formatted from binary64 TDR data cannot
restore digits that were absent from the source. Mobility, band-edge models
and contact population targets still use the existing concentration view.

### regions[] entries

Common fields:
- name: string
- fixed_charge_m3: number (optional)

Notes:
- Region `fixed_charge_m3` is additive source charge for Poisson and for drift-diffusion Poisson substeps.
- Duplicate definitions for the same region are rejected.
- With `scaling.mode: "unit_scaling"`, `fixed_charge_m3` numeric inputs are
  read and kept internally as `cm^-3`.

### interfaces[] entries

Region pair selectors:
- regions: [regionA, regionB] (preferred)
- or region0 + region1 (legacy-compatible form)

Charge fields (all optional, each in m^-2 and in units of elementary charge):
- sheet_charge_m2
- fixed_charge_m2
- trap_density_m2
- trap_occupancy (0..1; requires `trap_density_m2`)

Notes:
- Effective interface sheet charge is:
  `sheet_charge_m2 + fixed_charge_m2 + trap_density_m2 * trap_occupancy`
- Charge is distributed edge-by-edge across the requested region pair.
- `interfaces[]` is consumed by the standalone Poisson driver and by drift-diffusion DC sweep Poisson substeps (`solver.method: gummel` and `solver.method: newton`).
- `trap_density_m2` is signed in this prototype. Positive occupied traps contribute positive sheet charge; negative values contribute negative sheet charge. Use `fixed_charge_m2` for bias-independent fixed charge when you do not want the value scaled by `trap_occupancy`.
- The trap model is a quasi-static prototype: `trap_occupancy` is a fixed user-supplied constant for the whole run/sweep. Bias-dependent trap occupancy, frequency dispersion, and trap statistics are not implemented.
- With `scaling.mode: "unit_scaling"`, `sheet_charge_m2`,
  `fixed_charge_m2`, and `trap_density_m2` numeric inputs are read as
  `cm^-2` and kept internally as `cm^-2`.

## contacts[]

Required per entry:
- name: string
- bias: number (V)

Optional fields:
- type: string
- flatband_voltage: number (V)
- work_function_eV: number (eV)
- barrier_eV: number (eV)
- electron_barrier_eV: number (eV)
- hole_barrier_eV: number (eV)
- surface_recombination_velocity_m_per_s: number (m/s)
- surface_recombination_velocity: number (legacy alias)
- electron_surface_recombination_velocity_m_per_s: number
- hole_surface_recombination_velocity_m_per_s: number
- emission_model: string (prototype Schottky selector)

Recognized type values (case-insensitive, `-`/`_` normalized):
- ohmic
- dirichlet
- metal_gate (alias: gate)
- schottky
- floating

Compatibility policy:
- If `type` is omitted, behavior is exactly legacy-compatible and treated as `ohmic`.

Current implementation status:
- Poisson driver: `ohmic`, `dirichlet`, `metal_gate`, `schottky` supported; `floating` rejected.
- DC sweep (Gummel): `ohmic` and `schottky` with
  `emission_model: "dirichlet_barrier"` supported.
- DC sweep (Newton): `schottky` is supported only with
  `emission_model: "thermionic_robin"`, a positive electron velocity, and a
  positive hole velocity. Poisson-only and Gummel-Newton sweeps reject
  thermionic contacts.

Prototype notes:
- The Gummel `dirichlet_barrier` path fixes carrier densities and remains an
  engineering smoke model.
- The Newton `thermionic_robin` path pins only electrostatic potential at the
  Schottky contact. Carrier quasi-Fermi potentials remain unknown, and the
  continuity rows receive `sn*(n-n0)` and `-sp*(p-p0)` boundary fluxes.
- In `unit_scaling` mode, all three velocity fields are interpreted in cm/s;
  without it they are interpreted in m/s. The `_m_per_s` spelling is retained
  for schema compatibility.

## boundaries[]

Required per entry:
- name: string
- type: string
- node_ids: array of node ids, length >= 2

Optional per entry:
- normal_displacement_C_per_m2: number (for `type: neumann`)

Recognized boundary type values:
- neumann
- insulating
- symmetry
- dirichlet (parsed but currently rejected in boundaries path; use contacts for Dirichlet)

Behavior:
- `neumann` uses `normal_displacement_C_per_m2` (default 0.0 if omitted).
- `insulating` and `symmetry` are equivalent to zero-Neumann.
- Unknown type, short `node_ids`, or non-finite Neumann values are rejected.

Unit interpretation:
- The field name `normal_displacement_C_per_m2` is stable for compatibility.
- Legacy mode reads numeric values as SI displacement (`C/m^2`).
- `unit_scaling` mode reads numeric values in common TCAD area units
  (`C/cm^2`) and assembles it with the active internal area convention.

About node_ids or edge_node_pairs:
- `node_ids` polyline is the implemented path.
- `edge_node_pairs` is not currently parsed by the C++ drivers and should be treated as reserved for future schema extensions.

About normal_displacement_C_per_m2 / normal_electric_field_V_per_m:
- `normal_displacement_C_per_m2` is implemented for Neumann Poisson boundaries.
- `normal_electric_field_V_per_m` is not currently parsed; treat it as a future extension placeholder.

## runtime_log

The optional `runtime_log` object controls file-based runtime logging:

```json
"runtime_log": {
  "enabled": true,
  "file": "run.log",
  "level": "info",
  "flush_level": "error",
  "append": false,
  "profile": "default"
}
```

Fields:
- `enabled`: boolean. If omitted, runtime log is enabled by default.
- `file`: output file path. Relative paths are resolved from the config file directory.
- `level`: one of `trace|debug|info|warn|error|critical|off`.
- `flush_level`: flush threshold with the same value set as `level`.
- `append`: append mode when true; overwrite mode when false.
- `profile`: one of `minimal|default|debug`, controls output detail density.

## solver

The solver object is used by DD sweep and Newton solve paths.

Method selection:
- method: `gummel`, `newton`, `poisson_only`, `frozen_state`, or `gummel_newton`
- type: alias for method (legacy compatibility)
- `frozen_state` (alias `diagnostic_state_replay`) is diagnostic-only. It
  requires `sweep.initial_state_file`, performs no nonlinear solve, preserves
  all supplied state fields, and reports zero terminal current by default.
  `sweep.frozen_state_compute_current=true` evaluates terminal currents on the
  supplied state. Replay success and zero residual placeholders do not imply
  that the state satisfies the nonlinear equations or any acceptance gates.

Commonly used controls:
- max_iter
- reltol
- abstol
- temperature_K
- mobility
- recombination
- auger_cn_m6_per_s
- auger_cp_m6_per_s
- auger_excess_product
- auger_with_generation
- auger_density_dependence
- recombination_materials
- impact_ionization

Gummel-specific keys:
- damping_psi
- carrier_floor_m3
- `gummel_carrier_convergence_norm`: `l2` (default) or `componentwise`.
  `componentwise` requires every node's carrier update divided by the maximum
  of old density, new density, and the solve carrier floor to meet `reltol`;
  it prevents a large majority-carrier norm from hiding a local minority-
  carrier update.
- taun
- taup
- bandgap_narrowing

Newton-specific keys:
- damping_factor
- max_update
- quasi_fermi_update_limit_V
- quasi_fermi_update_limit_minority_V
- quasi_fermi_update_limit_mode
- quasi_fermi_trust_region_growth_factor
- quasi_fermi_trust_region_max_multiplier
- quasi_fermi_trust_region_expansion_threshold
- quasi_fermi_trust_region_shrink_factor
- quasi_fermi_trust_region_min_multiplier
- electron_quantum_potential
- carrier_regularization_scale
- line_search
- bank_rose (explicit alternative to backtracking, default disabled)
- warm_start
- contact_boundary_reconstruction
- contact_boundary_minority_electron_relaxation
- contact_boundary_minority_electron_relaxation_bias_threshold_V
- contact_boundary_minority_electron_relaxation_two_terminal_only
- contact_boundary_minority_electron_relaxation_contact_side
- verbose
- diagnostics / diagnostic_history
- jacobian (`analytic` or `finite_difference`)
- finite_difference_step
- quasi_fermi_reference (`none`, `contact_majority`, `contact_basin`, or `contact_nearest_state`)
- carrier_row_convergence
- continuity_row_scaling
- linear_equilibration (`off` or `l2_row_column`)
- global_continuity_closure
- residual_norm (`block` or `l2`)
- residual_weights
- residual_scales
- taun
- taup
- bandgap_narrowing

Hybrid Gummel-Newton keys:
- `handoff.fallback`: `none` or `gummel_on_newton_failure`
- `handoff.require_gummel_convergence`: boolean, default `true`
- `handoff.gummel_max_iter`: optional non-negative integer overriding only the
  Gummel initializer iteration limit
- `handoff.newton_max_iter`: optional non-negative integer overriding only the
  Newton handoff stage iteration limit

Notes:
- `line_search` and `damping_factor` apply to Newton config.
- `bank_rose` selects the published Bank--Rose **Algorithm Global** damping
  recurrence, not a replica of proprietary Sentaurus internals. Use
  `line_search: false`, `line_search_mode: "merit"`, `damping_factor: 1`, and
  `stable_merit_comparison: false` with `bank_rose.enabled: true`; conflicting
  selections are rejected by both JSON parsing and the C++ constructor.
  Options and defaults are `not_damped: 25`, `omega: 1e-4`, `initial_K: 0`,
  `seed_K: 1`, `growth: 10`, `decay: 0.1`, `min_damping: 1e-8`,
  `max_trials: 32`. `not_damped` must be a non-negative integer, `max_trials`
  an integer in `[1,256]`, omega and decay in `(0,1)`, initial_K non-negative,
  seed_K positive, growth greater than one and min_damping in `(0,1]`;
  all real parameters must be finite.
  For one-based iteration `i <= not_damped`, take one full, finite,
  carrier-valid trial without requiring norm descent. Afterwards use
  `alpha = 1/(1 + K*||F||)` and require
  `||F_trial|| <= (1-omega*alpha)*||F||` with strict norm reduction.
  Rejecting a trial changes `K=0` to seed_K, otherwise multiplies K by growth;
  accepting it carries `decay*K` into the next iteration. K has inverse
  residual-norm units: changing norm scaling changes the meaning of initial_K
  and seed_K. K resets at each solve, and rejected searches do not commit it.
  No acceptable trial, nonfinite input, K overflow or the damping floor returns
  failure without installing a rejected state. A zero norm does not authorize
  a step or bypass device qualification. `max_trials` includes the first trial.
  A typed `PopulationDomainError` from trial construction or residual assembly
  (currently the SG-source finite-positive-population check) rejects that trial
  as `bank_rose_population_domain`; active damping increases K and retries within
  the same trial/floor bounds. No density clamp or residual is substituted.
  Before activation, the single invalid full trial still fails. Invalid parent
  evaluation and unrelated exceptions propagate normally; they are not retries.
  The selected block/L2 norm has fixed scales within a solve. Existing direction
  limits, Poisson recorrection, residual equations, and final qualification
  gates remain in effect. These constraints mean the paper's smooth exact-Newton
  convergence theorem is not automatically a guarantee for a device run.
  Applies to coupled Newton (including DC Newton solves) and its Poisson-only
  initializer, not the separate coupled-load-line/arclength correctors or
  Gummel/linear Poisson drivers. With `diagnostics: true`, trial history records
  `bank_rose_active`, `bank_rose_K`, alpha, norm, target and rejection reason;
  accepted-iterate progress JSON also includes the Bank--Rose trial history.
  Source: Bank and Rose, *Global Approximate Newton Methods* (1981), Section 3,
  [author-hosted paper](https://ccom.ucsd.edu/~reb/m270b/newton.pdf).
- `electron_quantum_potential` enables the electron density-gradient quantum
  correction.  The object accepts `coupling_mode` (`outer`, the default, or
  `frozen` for an imported restart potential), `transport_coupling`
  (`direct_band_edge`, the compatibility default corresponding to Sentaurus
  `DirectQuantumCorrection`, or experimental `sentaurus_exponential`), and
  `transport_coupling_weight` (finite in `[0,1]`, default `1`; an optional
  homotopy parameter for the experimental exponential mode). The exponential
  mode keeps the classical Fermi argument and multiplies the transported
  density by `exp(-weight*Lambda/Vt)`; it is not yet qualified for production
  self-consistent DG sweeps. The object also accepts `formulation`
  (`potential_based`, the Sentaurus default, or `density_based` for audit), `gamma`,
  `effective_mass_ratio`, optional `coefficient_mass_ratio`,
  `max_iterations`, `outer_max_iterations`,
  `relative_tolerance`, `absolute_tolerance_V`, optional
  `outer_absolute_tolerance_V`, `damping`, and
  `max_update_V`.  The outer fixed-point loop also accepts
  `outer_acceleration` (`none` or vector `aitken`), `outer_relaxation`,
  `outer_relaxation_min`, and `outer_relaxation_max`; every outer iteration is
  written to the runtime log with its raw/applied potential change.  `max_update_V`
  is the non-negative infinity-norm cap on an
  inner density-gradient Newton update; its default is `0.1 V`, and `0`
  disables the cap.
  `outer_absolute_tolerance_V` controls only the outer fixed-point Lambda
  change. When omitted, `absolute_tolerance_V` is used for both inner and
  outer convergence for backward compatibility.
  The qualified semiconductor-side interface model may set
  `interface_boundary: sentaurus_step` to apply the O-2018.06 Eq. 233
  inhomogeneous step boundary at semiconductor-to-unsolved-nonmetal edges.
  Its material controls are `theta` (default `0.5`),
  `insulator_effective_mass_ratio` (the unsolved-side DOS mass, default
  `0.42`), and `conduction_band_narrowing_fraction` (default `0.5`).  The
  generic default remains `homogeneous_neumann`; imported Sentaurus decks
  with `eQuantumPotential` select `sentaurus_step` explicitly.
  Material-resolved global interface experiments may instead set
  `include_insulators: true` together with
  `insulator_gamma` (default `1`) and `insulator_effective_mass_ratio`
  (default `0.42`). `insulator_coefficient_mass_ratio` is the optional
  coefficient-mass counterpart. When either coefficient-mass control is
  omitted, its DOS mass retains the legacy shared-mass behavior. This global
  mode is experimental.
  `gss_density_fitted` represents a continuous quantum potential Lambda while
  reconstructing the band/DOS potential and sqrt-density traces independently
  on each material side. The qualified default keeps semiconductor-side
  homogeneous Neumann treatment unless `sentaurus_step` is selected.
  For a fixed-reference audit of that experimental global mode,
  `residual_diagnostic_prefix` writes `<prefix>_cells.csv`, `_nodes.csv`,
  `_regions.csv`, and `_summary.txt` with the initial Eq. 231 stiffness,
  gradient-square, reaction, and total residual contributions. Set
  `residual_diagnostic_use_initial_state: true` to evaluate the drives from
  the imported checkpoint rather than the preceding frozen-quantum DD solve.
  Both diagnostic controls default to disabled and do not alter the nonlinear
  equation when omitted.
  `global_discretization` selects the all-material Eq. 231 spatial operator:
  `p1_direct` (default) follows the expanded potential-based equation used by
  the P1/box endpoint audit. `exponential_fitted` exactly preserves the
  `theta=0.5` identity with `u=exp(w/2)` for manufactured-solution checks, but
  is opt-in because steep imported oxide quantum-potential gradients produce
  exponentially ill-conditioned element ratios.
  `cvfem_full` is an experimental median-dual control-volume finite-element
  integration of the complete expanded Laplacian and quadratic-field terms,
  intended for comparison with Charon's density-gradient CVFEM architecture.
  `p1_lambda_direct` is a diagnostic two-factor control: it retains the
  expanded P1 volume operator but uses continuous Lambda with material-side
  band/DOS traces. It isolates the interface unknown from the fitted flux.
  `gss_potentiallike_fitted` is the orthogonal diagnostic: it retains the
  continuous potential-like state and changes only the spatial operator to
  the GSS Eq. 9.128 fitted flux.
  `sentaurus_box` is an experimental interface-trace prototype. It combines
  that fitted edge operator with signed AverageBox circumcentric
  element-vertex measures and
  reconstructs side-specific fitted and reaction traces at a shared vertex.
  It is not a qualified default and must pass the fixed-state
  and endpoint gates before use in a sweep.
  `sentaurus_interface_insulator_half_jump_offset`,
  `sentaurus_interface_silicon_half_jump_offset`, and
  `sentaurus_interface_polysilicon_half_jump_offset` independently shift the
  fitted dimensionless half jump on each region side of a shared material
  vertex. They default to `0`. The SingleDevice O-2018.06 direct
  residual/Jacobian oracle selects about `0.02012` on the insulator side,
  equivalent to a potential-like region trace shift of about `-1.03 mV` at
  300 K. This calibration is local to the experimental `sentaurus_box` mode;
  it is not applied to other discretizations.
  Four side-specific reaction weights replace the invalid single interface
  multiplier: `sentaurus_interface_silicon_reaction_weight`,
  `sentaurus_interface_polysilicon_reaction_weight`,
  `sentaurus_interface_insulator_at_silicon_reaction_weight`, and
  `sentaurus_interface_insulator_at_polysilicon_reaction_weight`. Each scales
  both the material-side reaction and its diagonal Jacobian and defaults to
  `1`. They represent an algebraically eliminated region-side trace. The
  SingleDevice fixed-state oracle calibrates them explicitly; they are not
  global material parameters or qualified defaults.
  Four matching `*_reaction_offset_V` controls add the constant term of each
  affine material-side trace. They default to `0`; the residual uses
  `weight * Lambda + offset` while the consistent diagonal Jacobian contains
  only `weight`. This distinguishes a true eliminated interface trace from a
  simple reaction-volume multiplier.
  `sentaurus_insulator_reentrant_corner_reaction_weight` applies only to a
  pure nontransport vertex with six incident triangles split 2:4 between two
  distinct insulator materials. It scales that corner's reaction and
  consistent diagonal Jacobian, defaults to `1`, and is ignored outside
  `sentaurus_box`. The earlier SingleDevice SiO2/Nitride corner fit selected
  about `0.9713` only as compensation for a positive mixed-area obtuse-cell
  fallback. After adopting the signed AverageBox measure, the qualified
  SingleDevice configuration restores this control to its neutral value `1`;
  ordinary nodes along the same interface are unchanged.
  `conservative_sqrt_fitted` retains the continuous potential-like state but
  assembles the exact theta=0.5 sqrt-density weak form, including the matching
  Lambda-times-sqrt-density reaction. A common fixed row scale prevents
  overflow without changing the nonlinear root.
  `gss_density_fitted` follows the GSS unexpanded sqrt-density control-volume
  flux. It uses the Eq. 9.128 exponential/second-order piecewise flux with an
  `expm1` small-jump evaluation and a globally continuous Lambda unknown; it
  is the experimental explicit-oxide path for Si/SiO2/PolySilicon structures.
  `oxide_boundary: devsim_wkb` adds the experimental DEVSIM/Garcia-Asenov
  oxide interface closure. Its defaults follow the public DEVSIM MOSCAP
  example: oxide quantum mass `oxide_quantum_mass_ratio=0.14`, WKB barrier
  mass `oxide_barrier_mass_ratio=0.4`, and electron barrier
  `oxide_barrier_height_V=3.15`. `oxide_boundary` defaults to `none`; these
  controls are diagnostic and do not change the qualified SingleDevice path.
  WKB is a truncation boundary for a mesh that omits the oxide quantum
  equation. It is therefore rejected when `include_insulators: true`; an
  explicit SiO2 solve must use `oxide_boundary: none`.
  A material JSON entry may override the solver-level parameters with
  `electron_quantum_gamma`, `electron_quantum_dos_mass_ratio`, and
  `electron_quantum_coefficient_mass_ratio`. The DOS mass controls only the
  logarithmic DOS drive; the coefficient mass controls
  `gamma*hbar^2/(6*m*q)`. If the latter is omitted it falls back to the DOS
  mass for backward compatibility. Eq. 231 evaluates its OldSlotboom band
  drive from total impurity only; this does not change the generic DD/BGN
  model used by the transport equations.
- `contact_boundary_reconstruction` controls only Ohmic-contact boundary
  quasi-Fermi/carrier reconstruction in the Newton path. Accepted values are
  `dominant_signed_contact_mean` (default, current behavior) and
  `legacy_node_local` (uses each node's local signed doping without
  contact-mean polarity alignment). This is a boundary reconstruction policy,
  not a mobility calibration knob. Use `legacy_node_local` for a physical
  multi-polarity terminal such as the Templates/LDMOS source metal that shorts
  a p+ body pickup to the n+ source. `dominant_signed_contact_mean` is intended
  for isolated compensated/tie-node sign outliers; applying it to a real p+/n+
  terminal replaces the minority segment with the terminal-average polarity.
  With `quasi_fermi_reference: "contact_basin"`, reference seeds follow the
  selected boundary doping rule per contact node. Thus a mixed source/body
  terminal in `legacy_node_local` mode supplies both electron and hole seeds;
  its mean polarity does not suppress one of the carrier reference basins.
- `contact_boundary_minority_electron_relaxation` controls whether Newton
  relaxes minority-electron quasi-Fermi pinning on p-side Ohmic contacts at
  higher applied bias. Default is `true` (current behavior).
- `contact_boundary_minority_electron_relaxation_bias_threshold_V` sets the
  minimum `|Vbias|` that activates this relaxation. It must be finite and
  non-negative. Default is `0.1`.
- `contact_boundary_minority_electron_relaxation_two_terminal_only` keeps this
  relaxation restricted to two-terminal decks when `true` (default). Set to
  `false` to allow the same relaxation policy in multi-terminal decks.
- `contact_boundary_minority_electron_relaxation_contact_side` selects which
  contact polarity relaxes minority-carrier pinning once relaxation is active.
  Accepted values are `p_contact_only` (default), `n_contact_only`, and
  `both_contacts`.
- `contact_boundary_minority_electron_relaxation_strength` controls how much of
  the selected minority pinning is released when relaxation is active. `1.0`
  matches the current full-relaxation behavior, while `0.0` leaves the
  minority contact pinned.
- `max_update` is an optional non-negative infinity-norm cap on one Newton
  update in solver unknown units; `0` disables the cap.
- `quasi_fermi_update_limit_V` is an optional non-negative physical-voltage
  cap applied only to Newton quasi-Fermi updates (`phin` and `phip`) after
  `max_update` and before line search; `0` disables the cap. In
  `unit_scaling` mode the value is converted to solver unknown units using the
  potential scale, so the accepted `phin`/`phip` delta remains capped in volts.
- `quasi_fermi_update_limit_minority_V` is an optional non-negative
  physical-voltage cap applied only to the minority-carrier quasi-Fermi update
  at each node, classified by local net doping (`netDoping < 0` p-type makes
  `phin` the minority update; `netDoping > 0` n-type makes `phip` the minority
  update). When greater than zero it tightens only the minority carrier to
  `min(quasi_fermi_update_limit_V, quasi_fermi_update_limit_minority_V)` while
  the majority carrier keeps the looser `quasi_fermi_update_limit_V`. `0`
  disables the minority-specific cap and reproduces the global uniform behavior
  (default). It must be finite and non-negative.
- `quasi_fermi_update_limit_mode` selects how those physical quasi-Fermi caps
  constrain a Newton step. `componentwise_poisson_recorrect` (default) clips
  each `phin`/`phip` component and then re-solves the Poisson block.
  `uniform_trust_region` computes the most restrictive active global/minority
  cap and uniformly scales the complete `(psi, phin, phip)` Newton vector. The
  latter preserves the linearized Newton direction and is intended for strongly
  nonlinear avalanche branch correction; it does not relax convergence or
  forcing tolerances.
- `poisson_only_quasi_fermi_initialization` controls the carrier state held
  fixed by `solver.method: "poisson_only"`. `contact_basin` (default) rebuilds
  the interior electron and hole quasi-Fermi fields from the nearest
  majority-carrier contact biases, matching the ABA initialization path.
  `preserve_initial` retains the supplied restart quasi-Fermi fields and only
  projects explicit contact boundary values. Use the latter when correcting
  the Poisson discretization of an imported biased TCAD state; it also prevents
  pure-insulator gauge rows from inheriting a nearby terminal bias.
  Consistent saved reference/increment pairs retain their small free-node
  coordinates, including increments below the absolute-QF rounding precision.
  If an absolute QF was explicitly edited and is inconsistent with its saved
  pair under the ordinary Newton restart check, the absolute value takes
  precedence. `contact_basin` continues to reset the QF fields intentionally.
- The explicit BVMethods CMD sidecar profile defaults both QF cap fields to
  `0`, overriding any implicit cap inherited from warm-start emission. Explicit
  sidecar cap values remain available for comparisons. This profile rule does
  not alter the separate load-line inner-voltage update control or the
  residual-based line search. See the [operator contract](bvmethods_operator_contract.md).
- In the coupled external-resistor corrector, an active QF cap uniformly
  scales the full device/inner-voltage direction before residual line search.
  Setting both caps to zero removes that scaling but retains raw-update
  diagnostics, line search and the circuit's separate voltage-update control.
- For `uniform_trust_region`, `quasi_fermi_trust_region_growth_factor` and
  `quasi_fermi_trust_region_max_multiplier` optionally expand the physical
  update radius after a boundary-limited step. Expansion occurs only when the
  line search accepts the complete step and the ratio of actual residual
  decrease to linearized predicted decrease is at least
  `quasi_fermi_trust_region_expansion_threshold`. Defaults are `1`, `1`, and
  `0.75`, so adaptive expansion is disabled unless explicitly requested.
  `quasi_fermi_trust_region_shrink_factor` contracts the multiplier after an
  accepted step required line-search backtracking, down to
  the positive `quasi_fermi_trust_region_min_multiplier`. Their defaults are
  both `1`, so contraction is likewise opt-in.
- `line_search_mode: block_filter` accepts a trial when at least one enabled
  carrier block achieves its configured sufficient decrease while every block
  remains inside `residual_filter_envelope_factor`. For a limited direction,
  `residual_filter_gamma` multiplies the block decrease predicted by
  `F + J*dx`, rather than incorrectly assuming that the limited step predicts
  a zero residual. When authoritative `block_absolute_convergence` is enforced,
  the sufficient-decrease test applies to blocks that currently exceed their
  absolute ceiling, and an already-qualified block uses that ceiling as the
  lower bound of its envelope. This lets a Poisson boundary update proceed
  without demanding artificial improvement from a carrier block already below
  its hard gate. This changes globalization only; the final Newton convergence
  tolerances are unchanged.
- `line_search_mode: block_closure_filter` extends that filter with electron
  and hole global-continuity blocks. A trial may improve either the enabled
  carrier/device blocks or a global-closure block only while every device and
  closure block remains inside its configured residual envelope. The closure
  envelope uses the stricter applicable relative or low-source absolute target;
  final convergence still requires `global_continuity_closure` itself to pass.
- `local_update_diagnostics` is an opt-in per-node Newton transformation trace.
  In addition to raw, capped, applied, and selected node updates, it records
  aggregate current/raw-predicted/capped-predicted/selected block norms, the
  number and largest physical magnitude of quasi-Fermi cap changes, and the
  signed free electron/hole residual sums. This allows a rejected lift to be
  classified without making those diagnostic quantities part of acceptance.
- `progress_directory` (default empty) opts the classical Newton solver into
  diagnostic HDF5/JSON snapshots at the initial state and each accepted update.
  Use a new absolute directory for every solve. Existing iteration files are
  never overwritten. JSON retains full double precision, block residuals, row
  violations, closure ratios, raw/capped direction norms, actual packed update
  maxima/nodes, cap counts, rounded-away update counts and elapsed wall time.
  `potential_scale_V` converts packed updates to volts. Snapshot I/O cost is
  recorded as `newton.progress_snapshot`. These files always carry
  `qualified: false`; a complete normal
  Newton result and its original gates are still required. Ignore an incomplete
  final file after a hard process kill; earlier unique snapshots remain usable.
- `capture_linear_system` (default false) adds in-memory linear diagnostics to
  a nonempty `progress_directory`, for classical non-split solves only. Each
  assembly attempt gets a unique `linear_N` directory with iteration identity,
  physical/regularized/closure-replaced matrices and packed state, residual,
  reference volts, linear/raw/processed/selected directions. It does not
  reassemble the physical Jacobian or change acceptance. `assembled.csv` uses
  its direction column for the closure-system RHS; other files use it for the
  named direction (selected: actual accepted delta). `processed` precedes any
  subsequent line-search restoration; `selected` is written only on acceptance.
  Matrices use `VELANJ01`, three little-endian uint64 dimensions (rows, columns,
  nnz), int32 CSC outer/inner indices and binary64 values. CSV uses 17 digits.
  Partial files after interruption are not complete evidence. Capture adds I/O
  cost and changes no qualification gate; output paths cannot be reused.
- `soft_time_limit_seconds` (default zero/off) stops classical Newton with
  `converged=false`, reason `soft_time_limit`, at an iteration boundary after
  this wall-clock budget. It does not interrupt an assembly or line search and
  therefore does not replace the caller's hard child-process timeout. Both new
  options reject quantum or recursive carrier-row-recovery solves to avoid
    ambiguous nested budgets/output ownership. They change no convergence gate.
  - `diagnostic_stop_policy` defaults to `off`. Opt-in `row_trend_v1` is a
    failure-only stop monitor for classical packed Newton with enforced carrier
    rows, ordinary residual tolerances and a nonempty `progress_directory`.
    It does not modify residuals, Jacobians, line search or success gates.
    Per carrier it records N/E/max, where E is the sum of positive
    log10(row_ratio/eps_row); it also records the actual primary residual,
    closure (including absolute branch) and configured contact-QF gate ratios.
    Beginning at update 12, compare median samples k-5..k-3 with k-2..k.
    Three consecutive windows with <1% improvement in every row metric and
    every still-failing global gate stop only if the raw pre-cap direction
    shrinks <10%, or at least two of the last three alphas are <=1/128.
    Three consecutive alphas <=1/1024 with <1% improvement at each transition
    can stop earlier. Reasons are `stagnation_detected`, `tiny_step_stagnation`,
    or `nonfinite_trend_metric`. Existing success checks take priority.
    With this policy enabled, reaching `max_iter` returns `iteration_limit`
    after those checks, without post-budget restoration. Full state and metrics
    are saved before stopping. This is an empirical diagnostic policy, not a
    proof of nonexistence or a relaxed acceptance criterion.
- Diagnostic `newton_block_step_probe` additionally accepts
  `block_step_damping` in `(0,1]` (default 1) and mode `carrier_local`.
  That mode requires `local_carrier_nodes`, an array of zero-based
  `{node_id, carrier: "electron" | "hole"}` records, 1 to 64 distinct active
  noncontact rows. It solves the restricted QF Jacobian, applies existing caps
  then damping, and reconstructs densities with the configured assembler,
  including Fermi statistics. Complementary unknowns are frozen; neighboring
  residuals are not. The output is an unqualified diagnostic seed, requiring
  all-row/global audits and full Newton validation. The older `carrier_density`
  mode retains its Boltzmann restriction and its own backtracking.
- `return_best_feasible_on_failure` (default `false`) keeps the failure status
  and diagnostics unchanged but exposes, as the rejected restart state, the
  lowest-residual accepted iterate that satisfied both enforced local carrier
  rows and global continuity closure. This is intended for deterministic
  restart/corrector audits and does not turn a failed solve into convergence.
- `poisson_line_search_stall_contact_majority_qf_drop_limit_V` is an optional
  non-negative physical-voltage guard for the `poisson_line_search_stall_floor`
  acceptance path. When greater than zero, a Poisson-block line-search stall is
  accepted only if the maximum majority-carrier quasi-Fermi drop across contact-to-interior edges
  is below this value. The default is `5e-11` V; `0` disables this guard.
  The species is selected from the net doping at each edge's contact endpoint:
  donors select electrons, acceptors select holes, and zero net doping checks
  both. A mixed n+/p+ electrode is not classified by its mean doping. Edges
  between contacts and species with no positive finite transport density are
  excluded. Referenced QF differences retain their reference/increment split.
- `block_absolute_convergence` is an optional authoritative per-block Newton
  convergence contract. With `mode: enforce`, scalar `reltol`, `abstol`, and
  numerical-floor acceptance no longer bypass the block contract: the raw
  `psi`, electron-continuity, and hole-continuity L2 norms must be at or below
  `psi_residual_ceiling`, `electron_residual_ceiling`, and
  `hole_residual_ceiling`. All three ceilings must be positive. The default
  mode is `off`, preserving existing solver behavior.
- With `line_search_mode: block_filter` and both block and carrier-row
  convergence enforced, an iterate that already meets every block ceiling
  may still correct unbalanced local carrier rows. Such a trial must remain
  within every original block ceiling and decrease the maximum qualified
  carrier-row residual ratio by the configured filter sufficient-decrease
  factor. Carrier-row `off` and `report` modes retain the usual block filter.
- `quasi_fermi_recenter_on_stall` (default `false`) enables bounded in-solve
  quasi-Fermi coordinate repartition in classical coupled Newton with enforced
  absolute block convergence. After the Poisson block qualifies, an unresolved
  carrier row whose nonzero line-search correction rounds away can trigger a
  change of reference. Reference-plus-increment low bits are retained; residuals
  and coordinate-dependent caches are refreshed, and the direction is recomputed.
  At most eight retries are allowed per solve. A retry is not a Newton update,
  does not reset the iteration budget or residual normalization, and does not
  bypass any original stopping gate. Profiling counts `newton.qf_recenters` and
  `newton.qf_recenter_discarded_line_search_trials`. Enabling local-update
  diagnostics also records carrier/residual/Jacobian changes at repartition.
- `local_update_diagnostics` optionally writes selected carrier rows and linear
  solve diagnostics (`enabled`, `csv_file`, `nodes`, `first_iterations`,
  `every_iterations`). The CSV includes the uncapped raw system's infinity-norm
  relative residual, normwise backward error
  `||J dx + r||inf / (||J||inf ||dx||inf + ||r||inf)`, and maximum componentwise
  backward error `max_i |J dx+r|_i / (|J| |dx|+|r|)_i`. These are measured before
  continuity row weighting and equilibration, after mapping the linear solution
  back to the original coordinates. Profiling records this diagnostic's cost
  as `newton.local_update_diagnostics`; it does not change Newton acceptance.
  The appended coordinate-writeback columns distinguish the selected high
  coordinate change from the full split high/low change and from the nominal
  damped step. Rejected split candidates whose low parts were discarded have
  blank full-change/error fields and `complete_coordinate_change_available=0`.
  These selected-candidate records must not be interpreted as accepted updates
  unless `line_search_accepted=1`. BVDS has no concentration finite-step
  transformation in this path; that stage belongs to `slot_population_memory`.
- For `simulation_type: newton_jvp_probe`, top-level
  `freeze_transport_mobility: true` holds each edge's mobility at the supplied
  base state for the positive/negative perturbations. This diagnostic requires
  `solver.mobility.jacobian_field_derivatives: false`, edge SG transport, and
  no coupled avalanche. The default remains a finite difference of the live
  residual. This option permits a matched audit of a deliberately lagged
  mobility Jacobian without changing the production residual or its root.
- `contact_majority_qf_branch_drop_limit_V` is an optional non-negative guard
  applied to every formal Newton convergence path and to best-iterate
  eligibility. `contact_majority_qf_branch_guard_contacts` may restrict the
  maximum contact-to-interior majority-carrier QF drop to named transport
  contacts; an empty list checks all non-metal-gate contacts. Unknown contact
  names are rejected. The default limit is `0`, so this guard is opt-in.
- A failed Newton attempt retains the lowest-residual accepted iterate that
  passes the configured QF branch guard. When rejected-state capture is
  enabled, `newton_attempts.csv` records its iteration, residual, and QF drop,
  and writes an additional `*_best.csv` full-precision restart beside the
  parent, initial, and final snapshots.
- `carrier_regularization_scale` is an experimental non-negative Newton
  stabilization knob. When greater than zero, Vela adds
  `sign(diagonal) * scale * carrier_row_abs_sum` to each carrier continuity
  diagonal of each free carrier row before solving the coupled Newton step;
  carrier Dirichlet identity rows are excluded, so this option can be used with
  `exact_dirichlet_updates`. The carrier-only diagnostic uses the same boundary
  exclusion. The row sum includes the
  full coupled carrier row across `psi`, `phin`, and `phip` columns. `0`
  disables the regularization. This is intended for BV branch-stability
  diagnostics and is not enabled by default.
- `carrier_floor_m3` is an optional non-negative Gummel carrier floor used to
  keep reconstructed quasi-Fermi potentials consistent with solved carrier
  densities.
- In `dc_sweep`, when `solver.method` is `newton`, `poisson_only`, or
  `gummel_newton` and
  `solver.diagnostics: true`, the CSV appends opt-in recombination diagnostics:
  `recombination_max_abs_rate_m3_per_s`,
  `recombination_mean_abs_rate_m3_per_s`, and
  `carrier_product_max_np_over_ni2`.
  These columns are disabled by default.
- `taun` and `taup` override the SRH lifetime defaults (`1e-5 s` and `3e-6 s`).
- `auger_cn_m6_per_s` and `auger_cp_m6_per_s` override the Auger coefficients
  passed into the recombination model (`2.90e-43` and `1.028e-43 m^6/s`
  Sentaurus 2018 silicon-at-300 K defaults).
  Negative values are rejected by the recombination model validation.
  With `scaling.mode: "unit_scaling"` the compiled defaults are expressed in
  the internal unit system (`2.90e-31` and `1.028e-31 cm^6/s`), and an explicit
  deck value is read and kept internally as `cm^6/s`.
- `auger_excess_product` selects the Auger excess-carrier product under Fermi
  statistics. `generalized_fermi` is the compatibility default and reuses the
  generalized SRH product. `classical_np` uses `n*p-ni_eff^2`; its analytic
  Jacobian differentiates the transport carrier densities while treating the
  node material/BGN `ni_eff` as fixed during one Newton assembly.
- `auger_with_generation` is a boolean (shared isothermal compatibility default
  `true`). Set `false` to suppress nonpositive Auger excess products and their
  derivatives, including the inactive-side derivative at equilibrium. This
  matches plain Sentaurus `Auger`; `true` corresponds to `Auger(WithGeneration)`.
  SRH generation remains enabled. Newton, Gummel and diagnostics share this
  switch. The independent audited thermal interface defaults to `false`.
- `auger_density_dependence` optionally enables
  `Cn(n)=Cn0*(1+Hn*exp(-n/N0n))` and the corresponding hole coefficient.
  Omission preserves constant coefficients. The object accepts `enabled`
  (default true when present), plus `electron` and `hole` objects containing
  `enhancement` (nonnegative H, default 0) and `reference_density_m3`
  (positive N0, default 1e24 m^-3). As with other concentration input fields,
  explicit N0 is in cm^-3 under `unit_scaling`, so native N0=1e18 cm^-3
  is entered as `1e18` in that mode. Newton, Gummel recovery and diagnostics
  use the same rate; the analytic Jacobian includes the coefficient's density
  derivative. Cn0/Cp0 still come from the explicit Auger coefficient fields;
  this option does not imply a temperature-dependent A/B/C polynomial.
- `recombination_materials` optionally restricts nodal SRH/Auger to a list of
  exact material names, for example `["Si"]`. Omission or an empty list preserves
  all-material behavior. Newton, Gummel, source Jacobian audits and isothermal
  exports use the same mask. Source volumes are independently selected by
  `recombination_volume_policy`. Unknown/duplicate
  names and nodes shared by included and excluded semiconductor materials are
  rejected. Insulator interfaces are allowed. This option does not restrict
  BTBT or avalanche and does not define an electrothermal material policy.
- `recombination_volume_policy`: `legacy_node` (default) preserves the existing
  nodal SRH/Auger measure. `material_local_mixed_voronoi` sums mixed-Voronoi
  vertex shares from semiconductor Tri3 cells only, excluding dielectric area
  at shared nodes. `material_local_circumcentric` instead sums signed
  circumcentric shares and requires positive finite assembled transport-node
  measures. Residuals, source derivatives, Gummel linearization and
  continuity diagnostics use the same measure. Transport, Poisson, BTBT and
  avalanche retain their own measures. This is an explicit isothermal control,
  not a change to template or C++ defaults.
- Both Gummel/Newton parse `mobility`, `recombination`, `impact_ionization`, `temperature_K`.
- With `scaling.mode: "unit_scaling"`, `bandgap_narrowing.reference_doping_m3`
  is read and kept internally as `cm^-3`.

`gummel_newton` runs the configured Gummel solve first, validates that solution,
then runs coupled Newton with `warm_start=true` from the Gummel state. The
default fallback policy is strict: a Newton failure fails the sweep point. Use
`gummel_on_newton_failure` only for diagnostic curves where a finite Gummel
result is preferable to aborting the sweep.

Reference-import configs may also use `vela_step` and `vela_stop` on an
individual simulation entry to override the generated Vela sweep range while
preserving the full imported reference curve. `vela_current_contact` may be set
when the Vela terminal current to compare differs from the swept bias contact.
A simulation `comparison` block can pass curve gate options such as
`candidate_scale`, `bias_min`, `bias_max`, `reference_column`,
`candidate_column`, `interpolation`, `max_orders_of_magnitude`,
`max_relative_error`, `min_points`, and `require_trend_match` to the comparison
report. `interpolation` defaults to `linear`; use `log_current` for exponential
current-magnitude comparisons when adjacent current samples keep the same
nonzero sign after scaling.
`runtime_diagnostic` is an
optional simulation block:

```json
"runtime_diagnostic": {
  "enabled": true,
  "doping_scale": 0.0001,
  "step": 0.1
}
```

When enabled, it creates an additional conservative runtime deck using
region-average scaled doping. When disabled or omitted, the faithful deck is
the executable comparison path.

Reference import configs may include:

```json
"tdr_doping": {
  "compensated_node_policy": "reported"
}
```

Supported policies:
- `reported`: preserve `doping.csv` exactly as merged from region-local TDR fields and report compensated nodes in `doping_metadata.json`.
- `dominant_signed_region`: when a global node receives equal donor and acceptor active concentrations, use the dominant signed `DopingConcentration` field to choose a single majority dopant for the node, and record the rewrite in `doping_metadata.json`. If equal positive and negative signed values tie, Vela preserves the existing first-region tie-break but records `signed_doping_tie: true` with `resolution_source: "signed_aggregate_tie_first_region"`.

For compensated nodes, `doping_metadata.json` records original donor/acceptor
values, resolved donor/acceptor values, whether a rewrite occurred, and
`resolution_source` (`reported`, `signed_aggregate_doping`,
`signed_aggregate_tie_first_region`, `neighbour_region_sign`, or
`unresolved_tie`).

### carrier_statistics

The coupled Newton drift-diffusion solver accepts either form below:

```json
"carrier_statistics": "fermi_dirac"
```

```json
"carrier_statistics": { "model": "fermi_dirac" }
```

Supported values are `boltzmann` (the compatibility default),
`fermi_dirac`, and the explicitly restricted `fermi_dirac_joyce_dixon`.
The last option requires `solver.coupled_operator: "bvmethods_tri3_qf"`;
see [the operator contract](bvmethods_operator_contract.md) for its supported
models, coordinates and diagnostics. `solver.coupled_operator` otherwise
defaults to `legacy`. The defining-integral Fermi-Dirac selection is one atomic physics path: carrier
densities use the normalized complete integral `F_{1/2}`, ideal Ohmic contacts
solve charge neutrality with the same statistics, continuity edges use the
generalized Einstein/Scharfetter-Gummel operator, and terminal-current
post-processing uses that identical edge operator. The density-form Gummel
solver currently rejects `fermi_dirac` explicitly; use `method: "newton"`.

An independent, opt-in transport discretization is available as
`"carrier_statistics": {"model":"fermi_dirac", "transport":"effective_band_sg"}`.
The default `transport: "generalized_einstein"` preserves the existing operator.
Both choices retain the same defining-integral Fermi populations, contacts and
BGN; this option changes the edge transport operator, not merely its Jacobian.
For either carrier let `d = (QF1-QF0)/Vt` for electrons and
`d = -(QF1-QF0)/Vt` for holes, `a = log(c1/c0)+d`, and
`K = mobility*Vt*edge_couple/edge_length` including the active unit conversion.
The effective-band particle continuity flux is `K*c1*B(a)*expm1(d)`.
It uses thermal SG weights without the generalized Einstein multiplier, cancels
flat QF exactly, and changes sign on edge reversal. Residuals, transport probes,
Jacobian evaluation and terminal currents use the selected operator together.
The Jacobian retains the full population and mobility dependence; existing
Fermi endpoint partials use centered finite differences, with analytic mobility
chains and the selected source Jacobian. It is not a wholly symbolic derivative.
This mode currently supports classical coupled Newton only: Boltzmann statistics,
quantum correction and `element_qf_gradient` transport are rejected. No material
constants, avalanche source profile or acceptance tolerances are changed by this
option. Previously qualified states require requalification under a new operator.

The Fermi integral now uses a high-accuracy piecewise Chebyshev evaluation,
with analytic derivatives and matched nondegenerate/degenerate expansions.
The inverse uses relative density error, including dilute carrier ratios.
The apparent BGN correction calls this same inverse; there is no separate
lower-accuracy statistics approximation inside BGN. This replaces the former
Bednarczyk approximation, so archived Fermi states must be requalified under
the new executable. Historical frozen executables remain the replay baseline.
Coefficient generation is a development-only operation using
`scripts/generate_fermi_half_coefficients.py` (NumPy/SciPy); normal builds use
the checked-in table and require neither package.

### bandgap_narrowing

`solver.bandgap_narrowing` accepts either a string or an object.

String form:

```json
"bandgap_narrowing": "slotboom"
```

Object form:

```json
"bandgap_narrowing": {
  "model": "slotboom",
  "reference_doping_m3": 1.0e23,
  "coefficient_eV": 0.009,
  "smoothing": 0.5,
  "offset_eV": 0.0,
  "fermi_statistics_correction": false
}
```

Supported `model` values:
- `none`
- `slotboom`
- `old_slotboom`
- `bennett_wilson`

`bennett_wilson` implements Sentaurus Bennett-Wilson (User Guide, Eq. 165):
`DeltaEg = coefficient_eV * ln(Ntotal/Nref)^2` for `Ntotal > Nref`, otherwise
zero. `Ntotal` is donors **plus** acceptors, not net doping or free-carrier
density. Defaults are `coefficient_eV=0.00684` and `Nref=3.162e18 cm^-3`
(`3.162e24 m^-3`). Both string and object forms work in Newton and Gummel.
In an object, `reference_doping_m3` follows the existing unit convention:
use `3.162e18` with `unit_scaling`, `3.162e24` with `legacy_si`.
`smoothing` and `offset_eV` must be zero for this model; a material's native
`Bandgap.dEg0(Bennett)` belongs in its base intrinsic density, not in the
positive narrowing term. `fermi_statistics_correction` remains false by default.
The SLOT IDL mapping uses Boltzmann statistics and symmetric band-edge splitting
(`Bgn2Chi=0.5`); arbitrary material-dependent splitting is not added here.

The text importer recognizes explicit `BennettWilson`. It does not silently
infer a BGN choice from an omitted Physics option; resolve native defaults from
the run log and material parameters before a comparison. Conflicting scoped BGN
choices and Bennett-Wilson/Fermi imports fail closed. General `.par` reference
temperature/material emission is a separate contract; this addition does not
claim arbitrary parameter-library import coverage.

The low-level Slotboom and OldSlotboom prototypes compute the positive
bandgap-narrowing term from the maximum of the supplied impurity concentration
and supplied carrier densities. Coupled Newton builds its effective intrinsic
density from **total impurity** (donors plus acceptors), passing zero carrier
arguments. Thus its BGN `ni_eff` varies between nodes but is held fixed during
one state/Jacobian evaluation. For Sentaurus `OldSlotboom` parity, the
`models.par` `Bandgap.dEg0(OldSlotboom) = -1.595e-2 eV` term is handled by the
material intrinsic-density override, while the `old_slotboom` BGN term uses
`Ebgn = 9e-3 eV`, `Nref = 1e17 cm^-3`, and `C = 0.5`. This is implemented in
Gummel and Newton configurations. With
`scaling.mode: "unit_scaling"`, `reference_doping_m3` numeric input is read as
`cm^-3` and kept internally as `cm^-3`. The compiled `Nref` default of both
`slotboom` and `old_slotboom` is likewise expressed in the internal unit
system (`1e23 m^-3` in `legacy_si`, `1e17 cm^-3` in `unit_scaling`).

For state-independent-mobility Boltzmann transport, the coupled Jacobian's BGN potential
partials differentiate the same stable, variable-ni SG flux used by the
residual. This avoids subtracting large drift and diffusion derivatives near
flat quasi-Fermi equilibrium. It applies with either value of the separate
`equal_ni_flux_evaluation` control, which selects the no-BGN equal-ni residual
evaluation (`legacy_factor_difference` or `compensated_log_expm1`). The stable
partial requires both population exponents inside the existing ±500 clamps
and a positive effective Bernoulli factor; other branches retain their
existing treatment. Mobility chain derivatives and Fermi-statistics paths
remain separate. Field/surface models with deliberately disabled mobility
chain derivatives retain their existing approximate-Jacobian protocol until
separately qualified. This change does not alter ni, the flux, or convergence gates.

Set `fermi_statistics_correction: true` to reproduce the additional bandgap-
narrowing correction that Sentaurus applies by default when `Fermi` and
`EffectiveIntrinsicDensity(OldSlotboom)` are active together.  The correction
is evaluated from donor/acceptor concentrations and the 300 K density of states.
Leave it `false` to reproduce Sentaurus `EffectiveIntrinsicDensity(NoFermi)` or
the historical Vela behavior.

The SDevice frontend emits this flag as `true` for OldSlotboom with `Fermi`
and as `false` for explicit `EffectiveIntrinsicDensity(NoFermi)`. Without an
explicit Fermi/NoFermi selection, it omits the flag and retains the C++ default `false`. `NoFermi` does not disable Fermi carrier
statistics. Direct JSON configuration retains its existing opt-in default.
Results produced by the older frontend with Fermi + OldSlotboom require a
new physical baseline; parser/deck golden updates alone do not qualify them.

### Newton diagnostics and residual options

Newton configs can opt into diagnostic history with either
`"diagnostics": true` or `"diagnostic_history": true`. The solver also accepts:

- `jacobian`: `analytic` or `finite_difference`
- `finite_difference_step`
- `quasi_fermi_reference`: `none` (default), `contact_majority`, `contact_basin`, or `contact_nearest_state`
- `residual_norm`: `block` or `l2`
- `residual_weights`: object with `psi`, `phin`, and `phip`
- `residual_scales`: object with `psi`, `phin`, and `phip`

`quasi_fermi_reference: "contact_majority"` keeps the external solution and
boundary-condition values as absolute volts, but stores the internal electron
quasi-Fermi unknown relative to the most strongly n-type biased contact and the
hole unknown relative to the most strongly p-type biased contact. Because each
carrier uses one constant reference, edge differences and the physical model
are unchanged. This avoids losing sub-femtovolt current-carrying increments
when an absolute quasi-Fermi potential is tens of volts from zero.

`contact_basin` partitions the transport graph by distance to majority-carrier
contact seeds. For mixed source/body terminals, seed classification follows
`contact_boundary_reconstruction`.

`contact_nearest_state` instead chooses, independently for each carrier and
node, the non-metal terminal voltage closest to the supplied quasi-Fermi state.
This is an explicit numerical-coordinate option; boundary conditions and
physical models do not change. References stay fixed during each solve. A
subsequent warm restart or diagnostic may choose a new partition from its input
state. Without an input state, initialization falls back to geometric contact
basins. Saved reference/increment pairs preserve sub-ULP changes when rebasing.

For the explicit `split_dd_state` path, Newton warm starts transform saved QF
increments into the new contact-reference frame using multiprecision arithmetic
before splitting the result into binary64 high/low components. The mesh identity
and potential scale must match. Unchanged references retain their original pairs;
changed references incur the rounding of the new pair, not a lossless guarantee.
Ordinary diagnostic checkpoint restoration still requires identical references.
The DC linear predictor extrapolates the complete split state and refreshes its
projected potentials, increments and 300 K Boltzmann initial-guess densities.
New contact constraints are applied by the solver; rejected-step retries use the
previous accepted state without extrapolation. This does not change convergence
gates or enable split arithmetic by default.
With explicit `bias_points`, accepted intermediate steps also advance the
predictor's state and bias history together; selecting fewer output points does
not change which accepted states define the next extrapolation.

The source-aware convergence controls are:

The runner's `carrier_row_convergence` result separates selected continuity rows
into diagnostic `scale_groups.floor_limited` and `non_floor_limited` groups.
The former means a positive `scale_floor` is at least as large as
`max(abs(fluxAbsSum), abs(recombination), abs(impact))`, including equality.
Each group reports selected row count, violation count, maximum ratio and
maximum absolute residual in internal equation units. Ignored/inactive rows
are not included. Violations also expose `physical_scale` (before the floor)
and `floor_limited`. `floor_absolute_limit = eps_row * scale_floor` states the
effective absolute residual limit for floor-limited rows. These fields are
diagnostic only: no row is exempted, no threshold or `satisfied`/`enforced`
semantics change, and a floor-limited ratio is not a port-current error estimate.

An independent explicit policy `carrier_row_convergence.floor_limited_action`
can be `enforce` (default, historical behavior) or `report`. With `report`,
floor-limited violations remain in the output but have `blocks_acceptance:false`;
`blocking_violation_count` excludes them and `satisfied` means no blocking
violations remain. Non-floor-limited violations still block when mode is
`enforce`. Local carrier-row recovery receives blocking rows only. This policy
does not alter residuals, scales, Newton directions, global closure, boundary or
numerical-validity gates. It is an acceptance-policy change, not a claim that
weak-row physical error is bounded. The SLOT validation batch opts in explicitly;
no general/default template is changed.
The carrier-row ratio is an additional local balance target, not an
established Sentaurus-equivalent termination criterion. In the BVmethods
no-BTBT control, all of the first five native accepted nonzero states exceed the
`0.001` row target under the same source eligibility profile. See the
[native row-gate audit](validation/bvmethods_native_row_gate_2026-10-21.md).
The audit's first revision omitted the TCAD current-density/length conversion;
its source eligibility counts have been corrected. See the
[full convergence-rule audit](validation/bvmethods_native_convergence_2026-10-21.md).
Keep reference-state diagnostics distinct from solver acceptance; this finding
does not change the production defaults or authorize a global threshold increase.

```json
"carrier_row_convergence": {
  "mode": "enforce",
  "eps_row": 0.001,
  "scale_floor": 1e-30,
  "min_source_scale_fraction": 0.0,
  "min_source_scale": 1e-18
},
"continuity_row_scaling": {
  "enabled": true,
  "flux_fraction": 0.0,
  "scale_floor": 1e-30,
  "min_source_scale": 1e-18,
  "min_weight": 1e-12,
  "max_weight": 1e18
},
"linear_equilibration": {
  "mode": "l2_row_column"
},
"global_continuity_closure": {
  "mode": "enforce",
  "tolerance": 0.01,
  "source_floor": 1e-14,
  "absolute_tolerance": 1e-12,
  "linear_system_stabilization": "off",
  "feasibility_restoration": "off",
  "feasibility_restoration_factor": 1.0
}
```

`continuity_row_scaling` left-scales the Newton continuity rows using the
larger of their configured flux/source measures. It changes conditioning only,
not the nonlinear equations. `min_source_scale` supplies an absolute scale for
low-current problems where a fraction of the largest source would suppress the
small rows of interest.

`linear_equilibration: {"mode": "l2_row_column"}` applies one pass of L2 row
and column scaling to the complete linearized Newton system after any
continuity-row scaling. The recovered physical Newton increment solves the
same equations; this option only improves finite-precision sparse
factorization on matrices with large row/column norm spreads. The default is
`off`. Boolean `true` is accepted as shorthand for `l2_row_column`.

`global_continuity_closure` independently sums the free-node SRH/avalanche
source and all contact fluxes for electrons and holes. `report` appends the
metrics without changing the ordinary merit line-search trajectory or
rejecting a solution; `enforce` preserves that same trajectory and prevents
final convergence unless each carrier with
`max(|contact flux|, |integrated source|) >= source_floor` has relative
mismatch no greater than `tolerance`. Qualifying on both sides prevents a
large contact-current/small-source state from bypassing the gate.
The denominator is the largest of the contact-flux magnitude,
integrated-source magnitude, and `source_floor`. Below that floor the relative
ratio remains diagnostic and the absolute mismatch must instead be no greater
than `absolute_tolerance` (default `1e-12`). The IV CSV records each carrier's
qualification flag, selected `absolute`/`relative` criterion, signed mismatch,
and relative ratio.

Global closure shapes intermediate line-search acceptance only when
`line_search_mode: "block_closure_filter"` is selected explicitly. Ordinary
`merit` and `block_filter` modes use the device residual alone; this keeps the
conservation identity as a final solution invariant instead of imposing it as
an implicit continuation-path constraint.

In `block_closure_filter`, the electron and hole closure normalizers are
frozen from the bias-point restart state, with `source_floor` as their lower
bound. They are not recomputed from intermediate iterates, so closure ratios
remain comparable while terminal current is rebalanced.

`linear_system_stabilization: "replace_free_row_sum"` is an opt-in
finite-precision stabilization. For each carrier it replaces one free-node
continuity row by the sum of all free-node continuity rows. This is an
invertible row transformation, so the exact Newton equations are unchanged,
but the solve must directly resolve the global conservation mode. The
implementation applies the rank-one/rank-two row replacements through a small
Woodbury border and reuses the original sparse Jacobian factorization; it does
not factor a matrix containing dense global rows. It is mutually exclusive
with `continuity_row_scaling`; the default is `"off"`.

For isolated SLOT comparisons, `scripts/slot_newton_stabilization.py` exposes
independent `--carrier-regularization on|off` and
`--closure-replacement on|off` controls. It writes a new configuration and
refuses to overwrite an existing file. Regularization `on` uses `1e-8` unless
`--regularization-scale` is supplied; `off` writes zero. Closure replacement
`off` disables only the linear row transformation, **not** global continuity
acceptance. Physics, boundary settings, and acceptance gates are preserved.
For example, from the worktree root:

```powershell
python -B scripts/slot_newton_stabilization.py input.json comparison.json --carrier-regularization off --closure-replacement off
```

`feasibility_restoration: "linear_qf_projection"` is a separate opt-in
globalization aid for componentwise-limited Newton steps. After the normal
update caps and Poisson recorrection, it forms the two sums of free electron
and hole Jacobian rows and applies the minimum-norm correction over only the
unconstrained electron/hole quasi-Fermi unknowns. The projected direction
satisfies both global continuity rows to first order without moving contact
Dirichlet unknowns. The normal QF caps are then reapplied without a second
Poisson recorrection, so the bounded direction is deliberately an inexact
feasibility restoration; nonlinear trial acceptance remains responsible for deciding
whether it is useful. Final local/global convergence gates remain unchanged.
Rank-deficient or non-finite projections are skipped; the default is `"off"`.
`feasibility_restoration_factor` selects a fraction in `(0, 1]` of the
projected correction before QF recapping; it defaults to `1`.

`global_continuity_closure.audit_jsonl_file` optionally appends JSON Lines
records for the selected replacement rows (initial state and every accepted
Newton iteration), final-restoration triggers, and every nonlinear restoration
trial, including rejected trials. Records distinguish physical contact/source
closure from assembled free-row residual sums. They include normalized-row SVD
conditioning, projection multipliers, predicted and reassembled sums, raw/
projected/capped/nominal/applied step norms, cap locations, carrier qualification,
and rejection reasons (`projection_failed`, `nonfinite_residual`, `device_floor`,
`carrier_rows`, `global_closure`, or `accepted`). Nonfinite JSON numbers are null.
The ideal alpha lower bound freezes the initial closure mode and assumes an exact
linear correction; `ideal_linear_infeasible` labels trials without skipping them
or changing the existing halving schedule. The final corrector uses QF caps only;
its audit explicitly reports that Poisson recorrection was not applied.

`global_continuity_closure.restoration_jacobian` defaults to `"stabilized"`.
The opt-in `"physical"` mode excludes both the assembler carrier diagonal floor
and solver carrier-row regularization from the **final** closure sensitivity.
The main Newton solve retains both configured stabilizers. With an audit file,
physical-Jacobian predictions are also recorded for stabilized restoration trials.
This is an A/B diagnostic option, not a change to default convergence or physics.

`global_continuity_closure.restoration_physical_rows_only` defaults to `false`.
When enabled, the **final** QF projection excludes inactive carrier gauge rows
from its constraint sums and excludes their unknowns from its correction.
Contact Dirichlet rows/unknowns remain excluded. Gauge equations enforce
auxiliary coordinates in non-transport regions; their residuals must not cancel
physical contact/source imbalance. The main Newton equations and stabilization
are unchanged. This can be combined with `quasi_fermi_reference: "contact_basin"`
to retain small corrections that absolute biased QF coordinates cannot represent.

`global_continuity_closure.roundoff_audit` defaults to `false` and requires
`audit_jsonl_file` to emit results. It adds expensive, read-only diagnostics at
the initial state and the first final-restoration trial if finite: fixed-residual
reductions (ordinary, pairwise, Kahan, Neumaier, long double, magnitude-ordered,
permuted and Boost quad reference), actual forward/reverse/permuted edge
assembly with ordinary or compensated per-node accumulation, and bitwise
repeatability. Records distinguish physical term sums from inactive gauge sums.
Audit variants never replace the residual used by Newton or its acceptance gates.

`feasibility_restoration: "final_linear_qf_projection"` uses the same
minimum-norm QF construction only after the accepted device residual is at or
below `stall_residual_floor`, carrier-row convergence passes, and the enforced
final closure gate fails. It recomputes the Jacobian at that accepted state,
applies the configured QF bounds, and backtracks the correction. A trial is
accepted only if its nonlinear device residual remains below the stall floor
and both carrier-row and global-closure gates pass. This mode therefore cannot
turn an unconverged device iterate into an accepted solution and does not
constrain intermediate Newton steps.

`warm_start: true` preserves supplied quasi-Fermi potentials when continuing
from a previous solution. The default `false` keeps the conservative
cold-start behavior.

### band_to_band

NonLocalPath is not a supported public model yet. The experimental C++
contracts, geometry, phonon kernel and private SI serializer are separate from `solver.band_to_band`;
they do not enable a nonlocal solver or import native parameter files. See the
[N1 interface record](validation/bvmethods_nonlocalpath_n1_2026-10-16.md) and
[N2 geometry record](validation/bvmethods_nonlocalpath_n2_2026-10-16.md), plus the
[N3 private kernel record](validation/bvmethods_nonlocalpath_n3_2026-10-16.md).
The kernel requires explicit reservoir windows. Private paired-source assembly
and frozen actual-residual diagnostics now exist, but native source/window
fidelity and derivatives remain unqualified; see the
[N4 investigation](validation/bvmethods_nonlocalpath_n4_2026-10-16.md).

`solver.band_to_band` enables a local electron-hole pair-generation source in
both Newton and Gummel continuity equations.  The `e2` model implements the
Sentaurus Device O-2018.06 simple silicon expression
`G = A |F|^2 exp(-B/|F|)`. Vela evaluates `F` in each semiconductor triangle,
integrates the local rate over that cell, and mass-lumps the paired source to
its vertices; insulator cells are excluded even at shared nodes. The default
parameters are `A=3.4e23 m^-1 s^-1 V^-2` and `B=2.26e9 V/m`, equivalent
to the Sentaurus parameter-file units `3.4e21 cm^-1 s^-1 V^-2` and
`22.6e6 V/cm`.

```json
"band_to_band": {
  "model": "e2",
  "A_cm_inv_s_inv_V_inv2": 3.4e21,
  "B_V_per_cm": 22.6e6,
  "minimum_field_V_per_cm": 0.0,
  "jacobian": "frozen_field"
}
```

The canonical SI keys are `A_m_inv_s_inv_V_inv2`, `B_V_per_m`, and
`minimum_field_V_per_m`.  The `cm` aliases above are provided for direct
transcription of Sentaurus parameters.  `jacobian: "frozen_field"` (default)
re-evaluates `G` from the current potential at every nonlinear residual but
omits `dG/dpsi` from that linear solve, matching the solver's controlled
high-field quasi-Newton strategy.  `potential_finite_difference` includes the
full derivative of the semiconductor-cell field/source integral and is intended
for focused small-mesh verification because its cost grows with the number of
potential unknowns.  `model: "none"` disables the source.

### impact_ionization

`solver.impact_ionization` accepts either a legacy model string or an object.

String form:

```json
"impact_ionization": "none"
```

Object form:

```json
"impact_ionization": {
  "model": "selberherr",
  "coupling_mode": "self_consistent",
  "parameter_set": "default",
  "driving_force": "electric_field",
  "eparallel_field_recovery": "edge_adjacent_cells",
  "quasi_fermi_gradient_discretization": "edge_difference",
  "generation": "carrier_density",
  "debug_raw_vanoverstraeten": false,
  "A_scale": 1.0,
  "electron_A_m_inv": 7.03e7,
  "electron_B_V_m": 1.231e8,
  "hole_A_m_inv": 1.582e8,
  "hole_B_V_m": 2.036e8,
  "carrier_velocity_m_s": 1.0e5
}
```

Supported `model` values:
- `none`
- `selberherr`
- `van_overstraeten`

Supported `coupling_mode` values:
- `self_consistent` (default): include avalanche generation in the continuity
  residual and Jacobian.
- `postprocess_only`: evaluate avalanche coefficients, driving fields,
  current/source diagnostics, and VTK output on the solved state, but exclude
  avalanche generation from the continuity residual and Jacobian. This is an
  observation mode; it does not change the default.

Field meanings (Selberherr prototype):
- `electron_A_m_inv` (1/m): electron ionization prefactor.
- `electron_B_V_m` (V/m): electron critical field.
- `hole_A_m_inv` (1/m): hole ionization prefactor.
- `hole_B_V_m` (V/m): hole critical field.
- `carrier_velocity_m_s` (m/s): effective saturated carrier speed used by the
  generation-rate proxy.
- `driving_force`: `electric_field` (default), `quasi_fermi_gradient`,
  `grad_potential_parallel_j`, `effective_field_parallel_j`, or `eparallel`.
  `eparallel` implements the Sentaurus carrier-specific vector projection
  `max(E dot Jn/|Jn|,0)` and `max(E dot Jp/|Jp|,0)` using a shared
  electrostatic-field geometry and cell-reconstructed SG current vectors.
  Sentaurus
  `Avalanche(VanOverstraeten)` decks use `quasi_fermi_gradient`; the
  current-aligned options are Charon-style SG edge-current diagnostics and
  require `generation: "current_density"`.
- `eparallel_field_recovery`: electric-field vector recovery used only with
  `driving_force: "eparallel"`. `edge_adjacent_cells` (default) preserves the
  existing area-weighted average of the cells sharing each edge.
  `nodal_vertex_star` first recovers a field at each endpoint from the complete
  transport-material triangle star, then averages the two endpoint vectors;
  this matches the Sentaurus vertex-field interpolation observed in the
  BVmethods NMOS reference without changing the compatibility default.
- `generation`: `carrier_density` (legacy `alpha*v*n/p` proxy) or
  `current_density` (`alpha_n*mu_n*n*|grad(phin)| + alpha_p*mu_p*p*|grad(phip)|`).
- `current_approximation` (with `generation: "current_density"`):
  `density_gradient`/`grad_qf` use the per-carrier SG continuity flux,
  `cell_reconstructed` uses `mu * n_mid * |driving field|` with a
  Bernoulli-weighted edge-midpoint density, `cell_current_reconstructed` and
  `cell_vector_current_reconstructed` use cell-smoothed SG flux magnitudes, and
  `conserved_total_current` feeds both carriers the conserved total-current
  magnitude `|F_p - F_n|` (divergence-free in the converged state, so the
  avalanche source does not collapse on the depleted side of a reverse-biased
  junction).
- `edge_source_partition`: SG edge-current source split. The default
  `symmetric` assigns each carrier source 50/50 to the edge endpoints. Use
  `qf_gradient` only for `grad_qf` or explicit experimental probes that need
  quasi-Fermi-gradient directional endpoint weights.
- `quasi_fermi_gradient_discretization`: `edge_difference` (default) preserves
  the existing Vela GradQf behavior. `cell_gradient` is Genius-compatible for
  `II.Force=GradQf`: Vela rebuilds electron/hole quasi-Fermi potentials from
  `psi`, carrier density, `ni`, and `Vt`, then uses the area-weighted adjacent
  cell-gradient magnitude as the ionization-coefficient driving field. This
  switch only affects `solver.impact_ionization`; it does not change
  `solver.mobility.high_field_driving_force`.
- `quasi_fermi_carrier_truncation`: default-off, non-negative low-density
  feedback-support diagnostic. When greater than zero, the quasi-Fermi
  potentials used only by the avalanche driving field are rebuilt with
  `n_eff=max(n,value*ni)` and `p_eff=max(p,value*ni)`. The physical carrier
  state, continuity transport current, contacts, and mobility driving field are
  unchanged. It is accepted with `triangle_gss_gradqf_truncated` so a
  controlled AvalDens-like support experiment can change this axis without
  changing source geometry or current reconstruction. This is not declared to
  be a Sentaurus-equivalent default.
- `source_geometry_scale`: positive finite diagnostic multiplier for the
  Scharfetter-Gummel edge-current avalanche source geometry. The default is
  `1.0`; values other than `1.0` are intended only for BV parity probes.
- `source_volume_policy`: SG edge-current avalanche source support preset.
  `edge_half_box` keeps the default `0.5 * h * edge.couple` ownership;
  `edge_box` uses `1.0 * h * edge.couple` for focused source-ownership probes.
- `source_volume_factor`: default-off diagnostic override for SG edge-current
  avalanche source support. `0` uses `source_volume_policy`; finite values in
  `[0.5, 1.0]` directly set the source area factor in
  `factor * h * edge.couple`.
- `minimum_field_V_m`: non-negative Charon-style cutoff for avalanche
  coefficients. `0` disables the cutoff; Charon van Overstraeten defaults are
  commonly probed with `5.0e6 V/m` (`5.0e4 V/cm`).
  - `debug_raw_vanoverstraeten`: default-off diagnostic switch for
    `model: "van_overstraeten"`. When `true`, Vela forces the avalanche
    coefficient drive to `|grad(eQuasiFermi)|` and `|grad(hQuasiFermi)|` and
  bypasses `minimum_field_V_m`, `driving_force_interpolation` RefDens blending,
  quasi-Fermi carrier truncation, contact-element electric-field fallback, and
    current-aligned driving-force projection. This is intended only to isolate
    cutoff/smoothing/RefDens suppression and should not be used as a production
    calibration knob.
  - `parameter_set`: explicit Van Overstraeten/de Man diagnostic parameter
    override. Accepted values are `default`, `sentaurus_fit_A_only`,
    `sentaurus_fit_A_B`, and `sentaurus_fit_A_B_switch`. The default leaves the
    built-in silicon parameters unchanged. The Sentaurus-fit options use the
    node-output effective fit `van_overstraeten_sentaurus_effective_fit_v1` only
    when requested; they do not replace the material-library defaults. The
    `sentaurus_fit_A_B_switch` option also sets the Van Overstraeten
    `switch_field_V_m` equivalent to `2.5e5 V/cm`.
  - `A_scale`: positive finite Van Overstraeten/de Man diagnostic multiplier
    for the four `A` prefactors only (`electron_a_low_m_inv`,
    `electron_a_high_m_inv`, `hole_a_low_m_inv`, `hole_a_high_m_inv`). The
    default `1.0` leaves the default path unchanged. It does not modify `B`,
    `switch_field_V_m`, low/high branch selection, cutoff, smoothing, or
    RefDens-style driving-force blending.
  - `B_scale`: positive finite Van Overstraeten/de Man diagnostic multiplier
    for the four `B` critical-field values only (`electron_b_low_V_m`,
    `electron_b_high_V_m`, `hole_b_low_V_m`, `hole_b_high_V_m`). The default
    `1.0` leaves the default path unchanged.

Validation:
- `electron_A_m_inv`, `hole_A_m_inv`, and `carrier_velocity_m_s` must be non-negative.
- `electron_B_V_m` and `hole_B_V_m` must be positive.
- `driving_force` must be `electric_field`, `quasi_fermi_gradient`,
  `grad_potential_parallel_j`, `effective_field_parallel_j`, or `eparallel`.
- Current-aligned driving forces require `generation: "current_density"` and
  `current_approximation: "density_gradient"` or `"grad_qf"`.
- `quasi_fermi_gradient_discretization` must be `edge_difference` or
  `cell_gradient`; `cell_gradient` is valid only with
  `driving_force: "quasi_fermi_gradient"`.
- `generation` must be `carrier_density` or `current_density`.
- `source_volume_policy` must be `edge_half_box` or `edge_box`.
- `source_volume_factor` must be `0` or finite within `[0.5, 1.0]`.
- `quasi_fermi_carrier_truncation` must be non-negative and finite.
  - `minimum_field_V_m` must be non-negative and finite.
  - `debug_raw_vanoverstraeten` requires `model: "van_overstraeten"`.
  - Non-default `parameter_set` values require `model: "van_overstraeten"`.
  - `A_scale` must be positive and finite; values other than `1.0` require
    `model: "van_overstraeten"`.
  - `B_scale` must be positive and finite; values other than `1.0` require
    `model: "van_overstraeten"`.

Scaling:
- With `scaling.mode: "unit_scaling"`, `electron_A_m_inv` and
  `hole_A_m_inv` numeric inputs are read as `cm^-1`, while
  `electron_B_V_m` and `hole_B_V_m` are read as `V/cm`. They are normalized
  as internal `cm^-1` and `V/cm` values before the impact-ionization model sees them.
- The Van Overstraeten model also accepts split low/high-field parameters:
  `electron_a_low_m_inv`, `electron_a_high_m_inv`, `electron_b_low_V_m`,
  `electron_b_high_V_m`, `hole_a_low_m_inv`, `hole_a_high_m_inv`,
  `hole_b_low_V_m`, `hole_b_high_V_m`, `switch_field_V_m`,
  `phonon_energy_eV`, `reference_temperature_K`, and `temperature_K`.
  Its defaults match the Sentaurus 2018 silicon `vanOverstraetendeMan`
  parameters at 300 K.
- `carrier_velocity_m_s` remains `m/s`.

Prototype note:
- The default `carrier_density` path preserves the original engineering source
  term for smoke diagnostics. For Genius-style BV decks, use
  `driving_force: "quasi_fermi_gradient"`, `generation: "current_density"`,
  `current_approximation: "density_gradient"` or `"grad_qf"`,
  `quasi_fermi_gradient_discretization: "cell_gradient"`, and usually
  `quasi_fermi_carrier_truncation: 1.0e-2`.

### mobility

The experimental object model `"ialmob"` requires
`carrier_current_discretization: "scharfetter_gummel_edge"` and
`high_field_driving_force: "quasi_fermi_gradient"`. Its additional `ialmob`
object contains:

- `geometry_file`: an absolute path to `vela.ialmob.transport_geometry.v1` JSON.
  The header records `node_count` and `cell_count`; every Si cell supplies
  `cell_id`, mesh-order `node_ids`, `vertex_measure_m2` in m2 and nonnegative
  dimensionless `edge_coefficients` in `(v0,v1),(v1,v2),(v2,v0)` order.
  The summed coefficients must reproduce the active transport coupling.
- `effective_electrodes`: explicit mesh contact names whose adjacent cells use
  electric-field fallback. All contact edges are excluded from ordinary
  exterior boundary projection, including contacts absent from this list.
- `avalanche_geometry_file`: optional absolute path to the controlled
  `vela.ialmob.observed_avalanche_geometry.v1` source geometry. The header
  supplies integer `node_count`, `cell_count` and `cells`; optional `scope` is
  descriptive provenance. Every Si triangle appears exactly once, with
  `cell_id`, mesh-order `node_ids`, three `coordinates_m` pairs,
  `mobility_vertex` (local index 0--2), and three `source_measure_m2` shares.
  Other materials are rejected. Coordinates must match the actual mesh;
  shares must be finite, nonnegative, have positive sum, and cannot exceed
  the triangle area (relative export-roundoff allowance `1e-10`). This is
  independent of the transport geometry and must be regenerated/revalidated
  after changing the mesh. It is measured geometry, not a carrier-state table.
- `crystal_x`, `crystal_y`: three-component crystal directions corresponding to
  the two mesh coordinate axes.
- `electron_parameters_cm`, `hole_parameters_cm`: objects keyed by orientation
  family (`"100"`, `"110"`, `"111"`). Fields use the native centimeter-based
  IALMob parameter convention in both scaling modes. Every family encountered
  by the mesh must have parameters; unknown parameter fields are rejected.
- `high_field`: default `true`; `reference_density_m3`: SI m^-3, default `1e18`,
  in both scaling modes. Outer electron/hole field parameters retain the normal
  active-unit interpretation.
- `constant_material_mobility_cm2_per_Vs`: optional mapping of explicit non-Si
  material names to positive finite `electron` and `hole` mobilities. For
  example, `{"PolySilicon":{"electron":1417,"hole":470.5}}` restores constant
  poly transport alongside Si IALMob. Values always use cm2/(V s), in both
  scaling modes. Names must exist on the mesh; an edge shared with Si IALMob
  is rejected. Constant mobility has zero field/density derivatives, while
  SG flux derivatives remain active. Material `ni`, doping and contacts must
  also enable transport; this mapping alone does not create carriers.
  The isolated SLOT poly extension is qualified here only at 300 K without
  avalanche. Its material energy reference, thermal behavior and avalanche
  coupling require separate validation. Use `recombination_materials: ["Si"]`
  for the original SLOT deck's Silicon-only SRH/Auger scope.

The live coupled path supports 300 K, Si/SiO2 Tri3 geometry, Fermi or Boltzmann
carrier statistics, without quantum potential. Avalanche additionally requires
Boltzmann statistics, a supported element source and
`source_jacobian: local_ad`. Its cell SG currents use the live element IALMob;
all nine potential columns include density, normal-field and HFS feedback.
The independent source audit recomputes mobility at each perturbed state;
multiprecision source audits are rejected because the IALMob kernel uses double.
This implementation does not establish SLOT high-voltage or breakdown qualification.
The additional controlled source is selected explicitly with
`impact_ionization.current_approximation: "ialmob_element_density_grad_qf"`.
It requires the observed avalanche geometry above, `van_overstraeten`,
`generation: "current_density"`, `driving_force: "quasi_fermi_gradient"`,
`quasi_fermi_gradient_discretization: "cell_gradient"`,
`source_mapping_mode: "element_vertex_box_measure"`, and `source_jacobian: "local_ad"`.
Use fallback scope `contact_node_cell`, fallback mode `cell_gradient_magnitude`,
and `ialmob_boundary_field_projection: true`. The solver-level
`grad_qf_contact_field` selects the contact field.
The IALMob reference density and avalanche interpolation/truncation must be zero;
complete mobility derivatives are required. The identified HFS parameters are
electron/hole saturation velocities `1.07e5`/`8.37e4 m/s` and beta `1.109`/`1.213`,
at 300 K. Source-volume scaling overrides, directional partition, debug raw
alpha and alternative current-magnitude policies are unsupported.
The source is restricted to Si and uses extrema-based harmonic carrier density,
the observed vertex mobility, full contact QF current gradient, and tangential
non-contact interface gradients. Contact alpha uses the electric field by default. All
nine local potential derivatives include carrier density and live IALMob/HFS.
The production mobility retains the mathematical screening minimum; inferred
native numerical clipping constants are diagnostic evidence only. This source
does not establish a universal native near-angle selector algorithm or qualify
an original-model BVDS curve. See the
[controlled source validation](validation/slot_ldmos_observed_source_integration_2026-09-23.md).
For a separate boundary-field control, `impact_ionization.ialmob_boundary_field_projection`
defaults to `false`. When enabled, the canonical element-edge SG/GSS-Laux
current-density source with quasi-Fermi driving force projects the avalanche
quasi-Fermi gradient onto the insulating-boundary tangent of prepared element
IALMob geometry. Bulk cells retain the full gradient; configured contact electric
fallback retains precedence. Scalar source and local AD use the same projection,
including its potential derivatives. This option changes the ionization field
only: SG source currents and IALMob transport remain unchanged. It does not
implement Sentaurus `AvalDensGradQF` density averaging or establish native parity.
Material interfaces
take precedence over ordinary exterior faces at interface terminations;
multiple nonparallel interfaces of equal priority are rejected. Distance,
orientation, unnormalized distance gradient, vertex HFS and Measure averaging
are evaluated on the live state. `jacobian_field_derivatives` controls the full
cell mobility feedback; qualification uses `true`. Residuals, SG terminal
currents and density recovery share this transport interface. External geometry
weights are evidence dependencies, not regenerated reference solution values.
This is an explicit experimental profile; consult
[current LDMOS validation](validation/templates_ldmos_current_status.md) before
claiming D4 curve qualification.

The isotropic SLOT profile supplies the original unnamed parameter block to
all three orientation families and explicitly sets `reference_density_m3: 0`.
Its original deck does not activate low-density HFS interpolation. The D4
default `1e18 m^-3` must not be transferred to SLOT. The standalone
`ialmob_element_probe` accepts `reference_density_m3` and `high_field` for the
same controlled comparison.

`solver.mobility` accepts either the legacy string form or an object. String decks remain compatible:

```json
"mobility": "caughey_thomas"
```

Object form:

```json
"mobility": {
  "model": "caughey_thomas_field_surface",
  "carrier_current_discretization": "scharfetter_gummel_edge",
  "high_field_driving_force": "electric_field",
  "high_field_gradient_discretization": "edge_projection",
  "electron_mu_min_m2_V_s": 0.00522,
  "electron_nref_m3": 9.68e22,
  "electron_alpha": 0.68,
  "hole_mu_min_m2_V_s": 0.00449,
  "hole_nref_m3": 2.23e23,
  "hole_alpha": 0.70,
  "electron_saturation_velocity_m_s": 1.07e5,
  "electron_field_beta": 1.109,
  "hole_saturation_velocity_m_s": 8.37e4,
  "hole_field_beta": 1.213,
  "surface": {
    "theta_electron_m_per_V": 2.0e-8,
    "theta_hole_m_per_V": 1.0e-8,
    "beta": 1.0,
    "reference_field_V_per_m": 0.0,
    "min_factor": 0.05,
    "max_factor": 1.0,
    "surface_region": "p_body",
    "surface_interface": ["p_body", "gate_oxide"]
  }
}
```

`carrier_current_discretization` selects the continuity-equation current
support. `scharfetter_gummel_edge` is the default and preserves existing Vela
behavior. `element_qf_gradient` is an explicit experimental Tri3 profile: it
reconstructs the P1 quasi-Fermi gradient in each semiconductor cell, evaluates
the cell current from area-averaged density and mobility, and assembles the
conservative P1 divergence. It is currently restricted to bulk mobility with
self-consistent impact ionization disabled. It must not be used as a silent
replacement for the SG profile; Templates/LDMOS fixed-state qualification on
2026-08-29 did not show an accuracy improvement.

`jacobian_field_derivatives` defaults to `true`. Setting it to `false`
holds the base-state mobility fixed only while forming the transport
Jacobian; the nonlinear residual still evaluates the live mobility model.
For bulk edge-projected quasi-Fermi HFS away from contact-field fallback,
the Jacobian applies the analytic mobility chain rule to the frozen-mobility
SG flux derivative. The cell-vector path includes the mobility feedback from
both edge endpoints and all vertices of adjacent transport cells. Contact
electric-field fallback replaces that QF feedback with its electric-potential
stencil. `constant_field` participates in the
field-derivative path even when contact fallback is disabled.

Supported `model` values are `constant`, `constant_field`, `caughey_thomas`,
`caughey_thomas_field`, `caughey_thomas_surface`,
`caughey_thomas_field_surface`, `masetti`, `masetti_field`,
`masetti_surface`, `masetti_field_surface`, and the explicit experimental
`ialmob` object described above.
For bulk/prototype field-saturation models, `high_field_driving_force` is `electric_field`
by default and may be set to `quasi_fermi_gradient` to match Sentaurus
`HighFieldSaturation`; electrons use `|grad(phin)|` and holes use
`|grad(phip)|`.
The solver-level `grad_qf_contact_field` applies to quasi-Fermi-gradient
driving for both high-field mobility and impact ionization. All JSON decks
default to `"electric_field"`. With `"quasi_fermi_gradient"` the contact
stencil uses the quasi-Fermi field. The Sentaurus command frontend emits the
override for `Math { ComputeGradQuasiFermiAtContacts = UseQuasiFermi }`.
This setting takes precedence over the legacy per-model
`contact_electric_field_fallback` booleans in JSON. It selects the driving
field without replacing the carrier current used to evaluate ionization.

The BVMethods Tri3 and SLOT IALMob element kernels classify a cell as touching
an electrode when at least one vertex belongs to an applicable contact.
BVMethods uses mesh contacts and SLOT uses the configured effective electrodes.
Each cell selects its own field, including cells that share an edge. Legacy
edge-assembled profiles retain their edge field reconstruction and weighting.

When quasi-Fermi driving is selected, `high_field_gradient_discretization`
controls how that magnitude is recovered. `edge_projection` (default) uses the
potential difference projected onto each edge and preserves legacy Vela
behavior. `transport_cell_vector` reconstructs the two-dimensional P1 gradient
inside adjacent semiconductor cells, area-averages the gradient vectors, then
takes the magnitude on each edge; this matches the orientation-independent `GradQuasiFermi` semantics used
by Sentaurus on unstructured 2-D meshes.

For referenced quasi-Fermi states, the cell-vector mobility path forms each
cell gradient from reference differences and explicit increments before
combining them. Residuals, analytic mobility feedback, edge diagnostics and
terminal currents share this arithmetic so a large common potential offset
does not erase sub-ULP variations. Restart packing subtracts the new reference
before adding the saved increment. Checkpoint frame translations must also
retain reference-rounding remainders in the increment; the Python
`scripts/translate_dd_state.py` utility provides this operation with explicit
transport node IDs and refuses to overwrite an existing checkpoint. A translated
state still requires the original convergence and state/current equivalence
checks before a sweep may accept it.

For SG current-density avalanche integration,
`impact_ionization.source_volume_policy=genius_conservative` normalizes the
three truncated edge-box pieces of every semiconductor triangle to its exact
area and excludes adjacent insulator cells. This is the conservative
counterpart to legacy `genius_truncated` and is intended for comparisons with
Sentaurus `IntegrSemiconductor AvalancheGeneration`.
For frozen-state diagnostics,
`impact_ionization.source_mapping_mode=nodal_eparallel_p1` reconstructs the
electric field and conventional electron/hole current vectors at transport
nodes, evaluates the two carrier ionization rates there, and integrates the
nodal generation with the exact semiconductor P1 measure (one third of each
incident triangle area). It is restricted to `postprocess_only`, `eparallel`,
`current_density`, `nodal_vector_current_reconstructed`, and
`eparallel_field_recovery=nodal_vertex_star`. Consequently the electric field,
carrier current, ionization coefficient, and P1 generation integral are all
evaluated at the same transport node. It is not a self-consistent
continuity-equation source mapping.
The `masetti` models implement a Masetti-style silicon doping-dependent
mobility shape with configurable electron/hole fields:
`*_mu_const_m2_V_s`, `*_mumin1_m2_V_s`, `*_mumin2_m2_V_s`,
`*_mu1_m2_V_s`, `*_pc_m3`, `*_cr_m3`, `*_cs_m3`,
`*_masetti_alpha`, and `*_masetti_beta`, where `*` is `electron` or
`hole`. The `surface` block is a MOS prototype for Si/SiO2-like channel
mobility degradation, not a calibrated Lombardi model. It applies a
vertical-field factor `mu_eff = mu_bulk / (1 + (theta * max(|E_normal| -
reference_field, 0))^beta)^(1/beta)`, optionally clamped by
`min_factor`/`max_factor`.

The first implementation estimates `E_normal` with the local edge electric-field magnitude on edges that match `surface_region` and/or the two-name `surface_interface`; this is sufficient for trend regressions but should not be interpreted as a calibrated normal-field extraction. If no matching surface edge is found for a mobility evaluation, surface degradation is disabled and the existing low-field or velocity-saturation behavior is used.

For the BVmethods full-physics comparison, `masetti_lombardi` and
`masetti_field_lombardi` select the Sentaurus O-2018.06 default Enhanced
Lombardi `Enormal` coefficients for silicon. The implementation evaluates the
acoustic-phonon and surface-roughness contributions from Eqs. 278--281 of the
Sentaurus Device User Guide, applies the `exp(-distance/l_crit)` interface
damping, and combines them with the bulk mobility using Matthiessen's rule.
When surface degradation and high-field saturation are both enabled, the
surface-limited **low-field** mobility is the input to Canali saturation.
This ordering applies to electrons and holes in the constant, Caughey--Thomas
and Masetti combinations; the PhuMob combinations already use this order.
Live field derivatives use the same composition as the residual. The
2026-10-10 correction changes the constitutive result, so earlier numerical
trajectory qualifications must be rerun with the corrected core.
`surface.surface_interface` should identify the semiconductor/insulator pair;
the normal field is the electric-field projection along that interface normal.

When a Lombardi carrier-dependent roughness exponent is enabled (`alpha != 0`
with nonzero roughness weight), edge transport evaluates mobility separately
at the two endpoint carrier densities and averages the two mobilities. The
edge doping and surface-field conventions stay the same. Residuals, terminal
currents, and the coupled Jacobian use the same density fields. With
`jacobian_field_derivatives=true`, local central differences include both
carrier quasi-Fermi potentials, the endpoint electrostatic potentials, and
all adjacent-cell field reconstruction dependencies. In particular,
`transport_cell_vector` freezes only its GradQF drive in the endpoint part;
it does not freeze the surface field or carrier densities. These are numerical
derivatives of the implemented flux, not symbolic derivatives of every model.
The default Lombardi `alpha=0` and bulk Masetti behavior are unchanged.

`slot_ialmob` and `ialmob_field` select the legacy SLOT isotropic IALMob model
(T-2022.03 equations 290--311). The latter also applies high-field saturation
to the complete IALMob low-field mobility. These models require separate
nodal donors and acceptors and `doping_concentration_basis: "total_impurity"`
(their default). Both carrier densities affect screening and scattering, even
when `alpha_sr=0`; the full surface/density Jacobian dependency path applies.
Use `surface.surface_region` and `surface.surface_interface` to identify Si
and its insulating interface.

The low-field SLOT name was `ialmob` before the September 2026 main integration;
use `slot_ialmob` to reproduce that endpoint model. The name `ialmob` now refers
to the element transport profile documented above and requires its own geometry
and orientation parameter block. Existing SLOT `ialmob_field` inputs retain their
original meaning. Neither model's validation transfers to the other.

The optional `mobility.ialmob` object contains `electron` and `hole` objects
whose coefficient names match the unnamed Sentaurus parameter block:
`mumax`, `mumin`, `theta`, `n_ref`, `alpha`, `ndop_ref`, `nsc_ref`, `S`, `p`,
`l_crit_c`, `B`, `C`, `lambda`, `k`, `alpha_ph2d_A/D`, `lambda_ph2d_A/D`,
`delta`, `lambda_sr`, `A`, `alpha_sr`, `nu`, `eta`, `N1`, `N2`,
`alpha_sr_A/D`, `lambda_sr_A/D`, `l_crit`, `EnormMinimum`, and the
`D1`, `D2`, `nu0`, `nu1`, `nu2`, `alpha1`, `alpha2` families with `_inv`
and `_acc` suffixes. Here `A/D` denotes two separate keys.
Shared keys are `nref_D/A`, `cref_D/A`, `me_over_m0`, `mh_over_m0`,
`full_phumob` (default true), `clustering_everywhere` (false), and
`phonon_combination` (0, 1, or 2; default 1). Unknown keys are rejected.

All coefficients inside `ialmob` always retain the source model's **cm-based
units**, including concentrations in `cm^-3`, lengths in `cm`, mobility in
`cm^2/(V s)`, and the dimensional fractional-power coefficients. They do not
change under `scaling` or `format_version: 2`. Runtime state quantities are
converted separately. Defaults are listed in
`include/vela/physics/IALMobilityParameters.h`. The implementation excludes
ThinLayer, stress, AutoOrientation, and low-temperature phonon corrections;
temperatures below 50 K are rejected. Current high-field parameter alignment
is validated at 300 K. Model/spatial comparison status is recorded in
[the IALMob validation report](validation/slot_ldmos_ialmob_alignment_2026-09-07.md).

The element-box mobility supports `solver.mobility.constant_materials`, an
optional list of material names for `element_box_phumob`. It returns each
selected material's own electron/hole
mobility before density, Lombardi and HFS corrections. Unlisted materials keep
the selected model. The SLOT IdVd package selects only `PolySilicon`; recombination
material scope is configured independently with `recombination_materials`.

The explicit surface configuration uses
`model: "phumob_lombardi"`, `edge_averaging: "element_box_phumob"`,
`doping_concentration_basis: "total_impurity"`, and
`surface.discretization: "element_distance_gradient"`. It uses true vertex
distances to the default semiconductor/insulator interface and the **nonunit**
gradient of their linear triangle interpolant. The cell field is shared by the
cell's vertices; each vertex retains its own distance, impurities and live
carrier densities before box averaging. The option requires element-box
transport and cell-material Poisson geometry. Custom interface selectors,
nonzero Lombardi `alpha`, and high-field saturation are outside this candidate.
Its numerical/self-consistent qualification is tracked in the SimpleMOS branch
status; selecting it does not assert that those checks passed.

`surface.acoustic_factor` and `surface.roughness_factor` multiply the respective
inverse-mobility terms for both carriers. They default to 1 and must be finite
and nonnegative. Setting both to zero recovers the bulk model. To reproduce a
common factor `f` in native Lombardi without stress, divide B, C, delta and eta
by `f`. Native `a_ac`/`a_sr` are stress-enhancement parameters and are inactive
without the corresponding stress model; they are not these direct multipliers.
Existing surface models retain
`surface.discretization: "legacy_cell_centroid"` by default. The explicit split
state runtime supports the candidate's 300 K, alpha-zero, Boltzmann/SRH scope;
potential and population chain derivatives include every adjacent cell vertex.

An additional explicit HFS candidate uses `model: "phumob_field_lombardi"`
with the same element-distance/box settings, plus
`high_field_driving_force: "quasi_fermi_gradient"`,
`high_field_gradient_discretization: "element_vertex_partial_layer"`, and
`jacobian_field_derivatives: true`. It currently requires 300 K and alpha zero.
Each vertex receives its Canali correction **before** box integration. The
drive is the linear cell quasi-Fermi gradient projected parallel to true
semiconductor boundary edges (PartialLayer); two nonparallel boundaries give
zero projection. Any contact vertex instead selects the full electric gradient.
Raw drives below 1 V/cm disable saturation on both branches. This cutoff is
discontinuous: the Jacobian is branch-local, not a smooth derivative across it.
The residual, split-state evaluator and analytic chain use the same geometry.
No global default changes. Numerical qualification is recorded separately in
the SimpleMOS branch status. A single edge drive cannot represent this nonlinear
mixture: split edge/port diagnostics report NaN for that field, and the legacy
single-drive factorization probe rejects this profile; use `element_box_probe`.

### SRH doping-dependent lifetime

`solver.srh_doping_dependence` enables the Sentaurus
`SRH(DopingDep)` Scharfetter lifetime law independently for electrons and
holes:

```text
tau(N) = tau_min + (tau_max - tau_min) / (1 + (N / Nref)^gamma)
```

Example matching the BVmethods `models.par` values at 300 K:

```json
"srh_doping_dependence": {
  "enabled": true,
  "concentration_basis": "total_impurity",
  "electron": {
    "tau_min_s": 0.0,
    "tau_max_s": 1.0e-7,
    "reference_doping_m3": 1.0e16,
    "gamma": 1.0
  },
  "hole": {
    "tau_min_s": 0.0,
    "tau_max_s": 1.0e-7,
    "reference_doping_m3": 1.0e16,
    "gamma": 1.0
  }
}
```

`concentration_basis` accepts `total_impurity` (`Nd + Na`) or `net_doping`
(`abs(Nd - Na)`). With `scaling.mode: "unit_scaling"`, numeric
`reference_doping_m3` values are read in the active TCAD concentration unit,
`cm^-3`, consistently with the existing mobility inputs. The compiled
`reference_doping_m3` default is likewise expressed in the internal unit
system (`1e22 m^-3` in `legacy_si`, `1e16 cm^-3` in `unit_scaling`). Omitting
this block preserves the legacy uniform `taun` and `taup` behavior.

For density-gradient runs, `solver.srh_density_coupling` selects the electron
density used only by the SRH denominator and generalized Fermi factors:

```json
"srh_density_coupling": "sentaurus_default"
```

`quantum` (the default) preserves the legacy Vela behavior and uses the
quantum-corrected transport density. `sentaurus_default` follows the Sentaurus
default density-gradient SRH convention: SRH uses the classical electron
density obtained without the electron quantum potential, while Poisson,
transport, contact current, and Auger recombination continue to use the
quantum-corrected density. Both modes use an analytic generalized Fermi-SRH
Jacobian.

With `scaling.mode: "unit_scaling"`, Caughey-Thomas and Masetti mobility
values are read as `cm^2/(V s)`, reference dopings as `cm^-3`, saturation
velocities as `cm/s`, surface reference fields as `V/cm`, and surface theta
coefficients as `cm/V`. They are kept internally in those TCAD units before
mobility evaluation.

## sweep

Required core fields:
- mode: `iv`, `cv_quasistatic`, or `bv_reverse` (aliases: `cv`, `bv`, `reverse_breakdown`)
- contact: swept contact name
- start: number
- stop: number
- step: non-zero nominal target spacing; its sign must move `start` toward
  `stop`
- initial_step: optional positive magnitude of the first internal voltage
  step. It must lie within `[min_step, max_step]`. When omitted it defaults to
  `abs(step)`, preserving legacy sweep behavior.
  - min_step: lower bound for adaptive retries. If shrinking a failed step
    crosses this bound, the controller tries exactly `min_step` once, provided
    `max_retries` permits another attempt. Failure at the bound terminates
    retries. A final or mandatory interval shorter than `min_step` may be
    attempted, but is never enlarged past its target or retried below the bound.
- bias_points: optional array of explicit biases. When present, the adaptive
  `start`/`stop`/`step` stepping loop is bypassed and Vela solves exactly these
  biases in order. Use a one-element array for a single-bias restart.

Output and current fields:
- current_contact
- write_vtk
- vtk_prefix
- vtk_format: `legacy` (default) writes the ASCII `.vtk` file with every
  diagnostic array under Vela's historical names; `vtu` writes an XML
  UnstructuredGrid `.vtu` file (Float64 arrays, `region_id` Int32 cell data,
  `FieldData` with `VelaVersion`, `SweepContact`, `Bias_V`, `PointIndex`,
  `Converged`, per-contact `<name>_TotalCurrent_A_per_um` when
  `csv_all_contacts` is on, and one `Unit:<array>` string per array). Array
  names and units follow the Sentaurus TDR data names of the manual's
  Appendix F (`ElectrostaticPotential` V, `eDensity`/`hDensity` cm^-3,
  `DopingConcentration`/`DonorConcentration`/`AcceptorConcentration`/
  `SpaceCharge` cm^-3, `ElectricField` V/cm, `eCurrentDensity`/
  `hCurrentDensity`/`TotalCurrentDensity`/`ConductionCurrentDensity` A/cm^2,
  `eMobility`/`hMobility` cm^2/(V*s), `eVelocity`/`hVelocity` cm/s,
  `srhRecombination`/`AugerRecombination`/`Band2BandGeneration`/
  `AvalancheGeneration`/`TotalRecombination` cm^-3*s^-1, `eAlphaAvalanche`/
  `hAlphaAvalanche` cm^-1, `eIonIntegral`/`hIonIntegral`/`MeanIonIntegral`,
  `BandGap`/`BandgapNarrowing`/`EffectiveBandGap`/`ElectronAffinity`/
  `ConductionBandEnergy`/`ValenceBandEnergy` eV, `EffectiveIntrinsicDensity`/
  `eEffectiveStateDensity`/`hEffectiveStateDensity` cm^-3, `eEparallel`/
  `hEparallel`/`eEnormal`/`hEnormal` V/cm, `eGradQuasiFermi`/`hGradQuasiFermi`
  V/cm). The authoritative list is `tdrFieldCatalogue()` in
  `src/io/VtuWriter.cpp` (`VtuWriter.h` declares the catalogue API).
- vtk_fields: optional array of TDR names for `vtk_format: "vtu"`. `Name`
  requests the scalar form (the magnitude for vector quantities),
  `Name/Vector` the three-component array; a vector request supersedes a
  scalar request for the same name. An empty or absent list writes the whole
  catalogue. Requested names the build cannot provide are not written; they
  are listed per point in `DCSweepPoint::skippedVtkFields` and logged as
  `vtk_fields_skipped`.
- vtk_binary: `false` (default, ASCII arrays) or `true` (base64 inline binary
  with a UInt32 header) for `vtk_format: "vtu"`.
- vtk_field_data: optional object of string or finite numeric values, for
  example `{"CmdSHA256":"abc","Stage":1,"t":0.25}`. Keys must match
  `[A-Za-z_][A-Za-z0-9_]*`; booleans, arrays and nested objects are rejected.
  Built-in names cannot be overridden: `VelaVersion`, `SweepContact`, `Bias_V`,
  `PointIndex`, `Converged`, `GoalReached`, `goal_reached`, swept-contact
  `<name>_OuterVoltage_V`/`<name>_InnerVoltage_V`, and any
  `<name>_TotalCurrent_A_per_um`. FieldData stays readable ASCII in binary VTU.
  Presence of this object enables per-point `GoalReached`/`goal_reached` (0/1)
  and swept-contact inner/outer voltages. If `t` is not supplied, it is computed
  from the voltage interval, or from the measured initial current and final
  current target for voltage-to-current stages. Goal status is local to this
  stage, using the existing voltage/current tolerance; it does not change
  convergence or stopping decisions. Omitting the object preserves the
  previous VTU bytes and metadata set.
- vtk_bias_points: optional array of biases the adaptive sweep must land on
  exactly (to 1e-12 V). Each becomes an additional accepted point with a VTK
  file even when `write_vtk` is false; the adaptive step size is not changed,
  so all other visited biases are identical to a run without the list.
  Points outside `[start, stop]` or equal to a nominal target are ignored.
  Explicit `bias_points` sweeps do not use this key.
- csv_all_contacts: when `true`, the sweep CSV gains
  `contact_<name>_electron_A_per_um`, `contact_<name>_hole_A_per_um` and
  `contact_<name>_total_A_per_um` for every mesh contact (the monitored
  `current_contact` columns are unchanged).
- csv_file
- initial_state_file: optional `vela.state/2` HDF5 file (`.h5`) used as the initial
  `DDSolution` for the first solved bias.
- write_state_file: optional restart-state file overwritten after every
  converged point with the latest `DDSolution`.
  `state_format` defaults to `hdf5` and accepts only that value. The required
  `VELA_ENABLE_HDF5_STATE=ON` build uses HDF5/HighFive directly, without
  FlatBuffers. CSV, VDS1, NPY and the old `vela.ddstate/1` schema are not runtime
  restart formats. Current cases and explicit paths must use `.h5`.
  Values are binary64, potentials in V and densities in m^-3. Mesh identity
  includes numbering/order, coordinates, connectivity, regions, contacts and
  length units. Loading checks the mesh and potential origin; reference changes
  require an explicit seed translation. QF references and increments are stored
  independently, as are enabled quantum fields. Writes close a same-directory
  temporary file before replacement. The shared four-equation representation
  also requires explicit lattice temperature in K; it never supplies missing
  temperature from a default. See the [storage interface](state_archive.md),
  [production migration contract](validation/templates_ldmos_hdf5_production_migration_plan_2026-09-23.md)
  and [current qualification scope](validation/templates_ldmos_hdf5_migration_execution_2026-09-23.md).
  Explicit split-DD checkpoints additionally preserve the three normalized
  high/low coordinate blocks, their potential scale and split mesh identity.
  These optional DD-only bundles are validated against the physical state;
  partial or inconsistent bundles are rejected. Historical CSV migration is a
  separate utility and does not restore a production CSV fallback.
- `initialization.mode`: optional first-point initialization mode. `none`
  preserves the baseline cold-start path; `poisson_block` runs one Newton
  Poisson block solve before the first coupled Newton solve and uses that state
  as the first-point handoff.
- `initialization.poisson_max_iter`: non-negative iteration budget, default `0`.
  A positive value requires `poisson_block` and solves nonlinear Poisson from
  the cold state with fixed contact-basin quasi-Fermi potentials before handoff.
  It retains the configured solver tolerances and aborts if initialization fails.
  `<sweep csv stem>_poisson_initialization.json` records the budget, iterations,
  residuals and result, including failures. The command frontend preserves the
  preceding `Coupled { Poisson }` budget here when merging equilibrium stages.
  The following frontend Coupled stage uses warm-start Newton to preserve this
  converged handoff instead of performing another Gummel pass.
  Zero retains the historical one-block-update behavior for existing JSON decks.
- `initialization.diagnostic_csv`: optional CSV written only when
  `initialization.mode` is `poisson_block`. It records the first-point
  initialization bias as `bias_V`, plus the cold-state and Poisson-block
  residual norms.
- `initialization.write_state_file`: optional HDF5 restart state written only
  when `initialization.mode` is `poisson_block`. It captures the Poisson-block
  handoff state for the first bias point before ordinary sweep continuation
  begins.

Initialization semantics:
- `sweep.initialization` applies only to the first solved bias point. Later
  points continue from the last accepted sweep state as usual.
- `initialization.mode="poisson_block"` cannot be combined with
  `initial_state_file`; use one first-point initialization source or the other.

### Versioned simulation templates

Qualified PN2D production starting points are stored separately for forward
and reverse operation:

- `configs/templates/pn2d_iv.template.json`: forward 0 to +10 V IV, low-field
  Masetti mobility, `total_impurity`, SRH plus Old
  Slotboom, and impact ionization disabled.
- `configs/templates/pn2d_bv.template.json`: reverse 0---20 V BV,
  `masetti_field`, `total_impurity` (`Nd + Na`), SRH plus Old Slotboom, and
  Van Overstraeten impact ionization.

IV v2 and BV v5 default to compensated potential and concentration with unit
scaling. They use Newton iteration counts for step growth, linear prediction,
and preserve the proposed step after landing on a requested output bias.
The default grids have 201 IV and 401 BV output biases. The numerical defaults
use the CODATA2018 reference profile. For the previous solver and sweep
defaults, render with `--set solver_profile=legacy`.

PN2D BV retains the version-3 atomic
`element_edge_sg_gss_laux` profile. Rendering it binds element-edge
SG/GSS-Laux current support, element-vertex box source mapping, Bernoulli
midpoint density, mixed-Voronoi node volumes, and
`mesh_geometry.require_non_obtuse=true`. Use the
`legacy_cell_reconstructed` profile for the complete cell-reconstructed and
barycentric rollback. The generator rejects partially mixed profiles. This
template policy applies to the BV geometry and avalanche settings.
See the PN2D validation reference for the scope of each template version.

For audited imported meshes that contain obtuse triangles but are Delaunay at
the assembled interior-edge level, `mesh_geometry.edge_coupling_policy` may be
set to `assembled_delaunay_fallback`. This sums the two raw cotangent
contributions before testing their sign, preserving a non-negative assembled
Voronoi dual face when one local contribution is negative. A genuinely
negative assembled face uses the existing positive barycentric fallback.
`local_positive_fallback` remains the parser default and preserves legacy
behavior. The SLOT-LDMOS preparation tool exposes this only through the
independent atomic `element_edge_sg_gss_laux_delaunay` profile, which also
binds mixed-Voronoi node volumes, SG/GSS-Laux current support, element-vertex
box source mapping, and Bernoulli midpoint density. It does not modify the
existing non-obtuse profile or any global default.

The Sentaurus BVmethods NMOS boundary-control equivalents are stored as
separate templates because their sweep voltage has different semantics from a
plain voltage sweep:

- `configs/templates/bvmethods_nmos_external_resistor.template.json`:
  self-consistent avalanche with a `1e7 ohm*um` series-resistor load line.
  `sweep.start/stop/step` are outer circuit voltages; `bias_V` remains the
  solved device/inner voltage.
- `configs/templates/bvmethods_nmos_voltage_to_current.template.json`:
  voltage continuation through `6.0 V`, followed by drain-current targets up
  to `1e-4 A/um`.

Both BVmethods templates freeze the validated NMOS physics combination:
Fermi-Dirac statistics, Old Slotboom with the Fermi correction, Masetti
high-field mobility, SRH, E2 band-to-band generation, and self-consistent Van
Overstraeten avalanche. They do not scale mobility or avalanche parameters.
Their defaults assume a restart state at `5.9 V`; override
`initial_state_file` when materializing a run directory.

Render a runnable configuration and its separate reproducibility manifest with:

```text
python scripts/generate_pn2d_config.py \
  --template pn2d_bv \
  --output runs/pn2d_bv/simulation.json \
  --set mesh_file="inputs/mesh.json" \
  --set node_doping_file="inputs/doping.csv" \
  --set materials_file="inputs/materials.json"
```

For example, render the external-resistor method with:

```text
python scripts/generate_pn2d_config.py \
  --template bvmethods_nmos_external_resistor \
  --output runs/bvmethods_external_resistor/simulation.json \
  --set mesh_file="inputs/mesh.json" \
  --set node_doping_file="inputs/doping.csv" \
  --set materials_file="inputs/materials.json" \
  --set initial_state_file="inputs/bv_prebias_5p9V_state.h5"
```

`--set` accepts only declared template parameters and parses its value as JSON
when possible. Vela resolves relative input and output paths against the
directory containing the generated configuration. Paths are relative by
default so a deck can be moved together with its run directory.
`--allow-absolute-paths` is reserved for legacy or external workflows. The
rendered deck is checked for parameter types, sweep direction, step bounds, and
the template's qualified IV/BV physics combination. The machine-readable final-config schema is
`configs/schema/vela-simulation.schema.json`.

Restart-state CSV files use this exact header:

```text
node_id,psi,phin,phip,electrons_m3,holes_m3
```

Density-gradient checkpoints may append
`electron_quantum_potential_V`. All-material Eq. 231 checkpoints append both
`electron_quantum_potential_V` and
`electron_quantum_potential_like_V`; the latter is the continuous primary
unknown used to reconstruct the material-side quantum correction after a
restart. Readers remain compatible with the legacy six-column and
Lambda-only seven-column forms.

Rows must cover every mesh node exactly once by `node_id`. Restart potentials
are in V and the `electrons_m3`/`holes_m3` fields are always in `m^-3`.
The reader converts densities to the active internal concentration unit.
For checkpoint-style BV runs, set `write_state_file` during the first run, then
resume from the latest converged point with `initial_state_file` plus either a
new adaptive range or `bias_points` for selected target biases. Sweep CSV files
are still opened as fresh outputs; use distinct CSV names for resumed segments.

Diagnostics fields:

`simulation_type: "ialmob_point_probe"` evaluates the production IALMob
formula without a mesh solve. It requires top-level `mobility`, `points_file`
(a nonempty JSON array), and `output_csv`; `temperature_K` defaults to 300.
It accepts explicitly SI inputs and rejects a `scaling` block. Each point
has a unique nonnegative `node_id`, `donors_m3`, `acceptors_m3`,
`electrons_m3`, `holes_m3`, `interface_distance_m`, and
`electron_driving_field_V_per_m`, `hole_driving_field_V_per_m`,
`electron_normal_field_V_per_m`, `hole_normal_field_V_per_m`.
CSV output contains both carriers' final mobility and scattering components
in `cm^2/(V s)`, together with screening and interface weights. Infinite
component mobilities denote a disabled scattering mechanism; the final
mobility must remain finite. Probe success does not qualify an I-V state.
The preparation and comparison scripts in `scripts/*slot_ldmos_ialmob_probe.py`
retain source parameters, field units, node support and averaging limitations.

The read-only `simulation_type: "newton_carrier_term_probe"` accepts a
`state_file` and `carrier_term_probe: {"solved_equation_terms": true}` to
export the actual carrier equation terms at a saved state. Its CSV includes
`electron_continuity_active` and `hole_continuity_active`: exclude rows marked
zero (contacts and inactive gauges) when assessing local continuity. The
`electron_density`/`hole_density` columns use active internal concentration
units. Term/residual columns retain assembly units. A successful probe only
means the diagnostic completed; it does not qualify the state as converged.
`scripts/audit_slot_ldmos_local_continuity.py` summarizes all active rows
without the production source-size filter, including absolute and signed sums
to expose cancellation. Its relative-error counts are diagnostic and do not
replace the configured solver gates.

`simulation_type: "newton_block_step_probe"` also accepts an optional
`output_state_prefix`. It saves each trial as `<prefix>_<block_mode>.h5` and
returns its `trial_state_file` in the block summary. Saved states include the
reference/increment coordinates and SI carrier densities. They are diagnostic
trials and must pass a full coupled validation before use as checkpoints.

The explicit `block_modes: ["carrier_density"]` candidate requires Boltzmann
statistics, no quantum correction, and already projected carrier contact
values. It transforms the full carrier Jacobian into density coordinates,
limits each depletion to 90 percent, then maps the positive candidate back to
quasi-Fermi coordinates. Backtracking requires a decrease of the carrier
residual norm with frozen physical row scales. It is a diagnostic candidate
method, not a convergence criterion or a new default Newton formulation.
When all 20 density backtracks are rejected, the CLI writes
`numerical_failure.json` beside the probe configuration and exits with code 3.
The `vela.block_candidate_failure.v1` report binds the exact configuration,
state and physical input file hashes, and records the initial merit plus each
backtrack's alpha, finite flag, trial merit (null when nonfinite), depletion
clip count and QF cap count. It does not qualify a state. Existing failure
reports are never overwritten. Input/runtime errors retain their other exit
semantics. The SLOT load-line adapter treats code 3 as retryable only after
validating this report and unchanged inputs; missing, stale or inconsistent
reports remain fatal evidence errors. No QF phase switch is implied by this
classification; the existing bounded voltage recovery handles the rejection.
The SLOT load-line controller writes checkpoints by flushing a temporary file
and atomically replacing the destination. Windows replacement errors 5, 32 and
33 receive up to six short retries (1.575 seconds total waiting); the previous
checkpoint is never deleted first. Persistent access errors remain environment
errors, not numerical rejections. Other write/replace errors propagate without
retry. Read-only Windows monitors should permit `FILE_SHARE_DELETE` so they do
not obstruct atomic replacement. Resume still revalidates persisted device and
circuit evidence before advancement.
The load-line CLI additionally offers an opt-in state predictor through
`--state-predictor-max-ratio R`. It extrapolates potential and split quasi-Fermi
references/increments from two already qualified device states. A sufficiently
spaced parent is preferred over a nearly coincident load-line correction point;
the extrapolation ratio is bounded by `R`. Predictions are unqualified initial
guesses, with parent hashes recorded, and must pass the unchanged full coupled,
all-active-row, global-closure and circuit-confirmation checks. A numerical
prediction failure consumes its own attempt and retries the same voltage from
the original qualified parent before the existing voltage-halving recovery.
Environment, evidence and total-budget errors do not enter this fallback.

`--adaptive-outer-initial-step H` enables success-based outer stepping for one
increasing final `--outer-targets` value and an explicit `--initial-outer-bias`.
`--adaptive-outer-growth` defaults to 1.6 and `--adaptive-outer-max-step` to 0.2 V.
Only a fully accepted device/circuit point advances the step policy; numerical
rejections since the last accepted outer point inhibit growth. The existing
inner and outer halving paths remain active, with inner proposals also bounded
by the current outer step and `--max-inner-step`. Step state, prediction options
and the final target are persisted and checked on resume. Both prediction and
adaptive growth are disabled unless explicitly requested; physical defaults and
qualification thresholds are unchanged.
The full SLOT BV campaign (`run_slot_ldmos_bv_campaign.py`) accepts
`--outer-goal` (physical external voltage, default 200000 V) and `--criterion`
(absolute drain current in A/um, default 1e-7). It stops at the first qualified
device/circuit record satisfying either condition and records the stop reason.
Proposals are clipped to the remaining external voltage; an independently
reclosed point exceeding the goal by more than the load-line tolerance is
rejected. These are not internal-voltage targets.
`--initial-loadline-checkpoint` imports a completed, physically matching scalar
campaign only after revalidating its device and circuit evidence. The imported
checkpoint and continuation policy are frozen for resume.
`--state-predictor-max-ratio` enables same-reference-frame qualified-parent
prediction in the low-voltage scalar and high-voltage coupled paths. High-voltage
prediction uses the external-voltage secant ratio; numerical rejection retries
the original qualified parent, while environment/evidence errors propagate.
`--adaptive-initial-step`, `--adaptive-growth` (default 1.4) and
`--low-voltage-max-step` (default 0.1 V) enable success-based growth. At physical
internal voltage 35 V and above, the existing resistor/current-based step bound
also applies. Growth follows independent circuit confirmation; recovered scalar
or high-voltage prediction failures suppress growth. Reference-frame changes
are independently requalified and never use cross-frame predictor parents.
The opt-in campaign `--port-step-control` replaces that legacy cap policy with
three independent limits: `--outer-max-step` (default 0.2 V),
`--accepted-inner-max-step` (0.2 V), and `--trial-inner-max-step` (0.1 V).
The last limit controls scalar device trials and, on the coupled branch, the
internal-voltage update cap; it is not the accepted-state increment limit.
Without this flag, legacy behavior including the 0.02 V coupled update cap is
unchanged. These are campaign CLI options, not new physical JSON deck fields.
The newest three separated qualified points in one reference frame estimate
dI/dVinner by weighted linear fitting (weights 1,2,3). Negative slopes, spacing
below 1e-6 V, insufficient history or adjacent slope disagreement above a factor
of two disable resistor amplification; this is not an NDR continuation method.
External proposals are bounded by the external cap, estimated internal motion,
estimated current growth, next-span feedback and the remaining voltage goal.
The optional `--proposal-max-step` additionally caps each proposed external
increment, including proposals after success feedback or retry. It defaults to
unset, must be finite and positive, and is frozen as `proposal_max_step_V` in
the continuation identity. It does not replace `--outer-max-step`: the latter
still checks the actual independently reclosed external increment. This permits
bounded small-step experiments without rejecting a valid strict correction
merely because it exceeds the smaller proposal cap.
`--step-current-fraction` (0.25) and `--step-current-floor` (1e-16 A/um) bound
absolute current change by fraction * max(abs(previous current), floor).
Actual corrected internal/external/current increments are checked before ledger
acceptance, in addition to every existing device/circuit qualification gate.
Port prediction error and all saved nonlinear history rows (including rejected
trials) plus elapsed time limit growth. `--step-target-iterations` (20) is a
history-row work target, not a solver iteration limit; `--step-target-seconds`
(120) is a time target, not a timeout. Recovery suppresses growth and halves the
next proposal. The feedback factor is limited to [0.5, adaptive-growth].
Limits and controller source hashes are frozen in the campaign identity;
`--accept-driver-update` cannot silently migrate a saved run to a new policy.
Predictions are heuristics, not physical-error bounds or convergence guarantees.
The optional `--actual-outer-continuation` requires port step control. Below
38.4 V internal bias it tries at most one extrapolated device state and one
qualified-parent fallback, then independently confirms the actual supply
coordinate Vinner+R*Id. It does not claim that the originally proposed supply
was solved: both coordinates and their difference are retained. All actual
increment/current limits, goal overshoot guards, device and circuit gates still
apply. Same-seed failed trials within 1e-9 V are skipped within this campaign.
`--candidate-max-seconds` bounds the whole candidate call (including nested
corrections); the global campaign deadline still takes precedence. Optional
`--corrector-stagnation-passes N` (zero/off or N>=2) stops unsuccessful recovery
when N full validations reject line-search decrease, all carrier rows and global
closure pass, and Poisson residuals vary by at most 1e-4 relative. This returns
failure, never qualification, and leaves solver acceptance thresholds unchanged.
These options are frozen in the continuation identity; start a new audited run
when changing them. Old partial ledgers require explicit independent validation
and a provenance-bearing seed export before reuse.
`--reference-origin-floor V` optionally keeps the scheduled numerical origin at
least V after an audited seed-window import. It does not alter physical bias;
all translated states and circuit proofs must be requalified. The option and
its value are frozen in identity. It is not available for an unseeded cold run.
Actual-coordinate candidates also save `qualified_device.json` before outer
acceptance so that a later rejection does not erase their device proof; such a
file alone is not an accepted circuit state.
With `--port-step-control --port-predictor limited_rational`, the baseline is
a limited three-point endpoint tangent of Vinner(Vouter), not an I(Vinner)
regression. Up to six qualified, same-frame points supply three causal one-step
holdouts. A three-point rational candidate is selected only if all holdouts are
eligible, its median normalized error is below 80% of the baseline (above a
1e-10 normalized-error tie floor), and its latest error is no worse. Selection
uses completed steps only. The endpoint slope is limited to [0,min(3*s_last,1)]
on a nonnegative historical branch. Observed negative internal or current
increments break the history instead of being clipped. Insufficient history
uses conservative limits without claiming a measured monotone law.
The rational form is deltaV=a*H/(1+b*H). It must have finite fit parameters,
scaled infinity-norm condition estimate <=1e8, denominator >=1e-3 and derivative
in [0,1] across the fitting and proposed intervals; otherwise the baseline is
used. Bounded bisection enforces predicted internal/current caps. Existing
actual-increment, device, all-row and circuit checks remain mandatory. The
model, eligibility, historical errors and decision are saved in each proposal.
`--port-predictor legacy` remains the default; predictor identity is frozen
with the rest of the policy and cannot be changed by ordinary resume.
The optional `--port-feedback-mode matched_actual_outer` requires
`--port-step-control --port-predictor limited_rational`. The default
`proposed_outer` retains the previous feedback calculation. The new mode
evaluates the already selected frozen model at the strictly accepted actual
outer voltage, then compares its inner voltage and load-line-derived current
with the strict result. Normalizers, cost/recovery factors and all final
qualification/increment gates are unchanged. No new model is selected using
the accepted point. Saved past samples revalidate the same coefficients and
the extended domain; insufficient history, an invalid domain or a reference
frame mismatch falls back to the original error, caps the feedback factor at
1 and prohibits growth-stage promotion (the growth controller backs off).
The proposal records resistance and reference origin. Feedback retains the
original normalized error, proposed/actual supply, supply translation, frozen
model, checked span, matched predictions and promotion eligibility. Candidate
current error and conservation remain separate growth-quality checks; candidate
current never substitutes for strict current. The mode and fallback protocol
are bound as `continuation_policy.port_feedback`, version
`matched_actual_outer_v1`. Changing mode requires a new audited run, not an
in-place resume or `--accept-driver-update` override.
The optional `--row-quality-repair` requires `--current-growth-control` and is
off by default. After strict device qualification, it permits at most one
fixed-potential `carrier_only` correction when the local row ratio exceeds the
growth-promotion threshold, while candidate error/conservation, resolved
progress, actual increments and predictor quality permit promotion. Recovered
or retried candidates are ineligible. This is a growth-quality optimization,
not a relaxation of the final row or convergence gates.
The correction must pass full strict device validation at the same inner
voltage, improve the row ratio and retain the candidate-current quality limit.
A numerical rejection retains the original qualified device point and blocks
promotion; budget, input, environment and evidence errors propagate normally.
Nested multi-pass correction is disabled during this strict revalidation.
Actual supply voltage and current increment checks, growth metrics and circuit
confirmation are recomputed for the selected final point. The complete repair
and revalidation cost, including rejected repairs, enters attempt feedback.
`quality_repair/summary.json` binds the original state and records eligibility,
before/after metrics, status and elapsed time. The option is bound in
`continuation_policy.row_quality_repair` as `one_pass_row_quality_reuse_v2`; changing
it requires a new audited run rather than ordinary resume.
The correction helper's completed full Newton validation and final all-active-row
probe are bound by `SolverBackend.adopt_validation`; a second identical solve
and probe are not run. Reuse checks the physical/qualification configuration,
actual bias, input and output files, full-solve lineage and active-row coverage.
Only iteration budgets and output paths may differ from the base contract.
The resulting proof records the runner and physical input identity.
Successful device and circuit-confirmation checks may be memoized in one process
(at most 128 entries). Every lookup still hashes all bound evidence and physical
resources. Changed contents, metadata, base configuration or input identity
invalidate reuse; missing dependencies are rejected. No persistent cache or
candidate-only success flag can qualify a point. Circuit confirmation remains
separate from full device qualification.

The default-off `--current-proposal-mode current_target` requires the existing
current-growth control and a positive qualified parent current. It is intended
for the near-vertical positive-current branch. Instead of feeding the achieved
voltage span into the next proposal, it proposes `f * R * abs(I_parent)`, where
`f` is the current growth stage times its backoff multiplier. Feedback supplies
only a bounded factor in `[0.5, 1]`. The remaining supply range, explicit retry
and proposal caps, current-scaled outer cap, and unchanged empirical raw
candidate-error cap all still limit this proposal. The final current guard
remains independent (25% in the SLOT experiment); the formula is a constant-inner
voltage forecast, not a guarantee of achieving that current growth.

This mode does not require a monotone inner-voltage port fit. Secant restart
history is selected by increasing outer voltage and current within the same
reference frame; missing history is recorded. The optional
`--accepted-inner-decrease-max-step` (default `0 V`) permits a bounded decrease
of accepted inner voltage, in addition to the existing numerical tolerance.
It must be nonnegative, no larger than the accepted-inner cap, and requires
`current_target`. External voltage and positive-branch current must still
increase. These settings are identity-bound as `current_target_fold_v1`,
`outer_current_v1`, and the nonzero decrease cap; changing them requires a new
audited run. They do not change device/row/closure/load-line gates, candidate
promotion metrics, or the requirement for strict device reclosure.

When `coupled_linearization_audit` is enabled and an output state is configured,
the usual seed diagnostics remain available. Additionally,
`sweep.external_circuit.coupled_diagnostic_min_updates` (default `0`, range
`0..3`, no larger than `max_iterations`) can suppress early convergence for a
bounded diagnostic. It retains the original residual tolerances, update limits
and line search, and may fail before reaching the requested update count.
When nonzero, every accepted update is saved beside `write_state_file` as
`.coupled_update_<outer>_<retry>_<update>.h5`, with diagnostic-only state role.
Every successfully evaluated line-search trial is also saved as
`.coupled_trial_<outer>_<retry>_<iteration>_<line_search>.h5`. Its log records
alpha, in-memory equation norm/index, circuit residual, merit and acceptance.
Rejected trials retain a diagnostic-only role and are never promoted. Export
errors propagate rather than being treated as numerical backtracking failures.
These snapshots grant no device qualification. The default performs no extra
updates and writes no update or trial snapshots.

When `coupled_linearization_audit` is enabled and an output state is configured,
the diagnostic runner saves the initial coupled seed separately, logs its
current and load-line residual, and logs successful as well as failed coupled
iteration traces. The saved seed is not qualified. The frozen
`newton_candidate_conservation_probe` additionally accepts
`candidate_port_diagnostics: true` to report all contact residual-functional
currents (A/um, signed as in the circuit functional) and frozen residual blocks.
The contact functional supports split DD states and restores their referenced
QF and packed low parts before evaluating the unconstrained physical residual.
These are state diagnostics, not an alternative convergence or acceptance gate.

The continuation driver's `--carrier-corrector-mode carrier_density` uses this
candidate and retains the full coupled, global closure and all-active-row gates.

The continuation runner's optional `--audit-local-continuity` adds an all-active
row audit at each accepted candidate, using the deck's `eps_row` and scale floor.
Failed audits stop advancement; the original solver and closure gates still
apply. This option requires per-point continuation and cannot be combined with
`--batch-predictor`. Diagnostic outputs are isolated per attempt. Secant restart
prediction extrapolates saved QF references and increments separately so tiny
leakage-current signals survive a moving bias reference.
`--secant-predictor` applies this prediction to each subsequent point once two
qualified states are available. `--carrier-corrector-passes N` optionally tries
up to N fixed-potential carrier Newton steps after a failed coupled attempt.
These steps use the physical carrier Jacobian without its diagonal floor. The
result must then pass the full original coupled deck and a fresh all-active row
audit. The attempt manifest reports coupled Newton iterations and carrier
correction steps separately; a zero-step validation is not counted as a fresh
coupled Newton iteration. This correction requires `--audit-local-continuity`.
If a coupled validation fails and saves a rejected state, the remaining
correction budget can continue from that state and revalidate again. Each
validation has an isolated directory, and its Newton iterations are included
in the attempt total. Only a complete coupled solve followed by the local
audit advances the accepted state.

The historically named `run_slot_ldmos_noavalanche_continuation.py` preserves
the base deck's impact-ionization configuration. Enabled self-consistent
avalanche uses separate `02_delaunay_avalanche_adaptive` paths and
`adaptive_avalanche_summary.json`; avalanche-off/postprocess-only retains the
existing no-avalanche names. Summaries record the actual physics branch and
source configuration. An enabled external circuit is rejected because the
local audit needs the solved inner voltage, not the outer supply voltage;
retained circuit settings with `enabled: false` are allowed.
Known Windows loader errors (missing/incompatible runtime DLLs) stop the
driver immediately and write `startup_failure.json`. They cannot be repaired
by voltage bisection. Launch Windows runs with the matching UCRT64 `bin` on
`PATH`, including when using a frozen copy of the Release executable.

- `diagnostics.transport.enabled`: appends aggregate mobility and high-field
  drive columns to the main sweep CSV.
- `diagnostics.contact_edge.enabled`: writes a separate contact-edge transport
  CSV. Optional `contacts` selects contacts and `csv_file` overrides the default
  `<sweep csv stem>_contact_edges.csv`.
- `diagnostics.continuity_balance.enabled`: writes contact-adjacent continuity
  balance rows. Optional `contacts` selects contacts and `csv_file` overrides
  the default `<sweep csv stem>_continuity_balance.csv`.
- `diagnostics.terminal_balance.enabled`: writes per-contact terminal current
  balance rows. Optional `contacts` selects contacts and `csv_file` overrides
  the default `<sweep csv stem>_terminal_balance.csv`.
- `diagnostics.poisson_dirichlet_reaction.enabled`: writes the physical
  unconstrained Poisson reaction at every node of each selected electrostatic
  Dirichlet contact, together with the summed contact charge in `C/m` of
  out-of-plane depth. Positive reaction denotes positive electrode charge.
  This reads the production mixed-material Poisson row before Dirichlet row
  replacement; it does not reconstruct a P1 boundary gradient. Optional
  `contacts` defaults to the swept contact and `csv_file` overrides
  `<sweep csv stem>_poisson_dirichlet_reaction.csv`.
- `diagnostics.srh_balance.enabled`: integrates the net, generation, and
  recombination parts of SRH over triangular cells in `material` (default
  `"Si"`). It compares the generated current with the electron current at
  `drain_contact` and the hole current at `substrate_contact`, and sums
  `kcl_contacts` for a port-conservation residual. A point is written as
  `numerically_unresolved` unless the magnitude of the drain total current is
  at least `resolution_margin_ratio` (default `10`) times the four-terminal
  KCL residual. Optional `csv_file` overrides the default
  `<sweep csv stem>_srh_balance.csv`.
- `diagnostics.sg_avalanche_edges.enabled`: writes C++ assembled SG
  edge-current avalanche source rows for `impact_ionization.generation:
  "current_density"` with an SG edge-current `current_approximation`
  (`"density_gradient"`, `"grad_qf"`, `"cell_reconstructed"`,
  `"cell_current_reconstructed"`, `"cell_vector_current_reconstructed"`, or
  `"conserved_total_current"`).
  Optional `csv_file` overrides the default
  `<sweep csv stem>_sg_avalanche_edges.csv`.
- `diagnostics.path_ionization_integrals.enabled`: traces monotone-potential
  field paths through local electric-field maxima and evaluates the
  Sentaurus Device Eq. 469/470 electron- and hole-injection integrals. Paths
  are ranked by `mean_ionization_integral`. `max_paths` controls the number of
  reported paths (`0` writes all); `break_rank: 3` and `break_value: 1.0`
  reproduce `BreakAtIonIntegral(3 1.)` stopping semantics. The path diagnostic
  has an independent `driving_force`, defaulting to `solver` so the path uses
  the avalanche driving force of the solved deck. Set it explicitly to
  `electric_field` for the field-only approximate-breakdown analysis
  recommended by the Sentaurus manual; this does not change the avalanche
  source driving force used by the continuity equations. Optional
  `stop_field_V_per_m` terminates both sides of a traced path when the
  electrostatic field falls below the specified non-negative SI threshold;
  the default `0` preserves the full monotone path. This is an explicit Vela
  depletion-support diagnostic because Sentaurus does not expose a named
  stop-field keyword in the 2018.06 command syntax. Optional
  `electron_stop_field_V_per_m` and `hole_stop_field_V_per_m` independently
  restrict the contiguous carrier support around each carrier's strongest
  driving field on that same geometric path. The coupled electron/hole alpha
  pair is retained inside each interval when evaluating Eq. 469 or Eq. 470;
  both defaults are `0`, which preserves the shared path interval.
  `mean_definition` is `carrier_integral_arithmetic` by default. The optional
  `carrier_alpha_length_arithmetic` diagnostic instead ranks
  `0.5*(integral(alpha_n ds)+integral(alpha_p ds))` on the independently
  supported intervals. It is exposed for Sentaurus `MeanIonIntegral` audits,
  not asserted to be Synopsys's undocumented strong-ionization formula.
  `break_ordering` defaults to `path_mean`. The alternative
  `carrier_integrals` flattens electron and hole values from every path before
  applying `break_rank`; it is an explicit ordering experiment and is not the
  validated default for the BVmethods deck.
  Optional
  `tracing_mode: "continuous_cell"` follows the barycentrically interpolated
  nodal electric field through triangle interiors and samples electron/hole
  coefficients at each triangle's largest-angle (Sentaurus best-vertex)
  node. Nodal electric-field and quasi-Fermi current directions use
  element-vertex-angle recovery. `seed_mode` defaults to
  `nodal_local_maxima`, using maxima of the recovered P1 field; the diagnostic
  alternative `cell_local_maxima` uses raw element maxima. `path_retention`
  defaults to `distinct_local_maxima`: each local-maximum seed contributes one
  strongest trajectory. `numbered_peak_groups` is the BVmethods WriteAll
  policy: nodal maxima separated by at most two mesh edges remain separate
  path numbers but share the strongest trajectory of that P1 peak group.
  `all_seed_trajectories` is an audit mode that retains every distinct
  incident-cell launch from a nodal seed. The compatibility value
  `corridor_deduplicated` also merges adjacent/two-ring maxima sharing a path
  corridor.
  The legacy `edge_graph` mode remains the default. `seed_field_V_per_m` can
  qualify local nodal-field maxima independently of the path stopping field.
  For continuous-cell paths, `tracing_vector` may be `electric_field`
  (default), `electron_current`, `hole_current`, `electron_qf_gradient`,
  `hole_qf_gradient`, `cell_electric_field`, `electric_field_rk4`, or
  `sentaurus_eparallel_adaptive`. The `cell_electric_field` policy traces the
  element-constant P1 potential gradient while retaining nodal local maxima
  for seed discovery. `electric_field_rk4` follows the continuous nodal P1
  electric field with subcell RK4 integration instead of a single straight
  chord per triangle. The adaptive policy
  follows majority-carrier current for interior peaks and the minority-carrier
  quasi-Fermi gradient for peaks on a transport boundary. This preserves bulk
  current paths while preventing low-current surface paths from falling back
  to a short electrostatic-field tangent. Its
  `tracing_qf_relative_floor` (default `5.1e-3`) requires the selected
  carrier's nodal quasi-Fermi gradient to be observable relative to its global
  maximum; otherwise that surface seed follows the electric field. The fixed
  current choices use
  reconstructed SG current; quasi-Fermi-gradient choices are explicit
  diagnostics. Seed and stop
  thresholds continue to use electrostatic-field magnitude. For SG-current
  tracing, `tracing_current_relative_floor` defaults to `1e-8`; nodal currents
  below that fraction of the global current maximum use the electric-field
  direction. The adaptive policy applies that current rule to interior peaks
  and directly selects a minority-carrier quasi-Fermi family for boundary
  peaks and their two-ring P1 aliases.
  `tracing_direction` is `bidirectional` by default; carrier-current path
  searches may select `along_vector` or `opposite_vector` for one-way carrier
  injection trajectories.
  Continuous-cell tracing supports `driving_force: "electric_field"` and
  `driving_force: "eparallel"`.
  Optional
  `csv_file` overrides
  `<sweep csv stem>_path_ionization_integrals.csv`.
- `diagnostics.path_ionization_integrals.segments_csv_file`: optional ordered
  one-row-per-segment trace containing path coordinates, endpoint potential,
  local electric field, electron/hole alpha, alpha-ds, cumulative alpha-ds,
  and prefix/full-path ionization integrals.  If omitted it defaults beside
  `csv_file` with a `_segments.csv` suffix.
- `diagnostics.bv_process_probe.enabled`: writes a normalized solver-used
  avalanche process record for every active edge/cell support and carrier.
  Records include endpoint state, electric/QF-gradient vectors, low- and
  high-field mobility stages, native or reconstructed current provenance,
  alpha, carrier-split generation, source measure, qG, residual scatter,
  every selected branch flag, and deterministic configuration/branch
  fingerprints. The writer checks its summed source against the production
  assembled source to `1e-12` before emitting rows. Optional `csv_file`
  overrides `<sweep csv stem>_bv_process_probe.csv`.
- `diagnostics.avalanche_internal_source_current_audit.enabled`: writes the
  internal SG edge-current avalanche source terms used by assembly, including
  `Fn/Fp` in `V/cm`, `alpha_n/alpha_p` in `cm^-1`, `Jn/Jp` in `A/cm^2`,
  `Gava` in `cm^-3 s^-1`, 2-D contribution area in `cm^2`, and qG contribution
  in `A/um`. Optional `csv_file` and `summary_file` override the default
  `avalanche_internal_source_current_audit.csv` and
  `avalanche_internal_source_current_audit_summary.md` next to the sweep CSV.
  This diagnostic is an internal source audit; exported node current density is
  a separate post-processing quantity.
- `diagnostics.release_bv_config_audit.enabled`: writes a per-bias BV parity
  metadata CSV and a Markdown summary without changing the solve. It records
  the resolved avalanche model, `driving_force`, `parameter_set`, `A_scale`,
  `B_scale`, `switchField`, cutoff/minimum field, RefDens values,
  `source_mapping_mode`, source support/lambda description, 2-D current
  normalization, qG full/junction-window integrals in `A/um`, terminal current
  in `A/um`, max electric field in `V/cm`, max `Gava` in `cm^-3 s^-1`, and
  convergence. Optional `csv_file` and `summary_file` override the default
  `release_bv_config_audit.csv` and `release_bv_config_audit_summary.md` next
  to the sweep CSV. Optional `diagnostic_reference_A_scale`,
  `diagnostic_reference_B_scale`,
  `diagnostic_reference_source_mapping_mode`,
  `diagnostic_reference_qG_full_A_per_um`, and
  `diagnostic_reference_qG_junction_A_per_um` let the summary compare the
  current release run against a known A2/B1.05 diagnostic reference.
- `diagnostics.newton_history.enabled`: writes a complete nonlinear transition
  trace without changing the solve. The compatibility `csv_file` contains
  accepted Newton iterations for converged points. `attempts_csv_file`
  contains every accepted or rejected continuation/Newton attempt, including
  parent bias/state hash, requested and actual target, retry number, typed
  reason, solver/handoff stage, residuals, clamp/damping summary, and final
  state hash. `iterations_csv_file` contains initial, accepted, and rejected
  Newton rows with raw and row-scaled equation-block residuals, update and
  line-search metrics, carrier-row convergence, top residual node per
  equation, and the solver-used source/Jacobian active-branch fingerprint.
  Defaults are `<sweep csv stem>_newton_history.csv`,
  `<sweep csv stem>_newton_attempts.csv`, and
  `<sweep csv stem>_newton_iterations.csv`.
- `write_state_every_point_prefix`: writes every accepted continuation state
  as `<prefix>_bias_<encoded bias>.csv`. These states are restartable with
  `initial_state_file`; use this together with the nonlinear trace to seal
  exact parent states for failed-transition reproduction.

For deterministic avalanche-source diagnostics,
`solver.impact_ionization.source_jacobian` accepts `local_ad`,
`finite_difference` (default), or `frozen`. For the canonical legacy
`triangle_gss_gradqf_truncated` profile, `local_ad` propagates the nine local
Tri3 potential/quasi-Fermi derivatives through the impact coefficient,
midpoint density, high-field mobility, edge GradQF magnitude, truncated source
volume, and symmetric source partition. `finite_difference` remains an audit
oracle. With surface mobility, the cell source uses the complete local
central-difference path even when `local_ad` is selected, because the current
AD kernel does not propagate the surface-field mobility response. Trial
potentials refresh the normal field; trial densities follow the selected
carrier statistics. Edge sources also include adjacent-vertex electrostatic
columns and cross-carrier density dependence. `frozen` keeps the
self-consistent avalanche residual but omits only
its Jacobian contribution; it is an ablation control and does not match the
default Sentaurus avalanche-derivative semantics.

### External circuit and voltage-to-current boundary control

`sweep.external_circuit` enables a Sentaurus-compatible 2-D series-resistor
load line. The supported mode is `series_resistor`. `resistance_ohm_um` uses
`ohm*um`, terminal current uses `A/um`, and the closed equation is

```text
OuterVoltage = InnerVoltage + current_direction * Iterminal * Rseries
```

`bias_points` (or `start/stop/step`) are outer/source voltages. The
targets may increase or decrease monotonically, including negative PN2D anode
biases. A reversal of sweep direction within one target list is rejected.
The nonlinear
device boundary remains the inner voltage and is solved repeatedly until
`load_line_residual_V` meets `residual_tolerance_V`. `current_direction` is
`+1` for the BVmethods NMOS drain convention and may be `-1` for an oppositely
oriented terminal. `initial_inner_voltage_V`, `max_inner_voltage_step_V`,
`voltage_tolerance_V`, `max_bracket_steps`, and `max_iterations` control only
the scalar boundary solve.

For `nested_scalar`, the outer targets are enumerated from `bias_points` or
`start/stop/step`. By default they retain the direct target schedule, and the
first bracket probe can move by `max_inner_voltage_step_V` even when the
requested outer-voltage increment is smaller.
`boundary_control.adaptive_outer_continuation=true` inserts adaptive outer
steps using `initial_step`, `min_step`, `max_step`, `growth_factor`,
`shrink_factor`, and `max_retries`. Growth uses the fixed multiplier.
The first inner prediction is bounded by both the outer increment and the
maximum inner step; subsequent predictions use accepted load-line history.
Only a device-converged, load-line-closed root advances the accepted state.
Failed trials shrink the outer step; exhaustion reports the last failure.
The terminal-current stop criterion applies to both schedules. In adaptive
outer mode the main point's attempted/accepted step diagnostics are outer
increments; `bias_V` and the nonlinear attempt trace remain inner voltages.
The flag defaults to false in JSON and is rejected outside a `nested_scalar`
external circuit. The SDevice frontend enables it for resistor QS decks.
`boundary_control.adaptive_device_continuation` separately controls device
continuation inside a bracket evaluation.
Use `boundary_control.evaluation_csv` and `solver.local_update_diagnostics`
to inspect bracket evaluations. `sweep.diagnostics.newton_history` now also
records these device trials, including failed predictors before fallback,
outer retry numbers, and rejected parent/initial/final states when configured.
Successful bracket evaluations are device trials, not independently accepted
outer/load-line points. Earlier behavior and matched-bias probes are documented
in the [BVmethods controller comparison](validation/bvmethods_loadline_biased_derivatives_2026-10-08.md).

`external_circuit.solver` defaults to `nested_scalar`. The opt-in
`coupled_newton` mode augments the device unknown vector with the inner
terminal voltage and solves the device equations together with
`InnerVoltage + Rseries * Iterminal - OuterVoltage = 0`. It supports
`coupled_equation_tolerance`, `current_directional_step`,
`coupled_damping_factor`, `coupled_max_line_search_steps`, and
`coupled_apply_device_update_limit`. `coupled_line_search_mode` is `merit`
(strict maximum normalized residual decrease) or `residual_filter` (one
equation block may decrease while the other remains inside a bounded
non-divergence envelope); both modes retain the same final device and
load-line convergence tolerances. The filter controls are
`coupled_filter_gamma` and `coupled_filter_envelope_factor`. The envelope is
applied independently to the device-equation and load-line residual blocks;
a large load-line residual cannot authorize growth of the device residual.

Coupled circuit and pseudo-arclength restarts use their actual initial contact
bias and one seed-aware QF reference coordinate system for the device, current,
and state conversions. When referenced QF coordinates are enabled, packing
uses the saved reference/increment fields instead of rounded absolute QFs.
References stay fixed during that continuation context;
this does not guarantee sub-ULP carrier resolution after a large bias change.

`coupled_positive_density_updates` (default `false`, `coupled_newton` only)
uses a nonlinear positive-density trial path for classical Boltzmann carriers.
The full device and current Jacobians are unchanged. The linearized relative
density step is bounded below by -0.9, then mapped back with `log1p`; configured
quasi-Fermi caps apply to the mapped components instead of uniformly scaling
the circuit voltage with the largest raw quasi-Fermi step. This trial path
retains the original Newton tangent at zero step length. Each trial still
passes the selected residual line search and final equation/load-line tests.
Fermi-Dirac statistics and electron quantum potential are rejected. This
experimental globalization option does not replace carrier-row or global
current-closure qualification of a simulation checkpoint.
The opt-in `report_only` mode evaluates the `merit` backtracking sequence
starting at alpha=1 and records `reported_alpha` (zero if no trial qualifies),
but applies the full alpha=1 step even when merit increases. It ignores
`coupled_damping_factor`; separate voltage and device update limits still
apply. A nonfinite or model-invalid full step fails the corrector without
substituting a damped trial. Final residual, load and physical acceptance
gates are unchanged. The coupled update diagnostic CSV also records the full
step merit and raw/limited update sizes. Shadow trial evaluations add cost;
this mode is intended for controlled nonlinear investigations, not a new default.

`full_step` evaluates and applies only alpha=1, without shadow backtracking.
It retains the same independent update limits and invalid-state/terminal gates.
The BVMethods CMD profile now defaults to this mode for Sentaurus
`LineSearchDamping=1` (also the vendor default) and a `NotDamped` threshold
at least as large as the effective iteration budget. Math settings are parsed
and local Coupled settings take precedence. Enabled vendor line search or
Bank–Rose damping is rejected by this profile because those algorithms are
not implemented; it is not silently mapped to Vela's merit rule. Explicit
sidecar modes remain diagnostic overrides. The global C++ default is unchanged.

Scalar coupled convergence is followed by a physical acceptance check at the
actual packed state and inner bias. The configured `carrier_row_convergence`
and `global_continuity_closure` modes apply here as in the standard Newton path;
an enforced failure rejects the point before it advances continuation or writes
an accepted checkpoint. All terminal currents must be finite and satisfy
`abs(sum(I)) <= max(1e-12, 1e-5 * abs(Id))` in A/um. The circuit residual is also
recomputed from the actual extracted terminal current with the configured
`residual_tolerance_V`. For `coupled_newton`, the explicit
`load_line_convergence_mode` selects `enforce` (default) or `report`. Report
keeps the scalar load-line equation, its Jacobian, merit contribution and
measured residual. The main residual branch always checks the complete
bordered residual: `max(device_norm / effective_device_tolerance,
abs(load_residual_V) / residual_tolerance_V) <= 1`. Report removes only the
additional independent load tolerance gate, including final extracted-current
acceptance; it does not remove the circuit row from the main residual branch.
The OR update branch still requires a valid update actually applied during the
current solve. Its density RMS includes the active circuit voltage through the
contact aliases, scaled by `contact_reference_V`; no extra sample is added.
A new target may accept zero updates only if its complete coupled residual is
already satisfied. Nonfinite load residuals, terminal
KCL and configured enforced carrier/closure checks still reject the point.
`nested_scalar` requires `enforce`. No scalar or inexact device tolerance
bypasses the remaining checks. The point's carrier/closure and block residual fields contain the actual
measurements; `validation_diagnostics` includes the JSON
`coupled_acceptance_checked` flag, enforcement status and KCL/load results.
Global closure retains the standard source-floor qualification rule; the JSON
records each carrier's `global_*_qualified` flag. An unqualified carrier does
not establish a strict independent ratio pass. Commas and quotes in diagnostics
are escaped as CSV fields so subsequent columns remain aligned.

`sweep.external_circuit.coupled_step_control` selects the coupled resistor's
outer-voltage stepping policy:

- `quasistationary`: use outer voltage only. The controller takes the ordinary
  sweep's `initial_step`, `min_step`, `max_step`, `growth_factor`, and
  `shrink_factor` (all steps in outer volts), unless explicit `coupled_*outer*`
  controls are supplied. Growth uses the accepted augmented Newton update count
  and `external_circuit.max_iterations`; failure divides the step by Decrement
  (`shrink_factor` is its reciprocal). If the divided step is below MinStep the
  sweep exits, rather than inventing an extra trial at the floor. There is no
  fixed Newton inner-voltage limit, predicted-inner limit, or extra retry at the
  same outer target from the accepted state. The default retry budget is
  unbounded by a separate count; an explicit `coupled_max_step_retries` remains
  an optional work safeguard. Physical/linear/nonlinear acceptance is unchanged.
- `legacy_inner_guard`: retain the existing explicit inner-voltage corrector,
  predictor and retry controls described below.

Unguarded coupled JSON decks default to `quasistationary`. To preserve historical
JSON trajectories, a missing selector resolves to `legacy_inner_guard` when
`max_inner_voltage_step_V` is present or `coupled_continuation.enabled` is true.
An explicitly selected `quasistationary` rejects these conflicting settings,
including `inner_step_control` and enabled `retry_from_accepted`; remove them
rather than silently ignoring them. `nested_scalar` is unchanged and cannot
select a coupled policy. This selector does not change device physics, state
coordinates, density protection, or convergence tolerances.

The qualified BVMethods CMD profile emits `quasistationary` by default. Its
Quasistationary dimensionless steps are converted by the **external** voltage
span, and Decrement is converted to its reciprocal. For a CMD sweep using
`coupled_newton`, the sweep's local `Coupled(Iterations=...)` budget takes
precedence over `Math.Iterations` and supplies `external_circuit.max_iterations`.
An explicit sidecar circuit budget takes precedence over both. An existing sidecar that
requests inner guards must explicitly select `legacy_inner_guard`. This is a
compatibility contract for documented stepping, not a promise of identical
Sentaurus nonlinear trajectories or proprietary solver internals.

For example, to migrate a coupled deck:

```json
"external_circuit": {
  "solver": "coupled_newton",
  "resistance_ohm_um": 1e7,
  "coupled_step_control": "quasistationary"
}
```

Keep the existing physical/acceptance settings, and remove
`max_inner_voltage_step_V` and the old `coupled_continuation` object. Each
accepted coupled point records `external_step_control` and
`newton_inner_voltage_limit_enabled`; a disabled predicted limit is recorded
as zero, not as an active zero-volt constraint. Sweep steps are separate from
Newton line search and from the update RMS stopping criterion.

`sweep.external_circuit.coupled_continuation` explicitly enables a numerical
continuation option for `solver: coupled_newton`:

```json
"coupled_continuation": {
  "enabled": true,
  "retry_from_accepted": true,
  "inner_step_control": {"max_step_V": 0.02}
}
```

After a failed secant seed, `retry_from_accepted` retries the same external
voltage from the last accepted state before shrinking the external step.
Every attempt retains the same Newton, load-line and physical acceptance gates.
The accepted inner/outer voltage secant bounds the predicted inner increment.
The cap is resolved from `inner_step_control.max_step_V`, or a normalized
`max_step_fraction` with explicit device-voltage endpoints, then an explicitly
provided `max_inner_voltage_step_V`. A resistor sweep's external `MaxStep`
does not uniquely specify a device-voltage cap. When no applicable value is
available, the fallback is 0.02 V; this is an empirical default, not a universal
device scale. Its source and fallback reason are recorded in each accepted
point's `validation_diagnostics`. See the
[resolution contract](bvmethods_operator_contract.md) for details. A cap
requiring a step below the configured external minimum terminates the advance;
it is not silently exceeded. This option is disabled by default.

`sweep.external_circuit.coupled_update_convergence` is an opt-in classical
isothermal DD stopping experiment for `solver: coupled_newton`:

```json
"coupled_update_convergence": {
  "mode": "report",
  "digits": 5,
  "potential_reference_V": 0.025852,
  "electron_reference_cm3": 1e10,
  "hole_reference_cm3": 1e10,
  "contact_reference_V": 0.025852,
  "diagnostic_csv": "density_update.csv"
}
```

Modes are `off` (default), `report` (observations without changing stopping),
`enforce` (device residual **or** normalized density-update RMS < 1), and
`enforce_joint` (device residual **and** normalized density-update RMS < 1).
Joint mode permits a residual-only zero-voltage initial anchor before any update;
all nonzero accepted points require a valid raw-update report. This additional
explicit mode retains the historical `enforce` semantics. In particular, a
residual-only pass with a large reported RMS is not a joint qualification pass.
Both branches retain the configured load-line mode and terminal KCL checks;
configured enforced carrier-row/global-closure checks also still apply.
No conservation mode is automatically changed. To study Newton stopping
separately, explicitly configure the conservation checks as `report` in the
experimental input; this does not qualify a formerly rejected state under
the old enforced contract.

The RMS numerator uses the **raw linear direction before all update
limits and backtracking**. QF mode transforms it analytically through the
assembler's actual carrier statistics; `potential_density` uses δn/δp directly,
without a QF transformation or QF voltage cap. Its denominator uses the accepted
physical state plus the declared reference. Potential entries include every
node, carrier entries include free semiconductor nodes, and the active
contact contributes one voltage entry per contact-node alias. These are metric
entries, not independent linear unknowns. Density references are always cm^-3;
potential/contact references are always V. `digits` must be in [0,15], with
finite positive references. Split DD and quantum formulations are unsupported.
No accepted update, a failed evaluator or nonfinite operands cannot authorize
the update branch; tiny applied steps do not shrink its raw numerator.

This follows the observed BVmethods factor-one density/potential error
composition. The error metric alone does **not** implement native iteration,
protection, damping or RHS scaling; primary coordinates and population protection
are selected separately by `state_variables` and `population_update`.
The residual branch remains Vela's device infinity norm and
`coupled_equation_tolerance`, not Sentaurus `RhsMin`. Accepted points identify
`coupled_density_update` versus `coupled_residuals`; the diagnostic CSV records
per-block RMS, counts and residuals for all attempts. Relative CSV paths resolve
from the config directory. See the [observed native contract](validation/contracts/bvmethods_update_error_s0.v1.json).
In that CSV, `iteration=k` identifies the current accepted state (zero is the
seed). For k>0 its RMS uses the raw direction computed at k-1 and the state at
k; `accepted_alpha` describes the next transition k to k+1. It is not the alpha
of the RMS numerator. `convergence_reason` is the final stopping reason of the
whole attempt, repeated on its rows, not a per-row declaration of convergence.

With global closure reporting enabled, the point CSV additionally records
`global_{electron,hole}_{absolute_mismatch,net_scale,contact_node_flux_abs_sum,source_abs_sum}`.
These retain the assembler's continuity-row units. The contact-node sum is not
the sum of absolute terminal currents. Coupled `validation_diagnostics` separately
records the latter as `terminal_{electron,hole}_current_abs_sum_A_per_um`.
The original net-source-normalized ratio and its qualification remain unchanged.
See the [BVMethods stopping control](validation/bvmethods_stopping_control_2026-10-24.md)
for the same-state residual-floor test, physical continuity-unit reconstruction,
and bounded Release scans. These do not establish full native RHS equivalence.
With `newton_history.rejected_state_directory`, physically rejected coupled
candidates are saved as `coupled_physical_<attempt>.h5`; these are never seeds
for accepted continuation. Earlier CSVs without the checked flag do not prove
carrier or global closure qualification.

`coupled_inexact_device_forcing` is an optional object for locator or branch-
tracking runs. It is disabled by default. When enabled, the device tolerance
remains `coupled_equation_tolerance` while the scalar row is far from closure.
Inside `load_activation_ratio * residual_tolerance_V`, it grows geometrically
toward `max_equation_tolerance` as the scalar residual closes. The same dynamic
tolerance is used by convergence, normalized merit, and residual-filter
globalization. A production endpoint should normally be reclosed with forcing
disabled and the strict device tolerance. Example:

```json
"coupled_inexact_device_forcing": {
  "enabled": true,
  "max_equation_tolerance": 5e-6,
  "load_activation_ratio": 100.0,
  "zero_converged_residual": true
}
```

When `zero_converged_residual` is true and the device block is already inside
the current effective inexact tolerance, the bordered tangent solve sets the device right-
hand side to zero while retaining the full device Jacobian. This lets the
scalar circuit row move along the local device manifold without repeatedly
polishing a device residual that already meets the locator tolerance. If a
trial leaves that active tolerance, the next iteration restores the device
right-hand side instead of suppressing the normal correction.

The SLOT campaign driver `scripts/run_slot_ldmos_bv_campaign.py` exposes
`--candidate-zero-converged-device-rhs` as a default-off candidate-only option.
In its high-voltage coupled branch (physical internal voltage at least 38.4 V),
this sets the maximum forcing tolerance to the original `1e-8` and enables
`zero_converged_residual`, without relaxing the tolerance. It also applies to
the qualified-parent fallback after a failed predicted candidate. Lower-voltage
scalar/device continuation, final independent device reclosure, and circuit
confirmation retain their original settings. This option is bound to the
campaign continuation identity: changing it requires a fresh output directory
with a revalidated seed window, not `--accept-driver-update` on the old ledger.

`coupled_linear_solver` accepts `direct_bordered` and
`direct_bordered_qr`. The QR spelling forces the rank-revealing sparse QR
backend first and falls back only if that factorization fails; the non-QR
spelling uses sparse LU first. `coupled_linearization_audit` is an opt-in,
diagnostic-only check that reports device and circuit linear solve residuals,
the analytic versus central-finite-difference terminal-current directional
derivative, and the full device JVP error (including the largest-error row).
It is intended for bounded developer probes because each audited Newton step
evaluates extra residuals.

`coupled_state_predictor="tangent"` optionally builds a full-state supply
tangent at the restart/current state using the bordered Jacobian and RHS
`[0; target_outer-parent_outer]`; the default is `"none"` (existing behavior).
The parent must satisfy the candidate device tolerance. The whole increment
is limited, and a bounded merit check falls back to the unchanged parent if
prediction fails. The predictor runs even for a single direct target. It is
only an initial guess, not device or load-line qualification. Diagnostics are
logged as `external_circuit: tangent_predictor`. The SLOT campaign exposes
`--candidate-state-predictor tangent`, bypasses Python secant extrapolation
for that candidate and binds `bordered_tangent_v1` into continuation identity.
Use a fresh evidence directory and revalidated parents when changing this policy.
The C++ developer API additionally provides
`CoupledLoadLineConfig::directRefinementSteps` (default `0`, range `0..8`).
It is not a JSON/frontend option. An explicit direct solver is required when
enabled. Corrections use a multiprecision defect of the original unscaled
bordered system and the two-sided scaled LU solve. Each retained correction
must reduce the worst componentwise backward error; the result must satisfy
`1e-10` and the existing scaled linear-residual bound. Failure is reported,
and nonlinear convergence, line search, physical gates and default profiles
remain unchanged. This control does not certify the device Jacobian or a BV
trajectory.

`sweep.external_circuit.coupled_state_predictor: "constant"` explicitly starts
each coupled target from the last accepted state, suppressing the built-in
secant extrapolation. The default `"none"` retains the existing secant policy;
the unrelated top-level sweep predictor does not control this coupled branch.
`"tangent"` retains its existing guarded prediction and accepted-state fallback.
When resuming, bind `initial_inner_voltage_V` to the saved device voltage,
separately from `sweep.start`, which is an outer supply voltage.

`sweep.external_circuit.coupled_preserve_contact_potential` (boolean, default
`false`) retains the contact potential rows and columns in direct bordered
`coupled_newton` solves with `solver.state_variables: potential_density`.
The carrier constraints still undergo exact affine substitution, including
their coupling to the retained potential and external circuit voltage.
This option can match a reference solver's system layout in a matrix audit.

For bounded backend diagnostics, set `VELA_PARDISO_INPUT_CAPTURE_DIR` to an
existing, empty output directory. The MKL adapter records the first four
solves of dimension at least 30,000 per process before it calls phase 33.
Each `VELAPD01` file contains the two little-endian 64-bit integers
`dimension` and `stored_nonzeros`, followed by one-based 32-bit CSR row
pointers and columns, binary64 coefficients, and the binary64 RHS.

`sweep.external_circuit.coupled_extended_precision_fallback` (boolean, default
`false`) is an independent direct-bordered recovery option. If the ordinary
solve and refinement fail the existing original-system linear gates, it tries
`sparselu_long_double` with up to eight corrections. It requires
`coupled_linear_solver: "direct_bordered"` and a platform with a wider
`long double`. Only a direction improving the original-system defect score is
retained; the `1e-10` componentwise gate and scaled absolute gate are unchanged.
Assembly, affine constraint substitution, public directions and nonlinear
states remain binary64. This is not full extended-precision TCAD. Update
diagnostics record whether recovery was selected and its correction count.
It adds no new nonlinear stopping rule and changes no default profile.
The persistent `sparselu_long_double` solver reuses factors for bitwise
identical matrices with different RHS. Changed values invalidate numeric
factors; changed sparsity invalidates symbolic analysis. This also applies
after the common prescribed-row substitution, and does not cache physical
evaluations or carry factors between different Newton matrices.

For isolated failure investigations, `VELA_BORDERED_FAILURE_DUMP` may name a
private output directory. A rejected componentwise solve then writes its raw
augmented matrix, RHS, direction, scales and row defects in a fresh child
directory. Export is disabled when the variable is unset; it cannot accept a
failed direction or alter thresholds. Keep generated systems outside Git.

Coupled outer-voltage continuation is controlled by
`coupled_initial_outer_step_V`, `coupled_min_outer_step_V`,
`coupled_max_outer_step_V`, `coupled_outer_growth_factor`,
`coupled_outer_shrink_factor`, and `coupled_max_step_retries`. Failed trials
do not replace the last accepted coupled state. Every accepted adaptive
intermediate point is recorded and is eligible for the normal DCSweep state
checkpoint/restart outputs.
With `sweep.step_growth_mode: "newton_iterations"` and `solver.method: "newton"`,
the coupled controller uses its actual accepted update count and
`external_circuit.max_iterations` as the work budget. Its increment is
`coupled_outer_growth_factor`; all step magnitudes remain outer volts.
The ordinary device seed solver's iteration budget is not used for growth.
The ordinary sweep/legacy default remains `fixed`; the `quasistationary` policy
selects `newton_iterations`. Current break criteria stop coupled continuation
after an accepted point, including intermediate points before an explicit target.
This does not change carrier-row or global-closure acceptance.

`sweep.voltage_to_current` first solves the voltage points through the final
point, which must equal `switch_voltage_V`, and then replaces the voltage
target with each entry in `current_points_A_per_um`. The accepted equation is
`current_direction * Iterminal = Itarget`; tolerances and voltage-step limits
have the same scalar-control meaning as for the resistor mode. The two modes
are mutually exclusive.

When either mode is active, the main sweep CSV appends
`boundary_control_mode`, `inner_voltage_V`, `outer_voltage_V`,
`series_resistance_ohm_um`, `load_line_residual_V`,
`target_current_A_per_um`, `current_boundary_residual_A_per_um`, and
`boundary_control_evaluations`. `bias_V` remains the device/inner voltage so
existing BV extraction and field diagnostics keep their physical meaning.

`sweep.boundary_control` configures persistence and prediction shared by both
boundary modes. `evaluation_csv` records every attempted device solve and its
load-line/current residual. Each successful evaluation is written below
`checkpoint_directory`; `resume=true` reloads the best matching checkpoint for
the current target. `predictor_max_step_factor` limits guarded secant
extrapolation relative to the normal inner-voltage step.
`preferred_max_evaluations` is a performance target: exceeding it emits a
warning but does not weaken the residual acceptance test.
- `contact_current_reporting.endpoint_qf_floor.enabled`: opt-in reporting-only
  policy for Sentaurus restart parity. When enabled, DCSweep captures tiny
  contact-edge hole quasi-Fermi endpoint drops from the external
  `initial_state_file` before contact projection and passes those drops only
  to terminal/contact-current extraction. Optional `contacts` selects contacts.
  The nonlinear residual, Jacobian, solution state, VTK output, and ordinary
  continuation states are unchanged; continuation from a previously solved
  Vela point is not used as a QF-floor source.
- `diagnostics.contact_current_qf_floor.enabled`: compatibility alias for the
  same reporting-only policy. Prefer
  `contact_current_reporting.endpoint_qf_floor` in new BV configurations.

VTK node current and mobility fields are reconstructed post-processing
quantities. `NodeReconstructedElectronCurrentDensityVector`,
`NodeReconstructedHoleCurrentDensityVector`,
`NodeReconstructedTotalCurrentDensityVector`,
`NodeReconstructedElectronMobility`, and
`NodeReconstructedHoleMobility` are explicit provenance aliases; the older
unqualified names remain for compatibility. `ElectronIonIntegral`,
`HoleIonIntegral`, and `MeanIonIntegral` are local alpha-length accumulations,
not Sentaurus path ionization integrals. Prefer the equivalent
`LocalElectronAlphaLengthProxy`, `LocalHoleAlphaLengthProxy`, and
`LocalMeanAlphaLengthProxy` names in new comparisons. These local products are
formed after converting length and inverse length to compatible SI units; in
`unit_scaling` this is essential because coordinates use micrometers while
ionization coefficients use inverse centimeters. The separate
`ElectronPathIonIntegral`, `HolePathIonIntegral`, and `MeanPathIonIntegral`
VTK fields contain the ranked field-path result rather than the local proxy.

When the full VTK writer receives exact coupled-DD edge-flux diagnostics, it
also writes `DualFaceSgElectronCurrentDensityVector`,
`DualFaceSgHoleCurrentDensityVector`, and
`DualFaceSgTotalCurrentDensityVector` in A/cm^2. For each carrier, the signed
production SG particle line flux is divided by the edge's finite-volume
dual-face length and converted to conventional current density. Incident edge
normal components are then recovered at each node by a dual-face-length-
weighted least-squares fit. These node vectors are spatial-comparison
diagnostics; the original SG line flux remains the conservative transport
authority. The fields are omitted when exact edge diagnostics are unavailable.

For Sentaurus BV parity work, compare both the `.plt` terminal current and the
TDR-exported `ContactCurrentFlux` when judging the remaining terminal-current
gap. Near high reverse bias these two Sentaurus outputs can differ by
solver/output convention, so do not tune Vela mobility, SG flux, or avalanche
transport solely against the `.plt` value before checking the TDR contact flux.
The helper script
`scripts/verify_pn2d_sentaurus_terminal_current_crosscheck.py` performs this
cross-check and reports a JSON `sentaurus_plt_contact_flux_mismatch` status
when the two Sentaurus outputs exceed the configured relative tolerance.

For 2-D devices, currents and terminal charges are per-depth quantities by
default. Legacy CSV current values are per meter of device depth, and
`charge_C_per_m` / `capacitance_F_per_m` are also per meter. In
`unit_scaling` mode the CSV keeps those legacy column names but their numeric values are active internal units; it also appends per-micron
display columns (`*_A_per_um`, `charge_C_per_um`, `capacitance_F_per_um`) by
dividing per-meter values by `1e6`.

Step control fields:
- initial_step
- min_step
- max_step
- growth_factor
- step_growth_mode
- shrink_factor
- max_retries
- stop_on_failure

`step` and `initial_step` have distinct roles. For example, the following
configuration records nominal targets every `0.05 V`, starts continuation at
`1e-4 V`, grows successful internal steps by `1.2`, and never exceeds
`0.05 V`:

```json
{
  "step": -0.05,
  "initial_step": 1e-4,
  "min_step": 1e-10,
  "max_step": 0.05,
  "growth_factor": 1.2,
  "shrink_factor": 0.5
}
```

For explicit `bias_points`, `initial_step` is used only for the first interval.
Later intervals inherit the adaptive step magnitude from the preceding
accepted interval, capped by the next target distance and `max_step`.

`step_growth_mode` defaults to `fixed`, retaining the existing successful-step
multiplier. Opt in to `newton_iterations` to adjust growth using the number of
Newton updates in the last accepted voltage step, following the steady-state
strategy in the [Sentaurus training, section 6.3](https://ghzphy.github.io/Sentaurus_Training/sd/sd_6.html#sec3):

```text
G(N) = 1 + (growth_factor - 1) * max(0, 1 - max(N - 1, 0) / (0.75 * M))
next_step = clamp(abs(accepted_step) * G(N), min_step, max_step)
N = accepted step's Newton update count; M = solver.max_iter
```

`growth_factor` must be finite and at least `1`; it remains the upper bound on
growth. A zero-update acceptance uses the same multiplier as one update. Hard
accepted steps keep their step magnitude instead of growing; failed steps
still use the existing `shrink_factor`, rollback, and `max_retries` rules.
Acceptance involving QF-bounds or carrier-row recovery also suppresses growth,
because the final recovery solve's iteration count understates total work.
For example, `growth_factor: 1.35`, `solver.max_iter: 25`, and a `10 mV` step
requiring five Newton updates give a next proposal of about `12.7533 mV`.

This mode supports ordinary voltage sweeps with `solver.method: "newton"`
and a positive `solver.max_iter`, for both range and explicit `bias_points`
inputs. It also supports series resistors with `solver: "coupled_newton"`
as described above. Other solvers, nested-scalar/current-controlled sweeps, and
pseudo-arclength continuation are rejected explicitly. It does not change the
Newton stopping criteria, physical models, or predictor configuration.

Exact nominal/reference targets still clip steps, possibly below `min_step`.
In `newton_iterations` mode the next proposal grows from that actual clipped
step and is bounded again; the default `fixed` mode retains its previous
unclipped-proposal growth behavior. Runtime logs record the accepted bias,
Newton update count, effective multiplier, and next proposed magnitude on
`step_control: mode=newton_iterations` lines. Proposals may be clipped again
at the next target. These log lines include accepted internal steps; explicit
`bias_points` result/CSV rows still contain only the requested reference points.
External checkpoint drivers that launch separate sweeps
do not automatically inherit this controller's history.

Sweep continuation fields:
- `sweep.continuation.predictor.mode`: `none`, `constant`, `linear`, or `secant`.
- `sweep.continuation.predictor.fields`: optional subset of `psi`, `phin`, and
  `phip`; omitted non-`none` predictors default to all three fields.
- `sweep.continuation.predictor.max_extrapolation_ratio`: finite value at least `1`.
- `sweep.continuation.branch_acceptance.terminal_current_consistency`: when true,
  rejects a converged sweep point if terminal current cancellation is too large.
- `sweep.continuation.branch_acceptance.min_terminal_current_ratio`: non-negative
  threshold used by `terminal_current_consistency`.
- `sweep.continuation.branch_acceptance.psi_phin_jump`: when true, compares
  `psi - phin` against the previous accepted sweep state before accepting a
  new point. This is a diagnostic guard for reverse-bias branch jumps.
- `sweep.continuation.branch_acceptance.max_psi_phin_jump_V`: finite non-negative
  maximum allowed absolute nodewise jump in `psi - phin` when
  `psi_phin_jump` is enabled.
- `sweep.continuation.branch_acceptance.carrier_density_jump`: when true,
  compares electron density against the previous accepted sweep state in log10
  space before accepting a new point. This is a diagnostic guard for abrupt
  carrier-density branch jumps during reverse-bias continuation.
- `sweep.continuation.branch_acceptance.max_electron_density_jump_dex`: finite
  non-negative maximum allowed absolute nodewise electron-density jump in dex
  when `carrier_density_jump` is enabled.
- `sweep.continuation.branch_acceptance.max_electron_density_jump_p95_abs_dex`:
  optional finite non-negative maximum allowed p95 absolute electron-density
  jump in dex when `carrier_density_jump` is enabled. This is usually less
  sensitive to isolated node spikes than the nodewise maximum.

Pseudo-arclength continuation fields (`sweep.continuation.arclength`, default
disabled; the bias on `sweep.contact` is the continuation parameter):
- `sweep.continuation.arclength.enabled`: when `false` (default) the sweep keeps
  its voltage-parameterized stepping unchanged. Bounds below are only validated
  when this is `true`.
- `sweep.continuation.arclength.predictor`: must be `tangent` (the only currently
  supported predictor).
- `sweep.continuation.arclength.initial_step`: positive initial arclength step
  length; must lie within `[min_step, max_step]`.
- `sweep.continuation.arclength.min_step` / `max_step`: positive adaptive
  arclength step bounds; `min_step` must not exceed `max_step`.
- `sweep.continuation.arclength.growth_factor`: finite multiplier at least `1`
  applied after a clean corrector convergence.
- `sweep.continuation.arclength.shrink_factor`: multiplier in `(0, 1)` applied
  after a corrector failure.
- `sweep.continuation.arclength.max_corrector_iterations`: positive cap on the
  bordered-Newton corrector iterations per attempted step.
- `sweep.continuation.arclength.corrector_tolerance`: finite positive convergence
  tolerance on `max(||F||_inf, |arclength residual|)`.
- `sweep.continuation.arclength.max_step_retries`: non-negative number of shrink
  retries before a step is abandoned.
- `sweep.continuation.arclength.parameter_scale`: positive weight `theta` of the
  bias component in the arclength norm.
- `sweep.continuation.arclength.state_weight`: finite non-negative weight for
  state-space inner products. The default `0` selects the mesh-size-independent
  value `1 / state_dimension`, so the norm is
  `state_weight * ||x_dot||^2 + theta^2 * lambda_dot^2 = 1`.
- `sweep.continuation.arclength.damping_factor`: finite initial corrector
  line-search damping factor in `(0, 1]`.
- `sweep.continuation.arclength.max_line_search_steps`: non-negative number of
  corrector backtracking halvings allowed per bordered-Newton update.
- `sweep.continuation.arclength.bias_finite_difference_step_V`: positive finite
  voltage step used to estimate `dF/dV` for the bordered system.

When any continuation diagnostic is enabled, the sweep CSV appends
`predictor_mode`, `predicted_initial_state`, `branch_acceptance_status`,
`branch_acceptance_reason`, `terminal_current_consistency_ratio`, and
`psi_phin_max_jump_V`, plus carrier-density jump columns
`electron_density_jump_median_dex`, `electron_density_jump_p95_abs_dex`,
`electron_density_jump_max_abs_dex`, and `electron_density_jump_max_node`.

terminal_charge (for legacy single-terminal CV):
- terminal_charge.contact
- terminal_charge.regions
- terminal_charge.contact_radius
- terminal_charge.include_mobile_charge
- terminal_charge.include_ionized_dopants
- terminal_charge.per_meter
- terminal_charge.depth_m

terminal_charges (for multi-terminal quasi-static CV prototype):
- `terminal_charges` is an optional array of terminal-charge objects. When present,
  each entry is computed independently while a single sweep contact is varied.
- Each entry accepts `name`, `contact`, `regions`, `contact_radius`,
  `include_mobile_charge`, `include_ionized_dopants`, `per_meter`, and `depth_m`.
- `name` is sanitized to lowercase alphanumeric/underscore for CSV columns. If it
  is omitted, the contact name is used.
- The implementation is a quasi-static finite-difference prototype: for a sweep of
  contact `gate`, `capacitance_Cgate_drain_F_per_m` means `dQ_drain / dV_gate`. It is not
  an AC small-signal matrix solve or matrix inversion.
- CV CSV output always retains legacy `charge_C_per_m` / `capacitance_F_per_m`
  (or total-charge `charge_C` / `capacitance_F`) for compatibility, populated
  from the first configured terminal charge. With `terminal_charges`, additional
  columns are emitted as `charge_<name>_C_per_m` (or `_C`) and
  `capacitance_C<swept contact>_<name>_F_per_m`
  (or `_F`), for example `charge_gate_C_per_m`, `charge_drain_C_per_m`,
  `capacitance_Cgate_gate_F_per_m`, `capacitance_Cgate_drain_F_per_m`,
  `capacitance_Cgate_source_F_per_m`, and `capacitance_Cgate_body_F_per_m`
  for a gate sweep. Full sanitized names are used rather than initials so
  terminals such as `source` and `substrate` cannot collide.
- With `scaling.mode: "unit_scaling"` and `per_meter: true`, CV CSV output also
  appends `charge_C_per_um` and `capacitance_F_per_um` for the compatibility
  terminal charge columns.

stored_charge (optional IV/CV/BV mobile-charge proxy):
- `stored_charge` is an optional object under `sweep` for IV (`mode: "iv"`), quasi-static CV (`mode: "cv_quasistatic"`), or BV reverse (`mode: "bv_reverse"`) decks.
- Fields: `regions` (array of region names), `per_meter` (bool, default true),
  and `depth_m` (required > 0 when `per_meter` is false).
- When enabled, CSV adds `stored_charge_C_per_m` (or
  `stored_charge_C`) computed as a coarse proxy `q * integral(n + p) dV` over
  selected regions.
- This is a smoke-level stored-charge indicator, not a calibrated dynamic
  charge model.

Legacy aliases still accepted:
- charge_contact
- charge_regions
- charge_contact_radius
- charge_per_meter
- charge_depth_m


impact_ionization (optional BV source prototype):
- `impact_ionization` is an optional object under `sweep` for `mode: "bv_reverse"`.
- Fields:
  - `model`: `none` (default) or `selberherr`.
  - `a_n_per_m`, `b_n_V_per_m`: electron Selberherr coefficients.
  - `a_p_per_m`, `b_p_V_per_m`: hole Selberherr coefficients.
  - `min_field_V_per_m`: floor below which generation is forced to zero.
  - `max_generation_rate_per_m3s`: optional clamp for smoke-level stability.
- The current implementation is an engineering prototype used by BV trend
  decks (for example IGBT `simulation_bv_ii.json`); coefficients are not
  process-calibrated by default.

breakdown (for BV reverse):
- breakdown.max_electric_field_V_per_m
- breakdown.current_jump_ratio
- breakdown.non_convergence

Terminal-current stop criteria (every sweep mode; Sentaurus
`BreakCriteria{ Current(...) }`):
- breakdown.current_abs_A_per_um: stop once `|I| >= value` (positive number).
- breakdown.current_min_A_per_um: stop once the signed current `I <= value`.
- breakdown.current_max_A_per_um: stop once the signed current `I >= value`.
- breakdown.current_contact: contact whose total current is tested; defaults
  to `sweep.current_contact`.

`I` is the total current of that contact in A/um for a 2-D device. The
criterion is evaluated after each accepted point; the point that triggers it
is recorded (`breakdown_detected=1`, `criterion` one of `current_abs`,
`current_min`, `current_max`) and the sweep ends there. When any of these
criteria is configured the CSV gains `stop_criterion` (empty until the stop
point) and `stop_current_A_per_um` columns.

In `unit_scaling` mode, BV CSV output keeps `max_electric_field_V_per_m` as a compatibility column name using active internal electric-field values, and
also appends `max_electric_field_V_per_cm` using the active unit system conversion.

Legacy aliases still accepted:
- breakdown_max_electric_field_V_per_m
- breakdown_current_jump_ratio
- breakdown_on_non_convergence

## Mixed-material MOS DD deck example

A compact Si/SiO2 NMOS drift-diffusion prototype can use a semiconductor-only
source/body/drain set of ohmic contacts plus a `metal_gate` contact on the oxide
region. The continuity equations are intended to carry transport only in the Si
regions; SiO2 uses zero `ni`, `mun`, and `mup` from the built-in material
database so oxide carrier rows are pinned internally rather than treated as
transport unknowns.

Minimal contact and sweep fragments:

```json
{
  "mesh_file": "mesh.json",
  "doping": [
    { "region": "p_body", "donors": 0.0, "acceptors": 1e21 },
    { "region": "n_source", "donors": 5e21, "acceptors": 0.0 },
    { "region": "n_drain", "donors": 5e21, "acceptors": 0.0 },
    { "region": "gate_oxide", "donors": 0.0, "acceptors": 0.0 }
  ],
  "contacts": [
    { "name": "body", "type": "ohmic", "bias": 0.0 },
    { "name": "source", "type": "ohmic", "bias": 0.0 },
    { "name": "drain", "type": "ohmic", "bias": 0.1 },
    { "name": "gate", "type": "metal_gate", "bias": 0.1,
      "flatband_voltage": 0.0 }
  ],
  "sweep": {
    "mode": "cv_quasistatic",
    "contact": "gate",
    "start": 0.0,
    "stop": 0.1,
    "step": 0.05,
    "current_contact": "drain",
    "terminal_charge": {
      "contact": "gate",
      "regions": ["p_body", "n_source", "n_drain"],
      "per_meter": true,
      "contact_radius": 1e-6
    },
    "terminal_charges": [
      { "name": "gate", "contact": "gate", "regions": ["p_body", "n_source", "n_drain"],
        "per_meter": true, "contact_radius": 1e-6 },
      { "name": "drain", "contact": "drain", "regions": ["n_drain"],
        "per_meter": true, "contact_radius": 1e-6 },
      { "name": "source", "contact": "source", "regions": ["n_source"],
        "per_meter": true, "contact_radius": 1e-6 },
      { "name": "body", "contact": "body", "regions": ["p_body"],
        "per_meter": true, "contact_radius": 1e-6 }
    ]
  }
}
```

For off-state high-field diagnostics, set `sweep.mode` to `bv_reverse` and add
`breakdown.max_electric_field_V_per_m`, `breakdown.current_jump_ratio`, and
`breakdown.non_convergence` under the sweep block. The corresponding fields are
exercised by focused unit tests. This configuration family is an engineering
prototype and is not a calibrated MOSFET model.

## Experimental P1 convergence shadow capture

`VELA_P1_SHADOW_CAPTURE=1` enables an experimental observer in the classical
isothermal Newton progress output. It is off by default and does not select a
new convergence policy or change legacy acceptance. Existing progress output
must be configured; the observer adds `p1_shadow` to its JSON records.

The capture binds the actual advancement direction, processed update, damping,
constraint identities and residual to each saved iteration. Its compensated
linear diagnostics separately identify the pre-closure, post-constraint
matrix and the returned-coordinate solve equations after configured closure
replacement or continuity row scaling. Neither claims to capture the internal
eliminated/equilibrated PARDISO CSR or to certify forward accuracy.
The observer reuses existing row/closure/branch verdicts where available;
absent guard results remain `not_evaluated`. Physical residuals are reported
in C/m (free Poisson rows), A/µm (free continuity rows), and V (constraints).
Unit-valued reporting references are fixed and independent of the starting
residual; they are not calibrated device-accuracy tolerances.
`scripts/analysis/slot_shadow_convergence.py` replays adjacent archived states
using their exact binary64 reference/high/low/scale composition. New numeric
metrics are `not_calibrated`; missing inputs are `not_evaluated`. Neither is a
pass or a failure, and neither controls continuation.

This experiment does not implement `newton.convergence_policy`, calibrate
physical residual scales, or replace existing guards. See the
[P1 validation report](validation/slot_ldmos_p1_shadow_convergence_2026-10-20.md)
for coverage and missing evidence.

## Independent prescribed-source lattice heat diagnostic

`lattice_heat_probe INPUT.json OUTPUT.json` is a separate diagnostic executable.
It does not enable self-heating in `vela_example_runner` or DC sweeps. Build it
with the matching CMake preset and target `lattice_heat_probe`.

| Input field | Meaning |
|---|---|
| `mesh_file` | Existing Tri3 JSON mesh; use an absolute path for reproducibility. |
| `coordinate_to_metres` | Required positive conversion from the mesh's numeric coordinates to metres; e.g. `1e-6` for um. No electrical `scaling` setting is inferred. |
| `temperature_K` | One positive finite initial temperature per mesh node, in contiguous node-ID order. |
| `cell_source_W_per_m3` | One finite prescribed constant heat source per triangle, in contiguous cell-ID order. Negative cooling sources are permitted. |
| `region_conductivity` | One explicit law per region: `region_id`, `model: "constant"`, `value_W_per_m_K`; or `model: "inverse_quadratic"`, `numerator`, `denominator: [a,b,c]`, giving k(T)=numerator/(a+bT+cT²). |
| `thermodes` | Exterior-edge objects with `nodes: [i,j]`, `ambient_K`, `conductance_W_per_m2_K`. Duplicate and interior edges are rejected; unspecified outer boundaries insulate. |
| `solve_prescribed_source` | Default `false`: evaluate residual/Jacobian only. `true` applies safeguarded thermal Newton updates with the prescribed heat held fixed. |
| `absolute_residual_W_per_m` | Positive finite diagnostic infinity-norm tolerance, default `1e-8`; this is separate from device-level thermal acceptance gates. |

Temperature-dependent conductivity uses the element mean temperature and its
analytic derivative. Cell heat is deposited with area/3 weights. Surface heat
uses the consistent two-node Robin matrix. The output contains ordered node IDs,
temperatures, nodal areas, residuals in W/m, integrated source and outward thermal
boundary power in W/m, and thermal iteration history. The output path must not
already exist. `solved: false` is never a self-consistent thermal qualification.

`scripts/analyze_templates_ldmos_thermal.py` scores aligned reference and candidate
temperature JSON with an explicit contract. The approved LDMOS thermal contract
is in `reference_tcad/templates_ldmos_sentaurus2022/thermal/d0_thermal_acceptance.json`,
separate from the frozen generic `contracts/` bundle. A thermal-only pass does not
include electrical acceptance; zero power requires a separate absolute/equilibrium
check and is not reported as a relative heat-balance pass.


### Independent IALMob temperature diagnostic

`ialmob_probe INPUT.json` accepts a global `temperature_K` default and optional
`temperature_K` in each `states_SI` row (T >= 50 K). Existing six mobility partials
retain their order and units. With `derivatives=true`, the separate
`mobility_temperature_derivative_m2_per_Vs_K` holds densities, field and distance
fixed. Low-temperature phonon-exponent corrections are not supported.

An optional `high_field` object requires all four explicit keys:
`vsat300_m_per_s`, `beta300`, `vsat_exponent`, `beta_exponent`. Every state then
requires `drivingField_V_per_m`. This is the alpha=0, Vsat_Formula=1 local
primitive, with no automatic interface/contact driving-field selection.
`high_field_mobility_m2_per_Vs` is additional to the existing low-field output;
its temperature derivative includes the low-field mobility chain at fixed
carrier densities and driving field. A drive without a high-field configuration
is rejected. This diagnostic does not enable temperature feedback in DC sweeps.
See [native thermal mobility comparison](validation/templates_ldmos_d0_mobility_temperature_2026-09-12.md).

### Experimental four-equation electrothermal input

`electrothermal_probe INPUT.json OUTPUT.json` uses a separate explicit SI input;
it does not interpret production DC `scaling` or enable self-heating in the
production runner. `scripts/prepare_templates_ldmos_electrothermal.py` prepares
the audited LDMOS profile and input hashes. See the
[D0 report](validation/templates_ldmos_d0_electrothermal_2026-09-12.md) for its
model and qualification boundaries.

The mesh, coordinate conversion, conductivity laws and thermodes have the same
meaning as in the independent lattice-heat tool. Additional fields are:

| Input field | Meaning |
|---|---|
| `state_interleaved` | Contiguous per-node `[psi_V, fn_V, fp_V, temperature_K]`; T >= 50 K. |
| `donors_m3`, `acceptors_m3` | Per-node SI dopant concentrations. |
| `silicon_area_m2` | Silicon charge-control areas; zero marks non-silicon nodes. |
| `recombination_area_m2` | Continuity source-control areas; defaults to silicon areas. The qualified LDMOS adapter explicitly supplies global barycentric node areas at silicon nodes. |
| `fixed_charge_C_per_m` | Explicit per-node fixed charge. |
| `edge_geometry` | Complete edge list: `nodes`, `poisson_F_per_m`, nonnegative dimensionless `transport_weight`. |
| `mobility_SI` | Explicit constant or IALMob configuration in SI; the preparation adapter converts the audited TCAD saturation velocities. |
| `boundaries` | Objects with `node`, `kind` (`psi`, `fn`, `fp`, `temperature`, `neutral_contact`) and `value`. Neutral contacts constrain all three electrical variables using local charge neutrality and temperature. |
| `referenced_state_interleaved`, `electron_qf_reference_V`, `hole_qf_reference_V` | Optional complete restart pack. The fn/fp slots hold increments relative to the corresponding per-node references; preserve all three arrays together. |
| `use_qf_references` | Defaults to true. Without a restart pack, reference coordinates are prepared from the supplied absolute state. |
| `potential_origin_V` | Provenance offset for conversion back to the absolute external electrical gauge; does not itself alter the supplied state or boundaries. |
| `diagnostic_newton_max_iterations` | Default zero for frozen evaluation; a positive value enables bounded four-equation Newton updates. |
| `electrical_current_scale_A_per_m` | Explicit scale used to interpret existing numerical electrical gates. |
| `electrical_gate_solver` | Existing `carrier_row_convergence`, `continuity_row_scaling` and `block_absolute_convergence` configuration. Device qualification requires explicit original settings. |

The probe outputs physical and referenced states, original block/row gate
results, terminal currents, nodal temperatures/areas and integrated discrete
lattice source/boundary heat. A successful process exit can describe a failed
diagnostic solve; inspect the gate results. The continuation script additionally
checks KCL and heat balance before accepting a state. Full qualification uses
`scripts/analyze_templates_ldmos_d0.py` and both 31-point native curves and all
62 native temperature fields. No prescribed native heat source enters this
four-equation solve.

The explicit Newton-only diagnostic
`solver.region_resolved_interface_assembly.srh_signed_transport_volume_fraction`
defaults to 0 (unchanged SRH measure). It blends the inherited SRH area with
signed transport-material area: `V_SRH = V_old + f*(V_signed - V_old)`;
`f=1` uses the signed area directly. Fractions in [-1,1] permit independent
small positive/negative response calibration. Nonfinite values, negative
assembled areas, global continuity-volume flags, disabled SRH, and enabled
Auger are rejected. This option changes SRH residuals, Jacobian entries and
split-state source evaluation together; mesh volumes, Poisson charge, edge
transport, BTBT and avalanche measures are unchanged. It is not a global default.

### Explicit PhuMob element avalanche for SLOT IdVd

`solver.impact_ionization.current_approximation: "phumob_element_sg_best_vertex"`
selects an independent Si-only source profile. It uses the same 300 K PhuMob /
element Lombardi / full-derivative guards as the density-GradQF profile below;
it does not change that profile or the transport current discretization.
For each triangle, the vertex with the largest element box measure is selected
(ties use the smallest global node ID). Two incident SG edge drops reconstruct
one current vector. Mobility is evaluated at that vertex, and the constant
element generation rate is distributed with the existing vertex box measures.
The effective-band SG formula uses live positive populations and QF differences:
`a = log(c_j)-log(c_i) + s*(phi_j-phi_i)/Vt`,
`f = Vt*c_j*B(a)*expm1(s*(phi_j-phi_i)/Vt)`, with `s=+1` for
electrons and `s=-1` for holes; the physical particle-current drop is `-s*f`.
The statistical band correction is represented by the population/QF relation,
not by adding a generalized-Einstein multiplier to this source formula.
All live-population, mobility, field and source derivatives enter local AD.
Both the reconstructed source current and the avalanche GradQF force follow
the element partial-layer projection;
the existing explicit contact electric-field fallback has precedence.
Exactly zero current uses the existing zero subgradient of its Euclidean norm;
this is not a differentiability claim at zero. Modern Vela constants and the
definition-integral Fermi mapping are retained. Observed agreement with a native
kernel is a version/state-specific validation, not a general equivalence guarantee.
See `configs/regression/slot_ldmos_idvd_sg_source_control.json` for the opt-in
AverageBox configuration; it retains all device acceptance conditions.

`solver.impact_ionization.current_approximation: "phumob_element_density_grad_qf"`
is an opt-in, Si-only classical Newton profile. It requires `van_overstraeten`,
`generation: current_density`, `driving_force: quasi_fermi_gradient`,
`quasi_fermi_gradient_discretization: cell_gradient`,
`source_mapping_mode: element_vertex_box_measure`, `source_jacobian: local_ad`,
and the complete `phumob_field_lombardi` element-distance/HFS mobility profile
with `jacobian_field_derivatives: true`. It supports live Fermi or Boltzmann
populations at 300 K. RefDens interpolation, quantum corrections, IALMob boundary
projection, extra source/geometry factors and a constant-mobility override for
Si are rejected. PolySilicon and insulators contribute no avalanche source.

The source at vertex k of a Si triangle is
`M_k [alpha_n mu_n n_k |grad(phi_n)| + alpha_p mu_p p_k |grad(phi_p)|]`,
with explicit internal-to-SI factors. `M_k` is the existing cell box measure.
The vertex mobility includes PhuMob, Lombardi and the transport HFS field
stencil; the avalanche driving field is the configured cell quasi-Fermi
gradient (or the explicitly selected contact electric-field fallback).
The same local evaluator supplies residuals, source Jacobian, global source
audit and process-probe records. Local AD propagates the live Fermi density
response and analytic mobility partials; it does not differentiate through
finite differences of the complete source.

This profile uses `mesh_default` with `transport_edge_geometry: element_box`
and cell-material Poisson geometry. Imported external AverageBox edge couples
are not used by this element transport assembly; the external-profile avalanche
guard remains in force.

This is an explicit discretization choice, not a claim of numerical identity
with Sentaurus or completion of the SLOT 0--60 V validation. It does not change
the PN2D SG/GSS-Laux bundle, historical IALMob BVDS sources, split-state support,
or any convergence gate. The separate IdVd configurations are
  `configs/regression/slot_ldmos_idvd_fermi_{control,materials}.json`.

### Opt-in global population step

`solver.population_step_control` defaults to `{"mode":"legacy"}`. An explicit
`{"mode":"fraction_to_boundary","fraction":0.99}` selects a common positive
population step bound. The fraction must be finite and strictly between 0 and 1.
No concentration cutoff is used. Free electron and hole population tangents limit
one factor for the entire update; constrained and eliminated carrier slots do not
introduce artificial population bounds. A bordered system scales its contact
voltage direction by the same factor. Finite population changes are mapped back
through the configured Boltzmann/Fermi statistics, preserving split-coordinate
writeback where enabled. Subsequent line search and all existing acceptance,
finite-value, boundary, inverse-mapping and circuit protections remain active.

The SLOT `slot_population_memory` runner uses the same bound and records
`global_population_step` per update. This opt-in path replaces its componentwise
Poisson damping and nonpositive-population projection with a common step. Its
raw-direction U and existing stop conditions remain unchanged. This option is an
experimental numerical policy, not a validated replacement for legacy defaults.

### Opt-in population lower bound and predictor cutoff

For `solver.state_variables: potential_density`, the following controls are
independent. BVMethods population candidate floor defaults to 1e-100 cm^-3; legacy population mode defaults to no floor; predictor cutoff remains zero. Explicit zero floor restores the disabled-floor policy:

- `solver.population_update.minimum_density_cm3`: finite nonnegative lower
  bound applied after each free-carrier candidate update and its existing
  nonpositive protection. Convert cm^-3 to the configured internal concentration
  units at the solver callback. Contact prescribed slots and inactive oxide slots
  are not clipped. It does not change the original update RMS or residual gates.
- `solver.population_predictor_low_density_limit_cm3`: finite nonnegative
  cutoff for logarithmic population extrapolation. At a node whose most recent
  accepted concentration is at or below the cutoff, retain that concentration;
  do not raise it to the cutoff. Other nodes retain the existing log predictor.
  Potential and circuit-voltage prediction are unchanged.

These explicit experimental controls support separate Sentaurus compatibility
validation. T-2022.03 documents default `DensLowLimit=1e-100 cm^-3` for Newton
updates and `Extrapolate(LowDensityLimit=1e-100)` for extrapolation. The candidate floor does not establish native update equivalence; the predictor
cutoff remains independently disabled. The exact native mixed
history/equality boundary and boundary-slot behavior require separate evidence;
using these fields alone is not a claim of complete native equivalence.
Candidate floor and predictor-cutoff events are recorded in the runtime log;
these include rejected candidates and are not counts of accepted updates.

The standard deck canonicalizer converts the physical `*_cm3` keys above to
`minimum_density_m3` and `population_predictor_low_density_limit_m3`. These are parser spellings, not a numeric SI conversion: format v2 keeps
cm^-3 values and its field parser uses the supplied scaling context. Legacy
SI decks use m^-3 values (1 cm^-3 = 1e6 m^-3). Specifying both spellings for
the same control is rejected. Internal callback conversion remains separate.
