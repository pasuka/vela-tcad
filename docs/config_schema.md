# Config Schema Reference

This document is the implementation-aligned reference for JSON config files used
by Vela. It describes fields currently parsed by the C++ Poisson, DC sweep, and
single-bias Newton paths.

Use this file as the field-level reference. Use
[examples.md](examples.md) for deck support status and
[architecture.md](architecture.md) for solver path boundaries.

Scope and conventions:
- No `scaling` field keeps the legacy SI input behavior used by existing decks.
- In legacy SI mode, all numeric values are SI unless a field name includes a unit suffix such as `_eV`.
- The optional public input-unit mode is `scaling.mode: "unit_scaling"`; it is a unit interpretation and numeric scaling foundation only, not a calibration feature.
- Relative paths are resolved from the directory of the config JSON file.
- Legacy decks remain supported where noted.
- Prototype features are marked explicitly.
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
| output_state_file | string | Optional | Restart-state CSV written by `newton_solve_from_state`. Uses Vela restart format: `node_id,psi,phin,phip,electrons_m3,holes_m3`. |
| state_fields_dir | string | Required for `newton_solve_from_state` and probe-style external-state tools | Directory containing Sentaurus-export-style scalar field CSVs. |
| runtime_log | object | Optional | Runtime log control. Default is enabled and writes `<config stem>.log` next to the config file. |
| scaling | object | Optional | Input unit interpretation. Omit for legacy SI behavior, or set `mode` to `unit_scaling`; see below. |
| mesh_geometry | object | Optional | Mesh box-geometry options; see below. |
| doping | array | Yes | Region doping definitions; see below. |
| regions | array | Optional | Region-level fixed charge definitions; see below. |
| interfaces | array | Optional | Interface sheet/fixed/trap charge definitions; see below. |
| contacts | array | Yes | Contact bias and type definitions; see below. |
| boundaries | array | Optional | Explicit Poisson boundary segments (Neumann/insulating/symmetry); see below. |
| solver | object | Optional | Gummel/Newton settings for DD sweep and Newton solve. |
| sweep | object | Required for dc_sweep | Sweep mode, range, outputs, and diagnostics controls. |
| regression | object | Optional | Regression assertions consumed by `scripts/run_regression.py`. |

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
  "node_volume_policy": "mixed_voronoi"
}
```

Supported `node_volume_policy` values:
- `barycentric` (default): each Tri3 cell contributes one third of its area to
  each adjacent node.
- `mixed_voronoi`: acute and right triangles use the cotangent Voronoi vertex
  area; obtuse triangles assign half the area to the obtuse vertex and one
  quarter to each other vertex. Edge couplings still use the existing
  cotangent box coefficients and negative-cotangent fallback policy.

This option changes the mesh node control volumes used by Poisson charge,
carrier continuity volume terms, stored/terminal charge, and SG avalanche
source density normalization. It does not change donor/acceptor values used for
Ohmic contact boundary conditions.

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
- DC sweep (Gummel): `ohmic` and prototype `schottky` supported.
- DC sweep (Newton): Schottky currently rejected with a clear error.

Prototype notes:
- Schottky is an engineering prototype (Dirichlet-barrier style), not a calibrated thermionic-emission model.

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
  all supplied state fields, and reports zero terminal current.

Commonly used controls:
- max_iter
- reltol
- abstol
- temperature_K
- mobility
- recombination
- auger_cn_m6_per_s
- auger_cp_m6_per_s
- impact_ionization

Gummel-specific keys:
- damping_psi
- carrier_floor_m3
- taun
- taup
- bandgap_narrowing

Newton-specific keys:
- damping_factor
- max_update
- quasi_fermi_update_limit_V
- quasi_fermi_update_limit_minority_V
- electron_quantum_potential
- carrier_regularization_scale
- line_search
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
- quasi_fermi_reference (`none` or `contact_majority`)
- carrier_row_convergence
- continuity_row_scaling
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
- `electron_quantum_potential` enables the electron density-gradient quantum
  correction.  The object accepts `coupling_mode` (`outer`, the default, or
  `frozen` for an imported restart potential), `formulation`
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
  not a mobility calibration knob.
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
- `poisson_line_search_stall_contact_majority_qf_drop_limit_V` is an optional
  non-negative physical-voltage guard for the `poisson_line_search_stall_floor`
  acceptance path. When greater than zero, a Poisson-block line-search stall is
  accepted only if the maximum majority-carrier quasi-Fermi drop across contact-to-interior edges
  is below this value. The default is `5e-11` V; `0` disables this guard.
- `carrier_regularization_scale` is an experimental non-negative Newton
  stabilization knob. When greater than zero, Vela adds
  `sign(diagonal) * scale * carrier_row_abs_sum` to each carrier continuity
  diagonal before solving the coupled Newton step; the row sum includes the
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

Supported values are `boltzmann` (the compatibility default) and
`fermi_dirac`. The Fermi-Dirac selection is one atomic physics path: carrier
densities use the normalized complete integral `F_{1/2}`, ideal Ohmic contacts
solve charge neutrality with the same statistics, continuity edges use the
generalized Einstein/Scharfetter-Gummel operator, and terminal-current
post-processing uses that identical edge operator. The density-form Gummel
solver currently rejects `fermi_dirac` explicitly; use `method: "newton"`.

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

The Slotboom and OldSlotboom prototypes compute the positive effective
bandgap-narrowing term from the maximum of absolute net doping and local
carrier densities, then feed the resulting effective intrinsic density into
the drift-diffusion statistics path. For Sentaurus `OldSlotboom` parity, the
`models.par` `Bandgap.dEg0(OldSlotboom) = -1.595e-2 eV` term is handled by the
material intrinsic-density override, while the `old_slotboom` BGN term uses
`Ebgn = 9e-3 eV`, `Nref = 1e17 cm^-3`, and `C = 0.5`. This is implemented in
Gummel and Newton configurations. With
`scaling.mode: "unit_scaling"`, `reference_doping_m3` numeric input is read as
`cm^-3` and kept internally as `cm^-3`. The compiled `Nref` default of both
`slotboom` and `old_slotboom` is likewise expressed in the internal unit
system (`1e23 m^-3` in `legacy_si`, `1e17 cm^-3` in `unit_scaling`).

Set `fermi_statistics_correction: true` to reproduce the additional bandgap-
narrowing correction that Sentaurus applies by default when `Fermi` and
`EffectiveIntrinsicDensity(OldSlotboom)` are active together.  The correction
is evaluated from donor/acceptor concentrations and the 300 K density of states.
Leave it `false` to reproduce Sentaurus `EffectiveIntrinsicDensity(NoFermi)` or
the historical Vela behavior.

### Newton diagnostics and residual options

Newton configs can opt into diagnostic history with either
`"diagnostics": true` or `"diagnostic_history": true`. The solver also accepts:

- `jacobian`: `analytic` or `finite_difference`
- `finite_difference_step`
- `quasi_fermi_reference`: `none` (default) or `contact_majority`
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

The source-aware convergence controls are:

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
"global_continuity_closure": {
  "mode": "enforce",
  "tolerance": 0.01,
  "source_floor": 1e-14
}
```

`continuity_row_scaling` left-scales the Newton continuity rows using the
larger of their configured flux/source measures. It changes conditioning only,
not the nonlinear equations. `min_source_scale` supplies an absolute scale for
low-current problems where a fraction of the largest source would suppress the
small rows of interest.

`global_continuity_closure` independently sums the free-node SRH/avalanche
source and all contact fluxes for electrons and holes. `report` appends the
metrics without rejecting a solution; `enforce` prevents convergence unless
each carrier with `|integrated source| >= source_floor` has relative mismatch
no greater than `tolerance`. The denominator is the largest of the contact
flux magnitude, integrated-source magnitude, and `source_floor`.

`warm_start: true` preserves supplied quasi-Fermi potentials when continuing
from a previous solution. The default `false` keeps the conservative
cold-start behavior.

### band_to_band

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

`solver.mobility` accepts either the legacy string form or an object. String decks remain compatible:

```json
"mobility": "caughey_thomas"
```

Object form:

```json
"mobility": {
  "model": "caughey_thomas_field_surface",
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

Supported `model` values are `constant`, `caughey_thomas`,
`caughey_thomas_field`, `caughey_thomas_surface`,
`caughey_thomas_field_surface`, `masetti`, `masetti_field`,
`masetti_surface`, and `masetti_field_surface`.
For field-saturation models, `high_field_driving_force` is `electric_field`
by default and may be set to `quasi_fermi_gradient` to match Sentaurus
`HighFieldSaturation`; electrons use `|grad(phin)|` and holes use
`|grad(phip)|`.
When quasi-Fermi driving is selected, `high_field_gradient_discretization`
controls how that magnitude is recovered. `edge_projection` (default) uses the
potential difference projected onto each edge and preserves legacy Vela
behavior. `transport_cell_vector` reconstructs the two-dimensional P1 gradient
inside adjacent semiconductor cells and area-averages its magnitude onto each
edge; this matches the orientation-independent `GradQuasiFermi` semantics used
by Sentaurus on unstructured 2-D meshes.

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
`surface.surface_interface` should identify the semiconductor/insulator pair;
the normal field is the electric-field projection along that interface normal.

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
- bias_points: optional array of explicit biases. When present, the adaptive
  `start`/`stop`/`step` stepping loop is bypassed and Vela solves exactly these
  biases in order. Use a one-element array for a single-bias restart.

Output and current fields:
- current_contact
- write_vtk
- vtk_prefix
- csv_file
- initial_state_file: optional restart-state CSV used as the initial
  `DDSolution` for the first solved bias.
- write_state_file: optional restart-state CSV overwritten after every
  converged point with the latest `DDSolution`.
- `initialization.mode`: optional first-point initialization mode. `none`
  preserves the baseline cold-start path; `poisson_block` runs one Newton
  Poisson block solve before the first coupled Newton solve and uses that state
  as the first-point handoff.
- `initialization.diagnostic_csv`: optional CSV written only when
  `initialization.mode` is `poisson_block`. It records the first-point
  initialization bias as `bias_V`, plus the cold-state and Poisson-block
  residual norms.
- `initialization.write_state_file`: optional restart-state CSV written only
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

- `configs/templates/pn2d_iv.template.json`: forward 0--20 V IV, low-field
  Masetti mobility, `cell_reconstructed_total_impurity`, SRH plus Old
  Slotboom, and impact ionization disabled.
- `configs/templates/pn2d_bv.template.json`: reverse 0---20 V BV,
  `masetti_field`, `net_doping`, SRH plus Old Slotboom, and
  Van Overstraeten impact ionization.

The PN2D BV version 3 default is the qualified atomic
`element_edge_sg_gss_laux` profile. Rendering it binds element-edge
SG/GSS-Laux current support, element-vertex box source mapping, Bernoulli
midpoint density, mixed-Voronoi node volumes, and
`mesh_geometry.require_non_obtuse=true`. Use the
`legacy_cell_reconstructed` profile for the complete cell-reconstructed and
barycentric rollback. The generator rejects partially mixed profiles. This
template policy does not change the global barycentric parser default or the
PN2D IV template.

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
  --set initial_state_file="inputs/bv_prebias_5p9V_state.csv"
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

Rows must cover every mesh node exactly once by `node_id`. All state values are
stored in active internal solver units: potentials in V and carrier densities in legacy `m^-3` or `unit_scaling` `cm^-3`.
For checkpoint-style BV runs, set `write_state_file` during the first run, then
resume from the latest converged point with `initial_state_file` plus either a
new adaptive range or `bias_points` for selected target biases. Sweep CSV files
are still opened as fresh outputs; use distinct CSV names for resumed segments.

Diagnostics fields:
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

### External circuit and voltage-to-current boundary control

`sweep.external_circuit` enables a Sentaurus-compatible 2-D series-resistor
load line. The supported mode is `series_resistor`. `resistance_ohm_um` uses
`ohm*um`, terminal current uses `A/um`, and the closed equation is

```text
OuterVoltage = InnerVoltage + current_direction * Iterminal * Rseries
```

`bias_points` (or `start/stop/step`) are outer/source voltages. The nonlinear
device boundary remains the inner voltage and is solved repeatedly until
`load_line_residual_V` meets `residual_tolerance_V`. `current_direction` is
`+1` for the BVmethods NMOS drain convention and may be `-1` for an oppositely
oriented terminal. `initial_inner_voltage_V`, `max_inner_voltage_step_V`,
`voltage_tolerance_V`, `max_bracket_steps`, and `max_iterations` control only
the scalar boundary solve.

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
`breakdown.non_convergence` under the sweep block. See
`examples/nmos2d_mos_dd/simulation_bv.json` for the CI smoke deck. This example
family is an engineering prototype and is not a calibrated MOSFET model.

## regression block

Regression fields are optional and consumed by the regression runner.

### regression.dc_sweep

Supported fields include:
- expected_rows
- max_abs_attempted_step
- max_abs_accepted_step
- max_retry_count
- require_monotone_abs_current
- require_monotone_max_field
- min_converged_rows

Also supported:
- allow_nonconverged_final_bv_point
- current_monotone_abs_tolerance
- current_monotone_rel_tolerance
- max_field_monotone_abs_tolerance
- max_field_monotone_rel_tolerance
- min_max_electric_field_V_per_m
- max_max_electric_field_V_per_m
- allow_zero_capacitance
- expected_zero_capacitance_rows
- min_nonzero_capacitance_rows

### regression (top-level)

Common fields used by examples:
- declared_converged
- dc_sweep: { ... }
- example-specific keys used by dedicated checks (for example MOS interface probes)
- ldmos_iv: optional regression-runner settings for the LDMOS DD-IV smoke check,
  including `drain_current_sign`, `current_monotone_abs_tolerance`, and
  `current_monotone_rel_tolerance`
- mos: optional Id-Vd / generated Id-Vg trend settings for MOS examples,
  including `device`, `drain_current_sign`, and nested `idvg` sweep controls.
- surface_mobility: optional comparison block for a surface-mobility variant
  against a baseline Id-Vg deck. Fields include `baseline_config`,
  `baseline_csv`, and `current_ratio_tolerance`.
- schottky_iv: optional Schottky IV trend block with current sign and monotonic
  tolerance fields.
- ldmos_fieldplate_trend: optional LDMOS field-plate comparison block with
  `baseline_config`, optional baseline/variant field columns, and
  `max_field_ratio_limit`.
- igbt_high_injection: optional high-injection IV comparison block with
  baseline CSV/config fields and stored-charge monotonicity settings.
- igbt_charge_cv: optional stored-charge and multi-terminal CV trend checks.
- igbt_bv: optional BV/impact-ionization comparison block with baseline config,
  bias-match tolerance, and current multiplier tolerance.

## Minimal examples

Poisson with explicit boundary/contact types:

```json
{
  "simulation_type": "poisson",
  "mesh_file": "mesh.json",
  "output_vtk": "outputs/result.vtk",
  "doping": [
    { "region": "silicon", "donors": 1e21, "acceptors": 0.0 }
  ],
  "contacts": [
    { "name": "anode", "type": "ohmic", "bias": 0.0 },
    { "name": "gate", "type": "metal_gate", "bias": 1.0 }
  ],
  "boundaries": [
    { "name": "left", "type": "symmetry", "node_ids": [0, 3, 6] },
    { "name": "right", "type": "insulating", "node_ids": [2, 5, 8] }
  ]
}
```

Newton solve initialized from external scalar fields:

```json
{
  "simulation_type": "newton_solve_from_state",
  "mesh_file": "mesh.json",
  "state_fields_dir": "path/to/state_fields",
  "output_state_file": "outputs/minus20_from_state.csv",
  "output_vtk": "outputs/minus20_from_state.vtk",
  "doping": [
    { "region": "n_region", "donors": 1e23, "acceptors": 0.0 },
    { "region": "p_region", "donors": 0.0, "acceptors": 1e23 }
  ],
  "contacts": [
    { "name": "anode", "type": "ohmic", "bias": -20.0 },
    { "name": "cathode", "type": "ohmic", "bias": 0.0 }
  ],
  "solver": {
    "method": "newton",
    "max_iter": 40,
    "warm_start": true,
    "line_search": true
  }
}
```

DC sweep with Gummel:

```json
{
  "simulation_type": "dc_sweep",
  "mesh_file": "mesh.json",
  "output_csv": "outputs/iv.csv",
  "doping": [
    { "region": "n_region", "donors": 1e23, "acceptors": 0.0 },
    { "region": "p_region", "donors": 0.0, "acceptors": 1e23 }
  ],
  "contacts": [
    { "name": "anode", "type": "ohmic", "bias": 0.0 },
    { "name": "cathode", "type": "ohmic", "bias": 0.0 }
  ],
  "solver": {
    "method": "gummel",
    "max_iter": 80,
    "reltol": 1e-5,
    "damping_psi": 0.5,
    "temperature_K": 300.0
  },
  "sweep": {
    "mode": "iv",
    "contact": "anode",
    "start": 0.0,
    "stop": 0.5,
    "step": 0.25,
    "current_contact": "anode",
    "write_vtk": true,
    "vtk_prefix": "outputs/iv"
  }
}
```


## power-device regression block examples

LDMOS/IGBT decks can combine `regression.dc_sweep` with device-specific trend
blocks while staying explicitly prototype-level:

```json
"regression": {
  "dc_sweep": {
    "expected_rows": 7,
    "min_converged_rows": 6,
    "require_monotone_abs_current": true,
    "require_monotone_max_field": true
  },
  "ldmos_fieldplate_trend": {
    "baseline_config": "simulation_bv.json",
    "max_field_ratio_limit": 1.20
  },
  "igbt_high_injection": {
    "baseline_config": "simulation_iv.json",
    "baseline_csv": "outputs/igbt2d_iv_baseline.csv",
    "baseline_final_current_min_ratio": 1.0,
    "require_stored_charge_monotone": true,
    "stored_charge_monotone_direction": "either",
    "stored_charge_monotone_abs_tolerance": 1e-24,
    "stored_charge_monotone_rel_tolerance": 1e-8
  },
  "igbt_charge_cv": {
    "require_stored_charge_monotone": true,
    "stored_charge_monotone_direction": "either",
    "stored_charge_monotone_abs_tolerance": 1e-24,
    "stored_charge_monotone_rel_tolerance": 1e-8
  }
}
```

These checks are trend validation guards (finite outputs + directional checks),
not calibrated silicon sign-off criteria. `stored_charge_monotone_direction` accepts
`"nondecreasing"`, `"nonincreasing"`, or `"either"`.
