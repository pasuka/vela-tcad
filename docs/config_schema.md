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
| simulation_type | string | Optional | Common values: `poisson`, `dc_sweep`, `newton`. The selected CLI/tool path determines which fields are consumed. |
| mesh_file | string | Yes | Input mesh JSON path. |
| materials_file | string | Optional | Optional material override file. Supported shapes: top-level array, object with `materials` array, or object map keyed by material name. |
| output_csv | string | Optional | Default CSV output path for DC sweep. Can be overridden by `sweep.csv_file`. |
| output_vtk | string | Poisson: Yes | Poisson VTK output file path. |
| output_vtk_prefix | string | Optional | Default VTK prefix for DC sweep point outputs. Can be overridden by `sweep.vtk_prefix`. |
| scaling | object | Optional | Input unit interpretation. Omit for legacy SI behavior, or set `mode` to `unit_scaling`; see below. |
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

For `dc_sweep`, omitting `solver.method` keeps the default Gummel path. Set
`solver.method` (or legacy alias `solver.type`) to `newton` to use the coupled
Newton sweep path where supported.

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
| electric field | V/m | V/cm |
| sheet density | m^-2 | cm^-2 |
| voltage | V | V |
| temperature | K | K |
| energy | eV | eV |

`scaling.mode` does not accept `si` or vendor-specific aliases.
Field names with explicit unit suffixes remain the stable public names; this
mode defines how their numeric values are read at the schema boundary.
The conversion is applied while reading mesh coordinates, material override
files, doping, region fixed charge, interface sheet/trap charge, mobility
settings, and electric-field-related solver settings. The Poisson and
drift-diffusion assemblers continue to receive SI values.

Poisson driver note:
- In legacy mode (no `scaling`), Poisson uses the historical SI assembly path.
- In `unit_scaling` mode, Poisson uses a scaled (dimensionless) assembly path
  and restores physical potential before producing outputs.
- VTK `potential_V` remains the physical potential in volts in both modes.

CSV output keeps the legacy SI columns for compatibility. In `unit_scaling`
mode, DC sweep CSV output also appends convenience columns using common
external TCAD display units, such as `current_total_A_per_um`,
`current_electron_A_per_um`, `current_hole_A_per_um`, `charge_C_per_um`,
`capacitance_F_per_um`, and `max_electric_field_V_per_cm`.

DC sweep CSVs include solver provenance columns:
- `solver_method`: selected nonlinear path, such as `gummel`, `newton`, or `gummel_newton`
- `gummel_iterations`: iterations used by the Gummel stage for this bias point
- `newton_iterations`: iterations used by the coupled Newton stage for this bias point
- `handoff_stage`: final accepted stage or failure stage, such as `newton`, `gummel_failed`, `newton_failed`, or `gummel_fallback`

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
  `fixed_charge_m3` numeric inputs are read as `cm^-3` and normalized to
  `m^-3` internally.

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
`cm^-3`, normalized to `m^-3` internally. In legacy SI mode, the same values are
read through the legacy concentration path.

### regions[] entries

Common fields:
- name: string
- fixed_charge_m3: number (optional)

Notes:
- Region `fixed_charge_m3` is additive source charge for Poisson and for drift-diffusion Poisson substeps.
- Duplicate definitions for the same region are rejected.
- With `scaling.mode: "unit_scaling"`, `fixed_charge_m3` numeric inputs are
  read as `cm^-3` and normalized to `m^-3` internally.

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
  `cm^-2` and normalized to `m^-2` internally.

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
  (`C/cm^2`) and normalizes to SI before assembly.

About node_ids or edge_node_pairs:
- `node_ids` polyline is the implemented path.
- `edge_node_pairs` is not currently parsed by the C++ drivers and should be treated as reserved for future schema extensions.

About normal_displacement_C_per_m2 / normal_electric_field_V_per_m:
- `normal_displacement_C_per_m2` is implemented for Neumann Poisson boundaries.
- `normal_electric_field_V_per_m` is not currently parsed; treat it as a future extension placeholder.

## solver

The solver object is used by DD sweep and Newton solve paths.

Method selection:
- method: `gummel`, `newton`, or `gummel_newton`
- type: alias for method (legacy compatibility)

Commonly used controls:
- max_iter
- reltol
- abstol
- temperature_K
- mobility
- recombination
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
- line_search
- warm_start
- verbose
- diagnostics / diagnostic_history
- jacobian (`analytic` or `finite_difference`)
- finite_difference_step
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
- `max_update` is an optional non-negative infinity-norm cap on one Newton
  update in solver unknown units; `0` disables the cap.
- `carrier_floor_m3` is an optional non-negative Gummel carrier floor used to
  keep reconstructed quasi-Fermi potentials consistent with solved carrier
  densities.
- Both Gummel/Newton parse `mobility`, `recombination`, `impact_ionization`, `temperature_K`.
- With `scaling.mode: "unit_scaling"`, `bandgap_narrowing.reference_doping_m3`
  is read as `cm^-3` and normalized to `m^-3`.

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
`candidate_column`, `max_orders_of_magnitude`, `max_relative_error`,
`min_points`, and `require_trend_match` to the comparison report.
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
- `dominant_signed_region`: when a global node receives equal donor and acceptor active concentrations, use the signed `DopingConcentration` field with the largest magnitude to choose a single majority dopant for the node, and record the rewrite in `doping_metadata.json`.

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
  "smoothing": 0.5
}
```

Supported `model` values:
- `none`
- `slotboom`

The Slotboom prototype computes an effective bandgap narrowing from the maximum
of absolute net doping and local carrier densities, then feeds the resulting
effective intrinsic density into the drift-diffusion statistics path. This is
implemented in Gummel and Newton configurations. With
`scaling.mode: "unit_scaling"`, `reference_doping_m3` numeric input is read as
`cm^-3` and normalized to `m^-3`.

### Newton diagnostics and residual options

Newton configs can opt into diagnostic history with either
`"diagnostics": true` or `"diagnostic_history": true`. The solver also accepts:

- `jacobian`: `analytic` or `finite_difference`
- `finite_difference_step`
- `residual_norm`: `block` or `l2`
- `residual_weights`: object with `psi`, `phin`, and `phip`
- `residual_scales`: object with `psi`, `phin`, and `phip`

`warm_start: true` preserves supplied quasi-Fermi potentials when continuing
from a previous solution. The default `false` keeps the conservative
cold-start behavior.

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

Field meanings (Selberherr prototype):
- `electron_A_m_inv` (1/m): electron ionization prefactor.
- `electron_B_V_m` (V/m): electron critical field.
- `hole_A_m_inv` (1/m): hole ionization prefactor.
- `hole_B_V_m` (V/m): hole critical field.
- `carrier_velocity_m_s` (m/s): effective saturated carrier speed used by the
  generation-rate proxy.

Validation:
- `electron_A_m_inv`, `hole_A_m_inv`, and `carrier_velocity_m_s` must be non-negative.
- `electron_B_V_m` and `hole_B_V_m` must be positive.

Scaling:
- With `scaling.mode: "unit_scaling"`, `electron_A_m_inv` and
  `hole_A_m_inv` numeric inputs are read as `cm^-1`, while
  `electron_B_V_m` and `hole_B_V_m` are read as `V/cm`. They are normalized
  to `1/m` and `V/m` before the impact-ionization model sees them.
- `carrier_velocity_m_s` remains `m/s`.

Prototype note:
- This is an engineering impact-ionization source term for smoke diagnostics; it
  is not a calibrated avalanche-breakdown prediction model.

### mobility

`solver.mobility` accepts either the legacy string form or an object. String decks remain compatible:

```json
"mobility": "caughey_thomas"
```

Object form:

```json
"mobility": {
  "model": "caughey_thomas_field_surface",
  "electron_mu_min_m2_V_s": 0.00522,
  "electron_nref_m3": 9.68e22,
  "electron_alpha": 0.68,
  "hole_mu_min_m2_V_s": 0.00449,
  "hole_nref_m3": 2.23e23,
  "hole_alpha": 0.70,
  "electron_saturation_velocity_m_s": 1.0e5,
  "electron_field_beta": 2.0,
  "hole_saturation_velocity_m_s": 1.0e5,
  "hole_field_beta": 2.0,
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

Supported `model` values are `constant`, `caughey_thomas`, `caughey_thomas_field`, `caughey_thomas_surface`, and `caughey_thomas_field_surface`. The `surface` block is a MOS prototype for Si/SiO2-like channel mobility degradation, not a calibrated Lombardi model. It applies a vertical-field factor `mu_eff = mu_bulk / (1 + (theta * max(|E_normal| - reference_field, 0))^beta)^(1/beta)`, optionally clamped by `min_factor`/`max_factor`.

The first implementation estimates `E_normal` with the local edge electric-field magnitude on edges that match `surface_region` and/or the two-name `surface_interface`; this is sufficient for trend regressions but should not be interpreted as a calibrated normal-field extraction. If no matching surface edge is found for a mobility evaluation, surface degradation is disabled and the existing low-field or velocity-saturation behavior is used.

With `scaling.mode: "unit_scaling"`, Caughey-Thomas mobility floors are read
as `cm^2/(V s)`, reference dopings as `cm^-3`, surface reference fields as
`V/cm`, and surface theta coefficients as `cm/V`. They are normalized to
`m^2/(V s)`, `m^-3`, `V/m`, and `m/V` before mobility evaluation.

## sweep

Required core fields:
- mode: `iv`, `cv_quasistatic`, or `bv_reverse` (aliases: `cv`, `bv`, `reverse_breakdown`)
- contact: swept contact name
- start: number
- stop: number
- step: non-zero number, sign must move start toward stop

Output and current fields:
- current_contact
- write_vtk
- vtk_prefix
- csv_file

For 2-D devices, currents and terminal charges are per-depth quantities by
default. Legacy CSV current values are per meter of device depth, and
`charge_C_per_m` / `capacitance_F_per_m` are also per meter. In
`unit_scaling` mode the CSV keeps those legacy columns and appends per-micron
display columns (`*_A_per_um`, `charge_C_per_um`, `capacitance_F_per_um`) by
dividing per-meter values by `1e6`.

Step control fields:
- min_step
- max_step
- growth_factor
- shrink_factor
- max_retries
- stop_on_failure

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

In `unit_scaling` mode, BV CSV output keeps `max_electric_field_V_per_m` and
also appends `max_electric_field_V_per_cm` by dividing the SI value by `100`.

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
