# Config Schema Reference

This document is the implementation-aligned reference for JSON config files used by Vela.
It describes fields currently parsed by the C++ code paths in Poisson and DC sweep drivers.

Scope and conventions:
- All numeric values are SI unless a field name includes a unit suffix such as `_eV`.
- Relative paths are resolved from the directory of the config JSON file.
- Legacy decks remain supported where noted.
- Prototype features are marked explicitly.

## Top-level fields

| Field | Type | Required | Notes |
| --- | --- | --- | --- |
| simulation_type | string | Optional | Common values: `poisson`, `dc_sweep`, `newton`. The selected CLI/tool path determines which fields are consumed. |
| mesh_file | string | Yes | Input mesh JSON path. |
| materials_file | string | Optional | Optional material override file. Supported shapes: top-level array, object with `materials` array, or object map keyed by material name. |
| output_csv | string | Optional | Default CSV output path for DC sweep. Can be overridden by `sweep.csv_file`. |
| output_vtk | string | Poisson: Yes | Poisson VTK output file path. |
| output_vtk_prefix | string | Optional | Default VTK prefix for DC sweep point outputs. Can be overridden by `sweep.vtk_prefix`. |
| doping | array | Yes | Region doping definitions; see below. |
| regions | array | Optional | Region-level fixed charge definitions; see below. |
| interfaces | array | Optional | Interface sheet/fixed/trap charge definitions; see below. |
| contacts | array | Yes | Contact bias and type definitions; see below. |
| boundaries | array | Optional | Explicit Poisson boundary segments (Neumann/insulating/symmetry); see below. |
| solver | object | Optional | Gummel/Newton settings for DD sweep and Newton solve. |
| sweep | object | Required for dc_sweep | Sweep mode, range, outputs, and diagnostics controls. |
| regression | object | Optional | Regression assertions consumed by `scripts/run_regression.py`. |

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

### regions[] entries

Common fields:
- name: string
- fixed_charge_m3: number (optional)

Notes:
- Region `fixed_charge_m3` is additive source charge for Poisson and for drift-diffusion Poisson substeps.
- Duplicate definitions for the same region are rejected.

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

About node_ids or edge_node_pairs:
- `node_ids` polyline is the implemented path.
- `edge_node_pairs` is not currently parsed by the C++ drivers and should be treated as reserved for future schema extensions.

About normal_displacement_C_per_m2 / normal_electric_field_V_per_m:
- `normal_displacement_C_per_m2` is implemented for Neumann Poisson boundaries.
- `normal_electric_field_V_per_m` is not currently parsed; treat it as a future extension placeholder.

## solver

The solver object is used by DD sweep and Newton solve paths.

Method selection:
- method: `gummel` or `newton`
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
- taun
- taup
- bandgap_narrowing

Newton-specific keys:
- damping_factor
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

Notes:
- `line_search` and `damping_factor` apply to Newton config.
- Both Gummel/Newton parse `mobility`, `recombination`, `impact_ionization`, `temperature_K`.

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

Legacy aliases still accepted:
- charge_contact
- charge_regions
- charge_contact_radius
- charge_per_meter
- charge_depth_m

breakdown (for BV reverse):
- breakdown.max_electric_field_V_per_m
- breakdown.current_jump_ratio
- breakdown.non_convergence

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
  including `drain_current_sign` and `current_monotone_abs_tolerance`

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
