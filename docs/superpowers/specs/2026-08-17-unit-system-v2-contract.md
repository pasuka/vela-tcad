# Deck `format_version: 2` Unit and Format Contract

Status: normative contract for the legacy-SI removal migration.
Supersedes the unit-mode portions of
`2026-07-08-tcad-internal-unit-system-design.md`.

## 1. Purpose

Vela currently ships two internal physical unit systems selected by
`scaling.mode`:

- `legacy_si` (default when `scaling` is omitted): internal units are SI.
- `unit_scaling`: internal units are TCAD units.

The goal of the migration is a **single** internal unit system (TCAD) and the
complete removal of `UnitScalingMode`, `PhysicalUnitSystem::legacySI()`,
`UnitScalingConfig::isUnitScaling()`, and every `scaling.enabled` branch.

Two properties of the current state make a naive removal unsafe, and this
document exists to close both:

1. **A deck without `scaling` is ambiguous.** Old decks omit `scaling` and mean
   SI. New decks would omit `scaling` and mean TCAD. The numbers differ by
   6-12 orders of magnitude and both parse successfully. Silent misreading is
   the single largest risk in this migration.
2. **Many deck keys carry SI-suffixed names but TCAD-interpreted values.** In
   `unit_scaling` mode `reference_doping_m3: 1.0e16` means `1.0e16 cm^-3`, not
   `m^-3`. Keeping those names after unification makes the trap permanent.

## 2. Canonical internal unit system

After the migration there is exactly one internal unit system.

| Quantity | Internal unit |
| --- | --- |
| length | um |
| area | um^2 |
| volume | um^3 |
| concentration | cm^-3 |
| sheet density | cm^-2 |
| mobility | cm^2/(V s) |
| diffusivity | cm^2/s |
| velocity | cm/s |
| electric field | V/cm |
| inverse length | cm^-1 |
| surface field coefficient | cm/V |
| current density | A/cm^2 |
| Auger coefficient | cm^6/s |
| terminal current per device depth | A/um |
| potential | V |
| temperature | K |
| time | s |
| energy | eV |

Composite factors that the assemblers rely on are fixed by this table and must
not drift. The authoritative regression values live in `tests/test_scaling.cpp`:

| Factor | Value | Meaning |
| --- | --- | --- |
| `chargeAreaFactor` | `1e-6` | 2-D volumetric charge term in the Poisson RHS |
| `chargeLineFactor` | `1e-2` | 2-D sheet charge term in the Poisson RHS |
| `chargeVolumeFactor` | `1e-12` | 3-D volumetric term (not used by the 2-D Poisson RHS) |
| `chargeSheetFactor` | `1e-8` | 3-D sheet term (not used by the 2-D Poisson RHS) |
| `fieldFromCoordinateDeltaFactor` | `1e4` | V per um of coordinate delta expressed in V/cm |
| `continuitySourceIntegralFactor` | `1e-4` | continuity source integrated over a 2-D control volume |
| `currentPerInternalDepthFactor` | `1e-6` | current density times length to A/um |

`currentDensity x length = 1e-2` is the display-side identity implied by the
last row and must be asserted alongside it.

## 3. Deck format version

`format_version` becomes a required top-level integer key.

```json
{
  "format_version": 2
}
```

Rules, in final state:

- `format_version: 2` — the deck is interpreted in the canonical internal units
  of section 2. This is the only accepted value.
- `format_version` missing — hard error. The message must name the migration
  tool and must not guess a unit system. Silently defaulting to either
  interpretation is prohibited.
- `format_version: 1` — hard error directing to the migration tool.
- Any `scaling` block present — hard error stating that the key was removed and
  pointing at `solver.normalization`.

During the transition (before the final phase) the rules are relaxed exactly
once: a deck may carry `format_version: 2` *or* the legacy pair
`scaling.mode: "unit_scaling"`, and the two must be provably equivalent. The
transitional equivalence test is what licenses the final deletion; it is
deleted together with the legacy path.

## 4. `scaling` is not deleted, it is split

`scaling` currently conflates two unrelated concerns:

1. **Input unit mode** — `scaling.mode`. This is what the migration removes.
2. **Equation non-dimensionalisation references** — `characteristic_length_um`,
   `reference_concentration_cm3`, `reference_mobility_cm2_V_s`, parsed by
   `parseUnitScalingReferenceConfig` in `src/core/UnitScalingSystem.cpp` and
   consumed by `UnitScalingSystem::fromInputs`. These remain fully functional
   solver controls and must survive.

Concern 2 moves to:

```json
"solver": {
  "normalization": {
    "characteristic_length_um": "auto",
    "reference_concentration_cm3": "auto",
    "reference_mobility_cm2_V_s": "auto"
  }
}
```

`"auto"` keeps the existing `std::nullopt` semantics. Key names are already
correct TCAD units and do not change.

`UnitScalingSystem` is renamed to `EquationScaling` so that the type no longer
reads as an input-unit mode. `UnitScalingReferenceConfig` becomes
`EquationScalingReferenceConfig`.

## 5. Field-level unit contract

The migration is **not** "multiply anything that looks like a concentration".
Every deck key falls into exactly one of five classes.

### Class A — TCAD-native input, renamed in v2

The value is already interpreted in TCAD units under `unit_scaling` (the parse
site wraps it in an identity `scaling.xxxToInternal(...)` boundary method). Only
the key **name** is wrong. v2 renames it; the v1 name is **rejected**, not
aliased.

Mesh and geometry:

| v1 key | v2 key | unit |
| --- | --- | --- |
| `mesh.nodes[].x` / `.y` | unchanged | um |
| `sweep.terminal_charge.depth_m` | `depth_um` | um |
| `sweep.terminal_charge.contact_radius` | `contact_radius_um` | um |
| `sweep.stored_charge.depth_m` | `depth_um` | um |

Materials and doping:

| v1 key | v2 key | unit |
| --- | --- | --- |
| `materials[].ni` | `ni_cm3` | cm^-3 |
| `materials[].mun` / `.mup` | `mun_cm2_V_s` / `mup_cm2_V_s` | cm^2/(V s) |
| `doping[].donors` / `.acceptors` | `donors_cm3` / `acceptors_cm3` | cm^-3 |
| `doping[].fixed_charge_m3` | `fixed_charge_cm3` | cm^-3 |
| `regions[].fixed_charge_m3` | `fixed_charge_cm3` | cm^-3 |

Interfaces (sheet quantities are per cm^2, not per um^2):

| v1 key | v2 key | unit |
| --- | --- | --- |
| `interfaces[].sheet_charge_m2` | `sheet_charge_cm2` | cm^-2 |
| `interfaces[].fixed_charge_m2` | `fixed_charge_cm2` | cm^-2 |
| `interfaces[].trap_density_m2` | `trap_density_cm2` | cm^-2 |

Solver core:

| v1 key | v2 key | unit |
| --- | --- | --- |
| `solver.auger_cn_m6_per_s` | `auger_cn_cm6_per_s` | cm^6/s |
| `solver.auger_cp_m6_per_s` | `auger_cp_cm6_per_s` | cm^6/s |
| `solver.carrier_floor_m3` | `carrier_floor_cm3` | cm^-3 |

Mobility (`solver.mobility.*`):

| v1 suffix | v2 suffix | unit |
| --- | --- | --- |
| `*_mu_min_m2_V_s`, `*_mu_const_m2_V_s`, `*_mumin1_m2_V_s`, `*_mumin2_m2_V_s`, `*_mu1_m2_V_s` | `*_cm2_V_s` | cm^2/(V s) |
| `*_nref_m3`, `*_pc_m3`, `*_cr_m3`, `*_cs_m3` | `*_cm3` | cm^-3 |
| `*_saturation_velocity_m_s` | `*_saturation_velocity_cm_s` | cm/s |
| `surface.theta_electron_m_per_V`, `surface.theta_hole_m_per_V` | `*_cm_per_V` | cm/V |
| `surface.reference_field_V_per_m` | `reference_field_V_per_cm` | V/cm |

Recombination and bandgap narrowing:

| v1 key | v2 key | unit |
| --- | --- | --- |
| `solver.bandgap_narrowing.reference_doping_m3` | `reference_doping_cm3` | cm^-3 |
| `solver.srh_doping_dependence.{electron,hole}.reference_doping_m3` | `reference_doping_cm3` | cm^-3 |
| `contacts[].surface_recombination_velocity_m_per_s` | `surface_recombination_velocity_cm_per_s` | cm/s |

`contacts[].surface_recombination_velocity_m_per_s` is a special case: it is
parsed with a raw `get<Real>()` and no unit boundary at all
(`src/boundary/BoundaryCondition.cpp:143-145`), and it has an unsuffixed alias
`surface_recombination_velocity`. Both v1 spellings are rejected in v2 in
favour of the single `_cm_per_s` name, and the parse site gains the missing
unit boundary.

Impact ionization (`solver.impact_ionization.*`):

| v1 key | v2 key | unit |
| --- | --- | --- |
| `minimum_field_V_m`, `switch_field_V_m` | `*_V_per_cm` | V/cm |
| `{electron,hole}_A_m_inv`, `{electron,hole}_a_low_m_inv`, `{electron,hole}_a_high_m_inv` | `*_cm_inv` | cm^-1 |
| `{electron,hole}_B_V_m`, `{electron,hole}_b_low_V_m`, `{electron,hole}_b_high_V_m` | `*_V_per_cm` | V/cm |
| `carrier_velocity_m_s` | `carrier_velocity_cm_s` | cm/s |
| `driving_force_interpolation.{electron,hole}_ref_density_m3` | `*_ref_density_cm3` | cm^-3 |

Diagnostics that take a field threshold:

| v1 key | v2 key | unit |
| --- | --- | --- |
| `sweep.diagnostics.bv_process_probe.max_electric_field_V_per_m` | `max_electric_field_V_per_cm` | V/cm |
| `sweep.diagnostics.avalanche_tracing.seed_field_V_per_m` | `seed_field_V_per_cm` | V/cm |
| `sweep.diagnostics.avalanche_tracing.stop_field_V_per_m` | `stop_field_V_per_cm` | V/cm |

### Class B — compiled SI defaults

The struct default is an SI literal and the parse site converts it with a real
`unitSystem().xxxToInternal(...)` call before the deck override is applied.
These are invisible in the deck but must be re-expressed as TCAD literals when
`legacySI()` is deleted, because the conversion helper disappears with it.

Known members: mobility defaults (`convertMobilityDefaultsToInternal`), impact
ionization defaults (`convertImpactDefaultsToInternal`), Slotboom and
OldSlotboom `Nref`, SRH `referenceDoping`, Auger `Cn`/`Cp`, material database
defaults, `GummelConfig::carrierFloor`, density-gradient defaults.

### Class C — genuinely SI, must stay SI

`solver.band_to_band.*`. `BandToBandTunnelingModel` documents that its stored
parameters are always SI and the model evaluates in SI. Do **not** internalise
it. It keeps its SI key names, keeps its existing `*_cm_*` input aliases that
convert to SI, and gains an explicit unit boundary at the call site where the
internal V/cm field is converted to V/m and the SI generation rate is converted
back to internal cm^-3 s^-1.

Affected keys: `A_m_inv_s_inv_V_inv2`, `A_cm_inv_s_inv_V_inv2`, `B_V_per_m`,
`B_V_per_cm`, `minimum_field_V_per_m`, `minimum_field_V_per_cm`.

### Class D — never scaled

Potentials (`*_V`), temperatures (`*_K`), times (`*_s`), energies (`*_eV`),
and dimensionless coefficients are identical in both unit systems and are not
touched. This includes `taun`, `taup`, `tau_min_s`, `tau_max_s`,
`temperature_K`, `reference_temperature_K`, `phonon_energy_eV`,
`coefficient_eV`, `offset_eV`, all sweep voltages and voltage tolerances,
`eps_r`, `gamma`, `beta`, `temperature_exponent`, and enum-valued keys.

### Class E — device-depth compound units

`sweep.external_circuit.resistance_ohm_um`,
`sweep.voltage_to_current.current_points_A_per_um`, and
`sweep.voltage_to_current.current_tolerance_A_per_um` are already expressed per
micron of device depth. They keep their names and values.

### Rejection rules

The v2 parser rejects, with a message naming the correct key:

- any v1 Class A key name;
- any unknown key ending in a unit suffix (`_m`, `_m2`, `_m3`, `_m_s`,
  `_m2_V_s`, `_V_per_m`, `_m_inv`, `_m6_per_s`, `_m_per_V`, `_m_per_s`),
  except the Class C allow-list.

This turns the ambiguity into a compile-time-like failure rather than a silent
numeric error.

## 6. Output contract

### 6.1 Sweep CSV

Today `src/simulation/DCSweep.cpp` emits an unsuffixed `current_total` column
(plus `current_electron`, `current_hole`, and their drift/diffusion splits) in
raw internal units, and appends the `*_A_per_um` columns **only** when
`isUnitScaling()` is true.

In v2:

- The `*_A_per_um` columns are emitted **unconditionally** and are the only
  canonical current columns.
- The unsuffixed `current_total`, `current_electron`, `current_hole`,
  `current_electron_drift`, `current_electron_diffusion`,
  `current_hole_drift`, and `current_hole_diffusion` columns are **deleted**.
  They are not kept as deprecated aliases; a column whose name claims no unit
  but whose value silently changed meaning is exactly the failure mode this
  migration exists to remove.
- Every consumer must be migrated before the deletion lands. Known consumers:
  `scripts/audit_pn2d_forward_lowfield_parity.py`,
  `scripts/summarize_pn2d_iv_full_range_debug.py`,
  `scripts/diagnose_pn2d_qf_cap_warmstart_steps.py`,
  `scripts/summarize_pn2d_bv_grid_scan.ps1`,
  `scripts/analyze_bvmethods_nmos_iic_postprocess.py`,
  `scripts/run_pn2d_forward_node_volume_policy_acceptance.py`.
  Note that some of these read `current_total` as the *Sentaurus reference*
  column, which is an independent external format and is unaffected.
- Any column that intentionally remains SI must be renamed to say so
  (`*_A_per_m`, `*_V_per_m`). Emitting an SI-named column holding a TCAD value
  is prohibited.

### 6.2 Restart CSV

`src/io/DDSolutionCsv.cpp` already writes SI to disk: densities are converted
out with the inverse of `m3ToInternalConcentration` on write and back in on
read. The on-disk header `node_id,psi,phin,phip,electrons_m3,holes_m3` is
therefore correct and **does not change**.

The migration requirement is a compatibility regression test proving that a
file written by the pre-migration writer is read identically by the
post-migration reader. There is no restart format break to announce.

### 6.3 VTK

Two gaps must be closed before the legacy path is deleted, because they
currently write raw internal values under SI-suffixed or unit-free names:

- `src/simulation/PoissonSimulation.cpp` writes `net_doping_m3` from the raw
  internal doping vector. Rename to `net_doping_cm3`.
- `src/solver/GummelSolver.cpp` writes `Electrons`, `Holes`, and `NetDoping`
  as raw internal values under unit-free names. Either suffix them
  (`Electrons_cm3`, ...) or document the array unit in the same place the field
  arrays are documented.

`writeRecoveredElectricFields` is already correct: it converts explicitly to
V/cm at the writer boundary via `fieldFromCoordinateDeltaFactor` and
`electricFieldVPerMPerInternal`. That pattern is the model for the fixes above.

Mesh coordinates written into VTK are internal um and stay internal um; the
contract is that VTK geometry is in the internal length unit.

### 6.4 JSON diagnostics

Diagnostic keys are out of scope for the renaming pass **except** where the key
name asserts a unit that the value does not have. Those are corrected in the
same commit that touches the producing code, not in a separate sweep.

## 7. Migration tool requirements

A single tool performs the v1 to v2 deck rewrite. It must:

- accept `--dry-run` and print the full diff without writing;
- be idempotent: running it on a v2 deck is a no-op that exits successfully;
- refuse to run on a deck it cannot fully classify, listing the unrecognised
  unit-suffixed keys rather than passing them through;
- rewrite `scaling.mode` decks and `scaling`-less decks differently, using the
  presence of `scaling.mode: "unit_scaling"` to decide whether numeric values
  need a unit conversion or only a rename;
- move `scaling.{characteristic_length_um,reference_concentration_cm3,reference_mobility_cm2_V_s}`
  into `solver.normalization`, then drop the `scaling` block;
- insert `format_version: 2`;
- print a per-key conversion summary (key, old value, new key, new value,
  applied factor);
- cover `configs/templates/*.template.json` and
  `configs/schema/vela-simulation.schema.json` in addition to `examples/` and
  test fixtures.

## 8. Acceptance

The migration is complete when all of the following hold:

- `grep -r "legacySI\|isUnitScaling\|UnitScalingMode\|scaling.enabled" src include tests` returns nothing.
- No assembler or spec struct has a default argument that selects a unit system;
  all specs are built through a validated factory.
- Every deck under `examples/`, `configs/`, and `tests/` carries
  `format_version: 2` and no `scaling` block.
- The full `ctest` suite is green, including the restart backward-compatibility
  test and the composite-factor regression in `tests/test_scaling.cpp`.
- `docs/config_schema.md` documents a unit for every numeric deck key.
