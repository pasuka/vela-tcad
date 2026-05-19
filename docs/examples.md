# Engineering Examples and Regression Suite

Vela ships small engineering examples that are intentionally coarse so they can
run quickly in CI while still exercising the Poisson, drift-diffusion,
post-processing, and file-output paths. Device-oriented examples live under
`examples/<device>/` with a shared `mesh.json` and named simulation decks such
as `simulation_iv.json`, `simulation_cv.json`, and `simulation_bv.json` when the
corresponding capability is implemented for that device.

## Current Support Matrix

The table below is the source of truth for the example decks. It describes what
is currently present in the repository, not a claim of calibrated production
TCAD coverage.

| Device example | Poisson-only | DD-IV | CV | BV |
| --- | --- | --- | --- | --- |
| `examples/pn_diode` | Implicit equilibrium in DD decks | `simulation_iv.json` forward IV | `simulation_cv.json` quasi-static junction C-V | `simulation_bv.json` reverse-bias diagnostics |
| `examples/nmos2d_dd` | Implicit equilibrium in DD decks | `simulation_iv.json` Id-Vd and `simulation_idvg.json` Id-Vg | `simulation_cv.json` gate/body quasi-static CV | Not yet |
| `examples/nmos2d_mos_dd` | Implicit equilibrium in DD decks on a mixed Si/SiO2 mesh | `simulation_iv.json` Id-Vd and `simulation_idvg.json` Id-Vg Si/SiO2 MOS DD prototype smoke sweeps | `simulation_cv.json` multi-terminal quasi-static MOS CV prototype | `simulation_bv.json` off-state drain high-field diagnostic smoke sweep |
| `examples/pmos2d_dd` | Implicit equilibrium in DD decks | `simulation_iv.json` Id-Vd and `simulation_idvg.json` Id-Vg | `simulation_cv.json` gate/body quasi-static CV | Not yet |
| `examples/pmos2d_mos_dd` | Implicit equilibrium in DD decks on a mixed Si/SiO2 mesh | `simulation_iv.json` negative-bias Id-Vd and `simulation_idvg.json` negative-gate Id-Vg Si/SiO2 MOS DD prototype smoke sweeps | `simulation_cv.json` multi-terminal quasi-static MOS CV prototype | `simulation_bv.json` off-state negative-drain high-field diagnostic smoke sweep |
| `examples/ldmos2d` | `simulation_iv.json` remains a reverse-biased Poisson field-distribution deck on the mixed Si/SiO2 mesh | `simulation_dd_iv.json` low-voltage on-state Id-Vd prototype smoke sweep (**prototype / not calibrated**) | Not yet | `simulation_bv.json` baseline off-state diagnostics (**prototype / not calibrated**) and `simulation_bv_fieldplate.json` RESURF/field-plate trend variant (**prototype / not calibrated**) |
| `examples/igbt2d` | `simulation_poisson.json` PNPN/drift off-state baseline | `simulation_iv.json` low-current collector DD prototype and `simulation_high_injection_iv.json` high-injection IV prototype smoke deck (**prototype / not calibrated**) | `simulation_charge_cv.json` quasi-static gate/collector/emitter charge + stored-charge proxy smoke deck (**prototype / not calibrated**) | `simulation_bv.json` off-state collector reverse-bias diagnostics prototype and `simulation_bv_ii.json` impact-ionization smoke variant (**prototype / not calibrated**) |
| `examples/moscap` | `simulation.json` Si/SiO2 Poisson interface check | Not applicable | Not yet | Not applicable |
| `examples/nmos2d` | `simulation.json` legacy simplified NMOS Poisson check | Not yet | Not yet | Not yet |
| `examples/schottky_diode_2d` | Implicit equilibrium in DD deck | `simulation_iv.json` prototype Schottky-anode IV smoke (M1.3, not a calibrated Schottky model) | Not yet | Not yet |


Planned extensions should keep this table conservative. The
`examples/ldmos2d/simulation_dd_iv.json` deck is a low-voltage,
multi-terminal mixed Si/SiO2 drift-diffusion on-state IV prototype with
regression coverage for finite convergence and current trend only; it is not a
calibrated LDMOS model and does not claim calibrated breakdown-voltage
coverage. The LDMOS BV decks are limited to executable off-state diagnostics and
trend checks;
IGBT now includes low-current plus high-injection IV prototype smoke decks with recombination/lifetime controls and trend checks; this remains an engineering prototype and is not a calibrated IGBT or calibrated on-state voltage model.
The `examples/ldmos2d/simulation_bv.json` and
`examples/ldmos2d/simulation_bv_fieldplate.json` pair form a coarse
RESURF/field-plate **trend smoke** (prototype / not calibrated): the variant extends the gate field-plate
overlap and lightens drift doping to compare 2-D peak-field tendencies
(`max_electric_field_V_per_m`) and final drain current under the same BV sweep.
This is intentionally a loose regression guard for field-distribution
directionality and finite diagnostics, not process optimization or calibrated
breakdown-voltage prediction.
The `examples/nmos2d_mos_dd` and
`examples/pmos2d_mos_dd` decks are true mixed-material Si/SiO2 MOS
drift-diffusion prototypes with metal-gate contacts over gate oxide and ohmic
source/drain/body contacts. They are intentionally coarse CI smoke decks for
polarity, multi-terminal charge/CV, and off-state drain reverse-bias
high-field/BV diagnostics, and must not be described as calibrated MOSFETs.
The mixed-material MOS DD CV decks use the
multi-terminal quasi-static CV prototype: a single gate-bias sweep is finite-
differenced to report terminal charges plus `Cgate_gate`, `Cgate_drain`,
`Cgate_source`, and `Cgate_body`-style
columns such as `capacitance_Cgate_gate_F_per_m` and
`capacitance_Cgate_drain_F_per_m`. These
columns are useful engineering smoke diagnostics, not AC small-signal admittance
or capacitance matrix results. The mixed-material MOS DD BV decks hold gate,
source, and body at off-state bias, sweep the drain in the reverse/high-field
direction with `sweep.mode: "bv_reverse"`, write per-point VTK snapshots, and
regress finite leakage current plus non-decreasing maximum edge-field columns;
these are smoke diagnostics, not calibrated breakdown-voltage predictions.

## Contact Schema Compatibility

Decks below use either the original `contacts[]` shape (`name`, `bias`, and an
optional `flatband_voltage` or `work_function_eV`) or the explicit typed form.
Starting with the unified
boundary parser, contact entries may also include a string `type` field. The
legacy untyped form is treated as `"type": "ohmic"`, so the example decks here
keep working without changes. The parser accepts case-insensitive names with
either `-` or `_` separators (for example `metal_gate` and `Metal-Gate` both
resolve to the same metal-gate model). Recognised values are `ohmic`,
`dirichlet`, `metal_gate`, `schottky`, and `floating`. The Poisson driver
currently maps `ohmic`, `dirichlet`, and `metal_gate` to an effective
Dirichlet potential; `schottky` contacts use a prototype Dirichlet-barrier
model that pins the surface carrier density via a Boltzmann factor
`n = Nc*exp(-phi_Bn/Vt)` and is supported by the Gummel solver and the Poisson
driver (the Newton solver rejects Schottky contacts with a clear error until a
future milestone). `floating` contacts are reserved by the schema and will
raise a clear error if a deck requests them today. See `docs/config_schema.md`
for the central field-level schema reference.


## Boundary Schema (Neumann / Insulating / Symmetry)

In addition to `contacts[]` for Dirichlet boundary conditions, decks may
include an optional top-level `boundaries[]` array to declare explicit
non-Dirichlet boundary conditions for the Poisson equation. Example:

```json
"boundaries": [
  {
    "name": "left_symmetry",
    "type": "symmetry",
    "node_ids": [0, 4, 8]
  },
  {
    "name": "top_neumann",
    "type": "neumann",
    "node_ids": [8, 9, 10],
    "normal_displacement_C_per_m2": 0.0
  }
]
```

Supported `type` values (case-insensitive, with `-` or `_` separators):

- `neumann`: Specifies the normal displacement on the segment via
  `normal_displacement_C_per_m2`, in units of C/m^2. Positive values describe
  flux pointing out of the simulation domain; the RHS contribution per polyline
  edge is `value * edge_length / 2` to each endpoint.
- `insulating`: Equivalent to zero normal displacement (`D_n = 0`). Accepted by
  parsers; does not change the matrix or RHS.
- `symmetry`: Same numerical meaning as `insulating`; intended for symmetry
  planes such as half-device decks.

Each `boundaries[]` entry must define `node_ids` as a polyline of at least two
existing mesh node ids; the segment length is taken from the Euclidean distance
between consecutive nodes. The unified parser rejects unknown `type` strings,
`node_ids` shorter than 2, and non-finite Neumann values. Naturally-bounded
edges that are not declared explicitly continue to behave as zero-flux, exactly
as before this milestone, so existing decks need no changes. Explicit
`boundaries[]` are currently consumed by the Poisson driver. DD sweeps keep
their historical natural zero-flux behavior and do not yet parse non-contact
boundary segments from deck JSON.

## Common Workflow

Build the project first:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
```

On Windows, run the same commands from the MSYS2 **UCRT64** shell. If MSYS2 is
installed at `D:\msys64`, the corresponding tools live under
`D:\msys64\ucrt64\bin`. Use `python` instead of `python3` below if the shell
does not provide a `python3` alias.

Run all examples through the regression harness:

```bash
python3 scripts/run_regression.py --runner build/vela_example_runner
```

CTest exposes the same entry point:

```bash
ctest --test-dir build --output-on-failure -R regression
```

The regression script copies examples into `build/regression_output/`, creates
per-example `outputs/` folders, runs `vela_example_runner --config <deck>`,
verifies expected files, checks convergence indicators, scans CSV/VTK text for
NaN/Inf values, and writes `build/regression_output/regression_summary.json`.

## Python API Example

The optional pybind11 API is a thin wrapper over the C++ core. Python examples
only document curve features that are already implemented and tested in C++
(DD-IV, quasi-static CV, and reverse-bias BV diagnostics for the decks listed
in this guide); experimental C++-only behavior should not be presented as a
Python-supported API until it has matching tests. Configure with
`VELA_ENABLE_PYTHON=ON` and make the generated package importable:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DVELA_ENABLE_PYTHON=ON
cmake --build build --parallel
PYTHONPATH=build/python/Debug python3 examples/python/run_device_curves.py
```

On Windows PowerShell, set `PYTHONPATH` with:

```powershell
$env:PYTHONPATH = "build\python\Debug"
python examples\python\run_device_curves.py
```

The script calls the lightweight helpers `vela.run_iv_curve()`,
`vela.run_cv_curve()`, and `vela.run_bv_curve()`. Those helpers still call the
C++ DC sweep binding internally and return the same point dictionaries as
`vela.run_dc_sweep()`, including `curve_type`, contact names, convergence
diagnostics, current/capacitance/breakdown fields, `output_csv`, and per-point
`output_vtk` paths when VTK output is enabled. The example covers PN IV/CV/BV
and the NMOS Id-Vd sweep, writing the usual CSV/VTK outputs under each
device's `outputs/` directory.

CTest registers the Python coverage as `python_api` when the extension is
enabled:

```bash
ctest --test-dir build --output-on-failure -R python_api
```

## PN Diode (`examples/pn_diode`)

**Physical meaning.** A tiny 2-D silicon p-n diode validates the
self-consistent drift-diffusion flow, contact-current extraction, quasi-static
terminal charge extraction, and reverse-bias breakdown diagnostics. The decks
are CI smoke cases rather than calibrated diode models.

**Contact schema note.** PN decks can omit `contacts[].type` (legacy-compatible
default `ohmic`) or set it explicitly to `ohmic`.

**Inputs.**

- `mesh.json`: 1 um x 1 um triangular silicon mesh split into `p_region` and
  `n_region`, with `anode` and `cathode` contacts.
- `simulation_iv.json`: forward anode IV from 0 V to 0.5 V.
- `simulation_cv.json`: reverse/anode quasi-static junction capacitance from
  0 V to -0.2 V using terminal-charge differencing.
- `simulation_bv.json`: reverse/anode sweep from 0 V to -0.5 V with maximum
  electric-field and current-jump diagnostic columns, plus
  last-stable/failed-bias diagnostics when non-convergence is treated as
  breakdown. The regression block requires the reported maximum edge-field
  diagnostic to be non-decreasing along the reverse-bias sweep, with a small
  floating-point tolerance.
- `newton_simulation.json`: single equilibrium Newton solve used by the CLI
  Newton smoke test.

**Run directly.**

```bash
mkdir -p examples/pn_diode/outputs
build/vela_example_runner --config examples/pn_diode/simulation_iv.json
build/vela_example_runner --config examples/pn_diode/simulation_cv.json
build/vela_example_runner --config examples/pn_diode/simulation_bv.json
```

The default PN diode IV/BV/CV sweeps use `solver.method: "gummel"`. To
exercise the coupled Newton CLI entry point for a single equilibrium bias point,
run:

```bash
build/vela_example_runner --config examples/pn_diode/newton_simulation.json
```

DC sweeps can also opt into Newton by setting `solver.method` (or
`solver.type`) to `"newton"` in the same `solver` block used for Newton
tolerances such as `max_iter`, `reltol`, `abstol`, `damping_factor`, and
`line_search`. Omitting `solver.method` / `solver.type` keeps the historical
Gummel sweep path.

**Outputs.**

- `examples/pn_diode/outputs/pn_iv.csv`: voltage, electron current, hole
  current, total current, convergence flag, and iteration count.
- `examples/pn_diode/outputs/pn_cv.csv`: quasi-static charge and capacitance
  columns in addition to the common sweep diagnostics.
- `examples/pn_diode/outputs/pn_bv.csv`: maximum electric field, current-jump
  ratio, breakdown flag, breakdown voltage, criterion, `last_stable_bias`,
  `failed_bias`, and `failure_reason` columns. The maximum electric field is an
  edge-difference engineering diagnostic (`max(|delta psi| / edge_length)`), not
  a high-order field reconstruction. The PN BV regression checks that this
  diagnostic does not decrease as reverse bias increases, within a small
  tolerance. If `breakdown.non_convergence` is enabled and a reverse-bias point
  cannot be continued, the diagnostic row reports
  `criterion=last_stable_before_nonconvergence` and uses the last converged bias
  as the reported breakdown voltage.
- `examples/pn_diode/outputs/pn_sweep_*.vtk` and
  `examples/pn_diode/outputs/pn_bv_sweep_*.vtk`: potential, quasi-Fermi
  levels, carrier densities, and net doping snapshots for converged points.

## NMOS Drift-Diffusion (`examples/nmos2d_dd`)

**Physical meaning.** A deliberately coarse, all-silicon NMOS-like
cross-section exercises the drift-diffusion Gummel path with four independently
biased terminals (`body`, `source`, `gate`, and `drain`). The source/drain
implant regions are abrupt n+ boxes next to a p-type body/channel.

**Current support level.** DD-IV and CV. This is a CI stability baseline, not a
calibrated transistor model. BV and oxide-interface MOS electrostatics are not
claimed for this deck.

**Contact schema note.** NMOS decks remain compatible with omitted
`contacts[].type` and may also use explicit `ohmic` / `metal_gate` entries.

**Inputs.**

- `mesh.json`: all-silicon body, source implant, and drain implant regions.
- `simulation_iv.json`: Id-Vd drain sweep at fixed small positive gate bias.
- `simulation_idvg.json`: Id-Vg gate sweep at fixed drain bias.
- `simulation_cv.json`: gate/body quasi-static CV sweep on the same DD mesh.
- `test_mos_solver_crosscheck`: CTest single-bias smoke coverage that solves
  the same low-bias NMOS point with Gummel and Newton, checks both solutions
  are finite and converged, and compares drain-current sign and order of
  magnitude.

**Run directly.**

```bash
mkdir -p examples/nmos2d_dd/outputs
build/vela_example_runner --config examples/nmos2d_dd/simulation_iv.json
build/vela_example_runner --config examples/nmos2d_dd/simulation_idvg.json
build/vela_example_runner --config examples/nmos2d_dd/simulation_cv.json
```

**Outputs and acceptance metrics.**

- `examples/nmos2d_dd/outputs/nmos2d_dd_iv.csv`: Id-Vd sweep rows with
  `converged=1` for each declared point.
- `examples/nmos2d_dd/outputs/nmos2d_dd_idvg.csv`: fixed-drain Id-Vg rows.
- `examples/nmos2d_dd/outputs/nmos2d_dd_cv.csv`: gate/body charge and
  capacitance columns.
- `examples/nmos2d_dd/outputs/nmos2d_dd_sweep_*.vtk`: finite potential,
  quasi-Fermi, carrier-density, and net-doping fields.

## PMOS Drift-Diffusion (`examples/pmos2d_dd`)

**Physical meaning.** The PMOS deck is the polarity-complement counterpart to
the NMOS stability case. It uses abrupt p+ source/drain implants in an n-type
body to exercise sign handling, descending DC sweep targets, and contact-bias
application in the drift-diffusion solver.

**Current support level.** DD-IV and CV. BV is not yet covered.

**Contact schema note.** PMOS decks remain compatible with omitted
`contacts[].type` and may also use explicit `ohmic` / `metal_gate` entries.

**Inputs.**

- `mesh.json`: all-silicon n-body, p+ source implant, and p+ drain implant
  regions.
- `simulation_iv.json`: Id-Vd drain sweep at fixed small negative gate bias.
- `simulation_idvg.json`: Id-Vg gate sweep at fixed drain bias.
- `simulation_cv.json`: gate/body quasi-static CV sweep on the same DD mesh.
- `test_mos_solver_crosscheck`: CTest single-bias smoke coverage that solves
  the same low-bias PMOS point with Gummel and Newton, checks both solutions
  are finite and converged, and compares drain-current sign and order of
  magnitude.

**Run directly.**

```bash
mkdir -p examples/pmos2d_dd/outputs
build/vela_example_runner --config examples/pmos2d_dd/simulation_iv.json
build/vela_example_runner --config examples/pmos2d_dd/simulation_idvg.json
build/vela_example_runner --config examples/pmos2d_dd/simulation_cv.json
```

**Outputs and acceptance metrics.**

- `examples/pmos2d_dd/outputs/pmos2d_dd_iv.csv`: Id-Vd sweep rows with
  `converged=1` for each declared point.
- `examples/pmos2d_dd/outputs/pmos2d_dd_idvg.csv`: fixed-drain Id-Vg rows.
- `examples/pmos2d_dd/outputs/pmos2d_dd_cv.csv`: gate/body charge and
  capacitance columns.
- `examples/pmos2d_dd/outputs/pmos2d_dd_sweep_*.vtk`: finite potential,
  quasi-Fermi, carrier-density, and net-doping fields.

## LDMOS (`examples/ldmos2d`)

**Physical meaning.** A coarse lateral DMOS-like mixed Si/SiO2 mesh covers a
p-body / n-drift lateral junction with n+ source/drain regions and a field-oxide
slab. The Poisson deck catches material-interface regressions under a small
reverse drain bias. The drift-diffusion decks are engineering smoke/prototype
cases for off-state diagnostics and low-bias current continuity; they are not
calibrated LDMOS performance models.

**Current support level.** Poisson baseline, low-bias DD-IV smoke sweep, and an
off-state drain reverse-bias BV diagnostic prototype. The BV deck reports
breakdown-oriented diagnostics from the generic `bv_reverse` DC sweep mode, but
`breakdown_voltage` is only the first configured diagnostic crossing (or
last-stable-before-nonconvergence marker) and must not be interpreted as a
calibrated avalanche breakdown prediction.

**Contact schema note.** LDMOS decks use explicit `ohmic` source/body/drain
contacts and a `metal_gate` contact. In the BV prototype, source/body/gate are
held at 0 V and the drain is swept on the positive high-voltage/off-state side,
which stays on one reverse-bias polarity side for `bv_reverse`.

**Inputs.**

- `mesh.json`: p-body/drift silicon, n+ source, n drain/drift, field oxide, and
  body/source/gate/drain contacts.
- `simulation_iv.json`: reverse-biased Poisson field-distribution deck; the
  filename is retained so every device directory has a stable IV deck slot, but
  this file itself is not the DD-IV sweep.
- `simulation_dd_iv.json`: low-bias drain DD-IV smoke sweep with Gummel.
- `simulation_bv.json`: LDMOS off-state drain `bv_reverse` sweep with
  source/body/gate at 0 V, `impact_ionization.model: "none"`, sweep VTK output,
  and BV diagnostic CSV columns: `max_electric_field_V_per_m`,
  `current_jump_ratio`, `breakdown_detected`, `breakdown_voltage`, `criterion`,
  `last_stable_bias`, `failed_bias`, and `failure_reason`. The max-field value
  is the edge-difference diagnostic (`max(|delta psi| / edge length)`) from the
  mesh potential, not a calibrated high-order electric-field reconstruction or
  calibrated BV threshold.

**Run directly.**

```bash
mkdir -p examples/ldmos2d/outputs
build/vela_example_runner --config examples/ldmos2d/simulation_iv.json
build/vela_example_runner --config examples/ldmos2d/simulation_dd_iv.json
build/vela_example_runner --config examples/ldmos2d/simulation_bv.json
```

**Outputs and acceptance metrics.**

- `examples/ldmos2d/outputs/ldmos2d_reverse_poisson.vtk`: electrostatic
  potential and net doping across silicon and field oxide.
- `examples/ldmos2d/outputs/ldmos2d_dd_iv.csv` and
  `examples/ldmos2d/outputs/ldmos2d_dd_iv_sweep_*.vtk`: finite low-bias DD-IV
  current rows with monotone absolute drain-current trend.
- `examples/ldmos2d/outputs/ldmos2d_bv.csv` and
  `examples/ldmos2d/outputs/ldmos2d_bv_sweep_*.vtk`: finite off-state BV
  diagnostic rows. Regression requires the BV schema to be present,
  `max_electric_field_V_per_m` to be finite and non-decreasing across converged
  rows, and any allowed final non-converged BV row to report valid
  last-stable/failed-bias diagnostics.


## IGBT (`examples/igbt2d`)

**Physical meaning.** A coarse vertical IGBT-like stack contains a p-type
collector, n buffer, lightly doped n drift, and top p-base/emitter region. The
Poisson deck targets off-state electrostatic stability; the IV deck is a
low-current drift-diffusion prototype, and BV decks are reverse-bias diagnostic
smoke sweeps only.

**Current support level.** Poisson baseline, low-current DD-IV, high-injection
IV smoke trends, and reverse-bias BV diagnostics prototype. These decks are not
calibrated IGBT or calibrated avalanche-breakdown predictions.

**Contact schema note.** IGBT decks remain compatible with omitted
`contacts[].type`; explicit `ohmic` / `metal_gate` typing is supported by the
shared parser.

**Inputs.**

- `mesh.json`: vertical PNPN/drift stack with collector, gate, and emitter
  contacts.
- `simulation_poisson.json`: off-state collector-bias Poisson baseline.
- `simulation_iv.json`: low-current collector sweep DD prototype.
- `simulation_bv.json`: off-state (`emitter=0 V`, `gate=0 V`) collector
  reverse-bias BV/leakage diagnostics with `impact_ionization.model: "none"`.
- `simulation_bv_ii.json`: optional Selberherr impact-ionization variant for a
  short off-state BV smoke comparison path.

**Run directly.**

```bash
mkdir -p examples/igbt2d/outputs
build/vela_example_runner --config examples/igbt2d/simulation_poisson.json
build/vela_example_runner --config examples/igbt2d/simulation_iv.json
build/vela_example_runner --config examples/igbt2d/simulation_bv.json
build/vela_example_runner --config examples/igbt2d/simulation_bv_ii.json
```

**Outputs and acceptance metrics.**

- `examples/igbt2d/outputs/igbt2d_poisson.vtk`: electrostatic potential and net
  doping through the vertical stack.
- `examples/igbt2d/outputs/igbt2d_iv.csv`: low-current collector sweep with
  convergence and step diagnostics.
- `examples/igbt2d/outputs/igbt2d_iv_sweep_*.vtk`: DD field snapshots for the
  low-current prototype.
- `examples/igbt2d/outputs/igbt2d_bv.csv` and
  `examples/igbt2d/outputs/igbt2d_bv_sweep_*.vtk`: off-state collector reverse-bias
  diagnostics that include `max_electric_field_V_per_m`,
  `current_jump_ratio`, `breakdown_detected`, `breakdown_voltage`, `criterion`,
  plus last-stable/failed-bias non-convergence diagnostics.
- `examples/igbt2d/outputs/igbt2d_bv_ii.csv` and
  `examples/igbt2d/outputs/igbt2d_bv_ii_sweep_*.vtk`: short Selberherr impact-ionization
  smoke variant for finite/trend comparison against the no-impact-ionization
  baseline. This is a diagnostic trend check, not a calibrated BV predictor.

## MOS Capacitor (`examples/moscap`)

**Physical meaning.** A two-material Si/SiO2 MOS capacitor validates Poisson
assembly across a dielectric interface. The body contact is grounded and the
gate is biased at 1 V. No mobile-carrier equation is solved in this first
regression; the case is Poisson-only.

**Inputs.**

- `mesh.json`: layered triangular mesh with a silicon slab, a silicon-dioxide
  slab, shared interface nodes at `y = 0.30 um`, and `body`/`gate` contacts.
- `simulation.json`: Poisson configuration with nonzero oxide
  `fixed_charge_m3` and a `regression` block describing interface probe
  locations and tolerances.

**Run directly.**

```bash
mkdir -p examples/moscap/outputs
build/vela_example_runner --config examples/moscap/simulation.json
```

**Outputs and checks.**

- `examples/moscap/outputs/moscap_poisson.vtk`: electrostatic potential and net
  doping.
- Regression confirms the shared-node potential jump at the Si/SiO2 interface
  is effectively zero and that the normal electric-displacement estimates on
  the silicon and oxide sides are approximately continuous.

## Legacy Simplified NMOS (`examples/nmos2d`)

**Physical meaning.** A coarse rectangular NMOS cross-section exercises a mixed
semiconductor/oxide mesh with four contacts: `source`, `drain`, `body`, and
`gate`. This legacy case is Poisson-only and uses rectangular n+ source/drain
implant regions plus a p-type body as a CI-friendly approximation of analytic
implant profiles.

**Inputs.**

- `mesh.json`: p-body silicon, n+ source/drain implant boxes, top gate oxide,
  and four contacts.
- `simulation.json`: fixed doping for each region and a single Poisson bias
  point (`Vg = 0.8 V`, `Vd = 0.1 V`, source/body grounded).

**Run directly.**

```bash
mkdir -p examples/nmos2d/outputs
build/vela_example_runner --config examples/nmos2d/simulation.json
```

**Outputs.**

- `examples/nmos2d/outputs/nmos_poisson.vtk`: potential and net doping for
  quick inspection in a VTK viewer.

## Adding Future Regression Cases

Use this minimal checklist when adding a new device regression case:

1. Create or update `examples/<device_name>/` with a small `mesh.json` and a
   named simulation deck such as `simulation_iv.json`, `simulation_cv.json`, or
   `simulation_bv.json` that exercises only implemented solver features.
2. Keep generated files under `examples/<device_name>/outputs/`; do not commit
   those outputs unless a task explicitly asks for golden artifacts.
3. Add the case to `EXAMPLES` in `scripts/run_regression.py` with every expected
   CSV/VTK file that the runner should verify.
4. Add the smallest physics-specific regression check that matches the current
   implementation, such as convergence flags, NaN/Inf scanning, interface
   continuity, current-direction sanity, maximum edge-field range bounds, or a
   monotonic trend already supported by the model.
5. Document the case in this file with physical meaning, inputs, direct run
   command, outputs, and any special checks; keep the support matrix from
   over-promising.
6. If the new case exposes solver, mesh, physics, or discretization behavior,
   add or update the closest Catch2 test in `tests/` instead of relying only on
   the example.
7. Build and run `ctest --test-dir build --output-on-failure -R regression`
   before committing; run the full CTest suite when core code changes are part
   of the regression update.
