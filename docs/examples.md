# Engineering Examples And Regression Suite

Vela examples are small engineering fixtures. They are intentionally coarse so
the full suite can run in CI while still exercising Poisson, drift-diffusion,
post-processing, file output, and regression assertions.

The examples are not calibrated production TCAD decks. When a device row says
"prototype" or "diagnostic", it means the deck is used for finite-output,
polarity, monotonicity, and trend checks only.

## How Examples Are Run

Build first:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
```

Run one deck:

```bash
build/vela_example_runner --config examples/pn_diode/simulation_iv.json
```

Run the checked-in engineering regression suite:

```bash
python scripts/run_regression.py --runner build/vela_example_runner
```

Python helper examples are available under `examples/python/`:

- `run_pn_diode.py`: minimal script-style entry for the PN diode decks.
- `run_device_curves.py`: helper entry for running multiple curve families.

These scripts are convenience wrappers for local workflows and are not part of
the core CTest regression matrix.

CTest exposes the same suite:

```bash
ctest --test-dir build --output-on-failure -R regression
```

The regression script copies each selected deck into
`build/regression_output/<case>/`, canonicalizes the selected deck as
`simulation.json`, runs `vela_example_runner`, verifies expected CSV/VTK files,
scans outputs for NaN/Inf, applies configured trend checks, and writes
`build/regression_output/regression_summary.json`.

## Support Matrix

This table mirrors the current example deck inventory and the regression
coverage in `scripts/run_regression.py`.

| Example | Poisson | IV / Id-Vd | Id-Vg | CV | BV / high field | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| `examples/pn_diode` | Implicit equilibrium in DD decks | `simulation_iv.json` | Not applicable | `simulation_cv.json` | `simulation_bv.json` | PN smoke deck plus `newton_simulation.json` single-bias Newton entry. |
| `examples/nmos2d_dd` | Implicit equilibrium in DD decks | `simulation_iv.json` | `simulation_idvg.json` | `simulation_cv.json` | Not present | All-silicon NMOS-like baseline. |
| `examples/pmos2d_dd` | Implicit equilibrium in DD decks | `simulation_iv.json` | `simulation_idvg.json` | `simulation_cv.json` | Not present | All-silicon PMOS-like baseline. |
| `examples/nmos2d_mos_dd` | Implicit equilibrium in DD decks | `simulation_iv.json` | `simulation_idvg.json`, `simulation_idvg_surface.json` | `simulation_cv.json` | `simulation_bv.json` | Mixed Si/SiO2 MOS prototype with surface-mobility smoke deck. |
| `examples/pmos2d_mos_dd` | Implicit equilibrium in DD decks | `simulation_iv.json` | `simulation_idvg.json` | `simulation_cv.json` | `simulation_bv.json` | Mixed Si/SiO2 PMOS prototype. |
| `examples/ldmos2d` | `simulation_iv.json` is a Poisson field-distribution deck | `simulation_dd_iv.json`, `simulation_dd_iv_unit_scaling.json` | Not present | Not present | `simulation_bv*.json` | Power-device trend prototype; regression uses unit-scaling DD-IV and BV decks. |
| `examples/igbt2d` | `simulation_poisson.json` | `simulation_iv*.json`, `simulation_high_injection_iv*.json` | Not present | `simulation_charge_cv*.json` | `simulation_bv*.json`, `simulation_bv_ii*.json` | Power-device trend prototype with stored-charge and impact-ionization smoke variants. |
| `examples/moscap` | `simulation.json` | Not applicable | Not applicable | `simulation_cv.json` exists as a deck artifact but is not in the main regression matrix | Not applicable | Si/SiO2 Poisson interface and displacement continuity check. |
| `examples/nmos2d` | `simulation.json` | Not present | Not present | Not present | Not present | Legacy simplified Poisson-only NMOS cross-section. |
| `examples/schottky_diode_2d` | Implicit equilibrium in DD deck | `simulation_iv.json` | Not applicable | Not present | Not present | Prototype Schottky-anode Gummel IV smoke deck. |

No `scaling` field keeps legacy SI input units. Decks with:

```json
"scaling": { "mode": "unit_scaling" }
```

use common external TCAD input units at the schema boundary. See
[config_schema.md](config_schema.md) for exact unit rules.

## Current Regression Cases

`scripts/run_regression.py` currently runs:

- PN diode IV, CV, and BV.
- MOS capacitor Poisson interface check.
- Legacy simplified NMOS Poisson check.
- All-silicon NMOS/PMOS DD IV, Id-Vg, and CV.
- Mixed-material NMOS MOS DD IV, Id-Vg, surface-mobility Id-Vg, CV, and BV.
- Mixed-material PMOS MOS DD IV, Id-Vg, CV, and BV.
- LDMOS Poisson baseline, unit-scaling low-bias DD-IV, unit-scaling BV, and
  unit-scaling field-plate BV variant.
- IGBT Poisson baseline, unit-scaling IV, high-injection IV, charge/CV,
  BV, and impact-ionization BV variant.
- Schottky diode IV prototype.

Regression checks include finite-output scanning, CSV convergence policy,
expected row counts, step-control diagnostics, monotone current or max-field
checks, surface-mobility current-ratio checks, LDMOS field-plate comparison,
IGBT high-injection/stored-charge checks, and Schottky current trend checks.

## Device Notes

### PN Diode

`examples/pn_diode` validates the end-to-end drift-diffusion flow on a small
abrupt junction. It emits IV, quasi-static CV, and reverse-bias/BV diagnostic
CSV files plus optional VTK snapshots. `examples/pn_diode/newton_simulation.json`
is the single-bias Newton CLI smoke deck.

Run:

```bash
build/vela_example_runner --config examples/pn_diode/simulation_iv.json
build/vela_example_runner --config examples/pn_diode/simulation_cv.json
build/vela_example_runner --config examples/pn_diode/simulation_bv.json
build/vela_example_runner --config examples/pn_diode/newton_simulation.json
```

### All-Silicon MOS DD

`examples/nmos2d_dd` and `examples/pmos2d_dd` are all-silicon MOS-like
stability fixtures. They exercise multi-terminal biasing, DD current
extraction, Id-Vg trend behavior, and quasi-static CV. They are not calibrated
MOSFET models.

### Mixed-Material MOS DD

`examples/nmos2d_mos_dd` and `examples/pmos2d_mos_dd` use Si/SiO2 meshes with
ohmic source/drain/body contacts and a `metal_gate` contact. Their CV decks use
the multi-terminal quasi-static terminal-charge prototype, and their BV decks
are off-state high-field diagnostics. The NMOS family also carries an
Id-Vg surface-mobility smoke deck using `caughey_thomas_field_surface`.

These are engineering prototypes. They should not be described as calibrated
MOSFET models or AC capacitance extractors.

### LDMOS

`examples/ldmos2d` contains:

- `simulation_iv.json`: historical Poisson field-distribution deck.
- `simulation_dd_iv.json`: legacy SI low-bias DD-IV smoke deck.
- `simulation_dd_iv_unit_scaling.json`: unit-scaling DD-IV deck used by the
  main regression.
- `simulation_bv.json`: legacy SI off-state diagnostic deck.
- `simulation_bv_unit_scaling.json`: unit-scaling BV deck used by the main
  regression.
- `simulation_bv_fieldplate*.json`: field-plate/RESURF trend variants.

The LDMOS decks validate finite outputs, monotone current/max-field trends, and
relative field-plate behavior. They do not claim calibrated breakdown voltage.

### IGBT

`examples/igbt2d` contains:

- `simulation_poisson.json`: PNPN/drift off-state Poisson baseline.
- `simulation_iv*.json`: low-current collector IV prototype.
- `simulation_high_injection_iv*.json`: high-injection trend smoke deck with
  stored-charge-related physics knobs.
- `simulation_charge_cv*.json`: quasi-static gate/collector/emitter charge and
  stored-charge proxy deck.
- `simulation_bv*.json`: off-state collector high-field diagnostics.
- `simulation_bv_ii*.json`: Selberherr impact-ionization smoke variant.

The IGBT decks are engineering trend validation fixtures only.

### Schottky Diode

`examples/schottky_diode_2d/simulation_iv.json` exercises the prototype
Schottky contact path through the Gummel sweep. Schottky is not available in
Newton sweeps.

## Contact And Boundary Compatibility

Contacts may omit `type`; omitted type is treated as `ohmic`. Recognized
contact types are `ohmic`, `dirichlet`, `metal_gate`, `schottky`, and
`floating`. `floating` is reserved and currently rejected.

Explicit `boundaries[]` segments are consumed by the Poisson driver for
`neumann`, `insulating`, and `symmetry` boundary segments. DD sweeps keep the
historical natural zero-flux behavior on non-contact boundaries.

## Adding A New Regression Case

1. Add a small `examples/<device>/mesh.json` plus one or more named
   `simulation_*.json` decks.
2. Keep generated outputs under `examples/<device>/outputs/`; do not commit
   them unless the task explicitly asks for golden artifacts.
3. Add the deck to `EXAMPLES` in `scripts/run_regression.py`.
4. Add the smallest check that matches implemented behavior: finite output,
   convergence policy, monotone current, max-field trend, interface continuity,
   or a device-specific trend.
5. Update this document and [config_schema.md](config_schema.md) if the deck
   introduces a new public field.
6. Add or update the nearest Catch2 test when solver, mesh, physics, or
   discretization behavior changes.
7. Run `ctest --test-dir build --output-on-failure -R regression`; run the full
   suite when core behavior changed.
