# Sentaurus TDR/HDF5 Import

This document describes the optional Sentaurus import path used to convert
Sentaurus TDR/HDF5 and related text outputs into Vela-friendly neutral fixtures.

TDR and the optional `vela.ddstate/1` HDF5 restart files share the HDF5 container,
but use different schemas. They are not interchangeable by renaming files.
See the [implementation and sample compatibility audit](validation/templates_ldmos_tdr_hdf5_compatibility_2026-09-22.md)
for supported geometry/state layouts, units, node mapping, and conversion limits.

## Scope

The repository provides two complementary import paths:

- `sentaurus_import` (C++ executable): reads `.tdr` and can emit
  inventory JSON and neutral exports.
- `scripts/sentaurus_import.py` (Python script): parses surrounding Sentaurus
  text artifacts such as `.plt` and selected `.cmd` metadata.

Use both when building reference fixture chains under `reference_tcad/`.

## Direct command execution: vela device

For a preprocessed single-device SDevice `.cmd`, use the native `vela device`
entry point. It converts the grid, imports the explicitly declared parameter
file, validates/ emits stage decks and orchestrates the solver without an SDE
input. `--check` stops before solving. TDR input requires the detected HDF5
importer; the DF-ISE `.grd`/`.dat` path is also supported.

The [current workflow](sentaurus_sdevice_workflow.md) contains the supported
subset, gates, exit codes, output units and an executable
[neutral mesh/sidecar example](sentaurus_sdevice_workflow.md#210-可执行中性示例与-sidecar).
Use its section 2.10 commands from the worktree root. Help can be read with:

```powershell
.\build-release\vela.exe device -h
.\build-release\sentaurus_import.exe --help
```

```bash
./build-linux26/vela device -h
./build-linux26/sentaurus_import --help
```

The [runtime packages](sdevice_packaging.md) target Windows and Ubuntu 26.04,
and include Python, libraries and notices. The C++ importer exports neutral
reference fixtures; it does not execute an SDevice deck. The older Python
`scripts/sentaurus_import.py device --sde ... --sdevice ...` is a separate
conversion workflow with its own substitution policy described below.
Its examples containing `path/to/...` and the `FILE` CLI synopsis are templates
to replace with real inputs, not complete runnable examples.

For the native frontend, `--allow-lossy` accepts only material-matrix parameters
classified as `approximated`; it does not bypass unsupported physics, thermal
or mixed-mode gates. Numerical qualification and its limits are recorded in
[WP10](validation/sdevice_wp10_closeout_2026-10-03.md), with current delivery/DoD
evidence in the [native package follow-up](validation/sdevice_native_package_e2e_2026-10-05.md).

## Build Requirements

Sentaurus import support is controlled by `VELA_ENABLE_HDF5`:

- `VELA_ENABLE_HDF5=ON` (default in shipped presets) requests HDF5 support.
- If CMake finds an HDF5 target, it builds:
  - `vela_sentaurus` library
  - `sentaurus_import` executable
  - HDF5-linked tests (`sentaurus_tdr`, `sentaurus_sample_integration`)
- If HDF5 is missing, core solver targets still build, but Sentaurus TDR/HDF5
  targets are skipped.

Typical dependencies:

- Ubuntu/Debian: `libhdf5-dev`
- MSYS2 UCRT64: `mingw-w64-ucrt-x86_64-hdf5`

## CLI: sentaurus_import

Basic usage:

```bash
sentaurus_import --tdr FILE [--inventory-json FILE] [--export-dir DIR] \
  [--compensated-doping-policy reported|dominant_signed_region] \
  [--coordinate-unit um|cm]
```

Arguments:

- `--tdr`: input Sentaurus TDR file (required)
- `--inventory-json`: write parsed inventory metadata as JSON
- `--export-dir`: write neutral fixture exports for conversion/comparison flows
- `--coordinate-unit`: explicitly override the source geometry unit with `um`
  or `cm`. If omitted, use the TDR geometry unit metadata. Neutral coordinates
  are exported in micrometres in either case.
- `--compensated-doping-policy`: policy for compensated doping handling;
  accepted values:
  - `reported` (default)
  - `dominant_signed_region`

Examples:

```bash
build/sentaurus_import --tdr path/to/sample.tdr
build/sentaurus_import --tdr path/to/sample.tdr --inventory-json build/sample_inventory.json
build/sentaurus_import --tdr path/to/sample.tdr --export-dir reference_tcad/sample
```

## Python Text Import Helpers

Use `scripts/sentaurus_import.py` to process Sentaurus text outputs that are not
part of the TDR binary path (for example `.plt` curves and selected `.cmd`
metadata summaries).

Example pattern:

```bash
python scripts/sentaurus_import.py --help
```

### Explicit lossy impact-ionization substitution

The `device` command fails closed when an SDevice deck selects an impact-
ionization coefficient model that Vela does not implement.  `--allow-lossy`
permits only the coefficient models covered by the importer's explicit
substitution policy and records both the unsupported source model and the
chosen Vela model:

```bash
python scripts/sentaurus_import.py device \
  --sde path/to/device_sde.cmd \
  --sdevice path/to/device_sdevice.cmd \
  --output-dir build/imported_device \
  --allow-lossy
```

This option is not a general physics-whitelist override. Models without a
downstream substitution, including hydrodynamic and carrier-temperature
transport, remain fatal. `execution_ir.json` retains the unsupported source
model and `run_manifest.json` records the substitution warning.

## End-to-End Workflow

1. Build with HDF5 enabled and verify `sentaurus_import` exists.
2. Export TDR inventory and/or neutral files with `sentaurus_import`.
3. Use `scripts/convert_tcad_export.py` to generate Vela `mesh.json` and
   `simulation_*.json` from neutral exports.
4. Use `scripts/compare_reference_curves.py` to compare Vela outputs against
   reference curves.
5. Run tests:

```bash
ctest --test-dir build --output-on-failure -R sentaurus
ctest --test-dir build --output-on-failure -R reference_tcad_regression
```

## Related Documents

- `README.md`
- `reference_tcad/README.md`
- `tests/regression/README.md`
