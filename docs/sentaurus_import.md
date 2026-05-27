# Sentaurus TDR/HDF5 Import

This document describes the optional Sentaurus import path used to convert
Sentaurus TDR/HDF5 and related text outputs into Vela-friendly neutral fixtures.

## Scope

The repository provides two complementary import paths:

- `sentaurus_import` (C++ executable): reads `.tdr` and can emit
  inventory JSON and neutral exports.
- `scripts/sentaurus_import.py` (Python script): parses surrounding Sentaurus
  text artifacts such as `.plt` and selected `.cmd` metadata.

Use both when building reference fixture chains under `reference_tcad/`.

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
  [--compensated-doping-policy reported|dominant_signed_region]
```

Arguments:

- `--tdr`: input Sentaurus TDR file (required)
- `--inventory-json`: write parsed inventory metadata as JSON
- `--export-dir`: write neutral fixture exports for conversion/comparison flows
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
- `docs/examples.md`
- `reference_tcad/README.md`
- `tests/regression/README.md`
