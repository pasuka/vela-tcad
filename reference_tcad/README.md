# Reference TCAD CSV Fixtures

`reference_tcad/` contains neutral, text-based reference exports used to
cross-check Vela `unit_scaling` decks. The directory intentionally uses generic
CSV formats and avoids proprietary binary formats or commercial-tool-specific
names.

These fixtures validate signs, trends, finite outputs, and rough orders of
magnitude. They do not make calibration claims.

## Directory Shape

Each device directory can contain:

- `nodes.csv`: `id,x_um,y_um`
- `elements.csv`: `id,node0,node1,node2,region,material`
- `contacts.csv`: `name,node_ids,region`
- `doping.csv`: `node_id,donors_cm3,acceptors_cm3`
- `reference_curves/*.csv`: neutral reference curve summaries
- `vela/mesh.json`: converted Vela mesh
- `vela/simulation_*.json`: Vela decks using `scaling.mode = "unit_scaling"`
- `vela/*.csv`: checked-in Vela candidate curve outputs
- `reports/*.json` and `reports/*.md`: generated comparison reports
- `<device>_reference.json`: metadata inventory for checked-in Vela/reference
  curve fixtures

Only the README files and validation notes are hand-written. Files under
`reports/` are generated comparison outputs and should be regenerated through
the tool workflow instead of manually edited.

Checked-in fixture inventories use schema
`vela.reference_tcad.checked_in.v1`. These configs list the mesh, Vela decks,
candidate CSVs, reference curves, comparison reports, and curve kind for each
reusable sample. They are metadata-only for checked-in CSV fixtures; generated
Sentaurus imports use `vela.reference_tcad.sentaurus_reference.v1`.

## Tools

Sentaurus import workflow (HDF5/TDR + text artifacts):

1. Use the C++ `sentaurus_import` executable to read a `.tdr` file and export
  neutral mesh/doping/contact CSV files.
2. Use `scripts/sentaurus_import.py` for text artifacts such as `.plt` curve
  extraction and `.cmd`-derived summaries.
3. Convert neutral exports to Vela decks with `scripts/convert_tcad_export.py`.
4. Compare candidate and reference curves with `scripts/compare_reference_curves.py`.

Generate inventory JSON and neutral exports from TDR:

```bash
build/sentaurus_import --tdr path/to/device.tdr --inventory-json build/device_inventory.json
build/sentaurus_import --tdr path/to/device.tdr --export-dir reference_tcad/sample
```

Add `--compensated-doping-policy dominant_signed_region` when the default
`reported` policy is not desired for compensated regions.

Convert neutral CSV exports into Vela decks:

```bash
python scripts/convert_tcad_export.py \
  --input-dir reference_tcad/pn_diode \
  --output-dir reference_tcad/pn_diode/vela \
  --device pn_diode \
  --simulation-types iv,cv,bv
```

Compare reference and candidate curves:

```bash
python scripts/compare_reference_curves.py \
  --reference reference_tcad/pn_diode/reference_curves/pn_diode_reference_summary.csv \
  --candidate reference_tcad/pn_diode/reports/pn_diode_vela_summary.csv \
  --output-json reference_tcad/pn_diode/reports/pn_diode_comparison.json \
  --output-md reference_tcad/pn_diode/reports/pn_diode_comparison.md
```

Run tool and fixture checks:

```bash
ctest --test-dir build --output-on-failure -R reference_tcad_regression
ctest --test-dir build --output-on-failure -R sentaurus
```

## Checked-In Validation Chains

- `pn_diode`: 2D abrupt silicon PN diode with forward IV, reverse
  quasi-static CV, and reverse BV/max-field diagnostics.
- `nmos2d`: mixed Si/SiO2 NMOS prototype with interface charge, surface
  mobility, Id-Vd, Id-Vg, multi-terminal CV, and BV diagnostics.
- `pmos2d`: mixed Si/SiO2 PMOS prototype with the complementary polarity
  checks.
- `ldmos2d`: mixed Si/SiO2 LDMOS-like engineering trend validation for
  low-bias DD-IV, BV diagnostics, and field-plate max-field comparison.
- `igbt2d`: IGBT-like engineering trend validation for low-current IV,
  high-injection IV, stored charge proxy CV, BV diagnostics, and
  impact-ionization smoke diagnostics.

See `docs/validation/` for the hand-written validation summaries.
