# Reference TCAD CSV Fixtures

`reference_tcad/` contains neutral, text-based reference exports used to
cross-check Vela `unit_scaling` decks. Reusable reference cases may also retain
text-only commercial-tool input decks when they are required to reproduce a
published comparison. Proprietary binary outputs remain generated artifacts
and are not checked in.

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
- `source/*.cmd` and `source/*.par`: optional text-only source decks used to
  regenerate a commercial-tool reference

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

- `bvmethods_sentaurus2018`: Sentaurus Training NMOS BV method inputs covering
  ABA, external resistor, voltage-to-current, continuation, and transient
  approaches, plus the corresponding supported Vela template mapping.
- `schottky_charon_sentaurus2018`: Charon-derived 2-D n-silicon Schottky
  diode translated to Sentaurus O-2018.06-SP2, with a full 0--1 V reference
  curve and a full 0--1 V two-stage Vela thermionic-Robin acceptance chain
  using voltage continuation followed by pseudo-arclength continuation.
- `transportmodels_sentaurus2022`: Sentaurus T-2022.03-SP2 50 nm NMOS
  TransportModels inputs for matched DD/electron-DG Id-Vg and Id-Vd runs,
  including SDE/SDevice sources, neutral mesh/doping inputs, frozen Vela
  contracts, and the strict 12-stage continuous-scan configurations.
- `genius_bjt_sentaurus2022`: the Genius TCAD two-dimensional NPN BJT rebuilt
  with Sentaurus SDE/SDevice T-2022.03-SP2 and imported on common meshes into
  Vela.  The reusable fixture covers 31-point terminal curves, common-mesh
  potential/carrier fields, conservative section fluxes, SRH/Auger sources,
  local mesh sensitivity, and diagnostic current-vector recovery.  See its
  `CASE_SUMMARY.md` for the current acceptance boundary.

See `docs/validation/` for the hand-written validation summaries.
