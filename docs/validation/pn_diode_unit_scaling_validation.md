# PN Diode Unit Scaling Validation

Date: 2026-05-20

This page records the first checked-in PN diode reference_tcad validation chain
for Vela `unit_scaling`. It is a trend and order-of-magnitude validation only;
there is no calibration claim.

## Fixture

- Source export: `reference_tcad/pn_diode/nodes.csv`,
  `elements.csv`, `contacts.csv`, and `doping.csv`.
- Structure: 2D silicon abrupt junction with `p_region`, `n_region`,
  `anode`, and `cathode`.
- Decks: `reference_tcad/pn_diode/vela/simulation_iv.json`,
  `simulation_cv.json`, and `simulation_bv.json`.
- Scaling: every generated deck uses `"scaling": {"mode": "unit_scaling"}`.

## Commands

```powershell
python scripts/convert_tcad_export.py --input-dir reference_tcad/pn_diode --output-dir reference_tcad/pn_diode/vela --device pn_diode --simulation-types iv,cv,bv
build\vela_example_runner.exe --config reference_tcad\pn_diode\vela\simulation_iv.json
build\vela_example_runner.exe --config reference_tcad\pn_diode\vela\simulation_cv.json
build\vela_example_runner.exe --config reference_tcad\pn_diode\vela\simulation_bv.json
python scripts/compare_reference_curves.py --reference reference_tcad/pn_diode/reference_curves/pn_diode_reference_summary.csv --candidate reference_tcad/pn_diode/reports/pn_diode_vela_summary.csv --output-json reference_tcad/pn_diode/reports/pn_diode_comparison.json --output-md reference_tcad/pn_diode/reports/pn_diode_comparison.md
```

## Results

| Check | Result |
| --- | --- |
| IV monotonic | Pass: forward current increases from `2.639e-18` to `3.176e-6` to `1.718e-2` A/m. |
| finite capacitance | Pass: reverse CV has finite nonzero points, including `3.818e-15` F/m and `3.049e-17` F/m. |
| max field non-decreasing | Pass: BV diagnostic rises from `3491.77` to `5991.93` to `8491.93` V/cm. |
| Reference trend match | Pass: IV, CV, and BV trend fields match in `pn_diode_comparison.json`. |

The comparison uses explicit CSV/text reference data and checked-in Vela outputs.
It validates signs, trends, and key orders of magnitude for this fixture.
