# Genius NPN BJT comparison figures

These figures use the accepted 5611-node, 10940-triangle common mesh and the
corrected M1 Vela run. Current is normalized to A/um. The spatial comparison is
at VBE=0.70 V and VCE=3.00 V; carrier densities are shown in cm^-3.

## Figures

- `genius_bjt_device_mesh.png`: full-device and top-junction mesh views over
  the signed net-doping field, with all three contacts identified.
- `genius_bjt_m1_vce3_field_comparison.png`: common-scale Sentaurus and Vela
  maps for potential, electron density, and hole density, plus Vela-minus-
  Sentaurus differences.
- `genius_bjt_terminal_curve_comparison.png`: M0/M1 Ic, Ib, and beta curves at
  all 31 requested collector biases. Current panels use magnitude and a
  logarithmic vertical scale.
- `genius_bjt_m1_parity_error.png`: M1 active-region absolute log10 errors and
  the pre-registered 0.05-decade gate.

`figure_manifest.json` records source and output hashes, mesh counts, coordinate
alignment, and spatial error metrics. The Sentaurus and Vela node coordinates
match exactly in the plotted data.

## Reproduction

From the repository root with the MSYS2 UCRT64 Python environment:

```powershell
D:\msys64\ucrt64\bin\python.exe scripts\plot_genius_bjt_sentaurus_vela.py
```

The field plot requires the ignored TDR export under
`build-release/reference_tcad/genius_bjt_sentaurus2022/m1_current_diagnosis/sentaurus_vce3`
and the corrected Vela state under
`build-release/reference_tcad/genius_bjt_sentaurus2022/vela_wp3_wp5`.
