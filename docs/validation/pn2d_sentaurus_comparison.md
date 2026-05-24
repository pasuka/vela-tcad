# pn2d Sentaurus Comparison

The `reference_tcad/pn2d` case is a 2-D abrupt PN diode imported from
Sentaurus. The structure is `L=2.0 um`, `H=0.5 um`, with junction position
`Xj=1.0 um`, `Na=Nd=1e17 cm^-3`, an `Anode` contact on the left boundary, and
a `Cathode` contact on the right boundary.

The import path preserves the Sentaurus TDR mesh, contacts, material regions,
reference IV/BV curves, command summaries, and full `doping.csv` node-level
donor/acceptor table. Generated faithful Vela decks keep `node_doping_file:
"doping.csv"` so the exact imported doping is available for solver
development and regression inspection.

The current executable comparison uses explicit runtime approximation decks:

- IV: `simulation_iv_runtime.json` uses the imported TDR mesh and Sentaurus
  sweep definition, but runs with region-average doping scaled by `1e-4` for
  Vela convergence diagnostics. The IV report is trend-gated with at least 10
  compared points.
- BV: `simulation_bv_runtime.json` uses the imported TDR mesh, Vela's
  Selberherr impact-ionization diagnostic in place of Sentaurus
  `Avalanche(OkutoCrowell)`, and the same `1e-4` region-average doping runtime
  scale. The BV report is diagnostic-only and does not require trend match.

Vela currently treats Sentaurus Fermi statistics as Boltzmann carrier
statistics. Sentaurus Okuto-Crowell avalanche is approximated by Vela
Selberherr impact ionization for BV diagnostics. These reports establish a
readable, runnable, finite, and comparable regression loop; they are not yet a
calibrated numerical match to Sentaurus.

## Hybrid Solver Status

Faithful pn2d decks now use `solver.method: "gummel_newton"` and preserve
`node_doping_file: "doping.csv"`. Each faithful bias point is configured to run
Gummel first and hand the validated state to coupled Newton with
`warm_start=true`.

The runtime approximation remains available as the executable comparison path:
it uses imported TDR geometry, region-average doping scaled by
`runtime_doping_scale`, and conservative `solver.method: "gummel"`. Current
gate:

- faithful IV/BV deck generation is required;
- faithful decks must preserve node-level doping and strict Newton handoff
  settings;
- runtime IV/BV execution must remain finite;
- strict Sentaurus numerical agreement is not yet required.

Useful local verification command:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\sentaurus_import.py reference --config reference_tcad\pn2d\pn2d_reference.json --source-dir reference_tcad\pn2d --output-dir build\reference_tcad\pn2d --tdr-importer build\sentaurus_import.exe --runner build\vela_example_runner.exe
```
