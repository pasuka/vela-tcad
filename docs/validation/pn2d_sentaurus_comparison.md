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

The current executable comparison uses the faithful node-level doping decks.
The IV deck uses `vela_stop: 0.3` and `vela_step: 0.1`; the local gate compares
the 0.2-0.3 V forward-bias window against the full imported Sentaurus IV
reference. The BV deck uses `vela_stop: 10.0` and `vela_step: 5.0`; Vela
currently disables the Sentaurus `Avalanche(OkutoCrowell)` approximation in
the strict handoff deck because the imported 0 V state is not yet robust with
impact-ionization Jacobian terms. These reports are diagnostic-only and do not
yet require trend match.

Vela currently treats Sentaurus Fermi statistics as Boltzmann carrier
statistics. Sentaurus Okuto-Crowell avalanche is approximated by Vela
Selberherr impact ionization for BV diagnostics. These reports establish a
readable, runnable, finite, and comparable regression loop; they are not yet a
calibrated numerical match to Sentaurus.

## Hybrid Solver Status

Faithful pn2d decks now use `solver.method: "gummel_newton"` with
`handoff.fallback: "none"` and preserve `node_doping_file: "doping.csv"`. Each
accepted faithful IV/BV row must end with `handoff_stage: "newton"` and
`newton_iterations > 0`; Gummel fallback is no longer part of the default pn2d
gate.

Current solver limits are intentionally explicit: Gummel is used as a one-step
initializer (`handoff.gummel_max_iter: 1`) and coupled Newton owns the accepted
state. This is a strict provenance gate, not a calibrated convergence claim:
the IV comparison currently allows a wide order-of-magnitude envelope while
the coupled Newton residual and physical current calibration are improved.

The old region-average `runtime_doping_scale` path is no longer required for
pn2d. It remains available through an opt-in `runtime_diagnostic` config block
for future debugging, but the default bundled pn2d reference runs the imported
node-level doping. Current gate:

- faithful IV/BV deck generation is required;
- faithful decks must preserve node-level doping and hybrid handoff settings;
- faithful IV/BV execution must remain finite and end in Newton handoff;
- comparison reports align by `bias_V` and use configured bias windows;
- strict Sentaurus numerical agreement is not yet required.

Useful local verification command:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\sentaurus_import.py reference --config reference_tcad\pn2d\pn2d_reference.json --source-dir reference_tcad\pn2d --output-dir build\reference_tcad\pn2d --tdr-importer build\sentaurus_import.exe --runner build\vela_example_runner.exe
```
