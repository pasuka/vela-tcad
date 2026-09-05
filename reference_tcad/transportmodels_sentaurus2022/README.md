# TransportModels Sentaurus 2022 DD/DG reference inputs

This directory freezes the text inputs used by the Vela versus Sentaurus
TransportModels MOS comparison.  The source case is the Sentaurus Applications
Library `GettingStarted/sdevice/TransportModels` example, executed with
Sentaurus T-2022.03-SP2.

Only source and neutral text inputs are checked in.  Sentaurus TDR/PLT files,
Vela state files, logs, plots, and other generated results remain under the
ignored `build*/reference_tcad/` trees.

## Directory layout

- `source/nMOS_dvs.cmd`: SDE geometry, contacts, doping, and mesh refinement.
- `source/n1_half_msh.cmd`: generated Sentaurus Mesh command retained to audit
  the exact mesh contract.
- `source/IdVgs_des.cmd` and `source/IdVds_des.cmd`: original Workbench SDevice
  templates.
- `source/pp6_des.*` and `source/pp7_des.*`: preprocessed DD and electron-DG
  Id-Vg decks.
- `source/pp12_des.*` and `source/pp13_des.*`: preprocessed DD and electron-DG
  Id-Vd decks.
- `source/Silicon.par` and `source/sdevice.par`: material parameter include
  chain retained from the Application Library case.
- `vela/mesh.json` and `vela/doping.csv`: neutral text inputs imported from
  `n1_msh.tdr`; no proprietary binary mesh is stored here.
- `vela/materials_sentaurus2022.json`: frozen Vela material model.
- `vela/contracts/`: frozen physical, continuous-scan, contact-basin, and
  material contracts.
- `vela/configs/00_*.json` through `11_*.json`: portable snapshots of the exact
  12-stage DD/DG Id-Vg and Id-Vd workflow used on 2026-08-24.

The Vela configuration paths are relative to each file under `vela/configs/`,
matching the runner's path-resolution rules.  Generated outputs are directed
to `build/reference_tcad/transportmodels_sentaurus2022/work` so running the
fixture does not dirty this checked-in directory.

## Sentaurus workflow

The SDE deck retains Sentaurus Workbench `@node@` placeholders.  Run it inside
Workbench, or preprocess the node placeholder consistently before standalone
execution.  It creates the half structure and reflects it with `tdx` to produce
`n1_msh.tdr`.

With the mesh and matching parameter file in one directory, the four
preprocessed reference simulations are:

```text
sdevice pp6_des.cmd   # DD Id-Vg, Vd=1.1 V, Vg=-1.0..2.2 V
sdevice pp7_des.cmd   # DG Id-Vg, Vd=1.1 V, Vg=-1.0..2.2 V
sdevice pp12_des.cmd  # DD Id-Vd, Vg=1.0 V, Vd=0.0..2.0 V
sdevice pp13_des.cmd  # DG Id-Vd, Vg=1.0 V, Vd=0.0..2.0 V
```

The DD/DG distinction is limited to `eQuantumPotential` and the corresponding
coupled equation set.  The classical mobility, Fermi statistics,
OldSlotboom BGN, and doping/temperature-dependent SRH definitions are shared.

## Vela workflow

Run from the repository root after building `vela_example_runner`:

```powershell
New-Item -ItemType Directory -Force `
  build/reference_tcad/transportmodels_sentaurus2022/work | Out-Null

$configs = @(
  "00_dd_idvg_equilibrium.json",
  "01_dd_idvg_drain_ramp.json",
  "02_dd_idvg_final_bias_relax.json",
  "03_dd_idvg_curve.json",
  "04_dd_idvd_equilibrium.json",
  "05_dd_idvd_curve.json",
  "06_dg_idvg_equilibrium.json",
  "07_dg_idvg_drain_ramp.json",
  "08_dg_idvg_final_bias_relax.json",
  "09_dg_idvg_curve.json",
  "10_dg_idvd_equilibrium.json",
  "11_dg_idvd_curve.json"
)

foreach ($config in $configs) {
  build-release/vela_example_runner.exe --config `
    "reference_tcad/transportmodels_sentaurus2022/vela/configs/$config"
  if ($LASTEXITCODE -ne 0) { throw "TransportModels stage failed: $config" }
}
```

Stages must remain in the listed order because each continuation stage consumes
only the immediately preceding accepted state.  The curve configurations also
enforce the deep-off terminal criterion `Id/abs(KCL) >= 10`.

Canonical orchestration and comparison scripts remain in the repository-level
`scripts/` directory:

- `run_transportmodels_dd_dg_continuous_baseline.py`
- `run_transportmodels_sentaurus_three_regime_spatial_oracles.py`
- `run_transportmodels_vela_three_regime_spatial_states.py`
- `compare_transportmodels_three_regime_spatial_fields.py`

See `transportmodels_sentaurus2022_reference.json` for the inventory and
SHA-256 fingerprints.
