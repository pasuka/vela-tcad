# BV Impact-Ionization Validation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a defensible PN2D BV validation workflow that refreshes Sentaurus reference artifacts, confirms the active avalanche model, documents Vela's impact-ionization implementation, and promotes only trend/methodology gates that match the current evidence.

**Architecture:** Treat `reference_tcad/pn2d_sentaurus2018/source/pn2d_bv_sdevice.cmd` as the Sentaurus source of truth and keep VM regeneration, import, Vela execution, and comparison as separate stages. Confirm `Avalanche(VanOverstraeten)` before any model work; do not implement Okuto-Crowell unless a freshly generated deck proves it is active. Gate the workflow on reproducibility, model identity, monotonic/trend behavior, and documented current-band windows rather than full -20 V absolute-current parity.

**Tech Stack:** C++20, CMake/Ninja, Catch2, Python standard library diagnostics, MSYS2 UCRT64 on Windows, existing `sentaurus_import`, `vela_example_runner`, `DCSweep`, `ImpactIonizationModel`, and PN2D Sentaurus2018 fixtures.

---

## Current Evidence To Preserve

- The checked-in Sentaurus BV deck uses `Recombination(SRH Avalanche(VanOverstraeten))` and targets `Goal { Name="Anode" Voltage=-20.0 }`.
- The BV deck writes `pn2d_bv_multibias` snapshots with `Intervals=200`; this produces 201 endpoint-inclusive TDR files, `pn2d_bv_multibias_0000_des.tdr` through `pn2d_bv_multibias_0200_des.tdr`.
- Vela already has `van_overstraeten` in `include/vela/physics/ImpactIonizationModel.h` and `src/physics/ImpactIonizationModel.cpp`.
- The committed PN2D Sentaurus2018 BV generated-reference config is still a low-reverse-bias smoke gate with `vela_stop: -0.05`; any -20 V Vela run in this plan must use an explicit derived deck or a separately reviewed promotion step.
- Prior high-bias evidence is mixed by window: local `-13.0..-13.2 V` current can be close enough for a bounded-band diagnostic, while existing knee-shape evidence still shows meaningful high-bias curve-shape mismatch. Do not collapse these into a single "BV parity passes" claim.

## File Responsibilities

- Create `docs/validation/bv_impact_ionization_theory.md`: theory and code mapping for Chynoweth, van Overstraeten-de Man, Okuto-Crowell as contrast, ionization integral, multiplication, driving force choices, and SG edge-current source.
- Create or update `docs/validation/pn2d_bv_validation.md`: method, artifacts, accepted gates, non-goals, and interpretation of known residuals.
- Modify `scripts/run_sentaurus_vm_reference.py` only if live-run ergonomics need additional explicit phases; otherwise use its current `pn2d --stages ... --dry-run` and live-run CLI.
- Modify `tests/regression/test_sentaurus_vm_reference_runner.py` if the VM runner behavior or manifest schema changes.
- Modify `scripts/sentaurus_import.py` and `tests/regression/test_sentaurus_import_tools.py` only if model provenance is not already captured strongly enough in generated decks/manifests.
- Modify `scripts/compare_pn2d_bv_multibias_fields.py` and `tests/regression/test_reference_tcad_tools.py` only if a new trend-summary or same-dimension avalanche-source metric is needed.
- Modify `scripts/run_regression.py` and `tests/regression/test_run_regression.py` only for a lightweight BV trend check that does not require Sentaurus or VM access.
- Modify `tests/test_dc_sweep.cpp` only for generic BV output semantics or current/magnitude-band helper behavior, not for VM/Sentaurus artifact availability.

---

### Task 1: Document The Correct BV Scope And Stale Assumptions

**Files:**
- Create: `docs/validation/bv_impact_ionization_theory.md`
- Modify: `docs/validation/pn2d_bv_validation.md` if it already exists; otherwise create it in Task 10.
- Reference only: `reference_tcad/pn2d_sentaurus2018/source/pn2d_bv_sdevice.cmd`
- Reference only: `include/vela/physics/ImpactIonizationModel.h`
- Reference only: `src/physics/ImpactIonizationModel.cpp`
- Reference only: `include/vela/equation/AssemblerUtils.h`

- [ ] **Step 1: Verify source facts before writing**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
rg -n "Avalanche|Goal|Intervals|AvalancheGeneration" reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_sdevice.cmd
rg -n "van_overstraeten|VanOverstraeten|electronCoefficient|holeCoefficient|switchField|phononEnergy" include\vela\physics src\physics tests\test_impact_ionization.cpp
rg -n "sgEdgeCurrentAvalancheSourceRecords|electronAlpha|holeAlpha|edgeSourceIntegral" include\vela\equation\AssemblerUtils.h src\simulation\DCSweep.cpp
```

Expected:

```text
reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_sdevice.cmd contains Avalanche(VanOverstraeten), Goal ... Voltage=-20.0, and Intervals=200.
Vela impact-ionization sources contain van_overstraeten and SG edge-current avalanche source records.
```

- [ ] **Step 2: Write the theory document**

Create `docs/validation/bv_impact_ionization_theory.md` with these exact sections:

````markdown
# BV Impact-Ionization Theory And Vela Mapping

## Scope

This document describes the PN2D Sentaurus2018 BV validation target. The checked-in BV deck uses `Avalanche(VanOverstraeten)`, not Okuto-Crowell, and sweeps the Anode to `-20.0 V`. Okuto-Crowell is included only as contrast unless a fresh Sentaurus run proves a different active model.

## Chynoweth Form

The common local ionization coefficient form is:

```math
\alpha(E) = A \exp(-B / |E|)
```

`A` is a prefactor in inverse length and `B` is a critical field. Vela evaluates coefficient functions through `ImpactIonizationModel::electronCoefficient` and `ImpactIonizationModel::holeCoefficient`.

## Van Overstraeten-de Man

Vela's `VanOverstraetenImpactIonization` uses low-field and high-field coefficient sets selected by `switchField`, plus a temperature factor:

```math
\gamma(T) =
\frac{\tanh(\hbar\omega / (2 k_B T_\mathrm{ref}))}
{\tanh(\hbar\omega / (2 k_B T))}
```

```math
\alpha(E, T) = \gamma(T) A_\mathrm{region} \exp(-\gamma(T) B_\mathrm{region} / |E|)
```

The code mapping is `ImpactIonizationModelConfig::{electronALow,electronAHigh,electronBLow,electronBHigh,holeALow,holeAHigh,holeBLow,holeBHigh,switchField,phononEnergy,referenceTemperature_K,temperature_K}`.

## Okuto-Crowell Contrast

Okuto-Crowell is commonly written as:

```math
\alpha(E) = a E^2 \exp(-(b/E)^2)
```

This is not the current PN2D BV fixture target. Do not implement `OkutoCrowellImpactIonization` as part of this validation unless the freshly run deck or exported parameter provenance contradicts the checked-in `Avalanche(VanOverstraeten)` source.

## Driving Force

Sentaurus isothermal avalanche defaults are interpreted as quasi-Fermi-gradient driven unless the deck enables a different option. Vela maps this through `impact_ionization.driving_force = "quasi_fermi_gradient"` and `current_approximation = "density_gradient"` for the Sentaurus-default SG edge-current avalanche path.

## SG Edge-Current Source

The same-dimension comparison target is not `alpha(E)` versus `AvalancheGeneration`. Compare coefficient to coefficient, or compare generation/source-integral to generation/source-integral. Vela's SG path records `electronAlpha`, `holeAlpha`, electron/hole flux proxies, and edge/node source integrals through `sgEdgeCurrentAvalancheSourceRecords`.

## Ionization Integral And Multiplication

The one-dimensional intuition is:

```math
M = \frac{1}{1 - \int \alpha \, dl}
```

This is useful as a diagnostic for field-line breakdown propensity. It is not the production acceptance path for this PN2D Sentaurus-default SG edge-current validation unless a later task explicitly adds a post-processing criterion.

## Current Validation Boundary

This validation promotes reproducibility, model identity, documented field/current trends, and windowed current diagnostics. It does not claim full absolute-current parity over the entire `0..-20 V` sweep and does not promote hidden scalar calibration knobs such as `source_geometry_scale`.
````

- [ ] **Step 3: Check the document for stale claims**

Run:

```powershell
rg -n "Okuto|Selberherr|source_geometry_scale|absolute-current parity|200 multibias|201" docs\validation\bv_impact_ionization_theory.md
```

Expected:

```text
Okuto appears only as contrast or a conditional future task.
No wording claims full absolute-current parity.
No wording says there are only 200 multibias TDR files.
```

### Task 2: Correct And Use The Sentaurus VM Runner Workflow

**Files:**
- Reference first: `scripts/run_sentaurus_vm_reference.py`
- Test: `tests/regression/test_sentaurus_vm_reference_runner.py`
- Optional modify: `scripts/run_sentaurus_vm_reference.py`
- Optional modify: `tests/regression/test_sentaurus_vm_reference_runner.py`

- [ ] **Step 1: Run the current dry-run manifest command**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts\run_sentaurus_vm_reference.py pn2d `
  --ssh-target sentaurus `
  --source-dir reference_tcad\pn2d_sentaurus2018\source `
  --local-output-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs `
  --remote-root "~/sentaurus_runs/vela_oracle" `
  --run-id pn2d_bv_validation_plan_dry_run `
  --stages 0v,iv,bv `
  --dry-run
```

Expected:

```text
The command exits 0 and writes build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_validation_plan_dry_run\sentaurus_vm_run_manifest.json.
The manifest contains stages ["0v", "iv", "bv"] and commands for sde plus all three sdevice decks.
```

- [ ] **Step 2: Decide whether to add explicit phase flags**

Keep the current CLI unless the dry-run/live ergonomics are insufficient. The current supported workflow is:

```text
dry-run planning: include --dry-run
live upload/run/fetch: omit --dry-run
```

If explicit `--plan`, `--upload`, `--launch`, or `--fetch` flags are added, update `tests/regression/test_sentaurus_vm_reference_runner.py` in the same task so both old and new behavior are unambiguous. Do not write documentation that mentions flags the script does not support.

- [ ] **Step 3: Run VM runner regression tests**

Run:

```powershell
python -m unittest tests.regression.test_sentaurus_vm_reference_runner -v
```

Expected:

```text
test_dry_run_writes_manifest_without_ssh ... ok
test_missing_required_deck_fails_before_ssh ... ok
```

### Task 3: Add A Fetched-Artifact Validator

**Files:**
- Create: `scripts/validate_pn2d_sentaurus_bv_artifacts.py`
- Test: `tests/regression/test_reference_tcad_tools.py`

- [ ] **Step 1: Write regression tests for artifact validation**

Add tests that create a temporary `source` directory with:

```text
pn2d_bv.plt
pn2d_bv.log
pn2d_bv_multibias_0000_des.tdr
pn2d_bv_multibias_0001_des.tdr
...
pn2d_bv_multibias_0200_des.tdr
```

The positive test must expect success with 201 TDR files. The negative test must remove `pn2d_bv_multibias_0200_des.tdr` and expect a non-zero exit plus text containing:

```text
expected 201 pn2d_bv_multibias TDR files
```

- [ ] **Step 2: Run the new tests and confirm they fail**

Run:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest -v
```

Expected before implementation:

```text
The new validator tests fail because scripts\validate_pn2d_sentaurus_bv_artifacts.py does not exist.
```

- [ ] **Step 3: Implement the validator**

Create `scripts/validate_pn2d_sentaurus_bv_artifacts.py` with behavior:

```text
arguments:
  --source-dir PATH
  --require-final-bias -20.0
  --expected-multibias-count 201

checks:
  pn2d_bv.plt exists
  a BV log artifact exists: `pn2d_bv.log`, `pn2d_bv.log_des.log`, or `run_pn2d_bv.out`
  no fatal/license/error marker appears in pn2d_bv.log, except harmless case-insensitive words inside paths are ignored only if tests cover them
  exactly 201 files match pn2d_bv_multibias_*_des.tdr by default
  first index is 0000 and last index is 0200
  if the PLT parser can identify Anode OuterVoltage, the minimum or final bias reaches -20.0 within 1e-6 V
```

Keep the parser permissive: if the PLT format cannot be parsed in a minimal test fixture, report a warning rather than failing solely on parser format.

- [ ] **Step 4: Run validator tests**

Run:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest -v
```

Expected:

```text
The new artifact validator tests pass.
```

### Task 4: Preserve Avalanche Model Provenance Through Import

**Files:**
- Reference first: `scripts/sentaurus_import.py`
- Test: `tests/regression/test_sentaurus_import_tools.py`
- Optional modify: `scripts/sentaurus_import.py`
- Optional modify: `tests/regression/test_sentaurus_import_tools.py`

- [ ] **Step 1: Inspect current model provenance behavior**

Run:

```powershell
python -m unittest tests.regression.test_sentaurus_import_tools -v
rg -n "VanOverstraeten|Okuto|runtime_approximation|impact_ionization|current_approximation" scripts\sentaurus_import.py tests\regression\test_sentaurus_import_tools.py reference_tcad\pn2d_sentaurus2018\pn2d_sentaurus2018_reference.json
```

Expected:

```text
Existing tests pass.
The PN2D Sentaurus2018 BV config maps to model "van_overstraeten", driving_force "quasi_fermi_gradient", generation "current_density", and current_approximation "density_gradient".
```

- [ ] **Step 2: Add or strengthen a provenance test only if needed**

If no test explicitly protects the BV model mapping, add one to `tests/regression/test_sentaurus_import_tools.py` that parses a miniature command block containing:

```text
Recombination(
  SRH
  Avalanche(VanOverstraeten)
)
```

Expected generated solver snippet:

```json
{
  "impact_ionization": {
    "model": "van_overstraeten",
    "driving_force": "quasi_fermi_gradient",
    "generation": "current_density",
    "current_approximation": "density_gradient"
  }
}
```

Also assert that no warning claims Okuto-Crowell approximation for this deck.

- [ ] **Step 3: Run import tests**

Run:

```powershell
python -m unittest tests.regression.test_sentaurus_import_tools -v
```

Expected:

```text
All sentaurus_import tool tests pass.
```

### Task 5: Define Same-Dimension Avalanche Comparison Outputs

**Files:**
- Reference first: `scripts/compare_pn2d_bv_multibias_fields.py`
- Reference first: `include/vela/equation/AssemblerUtils.h`
- Optional modify: `scripts/compare_pn2d_bv_multibias_fields.py`
- Optional modify: `tests/regression/test_reference_tcad_tools.py`
- Optional use: `scripts/compare_pn2d_bv_sg_edge_source_dump.py`

- [ ] **Step 1: Keep field comparison semantics explicit**

Do not describe the existing `avalanche_generation` field comparator as an alpha comparator. It compares Sentaurus `AvalancheGeneration` or `ImpactIonization` to Vela `AvalancheGeneration` with the existing unit scale.

Run:

```powershell
rg -n "avalanche_generation|AvalancheGeneration|ImpactIonization|electronAlpha|holeAlpha|edgeSourceIntegral" scripts tests include\vela\equation\AssemblerUtils.h
```

Expected:

```text
The codebase exposes both generation-field comparison and Vela SG edge alpha/source records.
```

- [ ] **Step 2: Add a summary metric only if the existing reports lack it**

If the existing report does not produce a compact trend summary, extend `scripts/compare_pn2d_bv_multibias_fields.py` to write `bv_trend_summary.json` containing:

```json
{
  "biases_compared": [-0.5, -2.0, -5.0, -10.0, -20.0],
  "max_field_monotonic": true,
  "current_windows": [
    {
      "name": "mid_bias_current_band",
      "bias_min": -13.2,
      "bias_max": -13.0,
      "min_ratio": 0.6,
      "max_ratio": 1.4,
      "status": "pass"
    }
  ],
  "high_bias_knee_shape": {
    "status": "diagnostic",
    "reason": "current evidence does not support full -20 V absolute-current parity"
  }
}
```

Use available curve CSVs for current-window calculations. If a requested window is absent from a curve, report `status: "not_evaluated"` instead of silently passing.

- [ ] **Step 3: Add tests for missing-window behavior**

In `tests/regression/test_reference_tcad_tools.py`, create a tiny reference/candidate curve pair that lacks `-13.0..-13.2 V`. The expected summary must contain:

```json
{
  "status": "not_evaluated"
}
```

This prevents accidental promotion of a current-band gate when the supporting bias points are absent.

- [ ] **Step 4: Run comparison tests**

Run:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools -v
```

Expected:

```text
All reference TCAD tool tests pass.
```

### Task 6: Run Fresh Sentaurus And Import Artifacts

**Files:**
- Generated only: `build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_vm_runs/<run-id>/`
- Generated only: `build-release/reference_tcad/pn2d_sentaurus2018/`
- Reference: `scripts/run_sentaurus_vm_reference.py`
- Reference: `scripts/validate_pn2d_sentaurus_bv_artifacts.py`
- Reference: `scripts/sentaurus_import.py`

- [ ] **Step 1: Build the `build-release` binaries used by Tasks 6-8**

Tasks 6, 7, and 8 invoke `build-release\sentaurus_import.exe` and `build-release\vela_example_runner.exe`. The plan must build that tree explicitly; do not assume it already exists, and do not silently fall back to the `build` (debug) tree, because the binary tree must be consistent across import, candidate run, and comparison.

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --preset windows-ucrt64-release
cmake --build build-release --parallel --target vela_example_runner sentaurus_import
```

Expected:

```text
build-release\vela_example_runner.exe and build-release\sentaurus_import.exe exist and are newly built.
If the windows-ucrt64-release preset is unavailable, switch ALL build-release paths in Tasks 6-8 to the build (debug) tree instead, so import, candidate run, and comparison share one tree.
```

- [ ] **Step 2: Confirm VM reachability**

Run only when the Sentaurus VM and license should be available:

```powershell
C:\Windows\System32\OpenSSH\ssh.exe sentaurus "hostname && which sdevice && which sde"
```

Expected:

```text
The SSH command exits 0 and prints paths for sdevice and sde.
```

- [ ] **Step 3: Launch the live VM run**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts\run_sentaurus_vm_reference.py pn2d `
  --ssh-target sentaurus `
  --source-dir reference_tcad\pn2d_sentaurus2018\source `
  --local-output-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs `
  --remote-root "~/sentaurus_runs/vela_oracle" `
  --run-id pn2d_bv_validation_refresh `
  --stages 0v,iv,bv
```

Expected:

```text
The command exits 0 and fetches artifacts under build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_validation_refresh\source.
Warnings are acceptable only for artifact globs that are genuinely absent and not required by the BV validation.
```

- [ ] **Step 4: Validate fetched BV artifacts**

Run:

```powershell
python scripts\validate_pn2d_sentaurus_bv_artifacts.py `
  --source-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_validation_refresh\source `
  --require-final-bias -20.0 `
  --expected-multibias-count 201
```

Expected:

```text
Validation exits 0.
The report confirms pn2d_bv.plt, a clean BV log artifact, and 201 endpoint-inclusive multibias TDR files.
```

- [ ] **Step 5: Import refreshed references**

Run:

```powershell
python scripts\sentaurus_import.py reference `
  --config reference_tcad\pn2d_sentaurus2018\pn2d_sentaurus2018_reference.json `
  --source-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_validation_refresh\source `
  --output-dir build-release\reference_tcad\pn2d_sentaurus2018 `
  --tdr-importer build-release\sentaurus_import.exe `
  --runner build-release\vela_example_runner.exe
```

Expected:

```text
The import exits 0.
Reference curves and generated Vela decks are refreshed under build-release\reference_tcad\pn2d_sentaurus2018.
The committed source fixture is not overwritten by this task.
```

### Task 7: Define The Vela -20 V Candidate Deck Explicitly

**Files:**
- Generated only: `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_validation_candidate/simulation_bv_minus20_validation.json`
- Generated only: `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_validation_candidate/`
- Reference: `reference_tcad/pn2d_sentaurus2018/pn2d_sentaurus2018_reference.json`
- Reference: `build-release/reference_tcad/pn2d_sentaurus2018/vela/simulation_bv.json`

- [ ] **Step 1: Inspect generated BV deck before deriving**

Run:

```powershell
@'
import json
from pathlib import Path
path = Path("build-release/reference_tcad/pn2d_sentaurus2018/vela/simulation_bv.json")
data = json.loads(path.read_text())
print(json.dumps({
    "stop": data.get("sweep", {}).get("stop"),
    "step": data.get("sweep", {}).get("step"),
    "impact_ionization": data.get("solver", {}).get("impact_ionization"),
    "recombination": data.get("solver", {}).get("recombination"),
}, indent=2))
'@ | python -
```

Expected:

```text
The generated deck uses van_overstraeten/current_density/density_gradient.
If the generated deck still stops at -0.05 V, continue with a derived candidate deck rather than editing the committed reference JSON in this task.
```

- [ ] **Step 2: Create a derived -20 V candidate deck**

Use a small Python one-off or a checked-in helper only if this process needs to be repeated. The derived deck must:

```json
{
  "sweep": {
    "stop": -20.0,
    "step": -0.05
  },
  "solver": {
    "impact_ionization": {
      "model": "van_overstraeten",
      "driving_force": "quasi_fermi_gradient",
      "generation": "current_density",
      "current_approximation": "density_gradient"
    },
    "recombination": ["srh"]
  }
}
```

Do not set `source_geometry_scale` to a value other than `1.0`.

- [ ] **Step 3: Confirm the option-A avalanche Jacobian is present before any -20 V claim**

Reaching -20 V depends on the already-landed option-A avalanche Jacobian behavior that includes the breakdown-critical field-dependent alpha loop-gain terms. Without it, the calibrated BV sweep historically stalls near `-13.2 V` with `line_search_non_decrease`, and the Task 7/8 "-20 V" comparison silently degrades to a ~-13 V diagnostic. Verify the enabling behavior through the focused impact-ionization tests compiled into the `build-release` binary before running the candidate.

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --parallel --target test_impact_ionization
ctest --test-dir build-release --output-on-failure -R "impact"
```

Expected:

```text
The test_impact_ionization case "Coupled DD SG edge-current avalanche Jacobian captures field-dependent alpha" passes.
If the focused impact-ionization test fails, STOP: the -20 V candidate cannot be trusted; treat any high-bias result as a <=-13 V diagnostic only.
```

- [ ] **Step 4: Run the derived candidate**

Run:

```powershell
build-release\vela_example_runner.exe --config build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_validation_candidate\simulation_bv_minus20_validation.json
```

Expected:

```text
The run exits 0 or records a controlled non-convergence breakdown row.
If it does not reach -20 V even with the option-A Jacobian present, keep the result diagnostic and do not promote a -20 V ctest gate.
```

### Task 8: Compare BV Curves And Multibias Fields

**Files:**
- Generated only: `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_validation_compare/`
- Reference: `scripts/compare_pn2d_bv_multibias_fields.py`

- [ ] **Step 1: Verify the ACTUAL generated artifact filenames before running the comparator**

The comparator invocation in Step 3 hard-codes curve/field paths (`reference_curves\pn2d_sentaurus2018_bv_reference.csv`, `reports\bv_validation_candidate\iv.csv`, `sentaurus_multibias\sentaurus_<bias>v`). These names are assumptions about `sentaurus_import.py` / candidate-run output layout and WILL break the comparator on any mismatch. Discover the real names first and substitute them into the later steps.

Run:

```powershell
Get-ChildItem -Recurse build-release\reference_tcad\pn2d_sentaurus2018\reference_curves\*bv*.csv
Get-ChildItem -Recurse build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_validation_candidate\*.csv
Get-ChildItem -Directory build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias\
python scripts\compare_pn2d_bv_multibias_fields.py --help
```

Expected:

```text
The actual BV reference-curve CSV, candidate IV CSV, and per-bias Sentaurus export directory names are listed.
The comparator --help confirms the real flag names (--curve-reference/--curve-candidate/--sentaurus-root/--vela-vtk-root or their actual equivalents).
Update the paths/flags in Step 2 and Step 3 to match before running; do not run with the placeholder names unmodified.
```

- [ ] **Step 2: Export selected Sentaurus multibias fields if needed**

If the import did not already create `sentaurus_multibias/sentaurus_<bias>v` folders, run the existing `sentaurus_import.exe --tdr ... --export-dir ...` pattern for these bias files:

```text
pn2d_bv_multibias_0000_des.tdr -> sentaurus_0v
pn2d_bv_multibias_0005_des.tdr -> sentaurus_-0.5v
pn2d_bv_multibias_0020_des.tdr -> sentaurus_-2v
pn2d_bv_multibias_0050_des.tdr -> sentaurus_-5v
pn2d_bv_multibias_0100_des.tdr -> sentaurus_-10v
pn2d_bv_multibias_0200_des.tdr -> sentaurus_-20v
```

Expected:

```text
Each selected export contains nodes.csv, elements.csv, and field CSVs for ElectricField, eDensity, hDensity, and AvalancheGeneration or ImpactIonization.
```

- [ ] **Step 3: Run the multibias comparator**

Run:

```powershell
python scripts\compare_pn2d_bv_multibias_fields.py `
  --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_validation_candidate\vtk `
  --curve-reference build-release\reference_tcad\pn2d_sentaurus2018\reference_curves\pn2d_sentaurus2018_bv_reference.csv `
  --curve-candidate build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_validation_candidate\iv.csv `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_validation_compare `
  --biases 0,-0.5,-2,-5,-10,-20
```

Expected:

```text
The comparator writes markdown/JSON/CSV reports.
Avalanche generation is interpreted as generation/source-density comparison, not as alpha comparison.
```

- [ ] **Step 4: Record the model decision**

Write a short decision block in the comparison report or `docs/validation/pn2d_bv_validation.md`:

```markdown
## Avalanche Model Decision

The refreshed Sentaurus BV source and generated artifacts use `Avalanche(VanOverstraeten)`. Vela already implements this coefficient family through `van_overstraeten`; therefore no Okuto-Crowell model is added in this validation pass.
```

If the refreshed run proves a different active model, stop here and write a separate implementation plan for the new model.

### Task 9: Add Lightweight Regression Gates

**Files:**
- Modify: `scripts/run_regression.py`
- Test: `tests/regression/test_run_regression.py`
- Optional modify: `tests/test_dc_sweep.cpp`

- [ ] **Step 1: Write tests for trend-only BV summary semantics**

In `tests/regression/test_run_regression.py`, add a fixture CSV with BV rows containing:

```text
bias_V,current_total_A_per_um,converged,max_electric_field_V_per_m,breakdown_detected,breakdown_voltage,criterion,last_stable_bias,failed_bias,breakdown_failure_reason
0,-1e-18,1,1.0e5,0,,,,,
-1,-2e-18,1,1.5e5,0,,,,,
-2,-4e-18,1,2.0e5,0,,,,,
```

Expected regression result:

```json
{
  "max_field_trend_checked": true,
  "breakdown_detected": false
}
```

Add a negative fixture where max electric field decreases with increasing reverse magnitude and assert that the check fails with text containing:

```text
max electric field is not monotonic
```

- [ ] **Step 2: Run the tests and confirm failure before implementation**

Run:

```powershell
python -m unittest tests.regression.test_run_regression -v
```

Expected before implementation:

```text
The new monotonic max-field test fails if the existing regression checker does not enforce the trend.
```

- [ ] **Step 3: Implement only the lightweight trend gate**

Update `scripts/run_regression.py` so BV examples check:

```text
finite outputs
accepted-row convergence semantics
breakdown diagnostic columns
max electric field non-decreasing with increasing reverse-bias magnitude, allowing a small numerical tolerance
```

Do not add a Sentaurus-dependent ctest and do not require -20 V absolute-current parity in `run_regression.py`.

- [ ] **Step 4: Run regression unit tests**

Run:

```powershell
python -m unittest tests.regression.test_run_regression -v
```

Expected:

```text
All run_regression unit tests pass.
```

### Task 10: Write The PN2D BV Methodology Document

**Files:**
- Create or modify: `docs/validation/pn2d_bv_validation.md`
- Reference: `docs/validation/bv_impact_ionization_theory.md`
- Reference: `docs/validation/pn2d_sentaurus_comparison.md`

- [ ] **Step 1: Write methodology sections**

Create or update `docs/validation/pn2d_bv_validation.md` with these sections:

````markdown
# PN2D BV Validation Methodology

## Source Of Truth

The source deck is `reference_tcad/pn2d_sentaurus2018/source/pn2d_bv_sdevice.cmd`. It uses `Avalanche(VanOverstraeten)` and sweeps Anode to `-20.0 V`.

## Artifact Refresh

Sentaurus artifacts are refreshed through `scripts/run_sentaurus_vm_reference.py pn2d --stages 0v,iv,bv`. Dry-run planning uses `--dry-run`; live upload/run/fetch omits `--dry-run`.

## Required Artifact Checks

The BV refresh must contain `pn2d_bv.plt`, a clean BV log artifact such as `pn2d_bv.log` or `pn2d_bv.log_des.log`, and 201 endpoint-inclusive `pn2d_bv_multibias_*_des.tdr` files.

## Model Decision

When the refreshed source remains `Avalanche(VanOverstraeten)`, Vela's existing `van_overstraeten` implementation is the target and Okuto-Crowell remains contrast-only.

## Comparison Layers

1. Curve comparison checks current trend and documented windows.
2. Field comparison checks potential, electric field, carrier density, mobility, and avalanche generation/source density at selected biases.
3. Coefficient checks compare alpha to alpha; generation checks compare generation/source integral to generation/source integral.

## Accepted Gates

The promoted automated gates are VM-free and lightweight: parser/import provenance, artifact validation with synthetic fixtures, BV max-field trend, and documented comparison summaries.

## Non-Goals

This pass does not claim full `0..-20 V` absolute-current parity, does not promote hidden scalar source calibration, does not rewrite SG flux divergence, and does not add LDMOS/IGBT/MOS BV validation.

## High-Bias Interpretation

Windowed current agreement and high-bias knee shape are reported separately. If knee-shape evidence remains divergent, the methodology records it as an open physics/parity limit rather than hiding it behind a broad current band.
````

- [ ] **Step 2: Verify docs do not contradict current CLI**

Run:

```powershell
rg -n -- "--plan|--upload|--launch|--fetch|200 `pn2d_bv_multibias|Okuto-Crowell.*target|absolute-current parity" docs\validation docs\superpowers\plans\2026-06-22-bv-impact-ionization-validation-revised.md
```

Expected:

```text
No new document says the VM runner supports --plan/--upload/--launch/--fetch unless Task 2 explicitly implemented those flags.
No new document says there are only 200 endpoint-inclusive BV multibias TDR files.
No new document treats Okuto-Crowell as the active PN2D BV target.
No new document claims full -20 V absolute-current parity.
```

### Task 11: Verification Before Completion

**Files:**
- No new files unless a previous task changed implementation.

- [ ] **Step 1: Build required binaries**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --preset windows-ucrt64-debug
cmake --build build --parallel --target vela_example_runner sentaurus_import test_impact_ionization test_dc_sweep
cmake --build build-release --parallel --target vela_example_runner sentaurus_import test_impact_ionization
```

Expected:

```text
All requested debug and release targets build successfully. The release targets are included because Tasks 6-8 use `build-release`.
```

- [ ] **Step 2: Run focused C++ tests**

Run:

```powershell
ctest --test-dir build --output-on-failure -R "impact|dc_sweep"
ctest --test-dir build-release --output-on-failure -R "impact"
```

Expected:

```text
All selected debug tests pass, and the release impact-ionization test used by the -20 V candidate path passes.
```

- [ ] **Step 3: Run focused Python regression tests**

Run:

```powershell
python -m unittest tests.regression.test_sentaurus_vm_reference_runner tests.regression.test_sentaurus_import_tools tests.regression.test_reference_tcad_tools tests.regression.test_run_regression -v
```

Expected:

```text
All selected regression tests pass.
```

- [ ] **Step 4: Run full regression only if the workspace has enough time**

Run:

```powershell
ctest --test-dir build --output-on-failure
```

Expected:

```text
The full suite passes, or any failures are documented with exact failing tests and why they are unrelated or blocked.
```

---

## Final Acceptance Criteria

- The revised docs name VanOverstraeten as the active PN2D BV target and keep Okuto-Crowell conditional.
- The VM workflow uses the actual current runner CLI or includes tests for any newly added CLI flags.
- The artifact validator expects 201 endpoint-inclusive BV multibias TDR files.
- Any Vela -20 V run is based on an explicit derived candidate deck or a separately reviewed reference promotion.
- Avalanche coefficient comparisons are not mixed with generation/source-density comparisons.
- Automated gates remain lightweight and reproducible without VM/license access.
- High-bias knee-shape mismatch, if still present, is reported as diagnostic evidence rather than hidden behind a broad current tolerance.
