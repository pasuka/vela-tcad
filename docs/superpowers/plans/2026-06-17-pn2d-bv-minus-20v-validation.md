# PN2D BV -20 V Validation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reproduce the pn2d Sentaurus2018 default BV avalanche semantics in Vela, extend the validation from the current low-reverse-bias gate to a 0 V to -20 V run only after the Sentaurus-default source path passes high-bias parity checks, then compare BV current and same-bias spatial fields against Sentaurus.

**Architecture:** Use the existing `reference_tcad/pn2d_sentaurus2018` fixture as the source of truth, with Sentaurus's default Scharfetter-Gummel edge-current avalanche discretization as the required Vela parity target. Keep Sentaurus field export and Vela execution separate, then join them through `scripts/compare_pn2d_bv_multibias_fields.py` so curve error, field error, and avalanche-source error are ranked at the same bias points. Treat `current_approximation = "density_gradient"` as the Sentaurus-default path; keep `mobility_density_gradient` only as an `AvalDensGradQF`-like control path. When the Sentaurus 2018 VM is reachable over SSH, use it as an opt-in oracle to regenerate BV reference artifacts and one-factor Sentaurus variants instead of relying only on previously captured `.tdr` files.

**Tech Stack:** C++20, CMake/Ninja, MSYS2 UCRT64 on Windows, Python standard library plus NumPy for analysis scripts, existing Vela runner and Sentaurus import tooling, Windows OpenSSH/`scp` to the `sentaurus` SSH host for opt-in Sentaurus 2018 VM runs.

---

## Current Branch Baseline

- Latest commit: `768ba5ee092dad231c145758d2709d8cbfe1b1d3` (`Align pn2d BV Sentaurus physics diagnostics`).
- Branch: `codex-pn2d-sentaurus2018-calibration`.
- Worktree status observed during planning: clean.
- The commit adds `scripts/compare_pn2d_bv_multibias_fields.py` and `scripts/diagnose_pn2d_bv_mobility.py`, expands impact-ionization support, and updates the pn2d Sentaurus2018 reference JSON.
- Sentaurus BV deck `reference_tcad/pn2d_sentaurus2018/source/pn2d_bv_sdevice.cmd` already targets `Goal { Name="Anode" Voltage=-20.0 }` and writes `pn2d_bv_multibias` snapshots over 200 intervals.
- New capability as of 2026-06-19: the local Windows machine can access a Sentaurus 2018 VM through SSH. Use the existing `sentaurus` SSH alias and `docs/sentaurus_vm_ssh_workflow.md` as the operational guide. Treat live Sentaurus execution as opt-in and keep fetched VM results under `build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_vm_runs/` until they are reviewed; do not overwrite `reference_tcad/pn2d_sentaurus2018/source/` during exploratory runs.
- Vela pn2d BV config in `reference_tcad/pn2d_sentaurus2018/pn2d_sentaurus2018_reference.json` still gates only `vela_stop: -0.05` with `vela_step: -0.05`.

## Impact-Ionization Formula Check

Vela now has two executable impact-ionization models:

```text
selberherr:
  alpha_n(E) = A_n * exp(-B_n / |E|)
  alpha_p(E) = A_p * exp(-B_p / |E|)
  G_density = v_sat * (alpha_n * max(n, 0) + alpha_p * max(p, 0))

van_overstraeten:
  gamma(T) = tanh(hbar_omega / (2 k T_ref)) / tanh(hbar_omega / (2 k T))
  alpha(E) = gamma * A_region * exp(-(gamma * B_region) / |E|)
  low/high region selected by |E| < switch_field
  G_current = alpha_n * mu_n * max(n, 0) * |F_n| +
              alpha_p * mu_p * max(p, 0) * |F_p|
```

Default coefficients in `include/vela/physics/ImpactIonizationModel.h` are SI:

```text
Selberherr defaults:
  electron_A = 7.03e7 1/m
  electron_B = 1.231e8 V/m
  hole_A     = 1.582e8 1/m
  hole_B     = 2.036e8 V/m
  carrier_velocity = 1.0e5 m/s

Van Overstraeten defaults:
  electron_a_low  = 7.03e7 1/m
  electron_a_high = 7.03e7 1/m
  electron_b_low  = 1.231e8 V/m
  electron_b_high = 1.231e8 V/m
  hole_a_low      = 1.582e8 1/m
  hole_a_high     = 6.71e7 1/m
  hole_b_low      = 2.036e8 V/m
  hole_b_high     = 1.693e8 V/m
  switch_field    = 4.0e7 V/m
  phonon_energy   = 0.063 eV
  T_ref = T = 300 K by default, so gamma = 1 at 300 K
```

The local Sentaurus2018 BV deck uses `Avalanche(VanOverstraeten)`, not Okuto-Crowell. If another pn2d source deck still uses Okuto-Crowell, handle it as a separate model-port task rather than mixing coefficients into this Van Overstraeten validation.

## Sentaurus SDevice Manual Cross-Check 2026-06-18

The Sentaurus Device User Guide P-2019.03 was re-read with PyMuPDF from `D:\工作\学习资料\TCAD软件手册\Sentaurus PDFManual 2019\data\sdevice_ug.pdf`. Relevant printed pages are `425` through `429` and `439` through `440`; the corresponding PDF one-based pages are `469` through `473` and `483` through `484`.

Manual conclusions to carry into the remaining work:

- Printed page `425`: `Avalanche` defaults to `vanOverstraeten`, and the isothermal default driving force is `GradQuasiFermi`; only hydrodynamic simulations default to `CarrierTempDrive`.
- Printed page `426`: Sentaurus default Eq. 431 computes `Jn` and `Jp` in each mesh element using the Scharfetter-Gummel approximation applied to each element edge. The alternative `Math { AvalDensGradQF }` approximation uses `Jn = -q mu_n n grad(Phi_n)` and `Jp = -q mu_p p grad(Phi_p)`, can slightly shift breakdown voltage, and is documented as having better convergence and stability properties.
- Printed page `427`: `Math { ElementVolumeAvalanche }` uses truncated element-vertex volumes for avalanche so obtuse elements do not exaggerate generated carriers. `AvalFlatElementExclusion` can exclude nearly flat elements from avalanche generation.
- Printed pages `428` and `429`: Vela's Van Overstraeten coefficient formula and Table 81 defaults match the manual after converting Sentaurus `cm` units to SI. Exact `E == E0` low/high selection is a negligible boundary convention unless a test lands exactly on the switch field.
- Printed page `439`: `GradQuasiFermi` and `Eparallel` driving forces are affected by `ParallelToInterfaceInBoundaryLayer`. Sentaurus also supports interpolation of avalanche driving force to `ElectricField` at low carrier concentrations through `RefDens_*GradQuasiFermi_ElectricField` and `RefDens_*Eparallel_ElectricField_Aval`.
- Printed page `440`: hydrodynamic `CarrierTempDrive` formulas do not apply to the current pn2d drift-diffusion deck.

Revised interpretation:

```text
Current pn2d source deck:
  Recombination(Avalanche(VanOverstraeten))
  no BandgapDependence
  no explicit AvalDensGradQF

Current Vela imported deck:
  model = van_overstraeten
  driving_force = quasi_fermi_gradient
  generation = current_density
  current_approximation must be density_gradient for Sentaurus-default BV parity

Parity target:
  Sentaurus default BV uses SG edge-current avalanche generation.
  Vela current_approximation = density_gradient is the required parity path.
  Vela current_approximation = mobility_density_gradient is closest to
  Sentaurus Math { AvalDensGradQF } and must remain a control path, not the
  acceptance path for Sentaurus-default BV.

Current -13.208 V failure priority:
  1. No-impact high-bias branch parity, because the carrier-density mismatch is
     already visible before avalanche feedback is enabled.
  2. Sentaurus-default SG edge-current avalanche source parity.
  3. Avalanche hotspot geometry/control-volume amplification.
  4. Analytic Jacobian completeness for avalanche driving-field and mobility terms.
  5. Low-density driving-force interpolation to ElectricField only if the
     Sentaurus deck enables RefDens_*Aval or a controlled probe requires it.
  6. BandgapDependence only if a future Sentaurus deck explicitly enables it.
```

## Files To Modify Or Use

- Modify: `reference_tcad/pn2d_sentaurus2018/pn2d_sentaurus2018_reference.json`
  - Keep the committed BV gate at low reverse bias until Sentaurus-default SG edge-current parity passes.
  - When the gate is promoted, extend BV run target to `vela_stop: -20.0` and require `solver.impact_ionization.current_approximation: "density_gradient"`.
  - Use a two-stage step policy through generated probe configs instead of immediately promoting a coarse `-0.1` or `-0.5` step.
- Use: `scripts/sentaurus_import.py`
  - Regenerate Vela decks and Sentaurus curve references.
- Use: `scripts/compare_pn2d_bv_multibias_fields.py`
  - Compare potential, electric field, electron density, hole density, SRH recombination, avalanche generation, and mobility fields.
- Use: `scripts/diagnose_pn2d_bv_mobility.py`
  - Localize mobility and high-field drive mismatches when the curve or field comparison fails.
- Create: `scripts/diagnose_pn2d_bv_avalanche_hotspots.py`
  - Read a Vela VTK at the last stable avalanche bias, rank avalanche generation by node, and join it to mesh geometry/control-volume proxies.
  - Report whether the avalanche source is concentrated in a few nodes/elements before changing physics.
- Modify: `include/vela/physics/ImpactIonizationModel.h`, `include/vela/equation/AssemblerUtils.h`, `src/equation/CoupledDDAssembler.cpp`, `src/equation/DDAssembler.cpp`, `src/solver/NewtonSolver.cpp`, and `src/solver/GummelSolver.cpp`
  - Only if the diagnostic confirms the need for new solver knobs: add explicit avalanche discretization and driving-force-interpolation configuration.
- Modify: `include/vela/simulation/DCSweep.h` and `src/simulation/DCSweep.cpp`
  - Add opt-in continuation predictor and branch-acceptance diagnostics only after the external predictor proxy is converted into TDD-covered behavior.
- Create: `scripts/run_sentaurus_vm_reference.py`
  - Add a deterministic, opt-in SSH/SCP runner for the Sentaurus 2018 VM. It must support `--dry-run` for unit tests, stage every live run under a timestamped or explicit run id, and copy artifacts into run-specific directories such as `build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_vm_runs/pn2d_bv_vm_smoke/source`.
- Create: `tests/regression/test_sentaurus_vm_reference_runner.py`
  - Cover command planning, dry-run manifests, source-file validation, and live-test opt-in gates without requiring the VM in ordinary test runs.
- Modify: `docs/sentaurus_vm_ssh_workflow.md`
  - Add a short machine-readable contract section for the runner: SSH alias `sentaurus`, remote root `~/sentaurus_runs/vela_oracle`, source deck directory, expected output suffixes, and cleanup policy.
- Use: `D:\code-repo\tcad-charon\src\evaluators\Charon_Avalanche_vanOverstraeten_impl.hpp`
  - Reference Charon's avalanche driving-force variants, low-density force damping, and `alpha * |J|` generation form when checking Vela source semantics.
- Use: `D:\code-repo\tcad-charon\src\Charon_CurrentConstraintModelEvaluator.hpp` and `D:\code-repo\tcad-charon\src\solver\Charon_Solver_SteadyStateConstraint.cpp`
  - Reference Charon's current-constraint continuation design before adding any Vela current-controlled BV path.
- Use: `D:\code-repo\devsim\python_packages\ramp.py`, `D:\code-repo\devsim\python_packages\simple_dd.py`, and `D:\code-repo\devsim\src\AutoEquation\ExprEquation.cc`
  - Reference DEVSIM's bias-ramp fallback, Scharfetter-Gummel expression models, and edge/node-volume assembly when designing lightweight Vela probes.
- Test: `tests/test_impact_ionization.cpp`
  - Add coefficient and generation checks only if the formula review finds a mismatch.
- Test: `tests/test_newton_solver.cpp` and `tests/test_dc_sweep.cpp`
  - Add parser/restart coverage for any new avalanche discretization or driving-force-interpolation knobs.
- Test: `tests/regression/test_reference_tcad_tools.py`
  - Add regression coverage if command-line behavior or generated report schema changes.

## Task 1: Reproduce The Imported Reference And Existing Gate

- [x] **Step 1: Configure the MSYS2 UCRT64 build**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --preset windows-ucrt64-release
```

Expected: CMake configures `build-release/` successfully and finds Ninja plus the UCRT64 compiler.

- [x] **Step 2: Build the required tools**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --parallel --target vela_example_runner sentaurus_import test_impact_ionization
```

Expected: `build-release/vela_example_runner.exe`, `build-release/sentaurus_import.exe`, and `build-release/test_impact_ionization.exe` exist.

- [x] **Step 3: Run the impact-ionization unit tests before changing BV reach**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
   build-release\test_impact_ionization.exe
```

Expected: `All tests passed`, including Van Overstraeten coefficient checks.

- [x] **Step 4: Regenerate the pn2d Sentaurus2018 imported reference tree**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts\sentaurus_import.py reference --config reference_tcad\pn2d_sentaurus2018\pn2d_sentaurus2018_reference.json --source-dir reference_tcad\pn2d_sentaurus2018\source --output-dir build-release\reference_tcad\pn2d_sentaurus2018 --tdr-importer build-release\sentaurus_import.exe --runner build-release\vela_example_runner.exe
```

Expected: current low-bias BV gate still runs, and `build-release/reference_tcad/pn2d_sentaurus2018/reference_curves/pn2d_sentaurus2018_bv_reference.csv` exists.

## Task 2: Import Sentaurus BV Multi-Bias Field Snapshots

- [x] **Step 1: Confirm the source snapshots exist**

Run:

```powershell
Get-ChildItem reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_multibias_*_des.tdr | Measure-Object
```

Expected: count is at least `201`, covering normalized times from 0 to 1.

- [x] **Step 2: Export the comparison bias points**

Run one command per bias:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\sentaurus_import.exe --tdr reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_multibias_0000_des.tdr --export-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias\sentaurus_0v
```

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\sentaurus_import.exe --tdr reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_multibias_0005_des.tdr --export-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias\sentaurus_-0.5v
```

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\sentaurus_import.exe --tdr reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_multibias_0020_des.tdr --export-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias\sentaurus_-2v
```

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\sentaurus_import.exe --tdr reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_multibias_0050_des.tdr --export-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias\sentaurus_-5v
```

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\sentaurus_import.exe --tdr reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_multibias_0100_des.tdr --export-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias\sentaurus_-10v
```

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\sentaurus_import.exe --tdr reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_multibias_0200_des.tdr --export-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias\sentaurus_-20v
```

Expected: each output directory contains `nodes.csv`, `elements.csv`, and `fields/*.csv`.

- [x] **Step 3: Verify the imported Sentaurus fields needed by the comparison**

Run:

```powershell
Get-ChildItem build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias\sentaurus_-20v\fields
```

Expected: files exist for `ElectrostaticPotential`, `ElectricField`, `eDensity`, `hDensity`, and `AvalancheGeneration` or `ImpactIonization`.

## Task 3: Create A Controlled Vela -20 V Probe Deck

- [x] **Step 1: Copy the generated BV deck to a probe file**

Run:

```powershell
Copy-Item build-release\reference_tcad\pn2d_sentaurus2018\vela\simulation_bv.json build-release\reference_tcad\pn2d_sentaurus2018\vela\simulation_bv_minus20_probe.json
```

Expected: the probe deck exists and the committed reference JSON remains unchanged during the first exploration.

- [x] **Step 2: Edit the probe deck for a conservative reach test**

Change only the BV sweep controls:

```json
"sweep": {
  "mode": "bv_reverse",
  "stop": -20.0,
  "step": -0.1,
  "write_vtk": true
}
```

Keep solver physics unchanged from the imported BV deck:

```json
"impact_ionization": {
  "model": "van_overstraeten",
  "driving_force": "quasi_fermi_gradient",
  "generation": "current_density"
}
```

Expected: the probe deck remains valid JSON and does not change IV settings.

- [x] **Step 3: Run the Vela -20 V probe**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\vela_example_runner.exe --config build-release\reference_tcad\pn2d_sentaurus2018\vela\simulation_bv_minus20_probe.json
```

Expected: the run either reaches -20 V with VTK files at comparison bias points, or stops with a clear last-stable bias and failure reason in the BV CSV.

- [x] **Step 4: If the probe fails before -20 V, rerun with staged continuation**

Edit only the probe deck step size and rerun:

```json
"step": -0.05
```

Expected: the last-stable bias improves or the failure mode stays identical, which separates continuation stiffness from step-size artifacts.

### Execution Note 2026-06-18

- Release configure/build succeeded with `windows-ucrt64-release`; `build-release/test_impact_ionization.exe` passed all 71 assertions in 8 test cases.
- Regenerating the imported reference tree into `build-release/reference_tcad/pn2d_sentaurus2018` succeeded.
- The `simulation_bv_minus20_probe.json` run with Van Overstraeten impact ionization did not reach -20 V. It stopped after 901 points with last stable bias `-13.208218617327727 V`, failed bias `-13.208218617474687 V`, attempted step `-1.4695977768042212e-10 V`, `retry_count=26`, and `failure_reason=line_search_non_decrease`.
- The Newton failure diagnostics show the residual is dominated by the electron quasi-Fermi block (`phin=7.18897253390126e-09`) while Poisson residuals remain small (`psi=9.925545369294046e-12`) and carrier densities remain positive and finite.
- A no-impact-ionization control deck with the same -20 V sweep converged to `-20.000000000000014 V`, so the observed failure requires avalanche coupling rather than the reverse-bias continuation, contacts, mesh, or mobility model alone.
- A restart from the last stable saved state with analytic Jacobian immediately reproduced `line_search_non_decrease`. A finite-difference Jacobian restart was too slow for this 1943-node mesh in the interactive run and timed out before producing a useful comparison.

### Execution Note 2026-06-19

- `cmake --preset windows-ucrt64-release` completed and `cmake --build build-release --parallel` completed after a resumed Ninja build.
- `build-release\test_impact_ionization.exe` passed `123 assertions in 13 test cases`.
- Sentaurus multibias exports were regenerated for `0`, `-0.5`, `-2`, `-5`, `-10`, and `-20 V`; each export contains `ElectrostaticPotential`, `ElectricField`, `eDensity`, `hDensity`, and `ImpactIonization` field CSVs.
- The direct `simulation_bv_minus20_probe.json` run reached the high-bias branch but hit the interactive timeout after the last complete CSV row at `-17.3 V`; the next row was partial because the process was terminated.
- A restart state was reconstructed from `dc_sweep_1116_-17.3V.vtk` into `restart_from_minus17p3_state.csv`, then `simulation_bv_minus20_resume_from_minus17p3.json` resumed the same BV path from `-17.3 V` to `-20 V`.
- The resume segment completed with `converged=true`, `points=55`, and last bias `-20.000000000000039 V`.
- A clean combined curve was written to `pn2d_sentaurus2018_bv_minus20_combined.csv`, omitting the timeout-truncated row and appending the complete resume segment.

## Task 4: Compare BV Curve And Spatial Fields

- [x] **Step 1: Run the multibias comparison script**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts\compare_pn2d_bv_multibias_fields.py --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\vela --curve-reference build-release\reference_tcad\pn2d_sentaurus2018\reference_curves\pn2d_sentaurus2018_bv_reference.csv --curve-candidate build-release\reference_tcad\pn2d_sentaurus2018\vela\pn2d_sentaurus2018_bv.csv --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_multibias --biases 0,-0.5,-2,-5,-10,-20
```

Expected: `curve_compare.csv`, `field_compare.csv`, `debug_ranking.json`, and `README.md` are written.

- [x] **Step 2: Read the ranked failure order**

Run:

```powershell
Get-Content build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_multibias\debug_ranking.json
```

Expected: the first ranked items identify whether the dominant mismatch is curve current, potential, electric field, carrier density, mobility, SRH recombination, or avalanche generation.

- [x] **Step 3: Promote only fields that are same-bias and present on both sides**

Acceptance:

```text
potential: compare RMS error after same-bias pairing
electric_field: compare relative p95 and max-field location
electron_density: compare log10 p95
hole_density: compare log10 p95
avalanche_generation: compare log10 p95 only when both p99 values are above the script avalanche floor
bv_curve: compare abs log10 current ratio at non-floor current points
```

Expected: no field comparison is called passing when its Vela VTK is missing, its Sentaurus export is missing, or the bias is interpolated from the wrong Sentaurus snapshot.

### Execution Note 2026-06-18

- Sentaurus snapshots for `0`, `-0.5`, `-2`, `-5`, `-10`, and `-20 V` were exported under `build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias`.
- The multibias comparison report was generated under `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_multibias`.
- `curve_compare.csv` reports candidate data at `0`, `-0.5`, `-2`, `-5`, and `-10 V`; the `-20 V` row is `missing_candidate` because the Vela avalanche run stopped at `-13.208218617474687 V`.
- Pre-failure curve current errors are about `0.18` to `0.22` decades at the non-floor comparison points.
- `debug_ranking.json` ranks `electron_density`, `electron_mobility`, `electric_field`, and `avalanche_generation_thresholded` as the leading follow-up fields, but the first numerical blocker is the avalanche-enabled electron-continuity Newton line search near `-13.208 V`.

### Execution Note 2026-06-19

- Re-ran `scripts\compare_pn2d_bv_multibias_fields.py` using the regenerated Sentaurus multibias exports, Vela VTKs including the resumed `-20 V` state, and `pn2d_sentaurus2018_bv_minus20_combined.csv`.
- Report outputs were regenerated under `build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_multibias`: `curve_compare.csv`, `field_compare.csv`, `debug_ranking.json`, `pn2d_bv_multibias_field_compare.csv`, `pn2d_bv_multibias_field_compare.json`, and `README.md`.
- `curve_compare.csv` now contains an `ok` row at `-20 V`; `debug_ranking.json` reports the largest curve error at `-20 V` with absolute log10 current error `1.3880802290078997`.
- The leading ranked field mismatch remains electron density at `-20 V`, with `log10_p95` error `4.009914625448848`. This preserves the high-bias branch-mismatch diagnosis while proving that the current code path can be continued to `-20 V` via a restart segment.

## Task 5: Localize Manual-Guided Avalanche Parity And Stability Gaps

- [x] **Step 1: Confirm the Sentaurus BV deck uses default avalanche options**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
rg -n "Avalanche|eAvalanche|hAvalanche|BandgapDependence|AvalDensGradQF|ElementVolumeAvalanche|RefDens_.*Aval|ParallelToInterfaceInBoundaryLayer" reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_sdevice.cmd reference_tcad\pn2d_sentaurus2018\source\models.par
```

Expected:

```text
reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_sdevice.cmd:21:    Avalanche(VanOverstraeten)
reference_tcad\pn2d_sentaurus2018\source\models.par:...: vanOverstraetendeMan * Impact Ionization:
```

There must be no active `BandgapDependence`, `AvalDensGradQF`, `ElementVolumeAvalanche`, or `RefDens_*Aval` match in `pn2d_bv_sdevice.cmd`. If any of those are present, update the Vela parity target before continuing.

- [x] **Step 2: Keep Van Overstraeten coefficient parity as a guardrail, not the first fix**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\test_impact_ionization.exe
```

Expected: `All tests passed`. If this fails, fix the coefficient/unit test first. If it passes, do not tune `a`, `b`, `E0`, or `hbarOmega` while investigating the -13.208 V Newton failure.

- [x] **Step 3: Add an avalanche hotspot and geometry diagnostic**

Create `scripts/diagnose_pn2d_bv_avalanche_hotspots.py` with this interface:

```text
required args:
  --vtk build-release\reference_tcad\pn2d_sentaurus2018\vela\dc_sweep_0899_-13.2082V.vtk
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\vela\mesh.json
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_avalanche_hotspots

required outputs:
  avalanche_hotspots.csv
  avalanche_hotspots_summary.json
```

The report must include, for each top avalanche node, at least:

```text
node_id
x_um
y_um
avalanche_generation_m3_s
electron_density_m3
hole_density_m3
electric_field_V_m
electron_high_field_drive_V_m
hole_high_field_drive_V_m
node_control_volume_m2
adjacent_element_count
min_adjacent_angle_degrees
max_adjacent_angle_degrees
obtuse_adjacent_element_count
```

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts\diagnose_pn2d_bv_avalanche_hotspots.py --vtk build-release\reference_tcad\pn2d_sentaurus2018\vela\dc_sweep_0899_-13.2082V.vtk --mesh build-release\reference_tcad\pn2d_sentaurus2018\vela\mesh.json --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_avalanche_hotspots
```

Expected: `avalanche_hotspots_summary.json` identifies whether the top 1, 5, and 20 nodes dominate avalanche generation and whether they sit on obtuse or high-control-volume mesh geometry.

- [x] **Step 4: If hotspots concentrate on obtuse/high-volume nodes, implement an ElementVolumeAvalanche-style probe**

Add a Vela-only probe option under `solver.impact_ionization`:

```json
"element_volume_policy": "geometric"
```

Keep the default as:

```json
"element_volume_policy": "node_control_volume"
```

Expected implementation behavior:

```text
node_control_volume:
  current Vela behavior; multiply nodal avalanche G by the existing node control volume.

geometric:
  distribute each adjacent triangle's geometric area equally to its three nodes for avalanche only,
  so an obtuse element cannot contribute more avalanche volume than its geometric area.
```

Run the -20 V probe with only this option changed:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\vela_example_runner.exe --config build-release\reference_tcad\pn2d_sentaurus2018\vela\simulation_bv_minus20_element_volume_probe.json
```

Expected: compare last-stable bias against `-13.208218617327727 V`. If the run moves substantially farther or current drops sharply at the hotspot, prioritize geometric avalanche volume parity before Jacobian tuning.

- [x] **Step 5: Add a Sentaurus-default SG edge-current avalanche probe before any -20 V promotion**

Add a Vela-only probe option under `solver.impact_ionization`:

```json
"current_approximation": "density_gradient"
```

Keep the current implementation available as:

```json
"current_approximation": "mobility_density_gradient"
```

Expected mapping:

```text
mobility_density_gradient:
  current Vela behavior; G = alpha_n mu_n n |grad(Phi_n)| + alpha_p mu_p p |grad(Phi_p)|.
  This is closest to Sentaurus Math { AvalDensGradQF }.

density_gradient:
  Required Sentaurus-default parity path; derive avalanche current contribution from the same SG edge-current approximation used by the drift-diffusion fluxes, then accumulate Eq. 431 by element/node.
```

Acceptance: do not promote `vela_stop: -20.0` until the `density_gradient` SG edge-current path is explicitly documented, the Sentaurus deck is confirmed to remain on its implicit default avalanche discretization, and the high-bias branch/source gates in Tasks 7-9 pass.

- [x] **Step 6: Add a GradQuasiFermi-to-ElectricField interpolation probe**

Add optional damping densities under `solver.impact_ionization`:

```json
"driving_force_interpolation": {
  "mode": "quasi_fermi_to_electric_field",
  "electron_ref_density_m3": 1.0e16,
  "hole_ref_density_m3": 1.0e16
}
```

The interpolation must follow the manual intent on printed page 439: at low carrier density, blend the avalanche driving force toward the plain electric field. Probe at least these densities:

```text
1.0e14 m^-3
1.0e16 m^-3
1.0e18 m^-3
```

Expected: if one interpolation setting crosses the -13.208 V failure without materially changing pre-failure BV current at `-10 V`, keep it as a solver-stability candidate and compare fields again.

- [x] **Step 7: Only after Steps 3-6, revisit analytic Jacobian sensitivity**

Create restart probes from `pn2d_bv_minus20_last_stable_state.csv` differing only in:

```json
"jacobian": "analytic"
```

and:

```json
"jacobian": "finite_difference"
```

Run both from the same saved state and target a small voltage window:

```text
start = -13.208218617327727
stop  = -13.25
step  = -0.005
```

Expected: if finite difference reaches farther after volume/current/interpolation probes are controlled, the analytic avalanche Jacobian must be extended to include driving-field and mobility derivatives before any Task 9 promotion gate can pass.

- [x] **Step 8: Run mobility decomposition as a secondary field-mismatch check**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts\diagnose_pn2d_bv_mobility.py --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\vela --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_mobility --biases 0,-0.5,-2,-5,-10
```

Expected: use this only after the avalanche numerical blocker is understood. Mobility remains relevant to field/current parity, but it is no longer the first branch for the -13.208 V Newton failure.

### Task 5 Execution Notes 2026-06-18

Release execution context:

- CMake preset/build directory: `windows-ucrt64-release` / `build-release`.
- Built targets: `vela_example_runner`, `sentaurus_import`, `test_impact_ionization`, and `test_newton_solver`.
- Guardrail test: `build-release\test_impact_ionization.exe` passed after adding the interpolation coverage.

Manual/parity checks:

- `pn2d_bv_sdevice.cmd` contains `Avalanche(VanOverstraeten)` and no active `BandgapDependence`, `AvalDensGradQF`, `ElementVolumeAvalanche`, or `RefDens_*Aval`.
- The original Vela `generation = "current_density"` node-local implementation is closest to the SDevice manual's `Math { AvalDensGradQF }` approximation, while the Sentaurus deck is using the default SG edge-current avalanche path. The newly added `current_approximation = "density_gradient"` implementation is therefore the required Sentaurus-default parity path.

Hotspot geometry evidence:

- Script added: `scripts/diagnose_pn2d_bv_avalanche_hotspots.py`.
- Test added: `tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_avalanche_hotspot_diagnostic_reports_geometry`.
- Report: `build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_avalanche_hotspots\avalanche_hotspots_summary.json`.
- Top hotspot: node `26`, `(x, y) = (0.03125 um, 0.1875 um)`, avalanche generation `6.52555e26 m^-3 s^-1`.
- Concentration: top 1 node `4.23%`, top 5 `21.13%`, top 20 `76.77%` of reported top-50 hotspot generation.
- Geometry: top node has `6` adjacent elements, max adjacent angle `90 deg`, `0` obtuse adjacent elements, and node control volume `9.765625e-16 m^2`, about the 80th percentile of node control volumes.
- Decision: Step 4 is closed as "not the first supported branch"; no `ElementVolumeAvalanche`-style implementation was added because the hotspot is not an obtuse/high-volume artifact.

Driving-force interpolation probes:

- Implemented optional JSON:

```json
"driving_force_interpolation": {
  "mode": "quasi_fermi_to_electric_field",
  "electron_ref_density_m3": 1.0e16,
  "hole_ref_density_m3": 1.0e16
}
```

- Implementation files: `include/vela/physics/ImpactIonizationModel.h`, `include/vela/equation/AssemblerUtils.h`, `src/solver/NewtonSolver.cpp`, `src/solver/GummelSolver.cpp`, `src/equation/CoupledDDAssembler.cpp`.
- Behavior test: `Quasi-Fermi avalanche interpolation falls back to electric field at low density`.
- JSON test: `JSON solver config selects impact ionization model`.
- Strong interpolation probe `simulation_bv_minus20_qf_interp_1e16.json` failed at `0 V` with line-search non-decrease, so it is too aggressive for this deck.
- Weak electron-only interpolation probe `simulation_bv_minus20_qf_interp_e1e2_h0.json` reached only last-stable `-13.208776212893026 V`; this is not a material improvement over the baseline `-13.208218617327727 V`.
- Decision: driving-force interpolation is not the root fix for this failure.

Continuation/Jacobian sensitivity probes:

- Local no-line-search restart from `pn2d_bv_minus20_last_stable_state.csv` still failed at the same state with `max_iterations`; residual stayed dominated by `phin ~= 7.19e-9`.
- `abstol = 1e-8` probe reached last-stable `-13.212059900665954 V`; failure residual was `phin ~= 1.0e-8`.
- `abstol = 1e-7` probe reached last-stable `-13.237598791741954 V`; failure residual was `phin ~= 1.0e-7`.
- These runs show that the strict absolute residual threshold is one visible stop trigger, but loosening it only advances slightly while Vela current grows rapidly.

Reference-current cross-check:

- Sentaurus reference current at `-13.2 V`: `-8.38472088807e-17 A`.
- Sentaurus reference current at `-20 V`: `-9.10455666344e-16 A`.
- Vela baseline near `-13.208 V`: about `-1.7e-12 A/um`.
- Vela `abstol = 1e-7` near `-13.2376 V`: about `-2.27e-11 A/um`.
- Conclusion: the Newton failure is coupled to an over-strong/too-early Vela avalanche source, not merely to an overly strict convergence gate.

SG edge-current execution notes:

- Implemented `solver.impact_ionization.current_approximation = "density_gradient"` as the Sentaurus-default SG edge-current path. The existing Vela node-local path remains available as `mobility_density_gradient` for `AvalDensGradQF`-like control probes.
- Added parser and behavior/Jacobian coverage in `tests/test_impact_ionization.cpp`:
  - `SG edge-current avalanche approximation cancels flat quasi-Fermi current`
  - `Coupled DD SG edge-current avalanche Jacobian matches carrier finite differences`
- The first implementation incorrectly multiplied edge current density by dual edge length. That treated a boundary flux measure as a volume source measure and caused low-reverse-bias runaway near `-0.255 V`.
- Corrected the source integration to use a diamond-area proxy, `0.5 * edge.length * edge.couple`, with half the edge source assigned to each endpoint. The matching analytic Jacobian source derivative uses the same area scaling.
- After the area fix, `simulation_bv_minus20_sg_edge_current_probe.json` reached `-20.000000000000014 V` with `1098` sweep rows.
- Key current comparison against Sentaurus reference:

```text
bias_V   Vela current A/um        Sentaurus current A      log10(|Vela/ref|)
-0.5     -3.14825115631e-18       -5.34603361759e-18      -0.22996
-2       -9.17718156549e-18       -1.41761303955e-17      -0.18885
-5       -1.82045736574e-17       -2.84265164766e-17      -0.19354
-10      -1.22111684697e-17       -5.45421425346e-17      -0.64997
-13.2    -1.57777607262e-14       -8.38472088807e-17       2.27456
-20      -2.22504825355e-14       -9.10455666344e-16       1.38808
```

SG edge-current multibias field comparison:

- Probe deck: `build-release\reference_tcad\pn2d_sentaurus2018\vela\simulation_bv_minus20_sg_edge_current_vtk_probe.json`.
- Vela VTK root: `build-release\reference_tcad\pn2d_sentaurus2018\vela\sg_edge_current_vtk`.
- Imported Sentaurus `-13.2 V` TDR: `reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_multibias_0132_des.tdr`.
- Report: `build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_sg_edge_current_multibias`.
- Compared fields: `electric_field`, `electron_density`, `hole_density`, and `avalanche_generation`.
- The VTK replay solved exactly the requested bias list: `0`, `-0.5`, `-2`, `-5`, `-10`, `-13.2`, and `-20 V`.
- Current ratios from this replay:

```text
bias_V   Vela current A/um        Sentaurus current A      log10(|Vela/ref|)
-0.5     -3.14825115631e-18       -5.34603361759e-18      -0.22996
-2       -9.17707480763e-18       -1.41761303955e-17      -0.18885
-5       -1.82061598401e-17       -2.84265164766e-17      -0.19351
-10      -2.42465851066e-17       -5.45421425346e-17      -0.35208
-13.2    -2.22562034152e-14       -8.38472088807e-17       2.42396
-20      -2.22504619592e-14       -9.10455666344e-16       1.38808
```

- At `-10 V`, the junction electric-field relative p95 error is `0.147`, electron-density log10 p95 error is `0.467`, hole-density log10 p95 error is `0.485`, and avalanche p99 is lower than Sentaurus (`3.72e14` versus `1.40e15 cm^-3 s^-1`).
- At `-13.2 V`, the junction electric-field relative p95 error remains small (`0.0887`), but electron-density and hole-density log10 p95 errors jump to `3.06` and `1.87`. Avalanche p99 jumps to `6.02e17` versus Sentaurus `3.19e15 cm^-3 s^-1`, matching the `2.42` decade current overshoot.
- At `-20 V`, the junction electric-field relative p95 error is `0.252`, electron-density and hole-density log10 p95 errors are `3.15` and `1.67`, and avalanche p99 is `2.48e18` versus Sentaurus `5.60e16 cm^-3 s^-1`.
- Low-field/contact `electric_field` relative p95 values are dominated by near-zero Sentaurus denominators and should not be interpreted as the main high-field mismatch. The junction-local electric-field metric is the relevant value for the avalanche branch.
- Avalanche peak locations:
  - `-10 V`: Sentaurus `(1.0 um, 0.015625 um)`, Vela `(1.0 um, 0.0 um)`.
  - `-13.2 V`: Sentaurus `(1.0 um, 0.015625 um)`, Vela `(1.0 um, 0.0 um)`.
  - `-20 V`: Sentaurus `(1.0078125 um, 0.046875 um)`, Vela `(0.0 um, 0.5 um)`.

SG source-decomposition diagnostic:

- Script added: `scripts/diagnose_pn2d_bv_sg_avalanche_edges.py`.
- Regression coverage added: `ReferenceTcadToolsTest.test_pn2d_bv_sg_avalanche_edge_diagnostic_decomposes_sources`.
- Report root: `build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_sg_edge_current_source_decomposition`.
- The diagnostic reconstructs SG edge source using:
  - Vela box policy `0.5 * edge.length * edge.couple`.
  - Van Overstraeten-de Man coefficients.
  - Quasi-Fermi-gradient driving force.
  - Quasi-Fermi variable-ni SG flux.
  - `doping.csv + old_slotboom` effective intrinsic density for the PN2D run.
  - VTK endpoint-average mobility; this makes the report a localization diagnostic, not a byte-for-byte replacement for C++ `edgeMobility`.
- Reconstruction against VTK node-integrated source:

```text
bias_V   reconstructed integral   VTK node integral       relative error
-10      1.999579e7               2.032035e7             -1.60%
-13.2    3.952609e10              4.025590e10            -1.81%
-20      3.951536e10              1.176184e11            -66.4%
```

- Edge-class source fractions from reconstructed SG edge source:

```text
bias_V   interior_bulk   boundary_noncontact   contact_adjacent/contact_boundary
-10      96.78%          3.22%                 ~0%
-13.2    96.65%          3.35%                 ~0%
-20      96.65%          3.35%                 ~0%
```

- Top reconstructed SG edge at all three biases: edge `2886`, nodes `351-986`, `(1.0 um, 0.015625 um)` to `(1.0078125 um, 0.015625 um)`, classified as `interior_bulk`.
- VTK node-integrated source fractions:

```text
bias_V   interior nodes   boundary nodes   contact nodes
-10      96.78%           3.22%            ~0%
-13.2    96.63%           3.37%            ~0%
-20      63.90%           3.21%            32.89%
```

- Interpretation of the `-20 V` mismatch: VTK `AvalancheGeneration * node_volume` reports a large contact-node contribution that the SG edge-source reconstruction does not reproduce. Contact-adjacent SG fluxes are near zero under the same quasi-Fermi variable-ni form used by the C++ source helper. Treat the `-20 V` contact-node VTK peak as a diagnostic/output inconsistency until a C++ edge-source dump confirms otherwise.
- Decision: the boundary/contact-edge hypothesis is not supported for the `-13.2 V` current jump. The physically relevant excess source is junction-local/interior, not contact-driven.

Local avalanche-factor comparison around top SG edge:

- Script added: `scripts/diagnose_pn2d_bv_local_avalanche_factors.py`.
- Regression coverage added: `ReferenceTcadToolsTest.test_pn2d_bv_local_avalanche_factor_diagnostic_writes_edge_summary`.
- Report: `build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_local_avalanche_factors`.
- Target edge: Vela SG edge `2886`, nodes `351-986`, midpoint `(1.00390625 um, 0.015625 um)`. Nearest Sentaurus node is `351`, distance `0.00390625 um`.
- The local-factor diagnostic compares Vela edge source density against Sentaurus nearest-node `ImpactIonization`, `ElectricField`, `e/hCurrentDensity`, carrier densities, and mobilities. Sentaurus weighted alpha is inferred as `G / (|Jn|/q + |Jp|/q)`.

```text
bias_V   log10(Vela G / Sentaurus G)   Vela E      Sentaurus E   Vela e-flux     Sentaurus e-flux   Vela eDensity   Sentaurus eDensity
-10      -0.259                        4.025e7     4.022e7       1.838e14        3.264e14           1.877e3         3.230e3
-13.2     2.614                        4.546e7     4.589e7       2.662e17        4.774e14           2.684e6         4.724e3
-20       1.374                        4.546e7     5.609e7       2.662e17        4.189e15           2.683e6         4.183e4
```

- Local alpha and mobility are not the primary `-13.2 V` mismatch:
  - `-13.2 V` Vela electron alpha is `4.67e6 m^-1`; Sentaurus alpha from local electric field is `4.81e6 m^-1`.
  - `-13.2 V` Vela electron mobility is `22.02 cm^2/V/s`; Sentaurus is `22.90 cm^2/V/s`.
  - Hole alpha and mobility are also close enough that they cannot explain a `2.61` decade local source-density error.
- The dominant `-13.2 V` local error is carrier/flux feedback:
  - Vela electron flux is about `5.6e2` times Sentaurus.
  - Vela electron density is about `5.7e2` times Sentaurus.
  - Vela hole flux is about `7.9e1` times Sentaurus.
  - Vela hole density is about `6.6e1` times Sentaurus.
- Decision: do not tune Van Overstraeten coefficients or high-field mobility first. The next root-cause branch is why the Vela continuity solution accumulates much larger minority/majority carrier densities and SG flux at the interior junction edge after `-10 V`.

Continuity-feedback localization around edge `2886`:

- Script added: `scripts/diagnose_pn2d_bv_continuity_feedback.py`.
- Regression coverage added: `ReferenceTcadToolsTest.test_pn2d_bv_continuity_feedback_diagnostic_writes_local_terms`.
- Report: `build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_continuity_feedback`.
- The diagnostic compares the focus edge and incident edges, writing:
  - `continuity_feedback_edges.csv`: same-edge Vela/Sentaurus potential and quasi-Fermi drops, SG flux proxy, source density, and density ratios.
  - `continuity_feedback_nodes.csv`: endpoint/neighbor `psi`, `phin`, `phip`, Boltzmann density reconstruction terms, inferred Sentaurus effective `ni`, Vela SG transport integrals, Vela node avalanche/SRH terms, and Sentaurus node generation/SRH terms.

Focus edge `2886`:

```text
bias_V   Vela dphin   Sentaurus dphin   Vela e-flux      Sentaurus e-flux   log10 G ratio   log10 n ratio   log10 p ratio
-10      0.313210     0.308373          1.838e14         3.374e14           -0.254          -0.248          -0.388
-13.2    0.354690     0.352456          2.662e17         4.954e14            2.616           2.740           1.832
```

- The focus-edge quasi-Fermi drops remain close at `-13.2 V`; the electron `dphin` differs by only `0.00223 V` across the edge. Therefore the local edge slope itself is not the first-order cause of the flux/source jump.
- The endpoint density exponent terms do jump:

```text
bias_V   node   delta(psi-phin) Vela-Sentaurus   delta(phip-psi) Vela-Sentaurus   log10 n ratio   log10 p ratio
-10      351     0.00276                         -0.00587                         -0.246          -0.391
-10      986    -0.00182                         -0.00984                         -0.250          -0.384
-13.2    351     0.18115                          0.12700                          2.750           1.841
-13.2    986     0.17555                          0.12158                          2.730           1.824
```

- Interpretation: by `-13.2 V`, Vela has moved the local Boltzmann density exponents (`psi - phin` for electrons and `phip - psi` for holes) onto a much higher carrier-density branch. Reconstructing Vela densities from its own `psi/phin/phip/ni_eff` matches the VTK densities, while reconstructing with Sentaurus `psi/phin/phip` and Vela `ni_eff` gives only `~2.4e3-3.0e3 cm^-3`, close to the Sentaurus electron density scale. The excessive SG flux is therefore downstream of the absolute local electrostatic/quasi-Fermi state, not a separate edge-current formula error.
- Next root-cause branch: localize why the coupled solve lets `psi-phin` and `phip-psi` move by `~0.12-0.18 V` between `-10 V` and `-13.2 V`. Compare Poisson charge balance and continuity residuals at nodes `351/986` and neighbors, ideally by running the C++ Newton residual probe on the Vela `-10 V` and `-13.2 V` states and then on the Sentaurus `-13.2 V` state through the same Vela residual evaluator.

Newton residual probe on Vela/Sentaurus local states:

- Script added: `scripts/diagnose_pn2d_bv_newton_residual_states.py`.
- Regression coverage added: `ReferenceTcadToolsTest.test_pn2d_bv_newton_residual_state_diagnostic_prepares_probe_inputs`.
- Report: `build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_newton_residual_states`.
- Evaluated these external states with the same release Vela Newton residual evaluator and the SG edge-current avalanche deck:
  - Vela VTK state at `-10 V`.
  - Vela VTK state at `-13.2 V`.
  - Sentaurus exported state at `-13.2 V`.

Block residuals:

```text
state              psi block      phin block      phip block      combined
Vela -10 V         8.9937e-2      4.4916e-15     1.2396e-13     8.9937e-2
Vela -13.2 V       2.4589e-1      1.6315e-9      7.5151e-10     2.4589e-1
Sentaurus -13.2 V  5.9905e1      5.3639e-10     1.6636e-12     5.9905e1
```

Focus nodes `351/986`:

```text
state              node   psi residual   phin residual   phip residual
Vela -10 V         351    1.160e-3      -3.987e-18       7.525e-19
Vela -10 V         986   -5.531e-4       6.072e-19      -2.936e-18
Vela -13.2 V       351    3.868e-4      -3.462e-15       2.375e-16
Vela -13.2 V       986   -3.597e-4      -4.322e-15       2.788e-16
Sentaurus -13.2 V  351   -2.853e-14     -1.296e-15      -1.964e-15
Sentaurus -13.2 V  986    5.999e-6       2.159e-15       3.982e-16
```

- Interpretation:
  - The over-dense Vela `-13.2 V` state is locally well-balanced by Vela's own continuity equations at the focus edge endpoints. The high avalanche source is not caused by an unbalanced local continuity residual at nodes `351/986`.
  - The Vela `-13.2 V` global continuity residual is real (`phin ~= 1.63e-9`, `phip ~= 7.52e-10`) but its largest `phin` residual nodes are around `x ~= 0.156-0.188 um`, not at the junction focus edge near `x ~= 1.0 um`.
  - The Sentaurus `-13.2 V` state is locally acceptable around `351/986` in the Vela continuity residual, but globally very incompatible with Vela's Poisson equation (`psi block ~= 59.9`), with top Poisson residual nodes on the right n-side/contact-side plateau around `x ~= 1.53-1.75 um`.
  - Therefore the current evidence favors a branch-selection/global electrostatic-state mismatch, not a local avalanche coefficient, mobility, edge-current, or focus-edge continuity-source imbalance.
- Next root-cause branch: construct hybrid external-state probes (`Vela psi + Sentaurus phin/phip`, `Sentaurus psi + Vela phin/phip`, and possibly Sentaurus-shifted potential gauges) to determine whether the high-density branch is primarily selected by Poisson/electrostatic potential, quasi-Fermi contact anchoring, or a gauge/reference offset between Vela and Sentaurus.

Hybrid external-state residual probes:

- Extended script: `scripts/diagnose_pn2d_bv_newton_residual_states.py` now supports:
  - `hybrid_vpsi_sqf`: Vela `psi` + Sentaurus `phin/phip`.
  - `hybrid_spsi_vqf`: Sentaurus `psi` + Vela `phin/phip`.
  - `hybrid_spsi_shift_vqf`: Sentaurus `psi` shifted to the selected-node Vela mean + Vela `phin/phip`.
- Regression coverage extended in `ReferenceTcadToolsTest.test_pn2d_bv_newton_residual_state_diagnostic_prepares_probe_inputs`.
- Report: `build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_hybrid_residual_states`.

Block residuals at `-13.2 V`:

```text
state                         psi block      phin block      phip block      combined
Vela                          2.4589e-1      1.6315e-9      7.5151e-10     2.4589e-1
Sentaurus                     5.9905e1       5.3639e-10     1.6636e-12     5.9905e1
Vela psi + Sentaurus qF       1.2005e2       8.2326e-10     1.9849e-12     1.2005e2
Sentaurus psi + Vela qF       2.5925e6       1.6365e-9      7.4408e-8      2.5925e6
Sentaurus psi shifted + Vela qF 1.8606e4     6.4908e-9      1.2604e-9      1.8606e4
```

Focus-node density exponent terms:

```text
state                         node   psi-phin      phip-psi
Vela                          351   -0.21282      -0.25607
Vela                          986   -0.21235      -0.25880
Sentaurus                     351   -0.39397      -0.38307
Sentaurus                     986   -0.38790      -0.38039
Vela psi + Sentaurus qF       351   -0.25108      -0.52596
Vela psi + Sentaurus qF       986   -0.24838      -0.51991
Sentaurus psi + Vela qF       351   -0.35571      -0.11318
Sentaurus psi + Vela qF       986   -0.35187      -0.11928
Sentaurus psi shifted + Vela qF 351 -0.21450      -0.25439
Sentaurus psi shifted + Vela qF 986 -0.21067      -0.26048
```

- `hybrid_spsi_shift_vqf` uses a selected-node `psi` shift of `+0.141206055 V`, which nearly reproduces the Vela high-density exponents at nodes `351/986`. However, its global Poisson residual is still very large (`psi block ~= 1.86e4`).
- `hybrid_spsi_vqf` is even more incompatible globally (`psi block ~= 2.59e6`) and gives strongly distorted focus-node exponents.
- `hybrid_vpsi_sqf` reduces the global mismatch relative to Sentaurus-only but still has a larger Poisson residual than either native state and introduces nontrivial local qF residuals at the focus edge.
- Interpretation: the remaining mismatch is not explained by a simple constant electrostatic-potential gauge offset. The high-density branch depends on a coupled, spatially nonuniform electrostatic/quasi-Fermi state. Vela's local high-density branch is internally consistent near the focus edge, but that branch is not Sentaurus-compatible globally.
- Next root-cause branch: compare the spatial potential shape, not just absolute gauge, especially along horizontal/vertical cuts from the left p-side/contact region through the junction to the right n-side/contact region. The top Vela `-13.2 V` continuity residuals occur around `x ~= 0.156-0.188 um`, while Sentaurus-in-Vela Poisson residuals peak near `x ~= 1.53-1.75 um`; the next diagnostic should connect these global electrostatic plateaus to the junction exponent shift.

Full-device potential/quasi-Fermi profile comparison:

- Script added: `scripts/diagnose_pn2d_bv_potential_profiles.py`.
- Regression coverage added: `ReferenceTcadToolsTest.test_pn2d_bv_potential_profile_diagnostic_writes_plateaus`.
- Report: `build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_potential_profiles`.
- The diagnostic writes:
  - `potential_profile_samples.csv`: horizontal/vertical profile samples of `psi`, `phin`, `phip`, `psi-phin`, and `phip-psi`.
  - `potential_plateau_offsets.csv`: median offsets in full-device x-bands.

Plateau median offsets, Vela minus Sentaurus:

```text
bias_V   band             dpsi       dphin      dphip      d(psi-phin)   d(phip-psi)
-10      left_p          -0.0130     0.00015   ~0         -0.0132        0.0130
-10      pre_junction_p  -0.00275   -0.00329   -0.0187     0.00027      -0.0112
-10      junction        ~0         -0.00196   -0.00560   -0.00057      -0.00828
-10      post_junction_n  0.00266    0.0168    -0.00348   -0.00455      -0.00614
-10      right_n          0.0130    ~0         -0.00050    0.0130       -0.0135
-13.2    left_p           0.261      0.269      0.274     -0.00827       0.0130
-13.2    pre_junction_p   0.208      0.00093    0.330      0.204         0.123
-13.2    junction         0.137     -0.0410     0.261      0.176         0.122
-13.2    post_junction_n  0.0665    -0.0887     0.135      0.162         0.0744
-13.2    right_n          0.0130    ~0         -0.00045    0.0130       -0.0135
```

Selected horizontal profile near the junction top edge (`y = 0.015625 um`):

```text
bias_V   x_um      dpsi      dphin     dphip     d(psi-phin)   d(phip-psi)
-10      0.890625  0.00515   0.00403  -0.0125     0.00112      -0.0176
-10      1.007812  0.00561   0.00744  -0.00423   -0.00182      -0.00984
-10      1.101562  0.0116    0.0127    0.00425   -0.00115      -0.00734
-13.2    0.890625  0.1936   -0.00549   0.312      0.199         0.119
-13.2    1.007812  0.1395   -0.0360    0.261      0.176         0.122
-13.2    1.101562  0.1021   -0.0660    0.195      0.168         0.0931
```

- Interpretation:
  - At `-10 V`, Vela/Sentaurus exponent offsets remain near zero around the junction; the current mismatch is still modest.
  - At `-13.2 V`, the left p plateau has a large positive `psi` offset (`~0.261 V`) but `phin/phip` shift with it, so `psi-phin` and `phip-psi` remain near Sentaurus there.
  - The mismatch is created spatially between the pre-junction p side and post-junction n side: `psi` remains elevated while electron qF shifts little or negative and hole qF shifts strongly positive. This produces `d(psi-phin) ~= 0.16-0.20 V` and `d(phip-psi) ~= 0.07-0.12 V`, matching the carrier-density exponent jump at edge `2886`.
  - The right n plateau returns to the same small `~0.013 V` `psi` offset already seen at `-10 V`, so the error is a nonuniform shape/transition-region issue, not a device-wide offset.
- Next root-cause branch: inspect why the high-bias coupled solve reshapes the pre-junction/junction electrostatic transition. Candidate checks are contact Dirichlet/effective Poisson boundary values versus Sentaurus at `-13.2 V`, depletion-region charge/Poisson balance across the p-to-n transition, and whether avalanche feedback shifts `psi` before the edge `2886` source spike.

Poisson boundary, doping, and charge-balance comparison:

- Script added: `scripts/diagnose_pn2d_bv_poisson_boundary_charge.py`.
- Regression coverage added: `ReferenceTcadToolsTest.test_pn2d_bv_poisson_boundary_charge_diagnostic_writes_doping_contact_charge`.
- Report: `build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_poisson_boundary_charge`.
- The diagnostic writes:
  - `doping_distribution_summary.csv`: Vela and Sentaurus donor/acceptor grouping by value and x range.
  - `contact_boundary_summary.csv`: contact-node `psi/phin/phip` medians versus the Vela ohmic built-in boundary expectation.
  - `charge_balance_bands.csv`: band medians for fixed doping, carrier densities, and Poisson charge density `p - n + Nd - Na`.

Doping distribution:

```text
source      donor          acceptor       net           count   x range
Vela        0              1e17           -1e17         955     0.0 to 0.9921875 um
Vela        1e17           0              +1e17         955     1.0078125 to 2.0 um
Vela        1e17           1e17            0            33      x = 1.0 um
Sentaurus   0              1e17           -1e17         955     0.0 to 0.9921875 um
Sentaurus   1e17           0              +1e17         955     1.0078125 to 2.0 um
Sentaurus   1e17           1e17            0            33      x = 1.0 um
```

- This confirms the user's observation: the current PN2D BV case has constant `1e17 cm^-3` p/n doping on the two sides, plus 33 compensated nodes at the junction. Vela and Sentaurus imported donor/acceptor fields match to floating-point roundoff.

Contact-boundary medians:

```text
bias    contact   expected Vela psi   Vela psi    Sentaurus psi   Vela-Sentaurus
-10     Cathode    0.416685           0.416685     0.403651       +0.013034
-10     Anode    -10.416514         -10.416700   -10.403651       -0.013049
-13.2   Cathode    0.416685           0.416685     0.403651       +0.013034
-13.2   Anode    -13.616514         -13.616700   -13.603651       -0.013049
```

- Vela contact nodes match its own ohmic built-in expectation nearly exactly. Sentaurus contact-node potentials differ by a stable `~13 mV` built-in convention offset at both biases. This offset is much smaller than the `-13.2 V` interior potential-shape mismatch (`~0.07-0.26 V`), so contact Dirichlet anchoring is not the first-order cause.

Charge/density bands:

```text
bias    band             log10(n V/S)   log10(p V/S)   median charge delta cm^-3
-10     junction         -0.260         -0.387         -1.20e3
-13.2   pre_junction_p   +3.200         +1.846         -1.17e6
-13.2   junction         +2.697         +1.794         -1.92e6
-13.2   post_junction_n  +2.501         +1.032         -2.75e6
```

- Fixed net-doping medians differ only by roundoff (`~16 cm^-3` on `1e17 cm^-3` doped bands, zero at the compensated junction band).
- The large `-13.2 V` carrier-density branch mismatch is reproduced in the charge table, but in doped pre/post-junction bands the Poisson charge remains fixed-doping dominated. At the compensated junction, Vela's median charge is about `-1.91e6 cm^-3` versus Sentaurus `+2.48e3 cm^-3`, matching the higher carrier-density branch but still far below the fixed `1e17 cm^-3` side doping scale.
- Interpretation: the current evidence does not support a wrong imported doping distribution or a contact boundary sign/value mistake. The next root-cause branch should compare Vela's Poisson discretization/control-volume flux balance against the same fixed doping and contact boundary conditions, preferably with a no-avalanche or frozen-carrier Poisson reconstruction to separate electrostatic depletion-shape error from avalanche feedback.

Poisson control-volume flux-balance comparison:

- Script added: `scripts/diagnose_pn2d_bv_poisson_flux_balance.py`.
- Regression coverage added: `ReferenceTcadToolsTest.test_pn2d_bv_poisson_flux_balance_diagnostic_writes_top_nodes`.
- Report: `build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_poisson_flux_balance`.
- The diagnostic reconstructs the Vela-style Poisson balance:

```text
edge flux term = sum_edges eps * couple / length * (psi_i - psi_j)
charge term    = q * (p - n + Nd - Na) * node_control_volume
residual       = edge flux term - charge term
```

- Dirichlet contact nodes are excluded from the default ranking and band aggregates because the C++ residual replaces those rows with boundary-condition rows.
- Main CSV outputs:
  - `poisson_flux_balance_top_nodes.csv`: top interior nodes by absolute reconstructed Poisson residual.
  - `poisson_flux_balance_bands.csv`: Vela/Sentaurus band medians for flux term, charge term, and residual.
  - `poisson_flux_balance_band_compare.csv`: Vela minus Sentaurus residual deltas by bias and x band.

Band median comparison using the Vela-style control volumes:

```text
bias    band             Vela median residual/eps   Sentaurus median residual/eps
-10     pre_junction_p   -7.01e-7                  -9.25e-7
-10     post_junction_n  +1.50e-7                  +9.25e-7
-13.2   left_p           +2.50e-6                  ~0
-13.2   pre_junction_p   -3.35e-6                  -1.24e-6
-13.2   junction         +3.51e-7                  ~0
-13.2   post_junction_n  +1.67e-7                  +1.24e-6
```

Top interior reconstructed residuals:

```text
state       bias    top region        top node/x                         residual/eps
Vela        -10     pre_junction_p    node 244, x=0.78125 um             -2.35e-4
Sentaurus   -10     post/pre side     nodes 1089 x=1.125, 955 x=0.875    +/-1.73e-2
Vela        -13.2   pre_junction_p    node 297, x=0.8984375 um           -3.21e-4
Sentaurus   -13.2   post/pre side     nodes 1089 x=1.125, 955 x=0.875    +/-1.73e-2
```

- Interpretation:
  - Vela's own `-13.2 V` state is internally close to balanced under the Vela-style Poisson control-volume reconstruction; its largest interior residual scale remains `~3e-4` after excluding Dirichlet contact nodes.
  - The Sentaurus electrostatic state is much less compatible with the Vela-style Poisson control volumes at the abrupt junction shoulders (`x ~= 0.875 um` and `x ~= 1.125 um`), with reconstructed residual/eps around `1.7e-2` at both `-10 V` and `-13.2 V`.
  - This agrees with the earlier C++ residual probe where the Sentaurus state had a large global Vela Poisson residual while the focus-edge continuity residuals were locally acceptable.
  - The remaining branch is therefore not a constant doping or contact-boundary error. It is more likely a Poisson discretization/control-volume parity issue near the abrupt junction, possibly involving Sentaurus box-volume/truncated-volume handling versus Vela's current exported mesh control volumes/couplings. Avalanche feedback then amplifies the resulting electrostatic-shape difference into the high-density branch near `-13.2 V`.
  - Next root-cause task: export or reconstruct the exact Vela C++ node volumes and edge couplings and compare them against Sentaurus/TDR element geometry around nodes `955`, `1089`, and the focus edge neighborhood. In parallel, run a no-avalanche or frozen-carrier Poisson reconstruction on the same control volumes to separate the base electrostatic discretization gap from avalanche-amplified carrier feedback.

Junction geometry and box-coefficient comparison:

- Script added: `scripts/diagnose_pn2d_bv_junction_geometry.py`.
- Regression coverage added: `ReferenceTcadToolsTest.test_pn2d_bv_junction_geometry_diagnostic_writes_box_terms`.
- Report: `build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_junction_geometry`.
- The diagnostic writes:
  - `junction_geometry_nodes.csv`: selected-node barycentric control volumes, incident edge counts/coupling sums, contact names, and net doping.
  - `junction_geometry_edges.csv`: selected/incident edge length, cotangent coupling, local contribution/fallback counts, contact and junction-touching classification.
  - `junction_geometry_cells.csv`: selected-cell area and triangle angle diagnostics.
  - `junction_geometry_bands.csv`: geometry aggregates over pre-shoulder, compensated junction, and post-shoulder x windows.

Geometry import check:

```text
TDR inventory vertices:       1943
Vela mesh nodes:              1943
TDR inventory triangles:      3680
Vela mesh cells:              3680
max coordinate delta:         0.0 um
triangle sets match:          true
selected negative cot edges:  0
selected fallback edges:      0
```

Selected band geometry:

```text
band             node_count   total volume m2    median node volume m2   median edge couple m   median couple/length
pre_shoulder     166          4.4921875e-14      1.220703125e-16         1.1048543456e-08       0.5
junction         429          5.0781250e-14      1.220703125e-16         7.8125e-09             0.5
post_shoulder    166          4.4921875e-14      1.220703125e-16         1.1048543456e-08       0.5
```

Focus nodes from the Sentaurus-in-Vela Poisson residual ranking:

```text
node   x_um    y_um   net doping cm^-3   volume um^2       adjacent cells   incident edges   incident couple sum um
955    0.875   0.5    -1e17              4.475911e-4      4                5                0.0501110
1089   1.125   0.0    +1e17              4.475911e-4      4                5                0.0501110
```

- Nodes `955` and `1089` are top/bottom boundary shoulder nodes, not interior focus-edge nodes. They are geometrically symmetric p/n shoulder nodes and are exactly where the Sentaurus state has the largest reconstructed residual under the Vela Poisson control volumes.
- Their incident-edge lists include one zero-coupled diagonal edge each, but this occurs without negative cotangent or fallback; it is the ordinary cotangent result for the local right-triangle pattern, not a bad/degenerate cell.
- The compensated junction-touching edge sequence around the focus edge is regular. For example edge `2886` connects node `351` at `x=1.0 um` to node `986` at `x=1.0078125 um`, length `7.8125e-09 m`, and is part of the zero-net junction to n-side transition.
- Interpretation:
  - The Vela mesh import did not distort the TDR coordinates or connectivity.
  - The Vela box geometry near the junction shoulders is regular and symmetric, with no obtuse fallback path.
  - The mismatch is now localized more specifically to boundary-shoulder electrostatic box balance under Sentaurus's potential field. That points toward boundary box/Neumann treatment or Sentaurus's internal box-volume convention near the abrupt junction shoulders, rather than a malformed mesh, wrong doping, or a Vela negative-cotangent artifact.
  - Next root-cause task: run a no-avalanche or frozen-carrier Poisson reconstruction using the Vela control volumes and the same contact boundaries. If the boundary-shoulder shape gap appears without avalanche, focus on Poisson/box boundary parity. If it appears only with avalanche-coupled carriers, focus on the feedback path from continuity/impact-ionization into the boundary shoulder depletion shape.

Frozen-carrier Poisson reconstruction:

- Script added: `scripts/diagnose_pn2d_bv_poisson_reconstruction.py`.
- Regression coverage added: `ReferenceTcadToolsTest.test_pn2d_bv_poisson_reconstruction_diagnostic_solves_frozen_states`.
- Report: `build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_poisson_reconstruction`.
- The diagnostic solves the Vela-style linear Poisson equation on the imported mesh with:
  - `depletion`: fixed donor/acceptor charge only.
  - `vela_frozen`: Vela VTK `n,p` plus fixed doping.
  - `sentaurus_frozen`: Sentaurus exported `n,p` plus fixed doping.
  - Contact BC sources `vela_expected` and `sentaurus_state`.

Important caution:

- `depletion` is not a useful physical baseline for this whole-device reverse-bias PN diode. Removing majority carriers from the quasi-neutral p/n regions creates enormous potential overshoot (`~5-16 V` away from the actual fields), so the meaningful reconstructions are the frozen-carrier cases.

Key `-13.2 V` band results:

```text
bc source        charge source       match target        typical median error
vela_expected    vela_frozen         Vela psi            < 1.5e-4 V
sentaurus_state  sentaurus_frozen    Sentaurus psi       ~0 to 2.2e-3 V
vela_expected    sentaurus_frozen    Sentaurus psi       ~0 to 1.2e-2 V
sentaurus_state  vela_frozen         Vela psi            ~0 to 1.1e-2 V
```

Selected `-13.2 V` focus nodes:

```text
charge source       node   reconstructed - Vela   reconstructed - Sentaurus
vela_frozen         955       +1.4e-4 V             +0.171 V   (Vela BC)
vela_frozen         1089      +7.1e-5 V             +0.103 V   (Vela BC)
sentaurus_frozen    955       -0.192 V              -0.020 V   (Sentaurus BC)
sentaurus_frozen    1089      -0.082 V              +0.021 V   (Sentaurus BC)
sentaurus_frozen    351       -0.136 V              +0.0067 V  (Sentaurus BC)
sentaurus_frozen    986       -0.133 V              +0.0066 V  (Sentaurus BC)
```

- Interpretation:
  - Vela frozen carriers reconstruct the Vela electrostatic branch under the Vela Poisson matrix to within numerical/reporting tolerance.
  - Sentaurus frozen carriers reconstruct the Sentaurus electrostatic branch under the same Vela-style Poisson matrix to within a few millivolts in band medians and about `~20 mV` at the shoulder focus nodes.
  - Switching only the contact BC source mostly adds the already-known `~13 mV` built-in offset; it does not move the solution between the Vela and Sentaurus high-bias branches.
  - Therefore the high-bias shape gap is not primarily caused by an inability of Vela's Poisson matrix/control volumes to represent Sentaurus's potential. The decisive variable is the frozen carrier charge distribution selected by the coupled continuity/avalanche solve.
  - This supersedes the earlier suspicion that boundary-shoulder Poisson box geometry was the first-order root cause. The shoulder residuals are useful localization markers, but frozen-carrier reconstruction shows they are downstream of the carrier branch, not proof of a malformed Poisson discretization.
  - Next root-cause task: compare the avalanche-enabled Vela carrier branch against a no-impact-ionization Vela branch and the Sentaurus branch at the same bias. If no-impact Vela stays close to Sentaurus, the avalanche generation/continuity feedback selects the high-density branch. If no-impact Vela already follows the high-density branch, the issue is pre-avalanche continuity or recombination/mobility coupling.

Three-branch avalanche/no-impact/Sentaurus profile comparison:

- Script added: `scripts/diagnose_pn2d_bv_branch_profiles.py`.
- Regression coverage added: `ReferenceTcadToolsTest.test_pn2d_bv_branch_profile_diagnostic_compares_three_states`.
- Report: `build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_branch_profiles`.
- Inputs:
  - Avalanche-enabled Vela: `build-release\reference_tcad\pn2d_sentaurus2018\vela\sg_edge_current_vtk`.
  - No-impact Vela: `build-release\reference_tcad\pn2d_sentaurus2018\vela` root VTK files from `simulation_bv_minus20_no_impact_probe.json`, with `impact_ionization.model = none`.
  - Sentaurus multibias exports: `build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias`.

Key `-13.2 V` median comparison:

```text
band             dpsi av-sent   dpsi noimp-sent   dpsi av-noimp   log10 e av/sent   log10 e noimp/sent   log10 h av/sent   log10 h noimp/sent
left_p           +0.26095       +0.26135          -0.00040        -0.36            -0.33                 +0.00            +0.00
pre_junction_p   +0.19487       +0.19537          -0.00050        +3.07            +3.11                 +1.80            +2.31
junction         +0.13709       +0.13744          -0.00035        +2.70            +2.84                 +1.79            +2.33
post_junction_n  +0.07930       +0.07949          -0.00017        +2.45            +2.67                 +1.02            +1.76
right_n          +0.01303       +0.01303          +0.00000        +0.00            +0.00                 -0.45            -0.45
```

Key `-13.2 V` exponent medians:

```text
band             branch      psi-phin V   phip-psi V   log10 e cm^-3   log10 h cm^-3
pre_junction_p   avalanche   -0.22065     -0.24360     6.29            5.91
pre_junction_p   noimpact    -0.21920     -0.21330     6.32            6.42
pre_junction_p   sentaurus   -0.42421     -0.36650     3.09            4.06
junction         avalanche   -0.21600     -0.25887     6.37            5.65
junction         noimpact    -0.20803     -0.22672     6.51            6.19
junction         sentaurus   -0.39162     -0.38041     3.67            3.86
post_junction_n  avalanche   -0.21150     -0.34182     6.45            4.26
post_junction_n  noimpact    -0.19877     -0.29801     6.66            4.99
post_junction_n  sentaurus   -0.37340     -0.41627     3.95            3.23
```

- Interpretation:
  - At `-10 V`, avalanche-enabled and no-impact Vela have indistinguishable electrostatic medians, and both remain close to Sentaurus except for the known contact/built-in offset scale.
  - At `-13.2 V`, avalanche-enabled and no-impact Vela are still almost the same electrostatic branch: median `avalanche - noimpact` potential is only `0.17-0.50 mV` in the junction-side bands and zero in the right contact plateau.
  - The large mismatch is already present without impact ionization. Both Vela branches sit `~0.08-0.20 V` above Sentaurus in the pre-junction/junction/post-junction potential shape and have electron medians `~2.5-3.1` decades above Sentaurus in the same bands.
  - Therefore impact-ionization feedback is not the first selector of the high-density branch. Avalanche changes the already-selected high-bias branch, especially hole density, but the branch is present in the no-impact continuity solve.
  - Next root-cause task: compare the no-impact Vela branch against Sentaurus with avalanche disabled or source terms ignored, focusing on continuity terms that remain active: SRH recombination, high-field mobility/quasi-Fermi driving force, and contact minority-carrier boundary treatment. The immediate test should run no-impact variants that toggle SRH and high-field mobility while preserving the same reverse-bias continuation, then compare the same branch-profile report.

No-impact continuity-coupling variant scan:

- Script added: `scripts/diagnose_pn2d_bv_noimpact_variant_scan.py`.
- Regression coverage added: `ReferenceTcadToolsTest.test_pn2d_bv_noimpact_variant_scan_prepares_isolation_configs`.
- Report: `build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_noimpact_variant_scan`.
- Compared variants:
  - `baseline`: no-impact, SRH, Masetti high-field mobility driven by quasi-Fermi gradient.
  - `low_field_masetti`: no high-field mobility limiting, keep Masetti doping dependence.
  - `electric_field_drive`: keep high-field mobility, but use electric field instead of quasi-Fermi gradient as the mobility driving force.
  - `contact_relax_n_0p10`: enable n-contact minority-electron relaxation above `0.10 V`.
  - `no_srh`: prepared and attempted, but did not finish the `-10 V` point within the initial long scan; treat it as a slow/nonconvergent branch requiring a dedicated solver-stability probe rather than mixing it with the completed variants.

Terminal current at `-13.2 V`:

```text
variant                current A/um
baseline               -1.3749e-14
low_field_masetti      -8.5167e-14
electric_field_drive   -1.3417e-14
contact_relax_n_0p10   -1.3749e-14
```

Key `-13.2 V` branch comparison:

```text
variant                band             dpsi cand-sent   dpsi cand-base   log10 e cand/sent   log10 e cand/base   log10 h cand/sent   log10 h cand/base
baseline               pre_junction_p   +0.19197         +0.00000         +3.01               +0.00               -0.60               +0.00
baseline               junction         +0.13511         +0.00000         +2.56               +0.00               -0.59               +0.00
baseline               post_junction_n  +0.07817         +0.00000         +2.25               +0.00               -0.41               +0.00
low_field_masetti      pre_junction_p   +0.19197         +0.00000         +2.51               -0.44               -1.40               -0.83
low_field_masetti      junction         +0.13511         +0.00000         +1.81               -0.69               -1.72               -1.13
low_field_masetti      post_junction_n  +0.07817         +0.00000         +1.79               -0.39               -1.25               -0.87
electric_field_drive   pre_junction_p   +0.19257         +0.00060         +3.00               -0.01               -0.60               +0.00
electric_field_drive   junction         +0.13554         +0.00043         +2.55               -0.01               -0.59               +0.00
electric_field_drive   post_junction_n  +0.07843         +0.00021         +2.23               -0.01               -0.41               +0.00
contact_relax_n_0p10   all three bands  identical to baseline within reported medians
```

Mobility decomposition secondary check:

- Script run: `scripts/diagnose_pn2d_bv_mobility.py`.
- Report: `build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_mobility`.
- At `-10 V`, the top density-error scale is still only `~0.56-0.63` decades, before the `-13.2 V` branch jump. The mobility top relative errors are significant (`~0.72` electron, `~0.56` hole), and the inferred Sentaurus high-field limiter from Vela low-field mobility is much stronger than Vela's limiter at the median (`electron ~= 0.065`, `hole ~= 0.114`).
- Interpretation:
  - Contact minority-electron relaxation is not implicated; the tested n-contact relaxation variant is byte-for-byte equivalent in the reported medians and terminal current.
  - Changing the high-field mobility driving force from quasi-Fermi gradient to electric field is also not a first-order selector; it changes `-13.2 V` density medians by only about `0.01` decade and potential medians by less than `1 mV`.
  - Removing high-field mobility limiting lowers the local carrier-density medians by `~0.4-0.7` decades for electrons and `~0.8-1.1` decades for holes, but the electrostatic potential branch remains unchanged and the terminal current becomes larger. Mobility therefore modulates the high-density branch after it is selected; it does not explain the nonuniform `psi` branch offset itself.
  - The remaining strongest branch is recombination/lifetime and continuity balance. The `no_srh` variant did not reach `-10 V` in the broad scan, so SRH cannot be dismissed, but its effect must be isolated with a dedicated continuation/stability probe rather than a direct all-variant batch run.

Dedicated SRH lifetime no-impact probe:

- Reused and extended `scripts/diagnose_pn2d_bv_noimpact_variant_scan.py` with explicit SRH lifetime variants.
- Regression coverage updated: `ReferenceTcadToolsTest.test_pn2d_bv_noimpact_variant_scan_prepares_isolation_configs` now checks `taun/taup` emission.
- Report: `build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_noimpact_srh_lifetime_scan`.
- Completed variants over `0 -> -10 -> -13.2 V`:
  - `baseline`: no explicit lifetime override, using solver defaults.
  - `srh_tau_asym_1e_m5_3e_m6`: explicit `taun=1e-5 s`, `taup=3e-6 s`, matching the current solver defaults.
  - `srh_tau_equal_1e_m5`: `taun=taup=1e-5 s`.
  - `srh_tau_equal_1e_m6`: `taun=taup=1e-6 s`.
- Partial/stability variant:
  - `srh_tau_equal_1e_m7` produced `0 V` and `-10 V` VTK states but did not complete the `-13.2 V` point in the long scan.

Terminal current at `-13.2 V` for completed lifetime variants:

```text
variant                      current A/um
baseline                     -1.3749126e-14
srh_tau_asym_1e_m5_3e_m6     -1.3749126e-14
srh_tau_equal_1e_m5          -1.3746413e-14
srh_tau_equal_1e_m6          -1.3791786e-14
```

Key `-13.2 V` lifetime comparison:

```text
variant                      band             dpsi cand-sent   dpsi cand-base   log10 e cand/sent   log10 e cand/base   log10 h cand/sent   log10 h cand/base
baseline                     pre_junction_p   +0.19197         +0.00000         +3.007              +0.000              -0.598              +0.000
baseline                     junction         +0.13511         +0.00000         +2.562              +0.000              -0.591              +0.000
baseline                     post_junction_n  +0.07817         +0.00000         +2.246              +0.000              -0.409              +0.000
srh_tau_equal_1e_m5          pre_junction_p   +0.19192         +0.00000         +3.007              -0.000              -0.784              -0.186
srh_tau_equal_1e_m5          junction         +0.13510         +0.00000         +2.561              -0.000              -0.778              -0.186
srh_tau_equal_1e_m5          post_junction_n  +0.07817         +0.00000         +2.246              -0.000              -0.576              -0.185
srh_tau_equal_1e_m6          pre_junction_p   +0.19202         +0.00010         +3.010              +0.002              +0.205              +0.811
srh_tau_equal_1e_m6          junction         +0.13515         +0.00004         +2.565              +0.004              +0.219              +0.810
srh_tau_equal_1e_m6          post_junction_n  +0.07820         +0.00002         +2.252              +0.005              +0.399              +0.806
```

Partial `srh_tau_equal_1e_m7` at `-10 V`:

```text
band             dpsi cand-sent   dpsi cand-base   log10 e cand/sent   log10 e cand/base   log10 h cand/sent   log10 h cand/base
pre_junction_p   -0.00950         +0.00000         +1.532              +1.796              +1.266              +1.801
junction         ~0.00000         +0.00000         +1.502              +1.808              +1.359              +1.808
post_junction_n  +0.00950         +0.00000         +1.371              +1.795              +1.437              +1.797
```

- Interpretation:
  - Making the solver-default lifetimes explicit reproduces the baseline exactly, so no hidden lifetime default/config mismatch is present in the current no-impact deck.
  - Reasonable lifetime changes from default to equal `1e-5 s` or equal `1e-6 s` do not move the nonuniform electrostatic branch at `-13.2 V`; potential medians remain within about `0.1 mV` of baseline.
  - Electron density, the main driver of the high branch, is nearly insensitive over the stable lifetime range (`<= ~0.005` decade versus baseline at `-13.2 V`). Hole density is more lifetime-sensitive, but that does not explain the electron high-density branch or the electrostatic-shape offset.
  - Very strong recombination (`taun=taup=1e-7 s`) is solver-expensive and already increases carrier medians by about `1.8` decades at `-10 V` while leaving the potential branch unchanged. It is not a route toward Sentaurus parity.
  - Therefore SRH lifetime tuning is not the first-order root cause. The no-impact high-density branch persists after excluding avalanche feedback, contact minority relaxation, mobility driving-force choice, high-field mobility limiting as the selector, and stable SRH lifetime differences.
  - Next root-cause task: inspect the no-impact continuity residual/Jacobian or quasi-Fermi boundary/reference handling between `-10 V` and `-13.2 V`. The remaining evidence points to how Vela's coupled no-impact continuity solve chooses quasi-Fermi absolute levels under the same electrostatic branch, not to local avalanche, mobility, SRH lifetime, or doping/contact-value mismatch.

Interpretation:

- SG edge-current source integration is now the first stable path that reaches -20 V with avalanche enabled.
- It is not ready for promotion because high-reverse-bias avalanche current is still too large by `~2.27` decades near `-13.2 V` and `~1.39` decades at `-20 V`.
- The remaining mismatch is no longer the Newton line-search blocker. It is avalanche source/current-distribution parity at high field. The `-20 V` VTK contact-node peak is now separated from the `-13.2 V` current jump and should be handled as a diagnostic-output/contact-node source consistency issue.

Current recommendation:

- Do not promote a relaxed `abstol` as the -20 V validation solution.
- Treat `density_gradient` as the only Sentaurus-default BV parity path, but do not promote the -20 V gate until the high-field current/source overshoot is localized with a field comparison or source-volume sensitivity probe.
- Treat `mobility_density_gradient` as an `AvalDensGradQF`-like control path. It can help isolate numerical stability and source sensitivity, but it cannot be the acceptance path for a Sentaurus-default BV claim.
- Next best task: construct a no-impact quasi-Fermi/continuity residual diagnostic across `-10 V` and `-13.2 V`, especially contact-adjacent and plateau quasi-Fermi absolute levels. The goal is to determine why the no-impact continuity solve selects much higher carrier-density exponents even when the electrostatic potential branch is unchanged by mobility, contact relaxation, and stable SRH lifetime changes.
- Separately add a C++ edge-source dump or VTK contact-node guard to resolve the `-20 V` contact-node diagnostic inconsistency.

## Task 6: Define The Sentaurus-Default Path And Keep -20 V Blocked

The stable SG edge-current run proves the required Sentaurus-default path is executable, not that BV parity has passed. Keep the committed BV reference gate at low reverse bias until both high-bias current parity and the no-impact branch mismatch are understood.

- [x] **Step 1: Keep the committed BV reference JSON on the low-bias gate**

Confirm `reference_tcad/pn2d_sentaurus2018/pn2d_sentaurus2018_reference.json` still contains the BV controls:

```json
"vela_stop": -0.05,
"vela_step": -0.05
```

Expected: do not set `vela_stop` to `-20.0` in this task. If the committed generated path records Sentaurus-default impact-ionization settings, it must use `current_approximation = "density_gradient"` and clearly remain a low-bias parity gate, not a promoted -20 V validation.

- [x] **Step 2: Mark `density_gradient` as the Sentaurus-default parity path**

Record in this plan, and later in validation docs only after fresh report regeneration, that:

```text
mobility_density_gradient:
  Vela legacy/default control path; closest to Sentaurus Math { AvalDensGradQF }.
  Useful for stability/source-sensitivity comparison.
  Must not be used as the acceptance path for Sentaurus-default BV.

density_gradient:
  Required Sentaurus-default SG edge-current parity path.
  Must be used in all Sentaurus-default BV probes and eventual promotion.
  Reaches -20 V in the current probe, but promotion is blocked by current/source overshoot:
    +2.27 decades near -13.2 V
    +1.39 decades at -20 V
```

Expected: no worker should treat the stable `density_gradient` run as sufficient evidence for promoting the -20 V gate, and no worker should fall back to `mobility_density_gradient` to claim Sentaurus-default BV parity.

- [x] **Step 3: Add a blocker note to `docs/validation/pn2d_sentaurus_comparison.md` only after regenerating release reports**

Append a short "PN2D BV -20 V blocked status" note with:

```text
Validation date:
Vela probe deck:
Sentaurus multibias export root:
Compared biases:
SG edge-current current errors:
No-impact branch mismatch summary:
Why -20 V promotion remains blocked:
Next required diagnostic:
```

Expected: the validation doc explains the blocked state without implying that -20 V BV is calibrated.

### Task 6 Execution Notes 2026-06-18

- Regenerated `build-release\reference_tcad\pn2d_sentaurus2018` with `scripts\sentaurus_import.py reference`.
- Updated `scripts/sentaurus_import.py` and `reference_tcad/pn2d_sentaurus2018/pn2d_sentaurus2018_reference.json` so the Sentaurus-default BV deck records:

```json
"impact_ionization": {
  "model": "van_overstraeten",
  "driving_force": "quasi_fermi_gradient",
  "generation": "current_density",
  "current_approximation": "density_gradient"
}
```

- The committed BV gate remains low-bias only: `vela_stop = -0.05`, `vela_step = -0.05`.
- Added regression assertions in `tests/regression/test_sentaurus_import_tools.py` and `tests/regression/test_reference_tcad_tools.py` so future Sentaurus-default BV imports cannot silently drop `current_approximation = "density_gradient"`.
- Added the blocked-status note to `docs/validation/pn2d_sentaurus_comparison.md`.

## Task 7: No-Impact Quasi-Fermi And Continuity Residual Diagnostic

- [x] **Step 1: Regenerate release no-impact and SG edge-current probe states**

Use the release build and the existing probe-generation scripts to produce same-bias states at `-10 V` and `-13.2 V` for:

```text
Vela no-impact branch
Vela SG edge-current avalanche branch
Sentaurus multibias exported branch
```

Expected: each branch has VTK or imported field data for `psi`, `phin`, `phip`, electron density, hole density, electric field, and enough mesh metadata to evaluate local continuity terms.

- [x] **Step 2: Extend the existing residual/feedback diagnostics only if the current output lacks the required columns**

Prefer reusing:

```text
scripts/diagnose_pn2d_bv_newton_residual_states.py
scripts/diagnose_pn2d_bv_continuity_feedback.py
```

Required report columns:

```text
bias_V
branch
node_id
x_um
y_um
psi_minus_phin_V
phip_minus_psi_V
electron_density_cm3
hole_density_cm3
electron_sg_flux_proxy
hole_sg_flux_proxy
electron_continuity_residual
hole_continuity_residual
contact_or_plateau_band
```

Expected: the report compares `-10 V` and `-13.2 V`, with special focus on contact-adjacent plateaus and the interior junction focus edge around nodes `351` and `986`.

- [x] **Step 3: Classify the no-impact branch selector**

Use the generated residual tables to choose one of these conclusions:

```text
quasi_fermi_anchoring:
  Contact or plateau quasi-Fermi absolute levels explain the density exponent shift.

continuity_residual_balance:
  Local SG flux, recombination, or continuity residual balance explains the branch shift.

spatial_electrostatic_shape:
  A nonuniform psi profile change explains the density exponent shift even when local qF slopes look similar.

unclassified:
  The current diagnostics are insufficient; add the smallest missing observable before proposing a physics fix.
```

Expected: do not tune avalanche, mobility, SRH lifetime, or Newton tolerances until this classification is complete.

### Task 7 Execution Notes 2026-06-18

- Generated diagnostic decks under `build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution`.
- A direct no-impact adaptive run timed out after 240 s at about `-7.75 V`, but it wrote a restart state. A restart deck using that state converged explicit no-impact bias points `-8 V`, `-10 V`, and `-13.2 V`.
- Exported Sentaurus multibias snapshots for `-10 V`, `-13.2 V`, and `-20 V`.
- Ran `scripts/diagnose_pn2d_bv_continuity_feedback.py` for the no-impact Vela branch at `-10 V` and `-13.2 V`.
- Classification: `quasi_fermi_anchoring` / high-bias quasi-Fermi absolute-level branch mismatch. The local focus-edge region is close at `-10 V`, then diverges before avalanche is needed:

```text
-10 V no-impact:
  median log10(Vela/Sentaurus electron density) = -0.297
  median log10(Vela/Sentaurus hole density)     = -0.449
  median delta(psi - phin)                      = -0.00421 V

-13.2 V no-impact:
  median log10(Vela/Sentaurus electron density) = +2.654
  median log10(Vela/Sentaurus hole density)     = -0.592
  median delta(psi - phin)                      = +0.172 V
```

- Therefore the `-13.2 V` high-density branch is already present without impact ionization. Avalanche source tuning is not the first fix; the next implementation task should inspect high-bias quasi-Fermi anchoring/state reconstruction between about `-10 V` and `-13.2 V`.

## Task 8: C++ SG Edge-Source Diagnostic

- [x] **Step 1: Add a runner-visible edge-source dump for the assembled SG avalanche source**

Expose the actual C++ edge-level source terms used by the `density_gradient` assembler path. The dump must include:

```text
bias_V
edge_id
node0
node1
x0_um
y0_um
x1_um
y1_um
edge_length_m
edge_couple_m
edge_area_proxy_m2
electron_alpha_m_inv
hole_alpha_m_inv
electron_flux_proxy
hole_flux_proxy
edge_source_integral
node0_source_integral
node1_source_integral
edge_class
```

Expected: the diagnostic is generated from the same helper path used by the solver, not from a separate Python reconstruction.

- [x] **Step 2: Compare C++ edge dump against the Python SG reconstruction**

Run `scripts/diagnose_pn2d_bv_sg_avalanche_edges.py` on the same state and compare:

```text
top edge id/order
total source integral
contact-edge source fraction
interior-bulk source fraction
node 351/986 source contribution
```

Expected: explain whether the `-20 V` VTK contact-node avalanche peak is an output normalization artifact, a node-volume/reporting artifact, or real assembled source behavior.

- [x] **Step 3: Add focused regression coverage**

Add tests that prove:

```text
flat quasi-Fermi edge current still produces zero SG avalanche source
edge-source dump total equals assembled nodal source integral
VTK AvalancheGeneration remains consistent with the chosen node-volume policy
```

Expected: future diagnostics cannot silently diverge from the solver source assembly.

### Task 8 Execution Notes 2026-06-18

- Added `detail::sgEdgeCurrentAvalancheSourceRecords(...)` in `include/vela/equation/AssemblerUtils.h`; the existing nodal `sgEdgeCurrentAvalancheSourceIntegrals(...)` now sums those records, so the diagnostic and assembler source share the same helper path.
- Added `sweep.diagnostics.sg_avalanche_edges` to `DCSweep`. When enabled with `impact_ionization.generation = "current_density"` and `current_approximation = "density_gradient"`, it writes a separate CSV containing edge id, endpoint coordinates, edge-area proxy, alpha values, mobility values, SG flux proxies, source integrals, node half-source integrals, and edge class.
- Added focused tests:
  - `SG edge-current avalanche records sum to assembled nodal source`
  - `DCSweep: SG avalanche edge diagnostics write assembled source rows`
  - Existing `SG edge-current avalanche approximation cancels flat quasi-Fermi current`
  - `VTK AvalancheGeneration uses SG edge nodal source over node volume`
- Added `scripts/compare_pn2d_bv_sg_edge_source_dump.py` plus regression coverage for comparing C++ edge dumps against Python reconstruction CSVs. The comparison reports top edge ids/order, total source integral, contact-edge fraction, interior-bulk fraction, selected node source contributions, and key deltas.
- Built the release runner with the new diagnostic and generated a 0 V Sentaurus-default SG edge-current C++ dump under `build-release\reference_tcad\pn2d_sentaurus2018\reports\sg_edge_source_cpp_compare\vela_run`.
- Ran `scripts/diagnose_pn2d_bv_sg_avalanche_edges.py` on the matching 0 V VTK and compared it with the C++ dump:
  - Summary JSON: `build-release\reference_tcad\pn2d_sentaurus2018\reports\sg_edge_source_cpp_compare\cpp_vs_python_0v_summary.json`
  - `log10(C++ total / Python total) = -5.3510817401023306e-05`
  - C++ interior-bulk fraction `0.9771051921163618`; Python interior-bulk fraction `0.9771030786581933`
  - C++ contact-edge fraction `4.404660969341172e-07`; Python contact-edge fraction `4.885308332737057e-07`
  - Top edge sets agree, with one near-tie order swap between edge ids `3474` and `4028`.
- High-bias comparison remains open. A direct `0 -> -13.2 V` density-gradient diagnostic deck wrote the 0 V dump/VTK and then exceeded a 300 s foreground timeout while solving the high-bias point. A sparse ramp deck (`0,-1,-2,-4,-8,-10,-12,-13.2`) similarly wrote the 0 V dump/VTK but did not finish the next point within the short probe timeout. Do not use this 0 V consistency result to classify the `-20 V` contact-node avalanche peak.
- Remaining Task 8 work: run a true continuation `density_gradient` pn2d high-bias probe long enough to emit C++ SG edge-source rows at `-13.2 V` and `-20 V`, then compare those rows to the Python reconstruction and classify the VTK contact-node peak as output normalization, node-volume/reporting, or real assembled source behavior.

### Task 8 Additional Execution Notes 2026-06-18

- Used the converged no-impact `-13.2 V` state as the initial state for a Sentaurus-default SG edge-current `-13.2 V` single-point solve. It converged and emitted C++ edge-source rows.
- Used the converged SG `-13.2 V` state as the initial state for a Sentaurus-default SG edge-current `-20 V` single-point solve. It converged and emitted C++ edge-source rows.
- Full-edge C++ vs Python reconstruction results:

```text
-13.2 V:
  log10(total C++ source / Python source) = +0.00804
  C++ interior-bulk source fraction       = 0.966401
  Python interior-bulk source fraction    = 0.966459
  C++ contact-edge source fraction        = 2.68e-6

-20 V:
  log10(total C++ source / Python source) = +0.34695
  C++ interior-bulk source fraction       = 0.442780
  Python interior-bulk source fraction    = 0.966522
  C++ contact-edge source fraction        = 0.541857
```

- Interpretation: the assembled SG source and independent Python reconstruction agree well at `-13.2 V`, including focus node source integrals near nodes `351` and `986`. At `-20 V`, the total source remains within a factor of about `2.2`, but the C++ diagnostic assigns over half the source to contact edges while the Python reconstruction classifies almost all source as interior bulk. The `-20 V` contact-edge/source-reporting path remains open and must be resolved before promotion.

## Task 9: High-Bias Branch Acceptance Gate

- [x] **Step 1: Add a local high-bias gate before any -20 V promotion**

At `-13.2 V`, require the focus-edge local source and carrier feedback to satisfy:

```text
abs(log10(Vela_G / Sentaurus_G)) < 0.5 decades
abs(log10(Vela_electron_density / Sentaurus_electron_density)) < 0.5 decades
abs(log10(Vela_electron_flux / Sentaurus_electron_flux)) < 0.5 decades
```

Use the edge near nodes `351-986` unless a regenerated diagnostic identifies a different dominant interior-bulk edge.

Expected: this is the first promotion gate. Stricter final tolerances can be set after the root cause is known.

- [x] **Step 2: Keep -20 V reference promotion behind the gate**

Only after Task 7 classification and Task 8 source consistency pass, update the committed BV config:

```json
"vela_stop": -20.0,
"vela_step": -0.1
```

The promoted Sentaurus-default impact-ionization config must include:

```json
"current_approximation": "density_gradient"
```

If SG edge-source parity is not confirmed, fix the SG edge-current path, source-volume policy, or high-bias branch selection before promotion. Do not switch the acceptance path to `mobility_density_gradient` to obtain a more stable or closer-looking curve.

Expected: the final validation must compare Sentaurus's implicit default avalanche run against Vela's documented SG edge-current implementation, not against a Vela `AvalDensGradQF`-like implementation.

### Task 9 Execution Notes 2026-06-18

- Ran `scripts/diagnose_pn2d_bv_continuity_feedback.py` on the Sentaurus-default SG VTK states at `-13.2 V` and `-20 V`.
- The promotion gate fails:

```text
-13.2 V SG branch:
  median log10(Vela/Sentaurus electron density) = +2.654
  focus-edge log10(Vela/Sentaurus generation)   = +2.511

-20 V SG branch:
  median log10(Vela/Sentaurus electron density) = +1.787
  median log10(Vela/Sentaurus hole density)     = +0.736
  focus-edge log10(Vela/Sentaurus generation)   = +1.382
```

- Required gate was `< 0.5 decades` for local source, electron density, and electron flux. Therefore `vela_stop = -20.0` remains blocked. Do not promote the generated reference JSON beyond the low-bias gate.

## Task 10: Ionization-Integral Diagnostic Scope Lock

Ionization-integral breakdown analysis is a useful research diagnostic from the Sentaurus, Silvaco, GSS, and Sze references, but it is not the acceptance path for this Sentaurus-default BV reproduction. Use it to understand field-line breakdown propensity after SG edge-current parity is under control.

- [x] **Step 1: Keep ionization integral out of the -20 V promotion gate**

Record in the validation notes that the active acceptance target remains:

```text
Sentaurus default:
  Recombination(Avalanche(VanOverstraeten))
  isothermal GradQuasiFermi driving force
  SG edge-current avalanche source
  no explicit AvalDensGradQF
  no RefDens_*Aval interpolation
  no ElementVolumeAvalanche unless the source deck enables it

Vela required counterpart:
  model = van_overstraeten
  driving_force = quasi_fermi_gradient
  generation = current_density
  current_approximation = density_gradient
  driving_force_interpolation = none
```

Expected: do not replace the self-consistent DD avalanche comparison with an ionization-integral `integral >= 1` criterion when claiming Sentaurus-default BV parity.

### Task 10 Execution Notes 2026-06-18

- The active acceptance path remains self-consistent drift-diffusion with Sentaurus-default SG edge-current avalanche.
- Ionization integral remains a future post-processing diagnostic only; it was not used to pass or fail the `-20 V` promotion gate in this execution.

- [x] **Step 2: Add ionization-integral output only as a post-processing diagnostic after Task 9**

If Task 9 still leaves a high-field source mismatch after the SG edge source is internally consistent, create `scripts/diagnose_pn2d_bv_ionization_integral.py` with these required outputs:

```text
bias_V
path_id
carrier
start_x_um
start_y_um
end_x_um
end_y_um
max_field_V_per_m
integral_alpha_dx
dominant_material
dominant_edge_or_cell_id
```

Expected: the script reads existing Vela and imported Sentaurus field snapshots, reports the largest electron and hole ionization integrals at each bias, and writes a summary JSON. It must not change solver behavior or sweep stopping criteria in this plan.

### Task 10 Additional Execution Notes 2026-06-18

- Added `scripts/diagnose_pn2d_bv_ionization_integral.py`.
- Added regression coverage:
  - `ReferenceTcadToolsTest.test_pn2d_bv_ionization_integral_diagnostic_writes_required_outputs`
- The diagnostic is an edge-local path proxy, not a field-line integrator. It reports the dominant mesh-edge `alpha * dx` values and uses `--field-scale 100` by default because Vela VTK high-field scalars are written in `V/cm` scale while the Van Overstraeten coefficients use `V/m`.
- Post-processed Sentaurus-default SG VTK states:

```text
-13.2 V:
  max edge-local integral = 0.08193
  electron max            = 0.08193
  hole max                = 0.02804

-20 V:
  max edge-local integral = 1.80364
  electron max            = 1.80364
  hole max                = 1.40371
```

- These values support using ionization integrals as a future explanatory diagnostic, but they do not override the failed self-consistent DD high-bias parity gate in Task 9.

## Task 11: No-Impact High-Bias Branch Probe

The next blocker after Task 10 is the no-impact high-bias branch mismatch. Test
continuation history, Gummel handoff, and material intrinsic-density alignment
before changing avalanche physics or promoting the `-20 V` gate.

- [x] **Step 1: Re-run no-impact from `-10 V` to `-13.2 V` with small continuation steps**

Generated a restart CSV from the converged no-impact `-10 V` VTK state and ran:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\noimpact_smallstep_from_minus10
```

The probe used `impact_ionization.model = none`, `step = -0.05 V`, no explicit
`bias_points`, and reached `-13.2 V` in 65 points.

Result: the high-density branch still appears. The focus-edge median
`log10(Vela/Sentaurus electron density)` at `-13.2 V` moves only from `+2.654`
in the large-jump restart to `+2.597` with small-step continuation. Therefore
the previous `-10 V -> -13.2 V` explicit-bias jump is not the first-order root
cause.

- [x] **Step 2: Test whether a real Gummel initializer avoids the high branch**

Generated a restart from the pre-jump `-12.7 V` no-impact state and ran:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\noimpact_gummel_handoff_from_minus12p7
```

The probe enabled `handoff.gummel_max_iter = 30` before Newton handoff.

Result: the run did not reach `-13.2 V`; it repeatedly shrank near
`-12.9466659 V` and then failed Newton. The completed points already show the
same high-density tendency near `-12.9 V`, so Gummel handoff changes the
failure mode but does not restore Sentaurus-like carrier densities.

- [x] **Step 3: Separate material `ni` parity from the high-bias branch selector**

Generated a no-impact small-step probe with a local material override:

```json
{"materials": [{"name": "Si", "ni": 1.6556153e10}]}
```

Because the deck uses `scaling.mode = unit_scaling`, this is interpreted as
`1.6556153e10 cm^-3`.

Result:

```text
-10 V default ni:
  median log10(Vela/Sentaurus electron density) = -0.297

-10 V ni override:
  median log10(Vela/Sentaurus electron density) = -0.072

-13.2 V default ni:
  median log10(Vela/Sentaurus electron density) = +2.597

-13.2 V ni override:
  median log10(Vela/Sentaurus electron density) = +2.995
```

The material `ni` override is a real low-bias calibration axis, but it does not
fix the high-bias no-impact branch. In this probe it worsens the high-bias
electron-density overshoot.

### Task 11 Execution Notes 2026-06-18

- The high-bias branch starts to depart around `-12.75 V` in the no-impact
  small-step run. The terminal current jumps from the `~1e-11 A/m` leakage
  scale into the `~1e-9` to `~1e-8 A/m` electron-current branch before the
  `-13.2 V` point.
- Contact Dirichlet quasi-Fermi values are not the immediate cause. At both
  `-10 V` and `-13.2 V`, contact `phin` and `phip` match the electrode biases
  to numerical precision; the known contact potential difference is the
  material-`ni` built-in offset of about `13 mV`.
- The `-13.2 V` focus-edge density error is mostly driven by the internal
  electrostatic and carrier branch: Vela `psi` is about `0.14 V` above
  Sentaurus near the junction, while `phin` differs by only about `0.03 V`.
- Updated root-cause priority: do not tune avalanche, SRH lifetime, high-field
  mobility, contact minority relaxation, or material `ni` as the next fix for
  the `-13.2 V` branch. The remaining target is the no-impact coupled
  continuity/Poisson branch selection around `-12.7 V` to `-13.0 V`, especially
  the electron-continuity residual/Jacobian balance that permits the
  high-density internal solution.
- The committed BV reference remains blocked at low bias. Do not promote
  `vela_stop = -20.0` until the no-impact branch gate at `-13.2 V` is below
  `0.5 decades`.

## Task 12: No-Impact Branch Residual Spectrum

The next step after Task 11 was to test whether the high-density branch is a
false convergence caused by residual scaling or a self-consistent alternate
DD/Poisson branch.

- [x] **Step 1: Allow Vela-only Newton residual probe preparation**

Updated `scripts/diagnose_pn2d_bv_newton_residual_states.py` so pure
`vela:<bias>` states do not require a same-bias Sentaurus field export. Hybrid
and Sentaurus states still require the Sentaurus export.

Regression coverage added:

```text
ReferenceTcadToolsTest.test_pn2d_bv_newton_residual_state_diagnostic_allows_vela_without_sentaurus
```

- [x] **Step 2: Run high-precision residual probes across the no-impact branch
  jump**

Generated high-precision target states by re-running the no-impact small-step
deck from the same `-10 V` restart to each target bias. This avoids using
truncated VTK values as Newton residual input.

Report root:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\noimpact_branch_high_precision_targets
```

Summary:

| bias (V) | electron current (A/m) | `psi` block | `phin` block | max `phin` node |
|---:|---:|---:|---:|---:|
| `-12.65` | `-1.785e-11` | `8.83e-12` | `5.86e-12` | `2` |
| `-12.70` | `-4.475e-11` | `3.28e-09` | `6.33e-12` | `2` |
| `-12.75` | `-6.424e-10` | `1.17e-09` | `1.01e-11` | `202` |
| `-12.80` | `-3.124e-09` | `1.01e-10` | `1.14e-11` | `1243` |
| `-12.85` | `-9.001e-09` | `2.24e-08` | `1.99e-11` | `201` |
| `-12.90` | `-1.389e-08` | `1.17e-11` | `6.72e-12` | `2` |
| `-12.95` | `-1.061e-09` | `3.82e-10` | `1.10e-10` | `0` |
| `-13.00` | `-1.375e-08` | `5.67e-11` | `6.80e-12` | `2` |
| `-13.20` | `-1.375e-08` | `1.37e-08` | `7.27e-12` | `2` |

The branch jump is visible in the state itself. At focus nodes `351/986`,
`log10(electrons_m^-3)` moves from about `9.88` at `-12.70 V` to about
`11.14` at `-12.75 V`, then saturates around `12.29` after `-12.90 V`.

- [x] **Step 3: Interpret the branch selector**

The high-density branch is not a residual-threshold artifact. In high-precision
states, the accepted no-impact branch has electron and hole continuity block
norms near `1e-11`, including after the current jump. The VTK-based probe had
shown `~0.1` Poisson block residuals, but converting the high-precision
`latest_state.csv` for `-13.2 V` gives `psi ~= 1.37e-8`, confirming that the
VTK residual was dominated by output precision.

Updated next research direction:

- Do not spend the next iteration on residual tolerances or local focus-edge
  source balancing; the high branch is internally residual-balanced.
- Add an accepted-step Newton history/Jacobian diagnostic around
  `-12.70 V -> -12.75 V`, recording per-iteration block residuals, line-search
  damping, update norms, and optionally the smallest/most singular coupled
  Jacobian directions.
- Test branch selection controls rather than physical knobs: smaller target
  steps below `0.05 V`, pseudo-transient/continuation damping, and a
  Sentaurus-like extrapolation policy for the coupled variables.
- If adding a solver change, gate it on recovering the low-density no-impact
  branch at `-13.2 V` without degrading the existing `-10 V` parity.

## Task 13: Accepted-Step Newton History At The No-Impact Branch Jump

Task 12 showed that the high-density no-impact branch is internally
residual-balanced. The next step was to inspect how accepted Newton iterations
enter that branch around `-12.70 V -> -12.75 V`.

- [x] **Step 1: Add accepted-step Newton history diagnostics**

Added `sweep.diagnostics.newton_history` to `DCSweep`. When enabled, it writes a
separate CSV for converged sweep points with:

```text
point_index
bias_V
bias_contact
solver_method
handoff_stage
iteration
residual_norm
relative_residual_norm
raw_step_norm
applied_step_norm
damping_factor
line_search_attempts
line_search_accepted
block_psi
block_phin
block_phip
block_combined
```

Also extended `NewtonIterationInfo` so each accepted Newton iteration records
the post-step `psi`/`phin`/`phip`/combined block residuals. The diagnostic is
off by default and does not change solver behavior.

Regression coverage added:

```text
DCSweep: Newton history diagnostic writes accepted iteration block residuals
```

- [x] **Step 2: Run the no-impact `-12.70 V -> -12.75 V` jump with history**

Report root:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\noimpact_newton_history_m12p70_to_m12p75
```

The `-0.05 V` step from the high-precision `-12.70 V` restart reached
`-12.75 V` with no retry and no line-search damping. At `-12.75 V`, Newton used
five full accepted iterations. The dominant residual was Poisson, not electron
continuity:

| iter | residual | raw step | damping | `psi` block | `phin` block |
|---:|---:|---:|---:|---:|---:|
| 1 | `5.00e-1` | `4.52e1` | `1` | `3.81e0` | `1.95e-10` |
| 2 | `2.63e-4` | `6.88e1` | `1` | `2.00e-3` | `9.49e-13` |
| 3 | `2.61e-4` | `5.00e0` | `1` | `1.99e-3` | `9.45e-13` |
| 4 | `1.07e-4` | `1.91e1` | `1` | `8.18e-4` | `5.50e-13` |
| 5 | `2.16e-9` | `7.22e0` | `1` | `1.64e-8` | `1.38e-13` |

Interpretation: the branch jump is not caused by line-search backtracking or a
large electron-continuity residual. Newton accepts full coupled updates, and
the accepted path is controlled by the electrostatic block.

- [x] **Step 3: Test smaller continuation step size**

Report roots:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\noimpact_newton_history_m12p70_to_m12p75_step005
build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\noimpact_newton_history_m12p70_to_m13p20_step005
```

Using `step = -0.005 V` from the same `-12.70 V` restart:

```text
-12.75 V: focus log10(electrons_m^-3) ~= 10.71
-12.77 V: first focus log10(electrons_m^-3) >= 11
-12.84 V: first focus log10(electrons_m^-3) >= 12
-13.20 V: focus log10(electrons_m^-3) ~= 12.34
```

The smaller step delays the high-density transition but does not recover the
Sentaurus-like low-density branch at `-13.2 V`. Across the `-0.005 V` run,
line-search damping remained `1.0`, maximum Newton iterations per point were
`3`, and the accepted final residuals were again dominated by Poisson block
residuals.

Updated next research direction:

- Smaller continuation steps are a useful diagnostic but not the final fix.
- The next solver experiment should target branch control in the electrostatic
  Newton update: pseudo-transient continuation, a trust-region/max-update
  policy scaled by block or physical field, or Sentaurus-like extrapolation
  controls for the coupled variables.
- Add any branch-control experiment behind an opt-in solver setting and gate it
  on reducing the no-impact `-13.2 V` focus electron-density error without
  degrading the existing `-10 V` parity.

## Task 14: Existing Trust-Region Proxy Scan With `max_update`

Task 13 showed that line search accepts full Newton updates and that smaller
continuation steps only delay the high-density transition. The next experiment
was to test whether the existing Newton `max_update` cap can act as a
trust-region proxy for branch control.

- [x] **Step 1: Scan `max_update` on the no-impact `-12.70 V -> -13.20 V`
  branch**

All variants used:

```text
initial_state_file = high-precision no-impact -12.70 V state
start              = -12.70 V
stop               = -13.20 V
step               = -0.005 V
impact_ionization  = none
sweep.diagnostics.newton_history.enabled = true
```

Report root:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\noimpact_max_update_branch_scan
```

Results:

| `max_update` | points | final focus `log10(electrons_m^-3)` | final electron current (A/m) | max Newton iters/point | min damping |
|---:|---:|---:|---:|---:|---:|
| `5.0` | 101 | `12.343` | `-1.568e-08` | 3 | 1 |
| `1.0` | 101 | `12.343` | `-1.568e-08` | 4 | 1 |
| `0.5` | 101 | `12.343` | `-1.568e-08` | 6 | 1 |
| `0.2` | 101 | `12.343` | `-1.568e-08` | 8 | 1 |
| `0.1` | 101 | `12.343` | `-1.568e-08` | 14 | 1 |
| `0.05` | 101 | `12.343` | `-1.568e-08` | 25 | 1 |

More aggressive caps were not useful:

```text
max_update = 0.02: fails at the -12.70 V start point with max_iterations
max_update = 0.01: fails at the -12.70 V start point with max_iterations
```

The `0.02` failure still has a small Poisson-dominated residual
(`psi ~= 1.13e-9`, `phin ~= 5.14e-15`) after 40 Newton iterations, which
indicates over-restriction rather than a recovered low-density branch.

- [x] **Step 2: Interpret branch-control result**

The existing `max_update` cap is not sufficient branch control for pn2d
Sentaurus-default BV reproduction. It only changes the iteration count and, if
made too strict, prevents convergence at the restart point. It does not reduce
the no-impact `-13.2 V` focus density overshoot.

Updated next research direction:

- Do not promote `max_update` tuning as the BV fix.
- The next opt-in solver experiment should be a true continuation-control
  change, not just a Newton step cap:
  - pseudo-transient/source-term homotopy on the electrostatic equation, or
  - explicit predictor/extrapolation control for coupled `psi/phin/phip`
    between bias points, or
  - a block-aware trust region that separately constrains electrostatic and
    continuity updates.
- Any implementation must be TDD-first and gated on the no-impact `-13.2 V`
  density target plus unchanged `-10 V` parity.

## Task 15: External Linear Predictor Proxy

Task 14 showed that plain Newton update capping is not enough. Before adding a
new continuation predictor to production code, run an external restart-based
proxy for Sentaurus-like extrapolation.

- [x] **Step 1: Build external predictor states from prior converged states**

The proxy constructs each target initial state from two previous converged
restart CSVs:

```text
x_pred(target) = x_curr + (target - bias_curr) / (bias_curr - bias_prev)
                 * (x_curr - x_prev)
```

Only `psi`, `phin`, and `phip` are extrapolated; carrier densities in the
restart CSV are kept positive placeholders because Newton recomputes densities
from the potentials.

- [x] **Step 2: Run coarse `-0.05 V` external predictor continuation**

Report root:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\noimpact_external_linear_predictor_m12p65_m13p20_step05
```

The coarse predictor used the high-precision no-impact `-12.65 V` and
`-12.70 V` states as the initial pair, then solved single target-bias decks
through `-13.20 V`.

Key result:

```text
-12.75 V: focus log10(electrons_m^-3) ~= 10.81
-12.90 V: focus log10(electrons_m^-3) ~= 12.30
-13.20 V: focus log10(electrons_m^-3) ~= 10.51
```

The predictor eventually leaves the high-density branch, but after about
`-12.95 V` the terminal current columns collapse to zero while drift/diffusion
subcolumns remain nonzero. This is not yet an acceptance-quality branch; it is
a branch-selection signal that needs classification.

- [x] **Step 3: Run fine `-0.005 V` external predictor continuation**

Report root:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\noimpact_external_linear_predictor_m12p65_m13p20_step005
```

Fine-step predictor summary:

```text
-12.75 V: focus log10(electrons_m^-3) ~= 10.73
-12.80 V: focus log10(electrons_m^-3) ~= 11.61
-12.85 V: focus log10(electrons_m^-3) ~= 12.10
-12.90 V: focus log10(electrons_m^-3) ~= 12.29
-12.95 V: focus log10(electrons_m^-3) ~= 11.96, terminal current columns = 0
-13.20 V: focus log10(electrons_m^-3) ~= 10.88, terminal current columns = 0
```

Interpretation:

- Linear predictor/extrapolation is a real branch-control lever; unlike
  `max_update`, it changes the selected high-bias branch.
- The current external predictor does not yet reproduce Sentaurus default BV:
  it improves the `-13.2 V` focus density by roughly `1.5-1.8 decades`, but it
  lands on a suspicious low-current branch with zero terminal current columns.
- Do not implement or promote predictor blindly. The next production change
  should first add an opt-in, TDD-covered continuation predictor with explicit
  diagnostics for predicted initial state quality, terminal current consistency,
  and branch acceptance gates.

## Task 16: Charon And DEVSIM Guided Continuation And Branch Diagnostics

Task 15 showed that a linear predictor changes the selected high-bias branch,
but the current external proxy can land on a suspicious branch with zero
terminal-current columns. Use the Charon and DEVSIM codebase review to turn
that observation into production-safe diagnostics before enabling any new
branch-control behavior by default.

Reference findings to preserve during this task:

```text
Charon:
  - van Overstraeten supports GradQuasiFermi, GradPotentialParallelJ,
    GradPotentialParallelJtot, EffectiveFieldParallelJ, and
    EffectiveFieldParallelJtot.
  - low-density driving-force damping uses n * F / (n + n0) and
    p * F / (p + p0).
  - current-density avalanche generation is alpha_n * |Jn| + alpha_p * |Jp|.
  - current constraints add contact voltages/current responses to a
    continuation-style extended nonlinear system.

DEVSIM:
  - rampbias restores the last converged bias and halves the step on
    convergence failure.
  - SG currents are model expressions over edge Bernoulli functions.
  - ExprEquation separates EdgeCouple, EdgeNodeVolume, NodeVolume, and
    element-volume assembly paths, which is useful as a parity checklist.

Vela current state:
  - already has van_overstraeten, quasi_fermi_gradient,
    current_density/density_gradient, ref-density interpolation, adaptive
    retry, max_update, line-search history, and Newton history diagnostics.
  - lacks an in-core continuation predictor, current-constraint continuation,
    and a branch acceptance gate that can reject the zero-terminal-current
    predictor branch.
```

- [x] **Step 1: Add a source-backed investigation note to this plan**

Append a short subsection under this task after the code changes are complete
with the exact Vela/Charon/DEVSIM files that were used to justify each solver
knob. The subsection must include these local paths:

```text
D:\code-repo\tcad-charon\src\evaluators\Charon_Avalanche_vanOverstraeten_impl.hpp
D:\code-repo\tcad-charon\src\Charon_CurrentConstraintModelEvaluator.hpp
D:\code-repo\tcad-charon\src\solver\Charon_Solver_SteadyStateConstraint.cpp
D:\code-repo\devsim\python_packages\ramp.py
D:\code-repo\devsim\python_packages\simple_dd.py
D:\code-repo\devsim\src\AutoEquation\ExprEquation.cc
include\vela\physics\ImpactIonizationModel.h
include\vela\equation\AssemblerUtils.h
src\simulation\DCSweep.cpp
src\solver\NewtonSolver.cpp
```

Expected: the plan records that Vela already contains the main Sentaurus-default
avalanche formula path, so the next change targets branch control and
diagnostics rather than new ionization coefficients.

- [x] **Step 2: Write failing parser tests for opt-in sweep predictor config**

Modify `tests/test_dc_sweep.cpp` with focused tests for a new optional config:

```json
"sweep": {
  "continuation": {
    "predictor": {
      "mode": "linear",
      "fields": ["psi", "phin", "phip"],
      "max_extrapolation_ratio": 2.0
    },
    "branch_acceptance": {
      "terminal_current_consistency": true,
      "min_terminal_current_ratio": 1.0e-6
    }
  }
}
```

Required assertions:

```text
mode accepts: none, constant, linear, secant
fields accepts only psi, phin, phip
max_extrapolation_ratio must be finite and >= 1
terminal_current_consistency defaults to false
min_terminal_current_ratio must be finite and non-negative
```

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\test_dc_sweep.exe "DCSweep: continuation predictor config"
```

Expected before implementation: the test fails because the config is not parsed.

- [x] **Step 3: Add predictor and branch-acceptance config structs**

Modify `include/vela/simulation/DCSweep.h`:

```cpp
struct SweepPredictorConfig {
    std::string mode = "none";
    std::vector<std::string> fields;
    Real maxExtrapolationRatio = 2.0;
};

struct SweepBranchAcceptanceConfig {
    bool terminalCurrentConsistency = false;
    Real minTerminalCurrentRatio = 0.0;
};

struct SweepContinuationConfig {
    SweepPredictorConfig predictor;
    SweepBranchAcceptanceConfig branchAcceptance;
};
```

Add `SweepContinuationConfig continuation;` to `DCSweepConfig`.

Expected: headers compile, but parsing still fails until Step 4.

- [x] **Step 4: Implement parser validation without changing default behavior**

Modify `src/simulation/DCSweep.cpp` to parse `sweep.continuation`.

Validation rules:

```text
predictor.mode:
  allowed = none, constant, linear, secant
predictor.fields:
  if empty, default to psi, phin, phip when mode != none
  every entry must be psi, phin, or phip
predictor.max_extrapolation_ratio:
  finite and >= 1.0
branch_acceptance.min_terminal_current_ratio:
  finite and >= 0.0
```

Expected: all new parser tests pass and existing sweep behavior is unchanged
when `sweep.continuation` is absent.

- [x] **Step 5: Write failing predictor state tests**

Add tests in `tests/test_dc_sweep.cpp` for a helper that predicts a restart
state from one or two accepted solutions:

```text
mode = none:
  returns the current accepted solution unchanged
mode = constant:
  returns the current accepted solution unchanged
mode = linear:
  x_pred = x_curr + clamp(ratio, -max_ratio, max_ratio) * (x_curr - x_prev)
mode = secant:
  same first implementation as linear, but keep the separate mode name for
  future LOCA-like predictor behavior
```

Expected before implementation: tests fail because the helper does not exist.

- [x] **Step 6: Implement predictor helper as internal DCSweep detail**

Implement the helper under `namespace detail` and keep it internal to DCSweep
continuation behavior.
Only extrapolate selected fields from `DDSolution`; leave derived carrier
densities to be recomputed by Newton from `psi`, `phin`, and `phip`.

Acceptance:

```text
predictor disabled by default
prediction is used only as the initial state for the next attempted bias
failed attempts do not update predictor history
accepted attempts update predictor history after branch acceptance passes
```

Expected: predictor unit tests pass without changing existing no-predictor
tests.

- [x] **Step 7: Add terminal-current branch acceptance diagnostics**

Extend `DCSweepPoint` and the sweep CSV with diagnostic columns:

```text
predictor_mode
predicted_initial_state
branch_acceptance_status
branch_acceptance_reason
terminal_current_consistency_ratio
```

Define the terminal-current consistency ratio as:

```text
abs(total_terminal_current) /
max(abs(electron_drift) + abs(electron_diffusion) +
    abs(hole_drift) + abs(hole_diffusion), floor)
```

Use a small fixed floor only to avoid division by zero in diagnostics. If
`terminal_current_consistency` is enabled and the ratio is below
`min_terminal_current_ratio`, reject the point with:

```text
failure_reason = branch_acceptance_failed
branch_acceptance_reason = terminal_current_inconsistent
```

Expected: the zero-terminal-current branch observed in Task 15 can be rejected
without affecting default runs.

- [x] **Step 8: Run the no-impact predictor experiment with the in-core option**

Generate a probe deck from the Task 15 setup with:

```json
"continuation": {
  "predictor": {
    "mode": "linear",
    "fields": ["psi", "phin", "phip"],
    "max_extrapolation_ratio": 2.0
  },
  "branch_acceptance": {
    "terminal_current_consistency": true,
    "min_terminal_current_ratio": 1.0e-6
  }
}
```

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\vela_example_runner.exe --config build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\noimpact_predictor_incore_m12p65_m13p20_step005\simulation.json
```

Expected: the run either rejects the suspicious low-current branch with a
clear branch-acceptance reason, or reaches `-13.2 V` with nonzero terminal
current columns and improved focus density. Do not treat a zero-terminal-current
result as progress.

- [x] **Step 9: Add SG edge-source parity report before changing avalanche physics**

Use the existing SG avalanche edge diagnostics and add missing fields only if
needed:

```text
bias_V
edge_id
node0
node1
edge_length_m
edge_couple
edge_area_proxy
electric_field_V_per_m
electron_impact_field_V_per_m
hole_impact_field_V_per_m
electron_alpha_1_per_m
hole_alpha_1_per_m
electron_flux_proxy
hole_flux_proxy
edge_source_integral
node0_source_integral
node1_source_integral
```

Compare this report against the Sentaurus manual expectations and the Charon
`alpha * |J|` integration-point form before adding more avalanche model knobs.

Expected: any future avalanche source change is backed by a specific edge or
volume-weight mismatch, not by curve fitting.

- [x] **Step 10: Defer current-constraint continuation until the predictor gate is classified**

Do not implement a current-controlled BV solver in this task. Instead, add a
short note under this task answering:

```text
Did in-core predictor recover a nonzero-current branch?
Did terminal-current branch acceptance reject the suspicious predictor branch?
Is the remaining mismatch still present in no-impact mode?
Would a Charon-style current constraint solve a confirmed branch-selection
problem, or would it mask a terminal-current assembly bug?
```

Expected: current constraint becomes a separate future plan only if the
branch-acceptance diagnostics show the solver is selecting a mathematically
valid but physically wrong voltage-controlled branch.

### Execution Note 2026-06-19: In-Core Predictor Gate

Source-backed investigation files used for this task:

```text
D:\code-repo\tcad-charon\src\evaluators\Charon_Avalanche_vanOverstraeten_impl.hpp
D:\code-repo\tcad-charon\src\Charon_CurrentConstraintModelEvaluator.hpp
D:\code-repo\tcad-charon\src\solver\Charon_Solver_SteadyStateConstraint.cpp
D:\code-repo\devsim\python_packages\ramp.py
D:\code-repo\devsim\python_packages\simple_dd.py
D:\code-repo\devsim\src\AutoEquation\ExprEquation.cc
include\vela\physics\ImpactIonizationModel.h
include\vela\equation\AssemblerUtils.h
src\simulation\DCSweep.cpp
src\solver\NewtonSolver.cpp
```

- Added opt-in `sweep.continuation` parsing with predictor modes
  `none|constant|linear|secant`, field validation for `psi|phin|phip`, bounded
  extrapolation ratio validation, and branch-acceptance threshold validation.
- Added `SweepPredictorConfig`, `SweepBranchAcceptanceConfig`, and
  `SweepContinuationConfig` to `DCSweepConfig`.
- Added the predictor helper in `include/vela/simulation/DCSweepPredictor.h`
  under `vela::detail` so it can be TDD-covered directly while remaining an
  internal DCSweep continuation helper.
- Wired the predictor into `DCSweep` as an initial-state generator only.
  Predictor history is advanced only after an accepted point passes branch
  acceptance.
- Added continuation CSV diagnostics:
  `predictor_mode`, `predicted_initial_state`,
  `branch_acceptance_status`, `branch_acceptance_reason`, and
  `terminal_current_consistency_ratio`.
- Added terminal-current branch acceptance. When enabled, a point whose
  consistency ratio falls below `min_terminal_current_ratio` is rejected with
  `failure_reason=branch_acceptance_failed` and
  `branch_acceptance_reason=terminal_current_inconsistent`.
- Focused verification passed:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build-release --parallel --target test_dc_sweep
build-release\test_dc_sweep.exe "DCSweep: continuation predictor config is validated"
build-release\test_dc_sweep.exe "DCSweep predictor: extrapolates selected coupled variables"
build-release\test_dc_sweep.exe "DCSweep: continuation predictor writes branch diagnostics"
```

Status: Step 8's no-impact in-core predictor probe is classified in the
following execution note before any Charon-style current constraint is split
into a follow-up plan.

### Execution Note 2026-06-19: Predictor Gate Classification And SG Source Parity

Step 8 no-impact in-core predictor probe:

- Generated `build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\noimpact_predictor_incore_m12p65_m13p20_step005\simulation.json` from the Task 15 no-impact predictor setup.
- Enabled `sweep.continuation.predictor.mode = "linear"`, fields `psi|phin|phip`, `max_extrapolation_ratio = 2.0`, and terminal-current branch acceptance with `min_terminal_current_ratio = 1.0e-6`.
- The run did not recover a nonzero-current branch through `-13.2 V`. It accepted 73 points, then rejected the next point at `-12.926207151135785 V`.
- Last accepted point: bias `-12.926207151032472 V`, `branch_acceptance_status=accepted`, `terminal_current_consistency_ratio=1.0000005340965768e-06`, and electron terminal current `-1.5589793636067414e-08 A`.
- Rejected point: bias `-12.926207151135785 V`, `branch_acceptance_status=rejected`, `branch_acceptance_reason=terminal_current_inconsistent`, `terminal_current_consistency_ratio=9.9999986099197646e-07`, and zero terminal-current columns.
- Conclusion: the in-core predictor confirms the suspicious low-current branch can be detected and rejected, but it does not solve the high-bias no-impact branch mismatch.

Step 9 SG edge-source parity report:

- Added TDD coverage to require SG avalanche edge diagnostics to include `electric_field_V_per_m`, `electron_impact_field_V_per_m`, and `hole_impact_field_V_per_m`.
- Extended the C++ SG edge diagnostics CSV writer with those fields, using the assembled `SgEdgeCurrentAvalancheSourceRecord` values.
- Regenerated a C++ single-point `-20 V` SG edge-source dump under `build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\sg_cpp_minus20_with_fields`.
- Compared it against the Python SG reconstruction from `sg_python_minus20_full\sg_avalanche_edges.csv` with `scripts\compare_pn2d_bv_sg_edge_source_dump.py`.
- C++ and Python both report `3830` SG edge rows at `-20 V`, but the total source integral differs by `log10(C++/Python)=0.47369639705083294`.
- The leading structural mismatch is source placement, not just source magnitude: C++ reports `contact_edge_source_fraction=0.6578333131556747`, while Python reports essentially zero contact-edge contribution and `interior_bulk_source_fraction=0.9665216700765034`.
- Hotspot node source integrals at nodes `351` and `986` are close between the two reports, so the next source-parity investigation should focus on contact-edge classification/exclusion and edge-volume weighting rather than ionization coefficients.

Step 10 current-constraint decision:

- Did in-core predictor recover a nonzero-current branch? No. It approached the branch-ratio floor and was rejected before `-13.2 V`.
- Did terminal-current branch acceptance reject the suspicious predictor branch? Yes, with `failure_reason=branch_acceptance_failed` and `branch_acceptance_reason=terminal_current_inconsistent`.
- Is the remaining mismatch still present in no-impact mode? Yes. The no-impact predictor probe still fails to reach the Sentaurus-like high-bias branch, so avalanche feedback is not the only branch-selection issue.
- Would a Charon-style current constraint solve a confirmed branch-selection problem, or mask a terminal-current assembly bug? It should be deferred. The current evidence first requires classifying voltage-controlled branch selection and terminal-current/edge-source assembly. A current-constraint continuation path may be useful only after the zero-current branch rejection and SG source-placement mismatch are resolved or isolated.

## Task 17: Add Sentaurus 2018 VM Oracle Runs Through SSH

The new SSH access changes the BV workflow from "analyze captured Sentaurus
artifacts" to "regenerate and perturb Sentaurus reference artifacts on demand".
Keep this capability opt-in because it depends on a licensed VM, but make the
local planning, manifests, and artifact staging testable without the VM.

- [x] **Step 1: Verify the Sentaurus VM contract from Windows**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
ssh sentaurus "hostname; pwd; whoami; command -v sde; command -v sdevice"
```

Expected:

```text
sentaurus
/home/tcad
tcad
```

The fourth and fifth output lines must be non-empty executable paths ending in
`sde` and `sdevice`.

If the Sentaurus commands are missing, update the VM shell startup file first
so non-interactive SSH can resolve `sde` and `sdevice`; do not encode a private
installation path in Vela source.

- [x] **Step 2: Write dry-run tests for the VM runner**

Create `tests/regression/test_sentaurus_vm_reference_runner.py` with tests that
exercise the runner without contacting the VM:

```python
#!/usr/bin/env python3

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
RUNNER = REPO / "scripts" / "run_sentaurus_vm_reference.py"


class SentaurusVmReferenceRunnerTest(unittest.TestCase):
    def test_dry_run_writes_manifest_without_ssh(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_vm_dry_") as tmp:
            root = Path(tmp)
            source = root / "source"
            source.mkdir()
            for name in [
                "pn2d_sde.cmd",
                "pn2d_bv_sdevice.cmd",
                "models.par",
            ]:
                (source / name).write_text(f"{name}\n")
            out = root / "runs"

            subprocess.run([
                sys.executable,
                str(RUNNER),
                "pn2d",
                "--ssh-target", "sentaurus",
                "--source-dir", str(source),
                "--local-output-dir", str(out),
                "--remote-root", "~/sentaurus_runs/vela_oracle",
                "--run-id", "pn2d_bv_vm_dry_run",
                "--stages", "bv",
                "--dry-run",
            ], check=True)

            manifest = json.loads(
                (out / "pn2d_bv_vm_dry_run" / "sentaurus_vm_run_manifest.json").read_text()
            )
            self.assertEqual(manifest["ssh_target"], "sentaurus")
            self.assertEqual(manifest["remote_source_dir"], "~/sentaurus_runs/vela_oracle/pn2d_bv_vm_dry_run/source")
            self.assertEqual(manifest["stages"], ["bv"])
            self.assertEqual(manifest["commands"], [
                "cd ~/sentaurus_runs/vela_oracle/pn2d_bv_vm_dry_run/source && sde -e -l pn2d_sde.cmd",
                "cd ~/sentaurus_runs/vela_oracle/pn2d_bv_vm_dry_run/source && sdevice pn2d_bv_sdevice.cmd > run_pn2d_bv.out 2>&1",
            ])

    def test_missing_required_deck_fails_before_ssh(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_vm_missing_") as tmp:
            root = Path(tmp)
            source = root / "source"
            source.mkdir()
            (source / "pn2d_sde.cmd").write_text("mesh\n")
            completed = subprocess.run([
                sys.executable,
                str(RUNNER),
                "pn2d",
                "--source-dir", str(source),
                "--local-output-dir", str(root / "runs"),
                "--run-id", "missing_bv",
                "--stages", "bv",
                "--dry-run",
            ], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("missing required source file", completed.stderr)


if __name__ == "__main__":
    unittest.main()
```

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python -m unittest tests.regression.test_sentaurus_vm_reference_runner
```

Expected before implementation: tests fail because
`scripts/run_sentaurus_vm_reference.py` does not exist.

- [x] **Step 3: Implement the dry-run capable SSH/SCP runner**

Create `scripts/run_sentaurus_vm_reference.py` with these command-line
contracts:

```text
python scripts\run_sentaurus_vm_reference.py pn2d \
  --ssh-target sentaurus \
  --source-dir reference_tcad\pn2d_sentaurus2018\source \
  --local-output-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs \
  --remote-root ~/sentaurus_runs/vela_oracle \
  --run-id pn2d_bv_vm_smoke \
  --stages bv \
  --dry-run
```

Implementation requirements:

```text
accepted device: pn2d
accepted stages: 0v, iv, bv
default ssh target: sentaurus
default remote root: ~/sentaurus_runs/vela_oracle
required common files: pn2d_sde.cmd, models.par
required stage files:
  0v -> pn2d_0v_sdevice.cmd
  iv -> pn2d_iv_sdevice.cmd
  bv -> pn2d_bv_sdevice.cmd
default example remote source dir: ~/sentaurus_runs/vela_oracle/pn2d_bv_vm_smoke/source
default example local source copy dir: build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_vm_runs/pn2d_bv_vm_smoke/source
default example manifest path: build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_vm_runs/pn2d_bv_vm_smoke/sentaurus_vm_run_manifest.json
```

For a different `--run-id`, replace only `pn2d_bv_vm_smoke` in the three
example paths with the exact run id passed on the command line.

For non-dry-run execution, use `subprocess.run(..., check=True)` with command
arrays rather than shell string concatenation:

```python
[
    ["ssh", "sentaurus", "mkdir -p ~/sentaurus_runs/vela_oracle/pn2d_bv_vm_smoke/source"],
    ["scp", "reference_tcad/pn2d_sentaurus2018/source/pn2d_sde.cmd", "sentaurus:~/sentaurus_runs/vela_oracle/pn2d_bv_vm_smoke/source/"],
    ["scp", "reference_tcad/pn2d_sentaurus2018/source/models.par", "sentaurus:~/sentaurus_runs/vela_oracle/pn2d_bv_vm_smoke/source/"],
    ["scp", "reference_tcad/pn2d_sentaurus2018/source/pn2d_bv_sdevice.cmd", "sentaurus:~/sentaurus_runs/vela_oracle/pn2d_bv_vm_smoke/source/"],
    ["ssh", "sentaurus", "cd ~/sentaurus_runs/vela_oracle/pn2d_bv_vm_smoke/source && sde -e -l pn2d_sde.cmd"],
    ["ssh", "sentaurus", "cd ~/sentaurus_runs/vela_oracle/pn2d_bv_vm_smoke/source && sdevice pn2d_bv_sdevice.cmd > run_pn2d_bv.out 2>&1"],
    ["scp", "sentaurus:~/sentaurus_runs/vela_oracle/pn2d_bv_vm_smoke/source/*.tdr", "build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_vm_runs/pn2d_bv_vm_smoke/source/"],
    ["scp", "sentaurus:~/sentaurus_runs/vela_oracle/pn2d_bv_vm_smoke/source/*.plt", "build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_vm_runs/pn2d_bv_vm_smoke/source/"],
    ["scp", "sentaurus:~/sentaurus_runs/vela_oracle/pn2d_bv_vm_smoke/source/*.log", "build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_vm_runs/pn2d_bv_vm_smoke/source/"],
    ["scp", "sentaurus:~/sentaurus_runs/vela_oracle/pn2d_bv_vm_smoke/source/*.grd", "build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_vm_runs/pn2d_bv_vm_smoke/source/"],
    ["scp", "sentaurus:~/sentaurus_runs/vela_oracle/pn2d_bv_vm_smoke/source/*.dat", "build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_vm_runs/pn2d_bv_vm_smoke/source/"],
    ["scp", "sentaurus:~/sentaurus_runs/vela_oracle/pn2d_bv_vm_smoke/source/run_pn2d_*.out", "build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_vm_runs/pn2d_bv_vm_smoke/source/"],
]
```

Handle missing glob classes as warnings after the run, not as run failure,
because some Sentaurus stages may not emit every suffix.

- [x] **Step 4: Document the VM runner contract**

Modify `docs/sentaurus_vm_ssh_workflow.md` by adding a section named
`## Vela Runner Contract` with this content:

```markdown
## Vela Runner Contract

The automated Vela runner uses:

- SSH target: `sentaurus`
- Remote root: `~/sentaurus_runs/vela_oracle`
- Local pn2d source: `reference_tcad/pn2d_sentaurus2018/source`
- Example staged output: `build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_vm_runs/pn2d_bv_vm_smoke/source`
- Required common files: `pn2d_sde.cmd`, `models.par`
- BV command: `sdevice pn2d_bv_sdevice.cmd > run_pn2d_bv.out 2>&1`

The runner must never overwrite `reference_tcad/pn2d_sentaurus2018/source`
directly. Review staged artifacts first, then copy selected files deliberately.
```

Expected: future workers can run the same workflow from the plan without
reverse-engineering the VM notes.

- [x] **Step 5: Run dry-run verification**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python -m unittest tests.regression.test_sentaurus_vm_reference_runner
python scripts\run_sentaurus_vm_reference.py pn2d --ssh-target sentaurus --source-dir reference_tcad\pn2d_sentaurus2018\source --local-output-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs --remote-root ~/sentaurus_runs/vela_oracle --run-id pn2d_bv_vm_dry_run --stages bv --dry-run
```

Expected:

```text
unittest passes
build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_vm_dry_run\sentaurus_vm_run_manifest.json exists
```

- [x] **Step 6: Run a live BV oracle smoke test**

Run only after Step 1 passes:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\run_sentaurus_vm_reference.py pn2d --ssh-target sentaurus --source-dir reference_tcad\pn2d_sentaurus2018\source --local-output-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs --remote-root ~/sentaurus_runs/vela_oracle --run-id pn2d_bv_vm_smoke --stages bv
```

Expected staged files:

```text
build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_vm_smoke\source\pn2d_bv.plt
build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_vm_smoke\source\pn2d_bv_des.tdr
build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_vm_smoke\source\pn2d_bv_multibias_0200_des.tdr
build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_vm_smoke\source\run_pn2d_bv.out
```

Expected log check:

```powershell
Get-Content build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_vm_smoke\source\run_pn2d_bv.out -Tail 80
```

The tail must show normal SDevice completion and no license, mesh, or physics
fatal error.

- [x] **Step 7: Import the VM-generated BV snapshots without overwriting local fixtures**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
build-release\sentaurus_import.exe --tdr build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_vm_smoke\source\pn2d_bv_multibias_0132_des.tdr --export-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_vm_smoke\exports\sentaurus_-13.2v
build-release\sentaurus_import.exe --tdr build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_vm_smoke\source\pn2d_bv_multibias_0200_des.tdr --export-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_vm_smoke\exports\sentaurus_-20v
```

Expected:

```text
exports\sentaurus_-13.2v\fields\ElectrostaticPotential_region0.csv
exports\sentaurus_-13.2v\fields\ElectricField_region0.csv
exports\sentaurus_-13.2v\fields\eDensity_region0.csv
exports\sentaurus_-13.2v\fields\hDensity_region0.csv
exports\sentaurus_-13.2v\fields\ImpactIonization_region0.csv
exports\sentaurus_-20v\fields\ElectrostaticPotential_region0.csv
exports\sentaurus_-20v\fields\ElectricField_region0.csv
exports\sentaurus_-20v\fields\eDensity_region0.csv
exports\sentaurus_-20v\fields\hDensity_region0.csv
exports\sentaurus_-20v\fields\ImpactIonization_region0.csv
```

- [x] **Step 8: Compare VM Sentaurus output against the captured local source**

Add or reuse a small comparison command that checks the VM-run `.plt` curve and
imported field snapshots against the current captured local artifacts. The
acceptance gate is:

```text
BV curve: max abs log10 current difference <= 1.0e-6 at matching bias rows
-13.2 V fields: electrostatic potential max abs diff <= 1.0e-9 V after import
-20 V fields: electrostatic potential max abs diff <= 1.0e-9 V after import
field row counts and node ids must match exactly
```

If this gate fails, do not use the VM run as the new oracle until the cause is
classified as version, command-deck, license-feature, randomization, or import
difference.

- [x] **Step 9: Generate Sentaurus one-factor BV variants for the known Vela gaps**

Use the VM runner to make staged copies of the BV deck and run one change at a
time:

```text
pn2d_bv_vm_default:
  original pn2d_bv_sdevice.cmd

pn2d_bv_vm_no_avalanche:
  remove Avalanche(VanOverstraeten) from Recombination

pn2d_bv_vm_avaldensgradqf:
  add Math { AvalDensGradQF }

pn2d_bv_vm_element_volume:
  add Math { ElementVolumeAvalanche }

pn2d_bv_vm_refdens_efield:
  add RefDens_eGradQuasiFermi_ElectricField=1e8
  add RefDens_hGradQuasiFermi_ElectricField=1e8
```

Each variant must have its own staged run directory under
`build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\`.

Expected: each variant produces a BV `.plt`, final `.tdr`, and multibias TDRs
or a clear SDevice failure log. Do not mix variants into the committed source
fixture.

- [x] **Step 10: Use VM variants to reprioritize the Vela BV roadmap**

After importing the staged variants, add a new execution note under this task
that references a generated CSV report:

```text
build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\variant_summary.csv
```

The CSV must contain exactly these columns:

```text
variant,bias_V,current_total_A,electron_density_log10_p50_vs_default,hole_density_log10_p50_vs_default,impact_ionization_log10_p99_vs_default,contact_edge_source_fraction,interior_bulk_source_fraction,run_status,run_log_tail
```

Required rows:

```text
pn2d_bv_vm_default at -13.2 V
pn2d_bv_vm_default at -20 V
pn2d_bv_vm_no_avalanche at -13.2 V
pn2d_bv_vm_no_avalanche at -20 V
pn2d_bv_vm_avaldensgradqf at -13.2 V
pn2d_bv_vm_avaldensgradqf at -20 V
pn2d_bv_vm_element_volume at -13.2 V
pn2d_bv_vm_element_volume at -20 V
pn2d_bv_vm_refdens_efield at -13.2 V
pn2d_bv_vm_refdens_efield at -20 V
```

Interpretation rules:

```text
If no-avalanche Sentaurus matches Vela no-impact high-density behavior:
  prioritize pre-avalanche continuity, mobility, SRH, and contact minority-carrier boundary parity.

If no-avalanche Sentaurus remains on the low-density branch:
  prioritize Vela no-impact branch selection before changing avalanche physics.

If AvalDensGradQF moves Sentaurus toward Vela mobility_density_gradient:
  keep mobility_density_gradient as a control path and document the expected breakdown shift.

If ElementVolumeAvalanche or RefDens changes source placement toward Vela:
  add a narrowly scoped Vela source-volume or driving-force interpolation task.

If none of the Sentaurus one-factor variants approach Vela:
  prioritize Vela terminal-current assembly, contact-edge exclusion, and edge-volume weighting.
```

Expected: the next implementation plan is based on live Sentaurus one-factor
evidence, not inferred from manual text alone.

### Execution Note 2026-06-19: Sentaurus VM Oracle And One-Factor Variants

VM contract:

- `ssh sentaurus "hostname; pwd; whoami; command -v sde; command -v sdevice"` passed through Windows OpenSSH.
- The VM reports `hostname=sentaurus`, `pwd=/home/tcad`, `whoami=tcad`, `sde=/usr/synopsys/sentaurus/O_2018.06-SP2/bin/sde`, and `sdevice=/usr/synopsys/sentaurus/O_2018.06-SP2/bin/sdevice`.

Runner and dry-run gate:

- Added `scripts/run_sentaurus_vm_reference.py`.
- Added `tests/regression/test_sentaurus_vm_reference_runner.py`.
- Added `## Vela Runner Contract` to `docs/sentaurus_vm_ssh_workflow.md`.
- `python -m unittest tests.regression.test_sentaurus_vm_reference_runner` passed.
- Dry-run manifest was written under `build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_vm_dry_run\sentaurus_vm_run_manifest.json`.

Live default BV oracle:

- Ran `pn2d_bv_vm_smoke` on the Sentaurus 2018 VM.
- The run completed normally. `run_pn2d_bv.out` ends with `Sentaurus Device simulation finished` and `Good Bye`.
- Fetched `pn2d_bv.plt`, `pn2d_bv_des.tdr`, `pn2d_bv_multibias_0000_des.tdr` through `pn2d_bv_multibias_0200_des.tdr`, and logs into `build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_vm_smoke\source`.
- Imported `pn2d_bv_multibias_0132_des.tdr` and `pn2d_bv_multibias_0200_des.tdr` under the smoke run's `exports` directory.
- Compared VM output against the captured local source in `vm_vs_captured_compare.json`: BV curve had `223` matching bias rows with `max_abs_log10_current_delta=0`, and both `-13.2 V` and `-20 V` electrostatic-potential fields had matching row counts, matching node ids, and `max_abs_diff=0 V`.

One-factor variants:

- Generated staged source directories under `build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_variant_sources`.
- Ran and fetched:
  - `pn2d_bv_vm_default`
  - `pn2d_bv_vm_no_avalanche`
  - `pn2d_bv_vm_avaldensgradqf`
  - `pn2d_bv_vm_element_volume`
  - `pn2d_bv_vm_refdens_efield`
- The first `pn2d_bv_vm_refdens_efield` run failed because Sentaurus 2018.06-SP2 rejects `RefDens_eGradQuasiFermi_ElectricField_Aval` and `RefDens_hGradQuasiFermi_ElectricField_Aval`. The corrected 2018-compatible GradQuasiFermi keywords are `RefDens_eGradQuasiFermi_ElectricField` and `RefDens_hGradQuasiFermi_ElectricField`; rerunning with those names passed.
- Imported `-13.2 V` and `-20 V` snapshots for all five variants.
- Wrote `build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\variant_summary.csv`.

Variant summary:

```text
variant                         bias    current_total_A       eDensity p50 log10 vs default   hDensity p50 log10 vs default   ImpactIonization p99 log10 vs default
pn2d_bv_vm_default              -13.2   -8.38472088807e-17    0                              0                              0
pn2d_bv_vm_default              -20     -9.10455666344e-16    0                              0                              0
pn2d_bv_vm_no_avalanche         -13.2   -4.87695104455e-17   -0.0264649426634              -0.0926398767196               0
pn2d_bv_vm_no_avalanche         -20     -6.16408127806e-17   -0.553372546463               -0.997417601305                0
pn2d_bv_vm_avaldensgradqf       -13.2   -8.50410488397e-17    0.00126082318652              0.000392508733719              0.165593033388
pn2d_bv_vm_avaldensgradqf       -20     -1.3424381206e-15     0.130430732558                0.163966360678                 0.18329092653
pn2d_bv_vm_element_volume       -13.2   -8.38472088807e-17    0                              0                              0
pn2d_bv_vm_element_volume       -20     -9.10455666344e-16    0                              0                              0
pn2d_bv_vm_refdens_efield       -13.2   -8.47663242962e-17    0.00498899901103              0.00568829570008               2.40533753541
pn2d_bv_vm_refdens_efield       -20     -1.04928016535e-15    0.0514114290454               0.0632378581813                1.38272628696
```

Interpretation:

- No-avalanche Sentaurus remains on the low-density branch and actually reduces high-bias carrier medians versus the default avalanche run. This confirms the Vela no-impact high-density branch is a Vela-side branch-selection/parity issue, not a Sentaurus default avalanche feedback feature.
- `AvalDensGradQF` moves Sentaurus only modestly at `-20 V` (`~0.13-0.16` density decades and current from `-9.10e-16 A` to `-1.34e-15 A`). Keep Vela `mobility_density_gradient` as a control path, not as the Sentaurus-default acceptance path.
- `ElementVolumeAvalanche` is numerically identical to the default for this PN2D deck at the reported metrics, so it is not the missing Vela parity switch.
- `RefDens_*GradQuasiFermi_ElectricField` strongly changes high-percentile impact-generation values but moves density/current only modestly; it does not reproduce Vela's multi-decade high-density branch.
- The next Vela work should prioritize the no-impact branch-selection mismatch and terminal/edge assembly localization: contact-edge exclusion, edge-volume weighting, and terminal-current consistency remain higher priority than adding more Sentaurus avalanche model knobs.

## Task 18: Compare VM No-Avalanche Sentaurus Against Vela No-Impact Branch

The Task 17 VM variants prove that Sentaurus without avalanche stays on the
low-density branch. Use that live VM oracle to compare the Vela no-impact
branch at the same `-13.2 V` bias, where Vela already shows the high-density
branch jump.

- [x] **Step 1: Reuse existing trusted artifacts**

Inputs:

```text
Sentaurus no-avalanche root:
build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_vm_no_avalanche\exports

Vela no-impact VTK root:
build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\bv_noimpact_resume_minus8_minus10_minus13p2\vtk

Vela avalanche-enabled reference VTK root:
build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\bv_sg_default_from_noimpact_minus13p2\vtk
```

The reliable Vela no-impact branch artifact reaches `-13.2 V`; no equally
trusted no-impact `-20 V` VTK exists because later predictor probes were rejected
around `-12.93 V` by the terminal-current consistency gate.

- [x] **Step 2: Run branch-profile comparison**

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\diagnose_pn2d_bv_branch_profiles.py `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\vela\mesh.json `
  --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_vm_no_avalanche\exports `
  --vela-avalanche-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\bv_sg_default_from_noimpact_minus13p2\vtk `
  --vela-no-impact-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\bv_noimpact_resume_minus8_minus10_minus13p2\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\vm_noavalanche_vs_vela_noimpact_branch_profiles `
  --biases -13.2
```

Output:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\vm_noavalanche_vs_vela_noimpact_branch_profiles
```

Key `Vela noimpact - Sentaurus no-avalanche` medians:

```text
band            d(psi) V   d(psi-phin) V   d(phip-psi) V   log10(n ratio)   log10(p ratio)
left_p          0.260851   -0.008371       0.013049        -0.359592        0.000109
pre_junction_p  0.194724    0.196708      -0.005648         3.084902       -0.313327
junction        0.137000    0.181523      -0.005449         2.819986       -0.310957
post_junction_n 0.079251    0.166323      -0.006608         2.575127       -0.329952
right_n         0.013034    0.013034      -0.013488         0.000000       -0.445544
```

Focus nodes:

```text
node  branch      psi-phin V   phip-psi V   log10(n cm^-3)   log10(p cm^-3)
955   noimpact    -0.231100   -0.386500    6.117656         3.506820
955   sentaurus   -0.419675   -0.383167    3.168737         3.782044
1089  noimpact    -0.214600   -0.424950    6.394935         2.861108
1089  sentaurus   -0.389688   -0.413625    3.672491         3.270373
351   noimpact    -0.217710   -0.401390    6.342555         3.256963
351   sentaurus   -0.405879   -0.399687    3.474296         3.578307
986   noimpact    -0.217720   -0.402250    6.342438         3.242527
986   sentaurus   -0.400753   -0.396340    3.486607         3.560747
```

- [x] **Step 3: Run continuity-feedback comparison near edge 2886**

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\diagnose_pn2d_bv_continuity_feedback.py `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\vela\mesh.json `
  --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\vela\doping.csv `
  --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_vm_no_avalanche\exports `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\bv_noimpact_resume_minus8_minus10_minus13p2\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\vm_noavalanche_vs_vela_noimpact_continuity_feedback `
  --biases -13.2 `
  --edge-id 2886
```

Output:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\vm_noavalanche_vs_vela_noimpact_continuity_feedback
```

Near nodes `351/986` and incident neighbors, `delta(psi-phin)` is consistently
`0.182-0.188 V`, while `delta(phip-psi)` is only about `-0.0015` to `-0.0060 V`.
The local electron-density ratio is therefore `~2.84-2.89` decades, with hole
density slightly lower in Vela by `~0.32` decades. Sentaurus no-avalanche has
zero impact-generation node integral in this focus region, as expected.

Top focus/incident edges show that local potential and electron quasi-Fermi
gradients are close:

```text
edge  relation           dpsi Vela/Sentaurus       dphin Vela/Sentaurus      log10(n ratio)
1085  incident_to_focus  0.355140 / 0.358528       0.355170 / 0.361893       2.877086
2886  focus              0.355160 / 0.358528       0.355170 / 0.353402       2.862002
2893  incident_to_focus  0.345750 / 0.349088       0.345760 / 0.348287       2.849000
```

Interpretation:

- The VM-backed no-avalanche comparison confirms that the Vela no-impact
  mismatch is not an avalanche model effect and not primarily an edge-gradient
  mismatch at the focus edge.
- The high electron density follows an absolute `psi-phin` offset of about
  `0.18 V` in the pre-junction and junction bands. Because the local `dpsi` and
  `dphin` gradients are close on edge `2886`, the next root-cause work should
  inspect electron quasi-Fermi absolute reference/contact anchoring and boundary
  equations, not just SG edge flux formulas.
- Hole quasi-Fermi alignment is much closer than electron quasi-Fermi alignment,
  so the next diagnostic should be electron-continuity specific: contact
  electron boundary values, Dirichlet/minority-carrier treatment, quasi-Fermi
  gauge/reference shifts, and the residual/Jacobian rows that fix `phin` under
  no-impact reverse bias.

## Task 19: Diagnose Electron Quasi-Fermi Anchor And Gauge Shift Hypotheses

Task 18 narrowed the no-impact mismatch to the absolute `psi-phin` exponent,
while showing that local edge gradients near edge `2886` are close. The next
diagnostic should separate three hypotheses:

1. Vela contact `phin` is mis-anchored.
2. A uniform electron quasi-Fermi gauge shift could reconcile the interior
   density branch.
3. The mismatch is an interior electron-continuity/electrostatic branch issue
   despite correct contact quasi-Fermi values.

- [x] **Step 1: Add a TDD regression for the qF-anchor diagnostic**

Added `ReferenceTcadToolsTest.test_pn2d_bv_qf_anchor_script_reports_internal_phin_shift`.

The synthetic fixture has three nodes:

```text
node 0: Anode contact
node 1: interior
node 2: Cathode contact
```

The Vela synthetic state sets only the interior `psi-phin` exponent `0.18 V`
above Sentaurus while keeping contact `phin` equal to Sentaurus. The expected
diagnostic result is:

```text
delta_psi_minus_phin_median_V = 0.18
uniform_phin_shift_to_match_electron_median_V = 0.18
contact_phin_violation_if_uniform_shift_V > 0.17
contact delta_phin_median_V ~= 0
```

RED verification:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_qf_anchor_script_reports_internal_phin_shift
```

The test failed because `scripts\diagnose_pn2d_bv_qf_anchor.py` did not exist.

- [x] **Step 2: Implement the qF-anchor diagnostic script**

Added `scripts/diagnose_pn2d_bv_qf_anchor.py`.

The script writes:

```text
qf_anchor_contact_summary.csv
qf_anchor_band_summary.csv
qf_anchor_focus_nodes.csv
qf_anchor_summary.json
```

For each band it reports:

```text
delta_psi_median_V
delta_phin_median_V
delta_phip_median_V
delta_psi_minus_phin_median_V
electron_log10_ratio_median
electron_log10_ratio_from_exponent_median
uniform_phin_shift_to_match_electron_median_V
contact_phin_violation_if_uniform_shift_V
```

GREEN verification:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_qf_anchor_script_reports_internal_phin_shift
```

Result: passed.

- [x] **Step 3: Run the diagnostic on VM no-avalanche versus Vela no-impact**

Command:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\diagnose_pn2d_bv_qf_anchor.py `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\vela\mesh.json `
  --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_vm_no_avalanche\exports `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\bv_noimpact_resume_minus8_minus10_minus13p2\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\vm_noavalanche_vs_vela_noimpact_qf_anchor `
  --biases -13.2
```

Output:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\vm_noavalanche_vs_vela_noimpact_qf_anchor
```

Contact summary:

```text
contact   dpsi V      dphin V       dphip V       d(psi-phin) V   log10(n ratio)   log10(p ratio)
Cathode   0.013034    2.30e-16      0             0.013034        ~0               -0.437921
Anode    -0.013049    8.88e-15      1.07e-14     -0.013049       -0.437921        ~0
```

The electron and hole quasi-Fermi contact values are pinned to Sentaurus within
roundoff. This rejects a simple contact `phin/phip` Dirichlet-value bug.

Band summary:

```text
band            dpsi V    dphin V    d(psi-phin) V   log10(n ratio)   log10(n from exponent)   contact violation if shifted
left_p          0.260851  0.267867  -0.008371       -0.359592        -0.140622                0.008371
pre_junction_p  0.194724 -0.008475   0.196708        3.084902         3.304553                0.196708
junction        0.137000 -0.045352   0.181523        2.819986         3.049447                0.181523
post_junction_n 0.079251 -0.077144   0.166323        2.575127         2.794111                0.166323
right_n         0.013034 ~0          0.013034        ~0               0.218962                0.013034
```

Focus nodes:

```text
node  d(psi-phin) V   log10(n ratio)   log10(n from exponent)   d(phip-psi) V
955   0.188575        2.948919         3.167913                 -0.003333
1089  0.175088        2.722444         2.941347                 -0.011325
351   0.188169        2.868259         3.161096                 -0.001703
986   0.183033        2.855832         3.074817                 -0.005910
```

Interpretation:

- A uniform `phin` shift of `~0.17-0.20 V` would be required to bring the
  pre-junction and junction electron exponent toward Sentaurus, but that shift
  would violate the correctly pinned contact `phin` rows by the same amount.
- The no-impact branch mismatch is therefore not a legal global quasi-Fermi
  gauge issue. It is a spatially varying interior branch/electrostatic-continuity
  issue: `psi` remains elevated while `phin` is pulled downward from the left-p
  plateau through the junction, producing the electron density branch.
- Next implementation/debug step: inspect the electron-continuity Jacobian and
  residual rows that connect the correctly pinned contact `phin` values to the
  pre-junction/junction interior, especially whether the no-impact branch is
  selected by an electron-continuity transport/recombination balance or by the
  Poisson-coupled interior `psi` shape. Do not add a global `phin` gauge knob.

## Task 20: Re-Evaluate No-Avalanche States With Vela Residual Probe

Task 19 showed that contact `phin` is correctly pinned and that a global `phin`
shift is illegal. Reuse the existing Newton residual-state diagnostic to test
which external state combinations Vela's no-impact residual evaluator accepts
or rejects.

- [x] **Step 1: Run residual-state comparison with the VM no-avalanche oracle**

Command:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\diagnose_pn2d_bv_newton_residual_states.py `
  --base-config build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\bv_noimpact_resume_minus8_minus10_minus13p2.json `
  --runner build-release\vela_example_runner.exe `
  --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_vm_no_avalanche\exports `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\bv_noimpact_resume_minus8_minus10_minus13p2\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\vm_noavalanche_vs_vela_noimpact_residual_states `
  --states vela:-13.2,sentaurus:-13.2,hybrid_vpsi_sqf:-13.2,hybrid_spsi_vqf:-13.2,hybrid_spsi_shift_vqf:-13.2 `
  --nodes 955,1089,351,986,349,352,987,988
```

Output:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\vm_noavalanche_vs_vela_noimpact_residual_states
```

- [x] **Step 2: Compare global residual blocks**

```text
state                         psi block        phin block       phip block       combined
vela_m13p2v                   1.799017e-1      1.634672e-9      1.214231e-13    1.799017e-1
sentaurus_m13p2v              5.990493e+1      5.363926e-10     1.647774e-12    5.990493e+1
hybrid_vpsi_sqf_m13p2v        1.200514e+2      8.230532e-10     1.640791e-12    1.200514e+2
hybrid_spsi_vqf_m13p2v        2.584038e+6      1.634448e-9      2.870442e-10    2.584038e+6
hybrid_spsi_shift_vqf_m13p2v  1.837405e+4      4.558133e-9      2.296583e-12    1.837405e+4
```

The Sentaurus no-avalanche state has a smaller Vela `phin` block than the Vela
no-impact state, but a much larger Poisson block. This means the low-density
Sentaurus electron quasi-Fermi field is not rejected by Vela's electron
continuity equation at `-13.2 V`; the dominant incompatibility remains the
electrostatic potential/charge balance.

- [x] **Step 3: Inspect focus-node residual rows**

At focus nodes `351/986/955/1089`:

```text
state                  node   psi residual       phin residual scale
vela_m13p2v            351    5.80226e-4         2.15019e-15
vela_m13p2v            986   -3.59691e-4        -1.24396e-15
sentaurus_m13p2v       351    6.45629e-14       -1.00421e-15
sentaurus_m13p2v       986    5.99874e-6         1.15578e-15
hybrid_vpsi_sqf_m13p2v 351    5.80226e-4        -4.49171e-13
hybrid_vpsi_sqf_m13p2v 986   -3.59691e-4         1.00816e-13
```

Sentaurus `phin/phip` on Vela `psi` increases local `phin` residuals but keeps
them far below the Poisson mismatch scale. Sentaurus `psi` with Vela qF fields
explodes the Poisson block, confirming that mixing the low-density electrostatic
state with Vela's high-density qF branch is not a viable solution path.

Interpretation:

- The branch mismatch should not be attacked as a contact `phin` pinning bug or
  global quasi-Fermi gauge issue.
- The immediate root-cause target is the coupled Poisson/electrostatic shape
  that Vela accepts under its no-impact carrier state. Vela's electron
  continuity residual also accepts the Sentaurus low-density `phin` field, so
  the next debug task should isolate why the coupled Newton solve moves `psi`
  away from the Sentaurus no-avalanche electrostatic shape after the `-12.7 V`
  to `-13.2 V` transition.
- Recommended next experiment: rerun or reconstruct the `-12.7 V -> -13.2 V`
  no-impact transition with a Poisson-state guard or frozen-carrier Poisson
  correction step, then compare whether keeping Sentaurus-like `psi` prevents
  the `psi-phin` exponent jump without violating contact `phin`.

## Task 21: Reconstruct No-Avalanche Poisson Branches With Frozen Carriers

Task 20 showed that Vela's continuity residual does not reject the Sentaurus
no-avalanche quasi-Fermi fields at `-13.2 V`, but that the Sentaurus state has a
large Poisson residual under Vela's coupled residual evaluator. Reuse the
frozen-carrier Poisson reconstruction diagnostic against the VM no-avalanche
oracle to decide whether this is a Poisson matrix/box-geometry representation
failure or a carrier-charge branch selection problem.

- [x] **Step 1: Run the frozen-carrier reconstruction for VM no-avalanche and Vela no-impact**

Command:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\diagnose_pn2d_bv_poisson_reconstruction.py `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\vela\mesh.json `
  --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\vela\doping.csv `
  --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_bv_vm_no_avalanche\exports `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\bv_noimpact_resume_minus8_minus10_minus13p2\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\vm_noavalanche_vs_vela_noimpact_poisson_reconstruction `
  --biases -13.2 `
  --focus-nodes 955,1089,351,986,349,352,987,988 `
  --bc-sources vela_expected,vela_state,sentaurus_state `
  --charge-sources depletion,vela_frozen,sentaurus_frozen
```

Output:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\vm_noavalanche_vs_vela_noimpact_poisson_reconstruction
```

- [x] **Step 2: Compare reconstructed branch medians**

At `-13.2 V`, using the same Vela Poisson matrix and replacing only the frozen
charge source:

```text
bc source        charge source       band            reconstructed - Vela   reconstructed - Sentaurus
vela_expected    vela_frozen         pre_junction_p  +9.72e-5 V            +0.19485 V
vela_expected    vela_frozen         junction        +8.52e-5 V            +0.13709 V
vela_expected    vela_frozen         post_junction_n +7.23e-5 V            +0.07933 V
vela_expected    sentaurus_frozen    pre_junction_p  -0.20437 V            -0.00465 V
vela_expected    sentaurus_frozen    junction        -0.13691 V            +8.55e-5 V
vela_expected    sentaurus_frozen    post_junction_n -0.06946 V            +0.00482 V
sentaurus_state  vela_frozen         pre_junction_p  +0.00200 V            +0.19658 V
sentaurus_state  vela_frozen         junction        -2.57e-6 V            +0.13700 V
sentaurus_state  vela_frozen         post_junction_n -0.00203 V            +0.07742 V
sentaurus_state  sentaurus_frozen    pre_junction_p  -0.20243 V            -0.00220 V
sentaurus_state  sentaurus_frozen    junction        -0.13700 V            ~0 V
sentaurus_state  sentaurus_frozen    post_junction_n -0.07156 V            +0.00219 V
```

The Vela frozen carrier charge reconstructs the Vela no-impact potential to
`~1e-4 V` with Vela expected contact BC and to `~2 mV` with Sentaurus-state
contact BC. The Sentaurus frozen carrier charge reconstructs the Sentaurus
no-avalanche potential to a few millivolts in the pre/post-junction band medians
and essentially exactly at the compensated junction median.

- [x] **Step 3: Inspect the shoulder focus nodes**

With Sentaurus-state contact BC:

```text
charge source       node   reconstructed - Vela   reconstructed - Sentaurus
vela_frozen         955    +0.00164 V             +0.17290 V
vela_frozen         1089   -0.00162 V             +0.10150 V
vela_frozen         351    -4.44e-6 V             +0.14280 V
vela_frozen         986    -0.00010 V             +0.13933 V
sentaurus_frozen    955    -0.19158 V             -0.02032 V
sentaurus_frozen    1089   -0.08239 V             +0.02073 V
sentaurus_frozen    351    -0.13613 V             +0.00667 V
sentaurus_frozen    986    -0.13278 V             +0.00665 V
```

The shoulder nodes `955/1089` still show about `20 mV` residual branch distance
when reconstructing with Sentaurus frozen carriers, but this is much smaller
than the `~0.10-0.17 V` high-density branch gap produced by Vela frozen
carriers. The residual shoulder localization is still useful, but it is not the
first-order explanation for the BV branch divergence.

Interpretation:

- Vela's Poisson matrix/control volumes can represent both the Vela high-density
  no-impact branch and the Sentaurus low-density no-avalanche branch when the
  corresponding carrier charge is frozen.
- The dominant difference is therefore the carrier charge distribution selected
  by the coupled nonlinear solve, not contact boundary pinning, not a global
  quasi-Fermi gauge, and not a pure Poisson box-geometry inability to represent
  the Sentaurus electrostatic shape.
- The next implementation/debug target is the `-12.7 V -> -13.2 V` no-impact
  transition mechanics: identify why Newton and the sweep predictor leave the
  Sentaurus-like low-density charge branch and settle on the Vela high-density
  charge branch.
- Candidate next code experiment: add a diagnostic-only continuation guard that
  computes a frozen-carrier Poisson reconstruction from the previous accepted
  carrier state before accepting a large reverse-bias step. If the guard detects
  a `psi-phin` exponent jump or a large frozen-charge Poisson displacement, rerun
  the step with a smaller voltage increment or with the frozen-carrier Poisson
  reconstruction as the initial `psi` predictor. Keep this behind a debug/config
  flag until it proves it moves Vela toward the Sentaurus default BV path.

## Task 22: Add A Diagnostic `psi-phin` Continuation Guard

Task 21 showed that the key branch signal is the electron exponent
`psi - phin`, not a legal global `phin` shift or a pure Poisson matrix
representation error. Add the first lightweight continuation guard for that
signal before attempting a heavier frozen-carrier Poisson predictor.

- [x] **Step 1: Add RED tests for the guard contract**

Added tests in `tests/test_dc_sweep.cpp`:

- `DCSweep branch acceptance: measures psi-phin exponent jumps`
- `DCSweep: psi-phin branch guard writes jump diagnostics`
- config validation for negative `max_psi_phin_jump_V`

RED verification:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build-release --parallel --target test_dc_sweep
```

Expected and observed failure before implementation:

```text
DCSweepPoint has no member named 'psiPhinMaxJump_V'
```

- [x] **Step 2: Implement the opt-in guard**

Modified:

- `include/vela/simulation/DCSweep.h`
- `include/vela/simulation/DCSweepPredictor.h`
- `src/simulation/DCSweep.cpp`
- `docs/config_schema.md`

New config:

```json
"continuation": {
  "branch_acceptance": {
    "psi_phin_jump": true,
    "max_psi_phin_jump_V": 0.05
  }
}
```

Runtime behavior:

- For each candidate point after an accepted previous state, compute the maximum
  nodewise absolute jump in `psi - phin`.
- If the jump exceeds `max_psi_phin_jump_V`, mark the attempt as
  `branch_acceptance_failed` with reason `psi_phin_jump_exceeded`.
- Keep the feature opt-in; default sweep behavior is unchanged.
- When continuation diagnostics are enabled, append `psi_phin_max_jump_V` to the
  main sweep CSV.

- [x] **Step 3: GREEN focused verification**

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build-release --parallel --target test_dc_sweep
build-release\test_dc_sweep.exe "[continuation]"
```

Result:

```text
All tests passed (57 assertions in 5 test cases)
```

Interpretation:

- This is not yet the full frozen-carrier Poisson correction from Task 21.
- It provides a narrow, test-covered tripwire for the observed BV branch
  symptom: large spatial jumps in `psi - phin`.
- Next experiment: enable this guard on the `-12.7 V -> -13.2 V` no-impact BV
  sweep with thresholds around `0.03`, `0.05`, and `0.10 V`. If it forces
  smaller accepted steps or prevents the high-density branch jump, promote the
  guard into the BV debug deck. If it only rejects without finding an alternate
  branch, add the frozen-carrier Poisson predictor as the next continuation
  initializer.

## Final Verification Commands

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\test_impact_ionization.exe
ctest --test-dir build-release --output-on-failure -R "impact_ionization|newton|dc_sweep|reference_tcad|sentaurus_import"
python -m unittest tests.regression.test_reference_tcad_tools tests.regression.test_sentaurus_import_tools
python -m unittest tests.regression.test_sentaurus_vm_reference_runner
```

Then run the full suite if the promoted JSON or C++ physics implementation changed:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
ctest --test-dir build-release --output-on-failure
```

Expected: focused tests pass before promotion; full suite passes before merging any C++ solver or physics change.

### Execution Note 2026-06-19: Final Verification

Final verification was run in the MSYS2 UCRT64 release build:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build-release --parallel --target test_dc_sweep test_impact_ionization vela_example_runner
build-release\test_impact_ionization.exe
ctest --test-dir build-release --output-on-failure -R "impact_ionization|newton|dc_sweep|reference_tcad|sentaurus_import"
python -m unittest tests.regression.test_reference_tcad_tools tests.regression.test_sentaurus_import_tools
ctest --test-dir build-release --output-on-failure
```

Results:

- `build-release\test_impact_ionization.exe`: passed `123 assertions in 13 test cases`.
- Focused CTest subset: passed `5/5` tests.
- Python regression unittest suite: passed `87` tests.
- Full CTest suite: passed `316/316` tests.

### Execution Note 2026-06-19: Task 21 Focused Verification

After adding the VM no-avalanche frozen-carrier reconstruction note, the focused
tooling checks were rerun in the MSYS2 UCRT64 environment:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python -m unittest tests.regression.test_reference_tcad_tools tests.regression.test_sentaurus_import_tools tests.regression.test_sentaurus_vm_reference_runner
ctest --test-dir build-release --output-on-failure -R "reference_tcad|sentaurus_import"
rg -n "^- \[ \]" docs\superpowers\plans\2026-06-17-pn2d-bv-minus-20v-validation.md
```

Results:

- Python regression unittest suite: passed `90` tests.
- Focused CTest subset: passed `2/2` tests.
- Unchecked plan item scan: no output, meaning no remaining `- [ ]` tasks.

### Execution Note 2026-06-19: Task 22 Focused Verification

After adding the opt-in `psi-phin` continuation guard, the focused C++ sweep
checks were rerun in the MSYS2 UCRT64 release build:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
build-release\test_dc_sweep.exe
ctest --test-dir build-release --output-on-failure -R "dc_sweep|reference_tcad|sentaurus_import"
rg -n "^- \[ \]" docs\superpowers\plans\2026-06-17-pn2d-bv-minus-20v-validation.md
```

Results:

- `build-release\test_dc_sweep.exe`: passed `819` assertions in `45` test cases.
- Focused CTest subset: passed `3/3` tests.
- Full CTest suite: passed `318/318` tests.
- Unchecked plan item scan: no output, meaning no remaining `- [ ]` tasks.

### Execution Note 2026-06-19: Task 23 Contact BGN Offset Experiment

The low-bias BV mismatch was narrowed to Anode contact `ni_eff` / built-in
potential parity. The contact-density diagnostic found:

- Anode minority-electron density is low by about `0.431 decade` at -0.5 to
  -2 V.
- About `0.219 decade` comes from lower Vela `ni_eff`.
- About `0.0126-0.0130 V` lower Vela electron `psi - phin` drive contributes
  another `0.212-0.219 decade`.

Added a sensitivity-matrix case:

```json
{
  "name": "bgn_contact_ni_match_offset_0p0197",
  "solver": {
    "bandgap_narrowing": {
      "model": "old_slotboom",
      "offset_eV": 0.0197
    }
  }
}
```

Generated and ran:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\run_pn2d_solver_sensitivity_matrix.py --base-config build-release\reference_tcad\pn2d_sentaurus2018\reports\low_bias_solver_sensitivity\baseline\simulation.json --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\low_bias_solver_sensitivity\bgn_contact_ni_offset_sensitivity --bias-window reverse_low
build-release\vela_example_runner.exe --config build-release\reference_tcad\pn2d_sentaurus2018\reports\low_bias_solver_sensitivity\bgn_contact_ni_offset_sensitivity\bgn_contact_ni_match_offset_0p0197\simulation.json
python scripts\compare_pn2d_bv_multibias_fields.py --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\low_bias_solver_sensitivity\bgn_contact_ni_offset_sensitivity\bgn_contact_ni_match_offset_0p0197\vtk --curve-reference build-release\reference_tcad\pn2d_sentaurus2018\reference_curves\pn2d_sentaurus2018_bv_reference.csv --curve-candidate build-release\reference_tcad\pn2d_sentaurus2018\reports\low_bias_solver_sensitivity\bgn_contact_ni_offset_sensitivity\bgn_contact_ni_match_offset_0p0197\iv.csv --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\low_bias_solver_sensitivity\bgn_contact_ni_offset_sensitivity\bgn_contact_ni_match_offset_0p0197_compare --biases 0,-0.5,-2 --quantities potential,electric_field,electron_density,hole_density,electron_mobility,hole_mobility,avalanche_generation
```

The full -5 V run timed out after 180 s but had already reached about -4.2 V,
so the planned 0, -0.5, and -2 V comparison points were available.

Result:

- IV improved from `-0.22996` to `0.02546` decade at -0.5 V.
- IV improved from `-0.18716` to `0.03885` decade at -2.0 V.
- 0 V potential RMS improved from `0.01149 V` to `1.0e-5 V`.
- -0.5 V electron-density p95 improved from `0.52358 decade` to
  `0.08034 decade`.
- -2.0 V electron-density p95 improved from `0.48880 decade` to
  `0.06016 decade`.

Interpretation:

- This is the strongest evidence so far that the low-bias BV mismatch is mainly
  a BGN/contact built-in-potential mismatch, not terminal current extraction,
  mobility, or nonlinear damping.
- `offset_eV = 0.0197` slightly overcorrects the IV, so the next action is a
  narrow offset scan, approximately `0.017-0.0197 eV`, limited to 0, -0.5, and
  -2 V.
- After selecting the best offset, rerun contact BGN density diagnostics and
  then return to avalanche/source current-density coupling.

### Execution Note 2026-06-19: Task 24 Sentaurus OldSlotboom Split Semantics

The offset scan was superseded by a direct `models.par` / Sentaurus-manual
formula verification. Sentaurus PN2D uses `EffectiveIntrinsicDensity(
OldSlotboom )` with no Fermi correction, and the observed `ni_eff` is matched
by a split interpretation:

- `Bandgap.dEg0(OldSlotboom) = -1.595e-2 eV` is handled through the material
  intrinsic density override, yielding `Si ni = 1.4638914958767616e10 cm^-3`.
- The runtime `old_slotboom` BGN term uses only the positive Slotboom term
  `Ebgn * (ln(N/Nref) + sqrt(ln(N/Nref)^2 + C))`.
- At `N = 1e17 cm^-3`, the positive term is `0.00636396103068 eV`, producing
  `ni_eff ~= 1.65563e10 cm^-3`, matching the Sentaurus contact/global inferred
  `ni_eff ~= 1.65562e10 cm^-3`.

Implemented:

- `src/physics/BandgapNarrowing.cpp`: `old_slotboom` default `offset` is now
  `0.0` instead of `-1.595e-2`.
- `reference_tcad/pn2d_sentaurus2018/pn2d_sentaurus2018_reference.json`: BV now
  uses `vela_materials_file = pn2d_sentaurus2018_iv_materials.json`, matching
  the IV material alignment.
- `scripts/diagnose_pn2d_bv_sg_avalanche_edges.py`: Python SG/avalanche
  diagnostics now use the same positive OldSlotboom BGN term.
- `scripts/diagnose_sentaurus_ni_bgn_formula.py`: added a reproducible
  formula-hypothesis report.

Verification:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure -R "bgn|BGN|Slotboom|OldSlotboom|bandgap"
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_sentaurus2018_iv_uses_models_par_alignment tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_sentaurus2018_bv_uses_full_sentaurus_physics tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_solver_sensitivity_matrix_writes_configs
python -m py_compile scripts/run_pn2d_solver_sensitivity_matrix.py scripts/diagnose_sentaurus_ni_bgn_formula.py scripts/diagnose_pn2d_bv_sg_avalanche_edges.py scripts/diagnose_pn2d_bv_local_avalanche_factors.py scripts/diagnose_pn2d_bv_continuity_feedback.py
python scripts\sentaurus_import.py reference --config reference_tcad\pn2d_sentaurus2018\pn2d_sentaurus2018_reference.json --source-dir reference_tcad\pn2d_sentaurus2018\source --output-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke --tdr-importer build-release\sentaurus_import.exe --runner build-release\vela_example_runner.exe --skip-vela-run
build-release\vela_example_runner.exe --config build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_low_bias\simulation.json
python scripts\compare_pn2d_bv_multibias_fields.py --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_low_bias\vtk --curve-reference build-release\reference_tcad\pn2d_sentaurus2018\reference_curves\pn2d_sentaurus2018_bv_reference.csv --curve-candidate build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_low_bias\iv.csv --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_low_bias_compare --biases 0,-0.5,-2,-5 --quantities potential,electric_field,electron_density,hole_density,electron_mobility,hole_mobility,avalanche_generation
```

Results:

- Full build succeeded.
- BGN-focused CTest subset passed `10/10`.
- Focused Python regression tests passed `3/3`.
- Importer smoke generated `vela/pn2d_sentaurus2018_iv_materials.json` and a
  BV deck with both `materials_file` and `old_slotboom`.
- Official low-bias BV run converged from `0` to `-5 V` with `101` points.
- IV log errors: `0.02544 dex` at `-0.5 V`, `0.03885 dex` at `-2 V`,
  `0.03135 dex` at `-5 V`.
- 0 V field parity is restored: potential RMS `1.09e-5 V`, electron-density
  log-p95 `0.000494 decade`, hole-density log-p95 `0.000495 decade`.

Decision:

- Do not promote the `offset_eV = 0.0197` hack or spend another iteration on
  offset scanning. It was an empirical stand-in for the now-implemented
  Sentaurus split OldSlotboom semantics.
- The low-bias/BGN branch is closed for now. The next debug branch returns to
  the high-reverse-bias transition: extend the official split-material BV deck
  through the `-10 V` to `-13.2 V` region and re-run the carrier-density,
  quasi-Fermi anchoring, and no-impact/branch-control diagnostics against this
  corrected material baseline.

### Execution Note 2026-06-19: Task 25 Corrected High-Reverse Transition

Extended the official split-material BV deck from `0 V` through `-13.2 V`:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
build-release\vela_example_runner.exe --config build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_high_transition\simulation.json
python scripts\compare_pn2d_bv_multibias_fields.py --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_high_transition\vtk --curve-reference build-release\reference_tcad\pn2d_sentaurus2018\reference_curves\pn2d_sentaurus2018_bv_reference.csv --curve-candidate build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_high_transition\iv.csv --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_high_transition_compare --biases -5,-10,-13.2 --quantities potential,electric_field,electron_density,hole_density,electron_mobility,hole_mobility,avalanche_generation
python scripts\diagnose_pn2d_bv_qf_anchor.py --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_high_transition\vtk --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_high_transition_qf_anchor --biases -10,-13.2 --focus-nodes 202,199,351,986
python scripts\diagnose_pn2d_bv_continuity_feedback.py --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\doping.csv --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_high_transition\vtk --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_high_transition_continuity_feedback --biases -10,-13.2 --edge-id 2886
```

Results:

- The corrected impact-enabled run converged with `288` points.
- IV parity remains good through `-10 V`: `0.03135 dex` at `-5 V` and
  `0.04966 dex` at `-10 V`.
- At `-13.2 V`, Vela jumps to `-5.562543e-14 A/um` while Sentaurus is
  `-8.384721e-17 A`, giving `2.82178 dex` current error.
- `debug_ranking.json` recommends the order `electron_density`,
  `electron_mobility`, `electric_field`, then thresholded
  `avalanche_generation`.
- At `-13.2 V`, field comparison reports electron-density log-p95 error
  `3.45809 dex`, hole-density log-p95 error `2.27106 dex`, potential RMS
  `0.16437 V`, and electric-field relative-p95 error dominated by near-zero
  reference regions.
- QF anchoring shows contact boundary values are still aligned. The Anode
  median offsets remain only `delta_psi ~= -4.90e-5 V` and
  `delta_phin ~= 9e-15 V` at `-13.2 V`.
- The high-bias error is internal: QF band summaries show
  `delta(psi-phin) ~= 0.2067 V` in `pre_junction_p`,
  `0.1846 V` in `junction`, and `0.1696 V` in `post_junction_n`, producing
  electron-density median ratios of `3.47`, `3.10`, and `2.85 dex`.
- Focus node `202` at `x=0.75 um, y=0` has
  `delta(psi-phin)=0.2942 V`, explaining a local electron-density ratio of
  `4.943 dex`.
- Around SG edge `2886`, Vela/Sentaurus edge drops remain close even at
  `-13.2 V` (`vela_dpsi=0.35482 V` versus
  `sentaurus_dpsi=0.35853 V`), so the current best hypothesis is an internal
  state/branch shift rather than a single-edge SG drop formula error.

Ran a corrected-material no-impact branch-control variant:

```powershell
build-release\vela_example_runner.exe --config build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_noimpact_high_transition\simulation.json
python scripts\compare_pn2d_bv_multibias_fields.py --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_noimpact_high_transition\vtk --curve-reference build-release\reference_tcad\pn2d_sentaurus2018\reference_curves\pn2d_sentaurus2018_bv_reference.csv --curve-candidate build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_noimpact_high_transition\iv.csv --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_noimpact_high_transition_compare --biases -10,-13.2 --quantities potential,electric_field,electron_density,hole_density,electron_mobility,hole_mobility
```

No-impact result:

- The no-impact run converged but required `1502` continuation points, confirming
  a difficult high-bias transition even without avalanche feedback.
- At `-13.2 V`, no-impact Vela is still high by `2.66838 dex` versus the
  Sentaurus default curve and still shows electron-density log-p95 error
  `3.44638 dex`.
- Its QF band signature is nearly the same as the impact-enabled run:
  `delta(psi-phin)` is `0.2060 V` in `pre_junction_p`, `0.1795 V` in
  `junction`, and `0.1607 V` in `post_junction_n`.

Decision:

- The next branch should prioritize continuation/nonlinear branch control over
  avalanche coefficient tuning. Avalanche feedback increases the final current
  error from `2.668 dex` to `2.822 dex`, but it does not create the core
  internal electron-density/QF branch shift.
- Focus the next minimum experiment on reproducing Sentaurus Bank-Rose-like
  guarded nonlinear progression: add diagnostics or a temporary option that
  rejects a step when internal `delta(psi-phin)` or electron-density log change
  crosses a threshold in the `pre_junction_p`, `junction`, or
  `post_junction_n` bands between `-10 V` and `-13.2 V`.
- Compare three branch-control variants from `-10 V` to `-13.2 V` on the
  corrected material baseline: stronger line-search/damping, explicit
  quasi-Fermi update cap, and frozen/Poisson-guarded carrier predictor. Measure
  whether they keep `delta(psi-phin)` below `~0.02 V` and electron-density
  log-p95 below `~0.2 dex` without degrading the already-good `-5 V` to
  `-10 V` IV parity.

### Execution Note 2026-06-19: Task 26 Branch-Control Window from -10 V

Converted the corrected impact-enabled `-10 V` VTK state into a restart CSV:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_control_from_m10\restart_from_official_split_m10.csv
```

Generated a focused `-10 V -> -13.2 V` branch-control window under:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_control_from_m10
```

Window cases:

- `restart_control`: no branch-control change, used to validate the VTK restart.
- `branch_guard_0p02`: `sweep.continuation.branch_acceptance.psi_phin_jump`
  with `max_psi_phin_jump_V = 0.02`.
- `qf_hard_limit_0p0259`: Newton `quasi_fermi_update_limit_V = 0.0259`.
- `damped_update_0p5`: `damping_factor = 0.5`, `max_update = 0.5`,
  `line_search = true`.

Results:

- `restart_control` converged with `88` points and reproduced the full-run
  high-bias mismatch. At `-13.2 V`, IV error is `2.8217847749 dex`;
  electron-density log-p95 is `3.4580883684 dex`, hole-density log-p95 is
  `2.2710588194 dex`, and potential RMS is `0.164372992 V`.
- `branch_guard_0p02` converged with `414` points. The guard reduced local
  accepted jumps (`max psi_phin_max_jump_V = 0.01915`,
  `p95 = 0.00776`) but did not change the final branch. At `-13.2 V`, IV and
  field metrics are numerically unchanged from `restart_control`.
- QF-anchor summaries for `branch_guard_0p02` are identical to the control:
  `delta(psi-phin)` remains `0.206665 V` in `pre_junction_p`, `0.184600 V`
  in `junction`, and `0.169580 V` in `post_junction_n`, producing electron
  density median ratios of `3.47178`, `3.10106`, and `2.84877 dex`.
- `qf_hard_limit_0p0259` reached only `-12.9379802219 V` within the first
  300 s window. Resuming from its `last_state.csv` immediately failed at the
  same bias with `max_iterations` after 40 Newton iterations. The failure
  diagnostics show the remaining residual is Poisson-dominated
  (`psi ~= 1.65e-8`, `phin ~= 2e-13`, `phip ~= 2.5e-10`), so the hard QF
  cap stabilizes carrier updates but stalls the coupled Poisson correction.
- `damped_update_0p5` did not reach the transition region in 300 s; it slowed
  to about `-10.3515 V` and left an unfinished row. This is too expensive for
  the current branch-localization loop and does not yet provide evidence that
  damping recovers the Sentaurus branch.

Decision:

- The existing adjacent-state `psi_phin_jump` guard is a diagnostic step-size
  limiter, not a branch selector. It can force many smaller accepted steps while
  preserving the same incorrect accumulated internal branch.
- QF hard limiting alone is not enough: it prevents progress before `-13.2 V`
  and leaves a Poisson-dominated residual.
- Plain damping/update caps are computationally unattractive and should not be
  the next primary path unless paired with a better acceptance target.
- Next minimum code experiment: add a branch diagnostic/acceptance metric that
  compares the current state to a reference predictor or reference band profile,
  not only to the previous accepted point. The most useful low-risk form is a
  diagnostic-only band monitor for `psi-phin` over `pre_junction_p`,
  `junction`, and `post_junction_n`, written into the sweep CSV. Then use it
  to evaluate a frozen-carrier or Poisson-guarded predictor that targets
  absolute internal branch shape rather than local step size.

### Execution Note 2026-06-19: Task 27 Intermediate-State Residual Localization

Imported the locally available Sentaurus 2018 intermediate TDR states
`pn2d_bv_multibias_0125_des.tdr` through `pn2d_bv_multibias_0132_des.tdr`
into:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports
```

Compared these actual Sentaurus intermediate fields with the corrected-material
Vela high-transition run:

```powershell
python scripts\diagnose_pn2d_bv_qf_anchor.py --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_high_transition\vtk --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_qf_anchor --biases -12.5,-12.6,-12.7,-12.8,-12.9,-13,-13.1,-13.2 --focus-nodes 202,199,351,986
python scripts\compare_pn2d_bv_multibias_fields.py --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_high_transition\vtk --curve-reference build-release\reference_tcad\pn2d_sentaurus2018\reference_curves\pn2d_sentaurus2018_bv_reference.csv --curve-candidate build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_high_transition\iv.csv --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_field_compare --biases -12.5,-12.6,-12.7,-12.8,-12.9,-13,-13.1,-13.2 --quantities potential,electric_field,electron_density,hole_density,electron_mobility,hole_mobility,avalanche_generation
```

Intermediate-state results:

- The IV transition begins between `-12.6 V` and `-12.8 V`. Vela current error
  is `0.03335 dex` at `-12.6 V`, `0.30509 dex` at `-12.7 V`,
  `1.98827 dex` at `-12.8 V`, and `2.77256 dex` at `-12.9 V`.
- Potential remains close while the carrier branch diverges: potential RMS is
  only about `0.00337 V` from `-12.5 V` through `-12.9 V`, but electron-density
  log-p95 grows from `0.3784 dex` at `-12.6 V` to `0.9708 dex` at
  `-12.7 V`, `2.7366 dex` at `-12.8 V`, and `3.4001 dex` at `-12.9 V`.
- QF band deltas show the internal branch shift first in the electron
  exponent. In `pre_junction_p`, `delta(psi-phin)` is `0.01669 V` at
  `-12.6 V`, `0.05478 V` at `-12.7 V`, `0.16366 V` at `-12.8 V`, and
  `0.20321 V` at `-12.9 V`.

Ran Sentaurus-seeded Vela single-bias solves at `-12.7 V` and `-12.8 V`:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_seeded_vela_single_bias
```

Results:

- Starting from the actual Sentaurus state does not preserve the Sentaurus
  branch. Vela converges at `-12.7 V` to `-1.564858814e-16 A/um`
  versus Sentaurus `-7.904911015e-17 A`, and at `-12.8 V` to
  `-7.268585431e-15 A/um` versus Sentaurus `-8.003381369e-17 A`.
- The seeded `-12.8 V` solve still has `delta(psi-phin)=0.15421 V` in
  `pre_junction_p`, `0.13289 V` in `junction`, and `0.11810 V` in
  `post_junction_n`.

Evaluated Vela Newton residuals for actual and hybrid external states:

```powershell
python scripts\diagnose_pn2d_bv_newton_residual_states.py --base-config build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_high_transition\simulation.json --runner build-release\vela_example_runner.exe --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_high_transition\vtk --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_residual_hybrids --states vela:-12.7,sentaurus:-12.7,hybrid_vpsi_sqf:-12.7,hybrid_spsi_vqf:-12.7,hybrid_spsi_shift_vqf:-12.7,vela:-12.8,sentaurus:-12.8,hybrid_vpsi_sqf:-12.8,hybrid_spsi_vqf:-12.8,hybrid_spsi_shift_vqf:-12.8 --nodes 202,351,986,955,1089
```

Residual localization:

- Vela's continuity residual accepts both Vela and Sentaurus quasi-Fermi fields:
  `phin/phip` block norms remain near `1e-11` or below.
- The residual is selected almost entirely by the electrostatic potential:
  `Vela psi + Sentaurus QF` has the same small Poisson block as Vela
  (`0.2218` at `-12.7 V`, `0.2273` at `-12.8 V`), while
  `Sentaurus psi + Vela QF` has the same large Poisson block as Sentaurus
  (`2.7215` at both biases).
- The Sentaurus-state Poisson residual hotspots are not at contacts. They are
  concentrated at the abrupt doping-step control volumes:
  node `955` (`x=0.875 um`, `y=0.5 um`, `pre_junction_p`) and node `1089`
  (`x=1.125 um`, `y=0`, `post_junction_n`), each with raw
  `abs_psi_residual = 0.66941`.

Ran Poisson flux-balance and frozen-charge reconstruction diagnostics:

```powershell
python scripts\diagnose_pn2d_bv_poisson_flux_balance.py --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\doping.csv --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_high_transition\vtk --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_poisson_flux --biases -12.7,-12.8 --include-contact-nodes --top 25
python scripts\diagnose_pn2d_bv_poisson_reconstruction.py --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\doping.csv --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_high_transition\vtk --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_poisson_reconstruction --biases -12.7,-12.8 --focus-nodes 202,351,986,955,1089 --bc-sources vela_expected,sentaurus_state --charge-sources vela_frozen,sentaurus_frozen
```

Poisson evidence:

- At `-12.8 V`, Sentaurus node `1089` has
  `flux_term = 5.378444745e-12 C/m` but Vela's nodal-charge RHS gives
  `charge_term = 7.171200754e-12 C/m`, leaving
  `residual = -1.792756009e-12 C/m`.
- The required net doping inferred from the Sentaurus flux is not the full
  nodal `1e17 cm^-3`: it is about `0.750006 * 1e17 cm^-3` at the y-boundary
  step nodes and `0.916674 * 1e17 cm^-3` for the adjacent interior step nodes.
- A `top=500` flux-balance scan shows a broader geometric pattern around the
  abrupt doping transitions: the `x=0.875/1.125 um` columns require ratios from
  about `0.75` to `1.1786`, while neighboring
  `x=0.890625/1.109375 um` columns require ratios from about `1.09` to `1.50`.
  This points to control-volume doping integration/projection, not a scalar
  physics parameter.

Inverse-Poisson doping proof:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\inverse_poisson_doping_probe
```

- Replacing the top-500 Sentaurus hotspot nodes with the net doping inferred
  from the Sentaurus Poisson flux reduces the `-12.8 V` Sentaurus-state Poisson
  block from `2.7215124` to `7.36225195e-4`.
- The same inverse-doping file makes the Vela `-12.8 V` state bad
  (`psi` block `2.73337172`), proving this is not a harmless tolerance or
  current-extraction artifact. The two branches are tied to different effective
  fixed-charge RHS values near the doping steps.

Decision:

- Stop prioritizing mobility, SG flux, terminal current extraction, and
  avalanche tuning as first-order causes for the `-12.7 V` to `-12.9 V`
  transition mismatch. The strongest current root cause is that Vela's Poisson
  RHS uses point-sampled nodal doping times control-volume area, while the
  Sentaurus branch behaves like a control-volume-integrated/projection-smoothed
  fixed-charge RHS around abrupt doping transitions.
- The next implementation branch should add a diagnostic and then an optional
  production path for integrated fixed doping charge. The first target is PN2D's
  piecewise-constant vertical doping transitions near `x=0.9 um` and
  `x=1.1 um`; integrate donor/acceptor density over each node control volume or
  box subvolume instead of assigning the full nodal value to the full volume.
- Verification gate for the next branch: without using inverse-Poisson
  bias-specific fitting, the corrected integrated-doping RHS should reduce the
  Sentaurus external-state Poisson block at `-12.7/-12.8 V`, keep `0 V` and
  `-5/-10 V` IV parity within the current good envelope, and delay or remove
  the `-12.8 V` carrier-density branch jump.

### Next Tasks

1. Add a reusable fixed-charge RHS diagnostic that writes, per node, the
   control-volume area, current nodal net doping, integrated net doping
   candidate, and resulting Poisson RHS delta for PN2D.
2. Implement an optional piecewise-constant doping integration mode for imported
   PN2D reference decks. Start diagnostic-only and wire it through
   `node_doping_file` or a new fixed-charge RHS input only after the generated
   effective-doping field matches the geometric ratios seen in Task 27.
3. Re-run the focused residual probes:
   `sentaurus:-12.7`, `sentaurus:-12.8`, `vela:-12.8`,
   `hybrid_vpsi_sqf:-12.8`, and `hybrid_spsi_vqf:-12.8`.
4. If Poisson residual improves without breaking low-bias parity, run a short
   corrected-material sweep from `-12.5 V` to `-12.9 V`, then compare IV,
   potential, electron/hole density, QF anchors, and Poisson flux balance.

### Execution Note 2026-06-19: Task 28 Fixed-Charge RHS And Box-Volume Probe

Added a reusable fixed-charge RHS diagnostic:

```text
scripts/diagnose_pn2d_bv_fixed_charge_rhs.py
```

It writes per-node barycentric control volume, nodal net doping, an
integrated-profile candidate, mixed/circumcentric Voronoi volume, and effective
mixed-volume doping fields. Regression coverage was added in
`tests/regression/test_reference_tcad_tools.py`.

Verification:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python -m py_compile scripts\diagnose_pn2d_bv_fixed_charge_rhs.py
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_fixed_charge_rhs_diagnostic_writes_integrated_doping
```

Results:

- Simple piecewise-constant doping integration from the imported nodal profile
  changed only the immediate junction columns and did not reduce the external
  Sentaurus-state Poisson residual. At `-12.8 V`, the Sentaurus `psi` block
  stayed large (`2.73272057` versus the original `2.7215124`).
- Reading the SDE deck showed that PN2D has a true doping step at
  `x = 1.0 um`, while the Poisson residual hotspots are at
  `x = 0.875 um` and `x = 1.125 um`, aligned with mesh refinement-transition
  box geometry rather than with the doping discontinuity itself.
- Mixed/circumcentric Voronoi node-volume ratios match the inverse-Poisson
  required-charge ratios from Task 27:
  `0.75` at the y-boundary hotspot nodes, `0.9166667` at adjacent interior
  nodes, `1.285714` on neighboring refinement columns, and `1.5` at transition
  endpoints. This identifies Vela's barycentric node volume as the Poisson RHS
  mismatch source around refinement transitions.
- A contact-preserved mixed-volume effective doping probe reduced the external
  Sentaurus-state Poisson residual at `-12.8 V` from `2.7215124` to
  `9.72147968e-4` without using bias-specific inverse fitting. The corresponding
  `hybrid_spsi_vqf_m12p8v` Poisson block was `9.7215866e-4`.
- A Sentaurus-seeded single-bias solve at `-12.8 V` with the contact-preserved
  mixed-volume RHS still converged to the Vela high-current branch:
  `-7.1841163107e-15 A/um` versus Sentaurus
  `-8.003381369e-17 A`, or `1.95310 dex` current error. Potential RMS was
  already small (`3.1927e-05 V`), but electron-density log-p95 remained
  `2.5777 dex`.
- A stronger bias-specific inverse-top500 RHS single-bias solve confirmed the
  same remaining branch issue. At `-12.8 V`, potential RMS was
  `2.31125e-05 V`, current error was `1.95315 dex`, electron-density log-p95
  was `2.57777 dex`, and QF anchoring still showed
  `delta(psi-phin)` medians of `0.15294 V` in `pre_junction_p`,
  `0.13256 V` in `junction`, and `0.11763 V` in `post_junction_n`.

Decision:

- Replace the previous "integrated nodal doping profile" hypothesis with a
  more precise geometry hypothesis: Sentaurus's box method behaves like a
  mixed/circumcentric Voronoi control volume for the fixed Poisson charge near
  refinement transitions, while Vela currently combines cotangent fluxes with
  barycentric node volumes.
- This mixed-volume correction is necessary for electrostatic parity but not
  sufficient for BV parity. Once `psi` is aligned, the remaining high-current
  error is in the coupled carrier branch: continuity/SG flux, impact source
  feedback, density-from-QF relation, or terminal-current extraction must be
  isolated with Poisson geometry fixed.

### Next Tasks After Task 28

1. Add a production, opt-in mesh/box geometry volume policy for mixed Voronoi
   node volumes instead of continuing to use `node_doping_file` as a Poisson RHS
   hack. The option must preserve contact doping for boundary conditions.
2. Add focused tests for the new volume policy on an acute refinement-transition
   mesh: barycentric remains the default, mixed Voronoi reproduces the diagnostic
   ratios, and contact-node doping values are unchanged.
3. Re-run the `-12.8 V` residual and Sentaurus-seeded single-bias probes with
   the production mixed-volume option, then compare against the temporary
   contact-preserved effective-doping probe to prove equivalence.
4. With Poisson geometry fixed, run carrier-branch isolation at `-12.8 V`:
   no-impact, frozen-impact/source-off, Sentaurus-psi fixed-state continuity
   residual, SG edge-current decomposition, mobility field comparison, and
   terminal-current extraction. Prioritize explaining the persistent
   `0.118-0.153 V` `psi-phin` band offset.
5. Only after the single-bias carrier branch is explained, run a short
   `-12.5 V` to `-12.9 V` transition sweep. Do not repeat the full
   `0 V` to high-bias sweep until the single-bias branch selector is understood.

### Execution Note 2026-06-20: Task 29 Production Mixed-Volume Policy

Implemented the mixed-volume correction as a production opt-in mesh geometry
policy instead of a `node_doping_file` workaround:

```json
"mesh_geometry": {
  "node_volume_policy": "mixed_voronoi"
}
```

Code path:

- `BoxGeometryBuilder::Options::nodeVolumePolicy` supports `barycentric`
  (default) and `mixed_voronoi`.
- `DeviceMesh::buildBoxGeometry(options)` rebuilds node volumes and edge
  couplings with the selected policy.
- `parseBoxGeometryOptions()` parses the JSON field for Poisson, DC sweep, and
  single-bias Newton/residual-probe runner paths.
- `mixed_voronoi` changes only mesh control volumes. It does not rewrite
  donor/acceptor values, so Ohmic contact doping reconstruction remains tied to
  the physical imported doping.

Verification:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build-release --target test_box_geometry test_dc_sweep vela_example_runner --parallel
build-release\test_box_geometry.exe
build-release\test_dc_sweep.exe "DCSweep: mesh_geometry node_volume_policy selects mixed Voronoi volumes"
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_newton_residual_state_diagnostic_prepares_probe_inputs
git diff --check -- include\vela\mesh\BoxGeometryBuilder.h include\vela\mesh\DeviceMesh.h src\mesh\BoxGeometryBuilder.cpp src\mesh\DeviceMesh.cpp include\vela\simulation\ConfigParsing.h src\simulation\ConfigParsing.cpp src\simulation\DCSweep.cpp src\simulation\PoissonSimulation.cpp src\tools\vela_example_runner.cpp tests\test_box_geometry.cpp tests\test_dc_sweep.cpp docs\config_schema.md
```

Production mixed-volume residual probe:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/production_mixed_volume_probe/residual_probe
```

- `sentaurus:-12.8` Poisson block: `9.721479677e-4`.
- `hybrid_spsi_vqf:-12.8` Poisson block: `9.721586590e-4`.
- `vela:-12.8` Poisson block: `2.733466280`.
- The focus nodes `955` and `1089` are no longer Sentaurus-state Poisson
  hotspots; their `abs_psi_residual` is about `1.65e-5`. The Vela high-current
  state becomes inconsistent with the corrected mixed-volume Poisson RHS at the
  same nodes (`~0.668-0.669`).

Production mixed-volume Sentaurus-seeded single-bias solve:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/production_mixed_volume_probe/sentaurus_seed_m12p8
```

- The solve converged in `7` Newton iterations, but still landed on the Vela
  high-current branch.
- Vela current at `-12.8 V`: `-7.275678063e-15 A/um`.
- Sentaurus current at `-12.8 V`: `-8.003381369e-17 A`.
- Current error: `1.95860 dex`.
- Potential RMS error: `2.38771e-5 V`.
- Electron-density log-p95 error: `2.58758 dex`; top error remains near node
  `202` at `x = 0.75 um`.
- QF band offsets remain large:
  `delta(psi-phin) = 0.15408 V` in `pre_junction_p`,
  `0.13289 V` in `junction`, and `0.11826 V` in `post_junction_n`.

Decision:

- The mixed/circumcentric Voronoi node-volume policy should be retained as the
  required Sentaurus-compatible electrostatic box geometry for PN2D BV.
- The remaining BV mismatch is no longer first-order Poisson RHS geometry. With
  potential aligned to `~2e-5 V`, the `~2.6 dex` electron-density error and
  `~0.12-0.15 V` `psi-phin` offset point to the carrier branch selected by the
  continuity equations and source/current discretization.

### Next Tasks After Task 29

1. Keep `barycentric` as Vela's default until broader regression impact is
   reviewed, but use `mesh_geometry.node_volume_policy = "mixed_voronoi"` for
   all PN2D Sentaurus2018 BV debug probes.
2. Add a carrier-branch residual diagnostic at fixed Sentaurus `psi` with
   mixed-volume geometry. Report electron and hole continuity residuals by SG
   edge, recombination/source term, and node volume contribution at nodes
   `202`, `351`, `955`, `986`, and `1089`.
3. Run three `-12.8 V` Sentaurus-seeded mixed-volume single-bias variants:
   impact disabled, avalanche source frozen/off, and mobility/high-field drive
   simplified. Accept or reject each hypothesis by whether the
   `pre_junction_p` `delta(psi-phin)` drops below `0.02 V`.
4. Compare terminal-current extraction on the same mixed-volume solved state:
   contact SG current, integrated continuity residual, and Sentaurus reference
   current. This should be treated as secondary unless it can explain the
   existing QF/density branch shift.
5. If none of the one-factor carrier variants recovers the Sentaurus branch,
   inspect the density-from-QF relation and contact quasi-Fermi anchoring with
   BGN/`ni_eff` frozen to Sentaurus-exported values at `-12.8 V`.

### Execution Note 2026-06-20: Task 30 Mixed-Volume Carrier One-Factor Probes

Ran Sentaurus-seeded `-12.8 V` single-bias probes with production mixed-volume
geometry from the same imported Sentaurus restart state:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/production_mixed_volume_probe/carrier_one_factor_m12p8
```

Cases:

- `impact_default`: production mixed-volume baseline with Van Overstraeten SG
  edge-current avalanche.
- `no_impact`: `solver.impact_ionization.model = "none"`.
- `impact_low_field_masetti`: keep impact enabled but set mobility model to
  low-field Masetti.
- `impact_electric_field_drive`: keep high-field mobility but drive it with
  electric field instead of quasi-Fermi gradient.
- `no_impact_low_field_masetti`: lower-bound combination of no impact and
  low-field Masetti.

Summary at `-12.8 V`:

```text
case,current_error_dex,electron_log10_p95,pre_delta_psi_phin_V,junction_delta_psi_phin_V,post_delta_psi_phin_V
impact_default,1.95860,2.58758,0.15408,0.13289,0.11826
no_impact,1.81393,2.58126,0.15373,0.12826,0.10981
impact_low_field_masetti,2.12185,1.67926,0.10170,0.05165,0.05392
impact_electric_field_drive,1.65524,2.28360,0.13593,0.12244,0.10361
no_impact_low_field_masetti,1.97667,1.63539,0.10170,0.04598,0.04551
```

Findings:

- Disabling avalanche does not recover the Sentaurus carrier branch. The
  `pre_junction_p` offset changes only from `0.15408 V` to `0.15373 V`, and
  electron-density log-p95 remains about `2.58 dex`.
- Changing the high-field mobility drive to electric field helps modestly
  (`pre_junction_p` to `0.13593 V`) but remains far outside the `0.02 V`
  branch-parity target.
- Low-field Masetti reduces the QF offset more strongly
  (`pre_junction_p` to `0.10170 V`, junction to `0.05165 V`) but terminal
  current becomes worse (`2.12 dex` for impact-enabled, `1.98 dex` for the
  no-impact combination). Mobility is therefore a branch modulator, not the
  root selector.
- The production mixed-volume residual probe already showed both the actual
  Sentaurus state and `Sentaurus psi + Vela QF` have very small carrier
  continuity block residuals. With fixed Sentaurus `psi`, the continuity block
  is not uniquely selecting the Sentaurus QF branch.
- The remaining branch is therefore likely in the coupled Poisson-carrier
  correction: a small mixed-volume Poisson residual correction is coupled to
  an exponential density-from-QF response that moves the Newton solve onto the
  high-density branch.

Rejected primary causes for the mixed-volume `-12.8 V` branch jump:

- Avalanche source feedback as the first-order selector.
- High-field mobility driving-force choice as the first-order selector.
- High-field mobility limiting as a sufficient fix.
- Terminal current extraction as the primary mismatch, because field/QF errors
  persist before any current-extraction interpretation.

### Next Tasks After Task 30

1. Add a diagnostic that evaluates the first Newton correction from the actual
   Sentaurus `-12.8 V` state under production mixed-volume geometry. Report
   block-wise update norms and band medians for `delta psi`, `delta phin`,
   `delta phip`, and the induced `delta(psi-phin)` before line search.
2. Repeat that first-step diagnostic with no impact, low-field Masetti, and
   no-impact low-field Masetti. If the first Newton step already points toward
   the high-density branch in all cases, the target becomes the coupled
   Jacobian/variable scaling rather than avalanche or mobility models.
3. Compare analytic and finite-difference Jacobian-vector products on the
   Sentaurus state for selected perturbation directions, not a full finite-
   difference matrix. Prioritize directions that preserve contacts and perturb
   `psi`, `phin`, and `psi-phin` in `pre_junction_p` and `junction`.
4. Add a "Poisson-only correction with QF frozen" probe and a "continuity-only
   correction with psi frozen" probe if the first-step diagnostic confirms the
   coupled update direction is the branch selector.
5. Only after the first-step direction is understood should a solver-side
   remedy be considered, such as Bank-Rose-style trust-region damping on the
   carrier-density exponent or a Sentaurus-like staged Poisson/carrier update.

### Execution Note 2026-06-20: Task 31 Sentaurus-State First Newton Step Probe

Added a one-step Newton diagnostic for external imported states:

- `NewtonSolver::evaluateStep()` evaluates the residual, assembles the same
  analytic Jacobian used by the coupled solve, solves one Newton correction,
  applies the configured update caps, and reports the trial residual without
  running line search or iteration.
- `vela_example_runner` now accepts `simulation_type = "newton_step_probe"` and
  writes per-node `delta_psi`, `delta_phin`, `delta_phip`,
  `delta(psi-phin)`, `delta(phip-psi)`, trial densities, and initial/trial
  residual blocks.

Ran the probe from the actual Sentaurus `-12.8 V` imported state with
production mixed-volume geometry:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/production_mixed_volume_probe/newton_first_step_m12p8
```

Generated summaries:

- `newton_step_case_summary.csv`
- `newton_step_band_summary.csv`
- `newton_step_focus_nodes.csv`

Case-level first-step behavior:

```text
case,initial_residual,trial_residual,residual_ratio,max_abs_delta_psi_V,max_abs_delta_phin_V,max_abs_delta_psi_minus_phin_V
impact_default,9.721479677e-4,9.710213886e-4,0.998841,1.62e-8,0.12926,0.12926
no_impact,9.721479677e-4,9.710213859e-4,0.998841,1.61e-8,0.12926,0.12926
impact_low_field_masetti,9.721479677e-4,9.699780121e-4,0.997768,3.16e-8,0.12926,0.12926
impact_electric_field_drive,9.721479677e-4,9.712368906e-4,0.999063,1.31e-8,0.12926,0.12926
no_impact_low_field_masetti,9.721479677e-4,9.699779958e-4,0.997768,3.09e-8,0.12926,0.12926
```

Band-median `delta(psi-phin)` induced by the first step:

```text
case,pre_junction_p,junction,post_junction_n
impact_default,+5.51066e-3 V,+2.42248e-3 V,+1.31128e-3 V
no_impact,+5.43398e-3 V,+2.02417e-3 V,+9.45495e-4 V
impact_low_field_masetti,+2.49360e-3 V,+3.64795e-4 V,+4.03164e-4 V
impact_electric_field_drive,+4.62185e-3 V,+2.05149e-3 V,+1.13243e-3 V
no_impact_low_field_masetti,+2.48058e-3 V,+2.84016e-4 V,+2.77930e-4 V
```

Interpretation:

- The first coupled Newton correction barely moves electrostatic potential
  (`|delta psi| <= 3.2e-8 V`), while it immediately moves electron QF in the
  high-density direction by increasing `psi-phin`.
- Default mobility/impact increases the p-side electron exponent by about
  `0.093 dex` in one capped step; low-field Masetti reduces this to about
  `0.042 dex`, matching Task 30's conclusion that mobility changes the branch
  slope but does not restore Sentaurus parity.
- Disabling avalanche has almost no effect on the first-step electron-QF
  direction. The first-step branch selector is therefore not avalanche source
  feedback.
- The trial residual improves only slightly (`0.1-0.2%`), so the accepted
  direction is a weak residual descent but a persistent exponential density
  drift. Repeated coupled steps can plausibly accumulate into the observed
  `0.10-0.15 V` `psi-phin` offset.

Decision:

- Treat the next root-cause branch as the coupled Newton linearization and
  scaling/conditioning of the carrier-density exponent, not as a one-factor
  avalanche or terminal-current extraction problem.
- Do not tune Bank-Rose/trust-region damping yet. Damping may improve
  stability and slow the drift, but it would mask whether the Vela Jacobian
  direction is physically correct.

### Next Tasks After Task 31

1. Add a selected-direction Jacobian-vector diagnostic around the Sentaurus
   `-12.8 V` state. Compare analytic `J*v` against finite-difference residual
   differences for contact-preserving perturbations in:
   `psi`, `phin`, `phip`, `psi-phin`, and `phip-psi` over
   `pre_junction_p`, `junction`, and `post_junction_n`.
2. Add a Poisson-only correction probe with QF frozen and a continuity-only
   correction probe with `psi` frozen. This should identify whether the
   high-density drift is created by off-diagonal Poisson-carrier coupling or
   by the carrier blocks alone.
3. Inspect scaling and cap interaction for `quasi_fermi_update_limit_V`.
   The raw first-step norm is enormous (`~9.3e3-2.3e4` scaled units) and the
   capped step is still dominated by `phin`; confirm whether this cap is
   shaping the direction or only limiting magnitude.
4. If the analytic JVP agrees with finite differences, prototype a
   Sentaurus-like nonlinear safeguard as an experiment only: Bank-Rose or
   trust-region damping on carrier exponent change, with acceptance criteria
   based on both residual decrease and bounded `delta(psi-phin)`.
5. If analytic JVP disagrees with finite differences, fix the Jacobian block
   first before running additional BV sweeps.

### Execution Note 2026-06-20: Task 32 Selected-Direction JVP And Step-Cap Scan

Added a selected-direction Jacobian-vector diagnostic:

- `NewtonSolver::evaluateDirectionalDerivative()` accepts a physical-voltage
  perturbation vector and compares analytic `J*v` against the central
  finite-difference residual change `(R(x+v)-R(x-v))/2`.
- `vela_example_runner` now accepts `simulation_type = "newton_jvp_probe"`.
  Probe directions support contact-preserving masks, `x_min_um/x_max_um`, and
  modes: `psi`, `phin`, `phip`, `psi_minus_phin`, and `phip_minus_psi`.

Verification tests:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
build-release\test_newton_solver.exe "NewtonSolver: evaluateDirectionalDerivative compares analytic and finite-difference Jv"
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_runner_writes_newton_jvp_probe_for_external_state
```

Ran the JVP probe from the actual Sentaurus `-12.8 V` imported state with
production mixed-volume geometry:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/production_mixed_volume_probe/jvp_m12p8
```

Directions:

- Bands: `pre_junction_p = 0.7-0.9 um`, `junction = 0.95-1.05 um`,
  `post_junction_n = 1.1-1.3 um`.
- Modes: `psi`, `phin`, `phip`, `psi_minus_phin`, `phip_minus_psi`.
- Amplitude: `1e-6 V`; contact nodes excluded.

Summary:

- All 15 selected directions match central finite differences.
- Maximum relative JVP error: `1.36167e-12`.
- This rejects an analytic-Jacobian implementation bug as the current primary
  cause of the `-12.8 V` branch drift.

Important directional stiffness observations:

```text
direction,selected_nodes,analytic_norm,analytic_psi_norm,analytic_phin_norm,analytic_phip_norm
pre_junction_p_phin,168,3.52e-16,1.74e-17,3.51e-16,1.57e-20
junction_phin,429,2.75e-17,2.67e-17,6.78e-18,2.31e-19
post_junction_n_phin,168,2.13e-4,2.07e-4,4.74e-5,8.75e-20
pre_junction_p_psi_minus_phin,168,3.64e-4,3.64e-4,1.83e-16,1.01e-17
junction_psi_minus_phin,429,4.34e-4,4.34e-4,6.50e-18,4.49e-18
```

Interpretation:

- In the p-side and junction bands, pure `phin` perturbations have nearly zero
  residual response under the Sentaurus state. The electron quasi-Fermi
  unknown is therefore close to a null/weakly constrained direction in exactly
  the region where the final Vela branch accumulates excessive electron
  density.
- `psi_minus_phin` perturbations in those same bands are dominated by the
  Poisson block because they include a `psi` component; the carrier blocks
  remain almost insensitive. This explains why the coupled linear solve can
  produce a large `phin` correction while barely moving `psi`.

Ran a first-step `max_update` cap scan:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/production_mixed_volume_probe/newton_step_cap_scan_m12p8
```

Summary:

```text
max_update,max_abs_delta_phin_V,pre_delta(psi-phin),junction_delta(psi-phin),post_delta(psi-phin),residual_ratio
0.5,0.012926,0.000551066,0.000242248,0.000131128,0.999884
1.0,0.025852,0.001102130,0.000484496,0.000262255,0.999768
2.0,0.051704,0.002204260,0.000968992,0.000524510,0.999536
5.0,0.129260,0.005510660,0.002422480,0.001311275,0.998841
10.0,0.258520,0.011021300,0.004844960,0.002622550,0.997682
```

Decision:

- `max_update` scales the high-density drift almost linearly. It is limiting
  magnitude, not changing the Newton direction.
- Bank-Rose/trust-region damping is still a plausible stability/safeguard
  feature, but it should be treated as a nonlinear globalization strategy, not
  as the root-cause fix for BV parity.
- The next root-cause target is the weakly constrained carrier-QF subspace:
  why Vela's coupled linear solve uses this near-null `phin` direction to
  reduce the Poisson-dominated residual, and whether Sentaurus's default
  Bank-Rose/nonlinear staging suppresses that branch or whether Vela is
  missing a physical continuity/current term that should constrain it.

### Next Tasks After Task 32

1. Add block-correction probes from the Sentaurus `-12.8 V` state:
   - Poisson-only correction with `phin/phip` frozen.
   - Carrier-only correction with `psi` frozen.
   - Off-diagonal-disabled variants if needed.
   Report whether the high-density `delta(psi-phin)` is produced by carrier
   blocks alone or by Poisson-carrier off-diagonal coupling.
2. Add a linear-system conditioning diagnostic for selected bands:
   row norms, diagonal entries, off-diagonal block norms, and the solved
   correction norm for `phin` in p-side/junction. The goal is to quantify the
   near-null electron-QF subspace instead of inferring it only from JVP.
3. Compare the carrier continuity residual terms for the weak bands:
   SG flux divergence, SRH/Auger terms, impact generation, and volume scaling.
   If all are physically tiny under Sentaurus, the branch difference is likely
   a continuation/globalization issue; if a term is missing or under-scaled,
   fix the discretization first.
4. Prototype Bank-Rose-style damping only after the block-correction diagnostic:
   acceptance should constrain both residual decrease and band
   `delta(psi-phin)`, and the experiment should be compared against Sentaurus
   `-12.8 V` QF/density, not just convergence stability.

### Execution Note 2026-06-20: Task 33 Block-Correction Probe

Added a block-step Newton diagnostic:

- `NewtonSolver::evaluateBlockStep(state, "poisson_only")` solves only the
  `psi` block `J_psi,psi * delta_psi = -R_psi` with `phin/phip` frozen.
- `NewtonSolver::evaluateBlockStep(state, "carrier_only")` solves only the
  carrier block `J_carrier,carrier * [delta_phin, delta_phip] = -R_carrier`
  with `psi` frozen.
- `vela_example_runner` now accepts
  `simulation_type = "newton_block_step_probe"` and writes per-node deltas for
  each block mode.

Verification tests:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
build-release\test_newton_solver.exe "NewtonSolver: evaluateBlockStep freezes complementary unknown blocks"
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_runner_writes_newton_block_step_probe_for_external_state
```

Ran the block-step probe from the actual Sentaurus `-12.8 V` imported state
with production mixed-volume geometry:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/production_mixed_volume_probe/block_step_m12p8
```

Global block results:

```text
mode,raw_step_norm,step_norm,initial_combined_residual,trial_combined_residual,trial_psi_residual
poisson_only,0.0121662,0.0121662,9.721479677e-4,4.086281331e-8,4.086272288e-8
carrier_only,18407.0004,21.3312631,9.721479677e-4,9.723224563e-4,9.723224563e-4
```

Band-median `delta(psi-phin)`:

```text
mode,pre_junction_p,junction,post_junction_n
poisson_only,+1.27361e-5 V,-5.21136e-16 V,-1.27361e-5 V
carrier_only,+5.51066e-3 V,+2.42248e-3 V,+1.31128e-3 V
```

Focus-node examples:

```text
mode,node,x_um,delta_psi,delta_phin,delta(psi-phin),trial_psi_residual
poisson_only,955,0.875,+1.35569e-5,0,+1.35569e-5,-4.79589e-14
poisson_only,1089,1.125,-1.35569e-5,0,-1.35569e-5,-1.20651e-15
carrier_only,955,0.875,0,-5.13079e-3,+5.13079e-3,-1.64965e-5
carrier_only,1089,1.125,0,-1.51566e-3,+1.51566e-3,+1.64965e-5
```

Interpretation:

- The Poisson-only correction almost completely removes the mixed-volume
  Poisson residual without generating the high-density branch. Its
  `delta(psi-phin)` is only `~1e-5 V` in the side bands and essentially zero
  in the junction, equivalent to `~2e-4 dex` electron-density change.
- The carrier-only correction, even with `psi` frozen, reproduces the same
  high-density `delta(psi-phin)` direction and magnitude as the full first
  Newton step from Task 31. It also leaves the Poisson residual unchanged or
  slightly worse.
- Therefore the `-12.8 V` Sentaurus-state drift is not created by
  off-diagonal Poisson-carrier coupling. It is already present inside the
  carrier continuity sub-system under fixed electrostatic potential.
- This sharpens the root-cause target from "coupled Newton linearization" to
  "carrier continuity block is weakly constrained or missing a Sentaurus
  constraint/source/normalization in the low-density p-side/junction region."

Decision:

- Prioritize carrier-block diagnostics over Poisson/global coupling:
  SG flux divergence, recombination, impact source, mobility/current
  discretization, contact/QF boundary anchoring, and row/diagonal scaling for
  the electron continuity equation.
- Do not use Poisson-only correction as a solver remedy. It proves Poisson
  geometry is now mostly aligned, but it does not address current/BV parity.
- Bank-Rose damping can still be tested later as a globalization method, but
  it is now secondary to explaining why the carrier-only linear solve has a
  near-null high-density direction.

### Next Tasks After Task 33

1. Add a carrier-block row diagnostic for the Sentaurus `-12.8 V` state:
   for nodes/bands `pre_junction_p`, `junction`, and `post_junction_n`, report
   electron-continuity row diagonal, row norm, off-diagonal sum, RHS residual,
   solved `delta_phin`, and SG/recombination/impact term magnitudes.
2. Add a carrier-only no-source matrix comparison:
   repeat the block-step probe with impact disabled, SRH disabled, low-field
   mobility, and electric-field mobility drive, but inspect row conditioning
   and solved `delta_phin`, not just final current.
3. Compare Vela carrier-continuity discretization against Sentaurus/DEVSIM/
   Charon expectations for low-density reverse-bias regions:
   whether continuity is formulated in current density, quasi-Fermi potential,
   density, log-density, or uses additional Bank-Rose/pseudo-transient
   stabilization that changes the effective carrier-block diagonal.
4. If the carrier block is physically under-constrained, prototype a
   Sentaurus-like nonlinear stabilization in a separate experiment:
   Bank-Rose damping or pseudo-transient carrier diagonal, with a pass/fail
   gate requiring `pre_junction_p delta(psi-phin) < 0.02 V` and current error
   improvement against Sentaurus at `-12.8 V`.

### Execution Note 2026-06-20: Task 34 Carrier-Block Row Stiffness Probe

Added a carrier-row diagnostic:

- `NewtonSolver::evaluateCarrierRowDiagnostics()` assembles the same analytic
  Jacobian/residual as the coupled Newton solve, solves the carrier-only
  sub-system, and reports per-node electron/hole continuity row metrics.
- `vela_example_runner` now accepts
  `simulation_type = "newton_carrier_row_probe"`.
- The CSV includes electron/hole residuals, diagonal, row absolute sum,
  off-diagonal absolute sum, row L2 norm, raw carrier-only QF update, and
  capped carrier-only QF update.

Verification tests:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
build-release\test_newton_solver.exe "NewtonSolver: evaluateCarrierRowDiagnostics reports carrier row stiffness"
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_runner_writes_newton_carrier_row_probe_for_external_state
```

Ran the probe from the actual Sentaurus `-12.8 V` imported state with
production mixed-volume geometry:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/production_mixed_volume_probe/carrier_rows_m12p8
```

Global result:

```text
raw_carrier_step_norm = 18407.0004237
capped_carrier_step_norm = 21.3312631
electron block residual norm = 8.595601994e-11
```

Band-median electron-row result:

```text
band,e_abs_residual,e_abs_diag,e_row_abs_sum,offdiag/diag,raw_delta_phin,capped_delta_phin
pre_junction_p,4.78133e-16,2.05870e-14,5.51996e-14,1.77,-4.755215 V,-5.51066e-3 V
junction,3.53834e-16,4.75468e-14,1.33516e-13,1.80,-2.090390 V,-2.42248e-3 V
post_junction_n,8.73531e-15,1.20363e-13,4.04080e-13,1.70,-1.131515 V,-1.31128e-3 V
```

Focus-node examples:

```text
node,x_um,e_residual,e_diag,e_row_abs_sum,raw_delta_phin,capped_delta_phin
202,0.750,-1.14539e-14,-1.16362e-14,4.62357e-14,-111.540 V,-0.129260 V
351,1.000,+5.56832e-16,-6.78935e-14,1.71566e-13,-2.09012 V,-0.002422 V
955,0.875,+1.38025e-15,-1.28743e-14,3.64629e-14,-4.42742 V,-0.005131 V
986,1.008,+4.03872e-16,-4.76190e-14,1.33687e-13,-1.99696 V,-0.002314 V
1089,1.125,-7.12652e-15,-6.01815e-14,2.02040e-13,-1.30789 V,-0.001516 V
```

Top raw `delta_phin` nodes are clustered at `x = 0.75 um`; several have
`raw_delta_phin = -111.54 V`, then are clipped by `max_update = 5` to
`-0.12926 V`. This is not a normal local correction from a large residual.
The residual and diagonal are both tiny, and the solved correction is dominated
by a near-null carrier-block mode.

Ran the same carrier-row probe for the one-factor physics variants:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/production_mixed_volume_probe/carrier_rows_one_factor_m12p8
```

Summary:

```text
case,raw_norm,pre_raw_dphin,junction_raw_dphin,post_raw_dphin,pre_e_diag,junction_e_diag
impact_default,18407.0,-4.755,-2.090,-1.132,2.06e-14,4.75e-14
no_impact,18359.5,-4.689,-1.747,-0.816,2.06e-14,4.77e-14
impact_low_field_masetti,9300.0,-1.117,-0.163,-0.181,1.79e-13,7.38e-13
impact_electric_field_drive,22554.1,-4.932,-2.189,-1.208,2.06e-14,4.74e-14
no_impact_low_field_masetti,9298.1,-1.111,-0.127,-0.125,1.79e-13,7.43e-13
```

Interpretation:

- Disabling impact ionization barely changes the p-side electron-row stiffness
  or raw `delta_phin`; avalanche remains rejected as the root selector.
- Low-field Masetti increases electron-row diagonal/row sum by roughly an
  order of magnitude and reduces raw `delta_phin`, especially near the
  junction. This explains why Task 30 showed smaller QF offsets under low-field
  mobility, but it still does not restore the Sentaurus branch.
- Electric-field high-field drive does not improve conditioning; the raw
  carrier step norm is even larger than default.
- The root cause is now localized to the electron continuity block in
  low-density reverse-biased p-side/junction regions: Vela's steady
  quasi-Fermi formulation has extremely small row scale and a near-null mode
  that the linear solve turns into large negative `phin` updates.

Decision:

- Treat the carrier-block near-null mode as the primary current BV mismatch
  source at `-12.8 V`.
- The next question is whether Sentaurus avoids this via a missing physical
  term/normalization in Vela's continuity discretization or via nonlinear
  stabilization such as Bank-Rose/pseudo-transient carrier diagonal.
- A solver-side damping experiment is now justified, but it must be framed as
  a test of Sentaurus-like globalization/regularization, not as a final fix,
  until SG/recombination/impact term decomposition confirms no missing term.

### Next Tasks After Task 34

1. Add term-level carrier continuity diagnostics for selected nodes/bands:
   decompose electron row residual into SG flux divergence, SRH/Auger
   recombination, impact generation, and volume/source contributions. Compare
   term magnitudes against row diagonal and raw `delta_phin`.
2. Add an experimental carrier pseudo-transient/Bank-Rose regularization mode
   for carrier-only Newton steps only. Sweep the added diagonal or damping
   parameter at `-12.8 V` and require:
   `pre_junction_p capped delta(psi-phin) < 0.02 V`,
   no residual explosion, and improved Sentaurus current/density parity.
3. Compare Vela's carrier variable choice with Sentaurus/DEVSIM/Charon:
   determine whether Sentaurus default BV solve effectively regularizes
   quasi-Fermi potential in low-carrier regions, or solves/log-transforms a
   different carrier unknown.
4. Once a regularization hypothesis passes the single-state carrier-only gate,
   run a Sentaurus-seeded `-12.8 V` full solve and then a short
   `-12.5 V` to `-12.9 V` reverse-bias transition sweep.

### Execution Note 2026-06-20: Task 35 Carrier Continuity Term Decomposition

Added a term-level carrier-continuity diagnostic:

- `CoupledDDAssembler::carrierContinuityTermDiagnostics()` decomposes the same
  carrier residual used by the coupled solver into:
  SG flux divergence, recombination, impact generation, gauge, and boundary
  contributions.
- `NewtonSolver::evaluateCarrierTermDiagnostics()` exposes the assembler
  diagnostic with the same state scaling and boundary conditions as the
  Newton residual evaluator.
- `vela_example_runner` now accepts
  `simulation_type = "newton_carrier_term_probe"`.

The diagnostic reports scaled residual-units terms. For every row,
`flux + recombination + impact + gauge + boundary == residual`; the `-12.8 V`
production probe confirmed zero closure error in the inspected bands.

Verification tests:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
build-release\test_newton_solver.exe "NewtonSolver: evaluateCarrierTermDiagnostics decomposes continuity residual"
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_runner_writes_newton_carrier_term_probe_for_external_state
```

Ran the term probe from the actual Sentaurus `-12.8 V` imported state with
production mixed-volume geometry:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/production_mixed_volume_probe/carrier_terms_m12p8
```

Band-median electron terms:

```text
band,abs_flux,abs_recombination,abs_impact,abs_residual
pre_junction_p,3.41852e-15,7.79544e-16,3.76411e-17,4.78133e-16
junction,1.22963e-15,4.45454e-16,4.30347e-16,3.53834e-16
post_junction_n,7.79783e-15,7.79544e-16,1.09899e-16,8.73531e-15
```

Focus nodes:

```text
node,x_um,e_flux,e_recombination,e_impact,e_residual
202,0.750,-1.14534e-14,-5.09405e-19,-3.21949e-23,-1.14539e-14
351,1.000,+1.61541e-15,-5.27960e-16,-5.30618e-16,+5.56832e-16
955,0.875,+2.65764e-15,-1.22500e-15,-5.23970e-17,+1.38025e-15
986,1.008,+1.37026e-15,-4.45454e-16,-5.20934e-16,+4.03872e-16
1089,1.125,-5.70479e-15,-1.22500e-15,-1.96726e-16,-7.12652e-15
```

Interpretation:

- At the largest raw-update nodes near `x = 0.75 um`, the electron residual is
  almost entirely SG flux divergence; recombination and impact are negligible.
- In the junction, flux, recombination, and impact partially cancel, but all
  terms remain tiny (`~1e-15` scaled residual units). The carrier block can
  still solve these tiny row values into volt-scale raw `phin` changes because
  the row stiffness is also tiny.
- Impact generation is not the missing dominant term. Its median contribution
  is much smaller than flux in the p-side and post-junction bands, and no-impact
  probes preserve the near-null carrier update.

Ran one-factor term probes:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/production_mixed_volume_probe/carrier_terms_one_factor_m12p8
```

Summary:

```text
case,pre_abs_flux,junction_abs_flux,post_abs_flux,pre_abs_impact,junction_abs_impact
impact_default,3.42e-15,1.23e-15,7.80e-15,3.76e-17,4.30e-16
no_impact,3.42e-15,1.23e-15,7.80e-15,0,0
impact_low_field_masetti,3.74e-14,3.96e-14,5.58e-14,5.69e-16,1.08e-14
impact_electric_field_drive,3.38e-15,1.18e-15,9.72e-15,3.74e-17,4.29e-16
no_impact_low_field_masetti,3.74e-14,3.96e-14,5.58e-14,0,0
```

Decision:

- The missing source-term hypothesis is rejected for the current `-12.8 V`
  Sentaurus-state branch drift. Recombination and impact are too small to
  explain the p-side volt-scale raw `phin` mode.
- The dominant residual/stiffness scale is SG transport under the configured
  high-field mobility. Low-field Masetti raises the SG flux and row stiffness
  by roughly one order of magnitude, reducing but not eliminating the raw QF
  drift.
- The next highest-value experiment is therefore controlled carrier-block
  regularization/globalization: add pseudo-transient or Bank-Rose-like
  diagonal/damping to the carrier equations and test whether it reproduces
  Sentaurus's branch selection without damaging Poisson parity.

### Next Tasks After Task 35

1. Prototype an opt-in carrier pseudo-transient diagonal for diagnostic use.
   Start with carrier-only block-step probes at Sentaurus `-12.8 V`, sweeping
   the added diagonal relative to row absolute sum or node volume scale.
2. In parallel, prototype a Bank-Rose-style carrier update limiter that limits
   band `delta(psi-phin)` / carrier exponent change during Newton acceptance.
   Treat this as a globalization experiment, not a final model change.
3. Pass/fail gate for either experiment:
   - carrier-only `pre_junction_p delta(psi-phin) < 0.02 V`;
   - no residual explosion in carrier blocks;
   - Sentaurus-seeded full `-12.8 V` solve reduces electron-density log-p95
     and current error relative to the current `2.59 dex` / `1.96 dex` state.
4. If regularization works, compare its effective diagonal/update rule against
   Sentaurus Bank-Rose log behavior and DEVSIM/Charon continuation strategies.
5. If regularization does not work, inspect the high-field mobility model and
   SG flux normalization more deeply, because Task 35 shows SG transport scale
   is the dominant carrier-block stiffness source.

### Execution Note 2026-06-20: Task 36 Regularized Carrier-Step Probe

Added a diagnostic-only regularized carrier-step probe:

- `NewtonSolver::evaluateRegularizedCarrierStep(state, scale)` solves only the
  carrier sub-block with `psi` frozen.
- For each carrier row, it adds
  `sign(diagonal) * scale * row_abs_sum` to the carrier-block diagonal before
  solving. `scale = 0` is exactly the existing `carrier_only` block step.
- `vela_example_runner` now accepts
  `simulation_type = "newton_regularized_carrier_step_probe"` and writes
  per-node CSV columns for the step, trial state, residuals, and
  `regularization_diagonal_norm`.

Verification tests:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
build-release\test_newton_solver.exe "NewtonSolver: evaluateRegularizedCarrierStep damps carrier-only correction"
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_runner_writes_newton_regularized_carrier_step_probe_for_external_state
```

Ran the scale sweep from the actual Sentaurus `-12.8 V` imported state:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/production_mixed_volume_probe/regularized_carrier_step_m12p8
```

The probe swept `scale = 0, 0.01, 0.03, 0.1, 0.3, 1, 3, 10`.
Band windows use the existing BV diagnostics convention:
`pre_junction_p = 0.7..0.9 um`, `junction = 0.95..1.05 um`,
`post_junction_n = 1.1..1.3 um`.

Summary:

```text
scale  raw_step_norm  step_norm  pre d(psi-phin)  junction d(psi-phin)  post d(psi-phin)
0      1.8407e4       21.3313    +5.511e-3        +2.422e-3             +1.311e-3
0.01   1.3611e3       21.1197    +6.775e-3        +1.017e-3             +2.06e-4
0.03   3.3374e2       21.0561    +5.497e-3        +9.10e-5              -5.2e-5
0.1    5.2043e1       21.4500    +3.020e-3        -3.61e-4              ~0
0.3    9.9652         9.9652     +8.43e-4         -3.11e-4              ~0
1      2.2487         2.2487     +3.24e-4         -9.5e-5               ~0
3      7.022e-1       7.022e-1   +6.9e-5          -3.2e-5               ~0
10     2.079e-1       2.079e-1   +1.8e-5          -1.0e-5               -2.0e-6
```

Trial residual behavior:

- The global `phin` trial residual norm stays near the initial value for all
  scales (`~8.6e-11` initial; `8.53e-11` at `scale = 10`).
- In the three junction bands, `phin` trial residual ratios stay order-one and
  do not explode. Some bands improve, but the diagnostic is not intended to
  minimize the nonlinear residual by itself.

Interpretation:

- The row-absolute-sum diagonal strongly suppresses the near-null carrier mode:
  raw carrier correction drops from `1.84e4` to `0.208` as scale increases from
  `0` to `10`.
- At `scale >= 0.3`, the pre-junction electron exponent update is below
  `1 mV`; at `scale = 10`, it is only `18 uV`.
- This supports the stabilization hypothesis: Sentaurus's Bank-Rose/nonlinear
  globalization could plausibly prevent Vela's carrier block from taking the
  high-density branch in low-density reverse-biased regions.
- This does not yet prove the full BV mismatch is solved, because the
  regularized diagonal is diagnostic-only and has not been integrated into
  the accepted coupled Newton step or DC continuation.

Decision:

- Keep this probe as a diagnostic instrument.
- The next implementation experiment should be the smallest opt-in solver path
  that applies the same idea during Newton acceptance, guarded by config and
  disabled by default.
- Prefer a pseudo-transient/Bank-Rose globalization formulation over changing
  physical SG, recombination, or impact terms, because Tasks 34-36 show the
  immediate branch selector is carrier-block conditioning rather than a missing
  source term.

### Next Tasks After Task 36

1. Add an opt-in experimental carrier regularization mode to the coupled Newton
   solve path. Start diagnostic-only in config, disabled by default, and reuse
   the Task 36 row-absolute-sum diagonal for carrier rows.
2. Run a Sentaurus-seeded full `-12.8 V` single-bias solve with candidate
   scales `0.3`, `1`, `3`, and `10`. Gate on:
   - no worse Poisson residual;
   - reduced electron-density log-p95 error versus current `2.59 dex`;
   - reduced Anode current error versus current `1.96 dex`;
   - no artificial freezing of the carrier equations.
3. If the diagonal path improves full solve parity, compare it with a
   Bank-Rose-style accepted-step limiter on `delta(psi-phin)` and
   `delta(phip-psi)`. Prefer the version that preserves residual reduction
   and does not merely cap all carrier motion.
4. If neither full-solve experiment improves density/current parity, return to
   high-field mobility and SG flux normalization, using Task 35 term probes and
   Task 36 regularized-step probes as the reference diagnostics.

### Execution Note 2026-06-20: Task 37 Coupled Regularization and SG Clamp Root Cause

Implemented the first opt-in coupled Newton regularization experiment:

- New config key: `solver.carrier_regularization_scale`, default `0`.
- In the coupled Newton solve path, when the scale is positive, Vela adds
  `sign(diagonal) * scale * carrier_row_abs_sum` to each carrier-row diagonal
  before solving the full coupled Newton system.
- This is disabled by default and documented as an experimental BV
  stabilization knob in `docs/config_schema.md`.

TDD coverage added:

```powershell
build-release\test_newton_solver.exe "NewtonSolver: carrier regularization damps coupled Newton carrier mode"
build-release\test_newton_solver.exe "NewtonSolver: defaults to analytic Jacobian"
build-release\test_newton_solver.exe "NewtonSolver: parses block residual norm controls"
```

Initial full-solve experiment before fixing SG flux:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/production_mixed_volume_probe/regularized_coupled_m12p8
```

At Sentaurus-seeded `-12.8 V`, scales `0.3`, `1`, `3`, and `10` all
converged in 2 Newton iterations and reduced junction density errors
dramatically:

```text
case      junction electron log10_p95  junction hole log10_p95  current abs log10 error
baseline  2.588                      1.397                    1.959
scale_0.3 0.0458                     0.0623                   2.502
scale_1   0.0238                     0.0254                   2.500
scale_3   0.0124                     0.0117                   2.500
scale_10  0.00430                    0.00456                  2.499
```

Interpretation at this stage:

- Coupled carrier regularization fixed the high-density branch but worsened
  terminal current error. That separated the problem into two layers:
  carrier-state branch selection and terminal/SG current scale.
- Contact-edge diagnostics showed Vela IV current equals the sum of selected
  Anode contact-edge currents to roundoff, so the IV aggregation path is
  internally consistent.
- The contact-edge state at `scale_10` matched Sentaurus contact-node
  potential, carrier densities, and mobilities closely, but Vela electron edge
  flux was still much too large.

Root-cause finding:

- The large electron current came from
  `sgElectronContinuityFluxFromQuasiFermiVariableNi`.
- The old implementation evaluated the balanced quasi-Fermi form as separate
  factors such as `exp(psi/Vt)` and `exp(-phin/Vt)`, each individually
  clamped by `limitedExp`.
- At the PN2D `-12.8 V` p-contact, `psi/Vt ~= -511` is clamped to `-500`,
  while `-phin/Vt ~= 495` is not symmetrically cancelled. This creates an
  artificial factor of roughly `exp(11) ~= 5e4` in the SG flux.
- A new regression test reproduces this exact scale:

```powershell
build-release\test_sg_flux.exe "SG variable-ni quasi-Fermi flux matches density form at large absolute bias"
```

The SG fix:

- `sgElectronContinuityFluxFromQuasiFermi()` now forms physical densities from
  the combined exponent `(psi - phin) / Vt` and calls the density SG form.
- `sgHoleContinuityFluxFromQuasiFermi()` does the same with
  `(phip - psi) / Vt`.
- The variable-`ni` forms now use finite physical densities and the effective
  SG potential jump:
  - electron `eta = (psi1 - psi0)/Vt + log(ni1/ni0)`;
  - hole `eta = (psi1 - psi0)/Vt + log(ni0/ni1)`.
- Exact flat-QF fast paths preserve the existing zero-flux invariant for
  BGN/effective-`ni` edges.

Verification:

```powershell
build-release\test_sg_flux.exe
build-release\test_newton_solver.exe "NewtonSolver: carrier regularization damps coupled Newton carrier mode"
build-release\test_newton_solver.exe "NewtonSolver: evaluateRegularizedCarrierStep damps carrier-only correction"
```

Post-SG-fix full-solve experiment:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/production_mixed_volume_probe/regularized_coupled_m12p8_sgfix
```

With the corrected SG flux, the previous Sentaurus-seeded `-12.8 V`
single-point configurations no longer converge under the strict default
Newton tolerance:

```text
case       failure                  iterations  residual      psi block     phin block   phip block
baseline   line_search_non_decrease 1           1.13e-7       4.09e-8      5.55e-14     1.05e-7
scale_0.3  line_search_non_decrease 1           6.65e-8       4.09e-8      2.62e-12     5.24e-8
scale_1    line_search_non_decrease 1           4.79e-8       4.09e-8      4.32e-12     2.49e-8
scale_3    line_search_non_decrease 15          4.07e-8       4.07e-8      5.34e-12     8.45e-10
scale_10   max_iterations           40          3.71e-8       3.70e-8      5.80e-12     2.54e-9
```

Increasing `scale_10` to 200 Newton iterations reduces residual only to
`2.43e-8`, still dominated by the Poisson block. This is now a residual-floor
or globalization/tolerance problem rather than the old carrier high-density
branch drift.

Diagnostic accepted-state run:

- For diagnostic output only, `scale_10` was rerun with `abstol = 1e-7`
  and `handoff.newton_max_iter = 200`.
- It accepts after 1 Newton iteration and produces:

```text
current_total_A_per_um = -2.4659e-14
electron current       = -5.47e-19 A/um
hole current           = +2.4658e-14 A/um
Sentaurus current      = -8.0034e-17
current abs log error  = 2.489 dex
```

Field parity for that diagnostic state:

```text
quantity          field_error   junction_error
potential         2.39e-5       1.38e-5
electric_field    2.15e2        0.1116
electron_density  0.00135 dex   0.00210 dex
hole_density      0.00195 dex   0.00243 dex
```

Interpretation after SG fix:

- The electron-flux over-amplification was a real Vela bug and is now fixed by
  test-covered SG code.
- Carrier densities now match Sentaurus extremely well in the accepted
  diagnostic state, so the original high-density branch source is largely
  controlled.
- The remaining large current mismatch is no longer an electron SG clamp issue.
  It is now dominated by the p-contact hole current / terminal-current scale.
- Contact-edge aggregation remains internally consistent, so the next
  mismatch branch is Sentaurus terminal-current definition and Vela contact
  normal-flux scaling/sign/unit parity, especially comparing Sentaurus
  `ContactCurrentFlux`, `ContactCurrentDensity`, and `.plt` `Anode TotalCurrent`
  against Vela's `q * flux * edge.couple` per-depth convention.

### Next Tasks After Task 37

1. Keep the SG clamp fix. It has a concrete RED/GREEN regression and directly
   explains the old `~5e4` electron contact-flux amplification.
2. Do not tune `carrier_regularization_scale` as the final BV answer yet.
   With corrected SG flux, strict Newton convergence now fails at a Poisson
   residual floor around `2e-8..4e-8`; this needs a separate convergence
   strategy.
3. Run a dedicated Sentaurus-current-unit diagnostic:
   - parse `pn2d_bv.plt` at `-12.8 V`;
   - compare `.plt` `Anode TotalCurrent`, exported `ContactCurrentFlux`, and
     `TotalCurrentDensity`;
   - compute the geometric conversion needed to map Sentaurus current fields
     to Vela's `A/um` per-depth convention.
4. Add a Vela contact-current cross-check that reports both:
   - normal finite-volume edge current `q * flux * edge.couple`;
   - a Sentaurus-like boundary current-density integral over contact length.
5. After current-unit/contact-flux parity is resolved, rerun:
   - strict `-12.8 V` Sentaurus-seeded solve with SG fix;
   - short transition sweep around `-12.5..-12.9 V`;
   - full BV current comparison to verify whether density parity now translates
     into Sentaurus BV curve parity.

### Execution Note 2026-06-20: Task 38 SG Jacobian Consistency Follow-up

Subagent review found that Task 37 fixed the SG quasi-Fermi residual path but
left the coupled analytic Jacobian on the old separately clamped exponent
factorization. That made the residual and Jacobian inconsistent at BV absolute
potential scale (`psi ~= -13 V`, `phin/phip ~= -12.8 V`) and explained why the
post-SG-fix `-12.8 V` strict Newton run still stalled.

Code changes:

- `src/equation/CoupledDDAssembler.cpp`
  - Removed the Jacobian precomputation of separately clamped
    `exp(psi/Vt)`, `exp(-psi/Vt)`, `exp(-phin/Vt)`, and `exp(phip/Vt)`.
  - Re-derived electron and hole SG continuity derivatives from the same
    combined-exponent density form used by the residual:
    `n_i = ni_i exp((psi_i - phin_i)/Vt)` and
    `p_i = ni_i exp((phip_i - psi_i)/Vt)`.
- `src/discretization/ScharfetterGummel.cpp`
  - Added flat quasi-Fermi short-circuiting to the equal-`ni` plain SG helpers
    as well, so a large absolute potential plus electric field cannot produce
    a clamp-induced fake flat-QF flux.
- `src/solver/NewtonSolver.cpp`
  - Tightened `carrier_regularization_scale` row-sum semantics to include the
    full coupled carrier row (`psi`, `phin`, and `phip` columns), matching the
    documented `carrier_row_abs_sum` meaning.

New regression tests:

```text
build-release\test_sg_flux.exe "SG quasi-Fermi fluxes cancel flat QF at large absolute bias with electric field"
build-release\test_newton_solver.exe "CoupledDDAssembler: analytic Jacobian matches finite differences at BV absolute potential scale"
```

Verification run:

```text
cmake --build build-release --target test_sg_flux test_newton_solver vela_example_runner --parallel
build-release\test_sg_flux.exe
build-release\test_newton_solver.exe "CoupledDDAssembler: analytic Jacobian matches finite differences*"
build-release\test_newton_solver.exe "NewtonSolver: carrier regularization damps coupled Newton carrier mode"
python scripts\compare_pn2d_bv_multibias_fields.py --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\production_mixed_volume_probe\regularized_coupled_m12p8_sgfix_jacfix\scale_10_max200_strict\vtk --curve-reference build-release\reference_tcad\pn2d_sentaurus2018\reference_curves\pn2d_sentaurus2018_bv_reference.csv --curve-candidate build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\production_mixed_volume_probe\regularized_coupled_m12p8_sgfix_jacfix\scale_10_max200_strict\iv.csv --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\production_mixed_volume_probe\regularized_coupled_m12p8_sgfix_jacfix\scale_10_max200_strict\compare --biases -12.8 --quantities potential,electric_field,electron_density,hole_density
```

Strict `-12.8 V` Sentaurus-seeded smoke result after rebuilding
`vela_example_runner`:

```text
output:
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/production_mixed_volume_probe/regularized_coupled_m12p8_sgfix_jacfix/scale_10_max200_strict

solver result: converged = true
Newton iterations: 2
current_total_A_per_um = -1.5058175384774977e-17
electron current       = -5.4735612753122913e-19 A/um
hole current           = +1.4510819257243748e-17 A/um
Sentaurus current      = -8.00338136924e-17
current abs log error  = 0.725501160834031 dex
```

Field parity for the strict converged state:

```text
quantity          field_error   junction_error   contact-local note
potential         2.39e-5       1.38e-5          p-contact 4.90e-5
electric_field    0.7707        0.1116           p-contact 0.8273
electron_density  0.00201 dex   0.00268 dex      contacts ~floor
hole_density      0.00218 dex   0.00288 dex      contacts ~floor
```

Interpretation:

- The SG analytic Jacobian inconsistency was the immediate convergence
  blocker after Task 37. With residual and Jacobian aligned, the same
  `scale_10_max200` strict point changes from `max_iterations` at residual
  `2.43e-8` to convergence in 2 Newton iterations.
- Carrier-density parity remains excellent. The current mismatch is now
  reduced from the diagnostic accepted-state `2.489 dex` to `0.726 dex`, but
  is still large enough to keep the next debug branch on terminal current and
  contact normal-flux definition.
- Electric-field discrepancy is strongly contact/bulk weighted while the
  junction error is much smaller. This reinforces the next check: compare
  Sentaurus terminal-current fields and Vela contact boundary flux/integration
  conventions before tuning avalanche or mobility.

### Next Tasks After Task 38

1. Add a contact-current parity diagnostic that reports, at `-12.8 V`, all of:
   Vela terminal current from `post::ContactCurrent`, Vela finite-volume
   boundary-normal carrier flux per contact edge, Sentaurus `.plt`
   `Anode TotalCurrent`, exported `ContactCurrentFlux`, and exported
   `TotalCurrentDensity`.
2. Use that diagnostic to decide whether the remaining `0.725 dex` current
   error is a unit/depth conversion, contact-edge integration, sign convention,
   or actual transport-field discrepancy.
3. Re-run the strict transition sweep `-12.5..-13.2 V` with the SG Jacobian fix
   before tuning solver damping or `carrier_regularization_scale`; the
   convergence landscape changed substantially.

### Execution Note 2026-06-20: Task 39 Contact-Current Parity Diagnostic

Task 39 executed the first item after Task 38: split the remaining `0.725 dex`
current mismatch into terminal extraction, Sentaurus current-definition, and
actual transport-field branches.

Code changes:

- `scripts/diagnose_pn2d_bv_contact_current_extraction.py`
  - Added optional `.plt` parsing via `--sentaurus-plt`,
    `--plt-bias-column`, and `--plt-current-column`.
  - Added `ContactCurrentFlux` extraction from Sentaurus
    `field_manifest.json` using `region_name == contact`.
  - Added endpoint QF/potential drop diagnostics for contact edges:
    `sentaurus_hole_qf_drop_V`, `vela_hole_qf_drop_V`,
    `sentaurus_over_vela_hole_qf_drop`, and electron/potential analogues.
- `tests/regression/test_reference_tcad_tools.py`
  - Extended the contact-current extraction tests to cover `.plt`,
    `ContactCurrentFlux`, manifest `fields`, and endpoint QF-drop output.

Diagnostic output:

```text
runner output:
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/production_mixed_volume_probe/regularized_coupled_m12p8_sgfix_jacfix/scale_10_max200_strict_contact_diag

Anode summary:
build-release/reference_tcad/pn2d_sentaurus2018/reports/cc_diag_jacfix_m12p8

Cathode control:
build-release/reference_tcad/pn2d_sentaurus2018/reports/cc_diag_jacfix_m12p8_cathode
```

Verification:

```text
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_contact_current_extraction_filters_edges_and_classifies_parity tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_contact_current_extraction_aligns_sentaurus_nearest_fields
```

Key `-12.8 V` Anode findings:

```text
Vela edge sum current       = -1.505817538477498e-17 A/um
Vela IV current             = -1.5058175384774977e-17 A/um
edge vs IV relative error   = 2.05e-16
Sentaurus reference curve   = -8.00338136924e-17 A
Sentaurus .plt Anode total  = -8.00338136923547e-17 A
Sentaurus TDR contact flux  = -8.472042312441417e-17 A
.plt vs TDR flux rel error  = 5.53e-2
classification              = terminal_extraction_consistent
current_gap_branch          = transport_or_field_mismatch
```

Sentaurus `.plt` component split at `-12.8 V`:

```text
Anode eCurrent       = -5.4733131428545704e-19
Anode hCurrent       = -7.9486482378069298e-17
Anode TotalCurrent   = -8.0033813692354699e-17
Cathode eCurrent     = +7.9783405536963002e-17
Cathode hCurrent     = +2.5040815539147098e-19
Cathode TotalCurrent = +8.0033813692354502e-17
```

Vela terminal-balance split for the same solve:

```text
Anode electron current = -5.4735612753122913e-19 A/um
Anode hole current     = +1.4510819257243748e-17 A/um
Anode total current    = -1.5058175384774977e-17 A/um
Cathode electron       = +8.310065358305952e-17 A/um
Cathode hole           = -2.5042191671514656e-19 A/um
Cathode total          = +8.3351075499774659e-17 A/um
```

Interpretation:

- Vela contact-edge aggregation and IV extraction are not the cause: Anode
  edge-sum and IV agree to roundoff.
- Sentaurus `.plt` and TDR `ContactCurrentFlux` are mutually consistent at the
  current scale that matters here: the Anode flux differs by only `5.5%`.
- The remaining `0.725 dex` Anode current mismatch is almost entirely the
  Anode hole-current component. Vela's Anode electron current matches
  Sentaurus, and Vela's Cathode electron current is close to Sentaurus.
- On Anode contact-adjacent edges:

```text
mean Vela hole QF drop      = -1.7763568394002505e-15 V
mean Sentaurus hole QF drop = -8.881784197001252e-15 V
Sentaurus/Vela drop ratio   = 5.0
mean Vela electron QF drop  = -1.094579956765429e-3 V
mean Sentaurus electron QF  = -1.094172877079913e-3 V
electron drop ratio         = 0.9996
```

This matches the observed carrier split: the electron branch is aligned, while
the Anode hole branch is controlled by a contact-adjacent QF drop at the
`1e-15..1e-14 V` numerical floor. The next source of BV current difference is
therefore not current extraction or units, but Anode hole contact-adjacent
quasi-Fermi-gradient parity / numerical floor behavior.

### Next Tasks After Task 39

1. Re-run the strict transition sweep `-12.5..-13.2 V` with the SG Jacobian fix
   and emit contact-edge diagnostics for Anode at each bias. Track:
   - curve error;
   - Anode `hCurrent` error;
   - `sentaurus_over_vela_hole_qf_drop`;
   - Cathode electron current parity.
2. Add a read-only sensitivity probe that perturbs only the Anode
   contact-adjacent hole QF drop in the converged Vela state to the Sentaurus
   endpoint drop and recomputes contact current. This should prove whether the
   5x QF-drop floor fully explains the remaining `0.725 dex`.
3. If the sensitivity probe closes the current gap, inspect why Vela pins the
   Anode hole QF drop closer to exact flatness than Sentaurus:
   - contact Dirichlet projection precision;
   - warm-start/state import rounding;
   - Newton convergence floor and Bank-Rose damping behavior;
   - contact boundary reconstruction for ohmic minority/majority carriers.
4. Only after the QF-drop sensitivity is resolved should we tune avalanche or
   mobility parameters; current evidence says those are not the primary source
   of the `-12.8 V` mismatch.

### Execution Note 2026-06-20: Task 40 Per-Bias Transition Contact Diagnostics

The first attempt at a multi-bias transition sweep used the `-12.8 V`
Sentaurus restart for the `-12.5 V` point and failed immediately with a large
`phip` residual (`combined ~= 2.42`, `phip ~= 2.42`). That run is invalid as a
physics comparison because the restart state and contact bias did not match.

The corrected diagnostic strategy generated one Vela restart CSV per bias from
the corresponding Sentaurus intermediate export, then ran strict single-point
Vela solves with contact-edge diagnostics:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag
```

All eight single-bias points converged:

```text
-12.5, -12.6, -12.7, -12.8, -12.9, -13.0, -13.1, -13.2 V
```

Aggregated Anode current/QF-drop output:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/anode_qf_drop_current_summary.csv
```

Summary:

```text
bias   abs log current err   Vela Anode I      Sentaurus I      hQF drop ratio   Cathode total
-12.5  0.708                 -1.506e-17        -7.688e-17       5.0              +7.114e-17
-12.6  2.155                 -5.457e-19        -7.801e-17       undefined        +8.384e-17
-12.7  2.160                 -5.466e-19        -7.905e-17       undefined        +8.384e-17
-12.8  0.726                 -1.506e-17        -8.003e-17       5.0              +8.335e-17
-12.9  0.734                 -1.495e-17        -8.099e-17       6.0              +8.335e-17
-13.0  0.739                 -1.495e-17        -8.193e-17       6.0              +8.335e-17
-13.1  0.744                 -1.495e-17        -8.288e-17       6.0              +8.335e-17
-13.2  0.774                 +1.410e-17        -8.385e-17       -5.0             +8.335e-17
```

Interpretation:

- The Anode current error is not isolated to `-12.8 V`; it repeats across the
  transition window.
- At `-12.6 V` and `-12.7 V`, Vela's Anode hole current is exactly zero in the
  contact-current diagnostic because the Anode hole QF drop is exactly zero,
  while Sentaurus still reports a finite Anode hole current around
  `-7.8e-17 A`.
- Where the Vela Anode hole QF drop is nonzero, the Sentaurus/Vela endpoint
  hole-QF-drop ratio is `5..6`, matching the current error scale.
- Vela Cathode total/electron current is already close to Sentaurus magnitude,
  so the remaining transition-window mismatch is a contact-side Anode hole
  QF-gradient/numerical-floor issue, not a global current unit or mobility
  amplitude issue.
- The `-13.2 V` TDR `ContactCurrentFlux` differs from `.plt` by `13.3%`, so
  that point should be treated cautiously for TDR-flux parity. The `.plt`
  curve remains the BV target for current comparison.

### Next Tasks After Task 40

1. Implement a read-only Anode hole-QF sensitivity probe:
   - load a converged Vela state and its contact-edge diagnostics;
   - replace only the Anode contact-adjacent `phip` endpoint drop with the
     Sentaurus endpoint drop while keeping `psi`, density, mobility, and
     electron QF fixed;
   - recompute the contact current and report whether it lands on Sentaurus
     `.plt hCurrent`.
2. If the sensitivity probe closes the gap, prototype a Sentaurus-like
   contact-current evaluation or contact-boundary projection policy that
   preserves the tiny Anode hole-QF numerical floor without corrupting the
   already-good Cathode electron branch.
3. If the sensitivity probe does not close the gap, inspect the hole SG
   coefficient path on Anode edges: effective `ni`, hole mobility, Bernoulli
   weights, edge-couple sum, and sign convention.

### Execution Note 2026-06-20: Task 41 Anode Hole-QF Sensitivity Closure

Task 41 implemented and ran the read-only Anode hole-QF sensitivity probe from
Task 40. The probe keeps the converged Vela contact-edge `psi`, `ni`, hole
mobility, edge geometry, electron current, and sign convention fixed, then
recomputes only the Anode hole SG contact-edge current after replacing the
Vela `phip` endpoint drop with the Sentaurus endpoint drop exported by the
contact-current diagnostic.

Code changes:

- `scripts/diagnose_pn2d_bv_anode_hole_qf_sensitivity.py`
  - Reads `contact_edges_filtered.csv` and `edge_summary.csv`.
  - Recomputes the Anode hole current with the combined-exponent SG hole flux
    formula and the Sentaurus contact-edge hole QF drop.
  - Emits per-edge CSV and summary JSON with baseline and sensitivity currents.
- `tests/regression/test_reference_tcad_tools.py`
  - Added a regression that verifies a tiny Sentaurus-like Anode hole QF drop
    changes the hole component and moves the total current toward the
    Sentaurus reference.

Verification:

```text
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_contact_current_extraction_filters_edges_and_classifies_parity tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_contact_current_extraction_aligns_sentaurus_nearest_fields tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_anode_hole_qf_sensitivity_recomputes_edge_current
```

Output:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/anode_hole_qf_sensitivity_summary.csv
```

Summary:

```text
bias   baseline err dex   sensitivity err dex   baseline I       sensitivity I    Sentaurus I
-12.5  0.708              0.021                 -1.506e-17      -7.328e-17      -7.688e-17
-12.6  2.155              0.029                 -5.457e-19      -7.303e-17      -7.801e-17
-12.7  2.160              0.035                 -5.466e-19      -7.289e-17      -7.905e-17
-12.8  0.726              0.038                 -1.506e-17      -7.328e-17      -8.003e-17
-12.9  0.734              0.033                 -1.495e-17      -8.744e-17      -8.099e-17
-13.0  0.739              0.028                 -1.495e-17      -8.744e-17      -8.193e-17
-13.1  0.744              0.023                 -1.495e-17      -8.744e-17      -8.288e-17
-13.2  0.774              0.060                 +1.410e-17      -7.303e-17      -8.385e-17
```

Interpretation:

- The sensitivity probe closes nearly all of the transition-window current
  mismatch. The log-current error falls from `0.708..2.160 dex` to
  `0.021..0.060 dex`.
- The remaining difference is only a few percent to about `15%`, comparable to
  the uncertainty of using endpoint field exports and, at `-13.2 V`, the
  already-observed `.plt` versus TDR `ContactCurrentFlux` spread.
- Therefore the main BV transition-window mismatch is now localized to the
  Anode contact-adjacent hole quasi-Fermi micro-gradient / numerical-floor
  policy. It is not a terminal extraction, current unit, mobility amplitude,
  or bulk SG flux amplitude problem.
- The pathological `-12.6 V` and `-12.7 V` cases are explained by Vela setting
  the Anode hole QF drop exactly to zero while Sentaurus retains a tiny finite
  drop that is sufficient to produce the observed minority-carrier contact
  current.

### Next Tasks After Task 41

1. Add a contact-boundary projection diagnostic that prints, for each ohmic
   contact edge and each Newton state import:
   - raw imported `phip` at the contact node and adjacent node;
   - projected/Dirichlet-applied `phip`;
   - residual row value before and after contact constraints;
   - whether equality was forced by boundary projection, restart rounding, or
     solver convergence.
2. Prototype a Sentaurus-like contact-current evaluation option that preserves
   the tiny contact-adjacent minority-carrier QF floor for current reporting
   only. This should be gated as a diagnostic mode first, not used to change the
   nonlinear solve until its effect on forward/reverse low-bias parity is known.
3. Inspect the ohmic contact Dirichlet implementation and restart importer for
   exact-value clamping of quasi-Fermi potentials. The next likely fix is to
   avoid erasing sub-ulp contact-adjacent minority-carrier gradients when
   reconstructing/reporting terminal current.
4. Keep avalanche and bulk mobility tuning frozen until the contact QF-floor
   branch is resolved; Task 41 shows those are downstream for this mismatch.

### Execution Note 2026-06-20: Task 42 Contact Boundary Projection Source

Task 42 traced the Anode hole QF micro-gradient one layer earlier than the
contact-current sensitivity probe. The goal was to determine whether the tiny
Sentaurus-like Anode hole `phip` drop is lost in:

- the Sentaurus-to-Vela restart CSV;
- the ohmic contact Dirichlet projection before Newton;
- Newton convergence;
- or contact-current CSV reporting.

Code changes:

- `scripts/diagnose_pn2d_bv_contact_boundary_projection.py`
  - Compares `sentaurus_state_restart.csv`, a simulated ohmic contact
    projection, `last_state.csv`, and `contact_edges.csv`.
  - Uses `outward_sign` to classify contact-side and interior-side edge nodes.
  - Reports per-edge `restart`, `projected`, `final`, and reported
    `phip` drops plus a source classification.
- `tests/regression/test_reference_tcad_tools.py`
  - Added a regression proving that a restart contact-node QF offset can be
    erased by projecting the contact node to the exact Vela bias.

Verification:

```text
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_contact_current_extraction_filters_edges_and_classifies_parity tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_contact_current_extraction_aligns_sentaurus_nearest_fields tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_anode_hole_qf_sensitivity_recomputes_edge_current tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_contact_boundary_projection_classifies_erased_qf_drop
```

Output:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/contact_boundary_projection_summary.csv
```

Summary:

```text
bias   restart max |dphip|   projected max |dphip|   final max |dphip|   dominant source
-12.5  8.882e-15           1.776e-15              1.776e-15          final_state_preserved_drop
-12.6  8.882e-15           0.000e+00              0.000e+00          contact_projection_erased_restart_drop
-12.7  8.882e-15           0.000e+00              0.000e+00          contact_projection_erased_restart_drop
-12.8  8.882e-15           1.776e-15              1.776e-15          final_state_preserved_drop
-12.9  1.066e-14           1.776e-15              1.776e-15          final_state_preserved_drop
-13.0  1.066e-14           1.776e-15              1.776e-15          final_state_preserved_drop
-13.1  1.066e-14           1.776e-15              1.776e-15          final_state_preserved_drop
-13.2  8.882e-15           1.776e-15              1.776e-15          final_state_preserved_drop
```

Interpretation:

- The Sentaurus-like Anode hole QF micro-gradient is present in every
  per-bias restart CSV before Newton. All 17 Anode contact edges per bias have
  a nonzero restart drop of about `5..6` double-precision ulps near
  `|V| ~= 13 V`.
- Vela's ohmic contact projection then overwrites the contact-node `phip` with
  the exact configured bias. That projection immediately reduces the restart
  drop to either `0` or one ulp (`1.776e-15 V`), depending on binary
  representation of the bias and adjacent interior value.
- Newton does not materially change the projected contact-edge hole QF drop:
  the final-state max equals the projected max for all eight points.
- Therefore the source of the Task 41 current mismatch is now localized to
  contact-boundary projection / contact-current endpoint policy. It is not
  caused by the restart importer losing Sentaurus precision, and it is not a
  Newton convergence drift after projection.
- The especially bad `-12.6 V` and `-12.7 V` current errors occur because the
  projection erases the restart drop completely at those biases. The
  `-13.2 V` sign anomaly is also consistent with projection changing the
  micro-gradient sign/magnitude relative to the restart value.

### Next Tasks After Task 42

1. Prototype a diagnostic-only contact-current reporting mode that evaluates
   Anode minority-hole contact SG flux using a preserved contact-side QF floor:
   - first use the raw restart contact-edge `phip` endpoint drop to confirm it
     reproduces Task 41's Sentaurus-drop closure;
   - then test a restart-free policy based on a fixed ulp floor (`5..6 ulp`)
     with sign chosen from the pre-projection/restart direction or from the
     terminal current convention.
2. Keep this reporting mode separate from nonlinear residual assembly until
   low-bias forward/reverse parity is rechecked. The residual currently enforces
   a mathematically clean Ohmic Dirichlet boundary; the Sentaurus parity issue
   is in default terminal-current endpoint numerics at the micro-gradient floor.
3. Add a focused comparison table for:
   - baseline Vela current;
   - Sentaurus-drop sensitivity current;
   - restart-drop sensitivity current;
   - ulp-floor reporting current;
   - Sentaurus `.plt` current.
4. If the restart-drop and ulp-floor variants both close the gap, implement the
   smallest guarded Vela option for BV validation, for example
   `contact_current_qf_floor = sentaurus_like_ulp`, and keep the default solver
   equations unchanged.

### Execution Note 2026-06-20: Task 43 Reporting-Only QF Floor Policy Probe

Task 43 implemented the diagnostic-only contact-current reporting probe from
Task 42. The probe does not change nonlinear residual assembly or the final
solution. It recomputes only the Anode minority-hole contact-edge SG current
with several endpoint-drop policies:

- baseline Vela final-state contact current;
- Sentaurus endpoint `hQuasiFermi` drop exported by the TDR field alignment;
- raw restart-state `phip(node1) - phip(node0)` before Vela contact
  projection;
- fixed `5 ulp` floor at the bias value, sign from restart;
- fixed `6 ulp` floor at the bias value, sign from restart.

Code changes:

- `scripts/diagnose_pn2d_bv_contact_qf_floor_reporting.py`
  - Reads `contact_edges_filtered.csv`, `edge_summary.csv`, and optional
    `sentaurus_state_restart.csv`.
  - Recomputes Anode hole SG current for Sentaurus-drop, restart-drop, and
    fixed-ulp floor policies while keeping electron current, `psi`, `ni`,
    mobility, and geometry fixed.
  - Emits per-edge CSV and summary JSON.
- `tests/regression/test_reference_tcad_tools.py`
  - Added a regression that checks restart-drop and fixed-ulp reporting
    policies move a flat Anode hole contact current toward the reference.

Verification:

```text
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_contact_current_extraction_filters_edges_and_classifies_parity tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_contact_current_extraction_aligns_sentaurus_nearest_fields tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_anode_hole_qf_sensitivity_recomputes_edge_current tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_contact_boundary_projection_classifies_erased_qf_drop tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_contact_qf_floor_reporting_compares_restart_and_ulp_policies
```

Output:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/contact_qf_floor_reporting_policy_comparison.csv
```

Summary:

```text
bias   baseline err   sentaurus/restart err   5 ulp err   6 ulp err
-12.5  0.708          0.021                   0.021       0.056
-12.6  2.155          0.029                   0.029       0.051
-12.7  2.160          0.035                   0.035       0.044
-12.8  0.726          0.038                   0.038       0.039
-12.9  0.734          0.033                   0.044       0.033
-13.0  0.739          0.028                   0.049       0.028
-13.1  0.744          0.023                   0.054       0.023
-13.2  0.774          0.060                   0.060       0.020
```

Restart-drop measured in ulps of the bias value:

```text
bias   restart edge-order hQF drop
-12.5  -5 ulp
-12.6  -5 ulp
-12.7  -5 ulp
-12.8  -5 ulp
-12.9  -6 ulp
-13.0  -6 ulp
-13.1  -6 ulp
-13.2  -5 ulp
```

Interpretation:

- The Sentaurus endpoint-drop and raw restart-drop reporting currents are
  identical across the full `-12.5..-13.2 V` window. This proves the
  Sentaurus-like contact current can be reproduced by preserving the
  pre-projection endpoint `phip` micro-gradient for reporting.
- The fixed `5 ulp` policy closes the gap for `-12.5..-12.8 V` and `-13.2 V`
  but underestimates `-12.9..-13.1 V`.
- The fixed `6 ulp` policy closes `-12.9..-13.1 V` and improves `-13.2 V`,
  but overestimates `-12.5..-12.7 V`.
- Therefore a single hard-coded fixed-ulp floor is not the right implementation
  target. The safer next implementation branch is a guarded reporting-only mode
  that carries the pre-contact-projection endpoint QF value (or endpoint drop)
  into terminal-current evaluation.
- If a restart-free option is still desired, it should be a secondary fallback
  with a dynamic ulp count, not a constant `5 ulp` knob.

### Next Tasks After Task 43

1. Implement a guarded C++ diagnostic/reporting option that preserves
   pre-projection contact-edge QF endpoint drops for terminal current only:
   - capture initial/restart contact-edge `phip` before applying ohmic contact
     Dirichlet projection;
   - pass the preserved edge-order drop into `ContactCurrent` only when the
     option is enabled;
   - leave CoupledDD residual/Jacobian boundary conditions unchanged.
2. Add C++ tests around the option using a tiny contact-edge `phip` drop:
   - default current remains the mathematically clean projected result;
   - enabled reporting mode reproduces the preserved-drop current;
   - electron branch and Cathode branch are unchanged when no preserved hole
     floor is provided.
3. Re-run forward/reverse low-bias and BV transition diagnostics with the
   reporting option enabled:
   - require no degradation in `0..5 V` forward and `0..-5 V` reverse windows;
   - require transition-window Anode current error to stay near the Task 43
     `0.02..0.06 dex` range.
4. Only after this reporting-only option is verified should we consider a
   restart-free dynamic-ulp fallback or any residual-level boundary-policy
   change.

### Execution Note 2026-06-20: Task 44 ContactCurrent Override Seam

Task 44 implemented the first C++ layer needed for the Task 43 reporting-only
fix. This step deliberately stopped at `ContactCurrent` and did not yet wire a
new `DCSweep` config option, so the nonlinear residual/Jacobian path remains
unchanged.

Code changes:

- `include/vela/post/ContactCurrent.h`
  - Added `ContactCurrentEdgeOverrides`.
  - Added `compute` / `computeDetailed` overloads accepting per-edge reporting
    overrides.
  - Added a diagnostic flag
    `ContactCurrentEdgeDiagnostic::holeQfDropOverrideApplied`.
- `src/post/ContactCurrent.cpp`
  - Applies `holeQuasiFermiDropByEdge[edgeId]` only to the hole SG reporting
    flux on the matching contact edge.
  - Leaves electron current, mobility evaluation, and default current reporting
    unchanged when no override is supplied.
- `tests/test_newton_solver.cpp`
  - Added a focused `ContactCurrent` regression: a flat final-state `phip`
    gives near-zero hole current by default; supplying a preserved edge-order
    `phip` drop changes only reporting current and marks the edge diagnostic.

Verification:

```text
cmake --build build-release --target test_newton_solver --parallel
build-release/test_newton_solver.exe "ContactCurrent: preserved edge hole QF drop changes reporting only"
```

Focused result:

```text
All tests passed (9 assertions in 1 test case)
```

Python diagnostic regression result:

```text
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_contact_current_extraction_filters_edges_and_classifies_parity tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_contact_current_extraction_aligns_sentaurus_nearest_fields tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_anode_hole_qf_sensitivity_recomputes_edge_current tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_contact_boundary_projection_classifies_erased_qf_drop tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_contact_qf_floor_reporting_compares_restart_and_ulp_policies
```

```text
Ran 5 tests in 1.375s
OK
```

Verification caveat:

```text
build-release/test_newton_solver.exe
```

currently fails two existing Newton convergence tests in this worktree:

```text
NewtonSolver: unit_scaling accepts a physical Gummel warm initial guess
NewtonSolver: ohmic contact BC resists compensated-node polarity flips
```

Those failures are outside the new `ContactCurrent` override path but must be
cleared before claiming full branch health or merging.

Interpretation:

- The C++ reporting seam needed by Task 43 now exists and is covered by a
  focused regression.
- The override is intentionally explicit and edge-local: it cannot affect
  CoupledDD residual assembly unless a caller opts in and passes overrides.
- This is the right insertion point for the Sentaurus-like BV terminal-current
  parity path: `DCSweep` can capture restart/pre-projection contact-edge drops
  and pass them only to terminal-current reporting.

### Next Tasks After Task 44

1. Design and implement the guarded `DCSweep` integration:
   - add a config field for reporting-only preserved contact QF drops;
   - capture contact-edge `phip(node1)-phip(node0)` from the supplied
     initial/restart state before Newton projects contact nodes;
   - pass captured drops into `ContactCurrent::computeDetailed` only for the
     intended contact/current reporting path;
   - keep branch-acceptance and residual/Jacobian behavior explicit.
2. Add `test_dc_sweep` coverage proving:
   - default output current is unchanged;
   - enabling the option changes only terminal current reporting for edges with
     preserved drops;
   - contact-edge diagnostics mark overridden hole QF drops.
3. Before broad validation, resolve or isolate the two currently failing
   `test_newton_solver` convergence cases noted above.
4. Re-run the PN2D per-bias BV transition window with the DCSweep option
   enabled and require the current error to match the Task 43 reporting probe
   (`0.02..0.06 dex`) before testing low-bias forward/reverse windows.

### Execution Note 2026-06-20: Task 45 DCSweep QF Floor Reporting Integration

Task 45 wired the Task 44 `ContactCurrent` override seam into `DCSweep` behind
an explicit diagnostics switch:

```json
"diagnostics": {
  "contact_current_qf_floor": {
    "enabled": true,
    "contacts": ["anode"]
  }
}
```

Implementation:

- `include/vela/simulation/DCSweep.h`
  - Added `ContactCurrentQfFloorDiagnosticsConfig`.
  - Added `SweepDiagnosticsConfig::contactCurrentQfFloor`.
- `src/simulation/DCSweep.cpp`
  - Parses `sweep.diagnostics.contact_current_qf_floor`.
  - Captures initial/restart contact-edge hole quasi-Fermi endpoint drops as
    `phip(node1)-phip(node0)` before the selected solver runs.
  - Passes the captured drops only to `ContactCurrent::computeDetailed` for
    reporting/current extraction.
  - Leaves Newton/Gummel residual assembly, Jacobian assembly, and branch
    acceptance on the default physical solution path.
  - Adds `hole_qf_drop_override_applied` to contact-edge diagnostics CSV.
- `tests/test_dc_sweep.cpp`
  - Added focused DCSweep integration coverage for the guarded option.
  - The regression writes default and enabled sweeps, verifies the default
    contact-edge report has no applied override, verifies the enabled report
    marks overridden hole QF drops, and checks the sweep point hole current is
    consistent with the detailed edge-current sum.

Verification:

```text
cmake --build build-release --target test_dc_sweep --parallel
build-release/test_dc_sweep.exe "DCSweep: contact current QF floor reporting uses initial edge drops only when enabled"
```

```text
All tests passed (22 assertions in 1 test case)
```

```text
build-release/test_dc_sweep.exe "[dc_sweep][diagnostics]"
```

```text
All tests passed (398 assertions in 9 test cases)
```

```text
cmake --build build-release --target test_newton_solver --parallel
build-release/test_newton_solver.exe "ContactCurrent: preserved edge hole QF drop changes reporting only"
```

```text
All tests passed (9 assertions in 1 test case)
```

```text
git diff --check -- include/vela/simulation/DCSweep.h src/simulation/DCSweep.cpp tests/test_dc_sweep.cpp include/vela/post/ContactCurrent.h src/post/ContactCurrent.cpp tests/test_newton_solver.cpp
```

```text
No whitespace errors.
```

Verification caveats:

- `python -m pytest tests/regression/test_reference_tcad_tools.py -q` could not
  run in the current UCRT64 shell because `D:\msys64\ucrt64\bin\python.exe`
  does not have `pytest` installed.
- The broader `test_newton_solver.exe` binary still has the two pre-existing
  convergence failures recorded in Task 44. They remain outside the
  reporting-only QF floor path but still gate full branch health.

Interpretation:

- Vela now has a guarded reporting-only mechanism matching the Task 43 probe:
  contact-edge hole QF micro-drops can be preserved from restart/pre-projection
  state for terminal-current extraction without changing the nonlinear solve.
- The option is intentionally diagnostic/compatibility-scoped. It should be
  enabled only for Sentaurus-parity BV investigations until the real PN2D
  transition-window validation proves it closes the `-12.5..-13.2 V` current
  gap without damaging low-bias windows.

### Next Tasks After Task 45

1. Run the real PN2D BV transition window with:
   - default Vela reporting;
   - `diagnostics.contact_current_qf_floor.enabled=true` for `anode`;
   - the same Sentaurus restart/reference alignment used by Tasks 41-43.
2. Require the enabled run to reproduce the Task 43 reporting probe target:
   transition-window Anode current error near `0.02..0.06 dex`, instead of the
   default `0.708..2.160 dex` mismatch.
3. Re-run the `0..5 V` forward and `0..-5 V` reverse windows with the option
   enabled to check that low-bias IV, potential, field, electron density, and
   hole density comparisons do not regress.
4. If the real sweep does not close the transition-window error, inspect the
   captured `hole_qf_drop_override_applied` contact-edge CSV rows first:
   missing flags imply restart/initial-state capture is not reaching DCSweep;
   present flags with wrong current imply edge orientation or terminal-current
   summation mismatch.

### Execution Note 2026-06-20: Task 46 Real Transition-Window QF Floor Validation

Task 46 executed the first real PN2D validation after the DCSweep integration:
reuse the existing per-bias Sentaurus-seeded transition-window configs from
Task 40, enable `diagnostics.contact_current_qf_floor` for `Anode`, and run the
eight strict single-point solves:

```text
-12.5, -12.6, -12.7, -12.8, -12.9, -13.0, -13.1, -13.2 V
```

Generated configs and outputs:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/m12p5_qf_floor_enabled
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/m12p6_qf_floor_enabled
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/m12p7_qf_floor_enabled
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/m12p8_qf_floor_enabled
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/m12p9_qf_floor_enabled
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/m13p0_qf_floor_enabled
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/m13p1_qf_floor_enabled
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/m13p2_qf_floor_enabled
```

Aggregate report:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/qf_floor_enabled_window_summary.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/qf_floor_enabled_window_summary.json
```

Verification commands:

```text
cmake --build build-release --target vela_example_runner --parallel
build-release/vela_example_runner.exe --config <per-bias qf_floor_enabled simulation.json>
```

All eight points converged with one point per run.

Result summary:

```text
bias    enabled err dex   restart-probe err dex   baseline err dex   override edges
-12.5   0.021             0.021                   0.708              17/17
-12.6   0.029             0.029                   2.155              17/17
-12.7   0.035             0.035                   2.160              17/17
-12.8   0.038             0.038                   0.726              17/17
-12.9   0.033             0.033                   0.734              17/17
-13.0   0.028             0.028                   0.739              17/17
-13.1   0.023             0.023                   0.744              17/17
-13.2   0.060             0.060                   0.774              17/17
```

The enabled DCSweep currents match the Task 43 restart-drop predictions to
roundoff. Example at `-12.7 V`:

```text
baseline Vela current:       -5.465540102009639e-19 A/um
enabled DCSweep current:     -7.288777764227919e-17 A/um
Task 43 restart prediction:  -7.288777764227923e-17 A/um
Sentaurus .plt target:       -7.904911015264140e-17 A
```

Interpretation:

- The production DCSweep reporting option now reproduces the read-only Task 43
  restart-drop policy on the real transition-window artifacts.
- The transition-window Anode current mismatch is reduced from
  `0.708..2.160 dex` to `0.021..0.060 dex`, matching the acceptance target from
  Task 45.
- Every Anode contact edge in the window reports
  `hole_qf_drop_override_applied=1`; Cathode rows remain unmodified by the
  Anode-only option.
- This confirms the dominant BV transition-window current mismatch source is
  the Anode minority-hole contact QF micro-gradient erased by contact
  projection, not bulk mobility, SG flux amplitude, current units, or terminal
  summation.

### Next Tasks After Task 46

1. Run low-bias risk checks with the option enabled:
   - `0..5 V` forward window;
   - `0..-5 V` reverse window;
   - compare IV plus potential/electric-field/electron-density/hole-density
     summaries against the existing `forward_reverse_windows` report.
2. Decide the production/default policy:
   - keep the option disabled by default and document it as
     Sentaurus-parity/reporting compatibility;
   - or enable it only in imported Sentaurus restart-validation configs.
3. Add schema documentation for
   `sweep.diagnostics.contact_current_qf_floor.enabled` and `contacts`.
4. After low-bias risk checks pass, consider a checked regression that runs one
   reduced per-bias PN2D restart point, probably `-12.7 V`, and asserts the
   enabled current stays near the Task 46 value without changing the nonlinear
   solve result.

### Execution Note 2026-06-20: Task 47 Low-Bias Risk Check And Capture-Source Tightening

Task 47 ran the low-bias risk checks requested after Task 46 and found one
important policy issue before finalizing the option:

- If `contact_current_qf_floor` captured from every DCSweep initial guess, then
  an ordinary continuation sweep could use the previous Vela solution as the
  reporting floor source.
- In the reverse low-bias `0..-5 V` sweep this changed the reported `-0.5 V`
  current by about `12%` relative to a current-code baseline.
- That behavior is not the intended Sentaurus restart-parity semantics. The
  intended source is an external restart/import state, not Vela's own previous
  continuation point.

Implementation tightening:

- `src/simulation/DCSweep.cpp`
  - `solvePoint` now receives an explicit
    `allowContactCurrentQfFloorCapture` flag.
  - The flag is true only for the first point initialized directly from
    `sweep.initial_state_file`.
  - Predictor states and ordinary continuation from `previousSolution` pass
    `false`.
  - This preserves the per-bias Sentaurus restart validation behavior while
    preventing low-bias continuation sweeps from changing current reporting.
- `tests/test_dc_sweep.cpp`
  - Added coverage that an ordinary continuation sweep with the option enabled
    does not mark contact-edge QF floor overrides.
- `docs/config_schema.md`
  - Documented `sweep.diagnostics.contact_current_qf_floor.enabled` and
    `contacts`.
  - Documented that the option is reporting-only, external-restart sourced, and
    does not change residual/Jacobian/solution/VTK or ordinary continuation.

Focused C++ verification:

```text
cmake --build build-release --target test_dc_sweep --parallel
build-release/test_dc_sweep.exe "DCSweep: contact current QF floor reporting uses initial edge drops only when enabled"
build-release/test_dc_sweep.exe "DCSweep: contact current QF floor reporting ignores continuation states"
```

Results:

```text
All tests passed (22 assertions in 1 test case)
All tests passed (8 assertions in 1 test case)
```

Transition-window recheck after tightening:

```text
bias    enabled err dex   restart-probe err dex   baseline err dex   override edges
-12.5   0.021             0.021                   0.708              17/17
-12.6   0.029             0.029                   2.155              17/17
-12.7   0.035             0.035                   2.160              17/17
-12.8   0.038             0.038                   0.726              17/17
-12.9   0.033             0.033                   0.734              17/17
-13.0   0.028             0.028                   0.739              17/17
-13.1   0.023             0.023                   0.744              17/17
-13.2   0.060             0.060                   0.774              17/17
```

Low-bias risk checks after tightening:

Reverse `0..-5 V`, current-code baseline versus qf-floor enabled:

```text
bias   baseline A/um      enabled A/um       relative delta
 0.0   -3.547878e-19      -3.547878e-19      0
-0.5   -5.662533e-18      -5.662533e-18      0
-2.0   -1.309965e-17      -1.309965e-17      0
-5.0   -2.956164e-17      -2.956164e-17      0
```

Forward `0..5 V`, current-code baseline versus qf-floor enabled:

```text
bias   baseline A/um      enabled A/um       relative delta
 0.0   -5.658605e-19      -5.658605e-19      0
 0.5   -1.554630e-10      -1.554630e-10      0
 1.0   -1.261239e-04      -1.261239e-04      0
 2.0   -2.399073e-03      -2.399073e-03      0
 5.0   -1.004373e-02      -1.004373e-02      0
```

Generated low-bias diff reports:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_low_bias_qf_floor_enabled/iv_delta_vs_current_baseline.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_forward_0_to_5_qf_floor_enabled/iv_delta_vs_current_baseline.csv
```

Interpretation:

- The option is now scoped correctly: it closes the Sentaurus-seeded BV
  transition-window current gap, but has zero IV effect on ordinary forward and
  reverse low-bias continuation sweeps.
- Because the option does not alter the nonlinear solution path, low-bias
  potential, electric field, electron density, and hole density are unchanged
  for ordinary sweeps. Existing low-bias field mismatch investigations remain
  valid and separate from the contact-current reporting parity fix.

### Next Tasks After Task 47

1. Add a compact regression or scripted verifier for one real PN2D restart
   point, preferably `-12.7 V`, so future changes cannot silently lose the
   restart-sourced contact-current parity behavior.
2. Decide whether checked-in reference configs should expose
   `contact_current_qf_floor` only in Sentaurus restart-validation decks, not in
   normal production sweeps.
3. Continue the remaining BV alignment work on the separate low-bias reverse
   density/electrostatics mismatch, which is not fixed by the reporting-only
   contact QF floor option.

### Execution Note 2026-06-20: Task 48 Compact QF Floor Restart-Point Verifier

Task 48 added a narrow verifier for the key real PN2D restart point discovered
in Tasks 43-47. The intent is to preserve the known Sentaurus-restart
contact-current parity behavior without requiring a full transition-window
rerun after every low-level transport or current-extraction edit.

New tool:

```text
scripts/verify_pn2d_bv_qf_floor_restart_point.py
```

Verifier checks:

- the requested `iv.csv` bias row exists and is converged;
- the enabled Vela `current_total_A_per_um` is within a configurable log10
  error threshold versus the Sentaurus reference current from the restart
  diagnostic summary;
- the enabled Vela current matches the standalone restart-drop current
  prediction within a tight absolute tolerance;
- the requested contact has qf-floor override flags in `contact_edges.csv`;
- the expected number of override edges is present;
- optional guard: no other contact receives override flags.

Regression coverage:

```text
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_qf_floor_restart_point_verifier_checks_enabled_current
```

Result:

```text
OK
```

Real artifact verification command:

```text
python scripts/verify_pn2d_bv_qf_floor_restart_point.py \
  --iv-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/m12p7_qf_floor_enabled/iv.csv \
  --contact-edge-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/m12p7_qf_floor_enabled/contact_edges.csv \
  --restart-summary build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/m12p7/contact_qf_floor_reporting/contact_qf_floor_reporting_summary.json \
  --bias -12.7 \
  --contact Anode \
  --expected-override-edges 17 \
  --out-json build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/m12p7_qf_floor_enabled/qf_floor_restart_point_verification.json
```

Result:

```text
qf-floor restart point verified: bias=-12.7 V contact=Anode log10_error=0.0352423 override_edges=17
```

Generated summary:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/per_bias_sentaurus_seed_contact_diag/m12p7_qf_floor_enabled/qf_floor_restart_point_verification.json
```

Key values:

```text
current_total_A_per_um:                 -7.288777764227919e-17
reference_current_A:                    -7.904911015264140e-17
restart_drop_total_current_A_per_um:    -7.288777764227923e-17
abs_log10_error_vs_reference:            0.03524227675351099
abs_delta_vs_restart_current_A_per_um:   3.697785493223493e-32
override_edge_count:                     17
other_contact_override_edge_count:       0
```

Interpretation:

- The qf-floor parity fix now has a one-command guard at the representative
  `-12.7 V` transition-window point.
- This verifier protects the reporting/current-extraction compatibility fix,
  but it deliberately does not claim to validate the bulk electrostatic or
  carrier-density branch.
- The remaining BV mismatch work should now move back to the separate
  low-bias and pre-breakdown branch differences: potential shape, electric
  field, carrier density, mobility, and continuity residual/Jacobian behavior.

### Next Tasks After Task 48

1. Promote the verifier into a lightweight regression target or documented
   developer command, ideally without committing large generated PN2D outputs.
2. Run the verifier after every change touching:
   - `ContactCurrent`;
   - SG flux evaluation;
   - contact boundary projection;
   - `DCSweep` restart/initial-state handling;
   - terminal-current extraction.
3. Resume BV physics alignment on the non-reporting path:
   - compare `0..-5 V` reverse low-bias fields and carrier densities against
     Sentaurus with the qf-floor option disabled;
   - localize whether the remaining IV error starts from mobility/SG flux,
     density branch selection, or terminal-current extraction;
   - use no-impact or frozen-carrier residual probes before adding any new
     physical model knobs.

### Execution Note 2026-06-20: Task 49 Reverse Low-Bias Signed Partition Recheck

Task 49 rechecked the reverse low-bias state using the current worktree and
current generated artifacts, because the older
`forward_reverse_windows/summary.md` snapshot still reported a stale
`30..40%` reverse-current deficit and `~0.5 dex` density error.

New diagnostic:

```text
scripts/diagnose_pn2d_bv_signed_field_partitions.py
```

Purpose:

- reuse the existing PN2D multibias field-compare mesh, mask, VTK, and nearest
  mapping helpers;
- report signed Vela/Sentaurus differences by spatial partition:
  `all`, `centerline`, `p_contact`, `n_contact`, `junction`, `bulk`;
- emit both per-partition summary rows and top signed-error node rows;
- distinguish "Vela higher" from "Vela lower", which the older absolute p95
  field summary could not do.

Regression coverage:

```text
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_signed_partition_diagnostic_reports_direction
```

Result:

```text
OK
```

Real-data inputs:

```text
Sentaurus -0.5 V: reference_tcad/pn2d_sentaurus2018/source/pn2d_bv_multibias_0005_des.tdr
Sentaurus -2.0 V: reference_tcad/pn2d_sentaurus2018/source/pn2d_bv_multibias_0020_des.tdr
Sentaurus -5.0 V: reference_tcad/pn2d_sentaurus2018/source/pn2d_bv_multibias_0050_des.tdr

Vela -0.5 V: build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_low_bias/vtk/official_split_low_bias_0010_-0.5V.vtk
Vela -2.0 V: build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_low_bias/vtk/official_split_low_bias_0040_-2V.vtk
Vela -5.0 V: build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_low_bias/vtk/official_split_low_bias_0100_-5V.vtk
```

Generated reports:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/reverse_low_bias_signed_partitions/m0p5/signed_field_partition_summary.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/reverse_low_bias_signed_partitions/m2p0/signed_field_partition_summary.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/reverse_low_bias_signed_partitions/m5p0/signed_field_partition_summary.csv
```

Current reverse IV cross-check, from
`official_split_low_bias_compare/curve_compare.csv`:

```text
bias    Sentaurus A       Vela A/um          abs log10 error
-0.5    -5.3460336e-18    -5.6685190e-18    0.02544
-2.0    -1.4176130e-17    -1.5502551e-17    0.03885
-5.0    -2.8426516e-17    -3.0554590e-17    0.03135
```

This contradicts the stale `forward_reverse_windows/summary.md` statement that
the current reverse window is still `30..40%` low. In the current artifacts,
reverse low-bias IV is already within about `0.025..0.039 dex`.

Signed partition highlights:

```text
bias   e density p95 dex   h density p95 dex   h density junction median dex   h mobility junction median rel
-0.5   0.080               0.133               -0.002                         +0.031
-2.0   0.061               0.126               -0.099                         +0.209
-5.0   0.055               0.131               -0.110                         +0.221
```

Additional observations:

- Potential agreement is now tight in the low-bias reverse window:
  `5.7e-5 V RMS` at `-0.5 V`, `1.6e-3 V RMS` at `-2 V`, and
  `2.3e-3 V RMS` at `-5 V`.
- Electron density is close in signed median and p95: the all-region median
  Vela/Sentaurus log ratio remains near zero and p95 stays below `0.081 dex`.
- Hole density is the remaining carrier-density outlier, especially in the
  junction partition at `-2..-5 V`, where Vela is lower by about `0.10 dex`
  median.
- Junction hole mobility is consistently higher in Vela than Sentaurus by
  about `0.21..0.22` relative median at `-2..-5 V`.
- Electric-field relative p95 is still large, but the current and density
  errors are already small enough that global E-field p95 alone is a poor BV
  decision metric.

Interpretation:

- The current worktree no longer supports the earlier hypothesis of a broad
  reverse low-bias density-branch failure.
- The low-bias residual difference is now localized: junction hole density,
  hole mobility, and electric-field/gradient metrics deserve priority.
- Terminal-current extraction is unlikely to be the dominant low-bias source:
  current error is already below `0.04 dex` at `-0.5`, `-2`, and `-5 V`.
- The next BV parity risk is whether the small low-bias junction mobility/field
  differences amplify into the high-bias pre-breakdown branch. That should be
  tested before adding new physical knobs.

### Next Tasks After Task 49

1. Refresh or supersede `forward_reverse_windows/summary.md` so it no longer
   reports stale `30..40%` reverse-current deficits.
2. Extend the signed partition diagnostic to the transition/pre-breakdown
   window, starting at `-10 V`, `-12.0 V`, and `-13.2 V`, using the same
   partitions and signed density/mobility/field summaries.
3. Add a focused mobility/field driver comparison around the junction top-error
   nodes:
   - Vela and Sentaurus `eMobility`, `hMobility`;
   - quasi-Fermi gradients;
   - density-gradient SG flux proxies;
   - local current density components.
4. If high-bias signed partitions show the same hole-mobility/junction-density
   trend growing with voltage, inspect the C++ mobility high-field driving
   force and SG edge interpolation first.
5. If high-bias partitions instead show a new density branch shift, return to
   the Newton residual/block-step probes before changing mobility formulas.

### Execution Note 2026-06-20: Task 50 High-Bias Signed Partition Branch Split

Task 50 extended the Task 49 signed partition diagnostic into the
transition/pre-breakdown window on the same continuous Vela branch:

```text
Vela branch:
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_control_from_m10/branch_guard_0p02

Sentaurus -10.0 V: reference_tcad/pn2d_sentaurus2018/source/pn2d_bv_multibias_0100_des.tdr
Sentaurus -12.0 V: reference_tcad/pn2d_sentaurus2018/source/pn2d_bv_multibias_0120_des.tdr
Sentaurus -13.2 V: reference_tcad/pn2d_sentaurus2018/source/pn2d_bv_multibias_0132_des.tdr
```

Generated reports:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/high_bias_signed_partitions/m10p0/signed_field_partition_summary.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/high_bias_signed_partitions/m12p0/signed_field_partition_summary.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/high_bias_signed_partitions/m13p2/signed_field_partition_summary.csv
```

Current comparison, using Vela `branch_guard_0p02/iv.csv` and Sentaurus TDR
`ContactCurrentFlux_region2.csv`:

```text
bias    Vela A/um        Sentaurus A       abs log10 error   |Vela/Sentaurus|
-10.0   -6.066045e-17    -4.864012e-17    0.0959            1.25
-12.0   -6.064936e-17    -7.268505e-17    0.0786            0.83
-13.2   -5.562543e-14    -7.269472e-17    2.8838            765
```

Signed partition highlights:

```text
bias    e density all median dex   h density all median dex   h density junction median dex   h mobility junction median rel
-10.0   -0.008                     -0.114                     -0.159                         +0.213
-12.0   -0.015                     -0.133                     -0.201                         +0.210
-13.2   +2.930                     +1.663                     +2.177                         +0.228
```

Other field/density observations:

- At `-10 V` and `-12 V`, the branch is still broadly Sentaurus-like:
  electron-density p95 is `0.072..0.105 dex`; hole-density p95 is
  `0.175..0.219 dex`; potential p95 is about `0.0065..0.0069 V`.
- The low-bias trend from Task 49 continues through `-12 V`: Vela hole density
  is lower than Sentaurus around the junction, while Vela junction hole
  mobility is about `+0.21` relative median higher than Sentaurus.
- At `-13.2 V`, this pattern changes qualitatively. Vela is no longer just a
  small-mobility/field variant of the Sentaurus branch; it has selected a
  high-density branch:
  - electron-density all-region median is `+2.93 dex`;
  - hole-density all-region median is `+1.66 dex`;
  - junction electron-density median is `+3.10 dex`;
  - junction hole-density median is `+2.18 dex`;
  - top electron-density errors occur near `x ~= 0.75 um` with Vela
    `~5e6 cm^-3` versus Sentaurus `~61 cm^-3`;
  - top hole-density errors occur near `x ~= 0.90..0.94 um` with Vela
    `~2e6 cm^-3` versus Sentaurus `~1e4 cm^-3`.
- The `-13.2 V` electric-field relative metric becomes numerically unstable in
  low-field bulk regions, but the density and current evidence are already
  sufficient to classify the point as a different carrier branch.

Branch-transition IV tail from `branch_guard_0p02/iv.csv`:

```text
bias       Vela current A/um
-12.4000   -6.069078e-17
-12.5000   -6.067467e-17
-12.5507   -6.190283e-17
-12.6000   -7.284957e-17
-12.6519   -1.148209e-16
-12.7000   -1.590840e-16
-12.7498   -8.602542e-16
-12.7873   -4.283388e-15
-12.8218   -1.455574e-14
-12.8645   -3.444874e-14
-12.9011   -4.831333e-14
-13.2000   -5.562543e-14
```

Interpretation:

- The remaining high-bias BV mismatch is now classified as a branch-selection
  problem in the continuous Vela solve between roughly `-12.55 V` and
  `-12.90 V`.
- Low-bias and pre-transition mobility/field differences are real, but they do
  not by themselves explain the `-13.2 V` current error; the decisive symptom is
  the abrupt carrier-density branch jump.
- The qf-floor reporting fix from Tasks 45-48 closes Sentaurus-seeded
  transition-window terminal-current extraction, but it does not change this
  continuous-sweep branch selection.

### Next Tasks After Task 50

1. Run a focused branch-transition diagnostic over `-12.5..-12.9 V`:
   - signed partitions at `-12.55`, `-12.65`, `-12.75`, `-12.85`;
   - current and density jump metrics between adjacent accepted states;
   - top signed-density nodes and their quasi-Fermi/electric-field values.
2. Use Newton residual/block-step probes on paired pre/post jump states:
   - last low-density state near `-12.55..-12.60 V`;
   - first high-density state near `-12.75..-12.85 V`;
   - Sentaurus state at matching bias.
3. Test one-variable hypotheses before changing physics:
   - no-impact transition sweep: does the branch jump still occur?
   - frozen-carrier Poisson reconstruction: does electrostatics trigger the
     jump, or only continuity feedback?
   - mobility/field-driver swap: does the `+0.21` junction hole-mobility
     offset move the jump onset?
4. Only if the residual probes show the nonlinear solve is unstable around the
   same state, evaluate Bank-Rose style damping/continuation controls. Do not
   tune damping as a substitute for explaining the density branch.

### Execution Note 2026-06-20: Task 51 Branch-Transition Jump Report

Task 51 made the transition-window check reproducible by adding:

```text
scripts/diagnose_pn2d_bv_branch_transition_jumps.py
tests/regression/test_reference_tcad_tools.py::test_pn2d_bv_branch_transition_jump_report_finds_density_step
```

The script reads existing signed partition reports plus adjacent Vela VTK
states, then writes:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/branch_transition_jump_report/branch_transition_summary.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/branch_transition_jump_report/branch_transition_adjacent_jumps.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/branch_transition_jump_report/branch_transition_summary.json
```

Input points:

```text
-12.5 V: m12p5, branch_guard_0p02_0275_-12.5V.vtk
-12.6 V: m12p6, branch_guard_0p02_0297_-12.6V.vtk
-12.7 V: m12p7, branch_guard_0p02_0317_-12.7V.vtk
-12.8 V: m12p8, branch_guard_0p02_0349_-12.8V.vtk
-12.9 V: m12p9, branch_guard_0p02_0361_-12.9V.vtk
```

Same-bias all-region summary:

```text
bias    current err dex   e-density median dex   e-density p95 dex   h-density median dex   h-density p95 dex
-12.5   -0.0786           -0.018                 0.120               -0.143                 0.236
-12.6   +0.0010           +0.025                 0.384               -0.137                 0.209
-12.7   +0.2736           +0.405                 0.968               -0.044                 0.110
-12.8   +1.9266           +2.070                 2.582               +0.823                 1.399
-12.9   +2.7451           +2.881                 3.399               +1.605                 2.226
```

Adjacent Vela-state all-region jumps:

```text
step           e-density median/p95 dex     e top node (x,y), jump        h-density median/p95 dex     h top node (x,y), jump
-12.5->-12.6   +0.118 / 0.387               219 (0.75,0.1875), +1.300   +0.009 / 0.035               225 (0.78125,0), -0.213
-12.6->-12.7   +0.408 / 0.629               951 (0.75,0.5), +0.896      +0.043 / 0.172               225 (0.78125,0), -0.215
-12.7->-12.8   +1.698 / 1.790               1236 (1.1875,0.03125), +1.962  +0.923 / 1.445           228 (0.8125,0), +1.611
-12.8->-12.9   +0.803 / 0.840               1236 (1.1875,0.03125), +0.946  +0.789 / 0.837           228 (0.8125,0), +0.910
```

Interpretation:

- The transition is not a gradual low-bias current-extraction error. Vela
  remains close at `-12.6 V`, is already drifting at `-12.7 V`, and has clearly
  selected the high-density branch by `-12.8 V`.
- The decisive carrier jump is `-12.7 -> -12.8 V`: electron density jumps by
  `+1.70 dex` median over the whole device and `+1.73 dex` around the junction;
  hole density jumps by `+0.92 dex` whole-device median and `+1.39 dex` around
  the junction.
- The top density-jump nodes are near the junction/contact transition:
  electron node `1236` at `(1.1875 um, 0.03125 um)` and hole node `228` at
  `(0.8125 um, 0 um)`.
- Potential changes at these top nodes are modest compared with the carrier
  jump (`~0.018 V` for the electron top node, `~0.082 V` for the hole top node
  in the decisive step), so the next investigation should target continuity
  residuals, quasi-Fermi gradients, SG fluxes, and impact-generation feedback.

### Next Tasks After Task 51

1. Add a paired-state residual/block-step probe for `-12.7 V` and `-12.8 V`:
   evaluate Poisson, electron continuity, and hole continuity residual norms
   on both Vela states, plus cross-evaluate each state at the other bias.
2. Dump local edge diagnostics around electron node `1236` and hole node `228`
   for the `-12.7 -> -12.8 V` step:
   - quasi-Fermi potentials and gradients;
   - SG electron/hole flux terms;
   - mobility and high-field driving force;
   - avalanche generation/current-density contribution.
3. Run the same transition-window report for no-impact and electric-field
   impact-drive branches if those VTKs are available. If not available, run the
   shortest sweep that captures `-12.6`, `-12.7`, and `-12.8 V`.
4. Use damping/Bank-Rose controls only after the residual/block-step report
   shows solver-step instability rather than a deterministic physics-feedback
   branch.

### Execution Note 2026-06-20: Task 52 Residual, Local Flux, and No-Impact Branch Probe

Task 52 executed the first three Task 51 follow-ups around the decisive
`-12.7 -> -12.8 V` branch jump.

Small tool compatibility fixes were added so transition exports named
`sentaurus_m12p7v` are accepted consistently:

```text
scripts/diagnose_pn2d_bv_newton_residual_states.py
scripts/diagnose_pn2d_bv_continuity_feedback.py
tests/regression/test_reference_tcad_tools.py::test_pn2d_bv_newton_residual_state_finds_m_token_sentaurus_exports
tests/regression/test_reference_tcad_tools.py::test_pn2d_bv_continuity_feedback_finds_m_token_sentaurus_exports
```

Residual probe:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/branch_transition_residual_states_m127_m128/newton_residual_state_summary.json
build-release/reference_tcad/pn2d_sentaurus2018/reports/branch_transition_residual_states_m127_m128/newton_residual_state_nodes.csv
```

States evaluated at `-12.7 V` and `-12.8 V`:

```text
vela
sentaurus
hybrid_vpsi_sqf        # Vela psi + Sentaurus quasi-Fermi
hybrid_spsi_vqf        # Sentaurus psi + Vela quasi-Fermi
hybrid_spsi_shift_vqf  # Sentaurus psi shifted to selected-node Vela mean + Vela quasi-Fermi
```

Block residual highlights:

```text
state/source                 bias    psi block     phin block      phip block      combined
vela                         -12.7   0.222         3.11e-13       3.84e-13       0.222
sentaurus                    -12.7   2.722         4.94e-12       2.56e-12       2.722
hybrid_vpsi_sqf              -12.7   0.222         4.95e-12       2.57e-12       0.222
hybrid_spsi_vqf              -12.7   2.722         3.44e-13       3.70e-13       2.722
vela                         -12.8   0.227         2.20e-11       4.34e-10       0.227
sentaurus                    -12.8   2.722         6.22e-12       2.92e-12       2.722
hybrid_vpsi_sqf              -12.8   0.227         6.25e-12       2.93e-12       0.227
hybrid_spsi_vqf              -12.8   2.722         2.38e-11       4.12e-10       2.722
```

Interpretation from residual/hybrid probe:

- The `newton_residual_probe` combined norm is dominated by the Poisson block,
  not continuity blocks. This means the raw block norm is not yet evidence that
  Bank-Rose damping is the primary missing mechanism.
- The Poisson block follows the electrostatic potential source:
  Vela-psi states are around `0.22`, Sentaurus-psi states around `2.72`.
- At `-12.8 V`, Vela quasi-Fermi states show larger continuity residuals than
  at `-12.7 V` (`phin ~2e-11`, `phip ~4e-10`), but the same trend appears when
  Vela QF is combined with Sentaurus psi. This points to QF/continuity feedback
  as the branch symptom rather than a pure Poisson trigger.

Local edge feedback probes:

```text
electron jump edge:
build-release/reference_tcad/pn2d_sentaurus2018/reports/branch_transition_continuity_feedback_edge3646

hole/contact jump edge:
build-release/reference_tcad/pn2d_sentaurus2018/reports/branch_transition_continuity_feedback_edge709
```

Focus edges:

```text
edge 3646: nodes 1236-1237, near electron top-jump node 1236
edge 709:  nodes 225-228, near hole/contact top-jump node 228
```

Key local evidence:

```text
edge 3646, focus electron side
bias    Vela/Sentaurus n avg dex   Vela electron flux abs     Sentaurus electron flux abs   Vela avalanche node integral at 1236
-12.7   +0.34                       8.80e11                   9.24e14                       4.68e4 s^-1
-12.8   +2.24                       6.18e13                   9.35e14                       1.98e6 s^-1

edge 709, contact-side feedback
bias    Vela/Sentaurus n avg dex   Vela electron flux abs     Sentaurus electron flux abs   Vela/Sentaurus generation dex
-12.7   +1.19                       5.62e14                   4.52e13                       +0.29
-12.8   +3.02                       3.93e16                   4.58e13                       +2.16
```

At focus endpoints:

```text
node 1236:
-12.7: electron density +0.32 dex, electron residual estimate -1.43e7 s^-1
-12.8: electron density +2.23 dex, electron residual estimate -4.25e8 s^-1

node 228:
-12.7: electron density +1.09 dex, hole density -0.13 dex
-12.8: electron density +2.91 dex, hole density +1.44 dex
```

No-impact transition contrast:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/noimpact_transition_signed_partitions
build-release/reference_tcad/pn2d_sentaurus2018/reports/noimpact_transition_jump_report
```

No-impact adjacent all-region jumps:

```text
step           e-density median/p95 dex     h-density median/p95 dex
-12.7->-12.8   +1.691 / 1.791              +0.0001 / 0.0077
-12.8->-12.9   +0.800 / 0.837              +0.0000 / 0.0084
```

Interpretation from no-impact contrast:

- Removing impact ionization does **not** remove the electron-density branch
  jump. The electron jump remains essentially the same as the impact-enabled
  branch (`+1.69 dex` vs `+1.70 dex` median for `-12.7 -> -12.8 V`).
- Removing impact ionization **does** remove the large hole-density branch jump
  (`+0.0001 dex` no-impact vs `+0.92 dex` impact-enabled median).
- Therefore the likely order of events is:
  1. electron continuity / quasi-Fermi / SG feedback selects the high-electron
     branch;
  2. impact ionization amplifies hole density and terminal current after the
     electron branch has moved;
  3. damping alone should not be tuned until this electron-branch trigger is
     explained.

### Next Tasks After Task 52

1. Reconstruct the electron SG/QF terms on edge `709` and neighboring edge
   `3646` using both Vela and Sentaurus quasi-Fermi drops at `-12.7` and
   `-12.8 V`, and compare:
   - variable-`ni` QF SG form;
   - density-form SG flux;
   - electric-field drive versus quasi-Fermi-gradient drive.
2. Run a no-impact minimal sweep with an electron-branch guard or predictor
   limiter around `-12.65..-12.85 V`. The acceptance criterion should be
   electron-density jump, not terminal current.
3. Compare Vela's accepted step path with Sentaurus log/step behavior in the
   same window. If Sentaurus uses Bank-Rose damping there, record the step
   reductions and residual trend, but treat damping as a continuation-control
   comparison rather than root-cause proof.
4. If SG/QF reconstruction shows Vela's electron flux is too large only in the
   QF variable-`ni` form, inspect the C++ implementation of
   `sgElectronContinuityFluxFromQuasiFermiVariableNi` and its edge `ni_eff`
   interpolation/sign conventions before changing solver damping.

### Execution Note 2026-06-20: Task 53 SG Flux-Form Decomposition

Task 53 added a focused edge flux-form diagnostic:

```text
scripts/diagnose_pn2d_bv_sg_flux_forms.py
tests/regression/test_reference_tcad_tools.py::test_pn2d_bv_sg_flux_form_diagnostic_compares_state_sources
```

Reports generated:

```text
impact-enabled branch:
build-release/reference_tcad/pn2d_sentaurus2018/reports/branch_transition_sg_flux_forms_edges709_3646/sg_flux_form_edges.csv

no-impact branch:
build-release/reference_tcad/pn2d_sentaurus2018/reports/noimpact_transition_sg_flux_forms_edges709_3646/sg_flux_form_edges.csv
```

The report compares, for edges `709` and `3646` at `-12.7`, `-12.8`, and
`-12.9 V`:

- density-form electron SG flux;
- variable-`ni` quasi-Fermi SG flux using Vela's OldSlotboom/doping `ni_model`;
- variable-`ni` quasi-Fermi SG flux using `ni` inferred from the local
  density/QF relation;
- QF-gradient and electric-field impact alpha values.

Key impact-enabled branch result:

```text
bias   edge   Vela/Sentaurus density flux dex   Vela/Sentaurus qf-model flux dex   Vela qf-model/density
-12.7  709    +1.256                          +1.255                              0.998
-12.8  709    +3.093                          +3.092                              0.999
-12.9  709    +3.933                          +3.933                              1.001
-12.7  3646   +4.712                          +4.712                              1.001
-12.8  3646   +6.554                          +6.554                              1.000
-12.9  3646   +7.446                          +7.445                              0.999
```

The no-impact branch gives the same electron-branch result:

```text
bias   edge   Vela/Sentaurus density flux dex   Vela qf-model/density
-12.7  709    +1.256                          0.997
-12.8  709    +3.094                          1.000
-12.9  709    +3.931                          1.000
-12.7  3646   +4.560                          1.000
-12.8  3646   +6.398                          0.999
-12.9  3646   +7.291                          0.999
```

State-variable decomposition on edge `709`:

```text
impact-enabled, edge 709
bias   endpoint electron density dex   mobility dex   QF-field dex   E-field dex
-12.7  +1.34 / +1.09                 -0.086         +0.017        +0.001
-12.8  +3.18 / +2.91                 -0.086         +0.017        +0.001
-12.9  +4.02 / +3.75                 -0.085         +0.017        +0.001
```

State-variable decomposition on edge `3646`:

```text
impact-enabled, edge 3646
bias   endpoint electron density dex   mobility dex   QF-field contrast
-12.7  +0.32 / +0.35                 -0.484         Vela dphin=-1.244e-3 V, Sentaurus dphin=+1.71e-8 V
-12.8  +2.23 / +2.25                 -0.427         Vela dphin=-9.42e-4 V, Sentaurus dphin=+1.71e-8 V
-12.9  +3.13 / +3.15                 -0.384         Vela dphin=-8.42e-4 V, Sentaurus dphin=+1.71e-8 V
```

Interpretation:

- The variable-`ni` QF SG form is internally consistent with the density-form
  SG flux for both Vela and Sentaurus states. `qf_model/density` and
  `qf_inferred/density` are essentially `1.0` on the focus edges.
- The high electron flux is therefore not caused by an isolated algebraic bug
  in `sgElectronContinuityFluxFromQuasiFermiVariableNi`, nor by the
  OldSlotboom `ni_model` on these edges.
- On edge `709`, Vela's flux excess is almost entirely explained by already
  inflated electron density. Mobility is lower in Vela and QF/electric-field
  differences are small.
- On edge `3646`, Sentaurus has nearly flat electron QF and zero electrostatic
  field across the edge, while Vela has a small but finite QF/potential slope.
  That makes the Sentaurus reference flux extremely small, but the Vela
  high-density branch is still the dominant branch-selection symptom.
- The no-impact branch reproduces the electron flux/density jump, confirming
  that impact ionization is not the first trigger for electron branch
  selection.

Updated root-cause hypothesis:

```text
The primary mismatch is not SG flux formula parity. It is Vela's nonlinear
continuation/accepted-state path allowing the electron quasi-Fermi/density
state to move onto a high-electron-density branch between -12.7 and -12.8 V.
Impact ionization then couples that electron branch into a large hole/current
branch.
```

### Next Tasks After Task 53

1. Instrument accepted-step evolution in the no-impact branch over
   `-12.65..-12.85 V` with electron-density branch metrics:
   - max and median electron-density log jump per accepted state;
   - edge `709` and `3646` electron density/QF drops;
   - current and residual/iteration count.
2. Add a temporary electron-branch acceptance guard in sweep diagnostics only:
   reject or flag accepted steps whose electron-density median jump exceeds a
   threshold such as `0.25 dex`, then rerun the no-impact transition window.
   This is a diagnostic guard, not a final physics fix.
3. Compare Sentaurus log step reductions and Bank-Rose nonlinear damping events
   over the same `-12.7..-12.8 V` region. The question is whether Sentaurus
   avoids the high-electron branch by smaller continuation steps/damping, or by
   a different equation/model term.
4. If the diagnostic guard keeps the no-impact branch low-density without
   changing equations, implement a production-quality continuation criterion
   based on carrier-density branch movement. If it does not, inspect contact
   QF anchors and boundary state update around edge `709`.

### Execution Note 2026-06-20: Task 54 No-Impact Accepted-State Evolution

Task 54 added a postprocess diagnostic for accepted VTK/IV sweep states:

```text
scripts/diagnose_pn2d_bv_accepted_state_evolution.py
tests/regression/test_reference_tcad_tools.py::test_pn2d_bv_accepted_state_evolution_reports_density_jump
```

Report generated:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/noimpact_transition_accepted_state_evolution_m1265_m1285/accepted_state_evolution.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/noimpact_transition_accepted_state_evolution_m1265_m1285/accepted_state_edge_evolution.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/noimpact_transition_accepted_state_evolution_m1265_m1285/accepted_state_evolution_summary.json
```

Input branch:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_noimpact_high_transition
window: -12.65..-12.85 V
focus edges: 709, 3646
```

Accepted-state evolution highlights:

```text
step   bias       I(A/um)       iter  accepted_step   median e-jump dex   p95 e-jump dex   max e-jump dex   top node
1458   -12.6500   -8.65e-17     2     -0.00964        --                  --               --               --
1459   -12.6616   -9.54e-17     2     -0.01156        0.052               0.074            0.091            202
1460   -12.6754   -1.05e-16     2     -0.01388        0.051               0.067            0.081            207
1461   -12.6921   -1.15e-16     2     -0.01665        0.044               0.056            0.065            206
1462   -12.7000   -1.20e-16     2     -0.00791        0.015               0.032            0.049            951
1463   -12.7095   -1.48e-16     2     -0.00949        0.112               0.151            0.183            799
1464   -12.7209   -2.09e-16     2     -0.01139        0.176               0.202            0.228            202
1465   -12.7345   -3.35e-16     2     -0.01366        0.235               0.257            0.272            202
1466   -12.7500   -6.04e-16     2     -0.01546        0.285               0.299            0.332            1236
1467   -12.7686   -1.50e-15     3     -0.01855        0.334               0.343            0.433            1236
1468   -12.7908   -3.82e-15     3     -0.02226        0.404               0.425            0.479            1236
1469   -12.8000   -4.99e-15     2     -0.00918        0.128               0.133            0.137            1236
1470   -12.8110   -7.33e-15     2     -0.01102        0.168               0.179            0.199            1236
1471   -12.8242   -1.10e-14     2     -0.01322        0.174               0.184            0.205            1236
1472   -12.8401   -1.62e-14     2     -0.01587        0.164               0.171            0.193            1236
1473   -12.8500   -1.98e-14     2     -0.00989        0.066               0.068            0.084            1238
```

Edge `709` evolution:

```text
bias       e-QF drop V   edge avg e-density cm^-3   avg e-density jump dex
-12.6500   0.4416        5.65e9                    --
-12.7000   0.4450        9.72e9                    0.018
-12.7345   0.4477        4.05e10                   0.266
-12.7500   0.4486        8.20e10                   0.306
-12.7686   0.4497        1.82e11                   0.345
-12.7908   0.4509        4.89e11                   0.430
-12.8000   0.4514        6.66e11                   0.134
-12.8500   0.4541        2.69e12                   0.068
```

Edge `3646` evolution:

```text
bias       e-QF drop V     edge avg e-density cm^-3   avg e-density jump dex
-12.6500   -0.001328       2.72e9                    --
-12.7000   -0.001266       4.51e9                    0.027
-12.7345   -0.001190       1.68e10                   0.269
-12.7500   -0.001135       3.60e10                   0.331
-12.7686   -0.001028       9.74e10                   0.432
-12.7908   -0.000971       2.93e11                   0.479
-12.8000   -0.000968       4.02e11                   0.137
-12.8500   -0.000894       1.93e12                   0.084
```

Interpretation:

- The no-impact high-electron branch develops through a sequence of accepted
  states, not a single failed Newton solve.
- A diagnostic carrier-density branch guard with threshold `0.25 dex` on
  accepted-state electron-density median jump would first flag the path at
  `-12.7345..-12.75 V`, before the larger current rise at `-12.7686` and
  `-12.7908 V`.
- The maximum single-step electron jump in this window is `0.479 dex` at
  `-12.7908 V`, top node `1236`.
- `iterations` stays low (`2..3`) and `retry_count` remains `0`, so the current
  solver accepts the high-electron path as numerically converged. That makes a
  continuation/branch-acceptance diagnostic more relevant than a raw Newton
  failure handler.
- Existing `DCSweep` already has a `sweep.continuation.branch_acceptance`
  framework for terminal-current and `psi-phin` checks. The carrier-density
  guard should attach to that framework rather than creating a separate retry
  mechanism.

### Next Tasks After Task 54

1. Add a diagnostic-only branch acceptance option:
   `sweep.continuation.branch_acceptance.carrier_density_jump`, with
   `max_electron_density_jump_dex` and/or median/p95 variants.
2. Implement it first as a rejection/flag in the existing branch acceptance
   attempt path, using current accepted-state candidate versus previous
   solution; include CSV columns for the measured jump and top node.
3. Rerun the no-impact `-12.65..-12.85 V` transition window with
   `max_electron_density_jump_dex = 0.25` and a small retry shrink factor. Pass
   criterion: the sweep should avoid accepting the high-electron path through
   at least `-12.8 V` or clearly show that smaller steps still enter it.
4. Only after the diagnostic guard result is known, compare Sentaurus
   Bank-Rose/log behavior in the same window to decide whether the final fix
   should be a continuation-control policy, a contact/QF boundary adjustment,
   or a physics/model correction.

### Execution Note 2026-06-20: Task 55 Carrier-Density Branch Guard

Task 55 added an in-solver carrier-density branch diagnostic/guard to the
existing DCSweep continuation branch-acceptance framework.

Code and tests:

```text
include/vela/simulation/DCSweep.h
include/vela/simulation/DCSweepPredictor.h
src/simulation/DCSweep.cpp
tests/test_dc_sweep.cpp
docs/config_schema.md
```

New continuation settings:

```json
"continuation": {
  "branch_acceptance": {
    "carrier_density_jump": true,
    "max_electron_density_jump_dex": 100.0,
    "max_electron_density_jump_p95_abs_dex": 0.25
  }
}
```

New sweep CSV columns when continuation diagnostics are enabled:

```text
electron_density_jump_median_dex
electron_density_jump_p95_abs_dex
electron_density_jump_max_abs_dex
electron_density_jump_max_node
```

Verification:

```text
build-release/test_dc_sweep.exe "[branch_acceptance]"
ctest --test-dir build-release --output-on-failure -R DCSweep
```

Result:

```text
branch_acceptance: 6/6 passed
DCSweep ctest: 50/50 passed
```

No-impact short-window guard probes:

```text
base restart:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/noimpact_branch_high_precision_targets/m12p65/m12p65_latest_state.csv

max guard 0.25:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/noimpact_carrier_branch_guard_m12p65_m12p85_step005

max guard 0.40:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/noimpact_carrier_branch_guard_0p40_m12p65_m12p85_step005

diagnostic high max 100:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/noimpact_carrier_branch_diagnostics_m12p65_m12p85_step005

p95 guard 0.25, max guard 100:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/noimpact_carrier_p95_guard_0p25_m12p65_m12p85_step005
```

Observed behavior:

- `max_electron_density_jump_dex = 0.25` fails almost immediately at
  `-12.650001220703125 V` after 12 retries. The rejected candidate has median
  jump `-9.3e-5 dex`, p95 jump `0.0198 dex`, but max-node jump
  `0.34065 dex` at node `202`.
- `max_electron_density_jump_dex = 0.40` advances only to about
  `-12.650037673 V`, then again rejects node `202` with max jump
  `0.83356 dex` while p95 remains only `0.0161 dex`.
- With max threshold raised to `100 dex`, the short window reaches
  `-12.85 V` in 41 points. Its largest p95 jump is `0.03557 dex` at
  `-12.66 V`; its largest max-node jump is `0.60442 dex` at `-12.655 V`.
- With `max_electron_density_jump_p95_abs_dex = 0.25` and max threshold
  `100 dex`, the same window also reaches `-12.85 V` with zero rejected rows.
  Summary:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/noimpact_carrier_p95_guard_0p25_m12p65_m12p85_step005/carrier_p95_guard_summary.json
points: 41
last_bias: -12.850000000000032 V
max p95: 0.035573107765012324 dex at -12.660000000000002 V
max node jump: 0.60441868968065648 dex at -12.655000000000001 V, node 218
rejected rows: 0
```

Interpretation:

- The nodewise maximum is useful as a spike diagnostic, but too sensitive as
  the primary continuation accept/reject criterion. It can be dominated by one
  local node near startup and force tiny step retries without addressing the
  BV branch question.
- The p95 statistic is much more stable and matches the Task 54 evidence: the
  old no-impact high-electron branch showed p95 jumps around `0.257..0.425 dex`
  through `-12.7345..-12.7908 V`, while the high-precision restart probe shows
  p95 below `0.036 dex`.
- This suggests the immediate high-electron branch is path-dependent: the
  high-precision `m12p65` restart plus small fixed steps stays on a low-jump
  path through `-12.85 V`; the older `official_split_noimpact_high_transition`
  path had already drifted onto a rising electron-density branch.

### Next Tasks After Task 55

1. Run a p95 guard replay on the older `official_split_noimpact_high_transition`
   path, ideally starting from its `-12.65 V` accepted state, to confirm that
   `max_electron_density_jump_p95_abs_dex = 0.25` rejects the historical
   high-electron transition near `-12.7345..-12.75 V`.
2. Export or reconstruct restart CSVs from the older accepted VTK states at
   `-12.65 V`, `-12.70 V`, and `-12.7345 V`, then replay short windows with:
   p95 guard `0.25`, max guard `100`, VTK on for accepted states.
3. Compare the high-precision restart path against Sentaurus at `-12.8 V` and
   `-12.85 V`. If the field/current match improves, prioritize restart/path
   control; if it only hides mismatch, return to contact/QF boundary and
   Sentaurus Bank-Rose damping differences.
4. Keep nodewise max-jump in CSV diagnostics, but do not use it alone as the
   production BV branch criterion unless paired with p95/median or a known
   contact/interior mask.

### Execution Note 2026-06-20: Task 56 Old Official No-Impact Replay

Task 56 reconstructed restart CSVs from the older
`official_split_noimpact_high_transition` VTK accepted states and replayed the
same old bias points with the new carrier-density jump diagnostics.

Generated restart/replay workspace:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/official_noimpact_replay_from_old_vtk
```

Restart states:

```text
restart_from_old_m12p65.csv    from official_split_noimpact_high_transition_1458_-12.65V.vtk
restart_from_old_m12p7.csv     from official_split_noimpact_high_transition_1462_-12.7V.vtk
restart_from_old_m12p7345.csv  from official_split_noimpact_high_transition_1465_-12.7345V.vtk
```

Replay cases:

```text
from_m12p65_diag_p95_100
from_m12p65_guard_p95_0p25
from_m12p7_diag_p95_100
from_m12p7_guard_p95_0p25
from_m12p7345_diag_p95_100
from_m12p7345_guard_p95_0p25
```

High-threshold diagnostic replay results:

```text
case                         points  last bias   max p95 jump                      max node jump
from_m12p65_diag_p95_100     16      -12.85 V    0.1764 dex at -12.6616 V          0.7376 dex at node 204
from_m12p7_diag_p95_100      12      -12.85 V    0.3149 dex at -12.7095 V          0.8106 dex at node 893
from_m12p7345_diag_p95_100   9       -12.85 V    0.5868 dex at -12.7500 V          1.1457 dex at node 893
```

p95 guard replay results with `max_electron_density_jump_p95_abs_dex = 0.25`
and `max_electron_density_jump_dex = 100`:

```text
from_m12p65_guard_p95_0p25    completed 16/16 points to -12.85 V
from_m12p7_guard_p95_0p25     rejected second point at -12.7000092664 V
from_m12p7345_guard_p95_0p25  rejected second point at -12.7345386192 V
```

Rejected rows:

```text
from_m12p7_guard_p95_0p25:
  bias = -12.700009266411524 V
  retry_count = 16
  reason = electron_density_p95_jump_exceeded
  median electron jump = -0.04510 dex
  p95 electron jump = 0.31508 dex
  max electron jump = 0.78298 dex
  max node = 951

from_m12p7345_guard_p95_0p25:
  bias = -12.734538619178075 V
  retry_count = 27
  reason = electron_density_p95_jump_exceeded
  median electron jump = -0.21624 dex
  p95 electron jump = 0.34580 dex
  max electron jump = 0.42888 dex
  max node = 951
```

Restart projection check:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/official_noimpact_replay_from_old_vtk/restart_projection_diff_summary.json
```

Re-solving the first replay point preserves electrostatic potential almost
exactly but changes electron quasi-Fermi and electron density significantly:

```text
old -12.7 V VTK -> replay accepted -12.7 V:
  Potential median/p95/max diff = 0/0/0 V
  ElectronQuasiFermi median diff = +0.02771 V, p95 abs = 0.04380 V
  Electrons median log10 ratio = -0.46555 dex, p95 abs = 0.73557 dex

old -12.7345 V VTK -> replay accepted -12.7345 V:
  Potential median/p95/max diff = 0/0/0 V
  ElectronQuasiFermi median diff = +0.04323 V, p95 abs = 0.04880 V
  Electrons median log10 ratio = -0.72630 dex, p95 abs = 0.81926 dex
```

Interpretation:

- The old accepted VTK states at `-12.7 V` and `-12.7345 V` contain a high
  electron-density branch component. Re-solving the same bias from those states
  projects electron quasi-Fermi/density back toward a lower-density branch while
  leaving electrostatic potential nearly unchanged.
- Once started from the projected `-12.7 V` or `-12.7345 V` state, any small
  reverse-bias continuation step still triggers p95 electron-density jumps over
  `0.25 dex`; the new p95 guard catches that as
  `electron_density_p95_jump_exceeded`.
- Starting from the older `-12.65 V` VTK is recoverable: after re-solving, the
  replay stays below the p95 threshold through `-12.85 V`. This explains why
  the high-precision `m12p65` restart path and the replay-from-old-`-12.65`
  path are both much better behaved than the historical accepted-state sequence.

Current IV comparison for the recovered low-jump path:

```text
compare output:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/official_noimpact_replay_from_old_vtk/from_m12p65_guard_p95_0p25_compare

at -12.8 V:
  Sentaurus current = -8.00338136924e-17 A
  Vela current      = -4.4145100719249983e-17 A/um
  abs log10 error   = 0.25839099993864356 dex
```

The field comparison rows were `missing_input` for `-12.8 V` in the current
Sentaurus multibias field-export directory, so this is an IV-only improvement
check until a matching `-12.8 V` Sentaurus field export is refreshed or located.

### Next Tasks After Task 56

1. Refresh or locate Sentaurus field exports for `-12.8 V` and `-12.85 V`, then
   compare recovered low-jump Vela VTKs against Sentaurus for potential,
   electric field, electron/hole density, and mobility.
2. Promote p95 carrier-density guard from diagnostic to recommended BV
   continuation policy for pn2d reverse sweeps:
   `carrier_density_jump=true`,
   `max_electron_density_jump_p95_abs_dex=0.25`,
   `max_electron_density_jump_dex=100`.
3. Re-run the full Sentaurus-default BV path with impact ionization enabled
   using recovered low-jump continuation controls and compare IV at
   `-0.5, -2, -5, -10, -12.8, -13.2, -20 V`.
4. If the recovered path remains within about `0.3 dex` current error but field
   mismatch persists, prioritize field/density physics parity. If current error
   grows again after impact ionization is enabled, inspect avalanche source and
   terminal current extraction on the recovered branch.

### Execution Note 2026-06-20: Task 57 -12.8 V Field Export and Old-vs-Recovered Comparison

Task 57 imported the existing Sentaurus `-12.8 V` TDR snapshot that had not yet
been exported under `sentaurus_multibias`:

```text
input TDR:
reference_tcad/pn2d_sentaurus2018/source/pn2d_bv_multibias_0128_des.tdr

export:
build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias/sentaurus_-12.8v
```

The import contains the expected fields:

```text
ElectrostaticPotential, ElectricField, eDensity, hDensity,
eQuasiFermiPotential, hQuasiFermiPotential,
eCurrentDensity, hCurrentDensity, TotalCurrentDensity,
ImpactIonization, eMobility, hMobility,
ContactCurrentFlux
```

Recovered p95-guard low-jump comparison:

```text
compare output:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/official_noimpact_replay_from_old_vtk/from_m12p65_guard_p95_0p25_compare

at -12.8 V:
  Sentaurus current = -8.00338136924e-17 A
  Vela current      = -4.4145100719249983e-17 A/um
  current error     = 0.25839099993864356 dex
  potential rms     = 0.003374130062853403 V
  electric field relative p95 = 0.7751949188611165
  electron density log10 p95  = 0.24735097344167065 dex
  hole density log10 p95      = 0.37799529690039074 dex
  electron mobility relative p95 = 0.14739025730442992
  hole mobility relative p95     = 0.21488593498642963
```

Old high-transition path using the same Sentaurus `-12.8 V` export:

```text
compare output:
build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_noimpact_high_transition_compare_m12p8

at -12.8 V:
  Sentaurus current = -8.00338136924e-17 A
  Vela current      = -4.994568288717292e-15 A/um
  current error     = 1.795224443479167 dex
  potential rms     = 0.003374130062853403 V
  electric field relative p95 = 26390.95979329882
  electron density log10 p95  = 2.576327022637035 dex
  hole density log10 p95      = 0.37803024266919766 dex
  electron mobility relative p95 = 0.2717860963379762
  hole mobility relative p95     = 0.21491853374088576
```

Summary artifact:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/official_noimpact_replay_from_old_vtk/m12p8_old_vs_recovered_compare_summary.json
```

Interpretation:

- The recovered p95-guard low-jump branch reduces `-12.8 V` IV error from
  `1.80 dex` to `0.258 dex`.
- Electron density field error drops from `2.58 dex` to `0.247 dex`.
- Electrostatic potential is essentially unchanged between old and recovered
  branches relative to Sentaurus (`~0.00337 V rms`). The main old-path failure
  is therefore not potential, but electron-density/electric-field/current
  branch amplification.
- Remaining recovered-branch gaps at `-12.8 V` are now moderate carrier/mobility
  parity plus electric-field diagnostic/scaling/shape parity, not the previous
  multi-decade electron branch error.

### Next Tasks After Task 57

1. Run an impact-enabled pn2d BV replay with p95 carrier-density guard enabled
   and compare against Sentaurus at `-12.8 V` and `-13.2 V`.
2. If the impact-enabled path still follows the recovered low-density branch,
   inspect avalanche source/current extraction as the next dominant mismatch.
3. If impact feedback forces p95 guard rejection near `-12.7..-12.8 V`, compare
   Sentaurus Bank-Rose damping/log step behavior over the same interval and
   decide whether Vela needs continuation step control, nonlinear damping, or a
   model correction.

### Execution Note 2026-06-20: Task 58 Impact Replay with p95 Carrier Guard

Task 58 tried the next impact-enabled BV replay after Task 57 proved that the
no-impact recovered branch is much closer to Sentaurus at `-12.8 V`.

First attempt:

```text
case:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_to_m13p2

configuration:
  impact ionization enabled with Sentaurus-default pn2d settings
  carrier_density_jump = true
  max_electron_density_jump_p95_abs_dex = 0.25
  max_electron_density_jump_dex = 100
  predictor = linear
  target = 0 V -> -13.2 V

result:
  timed out after 240 s around -4.2389 V
```

Because the full restart-free path was too slow for the current diagnostic
loop, the next run restarted from the recovered old `-12.65 V` state and swept
through the dense `-12.65..-12.85 V` window plus a final `-13.2 V` target:

```text
case:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_from_m12p65_to_m13p2

restart:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/official_noimpact_replay_from_old_vtk/restart_from_old_m12p65.csv

result:
  converged = false
  accepted points = 17
  accepted through = -12.85 V
  failed bias = -12.907819789995315 V
```

Accepted-point diagnostics:

```text
at -12.8 V:
  current = -5.8691526883089689e-17 A/um
  electron_density_jump_p95_abs = 0.003683 dex
  electron_density_jump_max_abs = 0.064134 dex
  terminal_current_consistency_ratio = 0.999858

at -12.85 V:
  current = -5.8585512472682647e-17 A/um
  electron_density_jump_p95_abs = 0.003083 dex
  electron_density_jump_max_abs = 0.055966 dex
  terminal_current_consistency_ratio = 0.666626
```

Failure diagnostics at `-12.907819789995315 V`:

```text
newton_failure_diagnostics:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_from_m12p65_to_m13p2/iv_newton_failure_diagnostics.json

failure:
  handoff_stage = newton_failed
  newton_failure_class = line_search_non_decrease
  failed_iteration = 7
  line_search_attempts = 13
  damping_factor = 0
  residual_norm = 1.4086971951753875
  step_norm = 5

block residuals:
  combined = 20.701807216039953
  psi = 20.50457026655225
  phip = 2.850863061268496
  phin = 5.0301141938103836e-11

carrier diagnostics:
  positive_finite = true
  nonpositive electron/hole counts = 0
  nonfinite electron/hole counts = 0

top Poisson residual nodes:
  node 521 at x=0.03125 um, y=0.34375 um, p-side, residual=-5.2081
  node 3   at x=0.03125 um, y=0.0625 um, p-side, residual=-5.2081
  node 24  at x=0.03125 um, y=0.15625 um, p-side, residual=-5.2081
```

This failure is not a p95 carrier-density guard rejection. The branch guard was
never checked for the failed point because Newton did not converge. The failure
is currently a nonlinear solve/globalization issue dominated by Poisson
residuals at the left p-side/contact-adjacent region after impact feedback is
enabled.

Impact-enabled comparison against Sentaurus at `-12.8 V`:

```text
compare output:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_from_m12p65_to_m13p2_compare

Sentaurus current = -8.00338136924e-17 A
Vela current      = -5.869152688308969e-17 A/um
current error     = 0.1346981039577853 dex

potential rms = 0.003374130062853403 V
electric field relative p95 = 0.7751961321362968
electron density log10 p95 = 0.1275468624032094 dex
hole density log10 p95 = 0.24586861574381308 dex
electron mobility relative p95 = 0.14742787692632722
hole mobility relative p95 = 0.2160682568087798
avalanche_generation log10 p95 = 12.80555318783047 dex
avalanche_generation junction error = 0.4605666768502667 dex
Sentaurus avalanche p99 = 2.9255776772023465e15
Vela avalanche p99 = 1.142798e15
```

Compared with no-impact recovered replay at the same bias, enabling impact
improves current error from `0.258 dex` to `0.135 dex` and electron-density
error from `0.247 dex` to `0.128 dex`. The accepted `-12.8 V` point is therefore
on the desired low-jump branch and closer to Sentaurus. The remaining high
avalanche-generation full-field log error is dominated by low-generation bulk
regions; the junction-region error is much smaller at about `0.46 dex`.

Interpretation:

- The p95 carrier-density guard is doing its job on the recovered branch; it is
  not the blocker for impact-enabled progress beyond `-12.85 V`.
- The next dominant blocker is Newton globalization near `-12.91 V`: Poisson
  residuals dominate, carrier densities remain finite/positive, and line search
  reports non-decrease after all damping attempts.
- Sentaurus' Bank-Rose nonlinear method is directly relevant here. The next
  solver-side comparison should focus on whether Bank-Rose-style damping or a
  residual component merit function can accept a stabilizing step where Vela's
  current line search rejects every trial.
- Physics-side work should continue, but after the nonlinear blocker: at
  `-12.8 V`, the accepted impact branch is already close enough that the next
  high-value field investigations are avalanche spatial support, terminal
  current consistency at `-12.85 V`, and electric-field diagnostic parity.

### Next Tasks After Task 58

1. Inspect the C++ nonlinear line-search/globalization implementation around
   `line_search_non_decrease`: determine which merit norm is used, whether
   Poisson residual dominates scaling, and whether the current step cap of `5`
   is interacting with impact-enabled Poisson residuals near the p-side
   boundary.
2. Run a minimal local continuation experiment from the accepted `-12.85 V`
   state to `-12.91 V` with smaller forced voltage increments and, if currently
   configurable, stronger damping/lower maximum Newton step. The goal is to
   distinguish "no converged solution along this branch" from "line search
   globalization too aggressive for impact feedback".
3. Add or run a residual-probe report for the accepted `-12.85 V` state and the
   failed `-12.9078 V` trial, focusing on p-side/contact-adjacent Poisson
   residuals, contact boundary values, avalanche source contribution, and hole
   quasi-Fermi changes.
4. Once `-12.91 V` is passable or the exact nonlinear blocker is isolated,
   continue the Sentaurus comparison at `-13.2 V`: IV, potential, electric
   field, carrier densities, mobility, avalanche generation, and terminal
   current consistency.

### Execution Note 2026-06-20: Task 59 Isolate the -12.9078 V Failure to Continuation Initial Guess

Task 59 inspected Vela's current Newton line-search implementation and then
reran the failed `-12.907819789995315 V` point with detailed diagnostics.

Implementation notes:

```text
LineSearchConfig defaults:
  initialDamping = solver.damping_factor
  minDamping = 1e-4
  reduction = 0.5
  sufficientDecrease = 1e-4
  maxBacktracks = 12

Newton residual norm:
  solver.residual_norm = "block" by default
  normalized block L2 over psi/phin/phip with configurable weights/scales

Existing damping/update knobs:
  damping_factor
  line_search
  max_update
  quasi_fermi_update_limit_V
  carrier_regularization_scale
  residual_norm/residual_weights/residual_scales
```

The original Task 58 failure did not record per-attempt line-search history
because `solver.diagnostics` was not enabled. A minimal diagnostic replay was
therefore generated:

```text
case:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_m12p85_failure_history

initial_state_file:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_from_m12p65_to_m13p2/latest_state.csv

bias points:
  -12.85 V
  -12.907819789995315 V
```

Result:

```text
converged = true
points = 2

at -12.907819789995315 V:
  Newton iterations = 4
  line_search_attempts per iteration = 1
  damping_factor per iteration = 1
  final residual_norm = 9.473149082055066e-13
  current = -5.8578191860830973e-17 A/um
  electron_density_jump_p95_abs = 0.008806 dex
  electron_density_jump_max_abs = 0.137270 dex
  terminal_current_consistency_ratio = 0.666598
```

The same bias that previously failed with `line_search_non_decrease` converged
cleanly when restarted from the accepted `-12.85 V` state. The previous failure
is therefore not evidence that stronger damping is required for this solution
point. It is more likely caused by the original sweep's long target step to
`-13.2 V`, its linear predictor/retry initial state, or retry-state reuse after
shrinking from the failed final target.

Task 59 then segmented the remaining high-bias interval:

```text
case:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_m12p9078_to_m13p2_segmented

initial_state_file:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_m12p85_failure_history/latest_state.csv

bias points:
  -12.907819789995315
  -12.95
  -13.0
  -13.05
  -13.1
  -13.15
  -13.2
```

Result:

```text
converged = true
points = 7

Newton behavior:
  all non-start points converged with damping = 1
  all line_search_attempts = 1
  no retries
  final residuals at accepted points ~= 8e-12 or lower

carrier-density branch guard:
  max p95 electron jump = 0.0273 dex at -13.0 V
  max node jump can still spike locally, up to 0.844 dex at node 207,
  but p95 remains well below the 0.25 dex branch threshold.

terminal current consistency:
  -12.95 .. -13.1 V: about 0.99959
  -13.15 .. -13.2 V: about 0.66664
```

Sentaurus comparison at `-13.2 V`:

```text
compare output:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_m12p9078_to_m13p2_segmented_compare

Sentaurus current = -8.38472088807e-17 A
Vela current      = -5.858857423553487e-17 A/um
current error     = 0.15567568082397917 dex

potential rms = 0.0033987627186822652 V
electric field relative p95 = 0.7751432169096609
electron density log10 p95 = 0.13972605613955957 dex
hole density log10 p95 = 0.2592378707831501 dex
electron mobility relative p95 = 0.14185824683969814
hole mobility relative p95 = 0.2157402858368836
avalanche_generation log10 p95 = 12.806032696041484 dex
avalanche_generation junction error = 0.47421820326183145 dex
Sentaurus avalanche p99 = 3.19365234469398e15
Vela avalanche p99 = 1.217017e15
```

Interpretation:

- The high-bias impact-enabled low-jump branch is now reproducible through
  `-13.2 V` with ordinary Newton full steps. The original `-12.9078 V`
  `line_search_non_decrease` is best treated as a continuation predictor/retry
  pathology, not as proof that fixed damping must be strengthened.
- The `-13.2 V` current error is only `0.156 dex`, with electron density
  `0.140 dex`, hole density `0.259 dex`, and potential `3.4 mV rms`. This is a
  large improvement over the historical high-density branch and is close enough
  for model/current-extraction debugging.
- The remaining dominant observables are:
  1. electric-field diagnostic parity (`~0.775` relative p95 despite small
     potential RMS);
  2. terminal-current consistency drops to `~0.667` at `-13.15/-13.2 V`;
  3. avalanche-generation spatial support differs strongly in low-generation
     bulk regions, though junction error is about `0.47 dex` and peak location
     is near the same junction/contact corner.

### Next Tasks After Task 59

1. Fix or constrain continuation predictor/retry behavior before adding new
   Bank-Rose-like damping:
   - When a large target step fails and the sweep shrinks the voltage step,
     retry from the last accepted physical state or recompute the predictor for
     the shrunken step.
   - Add a regression test where a failed large target step cannot poison the
     retry initial guess for the smaller step.
   - Keep p95 carrier-density guard enabled as the recommended BV branch guard.
2. Promote a segmented high-bias replay preset for Sentaurus-default pn2d BV:
   `-12.65 -> -12.85`, then `-12.9078 -> -13.2` with explicit points or
   bounded target-step growth, until the predictor retry fix is implemented.
3. Investigate terminal-current consistency at `-13.15/-13.2 V`:
   compare contact-integrated current, edge/face current, and volume generation
   balance at the Anode, because the IV error is now moderate but the internal
   consistency ratio toggles between `~1.0` and `~0.667`.
4. Investigate electric-field and avalanche-generation diagnostics on the
   accepted branch:
   potential RMS is only millivolts, so the large electric-field relative p95
   may be due to gradient reconstruction, mesh/control-volume interpolation, or
   comparison metric sensitivity in low-field regions. Avalanche comparison
   should prioritize thresholded high-generation support near the junction
   before low-generation bulk log errors.

### Execution Note 2026-06-20: Task 60 Retry Predictor Fallback Patch

Task 60 implemented the first code-side fix implied by Task 59: shrunken retry
attempts no longer reuse the linear/secant continuation predictor. They now
restart from the last accepted physical state.

TDD red step:

```text
test:
DCSweep predictor: extrapolates selected coupled variables
section:
linear predictor is disabled for shrunken retry attempts

initial failure:
tests/test_dc_sweep.cpp: too many arguments to predictDCSweepInitialState(...)
```

Implementation:

```text
include/vela/simulation/DCSweepPredictor.h
  predictDCSweepInitialState(..., retryCount = 0)
  returns the current accepted state when retryCount > 0

src/simulation/DCSweep.cpp
  solvePointWithContinuation(..., retryCount = 0)
  bypasses predictor when retryCount > 0
  passes DCSweepStepControl retryCount from both bias_points and regular sweep paths
```

Verification:

```text
cmake --build build-release --target test_dc_sweep --parallel
  passed

build-release/test_dc_sweep.exe "[dc_sweep][continuation][predictor]"
  all tests passed: 49 assertions in 3 test cases

ctest --test-dir build-release --output-on-failure -R DCSweep
  all tests passed: 50/50
```

Full pn2d replay status:

```text
command:
build-release/vela_example_runner.exe --config
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_from_m12p65_to_m13p2/simulation.json

after rebuilding vela_example_runner:
  timed out at 180 s
  output had only progressed to about -12.824 V
```

The full pn2d replay should not yet be counted as verified for the patch. The
targeted evidence remains:

- A clean `-12.85 -> -12.907819789995315 V` restart without predictor converges
  in 4 full Newton steps.
- Segmented `-12.9078 -> -13.2 V` replay converges in full steps and matches
  Sentaurus at `-13.2 V` to `0.156 dex` current error.
- Automated DCSweep tests now enforce retry predictor fallback.

### Next Tasks After Task 60

1. Run a bounded real-case retry-regression deck that seeds predictor history
   near `-12.85 V`, intentionally targets `-13.2 V`, and verifies that the
   first shrunken retry records `predicted_initial_state=0` and converges.
   This should be much smaller than the full `-12.65 -> -13.2` replay.
2. Add a CSV-level regression if feasible: a failed first attempt with
   `retry_count > 0` should write `predicted_initial_state=0` for the accepted
   retry point.
3. After the bounded deck passes, re-run the full impact+p95 pn2d replay with a
   longer timeout and compare `-13.2 V` against the segmented reference.
4. Continue physics/current work on the now-stable branch:
   terminal-current consistency around `-13.15/-13.2 V`, electric-field
   reconstruction parity, and thresholded avalanche-generation support.

### Execution Note 2026-06-20: Task 61 Bounded Retry Probe and Predictor-Branch Classification

Task 61 followed up Task 60 with bounded real-case probes instead of another
long full replay.

Bounded retry-regression attempt:

```text
output:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression

bias_points:
  -12.85 V
  -12.907819789995315 V
  -13.2 V

result:
  converged all 3 points
  no retry was triggered
  final current at -13.2 V = -7.323944307607064e-17 A/um
  terminal_current_consistency_ratio = 1.0001000041295318
  electron_density_jump_p95_abs at -13.2 V = 0.023251638410184405 dex
```

Reduced-Newton variants did not provide the intended accepted retry sample:

```text
impact_p95_guard_bounded_retry_regression_newton6
  failed before the final target

impact_p95_guard_bounded_retry_regression_newton8
  reached -13.2 V but failed by max_iterations

impact_p95_guard_regular_retry_regression_newton8
  converged a regular step sequence but did not exercise the long target
  retry path; all accepted rows kept retry_count=0
```

Therefore the requested CSV-level real-case proof that the first shrunken retry
writes `predicted_initial_state=0` is still missing. The unit test remains the
only verified coverage for retry predictor fallback, and a synthetic solver hook
or refactor would be needed to make this specific CSV behavior deterministic.

The bounded probe exposed a more important branch-selection clue. The direct
long-step predictor branch at `-13.2 V` is closer to Sentaurus terminal current
than the segmented p95 branch, even though their field errors are nearly the
same:

```text
direct long-step branch:
  Vela current = -7.323944307607064e-17 A/um
  Sentaurus current = -8.38472088807e-17 A
  current error = 0.05874357711800847 dex
  potential RMS = 0.0033987627186822652 V
  electric-field relative p95 = 0.7751432169096609
  electron-density log10 p95 = 0.140394299016605 dex
  hole-density log10 p95 = 0.25921021594805405 dex
  avalanche-generation log10 p95 = 12.851162178734057 dex

segmented branch:
  Vela current = -5.858857423553487e-17 A/um
  Sentaurus current = -8.38472088807e-17 A
  current error = 0.15567568082397917 dex
  potential RMS = 0.0033987627186822652 V
  electric-field relative p95 = 0.7751432169096609
  electron-density log10 p95 = 0.13972605613955957 dex
  hole-density log10 p95 = 0.2592378707831501 dex
  avalanche-generation log10 p95 = 12.806032696041484 dex
```

The main observed difference is not the global potential, density, mobility, or
field comparison. It is the terminal-current decomposition:

```text
direct long-step branch at -13.2 V:
  total = -7.323944307607064e-17 A/um
  electron = -5.789104627816346e-19 A/um
  hole = 7.266053261328900e-17 A/um
  hole drift = 7.265320909891848e-17 A/um
  hole diffusion = 0
  terminal consistency = ~1.0001

segmented branch at -13.2 V:
  total = -5.858857423553487e-17 A/um
  electron = -5.807759822499503e-19 A/um
  hole = 5.800779825328491e-17 A/um
  hole drift = 7.265320909891598e-17 A/um
  hole diffusion = -1.465273436000408e-17 A/um
  terminal consistency = ~0.6666
```

No existing `contact_edges.csv` or `terminal_balance.csv` files were present in
the two candidate branch output directories, so the next diagnostic needs to
generate them deliberately or reconstruct the contact-edge quantities offline
from the saved states.

### Next Tasks After Task 61

1. Treat the direct long-step branch as a candidate Sentaurus-like branch for
   terminal-current purposes; do not force the segmented branch solely because
   it uses smaller voltage steps.
2. Generate single-point contact-edge and terminal-balance diagnostics for both
   `impact_p95_guard_bounded_retry_regression` and
   `impact_p95_guard_m12p9078_to_m13p2_segmented`, preferably without changing
   the solved states. If `vela_example_runner` must be used, record that it is
   a re-solve/projection probe and may perturb the branch.
3. Compare contact-edge sums against IV current at `-13.2 V`, then inspect
   anode hole-current components and hole quasi-Fermi endpoint drops. The key
   discriminator is whether the segmented branch's `-1.47e-17 A/um` hole
   diffusion term is a real solved-state difference, a contact projection/QF
   floor artifact, or a terminal-current reporting convention.
4. Keep electric-field and avalanche work secondary until the terminal-current
   branch is classified: both candidate branches show similar field and density
   parity but different terminal-current consistency.

### Execution Note 2026-06-20: Task 62 Contact-Current QF-Floor Probe at -13.2 V

Task 62 generated canonical single-point diagnostics for the two `-13.2 V`
candidate branches. The probes write to separate directories and do not
overwrite the original branch outputs:

```text
default re-solve/projection probes:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_contact_probe
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_m12p9078_to_m13p2_segmented_contact_probe

qf-floor reporting probes:
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_contact_qf_floor_probe
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_m12p9078_to_m13p2_segmented_contact_qf_floor_probe
```

The default probes are not pure offline readers. They re-enter DCSweep at the
same bias from `latest_state.csv`, so contact-boundary projection and current
reporting are recomputed.

Default probe result:

```text
bounded retry latest_state -> default probe:
  original current at -13.2 V = -7.323944307607064e-17 A/um
  probe current at -13.2 V    = -5.858870183437891e-17 A/um
  Newton iterations = 1
  anode edge sum matches probe IV exactly
  edge-vs-IV relative error = 0

segmented latest_state -> default probe:
  original current at -13.2 V = -5.858857423553487e-17 A/um
  probe current at -13.2 V    = -5.858857423553487e-17 A/um
  Newton iterations = 0
  anode edge sum matches probe IV to 2.1e-16 relative error
```

This proves the long-step branch's Sentaurus-like reported current is not
preserved by a plain same-bias reload/reprojection probe. The saved bulk state
reloads onto the segmented-style contact-current report unless the original
contact-current micro-floor is preserved.

The `contact_current_qf_floor` probes isolate that missing state. They capture
the anode contact-edge hole-QF drop from the initial state and use it only for
terminal/contact-current reporting:

```text
bounded retry qf-floor probe:
  current = -7.324143619438301e-17 A/um
  current error vs Sentaurus = 0.058731758506743864 dex
  anode qf-floor override edges = 17/17
  edge sum matches IV to 1.68e-16 relative error

segmented qf-floor probe:
  current = -5.858857423553487e-17 A/um
  current error vs Sentaurus = 0.15567568082397917 dex
  anode qf-floor override edges = 17/17
  edge sum matches IV to 2.1e-16 relative error
```

The qf-floor reporting script gives the root discriminator:

```text
bounded retry default probe:
  baseline_total_current = -5.858870183437891e-17 A/um
  restart_drop_total_current = -7.324143619438298e-17 A/um
  baseline hole-QF drop = -7.105427357601002e-15 V
  restart hole-QF drop  = -8.881784197001252e-15 V
  baseline error = 0.155674734983936 dex
  restart-drop error = 0.058731758506640475 dex

segmented default probe:
  baseline_total_current = -5.858857423553486e-17 A/um
  restart_drop_total_current = -5.858857423553486e-17 A/um
  baseline hole-QF drop = -7.105427357601002e-15 V
  restart hole-QF drop  = -7.105427357601002e-15 V
  baseline/restart error = 0.15567568082387576 dex
```

Both qf-floor probe verification commands passed:

```text
bounded retry:
qf-floor restart point verified: bias=-13.2 V contact=Anode
log10_error=0.0587318 override_edges=17

segmented:
qf-floor restart point verified: bias=-13.2 V contact=Anode
log10_error=0.155676 override_edges=17
```

Task 62 conclusion:

- For the observed `-13.2 V` terminal-current branch discrepancy, the
  bounded/segmented anode contact-edge hole quasi-Fermi micro-drop difference
  is highly consistent with the current difference, and the qf-floor override
  reproduces the bounded branch's reported current. The measured drops differ
  at one floating-point ULP scale (`-8.88e-15 V` versus `-7.11e-15 V`).
- The earlier field comparison outputs, not the contact-current probes alone,
  show that the two candidate branches have nearly identical bulk
  potential/density/mobility metrics at `-13.2 V`.
- The canonical contact-edge sum is internally consistent with each reported IV
  current. This is not a terminal-edge aggregation bug.
- A plain reload/re-solve destroys the long-step branch's reporting state, while
  `contact_current_qf_floor` restores it for the bounded branch. For this
  `-13.2 V` terminal-current discrepancy, the next implementation focus should
  be preservation/definition of contact-current reporting endpoint floors
  across continuation; this does not globally rule out later mobility, SG-flux,
  or avalanche refinements for other observables.
- The segmented qf-floor probe currently reproduces the restart drop preserved
  in that segmented state. Its synthetic ULP-floor candidate would move the
  current toward the bounded value, but the enabled policy does not force that
  synthetic floor.

### Next Tasks After Task 62

1. Decide whether `contact_current_qf_floor` should graduate from diagnostic to
   an opt-in BV reporting policy, with a name that reflects its purpose:
   preserving contact-current endpoint QF micro-drops from the accepted
   continuation state.
2. Add a regression around restart/reporting stability:
   - direct `-13.2 V` accepted point writes a contact-current reporting floor or
     equivalent reproducible metadata;
   - reloading the state and reporting terminal current reproduces the original
     `-7.32e-17 A/um` current when the policy is enabled;
   - default reporting remains unchanged unless the policy is explicitly
     enabled.
3. Compare against Sentaurus terminal-current definition more directly:
   current error is now `~0.059 dex`, so the remaining gap is small enough to
   inspect Sentaurus `ContactCurrentFlux`/PLT convention and anode edge endpoint
   floating-point floor behavior before touching avalanche.
4. Keep the p95 carrier-density guard and retry predictor fallback as solver
   stability tools, but do not use them to choose between two terminal-current
   reporting policies without the contact-QF-floor metadata.

### Execution Note 2026-06-20: Task 63 Opt-in Contact-Current Reporting Policy Entry

Task 63 promoted the existing diagnostic-only QF-floor reporting hook to a
clear opt-in reporting-policy entry while preserving the old diagnostics alias.

New preferred config:

```json
"sweep": {
  "contact_current_reporting": {
    "endpoint_qf_floor": {
      "enabled": true,
      "contacts": ["Anode"]
    }
  }
}
```

Compatibility config retained:

```json
"sweep": {
  "diagnostics": {
    "contact_current_qf_floor": {
      "enabled": true,
      "contacts": ["Anode"]
    }
  }
}
```

Implementation:

```text
src/simulation/DCSweep.cpp
  parses sweep.contact_current_reporting.endpoint_qf_floor
  maps it to the existing reporting-only contact-current QF floor mechanism
  leaves the default disabled behavior unchanged

docs/config_schema.md
  documents contact_current_reporting.endpoint_qf_floor as the preferred
  opt-in policy and marks diagnostics.contact_current_qf_floor as a
  compatibility alias
```

TDD red step:

```text
test:
DCSweep: contact current reporting policy preserves initial endpoint QF drops

failure:
REQUIRE( sawOverride )
with expansion: false
```

Green verification:

```text
cmake --build build-release --target test_dc_sweep --parallel
build-release/test_dc_sweep.exe "DCSweep: contact current reporting policy preserves initial endpoint QF drops"

result:
All tests passed (8 assertions in 1 test case)

build-release/test_dc_sweep.exe "[dc_sweep][diagnostics][contact_current_qf_floor]"
result:
All tests passed (30 assertions in 2 test cases)

build-release/test_dc_sweep.exe "[dc_sweep][contact_current_reporting]"
result:
All tests passed (8 assertions in 1 test case)
```

Review fix:

```text
issue:
  With continuation predictor enabled, the new reporting policy initially
  skipped external initial_state_file QF-floor capture because the solve used
  the predicted state as the Newton initial guess and disabled reporting-floor
  capture.

test update:
  DCSweep: contact current reporting policy preserves initial endpoint QF drops
  now enables continuation.predictor.mode = "constant" and verifies that
  contact-edge override rows preserve phip1 - phip0 = 1e-6 V from the external
  initial state.

implementation update:
  solvePointWithContinuation still solves from the predicted state, but when
  endpoint_qf_floor capture is allowed it builds contact-current overrides from
  the external initial state and attaches them to the attempt before reporting.
```

Post-review verification:

```text
build-release/test_dc_sweep.exe "DCSweep: contact current reporting policy preserves initial endpoint QF drops"
  All tests passed (16 assertions in 1 test case)

build-release/test_dc_sweep.exe "[dc_sweep][diagnostics][contact_current_qf_floor]"
  All tests passed (30 assertions in 2 test cases)

build-release/test_dc_sweep.exe "[dc_sweep][continuation][predictor]"
  All tests passed (49 assertions in 3 test cases)

ctest --test-dir build-release --output-on-failure -R DCSweep
  51/51 tests passed
```

Task 63 scope caveat:

- This is not a default physics/model change. It only provides a production
  naming path for the reporting-only behavior already validated in Task 62.
- The policy still captures endpoint QF micro-drops only from an external
  `initial_state_file`, not from ordinary in-memory continuation states.

### Next Tasks After Task 63

1. Add an end-to-end restart/reporting stability regression using the real
   `-13.2 V` bounded branch artifact or a compact fixture derived from it:
   default reload should report the baseline current, while
   `contact_current_reporting.endpoint_qf_floor` should reproduce the preserved
   restart-drop current.
2. If that regression is too heavy for normal CTest, add a Python regression
   that consumes the generated probe outputs and verifies the Task 62
   invariants.
3. Only after the restart/reporting regression is in place, consider enabling
   the policy in the recommended pn2d BV validation deck.

### Execution Note 2026-06-20: Task 64 Restart/Reporting Stability Regression

Task 64 implemented the Task 63 follow-up regression gate without adding a
large generated artifact to ordinary tests.

New verifier:

```text
scripts/verify_pn2d_bv_qf_floor_reporting_stability.py
```

The verifier compares one default reload/projection probe and one
`endpoint_qf_floor` probe against a restart-drop summary. It requires:

```text
default IV current ~= baseline_total_current_A_per_um
qf-floor IV current ~= restart_drop_total_current_A_per_um
default contact edges have zero qf-floor overrides
qf-floor contact edges have the expected override count
contact-edge sums match their IV currents
qf-floor log-current error improves over default log-current error
```

TDD red step:

```text
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_qf_floor_reporting_stability_compares_default_and_policy

failure:
can't open file
scripts/verify_pn2d_bv_qf_floor_reporting_stability.py
```

Green verification with compact fixture:

```text
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_qf_floor_reporting_stability_compares_default_and_policy

result:
Ran 1 test in 0.305s
OK
```

Real Task 62 bounded-probe verification:

```powershell
python scripts/verify_pn2d_bv_qf_floor_reporting_stability.py `
  --default-iv-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_contact_probe/iv.csv `
  --default-contact-edge-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_contact_probe/contact_edges.csv `
  --qf-floor-iv-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_contact_qf_floor_probe/iv.csv `
  --qf-floor-contact-edge-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_contact_qf_floor_probe/contact_edges.csv `
  --restart-summary build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_contact_probe/contact_qf_floor_reporting/contact_qf_floor_reporting_summary.json `
  --bias -13.2 `
  --contact Anode `
  --expected-override-edges 17 `
  --out-json build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_contact_qf_floor_probe/qf_floor_reporting_stability_verification.json
```

Result:

```text
qf-floor reporting stability verified:
bias=-13.2 V contact=Anode
default_log10_error=0.155675
qf_floor_log10_error=0.0587318
```

Focused Python regression verification:

```text
python -m unittest \
  tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_contact_qf_floor_reporting_compares_restart_and_ulp_policies \
  tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_qf_floor_restart_point_verifier_checks_enabled_current \
  tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_qf_floor_reporting_stability_compares_default_and_policy

result:
Ran 3 tests in 1.045s
OK
```

Task 64 conclusion:

- The Task 62 bounded-branch observation is now covered by a reusable verifier
  and a compact hermetic regression fixture.
- Ordinary tests do not depend on the large `build-release` probe artifacts,
  while the same verifier can be run against those artifacts when available.
- The next implementation step can safely consider adding
  `contact_current_reporting.endpoint_qf_floor` to a recommended pn2d BV
  validation/probe deck, but it should remain opt-in until the Sentaurus
  terminal-current convention is compared directly.

### Next Tasks After Task 64

1. Add the reporting policy only to the recommended BV debug/validation deck,
   not to global defaults.
2. Re-run the bounded `-13.2 V` probe/deck with
   `contact_current_reporting.endpoint_qf_floor` in config form, then run the
   new stability verifier on the generated default-vs-policy outputs.
3. Compare Sentaurus `ContactCurrentFlux` and PLT terminal current at `-13.2 V`
   before treating the remaining `~0.059 dex` current gap as a Vela transport
   error.

### Execution Note 2026-06-20: Task 65 Config-form BV Reporting Policy Probe

Task 65 executed the Task 64 next step with the preferred production-style
configuration key instead of the legacy diagnostics alias.

Generated probe:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_contact_reporting_policy_probe
```

Config distinction:

```json
"contact_current_reporting": {
  "endpoint_qf_floor": {
    "enabled": true,
    "contacts": ["Anode"]
  }
}
```

The probe was generated from the default contact probe and deliberately removed
`diagnostics.contact_current_qf_floor`, so this validates the new config form.
It does not change any global default or production fixture.

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build-release --target vela_example_runner --parallel
build-release\vela_example_runner.exe --config build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_contact_reporting_policy_probe\simulation.json
```

Result:

```text
converged = true
points = 1

at -13.2 V:
  current_total_A_per_um = -7.324143619438301e-17
  current_hole_A_per_um = 7.266053261328900e-17
  current_hole_drift_A_per_um = 7.265320909891597e-17
  current_hole_diffusion_A_per_um = -1.465273436000408e-17
  Newton iterations = 1
  Anode qf-floor override edges = 17/17
  Cathode qf-floor override edges = 0/17
```

Contact-current extraction:

```text
edge_current_A_per_um = -7.324143619438300e-17
edge_vs_iv_relative_error = 1.682920527577476e-16
Sentaurus PLT current = -8.384720888068e-17 A
Sentaurus TDR ContactCurrentFlux = -7.269471996838386e-17 A
Vela policy current error vs PLT = -0.058731758506743864 dex
Vela policy current error vs TDR ContactCurrentFlux = 0.0032539838945788765 dex
Vela policy relative error vs TDR ContactCurrentFlux = 0.007464575442624556
Sentaurus TDR ContactCurrentFlux vs PLT relative difference = 0.1330096619932436
Sentaurus TDR ContactCurrentFlux vs PLT log10 difference = -0.061985742401219186 dex
```

Stability gate:

```powershell
python scripts/verify_pn2d_bv_qf_floor_reporting_stability.py `
  --default-iv-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_contact_probe/iv.csv `
  --default-contact-edge-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_contact_probe/contact_edges.csv `
  --qf-floor-iv-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_contact_reporting_policy_probe/iv.csv `
  --qf-floor-contact-edge-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_contact_reporting_policy_probe/contact_edges.csv `
  --restart-summary build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_contact_probe/contact_qf_floor_reporting/contact_qf_floor_reporting_summary.json `
  --bias -13.2 `
  --contact Anode `
  --expected-override-edges 17 `
  --out-json build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_contact_reporting_policy_probe/qf_floor_reporting_stability_verification.json
```

Result:

```text
qf-floor reporting stability verified:
bias=-13.2 V contact=Anode
default_log10_error=0.155675
qf_floor_log10_error=0.0587318
```

Task 65 conclusion:

- The preferred `contact_current_reporting.endpoint_qf_floor` config form
  reproduces the Task 62/64 qf-floor behavior without relying on the legacy
  diagnostics alias.
- The policy probe is internally terminal-extraction consistent: contact-edge
  sum and IV agree to numerical precision.
- The policy current is very close to Sentaurus TDR `ContactCurrentFlux`
  (`0.00325 dex`, `0.746%` relative), while the remaining PLT current gap is
  nearly the same size as Sentaurus's own TDR `ContactCurrentFlux` versus PLT
  discrepancy (`0.06199 dex`, `13.3%` relative). Therefore the remaining
  `~0.059 dex` gap should be treated as a Sentaurus terminal-output convention
  investigation before changing Vela mobility, SG flux, or avalanche transport.
- `docs/config_schema.md` now includes the corresponding BV workflow note:
  compare Sentaurus `.plt` terminal current and TDR `ContactCurrentFlux` before
  tuning Vela transport against the remaining terminal-current gap.

Final verification for this task:

```text
python -m unittest \
  tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_contact_qf_floor_reporting_compares_restart_and_ulp_policies \
  tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_qf_floor_restart_point_verifier_checks_enabled_current \
  tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_qf_floor_reporting_stability_compares_default_and_policy

result:
Ran 3 tests in 1.379s
OK

ctest --test-dir build-release --output-on-failure -R DCSweep
result:
51/51 tests passed
```

### Next Tasks After Task 65

1. Keep `contact_current_reporting.endpoint_qf_floor` opt-in and recommend it
   only for Sentaurus restart/BV terminal-current parity probes.
2. Add a small note to the user-facing BV workflow or config schema that
   Sentaurus PLT current and TDR `ContactCurrentFlux` can differ at the
   `~0.06 dex` level near `-13.2 V`; compare both before tuning transport.
3. Continue with thresholded avalanche-support and electric-field reconstruction
   diagnostics only after terminal-output convention is documented.

### Execution Note 2026-06-20: Task 66 Sentaurus Terminal-current Cross-check Gate

Task 66 turns the Task 65 terminal-output observation into a reusable checker.
This keeps later BV work from mistaking a Sentaurus output-convention mismatch
for a Vela transport error.

Added:

```text
scripts/verify_pn2d_sentaurus_terminal_current_crosscheck.py
```

The checker reads one Sentaurus TDR export directory, extracts the contact
`ContactCurrentFlux`, reads the matching `.plt` terminal current column, and
writes a JSON summary with:

```text
contact_current_flux_A
plt_current_A
relative_difference
log10_contact_flux_over_plt
status
```

Regression:

```text
tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.
test_pn2d_bv_sentaurus_terminal_current_crosscheck_reports_flux_plt_mismatch
```

Real `-13.2 V` PN2D result:

```powershell
python scripts/verify_pn2d_sentaurus_terminal_current_crosscheck.py `
  --sentaurus-root build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias `
  --sentaurus-plt build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_vm_runs/pn2d_bv_vm_default/source/pn2d_bv.plt `
  --bias -13.2 `
  --contact Anode `
  --max-relative-difference 0.01 `
  --out-json build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/sentaurus_terminal_current_crosscheck_m13p2.json
```

Output:

```text
status = sentaurus_plt_contact_flux_mismatch
Sentaurus TDR ContactCurrentFlux = -7.269471996838386e-17 A
Sentaurus PLT Anode TotalCurrent = -8.384720888068e-17 A
relative_difference = 0.1330096619932436
log10_contact_flux_over_plt = -0.061985742401219186 dex
```

Conclusion:

- The earlier `~0.059 dex` Vela policy-vs-PLT residual is inside the same
  scale as Sentaurus's own PLT-vs-TDR terminal discrepancy.
- Before tuning mobility, SG flux, thresholded avalanche support, or electric
  field reconstruction against PLT current, run this cross-check and compare
  Vela against both Sentaurus terminal definitions.
- Next BV debug step can now move to electric-field reconstruction and
  thresholded avalanche support with terminal-current convention documented.

### Execution Note 2026-06-20: Task 67 BV Debug Direction Summary

Task 67 combines the now-separated terminal-current evidence with the existing
field-ranking evidence, so the next debug step is selected from current data
instead of from the older pre-terminal-crosscheck ranking alone.

Added:

```text
scripts/summarize_pn2d_bv_debug_direction.py
```

Regression:

```text
tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.
test_pn2d_bv_debug_direction_summarizes_terminal_and_field_evidence
```

Real-data command:

```powershell
python scripts/summarize_pn2d_bv_debug_direction.py `
  --debug-ranking build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_compare/debug_ranking.json `
  --terminal-crosscheck build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/sentaurus_terminal_current_crosscheck_m13p2.json `
  --qf-floor-stability build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_contact_reporting_policy_probe/qf_floor_reporting_stability_verification.json `
  --out-json build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/bv_debug_direction_after_task66.json
```

Output summary:

```text
terminal_assessment = plt_gap_explained_by_sentaurus_terminal_convention
qf_floor_abs_log10_error = 0.05873175850664031
Sentaurus ContactCurrentFlux-vs-PLT log10 gap = -0.061985742401219186 dex
Sentaurus ContactCurrentFlux-vs-PLT relative gap = 0.1330096619932436

priority 1 = thresholded_avalanche_support
  avalanche field_error = 12.851162178734057 dex
  avalanche junction_error = 0.4742765511395509 dex
  avalanche_status = thresholded_peak

priority 2 = junction_electric_field_reconstruction
  electric_field field_error = 0.7751432169096609 relative p95
  electric_field junction_error = 0.10751399902576818 relative p95

priority 3 = carrier_density_mobility_followup
  secondary quantities = electron_density, electron_mobility, hole_density,
  hole_mobility
```

Task 67 conclusion:

- The remaining Vela-vs-PLT terminal-current gap is no larger than Sentaurus's
  own TDR `ContactCurrentFlux` vs PLT terminal-current gap, so it should not
  drive immediate transport tuning.
- The next physics-side discrepancy is not the global electric-field relative
  p95 alone: global field error is denominator-sensitive, while the junction
  field error is only `0.1075`.
- Avalanche generation is the sharper next discriminator. Its global/log p95
  mismatch is dominated by thresholded/near-floor support, while the junction
  mismatch is `0.4743 dex`. The next task should compare high-generation
  support sets around the junction, including Sentaurus/Vela p99 masks,
  peak-node neighborhoods, and SG edge source ownership.

### Next Tasks After Task 67

1. Add or run a thresholded avalanche-support comparator at `-13.2 V`:
   - define Sentaurus and Vela active masks from each side's p99 or a shared
     absolute generation threshold;
   - report overlap/Jaccard, false-positive support, false-negative support,
     peak coordinate separation, and junction-local source integrals;
   - use the existing `thresholded_peak` field-ranking evidence to avoid
     interpreting low-generation near-floor nodes as physical error.
2. Reconstruct the junction electric-field stencil only on the active avalanche
   support:
   - compare Sentaurus exported `ElectricField`, Vela VTK `ElectricField`,
     and finite-difference/element-gradient reconstruction from potential;
   - classify whether the `0.1075` junction field error is due to field export
     convention, interpolation, or actual potential-gradient mismatch.
3. Only after the active-support and junction-field checks, revisit
   carrier-density/mobility:
   - electron density and mobility remain secondary because their current
     junction errors are smaller than the thresholded avalanche support signal;
   - tune mobility/SG flux only if support and field reconstruction are already
     aligned.

### Execution Note 2026-06-20: Task 68 Thresholded Avalanche-support Comparator

Task 68 executes the first Task 67 follow-up: compare high-generation avalanche
support sets instead of relying on global log-p95 error dominated by near-floor
nodes.

Added:

```text
scripts/diagnose_pn2d_bv_thresholded_avalanche_support.py
```

Regression:

```text
tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.
test_pn2d_bv_thresholded_avalanche_support_reports_mask_mismatch
```

Important implementation note:

- Vela VTK writes ASCII coordinates such as `1.04688 um` where Sentaurus CSV has
  `1.046875 um`. The comparator therefore uses a default coordinate
  quantization tolerance of `1e-4 um`, still much smaller than the local grid
  spacing (`0.015625 um`), to avoid dropping nodes from the support comparison.

Real-data command:

```powershell
python scripts/diagnose_pn2d_bv_thresholded_avalanche_support.py `
  --sentaurus-dir build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias/sentaurus_-13.2v `
  --vela-vtk build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression/vtk/impact_p95_bounded_retry_0002_-13.2V.vtk `
  --out-dir build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/thresholded_avalanche_support_m13p2 `
  --percentile 99
```

Side-specific p99 support result:

```text
matched_points = 1943
missing_vela_points = 0
sentaurus_active_count = 20
vela_active_count = 20
overlap_count = 10
false_positive_count = 10
false_negative_count = 10
union_count = 30
jaccard = 0.3333333333333333
peak_separation_um = 0.015625
Sentaurus active source sum = 6.393206116376254e16 cm^-3 s^-1
Vela active source sum = 2.452106e16 cm^-3 s^-1
Vela/Sentaurus active source-sum ratio = 0.383559658
```

Shared absolute-threshold checks:

```text
threshold = Sentaurus p99 = 3.193652344693980e15 cm^-3 s^-1:
  Sentaurus active = 20
  Vela active = 0
  overlap = 0
  false_negative = 20
  jaccard = 0

threshold = Vela p99 = 1.217856600000000e15 cm^-3 s^-1:
  Sentaurus active = 660
  Vela active = 20
  overlap = 20
  false_negative = 640
  jaccard = 0.030303030303030304
```

Support geometry:

- The dominant high-generation support is a vertical high-field line near
  `x = 1.0 um`.
- Sentaurus peak: node `351`, `(x, y) = (1.0, 0.015625) um`,
  `ImpactIonization = 3.1978521474416165e15 cm^-3 s^-1`.
- Vela peak: node `352`, `(x, y) = (1.0, -0.0) um`,
  `AvalancheGeneration = 1.24696e15 cm^-3 s^-1`.
- The peak location differs by one grid step in `y` (`0.015625 um`), but both
  peak nodes are in the side-specific p99 overlap set.

Task 68 conclusion:

- The high-generation support mismatch is not just a low-source near-floor
  artifact. With side-specific p99 masks, only half of the top support overlaps.
- The larger discriminator is source magnitude: Vela's active p99 support sum
  is only about `38.36%` of Sentaurus's, and no Vela node reaches the Sentaurus
  p99 absolute threshold.
- This points the next debug step at local avalanche source formation on the
  high-field line: compare Sentaurus/Vela electric field, high-field driving
  force, electron/hole current-density contribution, mobility, and alpha values
  on the same p99 overlap/false-negative nodes before changing global mobility
  or terminal-current reporting.

### Next Tasks After Task 68

1. Build a p99-support local-factor table for the nodes/edges in
   `thresholded_avalanche_support_m13p2`:
   - for overlap, false-negative, and false-positive support nodes, join
     Sentaurus `ElectricField`, `eCurrentDensity`, `hCurrentDensity`,
     `eMobility`, `hMobility`, `eDensity`, `hDensity`, and Vela VTK
     `ElectricField`, high-field drive, mobility, density, and source;
   - compute Vela/Sentaurus ratios on source, field, mobility, current-density
     magnitude, and carrier densities.
2. Reconstruct the Vela avalanche source on these nodes from its stored local
   quantities:
   - evaluate whether the `~0.3836` source-sum ratio is explained by field/alpha
     magnitude, current-density magnitude, source-volume assignment, or a
     unit/convention mismatch;
   - check whether false-negative nodes correspond to alternating grid rows or
     a systematic boundary/junction interpolation offset.
3. If local factors show field/alpha is dominant, proceed to Task 67 priority 2:
   junction electric-field reconstruction on the same active-support set.

### Execution Note 2026-06-20: Task 69 P99-support Local-factor Table

Task 69 executes the first Task 68 follow-up: join the p99 support nodes with
local Sentaurus/Vela field, driving-force, mobility, density, current-density,
and source quantities.

Added:

```text
scripts/diagnose_pn2d_bv_support_local_factors.py
```

Regression:

```text
tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.
test_pn2d_bv_support_local_factors_summarizes_source_ratios
```

Real-data command:

```powershell
python scripts/diagnose_pn2d_bv_support_local_factors.py `
  --support-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/thresholded_avalanche_support_m13p2/thresholded_avalanche_support_nodes.csv `
  --sentaurus-dir build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias/sentaurus_-13.2v `
  --vela-vtk build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression/vtk/impact_p95_bounded_retry_0002_-13.2V.vtk `
  --out-dir build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/support_local_factors_m13p2
```

Output summary:

```text
row_count = 30

overlap:
  count = 10
  source_ratio_median = 0.38301973033403697
  electric_field_ratio_median = 0.9996629041167242
  electron_drive_over_sentaurus_field_median = 1.007672007061604
  hole_drive_over_sentaurus_field_median = 1.0080010435304079
  electron_mobility_ratio_median = 0.9441711172275149
  hole_mobility_ratio_median = 1.2127105440764687
  electron_density_ratio_median = 0.8186974018624151
  hole_density_ratio_median = 0.5660786382619549

false_negative:
  count = 10
  source_ratio_median = 0.378348244854814
  electric_field_ratio_median = 0.9996596355425514
  electron_drive_over_sentaurus_field_median = 1.0078278090981103
  hole_drive_over_sentaurus_field_median = 1.00796073111555
  electron_mobility_ratio_median = 0.9440854802394835
  hole_mobility_ratio_median = 1.2127612060049497
  electron_density_ratio_median = 0.8042820731527076
  hole_density_ratio_median = 0.5686951067977273

false_positive:
  count = 10
  source_ratio_median = 0.38273742388206267
  electric_field_ratio_median = 0.999657469779057
  electron_drive_over_sentaurus_field_median = 1.0076001051249843
  hole_drive_over_sentaurus_field_median = 1.0080500788421132
  electron_mobility_ratio_median = 0.9441656876596268
  hole_mobility_ratio_median = 1.2126450503873
  electron_density_ratio_median = 0.8175765435702607
  hole_density_ratio_median = 0.5692607352475433
```

Task 69 conclusion:

- The p99 support source deficit is not caused by local electric-field magnitude:
  Vela/Sentaurus `ElectricField` is `~0.99966` across overlap,
  false-negative, and false-positive nodes.
- It is also not caused by the high-field driving-force magnitude used by Vela:
  Vela electron/hole high-field drive is `~1.008` times the Sentaurus exported
  electric field on the same nodes.
- Mobility is not the leading explanation: electron mobility is only about
  `5.6%` lower, and hole mobility is about `21%` higher in Vela.
- Carrier density contributes but still does not fully explain the `~0.38`
  source ratio: electron density is `~0.81` and hole density is `~0.57` of
  Sentaurus on the same p99 support nodes.
- Therefore the next leading suspect is the `alpha * current` side of
  avalanche source assembly: Vela SG edge flux/current-density reconstruction,
  source-volume ownership, or a current-density unit/convention mismatch in the
  node-level comparison.

### Next Tasks After Task 69

1. Build an active-support SG edge/source ownership comparator:
   - select Vela SG avalanche edges adjacent to the p99 support nodes;
   - sum edge source contributions into the same p99 overlap/false-negative/
     false-positive node sets;
   - compare node VTK `AvalancheGeneration` against edge-source backprojection
     to rule out VTK/source-volume export mismatch.
2. Compare Vela SG flux-derived `alpha * |J|` against Sentaurus current-density
   fields on the same support nodes:
   - convert Sentaurus `eCurrentDensity`/`hCurrentDensity` to particle flux and
     infer weighted alpha from `ImpactIonization`;
   - compute analogous Vela electron/hole SG flux magnitudes from edge dumps;
   - determine whether the `~0.38` source ratio follows current/flux magnitude,
     alpha weighting, or source-volume assignment.
3. Only if SG current/source ownership is consistent, revisit carrier density
   reconstruction around the p99 line, because density ratios (`~0.81`
   electron, `~0.57` hole) may still affect the SG flux but are not sufficient
   to explain the full source deficit by themselves.

### Execution Note 2026-06-20: Task 70 Active-support SG Source Ownership

Task 70 executes the first Task 69 follow-up: backproject SG edge-source
contributions onto the same thresholded p99 support nodes and compare those
node source integrals with Vela VTK `AvalancheGeneration` and Sentaurus
`ImpactIonization`.

Added:

```text
scripts/diagnose_pn2d_bv_sg_source_ownership.py
```

Regression:

```text
tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.
test_pn2d_bv_sg_source_ownership_backprojects_edges_to_support_nodes
```

Real-data preparation:

```powershell
python scripts/diagnose_pn2d_bv_sg_avalanche_edges.py `
  --vtk build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression/vtk/impact_p95_bounded_retry_0002_-13.2V.vtk `
  --mesh build-release/reference_tcad/pn2d_sentaurus2018/vela/mesh.json `
  --doping-csv build-release/reference_tcad/pn2d_sentaurus2018/doping.csv `
  --bias -13.2 `
  --top 100000 `
  --out-dir build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/sg_python_bounded_retry_m13p2_full
```

Ownership command:

```powershell
python scripts/diagnose_pn2d_bv_sg_source_ownership.py `
  --support-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/thresholded_avalanche_support_m13p2/thresholded_avalanche_support_nodes.csv `
  --sg-edge-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/sg_python_bounded_retry_m13p2_full/sg_avalanche_edges.csv `
  --mesh build-release/reference_tcad/pn2d_sentaurus2018/vela/mesh.json `
  --out-dir build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/sg_source_ownership_m13p2
```

Full-device SG reconstruction summary:

```text
total_source_integral_reconstructed = 5.471857927158714e7
total_source_integral_vtk = 8.157142405465497e7
total_source_relative_error = -0.3291942625039343
node_source_log10_p95_error = 5.114015403485291
```

P99 support ownership summary:

```text
overlap:
  reconstructed_over_vtk = 0.6791355127730886
  vtk_over_sentaurus = 0.38448975411833053
  reconstructed_over_sentaurus = 0.26112064631915116
  electron source integral = 7.093639767943259e5
  hole source integral = 2.5863542628044728e5

false_negative:
  reconstructed_over_vtk = 0.6790990738825788
  vtk_over_sentaurus = 0.37729255530521927
  reconstructed_over_sentaurus = 0.2562190248905661

false_positive:
  reconstructed_over_vtk = 0.6785561211256687
  vtk_over_sentaurus = 0.3833841346511587
  reconstructed_over_sentaurus = 0.26014765131001133
```

Task 70 conclusion:

- The Vela VTK source itself remains about `0.38` of Sentaurus on the p99
  support, consistent with Tasks 68 and 69.
- The Python SG edge reconstruction is lower than the Vela VTK source by another
  `~32%` (`reconstructed_over_vtk ~= 0.679`) on both p99 support nodes and the
  full device.
- Therefore the current Python SG diagnostic is not yet an authoritative proxy
  for the C++ assembled source on this bounded-retry state. Do not tune physics
  from Python SG flux ratios until the same state is dumped from C++.
- The next source of truth should be Vela's C++ `sweep.diagnostics.sg_avalanche_edges`
  output for the bounded-retry `-13.2 V` state, then repeat the ownership
  backprojection against that C++ edge dump.

### Next Tasks After Task 70

1. Generate a bounded-retry C++ SG edge dump:
   - copy the Task 68/69 bounded-retry config into a probe directory;
   - enable `sweep.diagnostics.sg_avalanche_edges.enabled`;
   - write `sg_avalanche_edges.csv` for the same `-13.2 V` final state.
2. Run `diagnose_pn2d_bv_sg_source_ownership.py` against the C++ edge dump:
   - if C++ `reconstructed_over_vtk ~= 1`, continue with C++ SG flux/current
     versus Sentaurus current-density comparison;
   - if C++ also gives `~0.679`, inspect VTK avalanche source export volume or
     node-source accumulation semantics.
3. Only after C++ edge ownership is reconciled with Vela VTK, compare
   Sentaurus-inferred `alpha * |J|` with Vela edge flux/source terms on the p99
   support line.

### Execution Note 2026-06-20: Task 71 Bounded-retry C++ SG Edge Dump Ownership

Task 71 executes the first Task 70 follow-up: generate a C++ SG avalanche edge
dump for the same bounded-retry `-13.2 V` state used by Tasks 68-70 and run the
same ownership backprojection against that C++ edge dump.

Generated probe config:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe/simulation.json
```

Only the output paths and diagnostics were changed from the Task 68/69 bounded
retry config:

```json
"sweep": {
  "diagnostics": {
    "sg_avalanche_edges": {
      "enabled": true,
      "csv_file": ".../impact_p95_guard_bounded_retry_regression_cpp_sg_probe/sg_avalanche_edges.csv"
    }
  }
}
```

Run:

```powershell
cmake --build build-release --target vela_example_runner --parallel
build-release/vela_example_runner.exe `
  --config build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe/simulation.json
```

Result:

```text
converged = true
points = 3
C++ sg_avalanche_edges rows at -13.2 V = 3830
```

The ownership comparator was extended to support C++ edge dump columns
(`edge_source_integral`) and a `--bias` filter. Regression:

```text
tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.
test_pn2d_bv_sg_source_ownership_filters_cxx_edge_dump_bias
```

Ownership command:

```powershell
python scripts/diagnose_pn2d_bv_sg_source_ownership.py `
  --support-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/thresholded_avalanche_support_m13p2/thresholded_avalanche_support_nodes.csv `
  --sg-edge-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe/sg_avalanche_edges.csv `
  --mesh build-release/reference_tcad/pn2d_sentaurus2018/vela/mesh.json `
  --bias -13.2 `
  --out-dir build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/sg_source_ownership_cpp_m13p2
```

C++ ownership summary:

```text
overlap:
  reconstructed_over_vtk = 1.0000014451520778
  vtk_over_sentaurus = 0.38448975411833053
  reconstructed_over_sentaurus = 0.3844903097644976
  reconstructed_electron_source_integral = 1.0445762392524544e6
  reconstructed_hole_source_integral = 3.8076639675331186e5

false_negative:
  reconstructed_over_vtk = 0.9999996546937229
  vtk_over_sentaurus = 0.37729255530521927
  reconstructed_over_sentaurus = 0.3772924250237316

false_positive:
  reconstructed_over_vtk = 0.9999996159868827
  vtk_over_sentaurus = 0.3833841346511587
  reconstructed_over_sentaurus = 0.38338398742662205
```

Task 71 conclusion:

- C++ SG edge-source ownership exactly matches the Vela VTK
  `AvalancheGeneration` node source integral on the p99 support sets.
- The previous Python SG reconstruction deficit (`~0.679` of VTK) is a Python
  diagnostic reconstruction difference, not a C++ source ownership or VTK export
  inconsistency.
- The real remaining physical discrepancy is now isolated to C++ edge source
  magnitude versus Sentaurus: C++/VTK source is only `~0.38` of Sentaurus on the
  same p99 support nodes.
- Since Task 69 already showed local field and Vela high-field drive are
  `~1.0` of Sentaurus, the next discriminator is the current/flux side of
  `alpha * |J|`: C++ SG electron/hole flux proxy and alpha weighting versus
  Sentaurus `eCurrentDensity`, `hCurrentDensity`, and inferred weighted alpha.

### Next Tasks After Task 71

1. Build a C++ edge-flux/current comparator on the p99 support line:
   - select C++ SG edges incident to the overlap/false-negative/false-positive
     p99 nodes;
   - aggregate `electron_flux_proxy`, `hole_flux_proxy`,
     `electron_alpha_m_inv`, `hole_alpha_m_inv`, and source integrals;
   - compare the resulting node/edge-scale `alpha * flux` components with
     Sentaurus `ImpactIonization`, `eCurrentDensity`, and `hCurrentDensity`.
2. Infer Sentaurus weighted alpha from the same nodes:
   - convert current density to particle flux using `|J|/q`;
   - compute `ImpactIonization / (|Jn|/q + |Jp|/q)`;
   - compare against C++ electron/hole alpha and flux-proxy weighted alpha.
3. If weighted alpha matches but flux does not, focus on SG flux/current
   approximation and carrier-density reconstruction. If flux matches but alpha
   does not, focus on van Overstraeten coefficient/temperature/field convention.

### Execution Note 2026-06-20: Task 72 C++ SG Flux/Current Comparator

Task 72 added a C++ edge-flux/current comparator for the same `-13.2 V` p99
thresholded avalanche support nodes used in Tasks 68-71.

New tool:

```powershell
python scripts/diagnose_pn2d_bv_cxx_edge_flux_current.py `
  --support-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/thresholded_avalanche_support_m13p2/thresholded_avalanche_support_nodes.csv `
  --sg-edge-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe/sg_avalanche_edges.csv `
  --sentaurus-dir build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias/sentaurus_-13.2v `
  --bias -13.2 `
  --out-dir build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/cxx_edge_flux_current_m13p2
```

Output:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/cxx_edge_flux_current_m13p2/cxx_edge_flux_current_nodes.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/cxx_edge_flux_current_m13p2/cxx_edge_flux_current_summary.json
```

The comparator converts Sentaurus `eCurrentDensity` and `hCurrentDensity` from
`A/cm2` to particle flux in `m^-2 s^-1`, converts `ImpactIonization` from
`cm^-3 s^-1` to `m^-3 s^-1`, and infers:

```text
sentaurus_weighted_alpha_m_inv =
  ImpactIonization_m3_s / (|Jn|/q + |Jp|/q)

cxx_weighted_alpha_m_inv =
  sum(alpha_n * |electron_flux_proxy| + alpha_p * |hole_flux_proxy|)
  / sum(|electron_flux_proxy| + |hole_flux_proxy|)
```

Regression coverage:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_cxx_edge_flux_current_comparator_infers_weighted_alpha
```

Real `-13.2 V` p99-support result:

| support class | count | C++ flux / Sentaurus flux median | C++ weighted alpha / Sentaurus weighted alpha median |
| --- | ---: | ---: | ---: |
| overlap | 10 | `1.5116653936098987` | `1.0135069350112558` |
| false_negative | 10 | `1.4960966677953818` | `1.011569500894982` |
| false_positive | 10 | `1.5155457302553186` | `1.0096777726114192` |

Interpretation:

- The C++ SG edge flux proxy is not too small; it is about `1.50x` the
  Sentaurus current-density-derived particle flux on the same p99 nodes.
- The C++ alpha weighting matches the Sentaurus inferred weighted alpha to
  within about `1-2%`.
- This rules out the local van Overstraeten alpha, high-field drive, mobility,
  and SG flux magnitude as the direct cause of the Task 71 `~0.38` C++/VTK
  source versus Sentaurus source ratio.
- The remaining discrepancy is therefore concentrated in the mapping from
  edge-level `alpha * flux` to nodal source density/integral: edge-area proxy,
  node/control-volume normalization, Sentaurus box-volume convention, or
  endpoint/source ownership near the abrupt junction.

### Next Tasks After Task 72

1. Add a source-geometry comparator for the same p99 support nodes:
   - join C++ edge `edge_area_proxy_m2`, endpoint source integrals, and Vela
     mesh node control volumes;
   - compute the effective factor that maps `alpha * flux` to node source
     density;
   - compare this factor against the factor implied by Sentaurus
     `ImpactIonization / (alpha_weighted * particle_flux)`.
2. Check whether the `~0.38` source-density ratio follows:
   - C++ edge-area proxy divided by Vela node control volume;
   - endpoint half-source ownership;
   - junction-adjacent box-volume truncation or compensated-node ownership.
3. If the geometry factor explains the deficit, prototype a bounded diagnostic
   switch for the source-volume policy and rerun only the `-13.2 V` p99-support
   comparison before touching the production BV acceptance path.

### Execution Note 2026-06-20: Task 73 Source-Geometry Factor Comparator

Task 73 tested the Task 72 geometry hypothesis directly. It added a comparator
that joins the C++ SG edge dump, Vela mesh control volumes, Sentaurus
`ImpactIonization`, and the p99 support mask.

New tool:

```powershell
python scripts/diagnose_pn2d_bv_source_geometry.py `
  --support-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/thresholded_avalanche_support_m13p2/thresholded_avalanche_support_nodes.csv `
  --sg-edge-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe/sg_avalanche_edges.csv `
  --sentaurus-dir build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias/sentaurus_-13.2v `
  --mesh build-release/reference_tcad/pn2d_sentaurus2018/vela/mesh.json `
  --bias -13.2 `
  --out-dir build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/source_geometry_m13p2
```

Output:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/source_geometry_m13p2/source_geometry_nodes.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/source_geometry_m13p2/source_geometry_summary.json
```

The comparator decomposes the source ratio as:

```text
C++ source / Sentaurus source
  = (C++ endpoint edge-area sum / Vela node volume)
    * (C++ area-weighted source-density proxy / Sentaurus generation)
```

It also reports an active-edge-only density proxy, where active edges are
incident edges with `alpha * flux` above a relative threshold of the node-local
maximum. This separates transverse near-zero source faces from high-current
faces.

Regression coverage:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_source_geometry_decomposes_area_volume_factor
```

Real `-13.2 V` p99-support result:

| support class | source ratio summed | endpoint area / node volume median | active endpoint area fraction median | active source-density / Sentaurus generation median | total source-density / Sentaurus generation median |
| --- | ---: | ---: | ---: | ---: | ---: |
| overlap | `0.3844903097644976` | `1.0` | `0.5` | `0.7660388358079806` | `0.3830194179039903` |
| false_negative | `0.3772924250237316` | `1.0` | `0.5` | `0.7566957020338947` | `0.37834785101694735` |
| false_positive | `0.38338398742662205` | `1.0` | `0.5` | `0.7654757797985718` | `0.3827378898992859` |

Concrete node evidence:

- Node `351` has four incident interior edges at `-13.2 V`.
- Two edges have high area-density `alpha * flux` around `2.49e21 m^-3 s^-1`.
- Two perpendicular/transverse edges have effectively zero `alpha * flux`
  (`~1e-284`).
- The endpoint edge-area sum equals the Vela node volume, so the box/control
  volume geometry is not missing a factor.
- Because only half the endpoint area carries the high-current source, Vela's
  full face-area average is about half the active-edge density. The active-edge
  density is itself only about `0.76x` Sentaurus, giving the observed
  `0.38x` total source.

Interpretation:

- The Task 72 geometry hypothesis is partially falsified: Vela is not losing
  source through a wrong endpoint-area/node-volume normalization. The summed
  endpoint area equals the node volume on the p99 support.
- The dominant discrepancy is a discretization/convention mismatch in the
  nodal avalanche density definition. Vela's finite-volume edge-source density
  averages high-current faces together with transverse near-zero-current faces;
  Sentaurus's exported `ImpactIonization` at these nodes behaves closer to a
  directional or nodal current-density magnitude, not a full face-area average.
- The remaining `~0.76` active-edge density deficit should be investigated
  separately after the `0.5` active-area dilution convention is understood.

### Next Tasks After Task 73

1. Add an edge-direction/source-density classifier:
   - classify incident p99-support edges by orientation relative to the junction
     and electric-field/current direction;
   - report high-source edge density, transverse zero-source area fraction, and
     reconstructed source under alternative policies:
     `full_face_average`, `active_edge_average`, and `active_edge_sum`.
2. Compare these reconstructed policies directly against Sentaurus
   `ImpactIonization`:
   - if `active_edge_average` is near Sentaurus, prototype a diagnostic
     Sentaurus-default source-density convention before changing solver behavior;
   - if `active_edge_sum` is nearer, inspect whether Sentaurus nodal current
     density represents a directional face-pair sum rather than an averaged box
     source.
3. Before changing production BV behavior, look for Sentaurus manual/par-file
   switches that affect avalanche volume/source post-processing, especially
   options related to element/box volume avalanche, flat-element exclusion, and
   current-density based avalanche generation.

### Execution Note 2026-06-20: Task 74 Edge-Direction Source Policy Classifier

Task 74 added an edge-direction/source-density classifier to compare the
candidate nodal avalanche-density policies from Task 73:

- `full_face_average`: current Vela finite-volume behavior, total endpoint
  source divided by Vela node control volume;
- `active_edge_average`: source density averaged only over active incident
  edges whose `alpha * flux` is above the node-local relative threshold;
- `active_edge_density_sum`: sum of active edge source densities, useful to
  test whether Sentaurus behaves like a directional face-pair sum.

New tool:

```powershell
python scripts/diagnose_pn2d_bv_edge_direction_source_policy.py `
  --support-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/thresholded_avalanche_support_m13p2/thresholded_avalanche_support_nodes.csv `
  --sg-edge-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe/sg_avalanche_edges.csv `
  --sentaurus-dir build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias/sentaurus_-13.2v `
  --mesh build-release/reference_tcad/pn2d_sentaurus2018/vela/mesh.json `
  --bias -13.2 `
  --out-dir build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/edge_direction_source_policy_m13p2
```

Output:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/edge_direction_source_policy_m13p2/edge_direction_source_policy_nodes.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/edge_direction_source_policy_m13p2/edge_direction_source_policy_summary.json
```

Regression coverage:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_edge_direction_policy_compares_active_reconstructions
```

Real `-13.2 V` p99-support result:

| support class | full face / Sentaurus median | active edge average / Sentaurus median | active edge density sum / Sentaurus median | active area fraction median | closest policy |
| --- | ---: | ---: | ---: | ---: | --- |
| overlap | `0.3830194179039903` | `0.7660388358079806` | `1.5320776716159612` | `0.5` | `active_edge_average` |
| false_negative | `0.37834785101694735` | `0.7566957020338947` | `1.5133914040677894` | `0.5` | `active_edge_average` |
| false_positive | `0.3827378898992859` | `0.7654757797985718` | `1.5309515595971437` | `0.5` | `active_edge_average` |

Direction classification:

- The junction-normal axis is `x` for the pn2d fixture.
- The active high-source edges are `x`-aligned.
- The transverse `y`-aligned edges contribute half the endpoint area but no
  meaningful source (`junction_tangent_active_area_fraction = 0`).

Interpretation:

- Sentaurus exported `ImpactIonization` is not consistent with Vela's current
  full face-area average on the p99-support nodes.
- Sentaurus also does not look like a direct sum of the two active face
  densities; that policy overshoots to about `1.51-1.53x`.
- The closest current diagnostic policy is `active_edge_average`, still low by
  about `23-24%`. This is now the leading parity target to investigate before
  changing the production source assembly.
- The remaining gap is no longer a broad alpha, mobility, field, SG flux, or
  control-volume issue. It is specifically the mapping between Sentaurus's
  nodal current-density avalanche export and Vela's edge-to-node source-density
  reconstruction along the active junction-normal edges.

### Next Tasks After Task 74

1. Compare Sentaurus nodal current density against active-edge C++ flux on the
   same p99 nodes with orientation awareness:
   - use only `x`-aligned active edges;
   - compute active-edge average particle flux and alpha-weighted generation;
   - compare to Sentaurus `eCurrentDensity`, `hCurrentDensity`, and
     `TotalCurrentDensity`.
2. Add a diagnostic Sentaurus-default source-density policy prototype that is
   reporting-only at first:
   - `full_face_average` keeps current Vela behavior;
   - `active_edge_average` emits a parallel avalanche field/current estimate;
   - do not change the nonlinear source until the diagnostic field comparison
     proves the policy tracks Sentaurus over more than one bias point.
3. Check Sentaurus documentation or exported parameters for whether the
   `ImpactIonization` field is a box-averaged generation density, an element
   value interpolated to nodes, or a current-direction nodal value. This decides
   whether the Vela fix belongs in source assembly, VTK/export diagnostics, or
   the Sentaurus comparison post-processing.

### Execution Note 2026-06-20: Task 75 Active-Edge Current-Density Comparator

Task 75 compared Sentaurus nodal current-density fields against Vela's
junction-normal active SG edge fluxes. This continues the Task 74 result that
`active_edge_average` is the closest source-density policy, but still low by
about `23-24%`.

New tool:

```powershell
python scripts/diagnose_pn2d_bv_active_edge_current_density.py `
  --support-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/thresholded_avalanche_support_m13p2/thresholded_avalanche_support_nodes.csv `
  --sg-edge-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe/sg_avalanche_edges.csv `
  --sentaurus-dir build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias/sentaurus_-13.2v `
  --bias -13.2 `
  --out-dir build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/active_edge_current_density_m13p2
```

Output:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/active_edge_current_density_m13p2/active_edge_current_density_nodes.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/active_edge_current_density_m13p2/active_edge_current_density_summary.json
```

The comparator uses only `x`-aligned active edges and computes area-weighted
active-edge averages for electron flux, hole flux, total particle flux,
weighted alpha, and alpha-weighted generation. It converts Sentaurus
`eCurrentDensity`, `hCurrentDensity`, and `TotalCurrentDensity` from `A/cm2` to
particle-flux-equivalent `m^-2 s^-1`.

Regression coverage:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_active_edge_current_density_matches_sentaurus_currents
```

Real `-13.2 V` p99-support result:

| support class | electron flux / Sentaurus e-current median | hole flux / Sentaurus h-current median | particle flux / Sentaurus e+h median | generation / Sentaurus median | weighted alpha / Sentaurus median |
| --- | ---: | ---: | ---: | ---: | ---: |
| overlap | `0.7895051205521759` | `0.6990414589372504` | `0.7417144309689099` | `0.7660388358079806` | `1.0335011569183161` |
| false_negative | `0.7773341669306977` | `0.7012797040278108` | `0.7341809638074493` | `0.7566957020338947` | `1.0300480563566938` |
| false_positive | `0.7881212636047819` | `0.7015909044148467` | `0.7416939160101677` | `0.7654757797985718` | `1.0322478366616914` |

`TotalCurrentDensity` check:

- `active_particle_flux_over_sentaurus_total_current_equiv` matches
  `active_particle_flux_over_sentaurus_eh_flux` to roundoff in all three
  classes.
- Therefore the Sentaurus `TotalCurrentDensity` export is consistent with the
  separate electron/hole current-density magnitudes for this diagnostic; it is
  not the source of the remaining gap.

Interpretation:

- The remaining `~0.76x` active-edge generation ratio is not caused by the
  avalanche coefficient. Vela's active-edge weighted alpha is actually slightly
  high (`~1.03x`).
- The gap is dominated by active-edge current/flux magnitude: particle flux is
  only `~0.734-0.742x` of the Sentaurus current-density-derived particle flux.
- Electron flux is low by `~21-22%`; hole flux is low by about `30%`. The hole
  side remains the larger active-edge current-density mismatch.
- The current leading root-cause branch is now active-edge SG flux /
  edge-to-node current-density mapping, not source-volume normalization,
  high-field alpha, mobility, electric field, terminal current extraction, or
  Sentaurus `TotalCurrentDensity` post-processing.

### Next Tasks After Task 75

1. Decompose active-edge SG flux into the quantities that drive the Bernoulli
   current:
   - active-edge endpoint carrier densities;
   - quasi-Fermi/electrostatic potential drops;
   - Bernoulli factors and mobility;
   - electron versus hole contributions separately.
2. Reconstruct the same active-edge flux using Sentaurus node fields on the
   Vela mesh:
   - if Sentaurus-state active-edge SG flux matches Sentaurus current density,
     Vela's state branch/densities remain the source of the flux deficit;
   - if Sentaurus-state SG flux is also about `0.74x`, the mismatch is in the
     edge-to-node current-density convention, interpolation, or Sentaurus export
     definition.
3. Only after the active-edge flux convention is explained, prototype a
   reporting-only `active_edge_average` avalanche diagnostic field. Do not
   change nonlinear source assembly until the flux convention is reconciled over
   multiple biases.

### Execution Note 2026-06-20: Task 76 Active-Edge SG Flux Factor Decomposition

Task 76 decomposed the Task 75 active-edge current deficit into SG/Bernoulli
flux factors for Vela state versus Sentaurus state on the same `x`-aligned
active p99-support edges.

New tool:

```powershell
python scripts/diagnose_pn2d_bv_active_edge_flux_factors.py `
  --support-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/thresholded_avalanche_support_m13p2/thresholded_avalanche_support_nodes.csv `
  --sg-edge-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe/sg_avalanche_edges.csv `
  --mesh build-release/reference_tcad/pn2d_sentaurus2018/vela/mesh.json `
  --doping-csv build-release/reference_tcad/pn2d_sentaurus2018/vela/doping.csv `
  --sentaurus-dir build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias/sentaurus_-13.2v `
  --vela-vtk-root build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe/vtk `
  --bias -13.2 `
  --out-dir build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/active_edge_flux_factors_m13p2
```

Output:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/active_edge_flux_factors_m13p2/active_edge_flux_factors_nodes.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/active_edge_flux_factors_m13p2/active_edge_flux_factors_summary.json
```

Regression coverage:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_active_edge_flux_factors_reconstructs_matching_states
```

Real `-13.2 V` p99-support result:

| support class | C++ flux / Vela density-SG median | Vela density-SG / Sentaurus-state density-SG median | Sentaurus-state density-SG / Sentaurus current median | Vela QF-model / Sentaurus-state QF-model median |
| --- | ---: | ---: | ---: | ---: |
| overlap | `1.0056086845685668` | `0.7312303768731796` | `1.0086692627539615` | `0.7312233313074645` |
| false_negative | `1.0056254524835082` | `0.723801725265691` | `1.0086692627533154` | `0.7238110745631495` |
| false_positive | `1.0056341400376154` | `0.7320182381793234` | `1.0075487331485857` | `0.7320370845719382` |

Factor medians on overlap p99 nodes:

| factor | Vela / Sentaurus median |
| --- | ---: |
| electric field magnitude | `0.9996572373228313` |
| electron quasi-Fermi field magnitude | `1.0009451843475397` |
| hole quasi-Fermi field magnitude | `1.0003426627173937` |
| electron mobility | `0.9424148302400709` |
| hole mobility | `1.2101927812318916` |
| electron endpoint density geometric mean | `0.8192917935861858` |
| hole endpoint density geometric mean | `0.5668364003170232` |
| electron density-form flux | `0.7782580758839022` |
| hole density-form flux | `0.6892373761967474` |

Interpretation:

- The C++ SG edge dump is internally consistent with Vela's density-form SG
  reconstruction (`~1.006x`).
- Reconstructing SG flux from Sentaurus node fields on the same Vela active
  edges matches Sentaurus `eCurrentDensity + hCurrentDensity` to about `1%`.
  Therefore the remaining mismatch is not an edge-to-node current-density
  convention or a Sentaurus current-density export problem.
- Vela's active-edge SG flux is `~0.72-0.73x` Sentaurus-state SG flux, matching
  the Task 75 active-edge current deficit.
- The driving fields are already matched to `~0.1%`; alpha was already slightly
  high in Task 75. The remaining deficit is carrier-state driven:
  electron density is `~0.82x` and electron mobility is `~0.94x`, giving
  electron flux `~0.78x`; hole density is only `~0.57x`, partly offset by
  hole mobility `~1.21x`, giving hole flux `~0.69x`.
- This moves the leading root cause from SG flux/current convention to the
  active-edge carrier-density state, especially hole density.

### Next Tasks After Task 76

1. Decompose active-edge endpoint density ratios into:
   - `ni_eff` or inferred intrinsic-density ratio;
   - Boltzmann exponent offsets from `psi - phin` for electrons and
     `phip - psi` for holes;
   - endpoint asymmetry between p-side and n-side active nodes.
2. Reconstruct active-edge densities under mixed states:
   - Vela `ni_eff` with Sentaurus `psi/qf`;
   - Sentaurus inferred `ni_eff` with Vela `psi/qf`;
   - Vela/Sentaurus mobility swapped after density is fixed.
3. If mixed-state reconstruction shows density mismatch comes from quasi-Fermi
   offsets, return to continuity equation/state-branch debugging. If it comes
   from `ni_eff`/BGN, revisit Task 63+ BGN validation and Sentaurus parameter
   defaults specifically on the active p99 edges.

### Execution Note 2026-06-20: Task 77 Active-Edge Density Factor Decomposition

Task 77 decomposed the Task 76 active-edge carrier-density mismatch into
inferred intrinsic-density and Boltzmann exponent factors.

New tool:

```powershell
python scripts/diagnose_pn2d_bv_active_edge_density_factors.py `
  --support-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/thresholded_avalanche_support_m13p2/thresholded_avalanche_support_nodes.csv `
  --sg-edge-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe/sg_avalanche_edges.csv `
  --mesh build-release/reference_tcad/pn2d_sentaurus2018/vela/mesh.json `
  --sentaurus-dir build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias/sentaurus_-13.2v `
  --vela-vtk-root build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe/vtk `
  --bias -13.2 `
  --out-dir build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/active_edge_density_factors_m13p2
```

Output:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/active_edge_density_factors_m13p2/active_edge_density_factors_nodes.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/active_edge_density_factors_m13p2/active_edge_density_factors_summary.json
```

Regression coverage:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_active_edge_density_factors_split_ni_and_exponent
```

Real `-13.2 V` p99-support result:

| support class | e-density ratio | e-inferred-ni ratio | e-Boltzmann ratio | e `psi-phin` delta | h-density ratio | h-inferred-ni ratio | h-Boltzmann ratio | h `phip-psi` delta |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| overlap | `0.8192917935861858` | `0.9999641950239109` | `0.8193710650826915` | `-0.005143811533930065 V` | `0.5668364003170232` | `1.0000097048946204` | `0.5668314451200728` | `-0.014672433882886499 V` |
| false_negative | `0.8054381516182242` | `1.0000728843109974` | `0.8053506558227881` | `-0.00558881148324053 V` | `0.5687183067828594` | `1.0000310973165982` | `0.5687277598582157` | `-0.014586183884296888 V` |
| false_positive | `0.8182056625974807` | `0.9999809449810138` | `0.8181995356636942` | `-0.005180173608212046 V` | `0.5697687278303936` | `0.9999399464899278` | `0.5698322304224902` | `-0.0145358666986318 V` |

Endpoint asymmetry:

- Electron left/right endpoint density ratios are both about `0.81-0.83`, so
  the electron deficit is not a one-sided endpoint artifact.
- Hole left/right endpoint density ratios are both about `0.56-0.57`, so the
  larger hole deficit is also distributed across the active edge endpoints.

Interpretation:

- The active-edge density mismatch is not caused by `ni_eff`, BGN, or intrinsic
  density defaults. Inferred `ni` ratios are `~1.0` for both carriers and all
  p99 support classes.
- The density mismatch is almost exactly the Boltzmann exponent mismatch:
  Vela's active-edge electron `psi - phin` is about `5 mV` below Sentaurus,
  and Vela's active-edge hole `phip - psi` is about `14.5-14.7 mV` below
  Sentaurus.
- Therefore the current leading root cause is the coupled carrier continuity
  state/quasi-Fermi branch on the active junction-normal edges, especially the
  hole quasi-Fermi relation. BGN/`ni_eff`, alpha, mobility, field magnitude, SG
  current discretization, source-volume normalization, and terminal-current
  extraction are no longer the leading branches for this `-13.2 V` p99 source
  deficit.

### Next Tasks After Task 77

1. Trace the active-edge quasi-Fermi exponent offsets back to continuity
   residual balance:
   - compare electron/hole continuity residual terms on active p99 edges;
   - separate drift/diffusion edge flux imbalance, recombination, and avalanche
     source contribution;
   - identify whether the `phip-psi` offset is set by hole continuity, contact
     boundary propagation, or avalanche feedback.
2. Run a mixed-state density/flux replay:
   - Vela `psi` with Sentaurus `phin/phip`;
   - Sentaurus `psi` with Vela `phin/phip`;
   - Vela state with only hole quasi-Fermi shifted by the measured
     `+14.6 mV` on active edges;
   - report predicted active-edge flux/source and terminal current direction.
3. If the mixed-state replay shows a small quasi-Fermi correction closes the
   active-edge source gap, prototype a diagnostic continuation/state-selection
   gate that monitors active-edge `psi-phin` and `phip-psi` against Sentaurus
   references before changing physical models.

### Execution Note 2026-06-20: Task 78 Active-Edge Mixed-State Replay

Task 78 replayed active-edge SG flux/source under mixed Vela/Sentaurus states
to test whether the Task 77 quasi-Fermi exponent offsets are sufficient to
explain the active-edge source gap.

New tool:

```powershell
python scripts/diagnose_pn2d_bv_active_edge_mixed_state_replay.py `
  --support-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/thresholded_avalanche_support_m13p2/thresholded_avalanche_support_nodes.csv `
  --sg-edge-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe/sg_avalanche_edges.csv `
  --mesh build-release/reference_tcad/pn2d_sentaurus2018/vela/mesh.json `
  --sentaurus-dir build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias/sentaurus_-13.2v `
  --vela-vtk-root build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe/vtk `
  --bias -13.2 `
  --electron-qf-shift-v -0.005143811533930065 `
  --hole-qf-shift-v 0.014672433882886499 `
  --out-dir build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/active_edge_mixed_state_replay_m13p2
```

Output:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/active_edge_mixed_state_replay_m13p2/active_edge_mixed_state_replay_nodes.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/active_edge_mixed_state_replay_m13p2/active_edge_mixed_state_replay_summary.json
```

Regression coverage:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_active_edge_mixed_state_replay_applies_qf_shift
```

Replay variants:

- `vela_baseline`: current Vela state.
- `sentaurus_baseline`: Sentaurus exported state, used as denominator.
- `vela_psi_sentaurus_qf`: Vela `psi/ni/mobility` with Sentaurus
  `phin/phip`.
- `sentaurus_psi_vela_qf`: Sentaurus `psi/ni` with Vela `phin/phip` and Vela
  mobility.
- `vela_qf_shift`: Vela state with `phin` shifted by `-5.1438 mV` and `phip`
  shifted by `+14.6724 mV`.
- `vela_qf_shift_sentaurus_mobility`: same qf-shifted Vela density state, but
  with Sentaurus mobility.

Real `-13.2 V` overlap p99-support result:

| variant | e-density | h-density | e-flux | h-flux | particle flux | generation |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `vela_baseline` | `0.8192917935861858` | `0.5668364003170232` | `0.7782580758839022` | `0.6892373761967474` | `0.7312303768731795` | `0.754521916587773` |
| `vela_psi_sentaurus_qf` | `1.1831899339634622` | `0.8452731845493741` | `1.1176384632600957` | `1.025124652510192` | `1.0673228431468016` | `1.0904683750514947` |
| `sentaurus_psi_vela_qf` | `0.6925612141426314` | `0.673326247909682` | `0.6565290923364693` | `0.8169047487394019` | `0.7437528326245539` | `0.7055376962869636` |
| `vela_qf_shift` | `0.9996565945847643` | `0.9998705388332872` | `0.949589418491291` | `1.2157796258963827` | `1.0968112406862471` | `1.031142804388406` |
| `vela_qf_shift_sentaurus_mobility` | `0.9996565945847643` | `0.9998705388332872` | `1.0073714051607303` | `1.0045969025392862` | `1.0078816548407237` | `1.009922844770538` |

Class-level behavior:

- `vela_qf_shift_sentaurus_mobility` closes the active-edge particle flux and
  generation to within about `1%` for all p99 support classes:
  - overlap particle flux median `1.0078816548407237`, generation median
    `1.009922844770538`;
  - false negative particle flux median `0.9980277985897577`, generation
    median `0.9979857955800552`;
  - false positive particle flux median `1.0091925951580474`, generation
    median `1.0104192611597498`.

Interpretation:

- The Task 77 measured qf exponent offsets are sufficient to close the carrier
  density gap on active edges. Applying `phin -= 5.14 mV` and
  `phip += 14.67 mV` brings electron and hole active-edge densities to
  `~1.0x` Sentaurus.
- After density is fixed, Vela mobility becomes the remaining active-edge flux
  mismatch: electron flux stays low with Vela electron mobility, while hole
  flux overshoots with Vela hole mobility. Swapping in Sentaurus mobility closes
  both carrier fluxes and alpha-weighted generation to about `1%`.
- Therefore the high-bias p99 source deficit is now decomposed into:
  1. quasi-Fermi exponent/state offset as the carrier-density root cause;
  2. mobility-model mismatch as the residual flux scaling after density is
     corrected.
- This is not a source-volume, SG-current convention, `ni_eff`/BGN, alpha,
  electric-field, terminal-current, or Sentaurus export issue.

### Next Tasks After Task 78

1. Trace the quasi-Fermi exponent offset to continuity equation balance:
   - active-edge hole and electron flux divergence;
   - SRH recombination;
   - avalanche generation;
   - contact/boundary contribution along the path from contacts to the junction.
2. Compare Vela and Sentaurus mobility model inputs on active p99 edges:
   - low-field mobility;
   - high-field mobility reduction;
   - carrier-density and field dependence;
   - whether mobility should use Vela or Sentaurus driving force on the
     qf-corrected state.
3. Prototype only diagnostic replays first:
   - qf-shifted density replay;
   - qf-shifted plus Sentaurus-mobility replay;
   - report predicted active-edge source and terminal-current direction across
     `-12.85`, `-12.9078`, and `-13.2 V` before changing nonlinear solver
     behavior.

### Execution Note 2026-06-20: Task 79 Active-Support Continuity Balance

Task 79 added a support-wide continuity-balance diagnostic for the thresholded
`-13.2 V` avalanche p99 nodes. Unlike the earlier single-edge
`continuity_feedback` probe, this report summarizes all non-inactive support
nodes and records:

- active x-edge ownership count and IDs;
- Vela/Sentaurus quasi-Fermi exponent deltas;
- node electron/hole density ratios;
- Vela signed SG transport integrals;
- Vela edge-integrated avalanche, VTK avalanche, and SRH integrals;
- Sentaurus ImpactIonization and SRH integrals on nearest nodes;
- residual estimates `transport + SRH - edge_avalanche`.

Implementation:

```powershell
python scripts/diagnose_pn2d_bv_active_support_continuity_balance.py `
  --support-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/thresholded_avalanche_support_m13p2/thresholded_avalanche_support_nodes.csv `
  --sg-edge-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe/sg_avalanche_edges.csv `
  --mesh build-release/reference_tcad/pn2d_sentaurus2018/vela/mesh.json `
  --doping-csv build-release/reference_tcad/pn2d_sentaurus2018/vela/doping.csv `
  --sentaurus-dir build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias/sentaurus_-13.2v `
  --vela-vtk-root build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe/vtk `
  --out-dir build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/active_support_continuity_balance_m13p2 `
  --bias -13.2
```

Outputs:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/active_support_continuity_balance_m13p2/active_support_continuity_balance_nodes.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/active_support_continuity_balance_m13p2/active_support_continuity_balance_summary.json
```

Key `-13.2 V` median results:

| support class | count | d(psi-phin) | d(phip-psi) | n ratio | p ratio | edge G ratio | VTK G ratio | SRH ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| overlap | 10 | `-5.17 mV` | `-14.71 mV` | `0.819` | `0.566` | `0.381` | `0.383` | `1.000006` |
| false_negative | 10 | `-5.64 mV` | `-14.59 mV` | `0.804` | `0.569` | `0.376` | `0.378` | `1.000006` |
| false_positive | 10 | `-5.21 mV` | `-14.56 mV` | `0.818` | `0.569` | `0.380` | `0.383` | `1.000006` |
| all | 30 | `-5.33 mV` | `-14.62 mV` | `0.814` | `0.568` | `0.379` | `0.381` | `1.000006` |

Interpretation:

- SRH is not the source of the active-support density mismatch: Vela and
  Sentaurus SRH integrals are equal to about `6e-6` relative error on these
  nodes.
- The quasi-Fermi exponent offsets from Task 77 persist support-wide, not only
  on a selected edge.
- Full node-integrated Vela avalanche is about `0.38x` Sentaurus, matching the
  earlier full-face dilution result. The active-edge average still explains the
  closer `~0.76x` local source density, so the remaining physical gap is the
  carrier state and mobility, not SRH or nodal volume normalization.
- Vela residual estimates are mostly negative relative to edge avalanche on
  overlap nodes (`electron ~-0.395x`, `hole ~-0.510x` median). This is a
  diagnostic estimate from exported fields, not a direct nonlinear residual
  assertion, but it argues for inspecting carrier-continuity boundary/current
  balance rather than retuning SRH.
- Strict same-bias multi-point support comparison is currently blocked by
  available Sentaurus exports: the local export directory contains
  `sentaurus_-12.8v` and `sentaurus_-13.2v`, while the current Vela SG probe VTK
  sequence contains `-12.85`, `-12.9078`, and `-13.2 V`. Before drawing a
  transition trend, export matching Sentaurus fields for `-12.85 V` and
  `-12.9078 V`, or regenerate the Vela probe at `-12.8 V`.

Verification:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_active_support_continuity_balance_summarizes_terms
```

Result: `OK`.

### Next Tasks After Task 79

1. Generate same-bias Sentaurus field exports for `-12.85 V` and `-12.9078 V`
   or rerun the Vela SG probe at `-12.8 V`, then run the support-wide
   continuity-balance diagnostic across the transition window.
2. Compare mobility inputs on active p99 edges after qf-shift density correction:
   use low-field mobility, high-field limiter, carrier density, electric field,
   and quasi-Fermi drive columns to isolate why Sentaurus mobility closes the
   remaining flux/source gap.
3. Inspect carrier-continuity boundary/current balance from contacts to the
   p99 support region, prioritizing hole `phip-psi` because the support-wide
   `~-14.6 mV` offset dominates the density loss.

### Execution Note 2026-06-20: Task 80 Active-Edge Mobility Inputs

Task 80 decomposed the remaining active-edge mobility gap after Task 78 showed
that qf-shifted density plus Sentaurus mobility closes the p99 source mismatch.
The diagnostic compares, on the same active x-edges used by the p99 support
reports:

- Vela final electron/hole mobility;
- Vela low-field mobility;
- Vela high-field mobility limiter;
- `low_field * limiter` versus exported final mobility;
- Sentaurus final e/h mobility;
- edge electric field and qf/high-field drive parity.

Implementation:

```powershell
python scripts/diagnose_pn2d_bv_active_edge_mobility_inputs.py `
  --support-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/thresholded_avalanche_support_m13p2/thresholded_avalanche_support_nodes.csv `
  --sg-edge-csv build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe/sg_avalanche_edges.csv `
  --mesh build-release/reference_tcad/pn2d_sentaurus2018/vela/mesh.json `
  --doping-csv build-release/reference_tcad/pn2d_sentaurus2018/vela/doping.csv `
  --sentaurus-dir build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias/sentaurus_-13.2v `
  --vela-vtk-root build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe/vtk `
  --out-dir build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/active_edge_mobility_inputs_m13p2 `
  --bias -13.2
```

Outputs:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/active_edge_mobility_inputs_m13p2/active_edge_mobility_inputs_nodes.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/active_edge_mobility_inputs_m13p2/active_edge_mobility_inputs_summary.json
```

Important unit note:

- Vela VTK `ElectricField`, `ElectronHighFieldDrive`, and
  `HoleHighFieldDrive` are exported in `V/cm` by `writeDDSolutionVTK`.
- Edge fields and Sentaurus qf fields in the diagnostic are `V/m`, so the new
  script reports both raw `V/cm` and converted `V/m` columns before computing
  drive ratios.

Key `-13.2 V` median results:

| support class | count | e mobility / S | h mobility / S | e final / low-field | h final / low-field | E-field Vela/S | e drive/S qf | h drive/S qf |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| overlap | 10 | `0.942415` | `1.210193` | `0.020288` | `0.054995` | `0.999657` | `1.006918` | `1.006355` |
| false_negative | 10 | `0.942362` | `1.210213` | `0.020286` | `0.054996` | `0.999650` | `1.007059` | `1.006337` |
| false_positive | 10 | `0.942433` | `1.210135` | `0.020289` | `0.054995` | `0.999657` | `1.007014` | `1.006408` |
| all | 30 | `0.942386` | `1.210186` | `0.020287` | `0.054995` | `0.999657` | `1.006999` | `1.006366` |

Additional consistency checks:

- `electron_limited_mobility_product / electron_final_mobility` median is
  about `1.118`.
- `hole_limited_mobility_product / hole_final_mobility` median is about
  `1.039`.
- The product mismatch means the exported node-level limiter/low-field product
  is close but not identical to the edge final mobility used in SG flux. This
  can come from edge versus node averaging and should be checked against the
  exact C++ edge mobility path before treating it as a physics discrepancy.

Interpretation:

- Electric field and high-field driving force are not the remaining mobility
  root cause; after unit conversion they match Sentaurus within about `1%`.
- The residual mobility gap is carrier-specific and directionally important:
  Vela electron mobility is about `5.8%` lower than Sentaurus, while Vela hole
  mobility is about `21%` higher.
- This exactly matches Task 78's mixed-state replay: once qf-density is shifted,
  swapping in Sentaurus mobility removes the residual flux/source gap. Therefore
  the next mobility task is not field-drive debugging; it is parameter/model
  parity for the high-field mobility law and the edge-averaging path.

Verification:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_active_edge_mobility_inputs_decompose_limiter
python -m py_compile scripts/diagnose_pn2d_bv_active_edge_mobility_inputs.py
```

Result: both passed.

### Next Tasks After Task 80

1. Compare Vela high-field mobility parameters and formulas against Sentaurus
   defaults/par-file values for the active-edge field range (`~4.6e7 V/m`):
   saturation velocity, beta/exponent, carrier-specific CT parameters, and any
   Sentaurus parallel-field mobility default.
2. Add an exact edge-mobility replay that uses the same C++ edge mobility path
   as SG assembly and prints low-field, limiter, final mobility, and selected
   driving field per active edge. This will decide whether the
   `low_field * limiter / final` offset is only node/edge averaging or an
   implementation inconsistency.
3. Continue continuity/source-root work in parallel:
   generate same-bias Sentaurus exports for `-12.85 V` and `-12.9078 V`, then
   rerun Task 79 to see where the qf exponent offset appears in the transition
   window.

### Execution Note 2026-06-20: Task 81 Exact C++ Edge-Mobility Probe

Task 81 added an exact C++ edge-mobility probe to separate real SG assembly
mobility from node-export/edge-averaging artifacts. The probe is a new
`vela_example_runner` `simulation_type`:

```json
{
  "simulation_type": "edge_mobility_probe",
  "state_fields_dir": ".../fields",
  "output_csv": ".../edge_mobility_probe_edges.csv",
  "solver": {
    "mobility": {
      "model": "masetti_field",
      "high_field_driving_force": "quasi_fermi_gradient"
    }
  }
}
```

It reuses the same C++ helper used by SG assembly:

```cpp
vela::detail::edgeMobility(...)
```

and writes per-edge:

- electric field;
- electron/hole qf field;
- selected electron/hole mobility field;
- low-field mobility;
- final mobility;
- limiter `final / low-field`;
- edge average net doping and adjacent cell count.

Implementation:

```text
src/tools/vela_example_runner.cpp
tests/regression/test_reference_tcad_tools.py::test_runner_writes_edge_mobility_probe_for_external_state
```

Real `-13.2 V` run:

```powershell
build-release\vela_example_runner.exe --config `
  build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\edge_mobility_probe_m13p2\edge_mobility_probe_config.json
```

Outputs:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/edge_mobility_probe_m13p2/edge_mobility_probe_edges.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/edge_mobility_probe_m13p2/edge_mobility_probe_active_nodes.csv
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/edge_mobility_probe_m13p2/edge_mobility_probe_active_summary.json
```

The active-node summary projects the exact C++ edge probe back onto the same
p99 active edges used by Tasks 76-80.

Key `-13.2 V` median results:

| support class | count | exact e final / VTK node final | exact h final / VTK node final | exact e final / Sentaurus | exact h final / Sentaurus |
| --- | ---: | ---: | ---: | ---: | ---: |
| overlap | 10 | `1.005993` | `1.005870` | `0.948062` | `1.217310` |
| false_negative | 10 | `1.006031` | `1.005866` | `0.948036` | `1.217287` |
| false_positive | 10 | `1.006020` | `1.005909` | `0.948148` | `1.217254` |
| all | 30 | `1.006020` | `1.005874` | `0.948060` | `1.217286` |

Interpretation:

- The exact C++ edge final mobility agrees with the node-export/active-edge
  Task 80 mobility within about `0.6%`.
- Therefore the Task 80 carrier-specific mobility gap is real and is not caused
  by VTK node export, active-edge averaging, or diagnostic reconstruction.
- The low-field and limiter components differ more strongly when comparing
  exact edge mobility to node-export fields:
  - electron low-field exact/node median `~0.826`;
  - electron limiter exact/node median `~1.090`;
  - hole low-field exact/node median `~0.906`;
  - hole limiter exact/node median `~1.068`.
- Those component differences mostly cancel in final mobility, so they are
  useful for understanding edge-vs-node semantics but are not the leading BV
  mismatch.
- The remaining mobility discrepancy to Sentaurus is now best stated using
  exact SG edge mobility:
  - electron mobility is about `5.2%` lower than Sentaurus;
  - hole mobility is about `21.7%` higher than Sentaurus.

Verification:

```powershell
cmake --build build-release --target vela_example_runner --parallel
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_runner_writes_edge_mobility_probe_for_external_state
```

Both passed.

### Next Tasks After Task 81

1. Compare Sentaurus 2018 high-field mobility defaults/par-file values against
   Vela's `MobilityModelConfig` for `masetti_field`:
   - electron and hole saturation velocity;
   - high-field beta/exponent;
   - whether Sentaurus uses GradQuasiFermi or a parallel-field projection with
     a carrier-specific smoothing/refdensity.
2. Prototype a reporting-only Sentaurus-mobility calibration replay:
   - keep Vela qf-shift density from Task 78;
   - vary only electron/hole saturation velocities or high-field beta until
     exact edge mobility ratios approach `1.0`;
   - report predicted active-edge source/current, without changing default
     solver behavior yet.
3. Continue qf exponent root-cause localization:
   - export matching Sentaurus states for `-12.85 V` and `-12.9078 V`;
   - rerun Task 79 support-wide continuity balance across the transition
   window;
   - prioritize the hole `phip-psi` branch offset because it remains the
     dominant density loss.

### Execution Note 2026-06-20: Task 82 Sentaurus High-Field Mobility Defaults

Task 82 closed the mobility-parameter branch from Tasks 80-81.

Local Sentaurus evidence:

```text
reference_tcad/pn2d_sentaurus2018/source/models.par
```

`HighFieldDependence` gives the Sentaurus 2018 Silicon Caughey-Thomas
high-field defaults:

```text
beta0 = 1.109 , 1.213
vsat0 = 1.0700e+07 , 8.3700e+06 cm/s
```

At `300 K`, these convert to:

| carrier | saturation velocity | beta |
| --- | ---: | ---: |
| electron | `1.07e5 m/s` | `1.109` |
| hole | `8.37e4 m/s` | `1.213` |

Before Task 82, Vela used the generic high-field defaults:

```text
electron: vsat = 1.0e5 m/s, beta = 2.0
hole:     vsat = 1.0e5 m/s, beta = 2.0
```

This explained the Task 81 exact edge-mobility discrepancy direction:

- electron final mobility was too low (`~0.948x` Sentaurus);
- hole final mobility was too high (`~1.217x` Sentaurus).

Parameter sensitivity replay:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/edge_mobility_probe_m13p2_sentaurus_hfs
```

Using only Sentaurus high-field `vsat0/beta0` values, exact C++ edge mobility on
the same `-13.2 V` p99 active edges became:

| support class | count | exact e final / Sentaurus | exact h final / Sentaurus |
| --- | ---: | ---: | ---: |
| overlap | 10 | `0.998715` | `0.998476` |
| false_negative | 10 | `0.998688` | `0.998458` |
| false_positive | 10 | `0.998746` | `0.998429` |
| all | 30 | `0.998702` | `0.998456` |

Implementation:

```text
include/vela/physics/MobilityModel.h
docs/config_schema.md
tests/test_mobility.cpp
```

The new default `MobilityModelConfig` high-field parameters are:

```cpp
FieldMobilityParameters electronField{1.07e5, 1.109};
FieldMobilityParameters holeField{8.37e4, 1.213};
```

TDD evidence:

1. Added
   `High-field mobility defaults match Sentaurus 2018 Silicon parameters`.
2. Confirmed RED:
   `100000.0 == Approx(107000.0)` failed before the default update.
3. Updated defaults and confirmed GREEN:
   all `4` assertions in the new test passed.

Post-update real exact edge replay:

```powershell
build-release\vela_example_runner.exe --config `
  build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\edge_mobility_probe_m13p2\edge_mobility_probe_config.json
```

Because this config has no explicit high-field overrides, it now uses the new
defaults. The resulting active p99 mobility ratios match the sensitivity replay:

| support class | count | exact e final / Sentaurus | exact h final / Sentaurus |
| --- | ---: | ---: | ---: |
| overlap | 10 | `0.998715` | `0.998476` |
| false_negative | 10 | `0.998688` | `0.998458` |
| false_positive | 10 | `0.998746` | `0.998429` |
| all | 30 | `0.998702` | `0.998456` |

Interpretation:

- The high-field mobility parameter branch is now closed for the `-13.2 V`
  p99 active edges.
- The previous source/flux residual after qf-density correction should no
  longer require swapping in Sentaurus mobility; Vela's exact SG edge mobility
  now agrees with Sentaurus to about `0.15%` median.
- This does not by itself fix the full BV mismatch, because Tasks 77-79 showed
  the dominant remaining gap is the carrier-density/quasi-Fermi exponent state:
  approximately `-5.3 mV` electron exponent and `-14.6 mV` hole exponent on the
  active support.

Verification:

```powershell
cmake --build build-release --target test_mobility vela_example_runner --parallel
build-release\test_mobility.exe "High-field mobility defaults match Sentaurus 2018 Silicon parameters"
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_runner_writes_edge_mobility_probe_for_external_state tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_active_edge_mobility_inputs_decompose_limiter
```

All passed.

### Next Tasks After Task 82

1. Rerun the qf-shift/source replay with Vela's updated exact edge mobility:
   confirm the former `qf_shift + Sentaurus mobility` success case is now
   reproduced by Vela default mobility.
2. Move the primary root-cause branch back to qf exponent / carrier state:
   export matching Sentaurus fields for `-12.85 V` and `-12.9078 V`, then rerun
   the Task 79 continuity-balance diagnostic across the transition window.
3. After the qf branch is corrected, rerun the actual BV sweep with the updated
   high-field mobility defaults and compare IV/current/source support against
   Sentaurus default BV.

### Execution Note 2026-06-20: Task 83 Transition-Window QF Offset Trend

Task 83 extended the Task 79 active-support continuity-balance diagnostic from
the `-13.2 V` endpoint into the pre-BV transition window.

Inputs:

- Vela C++ SG/VTK probe:
  `build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe`
- Sentaurus intermediate exports:
  `build-release/reference_tcad/pn2d_sentaurus2018/reports/official_split_branch_drift_monitor/sentaurus_intermediate_exports`
- Mesh/doping source:
  `build-release/reference_tcad/pn2d_sentaurus2018/reports/import_split_semantics_smoke/vela`

Commands:

```powershell
python scripts\diagnose_pn2d_bv_thresholded_avalanche_support.py `
  --sentaurus-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports\sentaurus_-12.9v `
  --vela-vtk build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_cpp_sg_probe\vtk\impact_p95_bounded_retry_cpp_sg_0001_-12.9078V.vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\thresholded_avalanche_support_m12p9078_vs_s12p9 `
  --percentile 99

python scripts\diagnose_pn2d_bv_active_support_continuity_balance.py `
  --support-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\thresholded_avalanche_support_m12p9078_vs_s12p9\thresholded_avalanche_support_nodes.csv `
  --sg-edge-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_cpp_sg_probe\sg_avalanche_edges.csv `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json `
  --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\doping.csv `
  --sentaurus-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports\sentaurus_-12.9v `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_cpp_sg_probe\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\active_support_continuity_balance_m12p9078_vs_s12p9 `
  --bias -12.9078
```

Additional weak-evidence run:

```powershell
python scripts\diagnose_pn2d_bv_thresholded_avalanche_support.py `
  --sentaurus-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports\sentaurus_-12.8v `
  --vela-vtk build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_cpp_sg_probe\vtk\impact_p95_bounded_retry_cpp_sg_0000_-12.85V.vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\thresholded_avalanche_support_m12p85_vs_s12p8 `
  --percentile 99
```

The `-12.85 V` run compares against the nearest saved Sentaurus multibias TDR
(`-12.8 V`), so it is trend evidence only. The local Sentaurus TDR sequence has
0.1 V saved plot intervals; `-12.9078 V` versus `-12.9 V` is the stronger
transition-window comparison.

Summary:

| comparison | Jaccard | d(psi-phin) median | d(phip-psi) median | n ratio median | p ratio median | edge G/Sentaurus G median | SRH ratio median |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Vela `-12.85` vs Sentaurus `-12.8` | `0.3333` | `-4.705 mV` | `-13.758 mV` | `0.8335` | `0.5874` | `0.3919` | `1.00000635` |
| Vela `-12.9078` vs Sentaurus `-12.9` | `0.3333` | `-4.896 mV` | `-14.013 mV` | `0.8273` | `0.5816` | `0.3868` | `1.00000636` |
| Vela `-13.2` vs Sentaurus `-13.2` | `0.3333` | `-5.328 mV` | `-14.621 mV` | `0.8138` | `0.5680` | `0.3793` | `1.00000637` |

Observations:

- The qF exponent offsets are already present before the `-13.2 V` endpoint and
  change smoothly across the transition window. They are not a sudden endpoint
  artifact.
- The active support shape is stable: all three comparisons have Jaccard
  `1/3` with 10 overlap, 10 false-positive, and 10 false-negative nodes.
- The median carrier-density ratios track the exponent offsets directly:
  electron density remains about `0.81-0.83x` Sentaurus, while hole density
  remains about `0.57-0.59x`.
- Vela/Sentaurus SRH remains `~1.000006x`, so SRH is still ruled out as the
  active-support carrier-density selector.
- Vela edge/Vela VTK avalanche density remains about `0.38-0.39x` Sentaurus
  across the window, consistent with the carrier-state deficit rather than a
  new avalanche coefficient or mobility branch.

Interpretation:

The remaining dominant mismatch is a smooth carrier-continuity state offset:
Vela solves a quasi-Fermi relation with `psi - phin` about `5 mV` lower and
`phip - psi` about `14-15 mV` lower on the active support. The high-field
mobility parameter branch from Task 82 is closed, and this transition-window
trend moves the next root-cause branch toward carrier-continuity boundary/current
balance, contact minority-carrier anchoring, or quasi-Fermi reference handling.

### Next Tasks After Task 83

1. Run a contact-to-active-support carrier-current balance for `-12.9078 V`
   using the same C++ SG probe state:
   - compare electron/hole terminal current extraction against integrated active
     support transport, SRH, and avalanche terms;
   - report whether the `~5 mV` electron offset or the `~14 mV` hole offset is
     already visible at contact-adjacent nodes.
2. Add or reuse a qF anchor diagnostic that reports contact Dirichlet values,
   nearest interior qF values, and Sentaurus export values at `-12.85`,
   `-12.9078/-12.9`, and `-13.2 V`.
3. If the qF offset is contact-local, inspect minority-carrier boundary
   reconstruction and contact current extraction. If it grows only from the
   interior continuity balance, inspect the carrier continuity residual/Jacobian
   terms for SG flux, SRH sign, and avalanche coupling using a one-row
   finite-difference check at the active-support nodes.

### Execution Note 2026-06-20: Task 84 Transition-Window QF Anchor Check

Task 84 reused the existing qF-anchor diagnostic on the Task 83 transition
window. The Vela C++ SG probe saved exact states at `-12.85 V`, `-12.9078 V`,
and `-13.2 V`, while local Sentaurus TDR exports are saved on `0.1 V`
intervals. For this diagnostic only, a generated alias VTK directory maps:

| alias bias | Vela state |
| --- | --- |
| `-12.8 V` | `-12.85 V` |
| `-12.9 V` | `-12.9078 V` |
| `-13.2 V` | `-13.2 V` |

Generated alias directory:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_cpp_sg_probe_qf_anchor_alias_vtk
```

Command:

```powershell
python scripts\diagnose_pn2d_bv_qf_anchor.py `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json `
  --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_cpp_sg_probe_qf_anchor_alias_vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\qf_anchor_transition_cpp_sg_probe_alias `
  --biases -12.8,-12.9,-13.2 `
  --focus-nodes 351,352,986,987,202,199
```

Outputs:

```text
qf_anchor_contact_summary.csv
qf_anchor_band_summary.csv
qf_anchor_focus_nodes.csv
qf_anchor_summary.json
```

Contact-anchor result:

| bias | contact | delta psi median | delta phin median | note |
| --- | --- | ---: | ---: | --- |
| `-12.8` alias | Cathode | `~6e-9 V` | `~0 V` | aligned |
| `-12.8` alias | Anode | `-50.05 mV` | `-50.00 mV` | expected alias bias mismatch |
| `-12.9` alias | Cathode | `~6e-9 V` | `~0 V` | aligned |
| `-12.9` alias | Anode | `-7.85 mV` | `-7.80 mV` | expected alias bias mismatch |
| `-13.2` | Cathode | `~6e-9 V` | `~0 V` | aligned |
| `-13.2` | Anode | `-0.049 mV` | `~0 V` | aligned |

Band qF exponent result:

| bias | band | d(psi-phin) median | d(phip-psi) median | electron log10 ratio median |
| --- | --- | ---: | ---: | ---: |
| `-12.9` alias | junction | `-4.908 mV` | `-13.999 mV` | `-0.08247` |
| `-12.9` alias | post_junction_n | `-7.265 mV` | `-7.758 mV` | `-0.12208` |
| `-13.2` | junction | `-5.339 mV` | `-14.601 mV` | `-0.08976` |
| `-13.2` | post_junction_n | `-7.837 mV` | `-8.085 mV` | `-0.13169` |

Focus-node examples:

- At node `351`, `-13.2 V`: Vela/Sentaurus `psi` differ by about `+5.82 mV`,
  while `phin` differs by about `+10.31 mV` and `phip` by about `-9.01 mV`.
  This gives the same active-support exponent offsets found in Task 83.
- At node `202`, `-13.2 V`: the qF values are nearly aligned
  (`phin` differs by about `-0.24 mV`, `phip` by roundoff), so the offset is
  not a uniform contact or global-gauge error.

Interpretation:

- Contact Dirichlet/qF anchoring is not the direct source of the Task 83
  `~5 mV` electron and `~14-15 mV` hole active-support exponent offsets.
- The offset appears in the junction and post-junction bands, while the exact
  `-13.2 V` contact medians remain aligned to roundoff.
- The next root-cause branch should inspect interior carrier-continuity balance
  and Jacobian/residual terms around the active support, especially whether the
  SG flux divergence, avalanche source, and SRH terms place a small but
  systematic qF offset into the junction band.

### Next Tasks After Task 84

1. Build a one-row active-support carrier-continuity finite-difference check at
   representative nodes `351/352/986/987`:
   - perturb `phin` by `+/-1 mV` and `phip` by `+/-1 mV`;
   - report electron/hole residual sensitivity from SG flux divergence, SRH,
     and avalanche separately.
2. Compare the same rows against the analytic coupled Jacobian used by Vela. If
   the finite-difference and analytic rows disagree, fix the Jacobian first. If
   they agree, the issue is likely model/discretization parity rather than
   Newton linearization.
3. Keep contact-current extraction as a secondary branch. Task 84 makes a pure
   contact qF anchor bug unlikely for the exact `-13.2 V` state.

### Execution Note 2026-06-20: Task 85 Active-Support Sensitivity And JVP Check

Task 85 executed the active-support one-row finite-difference branch requested
after Task 84. The goal was to determine whether the active-support qF exponent
offset is selected locally by SRH, avalanche generation, or an incorrect
analytic Jacobian row.

Implementation:

```text
scripts/diagnose_pn2d_bv_active_support_sensitivity.py
tests/regression/test_reference_tcad_tools.py::test_pn2d_bv_active_support_sensitivity_reports_term_derivatives
```

The new script reconstructs the Vela VTK state, recomputes carrier densities
from `psi`, `phin`, `phip`, and `ni_eff`, then perturbs one quasi-Fermi variable
at selected active-support nodes by `+/-1 mV`. For each selected node/carrier it
reports central finite-difference derivatives of:

- SG transport flux divergence;
- SRH recombination;
- impact-ionization source;
- total continuity residual.

TDD evidence:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_active_support_sensitivity_reports_term_derivatives
```

The test was first run before the script existed and failed with the expected
missing-file error. After implementing the script, it passed.

Real transition-window runs:

```powershell
python scripts\diagnose_pn2d_bv_active_support_sensitivity.py `
  --support-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\thresholded_avalanche_support_m12p9078_vs_s12p9\thresholded_avalanche_support_nodes.csv `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json `
  --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\doping.csv `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_cpp_sg_probe\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\active_support_sensitivity_m12p9078 `
  --bias -12.9078 `
  --nodes 351,352,986,987 `
  --delta-v 1e-3

python scripts\diagnose_pn2d_bv_active_support_sensitivity.py `
  --support-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\thresholded_avalanche_support_m13p2\thresholded_avalanche_support_nodes.csv `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json `
  --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\doping.csv `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_cpp_sg_probe\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\active_support_sensitivity_m13p2 `
  --bias -13.2 `
  --nodes 351,352,986,987 `
  --delta-v 1e-3
```

Outputs:

```text
active_support_sensitivity_rows.csv
active_support_sensitivity_summary.json
```

Median sensitivity summary:

| bias | carrier | d flux / dqF | d impact / dqF | d SRH / dqF | d residual / dqF |
| --- | --- | ---: | ---: | ---: | ---: |
| `-12.9078 V` | electron | `-1.859826652e8` | `-1.554315e6` | `-2.506e-1` | `-1.844283505e8` |
| `-12.9078 V` | hole | `+1.801325494e8` | `+5.444241e5` | `+8.054e-1` | `+1.795725376e8` |
| `-13.2 V` | electron | `-1.884456291e8` | `-1.624096586e6` | `-2.540e-1` | `-1.868215328e8` |
| `-13.2 V` | hole | `+1.831252103e8` | `+5.772956e5` | `+8.189e-1` | `+1.825329342e8` |

Derivative units are `s^-1/V`.

At both biases, avalanche sensitivity is small relative to SG transport
sensitivity:

```text
electron impact/flux derivative ratio: about 0.0082-0.0087
hole impact/flux derivative ratio:     about 0.0029-0.0034
SRH/flux derivative ratio:             about 1e-9 to 4e-9
```

Representative exact `-13.2 V` rows:

| node | carrier | baseline flux | baseline SRH | baseline impact | baseline residual | d flux / dqF | d impact / dqF |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `351` | electron | `276765.3` | `-184259.25` | `151227.15` | `-58721.11` | `-2.466e8` | `-2.152e6` |
| `351` | hole | `272057.43` | `-184259.25` | `151227.15` | `-63428.98` | `+2.516e8` | `+7.507e5` |
| `986` | electron | `337816.07` | `-155464.20` | `147978.92` | `+34372.95` | `-2.603e8` | `-2.190e6` |
| `986` | hole | `412636.58` | `-155464.20` | `147978.92` | `+109193.46` | `+2.406e8` | `+8.097e5` |

The active rows are therefore controlled primarily by SG transport coupling.
SRH is numerically negligible as a qF selector, and avalanche feedback is too
small in the local row derivative to explain the full active-support qF offset.

Task 85 also generated external-state C++ probes for the same four nodes:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/newton_row_jvp_active_support_m12p9078
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/newton_row_jvp_active_support_m13p2
```

Carrier-row probe results:

| bias | combined block residual | phin/phip block residual | raw carrier step norm | capped carrier step norm |
| --- | ---: | ---: | ---: | ---: |
| `-12.9078 V` | `3565.691` | `2058.650` | `15073.684` | `123.843` |
| `-13.2 V` | `3646.409` | `2105.253` | `15293.487` | `123.189` |

Focus-node C++ rows at `-13.2 V` are locally near balanced in scaled residual
units but still have large carrier-only Newton updates:

| node | e residual | h residual | e diagonal | h diagonal | raw d phin |
| ---: | ---: | ---: | ---: | ---: | ---: |
| `351` | `4.95219e-17` | `-1.23119e-16` | `-6.01824e-14` | `2.90183e-14` | `0.287958 V` |
| `352` | near zero | near zero | `-3.01497e-14` | `1.44975e-14` | `0.288295 V` |
| `986` | near zero | near zero | `-4.18950e-14` | `2.31333e-14` | `0.280473 V` |
| `987` | near zero | near zero | `-2.09680e-14` | `1.15456e-14` | `0.280849 V` |

Selected-node analytic JVP check:

```text
selection box: x = 0.999-1.009 um, y = -0.001-0.016 um
selected nodes: 4
```

| bias | variable | JVP relative error |
| --- | --- | ---: |
| `-12.9078 V` | phin | `1.5625e-16` |
| `-12.9078 V` | phip | `1.1083e-16` |
| `-13.2 V` | phin | `1.52765e-16` |
| `-13.2 V` | phip | `1.02886e-16` |

Interpretation:

- The exact selected-node analytic Jacobian-vector products agree with finite
  differences to roundoff.
- This rejects a local analytic-Jacobian sign/scale bug for the active-support
  qF rows/directions tested here.
- The remaining mismatch should be treated as a global SG transport/continuity
  coupling problem: the active rows are locally balanced, but the coupled
  carrier system still selects qF levels that differ from Sentaurus by the
  Task 83/84 `~5 mV` electron and `~14-15 mV` hole exponent offsets.

Verification:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_active_support_sensitivity_reports_term_derivatives
python -m py_compile scripts/diagnose_pn2d_bv_active_support_sensitivity.py
```

Both passed during Task 85 execution.

### Next Tasks After Task 85

1. Trace the SG transport coupling path from the active-support nodes to the
   contacts and quasi-neutral plateaus:
   - compute carrier qF drops along high-conductance paths;
   - report edge conductance weights and row-coupling strengths;
   - compare Vela and Sentaurus on the same path at `-12.9078/-12.9 V` and
     exact `-13.2 V`.
2. Replay the same active rows and path with Sentaurus-state fields through
   Vela's SG transport evaluator:
   - if both Vela and Sentaurus states are locally residual-balanced, focus on
     global boundary/plateau coupling and row scaling/normalization;
   - if Sentaurus-state rows are not balanced under Vela SG transport, focus on
     SG flux semantics, control-volume ownership, or contact-adjacent carrier
     coupling.
3. Request an exact VM Sentaurus export at `-12.85 V` and `-12.9078 V` only if
   the path analysis needs exact same-bias references. The current `0.1 V`
   Sentaurus alias is sufficient for trend evidence, but not for final parity
   claims.

### Execution Note 2026-06-20: Task 86 Active-Support SG Coupling Paths

Task 86 implemented the first global SG transport-coupling diagnostic after
Task 85 showed that local active-support qF-row sensitivity is dominated by SG
flux divergence and that the analytic JVP agrees with finite differences.

Implementation:

```text
scripts/diagnose_pn2d_bv_sg_coupling_paths.py
tests/regression/test_reference_tcad_tools.py::test_pn2d_bv_sg_coupling_paths_reports_contact_path_drops
```

TDD evidence:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_sg_coupling_paths_reports_contact_path_drops
```

The test first failed because `scripts/diagnose_pn2d_bv_sg_coupling_paths.py`
did not exist. After implementation, it passed.

Important implementation correction:

- The first real run used `abs(edge flux integral)` as the path weight.
- That was insufficient because zero net flux can still have large qF
  conductance.
- The script was corrected to use:

```text
max_abs_central_difference_of_edge_flux_integral_wrt_endpoint_qf
```

with default `delta_v = 1e-4 V`. Reported coupling units are therefore
`s^-1/V`, not `s^-1`.

Real runs:

```powershell
python scripts\diagnose_pn2d_bv_sg_coupling_paths.py `
  --support-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\thresholded_avalanche_support_m12p9078_vs_s12p9\thresholded_avalanche_support_nodes.csv `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json `
  --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\doping.csv `
  --sentaurus-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports\sentaurus_-12.9v `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_cpp_sg_probe\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\sg_coupling_paths_m12p9078_to_anode `
  --bias -12.9078 --nodes 351,352,986,987 --target-contact Anode

python scripts\diagnose_pn2d_bv_sg_coupling_paths.py `
  --support-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\thresholded_avalanche_support_m12p9078_vs_s12p9\thresholded_avalanche_support_nodes.csv `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json `
  --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\doping.csv `
  --sentaurus-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports\sentaurus_-12.9v `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_cpp_sg_probe\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\sg_coupling_paths_m12p9078_to_cathode `
  --bias -12.9078 --nodes 351,352,986,987 --target-contact Cathode

python scripts\diagnose_pn2d_bv_sg_coupling_paths.py `
  --support-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\thresholded_avalanche_support_m13p2\thresholded_avalanche_support_nodes.csv `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json `
  --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\doping.csv `
  --sentaurus-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports\sentaurus_-13.2v `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_cpp_sg_probe\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\sg_coupling_paths_m13p2_to_cathode `
  --bias -13.2 --nodes 351,352,986,987 --target-contact Cathode
```

The exact `-13.2 V -> Anode` strong-coupling run produced `0/8` paths. To
separate topology from SG conductance, Task 86 also ran a structural fallback:

```powershell
python scripts\diagnose_pn2d_bv_sg_coupling_paths.py `
  --support-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\thresholded_avalanche_support_m13p2\thresholded_avalanche_support_nodes.csv `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json `
  --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\doping.csv `
  --sentaurus-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports\sentaurus_-13.2v `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_cpp_sg_probe\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\sg_coupling_paths_m13p2_to_anode_structural `
  --bias -13.2 --nodes 351,352,986,987 --target-contact Anode `
  --min-coupling-s-inv-per-v -1
```

Outputs:

```text
sg_coupling_path_edges.csv
sg_coupling_path_summary.csv
sg_coupling_path_summary.json
```

Median path summary over nodes `351/352/986/987`:

| case | carrier | found | median Vela-Sentaurus total qF drop | median min qF coupling |
| --- | --- | ---: | ---: | ---: |
| `-12.9078 -> Anode` | electron | `4/4` | `-13.94 mV` | `6.953e6 s^-1/V` |
| `-12.9078 -> Anode` | hole | `4/4` | `+4.45 mV` | `1.269e8 s^-1/V` |
| `-12.9078 -> Cathode` | electron | `4/4` | `-6.14 mV` | `1.380e8 s^-1/V` |
| `-12.9078 -> Cathode` | hole | `4/4` | `+12.25 mV` | `3.180e6 s^-1/V` |
| `-13.2 -> Anode` | electron | `0/4` | `NA` | `NA` |
| `-13.2 -> Anode` | hole | `0/4` | `NA` | `NA` |
| `-13.2 -> Cathode` | electron | `4/4` | `-10.45 mV` | `1.401e8 s^-1/V` |
| `-13.2 -> Cathode` | hole | `4/4` | `+8.93 mV` | `3.224e6 s^-1/V` |

Structural fallback for `-13.2 -> Anode` found all `8/8` paths, but every path
has at least one zero-qF-coupling edge:

| carrier | found | median Vela-Sentaurus total qF drop | min qF coupling |
| --- | ---: | ---: | ---: |
| electron | `4/4` | about `-10.45 mV` | `0` |
| hole | `4/4` | about `+8.93 mV` | `0` |

The zero-coupling segment is not at the active-support edge itself. On the
representative horizontal path at `y = 0.03125 um`:

- electron paths first hit zero coupling around `x ~= 0.71875 um` and continue
  toward the Anode-side boundary;
- hole paths first hit zero coupling around `x ~= 0.8125 um` and continue
  toward the Anode-side boundary.

Superseded interpretation:

The interpretation below was the immediate Task 86 reading before the Python
diagnostic SG variable-`ni` qF flux was compared against the C++ implementation.
Task 87 supersedes it. The apparent Anode-side zero-coupling segment was a
diagnostic artifact caused by the Python script using an old absolute-qF
exponential form that clipped both endpoint exponentials to `exp(500)`.

Original Task 86 interpretation:

- The mesh/contact topology is connected; the structural fallback reaches
  Anode.
- At exact `-13.2 V`, the Vela SG qF differential-coupling graph from the
  active support to Anode is numerically cut by a zero-coupling segment.
- Cathode-side paths remain strongly coupled at the same bias.
- The path qF-drop deltas reproduce the already observed active-support offsets
  in global form:
  - electron qF drop is about `10 mV` lower than Sentaurus at `-13.2 V`;
  - hole qF drop is about `9 mV` higher on the Cathode/structural-Anode path;
  - the remaining active-support `phip-psi` mismatch therefore likely involves
    how the carrier system couples through the Anode-side quasi-neutral/low
    conductance segment, not a local active-row Jacobian defect.

Verification:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_sg_coupling_paths_reports_contact_path_drops
```

Passed after the coupling-weight correction.

### Next Tasks After Task 86

1. Add a Sentaurus-state replay mode to the coupling-path diagnostic:
   - compute qF differential coupling using Sentaurus `psi/phin/phip` on the
     same Vela mesh/control volumes;
   - if Sentaurus-state Anode paths remain nonzero under Vela SG, the zero
     segment is caused by the Vela high-bias state;
   - if Sentaurus-state Anode paths also collapse under Vela SG, inspect SG
     derivative semantics, effective `ni`, and log-domain flux derivatives.
2. Localize the zero-coupling segment on the exact edge rows:
   - report `psi`, `phin`, `phip`, `ni_eff`, endpoint density exponents,
     Bernoulli argument, and finite-difference forward/backward flux for the
     first zero edge around `x ~= 0.72-0.81 um`;
   - check whether the zero comes from physical depletion, exponential
     underflow, or a non-log-domain SG derivative cancellation.
3. Only after this path-segment cause is classified, revisit solver damping or
   Bank-Rose-style globalization. Task 85 and Task 86 point first to the
   transport graph/state selection, not to a generic Newton linearization bug.

### Execution Note 2026-06-20: Task 87 Sentaurus Replay And Python SG Diagnostic Correction

Task 87 executed the first follow-up from Task 86 and found that the Task 86
Anode-side zero-coupling conclusion was a Python diagnostic artifact, not a
Vela C++ solver or physical transport result.

TDD additions:

```text
tests/regression/test_reference_tcad_tools.py::test_pn2d_bv_sg_coupling_paths_replays_sentaurus_state_for_coupling
tests/regression/test_reference_tcad_tools.py::test_python_sg_variable_ni_qf_flux_matches_density_form_at_large_bias
```

First RED:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_sg_coupling_paths_replays_sentaurus_state_for_coupling
```

failed because `diagnose_pn2d_bv_sg_coupling_paths.py` did not recognize:

```text
--coupling-state sentaurus
```

Second RED:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_python_sg_variable_ni_qf_flux_matches_density_form_at_large_bias
```

failed on the Task 86 zero-edge representative values:

```text
electron qF flux: 0.0
electron density-form flux: -2.328821677e8
```

Root cause:

- C++ `sgElectronContinuityFluxFromQuasiFermiVariableNi` and
  `sgHoleContinuityFluxFromQuasiFermiVariableNi` already use the stable
  density-index form:

```text
n0 = ni0 * exp((psi0 - phin0) / Vt)
n1 = ni1 * exp((psi1 - phin1) / Vt)
flux = coef * (B(-eta) * n0 - B(eta) * n1)

p0 = ni0 * exp((phip0 - psi0) / Vt)
p1 = ni1 * exp((phip1 - psi1) / Vt)
flux = coef * (B(eta) * p0 - B(-eta) * p1)
```

- Python `scripts/diagnose_pn2d_bv_sg_avalanche_edges.py` still used the older
  absolute-qF factorization. At `qF ~= -13 V`, terms such as
  `exp(-phin/Vt)` clip to `exp(500)` at both endpoints, so their difference
  becomes `0.0`.
- This created the false Task 86 zero-coupling segment in downstream Python
  diagnostics.

Fix:

```text
scripts/diagnose_pn2d_bv_sg_avalanche_edges.py
```

now matches the C++ density-index variable-`ni` qF flux implementation.

The coupling-path diagnostic also gained:

```text
--coupling-state vela|sentaurus
```

and per-edge local derivative diagnostics:

```text
coupling_psi_from_V
coupling_psi_to_V
coupling_qf_from_V
coupling_qf_to_V
coupling_ni_from_m3
coupling_ni_to_m3
coupling_density_exp_from
coupling_density_exp_to
coupling_eta
coupling_forward_flux_from_m2_s
coupling_backward_flux_from_m2_s
coupling_derivative_from_s_inv_per_V
coupling_forward_flux_to_m2_s
coupling_backward_flux_to_m2_s
coupling_derivative_to_s_inv_per_V
```

Real reruns after the Python SG correction:

```powershell
python scripts\diagnose_pn2d_bv_sg_coupling_paths.py `
  --support-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\thresholded_avalanche_support_m13p2\thresholded_avalanche_support_nodes.csv `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json `
  --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\doping.csv `
  --sentaurus-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports\sentaurus_-13.2v `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_cpp_sg_probe\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\sg_coupling_paths_m13p2_to_anode `
  --bias -13.2 --nodes 351,352,986,987 --target-contact Anode

python scripts\diagnose_pn2d_bv_sg_coupling_paths.py `
  --support-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\thresholded_avalanche_support_m13p2\thresholded_avalanche_support_nodes.csv `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json `
  --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\doping.csv `
  --sentaurus-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports\sentaurus_-13.2v `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_cpp_sg_probe\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\sg_coupling_paths_m13p2_to_anode_sentaurus_replay `
  --bias -13.2 --nodes 351,352,986,987 --target-contact Anode `
  --coupling-state sentaurus
```

Updated path summary:

| case | carrier | found | median Vela-Sentaurus total qF drop | median min qF coupling |
| --- | --- | ---: | ---: | ---: |
| `-13.2 -> Anode`, Vela coupling | electron | `4/4` | `-10.45 mV` | `7.541e6 s^-1/V` |
| `-13.2 -> Anode`, Vela coupling | hole | `4/4` | `+8.93 mV` | `1.291e8 s^-1/V` |
| `-13.2 -> Anode`, Sentaurus replay | electron | `4/4` | `-10.45 mV` | `6.790e6 s^-1/V` |
| `-13.2 -> Anode`, Sentaurus replay | hole | `4/4` | `+8.93 mV` | `2.676e8 s^-1/V` |
| `-13.2 -> Cathode`, Vela coupling | electron | `4/4` | `-10.45 mV` | `1.401e8 s^-1/V` |
| `-13.2 -> Cathode`, Vela coupling | hole | `4/4` | `+8.93 mV` | `3.079e6 s^-1/V` |
| `-13.2 -> Cathode`, Sentaurus replay | electron | `4/4` | `-10.45 mV` | `1.979e8 s^-1/V` |
| `-13.2 -> Cathode`, Sentaurus replay | hole | `4/4` | `+8.93 mV` | `4.230e6 s^-1/V` |

Refreshed active-support sensitivity:

```powershell
python scripts\diagnose_pn2d_bv_active_support_sensitivity.py `
  --support-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\thresholded_avalanche_support_m13p2\thresholded_avalanche_support_nodes.csv `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json `
  --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\doping.csv `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_cpp_sg_probe\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\active_support_sensitivity_m13p2 `
  --bias -13.2 --nodes 351,352,986,987 --delta-v 1e-3
```

The refreshed Task 85 medians are unchanged in substance:

| bias | carrier | d flux / dqF | d impact / dqF | d SRH / dqF | d residual / dqF |
| --- | --- | ---: | ---: | ---: | ---: |
| `-12.9078 V` | electron | `-1.859826652e8` | `-1.554314954e6` | `-2.505849e-1` | `-1.844283505e8` |
| `-12.9078 V` | hole | `+1.801325494e8` | `+5.444241e5` | `+8.053527e-1` | `+1.795725376e8` |
| `-13.2 V` | electron | `-1.884456291e8` | `-1.624096586e6` | `-2.539667e-1` | `-1.868215328e8` |
| `-13.2 V` | hole | `+1.831252103e8` | `+5.772956e5` | `+8.188836e-1` | `+1.825329342e8` |

Updated interpretation:

- The active-support-to-contact SG qF differential-coupling graph is connected
  for both Vela and Sentaurus replay states after the Python diagnostic is
  aligned with C++.
- The `-13.2 V` path qF-drop mismatch is not a graph disconnection or a local
  Jacobian zero.
- The same median path qF deltas appear for Vela-coupling and Sentaurus-replay
  path selection, so path selection is not the source.
- The qF-drop mismatch is a real state mismatch concentrated around the
  high-field transition near `x ~= 0.78-0.69 um`, not at the contacts.

For node `351`, `-13.2 V -> Anode`, the largest electron path contributions are:

| step | x range | Vela qF drop | Sentaurus qF drop | delta |
| ---: | --- | ---: | ---: | ---: |
| `21` | `0.71875 -> 0.6875 um` | `-76.50 mV` | `-54.32 mV` | `-22.18 mV` |
| `22` | `0.6875 -> 0.65625 um` | `+2.00 mV` | `-22.15 mV` | `+24.15 mV` |
| `19` | `0.78125 -> 0.75 um` | `-283.50 mV` | `-268.58 mV` | `-14.92 mV` |
| `20` | `0.75 -> 0.71875 um` | `-161.20 mV` | `-176.33 mV` | `+15.13 mV` |

The largest hole path contribution is also in the same neighborhood:

| step | x range | Vela qF drop | Sentaurus qF drop | delta |
| ---: | --- | ---: | ---: | ---: |
| `18` | `0.8125 -> 0.78125 um` | `-211.10 mV` | `-226.32 mV` | `+15.22 mV` |

### Next Tasks After Task 87

1. Build a high-field transition edge-state comparator for the path edges around
   `x ~= 0.78-0.69 um`:
   - compare Vela/Sentaurus `psi`, `phin`, `phip`, `psi-phin`,
     `phip-psi`, electric field, qF field, and inferred carrier density on the
     exact edge sequence;
   - rank whether the qF-drop mismatch follows electrostatic potential shape,
     electron qF shape, hole qF shape, or `ni_eff`/BGN differences.
2. Replay Sentaurus `psi` with Vela qF and Vela `psi` with Sentaurus qF along
   just those transition edges:
   - if Sentaurus `psi` removes most qF-drop mismatch, return to Poisson/
     depletion-shape parity;
   - if Sentaurus qF removes it while Vela `psi` remains, focus on continuity
     boundary/current balance through the transition.
3. Keep Bank-Rose/damping as a stability tool only. Task 87 shows the next
   root-cause branch is state-shape parity through the high-field transition,
   not a disconnected transport graph or an analytic Jacobian zero.

### Execution Note 2026-06-20: Task 88 Transition Edge-State Comparator

Task 88 implemented the high-field transition edge-state comparator requested
after Task 87. The goal was to determine whether the `-13.2 V` qF-drop mismatch
around `x ~= 0.78-0.69 um` follows electrostatic potential shape, electron qF
shape, hole qF shape, or `ni_eff`/carrier-density differences.

Implementation:

```text
scripts/diagnose_pn2d_bv_transition_edge_state.py
tests/regression/test_reference_tcad_tools.py::test_pn2d_bv_transition_edge_state_comparator_ranks_state_terms
```

TDD evidence:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_transition_edge_state_comparator_ranks_state_terms
```

First run failed because:

```text
scripts/diagnose_pn2d_bv_transition_edge_state.py: No such file or directory
```

After implementation, the test passed. A second RED/GREEN added hybrid
decomposition columns:

```text
hybrid_vela_psi_sentaurus_phin_delta_exp_avg_V
hybrid_sentaurus_psi_vela_phin_delta_exp_avg_V
hybrid_vela_phip_sentaurus_psi_delta_exp_avg_V
hybrid_sentaurus_phip_vela_psi_delta_exp_avg_V
```

Real diagnostic run:

```powershell
python scripts\diagnose_pn2d_bv_transition_edge_state.py `
  --path-edge-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\sg_coupling_paths_m13p2_to_anode\sg_coupling_path_edges.csv `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json `
  --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\doping.csv `
  --sentaurus-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports\sentaurus_-13.2v `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_cpp_sg_probe\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\transition_edge_state_m13p2_anode_all_support `
  --bias -13.2 --carrier all --target-contact Anode `
  --x-min-um 0.65 --x-max-um 0.82
```

Outputs:

```text
transition_edge_state_compare.csv
transition_edge_state_summary.json
```

All-support transition-window summary over path edges in `x = 0.65-0.82 um`:

| metric | median absolute mismatch |
| --- | ---: |
| electron qF drop | `9.117 mV` |
| electron qF field | `2.9175e5 V/m` |
| electron exponent average | `3.221 mV` |
| electrostatic potential drop | `0.776 mV` |
| electric field | `2.4839e4 V/m` |
| hole qF drop | `~6.7e-12 V` |
| hole exponent average | `0.779 mV` |

Hybrid decomposition:

| hybrid | median absolute exponent mismatch |
| --- | ---: |
| `Vela psi + Sentaurus phin` | `0.779 mV` |
| `Sentaurus psi + Vela phin` | `3.259 mV` |
| `Vela phip + Sentaurus psi` | `~4e-12 V` |
| `Sentaurus phip + Vela psi` | `0.779 mV` |

Interpretation:

- In the high-field transition window, the dominant mismatch is electron qF
  shape/gradient, not Poisson/electrostatic potential shape.
- Replacing Vela `psi` with Sentaurus `psi` would not remove the electron
  exponent mismatch if Vela `phin` is kept.
- Replacing Vela `phin` with Sentaurus `phin` collapses the electron exponent
  mismatch to the small residual `psi` mismatch (`~0.78 mV`).
- Hole qF drop is essentially aligned on these edges; the remaining hole
  exponent mismatch is mostly the same small `psi` residual.
- This redirects the next root-cause branch from Poisson/depletion-shape parity
  toward electron continuity/qF state selection through the transition.

Representative `start_node = 351`, `-13.2 V -> Anode` edge rows:

| step | edge mid x | delta psi drop | delta phin drop | delta phip drop | dominant |
| ---: | ---: | ---: | ---: | ---: | --- |
| `19` | `0.765625 um` | `-0.776 mV` | `-14.923 mV` | `~0.001 mV` | electron qF drop |
| `20` | `0.734375 um` | `-0.778 mV` | `+15.133 mV` | `~0 mV` | electron qF drop |
| `21` | `0.703125 um` | `-0.352 mV` | `-22.183 mV` | `~0 mV` | electron qF drop |
| `22` | `0.671875 um` | `+0.005 mV` | `+24.152 mV` | `~0 mV` | electron qF drop |

The electron qF mismatch has a local oscillatory/over-corrected shape through
the transition: adjacent edges alternate sign while producing the net
active-support qF exponent offset found in Tasks 83-87.

### Execution Note 2026-06-20: Task 89 Transition Electron Flux Replay

Task 89 added a focused electron SG flux replay diagnostic for the Task 88
transition path edges. It reuses Vela geometry, edge coupling, mobility, and
`ni_eff`, but swaps the accepted state between:

```text
vela
sentaurus
vela_psi_sentaurus_phin
sentaurus_psi_vela_phin
```

Implementation:

```text
scripts/diagnose_pn2d_bv_transition_flux_replay.py
tests/regression/test_reference_tcad_tools.py::test_pn2d_bv_transition_flux_replay_compares_hybrid_states
```

TDD evidence:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_transition_flux_replay_compares_hybrid_states
```

First run failed because the script did not exist. After implementation, the
test passed.

Real all-support run:

```powershell
python scripts\diagnose_pn2d_bv_transition_flux_replay.py `
  --path-edge-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\sg_coupling_paths_m13p2_to_anode\sg_coupling_path_edges.csv `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json `
  --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\doping.csv `
  --sentaurus-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports\sentaurus_-13.2v `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_cpp_sg_probe\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\transition_flux_replay_m13p2_anode_all_support `
  --bias -13.2 --target-contact Anode --x-min-um 0.65 --x-max-um 0.82
```

Outputs:

```text
transition_flux_replay_edges.csv
transition_flux_replay_summary.json
```

All-support median absolute electron flux-integral mismatch versus Sentaurus:

| state | median absolute mismatch |
| --- | ---: |
| `vela` | `9.0767e4 s^-1` |
| `sentaurus_psi_vela_phin` | `1.2147e5 s^-1` |
| `vela_psi_sentaurus_phin` | `3.9914e3 s^-1` |

Inner transition run with `x = 0.65-0.78 um`:

| state | median absolute mismatch |
| --- | ---: |
| `vela` | `1.3583e5 s^-1` |
| `sentaurus_psi_vela_phin` | `1.2147e5 s^-1` |
| `vela_psi_sentaurus_phin` | `1.0669e3 s^-1` |

Representative `start_node = 351` rows:

| step | edge mid x | Vela delta | Vela psi + Sentaurus phin delta |
| ---: | ---: | ---: | ---: |
| `19` | `0.765625 um` | `-1.3583e5 s^-1` | `-8.4392e3 s^-1` |
| `21` | `0.703125 um` | `-2.3403e5 s^-1` | `-1.0669e3 s^-1` |
| `22` | `0.671875 um` | `+2.4036e5 s^-1` | `-3.3149e2 s^-1` |

The diagnostic window uses edge x-interval overlap, not strict midpoint
clipping, so boundary-crossing edges can appear slightly outside the requested
midpoint range.

Interpretation:

- The same state decomposition from Task 88 appears in the actual SG flux
  replay: swapping in Sentaurus electron qF on Vela electrostatic potential
  nearly closes the transition-edge electron flux mismatch.
- Swapping only Sentaurus electrostatic potential while keeping Vela electron
  qF does not close the flux mismatch and can be worse than Vela.
- This makes the immediate root-cause target the electron continuity equation's
  accepted `phin` profile through the transition, not the Poisson field, SG
  formula, contact current extraction, or a disconnected coupling graph.

### Execution Note 2026-06-20: Task 90 Transition Row Residual Decomposition

Task 90 implemented the electron-continuity row residual/Jacobian decomposition
requested after Task 89. It aggregates selected transition-edge electron SG
flux integrals into node residual contributions using the existing continuity
sign convention:

```text
node_from += flux_integral
node_to   -= flux_integral
```

It also reports finite-difference derivatives with respect to each endpoint
`phin` and a local delta-Newton step estimate for the selected edge subset.

Implementation:

```text
scripts/diagnose_pn2d_bv_transition_row_residual.py
tests/regression/test_reference_tcad_tools.py::test_pn2d_bv_transition_row_residual_decomposition_groups_node_terms
```

TDD evidence:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_transition_row_residual_decomposition_groups_node_terms
```

First run failed because:

```text
scripts\diagnose_pn2d_bv_transition_row_residual.py: No such file or directory
```

After implementation, the test passed.

Real `-13.2 V -> Anode` transition row decomposition:

```powershell
python scripts\diagnose_pn2d_bv_transition_row_residual.py `
  --path-edge-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\sg_coupling_paths_m13p2_to_anode\sg_coupling_path_edges.csv `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json `
  --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\doping.csv `
  --sentaurus-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports\sentaurus_-13.2v `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_cpp_sg_probe\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\transition_row_residual_m13p2_anode_inner `
  --bias -13.2 --target-contact Anode --x-min-um 0.65 --x-max-um 0.78
```

Inner transition median absolute node residual mismatch versus Sentaurus:

| state | median absolute mismatch |
| --- | ---: |
| `vela` | `7.1677e5 s^-1` |
| `sentaurus_psi_vela_phin` | `6.9488e5 s^-1` |
| `vela_psi_sentaurus_phin` | `7.3198e3 s^-1` |

Representative `start_node = 351` inner path:

| node | x | Vela residual delta | Vela psi + Sentaurus phin residual delta |
| ---: | ---: | ---: | ---: |
| `188` | `0.65625 um` | `-2.1562e5 s^-1` | `-2.51e1 s^-1` |
| `191` | `0.6875 um` | `+4.7439e5 s^-1` | `+7.35e2 s^-1` |
| `198` | `0.71875 um` | `-2.4096e5 s^-1` | `+2.9245e3 s^-1` |
| `201` | `0.75 um` | `+1.4277e5 s^-1` | `+4.4478e3 s^-1` |
| `224` | `0.78125 um` | `-1.3583e5 s^-1` | `-8.4392e3 s^-1` |

This confirms the Task 88/89 conclusion at the node-row level: the high-bias
transition residual gap is carried by Vela's electron `phin` profile. Keeping
Vela electrostatic potential while replaying Sentaurus `phin` reduces selected
row residual mismatch by about two orders of magnitude on the representative
inner path.

Pre-strong-avalanche comparison (`Vela -12.9078 V` against Sentaurus `-12.9 V`):

```powershell
python scripts\diagnose_pn2d_bv_sg_coupling_paths.py `
  --support-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\thresholded_avalanche_support_m12p9078_vs_s12p9\thresholded_avalanche_support_nodes.csv `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json `
  --doping-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\doping.csv `
  --sentaurus-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports\sentaurus_-12.9v `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_cpp_sg_probe\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\sg_coupling_paths_m12p9078_to_anode `
  --bias -12.9078 --target-contact Anode
```

The pre-strong-avalanche Anode coupling graph is connected (`60/60` paths).
In the same `x = 0.65-0.82 um` window:

| diagnostic | best state | Vela mismatch | Vela psi + Sentaurus phin mismatch |
| --- | --- | ---: | ---: |
| edge electron qF drop | electron qF still dominant but small | `1.1046 mV` median drop | not preferred |
| electron flux replay | `vela` | `1.8288e4 s^-1` | `6.0269e4 s^-1` |
| row residual decomposition | `vela` | `3.6076e4 s^-1` | `6.1158e4 s^-1` |

In the inner `x = 0.65-0.78 um` window, edge flux replay still prefers Vela
(`1.6711e4 s^-1` versus `5.3208e4 s^-1` for `Vela psi + Sentaurus phin`).
The row residual diagnostic mildly prefers `Vela psi + Sentaurus phin`
(`1.3139e4 s^-1` versus `2.3632e4 s^-1`), but the gap is small compared with
the `-13.2 V` two-order-of-magnitude split.

Interpretation:

- The large oscillatory electron `phin`/row-residual mismatch is not already
  fully present at `-12.9078 V`.
- It appears or is strongly amplified between the accepted `-12.9078 V` Vela
  state and the `-13.2 V` endpoint.
- Therefore the next root-cause branch should focus on avalanche-coupled
  electron continuity feedback and nonlinear continuation between these
  voltages, rather than SG formula parity or low-bias contact current
  extraction.

### Execution Note 2026-06-20: Task 91 Source-Adjusted Transition Row Residual

Task 91 extended the Task 90 row residual decomposition with frozen source
terms:

- Vela source terms from VTK:
  - `AvalancheGeneration * node_volume`
  - `SRHRecombination * node_volume`
- Sentaurus source terms from nearest exported node:
  - `ImpactIonization * 1e6 * node_volume`
  - `srhRecombination * 1e6 * node_volume`

The diagnostic now reports both transport-only selected-edge residual and
source-adjusted local electron residual:

```text
electron_full_residual = selected_transport + srh - avalanche
```

Implementation update:

```text
scripts/diagnose_pn2d_bv_transition_row_residual.py
tests/regression/test_reference_tcad_tools.py::test_pn2d_bv_transition_row_residual_decomposition_groups_node_terms
```

TDD evidence:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_transition_row_residual_decomposition_groups_node_terms
```

The first extended test failed because source/full-residual fields were absent.
A second RED/GREEN added `median_abs_source_delta_vs_sentaurus_s_inv` to the
summary so source-term scale is reproducible from the JSON report.

Real `-13.2 V -> Anode`, inner transition window `x = 0.65-0.78 um`:

| state | transport median mismatch | source-adjusted median mismatch |
| --- | ---: | ---: |
| `vela` | `7.1677e5 s^-1` | `7.1674e5 s^-1` |
| `sentaurus_psi_vela_phin` | `6.9488e5 s^-1` | `6.9486e5 s^-1` |
| `vela_psi_sentaurus_phin` | `7.3198e3 s^-1` | `7.3198e3 s^-1` |

The source-delta median is only `2.11e-2 s^-1` in this inner window. The
largest sampled Vela-state source delta in the same window is about
`2.43e2 s^-1`, while transport residual mismatch reaches `1.90e6 s^-1`.

Real `-13.2 V -> Anode`, all-support transition window `x = 0.65-0.82 um`:

| state | transport median mismatch | source-adjusted median mismatch |
| --- | ---: | ---: |
| `vela` | `4.7507e5 s^-1` | `4.7517e5 s^-1` |
| `sentaurus_psi_vela_phin` | `7.3808e5 s^-1` | `7.5354e5 s^-1` |
| `vela_psi_sentaurus_phin` | `1.4745e4 s^-1` | `1.4718e4 s^-1` |

The source-delta median is `2.70e1 s^-1`, still far below the transport-driven
row residual mismatch.

Pre-strong-avalanche comparison, `Vela -12.9078 V` against Sentaurus `-12.9 V`:

| window | best source-adjusted state | Vela mismatch | Vela psi + Sentaurus phin mismatch | source-delta median |
| --- | --- | ---: | ---: | ---: |
| `x = 0.65-0.78 um` | `vela_psi_sentaurus_phin` | `2.4242e4 s^-1` | `1.3134e4 s^-1` | `4.94e-3 s^-1` |
| `x = 0.65-0.82 um` | `vela` | `3.1544e4 s^-1` | `6.1158e4 s^-1` | `1.12e0 s^-1` |

Interpretation:

- Adding frozen local source terms does not change the `-13.2 V` conclusion:
  the mismatch is still overwhelmingly a selected-edge electron transport /
  `phin` profile mismatch.
- Source terms are not large enough on these transition-row nodes to flip the
  required correction direction.
- Since `-12.9078 V` does not yet show the same large `phin` transport split,
  the remaining root-cause target is the nonlinear continuation/Newton update
  between `-12.9078 V` and `-13.2 V`, not the instantaneous avalanche/SRH source
  magnitude on the final selected rows.

### Execution Note 2026-06-20: Task 92 Transition Continuation Overshoot

Task 92 added a continuation overshoot diagnostic for the transition nodes
identified by Tasks 88-91. It compares available accepted Vela VTK states
against nearest Sentaurus intermediate exports and reports node-level:

```text
psi
phin
psi - phin
electron density ratio
phin mismatch class
```

Implementation:

```text
scripts/diagnose_pn2d_bv_transition_continuation_overshoot.py
tests/regression/test_reference_tcad_tools.py::test_pn2d_bv_transition_continuation_overshoot_flags_sparse_jump
```

TDD evidence:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_transition_continuation_overshoot_flags_sparse_jump
```

First run failed because:

```text
scripts\diagnose_pn2d_bv_transition_continuation_overshoot.py: No such file or directory
```

After implementation, the test passed. A second RED/GREEN added
`first_alternating_large_phin_mismatch_bias_V` so whole-window offset and
alternating transition-node oscillation are not conflated.

Available accepted Vela states in the current
`impact_p95_guard_bounded_retry_regression_cpp_sg_probe` run:

```text
-12.85 V
-12.9078 V
-13.2 V
```

Sentaurus intermediate exports are available at:

```text
-12.5, -12.6, -12.7, -12.8, -12.9, -13.0, -13.1, -13.2 V
```

Real focused run over `-12.9078 V -> -13.2 V`:

```powershell
python scripts\diagnose_pn2d_bv_transition_continuation_overshoot.py `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\impact_p95_guard_bounded_retry_regression_cpp_sg_probe\vtk `
  --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\transition_continuation_overshoot_nodes_165_351_m129078_m132 `
  --bias-min -13.2 --bias-max -12.9078 `
  --nodes 165,188,191,198,201,224,351 `
  --phin-threshold-v 0.01 --max-accepted-gap-v 0.1
```

Summary:

| accepted Vela state | nearest Sentaurus | max abs `delta phin` | sign alternations | median log10 electron-density ratio |
| --- | --- | ---: | ---: | ---: |
| `-12.9078 V` | `-12.9 V` | `7.42 mV` | `0` | `-0.058` |
| `-13.2 V` | `-13.2 V` | `20.79 mV` | `1` | `-0.016` |

The first large and alternating `phin` mismatch is therefore only localized to
the sparse accepted-state interval:

```text
(-12.9078 V, -13.2 V]
```

The accepted-state gap before the first large alternating mismatch is
`0.2922 V`, above the `0.1 V` diagnostic threshold, so the script reports:

```text
needs_intermediate_restart = true
```

Representative node deltas:

| bias | node | x | delta phin | delta(psi-phin) | log10 n ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| `-12.9078` | `191` | `0.6875 um` | `-3.07 mV` | `-4.71 mV` | `-0.079` |
| `-12.9078` | `201` | `0.75 um` | `-6.44 mV` | `+0.36 mV` | `+0.006` |
| `-13.2` | `191` | `0.6875 um` | `-20.79 mV` | `+20.83 mV` | `+0.350` |
| `-13.2` | `201` | `0.75 um` | `-13.74 mV` | `+14.91 mV` | `+0.250` |
| `-13.2` | `351` | `1.0 um` | `+10.31 mV` | `-4.49 mV` | `-0.075` |

Interpretation:

- The final `-13.2 V` transition mismatch is not simply a uniform electron-qF
  offset. It has a spatially alternating sign between transition nodes.
- The current accepted Vela artifacts do not contain enough intermediate
  states to say whether the oscillation appears near `-13.0 V`, near `-13.1 V`,
  or only in the final step to `-13.2 V`.
- The next task must generate intermediate accepted states or residual probes
  inside `(-12.9078 V, -13.2 V]` before changing physics or solver damping.

### Next Tasks After Task 92

1. Generate a focused restart/sweep from the accepted `-12.9078 V` state with
   output near `-13.0 V`, `-13.1 V`, and `-13.2 V`:
   - include VTK output for every accepted point;
   - include Newton step/residual diagnostics for transition nodes
     `165,188,191,198,201,224,351`;
   - keep solver physics identical to the Task 83/91 run.
2. Re-run Task 92 on the denser VTK sequence and identify the first accepted
   state where `phin_sign_alternations > 0` and `max_abs_delta_phin > 10 mV`.
3. Only then test damping/Bank-Rose style globalization. The acceptance target
   is not just convergence: the dense sequence must reduce the transition-node
   `phin` alternation and the Task 91 source-adjusted row residual mismatch.

### Execution Note 2026-06-20: Task 93 Focused Restart Through `-13.2 V`

Task 93 implemented the focused restart generator and used it to densify the
critical `-12.9078 V -> -13.2 V` interval identified by Task 92.

Implementation:

```text
scripts/prepare_pn2d_bv_focused_restart.py
tests/regression/test_reference_tcad_tools.py::test_prepare_pn2d_bv_focused_restart_writes_restart_and_config
scripts/diagnose_pn2d_bv_transition_edge_state.py
scripts/diagnose_pn2d_bv_transition_flux_replay.py
scripts/diagnose_pn2d_bv_transition_row_residual.py
```

TDD evidence:

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_prepare_pn2d_bv_focused_restart_writes_restart_and_config
```

The real focused restart was generated from the accepted VTK state:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/impact_p95_guard_bounded_retry_regression_cpp_sg_probe/vtk/impact_p95_bounded_retry_cpp_sg_0001_-12.9078V.vtk
```

Output workspace:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m12p9078_m13p2
```

Generated artifacts:

```text
simulation.json
restart_from_vtk.csv
iv.csv
vtk/focused_restart_0000_-12.9078V.vtk
vtk/focused_restart_0001_-13V.vtk
vtk/focused_restart_0002_-13.1V.vtk
vtk/focused_restart_0003_-13.2V.vtk
iv_compare_as_iv.json
iv_compare_as_iv.md
```

Focused sweep result:

```text
converged = true
points = 4
```

Dense accepted-state IV. The `-12.9078 V` row uses the nearest available
Sentaurus reference point (`-12.9 V`) for orientation; the exact official
comparator rows are `-13.0`, `-13.1`, and `-13.2 V`.

| bias | Vela `current_total_A_per_um` | Sentaurus reference current | relative error | log10 ratio |
| --- | ---: | ---: | ---: | ---: |
| `-12.9078 V` | `-5.8946e-17` | `-8.0989e-17` | `27.22%` | `-0.1380` |
| `-13.0 V` | `-5.8553e-17` | `-8.1935e-17` | `28.54%` | `-0.1459` |
| `-13.1 V` | `-5.8772e-17` | `-8.2884e-17` | `29.09%` | `-0.1493` |
| `-13.2 V` | `-5.8560e-17` | `-8.3847e-17` | `30.16%` | `-0.1559` |

The official curve comparator, using the BV reference current as an IV-style
current column comparison over the exact reference points inside the focused
window (`-13.0`, `-13.1`, `-13.2 V`), reports:

```text
points_compared = 3
candidate_column = current_total_A_per_um
candidate_scale = 1.0
orders_of_magnitude = 0.15588996512771738
max_relative_error = 0.30159066609658763
status = pass
```

Re-running Task 92 on the dense focused VTK sequence changed the continuation
classification:

| accepted Vela state | nearest Sentaurus | max abs `delta phin` | sign alternations | median log10 electron-density ratio |
| --- | --- | ---: | ---: | ---: |
| `-12.9078 V` | `-12.9 V` | `7.84 mV` | `0` | `+0.00006` |
| `-13.0 V` | `-13.0 V` | `11.43 mV` | `0` | `+0.02196` |
| `-13.1 V` | `-13.1 V` | `11.60 mV` | `0` | `-0.00387` |
| `-13.2 V` | `-13.2 V` | `11.75 mV` | `0` | `-0.00802` |

The dense sequence reports:

```text
first_large_phin_mismatch_bias_V = -13.0
first_alternating_large_phin_mismatch_bias_V = null
largest_accepted_gap_V = 0.1
needs_intermediate_restart = false
```

Focused residual and flux probes at `-13.2 V`, `x = 0.65-0.78 um`, show that
the smaller-step restart collapses the earlier high-field transport error:

| diagnostic | sparse `-13.2 V` Vela mismatch | focused `-13.2 V` Vela mismatch | focused best hybrid |
| --- | ---: | ---: | ---: |
| electron SG flux replay | `1.3583e5 s^-1` | `3.9191e3 s^-1` | `1.0197e3 s^-1` |
| source-adjusted row residual | `7.1674e5 s^-1` | `3.0400e4 s^-1` | `6.6655e3 s^-1` |

Focused edge-state comparison:

| metric | focused median abs mismatch |
| --- | ---: |
| electron qF drop | `0.9478 mV` |
| electron qF field | `3.0328e4 V/m` |
| electron SG exponent average | `1.2350 mV` |
| psi drop | `0.3523 mV` |
| electric field | `1.1274e4 V/m` |
| `Vela psi + Sentaurus phin` hybrid exponent mismatch | `0.2134 mV` |

Interpretation:

- The Task 92 sparse `-12.9078 V -> -13.2 V` jump was a continuation-path
  artifact: densifying the restart interval removes the transition-node
  `phin` sign alternation and reduces the electron SG flux mismatch by about
  `35x`.
- The remaining `-13.2 V` IV gap is no longer dominated by the large alternating
  electron-qF branch error. It is a smoother current-magnitude gap of about
  `30%` across `-13.0..-13.2 V`.
- Because the focused row residual still prefers `Vela psi + Sentaurus phin`,
  the next root-cause target remains electron quasi-Fermi transport, but now in
  a small-signal, stable-branch regime rather than a gross branch-jump regime.

Review hardening after Task 93:

- transition edge-state, flux-replay, and row-residual diagnostics now fail
  explicitly when the path-edge filters select no rows;
- transition diagnostics now validate required Vela VTK scalar arrays before
  indexing, with clear missing-scalar/short-length errors;
- regression coverage was added for empty edge selection and missing
  `ElectronMobility` in row-residual replay.

### Next Tasks After Task 93

1. Promote the focused restart result into a reusable guard:
   - add a production or diagnostic config option that caps accepted reverse-BV
     continuation gaps near high-field transition regions;
   - acceptance should use `carrier_density_jump` and/or transition-node
     `psi-phin` diagnostics, not convergence alone.
2. Localize the remaining stable-branch `~30%` current deficit:
   - compare focused `-13.0`, `-13.1`, and `-13.2 V` contact current
     decomposition with and without `contact_current_qf_floor`;
   - check whether terminal current consistency degradation at focused
     `-13.2 V` (`0.6665`) is a reporting/extraction issue or a real continuity
     imbalance.
3. Re-run focused field comparison for the dense restart VTKs:
   - quantities: potential, electric field, electron density, hole density,
     electron mobility, hole mobility, avalanche generation;
   - biases: `-13.0`, `-13.1`, `-13.2 V`;
   - use this to decide whether the remaining current error follows mobility,
     avalanche driving force, contact-current extraction, or SG flux assembly.
4. Only after Tasks 1-3 classify the residual error, test Bank-Rose style
   damping/globalization. The expected benefit is stability and path control;
   it should not be treated as a substitute for matching the stable-branch
   physics/current extraction.

### Execution Note 2026-06-20: Task 94 Focused Stable-Branch Field/Source Split

Task 94 ran the dense focused restart states through the existing multibias
field comparator and source-factor diagnostics to classify the remaining
stable-branch current deficit after Task 93 removed the gross continuation
branch jump.

Field/current comparison workspace:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13x_field_compare
```

Command:

```powershell
python scripts\compare_pn2d_bv_multibias_fields.py `
  --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\focused_restart_m12p9078_m13p2\vtk `
  --curve-reference build-release\reference_tcad\pn2d_sentaurus2018\reference_curves\pn2d_sentaurus2018_bv_reference.csv `
  --curve-candidate build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\focused_restart_m12p9078_m13p2\iv.csv `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\focused_restart_m13x_field_compare `
  --biases -13.0,-13.1,-13.2 `
  --quantities potential,electric_field,electron_density,hole_density,electron_mobility,hole_mobility,avalanche_generation
```

Curve-current result:

| bias | Vela/Sentaurus current magnitude | Sentaurus/Vela current magnitude | log10 ratio |
| --- | ---: | ---: | ---: |
| `-13.0 V` | `0.7146` | `1.3993` | `-0.1459` |
| `-13.1 V` | `0.7091` | `1.4103` | `-0.1493` |
| `-13.2 V` | `0.6984` | `1.4318` | `-0.1559` |

Signed field/source ratios were computed by reusing
`compare_pn2d_bv_multibias_fields.py` loaders and nearest-node mapping. The
important junction-local signed medians are:

| bias | electron density log10(V/S) | hole density log10(V/S) | electron mobility relative(V-S)/S | hole mobility relative(V-S)/S | avalanche log10(V/S) |
| --- | ---: | ---: | ---: | ---: | ---: |
| `-13.0 V` | `-0.1077` | `-0.1467` | `-0.0166` | `-0.0162` | `-0.4149` |
| `-13.1 V` | `-0.1102` | `-0.1497` | `-0.0166` | `-0.0161` | `-0.4177` |
| `-13.2 V` | `-0.1126` | `-0.1529` | `-0.0165` | `-0.0161` | `-0.4204` |

Interpretation:

- Potential is close (`~3.4 mV` RMS, junction `~4.6 mV`), and the high-field
  electric-field maximum/p95 is essentially aligned. The global electric-field
  relative p95 remains large only because low-field regions have small
  denominators; field-source summary shows peak field log10 ratios within
  `~1.3e-4`.
- Mobility is not the leading residual. Junction mobility is only about
  `1.6%` low, far smaller than the `0.146..0.156 dex` current deficit.
- Carrier densities are systematically low in the junction active region:
  electrons by `~0.11 dex` and holes by `~0.15 dex`.
- Avalanche generation is lower by `~0.42 dex` in the same active region, but
  this follows the carrier/SG-current deficit rather than a field or
  ionization-coefficient mismatch.

Field-source summary workspace:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13x_field_source_summary
```

At `-13.2 V`, the electric-field and avalanche source summaries report:

| quantity | Sentaurus max | Vela max | log10 V/S max | Sentaurus p95 | Vela p95 | log10 V/S p95 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| electric field (`V/cm`) | `4.58916e5` | `4.58776e5` | `-1.32e-4` | `4.52392e5` | `4.58737e5` | `+0.00605` |
| avalanche (`cm^-3 s^-1`) | `3.19785e15` | `1.24192e15` | `-0.4108` | `3.03143e15` | `1.14486e15` | `-0.4229` |

The top 1% avalanche support comparison at `-13.2 V` shows the same source
deficit without a gross spatial miss:

```text
sentaurus_active_count = 20
vela_active_count = 20
overlap_count = 10
jaccard = 0.3333333333333333
peak_separation_um = 0.015625
sentaurus_active_sum_cm3_s = 6.393206116376254e16
vela_active_sum_cm3_s = 2.445126e16
```

The Vela/Sentaurus active-source ratio is therefore about `0.3825`, close to
the junction avalanche `-0.420 dex` result.

Local avalanche factor workspace:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13x_local_avalanche_factors
```

The default selected edge remains edge `2886`, nodes `351-986`, near
`x = 1.00390625 um`, `y = 0.015625 um`. At `-13.2 V`:

| factor | Vela | Sentaurus | Vela/Sentaurus |
| --- | ---: | ---: | ---: |
| source density (`m^-3 s^-1`) | `1.6856e21` | `3.1979e21` | `0.527` |
| electric field (`V/m`) | `4.5878e7` | `4.5892e7` | `0.9997` |
| electron alpha (`m^-1`) | `4.5978e6` | `4.8084e6` from field | `0.956` |
| hole alpha (`m^-1`) | `1.7272e6` | `1.6771e6` from field | `1.030` |
| electron flux magnitude (`m^-2 s^-1`) | `2.6834e14` | `4.7739e14` | `0.562` |
| hole flux magnitude (`m^-2 s^-1`) | `2.6162e14` | `5.6930e14` | `0.460` |
| electron density (`cm^-3`) | `3860.02` | `4723.84` | `0.817` |
| hole density (`cm^-3`) | `4781.49` | `7202.61` | `0.664` |

Interpretation:

- The local ionization coefficients are close enough that they cannot explain
  the source deficit.
- The local SG particle fluxes are low because the high-field carrier/qF state
  is low; this is the same mechanism seen in the signed field comparison.
- Mobility differences are negligible in this local source budget.

Mixed-state replay and qF-shift sensitivity:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p2_active_edge_mixed_state_replay
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p2_active_edge_mixed_state_replay_qf_shift_m6p5_p9p6mv
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p0_active_edge_mixed_state_replay_qf_shift_m6p5_p9p6mv
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p1_active_edge_mixed_state_replay_qf_shift_m6p5_p9p6mv
```

At `-13.2 V`, overlap active-edge median ratios:

| variant | generation/Sentaurus | particle flux/Sentaurus | electron density/Sentaurus | hole density/Sentaurus |
| --- | ---: | ---: | ---: | ---: |
| `vela_baseline` | `0.7536` | `0.7307` | `0.7765` | `0.6907` |
| `sentaurus_psi_vela_qf` | `0.7041` | `0.7433` | `0.6564` | `0.8202` |
| `vela_psi_sentaurus_qf` | `1.0801` | `1.0014` | `1.1832` | `0.8452` |
| `vela_qf_shift` with `phin -= 6.5 mV`, `phip += 9.6 mV` | `1.0022` | `1.0012` | `0.9984` | `1.0012` |

The same uniform qF-shift sensitivity is stable across the focused window:

| bias | `vela_baseline` generation/S | shifted generation/S | shifted particle flux/S |
| --- | ---: | ---: | ---: |
| `-13.0 V` | `0.7632` | `1.0148` | `1.0150` |
| `-13.1 V` | `0.7585` | `1.0087` | `1.0082` |
| `-13.2 V` | `0.7536` | `1.0022` | `1.0012` |

Conclusion:

- The remaining focused-window current/source gap is now localized to a small,
  nearly uniform high-field quasi-Fermi/carrier-density offset:
  approximately `phin` too high by `6.5 mV` and `phip` too low by `9.6 mV`
  in the active avalanche support.
- This qF offset is sufficient to explain the local SG particle flux/source
  deficit, while potential, electric-field peak, mobility, and impact
  ionization coefficients are secondary.
- This is a diagnostic sensitivity, not a proposed hard-coded correction.

### Next Tasks After Task 94

1. Localize the origin of the `~6.5 mV / ~9.6 mV` high-field qF offset:
   - compare `psi-phin` and `phip-psi` along the overlap/false-negative active
     support for `-13.0`, `-13.1`, and `-13.2 V`;
   - determine whether the qF offset is uniform, contact-anchored, or created
     by a local continuity/source balance near edge `2886`.
2. Re-evaluate continuity residual terms on the active avalanche support using
   the qF-shifted state as a diagnostic target:
   - if the shifted state has lower residual under Vela equations, the solver is
     selecting the wrong nonlinear branch;
   - if the shifted state has higher residual under Vela equations but matches
     Sentaurus source/current, the residual equation/source placement differs
     from Sentaurus.
3. Inspect boundary/contact and intrinsic-density references only for their
   ability to create a `5..10 mV` qF offset:
   - contact qF floor/reporting should not be changed blindly because Task 94
     points to interior active support, not only terminal reporting;
   - BGN/`ni_eff` should be checked through its effect on
     `psi-phin`/`phip-psi` in the active support.
4. After the qF-offset origin is classified, decide the first code experiment:
   - continuation/branch control if the shifted qF state is a lower-residual
     Vela solution;
   - SG source/continuity residual parity if Vela equations reject the shifted
     state despite Sentaurus parity;
   - contact/intrinsic reference alignment only if the offset is shown to be a
     uniform boundary/reference shift.

### Execution Note 2026-06-20: Task 95 Active-Support qF Offset Origin

Task 95 followed the Task 94 qF-shift sensitivity by checking whether the
`~6.5 mV / ~9.6 mV` active-support qF offset is a global/contact reference
shift, a spatially local hotspot, or a Vela-equation residual preference.

qF anchor workspace:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13x_qf_anchor
```

Command:

```powershell
python scripts\diagnose_pn2d_bv_qf_anchor.py `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json `
  --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\focused_restart_m12p9078_m13p2\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\focused_restart_m13x_qf_anchor `
  --biases -13.0,-13.1,-13.2 `
  --focus-nodes 351,352,986,191,201,216,893
```

Contact qF anchor result:

| bias | contact | median `delta psi` | median `delta phin` | median `delta phip` |
| --- | --- | ---: | ---: | ---: |
| `-13.0 V` | Cathode | `+5.98e-9 V` | `~0 V` | `0 V` |
| `-13.0 V` | Anode | `-4.90e-5 V` | `~0 V` | `~0 V` |
| `-13.1 V` | Cathode | `+5.98e-9 V` | `~0 V` | `0 V` |
| `-13.1 V` | Anode | `-4.90e-5 V` | `~0 V` | `~0 V` |
| `-13.2 V` | Cathode | `+5.98e-9 V` | `~0 V` | `0 V` |
| `-13.2 V` | Anode | `-4.90e-5 V` | `~0 V` | `~0 V` |

Therefore the qF offset is not a simple global qF reference error: contact
Dirichlet qF values are already aligned to Sentaurus within numerical precision.
The qF anchor summary also reports that applying a uniform `phin` correction
would create a contact `phin` violation of `8.28..8.68 mV`.

Band-local qF offset:

| bias | band | median `delta psi` | median `delta phin` | median `delta phip` | median `delta(psi-phin)` |
| --- | --- | ---: | ---: | ---: | ---: |
| `-13.0 V` | junction | `~0 V` | `+6.05 mV` | `-8.71 mV` | `-6.44 mV` |
| `-13.1 V` | junction | `~0 V` | `+6.18 mV` | `-8.88 mV` | `-6.58 mV` |
| `-13.2 V` | junction | `~0 V` | `+6.32 mV` | `-9.10 mV` | `-6.73 mV` |
| `-13.2 V` | post-junction n | `+0.42 mV` | `+7.41 mV` | `-0.68 mV` | `-8.68 mV` |

This shows the active-support qF offset is an interior/junction transition
feature, not a contact-anchor mismatch.

Active-support continuity balance workspace:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p2_active_support_continuity_balance
```

At `-13.2 V`, active support summary:

| support class | count | median `delta(psi-phin)` | median `delta(phip-psi)` | Vela edge avalanche / Sentaurus generation | Vela electron density / Sentaurus | Vela hole density / Sentaurus |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| all | `30` | `-6.72 mV` | `-9.55 mV` | `0.3790` | `0.7711` | `0.6913` |
| overlap | `10` | `-6.58 mV` | `-9.62 mV` | `0.3803` | `0.7751` | `0.6894` |
| false-positive | `10` | `-6.61 mV` | `-9.51 mV` | `0.3800` | `0.7744` | `0.6922` |
| false-negative | `10` | `-6.99 mV` | `-9.50 mV` | `0.3763` | `0.7632` | `0.6925` |

The offset is spatially coherent across the active avalanche support classes,
not isolated to a single support mismatch category.

Active-support sensitivity workspace:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p2_active_support_sensitivity
```

Finite-difference sensitivity at `-13.2 V`, `delta = 1 mV`:

```text
electron dR/dphin median = -2.3622e8 s^-1/V
hole     dR/dphip median = +2.5246e8 s^-1/V
SRH derivative terms are negligible relative to transport/source derivatives.
```

Using the Task 94 parity qF shifts (`phin -= 6.5 mV`, `phip += 9.6 mV`),
a linearized residual estimate was generated under:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p2_active_support_qf_shift_residual_estimate
```

Summary:

| carrier | qF shift | median abs baseline residual | median abs predicted shifted residual | shifted/baseline | shifted/impact | predicted sign |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| electron | `-6.5 mV` | `6.9852e4 s^-1` | `1.4636e6 s^-1` | `20.49x` | `9.91x` | positive at all 30 nodes |
| hole | `+9.6 mV` | `7.4864e4 s^-1` | `2.3495e6 s^-1` | `31.17x` | `15.91x` | positive at all 30 nodes |

Interpretation:

- The Sentaurus-parity qF shift is not favored by the current Vela residual
  equations. It would substantially worsen the local active-support residual
  according to the finite-difference residual derivative.
- Therefore the remaining BV current/source discrepancy is unlikely to be fixed
  by Bank-Rose damping, line search, or another continuation-only change.
  Those controls may improve path stability, but the stable branch selected by
  Vela is already consistent with Vela's current equations.
- The likely root branch is now residual-equation/source-placement parity:
  Vela's continuity/source balance rejects the qF state that reproduces
  Sentaurus active-edge particle flux and avalanche source.

### Next Tasks After Task 95

1. Compare Vela active-support residual assembly against a Sentaurus-form
   residual proxy at the same active nodes:
   - reconstruct electron/hole SG transport integrals on active edges using the
     exact Sentaurus state;
   - combine with Sentaurus `ImpactIonization` and `srhRecombination`;
   - compare sign and magnitude to Vela's residual/source convention.
2. Split the residual parity check into one-factor substitutions:
   - Vela geometry + Sentaurus state + Sentaurus source;
   - Sentaurus-nearest geometry proxy + Sentaurus state + Sentaurus source;
   - Vela geometry + Vela state + Sentaurus source;
   - Vela geometry + shifted Vela qF + Vela/Sentaurus source.
3. If Sentaurus-form residual is balanced while Vela-form residual rejects the
   state, inspect source ownership/placement around active edge `2886`:
   - endpoint half-edge source assignment;
   - node-volume / endpoint-area conversion;
   - sign convention for impact generation in electron and hole continuity;
   - whether Sentaurus couples generated carriers to both carrier equations with
     the same local support as Vela.
4. Defer code changes until this residual parity split identifies one concrete
   term. The next code experiment should be one-factor and TDD-covered, not a
   global qF offset or damping knob.

### Execution Note 2026-06-20: Task 96 Active-Support Source/Residual Proxy Split

Task 96 continued the Task 95 residual-equation/source-placement branch. The
first check was whether the `~0.38x` active-support Vela/Sentaurus avalanche
source ratio came from a mechanical source ownership bug.

Source placement workspaces:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p2_sg_source_ownership
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p2_edge_direction_source_policy
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p2_source_geometry
```

Result:

- C++ edge source reconstructed back to support nodes matches the Vela VTK
  node avalanche source (`reconstructed_over_vtk ~= 1.0`).
- Endpoint area sums match node control volume on the active support
  (`cxx_endpoint_area_over_node_volume ~= 1.0`).
- The active endpoint area fraction is `0.5`, as expected for a single dominant
  junction-normal edge per active support node.
- The full node source remains `~0.38x` Sentaurus, while the active-edge average
  density is `~0.76x` Sentaurus.

Conclusion: no evidence of a half-edge, endpoint ownership, or node-volume
conversion bug. The source deficit is already present in the active-edge source
density/particle flux, and full-node averaging reduces it by the expected
`~0.5x` active-area fraction.

Residual proxy diagnostic added:

```text
scripts/diagnose_pn2d_bv_active_support_residual_proxy.py
```

The script compares active-support continuity residual proxies on the same Vela
mesh while substituting state and source policy one factor at a time:

- `vela_state_vela_edge_source`;
- `vela_state_sentaurus_node_source`;
- `sentaurus_state_sentaurus_node_source`;
- `sentaurus_state_replayed_edge_source`;
- `shifted_vela_qf_sentaurus_node_source`;
- `shifted_vela_qf_replayed_edge_source`.

Focused `-13.2 V` workspace:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p2_active_support_residual_proxy
```

Command:

```powershell
python scripts\diagnose_pn2d_bv_active_support_residual_proxy.py `
  --support-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\focused_restart_m13p2_thresholded_avalanche_support\thresholded_avalanche_support_nodes.csv `
  --sg-edge-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\focused_restart_m12p9078_m13p2\sg_avalanche_edges.csv `
  --mesh build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\mesh.json `
  --sentaurus-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\official_split_branch_drift_monitor\sentaurus_intermediate_exports\sentaurus_-13.2v `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\focused_restart_m12p9078_m13p2\vtk `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\sentaurus_default_bv_execution\focused_restart_m13p2_active_support_residual_proxy `
  --bias -13.2 `
  --electron-qf-shift-v -0.0065 `
  --hole-qf-shift-v 0.0096
```

Overlap active-support medians:

| variant | impact/Sentaurus node source | electron residual / impact | hole residual / impact | electron transport / Sentaurus source | hole transport / Sentaurus source |
| --- | ---: | ---: | ---: | ---: | ---: |
| `vela_state_vela_edge_source` | `0.3820` | `+0.5472` | `+0.4834` | `1.0681` | `1.0386` |
| `vela_state_sentaurus_node_source` | `1.0000` | `-0.4041` | `-0.4335` | `1.0681` | `1.0386` |
| `sentaurus_state_sentaurus_node_source` | `1.0000` | `+0.1435` | `+0.1713` | `1.6156` | `1.6435` |
| `sentaurus_state_replayed_edge_source` | `0.5047` | `+1.2659` | `+1.3210` | `1.6156` | `1.6435` |
| `shifted_vela_qf_sentaurus_node_source` | `1.0000` | `-0.0988` | `+0.0335` | `1.3735` | `1.5057` |
| `shifted_vela_qf_replayed_edge_source` | `0.5058` | `+0.7671` | `+1.0444` | `1.3735` | `1.5057` |

The false-positive and false-negative support classes show the same qualitative
ordering: Vela state plus Sentaurus node source is under-driven relative to the
source (`~ -0.42x` residual/impact), Vela edge source flips the residual positive
because the edge source is only `~0.38x` Sentaurus, and the qF-shifted Vela state
with Sentaurus node source is closest to balance.

Interpretation:

- Source placement/ownership is not the primary bug; it preserves the C++/VTK
  source and the expected endpoint area. It does, however, amplifies the
  local-source deficit because only half of the node control volume is on the
  dominant active edge.
- The qF-shifted Vela state remains the best Sentaurus-parity local state under
  a density-SG residual proxy, but the current Vela qF-model residual derivative
  from Task 95 rejects it. This narrows the discrepancy to the transport/residual
  formulation used by Vela's continuity equations, not merely the nonlinear
  damping path.
- The Sentaurus state on the Vela mesh with Sentaurus node source is still not
  perfectly balanced (`~+0.14..0.17` residual/impact), so a pure Vela-geometry
  density-SG proxy cannot be treated as an exact Sentaurus residual. It is still
  useful for one-factor ordering.

### Next Tasks After Task 96

1. Add a residual proxy mode that uses Vela's actual qF-variable-`ni_eff`
   transport formula with old-Slotboom `ni_eff`, alongside the density-SG proxy.
   This will directly compare:
   - density-SG transport;
   - qF-inferred-`ni_eff` transport;
   - qF-old-Slotboom-`ni_eff` transport currently used by Vela.
2. Run the new transport-mode split on `-13.0`, `-13.1`, and `-13.2 V` focused
   restarts. The target is to determine whether the qF-shifted state is rejected
   specifically by old-Slotboom `ni_eff` transport, by source coupling, or by a
   remaining geometry/proxy approximation.
3. If old-Slotboom qF transport is the rejecting term, compare Sentaurus exported
   or inferred `ni_eff` against Vela `ni_eff` on the active support and test a
   one-factor `ni_eff` substitution in the residual proxy before changing solver
   code.
4. If all qF transport modes reject the shifted state, inspect the residual sign
   and impact/SRH coupling convention in the C++ continuity assembly against
   DEVSIM/Charon-style residual assembly patterns.

### Execution Note 2026-06-20: Task 97 qF-Transport Residual Proxy Split

Task 97 implemented the first "Next Tasks After Task 96" item by extending:

```text
scripts/diagnose_pn2d_bv_active_support_residual_proxy.py
```

with explicit transport modes:

- `density_sg`;
- `qf_inferred_ni`;
- `qf_old_slotboom_ni`.

The default remains `density_sg` so older Task 96 reports keep their original
meaning. The qF modes are enabled explicitly with:

```powershell
--doping-csv build-release\reference_tcad\pn2d_sentaurus2018\reports\import_split_semantics_smoke\vela\doping.csv `
--transport-modes density_sg,qf_inferred_ni,qf_old_slotboom_ni
```

Focused workspaces:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p0_active_support_residual_proxy_transport_modes
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p1_active_support_residual_proxy_transport_modes
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p2_active_support_residual_proxy_transport_modes
```

Overlap active-support medians for the two most diagnostic variants:

| bias | variant | transport mode | electron residual / impact | hole residual / impact | electron transport / Sentaurus source | hole transport / Sentaurus source |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| `-13.0 V` | `vela_state_sentaurus_node_source` | `density_sg` | `-0.3885` | `-0.4147` | `1.1047` | `1.0783` |
| `-13.0 V` | `vela_state_sentaurus_node_source` | `qf_old_slotboom_ni` | `-0.7779` | `-0.8182` | `0.7152` | `0.6748` |
| `-13.0 V` | `shifted_vela_qf_sentaurus_node_source` | `density_sg` | `-0.0727` | `+0.0701` | `1.4205` | `1.5632` |
| `-13.0 V` | `shifted_vela_qf_sentaurus_node_source` | `qf_old_slotboom_ni` | `-0.5734` | `-0.5148` | `0.9197` | `0.9783` |
| `-13.1 V` | `vela_state_sentaurus_node_source` | `density_sg` | `-0.4012` | `-0.4191` | `1.0813` | `1.0635` |
| `-13.1 V` | `vela_state_sentaurus_node_source` | `qf_old_slotboom_ni` | `-0.7815` | `-0.8127` | `0.7010` | `0.6697` |
| `-13.1 V` | `shifted_vela_qf_sentaurus_node_source` | `density_sg` | `-0.0921` | `+0.0592` | `1.3905` | `1.5418` |
| `-13.1 V` | `shifted_vela_qf_sentaurus_node_source` | `qf_old_slotboom_ni` | `-0.5812` | `-0.5115` | `0.9014` | `0.9709` |
| `-13.2 V` | `vela_state_sentaurus_node_source` | `density_sg` | `-0.4041` | `-0.4335` | `1.0681` | `1.0386` |
| `-13.2 V` | `vela_state_sentaurus_node_source` | `qf_old_slotboom_ni` | `-0.7772` | `-0.8206` | `0.6950` | `0.6515` |
| `-13.2 V` | `shifted_vela_qf_sentaurus_node_source` | `density_sg` | `-0.0988` | `+0.0335` | `1.3735` | `1.5057` |
| `-13.2 V` | `shifted_vela_qf_sentaurus_node_source` | `qf_old_slotboom_ni` | `-0.5785` | `-0.5277` | `0.8937` | `0.9445` |

The `qf_inferred_ni` and `qf_old_slotboom_ni` results are nearly identical at
all three biases. Example at `-13.2 V`, overlap support:

```text
shifted_vela_qf_sentaurus_node_source:
  qf_inferred_ni     eT/S = 0.8942, hT/S = 0.9431
  qf_old_slotboom_ni eT/S = 0.8937, hT/S = 0.9445
```

Interpretation:

- The active-support residual split does not support `ni_eff`/OldSlotboom BGN
  as the primary high-field discrepancy. Old-Slotboom `ni_eff` and carrier-wise
  inferred `ni_eff` produce effectively the same qF transport residual proxy.
- The dominant split is now between the density-SG/current proxy used for
  avalanche-source parity and the qF-variable-`ni_eff` transport used by the
  continuity residual. The same qF-shifted state that nearly balances under the
  density-SG proxy is under-transported under qF transport:
  `electron transport / Sentaurus source ~= 0.89..0.92`,
  `hole transport / Sentaurus source ~= 0.94..0.98`.
- Therefore the immediate next root-cause question is not "which `ni_eff`
  formula is wrong?" but "does Vela use the same carrier current definition for
  impact generation that Sentaurus uses for its default BV avalanche source, and
  is that current definition consistent with the continuity residual?"

### Next Tasks After Task 97

1. Add an edge-level source-current consistency diagnostic for the active
   support:
   - for each active edge, output density-SG flux, qF-inferred-`ni` flux,
     qF-old-Slotboom-`ni` flux, C++ dumped source flux proxy, and Sentaurus
     exported `eCurrentDensity`/`hCurrentDensity` particle-flux proxy;
   - report which flux definition best predicts Sentaurus `ImpactIonization`
     and which definition best matches Vela continuity transport.
2. Inspect Vela C++ avalanche source assembly and continuity transport assembly
   together, not separately:
   - confirm whether `impact_ionization.current_approximation =
     "density_gradient"` is used only for source generation while qF transport is
     used in continuity residual;
   - compare this mixed-current setup with Sentaurus default BV semantics and
     Charon/DEVSIM residual/source coupling patterns.
3. If Sentaurus exported current aligns with qF transport while
   `ImpactIonization` aligns with density-SG source, document that Sentaurus
   likely uses different current proxies for avalanche and continuity; then
   Vela needs a consistent mixed-current calibration, not a single transport
   replacement.
4. If Sentaurus exported current aligns with density-SG instead, the next code
   experiment should be a TDD-covered option to use density-SG carrier transport
   in the high-field BV residual path, tested first on the focused
   `-13.0..-13.2 V` restart window before any full `-20 V` sweep.

### Execution Note 2026-06-20: Task 98 Active-Edge Source/Current Consistency

Task 98 implemented the first "Next Tasks After Task 97" item by adding:

```text
scripts/diagnose_pn2d_bv_edge_source_current_consistency.py
```

The diagnostic writes one row per active support endpoint edge and compares:

- C++ dumped SG avalanche edge flux/source proxies from `sg_avalanche_edges.csv`;
- Sentaurus exported `eCurrentDensity`/`hCurrentDensity` particle-current
  proxies on the same support nodes;
- Sentaurus `ImpactIonization` at the support nodes and active-edge endpoints;
- Vela and Sentaurus density-SG, qF-inferred-`ni`, and qF-old-Slotboom-`ni`
  reconstructed edge fluxes.

Focused workspaces:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p0_edge_source_current_consistency
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p1_edge_source_current_consistency
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p2_edge_source_current_consistency
```

Overlap active-edge endpoint medians:

| bias | C++ flux / Sentaurus current | Vela density flux / Sentaurus current | Vela qF-old-Slotboom flux / Sentaurus current | Sentaurus density flux / Sentaurus current | Sentaurus qF-old-Slotboom flux / Sentaurus current | C++ source / Sentaurus support source |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `-13.0 V` | `0.7503` | `0.7471` | `0.7468` | `1.0084` | `1.0077` | `0.7738` |
| `-13.1 V` | `0.7453` | `0.7420` | `0.7418` | `1.0084` | `1.0077` | `0.7689` |
| `-13.2 V` | `0.7401` | `0.7369` | `0.7367` | `1.0083` | `1.0076` | `0.7640` |

Additional overlap medians:

- `C++ flux / Vela density flux ~= 1.0043`;
- `C++ flux / Vela qF-old-Slotboom flux ~= 1.0047..1.0048`;
- `C++ source / Sentaurus edge-average source ~= 0.783..0.793`.

Interpretation:

- Sentaurus exported terminal/current-density fields are internally consistent
  with Sentaurus reconstructed density and qF edge flux on the active edges:
  both reconstructed ratios are `~1.008x` the exported-current particle proxy.
- Vela's C++ dumped edge flux is internally consistent with Vela reconstructed
  density and qF edge flux: the C++/diagnostic ratios are only `~1.004..1.005`.
- The cross-simulator deficit is therefore not a dump/proxy mismatch and not a
  density-vs-qF edge-current convention mismatch at these active edges. Vela's
  local active-edge carrier flux/current magnitude is `~0.737..0.747x`
  Sentaurus, and its active-edge source is `~0.764..0.774x` Sentaurus support
  generation.
- This refines Task 97: the earlier residual-proxy split is affected by the
  node-volume/source policy used in the residual proxy. At the edge-current
  level, each simulator is internally self-consistent; the remaining root is the
  local state/current magnitude selected by Vela on the active avalanche branch.

### Next Tasks After Task 98

1. Add an edge-level one-factor state-substitution diagnostic on the same active
   edges and same Sentaurus exported-current denominator:
   - Vela `psi` plus Sentaurus carrier/quasi-Fermi state;
   - Sentaurus `psi` plus Vela carrier/quasi-Fermi state;
   - shifted Vela quasi-Fermi state using the previously observed active-support
     offsets;
   - report edge-current ratios against Sentaurus exported current at
     `-13.0`, `-13.1`, and `-13.2 V`.
2. Decompose the remaining active-edge flux deficit into edge-local factors:
   carrier-density geometric mean, Bernoulli/electric-potential factor,
   mobility, and quasi-Fermi drop/driving field. Prefer extending the Task 98
   edge-current denominator so all factors compare to the same Sentaurus
   `eCurrentDensity`/`hCurrentDensity` reference.
3. If the qF-shifted Vela state reaches Sentaurus exported current but the
   qF-transport residual proxy rejects it, inspect continuity residual/source
   balance and source placement rather than changing the edge-current formula.
4. If Sentaurus state on the Vela mesh cannot reproduce Sentaurus exported
   current on these edges, inspect geometry/current projection before touching
   solver physics.

### Execution Note 2026-06-20: Task 99 Edge-State Substitution vs Exported Current

Task 99 implemented the first "Next Tasks After Task 98" item by adding:

```text
scripts/diagnose_pn2d_bv_edge_state_substitution_current.py
```

Unlike the older node-aggregated mixed-state replay, this diagnostic keeps the
Task 98 denominator fixed: Sentaurus exported `eCurrentDensity` plus
`hCurrentDensity` particle flux at the same active support node or active-edge
average. It writes one row per active support endpoint edge and per state
variant:

- `vela_baseline`;
- `sentaurus_baseline`;
- `vela_psi_sentaurus_qf`;
- `sentaurus_psi_vela_qf`;
- `vela_qf_shift`;
- `vela_qf_shift_sentaurus_mobility`.

Focused workspaces:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p0_edge_state_substitution_current
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p1_edge_state_substitution_current
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p2_edge_state_substitution_current
```

Overlap active-edge endpoint median
`density_particle_flux_over_sentaurus_support_current`:

| bias | Sentaurus state | Vela state | Vela `psi` + Sentaurus qF | Sentaurus `psi` + Vela qF | shifted Vela qF | shifted Vela qF + Sentaurus mobility |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `-13.0 V` | `1.0084` | `0.7471` | `1.0104` | `0.7597` | `1.0241` | `1.0300` |
| `-13.1 V` | `1.0084` | `0.7420` | `1.0104` | `0.7547` | `1.0173` | `1.0232` |
| `-13.2 V` | `1.0083` | `0.7369` | `1.0103` | `0.7494` | `1.0101` | `1.0159` |

The qF-old-Slotboom particle-flux ratios are essentially the same ordering:

| bias | Sentaurus state | Vela state | Vela `psi` + Sentaurus qF | Sentaurus `psi` + Vela qF | shifted Vela qF | shifted Vela qF + Sentaurus mobility |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `-13.0 V` | `1.0077` | `0.7468` | `1.0098` | `0.7595` | `1.0236` | `1.0294` |
| `-13.1 V` | `1.0077` | `0.7418` | `1.0099` | `0.7545` | `1.0168` | `1.0227` |
| `-13.2 V` | `1.0076` | `0.7367` | `1.0100` | `0.7492` | `1.0095` | `1.0153` |

Interpretation:

- The active-edge current deficit follows the quasi-Fermi/carrier branch, not
  the electrostatic potential field. Substituting Sentaurus qF/carrier state
  onto Vela `psi` restores the edge current to `~1.01x` Sentaurus exported
  current, while substituting Vela qF/carrier state onto Sentaurus `psi` remains
  low at `~0.75x`.
- The fixed `phin -= 6.5 mV`, `phip += 9.6 mV` Vela qF shift closes the edge
  current at `-13.2 V` (`~1.01x`) and slightly overshoots at `-13.0 V`
  (`~1.02x`). This is consistent with a bias-dependent qF branch offset rather
  than a constant missing mobility or geometry factor.
- Sentaurus mobility substitution only adds another `~0.6%` on these edge
  currents. Mobility is not the leading source of the remaining `~25%` Vela
  baseline deficit.
- Because the qF-shifted state now reproduces Sentaurus exported current but
  Task 97 showed the qF-transport residual proxy under-transports/rejects that
  state, the next debug branch should inspect the continuity residual/source
  balance that selects the Vela qF branch.

### Next Tasks After Task 99

1. Build a focused active-edge continuity balance for the qF-shifted state using
   the same exported-current denominator:
   - compare C++ continuity residual terms at `vela_baseline`,
     `vela_qf_shift`, and `vela_psi_sentaurus_qf`;
   - split electron and hole transport divergence, impact generation, SRH, and
     net residual at the active support nodes;
   - report whether the qF-shifted state is rejected because transport is too
     small, source is too large/small, or the source is placed on the wrong
     control volume.
2. Add a bias-dependent qF-offset table from `-13.0`, `-13.1`, and `-13.2 V`
   by solving for the qF shifts that make edge current exactly match Sentaurus
   exported current. Use this only as a diagnostic branch monitor, not as a
   solver calibration.
3. If the qF-shifted continuity balance is still rejected with the correct
   source/current pair, inspect the C++ continuity assembly signs and scaling
   against DEVSIM/Charon residual/source coupling, especially how avalanche
   generation enters electron and hole rows.

### Execution Note 2026-06-20: Task 100 qF-Shifted Continuity Balance

Task 100 implemented the first "Next Tasks After Task 99" item by adding:

```text
scripts/diagnose_pn2d_bv_qf_shift_continuity_balance.py
```

The diagnostic keeps the Task 99 exported-current denominator while adding
continuity-balance terms for three focused state variants:

- `vela_baseline`;
- `vela_qf_shift`;
- `vela_psi_sentaurus_qf`.

For each active support node it reports:

- active-edge particle flux over Sentaurus exported current;
- electron and hole transport divergence integrals;
- C++ dumped active-edge source;
- replayed active-edge source from the same state/current model;
- Sentaurus node source from `ImpactIonization * node_volume`;
- SRH term;
- residuals using C++ source, replayed edge source, and Sentaurus node source.

Focused workspaces:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p0_qf_shift_continuity_balance
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p1_qf_shift_continuity_balance
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p2_qf_shift_continuity_balance
```

Overlap active-support medians using `qf_old_slotboom_ni` transport:

| bias | variant | J/S current | C++ source/S source | replay source/S source | e transport/S source | h transport/S source | e residual/C++ source | h residual/C++ source | e residual/replay source | h residual/replay source | e residual/S source | h residual/S source |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `-13.0 V` | `vela_baseline` | `0.7469` | `0.3869` | `0.3852` | `0.7152` | `0.6748` | `-0.4353` | `-0.5364` | `-0.4327` | `-0.5343` | `-0.7779` | `-0.8182` |
| `-13.0 V` | `vela_qf_shift` | `1.0234` | `0.3869` | `0.5122` | `0.9197` | `0.9783` | `+0.0844` | `+0.2427` | `-0.1798` | `-0.0610` | `-0.5734` | `-0.5148` |
| `-13.0 V` | `vela_psi_sentaurus_qf` | `1.0092` | `0.3869` | `0.5454` | `1.3813` | `0.9652` | `+1.2786` | `+0.2111` | `+0.6097` | `-0.1443` | `-0.1119` | `-0.5280` |
| `-13.1 V` | `vela_baseline` | `0.7420` | `0.3846` | `0.3828` | `0.7010` | `0.6697` | `-0.4359` | `-0.5145` | `-0.4333` | `-0.5122` | `-0.7815` | `-0.8127` |
| `-13.1 V` | `vela_qf_shift` | `1.0166` | `0.3846` | `0.5091` | `0.9014` | `0.9709` | `+0.0792` | `+0.2648` | `-0.1839` | `-0.0440` | `-0.5812` | `-0.5115` |
| `-13.1 V` | `vela_psi_sentaurus_qf` | `1.0088` | `0.3846` | `0.5452` | `1.3770` | `0.9561` | `+1.3096` | `+0.2230` | `+0.6218` | `-0.1412` | `-0.1056` | `-0.5264` |
| `-13.2 V` | `vela_baseline` | `0.7367` | `0.3820` | `0.3803` | `0.6950` | `0.6515` | `-0.4263` | `-0.5307` | `-0.4236` | `-0.5286` | `-0.7772` | `-0.8206` |
| `-13.2 V` | `vela_qf_shift` | `1.0094` | `0.3820` | `0.5058` | `0.8937` | `0.9445` | `+0.0853` | `+0.2360` | `-0.1793` | `-0.0664` | `-0.5785` | `-0.5277` |
| `-13.2 V` | `vela_psi_sentaurus_qf` | `1.0084` | `0.3820` | `0.5451` | `1.3720` | `0.9523` | `+1.3383` | `+0.2476` | `+0.6321` | `-0.1291` | `-0.1002` | `-0.5200` |

Interpretation:

- Task 99's current-matching result holds: `vela_qf_shift` reaches
  `~1.01x` Sentaurus exported current at `-13.2 V`, and modestly overshoots at
  `-13.0..-13.1 V`.
- With Sentaurus node source, the same `vela_qf_shift` state is still strongly
  under-balanced in the residual proxy: overlap medians are about
  `e residual/S source = -0.58`, `h residual/S source = -0.52`.
- With replayed active-edge source, `vela_qf_shift` is much closer:
  electron residual is about `-0.18x` replayed source and hole residual is about
  `-0.04..-0.07x`. This suggests the next discrepancy is source placement /
  active-edge source ownership, not merely edge-current reconstruction.
- With the original C++ dumped source, `vela_qf_shift` flips to a small positive
  residual (`e ~= +0.08x`, `h ~= +0.24x` C++ source). This is consistent with
  the C++ active-edge source being only `~0.382..0.387x` Sentaurus node source
  in node-integral form, while the replayed shifted edge source is
  `~0.506..0.512x`.
- The accepted `vela_baseline` is not balanced by this Python proxy either
  (`~ -0.43..-0.53x` C++ source), so the proxy must not be treated as the exact
  C++ residual. It is useful for bracketing source-policy effects, but the next
  decisive check needs the C++ residual evaluator on external shifted states.

### Next Tasks After Task 100

1. Add or extend an exact C++ external-state residual probe for the three Task
   100 variants:
   - `vela_baseline`;
   - `vela_qf_shift`;
   - `vela_psi_sentaurus_qf`;
   - report exact electron/hole continuity residual terms on the same overlap
     active-support nodes.
2. Use the exact C++ probe to answer one binary question: does the qF-shifted
   state have a near-zero C++ continuity residual when paired with C++ active
   edge source, or does exact assembly still reject it?
3. If exact C++ residual is near zero for qF-shifted state, the remaining BV
   mismatch is branch selection/continuation: add a diagnostic branch guard or
   predictor that steers toward the Sentaurus qF branch.
4. If exact C++ residual rejects qF-shifted state, inspect source assembly
   signs/scaling and edge-to-node ownership against Charon/DEVSIM patterns
   before changing solver damping or mobility.

### Execution Note 2026-06-20: Task 101 Exact C++ Carrier-Term Probe

Task 101 implemented the first two "Next Tasks After Task 100" items by adding:

```text
scripts/diagnose_pn2d_bv_exact_carrier_term_states.py
```

This wrapper uses the existing C++ `vela_example_runner` capability:

```text
simulation_type = "newton_carrier_term_probe"
```

It prepares external state fields for:

- `vela_baseline`;
- `vela_qf_shift`;
- `vela_psi_sentaurus_qf`;

then runs the exact C++ carrier-term diagnostic and extracts the same active
support rows used in Tasks 98-100. Unlike the Python residual proxy, this is the
actual C++ continuity assembly: SG flux, recombination, impact generation,
gauge/boundary terms, and row residual.

Focused workspaces:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p0_exact_carrier_term_states
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p1_exact_carrier_term_states
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p2_exact_carrier_term_states
```

Overlap active-support medians:

| bias | variant | phin block | phip block | e residual / abs(impact) | h residual / abs(impact) | e flux / abs(impact) | h flux / abs(impact) | recomb / abs(impact) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `-13.0 V` | `vela_baseline` | `4.20e-14` | `2.67e-13` | `+0.0152` | `+0.0059` | `2.2967` | `2.2774` | `-1.2784` |
| `-13.0 V` | `vela_qf_shift` | `1.4661` | `2.1653` | `+0.2565` | `+0.5187` | `2.2209` | `2.4828` | `-0.9614` |
| `-13.0 V` | `vela_psi_sentaurus_qf` | `1.02e-11` | `4.79e-12` | `+1.0215` | `+0.3500` | `2.9836` | `2.3421` | `-0.9432` |
| `-13.1 V` | `vela_baseline` | `4.77e-14` | `2.69e-13` | `-0.0005` | `+0.0175` | `2.2628` | `2.2855` | `-1.2585` |
| `-13.1 V` | `vela_qf_shift` | `1.4661` | `2.1653` | `+0.2371` | `+0.5360` | `2.1879` | `2.4913` | `-0.9463` |
| `-13.1 V` | `vela_psi_sentaurus_qf` | `1.29e-11` | `5.76e-12` | `+1.0267` | `+0.3509` | `2.9685` | `2.3227` | `-0.9231` |
| `-13.2 V` | `vela_baseline` | `4.85e-14` | `2.70e-13` | `+0.0038` | `-0.0018` | `2.2370` | `2.2457` | `-1.2397` |
| `-13.2 V` | `vela_qf_shift` | `1.4661` | `2.1653` | `+0.2314` | `+0.5082` | `2.1626` | `2.4476` | `-0.9321` |
| `-13.2 V` | `vela_psi_sentaurus_qf` | `1.58e-11` | `7.14e-12` | `+1.0211` | `+0.3578` | `2.9433` | `2.3098` | `-0.9035` |

Interpretation:

- The exact C++ probe overturns the branch-selection-only hypothesis. The
  accepted `vela_baseline` has essentially zero active-support carrier residual
  in exact C++ assembly (`phin/phip block ~= 1e-13`, residual/impact near zero).
- The Sentaurus-current-matching `vela_qf_shift` state is clearly rejected by
  exact C++ assembly: active-support residuals are about
  `+0.23..0.26x` electron impact and `+0.51..0.54x` hole impact, with global
  `phin/phip` block residuals `1.466/2.165`.
- `vela_psi_sentaurus_qf` is even more strongly rejected, especially in the
  electron row (`~+1.02x` impact), despite matching the exported active-edge
  current in Task 99.
- Therefore the current root is not "Newton damping failed to select a valid
  Sentaurus-like Vela branch." The Sentaurus-like qF/current branch does not
  satisfy Vela's exact continuity equations under the current recombination,
  impact-source, and carrier transport assembly.
- The strongest numerical discriminator is recombination/impact/flux balance:
  baseline needs recombination about `-1.24..-1.28x` impact to cancel flux,
  whereas qF-shifted/Sentaurus-like states have recombination only
  `-0.90..-0.96x` impact and leave positive residuals.

### Next Tasks After Task 101

1. Decompose why exact C++ recombination changes so much between
   `vela_baseline` and Sentaurus-like qF states:
   - export active-support electron/hole densities and SRH denominators from
     the exact C++ state, or reproduce the exact SRH formula in a Python
     comparator using the generated external state fields;
   - compare `n`, `p`, `ni_eff`, `n1`, `p1`, and lifetime terms on overlap
     active nodes.
2. Run a controlled exact C++ carrier-term probe with SRH disabled or lifetime
   scaled on the same three external states. If qF-shift residual drops toward
   zero when SRH is removed/scaled, the BV mismatch is tied to Sentaurus SRH
   parity under high-bias qF states rather than avalanche current alone.
3. Compare Sentaurus exported `srhRecombination` against Vela exact
   recombination for the same `vela_psi_sentaurus_qf` and `vela_qf_shift`
   states if a Sentaurus variant can be generated through the VM. Without a
   Sentaurus external-state residual, use the manual/par-file SRH defaults to
   reproduce Sentaurus SRH on the exported state.
4. Only after recombination parity is understood, revisit impact-source
   ownership. The exact probe shows impact/flux current matching is insufficient
   if recombination balance selects a different qF/carrier branch.

### Execution Note 2026-06-20: Task 102 Exact SRH Decomposition and SRH-Off Control

Task 102 extended:

```text
scripts/diagnose_pn2d_bv_exact_carrier_term_states.py
```

with an exact SRH decomposition for the active-support rows:

- `n = ni_eff * exp((psi - phin) / Vt)`;
- `p = ni_eff * exp((phip - psi) / Vt)`;
- `np - ni_eff^2 = ni_eff^2 * expm1((phip - phin) / Vt)`;
- `R_SRH = (np - ni_eff^2) / (taup * (n + ni_eff) + taun * (p + ni_eff))`.

The decomposition also reports `n/ni_eff`, `p/ni_eff`, the two lifetime
denominator terms, `R * node_volume`, and the inferred C++ continuity scaling
from the exact recombination column. The focused exact C++ outputs were
regenerated for `-13.0/-13.1/-13.2 V`.

Overlap active-support medians:

| bias | variant | n / ni | p / ni | excess / ni^2 | taup term / baseline | taun term / baseline | SRH rate / baseline | recomb / abs(impact) | e residual / abs(impact) | h residual / abs(impact) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `-13.0 V` | `vela_baseline` | `1.85e-7` | `2.51e-7` | `-1.0000` | `1.0000` | `1.0000` | `1.0000` | `-1.2742` | `+0.0152` | `+0.0059` |
| `-13.0 V` | `vela_qf_shift` | `2.38e-7` | `3.63e-7` | `-1.0000` | `1.0000` | `1.0000` | `1.0000` | `-0.9583` | `+0.2565` | `+0.5187` |
| `-13.0 V` | `vela_psi_sentaurus_qf` | `2.79e-7` | `3.02e-7` | `-1.0000` | `1.0000` | `1.0000` | `1.0000` | `-0.8994` | `+1.0215` | `+0.3500` |
| `-13.1 V` | `vela_baseline` | `1.86e-7` | `2.52e-7` | `-1.0000` | `1.0000` | `1.0000` | `1.0000` | `-1.2544` | `-0.0005` | `+0.0175` |
| `-13.1 V` | `vela_qf_shift` | `2.39e-7` | `3.65e-7` | `-1.0000` | `1.0000` | `1.0000` | `1.0000` | `-0.9433` | `+0.2371` | `+0.5360` |
| `-13.1 V` | `vela_psi_sentaurus_qf` | `2.82e-7` | `3.06e-7` | `-1.0000` | `1.0000` | `1.0000` | `1.0000` | `-0.8803` | `+1.0267` | `+0.3509` |
| `-13.2 V` | `vela_baseline` | `1.87e-7` | `2.53e-7` | `-1.0000` | `1.0000` | `1.0000` | `1.0000` | `-1.2358` | `+0.0038` | `-0.0018` |
| `-13.2 V` | `vela_qf_shift` | `2.40e-7` | `3.67e-7` | `-1.0000` | `1.0000` | `1.0000` | `1.0000` | `-0.9292` | `+0.2314` | `+0.5082` |
| `-13.2 V` | `vela_psi_sentaurus_qf` | `2.85e-7` | `3.10e-7` | `-1.0000` | `1.0000` | `1.0000` | `1.0000` | `-0.8618` | `+1.0211` | `+0.3578` |

Interpretation:

- The active overlap nodes are deep in depletion: `n/ni_eff` and `p/ni_eff`
  are only `~1e-7`.
- Because `np << ni_eff^2`, SRH is pinned near the net-generation plateau:
  `np - ni_eff^2 ~= -ni_eff^2`, and the denominator is dominated by the
  `+ni_eff` trap terms rather than by `n` or `p`.
- This explains why the qF-shift/Sentaurus-like variants change carrier
  densities but do not materially change the exact SRH term. The weaker
  `recomb / abs(impact)` ratio in Task 101 is primarily caused by larger
  impact/source terms in those variants, not by reduced SRH.
- The inferred continuity scaling from `R * node_volume / exact_recombination`
  is stable at `3.490027e20`, confirming that the Python SRH decomposition
  reproduces the exact C++ recombination column up to the known scaled residual
  units.

A controlled SRH-off exact probe was then run with a derived config:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m12p9078_m13p2/simulation_recomb_none.json
```

and outputs:

```text
focused_restart_m13p0_exact_carrier_term_states_recomb_none
focused_restart_m13p1_exact_carrier_term_states_recomb_none
focused_restart_m13p2_exact_carrier_term_states_recomb_none
```

SRH-off overlap residuals:

| bias | variant | recomb model | e residual / abs(impact) | h residual / abs(impact) |
| --- | --- | --- | ---: | ---: |
| `-13.0 V` | `vela_baseline` | `srh` | `+0.015` | `+0.006` |
| `-13.0 V` | `vela_baseline` | `none` | `+1.29` | `+1.27` |
| `-13.0 V` | `vela_qf_shift` | `srh` | `+0.257` | `+0.519` |
| `-13.0 V` | `vela_qf_shift` | `none` | `+1.22` | `+1.48` |
| `-13.1 V` | `vela_baseline` | `srh` | `-0.0005` | `+0.018` |
| `-13.1 V` | `vela_baseline` | `none` | `+1.26` | `+1.28` |
| `-13.1 V` | `vela_qf_shift` | `srh` | `+0.237` | `+0.536` |
| `-13.1 V` | `vela_qf_shift` | `none` | `+1.18` | `+1.48` |
| `-13.2 V` | `vela_baseline` | `srh` | `+0.004` | `-0.002` |
| `-13.2 V` | `vela_baseline` | `none` | `+1.24` | `+1.24` |
| `-13.2 V` | `vela_qf_shift` | `srh` | `+0.231` | `+0.508` |
| `-13.2 V` | `vela_qf_shift` | `none` | `+1.16` | `+1.44` |

This rejects the hypothesis that SRH parity is the cause of the qF-shift state
being rejected by Vela. SRH removal makes both the accepted Vela baseline and
the qF-shift/Sentaurus-like states less balanced. SRH is a required cancelling
term in Vela's exact continuity balance; it is not the missing mechanism that
would accept the Sentaurus-current-matching qF branch.

### Next Tasks After Task 102

1. Return to impact-source and SG-flux balance. For the exact C++ probe rows,
   decompose the impact source into edge-current magnitude, driving force,
   alpha, endpoint/node ownership, and electron/hole row contribution.
2. Compare the qF-shift state against Sentaurus exported avalanche generation
   and current density on the same active overlap nodes. Quantify the source
   multiplier needed to close the exact C++ residual separately for electron
   and hole rows.
3. Add a controlled exact carrier-term probe variant that scales only the
   impact source or replays a Sentaurus/edge-derived source, without changing
   SG flux or SRH. This will answer whether the remaining residual can be
   closed by source magnitude/ownership alone.
4. Re-read Charon/DEVSIM source assembly patterns for avalanche generation
   ownership after the above numeric multiplier is known. The key question is
   no longer SRH; it is whether Vela assigns density-gradient impact generation
   to nodes/carrier equations in the same way as Sentaurus.

### Execution Note 2026-06-20: Task 103 Exact Impact Closure Multiplier

Task 103 extended the exact carrier-term diagnostic with impact-closure
columns:

```text
electron_required_impact_multiplier
hole_required_impact_multiplier
electron_impact_multiplier_delta
hole_impact_multiplier_delta
electron_residual_over_abs_impact
hole_residual_over_abs_impact
electron_closed_residual
hole_closed_residual
```

For each carrier row, the required multiplier is computed from the exact C++
terms:

```text
k_required = -(flux + recombination) / impact
```

so that:

```text
flux + recombination + k_required * impact = 0
```

This is a controlled diagnostic only: it answers how much the existing exact
impact term would have to change to close the row if SG flux and SRH were held
fixed.

Focused exact outputs were regenerated for `-13.0/-13.1/-13.2 V`.

Overlap active-support medians:

| bias | variant | electron k_required | hole k_required | electron delta | hole delta | hole/electron k ratio |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `-13.0 V` | `vela_baseline` | `1.015` | `1.009` | `+0.015` | `+0.009` | `0.993` |
| `-13.0 V` | `vela_qf_shift` | `1.256` | `1.528` | `+0.256` | `+0.528` | `1.217` |
| `-13.0 V` | `vela_psi_sentaurus_qf` | `1.999` | `1.334` | `+0.999` | `+0.334` | `0.667` |
| `-13.1 V` | `vela_baseline` | `1.000` | `1.021` | `-0.0005` | `+0.021` | `1.021` |
| `-13.1 V` | `vela_qf_shift` | `1.235` | `1.534` | `+0.235` | `+0.534` | `1.242` |
| `-13.1 V` | `vela_psi_sentaurus_qf` | `2.004` | `1.335` | `+1.004` | `+0.335` | `0.666` |
| `-13.2 V` | `vela_baseline` | `1.004` | `0.998` | `+0.004` | `-0.002` | `0.994` |
| `-13.2 V` | `vela_qf_shift` | `1.236` | `1.505` | `+0.236` | `+0.505` | `1.217` |
| `-13.2 V` | `vela_psi_sentaurus_qf` | `1.999` | `1.342` | `+0.999` | `+0.342` | `0.671` |

Comparison against the Task 100 source-bracketing table for `vela_qf_shift`:

| bias | exact electron k_required | exact hole k_required | C++ edge source / Sentaurus node source | replayed edge source / Sentaurus node source |
| --- | ---: | ---: | ---: | ---: |
| `-13.0 V` | `1.256` | `1.528` | `0.3869` | `0.5122` |
| `-13.1 V` | `1.235` | `1.534` | `0.3846` | `0.5091` |
| `-13.2 V` | `1.236` | `1.505` | `0.3820` | `0.5058` |

Interpretation:

- The accepted Vela baseline is self-consistent: `k_required ~= 1` for both
  carrier rows.
- The Sentaurus-current-matching `vela_qf_shift` state cannot be fixed by a
  single global impact magnitude multiplier. Electron rows require about
  `1.24..1.26x` the current exact impact term, while hole rows require about
  `1.50..1.53x`.
- `vela_psi_sentaurus_qf` is even more split: electron rows need nearly `2x`
  impact, while hole rows need only `~1.34x`.
- Therefore the leading suspect is not just avalanche coefficient magnitude.
  The mismatch now points to carrier-row split / edge-to-node ownership / the
  current branch used to form generation. A uniform scalar applied to the Vela
  impact term would leave one carrier row under- or over-balanced.

### Next Tasks After Task 103

1. Add an exact C++ carrier-term probe variant that can evaluate hypothetical
   impact source policies without changing SG flux or SRH:
   - separate electron and hole impact multipliers;
   - optional edge-to-node ownership policy (`half_edge`, `node0_node1_dump`,
     or replayed explicit node source);
   - emit the same carrier-term rows and block residuals.
2. First run the minimal two-multiplier control using the medians above:
   `electron_impact_scale ~= 1.24`, `hole_impact_scale ~= 1.51` for
   `vela_qf_shift` at `-13.2 V`. If this closes active-support rows but not
   global `phin/phip`, inspect where non-active nodes disagree.
3. Then test a source-ownership policy rather than independent carrier scales.
   A physically plausible Sentaurus parity fix should explain the electron/hole
   split through edge current branch ownership or node ownership, not through
   arbitrary carrier-specific constants.
4. After the source-policy control identifies a candidate, compare the exact
   implementation against Charon/DEVSIM avalanche assembly patterns before
   promoting any solver change.

### Execution Note 2026-06-20: Task 104 Contact-Preserving qF Shift Control

Task 104 added two diagnostic capabilities:

- `newton_carrier_term_probe` now accepts diagnostic-only impact scales:

```json
"carrier_term_probe": {
  "electron_impact_scale": 1.23627,
  "hole_impact_scale": 1.50455
}
```

The probe preserves the original exact terms and adds adjusted impact/residual
columns plus `adjusted_block_residuals`. This does not change the production
solver; it is a carrier-term what-if probe.

- `scripts/diagnose_pn2d_bv_exact_carrier_term_states.py` now forwards those
  scales and can preserve contact quasi-Fermi values during `vela_qf_shift`:

```text
--preserve-contact-qf-on-shift
```

The immediate reason was a boundary-condition artifact in Task 103: the
uniform qF shift was applied to contact nodes too, so `vela_qf_shift` violated
the Dirichlet quasi-Fermi contact rows. For `-13.2 V`, preserving contact qF
reduced the global carrier blocks:

| state | phin block | phip block |
| --- | ---: | ---: |
| uniform qF shift, contacts shifted | `1.466` | `2.165` |
| qF shift with contacts preserved | `0.550` | `0.395` |

Then the Task 103 active-support impact scales were applied:

```text
electron_impact_scale = 1.23627
hole_impact_scale     = 1.50455
```

For `vela_qf_shift` at `-13.2 V`, overlap active-support medians were driven
essentially to zero:

| row | raw residual median | adjusted residual median |
| --- | ---: | ---: |
| electron overlap | `+1.31e-16` | `+2.02e-20` |
| hole overlap | `+2.88e-16` | `-3.03e-20` |

However the global carrier blocks remained:

| state | adjusted phin block | adjusted phip block |
| --- | ---: | ---: |
| contact-preserved qF shift + impact scales | `0.550` | `0.395` |

Top adjusted residual rows show why:

- electron residual is concentrated at the right-contact-adjacent column
  (`x ~= 1.96875e-6 m`) with `~1.40e-1` scaled residual, dominated by SG flux
  and with negligible impact/recombination;
- hole residual is concentrated at the left-contact-adjacent column
  (`x ~= 3.125e-8 m`) with `~1.00e-1` scaled residual, also dominated by SG
  flux and with negligible impact/recombination.

Interpretation:

- Preserving contacts confirms that a large part of the old qF-shift global
  residual was an artificial contact-BC mismatch.
- But a uniform interior qF shift with fixed contacts is also not a physically
  valid branch: it creates an abrupt quasi-Fermi boundary layer next to the
  contacts and large SG flux residuals outside the avalanche support.
- Therefore the qF-shift current-matching experiment is useful only as a local
  active-support sensitivity probe. It should not be treated as a globally
  valid Sentaurus-like state.
- The next branch-localization experiment must use a spatially consistent qF
  shape: either the actual Sentaurus qF field, or a smooth/contact-anchored qF
  perturbation fitted to Sentaurus-Vela qF differences, not a uniform offset.

### Next Tasks After Task 104

1. Decompose `vela_psi_sentaurus_qf` more carefully. Its global carrier block
   is small, but active-support residual/impact is large; this suggests the
   mismatch is localized and physically meaningful, unlike the uniform qF
   shift boundary artifact.
2. Build a smooth/contact-anchored qF perturbation:
   - zero shift on contact nodes;
   - fit the interior qF difference between Vela and Sentaurus, or solve a
     Laplace/smoothed extension from active-support target shifts;
   - rerun exact carrier-term probe and impact-scale closure.
3. For `vela_psi_sentaurus_qf`, run the same impact closure with the exact
   required multipliers (`~2.0x` electron, `~1.34x` hole at `-13.2 V`) and
   inspect whether active-support rows close without creating contact-adjacent
   SG flux residuals.
4. Only if a spatially consistent qF state can be made exact-balanced with a
   plausible source policy should solver behavior be changed. Otherwise the
   remaining difference is a transport/current branch mismatch, not merely
   avalanche source magnitude.

### Execution Note 2026-06-20: Task 105 Sentaurus-qF Exact Carrier-Term Control

Task 105 reran the exact C++ carrier-term probe on the spatially consistent
`vela_psi_sentaurus_qf` state for `-13.0`, `-13.1`, and `-13.2 V`, using the
Task 103 closure multipliers:

| bias | electron impact scale | hole impact scale |
| ---: | ---: | ---: |
| `-13.0 V` | `1.99927` | `1.33405` |
| `-13.1 V` | `2.00418` | `1.33498` |
| `-13.2 V` | `1.99930` | `1.34160` |

For the actual Sentaurus qF shape, active overlap residuals close without the
contact-boundary artifact seen in the uniform qF-shift probe:

| bias | raw phin block | raw phip block | overlap e raw | overlap e adjusted | overlap h raw | overlap h adjusted |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `-13.0 V` | `1.024e-11` | `4.790e-12` | `5.718e-16` | `1.921e-18` | `1.959e-16` | `-2.089e-19` |
| `-13.1 V` | `1.289e-11` | `5.762e-12` | `5.872e-16` | `2.124e-18` | `2.007e-16` | `-2.129e-19` |
| `-13.2 V` | `1.585e-11` | `7.135e-12` | `5.967e-16` | `1.858e-18` | `2.091e-16` | `-2.192e-19` |

Interpretation update:

- The earlier statement that `vela_psi_sentaurus_qf` is "strongly rejected"
  must be read locally, as residual divided by the local impact source. In
  absolute global carrier-block terms it is already small (`~1e-11`), unlike
  the uniform qF-shift state.
- The remaining mismatch is therefore localized in the avalanche-support
  carrier rows, not a global quasi-Fermi shape or contact-BC failure.
- The relevant unknown is the source policy/current branch used in the carrier
  rows, because the exact Sentaurus qF field is globally plausible but needs
  different effective impact strengths in the electron and hole equations.

### Execution Note 2026-06-20: Task 106 Impact Source Component Probe

Task 106 added diagnostic-only source decomposition:

- `detail::sgEdgeCurrentAvalancheSourceComponentIntegrals(...)` returns nodal
  SG avalanche source components:
  - `electron`: from `alpha_n * |electron SG flux|`;
  - `hole`: from `alpha_p * |hole SG flux|`;
  - `combined`: the existing `electron + hole` source.
- `CoupledDDCarrierTermDiagnostic` now carries
  `impactElectronSource`, `impactHoleSource`, and `impactCombinedSource`.
- `newton_carrier_term_probe` CSV now emits:
  `impact_electron_source`, `impact_hole_source`, and
  `impact_combined_source`.
- `scripts/diagnose_pn2d_bv_exact_carrier_term_states.py` includes these
  fields in `exact_carrier_term_state_nodes.csv`.

Validation:

```text
cmake --build build-release --target test_impact_ionization vela_example_runner -j 4
build-release/test_impact_ionization.exe
python -m unittest \
  tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_exact_carrier_term_states_exports_impact_source_components \
  tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_runner_writes_newton_carrier_term_probe_for_external_state
```

All passed. `ctest -R impact_ionization` did not match a registered test name,
so the executable was run directly.

For `vela_psi_sentaurus_qf` at `-13.2 V`, the active overlap support has:

| class | count | median electron source fraction | median hole source fraction | median combined source |
| --- | ---: | ---: | ---: | ---: |
| `overlap` | `10` | `0.77141` | `0.22859` | `5.84348e-16` |

Policy residual estimates on the same overlap nodes, normalized by the current
combined source:

| hypothetical policy | electron median residual | hole median residual |
| --- | ---: | ---: |
| current Vela: `combined` to both rows | `+0.9993` | `+0.3416` |
| `2 * combined` to both rows | `-0.0007` | `-0.6584` |
| component-only matching row | `+1.2548` | `+1.1130` |
| hole row: `combined + hole_component` | n/a | `+0.1130` |
| hole row: `combined + electron_component` | n/a | `-0.4298` |

This rules out the simplest "electron component goes only to electron row,
hole component goes only to hole row" hypothesis: it would make both carrier
rows worse. The electron row behaves as if it is missing nearly one additional
copy of the combined impact source. The hole row behaves differently; it is
closer to `combined + hole_component` but still under-balanced by about
`0.11 * combined` at the median.

### Next Tasks After Task 106

1. Add a probe-only source-policy evaluator using the three component columns,
   without changing production residual assembly:
   - `electron_policy`: `combined`, `double_combined`, `combined_plus_electron`,
     `combined_plus_hole`, `electron_only`, `hole_only`;
   - `hole_policy`: same choices;
   - report active-support and global block residuals for each policy.
2. Run the policy matrix on `vela_psi_sentaurus_qf` for `-13.0`, `-13.1`, and
   `-13.2 V`. A viable Sentaurus-parity policy should keep global blocks small
   while reducing overlap residuals across all three biases.
3. Compare the winning policy against Sentaurus/Charon/DEVSIM semantics:
   - whether carrier rows receive total pair-generation or branch-specific
     ionization source;
   - whether edge direction, contact ownership, or control-volume ownership
     changes the source seen by each row;
   - whether the electron/hole asymmetry tracks current direction rather than
     carrier type.
4. Only after a physically explainable policy is identified should production
   assembly be changed. The current evidence does not support a simple scalar
   avalanche multiplier or component-only split.

### Execution Note 2026-06-20: Task 107 Carrier-Row Source Policy Matrix

Task 107 added a probe-only postprocessor:

```text
scripts/diagnose_pn2d_bv_source_policy_matrix.py
```

Inputs:

- `exact_carrier_term_state_nodes.csv` for the selected support nodes;
- `carrier_terms/vela_psi_sentaurus_qf_carrier_terms.csv` for full carrier rows
  and, when needed, source-component columns backfilled by `node_id`;
- `variant = vela_psi_sentaurus_qf`;
- bias filter.

The policy set is:

```text
combined
double_combined
combined_plus_electron
combined_plus_hole
electron_only
hole_only
zero
```

For each electron/hole policy pair, the script recomputes:

```text
carrier residual = flux + recombination + gauge + boundary - policy_source
```

and reports support-class and full-row block metrics. The script includes a
Windows long-path write helper because the focused report path reached the
classic 260-character limit.

Real runs:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p0_exact_carrier_term_states_sentaurus_qf_impact_scaled/source_policy_matrix_vela_psi_sentaurus_qf
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p1_exact_carrier_term_states_sentaurus_qf_impact_scaled/source_policy_matrix_vela_psi_sentaurus_qf
build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/focused_restart_m13p2_exact_carrier_term_states_sentaurus_qf_impact_scaled/source_policy_matrix_vela_psi_sentaurus_qf
```

Active overlap result, current Vela policy versus best matrix policy:

| bias | current e median / combined | current h median / combined | current carrier L2 | best e policy | best h policy | best e signed median | best h signed median | best carrier L2 | best e abs median | best h abs median |
| ---: | ---: | ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| `-13.0 V` | `+0.9993` | `+0.3340` | `2.3516e-15` | `double_combined` | `combined_plus_hole` | `-0.0007` | `+0.1078` | `1.0196e-15` | `0.2498` | `0.1488` |
| `-13.1 V` | `+1.0042` | `+0.3350` | `2.4115e-15` | `double_combined` | `combined_plus_hole` | `+0.0042` | `+0.1076` | `1.0452e-15` | `0.2505` | `0.1432` |
| `-13.2 V` | `+0.9993` | `+0.3416` | `2.4620e-15` | `double_combined` | `combined_plus_hole` | `-0.0007` | `+0.1130` | `1.0612e-15` | `0.2388` | `0.1442` |

For `-13.2 V`, the top policy ranks in active overlap are:

| rank | electron policy | hole policy | carrier L2 | e signed median | h signed median | e abs median | h abs median |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | `double_combined` | `combined_plus_hole` | `1.0612e-15` | `-0.0007` | `+0.1130` | `0.2388` | `0.1442` |
| 2 | `double_combined` | `combined_plus_electron` | `1.1042e-15` | `-0.0007` | `-0.4298` | `0.2388` | `0.5466` |
| 3 | `combined_plus_electron` | `combined_plus_hole` | `1.2171e-15` | `+0.2548` | `+0.1130` | `0.2548` | `0.1442` |
| 4 | `double_combined` | `double_combined` | `1.2512e-15` | `-0.0007` | `-0.6584` | `0.2388` | `0.6584` |

Interpretation:

- The best active-overlap policy is stable across `-13.0/-13.1/-13.2 V`:
  electron row wants `double_combined`; hole row wants
  `combined_plus_hole`.
- This improves the active-overlap carrier L2 by about `2.3x`, but it does not
  close the rows. The best electron absolute median remains
  `~0.24-0.25 * combined`, and the best hole absolute median remains
  `~0.14-0.15 * combined`.
- Therefore a fixed carrier-row source policy is not sufficient. It explains
  the sign and average split, but not the node-to-node spatial variation.
- The remaining mismatch is more likely an edge-local ownership/current-branch
  weighting problem: the effective source seen by a row varies across active
  nodes, whereas the row-policy matrix applies the same algebraic combination
  everywhere.
- The full global matrix should not be used as the main discriminator here:
  global blocks stay at `~1e-11` and are dominated by tiny flux residuals far
  outside the avalanche support. Normalized global ratios are ill-conditioned
  where the local combined source is near zero.

### Next Tasks After Task 107

1. Build a node-local coefficient diagnostic on the active-overlap rows:
   fit `required_source = a * electron_component + b * hole_component` per
   carrier row and report `a,b` per node, per bias.
2. Correlate fitted `a,b` against incident active-edge direction, edge axis,
   current direction, contact side, and existing `edge_direction_source_policy`
   features. The goal is to determine whether the coefficient variation is
   geometric/ownership-driven rather than arbitrary.
3. If the coefficient map clusters by edge direction or support side, add an
   edge-local replay that assigns source to carrier rows before nodal lumping
   instead of applying a nodal fixed policy after lumping.
4. If the coefficient map does not cluster, revisit SG flux/mobility/current
   branch formation at the active edges; row-policy changes alone should not
   be promoted.

### Execution Note 2026-06-20: Task 108 Node-Local Source Component Coefficients

Task 108 added:

```text
scripts/diagnose_pn2d_bv_source_component_coefficients.py
```

The script computes, per selected support node:

```text
required_source = flux + recombination + gauge + boundary
```

for electron and hole rows, then reports:

- `required_source / combined_source`;
- `required_source / matching_component`;
- minimum-norm coefficients `(a,b)` satisfying:

```text
required_source = a * electron_component + b * hole_component
```

It can also merge node-level edge-direction features by `node_id`.

Real outputs:

```text
focused_restart_m13p0_exact_carrier_term_states_sentaurus_qf_impact_scaled/source_component_coefficients_vela_psi_sentaurus_qf
focused_restart_m13p1_exact_carrier_term_states_sentaurus_qf_impact_scaled/source_component_coefficients_vela_psi_sentaurus_qf
focused_restart_m13p2_exact_carrier_term_states_sentaurus_qf_impact_scaled/source_component_coefficients_vela_psi_sentaurus_qf
```

For `-13.2 V`, edge features were merged from:

```text
focused_restart_m13p2_edge_direction_source_policy/edge_direction_source_policy_nodes.csv
```

Active-overlap coefficient summary:

| bias | e required / combined median | h required / combined median | corr(e required/C, y) | corr(h required/C, y) | corr(e required/C, e source frac) | corr(h required/C, e source frac) | e source frac median | h source frac median |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `-13.0 V` | `1.99927` | `1.33405` | `-0.649` | `+0.992` | `+0.695` | `-0.997` | `0.77374` | `0.22626` |
| `-13.1 V` | `2.00418` | `1.33498` | `-0.648` | `+0.992` | `+0.694` | `-0.997` | `0.77258` | `0.22742` |
| `-13.2 V` | `1.99930` | `1.34160` | `-0.642` | `+0.992` | `+0.688` | `-0.997` | `0.77141` | `0.22859` |

For `-13.2 V`, all 30 focused support nodes in the merged edge-feature file
have:

```text
dominant_active_axis = x
active_endpoint_area_fraction = 0.5
junction_normal_active_area_fraction = 0.5
junction_tangent_active_area_fraction = 0.0
```

So the existing coarse edge-direction features do not explain the node-to-node
coefficient variation. The variation is continuous along the active column and
tracks component source fractions:

- `impact_electron_fraction` has correlation `-0.996` with `y`;
- `hole_required/combined` has correlation `+0.992` with `y`;
- `hole_required/combined` has correlation `-0.997` with
  `impact_electron_fraction`.

At `-13.2 V`, fitting `required/combined` versus `hole_source_fraction` over
active-overlap nodes gives:

| row | linear fit | median abs residual | max abs residual |
| --- | --- | ---: | ---: |
| electron | `2.9754 - 3.4896 * hole_fraction` | `0.1483` | `0.4665` |
| hole | `-0.6638 + 8.7150 * hole_fraction` | `0.0148` | `0.1031` |

Interpretation:

- The electron row consistently wants about `2 * combined` at the active
  overlap median, but its node-local scatter is not well explained by the
  component fraction alone.
- The hole row is much more tightly tied to the component mix, especially the
  hole-source fraction. This is a strong hint that Sentaurus' effective hole
  row source is not a fixed nodal policy, but a current-branch/source-fraction
  weighting that changes along the active column.
- Since all coarse active-edge geometry labels are identical in the current
  focused feature file, the next discriminator must move one level deeper:
  edge-local source assignment before nodal lumping, with current direction and
  endpoint ownership retained per incident edge.

### Next Tasks After Task 108

1. Build an edge-local replay for the `vela_psi_sentaurus_qf` active-overlap
   nodes:
   - keep each incident edge's electron/hole source component;
   - preserve edge axis, endpoint side, and SG flux/current sign;
   - assign candidate row sources before summing to nodes.
2. Test whether a direction-aware edge assignment can reproduce the local
   `required/combined` map better than the nodal policy matrix:
   - especially the hole-row near-linear dependence on hole-source fraction;
   - and the electron-row residual scatter around `2 * combined`.
3. If an edge-local assignment closes the active rows, promote it only as an
   opt-in diagnostic production candidate and compare against Charon/DEVSIM
   impact source semantics before changing defaults.
4. If edge-local assignment still fails, inspect the SG flux branch entering
   the impact source itself: density interpolation, mobility field, and current
   sign/absolute-value policy are then more likely than nodal carrier-row
   policy.

### Execution Note 2026-06-20: Task 109 Edge-Local Direction Replay

Task 109 added:

```text
scripts/diagnose_pn2d_bv_edge_local_source_replay.py
```

The script reads C++ `sg_avalanche_edges.csv`, keeps each incident edge's
electron/hole source components, classifies the neighbor direction
(`left/right/up/down`) and axis (`x/y`), then aggregates source before testing
carrier-row replay policies. Because `sg_avalanche_edges.csv` stores
unscaled `s^-1` source integrals while the exact carrier-term rows are scaled
continuity residuals, the script accepts:

```text
--carrier-term-csv
```

and uses each node's `impact_combined_source` to calibrate a node-local source
scale.

Real outputs:

```text
focused_restart_m13p0_exact_carrier_term_states_sentaurus_qf_impact_scaled/edge_local_source_replay_vela_psi_sentaurus_qf
focused_restart_m13p1_exact_carrier_term_states_sentaurus_qf_impact_scaled/edge_local_source_replay_vela_psi_sentaurus_qf
focused_restart_m13p2_exact_carrier_term_states_sentaurus_qf_impact_scaled/edge_local_source_replay_vela_psi_sentaurus_qf
```

Active-overlap replay summary:

| bias | all source median | x source median | y source median | e required / all median | h required / all median | e `double_all` residual / all median | h `combined_plus_hole` residual / all median | left source median | right source median |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `-13.0 V` | `5.59728e-16` | `5.59728e-16` | `3.70e-320` | `1.99927` | `1.33405` | `-0.00073` | `+0.06651` | `2.79288e-16` | `2.80440e-16` |
| `-13.1 V` | `5.71933e-16` | `5.71933e-16` | `3.77e-320` | `2.00418` | `1.33498` | `+0.00418` | `+0.06653` | `2.85429e-16` | `2.86504e-16` |
| `-13.2 V` | `5.84348e-16` | `5.84348e-16` | `3.81e-320` | `1.99930` | `1.34160` | `-0.00070` | `+0.07216` | `2.91523e-16` | `2.92825e-16` |

Additional `-13.2 V` endpoint analysis:

- `x_combined_source` equals `all_combined_source` in active overlap;
- `y_combined_source` is numerical underflow and physically zero in this
  support;
- left/right sources are almost exactly half of all source:
  `left_frac ~= 0.499`, `right_frac ~= 0.501`;
- left/right imbalance is only `0.0013..0.0028`, far too small to explain the
  node-local required-source variation.

Interpretation:

- A coarse edge-local replay based only on incident-edge axis or endpoint side
  collapses to the same information as nodal source lumping for this support:
  all meaningful source is on the two horizontal incident edges, split almost
  symmetrically left/right.
- This rules out a simple endpoint-half ownership bug as the source of the
  active-row mismatch.
- The remaining evidence points one level deeper: the source magnitude/current
  branch being assigned to those horizontal edges differs from Sentaurus in a
  way that varies along the active column. The next replay must preserve or
  reconstruct current sign and branch-specific source density, not just edge
  ownership.

### Next Tasks After Task 109

1. Extend the edge-local replay with `edge_source_current_consistency_edges.csv`
   fields:
   - `cxx_generation_density_m3_s`;
   - `sentaurus_support_generation_m3_s`;
   - `sentaurus_edgeavg_generation_m3_s`;
   - Vela/Sentaurus electron and hole current flux proxies.
2. Build candidate source magnitudes per incident edge using Sentaurus support
   and edge-average generation/current data, then scale them into carrier-term
   units and compare required row residuals.
3. If Sentaurus edge-average/support generation closes the residual trend,
   the root cause is source-magnitude/current-branch reconstruction rather than
   carrier-row ownership.
4. If Sentaurus current-derived magnitudes still do not close, inspect SG flux
   sign and density interpolation directly in the active horizontal edges,
   especially electron-row scatter around `2 * combined`.

### Execution Note 2026-06-20: Task 110 Sentaurus Generation/Current Source Replay

Task 110 extended:

```text
scripts/diagnose_pn2d_bv_edge_local_source_replay.py
```

with an optional:

```text
--edge-current-csv edge_source_current_consistency_edges.csv
```

input.  The replay now adds, per support node:

- `current_cxx_endpoint_source`;
- `sentaurus_support_generation_source`;
- `sentaurus_edgeavg_generation_source`;
- `sentaurus_support_current_scaled_source`;
- `sentaurus_edgeavg_current_scaled_source`.

Each consistency-derived candidate is scaled with the same node-local
`impact_combined_source / sg_avalanche_edges_source_sum` factor used by the
existing carrier-term replay.  This keeps the comparison in the same
continuity-residual units as the exact carrier-term rows.

Active-overlap summary:

| bias | Vela all source median | e required / Vela all | h required / Vela all | e required / Sentaurus support generation | h required / Sentaurus support generation | e required / Sentaurus edgeavg current-scaled | h required / Sentaurus edgeavg current-scaled |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `-13.0 V` | `5.59728e-16` | `1.99927` | `1.33405` | `1.55533` | `1.03229` | `1.50147` | `1.00063` |
| `-13.1 V` | `5.71933e-16` | `2.00418` | `1.33498` | `1.54943` | `1.02672` | `1.49525` | `0.99461` |
| `-13.2 V` | `5.84348e-16` | `1.99930` | `1.34160` | `1.53588` | `1.02502` | `1.48131` | `0.99251` |

Interpretation:

- Replacing the Vela endpoint source magnitude with Sentaurus edge-average
  current-scaled generation almost closes the hole active-overlap row:
  `h required / candidate ~= 0.993..1.001`.
- The same candidate still leaves the electron row short by a stable factor:
  `e required / candidate ~= 1.48..1.50`.
- Therefore the remaining dominant mismatch is not a simple source area or
  endpoint ownership issue.  The evidence now points to carrier-branch
  continuity semantics: the hole row wants Sentaurus's larger current-derived
  generation magnitude, while the electron row still wants roughly another
  half-source contribution beyond that magnitude.

### Next Tasks After Task 110

1. Add a probe-only row-source replay variant that substitutes
   `sentaurus_edgeavg_current_scaled_source` for the hole row and leaves the
   electron row as `2 * current_cxx_endpoint_source`; compare residual L2
   against the previous best policy matrix.
2. Decompose the electron-row missing term using active-edge current branches:
   compare candidate terms proportional to Sentaurus electron current,
   Sentaurus hole current, particle current, and their absolute values.
3. If the hybrid row replay closes carrier residuals, test it as an opt-in
   BV-only diagnostic source policy over the focused `-13.0/-13.2 V` restart
   window before touching default solver behavior.
4. If the electron residual remains near `0.5 * sentaurus_edgeavg_current`
   after branch decomposition, inspect the continuity equation assembly signs
   against Charon/DEVSIM carrier residual conventions and Sentaurus log
   equation scaling.

### Execution Note 2026-06-20: Task 111 Hybrid Row-Source Replay

Task 111 extended the edge-local replay with a probe-only hybrid residual:

```text
electron row source = 2 * current_cxx_endpoint_source
hole row source     = sentaurus_edgeavg_current_scaled_source
```

This tests the direct consequence of Task 110: keep the electron-row policy
that already gives near-zero signed median in active overlap, but replace the
hole-row source magnitude with the Sentaurus edge-average current-scaled
generation that closed the hole required/source ratio.

Active-overlap comparison against the Task 107 policy-matrix best row:

| bias | policy-matrix best L2 | hybrid L2 | hybrid / best | hybrid e residual / C median | hybrid h residual / C median |
| ---: | ---: | ---: | ---: | ---: | ---: |
| `-13.0 V` | `1.01958e-15` | `1.01658e-15` | `0.99705` | `-0.00073` | `+0.00090` |
| `-13.1 V` | `1.04524e-15` | `1.03630e-15` | `0.99145` | `+0.00418` | `-0.00717` |
| `-13.2 V` | `1.06116e-15` | `1.04567e-15` | `0.98540` | `-0.00070` | `-0.01007` |

Interpretation:

- The hybrid replay closes the signed active-overlap medians for both carrier
  rows, confirming that the hole-row bulk offset is mostly a source-magnitude
  / current-branch issue.
- The total L2 improves only slightly (`0.3%..1.5%`) because a few spatially
  structured outliers dominate the norm.
- At `-13.2 V`, the largest hybrid residual nodes are not random:
  high-y active-column nodes (`748`, `740`, `716`, `708`) have large positive
  hole residual, while low-y nodes (`351`, `360`, `384`) have positive electron
  residual and negative/small hole residual.  This points to a current-branch
  or direction-dependent source split, not a global source factor.

### Next Tasks After Task 111

1. Extend the edge-local replay to decompose Sentaurus edge-average current
   source into electron-current and hole-current branch candidates, including
   signed and absolute-value variants.
2. For active-overlap outliers, test whether the remaining hybrid residual
   aligns with the electron-current branch, hole-current branch, particle
   branch, or a sign flip across the active column.
3. If branch decomposition explains the outliers, implement an opt-in
   diagnostic source policy and run the focused BV restart window to see
   whether IV/current agreement improves.
4. If branch decomposition does not explain the outliers, move from source
   policy to carrier residual assembly convention: compare Vela carrier-row
   signs/scaling against Charon, DEVSIM, and Sentaurus equation logs.

### Execution Note 2026-06-20: Task 112 Edge-Average Current Branch Decomposition

Task 112 extended the same edge-local replay with Sentaurus edge-average
current branch candidates:

```text
sentaurus_edgeavg_electron_current_scaled_source
sentaurus_edgeavg_hole_current_scaled_source
sentaurus_edgeavg_abs_electron_current_scaled_source
sentaurus_edgeavg_abs_hole_current_scaled_source
```

The branch candidates use the same scale convention as the particle-current
candidate:

```text
cxx_endpoint_source * sentaurus_edgeavg_branch_current_flux / cxx_particle_flux
```

with absolute-value variants for sign-flip checks.

At `-13.2 V`, active-overlap nodes sorted by y show:

| node | y | hybrid e residual / C | hybrid h residual / C | e branch / particle | h branch / particle |
| ---: | ---: | ---: | ---: | ---: | ---: |
| `352` | `0` | `-0.21773` | `-0.16278` | `0.45635` | `0.54365` |
| `351` | `1.5625e-08` | `+0.49218` | `-0.16899` | `0.45630` | `0.54370` |
| `354` | `3.125e-08` | `+0.05413` | `-0.27252` | `0.45637` | `0.54363` |
| `360` | `4.6875e-08` | `+0.38490` | `-0.10775` | `0.45631` | `0.54369` |
| `384` | `7.8125e-08` | `+0.34985` | `-0.04894` | `0.45631` | `0.54369` |
| `392` | `1.09375e-07` | `+0.30229` | `+0.02880` | `0.45631` | `0.54369` |
| `708` | `2.65625e-07` | `-0.05554` | `+0.58668` | `0.45631` | `0.54369` |
| `716` | `2.96875e-07` | `-0.12607` | `+0.71085` | `0.45631` | `0.54369` |
| `740` | `3.28125e-07` | `-0.19776` | `+0.85976` | `0.45631` | `0.54369` |
| `748` | `3.59375e-07` | `-0.25993` | `+0.97756` | `0.45631` | `0.54369` |

Interpretation:

- The Sentaurus edge-average electron/hole branch fraction is essentially
  constant across the active column: about `45.63% / 54.37%`.
- The remaining hybrid residual is not constant: the hole residual grows
  strongly positive at high y, while the electron residual is largest positive
  near the lower active column and turns negative at high y.
- Therefore the remaining outliers are not explained by branch split or a
  simple current sign flip.  Continuing to tune electron-vs-hole branch
  weights would be a dead end.

### Next Tasks After Task 112

1. Reconstruct the active-overlap carrier residual using Sentaurus state on
   both SG flux and source sides, but vary only one interpolation/scaling choice
   at a time:
   - support-node generation vs edge-average generation;
   - endpoint vs edge-average density/current state;
   - old-Slotboom vs density SG flux form;
   - node-volume/continuity scaling.
2. Build an outlier-focused table for the ten `-13.2 V` overlap nodes that
   includes potential, quasi-Fermi potentials, carrier densities, mobility,
   SG flux, generation source, and residual pieces in y order.
3. If the residual trend tracks state interpolation, implement an opt-in
   diagnostic residual replay using the matching Sentaurus interpolation.
4. If the residual trend persists even with Sentaurus state/interpolation,
   inspect carrier-row assembly signs/scales against Charon/DEVSIM and the
   Sentaurus nonlinear equation scaling/log output.

### Execution Note 2026-06-20: Task 113 Active-Overlap Outlier Table

Task 113 added:

```text
scripts/diagnose_pn2d_bv_active_overlap_outlier_table.py
```

The script merges, by `node_id`, the focused active-overlap diagnostics from:

- exact carrier-term state rows: potential, quasi-Fermi potentials, carrier
  densities, `ni_eff`, node volume, and non-impact carrier terms;
- edge-local source replay rows: hybrid residuals and source candidates;
- edge source/current consistency rows: C++ and Sentaurus current/generation
  factors averaged over active incident edges;
- active-edge mobility inputs: mobility, electric-field, and quasi-Fermi drive
  ratios;
- active-edge mixed-state replay: density, particle-flux, and generation ratios;
- active-support residual proxy: selected Sentaurus-state and shifted-QF
  residual proxy variants.

Real output for `-13.2 V`:

```text
focused_restart_m13p2_exact_carrier_term_states_sentaurus_qf_impact_scaled/active_overlap_outlier_table_vela_psi_sentaurus_qf
```

Summary:

```text
row_count = 10
max_hybrid_residual_node_id = 748
max_hybrid_residual_l2 = 5.47594e-16
```

Key active-overlap trend table:

| node | y | e residual / C | h residual / C | mixed particle flux / Sentaurus | mixed generation / Sentaurus | e density / Sentaurus | h density / Sentaurus | e mobility / Sentaurus | h mobility / Sentaurus |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `352` | `0` | `-0.21773` | `-0.16278` | `1.06638` | `1.12084` | `1.25538` | `0.79650` | `0.94156` | `1.21143` |
| `351` | `1.5625e-08` | `+0.49218` | `-0.16899` | `1.06619` | `1.12011` | `1.25280` | `0.79813` | `0.94245` | `1.21030` |
| `354` | `3.125e-08` | `+0.05413` | `-0.27252` | `1.06609` | `1.11684` | `1.24541` | `0.80284` | `0.94249` | `1.21019` |
| `360` | `4.6875e-08` | `+0.38490` | `-0.10775` | `1.06607` | `1.11189` | `1.23394` | `0.81018` | `0.94254` | `1.21026` |
| `384` | `7.8125e-08` | `+0.34985` | `-0.04894` | `1.06656` | `1.09826` | `1.20213` | `0.83178` | `0.94255` | `1.21013` |
| `392` | `1.09375e-07` | `+0.30229` | `+0.02880` | `1.06808` | `1.08267` | `1.16425` | `0.85877` | `0.94255` | `1.21011` |
| `708` | `2.65625e-07` | `-0.05554` | `+0.58668` | `1.09427` | `1.01804` | `0.98349` | `1.01667` | `0.94238` | `1.21025` |
| `716` | `2.96875e-07` | `-0.12607` | `+0.71085` | `1.10307` | `1.00875` | `0.95124` | `1.05111` | `0.94237` | `1.21020` |
| `740` | `3.28125e-07` | `-0.19776` | `+0.85976` | `1.11328` | `1.00042` | `0.91947` | `1.08732` | `0.94235` | `1.21013` |
| `748` | `3.59375e-07` | `-0.25993` | `+0.97756` | `1.12493` | `0.99321` | `0.88842` | `1.12534` | `0.94230` | `1.21009` |

Correlations across these ten active-overlap nodes:

| target | factor | correlation |
| --- | --- | ---: |
| hole residual / C | y | `+0.99154` |
| hole residual / C | mixed particle flux / Sentaurus | `+0.98145` |
| hole residual / C | mixed generation / Sentaurus | `-0.99064` |
| hole residual / C | mixed electron density / Sentaurus | `-0.99503` |
| hole residual / C | mixed hole density / Sentaurus | `+0.99677` |
| hole residual / C | electron mobility / Sentaurus | `+0.02100` |
| electron residual / C | mixed particle flux / Sentaurus | `-0.72144` |
| electron residual / C | mixed electron density / Sentaurus | `+0.67249` |
| electron residual / C | mixed hole density / Sentaurus | `-0.68860` |

Interpretation:

- Mobility and high-field drive ratios are nearly flat along the active column;
  they do not explain the remaining structured hybrid residual.
- The hole residual is almost perfectly aligned with carrier-density and mixed
  SG flux state variation: as Vela/Sentaurus hole density and particle flux
  ratios increase with y, the hole residual grows positive.
- The source-generation ratio moves in the opposite direction and the current
  branch split was already shown to be constant, so the next likely source is
  the carrier-density / SG-flux interpolation state used in continuity.

### Next Tasks After Task 113

1. Build a focused state-interpolation residual replay for the same ten
   `-13.2 V` overlap nodes:
   - Vela density + Vela mobility + Vela source;
   - Sentaurus density + Vela mobility + Vela source;
   - Sentaurus density + Sentaurus mobility + Vela source;
   - Sentaurus density + Sentaurus mobility + Sentaurus edgeavg-current source.
2. Keep the source policy fixed while swapping density/flux state first.  The
   expected discriminator is whether the y-trending hole residual collapses
   when Sentaurus density or mixed SG state is used.
3. If Sentaurus density/flux state collapses the y trend, trace the specific SG
   flux formula and density averaging difference against Charon/DEVSIM before
   changing Vela defaults.
4. If the y trend remains even with Sentaurus density/flux state, move to
   continuity scaling and residual assembly conventions.

### Execution Note 2026-06-20: Task 114 State-Interpolation Residual Replay

Task 114 extended:

```text
scripts/diagnose_pn2d_bv_active_support_residual_proxy.py
```

with a focused state/source probe matrix:

- `probe_vela_density_vela_mobility_vela_source`;
- `probe_sentaurus_density_vela_mobility_vela_source`;
- `probe_sentaurus_density_sentaurus_mobility_vela_source`;
- `probe_sentaurus_density_sentaurus_mobility_edgeavg_current_source`.

Implementation details:

- `--edge-local-source-csv` accepts
  `edge_local_source_replay_nodes.csv`;
- the edgeavg-current source is used as a ratio,
  `sentaurus_edgeavg_current_scaled_source / current_cxx_endpoint_source`,
  multiplied back onto the residual proxy's local Vela edge source.  This keeps
  the edge-local exact carrier-row quantity and the residual-proxy `s^-1`
  source on the same scale;
- the script now uses Windows long-path wrappers for CSV/summary I/O, matching
  the deeper focused BV report paths.

Real outputs:

```text
focused_restart_m13p2_exact_carrier_term_states_sentaurus_qf_impact_scaled/edge_local_source_replay_vela_psi_sentaurus_qf
focused_restart_m13p2_state_interpolation_residual_proxy
```

Overlap active-support medians for the four probe variants:

| transport mode | variant | impact/Sentaurus source | electron transport/Sentaurus source | hole transport/Sentaurus source | electron residual/impact | hole residual/impact | hole residual vs y corr |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `density_sg` | `probe_vela_density_vela_mobility_vela_source` | `0.38205` | `1.06812` | `1.03861` | `+0.54718` | `+0.48343` | `+0.90917` |
| `density_sg` | `probe_sentaurus_density_vela_mobility_vela_source` | `0.38205` | `1.72919` | `1.78093` | `+2.29001` | `+2.43530` | `+0.12660` |
| `density_sg` | `probe_sentaurus_density_sentaurus_mobility_vela_source` | `0.38205` | `1.61565` | `1.64346` | `+2.00258` | `+2.06594` | `+0.66943` |
| `density_sg` | `probe_sentaurus_density_sentaurus_mobility_edgeavg_current_source` | `0.51680` | `1.61565` | `1.64346` | `+1.21410` | `+1.26650` | `+0.53893` |
| `qf_old_slotboom_ni` | `probe_vela_density_vela_mobility_vela_source` | `0.38205` | `0.69502` | `0.65151` | `-0.42627` | `-0.53073` | `+0.90089` |
| `qf_old_slotboom_ni` | `probe_sentaurus_density_vela_mobility_vela_source` | `0.38205` | `1.25969` | `1.22097` | `+1.06049` | `+0.96491` | `-0.13304` |
| `qf_old_slotboom_ni` | `probe_sentaurus_density_sentaurus_mobility_vela_source` | `0.38205` | `1.14326` | `1.08020` | `+0.76221` | `+0.59161` | `+0.57055` |
| `qf_old_slotboom_ni` | `probe_sentaurus_density_sentaurus_mobility_edgeavg_current_source` | `0.51680` | `1.14326` | `1.08020` | `+0.29946` | `+0.17660` | `+0.47194` |

Interpretation:

- Swapping from Vela to Sentaurus density/SG state does **not** close the
  active-overlap residual; it increases transport relative to Sentaurus source.
- Vela mobility with Sentaurus density makes the residual even larger, so the
  mobility ratio is not the missing damping term.
- Sentaurus mobility plus qF-old-Slotboom transport is closer than density-SG,
  but still leaves `~+0.30` electron and `~+0.18` hole residual/impact even
  with the edgeavg-current source.
- The y trend weakens when Sentaurus density/mobility/source are all used, but
  it does not collapse.  This points away from a pure density interpolation
  bug and toward continuity scaling/residual assembly conventions, including
  control-volume normalization, edge-to-node source ownership, or exact row
  sign/charge scaling.

### Next Tasks After Task 114

1. Build a focused continuity-scaling decomposition for the same ten overlap
   nodes:
   - compare residual proxy source units against exact C++ carrier-row units;
   - report per-node scale factors from `impact_source_s_inv` to exact
     `electron_required_source` / `hole_required_source`;
   - test whether a single control-volume, charge, or row-normalization factor
     explains the remaining `+0.30/+0.18` qF-old-Slotboom edgeavg-current
     residual.
2. If no scalar factor explains the gap, trace Vela continuity assembly row
   sign and source placement against Charon/DEVSIM-style residual assembly:
   transport divergence, recombination sign, impact source sign, and terminal
   current extraction convention.
3. Keep Sentaurus density/mobility/source as a diagnostic replay only.  The
   current evidence does not justify changing Vela defaults to Sentaurus-like
   density or mobility interpolation.

### Execution Note 2026-06-20: Task 115 Continuity-Scaling Decomposition

Task 115 added:

```text
scripts/diagnose_pn2d_bv_continuity_scaling_decomposition.py
```

The script joins the exact C++ carrier-row diagnostics, edge-local source
replay, and active-support residual proxy for the same ten `-13.2 V`
active-overlap nodes.  It reports scale factors between residual-proxy
`s^-1` quantities and exact carrier-row source units, plus carrier-specific
required-source scale ratios.

Real output:

```text
focused_restart_m13p2_continuity_scaling_decomposition/
  continuity_scaling_decomposition_nodes.csv
  continuity_scaling_decomposition_summary.json
```

Summary for variant `vela_psi_sentaurus_qf`,
transport model `qf_old_slotboom_ni`, and proxy variant
`probe_sentaurus_density_sentaurus_mobility_edgeavg_current_source`:

| metric | median | min | max |
| --- | ---: | ---: | ---: |
| `exact_edgeavg_current_over_proxy_impact_scale` | `4.10920e-21` | `3.64695e-21` | `4.17663e-21` |
| `exact_all_combined_over_proxy_impact_scale` | `3.04012e-21` | `2.68578e-21` | `3.10719e-21` |
| `exact_electron_required_over_proxy_required_scale` | `4.55228e-21` | `3.59655e-21` | `6.01787e-21` |
| `exact_hole_required_over_proxy_required_scale` | `3.47607e-21` | `2.83893e-21` | `5.32770e-21` |
| `electron_required_scale_over_edgeavg_scale` | `1.14396` | `0.98016` | `1.44084` |
| `hole_required_scale_over_edgeavg_scale` | `0.84285` | `0.68194` | `1.46086` |
| `edgeavg_scale_over_inverse_srh_scale` | `1.43412` | `1.27279` | `1.45765` |
| `proxy_electron_residual_over_impact` | `0.29946` | `0.28678` | `0.35448` |
| `proxy_hole_residual_over_impact` | `0.17660` | `0.07324` | `0.17915` |

Interpretation:

- A single scalar continuity-unit correction is rejected.  The proxy-to-exact
  edgeavg-current source scale is near `4.11e-21`, but the electron and hole
  required-source scales diverge by carrier and by y-position.
- The edgeavg-current scale is consistently larger than the inverse SRH-derived
  row scale by `~1.27x` to `~1.46x`, so a pure control-volume or row-unit
  normalization does not close the residual.
- `electron_required_scale_over_edgeavg_scale` and
  `hole_required_scale_over_edgeavg_scale` trend in opposite directions across
  the active column.  This points to row/source ownership or carrier-specific
  residual assembly conventions rather than a mobility, high-field drive, or
  one-factor scaling error.
- The relevant Vela C++ path currently subtracts the same SG avalanche
  `combined` source from both electron and hole continuity rows, and the
  carrier-term diagnostics divide flux, recombination, impact, and source
  components by `C0 * D0`.  Charon's CVFEM pattern separates flux scattering
  from scalar source integration, so the next discriminator should focus on
  source ownership and row-local source assembly.

### Next Tasks After Task 115

1. Instrument or export raw `sgAvalancheSourceIntegrals` and
   `sgAvalancheSourceComponents.{electron,hole,combined}` before
   `continuityScale` for the exact active rows.  Compare them directly with
   the Python edge-local and residual-proxy source units.
2. Add a carrier-specific source-ownership replay:
   - electron row uses electron branch only;
   - hole row uses hole branch only;
   - both rows use combined source;
   - both rows use Sentaurus edgeavg-current source.
   The discriminator is whether one ownership convention collapses the
   opposite y-trends in electron and hole required scale.
3. Check whether Sentaurus reports total avalanche generation while assembling
   branch-specific carrier-row sources.  Use the Charon density-damped
   avalanche and scalar-source integration patterns as references, but keep
   this as a diagnostic until Vela row-level evidence closes.
4. Only after row-source convention is identified, adjust Vela impact source
   assembly or add a Sentaurus-compatible diagnostic mode, then rerun the
   `-13.2 V` focused residual and the BV IV sweep.

### Execution Note 2026-06-20: Task 116 Source-Ownership Replay

Task 116 added:

```text
scripts/diagnose_pn2d_bv_source_ownership_replay.py
```

The script ranks carrier-row source policies directly from
`edge_local_source_replay_nodes.csv`, including both Vela component policies
and Sentaurus edge-averaged current-scaled total / electron / hole source
policies.

Real output:

```text
focused_restart_m13p2_source_ownership_replay/
  source_ownership_replay_summary.csv
  source_ownership_replay_summary.json
```

Top policies on the same ten `-13.2 V` active-overlap nodes:

| electron policy | hole policy | combined L2 | e median residual / combined | h median residual / combined | e corr(y) | h corr(y) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `vela_double_combined` | `vela_combined_plus_electron` | `1.02611e-15` | `-0.00070` | `-0.38895` | `-0.64188` | `+0.99242` |
| `vela_double_combined` | `sentaurus_edgeavg_current` | `1.04567e-15` | `-0.00070` | `-0.01007` | `-0.64188` | `+0.99154` |
| `vela_double_combined` | `vela_combined_plus_hole` | `1.11456e-15` | `-0.00070` | `+0.07216` | `-0.64188` | `+0.99149` |
| `sentaurus_edgeavg_current` | `sentaurus_edgeavg_current` | `1.73019e-15` | `+0.64939` | `-0.01007` | `-0.65590` | `+0.99154` |
| `sentaurus_edgeavg_electron_current` | `sentaurus_edgeavg_hole_current` | `3.24836e-15` | `+1.38328` | `+0.60672` | `-0.64828` | `+0.99174` |
| `vela_electron` | `vela_hole` | `3.56257e-15` | `+1.26577` | `+1.07216` | `-0.63327` | `+0.99149` |

Per-node discriminator:

| quantity | low-y range | high-y range | conclusion |
| --- | ---: | ---: | --- |
| `electron_required_source / all_combined_source` | `1.78` to `2.49` | `1.74` to `1.94` | electron row is closest to roughly `2 * combined`, but not exactly constant |
| `hole_required_source / all_combined_source` | `1.18` to `1.38` | `1.94` to `2.34` | hole row carries the dominant y trend |
| `hole_required_source / sentaurus_edgeavg_current_scaled_source` | `0.88` to `1.02` | `1.43` to `1.72` | Sentaurus total current-scaled source does not remove the high-y excess |
| `sentaurus_edgeavg_electron_current_scaled_source / all_combined_source` | `0.613` to `0.618` | `0.619` to `0.620` | Sentaurus branch split is nearly flat |
| `sentaurus_edgeavg_hole_current_scaled_source / all_combined_source` | `0.730` to `0.736` | `0.737` to `0.738` | branch split cannot explain the y-trending residual |

Interpretation:

- Carrier-specific Sentaurus branch-current sources are rejected as the missing
  carrier-row ownership rule.  They make both carrier rows worse than the
  total/combined-source candidates.
- Source ownership can reduce the electron median residual, but the hole row
  keeps a `~+0.99` y correlation under every tested policy.  The remaining
  mismatch is therefore not primarily the total-vs-branch avalanche source
  split.
- The required source is `flux + recombination + boundary/gauge`, so the next
  root-cause probe should decompose the non-impact transport row term rather
  than changing avalanche source ownership immediately.

### Next Tasks After Task 116

1. Build a focused non-impact row-term decomposition for the same ten overlap
   nodes:
   - electron/hole flux contribution by incident edge;
   - recombination contribution;
   - boundary/gauge contribution;
   - signed edge orientation and contact-adjacent edge class.
2. Compare `electron_required_source` and `hole_required_source` against the
   per-edge SG flux replay already used in the active-edge diagnostics.  The
   discriminator is whether the hole y trend follows a small set of incident
   edges, contact-edge treatment, or all local edges uniformly.
3. If the trend localizes to contact-adjacent or vertical active edges, inspect
   SG flux orientation and boundary projection before touching impact source
   assembly.
4. If the trend is uniform across all incident transport fluxes, compare Vela's
   qF-old-Slotboom SG expression against Sentaurus/Charon conventions for
   variable `ni_eff`, density placement, and row normalization.

### Execution Note 2026-06-20: Task 117 Non-Impact Row-Term and Edge-Flux Y-Trend

Task 117 added:

```text
scripts/diagnose_pn2d_bv_row_term_decomposition.py
scripts/diagnose_pn2d_bv_edge_flux_y_trend.py
```

Real outputs:

```text
focused_restart_m13p2_row_term_decomposition/
  row_term_decomposition_nodes.csv
  row_term_decomposition_summary.json

focused_restart_m13p2_edge_flux_y_trend/
  edge_flux_y_trend_nodes.csv
  edge_flux_y_trend_summary.json
```

The row-term decomposition used source policy
`sentaurus_edgeavg_current_scaled_source` on the same ten `-13.2 V`
active-overlap nodes.  Main correlations:

| term | corr(y) | median | min | max |
| --- | ---: | ---: | ---: | ---: |
| `hole_flux` | `+0.86967` | `1.34973e-15` | `6.36798e-16` | `1.79226e-15` |
| `hole_recombination` | `-0.37781` | `-5.27960e-16` | `-5.27960e-16` | `-2.63980e-16` |
| `hole_impact_source` | `-0.01121` | `5.84348e-16` | `3.16030e-16` | `6.32639e-16` |
| `hole_required_source` | `+0.94400` | `8.21770e-16` | `3.72818e-16` | `1.26430e-15` |
| `hole_required_over_source` | `+0.99158` | `0.99251` | `0.79711` | `1.71992` |
| `hole_residual_under_policy_over_source` | `+0.99158` | `-0.00749` | `-0.20289` | `+0.71992` |
| `electron_flux` | `-0.22751` | `1.71991e-15` | `8.27231e-16` | `2.10461e-15` |
| `electron_required_source` | `-0.37280` | `1.19195e-15` | `5.63251e-16` | `1.57665e-15` |
| `electron_residual_under_policy_over_source` | `-0.66104` | `+0.48131` | `+0.28147` | `+0.85405` |

Node-level source-normalized trend:

| node | y (um) | e flux/source | e required/source | h flux/source | h required/source | h residual/source |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `352` | `0` | `1.94982` | `1.32761` | `1.50096` | `0.87875` | `-0.12125` |
| `351` | `0.015625` | `2.47490` | `1.85405` | `1.49513` | `0.87428` | `-0.12572` |
| `354` | `0.03125` | `2.15401` | `1.52930` | `1.42182` | `0.79711` | `-0.20289` |
| `360` | `0.046875` | `2.39764` | `1.77211` | `1.54547` | `0.91993` | `-0.08007` |
| `384` | `0.078125` | `2.37405` | `1.74107` | `1.59672` | `0.96374` | `-0.03626` |
| `392` | `0.109375` | `2.34288` | `1.70076` | `1.66341` | `1.02128` | `+0.02128` |
| `708` | `0.265625` | `2.12664` | `1.43331` | `2.12579` | `1.43246` | `+0.43246` |
| `716` | `0.296875` | `2.08400` | `1.38158` | `2.22650` | `1.52409` | `+0.52409` |
| `740` | `0.328125` | `2.03962` | `1.32864` | `2.34481` | `1.63383` | `+0.63383` |
| `748` | `0.359375` | `1.99969` | `1.28147` | `2.43814` | `1.71992` | `+0.71992` |

The edge-flux localization used
`focused_restart_m13p2_edge_state_substitution_current/edge_state_substitution_current_edges.csv`
for variant `vela_psi_sentaurus_qf`.  It found:

```text
row_count = 10 support nodes
edge_count = 20 active edges
axis_counts = {"x": 20}
```

| node | y (um) | active axes | edges | e qF/Sentaurus edgeavg | h qF/Sentaurus edgeavg | particle qF/Sentaurus edgeavg |
| ---: | ---: | --- | ---: | ---: | ---: | ---: |
| `352` | `0` | `x` | `2` | `1.25864` | `0.79896` | `1.00848` |
| `351` | `0.015625` | `x` | `2` | `1.25746` | `0.80157` | `1.00935` |
| `354` | `0.03125` | `x` | `2` | `1.24904` | `0.80526` | `1.00754` |
| `360` | `0.046875` | `x` | `2` | `1.23876` | `0.81375` | `1.00744` |
| `384` | `0.078125` | `x` | `2` | `1.20688` | `0.83529` | `1.00461` |
| `392` | `0.109375` | `x` | `2` | `1.16894` | `0.86239` | `1.00203` |
| `708` | `0.265625` | `x` | `2` | `0.98721` | `1.02125` | `1.00549` |
| `716` | `0.296875` | `x` | `2` | `0.95490` | `1.05577` | `1.00951` |
| `740` | `0.328125` | `x` | `2` | `0.92294` | `1.09214` | `1.01471` |
| `748` | `0.359375` | `x` | `2` | `0.89184` | `1.13022` | `1.02123` |

Edge-level correlations:

| metric | corr(y) | median | min | max |
| --- | ---: | ---: | ---: | ---: |
| `electron_qf_over_sentaurus_edgeavg_median` | `-0.99810` | `1.18791` | `0.89184` | `1.25864` |
| `hole_qf_over_sentaurus_edgeavg_median` | `+0.99560` | `0.84884` | `0.79896` | `1.13022` |
| `particle_qf_over_sentaurus_edgeavg_median` | `+0.60237` | `1.00801` | `1.00203` | `1.02123` |

Interpretation:

- The remaining hole residual y trend is a transport-row problem.  It is
  dominated by `hole_flux`; SRH, impact source, boundary, and gauge terms do
  not explain it.
- The active-overlap edges are all x-directed and there are exactly two active
  edges per support node, so this is not a mixed-axis or isolated-edge artifact.
- Total particle flux already matches Sentaurus within roughly `0.2%` to
  `2.1%`, but the carrier split is wrong in opposite directions: Vela's
  qF-old-Slotboom electron flux is high at low y and low at high y, while hole
  flux is low at low y and high at high y.  This explains why total avalanche
  source can look close while the carrier continuity rows remain structured.
- The next likely source is the SG carrier split itself: quasi-Fermi / `ni_eff`
  endpoint placement, old-Slotboom density reconstruction, or Sentaurus current
  component convention.  The evidence no longer supports changing source
  ownership first.

### Next Tasks After Task 117

1. Build an edge-level SG carrier-split decomposition for the same 20 x-active
   overlap edges:
   - endpoint `psi`, `phin`, `phip`, `ni_eff`;
   - electron and hole quasi-Fermi drops;
   - Bernoulli arguments and endpoint density factors;
   - Vela qF-old-Slotboom electron/hole fluxes;
   - Sentaurus edgeavg electron/hole current-derived fluxes.
2. Test whether using Sentaurus `ni_eff` only, Sentaurus qF only, or Sentaurus
   density-only endpoint factors removes the opposite electron/hole y trends.
3. If total particle flux remains close but carrier split stays inverted, check
   Sentaurus current component sign/convention and whether `eCurrentDensity`
   / `hCurrentDensity` should be interpreted as carrier-current components or
   continuity particle-flux components for avalanche/source coupling.
4. Only after the edge carrier split is closed should Vela defaults be changed;
   until then keep all source ownership changes as diagnostic-only.

### Execution Note 2026-06-20: Task 118 SG Carrier-Split Decomposition

Task 118 added:

```text
scripts/diagnose_pn2d_bv_sg_carrier_split_decomposition.py
```

It first generated pure Vela/Sentaurus endpoint SG flux-form rows for the 20
x-active overlap edges:

```text
focused_restart_m13p2_sg_flux_form_active_edges/
  sg_flux_form_edges.csv
  sg_flux_form_summary.json
```

Then it merged those rows with the existing mixed-state current replay:

```text
focused_restart_m13p2_sg_carrier_split_decomposition/
  sg_carrier_split_decomposition_edges.csv
  sg_carrier_split_decomposition_summary.json
```

Real summary for variant `vela_psi_sentaurus_qf` at `-13.2 V`:

```text
row_count = 20
support_node_count = 10
edge_count = 20
edge_axes = ["x"]
```

Key correlations:

| metric | corr(y) | median | min | max |
| --- | ---: | ---: | ---: | ---: |
| `electron_fraction_delta` | `-0.99594` | `+0.08370` | `-0.06164` | `+0.11643` |
| `hole_fraction_delta` | `+0.99594` | `-0.08370` | `-0.11643` | `+0.06164` |
| `mixed_electron_over_sentaurus` | `-0.99670` | `1.18769` | `0.88607` | `1.26911` |
| `mixed_hole_over_sentaurus` | `+0.99446` | `0.84867` | `0.79491` | `1.13800` |
| `mixed_particle_over_sentaurus` | `+0.44764` | `1.00951` | `0.99710` | `1.02366` |
| `mixed_electron_density_factor_0` | `-0.97787` | `1.15194e12` | `8.44119e11` | `1.25401e12` |
| `mixed_electron_density_factor_1` | `-0.96767` | `1.08182e12` | `7.90144e11` | `1.18781e12` |
| `mixed_hole_density_factor_0` | `+0.95827` | `2.94620e7` | `2.63892e7` | `4.03073e7` |
| `mixed_hole_density_factor_1` | `+0.96981` | `3.14215e7` | `2.85640e7` | `4.30130e7` |
| `mixed_electron_density_factor_ratio_1_over_0` | `+0.00140` | `0.94213` | `0.93606` | `0.94903` |
| `mixed_hole_density_factor_ratio_1_over_0` | `+0.02600` | `1.07391` | `1.06293` | `1.08241` |
| `vela_electron_mobility_over_sentaurus` | `+0.02765` | `0.95426` | `0.94828` | `0.95840` |
| `vela_hole_mobility_over_sentaurus` | `-0.02404` | `1.23369` | `1.22898` | `1.23632` |
| `sentaurus_inferred_ni_electron_ratio_1_over_0` | `~0` | `1.01447` | `0.84373` | `1.18522` |
| `sentaurus_inferred_ni_hole_ratio_1_over_0` | `~0` | `1.01447` | `0.84373` | `1.18522` |

Sentaurus self-consistency check:

| metric | median | min | max |
| --- | ---: | ---: | ---: |
| `sentaurus_qf_electron_over_sentaurus` | `1.00744` | `0.99920` | `1.01504` |
| `sentaurus_qf_hole_over_sentaurus` | `1.00756` | `0.99726` | `1.01590` |
| `vela_qf_electron_over_sentaurus` | `927.89` | `821.24` | `1045.14` |
| `vela_qf_hole_over_sentaurus` | `140.71` | `118.97` | `166.04` |

Support-node medians:

| node | y (um) | e fraction delta | h fraction delta | mixed e/Sentaurus e | mixed h/Sentaurus h |
| ---: | ---: | ---: | ---: | ---: | ---: |
| `352` | `0` | `+0.11295` | `-0.11295` | `1.25864` | `0.79896` |
| `351` | `0.015625` | `+0.11193` | `-0.11193` | `1.25746` | `0.80157` |
| `354` | `0.03125` | `+0.10915` | `-0.10915` | `1.24904` | `0.80526` |
| `360` | `0.046875` | `+0.10454` | `-0.10454` | `1.23876` | `0.81375` |
| `384` | `0.078125` | `+0.09166` | `-0.09166` | `1.20688` | `0.83529` |
| `392` | `0.109375` | `+0.07581` | `-0.07581` | `1.16894` | `0.86239` |
| `708` | `0.265625` | `-0.00841` | `+0.00841` | `0.98721` | `1.02125` |
| `716` | `0.296875` | `-0.02478` | `+0.02478` | `0.95490` | `1.05577` |
| `740` | `0.328125` | `-0.04135` | `+0.04135` | `0.92294` | `1.09214` |
| `748` | `0.359375` | `-0.05787` | `+0.05787` | `0.89184` | `1.13022` |

Additional drop checks:

| quantity | corr(y) | median | min | max |
| --- | ---: | ---: | ---: | ---: |
| `mixed_psi_drop_V` | `-0.16853` | `0.35515` | `0.35514` | `0.35516` |
| `sentaurus_psi_drop_V` | `+0.50888` | `0.35853` | `0.35853` | `0.35853` |
| `mixed_phin_drop_V` | `-0.00032` | `0.35667` | `0.35241` | `0.36097` |
| `mixed_phip_drop_V` | `+0.00084` | `0.35696` | `0.35276` | `0.36121` |
| `mixed_electron_eta` | `-0.00027` | `13.73782` | `13.56750` | `13.90774` |
| `mixed_hole_eta` | `-0.00027` | `13.73762` | `13.56750` | `13.90813` |

Interpretation:

- Sentaurus current component convention is not the main error.  When pure
  Sentaurus endpoint state is passed through Vela's SG qF formula, electron and
  hole fluxes match Sentaurus current-derived fluxes within about `~1%`.
- The y-trending split appears only in the mixed state
  `Vela psi + Sentaurus qF`.  Total particle flux remains close, but electron
  and hole fractions shift in opposite directions.
- Endpoint `ni_eff` ratios, Bernoulli eta, and mobility ratios are essentially
  flat with y and cannot explain the carrier-split drift by themselves.
- The strongest correlated quantities are the absolute mixed density factors:
  electron factors decrease with y while hole factors increase with y.  This
  points to an electrostatic-potential / quasi-Fermi absolute alignment problem
  along the active column, not an avalanche source ownership problem and not a
  wrong SG formula.

### Next Tasks After Task 118

1. Build a focused potential/qF alignment replay for the same 20 active edges:
   - `Sentaurus psi + Sentaurus qF + Vela mobility`;
   - `Vela psi + Sentaurus qF + Sentaurus mobility`;
   - `Sentaurus psi + Sentaurus qF + Sentaurus mobility`;
   - `Vela psi shifted by best-fit constant + Sentaurus qF + Vela mobility`.
2. Determine whether a constant potential offset, a y-dependent potential
   offset, or mobility scaling is enough to collapse the carrier fraction
   deltas.  The first pass should fit only `psi` shifts, not qF shifts.
3. If a y-dependent `psi` mismatch is required, trace back to Poisson/fixed
   charge/contact potential alignment at the same active column before changing
   SG flux or impact source code.
4. Only if `Sentaurus psi + Sentaurus qF + Vela mobility` still fails should
   mobility interpolation or SG formula variants be promoted from diagnostic
   mode to solver changes.

### Execution Note 2026-06-20: Task 119 Potential/qF Alignment Replay

Task 119 added:

```text
scripts/diagnose_pn2d_bv_potential_qf_alignment_replay.py
```

The script replays the same 20 x-active overlap edges with controlled endpoint
substitution:

```text
focused_restart_m13p2_potential_qf_alignment_replay/
  potential_qf_alignment_replay_edges.csv
  potential_qf_alignment_replay_summary.json
```

Real summary for `-13.2 V`:

```text
row_count = 100
edge_count = 20
support_node_count = 10
```

| replay variant | psi shift (V) | fraction L2 | e/Sentaurus median | h/Sentaurus median | particle/Sentaurus median | e fraction delta median | h fraction delta median |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `sentaurus_psi_sentaurus_qf_sentaurus_mobility` | `0` | `0.00527` | `1.00744` | `1.00756` | `1.00711` | `+8.70e-5` | `-8.70e-5` |
| `sentaurus_psi_sentaurus_qf_vela_mobility` | `0` | `0.08887` | `0.96105` | `1.24107` | `1.11345` | `-0.06258` | `+0.06258` |
| `vela_psi_sentaurus_qf_sentaurus_mobility` | `0` | `0.76926` | `242.21` | `0.00469` | `114.61` | `+0.54369` | `-0.54369` |
| `vela_psi_sentaurus_qf_vela_mobility` | `0` | `0.76925` | `230.51` | `0.00580` | `109.08` | `+0.54368` | `-0.54368` |
| `vela_psi_shifted_const_sentaurus_qf_vela_mobility` | `-0.1365` | `0.10743` | `1.17380` | `1.13871` | `1.15531` | `+0.00757` | `-0.00757` |

Key y correlations:

| replay variant | e/Sentaurus corr(y) | h/Sentaurus corr(y) | e fraction delta corr(y) | h fraction delta corr(y) |
| --- | ---: | ---: | ---: | ---: |
| `sentaurus_psi_sentaurus_qf_sentaurus_mobility` | `+0.01468` | `+0.04764` | `-0.01760` | `+0.01760` |
| `sentaurus_psi_sentaurus_qf_vela_mobility` | `+0.01860` | `+0.02927` | `-0.00461` | `+0.00461` |
| `vela_psi_sentaurus_qf_vela_mobility` | `-0.86193` | `+0.87261` | `+0.00018` | `-0.00018` |
| `vela_psi_shifted_const_sentaurus_qf_vela_mobility` | `-0.86193` | `+0.87261` | `-0.86929` | `+0.86929` |

Interpretation:

- The pure Sentaurus endpoint state remains self-consistent: Vela's SG qF
  formula reproduces Sentaurus current-derived electron and hole components to
  about `~1%`.
- Mobility differences are visible but mostly flat with y.  With Sentaurus
  `psi/qF` and Vela mobility, carrier-fraction L2 is `0.08887` and the split
  bias is roughly `6.3%`, but the y correlation is nearly zero.
- Replacing only `psi` with Vela's electrostatic potential causes a very large
  absolute-density mismatch: electron current is amplified by `~230x` median
  while hole current falls to `~0.006x` median.  This confirms that absolute
  electrostatic potential / qF alignment, not the SG qF stencil itself, is the
  dominant mismatch.
- A best-fit constant Vela `psi` shift of `-0.1365 V` removes the gross
  magnitude error but leaves strong y-correlated component-ratio residuals.
  Therefore a single gauge offset is insufficient; the remaining problem is a
  spatially varying `psi` alignment or Poisson/contact-potential mismatch along
  the active column.

### Next Tasks After Task 119

1. Fit per-support-node or low-order y-dependent `psi` shifts on the same
   active edges and report the required shift versus y.  Compare that shift to
   direct `Vela psi - Sentaurus psi` endpoint offsets.
2. Trace the y-dependent `psi` offset back to Poisson inputs:
   - contact potential / built-in potential convention;
   - fixed charge and dopant ionization sign;
   - permittivity and mesh coupling weights;
   - charge density reconstruction under old-Slotboom BGN.
3. Keep SG carrier flux and avalanche source ownership unchanged until the
   electrostatic alignment check is closed.  Current evidence says SG qF is
   correct under Sentaurus endpoint state and mobility is secondary.

### Execution Note 2026-06-20: Task 120 Per-Node Psi Shift Scan

Task 120 extended:

```text
scripts/diagnose_pn2d_bv_potential_qf_alignment_replay.py
```

It now also writes:

```text
focused_restart_m13p2_potential_qf_alignment_replay/
  potential_qf_alignment_node_shift_scan.csv
```

The scan fits a separate constant Vela `psi` shift per active support node
using the same `Vela psi + Sentaurus qF + Vela mobility` replay.

Summary:

```text
node row_count = 10
best_psi_shift_V median = -0.138
best_psi_shift_V min/max = -0.1395 / -0.131
best_psi_shift_V corr(y) = +0.99693
direct Sentaurus-minus-Vela psi shift median = -0.14142
direct Sentaurus-minus-Vela psi shift min/max = -0.14294 / -0.13406
direct Sentaurus-minus-Vela psi shift corr(y) = +0.99727
per-node fraction L2 median = 0.05330
```

| node | y (um) | best psi shift (V) | direct Sentaurus-Vela psi median (V) | fraction L2 |
| ---: | ---: | ---: | ---: | ---: |
| `352` | `0` | `-0.1395` | `-0.14294` | `0.05331` |
| `351` | `0.015625` | `-0.1395` | `-0.14289` | `0.05320` |
| `354` | `0.03125` | `-0.1395` | `-0.14274` | `0.05312` |
| `360` | `0.046875` | `-0.1390` | `-0.14250` | `0.05344` |
| `384` | `0.078125` | `-0.1385` | `-0.14183` | `0.05317` |
| `392` | `0.109375` | `-0.1375` | `-0.14101` | `0.05347` |
| `708` | `0.265625` | `-0.1335` | `-0.13667` | `0.05335` |
| `716` | `0.296875` | `-0.1325` | `-0.13581` | `0.05316` |
| `740` | `0.328125` | `-0.1315` | `-0.13494` | `0.05328` |
| `748` | `0.359375` | `-0.1310` | `-0.13406` | `0.05363` |

Interpretation:

- The fitted per-node shift and the direct endpoint offset have nearly the
  same y trend.  The required correction is not an arbitrary flux-fitting
  artifact; it is visible directly in the electrostatic potential fields.
- Vela `psi` is higher than Sentaurus by about `0.134 V` to `0.143 V` on these
  active-edge endpoints.  The offset is not constant: it changes by about
  `8.9 mV` from low y to high y.
- Per-node fitting cuts the carrier-fraction L2 from the global constant-shift
  value `0.10743` to about `0.0533`, but it does not reach the pure Sentaurus
  endpoint L2 `0.00527`.  The remaining gap is consistent with the previously
  observed Vela/Sentaurus mobility split and possibly endpoint interpolation
  details, but the first-order defect is electrostatic.

### Next Tasks After Task 120

1. Compare raw Vela and Sentaurus `psi` fields on the same active column and
   neighboring nodes, not only on active edges:
   - direct `psi_sentaurus - psi_vela`;
   - residual after subtracting best constant gauge;
   - relation to doping, space charge, and contact adjacency.
2. Audit Poisson setup against Sentaurus defaults:
   - contact Dirichlet potential convention and built-in potential;
   - sign of ionized donors/acceptors in charge density;
   - old-Slotboom `ni_eff` use in equilibrium carrier initialization;
   - permittivity and finite-volume control-volume/coupling factors.
3. Run a Poisson-only replay or frozen-carrier Poisson residual check using
   Sentaurus carrier densities on the Vela mesh.  If the y-dependent `psi`
   residual is already present there, fix Poisson/contact alignment before any
   DD transport changes.

### Execution Note 2026-06-20: Task 121 Psi Field Alignment and BGN State Audit

Task 121 added:

```text
scripts/diagnose_pn2d_bv_psi_field_alignment.py
```

The script compares Vela and Sentaurus node fields around the active overlap
support column and its mesh-neighbor halo:

```text
focused_restart_m13p2_psi_field_alignment/
  psi_field_alignment_nodes.csv
  psi_field_alignment_summary.json

focused_restart_m13p2_psi_field_alignment_depth2/
  psi_field_alignment_nodes.csv
  psi_field_alignment_summary.json
```

Depth-2 real summary for the existing `-13.2 V` Vela VTK:

```text
row_count = 138
active_seed_psi_delta_median_V = -0.1394471128567183
psi_delta_V median/min/max = -0.138185 / -0.153194 / -0.123025
psi_delta_V corr(y) = +0.47921 over all selected nodes
```

Class summaries:

| node class | count | psi delta median (V) | min (V) | max (V) |
| --- | ---: | ---: | ---: | ---: |
| `active_support` | `10` | `-0.14142` | `-0.14294` | `-0.13406` |
| `active_edge_endpoint` | `20` | `-0.13928` | `-0.14633` | `-0.13067` |
| `neighbor_depth_1` | `57` | `-0.13719` | `-0.14975` | `-0.12685` |
| `neighbor_depth_2` | `51` | `-0.13656` | `-0.15319` | `-0.12303` |

Local plane fit over the depth-2 halo:

```text
psi_sentaurus - psi_vela ~= -0.14367185
                              + 0.43266340 * (x_um - 1.0)
                              + 0.02649348 * y_um
R^2 = 0.99916882
RMSE = 0.00021592 V
residual min/max = -0.000268 / +0.000732 V
```

This proves the active-region `psi` mismatch is a smooth local field mismatch,
not a noisy edge-local SG artifact.

The same report also revealed a stronger root cause in the existing Vela
`-13.2 V` state: it was not using the Sentaurus-aligned material `ni` / BGN
state.  Inferred intrinsic-density ratios from density + `psi/qF`:

| region/class | Vela inferred ni / model | Sentaurus inferred ni / model |
| --- | ---: | ---: |
| active edge endpoints | `~0.60401` | `~0.99999` |
| active support junction nodes | `~0.5096` | `~0.99999` |
| p/n high-doping halo nodes | `~0.6040` | `~0.99999` |

Examples:

- p-side neutral sample node `55`: old Vela `psi=-13.3426`, Sentaurus
  `psi=-13.60365`; Vela inferred `ni~1.0e16`, Sentaurus inferred
  `ni~1.6556e16`.
- active junction node `352`: old Vela inferred `ni/model~0.5096`, while
  Sentaurus inferred `ni/model~0.99999`.

Configuration audit:

- Current source reference config already sets
  `vela_materials_file = pn2d_sentaurus2018_iv_materials.json` for BV.
- The stale generated
  `build-release/reference_tcad/pn2d_sentaurus2018/vela/simulation_bv.json`
  did not contain top-level `materials_file`, so the old BV VTKs were produced
  with built-in default `Si.ni=1.0e16`.
- Re-running `sentaurus_import.py reference --skip-vela-run` refreshed the base
  generated BV deck and added:

```json
"materials_file": "pn2d_sentaurus2018_iv_materials.json"
```

The refresh command failed only at the final comparison step because it tried
to compare against a stale candidate CSV; the deck itself was rewritten.

### Execution Note 2026-06-20: Task 122 Materials/BGN Smoke Verification

Task 122 generated a short verification deck:

```text
materials_bgn_smoke/
  materials_bgn_smoke.json
  materials_bgn_smoke_0000_0V.vtk
  materials_bgn_smoke_0001_-0.05V.vtk
  materials_bgn_smoke.csv
```

Command:

```text
build-release/vela_example_runner.exe --config
  build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/materials_bgn_smoke/materials_bgn_smoke.json
```

Result:

```text
converged = true
points = 2
```

The refreshed deck plus current runner now gives BGN-consistent Vela state:

| sample node | ni_model | Vela inferred ni/model, electron | Vela inferred ni/model, hole |
| ---: | ---: | ---: | ---: |
| `55` p side | `1.65563e16` | `1.000006` | `0.999995` |
| `350` p-side edge | `1.65563e16` | `1.000016` | `0.999981` |
| `352` junction | `1.96229e16` | `0.999999` | `1.000002` |
| `987` n-side edge | `1.65563e16` | `1.000002` | `1.000000` |
| `1319` n side | `1.65563e16` | `0.999995` | `1.000005` |

Contact smoke values:

```text
Anode:   psi=-0.453651, phin=-0.05, phip=-0.05, minority n=2.74112e9
Cathode: psi=+0.403651, phin=0,     phip=0,     minority p=2.74112e9
```

These match the Sentaurus/BGN relation and explain why the stale `-13.2 V`
diagnostics showed a large `psi` and carrier-split mismatch.

### Next Tasks After Task 122

1. Treat old `-13.2 V` Vela VTKs and all diagnostics derived from them as
   contaminated by the missing `materials_file` unless their config explicitly
   contains top-level `materials_file`.
2. Regenerate the focused `-13.2 V` no-impact/default-BV branch from a
   materials-aligned deck.  Do this before changing SG flux, avalanche source,
   mobility, or solver damping.
3. Re-run the existing Task 119/121 diagnostics on the regenerated `-13.2 V`
   VTK:
   - `potential_qf_alignment_replay`;
   - `psi_field_alignment_depth2`;
   - carrier split / SG flux-form active edges.
4. Acceptance criteria for the next rerun:
   - Vela inferred `ni/model` should be near `1.0` on active support and
     high-doping halo nodes.
   - The old local plane `psi` residual (`~0.14 V` median, `~0.26 V` in the
     p-neutral plateau) should collapse substantially.
   - Only remaining component-current mismatch after this rerun should be
     revisited for mobility interpolation or current extraction.

### Execution Note 2026-06-20: Task 123 Materials-Aligned -13.2 V Rerun

Task 123 regenerated the focused BV branch from a materials-aligned deck:

```text
materials_aligned_to_m13p2/
  materials_aligned_to_m13p2.json
  materials_aligned_to_m13p2.csv
  materials_aligned_to_m13p2_0264_-13.2V.vtk
```

Command:

```text
build-release/vela_example_runner.exe --config
  build-release/reference_tcad/pn2d_sentaurus2018/reports/sentaurus_default_bv_execution/materials_aligned_to_m13p2/materials_aligned_to_m13p2.json
```

Result:

```text
converged = true
points = 265
final bias = -13.2 V
final Vela current_total = -5.855755130178548e-17 A/um
Sentaurus current_total = -8.38472088807e-17 A
Vela/Sentaurus terminal current ratio ~= 0.6985
```

The old missing-materials/BGN defect is resolved in this new state:

```text
active/halo Vela inferred ni/model median ~= 0.99998 to 1.00001
old stale run active-support psi_delta median ~= -0.14142 V
new materials-aligned active-support psi_delta median ~= -0.00435 V
new depth-2 halo psi_delta median ~= -0.000617 V
```

The `Vela psi + Sentaurus qF + Vela mobility` replay now shows that mobility is
not the main residual:

```text
sentaurus_psi_sentaurus_qf_vela_mobility fraction L2 = 0.00301
electron fraction median ratio = 1.00121
hole fraction median ratio = 1.00089
```

But the self-consistent Vela qF levels remain offset even though the qF
gradients are right:

```text
active-edge electron qF-field ratio median ~= 1.00078
active-edge hole qF-field ratio median ~= 1.00023
active-edge E-field ratio median ~= 0.99964
Sentaurus - Vela phin median ~= -10.8 mV
Sentaurus - Vela phip median ~= +5.1 mV
Vela/Sentaurus density ratio median ~= 0.79 for electrons and 0.69 for holes
```

Interpretation:

- The old `~0.14 V` electrostatic/BGN error has collapsed; do not tune BGN,
  material `ni`, or low-field mobility from the stale diagnostics.
- The remaining `~30%` terminal-current deficit is now controlled by absolute
  quasi-Fermi levels and the density/current source balance, not by qF gradient,
  SG edge gradient, material `ni`, or SRH lifetime.

### Execution Note 2026-06-20: Task 124 Source-Term Split At -13.2 V

Task 124 compared Vela and Sentaurus generation/recombination on the new
materials-aligned active support, active edge endpoints, and a depth-2 halo:

```text
materials_aligned_m13p2_generation_recombination_alignment/
  generation_recombination_alignment_nodes.csv
  generation_recombination_alignment_summary.json
```

Key active-region integral ratios:

| support set | Vela/Sentaurus avalanche source | Vela/Sentaurus SRH source | Vela/Sentaurus net source |
| --- | ---: | ---: | ---: |
| active support | `0.3833` | `1.000006` | `0.5811` |
| active edge endpoints | `0.3828` | `1.000005` | `0.5633` |
| active combined | `0.3830` | `1.000006` | `0.5695` |
| all selected depth-2 halo | `0.3816` | `1.000006` | `0.5734` |

This is the sharpest current debug signal so far:

- SRH is essentially identical after materials/BGN alignment.
- Vela avalanche generation is only about `0.38x` Sentaurus in the same
  high-field support.
- Vela net generation minus SRH is about `0.57x` Sentaurus, close enough to the
  `0.70x` terminal current ratio to be the next first-order feedback target.
- Therefore the next implementation/debug step should target impact-ionization
  source construction and continuity feedback, not SRH, BGN, or mobility.

Task 124 also regenerated a Python SG edge-source reconstruction:

```text
materials_aligned_m13p2_sg_avalanche_edges/
  sg_avalanche_edges.csv
  sg_avalanche_edge_summary.json
```

That reconstruction found:

```text
VTK AvalancheGeneration total source = 8.1421898e7
Python endpoint-averaged SG reconstruction = 5.4653781e7
relative difference = -32.9%
top source edges are the same x ~= 1.0 um active-column edges
```

This does not invalidate the VTK-vs-Sentaurus source comparison above, but it
means future edge-source diagnostics must use the C++ `SgEdgeCurrentAvalanche`
records directly or make the Python mobility/source reconstruction use exactly
the same `edgeMobility` path as C++.

### Next Tasks After Task 124

1. Add or enable a C++ runner diagnostic CSV for the materials-aligned `-13.2 V`
   run that exports raw `SgEdgeCurrentAvalancheSourceRecord` fields and nodal
   component sums:
   - `electronFluxProxy`, `holeFluxProxy`;
   - `electronAlpha`, `holeAlpha`;
   - electron/hole/combined nodal source integrals;
   - VTK nodal `AvalancheGeneration * node_volume`.
2. Re-run the active-support source comparison using those C++ records, not the
   Python endpoint-averaged reconstruction.  Acceptance criterion: C++ record
   nodal sums must reproduce VTK `AvalancheGeneration * node_volume`.
3. Compare C++ record-derived weighted alpha and particle flux against
   Sentaurus exported `ImpactIonization`, `eCurrentDensity`, and
   `hCurrentDensity` on the active support:
   - if Vela particle flux is low but weighted alpha matches Sentaurus, focus
     on continuity/qF absolute-level feedback;
   - if particle flux matches but weighted alpha is low, focus on
     VanOverstraeten coefficient/driving-force interpolation;
   - if both are low, separate density-caused current deficit from alpha-caused
     source deficit with a Sentaurus-qF replay that computes C++ edge records
     from imported Sentaurus fields.
4. Only after the C++ source-record audit, test one-factor probes:
   - Vela with an avalanche source scale chosen from the active-combined
     `0.383x` ratio, as a feedback sensitivity probe only;
   - Vela with Sentaurus-qF frozen/replayed edge source, to test whether the
     absolute qF offsets close when the source term is forced to Sentaurus.
5. Keep the acceptance target unchanged: Sentaurus default
   `Avalanche(VanOverstraeten)` with default SG edge-current avalanche, not the
   `AvalDensGradQF` control path.

### Execution Note 2026-06-20: Task 125 C++ SG Avalanche Record Audit

Task 125 did not add new code.  It used the existing DCSweep diagnostic:

```json
"sweep": {
  "diagnostics": {
    "sg_avalanche_edges": {
      "enabled": true
    }
  }
}
```

The materials-aligned `-13.2 V` sweep was rerun with VTK output disabled:

```text
materials_aligned_m13p2_cxx_sg_records/
  materials_aligned_m13p2_cxx_sg_records.json
  materials_aligned_m13p2_cxx_sg_records.csv
  materials_aligned_m13p2_cxx_sg_edges.csv
  materials_aligned_m13p2_cxx_sg_record_alignment_nodes.csv
  materials_aligned_m13p2_cxx_sg_record_alignment_summary.json
  materials_aligned_m13p2_cxx_sg_record_factorization_summary.json
```

Result:

```text
converged = true
points = 265
final point_index = 264
C++ edge total source = 8.142190381e7
VTK AvalancheGeneration * node_volume total = 8.142189798e7
C++/VTK total source = 1.00000007
```

This closes the previous Python-reconstruction discrepancy: the C++ assembled
SG avalanche source and the VTK nodal `AvalancheGeneration` field are
consistent.  The earlier `-32.9%` Python edge reconstruction came from using a
different endpoint-averaged mobility/source reconstruction and should not be
used for final source parity claims.

Active-region comparison against Sentaurus at `-13.2 V`:

| support set | C++/VTK source | VTK/Sentaurus source | C++/Sentaurus source | C++ particle flux / Sentaurus | C++ weighted alpha / Sentaurus |
| --- | ---: | ---: | ---: | ---: | ---: |
| active support | `0.9999993` | `0.3833353` | `0.3833350` | `1.5091` median | `1.0124` median |
| active edge endpoints | `1.0000000` | `0.3827837` | `0.3827837` | `1.4947` median | `1.0266` median |
| active combined | `0.9999997` | `0.3829721` | `0.3829720` | `1.4962` median | `1.0262` median |
| all selected depth-2 halo | `0.9999999` | `0.3815992` | `0.3815992` | `1.4970` median | `1.0198` median |

Factorization:

```text
C++/Sentaurus source
  ~= (C++ particle flux / Sentaurus particle flux)
    * (C++ weighted alpha / Sentaurus weighted alpha)
    * residual_geometry_factor

residual_geometry_factor median = 0.25
residual_geometry_factor min/max = 0.25 / 0.25 on active and halo sets
```

Interpretation:

- Vela is not under-generating avalanche because SG particle flux is too low;
  the C++ edge particle flux proxy is about `1.5x` Sentaurus on these supports.
- Vela is not under-generating because Van Overstraeten `alpha` is too low;
  the inferred weighted-alpha ratio is close to `1.0`.
- The entire `~0.383x` source deficit is explained by a stable `0.25`
  geometry/control-volume factor after accounting for flux and alpha.  This is
  now the highest-priority root-cause candidate.

### Next Tasks After Task 125

1. Audit the Vela SG avalanche source geometry against the Sentaurus manual's
   default element-edge current approximation:
   - current Vela record uses `edgeAreaProxy = 0.5 * edge.length * edge.couple`;
   - nodal source uses `0.5 * edgeSourceIntegral` per endpoint;
   - together this gives an effective `0.25 * edge.length * edge.couple`
     contribution per endpoint.
2. Build a one-factor diagnostic branch, not yet an accepted fix, that exposes
   an `sg_avalanche_geometry_scale` or equivalent source-record policy and tests
   scales `2x` and `4x` on the materials-aligned `-13.2 V` branch.
   Acceptance for this probe is diagnostic only:
   - `4x` should move the active avalanche source close to Sentaurus;
   - the terminal current should move toward the Sentaurus `-13.2 V` current
     without destabilizing the continuation.
3. In parallel, compare the Vela box geometry exported from `mesh.json` against
   Sentaurus element/control-volume geometry if the TDR export contains vertex
   box volumes or if Sentaurus can export them through the VM.  The goal is to
   decide whether Vela's edge-source volume policy is truly missing a factor,
   or whether Sentaurus `ImpactIonization` nodal CSV is a field interpolation
   that should not be multiplied by Vela node volume for source parity.
4. Do not change Van Overstraeten coefficients, SRH lifetime, mobility, BGN, or
   Bank-Rose damping based on Task 125.  Those have either matched or are not
   the dominant source deficit in the materials-aligned `-13.2 V` evidence.

### Execution Note 2026-06-20: Task 126 Alpha-Scale Probe Bracket

Task 126 tested whether simply scaling the Van Overstraeten prefactors can act
as a one-factor source-strength probe.  This is diagnostic only; it is not an
accepted physical fix for Sentaurus parity.

Generated:

```text
materials_aligned_m13p2_alpha_scale_probe/
  restart_from_materials_aligned_m13p2.csv
  alpha_scale_4/
  alpha_scale_2/
  alpha_scale_4_full/
  alpha_scale_2_full/
```

Single-point restarts at `-13.2 V` from the materials-aligned baseline state:

| probe | result | failure |
| --- | --- | --- |
| `alpha_scale_4` | failed at the target point | `line_search_non_decrease` after 3 Newton iterations |
| `alpha_scale_2` | failed at the target point | `line_search_non_decrease` |

Full continuation from `0 V`:

| probe | converged points | last converged bias | last current | failure |
| --- | ---: | ---: | ---: | --- |
| `alpha_scale_4_full` | `62` converged + one failure row | `-3.05 V` | `+6.84e-17 A/um` | `max_iterations` near `-3.050000000186 V` |
| `alpha_scale_2_full` | `121` converged + one failure row | `-3.745782673 V` | `+3.36e-15 A/um` | `line_search_non_decrease` near `-3.74578267336 V` |

Interpretation:

- Prefactor scaling does verify that increasing avalanche feedback strongly
  moves the branch: both `2x` and `4x` destabilize much earlier than the
  baseline `-13.2 V` run, with current sign flips and large jump ratios.
- Therefore the `0.25` geometry factor from Task 125 should not be patched by
  tuning Van Overstraeten coefficients or alpha prefactors.
- The next implementation should expose or correct the SG avalanche source
  geometry itself, then retest continuation stability.  A raw alpha/source
  multiplier is useful only as a sensitivity probe.

### Current Physical-Quantity Status After Task 126

Materials-aligned baseline IV error:

| bias | Vela current | Sentaurus current | Vela/Sentaurus |
| ---: | ---: | ---: | ---: |
| `-0.5 V` | `-5.5737e-18 A/um` | `-5.3460e-18 A` | `1.0426` |
| `-2 V` | `-1.4955e-17 A/um` | `-1.4176e-17 A` | `1.0550` |
| `-5 V` | `-2.9537e-17 A/um` | `-2.8427e-17 A` | `1.0391` |
| `-10 V` | `-4.4093e-17 A/um` | `-5.4542e-17 A` | `0.8084` |
| `-13.2 V` | `-5.8558e-17 A/um` | `-8.3847e-17 A` | `0.6984` |

This means the low-reverse-bias transport is now close after the materials/BGN
fix.  The large error appears with high-field avalanche feedback.

Field and carrier state at `-13.2 V`:

```text
max electric field, Vela = 4.5878e5 V/cm
active/halo psi_delta median = -0.000617 V
active-support psi_delta median = -0.00435 V
Vela inferred ni/model ~= 1.0
```

### Execution Note 2026-06-20: Task 127 SG Source Geometry Scale Probe

Task 127 added a default-off diagnostic `impact_ionization.source_geometry_scale`
knob to scale the Scharfetter-Gummel edge-current avalanche source geometry.
The default is `1.0`, so existing decks retain the old behavior unless they
explicitly opt in.

Implementation and tests:

- Added `ImpactIonizationModelConfig::sourceGeometryScale`.
- Parsed and validated `source_geometry_scale` in Gummel and Newton solver JSON.
- Applied the multiplier to the SG edge source geometry in
  `sgEdgeCurrentAvalancheSourceRecords`, which also feeds the nodal source
  integrals used by the continuity equations.
- Added tests that verify JSON parsing and that both edge records and nodal
  SG avalanche source integrals scale linearly.
- Verified:
  - `build-release\test_impact_ionization.exe "SG edge current avalanche source supports diagnostic geometry scale"`
  - `build-release\test_impact_ionization.exe "JSON solver config selects impact ionization model"`
  - `build-release\test_dc_sweep.exe "DCSweep: SG avalanche edge diagnostics write assembled source rows"`
  - `build-release\test_impact_ionization.exe "~[gummel]" --reporter compact`

Full `test_impact_ionization.exe` still exits `1` because the pre-existing
`[gummel]` case `Gummel reverse bias BV regression runs with impact ionization`
exits silently when run alone.  Excluding `[gummel]`, the remaining impact tests
pass.

Generated:

```text
materials_aligned_m13p2_source_geometry_scale_probe/
  geometry_scale_2_full/
  geometry_scale_4_full/
```

Both full sweeps reached `-13.2 V`.  The `4x` probe needed more accepted points
(`749`) than the baseline/`2x` runs (`265`), consistent with stronger avalanche
feedback and step-size contraction.

Current comparison:

| case | bias | Vela current | Sentaurus current | Vela/Sentaurus |
| --- | ---: | ---: | ---: | ---: |
| baseline | `-0.5 V` | `-5.5737e-18 A/um` | `-5.3460e-18 A` | `1.0426` |
| baseline | `-2 V` | `-1.4955e-17 A/um` | `-1.4176e-17 A` | `1.0550` |
| baseline | `-5 V` | `-2.9537e-17 A/um` | `-2.8427e-17 A` | `1.0390` |
| baseline | `-10 V` | `-4.4093e-17 A/um` | `-5.4542e-17 A` | `0.8084` |
| baseline | `-13.2 V` | `-5.8558e-17 A/um` | `-8.3847e-17 A` | `0.6984` |
| geometry scale `2x` | `-0.5 V` | `-5.6669e-18 A/um` | `-5.3460e-18 A` | `1.0600` |
| geometry scale `2x` | `-2 V` | `-1.4959e-17 A/um` | `-1.4176e-17 A` | `1.0552` |
| geometry scale `2x` | `-5 V` | `-2.9538e-17 A/um` | `-2.8427e-17 A` | `1.0391` |
| geometry scale `2x` | `-10 V` | `-5.8628e-17 A/um` | `-5.4542e-17 A` | `1.0749` |
| geometry scale `2x` | `-13.2 V` | `-8.7579e-17 A/um` | `-8.3847e-17 A` | `1.0445` |
| geometry scale `4x` | `-0.5 V` | `-5.6210e-18 A/um` | `-5.3460e-18 A` | `1.0514` |
| geometry scale `4x` | `-2 V` | `-1.4959e-17 A/um` | `-1.4176e-17 A` | `1.0552` |
| geometry scale `4x` | `-5 V` | `-2.9538e-17 A/um` | `-2.8427e-17 A` | `1.0391` |
| geometry scale `4x` | `-10 V` | `-7.3186e-17 A/um` | `-5.4542e-17 A` | `1.3418` |
| geometry scale `4x` | `-13.2 V` | `-2.1853e-16 A/um` | `-8.3847e-17 A` | `2.6063` |

Observed `max_electric_field_V_per_cm` stayed the same at the sampled biases
within CSV precision (`4.58776e5 V/cm` at `-13.2 V` for all three cases).  The
primary change is therefore in carrier continuity avalanche feedback and
terminal current, not in the Poisson/electric-field branch.

Interpretation:

- The old materials/BGN issue is resolved: low-reverse-bias current remains
  close to Sentaurus, and the remaining mismatch starts in high-field avalanche
  feedback.
- A `2x` SG source-geometry multiplier almost closes the `-13.2 V` current gap,
  while `4x` overshoots badly.  This strongly supports a remaining factor-of-two
  source-placement/control-volume convention mismatch, rather than a Van
  Overstraeten coefficient problem.
- The previous `2x`/`4x` alpha-prefactor probes failed early, but the source
  geometry probes converge.  This separates "coefficient scaling" from "where
  and how the SG edge generation is integrated into continuity".
- Next debug direction: replace the diagnostic multiplier with a physically
  named Sentaurus-parity source placement policy.  Check whether Sentaurus uses
  full edge box area, endpoint splitting, element-volume avalanche, or contact/
  boundary-layer edge treatment that makes Vela's current `0.5*h*couple` source
  geometry effectively too small by about `2x` on the active high-field support.

The replay diagnostics show that using Sentaurus `psi/qF` with Vela mobility is
already near-perfect:

```text
sentaurus_psi_sentaurus_qf_vela_mobility fraction L2 = 0.00301
electron fraction median ratio = 1.00121
hole fraction median ratio = 1.00089
```

With Vela self-consistent `psi` and Sentaurus qF, the remaining carrier split is
moderate:

```text
electron ratio median = 1.1877
hole ratio median = 0.8487
particle ratio median = 1.0095
```

Source-term status at `-13.2 V`:

```text
SRH source ratio ~= 1.000006
avalanche source ratio ~= 0.383 on active support
net generation minus SRH source ratio ~= 0.57
C++ SG edge records reproduce VTK AvalancheGeneration exactly
C++ particle flux / Sentaurus particle flux ~= 1.50
C++ weighted alpha / Sentaurus weighted alpha ~= 1.02
residual geometry factor ~= 0.25 exactly
```

Working cause hypothesis:

- BGN/material `ni`, SRH lifetime, mobility, qF gradient, and Van Overstraeten
  alpha formula are no longer first-order suspects.
- The leading candidate is an SG avalanche source geometry/control-volume
  policy mismatch: Vela currently computes an effective per-endpoint source
  contribution of `0.25 * edge.length * edge.couple` after applying
  `edgeAreaProxy = 0.5 * length * couple` and endpoint half-splitting.
- Because direct alpha scaling destabilizes the continuation, the next fix
  should target the source discretization geometry and then re-evaluate IV and
  source parity, not tune physical coefficients.

### Execution Note 2026-06-20: Branch Consolidation and Regression Gate

Before committing the accumulated BV debug branch, the local release build and
all regression unittests were refreshed:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build-release -j 4
python -m unittest discover -s tests/regression -p "test*.py"
```

Observed result:

```text
build-release build: exit 0
regression tests: Ran 184 tests in 112.883s, OK (skipped=1)
```

Current debug workflow from this checkpoint:

1. Keep `impact_ionization.source_geometry_scale` as a diagnostic knob only.
2. Promote the `2x` behavior into a named, TDD-covered source-geometry policy
   only after confirming whether Sentaurus parity is due to full edge box area,
   endpoint ownership, element-volume avalanche, or contact/boundary-layer edge
   treatment.
3. Re-run the materials-aligned `0 V` to `-13.2 V` BV sweep with VTK and SG
   edge diagnostics for the selected policy.
4. Compare IV, electrostatic potential, electric field, electron density, hole
   density, SRH, and avalanche source fields at `0`, `-0.5`, `-2`, `-5`, `-10`,
   and `-13.2 V` before extending back to `-20 V`.
