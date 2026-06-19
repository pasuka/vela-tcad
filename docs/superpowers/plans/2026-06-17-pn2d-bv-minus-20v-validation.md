# PN2D BV -20 V Validation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reproduce the pn2d Sentaurus2018 default BV avalanche semantics in Vela, extend the validation from the current low-reverse-bias gate to a 0 V to -20 V run only after the Sentaurus-default source path passes high-bias parity checks, then compare BV current and same-bias spatial fields against Sentaurus.

**Architecture:** Use the existing `reference_tcad/pn2d_sentaurus2018` fixture as the source of truth, with Sentaurus's default Scharfetter-Gummel edge-current avalanche discretization as the required Vela parity target. Keep Sentaurus field export and Vela execution separate, then join them through `scripts/compare_pn2d_bv_multibias_fields.py` so curve error, field error, and avalanche-source error are ranked at the same bias points. Treat `current_approximation = "density_gradient"` as the Sentaurus-default path; keep `mobility_density_gradient` only as an `AvalDensGradQF`-like control path.

**Tech Stack:** C++20, CMake/Ninja, MSYS2 UCRT64 on Windows, Python standard library plus NumPy for analysis scripts, existing Vela runner and Sentaurus import tooling.

---

## Current Branch Baseline

- Latest commit: `768ba5ee092dad231c145758d2709d8cbfe1b1d3` (`Align pn2d BV Sentaurus physics diagnostics`).
- Branch: `codex-pn2d-sentaurus2018-calibration`.
- Worktree status observed during planning: clean.
- The commit adds `scripts/compare_pn2d_bv_multibias_fields.py` and `scripts/diagnose_pn2d_bv_mobility.py`, expands impact-ionization support, and updates the pn2d Sentaurus2018 reference JSON.
- Sentaurus BV deck `reference_tcad/pn2d_sentaurus2018/source/pn2d_bv_sdevice.cmd` already targets `Goal { Name="Anode" Voltage=-20.0 }` and writes `pn2d_bv_multibias` snapshots over 200 intervals.
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

## Final Verification Commands

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\test_impact_ionization.exe
ctest --test-dir build-release --output-on-failure -R "impact_ionization|newton|dc_sweep|reference_tcad|sentaurus_import"
python -m unittest tests.regression.test_reference_tcad_tools tests.regression.test_sentaurus_import_tools
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
