# PN2D BV -20 V Validation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend the pn2d Sentaurus2018 BV validation from the current low-reverse-bias gate to a 0 V to -20 V Vela run, then compare BV current and same-bias spatial fields against Sentaurus.

**Architecture:** Use the existing `reference_tcad/pn2d_sentaurus2018` fixture as the source of truth. Keep Sentaurus field export and Vela execution separate, then join them through `scripts/compare_pn2d_bv_multibias_fields.py` so curve error and field error are ranked at the same bias points. Treat impact-ionization model parity as a first-class gate before tuning solver continuation.

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

Likely parity gap:
  Vela current_density generation is closer to Sentaurus Math { AvalDensGradQF }
  than to Sentaurus default SG edge-current avalanche discretization.

Current -13.208 V failure priority:
  1. Avalanche current-density discretization parity.
  2. Avalanche hotspot geometry/control-volume amplification.
  3. Low-density driving-force interpolation to ElectricField.
  4. Analytic Jacobian completeness for avalanche driving-field and mobility terms.
  5. BandgapDependence only if a future Sentaurus deck explicitly enables it.
```

## Files To Modify Or Use

- Modify: `reference_tcad/pn2d_sentaurus2018/pn2d_sentaurus2018_reference.json`
  - Extend BV run target to `vela_stop: -20.0`.
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
- Test: `tests/test_impact_ionization.cpp`
  - Add coefficient and generation checks only if the formula review finds a mismatch.
- Test: `tests/test_newton_solver.cpp` and `tests/test_dc_sweep.cpp`
  - Add parser/restart coverage for any new avalanche discretization or driving-force-interpolation knobs.
- Test: `tests/regression/test_reference_tcad_tools.py`
  - Add regression coverage if command-line behavior or generated report schema changes.

## Task 1: Reproduce The Imported Reference And Existing Gate

- [ ] **Step 1: Configure the MSYS2 UCRT64 build**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --preset windows-ucrt64-release
```

Expected: CMake configures `build-release/` successfully and finds Ninja plus the UCRT64 compiler.

- [ ] **Step 2: Build the required tools**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --parallel --target vela_example_runner sentaurus_import test_impact_ionization
```

Expected: `build-release/vela_example_runner.exe`, `build-release/sentaurus_import.exe`, and `build-release/test_impact_ionization.exe` exist.

- [ ] **Step 3: Run the impact-ionization unit tests before changing BV reach**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
   build-release\test_impact_ionization.exe
```

Expected: `All tests passed`, including Van Overstraeten coefficient checks.

- [ ] **Step 4: Regenerate the pn2d Sentaurus2018 imported reference tree**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts\sentaurus_import.py reference --config reference_tcad\pn2d_sentaurus2018\pn2d_sentaurus2018_reference.json --source-dir reference_tcad\pn2d_sentaurus2018\source --output-dir build-release\reference_tcad\pn2d_sentaurus2018 --tdr-importer build-release\sentaurus_import.exe --runner build-release\vela_example_runner.exe
```

Expected: current low-bias BV gate still runs, and `build-release/reference_tcad/pn2d_sentaurus2018/reference_curves/pn2d_sentaurus2018_bv_reference.csv` exists.

## Task 2: Import Sentaurus BV Multi-Bias Field Snapshots

- [ ] **Step 1: Confirm the source snapshots exist**

Run:

```powershell
Get-ChildItem reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_multibias_*_des.tdr | Measure-Object
```

Expected: count is at least `201`, covering normalized times from 0 to 1.

- [ ] **Step 2: Export the comparison bias points**

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

- [ ] **Step 3: Verify the imported Sentaurus fields needed by the comparison**

Run:

```powershell
Get-ChildItem build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias\sentaurus_-20v\fields
```

Expected: files exist for `ElectrostaticPotential`, `ElectricField`, `eDensity`, `hDensity`, and `AvalancheGeneration` or `ImpactIonization`.

## Task 3: Create A Controlled Vela -20 V Probe Deck

- [ ] **Step 1: Copy the generated BV deck to a probe file**

Run:

```powershell
Copy-Item build-release\reference_tcad\pn2d_sentaurus2018\vela\simulation_bv.json build-release\reference_tcad\pn2d_sentaurus2018\vela\simulation_bv_minus20_probe.json
```

Expected: the probe deck exists and the committed reference JSON remains unchanged during the first exploration.

- [ ] **Step 2: Edit the probe deck for a conservative reach test**

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

- [ ] **Step 3: Run the Vela -20 V probe**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\vela_example_runner.exe --config build-release\reference_tcad\pn2d_sentaurus2018\vela\simulation_bv_minus20_probe.json
```

Expected: the run either reaches -20 V with VTK files at comparison bias points, or stops with a clear last-stable bias and failure reason in the BV CSV.

- [ ] **Step 4: If the probe fails before -20 V, rerun with staged continuation**

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

## Task 4: Compare BV Curve And Spatial Fields

- [ ] **Step 1: Run the multibias comparison script**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts\compare_pn2d_bv_multibias_fields.py --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\vela --curve-reference build-release\reference_tcad\pn2d_sentaurus2018\reference_curves\pn2d_sentaurus2018_bv_reference.csv --curve-candidate build-release\reference_tcad\pn2d_sentaurus2018\vela\pn2d_sentaurus2018_bv.csv --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_multibias --biases 0,-0.5,-2,-5,-10,-20
```

Expected: `curve_compare.csv`, `field_compare.csv`, `debug_ranking.json`, and `README.md` are written.

- [ ] **Step 2: Read the ranked failure order**

Run:

```powershell
Get-Content build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_multibias\debug_ranking.json
```

Expected: the first ranked items identify whether the dominant mismatch is curve current, potential, electric field, carrier density, mobility, SRH recombination, or avalanche generation.

- [ ] **Step 3: Promote only fields that are same-bias and present on both sides**

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

- [x] **Step 5: Add a Sentaurus-default SG edge-current avalanche probe before promoting Vela's AvalDensGradQF-like path**

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
  Sentaurus-default parity target; derive avalanche current contribution from the same SG edge-current approximation used by the drift-diffusion fluxes, then accumulate Eq. 431 by element/node.
```

Acceptance: do not promote `vela_stop: -20.0` until the selected Vela path is explicitly documented as either `Sentaurus default SG avalanche` or `AvalDensGradQF-equivalent`, and the Sentaurus deck is made explicit if the latter is chosen.

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

Expected: if finite difference reaches farther after volume/current/interpolation probes are controlled, the analytic avalanche Jacobian must be extended to include driving-field and mobility derivatives before Task 6.

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
- The current Vela `generation = "current_density"` implementation is closest to the SDevice manual's `Math { AvalDensGradQF }` approximation, while the Sentaurus deck is using the default SG edge-current avalanche path. This remains the primary parity gap.

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

- Implemented `solver.impact_ionization.current_approximation = "density_gradient"` as the Sentaurus-default SG edge-current probe. The existing Vela node-local path remains the default `mobility_density_gradient`.
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
- Do not promote `density_gradient` directly until the high-field current overshoot is localized with a field comparison or source-volume sensitivity probe.
- Next best task: construct a no-impact quasi-Fermi/continuity residual diagnostic across `-10 V` and `-13.2 V`, especially contact-adjacent and plateau quasi-Fermi absolute levels. The goal is to determine why the no-impact continuity solve selects much higher carrier-density exponents even when the electrostatic potential branch is unchanged by mobility, contact relaxation, and stable SRH lifetime changes.
- Separately add a C++ edge-source dump or VTK contact-node guard to resolve the `-20 V` contact-node diagnostic inconsistency.

## Task 6: Promote The Stable -20 V Path

- [ ] **Step 1: Choose and document the avalanche parity target before promotion**

Promotion is blocked until one of these two paths is explicitly selected:

```text
Path A: Sentaurus-default SG avalanche parity
  Vela impact_ionization.current_approximation = "density_gradient"
  Sentaurus source deck remains Recombination(Avalanche(VanOverstraeten))

Path B: AvalDensGradQF-equivalent parity
  Vela impact_ionization.current_approximation = "mobility_density_gradient"
  Sentaurus source deck or validation note explicitly states Math { AvalDensGradQF }
```

Expected: the validation documentation must not compare an implicit Sentaurus-default avalanche run against an undocumented Vela AvalDensGradQF-like implementation.

- [ ] **Step 2: Update the committed reference JSON only after the selected probe is stable**

Modify `reference_tcad/pn2d_sentaurus2018/pn2d_sentaurus2018_reference.json` BV block:

```json
"vela_stop": -20.0,
"vela_step": -0.1
```

Use `-0.05` instead of `-0.1` only if the selected Task 5 probe proves `-0.1` is not stable. Include any selected non-default avalanche options in the generated BV deck rather than relying on manual probe files.

- [ ] **Step 3: Add or update regression assertions for the promoted path**

In `tests/regression/test_reference_tcad_tools.py`, assert that the generated BV deck includes:

```python
assert bv["solver"]["impact_ionization"]["model"] == "van_overstraeten"
assert bv["solver"]["impact_ionization"]["driving_force"] == "quasi_fermi_gradient"
assert bv["solver"]["impact_ionization"]["generation"] == "current_density"
assert bv["solver"]["impact_ionization"]["current_approximation"] in {
    "density_gradient",
    "mobility_density_gradient",
}
assert bv["sweep"]["stop"] == -20.0
```

Expected: future imports cannot silently fall back to a low-bias BV gate, a no-avalanche deck, or an undocumented avalanche-current approximation.

- [ ] **Step 4: Run focused regression**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\test_impact_ionization.exe
ctest --test-dir build-release --output-on-failure -R "reference_tcad|sentaurus_import"
```

Expected: focused tests pass.

- [ ] **Step 5: Update validation documentation with the final evidence**

Append a section to `docs/validation/pn2d_sentaurus_comparison.md` containing:

```text
PN2D BV -20 V validation date:
Vela run deck:
Sentaurus multibias export root:
Compared biases:
Curve max abs log10 error:
Dominant field mismatch:
Impact-ionization model:
Avalanche current approximation:
Avalanche volume policy:
Driving-force interpolation:
Accepted limitations:
```

Expected: the next worker can start from the generated report paths instead of rediscovering the run.

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
