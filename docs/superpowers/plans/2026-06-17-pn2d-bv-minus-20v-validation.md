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
- Test: `tests/test_impact_ionization.cpp`
  - Add coefficient and generation checks only if the formula review finds a mismatch.
- Test: `tests/regression/test_reference_tcad_tools.py`
  - Add regression coverage if command-line behavior or generated report schema changes.

## Task 1: Reproduce The Imported Reference And Existing Gate

- [ ] **Step 1: Configure the MSYS2 UCRT64 build**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --preset windows-ucrt64-debug
```

Expected: CMake configures `build/` successfully and finds Ninja plus the UCRT64 compiler.

- [ ] **Step 2: Build the required tools**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build --parallel --target vela_example_runner sentaurus_import test_impact_ionization
```

Expected: `build/vela_example_runner.exe`, `build/sentaurus_import.exe`, and `build/test_impact_ionization.exe` exist.

- [ ] **Step 3: Run the impact-ionization unit tests before changing BV reach**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
   build\test_impact_ionization.exe
```

Expected: `All tests passed`, including Van Overstraeten coefficient checks.

- [ ] **Step 4: Regenerate the pn2d Sentaurus2018 imported reference tree**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts\sentaurus_import.py reference --config reference_tcad\pn2d_sentaurus2018\pn2d_sentaurus2018_reference.json --source-dir reference_tcad\pn2d_sentaurus2018\source --output-dir build\reference_tcad\pn2d_sentaurus2018 --tdr-importer build\sentaurus_import.exe --runner build\vela_example_runner.exe
```

Expected: current low-bias BV gate still runs, and `build/reference_tcad/pn2d_sentaurus2018/reference_curves/pn2d_sentaurus2018_bv_reference.csv` exists.

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
build\sentaurus_import.exe --tdr reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_multibias_0000_des.tdr --export-dir build\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias\sentaurus_0v
```

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build\sentaurus_import.exe --tdr reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_multibias_0005_des.tdr --export-dir build\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias\sentaurus_-0.5v
```

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build\sentaurus_import.exe --tdr reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_multibias_0020_des.tdr --export-dir build\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias\sentaurus_-2v
```

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build\sentaurus_import.exe --tdr reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_multibias_0050_des.tdr --export-dir build\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias\sentaurus_-5v
```

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build\sentaurus_import.exe --tdr reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_multibias_0100_des.tdr --export-dir build\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias\sentaurus_-10v
```

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build\sentaurus_import.exe --tdr reference_tcad\pn2d_sentaurus2018\source\pn2d_bv_multibias_0200_des.tdr --export-dir build\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias\sentaurus_-20v
```

Expected: each output directory contains `nodes.csv`, `elements.csv`, and `fields/*.csv`.

- [ ] **Step 3: Verify the imported Sentaurus fields needed by the comparison**

Run:

```powershell
Get-ChildItem build\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias\sentaurus_-20v\fields
```

Expected: files exist for `ElectrostaticPotential`, `ElectricField`, `eDensity`, `hDensity`, and `AvalancheGeneration` or `ImpactIonization`.

## Task 3: Create A Controlled Vela -20 V Probe Deck

- [ ] **Step 1: Copy the generated BV deck to a probe file**

Run:

```powershell
Copy-Item build\reference_tcad\pn2d_sentaurus2018\vela\simulation_bv.json build\reference_tcad\pn2d_sentaurus2018\vela\simulation_bv_minus20_probe.json
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
build\vela_example_runner.exe --config build\reference_tcad\pn2d_sentaurus2018\vela\simulation_bv_minus20_probe.json
```

Expected: the run either reaches -20 V with VTK files at comparison bias points, or stops with a clear last-stable bias and failure reason in the BV CSV.

- [ ] **Step 4: If the probe fails before -20 V, rerun with staged continuation**

Edit only the probe deck step size and rerun:

```json
"step": -0.05
```

Expected: the last-stable bias improves or the failure mode stays identical, which separates continuation stiffness from step-size artifacts.

## Task 4: Compare BV Curve And Spatial Fields

- [ ] **Step 1: Run the multibias comparison script**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts\compare_pn2d_bv_multibias_fields.py --sentaurus-root build\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias --vela-vtk-root build\reference_tcad\pn2d_sentaurus2018\vela --curve-reference build\reference_tcad\pn2d_sentaurus2018\reference_curves\pn2d_sentaurus2018_bv_reference.csv --curve-candidate build\reference_tcad\pn2d_sentaurus2018\vela\pn2d_bv.csv --out-dir build\reference_tcad\pn2d_sentaurus2018\reports\bv_multibias --biases 0,-0.5,-2,-5,-10,-20
```

Expected: `curve_compare.csv`, `field_compare.csv`, `debug_ranking.json`, and `README.md` are written.

- [ ] **Step 2: Read the ranked failure order**

Run:

```powershell
Get-Content build\reference_tcad\pn2d_sentaurus2018\reports\bv_multibias\debug_ranking.json
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

## Task 5: Localize The First Failing Physics Lever

- [ ] **Step 1: Run mobility decomposition if mobility or field drive ranks high**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts\diagnose_pn2d_bv_mobility.py --sentaurus-root build\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias --vela-vtk-root build\reference_tcad\pn2d_sentaurus2018\vela --out-dir build\reference_tcad\pn2d_sentaurus2018\reports\bv_mobility --biases 0,-0.5,-2,-5,-10,-20
```

Expected: output ranks whether `masetti_field`, quasi-Fermi high-field drive, or local mobility saturation is the leading mismatch.

- [ ] **Step 2: If avalanche generation ranks high, isolate model parity before solver tuning**

Add a focused test in `tests/test_impact_ionization.cpp` only if the current defaults fail a hand calculation:

```cpp
TEST_CASE("Van Overstraeten temperature gamma adjusts ionization coefficients",
          "[impact][van_overstraeten]")
{
    ImpactIonizationModelConfig config = impactIonizationModelConfig("van_overstraeten");
    config.temperature_K = 400.0;
    VanOverstraetenImpactIonization model(config);

    constexpr Real kBoltzmann_eV_per_K = 8.617333262145e-5;
    const Real field = 5.0e7;
    const Real gamma = std::tanh(0.063 / (2.0 * kBoltzmann_eV_per_K * 300.0)) /
                       std::tanh(0.063 / (2.0 * kBoltzmann_eV_per_K * 400.0));
    const Real expected = gamma * 7.03e7 * std::exp(-(gamma * 1.231e8) / field);

    REQUIRE(model.electronCoefficient(field) == Catch::Approx(expected).epsilon(1.0e-12));
}
```

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build\test_impact_ionization.exe
```

Expected: the test passes before any tuning; otherwise fix formula or unit conversion first.

- [ ] **Step 3: If the run fails numerically before field mismatch analysis, isolate Jacobian sensitivity**

Create two probe decks differing only in Newton Jacobian mode:

```json
"jacobian": "analytic"
```

and

```json
"jacobian": "finite_difference"
```

Run both to the same target bias and compare last-stable bias.

Expected: if finite difference reaches farther, the analytic avalanche Jacobian is the first implementation target; if both fail at the same bias, continuation/physics stiffness is the first target.

## Task 6: Promote The Stable -20 V Path

- [ ] **Step 1: Update the committed reference JSON only after the probe is stable**

Modify `reference_tcad/pn2d_sentaurus2018/pn2d_sentaurus2018_reference.json` BV block:

```json
"vela_stop": -20.0,
"vela_step": -0.1
```

Use `-0.05` instead of `-0.1` only if Task 3 proves `-0.1` is not stable.

- [ ] **Step 2: Add or update regression assertions for the promoted path**

In `tests/regression/test_reference_tcad_tools.py`, assert that the generated BV deck includes:

```python
assert bv["solver"]["impact_ionization"]["model"] == "van_overstraeten"
assert bv["solver"]["impact_ionization"]["driving_force"] == "quasi_fermi_gradient"
assert bv["solver"]["impact_ionization"]["generation"] == "current_density"
assert bv["sweep"]["stop"] == -20.0
```

Expected: future imports cannot silently fall back to a low-bias BV gate or a no-avalanche deck.

- [ ] **Step 3: Run focused regression**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build\test_impact_ionization.exe
ctest --test-dir build --output-on-failure -R "reference_tcad|sentaurus_import"
```

Expected: focused tests pass.

- [ ] **Step 4: Update validation documentation with the final evidence**

Append a section to `docs/validation/pn2d_sentaurus_comparison.md` containing:

```text
PN2D BV -20 V validation date:
Vela run deck:
Sentaurus multibias export root:
Compared biases:
Curve max abs log10 error:
Dominant field mismatch:
Impact-ionization model:
Accepted limitations:
```

Expected: the next worker can start from the generated report paths instead of rediscovering the run.

## Final Verification Commands

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build\test_impact_ionization.exe
ctest --test-dir build --output-on-failure -R "reference_tcad|sentaurus_import"
```

Then run the full suite if the promoted JSON or C++ physics implementation changed:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
ctest --test-dir build --output-on-failure
```

Expected: focused tests pass before promotion; full suite passes before merging any C++ solver or physics change.
