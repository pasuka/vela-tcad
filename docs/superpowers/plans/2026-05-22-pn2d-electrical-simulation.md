# pn2d Electrical Simulation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make Vela TCAD run meaningful pn2d IV/BV electrical simulations from the committed Sentaurus reference artifacts and produce defensible comparison reports.

**Architecture:** Keep Sentaurus parsing in `scripts/sentaurus_import.py` and `src/io/SentaurusTdrReader.cpp`, convert neutral exports through `scripts/convert_tcad_export.py`, and extend the Vela deck/runtime path only where generic TCAD functionality is missing. The first milestone is correctness of the pn2d physical input, especially doping; later milestones improve physical equivalence with Sentaurus models.

**Tech Stack:** C++20, CMake/Ninja, HDF5 reader, nlohmann-json, Python standard library regression tools, Vela `dc_sweep` runner.

---

## File Structure

- Modify `src/io/SentaurusTdrReader.cpp`: map Sentaurus active dopant field aliases into neutral `doping.csv`.
- Modify `tests/test_sentaurus_tdr_reader.cpp`: synthetic HDF5 coverage for `BoronActiveConcentration` and `PhosphorusActiveConcentration`.
- Modify `scripts/convert_tcad_export.py`: preserve `doping.csv` as a node-level input file while retaining region averages for legacy decks.
- Modify `src/simulation/DCSweep.cpp`: resolve and load optional node doping files before falling back to region doping.
- Modify `src/simulation/ConfigParsing.cpp` and `include/vela/simulation/ConfigParsing.h`: expose reusable node-doping parsing if a header declaration exists; otherwise keep the helper local in `DCSweep.cpp`.
- Modify `scripts/sentaurus_import.py`: map SDevice physics models to Vela solver keys and classify approximations.
- Modify `tests/regression/test_reference_tcad_tools.py`: conversion tests for node doping and deck schema.
- Modify `tests/regression/test_sentaurus_import_tools.py`: reference config tests for solver physics mapping.
- Modify `tests/regression/test_sentaurus_sample_integration.py`: pn2d committed-artifact smoke test without requiring an environment variable.
- Optionally add `reference_tcad/pn2d/pn2d_reference.json`: reusable config for pn2d IV/BV import.

---

### Task 1: Fix pn2d Active Dopant Export

**Files:**
- Modify: `src/io/SentaurusTdrReader.cpp`
- Modify: `tests/test_sentaurus_tdr_reader.cpp`
- Test: `tests/test_sentaurus_tdr_reader.cpp`
- Test: `tests/regression/test_sentaurus_sample_integration.py`

- [ ] **Step 1: Write the failing synthetic TDR test**

In `tests/test_sentaurus_tdr_reader.cpp`, change the synthetic donor/acceptor dataset names from only `DonorConcentration` / `AcceptorConcentration` coverage to include Sentaurus active dopants:

```cpp
writeScalarField(state, "dataset_3", "PhosphorusActiveConcentration", 0, "cm^-3",
                 {1.0e17, 1.0e17, 1.0e17});
writeScalarField(state, "dataset_4", "BoronActiveConcentration", 1, "cm^-3",
                 {2.0e17, 2.0e17, 2.0e17});
```

Add assertions after neutral export:

```cpp
const auto doping = readCsvRows(outDir / "doping.csv");
REQUIRE(doping.at(1).at("donors_cm3") == "1e+17");
REQUIRE(doping.at(3).at("acceptors_cm3") == "2e+17");
```

- [ ] **Step 2: Verify the test fails**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build --target test_sentaurus_tdr_reader --parallel
ctest --test-dir build --output-on-failure -R sentaurus_tdr
```

Expected: FAIL because `doping.csv` still contains zero donors/acceptors for active dopant fields.

- [ ] **Step 3: Implement active dopant aliases**

In `src/io/SentaurusTdrReader.cpp`, replace the current name check around the doping export loop with:

```cpp
const auto donorLike = [](const std::string& name) {
    return name == "DonorConcentration" ||
           name == "PhosphorusActiveConcentration" ||
           name == "ArsenicActiveConcentration" ||
           name == "AntimonyActiveConcentration";
};
const auto acceptorLike = [](const std::string& name) {
    return name == "AcceptorConcentration" ||
           name == "BoronActiveConcentration" ||
           name == "AluminumActiveConcentration" ||
           name == "IndiumActiveConcentration";
};
if (!donorLike(field.name) && !acceptorLike(field.name)) {
    continue;
}
```

Use `donorLike(field.name)` in the assignment branch.

- [ ] **Step 4: Verify synthetic and real pn2d export**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build --target sentaurus_import test_sentaurus_tdr_reader --parallel
ctest --test-dir build --output-on-failure -R sentaurus_tdr
python tests\regression\test_sentaurus_sample_integration.py -k pn2d
```

Expected: PASS, and generated pn2d `doping.csv` has nonzero donor and acceptor populations.

- [ ] **Step 5: Commit**

```powershell
git add src/io/SentaurusTdrReader.cpp tests/test_sentaurus_tdr_reader.cpp tests/regression/test_sentaurus_sample_integration.py
git commit -m "Map Sentaurus active dopants into neutral exports"
```

---

### Task 2: Preserve Node-Level Doping in Vela Decks

**Files:**
- Modify: `scripts/convert_tcad_export.py`
- Modify: `src/simulation/DCSweep.cpp`
- Test: `tests/regression/test_reference_tcad_tools.py`
- Test: `tests/test_dc_sweep.cpp`

- [ ] **Step 1: Write failing converter test**

In `tests/regression/test_reference_tcad_tools.py`, extend the PN export conversion test to assert:

```python
self.assertEqual(iv["node_doping_file"], "doping.csv")
self.assertEqual(bv["node_doping_file"], "doping.csv")
self.assertEqual(iv["doping"][0]["region"], "p_region")
```

Expected behavior: deck keeps region averages for compatibility and also points the runtime to `doping.csv`.

- [ ] **Step 2: Write failing runtime test**

In `tests/test_dc_sweep.cpp`, add a small deck test that writes `doping.csv` beside the config:

```cpp
TEST_CASE("DCSweep reads node_doping_file before region averages", "[dc_sweep][doping]")
{
    // Create a two-region PN mesh, write region doping as zero, and write
    // node_doping_file with nonzero donor/acceptor values. Run a zero-bias
    // one-point sweep and assert it converges with finite current columns.
}
```

Use the existing JSON/temp-file helpers in `tests/test_dc_sweep.cpp`; keep the sweep `start=0`, `stop=0`, `step=0.1`.

- [ ] **Step 3: Verify tests fail**

Run:

```powershell
python tests\regression\test_reference_tcad_tools.py -k pn_export
ctest --test-dir build --output-on-failure -R dc_sweep
```

Expected: converter test fails because `node_doping_file` is absent; C++ test fails because runtime ignores the file.

- [ ] **Step 4: Add deck schema field in converter**

In `scripts/convert_tcad_export.py`, add to `base_deck(...)`:

```python
"node_doping_file": "doping.csv",
```

Keep the existing `doping` array unchanged.

- [ ] **Step 5: Add runtime node doping loader**

In `src/simulation/DCSweep.cpp`, change `dopingFromJson(...)` to receive `cfgDir`:

```cpp
DopingModel dopingFromJson(const DeviceMesh& mesh,
                           const nlohmann::json& cfg,
                           const std::filesystem::path& cfgDir,
                           UnitScalingConfig scaling)
```

If `cfg.contains("node_doping_file")`, resolve the path relative to `cfgDir`, read CSV columns `node_id,donors_cm3,acceptors_cm3`, convert concentrations with `scaling.concentrationToSI(...)`, and call `model.setNodeDoping(...)`. If the file has an invalid node id, throw `DCSweep: node_doping_file references missing node id N`.

- [ ] **Step 6: Verify**

Run:

```powershell
cmake --build build --parallel
python tests\regression\test_reference_tcad_tools.py -k pn_export
ctest --test-dir build --output-on-failure -R dc_sweep
```

Expected: all selected tests pass.

- [ ] **Step 7: Commit**

```powershell
git add scripts/convert_tcad_export.py src/simulation/DCSweep.cpp tests/regression/test_reference_tcad_tools.py tests/test_dc_sweep.cpp
git commit -m "Use node-level doping in generated TCAD decks"
```

---

### Task 3: Map pn2d Sentaurus Physics to Vela Solver Settings

**Files:**
- Modify: `scripts/sentaurus_import.py`
- Test: `tests/regression/test_sentaurus_import_tools.py`

- [ ] **Step 1: Write failing reference import assertions**

In `test_reference_import_config_generates_iv_bv_tree_and_reports`, assert:

```python
self.assertEqual(iv_deck["solver"]["mobility"]["model"], "caughey_thomas_field")
self.assertEqual(iv_deck["solver"]["recombination"], ["srh", "auger"])
self.assertEqual(iv_deck["solver"]["bandgap_narrowing"], "slotboom")
self.assertNotIn("impact_ionization", iv_deck["solver"])
self.assertEqual(bv_deck["solver"]["impact_ionization"]["model"], "selberherr")
self.assertIn("OkutoCrowell approximated by Selberherr", manifest["warnings"])
```

- [ ] **Step 2: Verify failure**

Run:

```powershell
python tests\regression\test_sentaurus_import_tools.py -k reference_import
```

Expected: FAIL because generated decks currently keep only generic Gummel settings.

- [ ] **Step 3: Implement physics mapper**

In `scripts/sentaurus_import.py`, add:

```python
def sentaurus_models(cmd_summary: dict[str, Any]) -> set[str]:
    models: set[str] = set()
    for physics in cmd_summary.get("physics", []):
        models.update(str(model) for model in physics.get("models", []))
    return models


def apply_solver_physics(deck: dict[str, Any],
                         cmd_summary: dict[str, Any],
                         sim: dict[str, Any]) -> list[str]:
    models = sentaurus_models(cmd_summary)
    solver = deck.setdefault("solver", {})
    warnings: list[str] = []
    if {"Mobility", "DopingDep"} & models:
        solver["mobility"] = {"model": "caughey_thomas"}
    if "HighFieldSaturation" in models:
        solver["mobility"] = {"model": "caughey_thomas_field"}
    recombination = []
    if "SRH" in models:
        recombination.append("srh")
    if "Auger" in models:
        recombination.append("auger")
    if recombination:
        solver["recombination"] = recombination
    if "OldSlotboom" in models:
        solver["bandgap_narrowing"] = "slotboom"
    if "Avalanche" in models:
        solver["impact_ionization"] = {"model": "selberherr"}
        if "OkutoCrowell" in models:
            warnings.append("OkutoCrowell approximated by Selberherr")
    if "Fermi" in models:
        warnings.append("Fermi statistics approximated by Boltzmann carrier statistics")
    return warnings
```

Call this from `patch_reference_deck(...)` and include returned warnings in the manifest warnings.

- [ ] **Step 4: Verify**

Run:

```powershell
python tests\regression\test_sentaurus_import_tools.py -k reference_import
```

Expected: PASS.

- [ ] **Step 5: Commit**

```powershell
git add scripts/sentaurus_import.py tests/regression/test_sentaurus_import_tools.py
git commit -m "Map Sentaurus pn2d physics into Vela decks"
```

---

### Task 4: Add Committed pn2d Reference Config and Smoke Gate

**Files:**
- Create: `reference_tcad/pn2d/pn2d_reference.json`
- Modify: `tests/regression/test_sentaurus_sample_integration.py`
- Test: `tests/regression/test_sentaurus_sample_integration.py`

- [ ] **Step 1: Add reusable config**

Create `reference_tcad/pn2d/pn2d_reference.json`:

```json
{
  "case": "pn2d",
  "device": "pn_diode",
  "mesh_tdr": "pn2d_msh.tdr",
  "sde_cmd": "pn2d_sde.cmd",
  "simulations": [
    {
      "name": "iv",
      "kind": "iv",
      "tdr": "pn2d_des.tdr",
      "cmd": "pn2d_sdevice.cmd",
      "plt": "pn2d_iv.plt",
      "bias_column": "Anode OuterVoltage",
      "current_column": "Anode TotalCurrent"
    },
    {
      "name": "bv",
      "kind": "bv",
      "tdr": "pn2d_bv_des.tdr",
      "cmd": "pn2d_bv_sdevice.cmd",
      "plt": "pn2d_bv.plt",
      "bias_column": "Cathode OuterVoltage",
      "current_column": "Cathode TotalCurrent"
    }
  ]
}
```

- [ ] **Step 2: Update smoke test to use committed config**

In `tests/regression/test_sentaurus_sample_integration.py`, change the pn2d test to use `REPO / "reference_tcad" / "pn2d"` when `VELA_SENTAURUS_PN2D_DIR` is not set, and load `pn2d_reference.json` instead of creating inline config.

- [ ] **Step 3: Add assertions for nonzero generated deck doping**

After import, assert:

```python
self.assertTrue(any(float(row["donors_cm3"]) > 0.0 for row in self._read_curve(out / "reference" / "doping.csv")))
self.assertTrue(any(float(row["acceptors_cm3"]) > 0.0 for row in self._read_curve(out / "reference" / "doping.csv")))
self.assertEqual(iv_deck["node_doping_file"], "doping.csv")
```

- [ ] **Step 4: Verify**

Run:

```powershell
python tests\regression\test_sentaurus_sample_integration.py -k pn2d
ctest --test-dir build --output-on-failure -R sentaurus_sample_integration
```

Expected: PASS without requiring `VELA_SENTAURUS_PN2D_DIR`.

- [ ] **Step 5: Commit**

```powershell
git add reference_tcad/pn2d/pn2d_reference.json tests/regression/test_sentaurus_sample_integration.py
git commit -m "Add committed pn2d Sentaurus smoke gate"
```

---

### Task 5: Run and Compare pn2d IV

**Files:**
- Modify: `scripts/sentaurus_import.py`
- Modify: `tests/regression/test_sentaurus_import_tools.py`

- [ ] **Step 1: Add IV execution test with fake runner plus generated solver fields**

Extend the fake runner test to verify the manifest includes:

```python
self.assertIn("reports/pn2d_iv_comparison.json", manifest["comparison_reports"])
self.assertIn("Fermi statistics approximated by Boltzmann carrier statistics", manifest["warnings"])
```

- [ ] **Step 2: Run real pn2d IV manually**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\sentaurus_import.py reference --config reference_tcad\pn2d\pn2d_reference.json --source-dir reference_tcad\pn2d --output-dir build\reference_tcad\pn2d --tdr-importer build\sentaurus_import.exe --runner build\vela_example_runner.exe
```

Expected: `build/reference_tcad/pn2d/vela/pn2d_iv.csv` exists and has finite values.

- [ ] **Step 3: If IV fails to converge, tune only generated deck settings**

Change `patch_reference_deck(...)` default solver for pn2d generated decks:

```python
solver.setdefault("max_iter", 150)
solver.setdefault("reltol", 1.0e-6)
solver.setdefault("damping_psi", 0.35)
```

Do not change core solver behavior in this task.

- [ ] **Step 4: Verify IV report**

Run:

```powershell
python scripts\compare_reference_curves.py --reference build\reference_tcad\pn2d\reference_curves\pn2d_iv_reference.csv --candidate build\reference_tcad\pn2d\vela\pn2d_iv.csv --output-json build\reference_tcad\pn2d\reports\pn2d_iv_comparison.json --output-md build\reference_tcad\pn2d\reports\pn2d_iv_comparison.md --kind iv --require-trend-match --min-points 10
```

Expected: PASS trend gate. Absolute error may remain large.

- [ ] **Step 5: Commit**

```powershell
git add scripts/sentaurus_import.py tests/regression/test_sentaurus_import_tools.py
git commit -m "Enable pn2d IV execution comparison"
```

---

### Task 6: Add BV Approximation Path

**Files:**
- Modify: `scripts/sentaurus_import.py`
- Modify: `scripts/compare_reference_curves.py` only if a non-current BV diagnostic comparison is needed.
- Test: `tests/regression/test_sentaurus_import_tools.py`

- [ ] **Step 1: Ensure BV deck uses approximate avalanche**

Assert in tests:

```python
self.assertEqual(bv_deck["solver"]["impact_ionization"]["model"], "selberherr")
self.assertEqual(bv_deck["sweep"]["mode"], "bv_reverse")
self.assertEqual(bv_deck["sweep"]["contact"], "Cathode")
```

- [ ] **Step 2: Keep BV comparison non-strict**

Set pn2d BV simulation config to include:

```json
"require_trend_match": false,
"comparison_kind": "iv"
```

This records magnitude and trend metadata without failing because Sentaurus uses Okuto-Crowell while Vela uses Selberherr.

- [ ] **Step 3: Run real BV manually**

Run:

```powershell
python scripts\sentaurus_import.py reference --config reference_tcad\pn2d\pn2d_reference.json --source-dir reference_tcad\pn2d --output-dir build\reference_tcad\pn2d --tdr-importer build\sentaurus_import.exe --runner build\vela_example_runner.exe
```

Expected: `pn2d_bv.csv` exists, all numeric values are finite, and comparison JSON is produced.

- [ ] **Step 4: Commit**

```powershell
git add scripts/sentaurus_import.py tests/regression/test_sentaurus_import_tools.py reference_tcad/pn2d/pn2d_reference.json
git commit -m "Add approximate pn2d BV comparison path"
```

---

### Task 7: Full Regression and Documentation

**Files:**
- Create: `docs/validation/pn2d_sentaurus_comparison.md`
- Modify: `docs/config_schema.md` if `node_doping_file` is introduced.

- [ ] **Step 1: Document current equivalence level**

Create `docs/validation/pn2d_sentaurus_comparison.md` with:

```markdown
# pn2d Sentaurus Comparison

The pn2d case is a 2-D abrupt PN diode with L=2.0 um, H=0.5 um,
Xj=1.0 um, Na=Nd=1e17 cm^-3, Anode on the left boundary, and Cathode
on the right boundary.

Vela currently treats Sentaurus Fermi statistics as Boltzmann carrier
statistics. Sentaurus Okuto-Crowell avalanche is approximated by Vela
Selberherr impact ionization for BV diagnostics. The IV comparison is
trend-gated; the BV comparison is diagnostic-only until Okuto-Crowell is
implemented or calibrated.
```

- [ ] **Step 2: Document `node_doping_file` schema**

Add to `docs/config_schema.md`:

```markdown
`node_doping_file` optionally points to a CSV with columns
`node_id,donors_cm3,acceptors_cm3`. When present, it overrides region-average
`doping` entries for drift-diffusion sweeps. With `scaling.mode =
unit_scaling`, concentrations use the same input concentration convention as
the rest of the deck.
```

- [ ] **Step 3: Run full verification**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Expected: all tests pass.

- [ ] **Step 4: Commit**

```powershell
git add docs/validation/pn2d_sentaurus_comparison.md docs/config_schema.md
git commit -m "Document pn2d Sentaurus comparison status"
```

---

## Self-Review

- Spec coverage: the plan covers active dopant import, node-level doping, physics mapping, pn2d config, IV/BV execution, comparison reports, and documentation.
- Placeholder scan: no deferred placeholders are used; approximation limits are explicit.
- Type consistency: `node_doping_file`, `solver.mobility`, `solver.recombination`, `solver.bandgap_narrowing`, and `solver.impact_ionization` match existing deck parsing conventions or are introduced with tests.

