# pn2d Calibration Next Phase Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the current pn2d Sentaurus import loop from a strict Newton smoke gate into a reproducible calibration workflow for IV/BV quantity, trend, and node-level doping fidelity.

**Architecture:** Keep the default pn2d reference deck conservative and runnable. Add isolated candidate scans and gates before promoting any calibration into `reference_tcad/pn2d/pn2d_reference.json`. Separate three concerns: Newton ownership, physical-model calibration, and TDR doping import policy.

**Tech Stack:** C++20, CMake/Ninja, Catch2, Python regression tests, PowerShell scan scripts, existing `scripts/sentaurus_import.py` reference workflow.

---

## Current State Summary

- Branch: `codex-sentaurus-import-v1`, currently ahead of origin by local commits.
- Default pn2d decks use `solver.method: "gummel_newton"` with `handoff.fallback: "none"`.
- Current generated IV/BV rows show `handoff_stage: "newton"` and `newton_iterations > 0`; Gummel fallback is not part of the default pn2d gate.
- Current IV comparison:
  - Window: `0.2-0.3 V`
  - Points compared: `6`
  - Trend: matched
  - Quantity delta: `orders_of_magnitude ~= 0.5048`
  - Max relative error: `~0.687`
- Current BV comparison:
  - Window: `0.05 V`
  - Points compared: `1`
  - Default BV has `recombination: ["none"]` to match the BV command physics block.
  - Quantity delta: `orders_of_magnitude ~= 0.3508`
  - Vela/reference ratio: `~0.4458`
- BV quick6 candidate:
  - Best candidate: `q_mu0p89_a0p89`
  - `ratio_vs_ref ~= 0.8628`
  - `orders ~= 0.0641`
  - This is reproducible with `scripts/scan_pn2d_bv_ct_quick6.ps1`, but it has not been promoted into the default reference gate.
- Main technical risk:
  - BV can be tuned very close at `0.05 V`, but promotion is unsafe until IV impact and multi-bias BV behavior are checked.

---

## File Structure

- Modify: `docs/validation/pn2d_sentaurus_comparison.md`
  - Keep current evidence table updated after each calibration experiment.
- Modify: `reference_tcad/pn2d/pn2d_reference.json`
  - Only promote calibration after candidate isolation passes.
- Modify/Create: `scripts/scan_pn2d_bv_ct_quick6.ps1`
  - Reproducible BV 0.05 V CT mobility quick refinement.
- Create: `scripts/scan_pn2d_candidate_iv_bv.ps1`
  - Run one named candidate against both IV and BV generated decks.
- Modify: `tests/regression/test_reference_tcad_tools.py`
  - Add coverage that calibration helper scripts exist and encode named candidate cases.
- Modify: `tests/regression/test_sentaurus_sample_integration.py`
  - Strengthen pn2d integration expectations after a candidate is promoted.
- Modify: `scripts/sentaurus_import.py`
  - If needed, support named reference-config overlays so calibration candidates can be generated without mutating the base config.
- Modify: `src/io/SentaurusTdrReader.cpp`
  - Improve compensated node handling only after isolating doping-policy evidence.
- Modify: `tests/test_sentaurus_tdr_reader.cpp`
  - Add fixture-style tests for compensated donor/acceptor node policy behavior.

---

## Task 1: Freeze Current pn2d Baseline Evidence

**Files:**
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`
- Test: `tests/regression/test_reference_tcad_tools.py`

- [ ] **Step 1: Run current reference import**

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
if (Test-Path build\pn2d_current_review) {
  Remove-Item -LiteralPath build\pn2d_current_review -Recurse -Force
}
python scripts\sentaurus_import.py reference `
  --config reference_tcad\pn2d\pn2d_reference.json `
  --source-dir reference_tcad\pn2d `
  --output-dir build\pn2d_current_review `
  --tdr-importer build\sentaurus_import.exe `
  --runner build\vela_example_runner.exe
```

Expected:
- Command exits `0`.
- `reports/pn2d_iv_comparison.json` has `orders_of_magnitude` near `0.5048`.
- `reports/pn2d_bv_comparison.json` has `orders_of_magnitude` near `0.3508`.

- [ ] **Step 2: Run BV quick6 reproduction**

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\scan_pn2d_bv_ct_quick6.ps1 `
  -BaseConfig build\pn2d_current_review\vela\simulation_bv.json `
  -ReferenceCsv build\pn2d_current_review\reference_curves\pn2d_bv_reference.csv `
  -OutputSummary build\pn2d_current_review\vela\pn2d_bv_ct_quick6_summary.csv
```

Expected:
- `q_mu0p89_a0p89` remains the best row.
- `orders` is near `0.0641`.

- [ ] **Step 3: Update the validation doc**

Add a dated short note under the quick6 section:

```markdown
Fresh baseline rerun on 2026-05-26:

- default IV: `orders ~= 0.5048`, trend matched;
- default BV: `orders ~= 0.3508`;
- quick6 best `q_mu0p89_a0p89`: `orders ~= 0.0641`.
```

- [ ] **Step 4: Run targeted tests**

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_quick_refinement_script_is_documented_and_reproducible -v
ctest --test-dir build --output-on-failure -R ascii_sources
```

Expected: both pass.

- [ ] **Step 5: Commit**

```powershell
git add docs/validation/pn2d_sentaurus_comparison.md scripts/scan_pn2d_bv_ct_quick6.ps1 tests/regression/test_reference_tcad_tools.py
git commit -m "Document reproducible pn2d BV quick refinement"
```

---

## Task 2: Add a Candidate IV/BV Isolation Runner

**Files:**
- Create: `scripts/scan_pn2d_candidate_iv_bv.ps1`
- Modify: `tests/regression/test_reference_tcad_tools.py`
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`

- [ ] **Step 1: Write failing test for candidate runner presence**

Add this test to `tests/regression/test_reference_tcad_tools.py`:

```python
def test_pn2d_candidate_iv_bv_script_records_best_bv_candidate(self) -> None:
    script = REPO / "scripts" / "scan_pn2d_candidate_iv_bv.ps1"
    self.assertTrue(script.is_file(), f"missing pn2d candidate IV/BV script: {script}")
    text = script.read_text()
    self.assertIn("q_mu0p89_a0p89", text)
    self.assertIn("pn2d_candidate_iv_bv_summary.csv", text)
    self.assertIn("simulation_iv", text)
    self.assertIn("simulation_bv", text)
```

- [ ] **Step 2: Run and confirm failure**

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_candidate_iv_bv_script_records_best_bv_candidate -v
```

Expected: fail because the script does not exist.

- [ ] **Step 3: Create candidate runner**

Create `scripts/scan_pn2d_candidate_iv_bv.ps1` that:
- accepts `-BaseDir build\pn2d_current_review`;
- reads `vela/simulation_iv.json` and `vela/simulation_bv.json`;
- writes candidate configs into `vela/`;
- applies candidate mobility:

```powershell
$cfg.solver.bandgap_narrowing = "none"
$cfg.solver.mobility = @{
    model = "caughey_thomas"
    electron_mu_min_m2_V_s = 52.2 * 0.89
    hole_mu_min_m2_V_s = 44.9 * 0.89
    electron_nref_m3 = 9.68e16
    hole_nref_m3 = 2.23e17
    electron_alpha = 0.68 * 0.89
    hole_alpha = 0.70 * 0.89
}
```

For IV output, use:
- `pn2d_iv_q_mu0p89_a0p89.csv`

For BV output, use:
- `pn2d_bv_q_mu0p89_a0p89.csv`

Summary CSV columns:

```text
case,kind,status,points,orders,ratio_vs_ref,total_A_per_um,csv_file,config
```

- [ ] **Step 4: Run candidate script**

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\scan_pn2d_candidate_iv_bv.ps1 `
  -BaseDir build\pn2d_current_review `
  -Case q_mu0p89_a0p89
```

Expected:
- BV remains near `orders ~= 0.0641`.
- IV result is present and finite.
- No default `simulation_iv.json` or `simulation_bv.json` is modified.

- [ ] **Step 5: Document IV impact**

Add a table to `docs/validation/pn2d_sentaurus_comparison.md`:

```markdown
## Candidate Mobility Impact: q_mu0p89_a0p89

| Metric | Default | Candidate |
| --- | ---: | ---: |
| IV orders | ... | ... |
| BV 0.05 V orders | 0.3508 | 0.0641 |
```

- [ ] **Step 6: Run tests**

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_candidate_iv_bv_script_records_best_bv_candidate -v
ctest --test-dir build --output-on-failure -R ascii_sources
```

Expected: pass.

- [ ] **Step 7: Commit**

```powershell
git add scripts/scan_pn2d_candidate_iv_bv.ps1 tests/regression/test_reference_tcad_tools.py docs/validation/pn2d_sentaurus_comparison.md
git commit -m "Add pn2d IV BV candidate calibration scan"
```

---

## Task 3: Decide Whether to Promote BV Candidate Calibration

**Files:**
- Modify: `reference_tcad/pn2d/pn2d_reference.json`
- Modify: `tests/regression/test_sentaurus_sample_integration.py`
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`

- [ ] **Step 1: Define promotion rule**

Use this rule unless new evidence contradicts it:

```text
Promote q_mu0p89_a0p89 only if:
- default strict Newton provenance remains unchanged;
- IV trend still matches;
- IV orders remains <= current default + 0.15;
- BV 0.05 V orders improves below 0.20;
- ascii_sources passes.
```

- [ ] **Step 2: Write failing integration expectation**

If promoting, update `tests/regression/test_sentaurus_sample_integration.py` so BV expects the promoted mobility object:

```python
self.assertEqual(bv_deck["solver"]["bandgap_narrowing"], "none")
self.assertEqual(bv_deck["solver"]["mobility"]["model"], "caughey_thomas")
self.assertAlmostEqual(bv_deck["solver"]["mobility"]["electron_mu_min_m2_V_s"], 46.458)
self.assertAlmostEqual(bv_deck["solver"]["mobility"]["hole_mu_min_m2_V_s"], 39.961)
self.assertAlmostEqual(bv_deck["solver"]["mobility"]["electron_alpha"], 0.6052)
self.assertAlmostEqual(bv_deck["solver"]["mobility"]["hole_alpha"], 0.623)
self.assertLess(bv_report["iv"]["orders_of_magnitude"], 0.20)
```

- [ ] **Step 3: Run and confirm failure**

```powershell
python -m unittest tests.regression.test_sentaurus_sample_integration.SentaurusSampleIntegrationTest.test_pn2d_reference_import_when_enabled -v
```

Expected: fail until the reference config is promoted.

- [ ] **Step 4: Promote only BV override**

In `reference_tcad/pn2d/pn2d_reference.json`, under the BV simulation `vela_solver`, set:

```json
"mobility": {
  "model": "caughey_thomas",
  "electron_mu_min_m2_V_s": 46.458,
  "hole_mu_min_m2_V_s": 39.961,
  "electron_nref_m3": 9.68e16,
  "hole_nref_m3": 2.23e17,
  "electron_alpha": 0.6052,
  "hole_alpha": 0.623
},
"bandgap_narrowing": "none"
```

Keep:

```json
"recombination": ["none"],
"impact_ionization": { "model": "none" }
```

- [ ] **Step 5: Verify promoted config**

```powershell
if (Test-Path build\pn2d_promoted_candidate) {
  Remove-Item -LiteralPath build\pn2d_promoted_candidate -Recurse -Force
}
python scripts\sentaurus_import.py reference `
  --config reference_tcad\pn2d\pn2d_reference.json `
  --source-dir reference_tcad\pn2d `
  --output-dir build\pn2d_promoted_candidate `
  --tdr-importer build\sentaurus_import.exe `
  --runner build\vela_example_runner.exe
python -m unittest tests.regression.test_sentaurus_sample_integration.SentaurusSampleIntegrationTest.test_pn2d_reference_import_when_enabled -v
ctest --test-dir build --output-on-failure -R ascii_sources
```

Expected: all pass.

- [ ] **Step 6: Commit**

```powershell
git add reference_tcad/pn2d/pn2d_reference.json tests/regression/test_sentaurus_sample_integration.py docs/validation/pn2d_sentaurus_comparison.md
git commit -m "Promote pn2d BV mobility calibration"
```

If IV impact violates the promotion rule, do not promote. Instead, document the candidate as BV-only evidence and continue Task 4.

---

## Task 4: Improve TDR Compensated Node Doping Policy

**Files:**
- Modify: `src/io/SentaurusTdrReader.cpp`
- Modify: `tests/test_sentaurus_tdr_reader.cpp`
- Modify: `docs/config_schema.md`
- Modify: `reference_tcad/pn2d/pn2d_reference.json`

- [ ] **Step 1: Inspect current compensated-node metadata**

Run:

```powershell
Get-Content build\pn2d_current_review\doping_metadata.json
Import-Csv build\pn2d_current_review\doping.csv | Where-Object {
  ([double]$_.donors_cm3 -gt 0) -and ([double]$_.acceptors_cm3 -gt 0)
} | Select-Object -First 20
```

Expected:
- Understand how many compensated nodes are present.
- Identify whether junction nodes are being assigned by region dominance, raw donor/acceptor presence, or reported net doping.

- [ ] **Step 2: Write a focused reader test**

Add a test in `tests/test_sentaurus_tdr_reader.cpp` for a synthetic node where:
- donor and acceptor are both nonzero;
- node is shared by p and n regions;
- policy is `dominant_signed_region`.

Expected behavior:
- p-side node keeps acceptor-dominant net sign;
- n-side node keeps donor-dominant net sign;
- metadata records compensation count and policy.

- [ ] **Step 3: Run and confirm failure**

```powershell
ctest --test-dir build --output-on-failure -R sentaurus_tdr_reader
```

Expected: fail until policy is improved.

- [ ] **Step 4: Implement node-level junction policy**

In `src/io/SentaurusTdrReader.cpp`, update compensated node handling so:
- exact donor/acceptor values remain exported;
- solver-facing signed compensation uses adjacent semiconductor region dominance when available;
- ambiguous ties are recorded in metadata instead of silently picking a side.

- [ ] **Step 5: Document policy**

In `docs/config_schema.md`, document:

```text
tdr_doping.compensated_node_policy:
- reported
- dominant_signed_region
```

Avoid vendor-specific terminology beyond neutral field names.

- [ ] **Step 6: Verify**

```powershell
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R "sentaurus_tdr_reader|sentaurus_sample|ascii_sources"
```

Expected: all pass.

- [ ] **Step 7: Commit**

```powershell
git add src/io/SentaurusTdrReader.cpp tests/test_sentaurus_tdr_reader.cpp docs/config_schema.md reference_tcad/pn2d/pn2d_reference.json
git commit -m "Improve compensated TDR doping node policy"
```

---

## Task 5: Add More Reference Device Configs

**Files:**
- Modify: `reference_tcad/README.md`
- Create or modify: `reference_tcad/igbt2d/*`
- Create or modify: `reference_tcad/ldmos2d/*`
- Modify: `tests/regression/test_reference_tcad_tools.py`

- [ ] **Step 1: Inventory available sample folders**

Run:

```powershell
Get-ChildItem reference_tcad -Directory
Get-ChildItem reference_tcad -Recurse -Filter '*reference.json'
```

Expected:
- List which devices already have reference configs.
- Identify missing IV/BV config coverage.

- [ ] **Step 2: Add regression expectation for config presence**

In `tests/regression/test_reference_tcad_tools.py`, add:

```python
def test_reference_tcad_device_configs_exist(self) -> None:
    expected = [
        REPO / "reference_tcad" / "pn2d" / "pn2d_reference.json",
        REPO / "reference_tcad" / "ldmos2d" / "ldmos2d_reference.json",
        REPO / "reference_tcad" / "igbt2d" / "igbt2d_reference.json",
    ]
    for path in expected:
        self.assertTrue(path.is_file(), f"missing reference config: {path}")
```

- [ ] **Step 3: Run and confirm failure**

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_reference_tcad_device_configs_exist -v
```

Expected: fail for missing configs.

- [ ] **Step 4: Add minimal config files**

For each new device config:
- include case/device names;
- point to available mesh/TDR/PLT/CMD artifacts;
- mark unsupported physics in docs;
- use conservative comparison gates first.

- [ ] **Step 5: Verify config discovery**

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_reference_tcad_device_configs_exist -v
ctest --test-dir build --output-on-failure -R ascii_sources
```

Expected: pass.

- [ ] **Step 6: Commit**

```powershell
git add reference_tcad/README.md reference_tcad/ldmos2d reference_tcad/igbt2d tests/regression/test_reference_tcad_tools.py
git commit -m "Add additional reference TCAD configs"
```

---

## Final Verification

Run after all promoted changes:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Expected:
- Full suite passes.
- pn2d sample integration passes when local Sentaurus sample artifacts are available.
- `ascii_sources` passes.

If the full suite is too slow during an intermediate checkpoint, run:

```powershell
ctest --test-dir build --output-on-failure -R "sentaurus|dc_sweep|newton|ascii_sources"
```

and explicitly record that full suite is still pending.
