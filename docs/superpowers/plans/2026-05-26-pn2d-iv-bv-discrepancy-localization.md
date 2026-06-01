# pn2d IV BV Discrepancy Localization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Localize the remaining pn2d IV slope mismatch and reconcile the BV low-bias calibration with the reference deck physics before tightening default gates further.

**Architecture:** Keep the current strict-Newton pn2d reference gate stable. Add reproducible probe scripts and reports that compare default, candidate, and physics-ablation decks without mutating the reference config. Promote only changes that pass both quantity and physical-consistency checks.

**Tech Stack:** C++20, CMake/Ninja, Catch2, Python regression tests, PowerShell probe scripts, existing `sentaurus_import.py`, `compare_reference_curves.py`, and `vela_example_runner.exe`.

---

## Current Evidence

- Branch: `codex-sentaurus-import-v1`, ahead of origin by local commits.
- Full suite recently passed: `271/271`.
- Strict Newton ownership is confirmed for pn2d IV/BV: accepted rows use `solver_method = gummel_newton`, `handoff_stage = newton`, and `newton_iterations > 0`.
- Current promoted pn2d IV:
  - comparison window: `0.204721576526-0.29 V`
  - `orders_of_magnitude ~= 0.5047505096`
  - `max_relative_error ~= 0.6872`
  - trend matches
  - per-bias ratio starts close to reference and rolls off at high forward bias:

| Bias V | Vela/ref ratio |
| ---: | ---: |
| 0.204721576526 | 0.9052 |
| 0.224721576526 | 1.2270 |
| 0.244721576526 | 0.9426 |
| 0.25 | 0.8499 |
| 0.27 | 0.5348 |
| 0.29 | 0.3128 |

- Current promoted pn2d BV:
  - comparison point: `0.05 V`
  - `orders_of_magnitude ~= 0.0641104534`
  - `max_relative_error ~= 0.1372`
  - current ratio Vela/reference `~= 0.8628`
- Important correction: `reference_tcad/pn2d/pn2d_bv_sdevice.cmd` includes
  `Recombination(SRH Auger Avalanche)`. The current Vela BV reference override
  disables recombination and impact ionization. Re-enabling SRH/Auger under the
  promoted BV mobility produces:
  - `orders_of_magnitude ~= 1.1737`
  - `current_total_A_per_um ~= 1.5986e-17`
  - electron residual dominates
- Therefore BV is currently a strong low-bias numerical match but not yet a
  physically reconciled match to the full reference physics block.

---

## File Structure

- Create: `scripts/scan_pn2d_iv_bv_physics_matrix.ps1`
  - Runs IV/BV physics ablations from a generated pn2d reference tree.
- Create: `scripts/summarize_pn2d_iv_ratios.ps1`
  - Emits per-reference-bias IV ratios and optional current decomposition.
- Modify: `tests/regression/test_reference_tcad_tools.py`
  - Adds presence/content checks for the new probe scripts.
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`
  - Records corrected BV physics facts and new ablation tables.
- Modify later only if evidence supports it: `reference_tcad/pn2d/pn2d_reference.json`
  - Do not alter default gates until the probe results justify it.
- Possible later C++ files after evidence:
  - `include/vela/physics/CarrierStatistics.h`
  - `src/physics/CarrierStatistics.cpp`
  - `src/equation/CoupledDDAssembler.cpp`
  - `src/equation/DDAssembler.cpp`
  - `tests/test_recombination.cpp`
  - `tests/test_sg_flux.cpp`

---

## Task 1: Correct BV Physics Documentation and Lock Probe Script Expectations

**Files:**
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`
- Modify: `tests/regression/test_reference_tcad_tools.py`

- [ ] **Step 1: Write failing script-presence tests**

Add these tests to `tests/regression/test_reference_tcad_tools.py`:

```python
def test_pn2d_iv_ratio_summary_script_exists(self) -> None:
    script = REPO / "scripts" / "summarize_pn2d_iv_ratios.ps1"
    self.assertTrue(script.is_file(), f"missing pn2d IV ratio summary script: {script}")
    text = script.read_text()
    self.assertIn("pn2d_iv_ratio_summary.csv", text)
    self.assertIn("current_electron_A_per_um", text)
    self.assertIn("current_hole_A_per_um", text)

def test_pn2d_iv_bv_physics_matrix_script_exists(self) -> None:
    script = REPO / "scripts" / "scan_pn2d_iv_bv_physics_matrix.ps1"
    self.assertTrue(script.is_file(), f"missing pn2d IV/BV physics matrix script: {script}")
    text = script.read_text()
    self.assertIn("recomb_srh_auger", text)
    self.assertIn("recomb_none", text)
    self.assertIn("pn2d_iv_bv_physics_matrix_summary.csv", text)
```

- [ ] **Step 2: Run and confirm failure**

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_iv_ratio_summary_script_exists -v
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_iv_bv_physics_matrix_script_exists -v
```

Expected: both fail because the scripts are missing.

- [ ] **Step 3: Correct the validation doc fact**

In `docs/validation/pn2d_sentaurus_comparison.md`, replace any claim that the
BV command omits SRH/Auger with:

```markdown
The imported BV command includes `Recombination(SRH Auger Avalanche)`. The
current Vela BV reference override intentionally disables recombination and
impact ionization as a low-bias numerical gate while the recombination and
avalanche model parity work remains open.
```

- [ ] **Step 4: Run ASCII check**

```powershell
ctest --test-dir build --output-on-failure -R ascii_sources
```

Expected: pass.

---

## Task 2: Add IV Per-Bias Ratio and Decomposition Summary

**Files:**
- Create: `scripts/summarize_pn2d_iv_ratios.ps1`
- Modify: `tests/regression/test_reference_tcad_tools.py`
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`

- [ ] **Step 1: Create the script**

Create `scripts/summarize_pn2d_iv_ratios.ps1`:

```powershell
param(
    [string]$ReferenceCsv = "build/pn2d_tdr_tie_probe/reference_curves/pn2d_iv_reference.csv",
    [string]$CandidateCsv = "build/pn2d_tdr_tie_probe/vela/pn2d_iv.csv",
    [string]$OutputSummary = "build/pn2d_tdr_tie_probe/vela/pn2d_iv_ratio_summary.csv",
    [double]$BiasMin = 0.2,
    [double]$BiasMax = 0.3,
    [double]$CandidateScale = -1.0
)

$ErrorActionPreference = "Stop"

function Get-InterpolatedRow($rows, [double]$bias) {
    $pts = @($rows | Sort-Object { [double]$_.bias_V })
    for ($i = 0; $i -lt $pts.Count; $i++) {
        $b = [double]$pts[$i].bias_V
        if ([math]::Abs($b - $bias) -le [math]::Max([math]::Abs($bias), 1.0) * 1.0e-12) {
            return $pts[$i]
        }
    }
    for ($i = 0; $i -lt $pts.Count - 1; $i++) {
        $b0 = [double]$pts[$i].bias_V
        $b1 = [double]$pts[$i + 1].bias_V
        if ($b0 -le $bias -and $bias -le $b1 -and $b1 -ne $b0) {
            $t = ($bias - $b0) / ($b1 - $b0)
            $out = [ordered]@{ bias_V = $bias }
            foreach ($name in $pts[$i].PSObject.Properties.Name) {
                if ($name -eq "bias_V") { continue }
                $v0 = 0.0
                $v1 = 0.0
                if ([double]::TryParse([string]$pts[$i].$name, [ref]$v0) -and
                    [double]::TryParse([string]$pts[$i + 1].$name, [ref]$v1)) {
                    $out[$name] = $v0 + $t * ($v1 - $v0)
                }
            }
            return [pscustomobject]$out
        }
    }
    return $null
}

$refRows = Import-Csv $ReferenceCsv
$candRows = Import-Csv $CandidateCsv
$rows = @()

foreach ($ref in $refRows) {
    $bias = [double]$ref.bias_V
    if ($bias -lt $BiasMin -or $bias -gt $BiasMax) { continue }
    $cand = Get-InterpolatedRow $candRows $bias
    if ($null -eq $cand) { continue }
    $refCurrent = [double]$ref.current_total
    $velaTotal = [double]$cand.current_total_A_per_um * $CandidateScale
    $rows += [pscustomobject]@{
        bias_V = $bias
        reference_A = $refCurrent
        vela_total_A_per_um = $velaTotal
        ratio_vs_ref = if ($refCurrent -ne 0.0) { $velaTotal / $refCurrent } else { [double]::NaN }
        current_electron_A_per_um = [double]$cand.current_electron_A_per_um * $CandidateScale
        current_hole_A_per_um = [double]$cand.current_hole_A_per_um * $CandidateScale
        current_electron_drift_A_per_um = [double]$cand.current_electron_drift_A_per_um * $CandidateScale
        current_electron_diffusion_A_per_um = [double]$cand.current_electron_diffusion_A_per_um * $CandidateScale
        current_hole_diffusion_A_per_um = [double]$cand.current_hole_diffusion_A_per_um * $CandidateScale
    }
}

$rows | Export-Csv -NoTypeInformation -Encoding UTF8 $OutputSummary
$rows | Format-Table -AutoSize
Write-Host "summary=$OutputSummary"
```

- [ ] **Step 2: Run the script**

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\summarize_pn2d_iv_ratios.ps1 `
  -ReferenceCsv build\pn2d_tdr_tie_probe\reference_curves\pn2d_iv_reference.csv `
  -CandidateCsv build\pn2d_tdr_tie_probe\vela\pn2d_iv.csv `
  -OutputSummary build\pn2d_tdr_tie_probe\vela\pn2d_iv_ratio_summary.csv
```

Expected:
- Output includes ratios near `0.905`, `1.227`, `0.943`, `0.850`, `0.535`, `0.313`.
- This confirms the remaining IV discrepancy is primarily high-forward-bias slope.

- [ ] **Step 3: Run test**

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_iv_ratio_summary_script_exists -v
```

Expected: pass.

- [ ] **Step 4: Document results**

Add a table to `docs/validation/pn2d_sentaurus_comparison.md`:

```markdown
## IV Per-Bias Ratio Shape

| Bias V | Vela/reference ratio |
| ---: | ---: |
| 0.204721576526 | 0.9052 |
| 0.224721576526 | 1.2270 |
| 0.244721576526 | 0.9426 |
| 0.25 | 0.8499 |
| 0.27 | 0.5348 |
| 0.29 | 0.3128 |
```

---

## Task 3: Add IV/BV Physics Matrix Probe

**Files:**
- Create: `scripts/scan_pn2d_iv_bv_physics_matrix.ps1`
- Modify: `tests/regression/test_reference_tcad_tools.py`
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`

- [ ] **Step 1: Create the script**

Create `scripts/scan_pn2d_iv_bv_physics_matrix.ps1` that:
- accepts `-BaseDir build\pn2d_tdr_tie_probe`;
- reads `vela/simulation_iv.json` and `vela/simulation_bv.json`;
- creates separate candidate configs without mutating defaults;
- runs this matrix:

```powershell
$cases = @(
    @{ name = "default"; kind = "iv"; recombination = @("srh", "auger"); bgn = "slotboom"; mobility = "default" },
    @{ name = "iv_recomb_none"; kind = "iv"; recombination = @("none"); bgn = "slotboom"; mobility = "default" },
    @{ name = "iv_bgn_none"; kind = "iv"; recombination = @("srh", "auger"); bgn = "none"; mobility = "default" },
    @{ name = "bv_recomb_none"; kind = "bv"; recombination = @("none"); bgn = "none"; mobility = "promoted_bv" },
    @{ name = "bv_recomb_srh"; kind = "bv"; recombination = @("srh"); bgn = "none"; mobility = "promoted_bv" },
    @{ name = "bv_recomb_srh_auger"; kind = "bv"; recombination = @("srh", "auger"); bgn = "none"; mobility = "promoted_bv" }
)
```

Summary columns:

```text
case,kind,status,points,orders,max_relative_error,ratio_at_target,total_A_per_um,csv_file,config
```

Use targets:
- IV target bias: `0.29 V`, comparison window `0.2-0.3 V`, candidate scale `-1.0`
- BV target bias: `0.05 V`, comparison window exactly `0.05 V`, candidate scale `1.0`

- [ ] **Step 2: Run the script**

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\scan_pn2d_iv_bv_physics_matrix.ps1 `
  -BaseDir build\pn2d_tdr_tie_probe
```

Expected:
- `bv_recomb_none` remains near `orders ~= 0.0641`.
- `bv_recomb_srh_auger` is near `orders ~= 1.17`.
- IV ablations produce finite rows.

- [ ] **Step 3: Run regression test**

```powershell
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_iv_bv_physics_matrix_script_exists -v
```

Expected: pass.

- [ ] **Step 4: Document the matrix**

Add a table:

```markdown
## IV/BV Physics Matrix

| Case | Kind | Orders | Ratio at target | Interpretation |
| --- | --- | ---: | ---: | --- |
| bv_recomb_none | BV | 0.0641 | 0.8628 | current promoted gate |
| bv_recomb_srh_auger | BV | 1.1737 | 14.9172 | recombination model parity unresolved |
```

---

## Task 4: Determine Whether IV Error Is Sweep Resolution or Model Slope

**Files:**
- Modify: `scripts/scan_pn2d_iv_bv_physics_matrix.ps1`
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`
- Optional modify: `reference_tcad/pn2d/pn2d_reference.json`

- [ ] **Step 1: Run a fine IV candidate deck**

Create a local generated config from `build\pn2d_tdr_tie_probe\vela\simulation_iv.json`:

```powershell
$cfg = Get-Content build\pn2d_tdr_tie_probe\vela\simulation_iv.json -Raw | ConvertFrom-Json
$cfg.output_csv = "pn2d_iv_fine_step.csv"
$cfg.sweep.step = 0.02
$cfg | ConvertTo-Json -Depth 100 | Set-Content -Encoding utf8 build\pn2d_tdr_tie_probe\vela\simulation_iv_fine_step.json
.\build\vela_example_runner.exe --config build\pn2d_tdr_tie_probe\vela\simulation_iv_fine_step.json
python scripts\compare_reference_curves.py `
  --reference build\pn2d_tdr_tie_probe\reference_curves\pn2d_iv_reference.csv `
  --candidate build\pn2d_tdr_tie_probe\vela\pn2d_iv_fine_step.csv `
  --output-json build\pn2d_tdr_tie_probe\reports\pn2d_iv_fine_step_comparison.json `
  --output-md build\pn2d_tdr_tie_probe\reports\pn2d_iv_fine_step_comparison.md `
  --kind iv `
  --candidate-column current_total_A_per_um `
  --candidate-scale -1.0 `
  --bias-min 0.2 `
  --bias-max 0.3 `
  --require-trend-match
```

Expected:
- If orders drops materially, interpolation/sweep resolution is part of the issue.
- If orders remains around `0.5`, IV model slope is the issue.

- [ ] **Step 2: Decide gate change**

If fine sweep improves IV without making runtime unacceptable:
- change `reference_tcad/pn2d/pn2d_reference.json` IV `vela_step` from `0.1` to `0.02`;
- update sample integration expected row count if needed;
- rerun pn2d sample integration.

If fine sweep does not improve:
- leave default gate unchanged;
- proceed to Task 5.

---

## Task 5: Recombination Model Parity Investigation

**Files:**
- Modify: `tests/test_recombination.cpp`
- Modify: `include/vela/physics/RecombinationModel.h`
- Modify: `src/physics/RecombinationModel.cpp`
- Modify: `src/solver/NewtonSolver.cpp`
- Modify: `src/solver/GummelSolver.cpp`
- Modify: `docs/config_schema.md`
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`

- [ ] **Step 1: Write failing config parse test for Auger coefficients**

In `tests/test_recombination.cpp`, add:

```cpp
TEST_CASE("RecombinationModelConfig accepts explicit Auger coefficients", "[recombination][json]")
{
    nlohmann::json j = {
        {"recombination", {"srh", "auger"}},
        {"taun", 1.0e-6},
        {"taup", 2.0e-6},
        {"auger_cn_m6_s", 1.0e-43},
        {"auger_cp_m6_s", 2.0e-43},
    };
    const auto cfg = newtonConfigFromJson(j);
    REQUIRE(cfg.taun == Catch::Approx(1.0e-6));
    REQUIRE(cfg.taup == Catch::Approx(2.0e-6));
    REQUIRE(cfg.augerCn == Catch::Approx(1.0e-43));
    REQUIRE(cfg.augerCp == Catch::Approx(2.0e-43));
}
```

Expected before implementation: compile or test failure because solver configs
do not expose Auger coefficient overrides.

- [ ] **Step 2: Implement explicit Auger coefficient config**

Add fields to `GummelConfig` and `NewtonConfig`:

```cpp
Real augerCn = 2.8e-43;
Real augerCp = 9.9e-44;
```

Parse JSON keys:

```cpp
cfg.augerCn = json.value("auger_cn_m6_s", cfg.augerCn);
cfg.augerCp = json.value("auger_cp_m6_s", cfg.augerCp);
```

Pass them into `RecombinationModelConfig` before constructing the model.

- [ ] **Step 3: Add BV recombination coefficient scan**

Extend `scripts/scan_pn2d_iv_bv_physics_matrix.ps1` with BV cases that vary:

```text
taun/taup: 1e-6, 1e-7, 1e-8
auger scale: 0.0, 0.1, 1.0
```

Expected:
- Identify whether SRH lifetime or Auger coefficient magnitude drives the
  `1.17` orders BV mismatch.

- [ ] **Step 4: Verify**

```powershell
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R "recombination|mobility|newton|sentaurus_sample|ascii_sources"
```

Expected: pass.

---

## Task 6: Carrier Statistics Triage for IV Slope

**Files:**
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`
- Optional future C++:
  - `include/vela/physics/CarrierStatistics.h`
  - `src/physics/CarrierStatistics.cpp`
  - `tests/test_recombination.cpp`
  - `tests/test_sg_flux.cpp`

- [ ] **Step 1: Document current limitation**

Record that both IV and BV command files request `Fermi`, while Vela currently
uses Boltzmann carrier statistics.

- [ ] **Step 2: Add a small analytic experiment**

Before changing solver equations, add a local document table estimating
degeneracy relevance at `1e17 cm^-3` using current effective density values.
Use this only as triage; do not alter transport equations from this estimate.

- [ ] **Step 3: Decide implementation path**

If physics matrix and fine sweep do not explain IV slope:
- write a separate plan for Fermi-Dirac carrier statistics support;
- keep it out of the current calibration branch unless scoped tests are ready.

---

## Final Verification

After any implementation tasks:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Expected:
- `100% tests passed`.
- pn2d IV remains strict Newton and trend-matched.
- pn2d BV remains strict Newton and within the promoted low-bias gate.

If full suite is too slow for an intermediate checkpoint, run:

```powershell
ctest --test-dir build --output-on-failure -R "SentaurusTdrReader|sentaurus_sample|reference_tcad|ascii_sources|recombination|newton|dc_sweep"
```

and explicitly record that full suite is pending.
