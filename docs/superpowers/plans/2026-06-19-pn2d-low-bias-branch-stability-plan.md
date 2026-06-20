# PN2D Low-Bias Branch And Stability Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Separate the PN2D low-reverse-bias physical mismatch from nonlinear damping/continuation stability effects, then decide whether Vela needs a Bank-Rose-like damping mode, a DEVSIM-style variable update limiter, or only model/boundary-condition fixes before continuing Sentaurus-default BV reproduction.

**Architecture:** Treat the local forward/reverse window report as the near-field error baseline, and use Sentaurus VM runs only to fill missing reference data. Keep solver-global damping experiments, variable-specific update limiting, and physical-model comparisons as separate branches so a convergence improvement cannot be mistaken for a physics fix. Use Charon as the reference for NOX/continuation-style Newton globalization and DEVSIM as the reference for variable-level `LOGDAMP`/`POSITIVE` update limiting.

**Tech Stack:** C++20, CMake/Ninja, MSYS2 UCRT64, existing Vela `DCSweep`/`NewtonSolver`, Python standard library plus NumPy diagnostics, Sentaurus 2018 VM over SSH, local references under `build-release/reference_tcad/pn2d_sentaurus2018`.

---

## Evidence Baseline

- Existing report: `build-release/reference_tcad/pn2d_sentaurus2018/reports/forward_reverse_windows/summary.md`.
- Forward local Sentaurus reference currently covers only `0` to `+2 V`; Vela diagnostic sweep already reaches `+5 V`.
- Forward current magnitude error:
  - `+0.5 V`: Vela lower by `21.8%`, `0.107` decades.
  - `+1.0 V`: Vela lower by `4.8%`, `0.021` decades.
  - `+1.5 V`: Vela lower by `1.4%`, `0.006` decades.
  - `+2.0 V`: Vela lower by `0.7%`, `0.003` decades.
- Reverse current magnitude error:
  - `-0.5 V`: Vela lower by `40.3%`, `0.224` decades.
  - `-2.0 V`: Vela lower by `35.3%`, `0.189` decades.
  - `-5.0 V`: Vela lower by `36.0%`, `0.194` decades.
- Reverse field mismatch:
  - potential RMS is `9.6` to `11.5 mV`;
  - electron/hole density p95 error is `0.47` to `0.58` decades, about a local factor of `3`;
  - junction electric-field relative p95 improves from `0.929` at `-0.5 V` to `0.424` at `-5 V`.
- Interpretation: the low reverse mismatch appears before avalanche turn-on, so do not tune avalanche coefficients or high-field generation to explain `-0.5` to `-5 V`.

## External-Code Lessons To Carry Forward

- Charon uses Trilinos NOX/LOCA rather than a local Bank-Rose implementation. The closest open-source analogue is `NOX + Newton + Line Search Based + continuation/rescue`, with examples enabling `Rescue Bad Newton Solve`.
- DEVSIM uses a direct Newton loop and stabilizes updates at the variable level:
  - `LOGDAMP` compresses large variable updates to a thermal-voltage-like scale.
  - `POSITIVE` prevents positive variables from crossing zero.
  - It also tracks consecutive divergence through `maximum_divergence`.
- Vela already has global backtracking line search, `damping_factor`, `max_update`, Gummel `damping_psi`, step retry/shrink, and an opt-in `psi-phin` branch guard. Missing pieces are:
  - a named Bank-Rose-like damping policy with recorded damping decisions;
  - variable-specific quasi-Fermi/log-carrier update limiting;
  - a controlled sensitivity matrix proving whether those mechanisms change only convergence or also the selected low-reverse-bias branch.

## Files To Modify Or Use

- Use: `build-release/reference_tcad/pn2d_sentaurus2018/reports/forward_reverse_windows/summary.md`
  - Source baseline for forward/reverse current and field errors.
- Use: `docs/superpowers/plans/2026-06-17-pn2d-bv-minus-20v-validation.md`
  - Historical BV branch diagnostics and already implemented `psi-phin` continuation guard.
- Use: `scripts/run_sentaurus_vm_reference.py`
  - Regenerate missing Sentaurus forward `+2` to `+5 V` multibias TDRs.
- Use or modify: `scripts/compare_pn2d_bv_multibias_fields.py`
  - Add a forward/reverse window mode only if current ad hoc outputs need repeatability.
- Create: `scripts/run_pn2d_solver_sensitivity_matrix.py`
  - Generate Vela configs for damping/line-search/step-size/branch-guard sweeps and summarize IV/field metrics.
- Test: `tests/regression/test_reference_tcad_tools.py`
  - Cover sensitivity-matrix config generation and summary parsing.
- Modify: `include/vela/solver/NewtonSolver.h`, `src/solver/NewtonSolver.cpp`, `include/vela/numerics/LineSearch.h`, `src/numerics/LineSearch.cpp`
  - Only after the sensitivity matrix proves global damping is worth productizing.
- Test: `tests/test_line_search_backtrack_failure.cpp`, `tests/test_newton_solver.cpp`
  - Cover any new Bank-Rose-like damping policy.
- Modify: `include/vela/equation/CoupledDDAssembler.h`, `src/equation/CoupledDDAssembler.cpp`, `src/solver/NewtonSolver.cpp`
  - Only if adding DEVSIM-style variable update limiting requires block-aware update transforms.
- Test: `tests/test_newton_solver.cpp`, `tests/test_dc_sweep.cpp`
  - Cover variable update limiting and DCSweep diagnostics.
- Modify: `docs/config_schema.md`
  - Document any promoted solver knobs.

## Task 1: Fill The Forward +2 V To +5 V Sentaurus Reference Gap

- [ ] **Step 1: Confirm local forward TDR coverage**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
(Get-ChildItem reference_tcad\pn2d_sentaurus2018\source -Filter 'pn2d_iv_multibias_*_des.tdr').Count
Get-ChildItem reference_tcad\pn2d_sentaurus2018\source -Filter 'pn2d_iv_multibias_*_des.tdr' | Select-Object -Last 5 -ExpandProperty Name
```

Expected: current local source shows `41` files and ends at `pn2d_iv_multibias_0040_des.tdr`.

- [ ] **Step 2: Run the VM forward deck without overwriting committed source**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\run_sentaurus_vm_reference.py `
  --run-id pn2d_iv_forward_0_to_5_refresh `
  --source-dir reference_tcad\pn2d_sentaurus2018\source `
  --deck pn2d_iv_sdevice.cmd `
  --dry-run
```

Expected: dry-run manifest lists the remote deck, expected `pn2d_iv_multibias_*_des.tdr`, `.plt`, and `.log` artifacts.

- [ ] **Step 3: Execute the live VM run after reviewing dry-run output**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\run_sentaurus_vm_reference.py `
  --run-id pn2d_iv_forward_0_to_5_refresh `
  --source-dir reference_tcad\pn2d_sentaurus2018\source `
  --deck pn2d_iv_sdevice.cmd
```

Expected: artifacts are copied under `build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_vm_runs/pn2d_iv_forward_0_to_5_refresh/source`.

- [ ] **Step 4: Import +2.5, +3, +4, and +5 V snapshots**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
build-release\sentaurus_import.exe --tdr build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_iv_forward_0_to_5_refresh\source\pn2d_iv_multibias_0050_des.tdr --export-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_iv_multibias\sentaurus_2.5v > build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_iv_multibias\sentaurus_2.5v_import.log
build-release\sentaurus_import.exe --tdr build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_iv_forward_0_to_5_refresh\source\pn2d_iv_multibias_0060_des.tdr --export-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_iv_multibias\sentaurus_3v > build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_iv_multibias\sentaurus_3v_import.log
build-release\sentaurus_import.exe --tdr build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_iv_forward_0_to_5_refresh\source\pn2d_iv_multibias_0080_des.tdr --export-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_iv_multibias\sentaurus_4v > build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_iv_multibias\sentaurus_4v_import.log
build-release\sentaurus_import.exe --tdr build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_vm_runs\pn2d_iv_forward_0_to_5_refresh\source\pn2d_iv_multibias_0100_des.tdr --export-dir build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_iv_multibias\sentaurus_5v > build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_iv_multibias\sentaurus_5v_import.log
```

Expected: each export contains `nodes.csv`, `elements.csv`, and `fields/ElectrostaticPotential_region0.csv`, `fields/ElectricField_region0.csv`, `fields/eDensity_region0.csv`, `fields/hDensity_region0.csv`.

## Task 2: Build A Repeatable Low-Bias Sensitivity Matrix

- [ ] **Step 1: Add a RED regression test for sensitivity config generation**

Add to `tests/regression/test_reference_tcad_tools.py`:

```python
def test_pn2d_solver_sensitivity_matrix_writes_configs(self):
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)
        cmd = [
            sys.executable,
            "scripts/run_pn2d_solver_sensitivity_matrix.py",
            "--base-config",
            "build-release/reference_tcad/pn2d_sentaurus2018/vela/simulation_bv.json",
            "--out-dir",
            str(out),
            "--dry-run",
            "--bias-window",
            "reverse_low",
        ]
        result = subprocess.run(cmd, cwd=REPO_ROOT, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        manifest = json.loads((out / "manifest.json").read_text())
        names = {case["name"] for case in manifest["cases"]}
        self.assertIn("baseline", names)
        self.assertIn("bank_rose_like_damped", names)
        self.assertIn("devsim_like_qf_limited", names)
```

- [ ] **Step 2: Run the RED test**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_solver_sensitivity_matrix_writes_configs
```

Expected: fails because `scripts/run_pn2d_solver_sensitivity_matrix.py` does not exist.

- [ ] **Step 3: Create `scripts/run_pn2d_solver_sensitivity_matrix.py`**

Implement a Python script with these dry-run cases:

```json
[
  {"name": "baseline", "solver": {}},
  {"name": "newton_damping_0p5", "solver": {"damping_factor": 0.5}},
  {"name": "newton_max_update_0p05", "solver": {"max_update": 0.05}},
  {"name": "branch_guard_0p05", "continuation": {"branch_acceptance": {"psi_phin_jump": true, "max_psi_phin_jump_V": 0.05}}},
  {"name": "bank_rose_like_damped", "solver": {"damping_factor": 0.5, "line_search": true, "max_update": 0.05}},
  {"name": "qf_hard_limit_0p0259", "solver": {"quasi_fermi_update_limit_V": 0.0259}}
]
```

The script must:

- copy the base JSON into case-specific config files;
- set reverse low-bias sweep to `start: 0`, `stop: -5`, `step: -0.05`, `write_vtk: true`;
- set `output_csv` and `vtk_prefix` inside each case directory;
- write `manifest.json`;
- in `--dry-run`, skip executing Vela.

- [ ] **Step 4: Run GREEN test**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_solver_sensitivity_matrix_writes_configs
```

Expected: passes.

## Task 3: Run Solver-Knob Sensitivity Before Adding New Solver Code

- [ ] **Step 1: Generate sensitivity configs**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\run_pn2d_solver_sensitivity_matrix.py `
  --base-config build-release\reference_tcad\pn2d_sentaurus2018\vela\simulation_bv.json `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\low_bias_solver_sensitivity `
  --bias-window reverse_low
```

Expected: case configs and a manifest are generated. Cases requiring not-yet-implemented `quasi_fermi_update_limit_V` are marked `unsupported` instead of run.

- [ ] **Step 2: Compare current and field metrics for supported cases**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\compare_pn2d_bv_multibias_fields.py `
  --sentaurus-root build-release\reference_tcad\pn2d_sentaurus2018\sentaurus_multibias `
  --vela-vtk-root build-release\reference_tcad\pn2d_sentaurus2018\reports\low_bias_solver_sensitivity\baseline\vtk `
  --curve-reference build-release\reference_tcad\pn2d_sentaurus2018\reference_curves\pn2d_sentaurus2018_bv_reference.csv `
  --curve-candidate build-release\reference_tcad\pn2d_sentaurus2018\reports\low_bias_solver_sensitivity\baseline\iv.csv `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018\reports\low_bias_solver_sensitivity\baseline_compare `
  --biases 0,-0.5,-2,-5 `
  --quantities potential,electric_field,electron_density,hole_density
```

Expected: baseline reproduces the existing low-bias mismatch within rounding.

- [ ] **Step 3: Classify damping sensitivity**

Pass criteria:

- If `newton_damping_0p5` and `newton_max_update_0p05` converge to the same reverse IV within `2%` at `-0.5`, `-2`, and `-5 V`, then low-bias mismatch is not a Newton damping artifact.
- If damping changes the converged current by more than `5%`, inspect the final field states; do not promote any damping knob until the final solution branch is understood.

### Task 3 Execution Note - 2026-06-19

- Generated sensitivity configs and fixed the generator so copied input paths resolve from `--base-config` and generated output paths are absolute.
- Fixed `branch_guard_0p05` generation so branch controls are written under `sweep.continuation`, where `DCSweep` reads them.
- Baseline ran through `-5 V`; comparison artifacts are under `build-release/reference_tcad/pn2d_sentaurus2018/reports/low_bias_solver_sensitivity/baseline_compare`.
- Fixed `newton_damping_0p5` was stopped around `-1.5 V`; it has 135 valid rows plus one truncated row. It did not satisfy Task 3 pass criteria.
- At `-0.5 V`, fixed damping changes Vela `current_total_A_per_um` from `-3.148e-18` baseline to `-1.134e-12`, while global potential/carrier field metrics remain close. Damping-specific comparison artifacts are under `build-release/reference_tcad/pn2d_sentaurus2018/reports/low_bias_solver_sensitivity/newton_damping_0p5_compare`.
- Bounded remaining runs show `branch_guard_0p05` completes to `-5 V` and matches the baseline-scale IV/field comparison, while `newton_max_update_0p05` and `bank_rose_like_damped` fail immediately at the initial 0 V Newton solve with `max_iterations`.
- Classification: Task 3 is diagnostic-only/incomplete. Do not conclude the mismatch is not a damping artifact, and do not promote fixed damping or strong global max-update limiting. Next step is to diagnose the 0 V failure mechanism or implement an opt-in quasi-Fermi/log-carrier update limiter that can be tested against the same matrix.

## Task 4: Add Opt-In Quasi-Fermi Update Limiting Only If Needed

Execution note: the first implemented limiter is a hard physical-voltage cap on
`phin`/`phip` Newton updates, not a full DEVSIM `LOGDAMP` transform. Keep the
case name and documentation explicit (`qf_hard_limit_0p0259`). A true smooth
log-damping transform remains a follow-up only if the hard cap improves the
diagnostic matrix without distorting the converged branch.

Execution update: `qf_hard_limit_0p0259` completes the 0 to -5 V low-bias BV
window and reproduces the same IV/field metrics as baseline. It is
non-regressive, but it does not reduce the Sentaurus low-bias current or field
errors. This points the next debug step back to global update-limit 0 V
stagnation and the remaining physics/current-extraction mismatch rather than to
promoting hard qF limiting as a calibration fix.

Max-update diagnostic update: the 0 V failure for `newton_max_update_0p05` and
`bank_rose_like_damped` is caused by a strong global infinity-norm cap applied
to the initial Poisson-dominated equilibrium solve. A one-point 0 V scan shows
`max_update <= 0.2` fails in 40 Newton iterations, `0.5` converges in 28
iterations, and `5` converges in 6 iterations. The residual remains Poisson
dominated with positive finite carriers. Do not use fixed global
`max_update = 0.05` as a Bank-Rose proxy; return to physics/current-extraction
debug unless a future adaptive, block-aware trust-region policy is explicitly
designed and tested.

BV full-quantity comparison update: `baseline_compare_full` confirms the low
reverse-bias current error is about `0.19` to `0.23` decade from `-0.5` to
`-5 V`, while electrostatic potential RMS error is only about `0.01 V`.
The largest stable physical-state discrepancies are carrier density
(`~0.48` to `0.57` log10 p95) plus mobility/impact-ionization derived fields.
Because branch guard and qF hard limiting are non-regressive but unchanged, the
next implementation/debug branch should inspect carrier-density reconstruction,
mobility/Einstein/SG flux coefficients, and terminal/contact-edge transport
drivers rather than additional fixed Newton damping.

- [ ] **Step 1: Add RED parser and behavior tests**

Add to `tests/test_newton_solver.cpp`:

```cpp
TEST_CASE("NewtonSolver: parses quasi-Fermi update limit", "[newton][config]")
{
    const auto cfg = newtonConfigFromJson(nlohmann::json{
        {"quasi_fermi_update_limit_V", 0.0259}
    });
    REQUIRE(cfg.quasiFermiUpdateLimit_V == Catch::Approx(0.0259));
}
```

Expected RED: `NewtonConfig` has no `quasiFermiUpdateLimit_V`.

- [ ] **Step 2: Implement config storage without changing defaults**

Add to `include/vela/solver/NewtonSolver.h`:

```cpp
Real quasiFermiUpdateLimit_V = 0.0; ///< 0 disables DEVSIM-style log damping.
```

Parse in `src/solver/NewtonSolver.cpp`:

```cpp
cfg.quasiFermiUpdateLimit_V =
    json.value("quasi_fermi_update_limit_V", cfg.quasiFermiUpdateLimit_V);
```

Validate:

```cpp
if (cfg.quasiFermiUpdateLimit_V < 0.0 || !std::isfinite(cfg.quasiFermiUpdateLimit_V))
    throw std::invalid_argument("newtonConfigFromJson: quasi_fermi_update_limit_V must be non-negative and finite.");
```

- [ ] **Step 3: Apply the limiter to Newton quasi-Fermi blocks before line search**

In `NewtonSolver::solve`, after solving `step = J^-1 (-r)` and after `max_update`, limit only the `phin` and `phip` blocks:

```cpp
if (cfg_.quasiFermiUpdateLimit_V > 0.0) {
    const Index nNodes = mesh_.numNodes();
    const Real limit = useScaledUnknowns ? cfg_.quasiFermiUpdateLimit_V / scaling.V0
                                         : cfg_.quasiFermiUpdateLimit_V;
    for (Index i = nNodes; i < 3 * nNodes; ++i) {
        const Real raw = step(static_cast<int>(i));
        if (std::abs(raw) > limit) {
            const Real sign = raw > 0.0 ? 1.0 : -1.0;
            step(static_cast<int>(i)) = sign * limit * std::log1p(std::abs(raw) / limit);
        }
    }
}
```

Use the actual local variable names from `NewtonSolver::solve`; if `useScaledUnknowns` or `scaling` is not in scope, derive them from the existing DD scaling spec already used by the solver.

- [ ] **Step 4: Run focused tests**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build-release --parallel --target test_newton_solver test_dc_sweep
build-release\test_newton_solver.exe "[newton][config]"
build-release\test_dc_sweep.exe "[dc_sweep]"
```

Expected: tests pass and default behavior remains unchanged.

## Task 5: Add A Named Bank-Rose-Like Damping Policy Only If Global Damping Helps

- [ ] **Step 1: Add RED line-search policy tests**

Add to `tests/test_line_search_backtrack_failure.cpp`:

```cpp
TEST_CASE("BacktrackingLineSearch: bank-rose policy starts from configured damping", "[line_search]")
{
    LineSearchConfig cfg;
    cfg.enabled = true;
    cfg.policy = "bank_rose";
    cfg.initialDamping = 0.25;
    cfg.recordHistory = true;
    BacktrackingLineSearch search(cfg);
    VectorXd x(1); x << 0.0;
    VectorXd step(1); step << 1.0;
    VectorXd r(1); r << 10.0;
    const auto result = search.search(
        x, step, r,
        [](const VectorXd& candidate) {
            VectorXd residual(1);
            residual << (candidate(0) < 0.3 ? 8.0 : 20.0);
            return residual;
        });
    REQUIRE(result.accepted);
    REQUIRE(result.damping == Catch::Approx(0.25));
}
```

Expected RED: `LineSearchConfig` has no `policy`.

- [ ] **Step 2: Implement policy as an alias first**

Add to `include/vela/numerics/LineSearch.h`:

```cpp
std::string policy = "backtracking"; ///< "backtracking" or "bank_rose".
```

In `BacktrackingLineSearch`, accept both values and initially make `bank_rose` share the same sufficient-decrease backtracking mechanics. This creates a stable configuration surface before adding Sentaurus-inspired damping heuristics.

- [ ] **Step 3: Add diagnostics columns**

Ensure Newton history records:

```text
damping_policy
damping_factor
line_search_attempts
line_search_failure_reason
```

Expected: existing diagnostic CSVs continue to parse and gain policy visibility for new runs.

## Task 6: Decide Physics Versus Stability Based On The Matrix

- [ ] **Step 1: Write a decision note**

Create:

```text
build-release/reference_tcad/pn2d_sentaurus2018/reports/low_bias_solver_sensitivity/decision.md
```

It must include:

- low-bias IV table for all supported cases;
- potential/electric-field/electron-density/hole-density field table at `-0.5`, `-2`, and `-5 V`;
- whether damping changes the final converged branch;
- whether quasi-Fermi limiting improves convergence only, current parity only, or both.

- [ ] **Step 2: Choose the next implementation branch**

Use these rules:

- If all damping variants preserve the same IV and field mismatch, open a physical-model branch: contact carrier statistics, intrinsic density/BGN, SRH lifetime defaults, and contact minority reconstruction.
- If damping changes the selected branch but worsens Sentaurus field parity, keep it diagnostic-only.
- If quasi-Fermi limiting improves field parity without harming forward `+1` to `+2 V`, promote it behind an opt-in config and add schema docs.
- If Bank-Rose-like damping only helps the avalanche-enabled `-13.2 V` failure, keep it under BV high-field plan and do not use it to explain low-bias leakage.

## Task 7: Final Verification

- [ ] **Step 1: Run focused tests**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python -m unittest tests.regression.test_reference_tcad_tools tests.regression.test_sentaurus_vm_reference_runner
ctest --test-dir build-release --output-on-failure -R "newton|dc_sweep|reference_tcad|sentaurus_import"
```

Expected: all focused tests pass.

- [ ] **Step 2: Run full CTest if C++ solver behavior changed**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
ctest --test-dir build-release --output-on-failure
```

Expected: full suite passes before merging any solver or schema change.

## Self-Review

- Spec coverage: plan covers forward `+2` to `+5 V` reference gap, reverse low-bias current/field mismatch, Charon NOX/continuation lessons, DEVSIM variable update limiting, Vela damping sensitivity, and final branch decision.
- Red-flag scan: no red-flag wording or unspecified implementation steps remain.
- Type consistency: new proposed config names are `quasi_fermi_update_limit_V` in JSON and `quasiFermiUpdateLimit_V` in C++; line-search policy uses `policy = "bank_rose"` as a named alias before heuristic changes.
