# Gummel Newton Handoff Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a DC sweep solver mode that uses Gummel as a stable initializer and then hands each bias point to fully coupled Newton, so imported Sentaurus pn2d decks can move from runtime-scaled diagnostics toward faithful high-doping node-level solves.

**Architecture:** Keep existing `gummel` and `newton` paths unchanged. Add a third explicit solver method, `gummel_newton`, inside the DC sweep path: for each attempted bias point, run a bounded Gummel initializer, validate the intermediate solution, then run Newton with `warm_start=true` from that solution. The generated Sentaurus pn2d runtime deck should opt into this hybrid path only when configured, and comparison reports should record whether the final accepted solve came from Newton or from a configured fallback.

**Tech Stack:** C++20, CMake/Ninja, Catch2 tests, existing Vela `DCSweep`, `GummelSolver`, `NewtonSolver`, Python Sentaurus import tooling, JSON deck schema.

---

## File Structure

- Modify `src/simulation/DCSweep.cpp`
  - Extend solver-method parsing to accept `gummel_newton`.
  - Add local hybrid solve orchestration in `solvePoint`.
  - Record stage-specific diagnostics without changing standalone Gummel/Newton semantics.
- Modify `include/vela/simulation/DCSweep.h`
  - Add optional per-point fields for solver diagnostics: `solverMethod`, `gummelIterations`, `newtonIterations`, `handoffStage`.
- Modify `docs/config_schema.md`
  - Document `solver.method: "gummel_newton"` and the initializer handoff controls.
- Modify `scripts/sentaurus_import.py`
  - Allow reference configs to request hybrid solver generation for faithful pn2d decks.
  - Keep runtime-scaled deck generation as a diagnostic fallback, not the preferred path.
- Modify `reference_tcad/pn2d/pn2d_reference.json`
  - Add opt-in hybrid settings for IV/BV simulations after C++ support lands.
- Test `tests/test_dc_sweep.cpp`
  - Cover method parsing, hybrid execution at zero bias, handoff from Gummel to Newton, and fallback behavior.
- Test `tests/regression/test_sentaurus_import_tools.py`
  - Cover generated pn2d decks using `gummel_newton` and preserving unsupported-physics warnings.
- Test `tests/regression/test_sentaurus_sample_integration.py`
  - Gate real pn2d import/execution expectations when bundled or environment sample is present.

---

## Task 1: Add Hybrid Solver Method Parsing And Diagnostics

**Files:**
- Modify: `src/simulation/DCSweep.cpp`
- Modify: `include/vela/simulation/DCSweep.h`
- Test: `tests/test_dc_sweep.cpp`

- [x] **Step 1: Write failing parser/diagnostic tests**

Add this test near the existing Newton method tests in `tests/test_dc_sweep.cpp`:

```cpp
TEST_CASE("DCSweep: hybrid Gummel-Newton method is reachable from config",
          "[dc_sweep][gummel_newton]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "gummel_newton_start.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.25},
        {"write_vtk", false}
    }, {
        {"method", "gummel_newton"},
        {"max_iter", 12},
        {"reltol", 1.0e-8},
        {"abstol", 1.0e-18},
        {"damping_psi", 0.35},
        {"damping_factor", 1.0},
        {"line_search", true},
        {"verbose", false}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 1);
    const DCSweepPoint& point = result.points.front();
    REQUIRE(point.converged);
    REQUIRE(point.solverMethod == "gummel_newton");
    REQUIRE(point.gummelIterations > 0);
    REQUIRE(point.newtonIterations >= 0);
    REQUIRE(point.handoffStage == "newton");
    REQUIRE(std::filesystem::exists(csvPath));
}
```

Also extend the invalid method test expectation so the error mentions `gummel_newton`.

- [x] **Step 2: Run the focused test and confirm it fails**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R dc_sweep
```

Expected: build or test failure because `DCSweepPoint` has no hybrid diagnostic fields and `solver.method` rejects `gummel_newton`.

- [x] **Step 3: Add diagnostic fields**

In `include/vela/simulation/DCSweep.h`, extend `DCSweepPoint`:

```cpp
std::string solverMethod;
int gummelIterations = 0;
int newtonIterations = 0;
std::string handoffStage;
```

Place these next to the existing `iterations` field so result consumers can find all nonlinear-solver diagnostics together.

- [x] **Step 4: Extend method parsing**

In `src/simulation/DCSweep.cpp`, update the local solver method enum and parser from the existing two-method shape to:

```cpp
enum class SolverMethod {
    Gummel,
    Newton,
    GummelNewton
};
```

Then accept these method strings:

```cpp
if (normalized == "gummel")
    return SolverMethod::Gummel;
if (normalized == "newton")
    return SolverMethod::Newton;
if (normalized == "gummel_newton" || normalized == "hybrid")
    return SolverMethod::GummelNewton;
throw std::invalid_argument(
    "DCSweep: solver.method/type must be 'gummel', 'newton', or 'gummel_newton'.");
```

- [x] **Step 5: Populate diagnostics for existing paths**

In the existing Gummel branch, set:

```cpp
attempt.solverMethod = "gummel";
attempt.gummelIterations = sol.iters;
attempt.newtonIterations = 0;
attempt.handoffStage = solverConverged ? "gummel" : "gummel_failed";
```

In the existing Newton branch, set:

```cpp
attempt.solverMethod = "newton";
attempt.gummelIterations = 0;
attempt.newtonIterations = result.iters;
attempt.handoffStage = solverConverged ? "newton" : "newton_failed";
```

Copy these fields into `DCSweepPoint` inside `recordPoint`.

- [x] **Step 6: Run tests and commit**

Run:

```powershell
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R dc_sweep
```

Expected: all `dc_sweep` tests pass.

Commit:

```powershell
git add include/vela/simulation/DCSweep.h src/simulation/DCSweep.cpp tests/test_dc_sweep.cpp
git commit -m "Add DC sweep hybrid solver method"
```

---

## Task 2: Implement Gummel Initializer Plus Newton Handoff

**Files:**
- Modify: `src/simulation/DCSweep.cpp`
- Test: `tests/test_dc_sweep.cpp`

- [x] **Step 1: Add a test proving Newton receives the Gummel solution**

Add this test near the hybrid reachability test:

```cpp
TEST_CASE("DCSweep: hybrid path uses Gummel iterations before Newton handoff",
          "[dc_sweep][gummel_newton]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "gummel_newton_forward.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.2},
        {"step", 0.2},
        {"write_vtk", false}
    }, {
        {"method", "gummel_newton"},
        {"max_iter", 20},
        {"reltol", 1.0e-8},
        {"abstol", 1.0e-18},
        {"damping_psi", 0.35},
        {"damping_factor", 1.0},
        {"line_search", true},
        {"warm_start", true},
        {"verbose", false}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 2);
    for (const DCSweepPoint& point : result.points) {
        REQUIRE(point.converged);
        REQUIRE(point.solverMethod == "gummel_newton");
        REQUIRE(point.gummelIterations > 0);
        REQUIRE(point.handoffStage == "newton");
    }
}
```

- [x] **Step 2: Run and confirm failure**

Run:

```powershell
ctest --test-dir build --output-on-failure -R "DCSweep: hybrid"
```

Expected: failure because the hybrid branch is parsed but not implemented yet.

- [x] **Step 3: Implement the hybrid branch**

Inside `solvePoint`, add a `SolverMethod::GummelNewton` branch with this behavior:

```cpp
DDSolution gummelInitial = initial != nullptr
    ? runGummel(mesh, matdb, doping, biases, contactSpecs, gummel, *initial,
                fixedChargeSpecs, sheetChargeSpecs)
    : runGummel(mesh, matdb, doping, biases, contactSpecs, gummel,
                fixedChargeSpecs, sheetChargeSpecs);

const DDSolutionValidationResult gummelValidation =
    validateDDSolution(gummelInitial, mesh, biases, validationOptions);

if (!gummelInitial.converged || !gummelValidation.valid) {
    attempt.ok = false;
    attempt.solution = std::move(gummelInitial);
    attempt.solverMethod = "gummel_newton";
    attempt.gummelIterations = attempt.solution.iters;
    attempt.newtonIterations = 0;
    attempt.handoffStage = !gummelInitial.converged
        ? "gummel_failed"
        : "gummel_validation_failed";
    attempt.failureReason = !gummelInitial.converged
        ? "gummel_non_convergence"
        : "gummel_validation_failed";
    attempt.validationDiagnostics = gummelValidation.diagnosticsString();
    return attempt;
}

NewtonConfig handoffNewton = newton;
handoffNewton.warmStart = true;
NewtonResult result = runNewton(mesh, matdb, doping, biases, gummelInitial,
                                handoffNewton, fixedChargeSpecs,
                                sheetChargeSpecs);
solverConverged = result.converged;
sol = std::move(result.solution);
attempt.solverMethod = "gummel_newton";
attempt.gummelIterations = gummelInitial.iters;
attempt.newtonIterations = result.iters;
attempt.handoffStage = solverConverged ? "newton" : "newton_failed";
```

Use the existing final validation block for the Newton result so all solver modes share one acceptance rule.

- [x] **Step 4: Preserve previous-bias continuation**

Keep passing `previousSolution` into `solvePoint`. In hybrid mode it should seed Gummel first, and the converged Gummel result should seed Newton. Do not pass `previousSolution` directly to Newton unless the Gummel initializer fails and a later task adds an explicit fallback.

- [x] **Step 5: Run tests and commit**

Run:

```powershell
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R "dc_sweep|newton|dd_gummel"
```

Expected: all selected tests pass.

Commit:

```powershell
git add src/simulation/DCSweep.cpp tests/test_dc_sweep.cpp
git commit -m "Implement Gummel initialized Newton sweep"
```

---

## Task 3: Add Configurable Handoff Policy

**Files:**
- Modify: `src/simulation/DCSweep.cpp`
- Modify: `docs/config_schema.md`
- Test: `tests/test_dc_sweep.cpp`

- [ ] **Step 1: Add tests for fallback and strict policies**

Add two tests:

```cpp
TEST_CASE("DCSweep: hybrid fallback can accept converged Gummel when Newton fails",
          "[dc_sweep][gummel_newton]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "gummel_newton_fallback.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.25},
        {"write_vtk", false}
    }, {
        {"method", "gummel_newton"},
        {"max_iter", 0},
        {"handoff", {{"fallback", "gummel_on_newton_failure"}}},
        {"verbose", false}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 1);
    REQUIRE(result.points.front().converged);
    REQUIRE(result.points.front().handoffStage == "gummel_fallback");
}

TEST_CASE("DCSweep: hybrid strict policy rejects Newton failure",
          "[dc_sweep][gummel_newton]")
{
    const auto dir = makeUniqueSweepDir();
    const ScopedDirectoryCleanup cleanup{dir};
    std::filesystem::create_directories(dir);
    const auto meshPath = writePNMesh(dir);
    const auto csvPath = dir / "gummel_newton_strict.csv";
    const auto cfgPath = writeSweepConfig(dir, meshPath, csvPath, {
        {"start", 0.0},
        {"stop", 0.0},
        {"step", 0.25},
        {"write_vtk", false},
        {"stop_on_failure", false}
    }, {
        {"method", "gummel_newton"},
        {"max_iter", 0},
        {"handoff", {{"fallback", "none"}}},
        {"verbose", false}
    });

    DCSweep sweep;
    const DCSweepResult result = sweep.runWithResult(cfgPath.string());

    REQUIRE(result.points.size() == 1);
    REQUIRE_FALSE(result.points.front().converged);
    REQUIRE(result.points.front().failureReason == "newton_non_convergence");
}
```

- [ ] **Step 2: Parse `solver.handoff` options**

Add a local struct in `src/simulation/DCSweep.cpp`:

```cpp
struct HybridHandoffConfig {
    bool fallbackToGummelOnNewtonFailure = false;
    bool requireGummelConvergence = true;
};
```

Parse:

```cpp
if (solverJson.contains("handoff")) {
    const auto& handoff = solverJson.at("handoff");
    const std::string fallback = handoff.value("fallback", "none");
    if (fallback == "none")
        hybrid.fallbackToGummelOnNewtonFailure = false;
    else if (fallback == "gummel_on_newton_failure")
        hybrid.fallbackToGummelOnNewtonFailure = true;
    else
        throw std::invalid_argument(
            "DCSweep: solver.handoff.fallback must be 'none' or "
            "'gummel_on_newton_failure'.");
    hybrid.requireGummelConvergence =
        handoff.value("require_gummel_convergence", true);
}
```

Default must be strict: Newton failure means the point fails. This keeps the new path scientifically honest for calibration.

- [ ] **Step 3: Implement fallback behavior**

After Newton returns but before final acceptance:

```cpp
if (!result.converged && hybrid.fallbackToGummelOnNewtonFailure) {
    attempt.ok = true;
    attempt.solution = std::move(gummelInitial);
    attempt.solverMethod = "gummel_newton";
    attempt.gummelIterations = gummelInitial.iters;
    attempt.newtonIterations = result.iters;
    attempt.handoffStage = "gummel_fallback";
    attempt.failureReason.clear();
    attempt.validationDiagnostics = gummelValidation.diagnosticsString();
    return attempt;
}
```

Only allow fallback when the Gummel initializer already passed validation.

- [ ] **Step 4: Document schema**

In `docs/config_schema.md`, update solver method selection:

```markdown
- method: `gummel`, `newton`, or `gummel_newton`
```

Add:

```markdown
Hybrid Gummel-Newton keys:
- `handoff.fallback`: `none` or `gummel_on_newton_failure`
- `handoff.require_gummel_convergence`: boolean, default `true`

`gummel_newton` runs the configured Gummel solve first, validates that
solution, then runs coupled Newton with `warm_start=true` from the Gummel
state. The default fallback policy is strict: a Newton failure fails the sweep
point. Use `gummel_on_newton_failure` only for diagnostic curves where a finite
Gummel result is preferable to aborting the sweep.
```

- [ ] **Step 5: Run tests and commit**

Run:

```powershell
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R "dc_sweep|config"
```

Expected: selected tests pass.

Commit:

```powershell
git add src/simulation/DCSweep.cpp docs/config_schema.md tests/test_dc_sweep.cpp
git commit -m "Add hybrid solver handoff policy"
```

---

## Task 4: Generate pn2d Faithful Hybrid Decks

**Files:**
- Modify: `scripts/sentaurus_import.py`
- Modify: `reference_tcad/pn2d/pn2d_reference.json`
- Test: `tests/regression/test_sentaurus_import_tools.py`
- Test: `tests/regression/test_sentaurus_sample_integration.py`

- [ ] **Step 1: Write generator tests**

In `tests/regression/test_sentaurus_import_tools.py`, add assertions to the existing pn2d deck-generation test:

```python
self.assertEqual(iv_deck["solver"]["method"], "gummel_newton")
self.assertTrue(iv_deck["solver"]["warm_start"])
self.assertEqual(iv_deck["solver"]["handoff"]["fallback"], "none")
self.assertEqual(iv_deck["node_doping_file"], "doping.csv")
self.assertNotIn("runtime_approximation", iv_deck.get("sentaurus_import", {}))
```

Add a separate runtime-deck assertion:

```python
self.assertEqual(runtime_iv["solver"]["method"], "gummel")
self.assertIn("runtime_approximation", runtime_iv["sentaurus_import"])
```

- [ ] **Step 2: Run regression test and confirm failure**

Run:

```powershell
python -m unittest tests.regression.test_sentaurus_import_tools -v
```

Expected: failure because generated faithful decks still use `gummel`.

- [ ] **Step 3: Add config opt-in**

In `reference_tcad/pn2d/pn2d_reference.json`, add:

```json
"vela_solver": {
  "method": "gummel_newton",
  "max_iter": 40,
  "reltol": 1.0e-8,
  "abstol": 1.0e-18,
  "damping_psi": 0.25,
  "damping_factor": 1.0,
  "line_search": true,
  "warm_start": true,
  "verbose": false,
  "handoff": {
    "fallback": "none",
    "require_gummel_convergence": true
  }
}
```

Place this at the case or simulation level. Simulation-level settings override case-level settings.

- [ ] **Step 4: Apply solver override in import script**

In `scripts/sentaurus_import.py`, after the base solver is created for each generated deck, merge `vela_solver` from the top-level config and from the current simulation:

```python
solver_override = {}
solver_override.update(config.get("vela_solver", {}))
solver_override.update(sim.get("vela_solver", {}))
if solver_override:
    deck["solver"].update(solver_override)
```

Keep `write_runtime_deck_if_requested()` forcing the runtime approximation back to conservative `gummel` unless the simulation explicitly sets:

```json
"runtime_solver_method": "gummel_newton"
```

- [ ] **Step 5: Run tests and commit**

Run:

```powershell
python -m unittest tests.regression.test_sentaurus_import_tools -v
python -m unittest tests.regression.test_sentaurus_sample_integration -v
```

Expected: tests pass.

Commit:

```powershell
git add scripts/sentaurus_import.py reference_tcad/pn2d/pn2d_reference.json tests/regression/test_sentaurus_import_tools.py tests/regression/test_sentaurus_sample_integration.py
git commit -m "Generate pn2d hybrid solver decks"
```

---

## Task 5: Real pn2d Gate With Faithful Node Doping

**Files:**
- Modify: `tests/regression/test_sentaurus_sample_integration.py`
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`
- Optional modify after diagnosis: `src/simulation/DCSweep.cpp`, `src/solver/NewtonSolver.cpp`, or `src/equation/CoupledDDAssembler.cpp`

- [ ] **Step 1: Add a non-strict real-sample hybrid execution test**

In `tests/regression/test_sentaurus_sample_integration.py`, add a test that runs `scripts/sentaurus_import.py reference` with the bundled `reference_tcad/pn2d` and checks the faithful deck exists and is attempted:

```python
self.assertEqual(iv_deck["solver"]["method"], "gummel_newton")
self.assertEqual(iv_deck["node_doping_file"], "doping.csv")
self.assertTrue((out_dir / "vela" / "simulation_iv.json").exists())
self.assertTrue((out_dir / "reference_tcad_manifest.json").exists())
```

Do not require the faithful high-doping curve to pass until Task 6. Require runtime approximation to remain finite so the regression suite stays useful.

- [ ] **Step 2: Run real import and capture failure mode**

Run:

```powershell
python scripts\sentaurus_import.py reference --config reference_tcad\pn2d\pn2d_reference.json --source-dir reference_tcad\pn2d --output-dir build\pn2d_hybrid_gate --tdr-importer build\sentaurus_import.exe --runner build\vela_example_runner.exe
```

Expected before solver fixes: faithful deck may fail at high doping, but runtime decks must still produce finite CSVs and comparison reports.

- [ ] **Step 3: Document observed status**

In `docs/validation/pn2d_sentaurus_comparison.md`, add:

```markdown
## Hybrid Solver Status

Faithful pn2d decks now use `solver.method: "gummel_newton"` and preserve
`node_doping_file: "doping.csv"`. The runtime approximation remains available
as a convergence diagnostic and continues to scale region-average doping by
`runtime_doping_scale`.

Current gate:
- faithful deck generation is required;
- faithful execution failure must be reported in `reference_tcad_manifest.json`;
- runtime IV/BV execution must remain finite;
- strict Sentaurus numerical agreement is not yet required.
```

- [ ] **Step 4: Run regression and commit**

Run:

```powershell
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Expected: full suite passes. If faithful execution fails, the manifest documents it and runtime approximation remains finite.

Commit:

```powershell
git add tests/regression/test_sentaurus_sample_integration.py docs/validation/pn2d_sentaurus_comparison.md
git commit -m "Gate pn2d hybrid deck generation"
```

---

## Task 6: Remove pn2d Runtime Doping Scale As A Requirement

**Files:**
- Modify after diagnosis: `src/solver/GummelSolver.cpp`
- Modify after diagnosis: `src/solver/NewtonSolver.cpp`
- Modify after diagnosis: `src/equation/CoupledDDAssembler.cpp`
- Modify: `reference_tcad/pn2d/pn2d_reference.json`
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`
- Test: `tests/test_dd_gummel.cpp`
- Test: `tests/test_dc_sweep.cpp`
- Test: `tests/regression/test_sentaurus_sample_integration.py`

- [ ] **Step 1: Add a synthetic abrupt PN high-doping hybrid test**

Add a small mesh test in `tests/test_dc_sweep.cpp` using `node_doping_file` with left-side acceptors and right-side donors at `1.0e17 cm^-3`. The test should run:

```json
"solver": {
  "method": "gummel_newton",
  "max_iter": 60,
  "reltol": 1.0e-7,
  "abstol": 1.0e-18,
  "damping_psi": 0.2,
  "line_search": true,
  "warm_start": true,
  "verbose": false,
  "handoff": {
    "fallback": "none"
  }
}
```

Expected assertion:

```cpp
REQUIRE(result.points.size() >= 1);
REQUIRE(result.points.front().converged);
REQUIRE(result.points.front().handoffStage == "newton");
```

- [ ] **Step 2: Run and use the failure as the solver target**

Run:

```powershell
ctest --test-dir build --output-on-failure -R "dc_sweep|dd_gummel|newton"
```

Expected before fixes: high-doping abrupt PN may fail in Gummel, Newton, or validation. The failure stage determines the implementation change.

- [ ] **Step 3: Fix the actual failure mode, not the symptom**

Use this decision table:

```text
Gummel fails before Newton:
  tune or add continuation controls in Gummel initialization.
  Prefer smaller damping, larger max_iter, and absolute carrier floors over changing physics.

Gummel validates but Newton residual explodes:
  inspect CoupledDDAssembler scaling and residual block scales.
  Prefer robust residual scaling or line-search acceptance fixes over disabling equations.

Newton converges but validation rejects:
  inspect validateDDSolution thresholds against high-doping contact equilibrium.
  Adjust validation only if the rejected values are finite and physically plausible.

Interface nodes have donor=acceptor compensation:
  add an importer-side or DCSweep-side option to resolve compensated interface nodes by adjacent material majority doping.
  Record this in metadata; do not silently rewrite `doping.csv`.
```

- [ ] **Step 4: Tighten pn2d sample gate**

In `tests/regression/test_sentaurus_sample_integration.py`, require faithful IV to produce a finite Vela CSV:

```python
faithful_iv = out_dir / "vela" / "pn2d_iv_vela.csv"
self.assertTrue(faithful_iv.exists())
self.assertGreaterEqual(len(read_csv_rows(faithful_iv)), 2)
self.assert_csv_has_finite_currents(faithful_iv)
```

- [ ] **Step 5: Retire required runtime scaling from pn2d config**

In `reference_tcad/pn2d/pn2d_reference.json`, remove `runtime_doping_scale` as a required path for IV. Keep runtime approximation only as an optional diagnostic block:

```json
"runtime_diagnostic": {
  "enabled": true,
  "doping_scale": 0.0001,
  "step": 0.1
}
```

Update `scripts/sentaurus_import.py` only if needed to support the renamed optional block while preserving backward compatibility with existing `runtime_doping_scale`.

- [ ] **Step 6: Full verification and commit**

Run:

```powershell
cmake --build build --parallel
ctest --test-dir build --output-on-failure
python scripts\sentaurus_import.py reference --config reference_tcad\pn2d\pn2d_reference.json --source-dir reference_tcad\pn2d --output-dir build\pn2d_hybrid_faithful --tdr-importer build\sentaurus_import.exe --runner build\vela_example_runner.exe
```

Expected:
- full CTest passes;
- faithful pn2d IV produces finite output;
- runtime diagnostic remains optional;
- manifest records hybrid solver settings and unsupported Sentaurus physics.

Commit:

```powershell
git add src include tests scripts reference_tcad/pn2d docs/validation/pn2d_sentaurus_comparison.md
git commit -m "Converge pn2d with hybrid Gummel Newton solver"
```

---

## Suggested Task Prompts

Use these prompts one at a time:

1. `请执行 Gummel Newton Handoff 计划的任务1：添加 gummel_newton 方法解析和 DCSweep 诊断字段，完成后运行 dc_sweep 测试并提交。`
2. `请执行 Gummel Newton Handoff 计划的任务2：实现每个偏置点先 Gummel 初始化、再 Newton warm_start 接管，完成后运行 dc_sweep/newton/dd_gummel 测试并提交。`
3. `请执行 Gummel Newton Handoff 计划的任务3：增加 handoff fallback/strict 策略和配置文档，完成后测试并提交。`
4. `请执行 Gummel Newton Handoff 计划的任务4：让 pn2d Sentaurus 导入生成 faithful gummel_newton deck，同时保留 runtime 诊断 deck，完成后测试并提交。`
5. `请执行 Gummel Newton Handoff 计划的任务5：添加真实 pn2d faithful hybrid deck 生成门控和状态文档，完成后跑完整回归并提交。`
6. `请执行 Gummel Newton Handoff 计划的任务6：定位并修复真实高掺杂 node_doping pn2d 的 faithful hybrid 收敛问题，逐步移除 runtime_doping_scale 作为必要路径，完成后完整验证并提交。`

---

## Self-Review

- Spec coverage: The plan covers explicit Gummel initialization, Newton full-coupled handoff, pn2d generated decks, strict/fallback policies, faithful node-doping execution, and regression/reporting gates.
- Placeholder scan: No task depends on an undefined future component; each task names files, test commands, expected failures, and implementation shape.
- Type consistency: `gummel_newton`, `handoff.fallback`, `gummelIterations`, `newtonIterations`, `handoffStage`, and `node_doping_file` are used consistently across C++, JSON, docs, and tests.
