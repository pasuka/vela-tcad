# PN2D BV Avalanche Branch Continuation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use `superpowers:subagent-driven-development` or `superpowers:executing-plans` before implementing this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the PN2D Sentaurus 2018 coarse 7x3 BV reverse sweep reach Vela's avalanche-multiplication branch at `-20 V`, while preserving the existing default low-bias and non-breakdown behavior.

**Architecture:** First add a default-off external-state Newton solve entry to prove the multiplication branch is reachable by Vela's own DD+avalanche discretization. If that gate passes, add a default-off branch-reaching continuation mode, with generation/avalanche-source homotopy as the recommended first production path and pseudo-arclength/current-boundary continuation retained as fallback options.

**Tech Stack:** C++20, CMake/Ninja, MSYS2 UCRT64 on Windows, Catch2, `vela_example_runner`, existing PN2D Sentaurus import artifacts, Python comparison scripts under `scripts/`.

---

## Authoritative Premise

Use the root-cause result from `/memories/repo/pn2d_iv_root_cause.md`, especially entry 97, as the controlling premise for this work:

- The missing `-19 V`/`-20 V` Sentaurus BV knee is not caused by ionization-coefficient parameters, electric-field reconstruction, mesh resolution, nonlocal avalanche effects, quasi-Fermi update caps, residual scaling, conditioning, tolerances, SRH, or thermal-generation seeding.
- Feeding the Sentaurus `-20 V` multiplication state into Vela gives electron and hole continuity residuals near `1.6e-12`, effectively zero.
- Therefore the multiplication branch is a real stable fixed point of Vela's discretized DD+avalanche equations. The blocker is branch selection in the bistable `G = alpha |J|` system: ordinary voltage stepping remains on the non-multiplication branch and does not cross the fold.

If future experiments contradict this premise, stop the implementation path, record the falsification explicitly, and do not keep tuning continuation on a false assumption.

Known conflicting historical note: current repository docs contain older/parallel BV notes that classify the issue as a physics-magnitude gap. For this plan, entry 97 supersedes those notes unless a fresh reproduction disproves it.

---

## Current Repo Touch Points

- `src/tools/vela_example_runner.cpp`
  - Already has `readExternalState(...)`, which reads `state_fields_dir` with `ElectrostaticPotential_region0.csv`, `eQuasiFermiPotential_region0.csv`, and `hQuasiFermiPotential_region0.csv`.
  - Currently uses that reader for probe-style simulation types, not full Newton solving.
- `include/vela/solver/NewtonSolver.h` and `src/solver/NewtonSolver.cpp`
  - Already expose `NewtonSolver::solve(const DDSolution& initial)`.
  - This is the core API for the external-state Newton gate.
- `src/simulation/DCSweep.cpp`
  - Already supports `initial_state_file`, `write_state_file`, `write_state_every_point_prefix`, VTK output, branch-acceptance diagnostics, and default-off pseudo-arclength BV continuation.
  - `initial_state_file` is Vela restart CSV format, not Sentaurus `state_fields_dir` format.
- `src/io/DDSolutionCsv.cpp`
  - Owns Vela restart CSV format: `node_id,psi,phin,phip,electrons_m3,holes_m3`.
- `scripts/compare_reference_curves.py`
  - Use for curve-level log-current comparison against `reference_curves/pn2d_sentaurus2018_bv_reference.csv`.
- `scripts/pn2d_bv_branch_discriminator.py`
  - Use for node7 density and `Jn_x` multiplication-profile checks.

---

## Phase 0: Reproduce Baseline And Locate Inputs

**Purpose:** Freeze the exact starting behavior and artifacts before adding any new branch-reaching path.

- [ ] Confirm Windows toolchain setup:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
Set-Location "D:\code-repo\vela-tcad"
cmake --preset windows-ucrt64-debug
cmake --build --preset windows-ucrt64-debug --parallel
```

- [ ] Run the existing clean BV deck with all new features disabled. Prefer the current coarse7x3 imported deck if present: `build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/imported_reference/vela/simulation_bv.json`.

Expected: the run reaches `-20 V` and remains on the flat leakage branch.

```powershell
$BaselineConfig = Resolve-Path "build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\imported_reference\vela\simulation_bv.json"
build\windows-ucrt64-debug\vela_example_runner.exe --config $BaselineConfig
```

- [ ] Save the baseline CSV/log outside version control, under `build/` or another ignored output directory.

- [ ] Verify the Sentaurus multiplication state directory exists and contains:

```text
ElectrostaticPotential_region0.csv
eQuasiFermiPotential_region0.csv
hQuasiFermiPotential_region0.csv
```

- [ ] Record exact paths for the baseline config, baseline CSV, Sentaurus `state_fields_dir`, and Sentaurus BV reference curve in the validation log.

---

## Phase 1: Add External-State Newton Solve Gate

**Decision:** Add a new gated runner type, recommended name `newton_solve_from_state`, instead of overloading `dc_sweep.initial_state_file`.

**Rationale:** `state_fields_dir` is a Sentaurus-export-style set of field CSVs. `initial_state_file` is a Vela restart CSV. Keeping them separate avoids format ambiguity and keeps default sweep behavior unchanged.

### Task 1.1: Add Full Newton Solve From `state_fields_dir`

**Files:**

- Modify: `src/tools/vela_example_runner.cpp`
- Test: `tests/test_newton_solver.cpp`; create `tests/test_external_state.cpp` and add it to `CMakeLists.txt` if the runner helper cannot be exercised cleanly from existing tests
- Docs: `docs/config_schema.md`

- [ ] Add a helper near existing probe helpers:

```cpp
nlohmann::json runNewtonSolveFromState(const std::string& configFile,
                                       const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    NewtonProblem problem = loadNewtonProblem(configFile, cfg);
    const vela::DDSolution initial =
        readExternalState(cfgDir, cfg, problem.mesh.numNodes());
    const vela::NewtonSolver solver(
        problem.mesh, problem.matdb, problem.doping, problem.biases, problem.newton);
    vela::NewtonResult result = solver.solve(initial);

    if (cfg.contains("output_state_file")) {
        vela::writeDDSolutionStateCsv(
            resolvePath(cfgDir, cfg.at("output_state_file").get<std::string>()),
            result.solution);
    }
    if (cfg.contains("output_vtk")) {
        vela::writeDDSolutionVTK(
            resolvePath(cfgDir, cfg.at("output_vtk").get<std::string>()),
            problem.mesh,
            problem.doping,
            result.solution);
    }

    return {
        {"nodes", problem.mesh.numNodes()},
        {"converged", result.converged},
        {"iterations", result.iters},
        {"initial_residual", result.initialResidualNorm},
        {"final_residual", result.finalResidualNorm},
    };
}
```

Use the same `NewtonResult` field names already consumed by `runNewtonConfig` in this file: `converged`, `iters`, `initialResidualNorm`, `finalResidualNorm`, and `solution`. Use the same `writeDDSolutionVTK` overload already used by `runNewtonConfig`.

- [ ] Register the new type in `main`:

```cpp
} else if (type == "newton_solve_from_state") {
    status.update(runNewtonSolveFromState(configFile, cfg));
```

- [ ] Include `"converged"` in the returned status so the process exit code fails on Newton failure.

- [ ] Add schema docs for:

```json
{
  "simulation_type": "newton_solve_from_state",
  "state_fields_dir": "path/to/state_fields",
  "output_state_file": "outputs/minus20_from_state.csv",
  "output_vtk": "outputs/minus20_from_state.vtk"
}
```

### Task 1.2: Add Tests For The Gate

**Files:**

- Modify: `tests/test_newton_solver.cpp`, `tests/test_csv_utils.cpp`, or create `tests/test_external_state.cpp`
- Modify: `CMakeLists.txt` only if a new test target is created

- [ ] Test that a tiny external field directory with all required rows is read and used as a Newton initial state.

Minimum behavior to assert:

- all node ids are present exactly once;
- the solver starts from the provided `psi`, `phin`, and `phip`;
- a converged external initial state remains converged or takes only a small number of Newton iterations.

- [ ] Test malformed input:

```text
missing eQuasiFermiPotential_region0.csv
missing node_id column
missing component0 column
duplicate node_id
out-of-range node_id
missing node row
```

Each case should throw a clear runtime error naming the bad file or column.

- [ ] Run the focused tests:

```powershell
ctest --preset windows-ucrt64-debug -R "newton|csv|external" --output-on-failure
```

Expected: PASS.

### Task 1.3: Run The Decisive `-20 V` Gate

**Files:**

- Create generated config under ignored `build/pn2d_bv_external_state_gate/`
- Do not commit generated CSV/VTK/log outputs

- [ ] Create a config equivalent to the clean `simulation_bv.json`, but with:

```json
{
  "simulation_type": "newton_solve_from_state",
  "state_fields_dir": "path/to/sentaurus_minus20_state_fields",
  "contacts": [
    {"name": "Anode", "bias": -20.0},
    {"name": "Cathode", "bias": 0.0}
  ],
  "output_state_file": "outputs/pn2d_minus20_from_sentaurus_state.csv",
  "output_vtk": "outputs/pn2d_minus20_from_sentaurus_state.vtk"
}
```

- [ ] Run:

```powershell
build\windows-ucrt64-debug\vela_example_runner.exe --config build\pn2d_bv_external_state_gate\minus20_from_state.json
```

- [ ] Gate result:

Pass if all are true:

- Newton converges;
- final electron/hole continuity residuals remain near the entry-97 residual scale;
- node7 electron density stays in the multiplication-state scale, approximately `1e3` to `1e4 cm^-3`;
- terminal total current is on the order of `1e-7 A/um`, not `1e-14 A/um`;
- `pn2d_bv_branch_discriminator.py` classifies the state as multiplication-like.

Fail if Newton collapses to `n7 ~= 0` or leakage-level terminal current. On failure, stop and update the root-cause note before attempting continuation.

---

## Phase 2: Implement Branch-Reaching Continuation

Start Phase 2 only if Phase 1 proves Vela can converge and remain on the multiplication branch from the Sentaurus `-20 V` initial state.

### Recommended Option: Generation / Avalanche-Source Homotopy

**Decision:** Implement this first.

**Rationale:** It directly addresses branch selection in a bistable `G = alpha |J|` system, is default-off, can be confined to BV scans, and avoids relying on pseudo-arclength progress through a high-dimensional state norm that has already shown shallow-bias step collapse in this repo.

Concept:

- introduce a homotopy parameter `lambda`;
- at large `lambda`, add enough generation or avalanche-source boost to force the multiplication branch;
- solve the deep reverse-bias state;
- reduce `lambda -> 0` while warm-starting from the prior multiplication solution;
- use the final `lambda = 0` solution as the real Vela physical state.

Candidate config surface:

```json
{
  "sweep": {
    "continuation": {
      "avalanche_branch_homotopy": {
        "enabled": true,
        "mode": "extra_generation",
        "start_lambda": 1.0,
        "stop_lambda": 0.0,
        "steps": 12,
        "activation_bias_V": -18.0,
        "extra_generation_m3_s": 1.0e30,
        "write_homotopy_csv": "outputs/pn2d_bv_homotopy.csv"
      }
    }
  }
}
```

Implementation notes:

- Keep this config default-off.
- Restrict it to explicit opt-in BV continuation.
- Do not change default impact-ionization physics or low-bias solver behavior.
- Prefer an additive generation homotopy over permanent coefficient retuning.
- Log every accepted `(bias, lambda, total_current, node7_n, residual)` pair.

### Fallback Option A: Current-Boundary Continuation

Use only if generation/source homotopy cannot robustly return to `lambda = 0`.

Concept:

- switch from voltage control to current or mixed V/I control near the knee;
- cross the fold in current space;
- recover the associated voltage and final DD state.

Trade-off:

- physically clean for fold crossing;
- larger boundary-condition and terminal-equation surface area than generation homotopy.

### Fallback Option B: Pseudo-Arclength Retuning

Use only if the homotopy path fails or if the project wants one generic fold-crossing engine.

Current repo status:

- `include/vela/simulation/PseudoArclength.h` exists;
- `NewtonSolver::makeArclengthSystem` exists;
- `DCSweep` already has default-off `sweep.continuation.arclength`;
- historical runs progressed only to shallow reverse bias and suffered step-size collapse.

If resumed, first change the state weighting/progress metric or add a fixed-bias corrector fallback so accepted arclength steps produce useful bias progress.

---

## Phase 3: Verification Matrix

All four gates must pass before claiming the task complete.

### Gate A: Baseline Unchanged

- Run the original baseline with all new features disabled.
- Expected:
  - reaches `-20 V`;
  - current remains flat/leakage-like;
  - CSV columns and default config semantics are unchanged;
  - no new output files are written unless explicitly configured.

### Gate B: Branch Reached

- Run the new feature-enabled BV flow to `-20 V`.
- Compare against `reference_curves/pn2d_sentaurus2018_bv_reference.csv`:

```powershell
python scripts\compare_reference_curves.py `
  --reference build-release\reference_tcad\pn2d_sentaurus2018\reference_curves\pn2d_sentaurus2018_bv_reference.csv `
  --candidate build\pn2d_bv_branch_continuation\pn2d_bv_branch.csv `
  --output-json build\pn2d_bv_branch_continuation\compare_bv.json `
  --output-md build\pn2d_bv_branch_continuation\compare_bv.md `
  --candidate-column current_total_A_per_um `
  --reference-column current_total `
  --candidate-scale -1.0 `
  --bias-min -20.0 `
  --bias-max -18.0 `
  --interpolation log_current `
  --min-points 3
```

Expected:

- `-19 V` and `-20 V` current growth markers appear;
- max absolute log10 current error over `-18 V` to `-20 V` decreases significantly from the baseline;
- `-20 V` terminal current is in the Sentaurus multiplication order of magnitude.

### Gate C: Physical Self-Consistency

Run:

```powershell
python scripts\pn2d_bv_branch_discriminator.py `
  --csv build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\reports\coarse_previous_full20_vector_current_20260630\coarse_node_field_compare_aligned.csv `
  --bias -20.0
```

Expected:

- node7 electron density remains multiplication-like;
- `Jn_x` along the active x direction shows qualitative monotonic multiplication consistent with Sentaurus;
- no isolated numerical spike is responsible for the terminal current.

### Gate D: Regression

Run focused then full tests:

```powershell
ctest --preset windows-ucrt64-debug -R "newton|dc_sweep|impact|pseudo_arclength" --output-on-failure
ctest --preset windows-ucrt64-debug --output-on-failure
```

Expected: all project tests pass. If an unrelated failure is already documented in the current workspace before this work starts, quote that exact failure and rerun enough focused tests to prove the new code path is not involved.

---

## Documentation And Memory Updates

- [ ] Update `docs/config_schema.md` for every new config key.
- [ ] Append validation results to `docs/validation/pn2d_bv_validation.md`.
- [ ] Append a concise entry to `/memories/repo/pn2d_iv_root_cause.md` after each completed phase:

```text
2026-07-01: Phase 1 external-state Newton gate. Config path: build/pn2d_bv_external_state_gate/minus20_from_state.json. Log path: build/pn2d_bv_external_state_gate/minus20_from_state.log. Outcome: pass/fail with final residual, node7 electron density, and terminal current.
```

- [ ] If a hypothesis is falsified, write the falsification in the same note and remove it from active assumptions.
- [ ] Do not commit generated simulation outputs unless explicitly requested.

---

## Stop Conditions

Stop and report instead of continuing if any condition occurs:

- external-state Newton from Sentaurus `-20 V` collapses to the leakage branch;
- the homotopy can reach multiplication only with nonzero artificial generation but cannot return to `lambda = 0`;
- feature-enabled runs improve current magnitude but fail the branch discriminator;
- default-off behavior changes baseline low-bias, IV, or non-breakdown examples;
- full regression failures point to solver/physics behavior changed outside the gated path.

---

## Recommended Execution Order

1. Implement `newton_solve_from_state`.
2. Add focused input-format and solve-path tests.
3. Run the Sentaurus `-20 V` external-state Newton gate.
4. If the gate passes, implement generation/source homotopy as the first branch-reaching continuation path.
5. Run A-D verification.
6. Update docs and memory with exact logs and conclusions.




