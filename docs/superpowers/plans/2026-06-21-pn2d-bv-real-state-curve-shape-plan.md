# PN2D BV Real-State Curve-Shape Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move PN2D BV work from fixture-level Jacobian confidence to real-state evidence, then decide whether the remaining Vela/Sentaurus curve-shape gap is a missing Jacobian term, a state-alignment issue, or an accepted discretization/model-parity limit.

**Architecture:** Preserve the current SG edge-current avalanche production path and avoid core SG flux-divergence rewrites. Reuse the existing restart CSV format as the single state handoff format, add real BV-state Jacobian block auditing through `NewtonSolver`, and gate any physics/config change by full `-10 V..-20 V` knee-shape movement rather than by one bias point.

**Tech Stack:** C++20, Catch2, CMake/Ninja, MSYS2 UCRT64, Python standard library diagnostics, existing `DCSweep`, `NewtonSolver`, `vela_example_runner`, and `build-release/reference_tcad/pn2d_sentaurus2018` artifacts.

---

## Current Baseline

Use these facts as the starting point for the work:

- The `avaljac` sweep reaches `-20 V` with `1136` converged points and no Newton failure.
- The current SG avalanche source Jacobian includes endpoint-density, alpha driving-field, and source mobility sensitivity on the `density_gradient/current_density` path.
- The fixture C++ Jacobian block probe reports finite analytic-vs-FD agreement, but it does not replay the full `pn2d_sentaurus2018` BV mesh/state.
- The knee-shape script reports Sentaurus first 1 V growth ratio `> 1.5` at `-19.0 V` and `> 2.0` at `-20.0 V`; Vela `avaljac` has no matching threshold in `-10..-20 V`, with max absolute log10 current error `0.891693` decades.
- Therefore convergence is solved for this branch; curve shape and real-state physics parity are still open.

## Files And Responsibilities

- Create `include/vela/io/DDSolutionCsv.h`: shared declarations for restart-state CSV read/write.
- Create `src/io/DDSolutionCsv.cpp`: shared implementation moved from `src/simulation/DCSweep.cpp`.
- Modify `src/simulation/DCSweep.cpp`: use shared restart-state I/O and optionally write one state CSV per accepted sweep point.
- Modify `include/vela/simulation/DCSweep.h`: add `writeStateEveryPointPrefix`.
- Modify `include/vela/solver/NewtonSolver.h`: expose a real-state Jacobian block-audit result API.
- Modify `src/solver/NewtonSolver.cpp`: compute analytic/FD Jacobian block norms on a provided solved or near-solved state.
- Modify `src/tools/vela_example_runner.cpp`: add `simulation_type: "newton_jacobian_block_probe"` using the shared restart CSV state.
- Modify `tests/test_dc_sweep.cpp`: cover per-point state CSV export.
- Modify `tests/test_newton_solver.cpp`: cover block audit rows on a small coupled-DD state.
- Modify `tests/regression/test_reference_tcad_tools.py`: cover runner config generation and CSV contract for the real-state probe.
- Modify `scripts/diagnose_pn2d_bv_jacobian_block_audit.py`: delegate to the runner real-state probe when state/config pairs are supplied; keep the fixture probe fallback.
- Create `scripts/run_pn2d_bv_real_state_jacobian_audit.py`: derive per-bias restart states from the `avaljac` deck, run the real-state probe, and write a combined audit report.
- Modify `scripts/diagnose_pn2d_bv_knee_shape.py`: optionally emit JSON for automated gates.
- Modify `docs/validation/pn2d_sentaurus_comparison.md`: resolve the current documentation conflict and record real-state evidence.

---

### Task 1: Resolve The BV Documentation Contract

**Files:**
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`

- [ ] **Step 1: Edit the BV acceptance wording**

Replace the conclusion in `docs/validation/pn2d_sentaurus_comparison.md` that says the high-bias BV gap is accepted with this narrower statement:

```markdown
### BV Acceptance Scope After `avaljac`

The `avaljac` branch demonstrates that the Sentaurus-faithful BV physics block can
be continued to `-20 V` without Newton failure. This is a convergence milestone,
not a final BV parity acceptance. The full-curve shape gate remains open because
the current Vela curve does not reproduce the Sentaurus one-volt growth knee in
the `-18 V..-20 V` region.

Accepted status is limited to:

- SG avalanche source Jacobian completeness for the current production path.
- SRH/Auger local derivative coverage in the coupled residual/Jacobian.
- End-to-end continuation robustness to `-20 V` for the current branch.

Open status remains:

- Real full-mesh BV-state Jacobian block replay.
- Curve-shape parity over `-10 V..-20 V`.
- The physical cause of the missing high-bias one-volt current-growth knee.
```

- [ ] **Step 2: Verify the document contains a single final BV status**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
Select-String -Path docs\validation\pn2d_sentaurus_comparison.md -Pattern "BV Acceptance Scope|PN2D BV Knee-Shape Acceptance Gate|irreducible difference|therefore accepted"
```

Expected: the old broad acceptance language is either removed or explicitly scoped below the new `BV Acceptance Scope After avaljac` heading.

- [ ] **Step 3: Commit**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git add docs/validation/pn2d_sentaurus_comparison.md
git commit -m "Clarify PN2D BV avaljac acceptance scope"
```

---

### Task 2: Extract Restart-State CSV I/O Into A Shared Module

**Files:**
- Create: `include/vela/io/DDSolutionCsv.h`
- Create: `src/io/DDSolutionCsv.cpp`
- Modify: `src/simulation/DCSweep.cpp`
- Modify: `CMakeLists.txt`
- Test: `tests/test_dc_sweep.cpp`

- [ ] **Step 1: Add the failing shared-I/O include test**

Add to `tests/test_dc_sweep.cpp` near the existing `write_state_file` tests:

```cpp
TEST_CASE("DDSolution CSV shared IO roundtrips restart state", "[dc_sweep]")
{
    const auto dir = std::filesystem::temp_directory_path() / "vela_ddsolution_csv_roundtrip";
    std::filesystem::create_directories(dir);
    const auto path = dir / "state.csv";

    vela::DDSolution solution;
    solution.psi = vela::VectorXd::LinSpaced(3, -0.1, 0.1);
    solution.phin = vela::VectorXd::LinSpaced(3, 0.2, 0.4);
    solution.phip = vela::VectorXd::LinSpaced(3, -0.4, -0.2);
    solution.n = vela::VectorXd::Constant(3, 1.0e16);
    solution.p = vela::VectorXd::Constant(3, 2.0e16);

    vela::writeDDSolutionStateCsv(path, solution);
    const vela::DDSolution loaded = vela::readDDSolutionStateCsv(path, 3);

    REQUIRE((loaded.psi - solution.psi).norm() == Catch::Approx(0.0));
    REQUIRE((loaded.phin - solution.phin).norm() == Catch::Approx(0.0));
    REQUIRE((loaded.phip - solution.phip).norm() == Catch::Approx(0.0));
    REQUIRE((loaded.n - solution.n).norm() == Catch::Approx(0.0));
    REQUIRE((loaded.p - solution.p).norm() == Catch::Approx(0.0));
}
```

- [ ] **Step 2: Run the test to confirm the missing header failure**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --target test_dc_sweep --parallel
```

Expected: build fails because `vela/io/DDSolutionCsv.h` does not exist.

- [ ] **Step 3: Create the shared header**

Create `include/vela/io/DDSolutionCsv.h`:

```cpp
#pragma once

#include "vela/core/Types.h"
#include "vela/solver/GummelSolver.h"

#include <filesystem>

namespace vela {

DDSolution readDDSolutionStateCsv(const std::filesystem::path& path,
                                  Index expectedNodeCount);

void writeDDSolutionStateCsv(const std::filesystem::path& path,
                             const DDSolution& solution);

} // namespace vela
```

- [ ] **Step 4: Move the implementation**

Create `src/io/DDSolutionCsv.cpp` by moving the existing `parseRestartStateReal`, `parseRestartStateNodeId`, `readDDSolutionStateCsv`, and `writeDDSolutionStateCsv` logic out of `src/simulation/DCSweep.cpp`. Keep the current CSV header exactly:

```text
node_id,psi,phin,phip,electrons_m3,holes_m3
```

Keep the current error text prefixes:

```text
DCSweep: initial_state_file
DCSweep: cannot write restart state
DCSweep: cannot open write_state_file
```

- [ ] **Step 5: Wire the module into the build**

Add `src/io/DDSolutionCsv.cpp` to the `vela_core` source list in `CMakeLists.txt`.

- [ ] **Step 6: Update DCSweep to use the shared header**

In `src/simulation/DCSweep.cpp`, include:

```cpp
#include "vela/io/DDSolutionCsv.h"
```

Delete the moved local function definitions from `DCSweep.cpp`. Leave the call sites at `readDDSolutionStateCsv(...)` and `writeDDSolutionStateCsv(...)` unchanged.

- [ ] **Step 7: Run focused tests**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --target test_dc_sweep --parallel
.\build-release\test_dc_sweep.exe "DDSolution CSV shared IO roundtrips restart state"
.\build-release\test_dc_sweep.exe "DCSweep: write_state_file stores latest converged restart state"
.\build-release\test_dc_sweep.exe "DCSweep: initial_state_file validates restart node coverage"
```

Expected: all three commands pass.

- [ ] **Step 8: Commit**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git add CMakeLists.txt include/vela/io/DDSolutionCsv.h src/io/DDSolutionCsv.cpp src/simulation/DCSweep.cpp tests/test_dc_sweep.cpp
git commit -m "Share DD restart state CSV IO"
```

---

### Task 3: Add Per-Accepted-Point State Snapshots To DCSweep

**Files:**
- Modify: `include/vela/simulation/DCSweep.h`
- Modify: `src/simulation/DCSweep.cpp`
- Test: `tests/test_dc_sweep.cpp`

- [ ] **Step 1: Add the failing per-point snapshot test**

Add to `tests/test_dc_sweep.cpp`:

```cpp
TEST_CASE("DCSweep: write_state_every_point_prefix stores accepted states", "[dc_sweep]")
{
    const auto dir = std::filesystem::temp_directory_path() / "vela_dc_sweep_point_states";
    std::filesystem::create_directories(dir);
    const auto configPath = dir / "sweep.json";
    const auto prefix = dir / "states" / "bv_state";

    writeSimpleDCSweepConfig(configPath, {
        {"start", 0.0},
        {"stop", -0.1},
        {"step", -0.05},
        {"write_state_every_point_prefix", prefix.string()}
    });

    vela::DCSweep sweep;
    const auto result = sweep.runWithResult(configPath.string());

    REQUIRE(result.points.size() == 3);
    REQUIRE(std::filesystem::exists(dir / "states" / "bv_state_bias_0p000000.csv"));
    REQUIRE(std::filesystem::exists(dir / "states" / "bv_state_bias_m0p050000.csv"));
    REQUIRE(std::filesystem::exists(dir / "states" / "bv_state_bias_m0p100000.csv"));
}
```

If the local helper is not named `writeSimpleDCSweepConfig`, use the existing helper in `tests/test_dc_sweep.cpp` that writes the minimal coupled-DD sweep JSON. Keep the expected file names unchanged.

- [ ] **Step 2: Run the test to see the missing config field failure**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --target test_dc_sweep --parallel
.\build-release\test_dc_sweep.exe "write_state_every_point_prefix"
```

Expected: fail because `write_state_every_point_prefix` is ignored or not parsed.

- [ ] **Step 3: Extend the config struct**

In `include/vela/simulation/DCSweep.h`, add:

```cpp
std::string writeStateEveryPointPrefix;
```

next to `writeStateFile`.

- [ ] **Step 4: Parse and resolve the new field**

In `src/simulation/DCSweep.cpp`, update `dcSweepConfigFromJson`:

```cpp
sweep.writeStateEveryPointPrefix =
    j.value("write_state_every_point_prefix", std::string{});
```

In the existing path-resolution block:

```cpp
if (!sweep.writeStateEveryPointPrefix.empty())
    sweep.writeStateEveryPointPrefix = resolve(sweep.writeStateEveryPointPrefix);
```

- [ ] **Step 5: Add deterministic bias-file naming**

Add a small helper in `src/simulation/DCSweep.cpp`:

```cpp
std::string biasToken(Real bias)
{
    std::ostringstream out;
    out << std::fixed << std::setprecision(6) << std::abs(bias);
    std::string token = out.str();
    std::replace(token.begin(), token.end(), '.', 'p');
    return (bias < 0.0 ? "m" : "") + token;
}
```

Include `<sstream>` if it is not already present.

- [ ] **Step 6: Write one state per converged accepted point**

At the same location that currently writes `sweep.writeStateFile` for a converged point, add:

```cpp
if (converged && !sweep.writeStateEveryPointPrefix.empty()) {
    const std::filesystem::path prefix(sweep.writeStateEveryPointPrefix);
    const std::filesystem::path path =
        prefix.parent_path() /
        (prefix.filename().string() + "_bias_" + biasToken(bias) + ".csv");
    writeDDSolutionStateCsv(path, sol);
}
```

Use the same `bias` variable that is recorded into the accepted `DCSweepPoint`.

- [ ] **Step 7: Run focused tests**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --target test_dc_sweep --parallel
.\build-release\test_dc_sweep.exe "write_state_every_point_prefix"
.\build-release\test_dc_sweep.exe "write_state_file stores latest"
```

Expected: both tests pass.

- [ ] **Step 8: Commit**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git add include/vela/simulation/DCSweep.h src/simulation/DCSweep.cpp tests/test_dc_sweep.cpp
git commit -m "Write DC sweep state snapshots per accepted point"
```

---

### Task 4: Add A Real-State Jacobian Block API To NewtonSolver

**Files:**
- Modify: `include/vela/solver/NewtonSolver.h`
- Modify: `src/solver/NewtonSolver.cpp`
- Test: `tests/test_newton_solver.cpp`

- [ ] **Step 1: Add the failing solver API test**

Add to `tests/test_newton_solver.cpp` near other coupled Newton diagnostic tests:

```cpp
TEST_CASE("NewtonSolver evaluates real-state Jacobian block audit rows", "[newton][coupled]")
{
    auto fixture = makeCoupledNewtonFixture();
    vela::NewtonSolver solver(
        fixture.mesh,
        fixture.matdb,
        fixture.doping,
        fixture.contactBiases,
        fixture.config);

    const vela::NewtonResult solved = solver.solve();
    REQUIRE(solved.converged);

    const auto rows = solver.evaluateJacobianBlockAudit(solved.solution, 1.0e-7);
    const auto hasBlock = [&](const std::string& name) {
        return std::any_of(rows.begin(), rows.end(), [&](const auto& row) {
            return row.block == name &&
                   std::isfinite(row.analyticNorm) &&
                   std::isfinite(row.fdNorm) &&
                   std::isfinite(row.relDiff);
        });
    };

    REQUIRE(hasBlock("poisson"));
    REQUIRE(hasBlock("transport"));
    REQUIRE(hasBlock("srh_auger"));
    REQUIRE(hasBlock("sg_avalanche"));
    REQUIRE(hasBlock("dirichlet_or_gauge"));
}
```

If `makeCoupledNewtonFixture()` does not exist, use the smallest existing helper in `tests/test_newton_solver.cpp` that creates a converged coupled-DD PN fixture. Keep the test assertion on the five block names.

- [ ] **Step 2: Run the test to confirm the missing API**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --target test_newton_solver --parallel
```

Expected: build fails because `evaluateJacobianBlockAudit` is not declared.

- [ ] **Step 3: Add public result types**

In `include/vela/solver/NewtonSolver.h`, add:

```cpp
struct NewtonJacobianBlockAuditRow {
    std::string block;
    Real analyticNorm = 0.0;
    Real fdNorm = 0.0;
    Real diffNorm = 0.0;
    Real relDiff = 0.0;
};
```

Add this public method:

```cpp
std::vector<NewtonJacobianBlockAuditRow> evaluateJacobianBlockAudit(
    const DDSolution& state,
    Real finiteDifferenceStep = 1.0e-7) const;
```

- [ ] **Step 4: Implement with the same scaling path as residual evaluation**

In `src/solver/NewtonSolver.cpp`, implement the method using the same assembler construction, scaling, `pack`, and boundary-condition path used by `evaluateResidual` and `evaluateStep`.

The method must:

```text
1. Build CoupledDDAssembler from mesh_, matdb_, doping_, cfg_.mobility, cfg_.recombination, cfg_.bandgapNarrowing, cfg_.impactIonization.
2. Build boundary conditions with buildBoundaryConditions(assembler).
3. Pack the physical DDSolution into the scaled CoupledDDState in the same way as evaluateResidual.
4. Compute assembler.assembleJacobian(x, bcs).
5. Compute assembler.finiteDifferenceJacobian(x, bcs, finiteDifferenceStep).
6. Return restricted-row norms for poisson, transport, srh_auger, sg_avalanche, and dirichlet_or_gauge.
```

Use row groups:

```text
poisson: rows [0, N)
transport: rows [N, 3N)
srh_auger: rows [N, 3N), with recombination enabled minus recombination none
sg_avalanche: rows [N, 3N), with impact enabled minus impact none
dirichlet_or_gauge: rows where boundary rows are imposed by CoupledDDBoundaryConditions
```

For `srh_auger` and `sg_avalanche`, construct comparison assemblers that differ by exactly one physics block. Keep mobility, BGN, temperature, doping, mesh, and boundary conditions identical.

- [ ] **Step 5: Run solver tests**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --target test_newton_solver --parallel
.\build-release\test_newton_solver.exe "[newton][coupled]"
```

Expected: all coupled Newton tests pass.

- [ ] **Step 6: Commit**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git add include/vela/solver/NewtonSolver.h src/solver/NewtonSolver.cpp tests/test_newton_solver.cpp
git commit -m "Add real-state Newton Jacobian block audit"
```

---

### Task 5: Expose The Real-State Probe In vela_example_runner

**Files:**
- Modify: `src/tools/vela_example_runner.cpp`
- Test: `tests/regression/test_reference_tcad_tools.py`

- [ ] **Step 1: Add a regression test for the runner config contract**

Add to `tests/regression/test_reference_tcad_tools.py`:

```python
def test_runner_real_state_jacobian_block_probe_config_contract(self) -> None:
    config = {
        "simulation_type": "newton_jacobian_block_probe",
        "mesh_file": "mesh.json",
        "node_doping_file": "doping.csv",
        "contacts": [
            {"name": "Anode", "bias": -13.2},
            {"name": "Cathode", "bias": 0.0},
        ],
        "solver": {
            "method": "gummel_newton",
            "recombination": ["srh"],
            "impact_ionization": {
                "model": "van_overstraeten",
                "driving_force": "quasi_fermi_gradient",
                "generation": "current_density",
                "current_approximation": "density_gradient",
            },
        },
        "state_file": "states/bv_state_bias_m13p200000.csv",
        "output_csv": "reports/jacobian_blocks_m13p2.csv",
        "finite_difference_step": 1.0e-7,
    }
    self.assertEqual(config["simulation_type"], "newton_jacobian_block_probe")
    self.assertEqual(config["state_file"], "states/bv_state_bias_m13p200000.csv")
    self.assertEqual(config["finite_difference_step"], 1.0e-7)
```

- [ ] **Step 2: Add runner implementation**

In `src/tools/vela_example_runner.cpp`, include:

```cpp
#include "vela/io/DDSolutionCsv.h"
```

Add a `runNewtonJacobianBlockProbe` helper:

```cpp
nlohmann::json runNewtonJacobianBlockProbe(const std::string& configFile,
                                           const nlohmann::json& cfg)
{
    const std::filesystem::path cfgDir = configDirectory(configFile);
    NewtonProblem problem = loadNewtonProblem(configFile, cfg);
    const std::filesystem::path statePath =
        resolvePath(cfgDir, cfg.at("state_file").get<std::string>());
    const vela::DDSolution state =
        vela::readDDSolutionStateCsv(statePath, problem.mesh.numNodes());
    const vela::Real fdStep = cfg.value("finite_difference_step", 1.0e-7);

    const vela::NewtonSolver solver(
        problem.mesh, problem.matdb, problem.doping, problem.biases, problem.newton);
    const auto rows = solver.evaluateJacobianBlockAudit(state, fdStep);

    const std::filesystem::path outputPath =
        resolvePath(cfgDir, cfg.at("output_csv").get<std::string>());
    if (!outputPath.parent_path().empty())
        std::filesystem::create_directories(outputPath.parent_path());
    std::ofstream out(outputPath);
    if (!out.is_open())
        throw std::runtime_error("Cannot write jacobian block probe CSV: " + outputPath.string());
    out << "block,analytic_norm,fd_norm,diff_norm,rel_diff\n";
    for (const auto& row : rows) {
        out << row.block << ','
            << row.analyticNorm << ','
            << row.fdNorm << ','
            << row.diffNorm << ','
            << row.relDiff << '\n';
    }

    return {
        {"nodes", problem.mesh.numNodes()},
        {"blocks", rows.size()},
        {"output_csv", outputPath.string()},
    };
}
```

Register it in `main`:

```cpp
} else if (type == "newton_jacobian_block_probe") {
    status.update(runNewtonJacobianBlockProbe(configFile, cfg));
```

- [ ] **Step 3: Run regression and build checks**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --target vela_example_runner test_newton_solver --parallel
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_runner_real_state_jacobian_block_probe_config_contract
```

Expected: build succeeds and the Python regression passes.

- [ ] **Step 4: Commit**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git add src/tools/vela_example_runner.cpp tests/regression/test_reference_tcad_tools.py
git commit -m "Expose real-state Jacobian block probe"
```

---

### Task 6: Automate Real BV-State Jacobian Audits

**Files:**
- Create: `scripts/run_pn2d_bv_real_state_jacobian_audit.py`
- Modify: `scripts/diagnose_pn2d_bv_jacobian_block_audit.py`
- Test: `tests/regression/test_reference_tcad_tools.py`

- [ ] **Step 1: Add a Python unit test for combined rows**

Add to `tests/regression/test_reference_tcad_tools.py`:

```python
def test_pn2d_bv_real_state_jacobian_audit_combines_bias_rows(self) -> None:
    module_path = REPO / "scripts" / "run_pn2d_bv_real_state_jacobian_audit.py"
    spec = importlib.util.spec_from_file_location("run_pn2d_bv_real_state_jacobian_audit", module_path)
    self.assertIsNotNone(spec)
    self.assertIsNotNone(spec.loader)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    with tempfile.TemporaryDirectory(prefix="vela_real_state_jacobian_audit_") as tmp:
        root = Path(tmp)
        block_csv = root / "block.csv"
        self._write_csv(block_csv, [
            "block", "analytic_norm", "fd_norm", "diff_norm", "rel_diff"
        ], [["sg_avalanche", 2.0, 2.0, 1.0e-8, 5.0e-9]])
        out = root / "combined.csv"
        module.write_combined_report(out, [module.BiasBlockReport(-13.2, block_csv)])
        rows = self._read_csv(out)

    self.assertEqual(rows[0]["bias_V"], "-13.2")
    self.assertEqual(rows[0]["block"], "sg_avalanche")
    self.assertLess(float(rows[0]["rel_diff"]), 1.0e-6)
```

- [ ] **Step 2: Create the audit runner**

Create `scripts/run_pn2d_bv_real_state_jacobian_audit.py` with:

```python
#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import json
import subprocess
from dataclasses import dataclass
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
DEFAULT_BASE_CONFIG = REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018" / "vela" / "simulation_bv_minus20_avaljac.json"
DEFAULT_RUNNER = REPO / "build-release" / "vela_example_runner.exe"
DEFAULT_OUT_DIR = REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018" / "reports" / "bv_real_state_jacobian_audit"


@dataclass(frozen=True)
class BiasBlockReport:
    bias: float
    block_csv: Path


def bias_token(bias: float) -> str:
    token = f"{abs(bias):.6f}".replace(".", "p")
    return ("m" if bias < 0 else "") + token


def write_combined_report(output: Path, reports: list[BiasBlockReport]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "bias_V", "block", "analytic_norm", "fd_norm", "diff_norm", "rel_diff"
        ])
        writer.writeheader()
        for report in reports:
            with report.block_csv.open(newline="", encoding="utf-8") as block_handle:
                for row in csv.DictReader(block_handle):
                    row = dict(row)
                    row["bias_V"] = f"{report.bias:g}"
                    writer.writerow(row)


def derive_snapshot_config(base_config: Path, out_dir: Path, biases: list[float]) -> Path:
    cfg = json.loads(base_config.read_text(encoding="utf-8"))
    sweep = cfg.setdefault("sweep", {})
    sweep["bias_points"] = biases
    sweep["csv_file"] = str(out_dir / "snapshot_sweep.csv")
    sweep["write_state_every_point_prefix"] = str(out_dir / "states" / "bv_state")
    cfg_path = out_dir / "snapshot_sweep.json"
    cfg_path.parent.mkdir(parents=True, exist_ok=True)
    cfg_path.write_text(json.dumps(cfg, indent=2), encoding="utf-8")
    return cfg_path


def derive_probe_config(base_config: Path, out_dir: Path, bias: float) -> Path:
    cfg = json.loads(base_config.read_text(encoding="utf-8"))
    cfg["simulation_type"] = "newton_jacobian_block_probe"
    cfg["state_file"] = str(out_dir / "states" / f"bv_state_bias_{bias_token(bias)}.csv")
    cfg["output_csv"] = str(out_dir / "blocks" / f"jacobian_blocks_{bias_token(bias)}.csv")
    cfg["finite_difference_step"] = 1.0e-7
    cfg["contacts"][0]["bias"] = bias
    cfg_path = out_dir / "configs" / f"jacobian_probe_{bias_token(bias)}.json"
    cfg_path.parent.mkdir(parents=True, exist_ok=True)
    cfg_path.write_text(json.dumps(cfg, indent=2), encoding="utf-8")
    return cfg_path


def run_runner(runner: Path, config: Path) -> None:
    subprocess.run([str(runner), "--config", str(config)], check=True, cwd=REPO)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-config", type=Path, default=DEFAULT_BASE_CONFIG)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--bias", action="append", type=float, default=[-11.5, -13.2, -18.0, -19.0, -20.0])
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    snapshot_config = derive_snapshot_config(args.base_config, args.out_dir, args.bias)
    run_runner(args.runner, snapshot_config)

    reports: list[BiasBlockReport] = []
    for bias in args.bias:
        probe_config = derive_probe_config(args.base_config, args.out_dir, bias)
        run_runner(args.runner, probe_config)
        reports.append(BiasBlockReport(
            bias,
            args.out_dir / "blocks" / f"jacobian_blocks_{bias_token(bias)}.csv",
        ))
    write_combined_report(args.out_dir / "jacobian_blocks_real_state.csv", reports)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 3: Run the new unit test**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_real_state_jacobian_audit_combines_bias_rows
```

Expected: pass.

- [ ] **Step 4: Run the real audit**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts\run_pn2d_bv_real_state_jacobian_audit.py
```

Expected:

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_real_state_jacobian_audit\jacobian_blocks_real_state.csv
```

exists with five block rows for each requested bias.

- [ ] **Step 5: Acceptance check**

Run:

```powershell
$csv = "build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_real_state_jacobian_audit\jacobian_blocks_real_state.csv"
Import-Csv $csv | Format-Table bias_V,block,rel_diff
```

Expected initial acceptance:

```text
sg_avalanche rel_diff <= 5e-5
srh_auger rel_diff <= 1e-4
transport rel_diff <= 1e-4
poisson finite and near the fixture-level scale
dirichlet_or_gauge finite and near the fixture-level scale
```

If `sg_avalanche` fails above `5e-5` on the real state, stop curve-shape tuning and inspect the SG avalanche source Jacobian around the failing bias.

- [ ] **Step 6: Commit**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git add scripts/run_pn2d_bv_real_state_jacobian_audit.py scripts/diagnose_pn2d_bv_jacobian_block_audit.py tests/regression/test_reference_tcad_tools.py
git commit -m "Add PN2D BV real-state Jacobian audit runner"
```

---

### Task 7: Make The Knee-Shape Gate Machine-Readable

**Files:**
- Modify: `scripts/diagnose_pn2d_bv_knee_shape.py`
- Test: `tests/regression/test_reference_tcad_tools.py`

- [ ] **Step 1: Add a JSON-output regression**

Extend `test_pn2d_bv_knee_shape_computes_thresholds_and_log_error` with:

```python
summary = module.build_summary(cand_points, ref_points, -12.0, -10.0)
self.assertEqual(summary["candidate"]["first_growth_over_2p0"], -11.0)
self.assertEqual(summary["reference"]["first_growth_over_1p5"], -12.0)
self.assertGreater(summary["max_abs_log10_current_error"], 0.0)
```

- [ ] **Step 2: Add `build_summary`**

In `scripts/diagnose_pn2d_bv_knee_shape.py`, add:

```python
def build_summary(candidate: list[tuple[float, float]],
                  reference: list[tuple[float, float]],
                  bias_min: float,
                  bias_max: float) -> dict[str, object]:
    reference_summary = summarize_curve(reference, bias_min, bias_max)
    candidate_summary = summarize_curve(candidate, bias_min, bias_max)
    return {
        "bias_min": bias_min,
        "bias_max": bias_max,
        "reference": {
            "first_growth_over_1p5": reference_summary.first_growth_over_1p5,
            "first_growth_over_2p0": reference_summary.first_growth_over_2p0,
        },
        "candidate": {
            "first_growth_over_1p5": candidate_summary.first_growth_over_1p5,
            "first_growth_over_2p0": candidate_summary.first_growth_over_2p0,
        },
        "max_abs_log10_current_error": max_abs_log10_error(
            candidate, reference, bias_min, bias_max),
    }
```

Add CLI option:

```python
parser.add_argument("--output-json", type=Path)
```

After computing `summary`, write it when requested:

```python
if args.output_json is not None:
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")
```

Import `json`.

- [ ] **Step 3: Run regression**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_bv_knee_shape_computes_thresholds_and_log_error
```

Expected: pass.

- [ ] **Step 4: Generate current gate JSON**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts\diagnose_pn2d_bv_knee_shape.py `
  --output-json build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_real_state_jacobian_audit\knee_shape_current.json
```

Expected JSON contains:

```json
{
  "reference": {
    "first_growth_over_1p5": -19.0,
    "first_growth_over_2p0": -20.0
  },
  "candidate": {
    "first_growth_over_1p5": null,
    "first_growth_over_2p0": null
  }
}
```

- [ ] **Step 5: Commit**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git add scripts/diagnose_pn2d_bv_knee_shape.py tests/regression/test_reference_tcad_tools.py
git commit -m "Emit machine-readable PN2D BV knee gate"
```

---

### Task 8: Run The Decision Matrix

**Files:**
- Read: `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_real_state_jacobian_audit/jacobian_blocks_real_state.csv`
- Read: `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_real_state_jacobian_audit/knee_shape_current.json`
- Read: `docs/validation/pn2d_sentaurus_comparison.md`
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`

- [ ] **Step 1: Classify the real-state Jacobian result**

Use this table:

```text
Case A: sg_avalanche real-state rel_diff <= 5e-5
Action: Treat current SG avalanche Jacobian as complete enough; continue to curve-shape physics/state experiments.

Case B: sg_avalanche real-state rel_diff > 5e-5 and fixture rel_diff <= 5e-5
Action: Inspect full-state-only terms: boundary rows, state scaling, low-density interpolation, high-field mobility driving-force input, and field floors.

Case C: transport real-state rel_diff > 1e-4 while sg_avalanche passes
Action: Do not tune avalanche first; isolate high-field mobility derivative in ordinary SG transport flux.

Case D: srh_auger real-state rel_diff > 1e-4
Action: inspect `totalRateDerivativesFromExcessProduct` and BGN/effective-ni derivative inputs before changing avalanche.
```

- [ ] **Step 2: Run existing curve-shape diagnostics before changing code**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts\diagnose_pn2d_bv_active_edge_flux_factors.py --help
python scripts\diagnose_pn2d_bv_source_policy_matrix.py --help
python scripts\diagnose_pn2d_bv_loop_gain_sensitivity.py --help
```

Expected: all three scripts are present and their help exits with code `0`.

- [ ] **Step 3: Choose exactly one next experiment**

Pick the first matching branch:

```text
1. Missing real-state Jacobian block: write a focused C++ test reproducing that block on a compact mesh, then fix the derivative.
2. Jacobian passes but Vela knee remains too late: run controlled alpha/model parameter parity checks and source-support semantics checks; do not alter SG flux divergence.
3. Jacobian passes and active-edge flux factors remain about 0.73x: focus on absolute quasi-Fermi/effective-ni alignment and source-support averaging, not current extraction.
4. Jacobian passes and no source/state diagnostic moves the knee: document a model/discretization limit and keep the `-20 V` run as convergence validation only.
```

- [ ] **Step 4: Record the decision**

Generate the compact decision note from the real reports:

```powershell
$audit = Import-Csv build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_real_state_jacobian_audit\jacobian_blocks_real_state.csv
$knee = Get-Content build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_real_state_jacobian_audit\knee_shape_current.json | ConvertFrom-Json
$maxSg = ($audit | Where-Object { $_.block -eq "sg_avalanche" } | ForEach-Object { [double]$_.rel_diff } | Measure-Object -Maximum).Maximum
$maxSrh = ($audit | Where-Object { $_.block -eq "srh_auger" } | ForEach-Object { [double]$_.rel_diff } | Measure-Object -Maximum).Maximum
$maxTransport = ($audit | Where-Object { $_.block -eq "transport" } | ForEach-Object { [double]$_.rel_diff } | Measure-Object -Maximum).Maximum
$branch = if ($maxSg -gt 5e-5) {
  "missing real-state SG avalanche Jacobian block"
} elseif ($maxTransport -gt 1e-4) {
  "ordinary transport high-field mobility Jacobian isolation"
} elseif ($maxSrh -gt 1e-4) {
  "SRH/Auger derivative isolation"
} elseif ($null -eq $knee.candidate.first_growth_over_1p5) {
  "curve-shape physics/state experiment with no production config change"
} else {
  "candidate curve-shape promotion review"
}
@"

### PN2D BV Real-State Decision

Real-state Jacobian block replay was run on `-11.5`, `-13.2`, `-18`, `-19`, and
`-20 V` states. The selected next branch is `$branch`.

Evidence:

- SG avalanche maximum relative Jacobian difference: `$maxSg`.
- SRH/Auger maximum relative Jacobian difference: `$maxSrh`.
- Transport maximum relative Jacobian difference: `$maxTransport`.
- Current knee gate: Sentaurus `>1.5` at `$($knee.reference.first_growth_over_1p5) V`; Vela candidate `$($knee.candidate.first_growth_over_1p5)`.
- Next production change: none until the selected branch produces a full-curve improvement.
"@ | Add-Content docs\validation\pn2d_sentaurus_comparison.md
```

- [ ] **Step 5: Commit**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git add docs/validation/pn2d_sentaurus_comparison.md
git commit -m "Record PN2D BV real-state decision"
```

---

### Task 9: Final Verification

**Files:**
- No source edits unless a previous task failed and required a targeted fix.

- [ ] **Step 1: Build all touched C++ targets**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --target vela_example_runner test_dc_sweep test_newton_solver test_impact_ionization pn2d_jacobian_block_audit --parallel
```

Expected: build succeeds.

- [ ] **Step 2: Run focused C++ tests**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
.\build-release\test_dc_sweep.exe "[dc_sweep]"
.\build-release\test_newton_solver.exe "[newton][coupled]"
.\build-release\test_impact_ionization.exe "[impact]~[gummel]"
```

Expected: all pass.

- [ ] **Step 3: Run Python regression coverage**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python -m unittest tests.regression.test_reference_tcad_tools
```

Expected: all tests pass. If the historical crash count appears in unrelated binaries, record it as pre-existing and do not conflate it with this plan.

- [ ] **Step 4: Run real BV diagnostics**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts\run_pn2d_bv_real_state_jacobian_audit.py
python scripts\diagnose_pn2d_bv_knee_shape.py `
  --output-json build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_real_state_jacobian_audit\knee_shape_current.json
```

Expected: both commands complete and write finite reports.

- [ ] **Step 5: Final status**

Run:

```powershell
$audit = Import-Csv build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_real_state_jacobian_audit\jacobian_blocks_real_state.csv
$knee = Get-Content build-release\reference_tcad\pn2d_sentaurus2018\reports\bv_real_state_jacobian_audit\knee_shape_current.json | ConvertFrom-Json
$maxSg = ($audit | Where-Object { $_.block -eq "sg_avalanche" } | ForEach-Object { [double]$_.rel_diff } | Measure-Object -Maximum).Maximum
$maxSrh = ($audit | Where-Object { $_.block -eq "srh_auger" } | ForEach-Object { [double]$_.rel_diff } | Measure-Object -Maximum).Maximum
$maxTransport = ($audit | Where-Object { $_.block -eq "transport" } | ForEach-Object { [double]$_.rel_diff } | Measure-Object -Maximum).Maximum
"real-state sg_avalanche max rel_diff = $maxSg"
"real-state srh_auger max rel_diff = $maxSrh"
"real-state transport max rel_diff = $maxTransport"
"Sentaurus knee >1.5 / >2.0 = $($knee.reference.first_growth_over_1p5) V / $($knee.reference.first_growth_over_2p0) V"
"Vela knee >1.5 / >2.0 = $($knee.candidate.first_growth_over_1p5) / $($knee.candidate.first_growth_over_2p0)"
```

---

## Non-Goals

- Do not rewrite the core SG flux-divergence discretization as part of this plan.
- Do not introduce `source_geometry_scale` or any hidden scalar calibration knob.
- Do not promote a production BV config change unless it moves the full `-10 V..-20 V` knee shape in the right direction.
- Do not remove SRH from the Sentaurus-faithful BV comparison.
- Do not treat `-20 V` convergence alone as BV parity.
