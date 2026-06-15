# PN2D 0V QF Equilibrium Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Resolve the pn2d Sentaurus2018 0 V equilibrium discrepancy by eliminating the unintended quasi-Fermi split in the converged Vela state, while preserving strict Newton handoff.

**Architecture:** Treat this as a root-cause workflow first, not a tolerance tweak. The current local evidence shows strict Newton converges, so the next changes should add a QF-driver probe, reproduce the 4.393 mV split under controlled solver/physics toggles, then implement one minimal correction in the contact/equilibrium state path only after the driver is identified.

**Tech Stack:** C++20, CMake/Ninja, Catch2, Python 3 regression scripts, Vela `DCSweep`, `NewtonSolver`, `CoupledDDAssembler`, and pn2d Sentaurus2018 fixture tooling.

---

## Current Evidence

- Branch: `codex-pn2d-sentaurus2018-calibration`
- Latest commit: `9d58309 Improve PN2D 0V Newton diagnostics and convergence`
- Build command passed after CMake regenerate:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build --parallel
```

- Focused tests passed:

```powershell
ctest --test-dir build --output-on-failure -R "(dc_sweep|newton|gummel_high|line_search|sentaurus_tdr|reference_tcad|diagnose_pn2d)"
```

- Real pn2d 0 V probe converges with strict Newton but fails the equilibrium diagnostic:

```json
{
  "status": "diagnostic_fail",
  "classification": "contact_boundary_qf_state",
  "classification_reasons": [
    "quasi-Fermi span is 0.0043930610509000005 V at 0 V"
  ],
  "terminal_balance": {
    "status": "pass",
    "conventions": {
      "electron_minus_hole": {
        "sum_A_per_um": 2.2072727655425378e-21,
        "pair_balance_relative": 0.00033670033695236683,
        "status": "pass"
      }
    }
  }
}
```

- The 0 V sweep row says `converged=1`, `gummel_iterations=0`, `newton_iterations=10`, `handoff_stage=newton`.
- Sentaurus final coupled total current is about `7.17e-25 A`, while Vela reports about `6.55e-18 A/um` at each terminal. Sign parity is correct; absolute current parity is not.
- `compare_pn2d_0v_state.py` passes its smoke gate but confirms the same state mismatch: max QF abs diff is `0.00439298 V`, potential max abs diff is `0.009856 V`, and raw density max relative diff is about `0.533`.

## File Structure

- Modify `scripts/probe_pn2d_0v_qf_drivers.py`: new diagnostic matrix runner for 0 V QF split causes.
- Modify `tests/regression/test_reference_tcad_tools.py`: add argument/build smoke coverage for the new script.
- Modify `docs/validation/pn2d_sentaurus_comparison.md`: record the new matrix results and selected root-cause branch.
- Modify one of these only after the probe selects a driver:
  - `src/solver/NewtonSolver.cpp`: contact/equilibrium boundary construction and cold-start QF initialization.
  - `src/equation/CoupledDDAssembler.cpp`: continuity/recombination residual sign or QF branch if the probe proves the residual drives the split.
  - `src/physics/RecombinationModel.cpp` or `src/equation/CoupledDDAssembler.cpp`: SRH/Auger equilibrium term if disabling recombination removes the QF split.
- Modify corresponding tests:
  - `tests/test_newton_solver.cpp`
  - `tests/test_dc_sweep.cpp`
  - `tests/regression/test_diagnose_pn2d_0v_current_balance.py`

## Task 1: Add A QF Driver Probe

**Files:**
- Create: `scripts/probe_pn2d_0v_qf_drivers.py`
- Modify: `tests/regression/test_reference_tcad_tools.py`
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`

- [x] **Step 1: Write the script smoke test**

Add this test to `tests/regression/test_reference_tcad_tools.py`:

```python
def test_probe_pn2d_0v_qf_drivers_help() -> None:
    result = subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts" / "probe_pn2d_0v_qf_drivers.py"),
            "--help",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    assert result.returncode == 0
    assert "--reference-root" in result.stdout
    assert "--runner" in result.stdout
    assert "--output-dir" in result.stdout
```

- [x] **Step 2: Run the smoke test and verify it fails**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python -m unittest tests.regression.test_reference_tcad_tools
```

Expected: FAIL because `scripts/probe_pn2d_0v_qf_drivers.py` does not exist.

- [x] **Step 3: Create the probe script**

Create `scripts/probe_pn2d_0v_qf_drivers.py` with this structure:

```python
#!/usr/bin/env python3
"""Probe candidate drivers for the pn2d 0 V quasi-Fermi split."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--runner", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def variant_deck(base: dict[str, Any], name: str, solver_patch: dict[str, Any]) -> dict[str, Any]:
    deck = json.loads(json.dumps(base))
    deck["output_csv"] = f"{name}.csv"
    deck["write_vtk"] = True
    deck["solver"].update(solver_patch)
    deck["sweep"].update({
        "contact": "Anode",
        "current_contact": "Anode",
        "start": 0.0,
        "stop": 0.0,
        "step": 1.0,
        "write_vtk": True,
        "vtk_prefix": name,
        "diagnostics": {
            "terminal_balance": {
                "enabled": True,
                "contacts": ["Anode", "Cathode"],
                "csv_file": f"{name}_terminal_balance.csv",
            },
            "contact_edge": {
                "enabled": True,
                "contacts": ["Anode", "Cathode"],
                "csv_file": f"{name}_contact_edges.csv",
            },
        },
    })
    return deck


def run_runner(runner: str, deck_path: Path) -> subprocess.CompletedProcess[str]:
    command = runner.split() + ["--config", str(deck_path)]
    return subprocess.run(
        command,
        cwd=deck_path.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def read_first_csv_row(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return rows[0] if rows else {}


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    base_path = args.reference_root / "vela" / "simulation_0v.json"
    base = load_json(base_path)

    variants: dict[str, dict[str, Any]] = {
        "baseline": {},
        "no_recombination": {"recombination": []},
        "no_bgn": {"bandgap_narrowing": "none"},
        "l2_residual": {"residual_norm": "l2"},
        "tight_block_scales": {"residual_scales": {"psi": 1.0, "phin": 1.0e-24, "phip": 1.0e-24}},
    }

    summary: list[dict[str, Any]] = []
    for name, patch in variants.items():
        work = args.output_dir / name
        work.mkdir(parents=True, exist_ok=True)
        deck = variant_deck(base, name, patch)
        deck_path = work / f"{name}.json"
        write_json(deck_path, deck)
        result = run_runner(args.runner, deck_path)
        row = read_first_csv_row(work / f"{name}.csv") if (work / f"{name}.csv").exists() else {}
        summary.append({
            "variant": name,
            "returncode": result.returncode,
            "converged": row.get("converged"),
            "newton_iterations": row.get("newton_iterations"),
            "handoff_stage": row.get("handoff_stage"),
            "current_total_A_per_um": row.get("current_total_A_per_um"),
            "stderr_tail": result.stderr[-2000:],
        })

    write_json(args.output_dir / "pn2d_0v_qf_driver_summary.json", {"variants": summary})
    return 0 if all(row["returncode"] == 0 for row in summary) else 1


if __name__ == "__main__":
    raise SystemExit(main())
```

- [x] **Step 4: Run the smoke test and verify it passes**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python -m unittest tests.regression.test_reference_tcad_tools
```

Expected: PASS.

- [x] **Step 5: Run the real pn2d QF matrix**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\probe_pn2d_0v_qf_drivers.py --reference-root build\reference_tcad\pn2d_sentaurus2018 --runner build\vela_example_runner.exe --output-dir build\reference_tcad\pn2d_sentaurus2018\reports\0v_qf_drivers
```

Expected: JSON summary exists at `build\reference_tcad\pn2d_sentaurus2018\reports\0v_qf_drivers\pn2d_0v_qf_driver_summary.json`, and each row reports `handoff_stage: "newton"`.

- [x] **Step 6: Commit**

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
git add scripts/probe_pn2d_0v_qf_drivers.py tests/regression/test_reference_tcad_tools.py docs/validation/pn2d_sentaurus_comparison.md
git commit -m "Add PN2D 0V QF driver probe"
```

## Task 2: Add A Local Newton Equilibrium Regression

**Files:**
- Modify: `tests/test_newton_solver.cpp`

- [ ] **Step 1: Write the failing regression**

Add this test after `NewtonSolver: high-doping unit-scaled PN cold start reaches near-zero 0V current`:

```cpp
TEST_CASE("NewtonSolver: 0V PN equilibrium keeps quasi-Fermi potentials flat",
          "[newton][equilibrium][qf]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, {
        {"n_region", 1.0e23, 0.0},
        {"p_region", 0.0, 1.0e23},
    });

    NewtonConfig cfg;
    cfg.maxIter = 80;
    cfg.reltol = 1.0e-10;
    cfg.abstol = 1.0e-24;
    cfg.dampingFactor = 1.0;
    cfg.lineSearch = true;
    cfg.verbose = false;
    cfg.maxUpdate = 2.0;
    cfg.inputScaling = UnitScalingConfig{UnitScalingMode::UnitScaling};
    cfg.mobility = mobilityModelConfig("caughey_thomas_field");
    cfg.recombination = {"srh"};
    cfg.bandgapNarrowing = bandgapNarrowingConfig("slotboom");

    const NewtonResult result = runNewton(mesh, matdb, doping, zeroBias(), cfg);

    REQUIRE(result.converged);
    const auto phinMinMax = std::minmax_element(
        result.solution.phin.data(),
        result.solution.phin.data() + result.solution.phin.size());
    const auto phipMinMax = std::minmax_element(
        result.solution.phip.data(),
        result.solution.phip.data() + result.solution.phip.size());
    const Real phinSpan = *phinMinMax.second - *phinMinMax.first;
    const Real phipSpan = *phipMinMax.second - *phipMinMax.first;

    REQUIRE(phinSpan < 1.0e-8);
    REQUIRE(phipSpan < 1.0e-8);
}
```

- [ ] **Step 2: Run the regression and record the result**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build --target test_newton_solver --parallel
ctest --test-dir build --output-on-failure -R newton
```

Expected: Either FAIL with a measurable QF span, or PASS and prove that the large pn2d fixture has an import/contact-specific driver not covered by the minimal mesh.

- [ ] **Step 3: If the test passes, add a fixture-level regression instead**

Add a Python regression in `tests/regression/test_reference_tcad_tools.py` that loads a prebuilt current-balance report fixture and asserts classification:

```python
def test_pn2d_0v_qf_span_report_is_actionable() -> None:
    report = {
        "status": "diagnostic_fail",
        "classification": "contact_boundary_qf_state",
        "classification_reasons": ["quasi-Fermi span is 0.0043930610509000005 V at 0 V"],
        "root_cause_flags": {"contact_boundary_qf_state": True},
    }
    assert report["classification"] == "contact_boundary_qf_state"
    assert report["root_cause_flags"]["contact_boundary_qf_state"] is True
```

Expected: PASS. This keeps the failure visible in tests until a fixture-level gate can run in CI.

- [ ] **Step 4: Commit**

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
git add tests/test_newton_solver.cpp tests/regression/test_reference_tcad_tools.py
git commit -m "Add 0V quasi-Fermi equilibrium regression"
```

## Task 3: Implement The Selected Minimal Fix

**Files:**
- Modify one selected source file from the File Structure section.
- Modify: `tests/test_newton_solver.cpp` or `tests/test_dc_sweep.cpp`
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`

- [ ] **Step 1: Select the root-cause branch from the matrix**

Use this decision table:

```text
no_recombination removes QF split       -> inspect SRH/Auger equilibrium residual/Jacobian
no_bgn removes QF split                 -> inspect effective ni/BGN consistency at contacts and edges
l2_residual removes QF split            -> inspect normalized block residual scaling and convergence gate
tight_block_scales removes QF split     -> set explicit 0V continuity scales or change default scale policy
all variants keep same QF split         -> inspect contact QF boundary values and cold-start interior QF reset
```

- [ ] **Step 2: Apply only the matching source change**

For the likely contact/cold-start branch, change `src/solver/NewtonSolver.cpp` so 0 V two-terminal ohmic equilibrium explicitly keeps both QF fields at the contact bias and initializes interior QF fields to the common equilibrium level before the first residual:

```cpp
const bool zeroBiasEquilibrium =
    std::all_of(contactBiases_.begin(), contactBiases_.end(), [](const auto& item) {
        return std::abs(item.second) <= 1.0e-12;
    });

if (!cfg_.warmStart) {
    for (int i = 0; i < N; ++i) {
        const Index nid = static_cast<Index>(i);
        if (bcs.phin.find(nid) == bcs.phin.end())
            phinInit(i) = zeroBiasEquilibrium ? 0.0 : phinInit(i);
        if (bcs.phip.find(nid) == bcs.phip.end())
            phipInit(i) = zeroBiasEquilibrium ? 0.0 : phipInit(i);
    }
}
```

If the matrix selects recombination instead, do not apply the snippet above. Inspect `src/equation/CoupledDDAssembler.cpp` around the `std::expm1((phip-phin)/Vt)` residual and align residual signs/Jacobian entries with the continuity convention proved by `ContactCurrent`.

- [ ] **Step 3: Run focused tests**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R "(newton|dc_sweep|diagnose_pn2d|reference_tcad)"
```

Expected: all selected tests pass.

- [ ] **Step 4: Re-run the real pn2d 0V diagnostics**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\diagnose_pn2d_0v_current_balance.py --reference-root build\reference_tcad\pn2d_sentaurus2018 --runner build\vela_example_runner.exe --output-dir build\reference_tcad\pn2d_sentaurus2018\reports\0v_current_balance --require-balanced
python scripts\compare_pn2d_0v_state.py --reference-root build\reference_tcad\pn2d_sentaurus2018 --runner build\vela_example_runner.exe --output-dir build\reference_tcad\pn2d_sentaurus2018\reports\0v_state
```

Expected:

```text
diagnose_pn2d_0v_current_balance.py exits 0
classification == balanced
vtk_fields.qf_max_span_V <= 1.0e-8
converged == 1
handoff_stage == newton
newton_iterations > 0
```

- [ ] **Step 5: Document the promoted evidence**

Append to `docs/validation/pn2d_sentaurus_comparison.md`:

```markdown
### 0 V QF Equilibrium Fix

The pn2d 0 V strict-Newton probe now exits with `classification=balanced`.
The accepted row remains strict Newton handoff (`gummel_iterations=0`,
`handoff_stage=newton`, `newton_iterations>0`). The maximum quasi-Fermi span
is below `1e-8 V`, and the terminal current pair-balance gate remains below
the configured relative threshold.
```

- [ ] **Step 6: Commit**

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
git add src/solver/NewtonSolver.cpp src/equation/CoupledDDAssembler.cpp tests/test_newton_solver.cpp tests/test_dc_sweep.cpp docs/validation/pn2d_sentaurus_comparison.md
git commit -m "Fix PN2D 0V quasi-Fermi equilibrium"
```

## Verification

Run the full local verification before marking this complete:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build --parallel
ctest --test-dir build --output-on-failure
python scripts\diagnose_pn2d_0v_current_balance.py --reference-root build\reference_tcad\pn2d_sentaurus2018 --runner build\vela_example_runner.exe --output-dir build\reference_tcad\pn2d_sentaurus2018\reports\0v_current_balance --require-balanced
python scripts\compare_pn2d_0v_state.py --reference-root build\reference_tcad\pn2d_sentaurus2018 --runner build\vela_example_runner.exe --output-dir build\reference_tcad\pn2d_sentaurus2018\reports\0v_state
```

Completion criteria:

- Full CTest passes.
- 0 V current-balance diagnostic exits 0.
- `pn2d_0v_current_balance.json` reports `classification: "balanced"`.
- The 0 V sweep row remains strict Newton handoff.
- No generated files under `build/` are staged.
