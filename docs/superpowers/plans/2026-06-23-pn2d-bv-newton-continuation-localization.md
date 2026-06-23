# PN2D BV Newton Continuation Localization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Localize and fix the PN2D BV continuation regression introduced between `c04edbf` and `51d8bf1`, where near-floor residual states produce raw Newton steps around `100+` and prevent continuation toward `-20 V`.

**Architecture:** Treat this as a solver/globalization regression first, not a Sentaurus artifact or avalanche-model parity task. Preserve the existing imported PN2D reference inputs and isolate changes through reproducible short BV decks, Newton history/failure diagnostics, and focused tests around Newton step capping, Poisson recorrection, and carrier-row behavior.

**Tech Stack:** C++20, CMake/Ninja with MSYS2 UCRT64, Catch2 tests, Python regression/diagnostic scripts, Vela `vela_example_runner` DC sweep JSON decks.

---

## Current Evidence

- `c04edbf` with `max_update=0` and `quasi_fermi_update_limit_V=0.1` reaches `-3 V` in 61 points; final raw step is about `4.50` and residual is about `4.07e-13`.
- `51d8bf1` with the same deck fails just past `-0.3 V` with `line_search_non_decrease`, residual about `8.25e-9`, raw step about `124.30`, positive finite carriers.
- `85924e5` with the same deck fails just past `-0.40625 V` with `line_search_non_decrease`, residual about `1.93e-8`, raw step about `107.59`, positive finite carriers.
- `51d8bf1` and `85924e5` with `impact_ionization.model = "none"` both fail identically just past `-0.8165 V` with `max_iterations`, residual about `8.58e-9`, raw step about `113.1`, positive finite carriers.
- Therefore the next target is the `51d8bf1` Newton/transport/globalization change set, not `85924e5` contact fallback and not Sentaurus artifact import.

## File Map

- Modify or create diagnostics under `scripts/`:
  - `scripts/compare_pn2d_bv_newton_history.py`: compare Newton history tails across two localization runs.
  - `scripts/prepare_pn2d_bv_localization_decks.py`: generate consistent short BV decks from one imported PN2D base config.
- Modify focused solver tests:
  - `tests/test_newton_solver.cpp`: add or extend tests for QF-only step limiting and Poisson recorrection behavior.
  - `tests/test_impact_ionization.cpp`: keep existing impact tests passing; do not use it as the primary continuation regression gate.
- Inspect before modifying:
  - `src/solver/NewtonSolver.cpp`
  - `include/vela/solver/NewtonSolver.h`
  - `src/equation/CoupledDDAssembler.cpp`
  - `include/vela/equation/AssemblerUtils.h`
  - `src/simulation/DCSweep.cpp`
- Generated artifacts stay under:
  - `build-release/bv_localization/`
  - Do not commit generated CSV/JSON/VTK outputs.

---

### Task 1: Freeze The Localization Harness

**Files:**
- Create: `scripts/prepare_pn2d_bv_localization_decks.py`
- Test: `tests/regression/test_reference_tcad_tools.py`

- [ ] **Step 1: Add a deck generator script**

Create `scripts/prepare_pn2d_bv_localization_decks.py`:

```python
#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-config", required=True, type=Path)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--out-root", required=True, type=Path)
    parser.add_argument("--case-name", required=True)
    parser.add_argument("--stop", required=True, type=float)
    parser.add_argument("--qf-limit", type=float)
    parser.add_argument("--max-update", type=float)
    parser.add_argument("--impact-model", choices=["keep", "none"], default="keep")
    parser.add_argument("--diagnostics", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    base = json.loads(args.base_config.read_text(encoding="utf-8-sig"))
    out_dir = args.out_root / args.case_name
    out_dir.mkdir(parents=True, exist_ok=True)

    base["mesh_file"] = str(args.reference_root / "vela" / "mesh.json")
    base["node_doping_file"] = str(args.reference_root / "doping.csv")
    base["materials_file"] = str(
        args.reference_root / "vela" / "pn2d_sentaurus2018_iv_materials.json"
    )
    base["output_csv"] = str(out_dir / "iv.csv")

    solver = base.setdefault("solver", {})
    if args.max_update is not None:
        solver["max_update"] = args.max_update
    if args.qf_limit is not None:
        solver["quasi_fermi_update_limit_V"] = args.qf_limit
    if args.impact_model == "none":
        solver["impact_ionization"] = {"model": "none"}
    if args.diagnostics:
        solver["diagnostics"] = True

    sweep = base.setdefault("sweep", {})
    sweep["stop"] = args.stop
    sweep["step"] = -0.05
    sweep["write_vtk"] = False
    sweep["max_step"] = 0.05
    sweep["min_step"] = 1.0e-10
    sweep["max_retries"] = 29
    if args.diagnostics:
        sweep["diagnostics"] = {
            "newton_history": {
                "enabled": True,
                "csv_file": str(out_dir / "newton_history.csv"),
            }
        }

    config_path = out_dir / "config.json"
    config_path.write_text(json.dumps(base, indent=2), encoding="utf-8")
    print(config_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Add a regression test for the deck generator**

Append a test to `tests/regression/test_reference_tcad_tools.py`:

```python
def test_prepare_pn2d_bv_localization_deck_writes_absolute_inputs(self) -> None:
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        reference_root = tmp / "reference"
        (reference_root / "vela").mkdir(parents=True)
        base = tmp / "base.json"
        base.write_text(json.dumps({
            "mesh_file": "mesh.json",
            "node_doping_file": "doping.csv",
            "materials_file": "materials.json",
            "output_csv": "old.csv",
            "solver": {
                "max_update": 5.0,
                "impact_ionization": {"model": "van_overstraeten"},
            },
            "sweep": {"stop": -0.05},
        }), encoding="utf-8")
        out_root = tmp / "out"
        run([
            sys.executable,
            str(REPO / "scripts" / "prepare_pn2d_bv_localization_decks.py"),
            "--base-config", str(base),
            "--reference-root", str(reference_root),
            "--out-root", str(out_root),
            "--case-name", "qflim",
            "--stop", "-3",
            "--qf-limit", "0.1",
            "--max-update", "0",
            "--diagnostics",
        ], check=True)
        deck = json.loads((out_root / "qflim" / "config.json").read_text())
        self.assertEqual(deck["sweep"]["stop"], -3.0)
        self.assertEqual(deck["solver"]["max_update"], 0.0)
        self.assertEqual(deck["solver"]["quasi_fermi_update_limit_V"], 0.1)
        self.assertTrue(Path(deck["mesh_file"]).is_absolute())
        self.assertTrue(deck["sweep"]["diagnostics"]["newton_history"]["enabled"])
```

- [ ] **Step 3: Run the new regression test**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_prepare_pn2d_bv_localization_deck_writes_absolute_inputs -v
```

Expected: `OK`.

- [ ] **Step 4: Generate the canonical decks**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts\prepare_pn2d_bv_localization_decks.py --base-config build-release\reference_tcad\pn2d_sentaurus2018\vela\simulation_bv.json --reference-root build-release\reference_tcad\pn2d_sentaurus2018 --out-root build-release\bv_localization\canonical --case-name qflim0p1_to_3V_diag --stop -3 --qf-limit 0.1 --max-update 0 --diagnostics
python scripts\prepare_pn2d_bv_localization_decks.py --base-config build-release\reference_tcad\pn2d_sentaurus2018\vela\simulation_bv.json --reference-root build-release\reference_tcad\pn2d_sentaurus2018 --out-root build-release\bv_localization\canonical --case-name qflim0p1_noimpact_to_3V_diag --stop -3 --qf-limit 0.1 --max-update 0 --impact-model none --diagnostics
```

Expected: both commands print a `config.json` path.

---

### Task 2: Compare Newton History Across The Regression Boundary

**Files:**
- Create: `scripts/compare_pn2d_bv_newton_history.py`
- Test: `tests/regression/test_reference_tcad_tools.py`

- [ ] **Step 1: Add the Newton history comparator**

Create `scripts/compare_pn2d_bv_newton_history.py`:

```python
#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def last_float(rows: list[dict[str, str]], column: str) -> float | None:
    for row in reversed(rows):
        value = row.get(column, "")
        if value:
            return float(value)
    return None


def summarize(path: Path) -> dict[str, object]:
    rows = read_rows(path)
    return {
        "path": str(path),
        "rows": len(rows),
        "last_bias_V": last_float(rows, "bias_V"),
        "last_residual_norm": last_float(rows, "residual_norm"),
        "last_raw_step_norm": last_float(rows, "raw_step_norm"),
        "last_block_psi": last_float(rows, "block_psi"),
        "last_block_phin": last_float(rows, "block_phin"),
        "last_block_phip": last_float(rows, "block_phip"),
        "max_raw_step_norm": max(
            (float(row["raw_step_norm"]) for row in rows if row.get("raw_step_norm")),
            default=None,
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--left", required=True, type=Path)
    parser.add_argument("--right", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    left = summarize(args.left)
    right = summarize(args.right)
    comparison = {
        "left": left,
        "right": right,
        "raw_step_ratio": (
            right["last_raw_step_norm"] / left["last_raw_step_norm"]
            if left["last_raw_step_norm"] and right["last_raw_step_norm"]
            else None
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(comparison, indent=2), encoding="utf-8")
    print(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 2: Add a comparator regression test**

Append a test to `tests/regression/test_reference_tcad_tools.py`:

```python
def test_compare_pn2d_bv_newton_history_reports_raw_step_ratio(self) -> None:
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td)
        left = tmp / "left.csv"
        right = tmp / "right.csv"
        header = "bias_V,iteration,residual_norm,raw_step_norm,block_psi,block_phin,block_phip\n"
        left.write_text(header + "-3,3,4e-13,4.5,3e-12,5e-14,5e-14\n", encoding="utf-8")
        right.write_text(header + "-0.3,6,7e-9,90,3e-8,1e-9,7e-9\n", encoding="utf-8")
        out = tmp / "summary.json"
        run([
            sys.executable,
            str(REPO / "scripts" / "compare_pn2d_bv_newton_history.py"),
            "--left", str(left),
            "--right", str(right),
            "--out", str(out),
        ], check=True)
        summary = json.loads(out.read_text())
        self.assertAlmostEqual(summary["raw_step_ratio"], 20.0)
        self.assertEqual(summary["left"]["rows"], 1)
        self.assertEqual(summary["right"]["rows"], 1)
```

- [ ] **Step 3: Run the comparator regression test**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_compare_pn2d_bv_newton_history_reports_raw_step_ratio -v
```

Expected: `OK`.

- [ ] **Step 4: Compare current localization histories**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts\compare_pn2d_bv_newton_history.py --left build-release\bv_localization\c04edbf\qflim0p1_to_3V_diag\newton_history.csv --right build-release\bv_localization\51d8bf1\qflim0p1_to_3V_diag\newton_history.csv --out build-release\bv_localization\reports\c04_vs_51_newton_history.json
python scripts\compare_pn2d_bv_newton_history.py --left build-release\bv_localization\c04edbf\qflim0p1_to_3V_diag\newton_history.csv --right build-release\bv_localization\85924e5\qflim0p1_to_3V_diag\newton_history.csv --out build-release\bv_localization\reports\c04_vs_859_newton_history.json
```

Expected: JSON summaries show `right.last_raw_step_norm` around `90..110+`, much larger than the `c04edbf` tail.

---

### Task 3: Audit The `51d8bf1` Step-Cap And Poisson Recorrection Change

**Files:**
- Inspect: `src/solver/NewtonSolver.cpp:397-430`
- Inspect: `src/solver/NewtonSolver.cpp:1803-1811`
- Modify: `tests/test_newton_solver.cpp`
- Possible modify: `src/solver/NewtonSolver.cpp`

- [ ] **Step 1: Read the relevant code blocks**

Run:

```powershell
Get-Content src\solver\NewtonSolver.cpp | Select-Object -Skip 380 -First 70
Get-Content src\solver\NewtonSolver.cpp | Select-Object -Skip 1788 -First 35
```

Expected: identify exactly how `applyConfiguredStepCapsAndPoissonRecorrection` applies global caps, QF caps, and Poisson recorrection.

- [ ] **Step 2: Add a focused unit test for QF cap not inflating psi step**

Add to `tests/test_newton_solver.cpp`:

```cpp
TEST_CASE("Newton step caps keep Poisson correction bounded after QF clipping",
          "[newton][bv][continuation]")
{
    const int N = 2;
    VectorXd step(3 * N);
    step << 1.0, -1.0, 100.0, -100.0, 100.0, -100.0;
    SparseMatrixd J(3 * N, 3 * N);
    J.setIdentity();
    VectorXd residual = VectorXd::Zero(3 * N);

    NewtonConfig cfg;
    cfg.maxUpdate = 0.0;
    cfg.quasiFermiUpdateLimit_V = 0.1;

    detail::applyConfiguredStepCapsAndPoissonRecorrectionForTest(
        step, J, residual, cfg, N, 1.0);

    REQUIRE(std::abs(step(2)) <= Catch::Approx(0.1));
    REQUIRE(std::abs(step(3)) <= Catch::Approx(0.1));
    REQUIRE(std::abs(step(4)) <= Catch::Approx(0.1));
    REQUIRE(std::abs(step(5)) <= Catch::Approx(0.1));
    REQUIRE(std::isfinite(step.norm()));
    REQUIRE(step.head(N).cwiseAbs().maxCoeff() <= 1.0);
}
```

If the helper is not externally visible, add a tiny test-only wrapper near the helper in `src/solver/NewtonSolver.cpp` and declare it in `include/vela/solver/NewtonSolver.h` under the local `detail` namespace.

- [ ] **Step 3: Run the focused test**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --parallel --target test_newton_solver
build-release\test_newton_solver.exe "[newton][bv][continuation]"
```

Expected before implementation: fail to compile if the wrapper is missing, or fail the bounded-step assertions if current recorrection inflates the psi block.

- [ ] **Step 4: Implement only the minimum code needed for the test**

If the test shows unbounded recorrection, change `applyConfiguredStepCapsAndPoissonRecorrection` so Poisson recorrection is either:

```cpp
if (cfg.quasiFermiUpdateLimit_V > 0.0) {
    const Real beforePsiMax = step.head(nodeCount).cwiseAbs().maxCoeff();
    recorrectPoissonStepForClippedQuasiFermi(step, jacobian, residual, nodeCount);
    const Real afterPsiMax = step.head(nodeCount).cwiseAbs().maxCoeff();
    if (beforePsiMax > 0.0 && afterPsiMax > beforePsiMax) {
        step.head(nodeCount) *= beforePsiMax / afterPsiMax;
    }
}
```

or skipped when the carrier block is already clipped and Poisson residual is below the numerical floor. Choose the smaller change that matches the observed failing test and does not mask residual errors.

- [ ] **Step 5: Re-run the focused unit test**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --parallel --target test_newton_solver
build-release\test_newton_solver.exe "[newton][bv][continuation]"
```

Expected: PASS.

---

### Task 4: Re-run The Minimal BV Gates After The Candidate Fix

**Files:**
- Use generated configs under `build-release/bv_localization/canonical/`
- No source file modifications unless Task 3 identified a solver fix.

- [ ] **Step 1: Rebuild the runner**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --parallel --target vela_example_runner
```

Expected: build succeeds.

- [ ] **Step 2: Run the impact-on QF-limit gate to `-3 V`**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\vela_example_runner.exe --config build-release\bv_localization\canonical\qflim0p1_to_3V_diag\config.json
```

Expected: `converged=true`, final CSV row at about `-3 V`.

- [ ] **Step 3: Run the no-impact QF-limit gate to `-3 V`**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\vela_example_runner.exe --config build-release\bv_localization\canonical\qflim0p1_noimpact_to_3V_diag\config.json
```

Expected: `converged=true`, final CSV row at about `-3 V`.

- [ ] **Step 4: If either `-3 V` gate fails, stop and report**

Do not continue to `-20 V`. Extract:

```powershell
python - <<'PY'
import csv, json, pathlib
for path in pathlib.Path("build-release/bv_localization/canonical").glob("*/iv.csv"):
    rows = list(csv.DictReader(path.open()))
    print(path, rows[-1])
PY
```

Expected: final report names the failed bias, failure reason, residual norm, raw step norm, and whether carriers are positive finite.

---

### Task 5: Only After `-3 V` Passes, Stage Toward `-20 V`

**Files:**
- Use `scripts/prepare_pn2d_bv_localization_decks.py`
- Generated outputs only under `build-release/bv_localization/canonical/`

- [ ] **Step 1: Generate staged decks**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts\prepare_pn2d_bv_localization_decks.py --base-config build-release\reference_tcad\pn2d_sentaurus2018\vela\simulation_bv.json --reference-root build-release\reference_tcad\pn2d_sentaurus2018 --out-root build-release\bv_localization\canonical --case-name qflim0p1_to_5V_diag --stop -5 --qf-limit 0.1 --max-update 0 --diagnostics
python scripts\prepare_pn2d_bv_localization_decks.py --base-config build-release\reference_tcad\pn2d_sentaurus2018\vela\simulation_bv.json --reference-root build-release\reference_tcad\pn2d_sentaurus2018 --out-root build-release\bv_localization\canonical --case-name qflim0p1_to_10V_diag --stop -10 --qf-limit 0.1 --max-update 0 --diagnostics
python scripts\prepare_pn2d_bv_localization_decks.py --base-config build-release\reference_tcad\pn2d_sentaurus2018\vela\simulation_bv.json --reference-root build-release\reference_tcad\pn2d_sentaurus2018 --out-root build-release\bv_localization\canonical --case-name qflim0p1_to_20V_diag --stop -20 --qf-limit 0.1 --max-update 0 --diagnostics
```

- [ ] **Step 2: Run the `-5 V` stage**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\vela_example_runner.exe --config build-release\bv_localization\canonical\qflim0p1_to_5V_diag\config.json
```

Expected: `converged=true`. If false, stop and extract failure diagnostics.

- [ ] **Step 3: Run the `-10 V` stage only if `-5 V` passes**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\vela_example_runner.exe --config build-release\bv_localization\canonical\qflim0p1_to_10V_diag\config.json
```

Expected: `converged=true`. If false, stop and extract failure diagnostics.

- [ ] **Step 4: Run the `-20 V` stage only if `-10 V` passes**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\vela_example_runner.exe --config build-release\bv_localization\canonical\qflim0p1_to_20V_diag\config.json
```

Expected: this is a diagnostic milestone, not final BV parity acceptance.

---

### Task 6: Verification And Documentation

**Files:**
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`
- Modify: `docs/validation/pn2d_bv_validation.md`

- [ ] **Step 1: Run focused C++ tests**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
build-release\test_newton_solver.exe "[newton]"
build-release\test_impact_ionization.exe
ctest --test-dir build-release --output-on-failure -R "impact|newton"
```

Expected: all selected tests pass.

- [ ] **Step 2: Run focused Python regression tests**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python -m unittest tests.regression.test_reference_tcad_tools -v
```

Expected: `OK`.

- [ ] **Step 3: Update validation docs with the new evidence**

In `docs/validation/pn2d_sentaurus_comparison.md`, add a subsection after the 2026-06-22 QF-limit staged sweep follow-up:

```markdown
### PN2D BV Newton Continuation Regression Localization (2026-06-23)

The regression boundary is `c04edbf -> 51d8bf1`. With the same imported PN2D BV
deck, `max_update=0`, and `quasi_fermi_update_limit_V=0.1`, `c04edbf` reaches
`-3 V`, while `51d8bf1` fails just beyond `-0.3 V` with a near-floor residual
and a raw Newton step above `100`. A no-impact control on `51d8bf1` and
`85924e5` fails identically just beyond `-0.8165 V`, which localizes the
remaining blocker to Newton/transport/globalization rather than the
`85924e5` contact avalanche driving-field fallback.
```

- [ ] **Step 4: Record final status**

Report:

- Which stages passed: `-3 V`, `-5 V`, `-10 V`, `-20 V`.
- The exact last stable and failed bias for any failed stage.
- The residual norm, raw step norm, line-search reason, and positive-finite carrier status.
- Whether full `-20 V` parity remains blocked or is ready for the existing multibias comparison gate.


### PN2D BV High-Bias Feedback Localization Update (2026-06-23)

Status after the frozen high-field mobility Jacobian candidate:

- Continuation reachability is unblocked: the frozen-Jacobian path reaches `-20 V`.
- BV parity is still open: the `-10..-20 V` Vela curve does not reproduce the Sentaurus one-volt current-growth knee.
- Field/source diagnostics show electric field is close at `-20 V`, while avalanche/source feedback is low.
- The continuity-feedback diagnostic now accepts `--material-ni-m3`; with the actual PN2D material `ni = 1.4638914958767616e16 m^-3`, Vela `ni_eff` matches Sentaurus-inferred `ni_eff` at the active focus endpoints.
- The remaining high-bias mismatch is localized to absolute `psi/phin/phip` state offsets: about `47..48 mV` in `psi-phin` and about `56 mV` in `phip-psi`, producing roughly `0.8..0.95` decades of carrier-density deficit and the low source/flux feedback.

Next task:

1. Build a minimal absolute-state feedback probe around the `-20 V` frozen-Jacobian VTK/state and Sentaurus active endpoints.
2. Test whether substituting Sentaurus-aligned `psi-phin` / `phip-psi` exponents into Vela source/flux diagnostics restores the missing source by the observed `0.8..0.95` decades.
3. If yes, inspect continuation/state initialization or branch policy for the quasi-Fermi carrier-density relation before touching production physics.
4. Do not promote alpha(E), material `ni`, BGN, or hidden source-scale changes from the current evidence.

### PN2D BV Absolute-State Feedback Probe Result (2026-06-23)

The minimal post-processing probe has been executed via
`scripts/diagnose_pn2d_bv_absolute_state_feedback.py` using the corrected
`continuity_feedback_material_ni` CSVs. For active `-20 V` edges, state-scaling
by Sentaurus absolute `psi/phin/phip` endpoint carrier-density factors recovers
`0.8827` decades of source gap; focus edge 2886 moves from `-0.8401` decades to
`+0.04295` decades relative to Sentaurus. At `-10 V`, the recovered gap is only
`0.0610` decades.

Updated next task:

1. Inspect the coupled high-bias solution branch that sets absolute `psi-phin`
   and `phip-psi` about `47..56 mV` too low at `-20 V`.
2. Compare state initialization, gauge/contact reference, and carrier-density
   reconstruction paths between Vela and Sentaurus exports.
3. Add a production-facing experiment only after identifying which branch/policy
   causes the absolute-state offset; do not promote the diagnostic source proxy
   itself.

### PN2D BV Absolute Branch Offset Probe Result (2026-06-23)

The full-node branch-offset probe has been executed via
`scripts/diagnose_pn2d_bv_absolute_branch_offsets.py` for `-2, -5, -10, -20 V`.
It shows that contact nodes are aligned at `-20 V`, while non-contact
impact-active nodes carry the high-bias branch offset: median `delta_psi` is
near zero, but median `delta_phin = +0.04609 V` and `delta_phip = -0.05398 V`,
which produces median `delta(psi-phin) = -0.04733 V` and
`delta(phip-psi) = -0.05547 V`. The active-support density deficits are
`-0.7949` decades for electrons and `-0.9318` decades for holes.

Updated next task:

1. Build a controlled high-field active-support replay/perturbation experiment
   that shifts interior quasi-Fermi branches while keeping contact Dirichlet
   values fixed.
2. Measure whether this moves the actual `-10..-20 V` curve knee or only the
   local source proxy.
3. Inspect carrier-continuity residual terms around the active support before
   proposing any production solver or boundary-policy change.
### PN2D BV Active-Support QF Shift Replay Result (2026-06-23)

The controlled active-support replay has been executed for the current frozen
high-field mobility-Jacobian `-20 V` state. Inputs were regenerated from the
current acceptance visual VTK under
`bv_frozen_mobility_jacobian_feedback/active_support_qf_shift/`: the
99th-percentile support has zero Sentaurus/Vela overlap (`20` false negatives,
`20` false positives, peak separation `0.04752 um`), and the SG edge replay uses
`quasi_fermi_gradient` with `quasi_fermi_variable_ni`.

The residual proxy now supports `--qf-shift-scope support_nodes`, backed by a
focused regression test. Three replay cases were run:

- baseline: no QF shift;
- `shift_all`: apply the active median branch offsets globally;
- `shift_support`: apply the same offsets only to thresholded active-support
  nodes.

On Sentaurus-only active nodes, baseline Vela-state transport/source ratios are
`0.0955` electron and `0.1185` hole, with residual/source medians `-0.927` and
`-0.904`. The global shift moves them to `0.596` and `1.013`, which confirms the
absolute QF/carrier-density branch is a real causal lever. The hard support-only
shift is not production-safe: it moves electron transport toward parity
(`0.844`) but overshoots hole transport to `14.20`; on Vela-only active nodes it
drives electron transport to `8.424` versus a Sentaurus-state reference near
`0.891`.

Updated next task:

1. Do not promote pointwise active-support QF shifts as a solver fix.
2. Build a smooth active-region or continuation-level branch-control experiment
   that preserves contact Dirichlet values and avoids discontinuous local QF
   jumps.
3. Judge the experiment by the actual `-10..-20 V` IV knee movement and terminal
   current parity, not only by local source/residual proxy recovery.
4. If the smooth branch-control experiment still overshoots or leaves the knee
   unchanged, inspect the carrier-continuity residual/Jacobian balance at the
   high-field active support before touching production physics.
### PN2D BV Smooth Branch-Control Backscan Result (2026-06-23)

The smooth active-region branch-control experiment has been executed. A new
helper, `scripts/prepare_pn2d_bv_smooth_branch_state.py`, prepares a
DCSweep-compatible `initial_state_file` by applying Gaussian QF shifts around
selected support nodes while preserving contact nodes at zero weight and
reconstructing carrier densities from the original Vela inferred `ni_eff`. The
focused regression test verifies the state CSV header, contact preservation,
unit-weight selected support, and sub-unit smooth neighbor weights.

The real `-20 V` run used false-negative support nodes, `decay_length_um = 0.05`,
`electron_qf_shift_v = -0.04732618818171197`, and
`hole_qf_shift_v = +0.05547180240410299`. It selected `20` support nodes, kept
`34` contact nodes unshifted, assigned nonzero smooth weights to `1008` interior
nodes, and converged a true DCSweep backscan over `-20, -19, ..., -10 V`.

Result: the curve-level knee does not move. The `-20 V` current is
`-1.1690341e-16 A/um`, only `0.00034` decades different from the frozen visual
baseline and still `-0.8914` decades below Sentaurus. The smooth `-20 -> -19 V`
growth ratio is `1.0007`, while Sentaurus is `2.204`. Residual probes on the
smooth and zero-shift `-20 V` states are Poisson dominated (`psi = 0.2443387`),
and false-negative support carrier residuals remain tiny (`~1e-14` after shift),
with local Poisson residual unchanged.

Updated next task:

1. Stop escalating QF shift strength/shape; smooth branch initialization alone
   is accepted by Newton but returns to the same low-current branch.
2. Build a mixed-state Poisson/space-charge consistency audit for the desired
   high-density active-support state: keep Vela potential, impose Sentaurus-like
   active carrier density factors, and measure Poisson residual/charge imbalance
   by support class, contact distance, and doping sign.
3. Use that audit to decide whether the missing branch is blocked by
   electrostatic charge consistency, carrier-density reconstruction policy, or a
   continuation predictor that must adjust potential and QF together.
### PN2D BV Mixed-State Charge Audit Result (2026-06-23)

The mixed-state Poisson/space-charge audit has been executed with a new helper,
`scripts/diagnose_pn2d_bv_mixed_state_charge_audit.py`. The helper keeps Vela
potential and fixed doping, replaces selected support-node carrier densities
with Sentaurus densities, integrates mobile/net charge changes over node control
volumes, and summarizes by support class, doping sign, and contact bucket.

At `-20 V`, replacing only false-negative active-support nodes changes selected
net charge by `1.4116e-23 C/m` versus `3.9116e-11 C/m` baseline net charge, a
ratio of `3.61e-13`. Replacing both false-negative and false-positive support
nodes changes `3.1196e-23 C/m`, still only `7.98e-13` of the baseline support
net charge. False-positive compensated nodes have a large relative change only
because their baseline net charge is nearly zero; the absolute charge remains
`~1.7e-23 C/m`.

Updated next task:

1. Stop treating Poisson/space-charge consistency as the leading blocker for the
   Sentaurus-like active carrier densities.
2. Inspect high-field active-edge carrier-continuity flux/Jacobian balance:
   compare SG current/flux sensitivity to the absolute QF density lever and
   explain why smooth QF initialization is accepted but returns to the same
   low-current branch.
3. Design the next experiment around coupled QF-gradient/current-density branch
   movement, not around stronger local QF shifts or Poisson charge correction.
## Active-Edge Flux Sensitivity And Restart Relaxation Probe

The next probe has been executed. Two small diagnostic improvements were added:

- `scripts/diagnose_pn2d_bv_active_edge_flux_factors.py` now accepts both old
  proxy SG-edge columns and the current generated columns
  (`electron_flux_abs`, `hole_flux_abs`, `edge_area_m2`, and optional
  `source_integral`). This fixed the false zero-active-edge replay result.
- `scripts/diagnose_pn2d_bv_restart_state_relaxation.py` compares a restart
  initial state, its converged final state, and a baseline state on the active
  support classes.

Real `-20 V` results:

- Active-edge mixed-state replay output:
  `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/active_edge_mixed_state_replay_m20_branch/`.
- Smooth restart single-point output:
  `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/smooth_branch_control/single_m20_decay0p05/`.
- On false-negative support, Vela baseline active-edge generation is `0.1380x`
  Sentaurus and particle flux is `0.1296x`; applying the uniform measured Vela QF
  branch shift recovers generation to `0.9618x` and particle flux to `0.9644x`.
  Therefore the SG source/flux path is sensitive to the absolute QF density
  lever.
- The one-point `-20 V` restart from `smooth_branch_state.csv` converges in two
  Newton iterations to `-1.169034088445e-16 A/um`, the same low-current branch as
  the backscan. False-negative support starts with median branch shifts
  `phin=-0.047326 V`, `phip=+0.055472 V` and density factors `n=6.238x`,
  `p=8.548x`, but the converged state retains only `phin=-0.009381 V`,
  `phip=+0.014140 V` and `n=1.438x`, `p=1.728x`. The retained absolute QF shift is
  `0.198x` for electrons and `0.255x` for holes.

Conclusion: the missing avalanche feedback is causally tied to the absolute QF
carrier-density branch, but local absolute-QF seeding is relaxed away by the
coupled carrier-continuity solve. The next task should target a coupled
QF-gradient/current-density branch movement or predictor around the active
high-field support, rather than stronger QF shifts, Poisson charge correction, or
source-scale calibration.

Next task:

1. Build a predictor/probe that perturbs the active high-field region using a
   coupled QF-gradient/current-density pattern, not a pure common-mode QF shift.
2. Measure whether the first Newton step preserves active-edge generation and
   terminal current growth or immediately relaxes the state back to the baseline
   QF branch.
3. If the coupled predictor is source-effective, promote it only as a diagnostic
   branch-control experiment and keep acceptance tied to the `-10..-20 V` knee
   shape gate.

## Coupled QF-Gradient/Current-Density Predictor Experiment

The selected predictor experiment has been executed. A new helper,
`scripts/prepare_pn2d_bv_coupled_qf_predictor_state.py`, constructs a diagnostic
restart state by selecting active SG edges around requested support classes,
blending both edge endpoints' `phin/phip` toward Sentaurus endpoint QF values,
and reconstructing `n/p` from Vela inferred-ni while leaving `psi` unchanged.
This is deliberately stricter than the previous common-mode smooth QF seed: it
moves the active-edge QF-gradient/current-density pattern as well as absolute
carrier density.

Additional replay support was added to
`scripts/diagnose_pn2d_bv_active_edge_mixed_state_replay.py` via
`--state-csv-variant NAME=CSV`, allowing explicit predictor restart states to be
measured directly at the active-edge flux/source level.

Real `-20 V` artifacts:

- False-negative predictor state:
  `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/state_m20_false_negative_blend1/`.
- False-negative single-point restart:
  `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/single_m20_false_negative_blend1/`.
- False-negative active-edge replay:
  `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/active_edge_replay_false_negative_blend1/`.
- All-support predictor state/restart/replay:
  `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/state_m20_all_support_blend1/`,
  `.../single_m20_all_support_blend1/`, and
  `.../active_edge_replay_all_support_blend1/`.

Findings:

1. False-negative-only predictor: `40` active edges, `60` endpoint nodes,
   density factors `n=7.13x`, `p=8.12x`, and QF-gradient deltas around
   `-0.81 mV` electron / `-0.33 mV` hole. Initial replay is source-effective:
   false-negative generation `1.0126x` Sentaurus and particle flux `0.9826x`.
2. The false-negative restart converges in two Newton iterations to
   `-1.169034089095e-16 A/um`, effectively identical to the previous low-current
   branch. Endpoint relaxation retains only `0.2216x` electron QF shift and
   `0.2462x` hole QF shift on false-negative support.
3. All-support predictor: `66` active edges, `92` endpoint nodes. Initial replay
   restores both support classes near Sentaurus (`false_negative` generation
   `1.0126x`, `false_positive` generation `1.0224x`), but the restart still
   lands at `-1.169034109498e-16 A/um`. Retained QF shift remains only about
   `0.22..0.26x` across support and edge endpoints.

Conclusion: moving the active-edge endpoint QF-gradient/current-density pattern
is source-effective before Newton, but still not a sufficient continuation
predictor. The coupled solve actively rolls back most of the desired QF-density
branch in the first two Newton iterations.

Next task:

1. Run a first-Newton-step audit from the coupled predictor state with residual
   block and update block decomposition.
2. Determine whether rollback is driven by carrier-continuity rows, Poisson
   coupling through `psi`, Dirichlet/gauge constraints, or line-search direction
   selection.
3. Only after that, test a predictor that constrains or co-moves the identified
   rollback direction; do not add stronger local density seeds.

## Predictor First-Newton-Step Audit

The first-step rollback audit has been executed with a new helper,
`scripts/diagnose_pn2d_bv_predictor_first_step_audit.py`. The helper can either
analyze existing probe CSVs or prepare/run `newton_step_probe` and
`newton_block_step_probe` from a predictor state CSV. It groups target nodes by
`support_class` and reports how much of the intended `psi-phin` and `phip-psi`
branch shift is rolled back by each step mode.

Artifacts:

- False-negative predictor first-step audit:
  `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/first_step_false_negative_blend1/`.
- All-support predictor first-step audit:
  `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/first_step_all_support_blend1/`.

Findings:

1. False-negative-only predictor: first full Newton step rolls back `0.438x` of
   the electron branch shift and `0.418x` of the hole branch shift on
   false-negative support. Carrier-only gives the same values, while Poisson-only
   is effectively zero rollback (`~2e-5`).
2. All-support predictor: same pattern. False-negative support full-step
   rollback is `0.438x/0.418x`; false-positive support rollback is
   `0.437x/0.425x`. Carrier-only matches full Newton, and Poisson-only remains
   near zero.
3. The predictor initial residual is Poisson dominated (`psi = 0.2443387`), but
   the QF rollback direction comes from the carrier block, not from Poisson/gauge
   motion. The carrier residuals are already small (`~1e-12`) and the carrier
   step still moves QF back toward the low-density branch to reduce them further.

Conclusion: the 75-80% two-iteration QF-density rollback begins with an
approximately 42-44% first carrier-block rollback. The blocker is now localized
to local carrier-continuity row/Jacobian signs and coefficients on active
high-field endpoints, not Poisson consistency, Dirichlet constraints, or a
missing source-effective predictor state.

Next task:

1. Probe local active endpoint carrier rows for the coupled predictor state:
   diagonal/off-diagonal Jacobian terms, RHS sign, and neighboring edge flux
   contributions.
2. Compare the same rows between baseline, predictor initial state, and first
   trial state to identify which term asks Newton to undo the QF-density branch.
3. Only after that, test a constrained diagnostic predictor that neutralizes the
   identified carrier-row rollback term.

## Predictor Carrier-Row/Term Audit

The active-endpoint carrier-row audit has been executed. A new helper,
`scripts/diagnose_pn2d_bv_predictor_carrier_row_audit.py`, can either analyze
existing carrier row/term CSVs or prepare and run `newton_carrier_row_probe` plus
`newton_carrier_term_probe` for labeled states. For the real `-20 V` predictor
cases it compared:

- Baseline Vela `-20 V` VTK state.
- Coupled predictor initial state CSV.
- First full-Newton trial state reconstructed from `newton_step_probe.csv`.

Artifacts:

- False-negative support audit:
  `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/carrier_row_audit_false_negative_blend1/`.
- All-support audit:
  `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/carrier_row_audit_all_support_blend1/`.

Findings:

1. False-negative-only predictor: predictor-minus-baseline median flux deltas
   are `7.09e-14` electron and `9.68e-14` hole, while the impact term only
   compensates `-8.24e-15` for each carrier. Median residual deltas therefore
   remain positive at `6.31e-14` electron and `8.85e-14` hole, and the raw row
   update rolls back `0.438x/0.418x` of the intended electron/hole branch shift.
2. The first trial state moves in the rollback direction and reduces the flux
   and impact magnitudes, but it does not return to baseline. False-negative
   residual deltas remain positive and the raw row update still represents
   `0.341x/0.336x` rollback.
3. All-support prediction repeats the same mechanism on both support classes.
   False-negative support has flux deltas `2.20e-14`/`2.07e-14`, impact
   compensation `-8.24e-15`, and residual deltas `1.32e-14`/`1.27e-14`.
   False-positive support has flux deltas `1.88e-14`/`1.54e-14`, impact
   compensation `-8.45e-15`, and residual deltas `1.04e-14`/`6.99e-15`.

Conclusion: the first-step rollback is localized to carrier-continuity
flux/source balance. The coupled predictor creates the missing active-edge
source, but it also increases the SG flux term more than the impact source term
cancels in Vela's carrier residual, so the local carrier row asks Newton to move
back toward the low-density branch. This is not a Poisson/gauge problem and not
a failure to seed the desired QF branch.

Next task:

1. Run an impact-source feedback sensitivity at the same row/term level:
   compare carrier residual and raw rollback under scaled electron/hole impact
   source terms for baseline, predictor, and first trial states.
2. If impact scaling can neutralize the rollback, inspect production impact
   sign, units, state averaging, and Jacobian coupling against Sentaurus/Charon
   conventions before changing predictor logic.
3. If impact scaling cannot neutralize it, move to SG flux-balance diagnostics
   on the active support edges and compare edge-level continuity terms against
   the Sentaurus reconstructed currents.

## Impact-Source Feedback Sensitivity

The next diagnostic has been executed. `scripts/diagnose_pn2d_bv_predictor_carrier_row_audit.py`
now accepts repeated `--impact-scale` values and writes
`predictor_carrier_impact_scale_nodes.csv` plus an `impact_scale_sensitivity`
summary. The calculation is local and non-invasive: it reuses the row/term probe
CSV, keeps SG flux and state fixed, and evaluates
`adjusted_residual = term_sum + (scale - 1) * impact` for each carrier row. It
also reports the exact per-row multiplier that would close each electron or hole
residual.

Artifacts:

- False-negative-only sensitivity:
  `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/impact_scale_sensitivity_false_negative_blend1/`.
- All-support sensitivity:
  `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/impact_scale_sensitivity_all_support_blend1/`.

Findings:

1. Baseline rows close at scale `~1.0`, as expected, which validates the
   sensitivity calculation against the already-converged Vela `-20 V` branch.
2. False-negative-only predictor rows need much stronger impact feedback:
   required median scale is `7.97x` electron and `10.26x` hole; the first trial
   still needs `8.47x`/`10.43x`. This single-support state restores local source
   magnitude but leaves the continuity rows badly over-fluxed.
3. All-support predictor rows are closer to balance but still under-coupled:
   predictor false-negative support needs `2.31x`/`2.37x`, and predictor
   false-positive support needs `2.07x`/`1.71x`. The first trial remains in the
   same broad range (`~1.61..2.16x`).

Conclusion: impact feedback scaling can neutralize the rollback at the local
carrier-row residual level, especially for the all-support predictor where the
needed multiplier is roughly two rather than ten. This does not justify a source
multiplier yet; it identifies the next root-cause axis as Vela's impact feedback
path relative to Sentaurus/Charon/Genius conventions.

Next task:

1. Audit impact source assembly semantics against references: carrier-current
   weighting, electron/hole source signs, finite-volume support, area/volume
   scaling, and whether source is tied to edge or node carrier continuity rows.
2. Compare Vela's active-edge impact integral and carrier-row source injection
   against Sentaurus reconstructed current/source terms for the all-support
   predictor, where the required multiplier is only about `2x`.
3. Only if the semantic audit supports it, test a focused production change in
   impact feedback/Jacobian coupling. Otherwise move to SG flux-balance support
   diagnostics.
## Impact Feedback Semantic Audit

The production-facing impact feedback audit has now been executed. The local
reference-code check found no sign reversal in the basic generation semantics:
Genius DDM forms `GII = alpha_n * |Jn| / q + alpha_p * |Jp| / q`, distributes it
with directional current weights, and injects the same generated pair source into
electron and hole continuity rows over the finite-volume truncated partial
volume. Charon computes `alpha_n * |Je| + alpha_p * |Jh|` and subtracts avalanche
from total recombination, i.e. treats it as generation. Vela likewise subtracts
the SG edge-current avalanche source from both carrier residual rows. The local
`devsim` tree did not expose an equivalent built-in impact-ionization assembly to
compare.

The concrete all-support predictor comparison is now:

- Active-edge replay: predictor false-negative generation is `1.01257x`
  Sentaurus and false-positive generation is `1.02244x` Sentaurus on the selected
  active x-edges.
- Source geometry replay, after updating
  `scripts/diagnose_pn2d_bv_source_geometry.py` to accept current C++
  `sg_avalanche_edges.csv` columns (`edge_area_m2`, `electron_flux_abs`,
  `hole_flux_abs`, `source_integral`), reports active endpoint area fraction
  `0.5` for both support classes.
- Multiplying those two facts gives effective active-support feedback of
  `0.506x` on false-negative and `0.511x` on false-positive nodes, matching the
  carrier-row sensitivity that needed roughly `2.31x/2.37x` and `2.07x/1.71x`
  electron/hole impact feedback to close the rows.

Artifacts:

- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/source_geometry_all_support_blend1/`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/active_edge_replay_all_support_blend1/`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/impact_scale_sensitivity_all_support_blend1/`

Conclusion: the all-support predictor's local active-edge avalanche physics is
already Sentaurus-sized; the remaining factor-of-two is explained by effective
finite-volume/support feedback into the carrier rows, not by a raw avalanche
coefficient sign or local QF branch-strength error. Do not introduce a hidden
source multiplier. The next task should inspect whether Vela's endpoint
`0.5 * edge_area` carrier-row injection should be compared to Genius/Sentaurus
truncated partial-volume ownership on active edges, and separately whether the
SG edge-current avalanche Jacobian should include source derivatives rather than
omitting them in the current Newton path.
## Impact Feedback Ownership Policy Probe

A follow-up ownership summary has been executed with the new helper
`scripts/summarize_pn2d_bv_impact_feedback_ownership.py`. The helper joins three
already generated all-support predictor diagnostics:

- active-edge replay generation ratio versus Sentaurus,
- source-geometry active endpoint area fraction,
- carrier-row impact scale needed to close the local residual.

Artifact:

- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_frozen_mobility_jacobian_feedback/coupled_qf_predictor/impact_feedback_ownership_all_support_blend1/`

Real `-20 V` all-support results:

| support class | active-edge generation / Sentaurus | active endpoint area fraction | endpoint feedback / Sentaurus | full active-edge feedback / Sentaurus | required e/h impact scale |
|---|---:|---:|---:|---:|---:|
| false_negative | `1.01257x` | `0.5` | `0.50628x` | `1.01257x` | `2.30979x / 2.36797x` |
| false_positive | `1.02244x` | `0.5` | `0.51122x` | `1.02244x` | `2.07195x / 1.71073x` |

This makes the factor-of-two failure mode concrete: the local active-edge source
strength is already approximately Sentaurus-sized, while endpoint-half ownership
leaves the carrier row with only about half of that active feedback. The product
`endpoint_feedback * required_scale` lands near unity (`1.17/1.20` for
false-negative and `1.06/0.87` for false-positive), so the previous `~2x`
impact-scale sensitivity is consistent with source-support ownership rather than
a raw ionization coefficient error.

Current C++ status: the SG edge-current avalanche residual path injects a source,
but the analytic Jacobian path still explicitly omits those nonlocal edge-source
derivatives. That is a separate Newton-coupling issue. The next production probe
should therefore be ordered as:

1. First test a focused source-ownership variant for SG edge-current avalanche,
   comparing Vela endpoint-half ownership against a directional/truncated-volume
   policy analogous to Genius/Sentaurus active-edge ownership. This must remain a
   gated experiment, not a hidden multiplier.
2. Then test source-derivative Jacobian completion for the SG edge-current path
   against finite-difference block probes.
3. Accept neither change without the curve-level `-10..-20 V` knee-shape gate and
   local carrier-row residual gate moving in the same direction.
## SG Edge-Box Source Volume Probe

Executed the focused production probe proposed above by adding an explicit gated
configuration knob, `impact_ionization.source_volume_policy`. The default remains
`edge_half_box`, preserving the previous SG edge-current avalanche source
support. The probe value `edge_box` changes only the SG edge source support from
`0.5 * h * edge.couple` to `1.0 * h * edge.couple`; it is not a hidden avalanche
coefficient multiplier and does not add the still-missing SG source-derivative
Jacobian terms.

Real `-20 V` all-support single-point deck:

- Baseline config: `single_m20_all_support_blend1/simulation_bv_coupled_qf_predictor_all_support_blend1_single_m20.json`
- Probe config: `single_m20_all_support_edge_box/simulation_bv_coupled_qf_predictor_all_support_edge_box_single_m20.json`
- Sentaurus reference current at `-20 V`: `-9.10455666344e-16 A`

| case | converged | Newton iterations | current_total_A_per_um | abs ratio vs baseline | abs ratio vs Sentaurus | decade error vs Sentaurus |
|---|---:|---:|---:|---:|---:|---:|
| endpoint-half baseline | `1` | `2` | `-1.16903410949798e-16` | `1.0000x` | `0.128401x` | `-0.891432` |
| `edge_box` probe | `1` | `2` | `-6.39910455999440e-16` | `5.47384x` | `0.702846x` | `-0.153140` |

The max field is effectively unchanged (`560748.547233267` to
`560748.547233097 V/cm`), so the terminal-current movement is carrier-continuity
feedback from the source ownership policy rather than a field-state movement. The
probe closes about `0.738` decades of the `0.891` decade gap and leaves a
remaining Sentaurus multiplier of `1.42279x`.

Conclusion: source ownership is now confirmed as a production-relevant direction,
but the `edge_box` probe alone is not yet an acceptance change. The next ordered
work is: run the same gated policy through the `-10..-20 V` knee-shape gate and
local carrier-row residual gate, then independently test SG edge-current
avalanche source-derivative Jacobian completion against finite-difference block
probes.

## SG Edge-Box Backscan Knee And Carrier-Row Gate

The `source_volume_policy=edge_box` experiment was extended from a single `-20 V`
point to a matched `-20 -> -10 V` all-support predictor backscan. Two decks were
generated from the same all-support restart family:

- endpoint-half baseline: `all_support_blend1_backscan/simulation_bv_coupled_qf_predictor_all_support_blend1_backscan_m20_to_m10.json`
- edge-box probe: `all_support_edge_box_backscan/simulation_bv_coupled_qf_predictor_all_support_edge_box_backscan_m20_to_m10.json`

Both backscans converged with `11` integer-bias points. The knee-shape gate now
reports:

| curve | first 1 V growth ratio > 1.5 | first 1 V growth ratio > 2.0 | max abs log10 current error |
|---|---:|---:|---:|
| Sentaurus | `-19.0 V` | `-20.0 V` | n/a |
| endpoint-half baseline backscan | none | none | `0.891432` decades |
| `edge_box` backscan | `-16.0 V` | `-19.0 V` | `0.519580` decades |

The probe therefore moves the curve-level knee in the right direction but does
not pass acceptance: the `>2.0` growth point is still one volt early, and the
`>1.5` threshold appears too early at `-16 V`. Pointwise, `edge_box` is close at
some mid-window biases but overshoots around `-19 V` (`3.308x` Sentaurus) while
remaining low at `-20 V` (`0.703x` Sentaurus).

A local carrier-row audit was also run under the `edge_box` policy with the
matched `-20 V` states:

- `carrier_row_audit_all_support_edge_box_policy/`

On the active support classes, the endpoint-half state evaluated with `edge_box`
already has required impact scales below unity (`false_negative` e/h
`0.8898/0.8774`, `false_positive` e/h `0.7784/0.7033`). The converged `edge_box`
state pushes the same rows further source-strong (`false_negative` e/h
`0.6452/0.5778`, `false_positive` e/h `0.6568/0.4142`) and lowers carrier raw
rollback magnitudes. This confirms source ownership is a real lever, but the
full-edge policy is too coarse as a production acceptance change.

Next ordered task: introduce and test an intermediate/truncated ownership factor
rather than binary `0.5` versus `1.0`, or move to the independent SG source-
derivative Jacobian probe if the goal is Newton coupling rather than curve-shape
calibration. Any intermediate policy must be gated by the same backscan knee
summary and carrier-row required-scale summary.
