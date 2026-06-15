# PN2D IV Residual Current Magnitude Debug Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Localize the remaining pn2d Sentaurus2018 IV current mismatch after the high-bias contact quasi-Fermi boundary artifact has been removed.

**Architecture:** Treat the remaining IV error as a curve-shape and current-magnitude calibration problem, not a contact-boundary failure. First freeze the no-relaxation baseline and quantify the error by bias region, then test one axis at a time: contact-current extraction/width convention, mobility and velocity saturation, recombination/effective-ni/BGN, and missing diagnostic observability.

**Tech Stack:** C++20, CMake/Ninja, Catch2, Python regression scripts, Pillow diagnostics, MSYS2 UCRT64 on Windows.

---

## Current Evidence

- The previous high-bias Anode `phin=0` failure was caused by Newton's p-contact minority-electron relaxation.
- `reference_tcad/pn2d_sentaurus2018/pn2d_sentaurus2018_reference.json` now disables IV `contact_boundary_minority_electron_relaxation`.
- Fixed 1.0 V probe converges through `1.0000000000000002 V`.
- Same-bias 1.0 V state fields are now close:
  - `electron_qf_V` max abs diff: `0.015289 V`.
  - `eDensity` mean abs diff: `1.42252e16 cm^-3`.
  - `hDensity` mean abs diff: `1.42250e16 cm^-3`.
  - `ElectricField` p95 abs diff: `332.808 V/cm`.
- IV current remains under-calibrated across the full curve:
  - `0.200 V`: Vela/Sentaurus `1.87215`.
  - `0.250 V`: Vela/Sentaurus `0.998334`.
  - `0.300 V`: Vela/Sentaurus `0.642462`.
  - `0.700 V`: Vela/Sentaurus `0.414665`.
  - `0.800 V`: Vela/Sentaurus `0.556912`.
  - `1.000 V`: Vela/Sentaurus `0.826676`.
- Interpretation: the residual error is not a constant scale factor. The ratio dips in the mid/high forward-bias region and partially recovers by 1 V.

### Baseline ratio buckets

- Diagnostic outputs:
  - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/baseline/iv_ratio_by_bias.csv`
  - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/baseline/iv_ratio_summary.json`
- `low_0p20_0p30`: 3 points, mean Vela/Sentaurus `1.17098`, min ratio `0.642462`, max abs relative error `0.872155`.
- `mid_0p30_0p80`: 10 points, mean Vela/Sentaurus `0.441964`, min ratio `0.405938`, max abs relative error `0.594062`.
- `high_0p80_1p00`: 4 points, mean Vela/Sentaurus `0.756881`, min ratio `0.666872`, max abs relative error `0.333128`.
- Worst relative point: `0.2 V`, Sentaurus `1.67742e-15 A`, Vela scaled `3.14039e-15 A/um`, ratio `1.87215`.
- Interpretation: the largest meaningful current deficit is the `0.3..0.8 V` trough; the `0.2 V` worst relative error is at the near-floor current scale.

### Contact extraction and width audit

- Diagnostic outputs:
  - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/baseline/contact_extraction_audit.csv`
  - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/baseline/contact_extraction_audit.json`
- Terminal-balance current and contact-edge summed current agree to roundoff: max abs edge-minus-terminal `1.35525e-20 A/um`.
- Contact edge count is stable at `17` for the audited terminal rows.
- `current_total / current_total_A_per_um = 1000000.0` at 1.0 V, so the A/m to A/um conversion is exact.
- Decision: do not modify `ContactCurrent` in this phase; the residual IV mismatch is not caused by contact-current extraction or width conversion.

### Mobility and velocity-saturation matrix

- Diagnostic outputs:
  - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/mobility_matrix_light/mobility_matrix_light_summary.csv`
  - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/mobility_matrix_light/mobility_matrix_light_summary.json`
- `baseline_ct_field`: ratios `0.25 V=0.998334`, `0.30 V=0.642462`, `0.80 V=0.556912`, `1.00 V=0.826676`.
- `no_mobility`: improves the mid-bias trough (`0.30 V=0.967282`, `0.80 V=1.00476`) but breaks the full shape (`0.25 V=1.32010`, `1.00 V=1.49128`, max abs relative error `1.18780`).
- `ct_silicon`: nearly identical to baseline (`0.30 V=0.644634`, `0.80 V=0.558488`, `1.00 V=0.828427`).
- `ct_bv_constants` and `ct_field_bv_constants`: slightly reduce the already-low mid/high ratios (`0.30 V=0.640575/0.638436`, `0.80 V=0.552825/0.551283`, `1.00 V=0.819985/0.818267`).
- `ct_field_silicon` is identical to `baseline_ct_field`.
- A script-level regression was added for `scripts/scan_pn2d_iv_mobility_candidates.py` so relative `materials_file` paths resolve correctly when candidate configs are copied to a temporary output directory.
- Decision: do not promote a mobility candidate. Removing mobility fixes part of the trough only by over-scaling the rest of the IV curve, while CT parameter variants do not explain the residual shape mismatch as a single clean axis.

### Recombination, effective-ni, and BGN matrix

- Diagnostic outputs:
  - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/recombination_matrix_light/recombination_matrix_light_summary.csv`
  - `build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/recombination_matrix_light/recombination_matrix_light_summary.json`
- `srh_default`: ratios `0.25 V=0.998334`, `0.30 V=0.642462`, `0.80 V=0.556912`, `1.00 V=0.826676`.
- `recomb_none`: destroys the low-bias near match (`0.25 V=0.396680`, `0.30 V=0.400454`) while leaving `0.80 V=0.556422` and `1.00 V=0.827723` essentially unchanged.
- `bgn_none`: worsens the full curve (`0.25 V=0.776103`, `0.30 V=0.501320`, `0.80 V=0.464690`, `1.00 V=0.782650`).
- `tau_1e-8`: improves neither full shape nor low-bias calibration; it overdrives low bias (`0.25 V=6.41314`, `0.30 V=2.82048`) and only nudges `0.80 V` to `0.561263`.
- `tau_3e-8`: also overdrives low bias (`0.25 V=2.40219`, `0.30 V=1.20714`) with little improvement at `0.80 V=0.558051`.
- `tau_1e-6`: damages low bias (`0.25 V=0.456844`, `0.30 V=0.424654`) and leaves high bias nearly unchanged (`0.80 V=0.556471`, `1.00 V=0.827618`).
- Decision: do not promote recombination, lifetime, effective-ni, or BGN tuning. These candidates either leave the mid/high trough essentially unchanged or break the already-good `0.25 V` point, so the plan stops here per the Task 4 stop condition.

## Files and Responsibilities

- `reference_tcad/pn2d_sentaurus2018/pn2d_sentaurus2018_reference.json`: faithful pn2d Sentaurus2018 deck overrides and comparison gate.
- `scripts/compare_reference_curves.py`: existing curve comparison gate; keep it as pass/fail infrastructure.
- `scripts/scan_pn2d_iv_mobility_candidates.py`: existing mobility scan entry point; reuse before writing new scanner code.
- `scripts/compare_pn2d_0v_current_related_quantities.py`: reusable field/current CSV and plotting helpers; borrow parsing utilities for IV residual reports.
- `src/post/ContactCurrent.cpp`: contact current extraction and detailed edge diagnostics.
- `src/simulation/DCSweep.cpp`: terminal-balance, contact-edge diagnostics, and optional per-bias diagnostic output.
- `src/physics/MobilityModel.cpp` and related headers: mobility/current-magnitude model axis.
- `src/physics/Recombination.cpp` and `src/physics/BandgapNarrowing.cpp`: recombination/effective-ni/BGN model axis.
- `tests/regression/test_reference_tcad_tools.py`: Python coverage for new scripts and reference config expectations.
- `tests/test_dc_sweep.cpp`, `tests/test_contact_current.cpp` if present, or `tests/test_sg_flux.cpp`: C++ coverage if diagnostics or current extraction behavior changes.
- Generated outputs under `build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/`: diagnostic artifacts only; do not commit.

---

### Task 1: Freeze the Fixed IV Baseline and Bias-Bucket Metrics

**Files:**
- Create: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/baseline/iv_ratio_by_bias.csv`
- Create: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/baseline/iv_ratio_summary.json`
- Read: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/fixed/iv_1v_fixed_probe.csv`
- Read: `build/reference_tcad/pn2d_sentaurus2018/reference_curves/pn2d_sentaurus2018_iv_reference.csv`

- [x] **Step 1: Generate a fixed-baseline ratio table**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
@'
import csv, json, math
from pathlib import Path

root = Path("build/reference_tcad/pn2d_sentaurus2018")
out = root / "reports" / "iv_residual_current" / "baseline"
out.mkdir(parents=True, exist_ok=True)
ref = [
    (float(r["bias_V"]), float(r["current_total"]))
    for r in csv.DictReader((root / "reference_curves" / "pn2d_sentaurus2018_iv_reference.csv").open())
]
ref.sort()
def interp(x):
    if x >= ref[-1][0]:
        return ref[-1][1]
    lo = None
    for b, c in ref:
        if b <= x:
            lo = (b, c)
        if b >= x:
            hi = (b, c)
            break
    if lo is None:
        return ref[0][1]
    if abs(hi[0] - lo[0]) < 1e-18:
        return lo[1]
    return lo[1] + (x - lo[0]) * (hi[1] - lo[1]) / (hi[0] - lo[0])

rows = []
candidate = root / "reports" / "iv_state" / "fixed" / "iv_1v_fixed_probe.csv"
for r in csv.DictReader(candidate.open()):
    if r["current_contact"] != "Cathode" or r["converged"] != "1":
        continue
    bias = float(r["bias_V"])
    if bias < 0.2:
        continue
    sent = interp(bias)
    vela = -float(r["current_total_A_per_um"])
    rows.append({
        "bias_V": bias,
        "sentaurus_A": sent,
        "vela_scaled_A_per_um": vela,
        "ratio_vela_to_sentaurus": vela / sent if sent else math.nan,
        "relative_error": (vela - sent) / sent if sent else math.nan,
    })

with (out / "iv_ratio_by_bias.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

buckets = {
    "low_0p20_0p30": [r for r in rows if 0.2 <= r["bias_V"] <= 0.3],
    "mid_0p30_0p80": [r for r in rows if 0.3 < r["bias_V"] <= 0.8],
    "high_0p80_1p00": [r for r in rows if 0.8 < r["bias_V"] <= 1.0 + 1e-12],
}
summary = {
    name: {
        "points": len(items),
        "mean_ratio": sum(r["ratio_vela_to_sentaurus"] for r in items) / len(items),
        "min_ratio": min(r["ratio_vela_to_sentaurus"] for r in items),
        "max_abs_rel_error": max(abs(r["relative_error"]) for r in items),
    }
    for name, items in buckets.items()
    if items
}
worst = max(rows, key=lambda r: abs(r["relative_error"]))
summary["worst"] = worst
(out / "iv_ratio_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(out / "iv_ratio_summary.json")
print(summary)
'@ | python -
```

Expected: a stable baseline shows the ratio dip is strongest in the mid/high region, while `0.25 V` is nearly matched.

- [x] **Step 2: Add the baseline numbers to the debug notes**

Append to `docs/superpowers/plans/2026-06-13-pn2d-iv-residual-current-magnitude-debug.md`:

```markdown
### Baseline ratio buckets

- low_0p20_0p30:
- mid_0p30_0p80:
- high_0p80_1p00:
- worst:
```

Expected: every later candidate is compared against these fixed bucket metrics.

---

### Task 2: Audit Contact-Current Extraction and Width Convention Across Bias

**Files:**
- Read: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/fixed/iv_1v_fixed_terminal_balance.csv`
- Read: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/fixed/iv_1v_fixed_contact_edges.csv`
- Create: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/baseline/contact_extraction_audit.csv`
- Create: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/baseline/contact_extraction_audit.json`

- [x] **Step 1: Compare terminal current and edge-sum current at every fixed bias**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
@'
import csv, json
from collections import defaultdict
from pathlib import Path

root = Path("build/reference_tcad/pn2d_sentaurus2018")
base = root / "reports" / "iv_state" / "fixed"
out = root / "reports" / "iv_residual_current" / "baseline"
out.mkdir(parents=True, exist_ok=True)

terminal = {}
for r in csv.DictReader((base / "iv_1v_fixed_terminal_balance.csv").open()):
    terminal[(r["bias_V"], r["contact"])] = float(r["current_total_A_per_um"])

edge_sum = defaultdict(float)
edge_abs = defaultdict(float)
edge_count = defaultdict(int)
for r in csv.DictReader((base / "iv_1v_fixed_contact_edges.csv").open()):
    key = (r["bias_V"], r["current_contact"])
    value = float(r["current_total_A_per_um"])
    edge_sum[key] += value
    edge_abs[key] += abs(value)
    edge_count[key] += 1

rows = []
for key, term in sorted(terminal.items(), key=lambda item: (float(item[0][0]), item[0][1])):
    esum = edge_sum.get(key, 0.0)
    rows.append({
        "bias_V": key[0],
        "contact": key[1],
        "terminal_total_A_per_um": term,
        "edge_sum_A_per_um": esum,
        "edge_abs_sum_A_per_um": edge_abs.get(key, 0.0),
        "edge_count": edge_count.get(key, 0),
        "edge_minus_terminal_A_per_um": esum - term,
        "edge_to_terminal_ratio": esum / term if term else "",
    })

with (out / "contact_extraction_audit.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)

max_abs_delta = max(abs(r["edge_minus_terminal_A_per_um"]) for r in rows)
summary = {
    "rows": len(rows),
    "max_abs_edge_minus_terminal_A_per_um": max_abs_delta,
    "max_edge_count": max(r["edge_count"] for r in rows),
}
(out / "contact_extraction_audit.json").write_text(json.dumps(summary, indent=2) + "\n")
print(summary)
'@ | python -
```

Expected: edge sums match terminal totals to numerical roundoff. If they do not, fix `ContactCurrent` extraction before any physics calibration.

- [x] **Step 2: Check the exact `A/m` to `A/um` conversion**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
@'
import csv
from pathlib import Path
p = Path("build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/fixed/iv_1v_fixed_probe.csv")
for r in csv.DictReader(p.open()):
    if r["current_contact"] == "Cathode" and r["converged"] == "1" and abs(float(r["bias_V"]) - 1.0) < 1e-9:
        per_m = float(r["current_total"])
        per_um = float(r["current_total_A_per_um"])
        print("current_total/current_total_A_per_um =", per_m / per_um)
'@ | python -
```

Expected: ratio is exactly `1e6`. If not, stop and debug unit conversion.

- [x] **Step 3: Decision gate**

Use this decision:

```text
If edge sums and A/m-to-A/um conversion are consistent, do not modify ContactCurrent in this phase.
If edge sums disagree with terminal totals, add a failing C++ test in tests/test_dc_sweep.cpp or tests/test_contact_current.cpp before changing ContactCurrent.
```

Expected: the remaining mismatch is either cleared as not extraction-related or pinned to a concrete current-extraction bug.

---

### Task 3: Run a No-Relaxation Mobility and Velocity-Saturation Matrix

**Files:**
- Read: `scripts/scan_pn2d_iv_mobility_candidates.py`
- Read: `reference_tcad/pn2d_sentaurus2018/pn2d_sentaurus2018_reference.json`
- Create: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/mobility_matrix/`

- [x] **Step 1: Inspect the existing mobility scan interface**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts/scan_pn2d_iv_mobility_candidates.py --help
```

Expected: the script exposes enough parameters to regenerate pn2d candidates. If it does not expose no-relaxation, add a solver override in the script before scanning.

- [x] **Step 2: Run the existing candidate matrix with the fixed no-relaxation IV config**

Run the scan into a fresh output directory:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts/scan_pn2d_iv_mobility_candidates.py --output-dir build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/mobility_matrix
```

Expected: a candidate report ranks mobility options by IV window and high-bias ratio without reintroducing contact relaxation.

- [x] **Step 3: If the existing scan cannot produce the needed fixed candidates, add a script-level test first**

Add this test to `tests/regression/test_reference_tcad_tools.py`:

```python
def test_pn2d_iv_mobility_scan_preserves_disabled_minority_relaxation(self) -> None:
    path = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "pn2d_sentaurus2018_reference.json"
    config = json.loads(path.read_text())
    iv = next(sim for sim in config["simulations"] if sim["name"] == "iv")
    self.assertIs(iv["vela_solver"]["contact_boundary_minority_electron_relaxation"], False)
```

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_pn2d_iv_mobility_scan_preserves_disabled_minority_relaxation
```

Expected: FAIL only if the test name is new and not yet added, then PASS after adding it.

- [x] **Step 4: Decision gate**

Use this decision:

```text
If one mobility candidate improves 0.3..1.0 V ratios without breaking 0.25 V, promote it to the IV reference config and add regression coverage.
If all mobility candidates shift the whole curve but do not fix the dip shape, keep current mobility and move to recombination/effective-ni.
```

Expected: mobility is either selected as the next fix axis or ruled out as a single-factor explanation.

---

### Task 4: Run Recombination, Effective-ni, and BGN Shape Matrix

**Files:**
- Create: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/recombination_matrix/`
- Candidate modify only if promoted: `reference_tcad/pn2d_sentaurus2018/pn2d_sentaurus2018_reference.json`
- Test if promoted: `tests/regression/test_reference_tcad_tools.py`

- [x] **Step 1: Create temporary IV configs for physics toggles**

Generate candidates from the fixed 1V probe config:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
@'
import json
from pathlib import Path

base = Path("build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/fixed/simulation_iv_1v_fixed_probe.json")
out = Path("build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/recombination_matrix")
out.mkdir(parents=True, exist_ok=True)
base_cfg = json.loads(base.read_text())
candidates = {
    "srh_default": {},
    "recomb_none": {"recombination": []},
    "bgn_none": {"bandgap_narrowing": "none"},
    "tau_1e-8": {"taun": 1.0e-8, "taup": 1.0e-8},
    "tau_3e-8": {"taun": 3.0e-8, "taup": 3.0e-8},
    "tau_1e-6": {"taun": 1.0e-6, "taup": 1.0e-6},
}
for name, solver_updates in candidates.items():
    cfg = json.loads(json.dumps(base_cfg))
    cfg["output_csv"] = str((out / f"{name}.csv").resolve())
    cfg["sweep"]["vtk_prefix"] = str((out / name).resolve())
    cfg["sweep"]["diagnostics"]["terminal_balance"]["csv"] = str((out / f"{name}_terminal_balance.csv").resolve())
    cfg["sweep"]["diagnostics"]["contact_edge"]["csv"] = str((out / f"{name}_contact_edges.csv").resolve())
    cfg.setdefault("solver", {}).update(solver_updates)
    (out / f"{name}.json").write_text(json.dumps(cfg, indent=2) + "\n")
    print(out / f"{name}.json")
'@ | python -
```

Expected: six candidate JSON configs are written, all preserving `contact_boundary_minority_electron_relaxation=false`.

- [x] **Step 2: Run every candidate**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
@'
from pathlib import Path
import subprocess

root = Path("build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/recombination_matrix")
runner = Path("build/vela_example_runner.exe")
for cfg in sorted(root.glob("*.json")):
    print("RUN", cfg)
    result = subprocess.run([str(runner), "--config", str(cfg)], text=True)
    print("RETURN", result.returncode)
'@ | python -
```

Expected: candidates either converge through 1 V or record a specific failure row. Do not promote a candidate that reintroduces validation failure.

- [x] **Step 3: Score candidates with the same ratio buckets**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
@'
import csv, json, math
from pathlib import Path

root = Path("build/reference_tcad/pn2d_sentaurus2018")
matrix = root / "reports" / "iv_residual_current" / "recombination_matrix"
ref = [
    (float(r["bias_V"]), float(r["current_total"]))
    for r in csv.DictReader((root / "reference_curves" / "pn2d_sentaurus2018_iv_reference.csv").open())
]
ref.sort()
def interp(x):
    if x >= ref[-1][0]:
        return ref[-1][1]
    lo = None
    for b, c in ref:
        if b <= x:
            lo = (b, c)
        if b >= x:
            hi = (b, c)
            break
    if abs(hi[0] - lo[0]) < 1e-18:
        return lo[1]
    return lo[1] + (x - lo[0]) * (hi[1] - lo[1]) / (hi[0] - lo[0])

summary = []
for csv_path in sorted(matrix.glob("*.csv")):
    if csv_path.name.endswith("_terminal_balance.csv") or csv_path.name.endswith("_contact_edges.csv"):
        continue
    rows = []
    for r in csv.DictReader(csv_path.open()):
        if r["current_contact"] == "Cathode" and r["converged"] == "1" and float(r["bias_V"]) >= 0.2:
            b = float(r["bias_V"])
            s = interp(b)
            v = -float(r["current_total_A_per_um"])
            rows.append((b, v / s, (v - s) / s))
    if not rows:
        summary.append({"candidate": csv_path.stem, "status": "no_converged_rows"})
        continue
    ratios = {round(b, 6): ratio for b, ratio, _ in rows}
    summary.append({
        "candidate": csv_path.stem,
        "points": len(rows),
        "ratio_0p3": ratios.get(0.3),
        "ratio_0p8": ratios.get(0.8),
        "ratio_1p0": ratios.get(1.0),
        "max_abs_rel_error": max(abs(rel) for _, _, rel in rows),
    })
(matrix / "candidate_score.json").write_text(json.dumps(summary, indent=2) + "\n")
for item in summary:
    print(item)
'@ | python -
```

Expected: identify whether recombination/lifetime/BGN changes flatten the ratio dip without destroying the 0.25 V match.

- [x] **Step 4: Decision gate**

Use this decision:

```text
If no recombination or BGN candidate improves both 0.3 V and 0.8 V ratios, do not tune SRH/ni/BGN in this phase.
If one candidate improves curve shape and keeps 1 V fields close, promote it with a config regression test.
```

Expected: recombination/effective-ni is either selected as a fix axis or ruled out.

---

### Task 5: Add Vela IV Transport Diagnostics if Mobility/Recombination Remain Ambiguous

**Files:**
- Modify: `src/simulation/DCSweep.cpp`
- Modify: `src/post/ContactCurrent.cpp`
- Modify: `include/vela/post/ContactCurrent.h`
- Test: `tests/test_dc_sweep.cpp`
- Create: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_residual_current/transport_diagnostics/`

- [ ] **Step 1: Write a failing test for optional transport diagnostics columns**

Add to `tests/test_dc_sweep.cpp`:

```cpp
TEST_CASE("DCSweep: transport diagnostics append mobility and current driver columns", "[dc_sweep][diagnostics]")
{
    auto config = makePnDiodeSweepConfig();
    config["solver"]["method"] = "gummel_newton";
    config["solver"]["handoff"] = {
        {"fallback", "none"},
        {"require_gummel_convergence", false},
        {"gummel_max_iter", 0},
        {"newton_max_iter", 4},
    };
    config["sweep"]["diagnostics"]["transport"] = {
        {"enabled", true},
    };

    const auto result = runDCSweepFromJson(config);

    REQUIRE(result.points.size() >= 1);
    const auto header = readCsvHeader(config["output_csv"].get<std::string>());
    REQUIRE(std::find(header.begin(), header.end(), "mean_electron_mobility_m2_V_s") != header.end());
    REQUIRE(std::find(header.begin(), header.end(), "mean_hole_mobility_m2_V_s") != header.end());
    REQUIRE(std::find(header.begin(), header.end(), "max_electric_field_V_per_cm") != header.end());
}
```

Expected before implementation: compile fails if helper names differ or test fails because columns are absent. Adjust helper names to existing local test helpers while preserving the column assertions.

- [ ] **Step 2: Implement minimal diagnostic columns**

Add opt-in columns only when `sweep.diagnostics.transport.enabled=true`:

```text
mean_electron_mobility_m2_V_s
mean_hole_mobility_m2_V_s
min_electron_mobility_m2_V_s
min_hole_mobility_m2_V_s
max_electric_field_V_per_cm
mean_srh_recombination_cm3_s
```

Expected: existing CSV schema is unchanged unless diagnostics are enabled.

- [ ] **Step 3: Run diagnostic-enabled fixed IV probe**

Run a derived fixed IV config with transport diagnostics enabled:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
build\vela_example_runner.exe --config build\reference_tcad\pn2d_sentaurus2018\reports\iv_residual_current\transport_diagnostics\simulation_iv_transport_probe.json
```

Expected: CSV contains mobility/current-driver columns that can be plotted against the ratio dip.

- [ ] **Step 4: Decision gate**

Use this decision:

```text
If mobility collapses in the same bias region where current ratio dips, investigate high-field mobility/velocity saturation.
If mobility is smooth but SRH/recombination changes sharply, investigate recombination/effective-ni.
If both are smooth, revisit current-density integration geometry or missing Fermi-Dirac statistics.
```

Expected: the next implementation target is selected from measured Vela transport drivers.

---

### Task 6: Promote One Fix Axis or Document No-Promotion

**Files:**
- Candidate modify: `reference_tcad/pn2d_sentaurus2018/pn2d_sentaurus2018_reference.json`
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`
- Test: `tests/regression/test_reference_tcad_tools.py`

- [ ] **Step 1: If a candidate is promoted, add a config regression test**

Add a test like:

```python
def test_pn2d_sentaurus2018_iv_promoted_physics_axis(self) -> None:
    path = REPO / "reference_tcad" / "pn2d_sentaurus2018" / "pn2d_sentaurus2018_reference.json"
    config = json.loads(path.read_text())
    iv = next(sim for sim in config["simulations"] if sim["name"] == "iv")
    solver = iv["vela_solver"]
    self.assertIs(solver["contact_boundary_minority_electron_relaxation"], False)
    self.assertEqual(solver["mobility"]["model"], "caughey_thomas_field")
```

Add extra assertions only for the promoted axis. For example, if a lifetime is promoted:

```python
self.assertEqual(solver["taun"], 3.0e-8)
self.assertEqual(solver["taup"], 3.0e-8)
```

Expected: test fails before config update and passes after config update.

- [ ] **Step 2: Update the validation report**

Update `docs/validation/pn2d_sentaurus_comparison.md` with:

```markdown
### IV residual current magnitude follow-up

- Baseline fixed no-relaxation ratios:
- Extraction audit:
- Mobility matrix result:
- Recombination/effective-ni matrix result:
- Promoted change or no-promotion decision:
- Remaining mismatch:
```

Expected: the report explains whether a new physics parameter was promoted or why no promotion was made.

- [ ] **Step 3: Run focused verification**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R "DCSweep|NewtonSolver|DDSolution validation|reference_tcad_regression|sentaurus_import_tools"
python -m unittest tests.regression.test_reference_tcad_tools tests.regression.test_sentaurus_import_tools
```

Expected: focused tests pass.

- [ ] **Step 4: Run full verification**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
ctest --test-dir build --output-on-failure
```

Expected: full suite passes.

---

## Stop Conditions

- Stop after Task 2 if current extraction or A/m-to-A/um conversion is inconsistent; that is a correctness bug and must be fixed before physics calibration.
- Stop after Task 3 if a mobility candidate clearly improves the full shape and does not regress 0.25 V; promote that single axis before testing recombination.
- Stop after Task 4 if recombination/BGN candidates all worsen the 0.25 V near-match; do not tune lifetimes to hide a transport-shape mismatch.
- Stop after Task 5 if transport diagnostics show missing observability in core solver output; add diagnostics first, then resume model tuning.
