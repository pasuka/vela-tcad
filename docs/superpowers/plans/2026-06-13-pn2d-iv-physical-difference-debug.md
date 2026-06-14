# PN2D IV Physical Difference Debug Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Localize and fix the high-bias IV physical-quantity mismatch between Sentaurus and Vela, starting from the Anode electron quasi-Fermi boundary anomaly.

**Architecture:** Treat the current IV mismatch as a data-flow bug until proven otherwise: sweep bias enters `DCSweep`, becomes Newton/Gummel contact boundary conditions, is assembled as Dirichlet rows, validated, written to VTK/CSV, and then compared against Sentaurus. First reproduce and instrument this path, then make the smallest fix supported by evidence.

**Tech Stack:** C++20, CMake/Ninja, Catch2, Python regression scripts, MSYS2 UCRT64 on Windows.

---

## Current Evidence

- Sentaurus IV field snapshot is at Anode `1.0 V`; Vela IV probe last converged at `0.8265625 V`.
- Vela next point `0.828125 V` failed validation with:

```text
contact 'Anode' node 2 phin=0 does not match bias=0.82812500000000011
```

- Vela terminal current is still balanced at the last converged point:

```text
Anode total  = +3.690157315e-6 A/um
Cathode total = -3.690157320e-6 A/um
```

- Largest physical differences are downstream of the same boundary anomaly:

```text
electron_qf_V max diff = 1.0 V at node 2, Sentaurus=1.0, Vela=0.0
eDensity mean abs diff = 1.1099e17 cm^-3
hDensity mean abs diff = 1.1010e17 cm^-3
ElectricField p95 diff = 4169 V/cm
Vela/Sentaurus current ratio at 0.8265625 V = 0.2279
```

### IV high-bias Anode boundary trace

- Last converged Vela bias: `0.82656250000000009 V`.
- Failed Vela bias: `0.82812500000000011 V`.
- First observed Anode minority-electron `phin` collapse in contact-edge diagnostics: already present at `0.82500000000000007 V`; Anode contact-side edge values show `phin1=0` while `phip1=bias`.
- Edge/node ids involved: example edge `node0=0`, `node1=2`, with `phin0=1.8414e-6`, `phin1=0`, `phip0=0.8246459`, `phip1=0.825`.
- Node 2 at the last converged field is not electron-majority: `eDensity=9.82989e16 cm^-3`, `hDensity=1.0e17 cm^-3`.
- Interpretation: the boundary relaxation/collapse exists before the failed point, but validation should still prefer `phip` at node 2 if it uses the last converged density majority. The failed point likely changes the solved `n/p` majority or validates a transient Newton result where `n > p`.

### IV minority-electron relaxation probe

- Active JSON only sets `contact_boundary_reconstruction=dominant_signed_contact_mean`; Newton defaults still enable `contactBoundaryMinorityElectronRelaxation=true`, threshold `0.1 V`, `p_contact_only`, strength `1.0`.
- Derived probe `simulation_iv_1v_no_minority_relax_probe.json` sets `contact_boundary_minority_electron_relaxation=false`.
- Result: no-relaxation probe converged through `1.0000000000000002 V` with 21 IV points.
- At Anode 1.0 V contact edge, contact-side `phin1=1.0000000000000002` and `phip1=1.0000000000000002`, so the boundary collapse is removed.
- 1.0 V current comparison: Sentaurus `1.25880187856e-4`, Vela scaled `1.0406217981174905e-4`, Vela/Sentaurus `0.826676`, relative error `-0.173324`.
- Interpretation: minority-electron relaxation is the immediate cause of the high-bias `phin=0` validation failure and a dominant contributor to the earlier high-bias current deficit.

### IV Anode node 2 density trace

- Diagnostic CSV: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/physical_quantity_compare/iv_anode_node_density_trace.csv`.
- At node 2 in the last converged relaxed run, Vela has `eDensity=9.82989e16 cm^-3`, `hDensity=1.0e17 cm^-3`, and `NetDoping=-1.0e17 cm^-3`.
- Therefore node 2 is still p-majority in the last converged output; the validator would check `phip`, not `phin`, for that state.
- Interpretation: the original failure is consistent with validating a failed next-step Newton state whose contact density became electron-majority while relaxed minority `phin` stayed at `0 V`.

### IV Sentaurus current-density and recombination drivers

- Diagnostic CSV: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/physical_quantity_compare/iv_sentaurus_current_recombination_summary.csv`.
- Sentaurus 1.0 V current density means: `eCurrentDensity=17494.4`, `hCurrentDensity=7681.58`, `TotalCurrentDensity=25175.98`.
- Sentaurus current split by mean density is roughly electron `69.5%`, hole `30.5%`.
- Vela no-relaxation 1.0 V terminal split at Cathode: electron `-7.11038e-5 A/um`, hole `+3.29584e-5 A/um`, total `-1.04062e-4 A/um`, giving a similar electron/hole split.
- Sentaurus SRH recombination at 1.0 V has mean `1.44817e22 cm^-3 s^-1` and max `2.18777e22 cm^-3 s^-1`.
- Interpretation: after removing the boundary relaxation artifact, the remaining ~17.3% current deficit is likely a whole-current magnitude calibration issue, not a single carrier branch sign/cancellation error. Next suspects are mobility/current-density magnitude, SRH lifetime/effective-ni calibration, or reference width/unit convention.

### Implemented minimal fix

- `reference_tcad/pn2d_sentaurus2018/pn2d_sentaurus2018_reference.json` now sets IV `vela_solver.contact_boundary_minority_electron_relaxation=false`.
- Added C++ tests for disabled-relaxation Newton contact QF boundaries and SolutionValidation majority-carrier contact checks.
- Added Python regression coverage that the pn2d Sentaurus2018 IV reference config keeps minority relaxation disabled.
- Targeted verification passed:
  - `ctest --test-dir build --output-on-failure -R "NewtonSolver|DDSolution validation|DCSweep"`: 58/58 passed.
  - `python -m unittest tests.regression.test_reference_tcad_tools`: 16/16 passed.

### Fixed 1.0 V IV physical-quantity comparison

- Fixed probe config: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/fixed/simulation_iv_1v_fixed_probe.json`.
- Fixed probe result: converged through `1.0000000000000002 V` with 21 points.
- Fixed curve report: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/fixed/iv_curve_compare.json`; status is still `fail` over `0.20547013066..1.0 V` because `max_relative_error=1.17508` and `orders_of_magnitude=0.337475`, but trend matches.
- High-bias current ratios after the fix:
  - `0.8 V`: Vela/Sentaurus `0.556912`.
  - `0.9 V`: Vela/Sentaurus `0.742684`.
  - `1.0 V`: Vela/Sentaurus `0.826676`.
- Same-bias 1.0 V field comparison output: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/fixed/physical_quantity_compare/iv_fixed_1v_physical_quantity_comparison.json`.
- Same-bias 1.0 V field error summary:
  - `potential_V`: mean abs `0.002572 V`, max abs `0.009852 V`.
  - `electron_qf_V`: mean abs `0.007607 V`, max abs `0.015289 V`; original diagnostic max was `1.0 V`.
  - `hole_qf_V`: mean abs `0.007690 V`, max abs `0.015057 V`.
  - `electron_density_cm3`: mean abs `1.42252e16`, max abs `2.2224e16`; original diagnostic mean was `1.10991e17`.
  - `hole_density_cm3`: mean abs `1.42250e16`, max abs `2.2142e16`; original diagnostic mean was `1.10105e17`.
  - `electric_field_V_per_cm`: mean abs `103.536`, p95 `332.808`, max abs `802.031`; original diagnostic p95 was `4169.25`.
- Interpretation: boundary relaxation removal fixes the dominant quasi-Fermi/field/density artifact. Remaining IV curve mismatch is now mostly a current magnitude calibration problem rather than a boundary-condition failure.

### Final verification

- Focused CTest passed: `ctest --test-dir build --output-on-failure -R "Scharfetter|NewtonSolver|DDSolution validation|DCSweep|reference_tcad_regression|sentaurus_import_tools|diagnose_pn2d_0v_current_balance_newton_failure"` -> 61/61 passed.
- Focused Python regression passed: `python -m unittest tests.regression.test_reference_tcad_tools tests.regression.test_sentaurus_import_tools tests.regression.test_run_regression` -> 67/67 passed.
- Full CTest passed: `ctest --test-dir build --output-on-failure` -> 288/288 passed.
- Post-documentation smoke passed: `ctest --test-dir build --output-on-failure -R "ascii_sources|reference_tcad_regression"` -> 2/2 passed.
- Post-documentation Python smoke passed: `python -m unittest tests.regression.test_reference_tcad_tools` -> 16/16 passed.

## Files and Responsibilities

- `src/simulation/DCSweep.cpp`: Applies sweep bias, calls Gummel/Newton, records per-point failure and contact-edge diagnostics.
- `src/solver/NewtonSolver.cpp`: Builds Newton contact boundary condition maps, including `contact_boundary_reconstruction` and minority-electron relaxation.
- `src/equation/CoupledDDAssembler.cpp`: Enforces Dirichlet rows for `psi`, `phin`, and `phip`.
- `src/solver/SolutionValidation.cpp`: Validates contact quasi-Fermi fields against contact bias using majority-carrier logic.
- `tests/test_newton_solver.cpp`: Unit-level tests for Newton contact boundary behavior.
- `tests/test_solution_validation.cpp`: Unit-level tests for validation behavior.
- `tests/test_dc_sweep.cpp`: Sweep-level regression tests for validation and diagnostics.
- `scripts/compare_pn2d_0v_current_related_quantities.py`: Reusable node-field/current comparison helper; extend or wrap for IV-specific output only if needed.
- `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/*`: Generated diagnostic inputs and outputs; do not commit.

---

### Task 1: Reproduce and Pin the Anode `phin=0` Failure

**Files:**
- Read: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/iv_1v_probe.csv`
- Read: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/iv_1v_probe_contact_edges.csv`
- Create: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/physical_quantity_compare/iv_anode_boundary_trace.csv`

- [ ] **Step 1: Extract the last two converged rows and the failed row**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python - <<'PY'
import csv
from pathlib import Path
p = Path("build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/iv_1v_probe.csv")
rows = list(csv.DictReader(p.open()))
for r in rows[-6:]:
    print(r["bias_V"], r["current_contact"], r["converged"], r["failure_reason"], r["validation_diagnostics"])
PY
```

Expected: the failed row at `0.82812500000000011` reports `validation_failed` and Anode `phin=0`.

- [ ] **Step 2: Extract Anode contact-edge quasi-Fermi values at `0.825` and `0.8265625`**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python - <<'PY'
import csv
from pathlib import Path
src = Path("build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/iv_1v_probe_contact_edges.csv")
out = Path("build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/physical_quantity_compare/iv_anode_boundary_trace.csv")
out.parent.mkdir(parents=True, exist_ok=True)
rows = []
for r in csv.DictReader(src.open()):
    if r.get("current_contact") == "Anode" and r.get("bias_V") in {"0.82500000000000007", "0.82656250000000009"}:
        rows.append({k: r.get(k, "") for k in ["bias_V", "contact", "node0", "node1", "psi0", "psi1", "phin0", "phin1", "phip0", "phip1"]})
with out.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
print(out)
print("rows", len(rows))
PY
```

Expected: Anode edge rows show whether `phin` is already collapsing before the failed point, or only at the failed point.

- [ ] **Step 3: Record the evidence in the debug notes**

Append a short note to the active debug summary:

```markdown
### IV high-bias Anode boundary trace

- Last converged Vela bias:
- Failed Vela bias:
- First observed Anode `phin` collapse:
- Edge/node ids involved:
- Interpretation:
```

Expected: the next task starts from observed boundary behavior, not from a guessed fix.

---

### Task 2: Determine Whether Minority-Electron Relaxation Is the Immediate Cause

**Files:**
- Read: `src/solver/NewtonSolver.cpp`
- Read: `include/vela/solver/NewtonSolver.h`
- Create: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/simulation_iv_1v_no_minority_relax_probe.json`

- [ ] **Step 1: Confirm the active Newton contact settings**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python - <<'PY'
import json
from pathlib import Path
p = Path("build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/simulation_iv_1v_probe.json")
cfg = json.loads(p.read_text())
solver = cfg["solver"]
for key in sorted(solver):
    if key.startswith("contact_boundary"):
        print(key, "=", solver[key])
PY
```

Expected: print the configured `contact_boundary_reconstruction` and any minority relaxation settings present in the IV probe deck.

- [ ] **Step 2: Run a one-variable probe with minority-electron relaxation disabled**

Create a derived config and run it:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python - <<'PY'
import json
from pathlib import Path
base = Path("build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/simulation_iv_1v_probe.json")
out = Path("build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/simulation_iv_1v_no_minority_relax_probe.json")
cfg = json.loads(base.read_text())
cfg["output_csv"] = str(out.with_suffix(".csv").resolve())
cfg["solver"]["contact_boundary_minority_electron_relaxation"] = False
cfg["sweep"]["vtk_prefix"] = str(out.with_suffix("").resolve())
cfg["sweep"]["diagnostics"]["terminal_balance"]["csv"] = str(out.with_name("iv_no_minority_relax_terminal_balance.csv").resolve())
cfg["sweep"]["diagnostics"]["contact_edge"]["csv"] = str(out.with_name("iv_no_minority_relax_contact_edges.csv").resolve())
out.write_text(json.dumps(cfg, indent=2) + "\n")
print(out)
PY
build\vela_example_runner.exe --config build\reference_tcad\pn2d_sentaurus2018\reports\iv_state\simulation_iv_1v_no_minority_relax_probe.json
```

Expected: one of these outcomes is observed:

```text
Outcome A: sweep reaches beyond 0.828125 V and Anode phin validation no longer fails.
Outcome B: same failure persists, so relaxation is not the immediate cause.
Outcome C: Newton fails earlier, so relaxation was masking a separate numerical instability.
```

- [ ] **Step 3: Compare terminal current at the last common bias**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python - <<'PY'
import csv
from pathlib import Path
for p in [
    Path("build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/iv_1v_probe.csv"),
    Path("build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/simulation_iv_1v_no_minority_relax_probe.csv"),
]:
    if not p.exists():
        continue
    rows = [r for r in csv.DictReader(p.open()) if r["current_contact"] == "Cathode" and r["converged"] == "1"]
    last = rows[-1]
    print(p.name, last["bias_V"], last["current_total_A_per_um"])
PY
```

Expected: decide whether disabling relaxation improves convergence only, or also moves the IV current toward Sentaurus.

---

### Task 3: Add a Unit Test for Contact Boundary Construction

**Files:**
- Modify: `tests/test_newton_solver.cpp`
- Read: `src/solver/NewtonSolver.cpp`

- [ ] **Step 1: Add a failing test that asserts ohmic contact majority and minority QFs at high bias**

Add this test near the existing Newton contact boundary tests:

```cpp
TEST_CASE("NewtonSolver: high-bias ohmic contacts keep quasi-Fermi boundary targets",
          "[newton][contacts]")
{
    DeviceMesh mesh;
    mesh.nodes = {
        {0, 0.0, 0.0},
        {1, 1.0e-6, 0.0},
    };
    mesh.elements = {{0, 1, 1}};
    mesh.contacts = {
        Contact{"Anode", {0}},
        Contact{"Cathode", {1}},
    };

    DopingProfile doping;
    doping.setNodeDoping(0, -1.0e23);
    doping.setNodeDoping(1, 1.0e23);

    NewtonConfig cfg;
    cfg.contactBoundaryMinorityElectronRelaxation = false;
    NewtonSolver solver(mesh, doping, cfg);

    std::unordered_map<std::string, Real> biases = {
        {"Anode", 0.828125},
        {"Cathode", 0.0},
    };
    auto result = solver.solve(biases);

    REQUIRE(result.solution.phin(0) == Catch::Approx(0.828125).margin(1.0e-10));
    REQUIRE(result.solution.phip(0) == Catch::Approx(0.828125).margin(1.0e-10));
    REQUIRE(result.solution.phin(1) == Catch::Approx(0.0).margin(1.0e-10));
    REQUIRE(result.solution.phip(1) == Catch::Approx(0.0).margin(1.0e-10));
}
```

Expected before any fix: if current defaults or code path collapse Anode `phin`, this fails and pins the contact-boundary bug.

- [ ] **Step 2: Run the targeted test**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R newton_solver
```

Expected: the new test result tells whether the issue reproduces at unit level or requires `DCSweep` warm-start/adaptive stepping.

---

### Task 4: Add a Validation Test for Majority-Carrier Contact Checks

**Files:**
- Modify: `tests/test_solution_validation.cpp`
- Read: `src/solver/SolutionValidation.cpp`

- [ ] **Step 1: Add a test for p-contact validation behavior**

Add:

```cpp
TEST_CASE("SolutionValidation: p-contact validates hole quasi-Fermi against bias",
          "[validation][contacts]")
{
    DeviceMesh mesh;
    mesh.nodes = {{0, 0.0, 0.0}};
    mesh.contacts = {Contact{"Anode", {0}}};

    DDSolution sol;
    sol.psi = VectorXd::Constant(1, 0.0);
    sol.phin = VectorXd::Constant(1, 0.0);
    sol.phip = VectorXd::Constant(1, 0.828125);
    sol.n = VectorXd::Constant(1, 1.0e10);
    sol.p = VectorXd::Constant(1, 1.0e17);
    sol.converged = true;

    std::unordered_map<std::string, Real> biases = {{"Anode", 0.828125}};
    DDSolutionValidationOptions options;
    options.checkContactQuasiFermiBias = true;

    auto result = validateDDSolution(sol, mesh, biases, options);
    REQUIRE(result.valid);
}
```

Expected: validation should accept a p-contact whose majority hole quasi-Fermi equals the bias, even if minority electron quasi-Fermi differs.

- [ ] **Step 2: Add the complementary ambiguous-contact test**

Add:

```cpp
TEST_CASE("SolutionValidation: ambiguous contact requires both quasi-Fermi fields",
          "[validation][contacts]")
{
    DeviceMesh mesh;
    mesh.nodes = {{0, 0.0, 0.0}};
    mesh.contacts = {Contact{"Anode", {0}}};

    DDSolution sol;
    sol.psi = VectorXd::Constant(1, 0.0);
    sol.phin = VectorXd::Constant(1, 0.0);
    sol.phip = VectorXd::Constant(1, 0.828125);
    sol.n = VectorXd::Constant(1, 1.0e12);
    sol.p = VectorXd::Constant(1, 1.0e12);
    sol.converged = true;

    std::unordered_map<std::string, Real> biases = {{"Anode", 0.828125}};
    DDSolutionValidationOptions options;
    options.checkContactQuasiFermiBias = true;

    auto result = validateDDSolution(sol, mesh, biases, options);
    REQUIRE_FALSE(result.valid);
}
```

Expected: the validation logic is documented by tests before changing it.

- [ ] **Step 3: Run validation tests**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R solution_validation
```

Expected: tests clarify whether the reported `phin` failure is a true boundary-condition bug or a validator choosing the wrong carrier due to density inversion.

---

### Task 5: Trace Density Inversion at the Failed Anode Node

**Files:**
- Read: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/iv_1v_probe_0018_0.826563V.vtk`
- Read: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/pn2d_iv_state_node_comparison.csv`
- Create: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/physical_quantity_compare/iv_anode_node_density_trace.csv`

- [ ] **Step 1: Extract node 2 and neighboring Anode nodes from the VTK and comparison CSV**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python - <<'PY'
import csv
from pathlib import Path
src = Path("build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/pn2d_iv_state_node_comparison.csv")
out = Path("build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/physical_quantity_compare/iv_anode_node_density_trace.csv")
keep = {"2"}
rows = []
for r in csv.DictReader(src.open()):
    if r["node_id"] in keep:
        rows.append(r)
with out.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
print(out)
print(rows[0])
PY
```

Expected: node 2 shows whether `n > p` at the Anode, which would make validation check `phin` instead of `phip`.

- [ ] **Step 2: Interpret validation failure against contact carrier majority**

Record:

```text
If node 2 has n > p, the validation failure is caused by contact density inversion.
If node 2 has p > n but validation checked phin, validation has a majority-detection bug.
If node 2 has phin=0 and p > n, the minority relaxation is probably acceptable but Sentaurus comparison still needs a policy decision.
```

Expected: this separates physical high-injection inversion from an implementation defect.

---

### Task 6: Compare Sentaurus Current-Density and Recombination Drivers

**Files:**
- Read: `build/reference_tcad/pn2d_sentaurus2018/sim_fields/iv/fields/eCurrentDensity_region0.csv`
- Read: `build/reference_tcad/pn2d_sentaurus2018/sim_fields/iv/fields/hCurrentDensity_region0.csv`
- Read: `build/reference_tcad/pn2d_sentaurus2018/sim_fields/iv/fields/TotalCurrentDensity_region0.csv`
- Read: `build/reference_tcad/pn2d_sentaurus2018/sim_fields/iv/fields/srhRecombination_region0.csv`
- Create: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/physical_quantity_compare/iv_sentaurus_current_recombination_summary.csv`

- [ ] **Step 1: Generate Sentaurus-only current/recombination stats by region group**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python - <<'PY'
import csv, statistics
from pathlib import Path
fields = Path("build/reference_tcad/pn2d_sentaurus2018/sim_fields/iv/fields")
out = Path("build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/physical_quantity_compare/iv_sentaurus_current_recombination_summary.csv")
rows = []
for name in ["eCurrentDensity", "hCurrentDensity", "TotalCurrentDensity", "srhRecombination"]:
    vals = [float(r["component0"]) for r in csv.DictReader((fields / f"{name}_region0.csv").open())]
    rows.append({
        "quantity": name,
        "min": min(vals),
        "mean": statistics.mean(vals),
        "median": statistics.median(vals),
        "max": max(vals),
    })
out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
print(out)
for r in rows:
    print(r)
PY
```

Expected: establish whether Sentaurus high current is mainly electron current, hole current, or recombination-driven.

- [ ] **Step 2: Map the dominant Sentaurus driver to Vela equations**

Use this mapping:

```text
eCurrentDensity mismatch -> inspect electron SG flux and electron mobility.
hCurrentDensity mismatch -> inspect hole SG flux and hole mobility.
srhRecombination mismatch -> inspect SRH lifetime/intrinsic-density/BGN inputs.
TotalCurrentDensity mismatch with balanced components -> inspect unit scaling and contact-current integration.
```

Expected: only one next physics subsystem is selected after boundary behavior is resolved.

---

### Task 7: Implement the Minimal Fix Supported by Evidence

**Files:**
- Candidate modify: `src/solver/NewtonSolver.cpp`
- Candidate modify: `src/solver/SolutionValidation.cpp`
- Candidate modify: `include/vela/solver/NewtonSolver.h`
- Candidate modify: `tests/test_newton_solver.cpp`
- Candidate modify: `tests/test_solution_validation.cpp`

- [ ] **Step 1: If Task 2 proves minority relaxation causes the collapse, disable it for ohmic Sentaurus-calibration IV**

Minimal config-only fix for the reference probe:

```json
{
  "solver": {
    "contact_boundary_minority_electron_relaxation": false
  }
}
```

Expected: Vela reaches beyond `0.828125 V` without the Anode `phin=0` validation failure.

- [ ] **Step 2: If Task 5 proves density inversion makes validation choose the wrong carrier, change validation only**

The intended behavior should be encoded as:

```cpp
// For ohmic contact validation, prefer configured/contact doping majority
// over solved n/p majority when high injection inverts the boundary node.
```

Expected: validation checks the physically prescribed majority for the contact, not a high-injection swapped carrier field.

- [ ] **Step 3: If both QFs should remain pinned for ohmic contacts, change boundary construction**

The invariant should be:

```cpp
bcs.phin[nid] = Vbias / potentialScale;
bcs.phip[nid] = Vbias / potentialScale;
```

Expected: both contact quasi-Fermi fields equal applied bias at ohmic contacts unless a specific non-ohmic contact model opts out.

- [ ] **Step 4: Re-run targeted tests**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R "newton_solver|solution_validation|dc_sweep"
```

Expected: targeted tests pass and no contact-boundary regression appears.

---

### Task 8: Re-run IV Physical Quantity Comparison After the Fix

**Files:**
- Read: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/simulation_iv_1v_probe.json`
- Create: `build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/fixed/*`

- [ ] **Step 1: Re-run the IV probe to 1.0 V**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
build\vela_example_runner.exe --config build\reference_tcad\pn2d_sentaurus2018\reports\iv_state\simulation_iv_1v_probe.json
```

Expected: the run either converges closer to `1.0 V`, or fails with a new non-boundary diagnostic.

- [ ] **Step 2: Recompute IV curve interpolation**

Run the existing comparison command:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts/compare_reference_curves.py --reference build/reference_tcad/pn2d_sentaurus2018/reference_curves/pn2d_sentaurus2018_iv_reference.csv --candidate build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/iv_1v_probe.csv --kind iv --output build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/fixed/iv_curve_compare.json
```

Expected: Vela/Sentaurus ratio at high bias improves from the current `0.228` at `0.8265625 V`.

- [ ] **Step 3: Recompute field comparison at the highest common bias**

Generate a new node-field summary and heatmaps using the same script pattern used for:

```text
build/reference_tcad/pn2d_sentaurus2018/reports/iv_state/physical_quantity_compare/
```

Expected:

```text
electron_qf_V max diff no longer contains Anode node 2 = 1.0 V vs 0.0 V.
eDensity/hDensity mean abs diff decreases from ~1.1e17 cm^-3.
ElectricField p95 diff decreases from ~4169 V/cm.
```

---

### Task 9: Full Verification

**Files:**
- Read/verify only.

- [ ] **Step 1: Run focused regression tests**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
ctest --test-dir build --output-on-failure -R "sg_flux|newton_solver|solution_validation|dc_sweep|reference_tcad_regression|sentaurus_import_tools"
```

Expected: all focused tests pass.

- [ ] **Step 2: Run the full suite if focused tests pass**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
ctest --test-dir build --output-on-failure
```

Expected: full suite passes.

- [ ] **Step 3: Update the validation report**

Update:

```text
docs/validation/pn2d_sentaurus_comparison.md
```

Include:

```markdown
### IV high-bias physical quantities

- Boundary anomaly fixed:
- Highest Vela converged bias:
- IV high-bias current ratio:
- Potential/ElectricField/eDensity/hDensity improvements:
- Remaining mismatch:
```

Expected: the report distinguishes fixed boundary behavior from remaining physics-model mismatch.

---

## Stop Conditions

- Stop after Task 2 if disabling minority relaxation does not change the failure; the next hypothesis should be DCSweep warm-start or validation majority selection.
- Stop after Task 5 if high-injection density inversion is confirmed; decide whether Vela should validate against contact doping majority or solved carrier majority before editing code.
- Stop after any fix if the high-bias IV current gets worse while boundary validation improves; that means the boundary fix exposed a separate mobility/recombination scaling mismatch.
