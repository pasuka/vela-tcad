# PN2D BV Next Stage Jacobian Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Determine whether the current PN2D BV mismatch is caused by missing Jacobian physics or by model/discretization parity, then improve the BV knee comparison without destabilizing the SG transport core.

**Architecture:** Keep the core SG flux-divergence discretization unchanged. Treat the current `density_gradient` SG edge-current avalanche path as the production BV path, add focused Jacobian diagnostics around real BV states, and only promote physics changes that improve the full IV/BV curve shape against Sentaurus.

**Tech Stack:** C++20, Catch2, CMake/Ninja, MSYS2 UCRT64, Python diagnostics, existing `build-release/reference_tcad/pn2d_sentaurus2018` artifacts.

---

### Task 1: Record Current Jacobian Completeness Baseline

**Files:**
- Read: `src/equation/CoupledDDAssembler.cpp`
- Read: `include/vela/equation/AssemblerUtils.h`
- Read: `src/physics/RecombinationModel.cpp`
- Read: `tests/test_impact_ionization.cpp`
- Read: `tests/test_newton_solver.cpp`
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`

- [ ] **Step 1: Run targeted Jacobian tests**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
.\build-release\test_impact_ionization.exe "[impact]~[gummel]"
.\build-release\test_newton_solver.exe "[newton][coupled]"
```

Expected:

```text
All tests passed
All tests passed
```

- [ ] **Step 2: Document the code-level baseline**

Append this exact status note to `docs/validation/pn2d_sentaurus_comparison.md`:

```markdown
### PN2D BV Jacobian Audit Baseline

Current BV production physics uses `van_overstraeten`, `driving_force:
quasi_fermi_gradient`, `generation: current_density`, and
`current_approximation: density_gradient`. In this SG edge-current avalanche
path, `CoupledDDAssembler::assembleJacobian` finite-differences the combined
edge avalanche source with respect to the six endpoint potentials, so the
matrix includes carrier-density, alpha driving-field, and local field-dependent
edge-mobility derivatives for that source discretization.

The non-SG node-local avalanche path remains intentionally approximate: it
includes local carrier-density derivatives but omits driving-field and mobility
derivatives. That path is not the current PN2D BV production path and must not
be used as evidence that the SG BV Jacobian is incomplete.
```

- [ ] **Step 3: Re-run the targeted tests after documentation**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
.\build-release\test_impact_ionization.exe "[impact]~[gummel]"
.\build-release\test_newton_solver.exe "[newton][coupled]"
```

Expected: both commands pass.

- [ ] **Step 4: Commit**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git add docs/validation/pn2d_sentaurus_comparison.md
git commit -m "Document PN2D BV Jacobian audit baseline"
```

### Task 2: Add Real-State Jacobian Block Audit

**Files:**
- Create: `scripts/diagnose_pn2d_bv_jacobian_block_audit.py`
- Modify: `tests/regression/test_reference_tcad_tools.py`

- [ ] **Step 1: Add a regression test for report parsing**

Append to `tests/regression/test_reference_tcad_tools.py`:

```python
def test_pn2d_bv_jacobian_block_audit_parses_fixture(tmp_path):
    csv = tmp_path / "audit.csv"
    csv.write_text(
        "bias_V,block,analytic_norm,fd_norm,diff_norm,rel_diff\n"
        "-13.2,sg_avalanche,2.0,2.0,1.0e-8,5.0e-9\n",
        encoding="utf-8",
    )
    rows = list(csv.read_text(encoding="utf-8").splitlines())
    assert rows[0].split(",") == [
        "bias_V",
        "block",
        "analytic_norm",
        "fd_norm",
        "diff_norm",
        "rel_diff",
    ]
    assert "-13.2,sg_avalanche" in rows[1]
```

- [ ] **Step 2: Run the failing test**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python -m pytest tests/regression/test_reference_tcad_tools.py -k jacobian_block_audit -q
```

Expected: pass, because this first test locks the CSV contract before the diagnostic exists.

- [ ] **Step 3: Create the diagnostic script**

Create `scripts/diagnose_pn2d_bv_jacobian_block_audit.py` with a command-line interface:

```python
#!/usr/bin/env python
import argparse
import csv
from pathlib import Path


FIELDS = [
    "bias_V",
    "block",
    "analytic_norm",
    "fd_norm",
    "diff_norm",
    "rel_diff",
]


def write_placeholder_report(output: Path, bias_values: list[float]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        for bias in bias_values:
            for block in ("srh_auger", "sg_avalanche", "transport", "poisson"):
                writer.writerow({
                    "bias_V": f"{bias:.6g}",
                    "block": block,
                    "analytic_norm": "nan",
                    "fd_norm": "nan",
                    "diff_norm": "nan",
                    "rel_diff": "nan",
                })


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--bias", action="append", type=float, required=True)
    args = parser.parse_args()
    write_placeholder_report(args.output, args.bias)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
```

- [ ] **Step 4: Replace placeholders with real C++/runner-backed extraction**

Use existing sweep states or restart configs from `build-release/reference_tcad/pn2d_sentaurus2018/vela/`. For each selected bias (`-11.5`, `-13.2`, `-18.0`, `-19.0`, `-20.0`), compute:

```text
rel_diff = ||J_analytic_block - J_fd_block|| / max(1, ||J_fd_block||)
```

Report blocks:

```text
poisson
transport
srh_auger
sg_avalanche
dirichlet_or_gauge
```

Acceptance:

```text
sg_avalanche rel_diff <= 5e-5 on small extracted submesh or replay fixture
srh_auger rel_diff <= 5e-5
transport rel_diff <= 1e-4
```

- [ ] **Step 5: Run the diagnostic**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts/diagnose_pn2d_bv_jacobian_block_audit.py `
  --bias -11.5 --bias -13.2 --bias -18 --bias -19 --bias -20 `
  --output build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_jacobian_block_audit/jacobian_blocks.csv
```

Expected: CSV exists and contains one row per block per bias.

### Task 3: Test Low-Density Avalanche Driving-Force Interpolation Derivatives

**Files:**
- Modify: `tests/test_impact_ionization.cpp`
- Modify: `src/equation/CoupledDDAssembler.cpp` only if the new test proves a missing derivative on the SG production path.

- [ ] **Step 1: Add a failing test for RefDens interpolation on SG path**

Add a Catch2 test next to `Coupled DD SG edge-current avalanche Jacobian captures field-dependent alpha`:

```cpp
TEST_CASE("Coupled DD SG avalanche Jacobian captures low-density driving-force interpolation",
          "[impact][newton]")
{
    DeviceMesh mesh = makePNMesh();
    MaterialDatabase matdb;
    const std::vector<RegionDopingSpec> specs = {
        {"n_region", 5.0e22, 0.0},
        {"p_region", 0.0, 5.0e22},
    };
    DopingModel doping = DopingModel::fromMeshAndRegions(mesh, specs);

    const Real Vt = 0.025852;
    const int N = static_cast<int>(mesh.numNodes());
    CoupledDDState state;
    state.psi = VectorXd::LinSpaced(N, -0.08, 0.08);
    state.phin = VectorXd::LinSpaced(N, 0.3, -0.3);
    state.phip = VectorXd::LinSpaced(N, -0.3, 0.3);

    ImpactIonizationModelConfig impactConfig;
    impactConfig.model = "van_overstraeten";
    impactConfig.drivingForce = "quasi_fermi_gradient";
    impactConfig.generation = "current_density";
    impactConfig.currentApproximation = "density_gradient";
    impactConfig.drivingForceInterpolation = "quasi_fermi_to_electric_field";
    impactConfig.electronDrivingForceRefDensity = 1.0e20;
    impactConfig.holeDrivingForceRefDensity = 1.0e20;

    CoupledDDAssembler assembler(
        mesh,
        matdb,
        doping,
        Vt,
        mobilityModelConfig("constant"),
        recombinationModelConfig({"none"}),
        BandgapNarrowingConfig{},
        impactConfig);

    const VectorXd x = assembler.pack(state);
    const CoupledDDBoundaryConditions bcs;
    const SparseMatrixd analytic = assembler.assembleJacobian(x, bcs);
    const SparseMatrixd finiteDifference = assembler.finiteDifferenceJacobian(x, bcs, 1.0e-7);
    const Eigen::MatrixXd diff = Eigen::MatrixXd(analytic - finiteDifference);
    const Eigen::MatrixXd ref = Eigen::MatrixXd(finiteDifference);
    REQUIRE(diff.norm() / std::max<Real>(1.0, ref.norm()) < 5.0e-5);
}
```

- [ ] **Step 2: Run only the new test**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
.\build-release\test_impact_ionization.exe "low-density driving-force interpolation"
```

Expected before any fix: if this fails, inspect the derivative of interpolation weight with respect to endpoint densities; if it passes, document interpolation is not the missing Jacobian term.

- [ ] **Step 3: Fix only if the test fails**

If the test fails, update `edgeAvalancheCombinedSource` in `src/equation/CoupledDDAssembler.cpp` so the finite-differenced source recomputes the endpoint densities and interpolated driving fields for each perturbed endpoint. Do not add a separate hand-derived interpolation Jacobian unless the finite-difference block is removed.

- [ ] **Step 4: Run impact tests**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
.\build-release\test_impact_ionization.exe "[impact]~[gummel]"
```

Expected: pass.

### Task 4: Quantify BV Knee Shape Against Sentaurus

**Files:**
- Create: `scripts/diagnose_pn2d_bv_knee_shape.py`
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`

- [ ] **Step 1: Create the knee-shape script**

Create a script that reads:

```text
build-release/reference_tcad/pn2d_sentaurus2018/vela/pn2d_sentaurus2018_bv_minus20_avaljac.csv
build-release/reference_tcad/pn2d_sentaurus2018/reference_curves/pn2d_sentaurus2018_bv_reference.csv
```

It should output:

```text
first bias where one-volt current growth ratio exceeds 1.5
first bias where one-volt current growth ratio exceeds 2.0
maximum absolute log10 current error over -10V to -20V
```

- [ ] **Step 2: Run the script**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts/diagnose_pn2d_bv_knee_shape.py
```

Expected current baseline:

```text
Sentaurus knee: approximately -18V to -19V
Vela current avaljac curve: earlier step-like turn-up near -11V to -12V
```

- [ ] **Step 3: Document acceptance gates**

Append:

```markdown
BV parity is not accepted solely because the -20V sweep converges. The next
acceptance gate requires the Vela knee to move toward the Sentaurus knee
window, approximately -18V to -19V, and requires the -10V to -20V curve to
avoid artificial plateaus or early step transitions near -11V to -12V.
```

### Task 5: Decide the Next Physics Change

**Files:**
- Read: `scripts/diagnose_pn2d_bv_loop_gain_sensitivity.py`
- Read: `scripts/diagnose_pn2d_bv_active_edge_flux_factors.py`
- Read: `scripts/diagnose_pn2d_bv_source_policy_matrix.py`
- Modify: `reference_tcad/pn2d_sentaurus2018/pn2d_sentaurus2018_reference.json` only after the diagnostic evidence supports a production configuration change.

- [ ] **Step 1: Compare three hypotheses**

Run diagnostics for:

```text
H1: low-density quasi-Fermi-to-electric-field interpolation changes the early knee
H2: avalanche source ownership/support changes the early knee
H3: ionization coefficient/model parameters change the early knee without harming 0V and low reverse bias
```

- [ ] **Step 2: Reject unsafe changes**

Reject any change that:

```text
rewrites core SG flux divergence
uses source_geometry_scale as a hidden calibration factor
improves only one bias point while worsening the knee location
requires disabling SRH when matching the Sentaurus BV physics block
```

- [ ] **Step 3: Promote a candidate only after full sweep**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
.\build-release\vela_example_runner.exe --config build-release\reference_tcad\pn2d_sentaurus2018\vela\simulation_bv_minus20_avaljac.json
```

Expected promotion criteria:

```text
converged=true
nonconverged points = 0
minimum bias <= -20V
Vela knee no earlier than about -16V unless a documented Sentaurus-physics reason explains the difference
```

