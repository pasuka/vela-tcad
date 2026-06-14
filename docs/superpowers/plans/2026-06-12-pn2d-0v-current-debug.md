# PN2D 0V Current Debug Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Localize why Vela's PN2D 0 V terminal current magnitude is about `9.13e6` times larger than Sentaurus while terminal signs, contact quasi-Fermi boundary values, and two-terminal conservation are correct.

**Architecture:** Treat the current mismatch as a measurement-chain and physics-chain problem. First prove whether the gap is caused by unit/current-definition conversion, then reconstruct terminal current from Sentaurus current-density fields, then compare the Vela edge-current ingredients against the Sentaurus state on the same contact nodes before attempting any solver fix.

**Tech Stack:** C++20, CMake/Ninja, Catch2, Python 3, Pillow for PNG summaries, Vela `ContactCurrent`, `DCSweep`, `NewtonSolver`, Sentaurus 2018 PN2D fixture exports.

---

## Current Evidence

- Vela terminal signs match Sentaurus:

```text
Anode:   Vela -6.5533928359887347e-18 A/um, Sentaurus -7.17389811693691e-25
Cathode: Vela  6.5556001087542772e-18 A/um, Sentaurus  7.17389811693687e-25
```

- Sentaurus/Vela absolute terminal-current ratio is about `1.094e-7`, so Vela is about `9.13e6` larger.
- Vela two-terminal current balance passes:

```text
electron_minus_hole_sum_A_per_um = 2.2072727655425378e-21
pair_balance_relative = 3.3670033695236683e-4
```

- Contact quasi-Fermi boundary values are not the source:

```text
Anode max eQF diff = 1.83689e-16 V
Anode max hQF diff = 4.59224e-17 V
Cathode max eQF diff = 2.29612e-16 V
Cathode max hQF diff = 0.0 V
```

- Node-state mismatch still exists in the body:

```text
Potential mean abs diff = 8.238815e-03 V
Electron QF mean abs diff = 7.461135e-05 V, max = 4.392980e-03 V
Hole QF mean abs diff = 7.461156e-05 V, max = 4.392980e-03 V
Electron density mean abs diff = 1.987844e14 cm^-3
Hole density mean abs diff = 1.989463e14 cm^-3
Electric field mean abs diff = 4.340789e02 V/cm
```

- Current-related comparison files already exist under:

```text
D:\code-repo\vela-tcad\build\reference_tcad\pn2d_sentaurus2018\reports\0v_current_related
```

## File Structure

- Modify `scripts/compare_pn2d_0v_current_related_quantities.py`: extend the report with stronger unit checks, boundary-integrated Sentaurus currents, and same-node current-driver diagnostics.
- Modify `tests/regression/test_reference_tcad_tools.py`: add smoke coverage for any new script arguments and output schema.
- Modify `docs/validation/pn2d_sentaurus_comparison.md`: record each debug stage and the selected root-cause hypothesis.
- Inspect only unless the evidence points here:
  - `src/post/ContactCurrent.cpp`
  - `include/vela/post/ContactCurrent.h`
  - `src/simulation/DCSweep.cpp`
  - `src/solver/NewtonSolver.cpp`
  - `src/equation/CoupledDDAssembler.cpp`

## Task 1: Lock Down Current Units And Definitions

**Files:**
- Modify: `scripts/compare_pn2d_0v_current_related_quantities.py`
- Modify: `tests/regression/test_reference_tcad_tools.py`

- [x] **Step 1: Add a schema smoke test for the current-related report**

Add this test to `tests/regression/test_reference_tcad_tools.py`:

```python
def test_compare_pn2d_0v_current_related_quantities_help() -> None:
    result = subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts" / "compare_pn2d_0v_current_related_quantities.py"),
            "--help",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    self.assertEqual(result.returncode, 0)
    self.assertIn("--reference-root", result.stdout)
    self.assertIn("--current-balance", result.stdout)
    self.assertIn("--edge-csv", result.stdout)
```

- [x] **Step 2: Run the smoke test**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python -m unittest tests.regression.test_reference_tcad_tools.ReferenceTcadToolsTest.test_compare_pn2d_0v_current_related_quantities_help
```

Expected: pass.

- [x] **Step 3: Add explicit current conversion candidates**

Extend the report with these candidates for each contact:

```python
conversion_candidates = {
    "sentaurus_raw_vs_vela_A_per_um": sentaurus_total / vela_total_a_per_um,
    "sentaurus_A_per_um_vs_vela_A_per_um": (sentaurus_total / 1.0) / vela_total_a_per_um,
    "sentaurus_A_per_cm_width_to_A_per_um": (sentaurus_total / 1.0e4) / vela_total_a_per_um,
    "sentaurus_A_per_m_width_to_A_per_um": (sentaurus_total / 1.0e6) / vela_total_a_per_um,
    "vela_A_per_m_vs_sentaurus": vela_total_a_per_m / sentaurus_total,
}
```

Keep these as diagnostic fields only; do not select a fix from this step.

- [x] **Step 4: Re-run the report**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\compare_pn2d_0v_current_related_quantities.py
```

Expected: report is regenerated, and the candidate ratios show whether a simple width conversion explains `9.13e6`.

## Task 2: Reconstruct Sentaurus Terminal Current From CurrentDensity

**Files:**
- Modify: `scripts/compare_pn2d_0v_current_related_quantities.py`

- [x] **Step 1: Build contact boundary edge sets from mesh triangles**

For each contact node set, find triangle edges where both nodes are contact nodes and the edge belongs to exactly one silicon triangle.

```python
def contact_boundary_edges(elements: list[list[int]], contact_nodes: set[int]) -> list[tuple[int, int]]:
    counts: dict[tuple[int, int], int] = {}
    for a, b, c in elements:
        for u, v in ((a, b), (b, c), (c, a)):
            edge = tuple(sorted((u, v)))
            counts[edge] = counts.get(edge, 0) + 1
    return [
        edge for edge, count in counts.items()
        if count == 1 and edge[0] in contact_nodes and edge[1] in contact_nodes
    ]
```

- [x] **Step 2: Integrate Sentaurus normal current-density magnitude along each contact**

Use the nodal average of `TotalCurrentDensity`, `eCurrentDensity`, and `hCurrentDensity`, edge length in cm, and unit depth assumptions:

```python
edge_current_A_per_cm_width = 0.5 * (j0 + j1) * edge_length_cm
```

Emit both signed and absolute sums for each contact. If Sentaurus scalar `CurrentDensity` lacks vector direction, mark the signed result as `direction_unresolved`.

- [x] **Step 3: Compare three Sentaurus currents**

Report for each contact:

```text
.plt TotalCurrent
TDR ContactCurrentFlux
CurrentDensity boundary integral, A/cm-width
```

Expected debug decision:
- If boundary integral matches `.plt` after one width factor, the mismatch is mostly Vela unit/current extraction.
- If boundary integral matches Vela scale, the mismatch is mostly Sentaurus `.plt`/TDR contact definition.
- If all three differ, continue to Task 3.

## Task 3: Compare Same-Node Current Drivers At Contacts

**Files:**
- Modify: `scripts/compare_pn2d_0v_current_related_quantities.py`
- Inspect: `src/post/ContactCurrent.cpp`

- [x] **Step 1: Add contact-node driver CSV**

Create `contact_current_driver_nodes.csv` with one row per Anode/Cathode node and these columns:

```text
contact,node_id,x_um,y_um,
sentaurus_potential_V,vela_potential_V,potential_diff_V,
sentaurus_eQF_V,vela_eQF_V,eQF_diff_V,
sentaurus_hQF_V,vela_hQF_V,hQF_diff_V,
sentaurus_eDensity_cm3,vela_eDensity_cm3,eDensity_ratio_vela_to_sentaurus,
sentaurus_hDensity_cm3,vela_hDensity_cm3,hDensity_ratio_vela_to_sentaurus,
sentaurus_eCurrentDensity_A_cm2,
sentaurus_hCurrentDensity_A_cm2,
sentaurus_totalCurrentDensity_A_cm2
```

- [x] **Step 2: Add contact-adjacent interior-node driver CSV**

For every Vela contact edge from `contact_edges.csv`, include both edge nodes and the adjacent interior node if it can be inferred from the triangle containing that edge. Emit:

```text
contact,edge_id,contact_node0,contact_node1,interior_node,
psi_contact_avg,psi_interior,
phin_contact_avg,phin_interior,
phip_contact_avg,phip_interior,
n_contact_avg,n_interior,
p_contact_avg,p_interior,
vela_current_electron_A_per_m,
vela_current_hole_A_per_m,
vela_current_total_A_per_m
```

- [x] **Step 3: Decide whether the contact gradient is the source**

Expected debug decision:
- If contact-node densities and QF agree but interior-node QF differs by millivolts, continue to Task 4.
- If contact-node densities differ strongly, inspect contact equilibrium density initialization in `src/solver/NewtonSolver.cpp`.
- If Vela edge currents are large despite tiny QF gradients, inspect `src/post/ContactCurrent.cpp` current formula and scaling.

## Task 4: Isolate BGN/Effective-Ni Impact On Current, Not Just QF Span

**Files:**
- Modify: `scripts/probe_pn2d_0v_qf_drivers.py`
- Modify: `scripts/compare_pn2d_0v_current_related_quantities.py`

- [x] **Step 1: Run the existing QF matrix again**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\probe_pn2d_0v_qf_drivers.py --reference-root build\reference_tcad\pn2d_sentaurus2018 --runner build\vela_example_runner.exe --output-dir build\reference_tcad\pn2d_sentaurus2018\reports\0v_qf_drivers
```

Expected baseline:

```text
baseline QF span ~= 0.004393061 V
no_bgn QF span ~= 8.39006e-13 V
```

- [x] **Step 2: Run current-related comparison on baseline and no_bgn VTKs**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
python scripts\compare_pn2d_0v_current_related_quantities.py --vtk build\reference_tcad\pn2d_sentaurus2018\reports\0v_qf_drivers\baseline\baseline_0000_0V.vtk --terminal-csv build\reference_tcad\pn2d_sentaurus2018\reports\0v_qf_drivers\baseline\baseline_terminal_balance.csv --edge-csv build\reference_tcad\pn2d_sentaurus2018\reports\0v_qf_drivers\baseline\baseline_contact_edges.csv --output-dir build\reference_tcad\pn2d_sentaurus2018\reports\0v_current_related\baseline
python scripts\compare_pn2d_0v_current_related_quantities.py --vtk build\reference_tcad\pn2d_sentaurus2018\reports\0v_qf_drivers\no_bgn\no_bgn_0000_0V.vtk --terminal-csv build\reference_tcad\pn2d_sentaurus2018\reports\0v_qf_drivers\no_bgn\no_bgn_terminal_balance.csv --edge-csv build\reference_tcad\pn2d_sentaurus2018\reports\0v_qf_drivers\no_bgn\no_bgn_contact_edges.csv --output-dir build\reference_tcad\pn2d_sentaurus2018\reports\0v_current_related\no_bgn
```

Expected debug decision:
- If `no_bgn` also collapses terminal current near Sentaurus scale, root cause is BGN/effective-ni consistency.
- If `no_bgn` fixes QF span but not current scale, root cause is current post-processing or unit convention.

## Task 5: Inspect And Test The Selected Root Cause

**Files:**
- Modify only one target after Tasks 1-4 select it:
  - `src/post/ContactCurrent.cpp`
  - `src/solver/NewtonSolver.cpp`
  - `src/equation/CoupledDDAssembler.cpp`
- Modify tests:
  - `tests/test_dc_sweep.cpp`
  - `tests/regression/test_reference_tcad_tools.py`

- [ ] **Step 1: Write the failing test for the selected behavior**

If current scaling is selected, add a focused unit test in `tests/test_dc_sweep.cpp` that asserts:

```text
current_total_A_per_um == current_total / 1.0e6
sum(contact_edges.current_total) == terminal_balance.current_total within 1e-12 relative
```

If BGN/effective-ni is selected, add a focused 0V equilibrium test that asserts:

```text
max(abs(ElectronQuasiFermi - contact_qf)) < 1.0e-6 V
max(abs(HoleQuasiFermi - contact_qf)) < 1.0e-6 V
```

- [ ] **Step 2: Run the failing focused test**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
ctest --test-dir build --output-on-failure -R "(dc_sweep|reference_tcad)"
```

Expected: the new assertion fails before the fix.

- [ ] **Step 3: Apply one minimal fix**

Only edit the selected file from this matrix:

```text
ContactCurrent.cpp      current formula, edge-length, width, A/m to A/um scaling
NewtonSolver.cpp        0V contact/equilibrium QF and density initialization
CoupledDDAssembler.cpp  continuity residual BGN/effective-ni consistency
```

- [ ] **Step 4: Verify focused and regression tests**

Run:

```powershell
$env:Path='D:\msys64\ucrt64\bin;D:\msys64\usr\bin;' + $env:Path
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R "(dc_sweep|newton|gummel_high|line_search|sentaurus_tdr|reference_tcad|diagnose_pn2d)"
```

Expected: focused test passes, existing PN2D diagnostics still pass, and the generated current-related report shows the terminal-current ratio has moved toward Sentaurus without breaking sign or conservation.

## Debug Stop Criteria

Stop and report instead of fixing if any of these are true:

- Sentaurus `.plt`, `ContactCurrentFlux`, and current-density boundary integral disagree by more than three orders of magnitude after all documented width/unit conversions.
- Vela edge-current sum disagrees with Vela terminal balance by more than `1e-6` relative.
- `no_bgn` removes QF split but does not materially change current magnitude.
- Three attempted minimal fixes fail to improve the same selected metric.

## Execution Result

- [x] Tasks 1-4 executed.
- [x] Task 5 stopped before solver edits because the first stop criterion is true: Sentaurus `.plt TotalCurrent`, TDR `ContactCurrentFlux`, and boundary `TotalCurrentDensity` integral differ by more than three orders after documented conversions.
- [x] Root-cause hypothesis selected for the next implementation phase: BGN/effective-ni consistency. Evidence: `no_bgn` collapses QF span from `0.0043930610509000005 V` to `8.39006e-13 V` and moves the Vela A/m current scale to roughly `1.2-1.5x` Sentaurus `.plt`.
- [x] No solver fix was applied in this debug pass.

## Follow-up Execution Result (2026-06-13)

- [x] Added failing SG and coupled-assembler tests for variable-`ni` flat-QF equilibrium.
- [x] Implemented variable-intrinsic-density quasi-Fermi SG fluxes.
- [x] Updated `CoupledDDAssembler` residual/Jacobian and `ContactCurrent` to use the variable-`ni` QF flux path instead of falling back to density SG on BGN edges.
- [x] Refreshed 0V diagnostics: current-balance status is now `pass`, classification is `balanced`, baseline QF span is `1.827779e-08 V`, and baseline Anode/Cathode terminal currents are near numerical zero.
