# P1 Triangle Electric Field Recovery Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add post-processing P1 Tri3 electric-field and quasi-Fermi gradient recovery diagnostics without changing conservative FVM/SG assembly or avalanche source defaults.

**Architecture:** Put recovery algorithms in `vela/post/ElectricFieldDiagnostics` as public post-processing APIs. Let VTK diagnostics and a standalone Python comparison script consume these fields explicitly. Keep assembler-side SG/FVM code paths unchanged.

**Tech Stack:** C++20, Eigen `Point2`/`Point3`, Catch2, existing VTK legacy writer, Python stdlib CSV/Markdown diagnostics.

---

### Task 1: Recovery API Tests

**Files:**
- Modify: `tests/test_electric_field_diagnostics.cpp`

- [ ] Add tests for Tri3 cell constant field and quasi-Fermi cell gradients.
- [ ] Add tests proving area average, LS 1/d, LS 1/d2, and SPR recover a linear potential to machine precision at interior nodes.
- [ ] Add tests for boundary/corner fallback and region-wise recovery that does not cross a material interface.
- [ ] Run `ctest --test-dir build --output-on-failure -R electric_field_diagnostics` and confirm the new tests fail before implementation.

### Task 2: Recovery Implementation

**Files:**
- Modify: `include/vela/post/ElectricFieldDiagnostics.h`
- Modify: `src/post/ElectricFieldDiagnostics.cpp`

- [ ] Add data structs for cell vectors and node recovery result vectors.
- [ ] Implement Tri3 scalar gradient and `CellElectricField = -grad(Potential)`.
- [ ] Implement cell quasi-Fermi gradients.
- [ ] Implement node area average, LS 1/d, LS 1/d2, and SPR with per-node region-aware patch selection.
- [ ] Keep corner/contact fallback deterministic: SPR falls back to LS 1/d, singular LS falls back to area average, empty patches return zero.
- [ ] Keep existing `maxEdgeElectricFieldMagnitude` behavior intact.

### Task 3: VTK Post-Processing Output

**Files:**
- Modify: `include/vela/io/VTKWriter.h`
- Modify: `src/io/VTKWriter.cpp`
- Modify: `src/solver/GummelSolver.cpp`

- [ ] Add cell scalar/vector output support to the VTK writer.
- [ ] Write cell vectors `CellElectricField`, `CellGradElectronQuasiFermi`, and `CellGradHoleQuasiFermi`.
- [ ] Write node vectors and magnitudes for `NodeElectricField_AreaAverage`, `NodeElectricField_LS_1overD`, `NodeElectricField_LS_1overD2`, and `NodeElectricField_SPR`.
- [ ] Keep existing `ElectricField` / `ElectricFieldVector` names compatible by continuing to use the existing LS 1/d convention.

### Task 4: Sentaurus Recovery Comparison Script

**Files:**
- Create: `scripts/compare_electric_field_recovery.py`

- [ ] Load Sentaurus neutral export `nodes.csv`, `elements.csv`, `contacts.csv`, and `fields/ElectricField_region*.csv`.
- [ ] Compute the same four node recovery methods from Sentaurus or Vela nodal potential input.
- [ ] Write `build/diagnostics/electric_field_recovery_compare.csv`.
- [ ] Write `build/diagnostics/electric_field_recovery_compare_summary.md` with interior, boundary, and contact statistics and the required post-processing caveat sentence.

### Task 5: Verification

**Files:**
- No production edits expected.

- [ ] Configure/build with the MSYS2 UCRT64 preset.
- [ ] Run `ctest --test-dir build --output-on-failure -R electric_field_diagnostics`.
- [ ] Run a broader relevant test slice if time allows: `ctest --test-dir build --output-on-failure -R "electric_field|poisson|sentaurus_tdr_reader"`.
- [ ] Inspect `git diff` to confirm no FVM/SG assembly or default avalanche source path changed.
