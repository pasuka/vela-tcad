# pn2d Recombination and IV Slope Follow-Up Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to execute this plan. Use superpowers:systematic-debugging for each numeric discrepancy before changing solver behavior. Track progress by checking off each step.

**Goal:** Turn the latest pn2d IV/BV discrepancy probes into targeted fixes or documented limits. Keep strict Newton ownership intact while separating IV sweep/slope error from BV recombination-physics error.

**Architecture:** Keep `reference_tcad/pn2d/pn2d_reference.json` stable until a candidate passes both numerical and physics-consistency checks. Add narrow scripts and tests first, then promote only the smallest config or solver change that is supported by the probe results.

**Tech Stack:** C++20, CMake/Ninja, Catch2, Python regression tests, PowerShell probe scripts, existing pn2d reference import and comparison tools.

---

## Current Evidence

- Latest reviewed commit: `0ce8e17 Add pn2d IV/BV physics matrix probe and findings`.
- The previous localization plan is partly complete:
  - `scripts/summarize_pn2d_iv_ratios.ps1` reports per-bias IV ratios.
  - `scripts/scan_pn2d_iv_bv_physics_matrix.ps1` reports IV/BV physics ablations.
  - `docs/validation/pn2d_sentaurus_comparison.md` records the new findings.
- Strict Newton handoff remains the desired invariant. Do not reintroduce Gummel fallback for reference gates.
- Current IV issue:
  - Vela/reference ratio is close around low forward bias, then rolls off at high forward bias.
  - At `0.29 V`, default IV ratio is about `0.3128`.
  - Turning recombination off or disabling BGN does not remove the slope error.
- Current BV issue:
  - The promoted BV gate is a strong low-bias numeric match with recombination disabled.
  - The reference BV command includes SRH, Auger, and avalanche physics.
  - Re-enabling SRH or SRH+Auger in Vela shifts BV by about `1.17` orders at the `0.05 V` check point.
- Config state:
  - `taun` and `taup` are already parsed by Gummel/Newton JSON config.
  - Auger coefficients exist in `RecombinationModelConfig`, but are not yet exposed through solver JSON config.

---

## Guardrails

- Do not weaken the strict Newton reference gate.
- Do not update `reference_tcad/pn2d/pn2d_reference.json` until a candidate passes the decision criteria in this plan.
- Do not mix IV and BV fixes unless one candidate improves both under the same physical assumptions.
- Keep generated probe output out of commits unless the task explicitly asks for a checked-in fixture.
- Keep documentation ASCII-only.

---

## Task 1: Isolate IV Sweep Resolution From True Slope Error

**Files:**
- Create: `scripts/scan_pn2d_iv_resolution.ps1`
- Modify: `tests/regression/test_reference_tcad_tools.py`
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`

- [x] Add a regression test that requires `scan_pn2d_iv_resolution.ps1` to exist and contain:
  - `pn2d_iv_resolution_summary.csv`
  - `0.02`
  - `0.01`
  - `current_total_A_per_um`
- [x] Implement the script to generate temporary pn2d IV configs from the imported reference tree with sweep steps:
  - existing promoted step
  - `0.02 V`
  - `0.01 V`
- [x] For each step, run strict `gummel_newton` and compare against the same reference IV window.
- [x] Report at least:
  - step size
  - accepted row count
  - orders of magnitude
  - max relative error
  - Vela/reference ratios at reference biases near `0.25 V`, `0.27 V`, and `0.29 V`
- [x] Document whether the high-bias rolloff is caused by coarse sweep/interpolation or by the solved current slope.

**Decision check:**
- If finer steps improve IV orders by at least `0.15`, consider promoting a smaller IV step only after runtime and strict-Newton acceptance are acceptable.
- If finer steps do not materially improve IV, keep the reference gate unchanged and continue to Tasks 2-4.

---

## Task 2: Sweep Existing SRH Parameters Before Adding New Knobs

**Files:**
- Modify: `scripts/scan_pn2d_iv_bv_physics_matrix.ps1`
- Modify: `tests/regression/test_reference_tcad_tools.py`
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`

- [x] Add script-content tests for named cases:
  - `bv_srh_tau1e-6`
  - `bv_srh_tau1e-8`
  - `iv_srh_tau1e-6`
  - `iv_srh_tau1e-8`
- [x] Extend the physics matrix to vary `taun` and `taup` while keeping all other knobs fixed.
- [x] Run BV cases with promoted BV mobility and recombination enabled.
- [x] Run IV cases with default IV mobility and recombination enabled.
- [x] Document whether a physically enabled SRH configuration can keep BV within `0.2` orders without worsening IV beyond the current gate.

**Decision check:**
- If SRH lifetime alone explains the BV jump, promote the lifetime only after adding a focused config test and updating the reference comparison doc.
- If no reasonable lifetime value reconciles BV, continue to Task 3.

---

## Task 3: Expose Auger Coefficients Through Solver JSON

**Files:**
- Modify: `include/vela/solver/GummelSolver.h`
- Modify: `include/vela/solver/NewtonSolver.h`
- Modify: `src/solver/GummelSolver.cpp`
- Modify: `src/solver/NewtonSolver.cpp`
- Modify: `docs/config_schema.md`
- Modify: `tests/test_mobility.cpp` or `tests/test_recombination.cpp`

- [x] Write a failing config-parsing test that verifies both Gummel and Newton accept:
  - `auger_cn_m6_per_s`
  - `auger_cp_m6_per_s`
- [x] Add the fields to Gummel/Newton config structs with defaults matching `RecombinationModelConfig`.
- [x] Pass the parsed coefficients into `RecombinationModelConfig` before constructing DD assemblers.
- [x] Validate that negative Auger coefficients still fail through the existing `RecombinationModel` checks.
- [x] Document the new keys under solver common controls.

**Decision check:**
- Keep this as a config-surface change unless a probe proves that default Auger coefficients are the primary BV or IV error.

---

## Task 4: Run Auger and Combined Recombination Parameter Matrix

**Files:**
- Modify: `scripts/scan_pn2d_iv_bv_physics_matrix.ps1`
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`

- [x] Add cases that scale Auger coefficients down and up while SRH is enabled.
- [x] Include at least one case with SRH disabled and Auger enabled, to separate Auger-only behavior from SRH behavior.
- [x] Compare IV and BV with the same summary columns used by the existing matrix.
- [x] Record the smallest parameter change that improves BV physics parity without damaging IV slope.

**Decision check:**
- If a combined recombination parameter set satisfies both IV and BV criteria, promote it as a reference candidate.
- If BV remains sensitive and IV remains shallow, document recombination as not sufficient and continue to Task 5.

---

## Task 5: Add a Narrow Recombination Diagnostic If Needed

**Files:**
- Prefer modifying an existing diagnostic path before adding a new public output contract.
- Possible files:
  - `src/equation/CoupledDDAssembler.cpp`
  - `include/vela/equation/CoupledDDAssembler.h`
  - `src/solver/NewtonSolver.cpp`
  - `tests/test_newton_solver.cpp`

- [x] Identify the smallest existing diagnostic hook that can report recombination magnitude or continuity residual contribution per terminal sweep point.
- [x] Add a test that the diagnostic is finite and opt-in.
- [x] Use the diagnostic to determine whether the BV SRH jump is caused by:
  - excessive carrier product near the junction
  - lifetime defaults
  - intrinsic-density or BGN interaction
  - boundary/contact carrier reconstruction
- [x] Keep diagnostic output disabled by default.

**Decision check:**
- If the diagnostic points to model mismatch, write a separate model-parity plan before implementing a broad physics change.
- If it points to a local bug, fix with a focused regression test.

---

## Task 6: Decide Whether Fermi Statistics Needs a Separate Implementation Plan

**Files:**
- Modify: `docs/validation/pn2d_sentaurus_comparison.md`
- Create only if justified: a new plan under `docs/superpowers/plans/`

- [x] Estimate degeneracy relevance from pn2d doping and carrier levels in the IV/BV windows.
- [x] Check whether the IV high-bias slope and BV SRH jump remain after Tasks 1-4.
- [x] If both remain unexplained, write a separate Fermi-statistics plan covering equations, tests, and reference impact. (Decision: not triggered in this stage; BV SRH jump is now localized to SRH sensitivity/model parity axis.)

**Decision check:**
- Do not start a broad Fermi implementation inside this plan unless the cheaper probes fail and the evidence specifically points there.

---

## Verification Commands

Run these before claiming the plan implementation is complete:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R "reference_tcad_tools|mobility|recombination|newton"
ctest --test-dir build --output-on-failure -R ascii_sources
```

For full closure, run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
ctest --test-dir build --output-on-failure
```
