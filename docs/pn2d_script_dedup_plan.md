# PN2D BV Contact-Tuning Script Dedup Plan

## Goal

Reduce maintenance overhead across the `scan_pn2d_bv_ct_*.ps1` scripts by
introducing a single parameterized base script and keeping compatibility through
thin wrappers.

This plan is intentionally non-functional: it should not change solver behavior,
only script organization and maintainability.

## Current Duplication Groups

Near-duplicate script families under `scripts/`:

- Grid/cartesian style:
  - `scan_pn2d_bv_ct_local.ps1`
  - `scan_pn2d_bv_ct_local_small.ps1`
- Explicit-case style:
  - `scan_pn2d_bv_ct_quick6.ps1`
  - `scan_pn2d_bv_ct_refine2.ps1`
  - `scan_pn2d_bv_ct_refine6.ps1`

Shared pipeline across these scripts:

1. Read base JSON config and reference CSV.
2. Build case list from mobility scaling dimensions.
3. Run `vela_example_runner.exe` for each case.
4. Parse convergence and sampled current rows.
5. Write ranked summary CSV.

Differences are mostly case generation strategy, timeout policy, and output
column granularity.

## Proposed Base Script

Create:

- `scripts/scan_pn2d_bv_ct_base.ps1`

Recommended interface:

- Input/output:
  - `-BaseConfig`
  - `-ReferenceCsv`
  - `-OutputSummary`
- Case definition:
  - `-CaseMode Grid|ExplicitList`
  - `-MuScales`
  - `-NrefScales`
  - `-AlphaScales`
  - `-ExplicitCases`
  - `-TagPrefix`
- Physics knobs:
  - `-BandgapNarrowing`
  - `-Recombination`
- Execution policy:
  - `-SecondsPerCase` (`0` means no timeout)
- Output policy:
  - `-ExtractComponentCurrents`
  - `-MaxResultsDisplay`

## Migration Stages

### Stage 0: Baseline

- Keep existing scripts unchanged.
- Capture baseline output CSV for each script variant.

### Stage 1: Introduce Base Script

- Add `scan_pn2d_bv_ct_base.ps1`.
- Validate output parity against at least:
  - `scan_pn2d_bv_ct_local.ps1`
  - `scan_pn2d_bv_ct_quick6.ps1`

### Stage 2: Migrate Grid Family

- Replace `local` and `local_small` script bodies with thin wrappers that call
  the base script with fixed parameter sets.
- Preserve existing script names and call signatures.

### Stage 3: Migrate Explicit Family

- Replace `quick6`, `refine2`, and `refine6` with thin wrappers.
- Keep old files as `.bak` during soak period if desired.

### Stage 4: Soak and Cleanup

- Run normal tuning workflow for one to two weeks.
- Remove archived `.bak` files once stable.

## Verification Checklist

For each migrated wrapper:

1. Exit code behavior matches previous script.
2. Summary CSV headers match expected schema.
3. Case count and tag naming match baseline.
4. Ranking and key ratios at probe bias (0.05 V) are unchanged or explainable.
5. Timeout behavior remains consistent (`SecondsPerCase` parity).

## Risk and Rollback

Primary risks:

- Timeout handling differences.
- CSV schema drift.
- Path handling assumptions in wrappers.

Rollback:

- Revert wrapper scripts to pre-migration bodies.
- Keep base script isolated so rollback does not affect solver code.

## Out of Scope

- Any drift-diffusion, Poisson, discretization, or model behavior change.
- Any regression baseline updates.
