# Sentaurus T-2022.03-SP2 Templates/LDMOS reference

This directory contains only neutral, reviewable metadata and small normalized
reference artifacts for the Synopsys Applications Library
`Templates/LDMOS` case. Proprietary source decks, TDR/PLT files, logs, and
archives must remain under the ignored top-level `reference_staging/` tree.

The executable plan is defined by:

- `docs/superpowers/plans/2026-08-26-templates-ldmos-sentaurus-vela-validation-plan.md`
- `docs/superpowers/plans/2026-08-26-templates-ldmos-phase-a-oracle-classical-validation-plan.md`

WP0, stages 0/1, stage 1.5 and WP1.75 are complete.  Stage 1.5 qualifies
state serialization, frozen replay and same-bias closure; WP1.75 defines the
strict material, solver-physics and phase-A discretization contracts.  Neither
result is a classic DD curve acceptance or authorization to treat pending
hRecVelocity, hQP, IALMob, thermal or Okuto physics as implemented.

WP2 and phase-2/3 execution status is recorded in
`docs/validation/templates_ldmos_phase23_execution_2026-08-27.md`. Generated
TDR/PLT and exact-mesh products remain under ignored `reference_staging`.

Expected persistent local layout:

```text
reference_staging/templates_ldmos_sentaurus2022/<run-id>/
  manifest/
  sentaurus_original/
  bundle/
  raw/
  normalized/
  imported_structure/
  comparisons/
  reports/
```

Before any selected neutral artifact is committed, its source manifest,
Sentaurus release, extraction method, units, current sign, and SHA-256 must be
reviewed.

The tracked entry points are:

The approved neutral execution budget is stored in `budget_freeze.json`;
run-local copies under `reference_staging/` remain the evidence-bound source
used by the summary generator.

```powershell
# WP0/stage 0: isolated, read-only-source VM run
python scripts/run_templates_ldmos_sentaurus_vm.py --live `
  --sentaurus-version T-2022.03-SP2 --run-id <unique-id> `
  --stages sprocess,idvg,idvd,bv

# Stage-0 fail-closed curve/log/state audit
python scripts/audit_templates_ldmos_oracle.py `
  --run-dir reference_staging/templates_ldmos_sentaurus2022/<unique-id> `
  --ssh-target sentaurus

# Output-only derivative state decks; these are never the official oracle
python scripts/prepare_templates_ldmos_state_decks.py `
  --bundle-dir reference_staging/templates_ldmos_sentaurus2022/<unique-id>/bundle `
  --idvg-curve <normalized-idvg-drain-curve.csv> `
  --bv-table <normalized-bv-full-table.csv> --output-dir <derived-dir>
python scripts/run_templates_ldmos_state_capture.py `
  --run-dir reference_staging/templates_ldmos_sentaurus2022/<unique-id> `
  --state-decks <derived-dir>

# A corrected or partial BV-only capture can be kept separate, then merged
# without overwriting either sealed source capture.
python scripts/run_templates_ldmos_state_capture.py `
  --run-dir reference_staging/templates_ldmos_sentaurus2022/<unique-id> `
  --state-decks <derived-dir> --stages bv --output-name representative_states_bv_v2
python scripts/consolidate_templates_ldmos_states.py `
  --base-states <representative-states> --bv-states <corrected-bv-states> `
  --state-decks <derived-dir> --output <representative-states-final>
python scripts/classify_templates_ldmos_states.py `
  --states-dir <representative-states-final>/raw `
  --importer build-release/sentaurus_import.exe `
  --output <representative-states-final>/state_inventory.json

# Stage 1: reported and dominant-signed doping imports, deterministic repeat,
# exact topology conversion, and structural audit
python scripts/analyze_templates_ldmos_structure.py --tdr <n1_fps.tdr> `
  --source-coordinate-unit cm `
  --output-dir reference_staging/templates_ldmos_sentaurus2022/<unique-id>/stage1

# Structural timing lower bound and unsigned review budget
python scripts/run_templates_ldmos_cost_probe.py `
  --stage1-dir reference_staging/templates_ldmos_sentaurus2022/<unique-id>/stage1 `
  --vela-runner build-release/vela_example_runner.exe `
  --vela-commit <commit> `
  --oracle-manifest reference_staging/templates_ldmos_sentaurus2022/<unique-id>/manifest/run_manifest.json

# Stage 1.5: generate exact-mesh qualification decks after exporting the
# three declared Sentaurus state CSV files, run the decks in dependency order,
# then score frozen replay and settled same-bias reclose.
python scripts/prepare_templates_ldmos_restart_qualification.py `
  --stage1-dir reference_staging/templates_ldmos_sentaurus2022/<unique-id>/stage1
python scripts/summarize_templates_ldmos_restart_qualification.py `
  --qualification-dir reference_staging/templates_ldmos_sentaurus2022/<unique-id>/stage1/qualification

# WP1.75: validate/materialize versioned contracts next to the run evidence
python scripts/prepare_templates_ldmos_wp175_contracts.py `
  --stage1-dir reference_staging/templates_ldmos_sentaurus2022/<unique-id>/stage1

# Assemble or refresh the review summary after all available prerequisite gates.
python scripts/finalize_templates_ldmos_phase01.py `
  --run-dir reference_staging/templates_ldmos_sentaurus2022/<unique-id> `
  --stage1-dir reference_staging/templates_ldmos_sentaurus2022/<unique-id>/stage1
```

`run_templates_ldmos_cost_probe.py` may then produce the draft budget from the
stage-1 exact mesh. Its Poisson-only placeholder materials are an explicit
runtime lower bound, not an accepted material or physics contract.
