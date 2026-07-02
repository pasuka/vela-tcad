# NMOS Id-Vd Sentaurus Reference Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a rectangular 2-D NMOS Sentaurus/Vela reference fixture for first-pass Id-Vd curve comparison.

**Architecture:** Add `reference_tcad/nmos2d_sentaurus2018` alongside the existing pn2d fixture. Keep Sentaurus source decks explicit and simple, import the TDR/PLT artifacts through the existing `scripts/sentaurus_import.py reference` path, and add one small import-tool hook so MOS contact metadata can mark the gate as `metal_gate` in the generated Vela deck.

**Tech Stack:** C++20 Vela runner, Python stdlib regression tests, Sentaurus SDE/SDevice command decks, JSON reference config.

---

### Task 1: Lock Reference Semantics With Tests

**Files:**
- Modify: `tests/regression/test_reference_tcad_tools.py`

- [ ] **Step 1: Add a failing inventory/config test**

Add tests that require `reference_tcad/nmos2d_sentaurus2018/nmos2d_sentaurus2018_reference.json`, validate `device == "nmos2d"`, require an `iv` simulation whose Sentaurus sweep is `Drain` current versus `Drain` voltage, and require the generated Vela current contact to be `Drain`.

- [ ] **Step 2: Add a failing source-deck content test**

Require the SDE deck to contain `R.Si`, `R.Ox`, `Source`, `Drain`, `Gate`, `Body`, p-body doping, and n+ source/drain doping. Require the SDevice IV deck to include pn2d-aligned physics: `Mobility(DopingDependence)`, `Recombination(SRH)`, and `EffectiveIntrinsicDensity(OldSlotboom)`.

- [ ] **Step 3: Add a failing import override test**

Import `scripts/sentaurus_import.py`, write a temporary generated Vela deck with contacts, call `patch_reference_deck()` with `vela_contact_overrides`, and assert that `Gate` becomes `type: "metal_gate"` with `flatband_voltage: 0.0` while `Drain` remains ohmic.

### Task 2: Add NMOS Source Fixture

**Files:**
- Create: `reference_tcad/nmos2d_sentaurus2018/source/nmos2d_sde.cmd`
- Create: `reference_tcad/nmos2d_sentaurus2018/source/nmos2d_idvd_sdevice.cmd`
- Create: `reference_tcad/nmos2d_sentaurus2018/nmos2d_sentaurus2018_reference.json`

- [ ] **Step 1: Write the SDE deck**

Use a rectangular 2.0 um by 0.4 um silicon body, a 1.3 um centered oxide gate from x=0.35 to 1.65 um, p-body acceptor doping `1e17 cm^-3`, and n+ source/drain donor doping `1e17 cm^-3`. Put contacts on Source, Drain, Gate, and Body.

- [ ] **Step 2: Write the SDevice Id-Vd deck**

Initialize all electrodes at 0 V, solve Poisson then coupled DD, hold Gate at 2 V, and sweep Drain from 0 V to 0.5 V with an initial step of `1e-3` and max step `0.02`. Export potential, carrier densities, quasi-Fermi potentials, current, doping, recombination, and mobility diagnostics.

- [ ] **Step 3: Write reference config**

Configure the imported simulation as `kind: "iv"`, `bias_column: "Drain OuterVoltage"`, `current_column: "Drain TotalCurrent"`, `vela_current_contact: "Drain"`, `vela_stop: 0.5`, `vela_step: 0.02`, and `comparison.candidate_column: "current_total_A_per_um"`.

### Task 3: Patch MOS Contact Overrides

**Files:**
- Modify: `scripts/sentaurus_import.py`

- [ ] **Step 1: Implement minimal override application**

Inside `patch_reference_deck()`, after contacts are reset to zero bias, read `sim.get("vela_contact_overrides", [])`. For each override, find the generated contact by name and update only the provided JSON fields.

- [ ] **Step 2: Keep behavior unchanged for pn2d**

If the override list is absent, no generated deck changes. Existing pn2d reference config and tests must remain valid.

### Task 4: Verify

**Files:**
- Test: `tests/regression/test_nmos2d_reference_fixture.py`

- [ ] **Step 1: Run targeted regression tests**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python -m unittest tests.regression.test_nmos2d_reference_fixture
```

Expected final state: all selected tests pass.

- [ ] **Step 2: Report runnable comparison command**

Provide the command:

```powershell
python scripts\sentaurus_import.py reference --config reference_tcad\nmos2d_sentaurus2018\nmos2d_sentaurus2018_reference.json --source-dir reference_tcad\nmos2d_sentaurus2018\source --output-dir build\reference_tcad\nmos2d_sentaurus2018 --tdr-importer build\sentaurus_import.exe --runner build\vela_example_runner.exe
```

Note that the command needs Sentaurus-generated `.tdr` and `.plt` artifacts in the fixture source directory before a live comparison can complete.
