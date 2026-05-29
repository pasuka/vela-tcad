# pn2d IV/BV — Future Development Plan

Status snapshot (2026-05-26). The Sentaurus pn2d IV slope gap and BV SRH
jump are still open. This document records the verified state of the
investigation and the next concrete steps so work can resume without
re-deriving context.

## Status update (2026-05-29, commit c6bb006)

Phases 0–2B are closed and committed:

- **Phase 2A (IV magnitude) — done.** The Sentaurus-faithful
  `caughey_thomas_field` mobility model was promoted into the IV deck via a
  per-simulation `vela_solver.mobility` override in
  `reference_tcad/pn2d/pn2d_reference.json`. The 0.2–0.3 V window orders drop
  from `0.5858` (default, fails the 0.50 gate) to `0.4214`, trend match holds,
  and the BV gate, 0.3 V terminal-current sum, and strict Newton handoff are
  unchanged. Scan: `scripts/scan_pn2d_iv_mobility_candidates.py`.
- **Phase 2B (IV local slope) — investigated, no promotion.** Re-running the
  contact-relaxation scan on the promoted base deck shows the n-contact-only QF
  relaxation meets the local-slope target (`Δ = 0.0702 ≤ 0.075`) but erases the
  Caughey-Thomas magnitude gain (IV orders rise to `0.4946`, only `0.0054`
  below the gate). The magnitude and local-slope axes trade against each other
  with these two levers, so the relaxation is **not** promoted; the residual
  local slope delta `~0.099` is carried as a documented known gap. Details in
  `docs/validation/pn2d_sentaurus_comparison.md`.
- **Phase 3 (cleanup) — done.** `ctest --preset windows-ucrt64-debug` is
  274/274 (no regression rebaseline needed); the validation doc is updated.

The remaining open item is **#2 the BV SRH/avalanche parity** below, which is a
separate physics-calibration track independent of the now-closed IV magnitude
axis.

## Phase 0 evidence refresh (2026-05-29, HEAD d05fc34)

Rebuilt (36/36 targets) and reran the IV probe; the older frozen
`|I_cath/ref| = 0.218 @ 0.30 V` figure is **stale**. On the current default
probe decks (which do not enable the minority-electron relaxation knob),
`scripts/probe_pn2d_iv_post_fix.py` now reports:

| bias V | \|I_cath/ref\| | two-contact sum (A/µm) |
| --- | --- | --- |
| 0.100 | 3.088 | −1.07e−17 |
| 0.200 | 0.912 | −3.58e−16 |
| 0.300 | 0.489 | −1.82e−14 |

The terminal-sum closure to ~4.3e−20 A/µm comes from the promoted iter2 deck
that enables the relaxation knob; these probe decks intentionally do not, so
the two figures are consistent.

Phase 1 discriminator (already on disk,
`build/pn2d_root_cause_probe/reports/taskI_qf_profile_compare_bias_aligned_1p0V_20260527.json`):
`psi` mean |diff| 0.0087 V, `phip` mean 0.0112 V (both agree), but `phin`
mean 0.0312 V with **max 0.9315 V Anode-localized**. So psi/phip agree and the
residual disagreement is the Anode electron-QF boundary node (handled by the
relaxation knob); the remaining magnitude gap is the mobility axis.

Mobility matrix `i0p29_ratio` (Vela/ref @ 0.29 V):
`default_field = 0.825`, `constant = 1.646`, `caughey_thomas = 0.997`,
`ct_promoted_bv = 0.971`. `caughey_thomas` nearly closes the 0.29 V point, so
Phase 2A targets the IV mobility model (the IV deck currently sets no explicit
mobility model while the BV deck already uses tuned `caughey_thomas`).

## Restart update after HEAD rerun

Fresh reruns on `3e84042` changed two conclusions that were true only for
older artifacts:

- The default Slotboom IV cathode current is already algebraically consistent
  with the coupled assembler residual. A direct residual probe at 0.30 V gives
  `assembler_residual = -1.897805e-14 A/um`, matching the CSV cathode total.
  The residual-based-contact-current hypothesis is therefore closed for the
  compared cathode terminal; the remaining IV error is in the solved
  transport/physics state or geometry/model calibration, not in cathode current
  extraction.
- The BGN-off/uniform-`ni` cases now take the QF branch introduced in the
  latest commit. For IV BGN-off at 0.30 V the density extraction would give
  `2.359720e-14 A/um`, but the QF/assembler residual gives
  `4.718129e-16 A/um`. The old "BGN-off essentially unchanged" note is stale.
- Sentaurus logs for both IV and BV report `no Lifetime file`,
  `no ModelParameters file`, and SRH with no field-, doping-, or
  temperature-dependent lifetimes. The BV SRH mismatch should not be described
  as missing Scharfetter doping-dependent lifetimes unless a later reference
  source proves hidden defaults. Current evidence points to constant-lifetime
  or SRH parameter parity.
- A constant equal-lifetime SRH scan cannot explain both IV and BV: IV is
  closest near `taun=taup=1e-8..3e-8 s`, while BV is closest near
  `taun=taup=3e-6 s`. Target ratios are recorded in
  `build/pn2d_tdr_tie_probe/vela/pn2d_taugrid_summary.csv`.
- A small asymmetric lifetime scan also does not find a joint explanation:
  the tested pairs improve neither metric together, with BV target ratios still
  `3.45x` or worse. Results are in
  `build/pn2d_tdr_tie_probe/vela/pn2d_taupair_summary.csv`.

## What is committed in this change

- `include/vela/post/ContactCurrent.h`, `src/post/ContactCurrent.cpp`:
  `ContactCurrent` now mirrors `CoupledDDAssembler::residual`'s
  Scharfetter–Gummel branch. When `ni_i == ni_j` on an edge it uses
  `sgElectronContinuityFluxFromQuasiFermi` / `sgHoleContinuityFluxFromQuasiFermi`
  (cancellation-free); otherwise it falls back to the density form.
  A `BandgapNarrowingConfig` parameter (default `{}`) was added to the
  constructor and to the static `compute(...)` overload so that the
  per-node `ni_` field is built with `buildValidatedEffectiveNodeNi`,
  exactly the same builder the assembler uses.
- `src/simulation/DCSweep.cpp`: forwards the previously-built
  `sweepBgnConfig` to `ContactCurrent`.
- `scripts/probe_pn2d_iv_post_fix.py`: minimal verification harness that
  re-runs the default (Cathode) and Anode-contact IV decks and prints
  per-bias `I_cath, I_anode, sum, |I_cath/ref|`.
- `scripts/probe_pn2d_iv_contact_decomposition.py`: read-only probe used
  during root-cause analysis.

These edits are correctness improvements (post-processor is now
algebraically consistent with the assembler for the uniform-`ni` edge
case) but, by design, do **not** change numerics on the default
`bandgap_narrowing: slotboom` deck, because `ni_i != ni_j` on every
edge there and both code paths use the density form.

## Verified findings

- Build: `cmake --build --preset windows-ucrt64-debug` — 35/35 targets.
- Tests: `ctest --preset windows-ucrt64-debug` — 272/273 pass. The only
  failure is `regression` (#270) because golden currents move with the
  consistency fix on decks that DO take the QF branch. Rebaseline only
  after the IV root cause is closed.
- Probe `scripts/probe_pn2d_iv_post_fix.py` on default Slotboom deck:

  | bias V | I_cath (A/µm) | I_anode (A/µm) |   sum    | \|I_cath/ref\| |
  |--------|---------------|----------------|----------|----------------|
  | 0.100  | −8.26e−17     | −7.69e−17      | −1.60e−16| 2.473          |
  | 0.200  | −8.50e−16     | −5.56e−16      | −1.41e−15| 0.566          |
  | 0.300  | −1.58e−14     | +4.16e−15      | −1.16e−14| 0.218          |

- BGN-off variant (where the new QF branch IS active in
  `ContactCurrent`) gives `I_cath ≈ −1.547e−14` A/µm at 0.30 V — same
  order as the Slotboom default. This **falsifies** the "post-processor
  catastrophic cancellation" hypothesis: the post-processor was already
  returning the correct discrete current; the slope shortfall and the
  two-contact non-conservation are properties of the discrete solution.

## Open issues (priority order)

### 1. pn2d IV slope shortfall + two-contact non-conservation (P0)

Symptoms: `|I_cath| ≈ 0.22 × reference` at 0.30 V; `I_cath + I_anode`
~1e−14 A/µm (same order as either contact); Newton converges in 2–3
iterations and tightening `reltol` from 1e-6 to 1e-10 does not move
the numbers; invariant under recombination on/off, BGN on/off, SRH
lifetime sweep, Auger sweep, and IV step refinement.

Next investigations (do in order, stop when conservation is restored):

- **F. Closed: residual-based cathode current check.** A Python recomputation
  of density, QF, and assembler-style residual currents shows the default
  Slotboom cathode CSV already equals the assembler-consistent residual at
  0.30 V. Do not implement a production residual-current path for this
  hypothesis unless a separate anode-oriented discrepancy is being targeted.

- **G. Geometry / unit-factor audit.** A bias-insensitive 0.2–0.3×
  shortfall is consistent with a per-µm vs per-cm or half-edge
  weighting bug. Compare at one bias (e.g. 0.30 V):
  - assembler edge weights and Voronoi areas used in the flux
    integral on contact-adjacent edges vs Sentaurus mesh dump,
  - the device thickness scaling (`scaling.thickness_um`) and how
    `ContactCurrent` maps the per-cell current to A/µm.

- **H. Mobility / SG-flux coefficient calibration.** Dump
  `mun, mup, n, p, psi` on contact-adjacent edges at 0.30 V from both
  Vela and Sentaurus and compare. Even with the same model family,
  Masetti / Caughey-Thomas constants can differ.

- **I. Quasi-Fermi profile check.** Plot Vela `phin, phip` vs Sentaurus
  `eQuasiFermi, hQuasiFermi` across the junction at 0.30 V. If `phin`
  and `phip` agree but currents don't, the gap is geometric or in the
  mobility model. If they disagree, the gap is in the assembler.

### 2. pn2d BV SRH jump (P1, independent of #1)

Constant `taun = taup = 1e-7 s` in
`include/vela/physics/RecombinationModel.h` does not match the observed BV
SRH sensitivity, but the Sentaurus logs explicitly say the run uses no
field-, doping-, or temperature-dependent lifetimes and no external lifetime
or model-parameter file. The next step is therefore to identify Sentaurus'
built-in constant SRH lifetime/parameter defaults for this model set, then
decide whether Vela needs configurable SRH trap parameters beyond `taun/taup`.

### 3. Cleanup once #1 is closed (P2)

- Rebaseline `tests/regression` golden currents (single failing test
  today: #270 in `ctest --preset windows-ucrt64-debug`). Do **not**
  rebaseline before #1 is closed.
- Update `docs/validation/pn2d_sentaurus_comparison.md` with the
  post-close IV ratio table and the final root-cause note.
- Consider fixing the diagnostic drift/diffusion split in
  `ContactCurrent` (the totals are correct, but the per-component
  midpoint formula still shows cancellation in the `electronDrift` /
  `electronDiffusion` columns).

## Reproducing the current state

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build --preset windows-ucrt64-debug
ctest --preset windows-ucrt64-debug
python scripts\probe_pn2d_iv_post_fix.py
```

Probe inputs (already generated in `build/pn2d_tdr_tie_probe/vela/`):
`simulation_iv_default.json`, `simulation_iv_anode_contact.json`,
`simulation_iv_iv_bgn_none.json`, plus `*_tighttol`,
`*_iv_recomb_none`, `*_iv_srh_tau{1e-6,1e-8}`, `*_iv_auger_only`,
`*_iv_srh_auger_{half,double}`, `*_resolution_step0p0{1,2}`,
`*_resolution_promoted`. Reference CSV:
`build/pn2d_tdr_tie_probe/reference_curves/pn2d_iv_reference.csv`.


## 2026-05-26 ROOT CAUSE: ContactCurrent totalCurrent sign-convention bug

### Finding

Tracing the SG flux sign conventions through `src/post/ContactCurrent.cpp`
and `src/discretization/ScharfetterGummel.cpp` reveals that the per-edge
accumulator stores:

- `result.electronCurrent = +I_n_outward`  (conventional electron current
  flowing out of the contact into the device)
- `result.holeCurrent     = -I_p_outward`  (NEGATIVE of conventional hole
  current flowing out of the contact)

The opposite sign on the hole accumulator arises because `holeFlux01 =
-pFlux` where `pFlux = sgHoleContinuityFlux` already represents J_p in the
0->1 direction (positive at anode in forward bias). The negation in
`holeFlux01 = -pFlux` therefore inverts the conventional hole sign.

Pre-fix line 168 was:

```cpp
result.totalCurrent = result.electronCurrent + result.holeCurrent;
//                  = I_n_outward - I_p_outward    <-- WRONG
```

This is not Kirchhoff-conserved. The correct total terminal current is

```cpp
result.totalCurrent = result.electronCurrent - result.holeCurrent;
//                  = I_n_outward + I_p_outward    <-- CORRECT
```

### Numerical evidence (pn2d default IV, V = 0.30 V)

| Contact | I_e (A/m)     | I_p (A/m)     | e+h (old)     | e-h (new)     |
|---------|---------------|---------------|---------------|---------------|
| Anode   | +2.135e-08    | -1.720e-08    | +4.158e-09    | +3.855e-08    |
| Cathode | -2.876e-08    | +9.787e-09    | -1.898e-08    | -3.855e-08    |
| SUM     |               |               | -1.482e-08    |  0.000e+00    |

Pre-fix sum-of-totals violates KCL by exactly 2*q*R_total = 1.48e-08 A/m
(volume recombination integral double-counted into terminal current).
Post-fix sum is identically zero at every bias point in `pn2d_iv_default.csv`
and `pn2d_iv_default_vtkprobe_anode_current.csv`.

### Impact on remaining gap to Sentaurus

At V=0.30 V Anode TotalCurrent:

- Pre-fix Vela / Sentaurus = 0.27x (~3.7x short)
- Post-fix Vela / Sentaurus = 0.55x (~1.8x short)

The fix closes roughly half of the gap. The residual ~1.8x is now
believed to be physics calibration: mobility model (Vela Caughey-Thomas
vs Sentaurus Masetti), BGN parameter parity, or the asymmetric mesh
(N-region has only 3 cathode nodes because the Sentaurus deck assigned
GlobalMesh to R.Si which became the P-half after R.NRegion overlay).

### Test fallout from the fix (must be addressed before merge)

The old formula coincidentally produced a positive number at MOS drains
under both NMOS and PMOS small-bias conditions (because
I_n_outward - I_p_outward yields the majority-carrier sign at the drain
in both cases). Several configs hard-coded `drain_current_sign = 1.0`
relying on that artifact:

- `examples/{nmos2d_dd,nmos2d_mos_dd,pmos2d_dd,pmos2d_mos_dd}/simulation_iv.json`
  - PMOS configs must be retuned: with the corrected convention
    (positive = current OUT of contact INTO device), a PMOS drain at
    Vds<0 reports a negative total. Update `drain_current_sign` to -1
    for the two PMOS decks, leave NMOS at +1.
- `tests/data/pn2d_iv_*` regression baselines may need refreshed
  reference magnitudes (now 2x larger by construction).
- `tests/regression/test_sentaurus_sample_integration.py` assertion
  `orders_of_magnitude < 0.20` is still violated (0.97 post-fix vs
  presumably worse pre-fix). The threshold is calibrated to the broken
  baseline; either tighten the physics gap further (see above) or relax
  the threshold to reflect known mobility-model differences.

Failing tests after fix (9 of 273):
- 122, 262: mos_solver_crosscheck (PMOS expectedSign now -1, not +1)
- 270: regression (PMOS drain polarity check)
- 273: sentaurus_sample_integration (orders_of_magnitude assertion)
- 249, 251, 252, 255, 272: SentaurusTdrReader.* (pre-existing,
  unrelated to ContactCurrent; observed failing on stashed baseline too)
