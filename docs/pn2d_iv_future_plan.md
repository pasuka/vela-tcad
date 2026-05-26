# pn2d IV/BV — Future Development Plan

Status snapshot (2026-05-26). The Sentaurus pn2d IV slope gap and BV SRH
jump are still open. This document records the verified state of the
investigation and the next concrete steps so work can resume without
re-deriving context.

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

- **F. Residual-based contact current.** Implement an alternate
  `ContactCurrent` path that sums the *un-replaced* electron + hole
  continuity residuals at Dirichlet contact nodes (the standard
  Sentaurus extraction). Compare with the SG-flux-integral result.
  If the residual-based sum collapses to ~Newton residual but is still
  off vs reference, the gap is physics/geometry, not post-processing.
  Touch points: `src/post/ContactCurrent.cpp`,
  `src/equation/CoupledDDAssembler.cpp` (need to expose the
  pre-replacement residual or recompute it on demand for contact nodes).

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
`include/vela/physics/RecombinationModel.h` does not match Sentaurus
Scharfetter doping-dependent lifetimes. Implement the Scharfetter
parameterization (`tau(N) = tau_max / (1 + (N / N_ref)^gamma)`) with
the standard `(tau_max, N_ref, gamma)` defaults from Sentaurus and
plumb it through `MaterialDB` so it can be overridden per-material in
the simulation JSON.

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
