# PN2D incremental variants

Eight unique silicon DD configurations are declared in [experiments.json](experiments.json).
Actual formulas, parameter values and contact conventions are documented in [MODEL_AUDIT.md](MODEL_AUDIT.md).
D0 and G0 are aliases of P0. The original case and existing templates are unchanged.
The prospective final gates are in [acceptance.json](acceptance.json); the initial,
more permissive current floor is retained in [acceptance_prebaseline.json](acceptance_prebaseline.json).
The floor was tightened using P0 equilibrium noise before inspecting new variant results.

Generated reviewable M0 SDE, SDevice, Vela and material inputs live under `inputs/<ID>`.
They are generated from the same manifest as the M1 refinement, not independently
maintained copies. `models.par` is copied unchanged from `../source` into each run.
Do not execute the review snapshots in the source tree: their common imported mesh
and node doping, all simulation outputs and both complete sweeps are materialized
in `build-release/pn2d_variants/<ID>/<mesh>`.

## Reproduce from this worktree

All commands run at the current worktree root in PowerShell:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --preset windows-ucrt64-release
cmake --build --preset windows-ucrt64-release --target vela_example_runner sentaurus_import test_sg_flux test_dd_gummel -j 4
python scripts/pn2d_variants.py generate --export-review-inputs

# Establish the same-version baseline and refinement controls first.
python scripts/pn2d_variants.py sentaurus --ids P0 --meshes original,M0,M1,M2,M3,E1 --workers 2
python scripts/pn2d_variants.py vela --ids P0 --meshes original,M0,M1,M2,M3,E1 --compact-spatial

# P1/P2/P3 copy each corresponding P0 TDR, including the original nodal doping.
python scripts/pn2d_variants.py sentaurus --ids P1,P2,P3 --meshes M0,M1,M2,M3 --workers 2
python scripts/pn2d_variants.py vela --ids P1,P2,P3 --meshes M0,M1,M2,M3 --compact-spatial
python scripts/compare_pn2d_variants.py

# These commands reject execution until all stage A gates pass.
python scripts/pn2d_variants.py sentaurus --ids D1,D2 --meshes M2,M3 --workers 4
python scripts/pn2d_variants.py vela --ids D1,D2 --meshes M2,M3 --compact-spatial
python scripts/compare_pn2d_variants.py
# C additionally requires completed, fingerprint-checked stage B qualification.
# E1 adds independent endpoint refinement; P0/E1 was solved above.
python scripts/pn2d_variants.py sentaurus --ids G1,G2 --meshes M2,M3,E1 --workers 4
python scripts/pn2d_variants.py vela --ids G1,G2 --meshes M2,M3,E1 --compact-spatial
python scripts/compare_pn2d_variants.py
python scripts/report_pn2d_variants.py --export-evidence

./build-release/test_sg_flux.exe
ctest --preset windows-ucrt64-release -R '^(dd|pn2d_variants|pn2d_config_templates|sentaurus_import_tools|reference_tcad_regression|pn2d_bv_m2_mixed_voronoi_self_consistent_control|pn2d_node_volume_policy_default_acceptance)$'
```

The default SSH alias is `sentaurus`; its binary banner must report
`T-2022.03-SP2`. The channel used for this campaign actually logs in as `root`,
despite the older environment document's `tcad` user entry. Job paths are isolated
under `/home/tcad/sentaurus_runs/vela_oracle_2022/pn2d_variants_<fingerprint>`.
Use `--license-server 27000@127.0.0.1` only for the same VM's local license daemon;
this changes the child process environment, not the VM configuration.
The current VM has four CPUs and 8 GB RAM; four single-thread native jobs fit
its measured memory budget. Lower `--workers` on smaller VMs. Native queues
start larger meshes first; `--workers` does not parallelize local Vela solves.

Run records contain process status, input fingerprints, tool version and the
actual Vela backend. A successful run is not acceptance: `results.json` contains
input, numerical, terminal, spatial, conservative-section, mesh and model-effect
gates. `report/report.md` and PNGs visualize the actual available results.
Missing evidence stays `not_run`/`incomplete`; failed numerical scans remain failures.

Independent modeling before stage A qualification is explicitly available via
`sentaurus-mesh --ids D1,D2,G1,G2 --meshes M0,M1 --workers 2`, followed by
`import-mesh --ids D1,D2,G1,G2 --meshes M0,M1`. The latter exports and audits exact
node doping and contact extents without launching a DD solve or bypassing its gate.

`--compact-spatial` solves every requested bias and retains each accepted Vela
restart. It suppresses all-point VTK writes, then reconstructs the five requested
full fields from those same independent states. `spatial_artifacts.json` records
state and output fingerprints; diagnostic frozen-state postprocessing is never
used as independent convergence evidence. Comparison with the original P0/M0
VTKs agreed to 2.12e-16 relative roundoff. Native Sentaurus extra Plot endpoints
remain on the VM; transfer archives include the five required full states.

M0/M1 did not meet the frozen P0 mesh thresholds. The failure remains in results.
[mesh_refinement.json](mesh_refinement.json) records continued refinement to M2
and M3, retaining the failed M1/M2 result (2.2456% P0 low-bias hole-current change).
Each halves junction x limits again while keeping M1 transverse/endpoint limits.
The final pair is M2/M3 for all eight configurations; the P0 original mesh remains
the same-version baseline control. Coarse B/C DD controls are optional, and their
mesh-only audits do not count as solved cases. M3 uses a full initial step for each
explicit SDevice Goal and saves only the required spatial states; adaptive retry,
Digits=8, iteration limits and all exact comparison biases remain unchanged.
This continuation is mesh qualification, not a relaxation of error tolerances.
E1 retains all M3 global/junction limits and halves only the four contact endpoint
window limits. P0/G1/G2 must additionally satisfy the same current/integral mesh
gates on M3/E1. This control was declared before any G1/G2 DD result was inspected.
Endpoint regions are circular neighborhoods of the actual Anode endpoints,
derived from imported contact nodes; unused segmentation vertices are excluded.
Regional RMS, p95, maxima and source integrals are recorded across refinements.

B/C native runs use `ExtendedPrecision(80)`, unchanged `Digits=8` and explicit
`RHSMin=1e-20`. Initial
D1/D2 double-precision runs failed small-current gates; an isolated EP128 probe
removed the native equilibrium residual and restored minority-current agreement.
Original results remain in each `initial_native_double` archive. The manifest's
`sentaurus_numerics_by_stage` keeps this repair separate from the physical factors.
Each B/C Goal starts with its full segment and permits adaptive retries; all exact
comparison voltages and required spatial states remain. A inputs, Vela inputs,
mesh scripts, material laws and thresholds are unchanged. An isolated D1/M3
EP80 probe also passes the zero/-0.1 V current gates (equilibrium 3.7e-23 A/um).
Strict-RHS EP80 controls in both D1/M3 and D2/M3 pass zero/0.02 V, so B/C use
this lower-cost precision. This choice precedes G results. Partial EP128 strict
cost controls are retained separately. Short probes do not qualify complete
sweeps; all original gates still apply.
Reproduce the isolated
checks with `python scripts/probe_pn2d_sentaurus_precision.py` and
`python scripts/diagnose_pn2d_equilibrium_roundoff.py`; add `--mesh M3 --precision 80`
to the native probe command for the EP80 control. These never qualify a sweep.
The first EP128 full sweep exposed a separate stop-condition issue at 0.02 V:
native update error was still 811 when default `RHSMin=1e-5` stopped Newton.
Strict-RHS probes (`--branch forward --rhs-min 1e-20`, with either precision)
reduce the port-component errors below 0.057%. Complete native sweeps are then
rerun; comparison points are not patched from the short probes. The preceding
complete/partial attempts remain under `initial_EP128_loose_rhs`.

## Numerical and physical conventions

- P0: low-field Masetti, fixed-lifetime SRH, OldSlotboom, Boltzmann, 300 K;
  no high-field saturation, avalanche, Auger or Fermi statistics.
- P1: constant 1417/470.5 cm2/(V s) mobility via the explicit parameter file.
- P2: SRH disabled in both decks. P3: BGN disabled, including its associated
  `dEg0(OldSlotboom)=-0.01595 eV` intrinsic-density convention.
- SRH lifetimes are 1e-5/3e-6 s. `SRH(DopingDep)` is not selected.
- Every branch independently starts from equilibrium; every comparison voltage
  is a solved SDevice Goal and an explicit Vela `bias_points` target.
- Vela uses the unchanged IV template as its starting point, with a campaign-local
  `continuity_row_scaling.flux_fraction=1`. P0 source-only scaling stalled at
  millivolt biases; flux-aware left scaling resolves that conditioning problem
  while retaining local and global conservation requirements.
- P3 exposed cancellation in the homogeneous-ni SG implementation. Its electron
  and hole helpers now reuse the existing expm1-factorized equal-ni limit, keeping
  sub-femtovolt conductive increments. This numerical fix changes no physical
  model/default or convergence threshold; the new small-signal conductance test
  fails before the fix and passes afterward.
- Anode current is positive under forward bias. Vela's conventional hole component
  is the negative of its raw `current_hole` column; total equals electron minus hole.
  AreaFactor and out-of-plane depth are one, giving A/um. No fitted current factors.
- Each binary SDevice mesh is imported directly. D1/D2 regenerate node doping;
  region-average doping overrides are never used. The junction-plane compensation
  remains exactly as exported (`reported` policy).
- M0 and M1 split the left boundary at y=0.125, 0.1875, 0.3125, 0.375 um;
  endpoint mesh windows are identical across contact variants. M1 halves length
  limits. Uncontacted segments remain natural insulating boundaries.
- Spatial node matching checks both unique IDs and coordinates. Carrier metrics
  use the reference-active population. Endpoint region statistics and area
  integrals accompany maxima; nodal reconstructed current vectors are diagnostic.
  The conservative authority is the exact production SG edge flux.
  `analyze_pn2d_sections.py` additionally integrates both tools' reconstructed
  nodal vectors across the same cuts, exported by the report command. These
  diagnostic line integrals are not native discrete face fluxes; their residuals
  are reported separately and are not labeled conservation passes.
- VTK densities are cm^-3. Restart CSV density columns are serialized in SI m^-3
  by the current writer and converted by its matching reader. Field comparison
  therefore uses VTK/TDR common units rather than reading restart columns as cm^-3.

Historical source limits (10 V), IV template defaults (20 V), old fixture acceptance
(0.2--0.3 V), and the old 1/2/5/10/15/20 V guard are distinct. This campaign does not
claim those high-voltage anchors have been requalified for Sentaurus 2022.

## Current-vector recovery diagnostic

The legacy `Sentaurus*CurrentDensityVector` output multiplies a recovered nodal
quasi-Fermi gradient by nodal conductivity. It is not the production SG flux.
For G1/G2 at 0.8 V, replaying the exported native state through the same recovery
reproduces the large nodal-vector discrepancy. Existing dual-face and cell-first
SG recoveries agree much more closely, without solving again or fitting currents.
The original metrics are retained. Reproduce on available saved states:

```powershell
python scripts/diagnose_pn2d_current_reconstruction.py --id G1 --mesh E1 --bias 0.8
python scripts/diagnose_pn2d_current_reconstruction.py --id G2 --mesh E1 --bias 0.8
```

These are frozen-state diagnostics, not independent-solve evidence or new gates.
Their JSON records, region statistics and shared-scale plots are exported with
`report_pn2d_variants.py --export-evidence`.
