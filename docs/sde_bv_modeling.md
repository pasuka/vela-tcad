# BVmethods source modeling (P5)

Corpus revision 4 models all eight unchanged primary SDE sources. BVmethods
now produces half and full `GRD`, `DAT`, and `MSH` files, with independent
native DF-ISE imports and analytic geometry/contact/doping checks. These
modeling gates do not qualify a BV sweep, NonLocalPath, or external-resistor R2/R3.

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts/sentaurus_gmsh_pipeline.py `
  --sde reference_tcad/bvmethods_sentaurus2018/source/bvmethods_nmos_sde.cmd `
  --parameter node=0 --output-dir "build-release/my BV mesh" `
  --runtime-dir build-release/gmsh-python-runtime
```

`device.*` is the full device; `half_device.*` and `half_mesh_field_bundle.json`
retain the authored half. The original `tdx -mtt -x ... drain=source` command
becomes a local x=0 reflection. Arbitrary shell commands remain unsupported.
All grids, TDR, and raw logs remain in ignored build directories.

## Geometry and lifecycle

`PolySi` explicitly maps to neutral `PolySilicon`; its original spelling stays
in Device IR. A surviving PolySi body is also tested through GRD/DAT export
and native import with its neutral PolySilicon identity. In this source the gate body is subsequently removed, so no
PolySi volume/material group is exported. Its constant-profile declaration
still references a known authored region and has no active volume. Bundle
provenance records it in `inactive_region_references` with the removal reason.

The bounded `(sdegeo:set-contact (find-body-id ...) name "remove")` form uses
creation-state ownership and an interior pick. It supports an unfilleted,
fully live Rectangle without later overlapping bodies at removal. Partially
overwritten/curved removals fail explicitly. Later bodies can fill a removed
vacancy; deleting a body never restores an older overwritten region.

Only surviving one-sided boundary curves become contacts. The BV gate is an L
in the half device and has total length 0.53 um after reflection. Its boundary
touches `R.Gateox`, `R.PolyReox`, and `R.PolyReox_mirrored`. Export/import retains
the complete adjacency list; legacy `region_id` is the smallest adjacent region
ID. This metadata convention does not qualify new solver boundary physics.

The source fillet radius remains 0.08 um. Maximum straight-edge chord error is
1e-5 um; each curved region's area budget is arc length times that error.
Geometry, contact endpoints, and mirror-axis merging use 1e-8 um tolerance.
The independent [BV contract](../scripts/sentaurus_bv_contract.py) checks source
areas, mirror topology/scalar fields, the absent gate, and nodal doping with
relative error <=1e-6. The corpus separately checks complete contact segments,
positive areas, connectivity, and per-region mesh spacing.

## Multibox and Erf

The supported multibox forms have six 2-D numerical size/ratio arguments and
a named Rectangle placement. Sizes must be finite, positive, and min<=max;
ratios are zero or have magnitude >=1. Unknown definitions/windows, invalid
ratios, and 3-D/non-Rectangle forms fail. A window clips the metric; overlapping
metrics use the finest directional constraint together with ordinary refinement.

For magnitude r>1, a positive ratio grades from the lower coordinate and a
negative ratio from the upper coordinate. The continuous local metric is
`min(h_max, h_min+(r-1)*distance)`. At geometric layer distances
`h_min*(r^n-1)/(r-1)` it gives `h_min*r^n`, until capped. Ratio zero imposes no
directional constraint. Ratio magnitude one uses **maximum size**, based on
live U-2022.12 controls: ratio=1 and max=max both gave 25 vertices/32 triangles;
min=min gave 576 vertices/1058 triangles. This resolves the conflicting
documentation for this implementation; no commercial topology identity is claimed.

The removed-gate multibox window has no surviving volume and is explicitly
inactive, including during adaptive splitting. The channel window remains
active. Bounded conforming refinement checks local x/y edge components against
1.35 times the metric and fails on resource exhaustion. Reports retain active
and inactive placements and their requested ratios/windows. Other constraints
can refine a graded zone further, so measured monotonic edge sizes are not
required where finer fields overlap. Mesh angles are reported; these meshes
do not have a non-obtuse numerical qualification.

Erf lateral diffusion is a continuous finite-baseline convolution. Its scale is
`Factor*sigma`, its endpoint weight is 1/2 for an isolated endpoint, and `Symm`
means both normal directions. The former outside-only Erf tail was incorrect
at endpoints and is corrected. Gaussian lateral profiles retain their existing
contract. See [the frozen gates](../configs/reference_bv_modeling_acceptance.json).

## Live oracle and verification

The user requested the compute VM's Sentaurus installation. Runs used the
existing private VM configuration and an isolated per-command license setting;
system/user profiles were unchanged. Private paths and job/process records stay
in ignored local logs. Commercial TDR, converted GRD/DAT, and logs are not Git
fixtures. Owned inputs are in [the control fixture directory](../tests/fixtures/sentaurus_p5).

The original BV source has an unbound `+z` when run directly by this SDE version.
A separate wrapper defines `(define +z "+z")` and loads the unchanged original.
That contextual run completed half meshing and x reflection. It does not prove
that the raw original executes standalone without preprocessing context.

The owned Erf control checked 306 native region-local values; BV checked 8,532.
Maximum peak-normalized errors were 4.352e-15 and 5.243e-15 respectively. Absolute
roundoff errors are recorded as well. Native region-local dataset ordering is
handled only by this oracle auditor, not newly enabled in the production grid
importer. Native circle tessellation is coarser than the approved chord budget;
its area is evidence, not the analytic area baseline for our finer mesh.

Reaudit fetched commercial outputs locally:

```powershell
python scripts/audit_reference_sde_p5_native.py `
  --input-dir "build-release/p5 native oracle" `
  --output "build-release/p5 native oracle/audit.json"
```

[Structured acceptance](validation/reference_sde_p5_acceptance_2026-10-30.json)
records source identities, numerical gates, per-run mesh hashes, native-oracle
statistics, and the actual regression/build outcomes. P1-P4 contracts and the
P2 eight-command checks remain required. P6 covers packaging/cross-platform
qualification separately.
