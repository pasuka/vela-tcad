# Sentaurus SDE to Gmsh M3 acceptance - 2026-08-25

## Outcome

M3 is complete. The public, license-free Gmsh preprocessing pipeline now accepts
the phase-4A refinement subset and preprocesses the complete checked-in PN2D SDE
deck. M3 does not run IV/BV or TDX; those remain M4.

## Frozen semantics

- Sentaurus Structure Editor T-2022.03 pages 554-561 define Line, Rectangle and
  closed counterclockwise Polygon Ref/Eval windows, refinement placements,
  refinement sizes and the `MaxLenInt` argument order.
- Sentaurus Mesh T-2022.03 pages 28-29 define the first interface name as the
  refined side, `DoubleSide` as both sides, `UseRegionNames` as regionwise
  matching, `Factor >= 1`, and an omitted factor as a one-layer default.
- Gmsh 4.15.2 `MathEvalAniso` maps `max_x/max_y` to a 2-D metric. `MinAniso`
  intersects overlapping metrics, so the finest directional constraint wins.
- The Gmsh 4.15.2 official t17 policy is used for anisotropic 2-D meshing:
  BAMG (`Mesh.Algorithm=7`), `Mesh.AnisoMax=1000`, one thread, first-order
  triangles, no recombination and no boundary-size interpolation.
- Sentaurus minimum spacing is recorded as a refinement stop target. It is not
  passed to Gmsh as a hard lower clamp.

## Implementation

- `scripts/sentaurus_device_ir.py` now preserves source locations for refinement
  definitions and placements, validates closed counterclockwise Polygon windows,
  and emits structured `MaxLenInt` fields.
- `scripts/sentaurus_gmsh_backend.py` adds M3 validation, boolean-owner-aware
  region/material conditions, analytic Rectangle/Polygon/Line window metrics,
  `Distance`-based `MaxLenInt`, a versioned field DAG and measured edge-quality
  gates.
- `scripts/sentaurus_gmsh_pipeline.py` defaults to M3 and writes
  `refinement_field_dag.json` together with the existing cross-format bundle.
- `MaxTransDiff` and `MaxGradient` are parsed with their dataset name but rejected
  as phase 4B. MaxInterval, multibox, offsetting and MaxLenInt-to-contact remain
  fail-closed with source-located reports.

## Complete PN2D evidence

The checked-in
`reference_tcad/pn2d_sentaurus2022_oracle/input/pn2d_sde.cmd` was processed with
Gmsh 4.15.2 and M3:

- representative run: 1044 nodes, 1954 Tri3 and 20 contact Line2 elements;
- one connected component; all triangles counterclockwise;
- global request 0.05 x 0.05 um: P95 edge components
  0.05018 x 0.04884 um;
- junction request 0.01 x 0.02 um: P95 edge components
  0.01164 x 0.02261 um;
- both placements passed the explicit 1.35x BAMG conformance envelope;
- MSH, VTU, VTK, DAT and CSV donor/acceptor/net fields had zero serialized
  difference in the cross-format report;
- MSH Physical Groups, Anode/Cathode line elements, DF-ISE regions and all five
  visualization fields round-tripped.

The M3 refinement fixture additionally passed material, Polygon and Line
placements, overlapping metric intersection and regionwise `MaxLenInt` with
`DoubleSide` and `UseRegionNames`. Its interface had 33 Line2 elements and a
maximum edge length of 0.03061 um for a 0.03 um request.

## Verification

- UCRT64 Release configure: passed.
- UCRT64 Release build: passed (38 relink steps).
- M3 Python regression: 6/6 passed.
- SDE/IR/M1/M2/M3 focused regression: 63/63 passed.
- Release CTest: 710/710 passed in 303.47 seconds.
- The licensed PN2D expected IR hash seal was updated only for newly retained
  refinement source locations; the source SDE, evaluated log and official mesh
  artifacts were not changed.

## M4 handoff

M4 must use this Gmsh-generated PN2D `mesh.json/doping.csv` in actual Vela IV/BV
contracts and perform the internal TDX/DF-ISE/TDR round-trip. It must not replace
the M3 mesh with tensor-legacy or official TDR results when claiming Gmsh E2E
acceptance.
