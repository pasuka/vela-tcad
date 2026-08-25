# Sentaurus SDE to Gmsh M1 acceptance

## Scope

M1 is the first public, offline preprocessing demonstration. It accepts the
checked-in two-region fixture with:

- two strictly separated `Rectangle` regions;
- one pick-point contact on the exterior boundary of each region;
- constant `NoReplace` doping placed by region or material;
- no refinement function.

It intentionally fails closed for touching or overlapping regions, ABA/BAB
ownership, inserted vertices, Gaussian profiles, replacement semantics and
refinement. Those capabilities remain assigned to M2 or M3. M1 is not a claim
that the full PN2D SDE input is supported.

## Public command

From an MSYS2 UCRT64 environment containing the pinned Gmsh 4.15.2 runtime:

```powershell
D:\msys64\ucrt64\bin\python.exe scripts\sentaurus_gmsh_pipeline.py `
  --sde tests\fixtures\sentaurus_gmsh_m1\two_regions_sde.cmd `
  --output-dir build-release\gmsh-m1 `
  --runtime-dir build-release\gmsh-m1-runtime `
  --mesh-size 0.25
```

The three runtime paths can also be fixed explicitly with `--gmsh-python`,
`--gmsh-library` and `--gmsh-executable`. Runtime discovery and staging reuse
the M0 contract and do not install packages or modify the host Python setup.

## Output contract

One successful run atomically publishes:

- `device.msh`: Gmsh MSH 4.1 ASCII with one mesh and five `NodeData` fields;
- `device.vtu` and `device.vtk`: ParaView/VisIt outputs with nodal fields and
  line/triangle physical ownership;
- `device.grd` and `device.dat`: generic DF-ISE grid and nodal datasets;
- `mesh.json` and `doping.csv`: Vela mesh and doping inputs;
- `device_ir.json`, `mesh_field_bundle.json`, `id_mapping.json`;
- `mesh_qualification.json`, `cross_format_report.json`, runtime manifest and
  final `run_manifest.json`.

`vela.mesh_field_bundle.v1` is the sole semantic source. Gmsh tags remain
positive and are preserved; Vela and generic DF-ISE identifiers are continuous
and zero-based, with explicit maps between the namespaces. Sentaurus material
name `Silicon` is mapped to Vela material name `Si` only in `mesh.json`.

The authoritative nodal fields are `DonorConcentration`,
`AcceptorConcentration`, `DopingConcentration`, `AbsDopingConcentration` and
`Log10AbsDopingConcentration`. Species-specific datasets are also written to
DF-ISE DAT where present.

## Acceptance gates

The registered CTest `sentaurus_gmsh_m1` checks:

1. the fixture stays inside the explicit M1 capability subset;
2. overlaps and touching interfaces fail closed with capability diagnostics;
3. Gmsh OCC produces linear, counter-clockwise Tri3 cells and both contacts;
4. MSH has one `$MeshFormat`, one `$Nodes`, one `$Elements` and five
   `$NodeData` sections;
5. MSH, VTU, VTK and DAT nodal fields round-trip to the bundle without numeric
   differences;
6. physical groups, continuous IDs, explicit maps, mesh quality and Vela
   material aliases satisfy the bundle contract;
7. the final manifest is written only after all exports and comparisons pass.
8. a second run with the same pinned runtime produces byte-identical semantic
   artifacts (the runtime and run manifests are path-dependent and excluded).

If the pinned runtime is absent, the three capability tests still run and the
real integration case is reported as skipped. A release acceptance run must be
performed on a host with Gmsh 4.15.2 and must not rely on that skip.

The final Windows UCRT64 Release acceptance used Gmsh 4.15.2, executed the real
integration case, and passed all 708 registered CTests in 209.44 seconds.
