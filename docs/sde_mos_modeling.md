# MOS source modeling (P4)

Both unchanged primary sources are supported by the offline Gmsh pipeline:
`singledevice_sentaurus2018/source/singledevice_sde.cmd` and
`transportmodels_sentaurus2022/source/nMOS_dvs.cmd`, under `reference_tcad/`.
Corpus revision 3 has seven generated cases and one expected BVmethods rejection.

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts/sentaurus_gmsh_pipeline.py `
  --sde reference_tcad/singledevice_sentaurus2018/source/singledevice_sde.cmd `
  --parameter node=0 `
  --output-dir "build-release/my MOS" `
  --runtime-dir build-release/gmsh-python-runtime
```

`device.grd/dat/msh` describe the complete reflected device. Independent half
products are `half_device.grd/dat/msh`, `half_mesh_field_bundle.json` and
`half_mesh_qualification.json`. All are generated locally and ignored by Git.
Both sets retain geometry, material and contact identities and nodal fields.
The corpus runner supplies its sealed `parameters.node=0` explicitly.
Unknown/unresolved placeholders, unsafe values and arbitrary shell commands fail.
The paired `sentaurus_import.py device` entry uses its existing `--template-var`
arguments for SDE placeholders as well as SDevice preprocessing.

## Geometry and refinement contract

The user approved [the numerical contract](../configs/reference_mos_acceptance.json)
before the first mesh: original radius 0.05 um, maximum chord error 1e-5 um,
pointwise doping relative error 1e-6, and axis/contact tolerance 1e-8 um.
Each curved region's area error budget is arc length times chord error. The
independent [oracle](../scripts/sentaurus_mos_contract.py) checks analytical
implant tails, incident semiconductor regions, curve geometry and reflection.

One unambiguous convex Rectangle corner per body can be rounded. The IR retains
the original radius, corner and source location. A deterministic chord chain has
enough segments that `r*(1-cos(theta/2)) <= 1e-5 um`; it is integrated into the
conformal OCC fragmentation before meshing. Unsupported/oversized fillets fail.
This is a certified polygonal representation of the circular boundary, not an
assertion that straight Tri3 edges are exact analytic arcs. Tensor legacy mode
rejects fillet/mirror inputs instead of silently retaining a rectangle.

The following bounded shortcut forms lower into ordinary, source-located IR:

| Source form | Lowering |
| --- | --- |
| refine-box name p0 p1 h, or hx hy hz | Rectangle window plus directional max-size placement; z must be zero; z spacing is parsed/validated but inactive in 2-D |
| refine-doping h | Minimum stop spacing h, net-doping MaxTransDiff Value 1, material Silicon |
| refine-interface region/material first second hlocal h, optional factor | MaxLenInt on the first side; factor defaults to 2; no implicit DoubleSide |

These defaults and scope conventions follow the official *Structure Editor User
Guide*, command entries 582–585 in N-2017.09:
[reference manual](https://kolegite.com/EE_library/books_and_lectures/%D0%90%D0%B2%D1%82%D0%BE%D0%BC%D0%B0%D1%82%D0%B8%D0%B7%D0%B0%D1%86%D0%B8%D1%8F%20%D0%BD%D0%B0%20%D0%9F%D1%80%D0%BE%D0%B5%D0%BA%D1%82%D0%B8%D1%80%D0%B0%D0%BD%D0%B5%D1%82%D0%BE%20%D0%B2%20%D0%95%D0%BB%D0%B5%D0%BA%D1%82%D1%80%D0%BE%D0%BD%D0%B8%D0%BA%D0%B0%D1%82%D0%B0/Sentaurus_manuals/PDFManual/data/sense_ug.pdf).
Optional forms outside this subset are rejected. A shortcut with no authored max
spacing uses a coarse nonbinding size, intersected with the existing base field
and all authored windows; it does not impose a uniform finest-size grid.

The [P3 adaptive contract](sde_bjt_modeling.md) now applies per CAD surface with
region-local doping, and per active region/material/window refinement. Shared
interface edges are split once, preserving all incident cells and CAD curve
classifications. MaxLenInt interface edges are also split until their target
length is met. The existing 1.35 sizing qualification envelope is unchanged.
Budgets still fail without publishing stale success artifacts.

## Local mesh reflection

Only the exact recognized `tdx -mtt -y -M 0 -S 0 -ren drain=source` operation with
the matching preceding half-mesh prefix is accepted. It becomes a typed local
operation and never launches TDX or a shell. It must be the final command.

The mirror copies mesh nodes and all fields, merges axis nodes, reverses triangle
winding, and rebuilds IDs and Physical Groups. Regions connected by an axis edge
keep their names; disconnected reox/spacer regions gain `_mirrored`. The mirrored
drain becomes source; gate and substrate contacts are continuous. Full and half
meshes are separately recorded. Refinement/field-DAG statistics refer to the
authored half before reflection; full geometry areas have a separate audit.

Both originals produce bounds [-0.1012,-0.185,1,0.185] um, seven regions and
source/drain/gate/substrate contacts. The independent oracle checks reflected
triangle identities, all scalar/species field parity, axis uniqueness, contact
spans, exact rectangular areas and bounded circular area error. Full GRD/DAT and
half GRD/DAT are independently imported through `vela_sdevice.mesh_import`.
The complete corpus also reads back MSH, VTU/VTK and CSV. Repeats require the same
contracts; BAMG byte determinism is not claimed and per-run hashes are saved.

## Consumer boundaries and qualification

The original SDevice physics and solve statements are preserved in Grid-only
local copies for `sd-sim --check`. Their current rejection boundaries are tested:

- `singledevice_deep_off_sdevice.cmd`: E403, because the original parameter input
  does not define PolySilicon numeric material constants.
- `pp6_des.cmd`: E105/E109 for unsupported Thermode and Solve System statements.

These are explicit command-input limitations alongside successful mesh imports.
No Silicon-to-PolySilicon parameter substitution or removal of physical options
was introduced. The existing eight P2 commands remain a separate positive gate.
No solver or TDR importer is invoked by the command boundary tests.

Thin dielectric cells may be obtuse: angles are reported and these meshes have
no PN2D mixed-Voronoi/non-obtuse or IV/BV qualification. Live commercial SDE was
attempted on the configured verification VM but its license was unavailable;
that attempt is not counted as validation. No T470p/Linux qualification or
numerical MOS solve is claimed. Raw transcripts and private job/path information
remain outside version control. No proprietary implementation source is included.

The next source-modeling scope is P5: PolySi aliases, creation-state body
contacts/removal and multibox for BVmethods. P4 does not unblock its numerical
NonLocalPath/external-resistor qualifications.

Measured results and per-run hashes: [P4 acceptance](validation/reference_sde_p4_acceptance_2026-10-30.json).
