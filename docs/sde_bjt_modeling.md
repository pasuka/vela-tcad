# BJT source modeling (P3)

This records the qualified BJT subset. [P4](sde_mos_modeling.md) extends adaptive
refinement to multiple CAD surfaces and region/material/window scopes; the BJT
numerical gates remain unchanged.

The two unmodified sources `genius_bjt_sentaurus2022/source/bjt_sde.cmd`
and `bjt_sde_local_refined.cmd` under `reference_tcad/` are supported by the
offline Gmsh preprocessor. Corpus revision 2 contains five generated cases and
three explicit rejections. This qualifies geometry, doping and exported mesh
contracts; numerical BJT IV results require their own validation.

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts/sentaurus_gmsh_pipeline.py `
  --sde reference_tcad/genius_bjt_sentaurus2022/source/bjt_sde.cmd `
  --output-dir "build-release/my BJT" `
  --runtime-dir build-release/gmsh-python-runtime
```

`device.grd`, `device.dat` and `device.msh` are the primary outputs. The
normalized bundle, ID maps, qualification and cross-format reports remain next
to them. For the independent BJT gates, use the
[corpus runner](reference_sde_cases.md), which also imports GRD/DAT through the
native consumer and compares contacts, materials, topology and fields.

## Source and IR semantics

- `(define silicon (sdegeo:create-rectangle ...))` executes the body creation
  at that source location and returns a typed `body:N` handle. Aliases preserve
  identity. Handles cannot enter scalar arithmetic. Creation order, boolean
  mode and source maps survive v1-to-v2 migration. Legacy numeric/string
  bindings keep their original representation.
- A named helper may contain only supported numeric operators and position
  construction, or calls to existing pure helpers. It has at most eight
  parameters, 32 levels of expression nesting, 16 levels of calls and at most
+  4096 total calls per source. Cross-kind rebinding fails closed. It cannot
  execute geometry mutations, I/O, control flow or arbitrary Scheme.
- The six-argument `sdedr:define-analytical-profile` supports `General` using
  device x/y/z coordinates (z=0 in 2-D). The initializer contains semicolon
  separated scalar assignments; the first assignment wins in source order.
  Supported expressions use finite constants, names, `+ - * / ^`, unary signs,
  `abs`, `exp`, `sqrt` and `log`. A whitelist AST has at most 256 nodes, 32 levels
  and 4096 source characters; exponent magnitude is bounded by 64. No Python
  eval/exec, attributes, indexing, imports or user-defined calls are allowed.
- General placements currently require a Line reference, `Both` and `Eval`;
  explicit zero-decay restrictions and ordered replacement remain available.
  Other modes fail closed. Arithmetic/domain failures use the explicitly
  authored default and are counted in `provenance.general_evaluation`.
  A negative concentration fails. The BJT contract requires zero defaults.
- Both Device IR schemas describe the tagged bindings and General AST. The
  backend revalidates the AST against its source expression before meshing.

## Bounded MaxTransDiff mapping

The supported subset is one Rectangle body, a region-scoped refinement,
`DopingConcentration`, a finite positive Value and both directional minima.
Other variables/scopes fail closed. The
[Sentaurus mesh training](https://ghzphy.github.io/Sentaurus_Training/smesh/smesh_02.html)
describes a hyperbolic-arcsine transform; Vela's explicit contract is
`abs(asinh(Nnet(a)/(1 cm^-3))-asinh(Nnet(b)/(1 cm^-3))) > Value` on each unique
triangle edge. The BJT source Value is 0.8. General device-coordinate functions
are described in the
[analytical-profile reference](https://ghzphy.github.io/Sentaurus_Training/smesh/smesh_04.html).

The initial BAMG mesh uses all original max-spacing placements. Subsequent
rounds split the failing edges at their midpoints with conforming red/green
triangles; adjacent triangles share exactly one new vertex per edge.
Longest-edge closure protects triangle shape. All CAD curves, endpoint
classifications and physical groups are preserved, including non-contact
boundary segments. Only affected triangles and the closure neighborhood split;
there is no global finest-size replacement.

An edge above Value stops only when **both** absolute x/y components are no
larger than the finest applicable min_x/min_y at its midpoint. Minima are stop
targets, never hard lower clamps. Such edges can still have a large transformed
difference across a junction: their counts and maximum defect are recorded
separately. Success means zero remaining refineable edges, not that every edge
meets Value. Default budgets are 12 refinement rounds and 200000 nodes:

```powershell
python scripts/sentaurus_gmsh_pipeline.py `
  --sde reference_tcad/genius_bjt_sentaurus2022/source/bjt_sde.cmd `
  --output-dir "build-release/budgeted BJT" `
  --maxtransdiff-max-iterations 12 --maxtransdiff-max-nodes 200000
```

The `maxtransdiff` section in the bundle provenance and mesh qualification
records every round's nodes, triangles, edge defects, marked/closure edges and
minimum-stop counts. Node growth is checked before a split. Budget exhaustion
fails and publishes no success outputs. This is an explicitly auditable local
mapping, without a claim of native snmesh triangulation equivalence or PN2D BV
non-obtuse qualification.

## Frozen acceptance

The user approved the numerical gates in
[`reference_bjt_acceptance.json`](../configs/reference_bjt_acceptance.json)
before the first mesh. The independent oracle
[`sentaurus_bjt_contract.py`](../scripts/sentaurus_bjt_contract.py) reads the
hash-sealed geometry contract and implements the outside-window-distance
Gaussian formula separately from the General AST evaluator.

- Region bounds [0,0,6,2] um, area 12 um², Silicon; base [1.25,2] and emitter
  [2.75,4.25] at y=0, collector [0,6] at y=2. Complete contact spans and lengths,
  one component, CCW positive triangles and original sizing P95 gates apply.
- Every node: donor/acceptor/net relative error <=1e-6, with a 1 cm^-3 absolute
  floor for near-zero values. At least 1000 total nodes and two per contact.
- Donor and acceptor integrals: exact product-erf integral versus the exact
  integral of the mesh's piecewise-linear nodal field, relative error <=1%.
  Reports use cm^-3·um²; multiplying by 1e-8 gives cm^-1. No extrusion width or
  numerical solver state is introduced.
- Vertical cuts x=1.625,2.6,3.5 um: all net-zero crossings in y=[0,2] must match
  in count and position. Analytic roots are bracketed independently, mesh roots
  use triangle-linear intersections. Error cannot exceed the source's finest
  local minimum y size (0.004 um in the base case; 0.002/0.0015 um in applicable
  locally refined windows).
- The local variant must reduce P95 edge length in all three authored local
  windows against the base mesh. Each run also passes cross-format and native
  import checks. BAMG repeats use `contract_equivalent`; each run's distinct
  hashes are retained and byte determinism is not claimed.

`sentaurus_bjt_modeling` registers AST/security, typed-binding/migration, oracle,
local split and budget-failure tests. `reference_sde_cases` exercises both real
scripts, repeats and independent imports. Downstream `sde_sdevice_checks` still
covers its eight named PN2D/gate/Schottky check-only inputs; added BJT mesh
support does not claim BJT command or numerical qualification.

TDR, raw logs and generated grids stay in ignored build directories. Public
evidence contains structured metrics and hashes. No live SDE/TDX or BJT solve
is required for these gates. The next development scope is P4 curved geometry
and complete SingleDevice/TransportModels modeling; their present rejection is
`fillet-2d` at source line 120. BVmethods remains blocked by its PolySi material
boundary until P5.

Measured results and per-run hashes: [P3 acceptance](validation/reference_sde_p3_acceptance_2026-10-29.json).
