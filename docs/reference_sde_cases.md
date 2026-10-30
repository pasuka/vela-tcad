# Offline SDE case regression

The P1 corpus is defined in
[`configs/reference_sde_cases.json`](../configs/reference_sde_cases.json) and its
[versioned schema](../schemas/vela.reference_sde_cases.v1.schema.json).
Corpus revision 3 covers eight current primary source scripts: seven generate meshes and one
is expected to fail at a specific source location. Archived variants, oracle
wrappers and retired cases are excluded.

## Run from the checkout root

Use the installed Gmsh Python API **4.15.2**, matching library and executable;
see [runtime discovery](gmsh_python_runtime.md). Missing dependencies fail the
batch and never count as an expected unsupported-script rejection.

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts/sentaurus_sde_cases.py `
  --output-dir "build-release/reference SDE cases" `
  --runtime-dir build-release/gmsh-python-runtime
```

On Linux the same command uses the selected Python interpreter and paths, for example:

```bash
python3 scripts/sentaurus_sde_cases.py \
  --output-dir 'build/reference SDE cases' --runtime-dir build/gmsh-python-runtime
```

The batch reuses the existing SDE-only CLI for each case. To generate a single
device without the corpus assertions:

```powershell
python scripts/sentaurus_gmsh_pipeline.py `
  --sde reference_tcad/gate_reference/source/gate_sde.cmd `
  --output-dir "build-release/my gate mesh" `
  --runtime-dir build-release/gmsh-python-runtime
```

The primary products are `device.grd`, `device.dat`, `device.msh`; helper files
include the normalized bundle, Vela mesh/doping, ID mapping and cross-format
report. The batch writes one case directory per ID and `case_summary.json`.
Exit 0 means the corpus met its declared expectations; **it means 7 generated
and 1 expected rejection**, not eight supported device scripts. Exit 2 means
a dependency, input identity, contract, determinism or rejection boundary failed.

## Contract and evidence

- Source identity is SHA-256 of UTF-8 source bytes with CRLF converted to LF;
  no other whitespace/content changes are normalized. This portable identity
  is explicit as `source_sha256_lf`. Every invocation also records the raw input
  hash and verifies that the run did not change the file.
- The material hash seals canonical JSON of the **SDE material-name to neutral
  material-name mapping** used by each case. It does not qualify numeric `.par`
  constants or a solver material library; those belong to P2 and numerical gates.
- Geometric contracts use analytic device bounds, each region's area, contact
  position and complete edge length. Coordinates and shared-window boundaries
  use `1e-12 um` tolerance. Nodal donor/acceptor/net values follow the declared
  analytic constant or PN window profiles, including both placements at the
  shared junction. These scalar comparisons use `rtol=1e-12, atol=1e-12`.
  The two BJT cases use the separately frozen [P3 contract](sde_bjt_modeling.md):
  original General expressions, pointwise relative error <= 1e-6, species
  integral error <= 1%, and cut junction error within the local minimum y size.
- Mesh quality additionally requires one connected component, positive triangle
  area/CCW orientation, and the per-region P95 absolute x/y components of unique
  edges within 1.35 times the source max spacing. This retains the existing M3
  BAMG conformance envelope and independently measures it; it does not qualify
  gate_reference for the PN2D non-obtuse solver profile.
- Existing independent GRD/DAT, MSH, VTU/VTK and CSV readers check cross-format
  consistency. The consumer `vela_sdevice.mesh_import` independently reads GRD/DAT;
  the batch compares coordinates/IDs, region-labeled topology, materials, full
  contact edges/nodes/ownership and donor/acceptor/net doping. Interface boundary
  metadata is checked separately against the original interface names and edges;
  it is not included in volume materials (see the P2 import correction).
- By default supported cases run twice (`--repeat 2`), including separate native
  imports. PN2D and Schottky require byte-identical GRD/DAT/MSH/mesh JSON/doping CSV.
  The gate_reference, BJT and MOS policy is `contract_equivalent`: each run must
  independently meet the geometry, material, contact, analytic doping and mesh
  quality contracts. The 4.15.2 BAMG wrapper traverses pointer-ordered vertices;
  seed/thread/Reproducible experiments did not give byte-identical results here.
  All per-run hashes, node/cell counts and measured sizes are retained, with an
  explicit `byte_identical_observed` flag; this flag is not a future guarantee.
  See the [pinned upstream implementation](https://raw.githubusercontent.com/gmsh-project/gmsh/gmsh_4_15_2/src/mesh/meshGFaceBamg.cpp).
  `--repeat 1` disables comparison between runs.
- Each case records parse, model, export, native-import and reference stages.
  Reference is `not_run`: no live SDE/TDX, `sd-sim --check` or physics solve is
  performed. This is P1, not the subsequent full SDevice integration qualification.
- The remaining rejection must retain exit 2, `sde_parse`, source basename, exact
  line/column and a diagnostic category. An unexpected generation also fails the
  corpus so that new support receives a reviewed positive contract.

## Safe repeated runs

Use a dedicated output directory per standalone device and do not run concurrent
writers into it. A prior valid pipeline manifest plus matching file hashes proves
ownership. Its generated files and success marker are removed before a rerun;
failure leaves no old generated grid pretending to be a new result. Recovery
after a failure regenerates outputs and removes the old failure report.

Unrelated files are kept. A reserved filename with no valid ownership manifest,
a modified previous output, a linked reserved file, or an input/output collision
fails with `output_safety` on stderr. The directory is left untouched in that
case, including its prior manifest; verify the manifest's hashes before using
any old result. Do not delete user files to make this guard pass: use a fresh
directory or explicitly preserve and relocate your edits. The existing combined
SDE/SDevice importer retains its separate output publication policy.

The batch invalidates its own previous summary before a new input check, so an
invalid manifest cannot leave a stale passing batch summary. Sources and
numerical thresholds are not changed by either CLI.

## Regression

```powershell
cmake --preset windows-ucrt64-release
ctest --test-dir build-release --output-on-failure -R '^reference_sde_cases$'
```

The test covers the entire corpus, path names with spaces, repeats,
success/failure/success recovery, input and user-file protection, malformed
ownership manifests, hash/schema drift and deliberate doping/contact corruption.
Developer machines without Gmsh skip the integration class explicitly; the batch
CLI itself fails. CI first runs the mandatory Gmsh probe before CTest.

Keep all generated grids, TDR and raw logs in ignored build directories. Only
source, small contracts, tests and structured acceptance summaries belong in Git.

The subsequent [P2 check-only workflow](sde_sdevice_checks.md) binds these outputs
to actual SDevice commands without running a numerical solver.

P4 adds the two complete MOS meshes with explicit `node=0` parameters, certified
fillet chords and independently imported half/full outputs. See [P4](sde_mos_modeling.md).
