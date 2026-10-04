# SimpleMOS engineering reference

Current qualified scope: **16 curves / 816 points**, n17–n24, Vd=0.05/1 V,
Vg=0:0.05:2.5 V, 300 K, existing process meshes. Both full cold startup and
independent M60 native-field initialization passed the frozen numerical gates.
Maximum absolute Id error versus M60 default port extraction is **0.186617394672%**
(n23, Vd=0.05 V, Vg=0.05 V). Later T470p Release regression passed **993/993**.
This directory packages existing evidence; its preparation did not rerun a solver.

## Configuration and meaning

- `engineering_contract_v1.json`: original engineering profile and frozen input
  hashes, unchanged. `inputs/n17` through `inputs/n24` contain byte-identical
  templates; `inputs/materials.json` contains their shared material definitions.
  These are driver templates, not standalone ready-to-run decks: mesh/doping
  paths are relative to the external `data/inputs` directory, and the driver
  sets simulation type, bias, state, and output paths.
- Physics: Boltzmann, OldSlotboom, PhuMob, Enormal, high-field saturation,
  doping-dependent SRH. Geometry: element-box Si transport, cell-material
  permittivity, signed Si Poisson and SRH source volumes, Delaunay transfer.
  Complete split state, binary128 Poisson residual, stable merit comparison,
  exact Dirichlet updates, and four linear refinement iterations are explicit.
  The executor selects **SparseLU**, with one linear/OMP/OpenBLAS thread.
- M60 reference: binary64, Digits=8, ErrRef(Electron/Hole)=100 cm^-3,
  default RhsMin=1e-5 and default terminal algorithm. EP128 is a separate
  calibration subset, not the reference for the entire matrix.
- Gates: Id error 2%; dual potential 1e-6 V, density relative 1e-4, current
  relative 1e-6; all-carrier-row ratio 1e-6, KCL/Id and port discrepancy 1e-8;
  density/SRH reconstruction 1e-12/1e-7. No gates or global defaults changed.
- Potential units are V, current A/um, exported density cm^-3 converted to
  state density m^-3. Legacy material concentration fields follow unit_scaling;
  do not infer their units from field names alone. Native field differences
  are descriptive, not additional acceptance thresholds. SRH weighted L1 is
  not a per-node relative error at zero crossings. Mobility field definitions
  still need separate coverage.

## Read-only entry (fresh checkout, no simulator needed)

From the repository root, with Python 3:

```text
python reference_tcad/simplemos_sentaurus2022/engineering/verify.py
```

This checks the byte manifest, 816 unique comparison keys, 4,896 field rows,
unchanged profile/gates, source ownership, recorded maxima and regression exits.
It verifies shipped summaries, not the absent large raw states or today's binary.
Optionally audit a restored input bundle as well:

```text
python reference_tcad/simplemos_sentaurus2022/engineering/verify.py --inputs <base>/data/inputs
```

## Numerical run entry (requires externally restored inputs)

Do not start from the earlier M3/M4 commands to reproduce this engineering run.
Use an isolated ignored `<base>` with this layout:

```text
<base>/data/engineering_contract.json  <- engineering_contract_v1.json
<base>/data/inputs/                    <- complete frozen input bundle
<base>/data/m60_joined.csv             <- m60_joined.csv
<base>/data/frozen_points.csv          <- frozen_points.csv (six-seed preflight)
<base>/data/seeds/*.h5                 <- six original complete split checkpoints
<base>/native/bundle/n*/              <- M60 field TDRs and IdVg PLT files
<base>/native/m60fields_<case>.exitcode
```

The complete input bundle needs mesh.json, doping.csv, geometry.json and two
original reference CSVs per device, plus materials.json, contract.json and
hashes.json. Its hashes must match `frozen_input_hashes`; hashes.json is that
mapping. The checked-in templates/materials/contract are exact copies, not
replacements for the missing meshes. The six preflight seeds and native export
archives must retain their original transfer manifests; check those before use.
For each case, the native bundle needs 51 `m60fields_<case>_state_*_des.tdr`
files and `IdVg_m60fields_<case>_des.plt`, with its zero exit record. These
licensed/generated inputs are intentionally not in Git. Do not fabricate seeds
or successful preflight records. Restore them from the authorized campaign
archive; without them, only the read-only check above is available.

The earlier portable input contract and its `_reference.csv` files retain
**original-M8** provenance. They support the cold-path diagnostic ledger;
the current engineering acceptance uses `m60_joined.csv`'s
`tight_default_Id_A_per_um` via the engineering matrix driver. Other columns
of that CSV, including proposed thresholds and earlier Vela currents, are
historical diagnostics, not new acceptance conditions.

Build `vela_example_runner` and `sentaurus_import` with actual HDF5/TDR support,
using the repository's platform build instructions. Python simulation drivers
need NumPy and h5py. Record source revision, binary/importer hashes and compiler
configuration for any new run; it is a new qualification, not the frozen run.
Then use the existing production validation entry points:

```text
python scripts/run_simplemos_t470p_preflight_20260929.py --base <base> --runner <absolute-runner>
python scripts/run_simplemos_engineering_matrix_20260929.py --base <base> --case n17_vd_1 --runner <absolute-runner> --importer <absolute-importer>
```

Despite its historical name, the preflight accepts a platform-specific runner.
Require all six preflight points before running cases. For the full matrix,
run each of the 16 unique case names in `inputs/contract.json`; preserve failure
records and stop dispatching on failure. Outputs go to `<base>/matrix/`.
Cold startup is nonlinear Poisson → coupled equilibrium → drain ramp → gate
sweep. The driver then independently solves from each native export and checks
dual states, all rows, ports, current and fields. A saved cold checkpoint alone
is not a passed dual-initialization point. This packaging task starts no jobs.

## Evidence and provenance

- `summary.json`: current combined closeout, implementation commit `cfcf13ff`.
- `evidence/joined/`: unchanged 816-point comparisons, six-field table, source
  ownership/hashes, 204-core-file hashes, original summary and byte seal.
- `evidence/regression/`: later 993/993 result, exits, source seal and retrieval
  verification. Logs and raw states remain outside Git.
- `evidence/cloud_source_verified.json`: Linux binary/importer identities;
  Windows runner identity is in the regression status. These are historical
  execution identities, not a claim that a new local build has those hashes.
- `sha256.json`: exact byte hashes of this small package, excluding itself.
  This protects consistency, not authenticity against a coordinated rewrite.

The original joined summary predates the provenance repair and still says
`full_regression_passed=false`. It is preserved as history; the later regression
record supersedes that status. The original n17 high-Vd retry exit was empty;
the independent cached replay passed without changing 162 states. This is not
recovery of the original process exit code. The older bad_alloc failure was not
reproduced in the successful retry. A later Windows Event 2004 was recovered
within 0.16 seconds of the failure record, with system committed memory near
its limit. This supports resource exhaustion; the largest process's owning
task has not been recovered.

The subsequent [mobility field audit](../../../docs/validation/simplemos_mobility_field_audit_2026-09-30.md)
covers all 816 existing states without new DC solves. Same-export cell mobility
maximum relative differences are 1.628364e-6 (electron) and 9.621782e-6 (hole).
The unchanged 1e-7 cell gate passes 805/816 and 0/816 states respectively; these
are separate from the qualified engineering current gates. Cell reconstruction
and full-split production edge mobility agree within 3.01e-13. Native node
display mobility is not interchangeable with conservative edge mobility; its
exact display projection remains unqualified. Raw audit outputs remain outside
this frozen numerical evidence package.

See the [matrix closeout](../../../docs/validation/simplemos_engineering_joined_closeout_2026-09-30.md)
and [subsequent provenance/commit review](../../../docs/validation/simplemos_provenance_and_commit_review_2026-09-30.md).
One historical source blob remains unrecovered as recorded there. Core source
qualification precedes a whitespace-only cleanup; binary identity is not
inferred from a configuration hash. No claims cover other meshes, temperatures,
Fermi statistics, quantum, avalanche, or arbitrary mobility combinations.

## Delivery checks

The package and all restored input hashes passed the read-only audit. The
reference-package, matrix, engineering-audit, join, Codespaces-controller and
exit-capture Python tests passed. A clean dependency copy without ignored input
data passed the package tests and loaded both documented run commands with
`--help`. The Release preset configured successfully and the newly registered
`simplemos_reference_package` CTest passed (1/1). This added check is not
retroactively included in the earlier 993-test run; no new full CTest or solver
simulation was run for packaging. The local driver import emitted an installed
h5py/HDF5 version warning; the standard-library package audit is independent of
HDF5, and this check does not qualify that local stack for numerical runs.

```text
python -m unittest tests.regression.test_simplemos_reference_package
ctest --preset windows-ucrt64-release -R "^simplemos_reference_package$" --output-on-failure
```
