# SDE offline modeling qualification (P6)

This is an independent modeling-tools distribution. Its products are
`device.grd`, `device.dat` and `device.msh`. The existing SDevice package and
reserved `sd-sim -d`, `-n`, `--max_threads` interface remain unchanged.
P6 local qualification passed on 2026-10-30; see the
[structured acceptance record](validation/reference_sde_p6_acceptance_2026-10-30.json).
The user explicitly selected local validation only: no push or actual remote CI
run was performed. The workflow definition remains available for later execution.

## Actual local results

Modeling-package source: `e1b83fd3e7cebf549d917aaf1c876db138bb0f0e`.
Both packages were assembled from a clean source tree. The final modeling
payload hashes still match the source modules and schemas in the worktree. Windows UCRT64
GCC 16.2.0 and Ubuntu 26.04.1 GCC 15.2.0 independently built the runtime targets;
Eigen 5.0.1, HDF5 and SuiteSparse UMFPACK/SPQR were enabled. The isolated
CPython versions are 3.12.10 (Windows) and 3.12.12 (Linux).

Each final package ran all eight primary sources twice: **32/32** modeling,
cross-format and independent GRD/DAT-import contracts passed. Each of those
32 meshes also passed the commercial geometry/contact/field/integral/cut/size
comparison. The eight native SDE+TDX source runs all returned zero. Different
unstructured runs retain their individual hashes, without byte-determinism claims.

| Case | Windows full-device nodes (two runs) | Ubuntu full-device nodes (two runs) | Worst donor/acceptor integral difference |
| --- | --- | --- | --- |
| pn2d | 1482 | 1482 | 0.1089% |
| schottky | 319 | 319 | 0.0000% |
| gate_reference | 163–165 | 163–164 | 0.0000% |
| bjt | 22413–22443 | 22366–22912 | 0.0295% |
| bjt_local | 42854–43130 | 44113–44117 | 0.0165% |
| singledevice | 60417–63048 | 58042–59839 | 0.8858% |
| transportmodels | 43061–47128 | 46073–47055 | 0.1633% |
| bvmethods | 43175–45257 | 44235–44571 | 0.3149% |

Across both platforms/runs, native pointwise relative error was at most
`2.65825e-12` and comparative integral error at most `0.885824%`, within
the frozen 1e-6 and 1% gates. Nonempty source-derived cut junctions passed their
local-spacing limits. Original radii, own chord limit and geometric/contact
tolerances were retained.

The final Windows second-run GRD/DAT pairs completed **8/8** licensed TDX
encode/decode loops. The mapped topology, materials, contacts, interfaces and
stored scalar/species fields passed; worst stored-field relative error was
`5.70112e-16`. Input and readback hashes identify these exact final meshes.
TDR and raw execution records remain private/ignored.

The final source revision passed **18 distinct CTest tests per platform**
(17/17 related checks plus P2 1/1). The final P6 selection is 3/3 CTest, covering 19 unit tests
on each platform. P2's eight real command checks and rejection tests passed on
both independently built consumers. One no-importer negative case conditionally
skips because the explicitly built importer is available; required Gmsh corpus
checks did not skip. The initial Linux old-package identity failure was corrected
and rechecked. An initial 1.069% integral failure motivated the Gaussian accuracy
control; a later BV iteration failure motivated local red closure. Neither
acceptance tolerances nor the 12-round/200000-node budget were relaxed.

Windows additionally generated the PN2D mesh with the package and outputs in
Unicode directories. Workflow YAML structure was checked locally; **actual CI
was not run**, following the user's no-push decision. Full C++ CTest and new
IV/BV/NonLocalPath/R2/R3 numerical qualification were not performed.

The final archives and checksum/manifest sidecars are local ignored products:
`build-release/p6/packages/vela-sde-p6-final-win64.zip` and
`build-release/p6/packages/vela-sde-p6-final-ubuntu26.04.tar.gz`.
Their exact hashes are in the structured record.


## Build and supply the package

Use Windows UCRT64 or Ubuntu 26.04 x86-64, following the repository build
instructions. `sd-sim`, `vela`, `vela_example_runner` and `sentaurus_import`
are independently built for runtime regressions. They are not included in the
modeling-tools archive. HDF5 and SuiteSparse availability must be read from the
actual CMake output, rather than inferred from options.

An explicit acquisition step downloads and verifies the pinned inputs in
[the offline input lock](../configs/sde_offline_inputs.json): official CPython
3.12 and six dependency wheels, plus matching upstream license texts.
The assembler only consumes these local inputs and does not install packages
or access the network. Linux additionally requires `patchelf`.

```bash
python3 scripts/prepare_sde_offline_assets.py --assets build-release/p6-assets
python3 scripts/package_sde_modeling.py --assets build-release/p6-assets \
  --output build-release/packages/vela-sde-ubuntu26.04.tar.gz
```

On Windows use UCRT64 `python` and a fresh `.zip` output name. The assembler
refuses to overwrite an archive. It writes SHA-256 and manifest sidecars;
`PACKAGE_MANIFEST.json` records every payload hash, input hashes, Python version,
source commit and dirty-tree state. A non-Git source snapshot must explicitly
supply `--source-commit <40-hex>` and `--source-tree-state clean|dirty`.
For example, Windows worktree `.git` paths are not automatically meaningful
inside WSL. A dirty snapshot is an explicitly identified prerelease artifact.

The archive includes Vela's Apache-2.0 license, upstream Python notices,
wheel licenses and notices for copied Linux libraries. **Gmsh is external**;
its wrapper, shared library and executable must all be supplied explicitly at
version 4.15.2. The runtime manifest verifies their versions and records paths
and hashes locally. Consult the [upstream Gmsh license](https://gmsh.info/doc/texinfo/#License)
for its terms. Neither a Gmsh binary nor proprietary Sentaurus software,
reference decks or generated simulation outputs are included in the archive.

For qualification, the separately acquired SDK/wheel is hash-locked in
[the external-runtime lock](../configs/sde_external_runtimes.json):

```bash
python3 scripts/prepare_sde_external_gmsh.py \
  --assets build-release/gmsh-inputs --runtime build-release/external-gmsh
```

Read `external-paths.json` to obtain the three paths. Windows source regressions
need both SDK `bin` and `lib` on `PATH` **when CMake configures the tests**.
The package entry sets its external-library search paths itself. Linux's Gmsh
wheel CLI needs its virtual environment's `bin` on `PATH`.

After relocating the extracted directory (spaces and Unicode paths are tested):

```bash
<package>/bin/python/bin/python3 -B -I -S <package>/run_sde.py \
  --sde <original-sde.cmd> --parameter node=0 --output-dir <fresh-output> \
  --gmsh-python <gmsh.py> --gmsh-library <libgmsh.so.4.15> \
  --gmsh-executable <gmsh>
```

On Windows substitute `<package>/bin/python/python.exe` and the SDK DLL/EXE.
Only supply parameters declared by the source. The entry verifies payload
hashes before execution and ignores injected `PYTHONPATH`/`PYTHONHOME`.
Missing or inconsistent runtime components fail explicitly.

## Integral convergence profile

The raw pipeline preserves its existing default (`--integral-relative-target 0`).
The offline package enables `0.00025` for Gaussian-containing unstructured
models; constant and General-only BJT models retain their existing paths.
The [profile](../configs/reference_sde_integral_control.json) controls separate
donor and acceptor integrals over semiconductor regions, including polysilicon.

Degree-five triangle quadrature estimates the local error of the linear nodal
field. Summed absolute cell errors prevent cancellation; bulk marking selects
the largest contributors covering 75% of the estimated error, followed by
conforming longest-edge closure. While integral control is active, parents
containing MaxTransDiff-marked edges receive red splits so their span halves;
this prevents a chain of green diagonals from deferring the doping stop past
the resource limit. The raw profile with integral control disabled is unchanged.
This is an estimator relative to quadrature,
not a mathematical bound on the unknown exact integral. Analytic Gaussian
convergence and resource-exhaustion tests check its behavior independently.
Source minima still stop MaxTransDiff; independent integral accuracy refinement
may produce smaller cells. The existing 12-round/200000-node budget remains
in force and exhaustion fails closed. The internal target does not replace or
relax the user-approved **1% comparative integral acceptance**.

## Qualification and commercial validation

[The frozen gates](../configs/reference_sde_p6_acceptance.json) retain geometry
and contact tolerance 1e-8 um, node doping relative tolerance 1e-6, own chord
error 1e-5 um and source radii 0.05/0.08 um. Each native/own mesh uses its own
measured chord budget for curved areas. Donor/acceptor whole-domain integrals
must differ by at most 1%. [Frozen cuts](../configs/reference_sde_p6_cuts.json)
must actually contain junctions; their position errors must not exceed the
applicable source directional minimum (or local graded multibox spacing).
Each mesh meets its own authored local-size contract. Different meshes need
not share node numbering or topology. Every repeated run retains its own hash.

```bash
python3 scripts/qualify_sde_package.py --package-root <extracted-package> \
  --output-dir <new-corpus-output> --gmsh-python <gmsh.py> \
  --gmsh-library <libgmsh> --gmsh-executable <gmsh>
python3 scripts/audit_sde_native_corpus.py --native-root <private-native-output> \
  --own-root <new-corpus-output> --run-index 1 --output <private-audit.json>
```

The package qualifier executes all eight frozen sources twice with isolated
Python, checks geometry/fields, cross-format consistency and independent
GRD/DAT imports, and preserves byte equality only for cases that require it.
Commercial native SDE runs use separately recorded `@node@=0` preprocessing
and the source's `+z` context; original and effective hashes are distinct.
The auditor checks native region/material identity, bounds, analytic areas,
contact positions/spans/lengths, interface identities/lengths, point fields,
integrals, nonempty cut roots and local spacing. These are modeling checks.

TDX validates generated GRD/DAT through a TDR and back to GRD/DAT. The returned
regional fields are mapped by coordinates for validation only, including all
stored scalar/species fields and exact mapped element/boundary identities.
Stored net doping is compared directly: subtracting rounded, nearly equal
carrier values would lose the compensated net field. The production GRD
reader's region-local boundary is unchanged. Ogawa TDR produces explicit E402;
the existing native HDF5 TDR reader does not read Ogawa.

```bash
python3 scripts/audit_sde_tdx_readback.py --original-root <corpus-output> \
  --readback-root <private-tdx-readback> --output <private-readback-audit.json>
```

[The CI workflow](../.github/workflows/sde-modeling.yml) independently builds
Windows/UCRT64 and Ubuntu 26.04 programs, assembles and relocates the isolated
package, runs the eight-source corpus twice and required regressions. It
retains package archives and structured summaries. Commercial tools and their
licenses are not dependencies of public CI. Actual workflow/run results must
be recorded separately from this workflow definition.

TDR, raw logs, runtime paths, remote jobs and generated meshes stay in ignored
local evidence. Public evidence contains statistics, source/output hashes and
CI links only. P6 does not qualify new IV/BV physics, NonLocalPath, external
resistor R2/R3 or the PN2D non-obtuse solver profile.

## Registered extensions

| Extension | Follow-up gate before implementation |
| --- | --- |
| SimpleMOS/LDMOS | User priority, frozen original-source inventory and physical contracts |
| Polygon body | Analytic area, contact and boolean ownership fixtures |
| Other boolean modes | Explicit mode semantics and overlap/interface oracle |
| Nonzero decay | Species/placement semantics and independent profile/integral oracle |

These entries are a backlog, not implemented capabilities. 3-D, ACIS and full
Elk Scheme remain outside the current scope.
