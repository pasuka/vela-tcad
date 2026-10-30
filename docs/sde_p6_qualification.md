# SDE offline modeling qualification (P6)

This is an independent modeling-tools distribution. Its products are
`device.grd`, `device.dat` and `device.msh`. The existing SDevice package and
reserved `sd-sim -d`, `-n`, `--max_threads` interface remain unchanged.
P6 qualification is in progress; the final evidence record must identify the
actual source revision, platform runs and CI conclusions before closure.

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
conforming longest-edge closure. This is an estimator relative to quadrature,
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
