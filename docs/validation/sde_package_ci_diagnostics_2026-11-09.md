# SDE package CI failure diagnostics and Codespaces verification (2026-11-09)

The failed [GitHub Actions run](https://github.com/pasuka/vela-tcad/actions/runs/42556107841)
reported a BVMethods Gmsh preprocessing failure, but did not retain its inner
`fail_report.json`. The original meshing cause therefore remains unidentified.
The changes here repair the missing diagnostic/evidence path; they do not claim
a correction to an unobserved meshing defect.

## Changes

Commit `a52df73d` makes the package qualification wrapper retain a failed
`qualification.json` and surface the child's structured stage/message.
Malformed, missing or unrelated reports fall back to the subprocess error.
A failed child remains a failure, and subsequent quality checks do not run.
The CI artifact now retains nested `package-corpus/**/fail_report.json` files.
New regression cases cover this failure propagation and stale-success receipt
replacement. Mesh algorithms, authored geometry, refinement settings, resource
budgets and acceptance tolerances are unchanged.

## Exact inputs and verification

- Original package source: `98f577e9362ac24eeae02b0a3ce1b8d384c21b9b`.
- Original package SHA-256: `ec53800ecdaa76f0600f91d079543795422433be6a3d98463a94f0829aeb0ad6`.
- Diagnostic wrapper/source checkout: `a52df73db9e268cf51add2a134991903c22073d5`, clean.
- The original archive and its modeling payload were not patched or repackaged.
- Local Ubuntu 26.04: one initial BVMethods run and eight additional independent
  BVMethods runs passed. The latter produced eight distinct mesh hashes and all
  finished refinement at round 11 under the unchanged 12-round / 200000-node cap.
- Linux package regression module: 13/13 passed, including actual extracted
  package relocation, corruption rejection and dependency boundaries; no skips.
- GitHub Codespaces Ubuntu 26.04.1, host Python 3.14.4, sealed package Python
  3.12.12 and hash-locked external Gmsh 4.15.2: 13/13 package regressions passed;
  all eight original sources passed twice, including the independent GRD/DAT
  import and each source's existing quality contract.

| Source | Full-device nodes, two runs | Triangles, two runs |
| --- | --- | --- |
| pn2d | 1482 / 1482 | 2800 / 2800 |
| schottky | 319 / 319 | 560 / 560 |
| gate_reference | 164 / 163 | 218 / 216 |
| bjt | 22469 / 22614 | 44673 / 44965 |
| bjt_local | 43545 / 43253 | 86783 / 86197 |
| singledevice | 61035 / 65291 | 121348 / 129834 |
| transportmodels | 42883 / 46162 | 85126 / 91682 |
| bvmethods | 51103 / 42571 | 101600 / 84348 |

The Codespaces corpus process exited 0 after 606.633 seconds.
Its source remained clean and no owned modeling processes remained at collection.
The evidence archive SHA-256 is `3c68b31965a38438e639886ebaca627a49eb7529f09fd577d9cd3b164bea0289`; all files listed in its
manifest were independently rehashed after download. Raw job records, private
connection details and generated meshes remain outside version control.

## Environment preparation and limits

The repository's original Codespaces container failed its unrelated C++ configure
step because installed Boost 1.83 did not satisfy the declared minimum 1.84.
For this investigation only, the new Codespace was rebuilt with the official
Ubuntu 26.04 base and the SDE prerequisites. Its first modeling attempt correctly
failed Gmsh runtime validation because graphics shared libraries were missing;
these were identified using `ldd` and installed before the fresh corpus run.
Neither environment failure reproduces the original CI failure: that CI had
already completed the first seven source cases. Temporary Codespace configuration
was not committed as a change to the project devcontainer.

Successful reruns do not identify or eliminate an intermittent meshing defect.
In particular, the observed refinement round counts do not establish that the
original failure exhausted its iteration budget. No retry, tolerance relaxation
or budget increase was introduced to hide a numerical failure. The next actual
CI failure can now retain the report needed for a targeted repair.

This is an archived-package verification, not qualification of a newly built
modeling archive, a rerun of C++ solver tests, or a new successful GitHub Actions
run. No GitHub branch was pushed for this work.
