# Optional oneMKL PARDISO integration into main

Date: 2026-10-23. Starting main: `179aa4e0`.

## Scope

Extracted the basic `pardiso_mkl` runtime adapter from the uncommitted SLOT
working tree based on `5e3a8cad`. This is a selective integration, not a merge
of the SLOT branch. The shared direct-backend interface, exact pattern/value
caches, error invalidation and default backend preference remain in place.

The adapter loads an explicitly configured absolute runtime path on Windows or
Linux, uses LP64 real nonsymmetric binary64 CSR, and reports vendor errors
without fallback. Calls set and restore the requested MKL thread-local setting.
Ordinary builds require no MKL headers or link-time MKL library. Linux links the
system dynamic-loader library. The LP64 size guard reserves room for one-based
CSR endpoints. No constraint elimination, SLOT population update, physics,
nonlinear acceptance or separate electrothermal backend was imported.

The configuration contract is in [the schema](../config_schema.md). Runtime
paths and dependent-library search paths belong in local environment settings.

## Verification

Local Windows UCRT64 Debug build used installed Eigen 5.0.1 and Boost 1.92.0.
Configuration confirmed HDF5/HighFive state storage, UMFPACK and SPQR enabled.
The real-runtime controls used installed oneMKL 2026.1.

| Check | Result |
| --- | --- |
| Configure and complete Debug build | Passed; final incremental build also passed |
| Full Debug CTest coverage, including targeted rechecks below | 1080 passed, 5 skipped, no unresolved failures out of 1085 |
| PARDISO contract with real MKL, requested 1/2/4 threads | Each: 3 cases, 46 assertions passed |
| Final CTest registration with real MKL | 3/3 passed |
| Poisson executable with `VELA_LINEAR_SOLVER=pardiso_mkl` | 23 cases, 89 assertions passed |
| Gummel executable with that override | 13 cases, 216 assertions passed |
| Newton executable with that override | 148 cases, 3823 assertions passed |
| DC sweep executable with that override | 120 cases, 4415 assertions passed |

The new tests cover missing/relative runtime paths, primary-variable precedence,
the legacy fallback, missing symbols, vendor failures during phases 11/22/33,
and subsequent recovery. Error fixtures are test-only shared libraries, never
linked into the production core. Real MKL tests check nonsymmetric equations,
non-unit affine constraints, multiple RHSs, numerical-factor reuse, forced
numeric invalidation, pattern changes, and recovery after invalid inputs.
Their numerical checks use solution and equation residuals, not implementation
identities or bitwise multithread agreement.

The first sandbox run could not create and clean nested temporary directories.
It was interrupted, its log retained, and the full suite was run with
`TEMP`, `TMP` and `TMPDIR` pointing into the ignored local build directory.

That full run completed all 1085 entries, with two follow-ups required:

- Catch2's absent-MKL skip was initially misclassified because the existing
  Windows PATH-property wrapper corrupted the generated skip metadata. The new
  target now uses direct Catch2 discovery, preserving inherited loader PATH and
  skip metadata; existing targets are unchanged. Without MKL, the real-runtime
  test is explicitly skipped. With MKL, all three tests pass through CTest.
- The existing `Element HFS matrix includes contact, projected and interior
  cell chains` exceeded the run's 300-second cap. Its unchanged high-precision
  finite-difference test passed in 320.39 seconds on a targeted CTest recheck
  with a 1200-second execution cap. No numerical tolerance or model changed.

The final cumulative full-suite result is 1080 passed and five conditional
skips: the optional real-MKL case in the default environment and four native
SDevice end-to-end cases. The real-MKL case was separately exercised successfully.

## Evidence and limits

Local evidence is retained under ignored `build/pardiso-assessment/`:
`configure.log`, `build.log`, `build_final_metadata.log`, `ctest_debug.log`,
`ctest_debug_recheck.log`, `ctest_mkl_final.log`, the three
`mkl_contract_threads_*.log` files and four `mkl_test_*.log` files.
Initial sandbox failures remain in `ctest_debug_initial_sandbox.log`.

Isolated snapshot jobs on T470p and C560 both reached compilation after source
hash verification. SSH connectivity was subsequently lost; their final build
and test results are unavailable and are not counted as passing evidence.
Remote paths, process IDs, manifests and connection details remain in ignored
local records. The remote snapshot predates the local CTest skip-metadata fix.

This qualifies the scoped Windows shared-solver integration and the exercised
regressions. It does not establish complete SLOT BVDS/IdVd convergence, a speed
ranking, or equivalence to Sentaurus's internal PARDISO. Fresh Linux execution
and Python processes with preinitialized NumPy/MKL remain unqualified here.
