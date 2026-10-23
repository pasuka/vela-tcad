# PARDISO main integration: T470p and C560 follow-up

Date: 2026-10-23. Tested source revision: `d8566e12`.

Both remote follow-up jobs completed successfully. This record supplies the
previously unavailable remote results from the
[main integration report](pardiso_main_integration_2026-10-23.md).
No production code, physical model, numerical tolerance, backend preference or
nonlinear acceptance rule changed during this follow-up.

## Source and environment identity

Each isolated remote source tree was checked against 5904 file hashes before
testing and checked again afterwards; both final checks passed. The only source
updates from the earlier remote snapshot were the committed CMake test-discovery
fix and the integration report. Existing Debug build directories were reused
with their original toolchains. Remote checkouts belonging to other work were
not changed.

| Environment | T470p | C560 |
| --- | --- | --- |
| Platform | Windows 10 x64, MSYS2 UCRT64 | Linux x86_64, glibc 2.43 |
| Compiler | GCC 16.2.0 | GCC 15.2.0 |
| CMake | 4.4.4 | 4.2.3 |
| Eigen | 5.0.1 | 5.0.1 |
| Boost headers | 1.92.0 | 1.90.0 |
| Selected MKL runtime installation | 2026.1 | 2025.2 |
| HDF5 state storage / UMFPACK / SPQR | Enabled | Enabled |

The MKL runtime binary hashes and complete source manifest are recorded in the
local evidence. Python API builds and preinitialized NumPy/MKL processes are
not covered by these executable-level checks.

## Results

| Check | T470p | C560 |
| --- | --- | --- |
| Initial complete Debug build | Passed | Passed |
| Final-source configure and incremental build | Passed | Passed |
| Full Debug CTest coverage, including failed-test rechecks | 1080 passed, 5 conditional skips, no unresolved failures | 1080 passed, 5 conditional skips, no unresolved failures |
| Final PARDISO CTest registration without MKL configured | 2 passed, 1 skipped | 2 passed, 1 skipped |
| Final PARDISO CTest registration with real MKL | 3/3 passed | 3/3 passed |
| Requested MKL threads 1 / 2 / 4 | Each: 3 cases, 46 assertions passed | Each: 3 cases, 46 assertions passed |
| PARDISO Poisson / Gummel / Newton / DC sweep executables | 304 cases, 8543 assertions passed | 304 cases, 8543 assertions passed |
| Source hashes unchanged after tests | Passed | Passed |

The four numerical executables comprise Poisson (23 cases, 89 assertions),
Gummel (13, 216), Newton (148, 3823) and DC sweep (120, 4415). They used the
explicit shared-solver override `VELA_LINEAR_SOLVER=pardiso_mkl`, requested one
MKL thread and the existing test configurations and tolerances. The separate
thread controls use equation/solution tolerances, not cross-thread bit equality.

The full-suite totals are cumulative: the original complete 1085-entry CTest
run plus targeted rechecks of its failures after fixing the environment.
Already passing tests were not rerun without cause. Five default-environment
skips are the optional real-MKL case and four native SDevice end-to-end cases;
the real-MKL case was subsequently exercised successfully on each host.

## Resolved execution problems

- DHCP addresses had changed. Connections were reestablished using saved SSH
  identities and strict host-key checks. Addresses and connection records remain
  outside version control.
- The unchanged HFS contact/projected/interior Jacobian test exceeded the
  original 300-second execution cap. It passed with a 1200-second recheck cap,
  taking approximately 305 seconds on T470p and 485 seconds on C560.
  No physical or numerical threshold was relaxed.
- C560 lacked `patchelf`. The distro package was downloaded and extracted into
  this task's tools directory; the packaging test then passed. No privileged
  system installation or upgrade was needed.
- T470p lacked `jsonschema`, and its package-mirror download did not complete.
  The existing local UCRT64 Python 3.14 packages and dependencies were copied
  into the task's private Python path, including package metadata. Import and
  schema validation passed before the SDE regression passed. The failed download
  left no running package manager or stale package lock.
- The Windows concurrency test required its existing `source/build-release`
  temporary-directory parent even for this Debug build. Creating that empty
  parent resolved the failure; the actual build remained in the Debug directory.

## Evidence and qualification boundary

Initial full-run logs, final logs, source manifests, runtime hashes and job
records are retained under ignored `build/pardiso-assessment/`. The final
downloaded archives were SHA-256 checked:

- T470p `t470p_followup_evidence.tgz`:
  `ae05a7c3f0d672fc123e37cc9eed71ff684a6418da65d2d05f920d4f934a0b3e`.
- C560 `c560_followup_evidence.tgz`:
  `00251064c0fd70ddb64513aa1d5964dbc2e3dead0cb32f7fee4f4966ddf44065`.

The remote Windows and Linux shared-solver qualification gap in the earlier
report is closed for these environments and tests. This does not qualify a
complete SLOT BVDS/IdVd trajectory, constraint-elimination variants, a performance
ranking, or equivalence to Sentaurus's internal solver.
