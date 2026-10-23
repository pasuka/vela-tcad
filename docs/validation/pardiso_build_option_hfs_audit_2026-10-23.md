# PARDISO build option and HFS test cost audit

Date: 2026-10-23. This follows the
[initial cross-host validation](pardiso_cross_host_validation_2026-10-23.md).

## Build behavior

`VELA_ENABLE_MKL_PARDISO` is a CMake option with default `OFF`. A disabled
build excludes `MklPardisoBackend.cpp`, the two mock MKL shared libraries and
the adapter's explicit dynamic-loader link dependency. The solver still
recognizes an explicit `pardiso_mkl` request and rejects it with the required
configure option, rather than silently changing backend. Its disabled-build
test runs without MKL and also checks that an explicit SparseLU solve works.

GitHub CI explicitly configures `-DVELA_ENABLE_MKL_PARDISO=OFF` and checks the
disabled configure message. The workflow needs no MKL installation. An enabled
Windows or Linux build retains the runtime-loading implementation, error-path
fixtures and real-runtime tests; the MKL SDK is not a compilation requirement.
Enabling the option on another platform fails at configure time.

For example, on the default Windows toolchain:

```powershell
cmake --preset windows-ucrt64-debug -DVELA_ENABLE_MKL_PARDISO=ON
cmake --build --preset windows-ucrt64-debug
```

Runtime selection additionally needs `VELA_LINEAR_SOLVER=pardiso_mkl`, an
explicit installed runtime path and its dependency search path; see the
[configuration reference](../config_schema.md). CMake caches retain explicit
selections; use `-DVELA_ENABLE_MKL_PARDISO=OFF` to turn a previous ON build off.
The current local Debug build is configured OFF.

## Why the HFS test takes minutes

The slow test is `Element HFS matrix includes contact, projected and interior
cell chains` in [test_element_box_transport.cpp](../../tests/test_element_box_transport.cpp).
It audits a complete Jacobian; it is not a device sweep or a linear solve.

- Its mesh has 20 nodes and 24 triangular cells (18 Si, 6 oxide), with 60
  potential/quasi-Fermi unknowns. Coordinates have a 10 nm grid spacing;
  doping is supplied in SI number density and temperature is 300 K.
- Two carrier-population states, 60 Jacobian columns and three perturbations
  (`1e-7`, `1e-12`, `1e-25` in packed coordinates) give 360 central differences.
  Each [symmetricDifference](../../src/equation/SplitDDRuntime.cpp) evaluates
  both perturbed residuals, giving at least 720 full residual evaluations,
  plus state and terminal-current consistency checks.
- [SplitDDState.h](../../include/vela/numerics/SplitDDState.h) selects
  `boost::multiprecision::cpp_bin_float_100`: 100 decimal digits, not 100 bits.
  This software arithmetic preserves sub-ULP responses when O(1 V) references
  coexist with O(1e-25 V) increments. Simply replacing it with double would
  change what the test qualifies.
- [SplitDDOperator.h](../../include/vela/equation/SplitDDOperator.h) recomputes
  full-mesh densities, PhuMob/Lombardi/Canali mobility, SRH and SG transport.
  These involve high-precision `exp`, `expm1`, `log` and noninteger powers.
  The runtime caches only the last exact high/low state; successive positive
  and negative column perturbations cannot reuse that evaluation. PhuMob's
  screening-minimum bisections are constructor work, not repeated per residual.
- The tested Debug configuration uses `-g` without an optimization flag.
  The column/step loops run serially inside this test, so CTest parallelism
  does not split its 360 differences across workers.

This establishes the work count and expensive arithmetic path from source;
it does not assign sampled per-function CPU percentages. The earlier Debug
runs took about 305 s on T470p and 485 s on C560. Those are historical test
times, not a controlled hardware or backend performance comparison.
PARDISO is not called by this test, so enabling MKL cannot accelerate its
reference residual differences. No HFS implementation, precision, test coverage
or numerical tolerance was changed for this task. A future optimization should
first profile a matching Release build and consider invariant precomputation
while preserving the existing numerical gates.

## Verification

Both isolated remote source snapshots were verified against 5905 file hashes
before and after testing. Their source base was `26bb5681`, plus the option
patch archived with SHA-256
`904f5a7757f1e13e008a51d97822b4a0120af442a6d63c1d424ef0642d60a8ac`.
The source/configuration/test changes are identical to the final local files;
this report and its documentation-index entry were added after the runs.
The existing environments and task-local regression dependencies from the
[cross-host record](pardiso_cross_host_validation_2026-10-23.md) were reused.

| Check | T470p, Windows/UCRT64 | C560, Linux |
| --- | --- | --- |
| Toolchain | GCC 16.2.0, CMake 4.4.4 | GCC 15.2.0, CMake 4.2.3 |
| Build type / Python API | Debug / OFF | Debug / OFF |
| HDF5 state storage / UMFPACK / SPQR | Enabled | Enabled |
| Default configure, without the new flag | OFF | OFF |
| OFF build graph excludes adapter and mock libraries | Passed | Passed |
| Complete OFF Debug build | Passed | Passed |
| Complete OFF CTest | 1079 passed, 4 conditional skips, 0 failed | 1079 passed, 4 conditional skips, 0 failed |
| Complete OFF CTest wall time | 395.93 s | 702.65 s |
| HFS contact/projected/interior test wall time in that run | 359.89 s, passed | 702.56 s, passed |
| ON configure and complete build | Passed | Passed |
| ON PARDISO CTest without a selected runtime | 2 passed, 1 skipped | 2 passed, 1 skipped |
| ON PARDISO CTest with real MKL | 3/3 passed | 3/3 passed |
| Requested MKL threads 1 / 2 / 4 | Each: 3 cases, 46 assertions passed | Each: 3 cases, 46 assertions passed |
| Real-MKL Poisson / Gummel / Newton / DC sweep | 304 cases, 8543 assertions passed | 304 cases, 8543 assertions passed |
| Source hashes unchanged afterwards | Passed | Passed |

Each OFF suite was a single complete 1083-entry run, with CTest parallelism 2
and a 1200 s per-test execution cap. The four conditional skips are the native
SDevice E2E tests; they are not failures or missing PARDISO coverage. The HFS
times above include competition from the other test worker and are not a
controlled performance benchmark. No new failure required a recheck.

The real-runtime checks selected installed MKL 2026.1 on T470p and 2025.2 on
C560, with the shared-solver override `VELA_LINEAR_SOLVER=pardiso_mkl` and
one requested thread for the physics executables. Those executables comprise
Poisson (23 cases, 89 assertions), Gummel (13, 216), Newton (148, 3823) and
DC sweep (120, 4415). Existing cases, meshes and numerical tolerances were
retained. The optional Python API and preinitialized NumPy/MKL processes are
outside this validation scope. The GitHub-hosted workflow itself was not run;
its OFF configuration is covered by these local/remote builds and tests.

The local Windows preset also completed its default OFF build. Its disabled
backend case passed 4 assertions and was discovered/passed by CTest; the
linear-solver executable passed 13 cases / 365 assertions.

## Evidence

Private job paths, source manifests, runtime binary hashes and full logs are
retained under ignored `build/pardiso-option/`. Both downloaded evidence
archives were SHA-256 checked:

- T470p `t470p_option_evidence.tgz`:
  `5ebc8f5b5fb390a5a719c4a8b0e4565602c72e5e5f8de0c25175f18772d00e1d`.
- C560 `c560_option_evidence.tgz`:
  `8fac5a44f4d79ac67d7d801bc2742d926a17e04bef2a70ea8bcd2abfd1d264d1`.

The base backend remains its own earlier commit (`d8566e12`); this build option
is a separate change. Constraint-elimination variants were not introduced.
