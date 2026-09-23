# Vela TCAD

Vela TCAD is a lightweight C++20 prototype for 2-D semiconductor device
simulation. It implements finite-volume Poisson and drift-diffusion solvers,
engineering device sweeps, post-processing diagnostics, and optional pybind11
Python bindings plus optional HDF5-based import tooling.

The project is intentionally transparent and small. It is useful for solver
development, regression experiments, and device-trend smoke tests. It is not a
calibrated commercial TCAD replacement.

## What Is Implemented

Core solver capabilities:

- 2-D triangular device meshes with regions, contacts, edge topology, and box
  geometry.
- Material database for Si and SiO2, with optional JSON material overrides.
- Region doping, fixed charge, interface sheet/fixed/trap charge, and explicit
  Poisson boundary segments.
- Poisson assembly and solve, including a scaled assembly path for
  `scaling.mode = "unit_scaling"`.
- Drift-diffusion assembly with Scharfetter-Gummel fluxes.
- Gummel and coupled Newton nonlinear solves.
- Adaptive DC sweeps for IV, quasi-static CV, and reverse-bias/BV diagnostics.
- Contact-current, terminal-charge, stored-charge, and electric-field
  diagnostics.
- Mobility, recombination, impact-ionization, and Slotboom bandgap-narrowing
  model hooks used by the validated reference cases and focused tests.
- CSV and VTK outputs for regression and visualization.
- Optional Python API for the implemented C++ paths.
- Optional HDF5 inventory and neutral export tooling.

Current device-level validation coverage is maintained under `reference_tcad/`.
It includes Sentaurus-backed PN, NMOS, Schottky, breakdown-method,
TransportModels, and Genius NPN BJT fixtures with explicit acceptance
boundaries.

Input unit modes:

- No `scaling` field keeps the legacy SI behavior used by older decks.
- `"scaling": {"mode": "unit_scaling"}` enables common external TCAD input
  units at the schema boundary.

Implementation boundaries:

- BV output is a diagnostic sweep with maximum edge-field and convergence
  indicators, not a calibrated avalanche breakdown prediction.
- Quasi-static CV is finite-difference terminal-charge extraction, not an AC
  small-signal matrix solve.
- Schottky support is a prototype barrier-style path; Newton sweeps reject
  Schottky contacts until a future implementation handles that model.

Current PN2D BV template policy:

- Template version 3 defaults to the qualified atomic
  `element_edge_sg_gss_laux` profile: element-edge SG/GSS-Laux current,
  element-vertex box source mapping, Bernoulli midpoint density,
  mixed-Voronoi node volumes, and non-obtuse mesh enforcement.
- `legacy_cell_reconstructed` is the complete rollback profile and restores
  the former cell-reconstructed/barycentric configuration.
- This qualification is limited to the validated non-obtuse PN2D BV M0/M2
  Tri3 mesh family. Global C++ defaults and the PN2D IV template are unchanged.

## Documentation Map

Start here:

- [docs/README.md](docs/README.md): documentation index and reading paths.
- [docs/architecture.md](docs/architecture.md): implementation and module map.
- [docs/config_schema.md](docs/config_schema.md): JSON configuration reference.
- [docs/validation/pn2d_bv_validation.md](docs/validation/pn2d_bv_validation.md):
  current PN2D BV validation decision, scope, and evidence map.
- [tests/regression/README.md](tests/regression/README.md): regression runner
  behavior and assertion fields.
- [reference_tcad/README.md](reference_tcad/README.md): neutral reference CSV
  fixture workflow.

Dated validation reports preserve point-in-time evidence. They are not the
source of truth for current defaults when a later decision document supersedes
them.

## Build

On Windows, this repository is developed primarily in MSYS2 UCRT64. If a tool or agent needs to build, test, or debug on Windows, it should assume `D:\msys64\ucrt64` is the default toolchain unless the user says otherwise.

Prerequisites:

- CMake 3.21 or newer
- C++20 compiler
- Boost headers (Boost.Multiprecision)
- Eigen3
- nlohmann/json
- Catch2 v3
- Python 3 interpreter for CTest regression orchestration
- HDF5 development package and HighFive 3 for restart storage. CMake uses an
  installed HighFive package or fetches the pinned v3.3.0 source. TDR import is
  separately controlled by `VELA_ENABLE_HDF5`; disabling import does not disable
  the required state archive library.
- Python `numpy` and `h5py` for state preparation and regression tests. Use the
  same Python interpreter selected by CMake; its h5py build and loaded HDF5
  runtime must be compatible.

Ubuntu/Debian:

```bash
sudo apt-get update
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y \
  build-essential \
  cmake \
  ninja-build \
  libboost-dev \
  libeigen3-dev \
  nlohmann-json3-dev \
  libspdlog-dev \
  catch2 \
  python3 python3-numpy python3-h5py \
  libhdf5-dev
```

Use these installation steps for initial setup or confirmed missing dependencies.
Verify that the distribution provides Catch2 v3; CMake requires version 3.
The Python interpreter is required even when the optional Python API is disabled.
If running as root without `sudo`, omit `sudo` from these commands.

Configure and build:

```bash
cmake --preset windows-ucrt64-debug
cmake --build --preset windows-ucrt64-debug
```

These presets require Ninja. For a Linux build without the Windows presets,
see the out-of-tree build commands in [AGENTS.md](AGENTS.md#environment-and-build).

The repository also ships `CMakePresets.json`, so CMake Tools and command-line
workflows can share the same Windows UCRT64 configuration.

Preset summary (from `CMakePresets.json`):

| Preset | Kind | Binary dir | Python | TDR import | Typical command |
| --- | --- | --- | --- | --- | --- |
| `windows-ucrt64-debug` | Configure / Build / Test | `build/` | OFF | ON | `cmake --preset windows-ucrt64-debug` |
| `windows-ucrt64-debug-python` | Configure / Build / Test | `build-python/` | ON | ON (inherited) | `cmake --preset windows-ucrt64-debug-python` |
| `windows-ucrt64-poisson` | Test | `build/` | OFF | ON | `ctest --preset windows-ucrt64-poisson` |

HDF5 and HighFive are required for production restart storage in every preset;
configuration fails if the state-storage dependencies are unavailable. The
optional TDR import tools can be disabled with `VELA_ENABLE_HDF5=OFF` without
disabling HDF5 state storage.

### Windows / MSYS2 UCRT64

The Windows development environment uses MSYS2 UCRT64, typically at
`D:\msys64`. The Windows presets prepend the required UCRT64 directories to
their configure, build, and test environments, so they also work from a plain
PowerShell session:

```powershell
D:\msys64\ucrt64\bin\cmake.exe --preset windows-ucrt64-debug
D:\msys64\ucrt64\bin\cmake.exe --build --preset windows-ucrt64-debug
```

Configured MinGW build trees also use the repository compiler and linker
launcher in `scripts/with_ucrt64_path.cmd`. Consequently, a direct command such
as `D:\msys64\ucrt64\bin\cmake.exe --build build` can start `cc1`/`cc1plus`,
`collect2`, and `ld` without manually editing `PATH`. This setting is
repository-local and does not change the Windows system or user environment.
Set the UCRT64 `PATH` manually only when invoking raw MSYS2 tools outside these
CMake workflows.

For initial setup or environment repair, install packages from a UCRT64 shell.
`pacman -Syu` is a system update step, not a prerequisite for each repository task:

```bash
pacman -Syu
pacman -S --needed \
  mingw-w64-ucrt-x86_64-toolchain \
  mingw-w64-ucrt-x86_64-cmake \
  mingw-w64-ucrt-x86_64-ninja \
  mingw-w64-ucrt-x86_64-boost \
  mingw-w64-ucrt-x86_64-eigen3 \
  mingw-w64-ucrt-x86_64-nlohmann-json \
  mingw-w64-ucrt-x86_64-spdlog \
  mingw-w64-ucrt-x86_64-hdf5 \
  mingw-w64-ucrt-x86_64-catch \
  mingw-w64-ucrt-x86_64-python \
  mingw-w64-ucrt-x86_64-python-numpy \
  mingw-w64-ucrt-x86_64-python-h5py \
  mingw-w64-ucrt-x86_64-gdb
```

Use `python` instead of `python3` if the UCRT64 shell does not provide a
`python3` alias.

## Test

Run the full suite:

```bash
ctest --preset windows-ucrt64-debug
```

Useful focused groups:

```bash
ctest --preset windows-ucrt64-poisson
ctest --test-dir build --output-on-failure -R "dd|newton|dc_sweep"
ctest --test-dir build --output-on-failure -R reference_tcad_regression
ctest --test-dir build --output-on-failure -R import
```

Named CTest targets exposed by `CMakeLists.txt`:

| Target | Always present | Notes |
| --- | --- | --- |
| `poisson` | Yes | Poisson-focused Catch2 entry |
| `dd` | Yes | DD/Gummel-focused Catch2 entry |
| `device_stability` | Yes | Device stability checks |
| `dc_sweep` | Yes | DC sweep behavior |
| `electric_field_diagnostics` | Yes | Electric-field diagnostics |
| `boundary` | Yes | Boundary-condition checks |
| `schottky` | Yes | Schottky contact prototype checks |
| `interface` | Yes | Interface charge checks |
| `reference_tcad_regression` | Yes | Reference TCAD conversion/comparison tools |
| `import_tools` | Yes | Python-level import-tool checks |
| `python_api` | Conditional | Present when `VELA_ENABLE_PYTHON=ON` |
| `import_tdr` | Conditional | Present when HDF5 target is found |
| `import_sample_integration` | Conditional | Present when HDF5 target is found |

The `reference_tcad_regression` target verifies the neutral CSV conversion and
comparison tools used by checked-in cross-TCAD fixtures.

For HDF5 import workflows and tool usage, see the documentation index and the
reference fixture workflow docs under `docs/` and `reference_tcad/`.

## Run Reference Cases

Validated device inputs, comparison reports, and reproduction entry points are
listed in [reference_tcad/README.md](reference_tcad/README.md).  Each fixture
documents the required run order because multi-stage continuation cases cannot
in general be launched from a single standalone deck.

## Optional Python API

The Python extension is disabled by default. To enable it, install Python
development headers and pybind11, then configure with `VELA_ENABLE_PYTHON=ON`:

```bash
# Ubuntu/Debian
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y python3-dev pybind11-dev

cmake --preset windows-ucrt64-debug-python
cmake --build --preset windows-ucrt64-debug-python
ctest --preset windows-ucrt64-debug-python -R python_api
```

The `windows-ucrt64-debug-python` preset still honors `VELA_ENABLE_HDF5=ON`
from the base preset, so the optional import tooling remains available when
HDF5 is installed.

On Windows/MSYS2:

```bash
pacman -S --needed mingw-w64-ucrt-x86_64-pybind11
```

The generated package is placed under `build-python/python/<config>/vela`, for example
`build-python/python/Debug/vela`. The CTest `python_api` target sets `PYTHONPATH`
automatically.

The Python API is intentionally thin and only documents behavior exercised by
the C++ core and tests.

## Debug

Use a Debug build and the debugger from the same toolchain as the build.

Windows/MSYS2 UCRT64:

```bash
gdb --args build/test_poisson.exe
```

For VS Code or another MI-compatible debugger, point `miDebuggerPath` to:

```text
D:\msys64\ucrt64\bin\gdb.exe
```

Avoid mixing MSYS2 UCRT64, CLANG64/MINGW64, and Visual Studio build outputs in
the same build directory.
