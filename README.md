# Vela TCAD

Vela TCAD is a lightweight C++20 prototype for 2-D semiconductor device
simulation. It implements finite-volume Poisson and drift-diffusion solvers,
engineering device sweeps, post-processing diagnostics, and optional pybind11
Python bindings.

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
  model hooks used by the current examples.
- CSV and VTK outputs for regression and visualization.
- Optional Python API for the implemented C++ paths.

Current example coverage:

- PN diode IV/CV/BV smoke decks.
- All-silicon NMOS/PMOS drift-diffusion IV, Id-Vg, and CV decks.
- Mixed Si/SiO2 NMOS/PMOS MOS DD prototypes with multi-terminal CV, BV
  diagnostics, and surface-mobility smoke coverage where present.
- Schottky diode IV prototype.
- LDMOS and IGBT power-device trend fixtures using `unit_scaling`.
- Reference TCAD CSV fixtures for PN, MOS, LDMOS, and IGBT trend comparison.

Input unit modes:

- No `scaling` field keeps the legacy SI behavior used by older decks.
- `"scaling": {"mode": "unit_scaling"}` enables common external TCAD input
  units at the schema boundary.

Prototype boundaries:

- Power-device decks are engineering trend validations, not calibrated
  LDMOS/IGBT models.
- BV output is a diagnostic sweep with maximum edge-field and convergence
  indicators, not a calibrated avalanche breakdown prediction.
- Quasi-static CV is finite-difference terminal-charge extraction, not an AC
  small-signal matrix solve.
- Schottky support is a prototype barrier-style path; Newton sweeps reject
  Schottky contacts until a future implementation handles that model.

## Documentation Map

Start here:

- [docs/README.md](docs/README.md): documentation index and reading paths.
- [docs/architecture.md](docs/architecture.md): implementation and module map.
- [docs/config_schema.md](docs/config_schema.md): JSON configuration reference.
- [docs/examples.md](docs/examples.md): supported example and regression matrix.
- [tests/regression/README.md](tests/regression/README.md): regression runner
  behavior and assertion fields.
- [reference_tcad/README.md](reference_tcad/README.md): neutral reference CSV
  fixture workflow.

Historical planning and handoff notes live under `docs/` and are retained as
archives. Current behavior should be verified against code, CMake targets, and
the schema/example documents above.

## Build

Prerequisites:

- CMake 3.20 or newer
- C++20 compiler
- Eigen3
- nlohmann/json
- Catch2 v3
- Python 3 interpreter for CTest regression orchestration

Ubuntu/Debian:

```bash
sudo apt-get update
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y \
  build-essential \
  cmake \
  ninja-build \
  libeigen3-dev \
  nlohmann-json3-dev \
  catch2
```

Configure and build:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
```

If Ninja is unavailable, omit `-G Ninja`.

### Windows / MSYS2 UCRT64

The Windows development environment uses MSYS2 UCRT64, typically at
`D:\msys64`. Keep UCRT64 tools first on `PATH`:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
```

Install packages from a UCRT64 shell:

```bash
pacman -Syu
pacman -S --needed \
  mingw-w64-ucrt-x86_64-toolchain \
  mingw-w64-ucrt-x86_64-cmake \
  mingw-w64-ucrt-x86_64-ninja \
  mingw-w64-ucrt-x86_64-eigen3 \
  mingw-w64-ucrt-x86_64-nlohmann-json \
  mingw-w64-ucrt-x86_64-catch \
  mingw-w64-ucrt-x86_64-python \
  mingw-w64-ucrt-x86_64-gdb
```

Use `python` instead of `python3` if the UCRT64 shell does not provide a
`python3` alias.

## Test

Run the full suite:

```bash
ctest --test-dir build --output-on-failure
```

Useful focused groups:

```bash
ctest --test-dir build --output-on-failure -R poisson
ctest --test-dir build --output-on-failure -R "dd|newton|dc_sweep"
ctest --test-dir build --output-on-failure -R regression
ctest --test-dir build --output-on-failure -R reference_tcad_regression
```

The `regression` CTest target runs `scripts/run_regression.py` against the
engineering examples with `vela_example_runner`. The
`reference_tcad_regression` target verifies the neutral CSV conversion and
comparison tools.

## Run Examples

After building:

```bash
build/vela_example_runner --config examples/pn_diode/simulation_iv.json
build/vela_example_runner --config examples/pn_diode/simulation_cv.json
build/vela_example_runner --config examples/pn_diode/simulation_bv.json
build/vela_example_runner --config examples/pn_diode/newton_simulation.json
```

On Windows, use `build\vela_example_runner.exe`.

Run the complete engineering example suite:

```bash
python scripts/run_regression.py --runner build/vela_example_runner
```

The script copies examples to `build/regression_output/`, runs each configured
deck, verifies CSV/VTK files, checks finite outputs and trend assertions, and
writes `build/regression_output/regression_summary.json`.

## Optional Python API

The Python extension is disabled by default. To enable it, install Python
development headers and pybind11, then configure with `VELA_ENABLE_PYTHON=ON`:

```bash
# Ubuntu/Debian
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y python3-dev pybind11-dev

cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DVELA_ENABLE_PYTHON=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R python_api
```

On Windows/MSYS2:

```bash
pacman -S --needed mingw-w64-ucrt-x86_64-pybind11
```

The generated package is placed under `build/python/<config>/vela`, for example
`build/python/Debug/vela`. The CTest `python_api` target sets `PYTHONPATH`
automatically.

Python examples:

```python
import vela

mesh = vela.load_mesh("examples/pn_diode/mesh.json")
potential = vela.run_poisson("examples/pn_poisson_2d.json")
iv_points = vela.run_iv_curve("examples/pn_diode/simulation_iv.json")
cv_points = vela.run_cv_curve("examples/pn_diode/simulation_cv.json")
bv_points = vela.run_bv_curve("examples/pn_diode/simulation_bv.json")
```

The Python API is intentionally thin and only documents behavior exercised by
the C++ core and tests.

## Debug

Use a Debug build and the debugger from the same toolchain as the build.

Windows/MSYS2 UCRT64:

```bash
gdb --args build/vela_example_runner.exe --config examples/pn_diode/simulation_iv.json
gdb --args build/test_poisson.exe
```

For VS Code or another MI-compatible debugger, point `miDebuggerPath` to:

```text
D:\msys64\ucrt64\bin\gdb.exe
```

Avoid mixing MSYS2 UCRT64, CLANG64/MINGW64, and Visual Studio build outputs in
the same build directory.
