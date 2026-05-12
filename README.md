# Vela TCAD

A lightweight **2-D semiconductor device drift-diffusion and Poisson solver**
based on the Finite-Volume Method (FVM/Box integration) and
Scharfetter-Gummel flux discretization.

---

## Project Goal

Vela provides an open, extensible C++20 prototype for simulating planar
semiconductor devices such as p-n diodes, MOS capacitors, and MOSFET
cross-sections. The core solver is written in C++ and now also exposes a small
optional Python API through pybind11.

---

## Current Stage - What Is Implemented

| Component | Status |
|-----------|--------|
| `Types.h` - `Real`, `Index`, `Point2`, `CellType` | done |
| `PhysicalConstants.h` - q, kb, eps0, h, m0, T0, Vt_300 | done |
| `ScalingSystem` - Debye scaling, scale/unscale helpers | done |
| `DeviceMesh` - nodes, cells, edges, regions, contacts | done |
| `JsonMeshReader` - reads JSON meshes | done |
| `VTKWriter` and `CSVWriter` - simulation output files | done |
| `MaterialDatabase` - Si and SiO2 built in | done |
| `DopingModel` - per-node Nd/Na from region JSON spec | done |
| `PoissonAssembler` - FVM/Box electrostatic Poisson | done |
| `ScharfetterGummel` and drift-diffusion assembly | done |
| `LinearSolver`, `GummelSolver`, and Newton helpers | done |
| `PoissonSimulation` - end-to-end Poisson driver | done |
| `DCSweep` - bias sweep with CSV/VTK outputs | done |
| Optional pybind11 Python API | done |
| Unit and regression tests | done |

---

## Build Instructions

**Prerequisites:** CMake >= 3.20, a C++20 compiler, Eigen3, nlohmann/json,
Catch2 v3, and Python 3 for test orchestration.

```bash
# Ubuntu/Debian
sudo apt-get install build-essential cmake ninja-build libeigen3-dev \
  nlohmann-json3-dev catch2 python3
```

Configure and build with Ninja:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
```

If Ninja is unavailable, omit `-G Ninja` and use the default CMake generator.

### Windows with MSYS2 UCRT64

The Windows development environment is based on **MSYS2 UCRT64**, typically
installed at `D:\msys64`. Use the **UCRT64** shell so that CMake, Ninja, GCC,
GDB, and Python all come from the same UCRT64 toolchain.

Install the required packages:

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

From the UCRT64 shell:

```bash
cd /d/code-repo/vela-tcad
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Use `python` instead of `python3` for local helper scripts if the UCRT64 shell
does not provide a `python3` alias.

From PowerShell, put the MSYS2 UCRT64 tools first on `PATH` before configuring:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

### Optional Python API

The Python extension is disabled by default. To build it, install pybind11 and
Python development headers, then enable `VELA_ENABLE_PYTHON`:

```bash
# Ubuntu/Debian
sudo apt-get install pybind11-dev python3-dev

cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DVELA_ENABLE_PYTHON=ON
cmake --build build --parallel
```

CMake places the importable package under `build/python/<config>/vela`, such as
`build/python/Debug/vela` for a Debug build. The CTest target sets `PYTHONPATH`
automatically. For ad hoc runs, point `PYTHONPATH` at the generated
`build/python/<config>` directory.

On Windows/MSYS2, install the optional binding dependency with:

```bash
pacman -S --needed mingw-w64-ucrt-x86_64-pybind11
```

---

## Running Tests

```bash
ctest --test-dir build --output-on-failure
```

The standard build runs the C++ Catch2 tests and the engineering regression
suite. When configured with `-DVELA_ENABLE_PYTHON=ON`, CTest also registers the
`python_api` test.

Run selected test groups:

```bash
ctest --test-dir build --output-on-failure -R poisson
ctest --test-dir build --output-on-failure -R regression
ctest --test-dir build --output-on-failure -R python_api
```

---

## Debugging

Debug builds include symbols through `-DCMAKE_BUILD_TYPE=Debug`. On
Windows/MSYS2 UCRT64, use the matching UCRT64 GDB:

```bash
gdb --args build/vela_example_runner.exe --config examples/pn_diode/simulation.json
```

Inside GDB, common commands are:

```text
break main
run
bt
```

Individual Catch2 tests can be debugged the same way, for example:

```bash
gdb --args build/test_poisson.exe
```

When launching from an editor on Windows, make sure the debugger and runtime
DLLs come from `D:\msys64\ucrt64\bin`; mixing MSYS2 environments or Visual
Studio toolchains with this build directory can produce confusing runtime
failures.

For VS Code or another MI-compatible debugger, set the debugger path to:

```text
D:\msys64\ucrt64\bin\gdb.exe
```

---

## Running the Poisson Solver

### Via CTest

The Poisson solver is exercised by the `test_poisson` test suite:

```bash
ctest --test-dir build --output-on-failure -R poisson
```

### Programmatically from C++

```cpp
#include "vela/simulation/PoissonSimulation.h"

vela::PoissonSimulation sim;
vela::VectorXd psi = sim.run("examples/pn_poisson_2d.json");
// psi(i) is the electrostatic potential [V] at node i
```

Integrations that need the solved mesh and net doping can use
`runWithResult()`:

```cpp
vela::PoissonResult result = sim.runWithResult("examples/pn_poisson_2d.json");
```

### From Python

After building with `-DVELA_ENABLE_PYTHON=ON` and setting `PYTHONPATH` to the
generated package directory:

```python
import vela

mesh = vela.load_mesh("examples/pn_diode_2d.json")
potential = vela.run_poisson("examples/pn_poisson_2d.json")
points = vela.run_dc_sweep("examples/pn_diode/simulation.json")
vela.write_vtk("examples/pn_diode/outputs/python_export.vtk")
```

The top-level Python API currently exposes `DeviceMesh`, `MaterialDatabase`,
`PoissonSimulation`, `DCSweep`, `load_mesh`, `run_poisson`, `run_dc_sweep`, and
`write_vtk`.

---

## Config JSON Schema

Poisson configurations use this shape:

```json
{
  "mesh_file": "pn_diode_2d.json",
  "output_vtk": "pn_poisson_2d.vtk",
  "doping": [
    { "region": "n_region", "donors": 1e23, "acceptors": 0.0 },
    { "region": "p_region", "donors": 0.0, "acceptors": 1e23 }
  ],
  "contacts": [
    { "name": "cathode", "bias": 0.0 },
    { "name": "anode", "bias": 0.0 }
  ]
}
```

Relative paths in `mesh_file`, `output_vtk`, and sweep output settings are
resolved relative to the directory containing the config JSON. VTK files can be
opened in ParaView to inspect electrostatic potential, net doping, and
drift-diffusion fields when available.

---

## Next Steps

- [x] Poisson FVM - assemble and solve the Poisson equation on the mesh
- [x] Scharfetter-Gummel flux - discretize carrier continuity fluxes
- [x] Gummel drift-diffusion sweep - run a PN diode IV curve
- [x] Python API - expose a compact pybind11 integration surface
- [ ] Voronoi control volumes - refine dual-cell areas and edge couplings
- [ ] Gmsh reader - replace `JsonMeshReader` with a real `.msh` importer
