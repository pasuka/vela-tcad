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
| `MaterialDatabase` - Si/SiO2 built in plus optional JSON overrides | done |
| `DopingModel` - per-node Nd/Na from region JSON spec | done |
| Boundary/contact schema parser (`contacts.type`, `boundaries`) | done |
| Explicit Poisson boundaries (Neumann/insulating/symmetry) | done |
| Schottky contact prototype (Poisson + Gummel sweep path) | done |
| `PoissonAssembler` - FVM/Box electrostatic Poisson | done |
| `ScharfetterGummel` and drift-diffusion assembly | done |
| `LinearSolver`, `GummelSolver`, and Newton helpers | done |
| `PoissonSimulation` - end-to-end Poisson driver | done |
| `DCSweep` - bias sweep with CSV/VTK outputs | done |
| Optional pybind11 Python API | done |
| Unit and regression tests | done |
| Power-device prototype examples (LDMOS/IGBT IV/CV/BV/stored-charge trend decks) | done (prototype / not calibrated) |

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
gdb --args build/vela_example_runner.exe --config examples/pn_diode/simulation_iv.json
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
iv_points = vela.run_iv_curve("examples/pn_diode/simulation_iv.json")
cv_points = vela.run_cv_curve("examples/pn_diode/simulation_cv.json")
bv_points = vela.run_bv_curve("examples/pn_diode/simulation_bv.json")
vela.write_vtk("examples/pn_diode/outputs/python_export.vtk")
```

The top-level Python API currently exposes `DeviceMesh`, `MaterialDatabase`,
`PoissonSimulation`, `DCSweep`, `load_mesh`, `run_poisson`, `run_dc_sweep`,
`run_iv_curve`, `run_cv_curve`, `run_bv_curve`, and `write_vtk`. Sweep points
are dictionaries with the C++ sweep data plus metadata such as `curve_type`,
`bias_contact`, `current_contact`, `convergence_diagnostics`, `output_csv`,
and per-point `output_vtk` paths. The Python API is a thin pybind11 wrapper
around documented C++ functionality; only capabilities already exercised and
tested in the C++ core are documented as supported from Python.

---

## Config JSON Schema

Full schema reference: [`docs/config_schema.md`](docs/config_schema.md).

Examples and support matrix: [`docs/examples.md`](docs/examples.md).

Poisson configurations use this shape:

```json
{
  "mesh_file": "pn_diode_2d.json",
  "materials_file": "materials.json",
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

Contact entries support an optional `type` field that selects the boundary
model. The legacy untyped form is equivalent to `"type": "ohmic"`, so existing
decks keep working unchanged. The parser accepts case-insensitive names with
either `-` or `_` separators (for example `metal_gate` and `Metal-Gate` both
resolve to `metal_gate`). Recognised values are `ohmic`, `dirichlet`,
`metal_gate`, `schottky`, and `floating`. The Poisson driver supports
`ohmic`, `dirichlet`, `metal_gate`, and a Schottky prototype path; `floating`
is still reserved and raises a clear error when requested. In DD sweeps,
Schottky is available on the Gummel path and rejected on the Newton path with
a clear error.

`materials_file` is optional. When present, it is loaded after the built-in
`Si` and `SiO2` materials, so entries with the same `name` override built-ins
and entries with new names can be referenced by mesh regions. The same optional
field is honored by Poisson runs, DC sweep runs, and the Python helpers that run
those same C++ config paths. A materials file may be either an object with a
`materials` array, a top-level array, or an object keyed by material name:

```json
{
  "materials": [
    {
      "name": "Si",
      "eps_r": 11.7,
      "ni": 1.0e16,
      "mun": 0.135,
      "mup": 0.048,
      "bandgap_eV": 1.12,
      "electron_affinity_eV": 4.05,
      "Nc_m3": 2.8e25,
      "Nv_m3": 1.04e25,
      "temperature_K": 300.0
    }
  ]
}
```

No `scaling` field keeps the legacy SI input behavior used by all existing
decks. In legacy SI mode, numeric config values use SI units unless the field
name explicitly states otherwise: lengths are meters, concentrations such as
`donors`, `acceptors`, `ni`, `Nc_m3`, and `Nv_m3` are in `m^-3`, mobilities are
in `m^2/(V s)`, electric fields are in `V/m`, interface sheet densities are in
`m^-2`, voltages are in volts, `temperature_K` is in kelvin, and band gap or
electron affinity fields use electron-volts as indicated by their `_eV` suffix.
To convert concentrations from `cm^-3` to `m^-3`, multiply by `1e6`; for
example, `1e17 cm^-3` is `1e23 m^-3`.

The public schema also defines an optional unit input mode:

```json
"scaling": { "mode": "unit_scaling" }
```

`unit_scaling` is an external TCAD common-unit input mode, not a calibration or
accuracy setting. It interprets lengths as `um`, concentrations as `cm^-3`,
mobilities as `cm^2/(V s)`, electric fields as `V/cm`, interface sheet
densities as `cm^-2`, voltages as `V`, temperatures as `K`, and energies as
`eV`. Vela normalizes those inputs to the existing SI solver kernel when it
reads mesh coordinates, material overrides, doping, fixed/interface charge,
mobility settings, and electric-field-related solver settings. The only
accepted `scaling.mode` value is `unit_scaling`; `si` and vendor-specific mode
aliases are intentionally not part of the schema.

Relative paths in `mesh_file`, optional `materials_file`, `output_vtk`, and
sweep output settings are resolved relative to the directory containing the
config JSON. VTK files can be opened in ParaView to inspect electrostatic
potential, net doping, and drift-diffusion fields when available.

Drift-diffusion sweeps use the Gummel solver unless the `solver` block
explicitly selects Newton:

```json
{
  "simulation_type": "dc_sweep",
  "solver": {
    "method": "newton",
    "max_iter": 10,
    "reltol": 1e-8,
    "abstol": 1e-18,
    "damping_factor": 1.0,
    "line_search": true,
    "warm_start": false,
    "verbose": false
  }
}
```

Use `"method": "gummel"` or omit the field to keep the default sweep path.
For Newton, `"warm_start": true` preserves supplied quasi-Fermi potentials
(`phin`/`phip`) when continuing from a previous solution; the default `false`
keeps the conservative cold-start behavior that resets interior quasi-Fermi
potentials before the coupled solve.

Gummel and Newton solver configs also accept `temperature_K` (default `300.0`)
which sets the thermal voltage used by contact built-in potentials, drift-
diffusion assembly, and DC-sweep contact-current post-processing. Gummel solver
configs also accept `abstol`; convergence is reported when either all relative
update checks pass or all absolute update norms are below that positive
threshold. A single-bias coupled Newton CLI run is available with
`"simulation_type": "newton"`; see
`examples/pn_diode/newton_simulation.json`.

---

## Development Notes

- 本周开发任务总结已同步到 [`docs/weekly_development_summary_2026-05-14.md`](docs/weekly_development_summary_2026-05-14.md)，按提交记录归纳了 2026-05-11 至 2026-05-14 的求解器、物理模型、性能、可观测性和测试进展。
- 2026-05-13 合并 Backlog 的 B01-B14 已全部完成；状态同步见 [`docs/multi_agent_merged_backlog_2026-05-13.md`](docs/multi_agent_merged_backlog_2026-05-13.md#31-2026-05-14-执行状态同步)。

---

## Next Steps

The roadmap below tracks future work only. Items already listed as `done` in
[Current Stage](#current-stage---what-is-implemented) are treated as implemented
baseline capabilities, while every milestone here remains planned until its
source changes and tests land.

### P0 - Numerical Stability

**Goal:** make the existing Poisson and drift-diffusion paths more robust on
coarse, highly doped, and multi-material meshes before expanding device
coverage.

- **Main source paths:** `src/discretization/`, `src/equation/`,
  `src/numerics/`, `src/solver/`, `src/mesh/BoxGeometryBuilder.cpp`, and the
  related public headers.
- **Test targets:** strengthen `tests/test_poisson.cpp`,
  `tests/test_dd_gummel.cpp`, `tests/test_sg_flux.cpp`,
  `tests/test_box_geometry.cpp`, `tests/test_newton_solver.cpp`, and the
  engineering regression suite (`ctest --test-dir build --output-on-failure -R regression`).
- **Planned milestones:** refine dual/control-volume geometry and edge
  couplings; add convergence guards for nonlinear solves; expand regression
  thresholds that catch NaN/Inf fields, non-converged sweeps, and current or
  potential sign regressions.

### P1 - MOS Physics

**Goal:** improve MOS capacitor and NMOS electrostatics/transport modeling
without marking Poisson-only examples as complete MOSFET simulation support.

- **Main source paths:** `src/material/`, `src/physics/`, `src/equation/`,
  `src/simulation/PoissonSimulation.cpp`, `examples/moscap/`, and
  `examples/nmos2d/`.
- **Test targets:** extend `tests/test_poisson.cpp`, `tests/test_mobility.cpp`,
  `tests/test_recombination.cpp`, and regression checks for `examples/moscap`
  and `examples/nmos2d`.
- **Planned milestones:** add explicit MOS interface/fixed-charge checks, improve
  oxide/semiconductor boundary validation, and introduce MOS-focused carrier
  physics tests only after the corresponding implementation exists.

### P2 - Device-Level Sweep

**Goal:** make bias sweeps reliable as device regression artifacts rather than
only single example runs.

- **Main source paths:** `src/simulation/DCSweep.cpp`,
  `src/tools/vela_example_runner.cpp`, `src/post/ContactCurrent.cpp`,
  `src/io/`, `scripts/run_regression.py`, and `examples/*/simulation_*.json`.
- **Test targets:** expand `tests/test_dc_sweep.cpp`,
  `tests/test_dd_gummel.cpp`, regression summary checks, and CSV/VTK output
  validation in `scripts/run_regression.py`.
- **Planned milestones:** add reusable sweep specifications for more devices,
  record stable per-bias convergence metadata, and tighten regression checks for
  monotonicity or physically expected trends where the implemented model
  supports them.

### P3 - Python/API/Documentation

**Goal:** keep the C++ core, optional Python bindings, examples, and docs aligned
with implemented behavior.

- **Main source paths:** `bindings/pyvela.cpp`, `python/vela/`,
  `examples/python/`, `docs/`, `README.md`, and public headers under
  `include/vela/`.
- **Test targets:** maintain `tests/python/test_python_api.py`, CTest's
  `python_api` target when `VELA_ENABLE_PYTHON=ON`, documentation examples, and
  the full `ctest --test-dir build --output-on-failure` suite.
- **Planned milestones:** document new regression-case patterns, keep Python
  wrappers limited to implemented C++ capabilities, and add API examples only
  when they are covered by tests.
