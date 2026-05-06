# Vela TCAD

A lightweight **2-D semiconductor device drift-diffusion solver** based on the
Finite-Volume Method (FVM/Box integration) and Scharfetter-Gummel flux
discretization.

---

## Project Goal

Vela aims to provide an open, extensible prototype for simulating planar
semiconductor devices (p-n diodes, MOSFETs, …) entirely in C++20.  
This first stage delivers the **engineering skeleton** only.

---

## Current Stage – What Is Implemented

| Component | Status |
|-----------|--------|
| `Types.h` – `Real`, `Index`, `Point2`, `CellType` | ✅ |
| `PhysicalConstants.h` – q, kb, eps0, h, m0, T0, Vt_300 | ✅ |
| `ScalingSystem` – Debye scaling, scale/unscale helpers | ✅ |
| `DeviceMesh` – nodes, cells, edges, regions, contacts | ✅ |
| `JsonMeshReader` – reads `examples/pn_diode_2d.json` | ✅ |
| `VTKWriter` – VTK legacy ASCII for ParaView | ✅ |
| `MaterialDatabase` – Si and SiO2 built-in | ✅ |
| `DopingModel` – per-node Nd/Na from region JSON spec | ✅ |
| `PoissonAssembler` – FVM/Box electrostatic Poisson | ✅ |
| `LinearSolver` – Eigen SparseLU direct solver | ✅ |
| `PoissonSimulation` – end-to-end Poisson driver | ✅ |
| Unit tests (Catch2) | ✅ |

---

## Build Instructions

**Prerequisites:** CMake ≥ 3.20, a C++20 compiler, and the following
system packages:

```bash
# Ubuntu/Debian
sudo apt-get install libeigen3-dev nlohmann-json3-dev catch2
```

```bash
git clone https://github.com/pasuka/vela-tcad.git
cd vela-tcad
mkdir build && cd build
cmake ..
cmake --build . --parallel
```

---

## Running Tests

```bash
cd build
ctest --output-on-failure
```

Expected output: all tests in `test_scaling` and `test_mesh` pass.

---

## Running the Poisson Solver

### Via the unit tests

The Poisson solver is exercised automatically by the `test_poisson` test suite:

```bash
cd build
ctest --output-on-failure -R poisson
```

### Programmatically from C++

```cpp
#include "vela/simulation/PoissonSimulation.h"

vela::PoissonSimulation sim;
vela::VectorXd psi = sim.run("examples/pn_poisson_2d.json");
// psi(i) is the electrostatic potential [V] at node i
```

### Config JSON schema (`examples/pn_poisson_2d.json`)

```json
{
  "mesh_file":  "pn_diode_2d.json",
  "output_vtk": "pn_poisson_2d.vtk",
  "doping": [
    { "region": "n_region", "donors": 1e23, "acceptors": 0.0  },
    { "region": "p_region", "donors": 0.0,  "acceptors": 1e23 }
  ],
  "contacts": [
    { "name": "cathode", "bias": 0.0 },
    { "name": "anode",   "bias": 0.0 }
  ]
}
```

Relative paths in `mesh_file` and `output_vtk` are resolved relative to
the directory containing the config JSON.  The output VTK file can be opened
in **ParaView** to visualise the electrostatic potential and net doping.

---

## Next Steps

- [x] **Poisson FVM** – assemble and solve the Poisson equation on the mesh
- [ ] **Voronoi control volumes** – compute dual-cell areas and edge couplings
- [ ] **Scharfetter-Gummel flux** – discretize the carrier continuity equations
- [ ] **Newton solver** – coupled Poisson + continuity iteration
- [ ] **Gmsh reader** – replace `JsonMeshReader` with a real `.msh` importer

