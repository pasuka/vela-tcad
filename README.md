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

## Next Steps

- [ ] **Poisson FVM** – assemble and solve the Poisson equation on the mesh
- [ ] **Voronoi control volumes** – compute dual-cell areas and edge couplings
- [ ] **Scharfetter-Gummel flux** – discretize the carrier continuity equations
- [ ] **Newton solver** – coupled Poisson + continuity iteration
- [ ] **Gmsh reader** – replace `JsonMeshReader` with a real `.msh` importer

