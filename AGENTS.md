# AGENTS.md

## Project overview

This repository contains a C++20 CMake project for Vela TCAD, a lightweight 2-D semiconductor device drift-diffusion and Poisson solver.

## Dependencies

This is not a Python or npm project. Do not install Python or npm dependencies unless a future task explicitly adds tooling that requires them.

For Ubuntu/Debian-based Codex environments, install the C++ build dependencies with:

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

If the environment runs setup commands as root and `sudo` is unavailable, remove `sudo` from the commands above.

## Build

Prefer an out-of-tree Ninja build:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
```

If Ninja is unavailable, omit `-G Ninja` and use the default CMake generator.

## Test

Run the full test suite with:

```bash
ctest --test-dir build --output-on-failure
```

To run only the Poisson tests:

```bash
ctest --test-dir build --output-on-failure -R poisson
```

## Code style and workflow

- Keep the code compatible with C++20.
- Prefer adding or updating Catch2 tests when changing solver, mesh, physics, or discretization behavior.
- Keep generated build artifacts inside `build/` or another ignored out-of-tree build directory.
- Do not commit generated simulation outputs unless a task explicitly asks for them.
