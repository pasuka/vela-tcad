# AGENTS.md

## Project overview

This repository contains a C++20 CMake project for Vela TCAD, a lightweight 2-D semiconductor device drift-diffusion and Poisson solver. The core remains C++, with an optional pybind11 Python API behind `VELA_ENABLE_PYTHON`.

## Dependencies

This is not an npm project. Do not install Python package dependencies unless a future task explicitly adds tooling that requires them. The optional Python API uses system Python development headers and pybind11 through CMake.

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

For tasks that explicitly build or test the Python API, also install:

```bash
sudo DEBIAN_FRONTEND=noninteractive apt-get install -y \
  python3-dev \
  pybind11-dev
```

If the environment runs setup commands as root and `sudo` is unavailable, remove `sudo` from the commands above.

On Windows, the project development environment is MSYS2 UCRT64, typically installed at `D:\msys64`. Use the UCRT64 shell or put `D:\msys64\ucrt64\bin` and `D:\msys64\usr\bin` first on `PATH`. Use `python` instead of `python3` if the UCRT64 shell does not provide a `python3` alias. Install dependencies with:

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

For Windows tasks that explicitly build or test the Python API, also install:

```bash
pacman -S --needed mingw-w64-ucrt-x86_64-pybind11
```

## Build

Prefer an out-of-tree Ninja build:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
```

If Ninja is unavailable, omit `-G Ninja` and use the default CMake generator.

To include the optional Python bindings:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DVELA_ENABLE_PYTHON=ON
cmake --build build --parallel
```

From PowerShell on Windows, initialize the UCRT64 toolchain first:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
```

## Test

Run the full test suite with:

```bash
ctest --test-dir build --output-on-failure
```

To run only the Poisson tests:

```bash
ctest --test-dir build --output-on-failure -R poisson
```

When configured with `-DVELA_ENABLE_PYTHON=ON`, run the Python API test with:

```bash
ctest --test-dir build --output-on-failure -R python_api
```

## Debug

For Windows debugging, use UCRT64 GDB from `D:\msys64\ucrt64\bin` against the Debug build:

```bash
gdb --args build/vela_example_runner.exe --config examples/pn_diode/simulation.json
gdb --args build/test_poisson.exe
```

Avoid mixing MSYS2 UCRT64, MSYS2 CLANG64/MINGW64, and Visual Studio build outputs in the same build directory.

For VS Code or another MI-compatible debugger, point `miDebuggerPath` at `D:\msys64\ucrt64\bin\gdb.exe`.

## Code style and workflow

- Keep the code compatible with C++20.
- Prefer adding or updating Catch2 tests when changing solver, mesh, physics, or discretization behavior.
- Keep generated build artifacts inside `build/` or another ignored out-of-tree build directory.
- Do not commit generated simulation outputs unless a task explicitly asks for them.
