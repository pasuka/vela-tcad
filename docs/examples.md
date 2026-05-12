# Engineering Examples and Regression Suite

Vela ships three small engineering examples that are intentionally coarse so
they can run quickly in CI while still exercising the Poisson,
drift-diffusion, post-processing, and file-output paths.

## Common Workflow

Build the project first:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
```

On Windows, run the same commands from the MSYS2 **UCRT64** shell. If MSYS2 is
installed at `D:\msys64`, the corresponding tools live under
`D:\msys64\ucrt64\bin`. Use `python` instead of `python3` below if the shell
does not provide a `python3` alias.

Run all examples through the regression harness:

```bash
python3 scripts/run_regression.py --runner build/vela_example_runner
```

CTest exposes the same entry point:

```bash
ctest --test-dir build --output-on-failure -R regression
```

The regression script copies examples into `build/regression_output/`, creates
per-example `outputs/` folders, runs `vela_example_runner --config
simulation.json`, verifies expected files, checks convergence indicators, scans
CSV/VTK text for NaN/Inf values, and writes
`build/regression_output/regression_summary.json`.

## Python API Example

The optional pybind11 API can run the same PN diode sweep from Python. Configure
with `VELA_ENABLE_PYTHON=ON` and make the generated package importable:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug -DVELA_ENABLE_PYTHON=ON
cmake --build build --parallel
PYTHONPATH=build/python/Debug python3 examples/python/run_pn_diode.py
```

On Windows PowerShell, set `PYTHONPATH` with:

```powershell
$env:PYTHONPATH = "build\python\Debug"
python examples\python\run_pn_diode.py
```

The script calls `vela.run_dc_sweep()` on `examples/pn_diode/simulation.json`
and writes the usual CSV/VTK outputs under `examples/pn_diode/outputs/`.

CTest registers the Python coverage as `python_api` when the extension is
enabled:

```bash
ctest --test-dir build --output-on-failure -R python_api
```

## PN Diode (`examples/pn_diode`)

**Physical meaning.** A tiny 2-D silicon p-n diode validates the
self-consistent drift-diffusion flow and contact-current extraction. The anode
is forward swept from 0 V to 0.5 V while the cathode is grounded. The
regression checks that every sweep point converges and that the forward-current
magnitude increases from the first to the final bias point.

**Inputs.**

- `mesh.json`: 1 um x 1 um triangular silicon mesh split into `p_region` and
  `n_region`, with `anode` and `cathode` contacts.
- `simulation.json`: region doping, contact biases, Gummel solver controls, and
  a DC sweep definition.

**Run directly.**

```bash
mkdir -p examples/pn_diode/outputs
build/vela_example_runner --config examples/pn_diode/simulation.json
```

**Outputs.**

- `examples/pn_diode/outputs/pn_iv.csv`: voltage, electron current, hole
  current, total current, convergence flag, and iteration count.
- `examples/pn_diode/outputs/pn_sweep_*.vtk`: potential, quasi-Fermi levels,
  carrier densities, and net doping snapshots for converged sweep points.

## MOS Capacitor (`examples/moscap`)

**Physical meaning.** A two-material Si/SiO2 MOS capacitor validates Poisson
assembly across a dielectric interface. The body contact is grounded and the
gate is biased at 1 V. No mobile-carrier equation is solved in this first
regression; the case is Poisson-only.

**Inputs.**

- `mesh.json`: layered triangular mesh with a silicon slab, a silicon-dioxide
  slab, shared interface nodes at `y = 0.30 um`, and `body`/`gate` contacts.
- `simulation.json`: Poisson configuration with zero fixed charge and a
  `regression` block describing interface probe locations and tolerances.

**Run directly.**

```bash
mkdir -p examples/moscap/outputs
build/vela_example_runner --config examples/moscap/simulation.json
```

**Outputs and checks.**

- `examples/moscap/outputs/moscap_poisson.vtk`: electrostatic potential and net
  doping.
- Regression confirms the shared-node potential jump at the Si/SiO2 interface
  is effectively zero and that the normal electric-displacement estimates on
  the silicon and oxide sides are approximately continuous.

## Simplified NMOS (`examples/nmos2d`)

**Physical meaning.** A coarse rectangular NMOS cross-section exercises a mixed
semiconductor/oxide mesh with four contacts: `source`, `drain`, `body`, and
`gate`. The first version is Poisson-only and uses rectangular n+ source/drain
implant regions plus a p-type body as a CI-friendly approximation of analytic
implant profiles.

**Inputs.**

- `mesh.json`: p-body silicon, n+ source/drain implant boxes, top gate oxide,
  and four contacts.
- `simulation.json`: fixed doping for each region and a single Poisson bias
  point (`Vg = 0.8 V`, `Vd = 0.1 V`, source/body grounded).

**Run directly.**

```bash
mkdir -p examples/nmos2d/outputs
build/vela_example_runner --config examples/nmos2d/simulation.json
```

**Outputs.**

- `examples/nmos2d/outputs/nmos_poisson.vtk`: potential and net doping for
  quick inspection in a VTK viewer.

## Adding Future Regression Cases

1. Add a new example directory with `mesh.json` and `simulation.json`.
2. Keep meshes small and outputs under an ignored `outputs/` subdirectory.
3. Extend `EXAMPLES` in `scripts/run_regression.py` with expected files and any
   physics-specific checks.
4. Run `ctest --test-dir build --output-on-failure -R regression` before
   committing.
