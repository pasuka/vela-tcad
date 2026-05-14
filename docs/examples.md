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

The default PN diode sweep uses `solver.method: "gummel"`. To exercise the
coupled Newton CLI entry point for a single equilibrium bias point, run:

```bash
build/vela_example_runner --config examples/pn_diode/newton_simulation.json
```

DC sweeps can also opt into Newton by setting `solver.method` (or
`solver.type`) to `"newton"` in the same `solver` block used for Newton
tolerances such as `max_iter`, `reltol`, `abstol`, `damping_factor`, and
`line_search`. Newton convergence and line-search merit checks use
`residual_norm: "block"` by default, which normalizes the Poisson, electron
continuity, and hole continuity residual blocks independently before combining
them. `residual_norm` only affects Newton solves; set it to `"l2"` for the
legacy raw L2 norm, or tune `residual_weights` / `residual_scales` with `psi`,
`phin`, and `phip` entries when a case needs explicit block weighting. Omitting
`solver.method` / `solver.type` keeps the historical Gummel sweep path.

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
- `simulation.json`: Poisson configuration with nonzero oxide `fixed_charge_m3` and a
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

## NMOS Drift-Diffusion Stability (`examples/nmos2d_dd`)

**Physical meaning.** A deliberately coarse, all-silicon NMOS-like cross-section
exercises the drift-diffusion Gummel path with four independently biased
terminals (`body`, `source`, `gate`, and `drain`). The source/drain implant
regions are represented as abrupt n+ boxes next to a p-type body/channel. This
is a CI stability baseline rather than a calibrated transistor model; it keeps
the mesh and voltage range small so that regressions catch solver instability
without requiring a production TCAD deck.

**Bias conditions.** Body and source are grounded, the gate contact is held at
0.05 V, and the drain is swept from 0 V to 0.10 V in 0.05 V steps.

**Run directly.**

```bash
mkdir -p examples/nmos2d_dd/outputs
build/vela_example_runner --config examples/nmos2d_dd/simulation.json
```

**Outputs and acceptance metrics.**

- `examples/nmos2d_dd/outputs/nmos2d_dd_iv.csv`: three drain-bias rows with
  `converged=1` for each declared point.
- `examples/nmos2d_dd/outputs/nmos2d_dd_sweep_*.vtk`: finite potential,
  quasi-Fermi, carrier-density, and net-doping fields.
- Regression requires all expected CSV/VTK files to exist, no NaN/Inf tokens in
  those outputs, no retries, and `attempted_step` / `accepted_step` magnitudes
  no larger than 0.05 V.

## PMOS Drift-Diffusion Stability (`examples/pmos2d_dd`)

**Physical meaning.** This is the polarity-complement counterpart to the NMOS
stability case. It uses a coarse all-silicon mesh with abrupt p+ source/drain
implants in an n-type body to exercise sign handling, descending DC sweep
targets, and contact-bias application in the drift-diffusion solver.

**Bias conditions.** Body and source are grounded, the gate contact is held at
-0.05 V, and the drain is swept from 0 V to -0.10 V in -0.05 V steps.

**Run directly.**

```bash
mkdir -p examples/pmos2d_dd/outputs
build/vela_example_runner --config examples/pmos2d_dd/simulation.json
```

**Outputs and acceptance metrics.**

- `examples/pmos2d_dd/outputs/pmos2d_dd_iv.csv`: three drain-bias rows with
  `converged=1` for each declared point.
- `examples/pmos2d_dd/outputs/pmos2d_dd_sweep_*.vtk`: finite potential,
  quasi-Fermi, carrier-density, and net-doping fields.
- Regression requires all expected CSV/VTK files to exist, no NaN/Inf tokens in
  those outputs, no retries, and `attempted_step` / `accepted_step` magnitudes
  no larger than 0.05 V.

## LDMOS Poisson Stability (`examples/ldmos2d_poisson`)

**Physical meaning.** A coarse lateral DMOS-like Poisson-only deck covers a
p-body / n-drift lateral junction with n+ source/drain regions and a field-oxide
slab. It is intended to catch Poisson assembly and material-interface
regressions under a small reverse drain bias before a full LDMOS
drift-diffusion model is introduced.

**Bias conditions.** Body and source are grounded, the gate/field-plate contact
is held at 0.2 V, and the drain contact is held at 1.0 V.

**Run directly.**

```bash
mkdir -p examples/ldmos2d_poisson/outputs
build/vela_example_runner --config examples/ldmos2d_poisson/simulation.json
```

**Outputs and acceptance metrics.**

- `examples/ldmos2d_poisson/outputs/ldmos2d_poisson.vtk`: electrostatic
  potential and net doping across silicon and field oxide.
- Regression requires the VTK file to exist, be nonempty, and contain no
  NaN/Inf tokens.

## IGBT Poisson Stability (`examples/igbt2d_poisson`)

**Physical meaning.** A coarse vertical IGBT-like Poisson baseline stacks a
p-type collector, n buffer, lightly doped n drift, and top p-base/emitter
region. It targets off-state electrostatic stability and high-gradient doping
assembly in a CI-friendly mesh.

**Bias conditions.** The collector is biased at 1.0 V while the top gate and
emitter contacts are grounded.

**Run directly.**

```bash
mkdir -p examples/igbt2d_poisson/outputs
build/vela_example_runner --config examples/igbt2d_poisson/simulation.json
```

**Outputs and acceptance metrics.**

- `examples/igbt2d_poisson/outputs/igbt2d_poisson.vtk`: electrostatic
  potential and net doping through the vertical stack.
- Regression requires the VTK file to exist, be nonempty, and contain no
  NaN/Inf tokens.

## Adding Future Regression Cases

Use this minimal checklist when adding a new device regression case:

1. Create `examples/<device_name>/` with a small `mesh.json` and
   `simulation.json` that exercise only implemented solver features.
2. Keep generated files under `examples/<device_name>/outputs/`; do not commit
   those outputs unless a task explicitly asks for golden artifacts.
3. Add the case to `EXAMPLES` in `scripts/run_regression.py` with every expected
   CSV/VTK file that the runner should verify.
4. Add the smallest physics-specific regression check that matches the current
   implementation, such as convergence flags, NaN/Inf scanning, interface
   continuity, current-direction sanity, or a monotonic trend already supported
   by the model.
5. Document the case in this file with physical meaning, inputs, direct run
   command, outputs, and any special checks.
6. If the new case exposes solver, mesh, physics, or discretization behavior,
   add or update the closest Catch2 test in `tests/` instead of relying only on
   the example.
7. Build and run `ctest --test-dir build --output-on-failure -R regression`
   before committing; run the full CTest suite when core code changes are part
   of the regression update.
