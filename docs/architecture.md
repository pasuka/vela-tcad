# Architecture

Vela is a C++20 CMake project centered on the `vela_core` static library.
`vela_example_runner` is the command-line executable used by examples and the
regression suite. Optional pybind11 bindings expose a thin Python API when
`VELA_ENABLE_PYTHON=ON`.

## Source Layout

| Area | Paths | Responsibility |
| --- | --- | --- |
| Core | `include/vela/core`, `src/core` | Scalar types, constants, SI scaling helpers, and `unit_scaling` input conversion. |
| Mesh | `include/vela/mesh`, `src/mesh` | Mesh entities, JSON mesh topology, edge construction, and box/control-volume geometry. |
| Material | `include/vela/material`, `src/material` | Built-in Si/SiO2 parameters and optional material override loading. |
| Boundary | `include/vela/boundary`, `src/boundary` | Contact and explicit Poisson boundary type parsing. |
| Physics | `include/vela/physics`, `src/physics` | Doping, carrier statistics, mobility, recombination, impact ionization, and bandgap narrowing. |
| Discretization | `include/vela/discretization`, `src/discretization` | Bernoulli and Scharfetter-Gummel flux kernels. |
| Equation | `include/vela/equation`, `src/equation` | Poisson, drift-diffusion, and coupled-DD assembly. |
| Numerics | `include/vela/numerics`, `src/numerics` | Residual norms and line search. |
| Solver | `include/vela/solver`, `src/solver` | Linear solve, Gummel iteration, coupled Newton solve, and solution validation. |
| Simulation | `include/vela/simulation`, `src/simulation` | Config parsing, Poisson driver, curve sweep helpers, and adaptive DC sweep orchestration. |
| Post | `include/vela/post`, `src/post` | Contact current, electric field, terminal charge, and stored charge diagnostics. |
| IO | `include/vela/io`, `src/io` | Mesh reading plus CSV/VTK output. |
| Tools | `src/tools` | `vela_example_runner` CLI. |
| Python | `bindings`, `python/vela` | Optional pybind11 wrapper and small Python helpers. |

## Solver Paths

Poisson-only:

1. Parse config and resolve relative paths from the config file directory.
2. Load mesh, material overrides, doping, fixed/interface charge, contacts, and
   optional explicit boundaries.
3. Assemble the finite-volume Poisson system.
4. Solve with the linear solver.
5. Write physical-voltage VTK output.

Drift-diffusion Gummel sweep:

1. Parse a `dc_sweep` deck.
2. For each accepted bias point, solve Poisson and electron/hole continuity
   subproblems with Scharfetter-Gummel fluxes.
3. Apply adaptive step control, retries, convergence checks, and optional VTK
   snapshots.
4. Write CSV rows with current, convergence, sweep, and requested diagnostic
   columns.

Coupled Newton:

1. Build the coupled Poisson/electron/hole residual.
2. Use analytic or finite-difference Jacobian options.
3. Apply damping and optional line search.
4. Use block or L2 residual norms with configurable weights/scales.
5. Optionally collect diagnostic history.

Newton is reachable through `simulation_type: "newton"` for a single configured
bias point and through `solver.method: "newton"` in DC sweeps. Gummel remains
the default DC-sweep method when no solver method is specified.

## Unit Scaling Boundary

Omitting `scaling` keeps legacy SI input behavior. Setting:

```json
"scaling": { "mode": "unit_scaling" }
```

uses common external TCAD input units at the schema boundary: um, cm^-3,
cm^2/(V s), V/cm, cm^-2, V, K, and eV. Inputs are normalized to the existing
solver kernel. The Poisson path also uses a scaled assembly and restores
physical voltage before returning or writing output.

## Implementation Boundaries

- Mixed-material MOS and power-device examples are engineering smoke tests.
- LDMOS and IGBT decks validate signs, trends, finite outputs, and diagnostic
  fields; they are not calibrated device models.
- BV sweeps report diagnostic max field/current jump/non-convergence markers,
  not calibrated breakdown voltages.
- CV sweeps use finite-difference terminal charge, not AC small-signal
  admittance extraction.
- Schottky support is available in the Poisson and Gummel prototype paths.
  Newton rejects Schottky contacts with a clear error.
