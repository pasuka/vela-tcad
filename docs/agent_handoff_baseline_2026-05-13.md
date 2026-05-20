# Vela TCAD Baseline Handoff

Original date: 2026-05-13

Status: archived and refreshed for current documentation navigation. This file
captures the working baseline for future agents. For current user-facing
behavior, prefer `README.md`, `docs/architecture.md`,
`docs/config_schema.md`, and `docs/examples.md`.

## Current Architecture

Vela builds one core C++ library and one command-line runner:

- `vela_core`: static library containing mesh, material, physics,
  discretization, equation assembly, nonlinear solvers, simulation drivers,
  post-processing, and IO.
- `vela_example_runner`: CLI entry point used by examples and regression tests.
- Optional `_core` Python extension when `VELA_ENABLE_PYTHON=ON`.

Primary source layers:

- `src/core`: scalar types, constants, SI scaling, and `unit_scaling` helpers.
- `src/mesh`: mesh entities, topology, and box geometry.
- `src/material`: built-in and JSON-overridden materials.
- `src/boundary`: contact and explicit Poisson boundary parsing.
- `src/physics`: doping, carrier statistics, mobility, recombination, impact
  ionization, and bandgap narrowing.
- `src/discretization`: Bernoulli and Scharfetter-Gummel kernels.
- `src/equation`: Poisson, DD, and coupled-DD assemblers.
- `src/numerics`: residual norms and line search.
- `src/solver`: linear, Gummel, Newton, and validation logic.
- `src/simulation`: config parsing, Poisson driver, curve sweeps, and adaptive
  DC sweeps.
- `src/post`: current, field, terminal-charge, and stored-charge diagnostics.
- `src/io`: mesh, CSV, and VTK IO.

## Current Implemented Capabilities

Implemented solver and infrastructure capabilities:

- Poisson FVM/box assembly and solve.
- Scaled Poisson assembly for `scaling.mode = "unit_scaling"`.
- Drift-diffusion assembly with Scharfetter-Gummel fluxes.
- Gummel nonlinear iteration.
- Coupled Newton single-bias solve and DC-sweep method selection.
- Adaptive DC sweep with retry, step diagnostics, and CSV/VTK outputs.
- Region fixed charge and interface sheet/fixed/trap charge.
- Explicit Poisson Neumann, insulating, and symmetry boundaries.
- Ohmic, Dirichlet, metal-gate, and prototype Schottky contact paths.
- Material override loading through `materials_file`.
- Mobility, recombination, impact-ionization, and Slotboom bandgap-narrowing
  options.
- Terminal charge, multi-terminal quasi-static CV, stored charge, electric
  field, and contact-current diagnostics.
- Optional Python bindings for implemented C++ paths.

Known boundaries:

- Schottky is supported in Poisson/Gummel prototype paths; Newton rejects it.
- LDMOS and IGBT decks are trend fixtures, not calibrated power-device models.
- BV diagnostics are max-field/current/convergence smoke checks, not calibrated
  breakdown-voltage predictions.
- Quasi-static CV is terminal-charge differencing, not AC small-signal matrix
  extraction.

## Test Inventory

CTest targets include:

- C++ unit tests: scaling, mesh, Poisson, linear solver, box geometry,
  Bernoulli, SG flux, Gummel, high-doping Gummel, device stability, LDMOS DD,
  Newton, MOS cross-check, mixed-material MOS, mobility, impact ionization,
  recombination, line-search failure, DC sweep, solution validation, stored
  charge, electric-field diagnostics, boundary conditions, Schottky contact,
  Poisson Neumann, and MOS interface charge.
- `regression`: engineering example suite driven by `scripts/run_regression.py`.
- `reference_tcad_regression`: neutral CSV conversion/comparison fixture tests.
- `python_api`: optional pybind11 test when `VELA_ENABLE_PYTHON=ON`.
- `ascii_sources`: ASCII check for source/test/example paths covered by the
  CMake script.

Useful focused commands:

```bash
ctest --test-dir build --output-on-failure -R "poisson|boundary|interface|scaling"
ctest --test-dir build --output-on-failure -R "dd|newton|line_search|dc_sweep"
ctest --test-dir build --output-on-failure -R regression
ctest --test-dir build --output-on-failure -R reference_tcad_regression
```

## Agent Handoff Guidance

When changing solver, mesh, physics, or discretization behavior:

1. Read the relevant implementation path before editing.
2. Add or update the nearest Catch2 test.
3. Update `docs/config_schema.md` if a public config field changes.
4. Update `docs/examples.md` if a checked-in deck or regression case changes.
5. Run a focused CTest target first, then the full suite when behavior changed.

When changing examples or reference fixtures:

1. Keep generated example outputs out of the repository unless the task asks for
   golden artifacts.
2. Add regression runner checks that match implemented behavior only.
3. Keep language conservative: "prototype", "diagnostic", and "trend
   validation" are intentional boundaries for MOS and power-device fixtures.

## Historical Note

The original 2026-05-13 handoff identified Newton entry exposure, Gummel
absolute tolerance, line-search failure coverage, residual normalization,
warm-start, temperature propagation, SG consistency, geometry/linear-solver
caching, BGN implementation, and diagnostics history as backlog items. Those
items have since been implemented and are documented in the current README,
schema, example matrix, and weekly summary archive.
