# Architecture

Vela is a C++20 CMake project centered on the `vela_core` static library.
`vela_example_runner` is the command-line executable used by JSON simulation
decks and reference-case workflows. Optional pybind11 bindings expose a thin Python API when
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

The experimental `vela::experimental` NonLocalPath contracts in
`physics/NonlocalPathBandToBand.h` define SI configuration, actual band-state
dependencies, path records, separate signed pair-source integrals and versioned
remote Jacobian support. `NonlocalPathGeometry` adds private autonomous Tri3
path traversal, endpoint roots and specular reflection. `NonlocalPathPhononKernel`
evaluates a private homogeneous-Si, isothermal P1 phonon-assisted rate with
explicit reservoir windows and a unique Keldysh crossover. It preserves signed
generation/recombination and logarithms for underflow diagnostics.
`NonlocalPathSourceAssembly` adds geometry-derived vertex-box quadrature,
separate endpoint scatter and contact ledgers. Its frozen
`CoupledDDAssembler::experimentalNonlocalResidual` diagnostic adds both sources
inside the actual residual loop, before gauge/scaling/contact replacement;
partial path sets cannot enter that diagnostic. Constitutive band identity and
native window/source collection remain unqualified; derivatives are unevaluated.
The named private `NonlocalPathCellSource` alternative adds centroid chords,
Delaunay cell quadrature, band-fraction integrals and first-crossing energy
collection. It preserves the canonical N3 entry point. Its explicit frozen
population residual audit retains the strict reconstructed-state gate and
checks physical temperature; neither audit qualifies production dispatch or
the full native RHS. See the [C++ N4 record](validation/bvmethods_nonlocalpath_n4_cpp_closure_2026-10-16.md).
Frozen-density audits now keep ordinary SRH on the actual substituted population,
while the explicit quantum SRH coupling remains separate. A dedicated Poisson
observer uses the frozen charge, retains absolute dielectric edge contributions,
and matches the actual endpoint scaling/subtraction order. Private original-statistics,
Gamma transport and element-source controls do not add production configuration
dispatch; see the [actual-assembly follow-up](validation/bvmethods_nonlocalpath_n4_actual_assembly_2026-10-16.md).
These APIs are not dispatched by simulation configuration. Private entry points
are `test_nonlocalpath_contract`, `test_nonlocalpath_geometry`,
`test_nonlocalpath_phonon`, `test_nonlocalpath_source` and the uninstalled
geometry/source probes. See the [N1 record](validation/bvmethods_nonlocalpath_n1_2026-10-16.md)
and [N2 geometry record](validation/bvmethods_nonlocalpath_n2_2026-10-16.md), plus
the [N3 kernel record](validation/bvmethods_nonlocalpath_n3_2026-10-16.md) and
[N4 source investigation](validation/bvmethods_nonlocalpath_n4_2026-10-16.md).
The public local E2 model remains in `BandToBandTunnelingModel`.

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

The experimental Newton report collector and independent `NewtonAuditScope`
are opt-in, bounded observers of existing classical-loop values. They copy
state, J/R, directions and line-search results without evaluating physics.
Reports and `ExactValueBits` comparison encoding run outside the solve; the
encoding preserves every aggregate field and floating-point bit pattern using
header-only Boost.PFR. It is not a restart format or a new convergence rule.
Split, quantum, recenter and recovery paths require separate qualification.

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

- Experimental `ialmob` transport is assembled through `IalTransport` and
  `IalElementMobility`: live element mobility and its nine potential derivatives
  feed the SG residual/Jacobian, terminal currents, and density recovery. The
  scalar `MobilityModel` API rejects this model to prevent a silent bulk fallback.
  It currently requires explicit native AverageBox geometry evidence, Si/SiO2
  Tri3 geometry, and a 300 K classical configuration without avalanche. See
  [current LDMOS validation](validation/templates_ldmos_current_status.md) for
  the distinction between local formula checks and self-consistent D4 curves.
- Device-level calibration and cross-tool claims require a checked-in
  `reference_tcad/` fixture with an explicit acceptance boundary.
- `LatticeHeatAssembler` and the `lattice_heat_probe` executable provide
  independent steady Tri3 conduction with temperature-dependent conductivity,
  prescribed cell heat sources and explicit Robin thermodes. Geometry conversion
  is explicit; thermal quantities are SI per unit width. This prescribed-source
  tool remains separate from electrical DC sweeps. See the
  [thermal operator report](validation/templates_ldmos_d0_thermal_operator_2026-09-12.md).
- `ElectrothermalAssembler` and the shared `solveElectrothermalPoint` service couple
  four SI unknowns per node: potential, electron/hole quasi-Fermi potentials and
  lattice temperature. `SiliconThermalPhysics`, temperature-aware IALMob and
  `ThermalSgCurrent` provide the analytic Jacobian. Default lattice heat uses the
  same discrete carrier currents and band-edge work; Thermodynamic, Peltier and
  RecGenHeat are not enabled. Explicit native geometry, ideal electron contacts
  and finite hole recombination boundaries are scoped to the audited LDMOS model.
  `simulation_type: "electrothermal_dc_sweep"` uses this service in process for
  initialization, Poisson gate prebias, four-equation drain continuation, failed
  step reduction, checkpoints and current/temperature/heat output. The separate
  `electrothermal_probe` CLI calls the same point solver for compatibility.
  Legacy `dc_sweep` remains isothermal. The configuration and integrated
  qualification boundaries are documented in the [schema](config_schema.md#explicit-electrothermal-dc-entry)
  and [current LDMOS status](validation/templates_ldmos_current_status.md).
- BV sweeps report diagnostic max field/current jump/non-convergence markers,
  not calibrated breakdown voltages.
- CV sweeps use finite-difference terminal charge, not AC small-signal
  admittance extraction.
- Schottky support is available in the Poisson and Gummel prototype paths.
  Newton rejects Schottky contacts with a clear error.
