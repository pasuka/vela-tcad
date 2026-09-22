# SLOT source-free poly equilibrium continuation plan

> **For agentic workers:** Use superpowers:subagent-driven-development for bounded implementation and review.

**Goal:** Make the already demonstrated, physically exact isolated-poly equilibrium restart reproducible without changing model equations or acceptance gates.

**Architecture:** A standalone Python helper validates the narrow SLOT profile and connected poly components, writes a candidate restart plus provenance, and never accepts a point. The existing continuation corrector invokes it only with a default-off CLI option, before the existing all-active audit and full coupled validation. Preserve every non-poly row and the input file. Any invalid applicability condition raises rather than silently projecting.

**Tech Stack:** Python stdlib, existing SLOT continuation and unittest regressions; UCRT64 Release actual-device verification.

The user already authorized stable initialization/continuation and all remaining validation work. The 15 V experiment establishes this specific numerical repair: 350 isolated Poly nodes, one zero-biased component, max correction 1.046e-21 V; original full/local gates then pass. This is a candidate initialization operation, not a density floor, acceptance relaxation, parameter fit, or new physical model.

- [x] Write failing tests for component/contact qualification, preservation of non-poly state, Boltzmann density transformation, split QF references/increments, invalid states/configurations.
- [x] Implement `scripts/slot_poly_equilibrium.py` with guarded `write_poly_equilibrium_restart(document, input_path, output_path, config_directory, bias)`; return provenance. Restrict to supported Boltzmann SLOT IALMob constant-poly profile, explicit recombination scope excluding Poly; avalanche absent or canonical current-driven GradQF profile. All component contacts must be ohmic and one common bias; reject floating or shared semiconductor nodes.
- [x] Add default-off `--restore-poly-equilibrium` to the existing continuation corrector, preserving full coupled + final all-active audit and update-history qualification. No acceptance shortcuts.
- [x] Run corresponding Python regressions and real Release continuation from the last qualified local point; verify independent 15/20 V points and Jacobian separately.
- [x] Review implementation and scope, update evidence and report remaining native-model mismatch explicitly.

Verification: helper 18 + continuation 20 Python tests pass; both registered targeted Release CTest entries pass. Independent review found and verified a versioned-material transport-key guard repair. Real 15/20 V helper output exactly matches previously fully validated candidates; new local continuation has qualified 10.5 V. Higher continuation and BVDS remain separate ongoing numerical tasks.
