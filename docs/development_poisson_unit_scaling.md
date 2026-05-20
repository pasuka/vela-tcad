# Poisson Unit-Scaling Developer Notes

This document describes the Poisson-specific scaled (dimensionless) solve path
introduced for public input mode:

```json
"scaling": { "mode": "unit_scaling" }
```

Scope:
- Poisson assembly and solve path only.
- DD/Gummel/Newton decks also use the shared `unit_scaling` input conversion
  at the schema boundary, but this note only describes the Poisson-specific
  scaled assembly path. Drift-diffusion assembly still receives normalized SI
  values.

## Goals

- Keep legacy behavior exactly when `scaling` is omitted.
- Support common TCAD input units at the schema boundary in `unit_scaling`.
- Use a true scaled Poisson assembly path in that mode.
- Restore physical potential before returning simulation results and writing VTK.

## Public compatibility rules

- Do not use `scaling.system` aliases.
- Supported value remains only `scaling.mode = "unit_scaling"`.
- Existing field names are preserved for compatibility, including
  `normal_displacement_C_per_m2`.

## Variables and reference scales

Scaled unknowns and coordinates:
- `psi_hat = psi / V0`
- `x_hat = x / L0`
- `C_hat = C / C0`

Where:
- `V0` is a potential scale from `UnitScalingSystem`.
- `L0` is a characteristic length used by the scaling system.
- `C0` is a concentration scale used by the scaling system.

Poisson scaled assembly uses a permittivity reference `eps_ref` so local
material contrast remains explicit through `eps/eps_ref`.

## Assembly strategy implemented

The existing SI Poisson assembly computes matrix `A_si` and RHS `b_si` from:
- edge flux terms with material `eps`
- volumetric charge terms (net doping and fixed charge)
- interface sheet/fixed/trap charge terms
- Neumann displacement boundary terms

The scaled path is implemented by transforming the assembled SI system:

- `A_hat = A_si / eps_ref`
- `b_hat = b_si / (eps_ref * V0)`

Then solve:
- `A_hat * psi_hat = b_hat`
- `psi = V0 * psi_hat`

This form is algebraically equivalent to writing the same box discretization in
scaled variables and keeps local multi-material variation because each edge
still carries its local `eps`, then is normalized by the same `eps_ref`.

## Source-term and boundary scaling coverage

In `unit_scaling` mode, config inputs are interpreted in common units and
normalized to SI before assembly. The scaled assembly then applies the global
`eps_ref` and `V0` normalization above.

Covered source terms:
- Region fixed charge (`fixed_charge_m3`)
- Interface sheet/fixed/trap charge (`sheet_charge_m2`, `fixed_charge_m2`,
  `trap_density_m2`, `trap_occupancy`)
- Neumann displacement boundaries (`normal_displacement_C_per_m2`)

Neumann note:
- The field name is stable for compatibility.
- Numeric interpretation is mode-dependent:
  - legacy mode: SI (`C/m^2`)
  - `unit_scaling` mode: common area units (`C/cm^2`), normalized to SI.

## PoissonSimulation path switch

`PoissonSimulation` behavior:
- legacy mode: construct `PoissonAssembler` without scaling spec and solve SI
  system as before.
- `unit_scaling` mode:
  - construct scaling references (`V0`, `eps_ref`) from mesh/material/doping
    context via `UnitScalingSystem`
  - pass `PoissonScalingSpec` to `PoissonAssembler`
  - solve scaled system
  - unscale potential to physical volts before output

Output behavior:
- `PoissonResult::potential` is physical volts.
- VTK `potential_V` is physical volts.
- Dimensionless `psi_hat` is internal only.

## Tests added for this feature

The following tests were added/updated in `tests/test_poisson.cpp`:
- scaled RHS maps back to SI RHS via `eps_ref * V0`
- scaled solve recovers legacy physical potential after unscaling
- `unit_scaling` Neumann displacement run matches legacy physical solution

Recommended focused test command:

```bash
ctest --test-dir build --output-on-failure -R "poisson|boundary|interface|scaling"
```

Full regression gate remains:

```bash
ctest --test-dir build --output-on-failure
```
