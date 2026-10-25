# PN2D BV Validation

Last updated: 2026-10-23

## Current decision

Template version 4 corrects the mobility impurity basis to `total_impurity`
(`Nd + Na`). Net doping does not describe impurity scattering at compensated
nodes. The user authorized this PN2D-specific correction; global mobility
defaults remain unchanged. The IV template already uses reconstructed total
impurity and remains unchanged.

Version 4 retains the atomic `element_edge_sg_gss_laux` profile on the
non-obtuse M0/M2 Tri3 mesh family. The historical full-curve acceptance below
belongs to version 3; it does not qualify version 4 automatically. Current
external-resistor controls and operator differences are recorded in the
[native matrix comparison](pn2d_external_resistor_native_matrix_2026-10-23.md).
The current investigation is restricted to zero bias. The
[staged Poisson-to-coupled comparison](pn2d_staged_zero_2026-10-23.md) gives
matching independent equilibrium fields with matched material inputs/constants.
The OldSlotboom electron DOS-mass conversion is now corrected; see the
[material conversion fix](sdevice_electron_dos_mass_fix_2026-10-23.md).
Finite-update/continuation qualification remains open; the equilibrium result
does not qualify a nonzero external-voltage sweep.
The [HFS and debug-precision follow-up](pn2d_zero_hfs_followup_2026-10-23.md)
aligns the finite RHS and one approximate Newton tangent, but later minority
carrier states, weak-row closure and applied-state acceptance remain unqualified.

The profile binds these settings as one unit:

- `impact_ionization.current_approximation=element_edge_sg_gss_laux`;
- `impact_ionization.source_mapping_mode=element_vertex_box_measure`;
- `impact_ionization.cell_reconstructed_midpoint_density=bernoulli`;
- `mesh_geometry.node_volume_policy=mixed_voronoi`; and
- `mesh_geometry.require_non_obtuse=true`.

The named `legacy_cell_reconstructed` profile is the complete rollback. It
restores cell-reconstructed current support, triangle GSS/QFP-gradient source
mapping, logistic midpoint density, barycentric node volumes, and no
non-obtuse qualification.

This is a template policy, not a global solver-default change. Configurations
that omit `mesh_geometry.node_volume_policy` still use barycentric volumes, and
the PN2D IV template remains unchanged with impact ionization disabled.

## Historical version-3 acceptance

The prospective atomic-default contract passed all M0, M2, forward-IV,
determinism, compatibility, geometry, and Release-regression gates. Both final
independent reviews returned `APPROVE`.

| Metric | M0 | M2 |
| --- | ---: | ---: |
| Exact reverse-bias lattice per branch/run | 29/29 | 29/29 |
| Avalanche-off current RMSE vs Sentaurus | 0.0000888 dex | 0.0000387 dex |
| Avalanche-on current RMSE | 0.001766 dex | 0.001923 dex |
| Avalanche-on maximum error | 0.003749 dex | 0.004647 dex |
| Knee RMSE | 0.002810 dex | 0.003051 dex |
| Absolute fitted breakdown-voltage error | 0.001 V | 0.001 V |
| Maximum avalanche-on continuity closure ratio | 0.002200 | 0.00000192 |

The two independent runs produced identical IV and state hashes. The
postprocess-only IIC branch produced exactly the same IV and states as the
corresponding avalanche-off branch at every requested bias, confirming that
postprocess-only impact ionization does not feed back into the continuity
equations.

The 201-point forward-IV guard also passed. At the declared Sentaurus anchors,
the mixed-Voronoi median relative error is 0.2700%, the maximum is 0.4066%, and
the maximum degradation relative to barycentric is `2.42e-14` absolute.

Direct Sentaurus box-measure evidence agrees with Vela mixed-Voronoi volumes to
floating-point precision on the actual non-obtuse M0/M2 meshes. The runtime
`require_non_obtuse=true` guard preserves this qualification boundary instead
of extending the result to arbitrary obtuse meshes.

The final acceptance Release suite passed 512/512 tests. The current `main`
also includes the CI portability fixes that declare Boost.Multiprecision,
remove the Linux `std::string`/`std::filesystem::path` conditional-expression
ambiguity, and keep `CoupledDDAssembler` test dependencies alive for the full
assembler lifetime.

## Model and branch semantics

The Sentaurus source deck remains
`reference_tcad/pn2d_sentaurus2018/source/pn2d_bv_sdevice.cmd` with
`Avalanche(VanOverstraeten)`. Vela's `van_overstraeten` implementation is the
comparison target.

The validation separates three branches:

1. avalanche off: establishes the transport and electrostatic baseline;
2. IIC/postprocess-only: evaluates impact ionization without adding it to the
   solved continuity residual or Jacobian; and
3. avalanche on/self-consistent: includes the source in both continuity
   equations and diagnostics.

Curve comparison, field/state comparison, source-process records, exact-state
hashes, continuity closure, and deterministic reruns are evaluated separately.
No acceptance result is inferred from a fitted scalar source factor.

## Qualification boundary

The accepted policy applies only to:

- the generated PN2D BV template version 3;
- the verified non-obtuse M0/M2 triangular mesh family;
- the documented nodal-doping, mobility, SG/GSS-Laux, and Van Overstraeten
  configuration; and
- the frozen reverse-bias and forward-IV contracts.

It does not establish:

- a global mixed-Voronoi or SG/GSS-Laux default;
- equivalence on arbitrary obtuse or non-Delaunay meshes;
- validation for other devices, materials, or contact models;
- a change to the PN2D IV template; or
- a calibrated commercial-grade breakdown prediction.

## Current sources of truth

- Template and profile implementation:
  `configs/templates/pn2d_bv.template.json` and
  `scripts/generate_pn2d_config.py`.
- Machine-readable configuration schema:
  `configs/schema/vela-simulation.schema.json`.
- Atomic acceptance contract:
  `docs/validation/contracts/pn2d_node_volume_policy_atomic_default_acceptance_v2.json`.
- Final acceptance report:
  `docs/validation/pn2d_node_volume_policy_atomic_default_acceptance_2026-08-01.md`.
- Final decision:
  `docs/validation/pn2d_node_volume_policy_atomic_default_decision_2026-08-01.md`.
- Independent reviews:
  `docs/validation/pn2d_node_volume_policy_atomic_default_independent_scientific_review_2026-08-01.md`
  and
  `docs/validation/pn2d_node_volume_policy_atomic_default_independent_code_review_2026-08-01.md`.
- Direct box-measure export:
  `docs/validation/sentaurus_box_measure_direct_export_2026-08-01.md`.

Earlier dated reports remain useful evidence for the diagnostic path that led
to this decision. Statements in those reports such as "production defaults
unchanged", "SG/Laux remains opt-in", or a recommended next experiment are
point-in-time conclusions and are superseded by the final atomic acceptance
and decision listed above.

## Reproduction entry points

Render the qualified default configuration:

```bash
python scripts/generate_pn2d_config.py \
  --template pn2d_bv \
  --output build/pn2d_bv.json
```

Run the generated deck and the repository regression suite:

```bash
build/vela_example_runner --config build/pn2d_bv.json
ctest --test-dir build --output-on-failure
```

Generated simulation and comparison outputs belong under `build/` or
`build-release/` and are not committed unless a task explicitly promotes a
small, stable validation fixture or contract.
