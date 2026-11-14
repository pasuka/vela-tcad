# PN2D BV Validation

Last updated: 2026-11-14

## Current PN2D regression and JSON defaults

CTest uses the GRD/DAT reference suite (`pn2d_ci_*`). The default IV v2 and
BV v5 JSON templates use the suite's compensated potential and concentration
coordinates, unit scaling, and Newton-controlled voltage steps. IV sweeps
to +10 V with 201 requested biases. BV sweeps to −20 V with 401 requested
biases and uses the atomic `element_edge_sg_gss_laux` profile. The standalone
`pn2d_ci_0v` equilibrium test remains part of this suite.

Use `--set solver_profile=legacy` with the template generator to reproduce
the earlier numerical defaults. The BV `legacy_cell_reconstructed` profile
also selects the legacy solver by default. Historical evidence below belongs
to the input and template versions recorded in those reports.

## Version 4 mobility and mesh policy

Template version 4 corrects the mobility impurity basis to `total_impurity`
(`Nd + Na`). Net doping does not describe impurity scattering at compensated
nodes. At the time, IV v1 used reconstructed total impurity.

Version 4 retains the atomic `element_edge_sg_gss_laux` profile on the
non-obtuse M0/M2 Tri3 mesh family. The historical full-curve acceptance below
belongs to version 3; it does not qualify version 4 automatically. Current
external-resistor controls and operator differences are recorded in the
[native matrix comparison](pn2d_external_resistor_native_matrix_2026-10-23.md).
The October external-resistor investigation was restricted to zero bias. The
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

## Sentaurus breakdown knee acceptance

The PN2D BV Sentaurus comparison passes the fitted breakdown knee check when
`abs(Vela knee - Sentaurus knee) / abs(Sentaurus knee)` is at most 0.01
(1%) and both curves have the same slope-crossing outcome. The report includes
the absolute voltage difference in volts and the dimensionless relative error.

## Compensated population frontend qualification

The PN2D SDevice frontend admits the explicit compensated population sidecar
on the original v2 non-obtuse mesh (1,943 nodes, 3,680 triangles), at 300 K.
Stage the [sidecar](../../configs/pn2d_bv_compensated.vela.json) beside
the BV deck as `pn2d_bv_sdevice.vela.json`, then run
`vela device pn2d_bv_sdevice.cmd --pn2d-population`. The acceptance contract pins the sidecar
and both emitted solver, geometry and contact policies. It reuses the
original reference hashes, ten bias anchors and field limits.

The Release frontend run reached −20 V with 458 accepted sweep points and
all 401 requested output biases. The maximum anode-current error was 0.199%,
maximum relative terminal imbalance was 0.000279%, and the fitted knee was
−19.404 V versus −19.405 V in Sentaurus. All ten reverse-bias field anchors
and the equilibrium field comparison passed.

The current reference regression is `pn2d_ci_bv` in the
[GRD/DAT suite](../../tests/regression/README.md#pn2d-release-regression).
Its reference archive supplies the pinned original mesh.

## Standalone equilibrium deck pending removal

The standalone reference deck
[`pn2d_0v_sdevice.cmd`](../../reference_tcad/pn2d_sentaurus2018/source/pn2d_0v_sdevice.cmd)
is marked for planned removal. The IV and BV decks both solve at 0 V
and their end-to-end acceptance compares that equilibrium field state.
Before removal, carry the standalone absolute-current and current-balance
checks into the sweep equilibrium acceptance. The standalone deck, archived
baseline, and CTest registration remain available during this transition.

## Vela-generated mesh comparison

The offline SDE path generated a non-obtuse PN2D mesh with 1,482 nodes and
2,800 triangles. Vela read its GRD/DAT directly. Sentaurus Device read a
TDX conversion of those same files. The comparison checked coordinate
bijection, triangle connectivity, contact node sets, and nodal donor and
acceptor densities before comparing fields.

- Standalone 0 V passed the equilibrium field limits. It exceeded the
  equilibrium current floor with ordinary population coordinates and
  passed with the compensated state enabled.
- IV reached +10 V with all 201 requested output biases. The maximum
  current error in the acceptance window (0.2–0.3 V) was 6.08%. The
  enforced terminal-imbalance maximum was 0.155%. Both field checks passed.
- Compensated BV reached −20 V with all 401 requested output biases.
  The maximum current error was 2.56%, terminal imbalance was 0.000165%,
  and the fitted knees were −19.371 V (Vela) and −19.381 V (Sentaurus).
  Those gates passed. All ten biased field checks exceeded the strict
  limits used for the original mesh. At −20 V the potential RMS error
  was 41.4 mV, the electron and hole density log RMS errors were 0.148
  and 0.0941 dex, and the global electric-field vector error was 2.18%.

A fixed-state replay traced the biased field discrepancy to 383 nominally
zero cotangents rounded slightly below zero. The geometry builder replaced
each with a positive barycentric edge contribution. With the Sentaurus
state frozen at −20 V, omitting these contributions reduced the maximum
Poisson residual from 2.65e−11 to 4.76e−24 C/m. A separate full self-consistent
run with `mesh_geometry.fallback_negative_cotangent=false` passed all ten
field anchors and reached −20 V. Its maximum current error was 0.235%, and
the fitted knee agreed at −19.381 V. This control used the same non-obtuse
mesh and the complete mixed-Voronoi SG/GSS-Laux profile.

### Unified near-right cotangent tolerance

A Release rerun used compensated potential-density solves with UMFPACK and
the default negative-cotangent fallback enabled on both meshes. The SDevice
precheck and box assembly used the same cotangent threshold. On the
generated mesh, 763 nonzero contributions were set to zero, including the
383 previously treated as negative. The fallback count was zero on both
meshes. The generated mesh's control areas summed to
0.9999999999999907 µm² against a geometric area of 1 µm².

| Result versus same-mesh Sentaurus | Original v2 mesh | Vela-generated mesh |
| --- | ---: | ---: |
| IV output biases through +10 V | 201/201 | 201/201 |
| IV maximum current error at 0.2–0.3 V | 0.00193% | 0.00933% |
| BV output biases through −20 V | 401/401 | 401/401 |
| BV maximum current error in the acceptance window | 0.199% | 0.235% |
| BV maximum relative terminal imbalance | 0.000279% | 0.000151% |
| BV fitted knee, Vela / Sentaurus | −19.404 / −19.405 V | −19.381 / −19.381 V |
| BV reverse-bias spatial anchors passed | 10/10 | 10/10 |

The equilibrium and +0.2 V IV field comparisons passed on both meshes.
At −20 V on the generated mesh, the potential RMS difference was
4.0e−15 V. The electron and hole density log RMS differences were
0.000168 and 0.000431 dex, and the global and junction electric-field vector
differences were 1.04% and 0.0891%. The formal compensated BV acceptance
also passed using the original v2 reference archive.

The formal frontend qualification continues to cover the pinned original
mesh. The generated-mesh comparisons use 300 K, the same material parameters,
potential in volts, concentrations in cm⁻³, and terminal currents in A/µm.

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

The historical version-3 acceptance applies only to:

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
