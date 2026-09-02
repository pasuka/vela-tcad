# Templates/LDMOS G3 material-local Poisson charge-volume qualification (2026-09-02)

## Outcome

The explicit, default-off `material_local` Poisson charge-volume policy closes
the remaining classical G3 maximum-gm failure without changing mobility,
flatband voltage, threshold parameters, HFS, carrier SG couples, or solver
acceptance thresholds. The exact 31-point, no-predictor curve passes all six
frozen Stage-3 gates.

The accepted diagnostic contract is:

- `contact_boundary_reconstruction = legacy_node_local`;
- carrier transport uses the qualified 16,237-edge external AverageBox profile;
- `discretization.poisson_charge_volume_policy = material_local`;
- predictor, IALMob, quantum potential, avalanche, and thermal coupling are off;
- Poisson carrier and ionized-dopant charge use the barycentric volume of
  adjacent transport cells only.

The global default remains `global`. Continuity sources, fixed/interface
charge, BTBT, SRH/Auger, impact ionization, stored charge, and transport SG
couples are unchanged.

The accepted combination is frozen as
`templates_ldmos_exact_topology_phase_a_classical_v2`; the versioned contract
selects both the archived external AverageBox carrier couple and the
material-local Poisson charge volume explicitly.

## P0: read-only evidence and frozen masks

The read-only audit reproduces the existing Vela Poisson operator before
changing only the charge volume. Input SHA256 values, interface/channel masks,
three channel windows, and output schema are frozen in
`vela.templates_ldmos_g3_interface_charge_volume_audit.v1`.

Six complete matched Sentaurus/Vela state pairs were available. The requested
`Vg=0` and `0.833333 V` pairs were not complete, so coverage is reported as
incomplete rather than silently filled or interpolated. This does not affect
the two endpoint implementation gate or the subsequent exact 31-point curve.

Frozen provenance:

| Object | SHA-256 |
|---|---|
| exact-topology mesh | `36f79a06a21550d69fb58ad221457b1c1ccc63fa3a929046ef6b90e185f09ff6` |
| nodal doping | `b9924f72bf36386d031e7a168adbc997d9f076db1662b7b8c1e3416ad68acb39` |
| Sentaurus Measure-derived node oracle | `8a6d78d9d129ac7905697013cd0da23014eb7bb33e682568f41572157f72c896` |
| free interface-node mask | `bcebfcbeb71a21baa5dcd274f7c6785354f27ee1b177a1b54b4880cef46b0902` |
| channel-interface mask | `8b88e746780ec0f35d79372d1d38c80f433412e0887f491a7e743798232546e6` |
| external AverageBox transport profile | `758161f73f570511ccd914628451ebcad58b5f8e58dd6b6351f00669e18b36cb` |

At `Vg=1.0/1.166667 V`, replacing global node volume with transport-cell
barycentric volume reduces the normalized channel-interface Poisson residual:

| Bias | Median reduction | P95 reduction | Maximum reduction |
|---:|---:|---:|---:|
| 1.0 V | 93.73% | 91.74% | 96.53% |
| 1.166667 V | 93.69% | 92.28% | 97.18% |

The three pre-registered channel windows independently show that the previous
global volume included oxide-side area at shared Si/SiO2 vertices. The free
interface median ratio `transport_volume/global_volume` is `0.49007`.

## P1: implementation contract

The new policy is selected under `discretization`, not `mesh_geometry`, because
it changes material ownership of Poisson charge without changing pure geometry.
The implementation:

1. accumulates `triangle_area/3` over transport cells;
2. applies the result to Poisson mobile-carrier and dopant residual/Jacobian
   terms together;
3. rejects composition with `mixed_voronoi` until independently qualified;
4. returns the existing node-volume vector for all-transport meshes, preserving
   the single-material path bitwise;
5. leaves the default policy and all unrelated source mappings unchanged.

## P2: endpoint reclose and exact 31-point qualification

Both Sentaurus-seeded endpoints converge in five Newton iterations in Debug and
Release, with identical currents and residuals. Relative KCL closure is of
order `1e-12`.

| Metric | Previous global volume | Material-local volume | Sentaurus |
|---|---:|---:|---:|
| Id at 1.0 V (A/um) | 1.0217975e-6 | 1.3722485e-6 | 1.2541552e-6 |
| Id at 1.166667 V (A/um) | 2.2130153e-6 | 2.9216905e-6 | 2.8279161e-6 |
| two-point gm / Sentaurus | 0.756924 | 0.984547 | 1.0 |

The production-style Release sweep uses all 31 exact reference biases, no
interpolation, no predictor, and no unscored warm-start points:

| Frozen metric | Result | Limit | Status |
|---|---:|---:|:---:|
| Exact resolved points | 31/31 | 31/31 | pass |
| Median absolute log-current error | 0.004492 dex | 0.10 dex | pass |
| P95 absolute log-current error | 0.099695 dex | 0.20 dex | pass |
| Strong-inversion endpoint error | 0.954% | 20% | pass |
| Diagnostic Vth error | 10.754 mV | 100 mV | pass |
| Maximum-gm error | 1.545% | 20% | pass |
| Worst resolved KCL error | 9.180e-10 | 1% | pass |

Maximum gm remains on the same `1.0--1.166667 V` segment in both engines.
Stage-3 classical G3 therefore passes without parameter calibration.

## P3: Poisson couple, cell permittivity, and Sentaurus Measure

With material-local charge volume held fixed, the independent region-local
AverageBox Poisson-couple A/B reconstructs 30,022 edges, including 760
heteromaterial edges, as

```text
sum_cell(eps_r(cell) * local_AverageBox_couple(cell))
```

Both endpoints converge, but the two-point gm ratio changes from `0.984547` to
`0.983021`. The per-cell permittivity/couple candidate therefore does not
improve the accepted charge-volume result and remains an explicit,
default-off diagnostic. Sentaurus `Measure` is used as the external volume
oracle; it is not copied into the production mesh or treated as a tunable
parameter.

## P4: regressions

- `test_mos_solver_crosscheck`: 458 assertions passed;
- `test_mos_mixed_material`: 1,452 assertions passed, including shared
  Si/SiO2 volume exclusion, unchanged continuity rows, and the all-silicon
  bitwise invariant;
- `test_mos_interface_charge`: 48 assertions passed;
- `test_ldmos_dd`: 357 assertions passed;
- `test_newton_solver` and `test_dc_sweep`: 4,936 assertions passed;
- relevant Python contract/runner/audit suites: 36 tests passed;
- full Debug build and the earlier Poisson/Newton/DC-sweep regression set pass.

## Decision and P5 handoff

1. The material-local policy is qualified for the Templates/LDMOS classical G3
   contract but remains default-off until the template contract explicitly
   selects it.
2. The former H2 Fermi/BGN SG-placement branch is closed as unnecessary for the
   maximum-gm failure: the frozen curve gate now passes under the unchanged SG
   physics contract. Its historical artifacts remain diagnostic evidence.
3. The known-difference ledger is not used to absorb the former gm error.
4. The next Phase-A work is the planned Stage-4 isothermal Id-Vd and WP3-D
   Sentaurus `hRecVelocity` influence decision. IALMob/hQP remain decision-gated;
   avalanche and self-heating remain Phase-B work.

## Reproducibility artifacts

- P0 audit: `scripts/audit_templates_ldmos_g3_interface_charge_volume.py`;
- P0 schema: `schemas/vela.templates_ldmos_g3_interface_charge_volume_audit.v1.schema.json`;
- endpoint runner: `scripts/run_templates_ldmos_g3_same_bias_reclose.py`;
- curve runner: `scripts/run_templates_ldmos_g3_averagebox_curve.py`;
- exact scorer: `scripts/analyze_templates_ldmos_g3_idvg.py`;
- P3 A/B: `scripts/audit_templates_ldmos_g3_side_local_poisson_ab.py`;
- ignored P0 evidence:
  `reference_staging/templates_ldmos_g3_interface_charge_volume_p0_20260902/`;
- ignored endpoint evidence:
  `reference_staging/templates_ldmos_g3_material_local_charge_reclose_corrected_20260902/`;
- ignored 31-point evidence:
  `reference_staging/templates_ldmos_g3_material_local_charge_curve_corrected_20260902/`;
- ignored P3 evidence:
  `reference_staging/templates_ldmos_g3_material_local_charge_side_local_poisson_ab_corrected_20260902/`.
