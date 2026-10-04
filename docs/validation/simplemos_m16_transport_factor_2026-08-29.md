# SimpleMOS M16 electron-transport factor attribution

## Scope

M16 is an SDevice-only diagnostic stage. It decomposes the M15 electron
transport response into mobility-state, mobility-drive, Bernoulli-weight, and
carrier-population interventions. No Sentaurus job is run and no production
default is changed.

The matrix contains the same 16 paired states used by M10-M15: `n17` and
`n21`, `Vd=0.05/1 V`, and `Vg=0/0.05/0.8/2.5 V`. M16 adds a diagnostic-only
production-equation probe because the older standalone `sg_edge_flux_probe`
did not reproduce the solved continuity-equation node terms closely enough in
strong inversion. Results from that rejected preliminary path are not frozen.

The M18 closure audit later found that standalone Enormal diagnostics depended
on whether a production residual had already initialized interface geometry.
M16 was regenerated after making that diagnostic initialization explicit and
call-order independent. The production residual and physical-model defaults
were not changed.

## Factorization

For each node, the new probe evaluates all 16 vertices of a four-factor game:

| Factor | Baseline to replacement intervention |
|---|---|
| Mobility state | state inputs used by Vela PhuMob/Enormal/HFS mobility evaluation |
| Mobility drive | Vela `transport_cell_vector` electron GradQuasiFermi field used by HFS |
| Bernoulli weights | variable-ni SG argument and Bernoulli weights |
| Carrier population | electron densities at the edge endpoints |

The empty and full vertices are independently evaluated solved-equation
electron-flux terms. The 14 mixed vertices use `cpp_dec_float_100` for the
Bernoulli-density product. Four-factor Shapley values are calculated per node
and weighted by the M14 terminal-current adjoint.

Bernoulli weights and carrier populations are physically and algebraically
coupled through the Boltzmann state relation. Their independent interventions
therefore generate a large, condition-sensitive cancellation. M16 also freezes
a three-factor grouped allocation in which they move together as the
`SG state kernel`. The grouped allocation is the primary interpretation view;
the four-factor view is retained to expose the cancellation rather than hide
it.

## Closure

| Check | Result | Acceptance |
|---|---:|---:|
| State count | 16 | 16 |
| Probe variants per state | 16 | 16 |
| Maximum endpoint closure | `8.67362e-19 A/um` | floating-point roundoff |
| Maximum absolute closure to M15 | `9.86732e-16 A/um` | `<1e-12 A/um` |
| Maximum relative closure, all states | `3.17681e-5` | `<1e-4` |
| Maximum relative closure, `Vg>=0.8 V` | `2.78340e-10` | `<1e-8` |

The all-state relative maximum is set by cancellation of much larger Shapley
terms at small transport responses. The absolute closure remains two orders of
magnitude below its acceptance ceiling.

## Key high-NWell, low-drain state

At `n21`, `Vd=0.05 V`, `Vg=0.8 V`, the grouped factors relative to Vela
baseline drain current are:

| Grouped factor | Signed response / Id | Absolute-magnitude share |
|---|---:|---:|
| SG state kernel | `+13.010104%` | `97.666256%` |
| Mobility state | `-0.052100%` | `0.391113%` |
| Mobility drive | `-0.258778%` | `1.942632%` |
| Sum | `+12.699226%` | — |

The sum reproduces the M15 electron-transport response. The spatial total is
`+12.949594%` in the channel and `-0.250369%` in the drain; source and body are
negligible.

The ungrouped four-factor result is intentionally not used as a root-cause
ranking: Bernoulli weights contribute `-2920.873379%` while carrier population
contributes `+2933.883659%`. Their nearly complete cancellation is the
mathematical consequence of breaking the normally coupled SG state relation.

Across all eight strong states (`Vg>=0.8 V`), the SG state kernel accounts for
`80.50%` to `98.68%` of grouped absolute contribution. Mobility drive accounts
for `0.214%` to `3.447%`. This moves the leading hypothesis away from the HFS
driving-field construction alone and toward the coupled SG state kernel:
electrostatic/variable-ni weights, carrier population, and quasi-Fermi state
consistency.

## Interpretation boundary

M16 supports:

- the SG state kernel as the dominant grouped Vela transport response in all
  strong states;
- a small direct mobility-drive share at the representative error state; and
- the M15 conclusion that the response is localized principally in the
  channel, with smaller opposing drain support.

M16 does not support:

- treating Bernoulli and carrier-population standalone values as independent
  physical error sources;
- claiming the Vela SG kernel is identical to an unavailable internal
  Sentaurus residual decomposition;
- claiming HFS has zero effect, because mobility-state/HFS interactions remain
  grouped inside the Vela mobility operator; or
- changing a default model from these diagnostics.

## Figures

- `m16_grouped_factor_response_n21_vg0p8.png`: primary grouped signed response.
- `m16_four_factor_cancellation_n21_vg0p8.png`: independent Bernoulli/population cancellation.
- `m16_grouped_factor_fraction_16_states.png`: grouped absolute shares across all states.
- `m16_region_response_n21_vg0p8.png`: signed response by M14 region.
- `m16_nodal_sg_kernel_response_n21_vg0p8.png`: nodal SG-kernel spatial support.
