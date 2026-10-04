# SimpleMOS M21 production-SG incident-edge audit

## Scope

M21 audits every edge incident on the first free-node layer next to the drain.
It evaluates the 16 paired frozen states through 32 read-only production SG
probes and records the endpoint electron QF, variable-ni Bernoulli terms,
electron density, production mobility, finite-volume geometry, and the exact
mobility drive. All 16 configurations explicitly use
`high_field_gradient_discretization=transport_cell_vector`.

The edge-local three-factor game from M17/M18 is reused to separate mobility,
SG secant conductance, and QF log imbalance. No Sentaurus job is rerun, no
state is re-solved, and no production model default is changed.

## Execution and closure

| Check | Result | Acceptance |
|---|---:|---:|
| Frozen states | 16 | 16 |
| Production SG probes | 32 | 16 x 2 |
| Selected state/edge comparisons | 352 | nonzero |
| Stable SG versus production-flux closure | exactly `0` | `<1e-12` |
| Maximum edge-factor closure | `6.93889e-18` | `<1e-12` |
| Maximum M20 node-ledger closure | `6.86950e-16` | `<1e-12` |
| Maximum geometry difference | exactly `0` | `0` |
| New Sentaurus execution | no | no |
| Default model changed | no | no |

## Key result

For `n21`, `Vd=0.05 V`, `Vg=0.8 V`, the eleven edges incident on nodes
991/992 have the following absolute response allocation:

| Factor | Absolute share |
|---|---:|
| Mobility | `0.002978%` |
| SG secant conductance | `0.003445%` |
| QF log imbalance | `99.993576%` |

The same result is not restricted to the key state. Across all 16 states, the
QF-imbalance share over the selected first-layer edges ranges from `99.7974%`
to `99.9997%`.

At node 991 the edge-factor sum gives `2.704239e-7`, matching the M20 flux
residual difference. The QF factor alone gives `2.714746e-7`; mobility and SG
conductance contribute `-4.873e-10` and `-5.634e-10`. At node 992 the QF term
gives `6.723670e-8`, while the complete result is `6.708855e-8`.

## Leading contact/internal edge pairs

| Edge | Class | Vela phin drop | Imported phin drop | Drop change | Flux change |
|---:|---|---:|---:|---:|---:|
| 2107 | drain cut, node 991 | `2.22748e-10 V` | `1.91758e-10 V` | `-13.912%` | `-1.72301e-4` |
| 2094 | internal into node 991 | `1.74658e-9 V` | `1.50323e-9 V` | `-13.933%` | `-1.69862e-4` |
| 2114 | drain cut, node 992 | `2.07245e-11 V` | `1.81064e-11 V` | `-12.633%` | `-3.31664e-5` |
| 2088 | internal into node 992 | `1.65070e-10 V` | `1.44175e-10 V` | `-12.658%` | `-3.32320e-5` |

With the node-oriented sign, each internal-edge change opposes its drain-cut
partner. Their relative QF-drop reductions are extremely close but not exactly
equal. The small differential between those paired reductions is the M20 row
residual.

## Input-by-input comparison

For the four leading edges:

- the endpoint QF drop changes by `12.63-13.93%`;
- the production `transport_cell_vector` drive changes by `12.63-13.95%`;
- production mobility changes by only `3.59-4.81 ppm`;
- endpoint electron density changes by at most `4.99-6.54 ppm`;
- Bernoulli weights change by at most `0.51-1.61 ppm`;
- edge length and finite-volume `couple` do not change.

The vector drive follows the spatial QF change, but these fields are far below
the velocity-saturation regime: a roughly 13% drive change produces only a few
ppm mobility response. Consequently HFS mobility is not the mechanism carrying
the paired edge-flux difference. Density, Bernoulli weights, and geometry are
also too stable; their combined effect is represented by the `0.003445%` SG
conductance share.

The raw SG left/right terms are themselves cancellation-conditioned. On edge
2107 their conditions are `2.32e8` for the Vela state and `2.70e8` for the
imported state. The production factorized SG form nevertheless closes exactly,
so this is physical/numerical sensitivity of the QF imbalance rather than an
implementation loss of precision in the audited Vela edge flux.

## Decision

M21 narrows the discrepancy from “contact versus internal electron transport”
to the spatial distribution of electron quasi-Fermi potential. It rules out
finite-volume geometry, Bernoulli-weight changes, carrier-density changes, and
HFS mobility response as material contributors for the key edges.

This remains a state-difference attribution, not its upstream cause. The next
useful stage is to trace how the coupled solve establishes the slightly
different `phin` slopes on the paired paths 996-991-984 and 993-992-983. That
requires a local Newton/continuity response audit along those paths, retaining
the production contact elimination and `transport_cell_vector` configuration.

## Figures

- `m21_key_edge_flux_balance.png`: node-oriented incident-edge flux changes.
- `m21_key_edge_input_changes.png`: QF, drive, mobility, density, and Bernoulli changes.
- `m21_key_factor_shares.png`: exact edge-local response allocation.
