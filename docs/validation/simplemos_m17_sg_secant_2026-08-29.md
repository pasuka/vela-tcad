# SimpleMOS M17 stable SG secant attribution

## Scope

M17 is an SDevice-only, read-only diagnostic stage. It replaces the
cancellation-dominated independent Bernoulli/population view from M16 with an
exact, stable reparameterization of the same Vela electron SG transport
operator. No Sentaurus job is run and no production model default is changed.

The state matrix remains `n17/n21`, `Vd=0.05/1 V`, and
`Vg=0/0.05/0.8/2.5 V`: 16 paired self-consistent Vela and mapped Sentaurus
states, with no interpolation.

The M18 edge-closure audit exposed call-order dependence in standalone Enormal
diagnostics. M17 was regenerated after the diagnostic path was made to
initialize the same interface geometry as the production residual. This
changes the allocation among the three factors but not their summed production
electron-transport response; no physical-model default is changed.

## Stable factorization

For one oriented edge, define the two positive SG state terms

`L = B(-eta) n0` and `R = B(eta) n1`.

The production electron flux can then be written exactly as

`F = mobility * K_secant * xi`,

where `xi = log(L/R)` and

`K_secant = Vt * edge_geometry * (L-R)/(log(L)-log(R))`.

`K_secant` is a positive logarithmic-mean conductance and remains finite as
`L -> R`. In the SimpleMOS variable-ni OldSlotboom configuration, and away
from the production exponent clamp, `xi=(phin1-phin0)/Vt`. M17 evaluates the
logarithmic mean in `cpp_dec_float_100`, evaluates all eight vertices of the
three-factor game, and applies exact node-level Shapley allocation with the
M14 production terminal-current adjoint.

| Factor | Physical content |
|---|---|
| Mobility | Complete Vela PhuMob + Enormal + HFS edge mobility |
| SG secant conductance | Coupled electrostatic, variable-ni, population and quasi-Fermi common-mode state weight |
| QF log imbalance | Signed edge quasi-Fermi imbalance that drives the electron SG flux |

## Closure

| Check | Result | Acceptance |
|---|---:|---:|
| State count | 16 | 16 |
| Probe variants per state | 8 | 8 |
| Maximum endpoint closure | `0 A/um` | exact production endpoints |
| Maximum absolute closure to M15 | `5.20417e-18 A/um` | `<1e-12 A/um` |
| Maximum relative closure, all states | `1.79734e-15` | `<1e-4` |
| Maximum relative closure, `Vg>=0.8 V` | `1.78639e-15` | `<1e-8` |

The stable factorization reduces the key-state absolute-factor amplification
from `461.057x` in the ungrouped M16 split to `1.16571x`, a `395.52x`
reduction. The M16 numerical cancellation is therefore closed rather than
merely hidden by grouping.

## Key high-NWell, low-drain state

At `n21`, `Vd=0.05 V`, `Vg=0.8 V`:

| Stable factor | Signed response / Id | Absolute-magnitude share |
|---|---:|---:|
| Mobility | `-0.310851%` | `2.09983%` |
| SG secant conductance | `-0.741354%` | `5.00792%` |
| QF log imbalance | `+13.751431%` | `92.89225%` |
| Sum | `+12.699226%` | — |

The spatial sum remains `+12.949594%` in the channel and `-0.250369%` in the
drain, with negligible source and body terms. Inside the channel, the QF
log-imbalance contribution is positive while mobility and conductance oppose
it. The residual at this representative error state is therefore controlled
by the quasi-Fermi imbalance, not by HFS mobility drive alone and not by the
positive SG conductance state.

## Bias-regime transition

The dominant factor is bias dependent:

- At `Vg=0` and `0.05 V`, QF log imbalance dominates every paired state.
- At `Vg=0.8 V`, QF log imbalance contributes `53.22%` for n17/0.05 V,
  `92.89%` for n21/0.05 V, and `82.44%` for n21/1 V; n17/1 V instead has a
  `78.18%` conductance share.
- At `Vg=2.5 V`, SG secant conductance dominates all four states with
  `75.16%` to `93.04%` of absolute contribution.
- Mobility remains secondary across all eight `Vg>=0.8 V` states, with
  `0.738%` to `16.074%` of absolute contribution.

This is a useful regime separation: the high-NWell weak-current discrepancy
points to the quasi-Fermi imbalance/continuity state, while the fully on-state
response is chiefly a conductance-state effect.

## Interpretation boundary

M17 supports prioritizing edge quasi-Fermi imbalance and its channel-localized
continuity-state construction for the n21 low-drain error region. It does not
identify whether that imbalance originates from boundary values, spatial
discretization, recombination balance, or self-consistent feedback. The
secant conductance is also a coupled factor and must not be relabeled as a
single electrostatic or carrier-density model.

The factors are Vela operator diagnostics, not internal Sentaurus residual
terms. Deep-off sub-fA relative responses remain conditioning-sensitive and
are retained as diagnostics rather than headline quantitative root-cause
evidence.

## Figures

- `m17_factor_response_n21_vg0p8.png`: signed key-state factors.
- `m17_factor_fraction_16_states.png`: factor shares for all 16 states.
- `m17_cancellation_reduction_n21_vg0p8.png`: M16-to-M17 cancellation reduction.
- `m17_region_response_n21_vg0p8.png`: signed response by M14 region.
- `m17_nodal_qf_imbalance_n21_vg0p8.png`: spatial QF-imbalance support.
