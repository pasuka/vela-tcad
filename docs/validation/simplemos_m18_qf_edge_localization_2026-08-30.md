# SimpleMOS M18 quasi-Fermi edge localization

## Scope

M18 is an SDevice-only, read-only diagnostic stage. It localizes the M17
electron SG factors on native finite-volume edges, then replaces the M15
electron-transport entry with the three stable M17 factors to form a complete
self-consistent feedback ledger. No Sentaurus job is run and no production
model default is changed.

The matrix remains `n17/n21`, `Vd=0.05/1 V`, and
`Vg=0/0.05/0.8/2.5 V`: 16 paired Vela/Sentaurus states without interpolation.

## Edge response definition

Each edge evaluates the complete eight-vertex M17 game:

`electron flux = mobility * SG secant conductance * QF log imbalance`.

The three edge Shapley values are multiplied by the production terminal-current
adjoint weight `lambda(node1)-lambda(node0)`. Constrained electron-continuity
rows have zero lambda. Contact cuts use the same finite-volume support as
terminal-current extraction: exactly one endpoint belongs to the named contact
node set. The source, drain, substrate and gate cuts plus the four internal
spatial buckets are mutually exclusive and exhaustive.

During M18 implementation, the first probe exposed call-order dependence in
standalone Enormal diagnostics: `residual()` initialized cached interface
normals/distances, while direct edge and factor probes did not. The diagnostic
paths now initialize the same surface geometry themselves. This is a
diagnostic-consistency fix; the production residual and model defaults are
unchanged. M17 was regenerated with the corrected call-independent path before
M18 was frozen.

## Closure

| Check | Result | Acceptance |
|---|---:|---:|
| State count | 16 | 16 |
| Edge probes / variants | 16 / 8 | 16 / 8 |
| Maximum edge-to-node absolute closure | `6.93889e-18 A/um` | `<1e-12 A/um` |
| Maximum edge-to-node relative closure | `1.05898e-11` | `<1e-8` |
| Maximum feedback-ledger closure | `9.97466e-18 A/um` | `<1e-12 A/um` |
| Maximum per-edge Shapley flux closure | `4.44089e-16` | floating-point roundoff |

## Key high-NWell, low-drain state

At `n21`, `Vd=0.05 V`, `Vg=0.8 V`, the M17 electron-transport response is:

| Stable factor | Signed response / Id | Absolute-factor share |
|---|---:|---:|
| Mobility | `-0.310851%` | `2.09983%` |
| SG secant conductance | `-0.741354%` | `5.00792%` |
| QF log imbalance | `+13.751431%` | `92.89225%` |
| Electron transport sum | `+12.699226%` | — |

The QF edge response is not primarily an undifferentiated internal-channel
effect:

| Edge partition | Signed response / Id | Absolute QF support |
|---|---:|---:|
| Drain contact cut | `+13.622821%` | `80.63628%` |
| Internal channel | `+0.128618%` | `19.36367%` |
| All other buckets | approximately zero | `<0.00005%` combined |

One drain-cut edge (`edge 2107`) supplies `66.59448%` of absolute QF support;
the two leading drain-cut edges (`2107`, `2114`) supply `79.41330%`. The leading
internal-channel edge is rank 3 and brings the cumulative support to
`80.92725%`.

This falsifies the strong form of the prior hypothesis that M17's QF response
is concentrated on low-terminal-sensitivity internal edges. For this state,
the controlling adjoint-weighted imbalance is localized mainly at the drain
finite-volume cut, with a smaller distributed internal-channel remainder.

## Self-consistent feedback ledger

| Component | Signed response / Id |
|---|---:|
| QF log imbalance | `+13.751431%` |
| Poisson feedback | `+0.575910%` |
| SG secant conductance | `-0.741354%` |
| Mobility | `-0.310851%` |
| Direct SRH | `+4.10e-9%` |
| Boundary/gauge plus other-carrier terms | negligible |
| Complete M15 response | `+13.275136%` |

Direct SRH is therefore not a material operator response at the key state.
This does not prove that SRH had no historical influence on the converged
carrier state; M18 only measures first-order response around that state.

Across the 16 states, contact cuts carry `55.0-80.6%` of absolute QF support at
`Vg=0.8 V` and less at the fully on state for three of four device/drain
combinations. The near-100% off-state fractions are retained for completeness
but are not headline evidence because sub-fA responses are ill-conditioned.

## Interpretation boundary

M18 localizes the Vela operator response; it is not a native Sentaurus residual
decomposition and does not identify a private Sentaurus formula. The drain-cut
result narrows the next investigation to the drain-side electron continuity
boundary/cut state, its quasi-Fermi interpolation and its self-consistent
coupling to the internal channel. It does not justify changing HFS, Enormal,
SRH, or any other default model.

## Figures

- `m18_qf_bucket_support_n21_vg0p8.png`: absolute QF support by FV partition.
- `m18_qf_edge_map_n21_vg0p8.png`: signed per-edge spatial response.
- `m18_qf_edge_pareto_n21_vg0p8.png`: response concentration by edge rank.
- `m18_feedback_ledger_n21_vg0p8.png`: material self-consistent components.
- `m18_qf_contact_fraction_16_states.png`: contact-cut share across all states.
