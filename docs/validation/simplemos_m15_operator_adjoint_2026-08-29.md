# SimpleMOS M15 operator-residual adjoint attribution

## Scope

M15 is an SDevice-only, read-only diagnostic stage. It evaluates the Vela
equation residual at the paired mapped Sentaurus state and weights the change
with the production terminal-current adjoint from M14. It does not run
Sentaurus, change C++, or modify a default physical model.

The 16-state matrix is `n17` and `n21`, `Vd=0.05/1 V`, and
`Vg=0/0.05/0.8/2.5 V`. These are the states for which both mapped Sentaurus
nodal fields and converged Vela states are already frozen. No interpolation is
used.

## Method

For every state, M15 forms

```text
delta F = F_vela(x_sentaurus_mapped) - F_vela(x_vela_converged)
delta I_operator = -lambda^T delta F
```

where `lambda` is the M14 solution of
`J^T lambda = dI_D/dx`. The residual difference is decomposed by equation term
and by the same source/channel/drain/body node partition used in M14. Two
`newton_carrier_term_probe` executions per state provide the solved-equation
terms at the Vela and mapped Sentaurus states.

This quantity is a first-order Vela operator response to the state mismatch.
It is not a Sentaurus residual decomposition because equivalent internal
Sentaurus equation terms are unavailable.

## Numerical closure

All 32 carrier-term probes completed. Summing all equation terms reconstructs
the M14 `field_all` adjoint relaxation in every state:

| Check | Result | Acceptance |
|---|---:|---:|
| State count | 16 | 16 |
| Carrier-term probes | 32 | 32 |
| Maximum absolute closure error | `8.61941e-18 A/um` | `<1e-12 A/um` |
| Maximum relative closure error | `4.25433e-14` | `<1e-8` |

The decomposition therefore adds no material numerical closure error.

## Representative high-NWell, low-drain state

At `n21`, `Vd=0.05 V`, `Vg=0.8 V`, the terminal signed error is
`+0.0639918 dex`. The adjoint-weighted operator response relative to the Vela
baseline current is:

| Equation component | Response / Id | Absolute-magnitude share |
|---|---:|---:|
| Electron transport | `+12.699226%` | `95.661739%` |
| Poisson | `+0.575910%` | `4.338261%` |
| Electron SRH | `+4.10e-9%` | `3.09e-8%` |
| Hole transport | `+2.63e-12%` | `1.98e-11%` |
| All boundary/gauge/impact terms | negligible | negligible |

The corresponding spatial decomposition is:

| Region | Response / Id |
|---|---:|
| Source | `+0.0000499%` |
| Channel | `+13.525400%` |
| Drain | `-0.250314%` |
| Body | `+5.60e-9%` |

Thus the non-cancelled response is principally an electron-continuity
transport residual in the channel. The drain region partially opposes it.
SRH is quantitatively excluded as a material contributor at this state.

At the `n21`, `Vd=0.05 V`, `Vg=0 V` off-state, electron transport contributes
`+109.699441%` of the baseline current response; channel and drain contribute
`+43.288412%` and `+66.411167%`, respectively. The `Vg=0.05 V` sub-fA states
remain ill-conditioned and are retained in the evidence but not used for
headline relative attribution.

## M11 mobility-factor cross-check

For the key `Vg=0.8 V` state, the M11 frozen log-current factorial main effects
are `-1.573149 dex` for PhuMob, `-3.65966e-5 dex` for Enormal, and
`-1.39004e-7 dex` for HFS. This cross-check shows that the M15 electron-
transport result cannot be narrowed to HFS alone. It identifies the broad
electron transport operator, which contains state gradients, mobility,
Scharfetter-Gummel flux construction, and their coupled response.

M15 therefore supports:

- electron transport, rather than SRH or a boundary term, as the dominant Vela
  operator reaction at the representative weak-inversion state;
- the channel as the dominant spatial support, with a smaller opposing drain
  response; and
- exact reconstruction of the M14 adjoint relaxation from solved-equation
  terms.

M15 does not support:

- claiming that HFS, PhuMob, electrostatics, or any private Sentaurus formula is
  the unique remaining cross-code root cause;
- interpreting the decomposition as a direct Sentaurus operator residual; or
- changing a production default based on this diagnostic.

## Figures

- `m15_component_response_n21_vg0p8.png`: signed equation-term response at the key state.
- `m15_region_response_n21_vg0p8.png`: signed spatial response at the key state.
- `m15_component_fraction_16_states.png`: absolute equation-family shares over all 16 states.
- `m15_m11_mobility_crosscheck_n21_vg0p8.png`: M11 frozen mobility-factor effects.
- `m15_nodal_response_n21_vg0p8.png`: nodal adjoint-weighted response map.
