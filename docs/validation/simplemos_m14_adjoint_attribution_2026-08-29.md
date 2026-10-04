# SimpleMOS M14 terminal-current adjoint attribution

## Scope

M14 is an SDevice-only, read-only diagnostic stage. It does not change the
SimpleMOS default physical model and does not run Sentaurus. It separates:

1. the frozen-state direct drain-current response to a Vela-to-Sentaurus state
   substitution; and
2. the first-order self-consistent response obtained from
   `J^T lambda = dI_D/dx` on the production Vela Jacobian.

The terminal-error spectrum covers all 80 high-NWell, low-drain weak-region
points (`n21`-`n24`, `Vd=0.05 V`, `Vg=0..0.95 V`). Spatial attribution is
limited to the 16 states that have paired Sentaurus nodal fields and
self-consistent Vela states. No interpolation is used.

## Frozen configuration

| Item | M14 setting |
|---|---|
| State variables | `psi`, `phin`, `phip` |
| Terminal functional | drain contact current from the no-BC continuity residual |
| Adjoint equation | `J^T lambda = dI_D/dx` |
| Spatial partition | body: `x>0.25 um`; otherwise source: `y<-0.3 um`, channel: `-0.3<=y<=0.3 um`, drain: `y>0.3 um` |
| Physics | unchanged M8 production configuration |
| Sentaurus execution | none |

Carrier densities are derived from electrostatic and quasi-Fermi potentials;
they are not independent Newton unknowns. Mobility and SRH are therefore not
invented as state-vector coordinates: M14 links their frozen evidence from
M10-M12 and confines the new adjoint calculation to the actual coupled-DD
unknowns.

## Results

The 16 adjoint solves close to a maximum relative residual of
`2.306e-16`, well below the `1e-12` acceptance ceiling.

All 80 weak-region points have positive signed error, so Vela predicts a larger
current than Sentaurus throughout this subset. The mean signed error is
`+0.0510963 dex`; the maximum is `+0.109419 dex`. Per-device linear fits have
positive slopes (`0.0135` to `0.0294 dex/V`), but the `Vg=0.05 V` spikes show
that a linear platform is not a complete model of the residual.

For the representative high-NWell low-drain states:

| State | Frozen full-field change / Id | Adjoint relaxation / Id | Two-layer remainder / Id |
|---|---:|---:|---:|
| n21, Vg=0 V | -109.700% | +109.700% | -0.0000379% |
| n21, Vg=0.8 V | -13.6220% | +13.2751% | -0.346868% |

At `Vg=0.8 V`, the field-by-field replacement identifies `phin` as the direct
terminal-current control coordinate: its frozen response is `-13.6228%`, while
`psi` is `+0.000878%` and `phip` is negligible. This does not mean that
electrostatics is irrelevant to the self-consistent solution; it means that the
frozen drain-current functional is locally controlled by electron quasi-Fermi
transport.

The mutually exclusive spatial replacement at the same state gives:

| Region | Frozen direct / Id | Adjoint relaxation / Id | Remainder / Id |
|---|---:|---:|---:|
| Source | 0.000000% | +0.000000105% | +0.000000105% |
| Channel | -11.2697% | +10.9228% | -0.346868% |
| Drain | -2.35233% | +2.35233% | +0.000000342% |
| Body | 0.000000% | negligible | negligible |

Thus, the drain-local frozen contribution is almost exactly cancelled by the
self-consistent response. The non-cancelled first-order state response is a
channel effect. This answers the M10 disconnect: improving mobility on many
active edges need not change terminal current when those edge changes are not
on the self-consistent channel control path, or when the coupled state relaxes
against the frozen direct response.

## Interpretation

M14 supports the following bounded conclusions:

- The production terminal-current adjoint is numerically well resolved.
- Frozen current is directly most sensitive to `phin`, not to an isolated
  `psi` substitution.
- At the representative `Vg=0.8 V` state, self-consistency cancels 97.45% of
  the full-field frozen current response; the residual response is localized
  to the channel partition.
- A mapped Sentaurus state inserted into the Vela operator is not a Vela
  self-consistent solution. The operator/model difference is required to
  sustain the alternative state and cannot be inferred from state replacement
  alone.

M14 does not support claiming that `phin`, HFS, source barrier, or SRH is the
unique remaining cross-code root cause. Deep-off `Vg=0.05 V` relative responses
remain dominated by sub-fA conditioning and are retained in raw evidence but
excluded from headline quantitative attribution.

## Figures

- `m14_weak_region_80_points.png`: complete weak-region terminal-error spectrum.
- `m14_field_cancellation.png`: frozen direct versus adjoint response by state field.
- `m14_region_response_n21_vg0p8.png`: signed spatial response at the key state.
- `m14_adjoint_quality.png`: transpose-system closure for all 16 states.
