# SimpleMOS M19 drain-contact cut audit

## Scope

M19 is an SDevice-only, read-only audit of the drain finite-volume contact cut.
It reuses the 16 exact-node paired states from M10/M14/M18 and runs nine Vela
SG probes per state: the eight vertices of a three-factor frozen endpoint-state
game plus one complete imported-Sentaurus-state control. No Sentaurus job is
run, no state is re-solved, and no production model default is changed.

The three frozen factors are:

1. all state fields on the constrained drain-contact nodes;
2. electron quasi-Fermi potential on the adjacent free interior nodes;
3. the remaining state fields on those adjacent free nodes.

The contact cut is exactly the terminal-current extraction support: an edge
with one and only one endpoint in the drain contact node set.

## Closure and qualification

| Check | Result | Acceptance |
|---|---:|---:|
| Paired state count | 16 | 16 |
| SG probes | 144 | 16 x 9 |
| Maximum drain-contact `phin` bias error | `6.66134e-16 V` | `<1e-12 V` |
| Maximum strong-state SG cut-current reconstruction error | `1.11082e-7` | `<1e-6` |
| Maximum frozen-game Shapley closure | `2.11758e-22 A/um` | `<1e-18 A/um` |
| New Sentaurus execution | no | no |
| Default model changed | no | no |

The maximum all-state current reconstruction discrepancy is retained as a
guarded diagnostic, not an acceptance metric. At `Vg=0/0.05 V`, mapped-state
sub-fA current is dominated by the previously documented cancellation and
state-serialization conditioning. M19 therefore qualifies SG cut-current
closure only for `Vg>=0.8 V`.

## Contact boundary condition

Vela and the imported Sentaurus states impose the same electron quasi-Fermi
value at every drain contact node. Across all 16 states, the maximum difference
from the configured drain bias is `6.66134e-16 V`.

At the key state `n21`, `Vd=0.05 V`, `Vg=0.8 V`, all three drain nodes have
`phin=0.05 V` in both states to floating-point precision. The terminal gap is
therefore not caused by a wrong applied drain voltage or a different contact
Dirichlet value.

## Unchanged-SG state replay

For the four stable `Vg=0.8 V` states, replaying the full imported Sentaurus
state through the unchanged Vela SG operator reduces the terminal error to at
most `0.0004435 dex`. The imported state closes `97.54-99.43%` of the original
Vela-to-Sentaurus current gap.

| State | Self-consistent error | Imported-state Vela-SG error | Gap closed |
|---|---:|---:|---:|
| n17, Vd=0.05 V, Vg=0.8 V | `0.0248444 dex` | `0.0003981 dex` | `98.44%` |
| n17, Vd=1 V, Vg=0.8 V | `0.0159426 dex` | `0.0003994 dex` | `97.54%` |
| n21, Vd=0.05 V, Vg=0.8 V | `0.0639918 dex` | `0.0003949 dex` | `99.43%` |
| n21, Vd=1 V, Vg=0.8 V | `0.0573070 dex` | `0.0003945 dex` | `99.36%` |

This is direct evidence that the Vela drain-cut SG flux construction can
reproduce the Sentaurus terminal current when supplied with the Sentaurus
state. It does not prove that every internal Sentaurus SG face formula is
identical; it specifically qualifies the observed drain terminal cut.

## Key-state frozen endpoint attribution

At `n21`, `Vd=0.05 V`, `Vg=0.8 V`:

| Frozen factor | Exact current change / Vela Id | Absolute-factor share |
|---|---:|---:|
| Drain contact endpoint state | `+0.000182%` | `0.00134%` |
| Adjacent free-node electron QF | `-13.622715%` | `99.99478%` |
| Adjacent free-node other state | `+0.000529%` | `0.00388%` |

The constrained contact endpoint is already correct. Nearly the complete
frozen drain-cut current correction comes from the electron quasi-Fermi value
on the five adjacent free interior nodes.

The two leading edges are:

| Edge | Vela `dphin` | Sentaurus `dphin` | Frozen adjacent-QF current change |
|---|---:|---:|---:|
| 2107 | `2.22748e-10 V` | `1.91758e-10 V` | `-3.26801e-9 A/um` |
| 2114 | `2.07245e-11 V` | `1.81064e-11 V` | `-6.29018e-10 A/um` |

The absolute voltage differences are only tens of picovolts, but the
low-drain weak-inversion current is proportional to these small contact-edge
QF drops, so their relative effect is material.

## Endpoint mapping and exported GradQuasiFermi

The imported Sentaurus state uses the original TDR node IDs; no spatial
interpolation is performed. Vela SG transport uses the direct endpoint QF
difference. Sentaurus's exported `eGradQuasiFermi` is a nodal vector diagnostic
and is not assumed to expose its private face interpolation.

For edge 2114, the imported-state endpoint secant (`0.00222852 V/cm`), exported
Sentaurus gradient (`0.00222930 V/cm`), and Vela transport-cell-vector value
(`0.00222854 V/cm`) agree closely. Edge 2107 is less coincident: the endpoint
secant is `0.0236014 V/cm`, the exported Sentaurus mean gradient is
`0.0284499 V/cm`, and the Vela transport-cell-vector value is
`0.0253760 V/cm`. Despite this vector-gradient difference, the complete
Sentaurus-state drain current replay remains within `0.000395 dex` at the key
state. This confirms the prior M10 result that the HFS drive mismatch is not
the controlling terminal-current error here.

## Frozen contribution versus self-consistent sensitivity

The two mathematical layers now agree in magnitude:

| Layer | Key-state response / Id |
|---|---:|
| Exact frozen replacement of adjacent interior `phin` | `-13.622715%` |
| M18 adjoint QF response on the drain contact cut | `+13.622821%` |
| Complete M18 self-consistent feedback response | `+13.275136%` |

The sign difference is expected: the frozen game replaces Vela by Sentaurus
and reduces current, while the M18 ledger attributes the opposite
Sentaurus-to-Vela operator/state discrepancy. Their magnitudes differ by only
about `0.000106 percentage point` for the drain-cut QF term.

## Decision

M19 excludes three candidate explanations for the stable weak-inversion gap:

- the drain contact voltage or constrained contact `phin` is not wrong;
- spatial interpolation of imported nodal fields is not involved;
- the observed drain-cut SG terminal-current construction is not the dominant
  error, because it reproduces Sentaurus current on the imported state.

The remaining discrepancy is a self-consistent electron-continuity state
difference at the first free drain-adjacent nodes. The next investigation
should compare the complete electron-continuity row balance at those nodes:
drain-cut flux, internal-edge fluxes, and source terms/Jacobian couplings. It
should not tune HFS, Enormal, SRH, or reference density globally.

## Figures

- `m19_vg0p8_gap_closure.png`: self-consistent error versus imported-state SG replay.
- `m19_key_state_frozen_factors.png`: exact three-factor frozen endpoint attribution.
- `m19_key_state_edge_concentration.png`: edge concentration of the adjacent-QF correction.
- `m19_contact_boundary_error.png`: contact `phin` boundary closure across 16 states.
