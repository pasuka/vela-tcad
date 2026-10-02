# SimpleMOS M20 drain-adjacent electron-row audit

## Scope

M20 is a read-only SDevice audit of the first free-node layer next to the drain
contact. It evaluates the 16 frozen M10/M14 states twice, using the converged
Vela state and the exact-node imported Sentaurus state in the unchanged Vela
operator. The target set is defined topologically: every free endpoint of an
edge whose other endpoint is a drain-contact node. For n21 this gives nodes
990, 991, 992, 995, and 997; for n17 the corresponding mesh-local IDs are 991,
992, 993, 996, and 998.

The audit separates raw electron continuity terms before row scaling and also
partitions every nonzero entry of the production Jacobian row by variable block
(`psi`, `phin`, `phip`) and by contact/free column. No Sentaurus job is rerun
and no production physics or default is changed.

## Execution and closure

| Check | Result | Acceptance |
|---|---:|---:|
| Frozen states | 16 | 16 |
| Read-only Vela probes | 96 | 16 x 2 x 3 |
| Maximum incident-edge versus raw flux closure | `9.68843e-16` | `<1e-10` |
| Maximum complete raw-term closure | `9.68843e-16` | `<1e-10` |
| Maximum exact Jacobian six-part closure | `2.03810e-16` | `<1e-12` |
| New Sentaurus execution | no | no |
| Default model changed | no | no |

The existing per-edge transport-Jacobian probe is retained only as a control.
It constructs an edge-projection mobility drive and therefore cannot reconstruct
the SimpleMOS production row, which uses `transport_cell_vector`. The headline
Jacobian result comes directly from exact production-matrix entries, not that
incompatible reconstruction.

## Key state: n21, Vd=0.05 V, Vg=0.8 V

The imported Sentaurus state is not a solution of the Vela equations, so its
Vela continuity residual identifies where the two self-consistent states cease
to agree. The raw residual difference is:

| Node | Contact-cut flux delta | Internal-edge flux delta | Residual delta | Delta cancellation |
|---:|---:|---:|---:|---:|
| 990 | `0` | `+2.18914e-7` | `+2.18914e-7` | `1.00x` |
| 991 | `-1.72301e-4` | `+1.72571e-4` | `+2.70424e-7` | `1275.3x` |
| 992 | `-3.31664e-5` | `+3.32335e-5` | `+6.70886e-8` | `989.7x` |
| 995 | `-2.86129e-6` | `+2.87450e-6` | `+1.32115e-8` | `434.2x` |
| 997 | `-3.03116e-7` | `+3.33586e-7` | `+3.04695e-8` | `20.90x` |

Node 990 is not directly connected to the drain-contact cut; its discrepancy
is entirely an internal-edge flux imbalance. Nodes 991 and 992 are the main
contact-adjacent rows. Their contact and internal SG changes are individually
about three orders of magnitude larger than the remaining residual and cancel
to 0.157% and 0.202%, respectively.

Across all five nodes, the absolute state-change ledger is 49.9282% contact
flux and 50.0718% internal flux. Their signed sums nearly cancel. The absolute
SRH share is `6.77e-19`; impact, gauge, and boundary source changes are zero on
these free rows. This excludes local SRH generation/recombination as the source
of the key-state residual platform.

## Production Jacobian coupling

At all five key-state rows, more than 99.9999993% of the absolute production
Jacobian row belongs to the electron quasi-Fermi (`phin`) block. Direct `psi`
shares range from `1.35e-11` to `6.14e-9`, while `phip` shares are below
`1.1e-24`.

| Node | Contact-column phin share | Free-column phin share | Imported-state flux condition |
|---:|---:|---:|---:|
| 990 | `0%` | `~100%` | `156.1x` |
| 991 | `44.4211%` | `55.5789%` | `7884.3x` |
| 992 | `44.4227%` | `55.5773%` | `6840.9x` |
| 995 | `44.4228%` | `55.5772%` | `3877.4x` |
| 997 | `0.000230%` | `99.999770%` | `121.1x` |

Thus the first-layer equations are not controlled by a local SRH derivative
or a direct Poisson-row coupling. They are stiff `phin` transport balances. At
the directly connected rows, about 44.42% of absolute stiffness comes through
contact columns and 55.58% through free interior columns.

## Decision

M20 localizes the stable weak-inversion mismatch to the balance of electron SG
fluxes across the drain interface and the immediately adjacent internal edges.
It does not support contact-voltage error (already excluded by M19), local SRH
source error, or direct `psi/phip` Jacobian coupling as the controlling cause.

The remaining distinction is narrower: determine why the imported Sentaurus
state produces a slightly different contact-versus-internal `phin` flux balance
under Vela, especially at nodes 991 and 992. Because those rows are conditioned
by cancellation factors near `10^3` (and absolute state flux conditions near
`10^4`), picovolt-level QF differences can be amplified into the observed row
residual. A follow-up should compare the individual incident-edge endpoint QF
drops, Bernoulli factors, carrier populations, mobility, and geometric couple
for these rows, while preserving the production `transport_cell_vector` drive.

## Figures

- `m20_key_state_flux_balance.png`: signed contact/internal flux changes and residual.
- `m20_key_state_flux_conditioning.png`: state-difference and imported-state cancellation factors.
- `m20_key_state_jacobian_partition.png`: exact production Jacobian contact/free and variable-block shares.
