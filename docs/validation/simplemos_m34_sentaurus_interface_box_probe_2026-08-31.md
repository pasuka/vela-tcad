# SimpleMOS M34 Sentaurus interface box probe

Date: 2026-08-31
Sentaurus Device: T-2022.03-SP2
Scope: n23 geometry oracle; one Poisson initialization; no production default change

## Outcome

The direct `MeasureCoefficients.debug` export closes the Si/Oxide edge geometry
question without requiring an explicit double-node implementation.

- The input grid has 1482 shared vertices and no coincident duplicate vertex
  groups.  The 20 Si/Oxide interface edges share 21 input vertex ids.
- SDevice reports 1564 internal vertices and 4309 internal edges, versus 1482
  input vertices and 4227 geometrical edges.  Both increases are 82, exactly
  the number of unique contact nodes.  `NumberOfDoubleEdges` is zero.
- All 2746 material triangles have complete `Measure` and `Coefficients`
  records.  Grid-element and device-element ids are identical; the local
  vertex slots use the fixed input-to-debug permutation `[1, 0, 2]`.

The evidence therefore supports shared Si/Oxide potential nodes with
region-local element contributions.  It does not support ordinary Si/Oxide
double points in this SimpleMOS run.  The extra internal vertices are fully
accounted for by contact-boundary expansion.

## Edge coefficient closure

After applying the documented local permutation, all 20 interface edges close:

| Quantity | Result |
|---|---:|
| Maximum absolute Si local-coefficient error | 2.43e-16 |
| Maximum absolute oxide local-coefficient error | 1.03e-15 |
| Vela region-local / Sentaurus Poisson coefficient | 0.999999999999995 to 1.000000000000006 |
| Legacy Vela / Sentaurus Poisson coefficient | 1.427269192633890 to 1.427269192633906 |

Thus the M33 region-local edge construction
`epsilon_Si * coefficient_Si + epsilon_oxide * coefficient_oxide` is the direct
Sentaurus geometry oracle on this interface.  The historical construction,
`average(epsilon) * total_coefficient`, overweights every interface edge by
42.7269%.

## Node-measure qualification

Sentaurus `AverageBoxMethod` node measures are region-local, but they are not
identical to the barycentric Si subvolume used by the first M33 diagnostic:

| Total/Si ratio over 21 interface nodes | Minimum | Median | Maximum |
|---|---:|---:|---:|
| Sentaurus direct `Measure` | 2.99204 | 4.76065 | 12.95225 |
| M33 barycentric diagnostic | 3.65606 | 4.83373 | 11.28750 |

Consequently, M33 proves the direction of the edge correction but its
`transport_node_volume` corner is not an exact Sentaurus node-volume oracle.
The 4.25% all-on current-gap closure must not be used to reject a coherent
region-local AverageBox implementation.

## Decision

Do not implement explicit Si/Oxide double nodes for SimpleMOS from the current
evidence.  The next discriminating task is an M35 self-consistent rerun using:

1. the now-validated region-local edge coefficients;
2. region-local signed AverageBox node measures instead of barycentric shares;
3. the existing shared electrostatic node and semiconductor-only carrier rows.

This will test the coherent Sentaurus geometry without changing the default
model and without introducing an unsupported interface-equation architecture.

## Artifacts

- `reference_tcad/simplemos_sentaurus2022/sentaurus_interface_box_probe/`
- `reference_tcad/simplemos_sentaurus2022/simplemos_m34_sentaurus_interface_box_probe_evidence.json`
- `scripts/run_simplemos_m34_sentaurus_interface_box_probe.py`
- `scripts/freeze_simplemos_m34_sentaurus_interface_box_probe_evidence.py`
