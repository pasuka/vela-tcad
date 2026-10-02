# SimpleMOS M36 boundary and contact measure audit

## Scope

M36 reuses the n23 `MeasureCoefficients.debug` oracle captured in M34 and
audits the signed `AverageBoxMethod` node measure on the complete 1,482-node,
2,746-triangle mesh.  It separates Si nodes into interior, Si/SiO2 material
interface, non-contact external boundary, direct-contact, and contact
first-layer classes.  No new Sentaurus execution, C++ change, or production
default change is involved.

The local `Measure` slots use the independently established M34 permutation
`[0, 2, 1]`.  Source and drain each contain three direct contact nodes and five
first-layer free nodes; the substrate contains 29 nodes in each class.

## Result

The assembled Si node measure closes against Sentaurus in every audited class:

| Node class | Si nodes | Median relative error | Maximum relative error |
|---|---:|---:|---:|
| Interior | 777 | 2.79e-16 | 3.99e-15 |
| Material interface | 39 | 1.73e-15 | 4.03e-15 |
| External boundary | 52 | 3.78e-16 | 2.41e-15 |
| Direct contact | 35 | 2.27e-16 | 1.78e-15 |
| Contact first layer | 39 | 2.02e-16 | 2.18e-15 |

The drain direct-contact measure sum is `1.3579358504782543e-05 um2` and its
first-layer sum is `9.887701219911767e-05 um2`; source is symmetric.  Both
agree with the raw signed circumcentric assembly to floating-point precision.

Individual element-local shares are not required to match slot by slot: their
largest relative mismatch is large when a Sentaurus local share is near zero.
The physically assembled node sums nevertheless close to a maximum relative
error of `4.03e-15`.  This shows that Sentaurus redistributes some local shares
while preserving the nodal `Measure` used by the assembled equations.

## Interpretation

M36 rejects the hypothesis that the remaining SimpleMOS deep-off mismatch is
caused by an unknown external-boundary or source/drain-contact node-measure
formula.  The signed formula reconstructed in Vela is already equivalent at
the assembled-node level across all relevant classes.

This does not mean the measure is numerically irrelevant: M35 showed that
enabling the signed node measure self-consistently changes the current by
about `-0.00513 dex`.  It means the remaining investigation should separate
where that known measure is applied from contact SG flux construction and
quasi-Fermi self-consistent feedback, rather than continue fitting the
measure formula itself.
