# SimpleMOS M45 post-QF rebaseline

M45 re-runs the frozen-state analyses affected by M43 SG-kernel unification and M44 quasi-Fermi reference/increment preservation.

## Stage result

| Stage | Re-run status | Impact | Historical contract |
|---|---|---|---|
| M31 | complete | numeric refresh; qualitative conclusion retained | retained |
| M32 | complete | numeric refresh; qualitative conclusion retained | retained |
| M37 | complete | bitwise unchanged control | retained |
| M41 | failed | M43 operator consistency now propagated | superseded; failure retained |
| M42 | failed | M43/M44 operator and QF packing closure propagated | superseded; failure retained |

## Key closure

- M42 maximum compensated cut/report gap: 2.314386e-15 dex
- M42 minimum legacy cut/report gap: 1.928655e-16 dex
- M42 compensated free-electron residual L2: 1.613741e-20
- Legacy/compensated residual separation: 1.207361e+09

M31/M32 retained their qualitative conclusions, M37 was bitwise unchanged, and M41/M42 deliberately retain failed historical checks. Those failures now mean that the old mixed-operator mismatch no longer reproduces; they are not execution failures.

The machine-readable metric table contains 13 rows.
