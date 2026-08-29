# Templates/LDMOS G3 Id-Vg exact-point qualification

Status: fail

Exact shared points: 31; resolved: 30; unresolved biases: [0.0].

| Metric | Value | L2 limit | Pass |
| --- | ---: | ---: | --- |
| Median abs log error (dex) | 0.0065760605 | 0.1 | True |
| P95 abs log error (dex) | 0.28721918 | 0.2 | False |
| Strong-inversion endpoint relative error | 0.01451373 | 0.2 | True |
| Fixed-current Vth absolute error (V) | 0.029626123 | 0.1 | True |
| Maximum gm relative error | 0.064249711 | 0.2 | True |
| Worst resolved KCL relative error | 0.06162799 | 0.01 | False |

Curve errors are evaluated only at exact shared CurrentPlot points; no curve interpolation is used.
Vth uses local log-current interpolation at 1.0e-08 A/um and remains diagnostic until that current level is formally frozen.
