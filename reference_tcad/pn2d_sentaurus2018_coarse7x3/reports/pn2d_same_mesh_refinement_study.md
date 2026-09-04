# PN 二极管三级同网格加密试验

三套网格均由原始网格嵌套细分，SDevice 与 Vela 在每一级使用完全相同的节点、三角形和掺杂数据。

## 结论

- 网格稀疏显著影响 SDevice 的反向端口漏电：7x3 结果仅为 25x9 的 0.196%，13x5 仍比 25x9 高 40.4%。因此原 7x3 网格不适合用作 -20 V 端口漏电的定量基准。
- 节点恢复顺序的差异随加密消失：13x5 和 25x9 的 direct/cell-first P95 基本相同。粗网格确实放大了恢复方法敏感性。
- 但 Vela 相对 SDevice 的绝对节点矢量 P95 没有随加密下降，典型节点 P50 也基本不变。剩余误差不能只归因于网格，输运模型、SG 边系数和 SDevice 节点矢量语义仍占主要部分。

## 端口与节点矢量

| 网格 | 节点/三角形 | SDevice |I| [A/um] | 相对25x9 | 电子 P50/P95 direct [dec] | 空穴 P50/P95 direct [dec] |
|---|---:|---:|---:|---:|---:|
| 7x3 | 21/24 | 6.986980e-19 | 0.00196411 | 0.238739/1.43699 | 0.1267/1.11722 |
| 13x5 | 65/96 | 4.994684e-16 | 1.40405 | 0.235262/1.69933 | 0.12908/1.35966 |
| 25x9 | 225/384 | 3.557333e-16 | 1 | 0.232151/1.69489 | 0.130057/1.35358 |

## 恢复方法 A/B

| 网格 | 电子 P95 direct/cell-first [dec] | 空穴 P95 direct/cell-first [dec] |
|---|---:|---:|
| 7x3 | 1.43699/1.46635 | 1.11722/1.13173 |
| 13x5 | 1.69933/1.69933 | 1.35966/1.35966 |
| 25x9 | 1.69489/1.69489 | 1.35358/1.35358 |

## 解释边界

- SDevice does not expose its directed internal edge flux; Vela SG section flux is checked for internal consistency and against SDevice terminal current.
- Vela section fluxes are recomputed on a frozen SDevice state, not a converged Vela state; their section spread diagnoses operator/model compatibility, not Vela self-consistent conservation.
- Absolute current-vector errors include transport-model differences as well as reconstruction differences.
- The 25x9 result is the finest level in this three-level experiment, not a proof of asymptotic mesh convergence.
