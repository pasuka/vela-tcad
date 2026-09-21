# Templates/LDMOS G3 全网格 AverageBox 只读回放审计（2026-08-31）

## 结论

在 IALMob 与 predictor 继续关闭、`Vd=0.1 V`、`Vg=0.5 V` 的同一 G3
Sentaurus 固定状态上，已使用完整 T-2022.03-SP2
`MeasureCoefficients.debug` 对全部 Silicon 输运单元和边完成只读回放。

仅把 Vela SG 边通量的 couple 替换为直接 AverageBox 系数后，5,674 个有效自由
电子 continuity 行的残差 L2 从 `2.42852e-6` 降到 `4.20316e-8`，即
`0.0173075x`；最大绝对残差从 `6.04064e-7` 降到 `2.61989e-8`，即
`0.0433710x`。冻结热点 3747、10233、4492、4538 全部超过 `2x` 改善门。

候选 Top-12 节点集合与基线 Top-12 不重合，说明残差排名发生变化；但候选最大值
只有原最大值的 `4.34%`，低于冻结的 `50%` 实质热点转移阈值，因此**没有继续
发生有实质幅值的热点转移**。全网格 external-profile 实现资格门通过。

把源项体积再由 Vela volume 替换为直接 AverageBox `Measure` 后，各项指标没有
可见改善，证明本固定状态的主导差异在输运 couple 系数，而非 SRH/impact 源项
体积。该结果只授权下一步可选诊断 profile 与一次 Newton A/B；尚未授权修改生产
默认、31 点曲线或 known-difference ledger 审批，ledger 继续保持 draft。

## 输入合同与覆盖范围

- exact SProcess topology：19,782 个 Tri3，其中 10,515 个 Silicon 输运单元；
- Silicon 输运图：5,723 个节点、16,237 条无向边；
- 电子 continuity：5,674 个有效自由行；47 个 carrier boundary 行不评分；
- fixed state：G3 `Vd=0.1 V, Vg=0.5 V` 的 SSS 状态；
- oracle：T-2022.03-SP2 `MeasureCoefficients.debug`；
- oracle SHA-256：
  `05ec1a43936e0ff71bea7cdfd956d80c7d692b228187ccae8625ec53ce58d49e`；
- debug 的每个 Tri3 均要求恰好 3 个 `Measure` 和 3 个 `Coefficients`，否则审计
  立即失败；
- 只累计 Silicon 单元，因而界面 carrier flux 不引入 oxide-side couple；
- 保持现有 SG Bernoulli、密度、迁移率和 HFS 语义，仅替换几何系数；
- production solver 未被修改。

为验证重装配无索引漂移，脚本先用生产 couple 重算全部边和节点通量：couple 最大
绝对复现误差为 `1.05879e-22 m`，节点 flux 最大绝对复现误差为
`2.51463e-22`。

## 全自由行结果

| 方案 | residual L1 | residual L2 | maximum absolute | L2 / baseline | max / baseline |
| --- | ---: | ---: | ---: | ---: | ---: |
| production baseline | `3.94960e-5` | `2.42852e-6` | `6.04064e-7` | `1` | `1` |
| AverageBox coefficient only | `2.74685e-7` | `4.20316e-8` | `2.61989e-8` | `0.0173075` | `0.0433710` |
| coefficient + Measure | `2.74685e-7` | `4.20316e-8` | `2.61989e-8` | `0.0173075` | `0.0433710` |

`coefficient + Measure` 与 `coefficient only` 的 L1 只相差约 `3.14e-14`，L2
差异低于显示精度。固定状态下源项 Measure 不是主导因素。

在 5,674 个有效行中，3,587 行绝对残差改善，1,303 行至少改善 `2x`，2,087
行持平或变差。后一个计数包含大量原本已处于低残差地板的行，不能单独作为热点
转移判据；本轮以全局 L2、最大残差和冻结热点三项联合门判定。

## 冻结热点与排名转移

| node | baseline residual | AverageBox residual | 绝对值比值 |
| ---: | ---: | ---: | ---: |
| 3747 | `6.04064e-7` | `-1.23418e-11` | `2.04e-5` |
| 10233 | `-5.77592e-7` | `3.09634e-11` | `5.36e-5` |
| 4492 | `-5.40895e-7` | `-1.31409e-10` | `2.43e-4` |
| 4538 | `-5.10737e-7` | `2.74112e-10` | `5.37e-4` |

候选最大行转为 node 4571，残差为 `2.61989e-8`；其自身基线残差为
`2.27884e-7`，仍改善约 `8.70x`。第二大候选 node 4541 从
`-8.21687e-8` 改善到 `-2.44316e-8`，约 `3.36x`。因此排名更换来自原热点
大幅关闭后低一级残差浮到前列，不是此前 node 4492 局部实验中那种近同量级转移。

冻结实现门为：全自由行 L2 与 maximum 均不高于 baseline 的 `0.5x`，且四个冻结
热点分别不高于 `0.5x`。三项全部通过。

## 解释与边界

1. node 4492 局部审计的“AverageBox 边系数语义”结论已被全网格证实，而其一环
   转移是假设只替换局部系数造成的边界不连续，不代表全网格候选失败。
2. 本轮直接系数同时关闭界面热点和非界面 node 4492，说明 region-local
   Si/SiO2 carrier coupling 与普通 Silicon 网格 couple 应由同一完整 AverageBox
   离散合同处理，而不是分别加入经验修正。
3. 这是 fixed-state residual replay，不包含候选 couple 对 Jacobian、line search、
   contact branch 或自洽状态反馈的影响；不能据此宣称 31 点 Id-Vg 已关闭。
4. 不继承 PN2D 私有 `element_edge_sg_gss_laux` 原子捆绑。LDMOS 候选需要独立
   external-profile 配置、Jacobian 测试和一次 Newton 资格。
5. IALMob 与 predictor 继续关闭；不通过下一资格门前不运行高成本 31 点曲线。

## 决策与后续结果

全网格只读回放任务完成，结论为“通过，无实质热点转移”。后续已实现
**显式选择、默认关闭**的 LDMOS external AverageBox 诊断 profile，并完成同一
fixed state 上的一次 Newton A/B 与同偏压 reclose。实现与接触 BC 两因素结果见
`templates_ldmos_g3_averagebox_newton_ab_2026-08-31.md`。

known-difference ledger 仍为 draft：本轮识别的是可关闭的输运 couple 与 LDMOS
source-short 接触合同差异，不是可批准的引擎固有地板；31 点曲线继续后置。

## 可复现工件

- 审计脚本：`scripts/audit_templates_ldmos_averagebox_full_mesh.py`；
- 回归测试：
  `tests/regression/test_audit_templates_ldmos_averagebox_full_mesh.py`；
- ignored 大型输出：
  `reference_staging/templates_ldmos_averagebox_full_mesh_20260831/`；
- 主要输出：`summary.json`、`cell_coefficients.csv`、`edge_replay.csv`、
  `node_residuals.csv`。
