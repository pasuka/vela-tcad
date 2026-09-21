# Templates/LDMOS G3 eContinuityRhs 量纲链与 AreaFactor A/B

## 结论

`AreaFactor=1 -> 2` 使 Sentaurus 一步 Newton 后的 drain electron、hole 和
conduction current 全部精确变为 `2.0x`，但相同 VSV 状态、相同 node 3721
`±1 uV` 扰动产生的 NewtonPlot `eContinuityRhs` 完全不变：正、负两个完整
5723-node CSV 分别逐字节同哈希，五节点一环中央差分的候选/基线最小二乘
比例为 `1.0`，拟合后相对 L2 为 `0.0`。

因此，NewtonPlot 中虽标注为 `A` 的 continuity RHS 是与器件面积因子解耦的
内部归一化量；前一轮测得的 `33.7806012x` 不能由二维默认厚度、`AreaFactor`
或 Vela 的 `A/um` 端口单位解释。该标量现在只剩 Sentaurus 内部 RHS
归一化或一个均匀 assembly coefficient 两类归属。它仍不是可批准的引擎
差异地板，ledger 保持 `draft`。

## 文档与 Vela 量纲链

T-2022.03-SP2 *Sentaurus Device User Guide* 给出三项直接约束：

- 第 51 页：二维器件默认第三维厚度为 `1 um`；
- 第 58 页：全局 `Physics.AreaFactor` 是电流和电荷的乘数，在一维/二维中通常
  表示剩余维度的延伸；
- 第 208--210 页：CNormPrint/NewtonPlot 输出是内部 Sentaurus Device 数据，
  属于 implementation dependent，解释可能随版本变化。手册没有给出
  `eContinuityRhs` 的行归一化公式。

Vela 的显式链为：TCAD 内部长度 `1 um`、浓度 `cm^-3`、迁移率
`cm2/(V s)`、电流密度 `A/cm2`；SG 审计得到 scaled continuity residual 到
particle line flux 的比例 `2.75651794697415e20 m^-1 s^-1`，再乘
`q * 1e-6` 得 `4.416428645843634e-5 A/um`。这条链没有自由宽度参数，且
与生产端口的 `A/um` 约定一致。默认 `1 um` 宽度只会给出数值比例 `1`，无法
产生 `33.7806`。

## 单因素合同

| 项目 | 基线 | 候选 |
| --- | ---: | ---: |
| 完整状态 | high-endpoint VSV | 同左 |
| QF 扰动 | node 3721, `±1 uV` | 同左 |
| 网格/物理/Math/Solve | 冻结 | 同左 |
| 全局 `Physics.AreaFactor` | `1.0` | `2.0` |
| production default | 不修改 | 不修改 |

日志分别确认 `DeviceAreaFactor = 1` 与 `2`。两个诊断均按预注册的一环
`[3714, 3720, 3721, 3722, 3723]` 评分。

## 结果

| 指标 | 结果 |
| --- | ---: |
| common silicon nodes | 5723 |
| active central-difference support | 5 nodes |
| minus RHS CSV，AF1 vs AF2 | SHA-256 完全相同 |
| plus RHS CSV，AF1 vs AF2 | SHA-256 完全相同 |
| AF2/AF1 derivative fitted scale | `1.0000000000000000` |
| cosine | `1.0000000000000002` |
| relative L2 after scale | `0.0` |
| drain electron current，minus | `1.626e-6 -> 3.252e-6 A` |
| drain hole current，minus | `-2.490e-14 -> -4.980e-14 A` |
| drain conduction current，minus | `1.626e-6 -> 3.252e-6 A` |

五个一环节点的中央差分 AF2/AF1 比例均逐位等于 `1.0`。正扰动态也给出
完全相同的 `2.0x` 端口缩放，因此结果不是扰动符号或单个接触表的偶然值。

## 解释边界与下一门

本实验关闭的是“`33.7806x` 来自二维宽度或 AreaFactor”的假设，并证明
NewtonPlot RHS 的 TDR unit 标签不能直接当作端口安培真值。它没有证明
`33.7806x` 本身是任意值：未扰动双偏压行模式和 node-3721 Jacobian 列仍
稳定落在同一标量附近，说明其来源依然是稳定的均匀内部合同。

该下一门已完成：Electron `ErrRef=1e10 -> 1e8 cm^-3` 只改变 update error，
不改变 RHS 或 Newton 增量；五节点 Measure 跨 `14.61x` 而行比例只跨
`1.005x`；iteration-0 NewtonPlot 又在第一次线性求解前写出。详见
`templates_ldmos_g3_internal_equation_scaling_2026-09-01.md`。现有公开合同仍不能
推导绝对标量，因此 NewtonPlot RHS 只能用于形状/相对变化比较，不能把其绝对
幅值作为跨引擎物理电流门。

## 可复跑工件

- 审计脚本：`scripts/audit_templates_ldmos_g3_rhs_dimension_chain.py`；
- 回归测试：`tests/regression/test_audit_templates_ldmos_g3_rhs_dimension_chain.py`；
- ignored oracle：`reference_staging/templates_ldmos_g3_rhs_dimension_chain_20260901/`；
- 汇总：上述目录的 `analysis/summary.json`；
- 逐节点表：上述目录的 `analysis/one_ring_rows.csv`。
