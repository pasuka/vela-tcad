# Templates/LDMOS G3 单模 QF 微扰审计（2026-09-01）

## 结论

已在同一个 `Vd=0.1 V, Vg=1.1666667 V` VSV 状态上，对非接触硅节点
`3721` 的电子准费米势施加 `±1 µV` 对称微扰，并分别执行 Sentaurus
T-2022.03-SP2 iteration-0 `NewtonPlot` 与 Vela fixed-state carrier-term
装配。比较集合在运行前固定为 node 3721 及与其共享 Silicon 三角形的完整一环：
`3714, 3720, 3721, 3722, 3723`。

Sentaurus 与 Vela 的中心差分电子 continuity 向量余弦为
`-0.9999998311`，最小二乘比例为 `-33.7806012`；负号来自双方 residual
符号约定。去掉单一标量后，相对 L2 仅 `0.0005812`，最大逐行剩余相对
Sentaurus L2 为 `0.0004085`。该单模比例相对未扰动 VSV 高端点的
`33.4599563` 仅漂移 `0.9583%`。

因此，前轮约 `33.5x` 并非“大基线残差偶然共线”或局部 QF-drive 空间形状
不同：它同样乘在受控的局部 Jacobian 列上。现有证据把主候选收窄为
Sentaurus NewtonPlot `eContinuityRhs` 归一化、二维宽度/电流约定，或双方统一
assembly coefficient 的标量差；不能据此调整迁移率、阈值或 SG 空间支撑。
ledger 继续保持 draft，直至该标量的单位归属被独立关闭。

## 实验合同

- 基准：G3 no-IALMob、predictor off、external AverageBox carrier couples、
  `legacy_node_local` 接触重构及已资格化 contact HFS fallback；
- 状态：Vela `psi` + Sentaurus `phin` + Vela `phip` 的高端点 VSV；
- 唯一自变量：node 3721 的 `phin`，幅值为 `±1e-6 V`；
- Sentaurus：两个复制的 loadable TDR，Load 后立即 Plot，再执行
  `Coupled(Iterations=1)` 并导出 iteration-0 RHS；
- Vela：同一对状态运行 `newton_carrier_term_probe`；
- 差分：`[R(+1 µV)-R(-1 µV)]/(2 µV)`；
- 生产默认、IALMob、predictor、迁移率与阈值参数均未改变。

两个 Sentaurus loaded 文件的 node 3721 `phin` 差为
`2.000000000002e-6 V`，误差 `2.00e-18 V`；全硅其余 `phin`、全部 `psi`
及 `phip` 的两支差值严格为零。每支 loaded 状态相对目标状态存在相同的全硅
`psi max=142.145 µV` 回读差，但由于两支逐点完全相同，它在中心差分中严格抵消，
不属于单模响应。

## 定量结果

| 指标 | 结果 |
| --- | ---: |
| 预注册一环节点数 | `5` |
| Sentaurus `dR/dphin` L2 | `2.8842581 A/V` |
| Vela 物理化 `dR/dphin` L2 | `0.08538207 A/(um V)` |
| 向量余弦 | `-0.9999998311` |
| signed Sentaurus/Vela 最小二乘标量 | `-33.7806012` |
| 去标量后相对 L2 | `0.05812%` |
| 相对未扰动高端点标量漂移 | `0.95829%` |
| Sentaurus even/odd L2 | `0.01869%` |
| Vela even/odd L2 | `0.02000%` |

两侧 even/odd 都远低于 `0.1%`，证明 `±1 µV` 位于线性区，中心差分没有被
二阶项污染。按全硅 `1e-8` 相对阈值自动识别的 active support 与预注册一环
完全一致，没有发生热点转移。

Vela continuity 物理化仍使用 SG edge probe 冻结的转换：particle scale
`2.75651794697415e20`，即每单位 scaled residual 为
`4.416428645843634e-5 A/um`；12,680 条非零边的最大相对离散为
`2.22e-16`。

## 判定与下一门

本实验支持：

1. 单模 Jacobian 列的空间方向已经跨引擎闭合；
2. `~33.5x` 是近似统一的行幅值/单位标量，不是局部 generalized-Einstein、
   QF edge average 或热点节点选择造成的形状差；
3. 不能把该比例直接解释为迁移率倍率，因为固定状态端口电流和 gm operator 已在
   先前实验中闭合到约 `0.02%`；
4. 31 点曲线、IALMob 和 predictor 继续暂缓，本轮不改生产 C++ 默认。

后续量纲链与 `AreaFactor=1 -> 2` A/B 已完成。端口电流精确变为 `2.0x`，
但正、负两个完整 NewtonPlot `eContinuityRhs` CSV 均逐字节不变，中央差分
候选/基线比例为 `1.0`。因此二维宽度与 AreaFactor 已被排除；NewtonPlot RHS
必须视为实现相关内部归一化量。详见
`templates_ldmos_g3_rhs_dimension_chain_2026-09-01.md`。ledger 仍保持 draft，
下一门转为 Sentaurus 方程缩放/控制盒 Measure/行预条件顺序审计。

## 可复现工件

- 状态/deck 准备脚本：
  `scripts/prepare_templates_ldmos_g3_qf_single_mode_perturbation.py`；
- 审计脚本：`scripts/audit_templates_ldmos_g3_qf_single_mode_perturbation.py`；
- 回归：
  `tests/regression/test_audit_templates_ldmos_g3_qf_single_mode_perturbation.py`；
- ignored TDR、日志、导出、逐行 CSV 与汇总：
  `reference_staging/templates_ldmos_g3_qf_single_mode_perturbation_20260901/`；
- 正式数值摘要：上述目录下 `analysis/summary.json` 和
  `analysis/one_ring_rows.csv`。
