# Templates/LDMOS G3 Sentaurus 电子方程余额对齐（2026-09-01）

## 结论

已在 Sentaurus T-2022.03-SP2 虚拟机上对 `Vd=0.1 V`、
`Vg=1.0/1.1666667 V` 两个最大 gm 端点执行 loadable-state 单步
`NewtonPlot(Residual Error Update)`，并将原生 `eContinuityRhs` 与 Vela 在完全相同
Sentaurus `psi/phin/phip`（SSS）状态上的 production SG/continuity 装配对齐。

七个预注册节点的坐标误差为零。Sentaurus 七节点原生电子 RHS L2 分别为
`3.50691e-15 A` 和 `4.41779e-15 A`，仅为端口电流的 `2.80e-9` 与
`1.56e-9`。Vela 在相同 SSS 状态上的物理化七节点余量分别为
`1.24147e-9 A/um` 与 `6.97873e-9 A/um`，相对端口电流为 `9.90e-4`
与 `2.47e-3`；高低端点向量余弦为 `0.990984`，L2 增长 `5.62135x`。

因此，同状态局部 continuity 算子差异得到原生 Sentaurus equation-balance 证实，
不是节点编号、坐标、状态所有权、端口电流幅值或 SRH 源项造成的假象。不过，收敛
equation-balance 只约束每行所有边贡献的总和；TDR 的 `eCurrentDensity` 是顶点
后处理场，`ElementEdgeCurrent` 只切换电流密度重构而不导出装配边通量。因此本轮
不能从接近零的 Sentaurus 行和唯一反演其逐边 QF-drive 或密度支撑公式。

结论等级为“同状态局部算子差已确认，具体逐边 QF-drive 规则仍开放”。不新增 C++
profile，不启动 31 点曲线，不调整迁移率、阈值、IALMob、predictor 或生产默认值；
known-difference ledger 继续保持 draft。

## Sentaurus 原生探针

两个有效探针均从对应的最终 loadable TDR 状态开始，物理段与 G3 生产参考一致：

- `Physics { Fermi }`；
- Silicon 中仅 `Mobility(HighFieldSaturation)`；
- `EffectiveIntrinsicDensity(OldSlotboom)`；
- `SRH(DopingDependence TempDependence)` 与 Auger；
- `Iterations=1`、`CNormPrint`、`NewtonPlot(Residual Error Update)`。

`Vg=1.0 V` 首行 electron CNorm 为 `6.450399e-11`，`Vg=1.1666667 V`
为 `1.920008e-14`。两次探针的漏端电子电流分别为约 `1.254e-6` 和
`2.828e-6 A/um`，与被冻结参考端点一致。初始化网格、缺少 SLP 的中间 Plot、错误
FilePrefix 等试运行未进入证据链。

Sentaurus 导出的全硅 RHS 包含接触/约束行，故全域 L2 不作为自由 continuity 行门槛；
本报告只对预注册的七个非欧姆接触热点节点作原生行余额判定。

## 同状态对齐

| 指标 | `Vg=1.0 V` | `Vg=1.1666667 V` |
| --- | ---: | ---: |
| 坐标最大误差 (um) | `0` | `0` |
| Sentaurus 七节点 RHS L2 (A) | `3.50691e-15` | `4.41779e-15` |
| Sentaurus 七节点 max abs (A) | `1.98709e-15` | `4.11889e-15` |
| Sentaurus 七节点 L2 / terminal | `2.79623e-9` | `1.56221e-9` |
| Vela SSS 七节点 residual L2 | `2.81103e-5` | `1.58018e-4` |
| Vela SSS 七节点 L2 (A/um) | `1.24147e-9` | `6.97873e-9` |
| Vela SSS 七节点 L2 / terminal | `9.89885e-4` | `2.46780e-3` |
| 七节点占 Vela 全硅 residual energy | `0.1284%` | `0.6277%` |
| Vela SSS terminal / Sentaurus | `1.00020215` | `1.00020219` |

Vela 的缩放余量使用同一次 SG 导出的
`electron_particle_line_flux_per_m_s / electron_flux` 反求 continuity particle scale。
两个端点都得到 `2.75651794697415e20`，12,680 条非零边的最大相对离散仅
`2.22e-16`；乘以基本电荷和 `1e-6 m/um` 后，转换因子为
`4.416428645843634e-5 A/um` 每单位缩放余量。

七节点 Vela 行的 recombination 项不超过 `4.52e-23`，远小于输运余量；行余量
相对入射边绝对通量和为约 `5.75e-5--2.56e-3`。这说明差异是多条大边通量相消后
留下的稳定小量，不能把任一顶点 `eCurrentDensity` 直接当成装配真值。

## 对 QF-drive 假设的判定

本轮证据关闭了以下问题：

1. Sentaurus 在这七个节点上的自身离散方程确实闭合到约 `1e-15 A`；
2. Vela 用同一 SSS 状态和已资格化 external AverageBox/HFS 合同不能逐行闭合；
3. Vela 七节点余量随偏压以几乎相同空间方向增长，而不是随机舍入噪声；
4. SSS 端口电流仍只差 `0.0202%`，所以端口积分匹配不能替代局部行等价。

本轮不能关闭“具体是哪一种逐边语义”问题。原生 Sentaurus 行和在零附近，无法把
40 条入射边的贡献唯一分解；已知十种 generalized-Einstein 因子替代全部失败，
所以剩余候选仍是完整 SG kernel 中的 QF-drive placement、密度支撑或未枚举的联合
约定，而不是一个可直接标定的单一因子。

下一决定性实验应在相同 VSV 状态上让 Sentaurus 输出 iteration-0 NewtonPlot，或
取得受支持的装配级 element-edge current。前者需要把 Vela 诊断状态写入复制的
loadable TDR；该副本已在本地完成字段回读验证，但向虚拟机传输具体 TDR 需要独立
数据传输授权，故没有绕过授权继续执行。

## 可复现工件

- 审计脚本：`scripts/audit_templates_ldmos_g3_sentaurus_equation_balance.py`；
- 回归：`tests/regression/test_audit_templates_ldmos_g3_sentaurus_equation_balance.py`；
- ignored 原生 deck、日志、NewtonPlot、导出与汇总：
  `reference_staging/templates_ldmos_g3_sentaurus_equation_balance_20260901/`；
- Vela SSS 输入：
  `reference_staging/templates_ldmos_g3_averagebox_state_feedback_vg1_20260831/` 与
  `reference_staging/templates_ldmos_g3_averagebox_state_feedback_vg1p166667_20260831/`。
