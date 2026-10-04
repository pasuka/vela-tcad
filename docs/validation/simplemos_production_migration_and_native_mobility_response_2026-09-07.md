# SimpleMOS 生产移植与原生有效迁移率响应验证

2026-09-07。**生产一致性修复已通过八点对照；原生有效迁移率的小扰动校准和有限替换已执行，但有限替换未通过全部候选放行条件。未扩展 16 工况或完整 0–1 V 曲线。**

承接[生产一致性及局部场审计](simplemos_production_consistency_and_local_fields_2026-09-07.md)。本轮完成 32 次生产对照 DC（含 8 次关闭线性修正的负控制）和 56 次迁移率候选 DC。正式模型默认值、原接受门槛、SRH 源体积和历史失败记录均保留。使用已有原生网格、单元迁移率和状态，无新增上传或 sdevice 仿真。

## 1. 正式实现及支持范围

- 装配器、普通端口提取器及迁移率探针共用 `computeEffectiveTransportEdgeCouplings`；显式 `transport_edge_coupling_ratios` 按 mesh edge id 指定已校准 Si 系数。长度、符号、有限性及与旧开关的冲突均检查。
- `poisson_charge_node_volume="signed_transport"` 只给电子、空穴和掺杂三项 Poisson 电荷赋予 signed Si 体积；对应电荷导数使用相同体积。介电矩阵 K、连续性/SRH 体积继续使用原定义。没有把方案实现成布尔总开关或全局 mesh 体积默认值。
- 六列状态先在物理电压中作参考差，再无量纲化；显式 reference/increment 的精度仍保留。
- equal-ni、Boltzmann、无 BGN 且指数未截断的 SG ψ 解析偏导改为稳定通量乘 Bernoulli 对数导数。其他分支保留原计算；本轮不声称任意截断态或所有物理模型的 Jacobian 已全部认证。
- `linear_refinement_iterations` 默认 0，当前验证显式设 4；以 100 位精度累加线性缺陷，使用同一个双精度矩阵、线性后端和行权重求修正。只移植耦合 Newton 的迭代路径，不改变接受条件，也不把独立 step 探针当成完整非线性迭代。

显式几何策略限定耦合 Newton；Gummel DD 直接收到这两项配置会明确报错。独立 Poisson、StoredCharge/TerminalCharge 等原有体积定义未自动改变。字段说明见[配置规范](../config_schema.md)。

## 2. 生产对照结果

24/24 个开启四次线性修正的状态通过原门槛（八点重放、Vela 初始化、原生初始化），8/8 双初始化通过。最大逐行比值为 9.422028e-7，门槛 1e-6；最大 KCL/Id 为 1.07033e-14，门槛 1e-8。关闭线性修正的 8/8 对照失败，失败记录未被成功结果覆盖。

与上轮隔离联合方案相比，最大 Id 相对变化 1.67368e-10；最大 ψ/φn/φp 变化分别为 8.13191e-11 / 2.62416e-09 / 1.72673e-08 V，最大密度相对变化 6.67928e-07。均满足原初始化不变性门槛。

24/24 项“合格联合态 / 原生六列态 / 原生 referenced 态”端口对照通过，最大提取器与残差函数相对差 8.77076e-15。48/48 项几何、边通量和迁移率探针一致性通过。合成接触边测试补上了原网格漏端接触未触及几何修改的覆盖缺口。

生产 carrier-term 探针未提供旧隔离程序的验收 JSON 字段，首个汇总脚本因此中止。原三点结果和错误日志保留；补充 32 次只读边通量导出，以全部 1814 条有效载流子行和接触边守恒和独立重算验收，没有调整门槛。全局闭合仍按 tolerance=1e-6、source_floor=1e-10 判断；不将低于 source_floor 的净源解读为已达到任意相对精度。

证据：[生产矩阵](../../reference_tcad/simplemos_sentaurus2022/production_migration_20260907/dc.csv)、[表示与端口](../../reference_tcad/simplemos_sentaurus2022/production_migration_20260907/ports.csv)、[旧新等价](../../reference_tcad/simplemos_sentaurus2022/production_migration_20260907/equivalence.csv)、[独立验收补充合同](../../reference_tcad/simplemos_sentaurus2022/production_migration_20260907/acceptance_addendum.json)。

## 3. 迁移率方向的独立校准

候选是固定的电子、空穴原生 Si box 有效迁移率：

`μ(α) = μVela × [1 + α × (μnative_eff / μVela − 1)]`。

几何、状态变量、材料、SRH 及边界定义保持联合基线。α=±0.001、±0.0005 各作自洽计算，用 α=0 状态重新求完整 DD 端口伴随，独立构造每条边的守恒源及端口直接项，再比较 `−λᵀδR + δId_direct` 与实际响应。不是使用重建电流显示值替代守恒端口电流。

首版候选只替换残差/探针的迁移率，漏掉解析 Jacobian 的缓存，导致 64 项有限替换 Jv 检查失败。零扰动生产 Jv 通过。该失败保留，未运行首版 DC；补齐候选缓存后重建至独立目录，受影响的 224 项 Jv 块检查全部通过，最大相对差 5.56092e-9（门槛 1e-4）。弱 SRH 交叉块保持相同解析值；不把其差分噪声作为新的弱列认证。

40/40 次小扰动 DC 合格，16/16 组幅度校准通过。最大预测相对误差 1.57848e-07、最大双幅度差 1.17466e-07、最大偶/奇分量比 0.000364103、最小信号/零漂移 313877。原门槛依次为 0.001、0.001、0.01、100。

证据：[首版失败](../../reference_tcad/simplemos_sentaurus2022/native_mobility_candidate_20260907/preflight.csv)、[缓存修复预检](../../reference_tcad/simplemos_sentaurus2022/native_mobility_candidate_20260907/cache_consistent/preflight.csv)、[同扰动校准](../../reference_tcad/simplemos_sentaurus2022/native_mobility_candidate_20260907/cache_consistent/response.csv)。

## 4. 有限替换及场差

表中误差为 `IdVela/IdSentaurus−1`，电流单位 A/μm；前值是合格联合几何基线，后值为 α=1 原生有效迁移率替换。状态资格见末列。

| NWell | Vg/V | Vd/V | 原生 Id | 替换前 | 替换后 | 双初始化资格 |
| --- | --- | --- | --- | --- | --- | --- |
| n19 | 0.8 | 0.05 | 5.77663231e-07 | -0.013368% | +0.100217% | True |
| n19 | 0.8 | 1.0 | 5.2568207e-06 | -0.287643% | +0.159276% | True |
| n23 | 0.8 | 0.05 | 2.0360598e-09 | +0.191855% | +0.255979% | True |
| n23 | 0.8 | 1.0 | 1.37763758e-08 | +0.351202% | +0.405286% | True |
| n19 | 1.0 | 0.05 | 6.48538281e-06 | -0.762563% | +0.021044% | True |
| n19 | 1.0 | 1.0 | 3.83983089e-05 | -1.013748% | +0.075233% | True |
| n23 | 1.0 | 0.05 | 2.08718283e-07 | +0.025813% | +0.146691% | True |
| n23 | 1.0 | 1.0 | 1.24108214e-06 | +0.073109% | +0.250034% | True |

有限替换 DC 资格为 16/16；双初始化资格为 8/8。候选还必须同时改善高 NWell 绝对误差、保持低 NWell 不恶化，并改善高低 NWell 的 log 电流配对；不是仅以状态收敛作为物理候选放行条件。

| Vg/V | Vd/V | 高 NWell 改善 | 低 NWell 不恶化 | 配对改善 | 通过 |
| --- | --- | --- | --- | --- | --- |
| 0.8 | 0.05 | False | False | True | False |
| 0.8 | 1.0 | False | True | True | False |
| 1.0 | 0.05 | False | True | True | False |
| 1.0 | 1.0 | False | True | True | False |

高 NWell、高 Vd 下节点 1000/1009 的带符号场差如下（势差单位 mV，SRH 为 Vela/原生−1）。

| Vg/V | 节点 | Δψ 前 → 后 | Δφn 前 → 后 | Δφp 前 → 后 | SRH 差前 → 后 |
| --- | --- | --- | --- | --- | --- |
| 1.0 | 1000 | +0.025321 → +0.032939 | -25.740675 → -0.026815 | -13.066517 → +0.030981 | -38.5516% → -0.0849% |
| 1.0 | 1009 | +0.022671 → +0.035281 | -23.433075 → -0.031768 | -20.446166 → +0.000662 | -59.4194% → -0.2566% |
| 0.8 | 1000 | +0.031835 → +0.031951 | -19.790833 → -0.046325 | -1.914368 → +0.046029 | -2.7144% → -0.0073% |
| 0.8 | 1009 | +0.030658 → +0.030848 | -21.762586 → -0.068077 | -17.414163 → -0.017775 | -49.2755% → -0.2806% |

原生 SRH 是节点 Plot 数据；Vela 速率由本轮原体积和装配源项还原。这里比较物理速率，不将原生节点速率乘体积冒充原生求解器导出的积分残差。所有八点的详细局部密度与场差见[局部场表](../../reference_tcad/simplemos_sentaurus2022/native_mobility_candidate_20260907/cache_consistent/local_fields.csv)。

在合格的 n23、Vg=Vd=1 V 状态，节点 1000 的 φn 差由 −25.7407 mV 降至 −0.026815 mV，φp 差由 −13.0665 mV 降至 +0.030981 mV；节点 1009 的 φn/φp 差降至 −0.031768/+0.000662 mV。这支持有效迁移率离散是这些局部准费米场差的重要来源。但同一点 Id 误差由 +0.0731% 增至 +0.2500%，ψ 差也从约 +0.0253/+0.0227 mV 变为 +0.0329/+0.0353 mV。剩余误差包含耦合与抵消，不能用局部准费米场接近原生来替代端口放行条件。

## 5. 验证、环境和后续边界

- Release 构建成功；新增 9 个 Catch2 用例、64 个断言全部通过，覆盖高精度 SG 差分、弱线性缺陷修正、端口几何、SRH/Jacobian 不变性、状态打包及配置冲突。
- 全套 CTest **754/770 通过**。16 项失败全部位于历史证据检查：13 项源码哈希/历史替代链不符，3 项引用已不存在的 `tests/test_mos_mixed_material.cpp`。没有修改旧哈希或删除失败测试以制造全绿；[逐项记录](../../reference_tcad/simplemos_sentaurus2022/production_migration_20260907/ctest_audit.json)。所以本轮不能报告“全套测试通过”。
- Windows MSYS2 UCRT64、C++20、Release；实际 DC 线性后端为 Eigen SparseLU/COLAMD。构建确认 HDF5/TDR、UMFPACK 和 SPQR 可用，不能把编译可用误写为本轮 DC 使用 UMFPACK。
- n19/n23：Boron 1e17/2e17 cm⁻³；节点 1480/1482、三角形 2742/2746；Si 942 节点、907 自由节点。300 K、Boltzmann/no-BGN、matched ni=1.0750038488844236e10 cm⁻³，Masetti 总杂质、掺杂相关 SRH；HFS、表面迁移率、Auger、雪崩、DG 均关闭。网格坐标 μm，状态密度 m⁻³，宽度 1 μm。
- max_iter=200、reltol=1e-7、abstol=1e-12、stall floor=1e-9；ψ/准费米上限 0.35/0.025 V、contact_basin、原标量线搜索及全有效载流子行 eps=1e-6 均保持不变。

**本轮不放行 16 工况和完整曲线。** 后续应先依据有限替换中未改善的控制点，分开检验迁移率候选与剩余 Poisson/介电及 SRH 离散约定的作用，并处理历史证据引用与生产回归接续。局部场改善不能替代绝对 Id 和高低配对条件。

生产与候选结果、首版失败、源码前像及本报告由[本轮最终证据](../../reference_tcad/simplemos_sentaurus2022/native_mobility_candidate_20260907/cache_consistent/final_evidence.json)统一关联。未提交或合并代码。
