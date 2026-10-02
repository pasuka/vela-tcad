# SimpleMOS 生产支路一致性与局部场差审计

2026-09-07。**已完成八点生产支路审计及节点 1000/1009 局部对照。现有生产开关不能直接复现已通过的联合方案；剩余场差已有明确的迁移率离散候选。** 另定位了原生六列状态回放中的准费米势表示舍入问题，并完成表示方式对照。正式求解器、模型默认值和收敛门槛均未修改。

本轮承接[联合方案八点验证](simplemos_joint_geometry_eight_point_validation_2026-09-07.md)。完成 136 次只读探针、八点原生 SG 端口校准、局部通量/Poisson/SRH 数据账本。没有新 DC、自洽候选或 sdevice 计算。结论不等同于生产移植完成，也不放行 M82/M83。

## 1. 生产路径的具体缺口

| 路径 | 当前源码行为 | 对联合方案的含义 |
| --- | --- | --- |
| 耦合 DD 输运 | `CoupledDDAssembler` 按 `transport_edge_coupling` 选择 Si 侧局部系数，并缓存至残差/Jacobian 的边核 | 开关确实生效，但局部 cotangent/fallback 系数仍与原生处理后系数有差别 |
| 普通端口提取 | `ContactCurrent` 构造器始终使用 `computeEdgeCouplings(mesh)` | 没有消费传入的输运几何开关；需要显式同步系数 |
| 端口残差函数及伴随 | `makeArclengthContactCurrentFunctional` 直接使用同一个耦合装配器的无边界替换残差/Jacobian | 会消费装配器系数；普通提取器与该函数并非天然共用几何 |
| Poisson 电荷体积 | 三个独立开关各自调用 `computeTransportNodeVolumes` | 这是 Si 三角形面积/3 的重心体积，不是原生 Si box 体积 |
| signed Si 体积 | `transport_signed_average_box_node_volume` 修改公共 `vol_` | 本网格数值上可得到原生 Si 体积，但同时改变连续性/SRH 体积；不能直接替代八点方案的 Poisson 专用修改 |
| Poisson Jacobian | ψ、φn、φp 电荷导数使用对应 `poissonElectronVol_`/`poissonHoleVol_` | 三个电荷体积与相应导数已有一致连接；净掺杂是固定源，没有状态导数 |
| SG/载流子行探针 | 使用装配器实际 `couple_` 和源项 | 可以核对真正进入方程的几何和通量 |
| 普通迁移率探针 | 迁移率数值来自同一物理函数，但 `couple_m` 直接输出 `edge.couple` | 此列代表原网格系数，不能当成修改后 SG 的实际系数 |
| 独立 Poisson、Gummel DD | `PoissonAssembler`、`DDAssembler` 使用网格原节点体积/边系数 | Newton 开关不会自动覆盖这些路径，生产化时需要明确支持范围 |
| 存储/端口体电荷诊断 | `StoredCharge`/`TerminalCharge` 默认读网格体积，选区域时使用面积/3 | 若用于联合方案的电荷对照，须明确其体积定义；不能假定已跟随 Poisson |

源码依据：[耦合装配器](../../src/equation/CoupledDDAssembler.cpp) 439–468、555–586、1633–1638、5384–5391 行；[端口提取器](../../src/post/ContactCurrent.cpp) 65、498–503 行；[端口残差函数及状态打包](../../src/solver/NewtonSolver.cpp) 34–79、3258–3326 行；[几何公共函数](../../include/vela/equation/AssemblerUtils.h) 251–331、404–420 行；[迁移率探针](../../src/tools/vela_example_runner.cpp) 的 `writeEdgeMobilityProbeCsv`；[Gummel DD](../../src/equation/DDAssembler.cpp)、[独立 Poisson](../../src/equation/PoissonAssembler.cpp)、[存储电荷](../../src/post/StoredCharge.cpp)、[端口电荷](../../src/post/TerminalCharge.cpp)。

本次实际量化结果：

- 每个网格有 1857 条原生 Si 系数为正的边。默认系数与原生不同的有 **54 条**；开启现有 Si 侧开关后仍有 **14 条**。最大 `abs(g_Vela/g_native−1)` 约 153.45/153.53；联合方案在相同支持上达到浮点精度一致。极大相对值对应小原生系数，不应直接解读为同量级端口电流误差。
- 942 个 Si 支持节点中，现有重心 Si 体积在 n19/n23 分别有 **852/844 个**与原生不同，最大相对差 **48.9646%**。独立 signed Si 体积在全部 942 节点与原生一致，最大相对差 **8.44e-15**；缺口在生产策略的连接方式，而非必须重新拟合体积。
- 迁移率非零的 2691 条边中，现有 Si 侧开关使 **40 条**边的实际 SG 系数与普通迁移率探针的同名列不同；联合方案为 **54 条**。迁移率数值本身在两种探针间完全相同。
- 八点的漏接触支持没有被这些几何差异改变。因此合格状态的 **24/24** 默认/开关/联合端口对照均通过 1e-8 门槛。这能证明本算例的当前漏端提取一致，但不能覆盖“接触支持本身跨材料或需改系数”的缺口，不能据此宣布 `ContactCurrent` 已正确实现该开关。

现有 `region_resolved_interface_assembly: true` 还会同时改变 Poisson 边系数及公共源体积，亦不等同于八点方案。本轮使用显式对象，仅打开输运边及三项 Poisson 体积开关作比较，未用布尔总开关冒充联合方案。

证据：[几何比较](../../reference_tcad/simplemos_sentaurus2022/production_consistency_20260907/geometry_comparison.csv)、[诊断支路](../../reference_tcad/simplemos_sentaurus2022/production_consistency_20260907/diagnostic_paths.csv)、[端口路径](../../reference_tcad/simplemos_sentaurus2022/production_consistency_20260907/port_paths.csv)。

## 2. 生产数值实现尚未包含隔离修正

本轮逐文本核对了当前生产源码与封存 runner：`CoupledDDAssembler::residualImpl` 完全相同；隔离版本另有构造器几何覆盖、关闭的电荷覆盖、只读诊断和 Jacobian 数值修正。运行时关闭电荷、向量链、迁移率步长等诊断覆盖；默认/现有开关对照同时关闭外部原生几何覆盖。

因此这些探针检查的是当前 Masetti **物理残差路径**，并不是重新编译后认证了完整生产求解器。稳定解析 SG ψ 导数、100 位残差累加的线性迭代修正仍位于隔离副本；本轮未把它们移入正式源码，也没有用隔离版 Jv 通过来声称生产 Jacobian 已通过。当前 Masetti 配置仅依赖固定掺杂，状态迁移率偏导为零；场相关或 PhuMob 分支仍需各自的导数和邻接列验证。

现有三项 Poisson 电荷体积改动与实际残差差分的独立复核 **16/16 通过**，最大相对差 1.55e-12；五种状态/系数组合的定向边通量与载流子行通量和 **40/40 通过**，使用预先冻结的 128 epsilon 绝对边流舍入界。该舍入界只用于回放恒等式，不是非线性接受门槛。

证据：[运行前合同和源码哈希](../../reference_tcad/simplemos_sentaurus2022/production_consistency_20260907/contract.json)、[算子回放检查](../../reference_tcad/simplemos_sentaurus2022/production_consistency_20260907/operator_replay_checks.csv)、[稳定 SG 隔离实现](../../scripts/diagnostics/simplemos_stable_sg_psi_derivative.hpp)。源码差分保存在忽略目录 `build-release/simplemos_production_consistency_20260907/`。

## 3. 原生状态回放的表示与求值差异

原生一致状态采用六列 ψ、φn、φp、n、p CSV。没有 reference/increment 时，状态打包先做 `φ/V0−reference/V0`；普通电流提取器直接使用物理势作差。缩放后相减会改变接近接触电压的微小准费米势差。

| 原生状态，n23 | 六列输入的端口相对差 | 同状态显式 reference/increment |
| --- | ---: | ---: |
| Vg=1、Vd=1 V | 1.73816e-8，失败 | 3.553e-15，通过 |
| Vg=0.8、Vd=0.05 V | 1.18840e-7，失败 | 2.442e-15，通过 |
| Vg=0.8、Vd=1 V | 1.99460e-6，失败 | 1.998e-15，通过 |

表示方式对照保持 ψ 和密度不变，使用已冻结 contact-basin 参考；按实际读入 double 计算，物理准费米势的最大表示变化仅 5.55e-17 V。八点均恢复到端口差 ≤6.00e-15。提取器电流不变到导出精度，仅一个点变化 −2.27e-23 A/μm；最大残差函数电流变化为 2.75e-14 A/μm。原失败记录保留。

这定位了一个状态输入路径的数值缺口，并给出已验证的回放格式；没有修改正式 `packReferencedSolution`。它不能解释毫伏级局部场差，也不会推翻原八点合格 reference/increment 终态的结论。

Poisson 回放还发现另一种不能混用的状态：CSV 中预存密度与装配器由 ψ/φ 重新求值的密度。直接用前者重构小残差时，八点相对余项为 3.19e-8–4.14e-8，均未通过原 1e-8 回放门槛。随后使用 SG 探针实际导出的求值态重构同一算子，八点全部通过，最大相对余项降至 **1.70e-11**。两种密度最大相对差仅 6.47e-13，ψ 打包差 ≤1.11e-16 V，但在大电荷相消后的残差中会放大。原失败未被覆盖；该结果证明的是 Vela 算子回放，不是原生 Poisson 行残差已获认证。

证据：[表示控制输入变化](../../reference_tcad/simplemos_sentaurus2022/production_consistency_20260907/representation_input_changes.csv)、[端口表示对照](../../reference_tcad/simplemos_sentaurus2022/production_consistency_20260907/representation_port_comparison.csv)、[实际求值态 Poisson 回放](../../reference_tcad/simplemos_sentaurus2022/production_consistency_20260907/packed_Poisson_replay.csv)、[原输入态 Poisson 分项](../../reference_tcad/simplemos_sentaurus2022/production_consistency_20260907/local_Poisson_ledger.csv)。

## 4. 节点 1000/1009 的剩余差异

以下为高 NWell、Vd=1 V、已合格联合终态相对原生的差。势差采用相同物理参考；密度取原生实际导出，SRH 原生值取节点 `srhRecombination`，Vela 值由实际连续性 SRH 积分除以其源体积得到。相对差为 `100×(Vela/native−1)`。

| Vg | 节点 | Δψ (mV) | Δφn (mV) | Δφp (mV) | Δn | Δp | SRH 相对差 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.8 | 1000 | +0.03184 | −19.79083 | −1.91437 | +115.28% | −7.25% | −2.71% |
| 0.8 | 1009 | +0.03066 | −21.76259 | −17.41416 | +132.33% | −49.07% | −49.28% |
| 1 | 1000 | +0.02532 | −25.74067 | −13.06652 | +170.93% | −39.73% | −38.55% |
| 1 | 1009 | +0.02267 | −23.43308 | −20.44617 | +147.76% | −54.70% | −59.42% |

节点坐标分别为 (0.04664907, 0.09890192) 和 (0.03360702, 0.09890192) μm。这些点的局部电势差只有约 23–32 μV，准费米势差和密度差却明显；不能用接近的 ψ 或端口 Id 代替输运场验证。原生和 Vela 的四个 SRH 值均为负，所以上表负的 SRH 相对差表示净生成幅值减小。原生节点 SRH 乘 box 体积只保留作诊断求积，未冒充原生积分源或原生行残差。

1000–1009 对应边 2152。原生 Si 几何已一致，但有效迁移率仍不同：

| 载流子 | Vela 边迁移率 (cm²/Vs) | 原生 box 加权有效迁移率 | Vela 相对差 |
| --- | ---: | ---: | ---: |
| 电子 | 621.13931 | 573.76561 | +8.25663% |
| 空穴 | 289.02708 | 269.64729 | +7.18709% |

数值在 Vg=0.8/1 V 相同，因为本配置迁移率只依赖固定掺杂。原生有效边迁移率由已校准的单元迁移率和 Si 单元边系数加权得到，不是直接拿节点迁移率图当内部边值。Vela 仍沿用现有边端点掺杂及其单元平均规则。

该边的 Vg=1 V 物理势差 `φ(j)−φ(i)`：电子 Vela/native 为 −25.68177/−27.98937 mV；空穴为 **−1.07273/+6.30692 mV，方向相反**。有效 SG 电导又受到密度和电势形成过程影响，因此不能把终态通量差只归结为 8% 的迁移率倍数。

证据：[局部节点物理量](../../reference_tcad/simplemos_sentaurus2022/production_consistency_20260907/local_node_ledger.csv)、[邻接边的迁移率、准费米差、SG 电导和通量分解](../../reference_tcad/simplemos_sentaurus2022/production_consistency_20260907/local_edge_ledger.csv)。同时保留了低 NWell、低 Vd、沟道 320/324/338 及漏接触附近 1091/1092 的对照。

## 5. 固定状态迁移率候选的检验

使用原生 ψ、n/p、准费米势和已校准单元迁移率/边几何重构 SG，八点端口都满足预先冻结的 1e-6 门槛，最大相对误差 **2.132e-7**。这校准了此次原生通量重构的符号、单位和端口积分，不能等同于导出了原生边未知量或逐行残差。

在显式参考形式的原生一致状态上，保持联合几何、状态及 SRH 源不变，仅把每条相邻边的 Vela 迁移率替换为原生有效迁移率，代数重算局部残差。每个边源在两端等量反号。以下是高 NWell、高 Vd 的 `abs(替换后行残差)/abs(替换前行残差)`，没有进行 Newton 重算：

| Vg | 节点 | 电子行剩余比例 | 空穴行剩余比例 |
| ---: | ---: | ---: | ---: |
| 0.8 | 1000 | 1.225e-7 | 4.861% |
| 0.8 | 1009 | 2.133e-7 | 4.367% |
| 1 | 1000 | 7.483e-7 | 3.790% |
| 1 | 1009 | 4.104e-7 | 16.981% |

这说明迁移率离散几乎解释了这两个节点在固定原生一致状态下的电子连续性不闭合，并解释了大部分空穴不闭合。剩余空穴差仍需要核对 SRH 求值、源体积及状态约定。**此结果支持把原生有效迁移率列为下一个候选，但尚未证明它会在自洽后改善高低 NWell 的 Id 配对或全部场差。** 局部代数源不能直接套用尚未独立校准的任意连续性权重。

证据：[八点原生 SG 端口校准](../../reference_tcad/simplemos_sentaurus2022/production_consistency_20260907/native_SG_port_calibration.csv)、[固定状态迁移率候选](../../reference_tcad/simplemos_sentaurus2022/production_consistency_20260907/local_mobility_counterfactual.csv)。电导/势差双因子代数分解的最大闭合相对余项为 3.0e-16；它只证明代数账本闭合，不是因果电流预测。

## 6. 验证边界与后续安排

仍使用 n19/n23 原网格（1480/1482 节点、2742/2746 三角形）、300 K、Boltzmann/no-BGN、matched ni、Masetti 总杂质和掺杂相关 SRH。网格坐标 μm、状态密度 m⁻³、迁移率表 cm²/Vs、端口 A/μm、宽度 1 μm；局部通量账本为粒子线通量 m⁻¹s⁻¹，Poisson 账本为 C/m。原 `terms.csv` 中部分带 `_m3` 的掺杂列仍是 unit_scaling 内部数值，本轮未将其当作物理 m⁻³ 使用。

运行工具为 Windows UCRT64 Python、已封存 C++20 Release 隔离 runner。继承的线性后端为 Eigen SparseLU，但本轮只读探针没有执行 DC 或线性迭代。五个新增 Python 脚本执行完成；未编译或修改 C++，未重跑全量 CTest。工作树原 runner 的 127 行新增/1 行删除及先前八点证据保持不变。首次汇总中的旧合同字段兼容和 NumPy 整数序列化问题已修正；已有仿真和表格未覆盖，恢复汇总时逐项校验了已写表格。

本轮审计任务完成，下一阶段应按证据分开推进：

1. **生产一致性修复**：统一装配器、端口提取器和探针的有效几何来源；为三项 Poisson 电荷提供独立 signed Si 体积策略，保持本轮未改的 SRH 体积；处理六列状态打包精度，分别验证稳定 SG 导数及线性修正的生产移植。不要用布尔总开关或全局 mesh 体积策略代替已验证联合方案。
2. **物理候选验证**：在合格联合状态上对原生有效迁移率做独立的小幅度正负自洽响应校准，覆盖两个 NWell 和两个 Vd，再做有限替换。保持原接受门槛，同时检查 Id、节点 1000/1009、SRH、端口守恒与双初始化。
3. 上述通过后才扩展 16 工况与完整 0–1 V 曲线。当前未执行这些新增生产修复或自洽计算。

复现入口：[准备与执行](../../scripts/prepare_simplemos_production_consistency_20260907.py)（`prepare → run`）、[状态表示控制](../../scripts/probe_simplemos_state_representation_20260907.py)、[主分析](../../scripts/analyze_simplemos_production_consistency_20260907.py)（`analyze`）、[实际求值态 Poisson 回放](../../scripts/audit_simplemos_packed_poisson_replay_20260907.py)、[局部迁移率候选账本](../../scripts/audit_simplemos_local_mobility_counterfactual_20260907.py)，最后运行主分析入口的 `seal`。从本验证工作树根目录使用 `D:/msys64/ucrt64/bin/python.exe`；已有目录封存，复跑须选择新输出目录并重新冻结输入。

汇总：[分析证据](../../reference_tcad/simplemos_sentaurus2022/production_consistency_20260907/analysis_evidence.json)、[实际求值态补充证据](../../reference_tcad/simplemos_sentaurus2022/production_consistency_20260907/packed_Poisson_evidence.json)、[局部候选补充证据](../../reference_tcad/simplemos_sentaurus2022/production_consistency_20260907/local_counterfactual_evidence.json)、[最终哈希](../../reference_tcad/simplemos_sentaurus2022/production_consistency_20260907/validation_evidence.json)。
