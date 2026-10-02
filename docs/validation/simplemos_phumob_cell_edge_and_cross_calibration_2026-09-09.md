# SimpleMOS PhuMob 单元、边及交叉导数校准

日期：2026-09-09。分支：`codex/simplemos-sdevice-validation`。

本轮完成了此前未收尾的 PhuMob 校准：16 个原生状态完成并取回，单元形成方式和 SG 端口重放已核验，生产装配器的弱交叉列已与独立高精度链式导数比较。结果支持“顶点本构 → box 单元平均 → 单元边系数加权”的形成方式，但**空穴 G(P) 截断值和生产差分导数仍有资格缺口，不能放行 PhuMob 的 element_box 生产分支或继续恢复 Enormal/HFS**。

## 1. 工况、求解器与证据

- 原生矩阵：n19/n23 × Vd=0.05/1 V × Vg=0.8/1 V，分别运行 PhuMob 和 Masetti 控制，共 16 点。保持 OldSlotboom、掺杂相关 SRH、300 K；只将 `Mobility(DopingDependence)` 改为 `Mobility(PhuMob)`，采用默认 arsenic 参数。Enormal、高场饱和保持关闭。
- 网格：沿用原始 SimpleMOS 二维网格；n19 为 1480 节点、2742 三角形，n23 为 1482 节点、2746 三角形，每个网格有 1750 个硅单元。坐标和原生 box 顶点面积分别为 μm、μm²；载流子为 cm⁻³，迁移率为 cm²/(V·s)，端口电流为 A/μm。
- 原生为 Sentaurus T-2022.03-SP2，`ExtendedPrecision(128) Method=Super RelErrControl Digits=12 ErrRef(Electron)=1e-2 ErrRef(Hole)=1e-2 RhsMin=1e-20 Iterations=40 ExitOnFailure CNormPrint`。16 点均通过退出码、运行版本、偏置和 KCL 检查；16 份 TDR 场导出完成。8 个 Masetti 控制的 Id 与前置参考相同。
- Vela 使用已冻结的 Release `vela_example_runner.exe` 和当前生产标量 PhuMob 内核。沿用已有 BGN 匹配材料、Delaunay box 转移、逐单元介电系数和 signed Si Poisson 电荷体积。实际交叉列检查仅装配残差和 Jacobian，不调用线性求解器，也不进行新的 PhuMob 自洽 DC；原 Masetti 合格状态只是固定状态输入。
- 9 月 8 日的 16:40 截止保护已结束。恢复时核实原生批次有 complete 标记及 16 个退出记录，没有待恢复的暂停组；本轮未继续运行旧截止保护。

原生合同及结果见 [native_contract.json](../../reference_tcad/simplemos_sentaurus2022/phumob_calibration_20260908/supported_export/native_contract.json)、[native_points.csv](../../reference_tcad/simplemos_sentaurus2022/phumob_calibration_20260908/supported_export/native_points.csv)。前置资格见 [BGN 电势导数修复](simplemos_bgn_psi_derivative_validation_2026-09-08.md)。

## 2. 单元形成方式与边通量

在同一个原生状态上，不拟合参数，对四种预先冻结的单元形成候选进行比较。下表为 8 个 PhuMob 状态全部硅单元的最大绝对相对差；每一项按“公式值/原生单元值−1”计算。

| 单元形成候选 | 电子 | 空穴 |
|---|---:|---:|
| 每个顶点使用本地 ND、NA、n、p 求 μ，再按原生 box 顶点面积平均 | 1.38265e-9 | 7.35772e-6 |
| 每个顶点使用本地掺杂及单元算术平均 n、p，再平均 μ | 1.22810 | 0.195626 |
| 每个顶点使用本地掺杂及单元几何平均 n、p，再平均 μ | 0.620457 | 0.208158 |
| 所有输入先按 box 面积平均，再求一次 μ | 0.491658 | 0.258440 |

“算术/几何平均”均使用原生 box 顶点面积权重，不是等权三顶点平均。独立 Masetti 控制的电子、空穴最大相对差均为 4.44e-16。原生 runtime 单元迁移率与 TDR 单元数据另行对照，避免把绘图节点平均值直接当成求解器本构。

边的有效迁移率系数按 `Σ_cell g_cell,edge μ_cell` 形成，保留原生有符号单元边系数。使用原生单元 μ、原生 n/p/ψ/准费米势/有效 ni 计算可变 ni 的 SG 边通量，再按接触边界求和，16 点端口重放全部通过 `1e-6` 相对门槛，最大相对差为 **3.59924e-7**。此重放使用已独立校准的原生常数约定，不修改 Vela 常数。

采用第一种本构候选时，32 个单元检查（16 状态 × 两种载流子）有 8 项失败，均为空穴 PhuMob；对应 16 个边检查有 8 项失败。`1e-7` 单元/边门槛保持原值。端口重放通过不能覆盖这项空穴本构失败，也不等于 Vela 自洽 Id 已通过。

明细：[cell_summary.csv](../../reference_tcad/simplemos_sentaurus2022/phumob_calibration_20260908/supported_export/analysis_20260909/cell_summary.csv)、[edge_comparison.csv](../../reference_tcad/simplemos_sentaurus2022/phumob_calibration_20260908/supported_export/analysis_20260909/edge_comparison.csv)、[port_replay.csv](../../reference_tcad/simplemos_sentaurus2022/phumob_calibration_20260908/supported_export/analysis_20260909/port_replay.csv)。

## 3. 剩余本构差定位到 G(P) 截断分支

依据本地 T-2022.03 手册第 400–402 页的公式 281–283及参数，独立计算 G(P) 最小值，并按每个单元中触发 `P < P_min` 的顶点数分类。两种载流子各有 14000 个单元观测。

| 载流子 | 无截断顶点的单元数 | 这些单元的最大相对差 | 含截断顶点的单元数 | 最大相对差 |
|---|---:|---:|---:|---:|
| 电子 | 10960 | 5.55e-16 | 3040 | 1.38265e-9 |
| 空穴 | 11040 | 5.55e-16 | 2960 | 7.35772e-6 |

这将已观察到的微小偏差定位到 G 截断分支，而不是一般单元平均或未截断本构式。进一步按冻结规则，只从 n19、Vd=0.05、Vg=0.8 的第一个全截断单元 148，分别反推一个有效 G 下限；其他全部单元作为留出对照：

| 载流子 | 当前公式精确最小值 | 原生数据反推的有效下限 |
|---|---:|---:|
| 电子 | 0.09810869654913677 | 0.09810965009295847 |
| 空穴 | 0.07994310999664854 | 0.07994487380175797 |

仅在原截断顶点替换该诊断下限后，全部留出单元的最大相对差降至 6.66e-16；独立公式与原 C++ 本构的差也不超过 6.66e-16。由此可说明剩余偏差具有“统一截断下限差”的形式。

**上述两个下限是数据反推值，不是已核实的 Sentaurus 内部参数。** 尚未识别原生最小值搜索、温度处理或近似形成算法；本轮没有把它们写进生产代码，也没有用拟合后的对照覆盖原 8 项失败。手册公式求精确极小值与原生有效下限的差异，不能直接定性为某一端的物理错误。

明细：[groups.csv](../../reference_tcad/simplemos_sentaurus2022/phumob_calibration_20260908/supported_export/analysis_20260909/screening_floor/groups.csv)、[floor_inference.csv](../../reference_tcad/simplemos_sentaurus2022/phumob_calibration_20260908/supported_export/analysis_20260909/screening_floor/floor_inference.csv)。

## 4. 标量与实际生产 Jacobian 交叉导数

先在 8 个此前合格状态的节点 792、1000、1009、1057 上，检查 32 组标量输入、128 个载流子偏导。独立 100 位链式导数与 100 位中心差分最大相对差为 1.66672e-51；C++ 标量值与高精度值最大相对差为 4.50679e-16。双精度 6 个对数密度步长的 768 次检查有 490 次失败、403 次返回零；其中 76 个偏导在全部 6 个步长都未通过。这里的步长是对数密度，不是电势伏特。

再将检查推进到生产装配器：取 Vg=0.8 的 n19/n23 × 两个 Vd 共 4 个固定状态，保留合格输运几何，迁移率明确使用现有 `phumob + edge_averaging=legacy`，关闭 SRH 以隔离迁移率交叉项。对四个节点的 φn/φp 使用 `1e-5`、`3e-6 V` 单列方向，导出实际 Jv 和残差差分。

独立参考在每条有效边上对 PhuMob 的平均载流子状态求 100 位链式导数，乘以同边导出的基态通量/迁移率，再按守恒符号装配至相邻行。单位是**装配器原始残差/物理伏特**。标量 μ 最大相对差 5.74e-16，基态边通量与行残差的误差除以绝对边通量和最大为 2.01e-16，均通过独立状态/单位预检。

| 交叉块 | 非零独立矩阵项 | 未过 1e-5 相对门槛 | 生产值为零 |
|---|---:|---:|---:|
| 电子行对 φp | 82 | 71 | 70 |
| 空穴行对 φn | 82 | 35 | 24 |
| 合计 | 164 | 106 | 94 |

两个外部步长合计 328 项检查，其中 212 失败、188 项生产值为零。独立参考绝对值范围为 2.70534e-29～3.69701e-16；最大绝对导数差为 3.69701e-16。按同输入列的所采样直接载流子块 L1 范数衡量，交叉误差最大比例约为 6.69130e-4。该比例只是量级参照，不能替代原导数门槛或自洽电流响应验证。

当前生产 PhuMob 分支确实包含电子/空穴互相依赖的逻辑，但通过**完整边通量的双精度中心差分**计算这些偏导，并非全部由符号求导后的解析式给出。对很弱的迁移率变化，通量差分发生消减，所以“源码存在这个交叉块”并不等于其数值完整性已验证。本轮没有证明这些小项是此前 Id 差异的主因。

现有 legacy 分支只依赖边两端的载流子状态。将来按单元顶点平均恢复 PhuMob 时，边相邻单元第三顶点的 n/p 也会影响迁移率；这些额外列需要显式装配和独立验证，不能仅删除现有 element_box 模型保护。

证据：[标量 summary.json](../../reference_tcad/simplemos_sentaurus2022/phumob_calibration_20260908/cross_derivatives/runtime_ready/summary.json)、[实际矩阵 summary.json](../../reference_tcad/simplemos_sentaurus2022/phumob_calibration_20260908/assembled_cross_20260909/active_transport_analysis/summary.json)、[cross_entries.csv](../../reference_tcad/simplemos_sentaurus2022/phumob_calibration_20260908/assembled_cross_20260909/active_transport_analysis/cross_entries.csv)、[column_scale_context.csv](../../reference_tcad/simplemos_sentaurus2022/phumob_calibration_20260908/assembled_cross_20260909/active_transport_analysis/column_scale_context.csv)。

## 5. 原生 PhuMob 对电流的影响

下表为原生 `Id_PhuMob/Id_Masetti−1`，两端使用相同 BGN、SRH 和偏置。**这不是 Vela–Sentaurus 的电流误差。**

| NWell | Vg | Vd=0.05 V | Vd=1 V |
|---|---:|---:|---:|
| n19 | 0.8 V | +15.67058% | +16.72184% |
| n19 | 1.0 V | +22.18661% | +21.36001% |
| n23 | 0.8 V | +21.94505% | +21.94611% |
| n23 | 1.0 V | +22.94281% | +23.97072% |

恢复 PhuMob 的整体影响远大于此次 G 截断的微小差异。旧 Masetti 或 BGN+Masetti 下的 Id 误差不能直接作为完整原始模型的验收结果。本轮未补做 PhuMob 自洽 A/B、0–1 V 全曲线或双初始化。

## 6. 失败记录与本轮验证边界

- 初始原生回调请求了不支持的 `Edge-RegionWise eMobility`；12 个已结束失败、4 个取消均保留。改用此前合格的单元迁移率/几何回调后，另建 supported_export 批次，未覆盖初始记录。
- 首次标量探针缺少运行库 PATH，启动失败已保留；在独立 runtime_ready 目录重跑。
- 本轮映射预检曾错误要求全部 runtime 顶点对 TDR 全局单射。第一个状态有 1562 个 runtime 顶点，但实际硅拓扑使用的 942 个顶点一一映射；改为核验硅顶点单射、全部硅单元双射及坐标精度，并保留原分析器和失败记录。没有放宽物理门槛。
- 实际矩阵分析首次对零迁移率的非输运边做除法而退出；原输入、原分析器和原始结果保持冻结。后续独立分析器只在同时确认零边通量时排除这些边，再通过基态残差闭合核验。
- 8 个本轮 Python 脚本语法检查通过；原生、场导出、单元/边、G 分支、标量导数、实际矩阵的 6 组结果身份复核通过。数值不合格项目按上表报告，未计作测试通过。本轮未改生产源码或重新构建求解器，未重跑完整 CTest；前次 766/782 是前置 BGN 阶段记录。

原失败证据：[初始运行账本](../../reference_tcad/simplemos_sentaurus2022/phumob_calibration_20260908/initial_attempt_ledger.csv)、[映射预检](../../reference_tcad/simplemos_sentaurus2022/phumob_calibration_20260908/supported_export/analysis_20260909/preflight_failure/failure.json)、[交叉分析预检](../../reference_tcad/simplemos_sentaurus2022/phumob_calibration_20260908/assembled_cross_20260909/analysis_preflight_failure.json)。所有生产默认、物理常数和收敛接受条件保持不变；全局 SRH 相对源门槛此前低于下限未激活，本轮也没有新增其独立响应资格。

## 7. 后续任务

1. 先为 PhuMob 标量载流子偏导提供避免完整通量相减的稳定链式求值，分别验证 G 截断分支、双载流子交叉块和装配后的全部相关列；以本轮 164 个矩阵项作为复核集，保留原失败结果。
2. 独立核实原生有效 G 下限的形成方式。反推值仅用于定位，不作为拟合生产参数；需要用额外温度或独立参数控制辨别原生最小值近似与其他处理。
3. 按已核验的顶点本构及 box 单元/边形成方式实现受保护的 PhuMob 候选，补齐相邻单元全部载流子列及端口/探针一致性。通过这些资格后，再进行两个 NWell、两个 Vd、两种 Vg 的自洽和双初始化对照，随后才恢复 Enormal → 高场饱和。

本阶段的校准和定位已完成；生产 PhuMob 恢复尚未放行。本轮保留在当前工作区，未另行提交或推送。
