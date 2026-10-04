# SimpleMOS PhuMob 单元平均候选与八点自洽验证

日期：2026-09-09。工作分支：`codex/simplemos-sdevice-validation`。前置报告为 [稳定链式导数修复](simplemos_phumob_stable_chain_repair_2026-09-09.md)。

**已实现显式 `element_box_phumob` 候选，补齐相邻单元全部载流子列，并完成八点自洽 A/B。** n19/n23 × Vd=0.05/1 V × Vg=0.8/1 V 的最大绝对相对电流误差，从 legacy 的 **1.130937% 降至 0.00652926%**，八点全部改善。两种离散各 16 次首次求解均合格，各八个双初始化比较通过；没有使用恢复重载或调整接受条件。

本轮验证的是单元/边迁移率形成及其 Jacobian 的显式候选。原生 G 截断内部算法仍未完全核实，未采用数据反推参数；该本构资格缺口、弱 SRH 相对源资格和完整 0–1 V 曲线仍独立保留。候选没有成为默认模型。

## 1. 配置与实现范围

对照保留同一 SimpleMOS Si/SiO2 网格和掺杂；n23 为 1482 节点、2746 三角形。两端均使用 plain PhuMob、OldSlotboom、匹配基础 ni 的独立材料副本、掺杂相关 SRH、300 K、Boltzmann 统计。几何仍为 `element_box` 输运、`signed_transport` Poisson 电荷体积、`cell_material` 介电装配及 `delaunay_transfer` 单元几何。两组仅切换迁移率的边平均方式。

Windows UCRT64 Release，全目标构建成功，线性求解后端为 UMFPACK。输入采用 `unit_scaling`；本报告电流为 A/μm、电势为 V、密度为 m⁻³。独立 PhuMob 高精度参考内部使用 cm⁻³ 和 cm²/(V·s)，比较时显式换算。没有新增 Sentaurus 仿真，沿用已合格的八个原生 PhuMob 点及其冻结输出。

新增 `solver.mobility.edge_averaging: element_box_phumob`，必须显式选择 plain `model: phumob`、`doping_concentration_basis: total_impurity` 和上述配套输运/介电几何。原 `element_box` 对 PhuMob 的保护仍保留；表面、高场模型不接受新候选，耦合 Newton 还要求 Boltzmann、关闭雪崩和未截断的密度指数区间。Gummel 继续拒绝该输运方案。

候选按以下顺序形成边迁移率：

1. 每个单元顶点使用本地 ND、NA、n、p 和材料温度计算 PhuMob。
2. 按实际 box 顶点面积平均为单元 μ，不将载流子或掺杂先做单元平均。
3. 按相邻输运单元的有符号局部边系数，形成 `Σ(g_cell,edge μ_cell)/Σg_cell,edge`。

稳定的 `dμ/dln(n)`、`dμ/dln(p)` 经同一权重传播至相邻单元每个顶点，再转换到 ψ、φn、φp 列，向边两端连续性行以相反符号装配。单元第三顶点的载流子列明确包含在内；稀疏矩阵模板同步扩展。G 的精确极小值规则、物理常数及 SRH 积分体积保持原值；迁移率平均方式改变了离散残差及自洽电流，这是本次 A/B 的预定变量。

相关实现见 [AssemblerUtils.h](../../include/vela/equation/AssemblerUtils.h)、[CoupledDDAssembler.cpp](../../src/equation/CoupledDDAssembler.cpp)、[MobilityModel.cpp](../../src/physics/MobilityModel.cpp)。[ContactCurrent.cpp](../../src/post/ContactCurrent.cpp)、[FixedStateOperatorAudit.cpp](../../src/equation/FixedStateOperatorAudit.cpp) 和 [vela_example_runner.cpp](../../src/tools/vela_example_runner.cpp) 同步传递完整的节点 n/p；缺少完整状态时拒绝候选计算，不退回端点均值。旧的端点 Jacobian 分解探针无法表示单元第三顶点列，因此显式拒绝候选，应使用完整 `newton_jvp_probe`。配置说明见 [config_schema.md](../config_schema.md)。

## 2. 首批输入资格失败及复核

首批固定态审计误用了模型恢复前的初始化文件。八个文件中的 n/p 仍对应旧基础 ni；当前配置下，Newton/SG 由 ψ/φn/φp 重算的 n/p 均为文件保存值的 **1.36174446517103 倍**。单元和独立迁移率探针直接读取保存的 n/p，因此两类探针实际使用了不同状态。

该批结果中，单元 μ 最大相对差为 0.31010922，探针一致性有 29649/67584 项失败、最大差 0.23670486；其边公式、864 个弱交叉项和 Jv 分块本身通过。整批资格仍判为失败，没有据这些局部通过结果进入 DC 对照。该差异是审计输入状态的语义不一致，不能解释为本轮的本构或离散误差。

随后另建 `qualified_state` 批次，固定态改用前置 BGN 导数修复阶段已合格、与当前材料一致的最终保存态。新增明确的同态预检：直接读取 n/p 与 SG 重算 n/p 的相对差需不超过 1e-10，再运行 Jv。八点预检最大差均为 **0**。该步骤仅修正固定态审计输入；DC 的两条初始化路径和接受门槛不变。

首批输入、输出、失败统计均保留，见 [初始固定态结果](../../reference_tcad/simplemos_sentaurus2022/phumob_box_20260909/fixed_summary.json)、[初始探针结果](../../reference_tcad/simplemos_sentaurus2022/phumob_box_20260909/probe_summary.json)、[输入失败账本](../../reference_tcad/simplemos_sentaurus2022/phumob_box_20260909/qualified_state/initial_state_failure_ledger.csv)。此前已完成的自洽结果不受该审计选种子错误影响，因为求解器本来就按当前材料和势变量重算密度。

## 3. 合格固定态的公式、全部邻接列及探针验证

在八个合格固定态上关闭 SRH，以独立 100 位 PhuMob 公式检查单元与边迁移率，并对节点 792、1000、1009、1057 的载流子输入列检查全部受影响的相邻行。弱交叉参考直接传播本构导数与实际边通量，不用完整残差相减充当弱项真值。

| 检查 | 数量 | 最大相对差 | 门槛 | 失败 |
|---|---:|---:|---:|---:|
| 单元迁移率 / 独立公式 | 28000 | 8.23680e-16 | 1e-12 | 0 |
| 边迁移率 / 独立公式 | 29680 | 8.03666e-16 | 1e-12 | 0 |
| 弱载流子交叉矩阵项，两步长 | 864 | 5.84836e-15 | 1e-5 | 0 |
| 基态边通量与行残差闭合 | 全部采样行 | 2.14331e-16 | 1e-10 | 0 |
| 迁移率探针 / SG 装配器 | 67584 | 4.02308e-16 | 1e-12 | 0 |

864 项对应 **432 个独立非零矩阵项**及两种物理电势步长 `1e-5`、`3e-6 V`；参考单位为原始装配残差/物理 V。参考账本另记录 **240 项非零的第三顶点到对边贡献**，避免只核对边端点而漏掉单元模板。

完整 Jv 使用每隔三个自由硅节点的 ψ/φn/φp 方向，四幅度为 1e-4、3e-5、1e-5、3e-6 V。相对范数去掉探针原有的 1 下限；小三幅度、排除弱交叉源分块后，原 1e-4 正式门槛下：

| 状态集 | 正式检查 | 失败 | 最大相对差 |
|---|---:|---:|---:|
| 八个固定态，关闭 SRH | 168 | 0 | 6.53973e-7 |
| 八个候选最终自洽态，恢复 SRH | 168 | 0 | 6.53973e-7 |

新增三个 Catch2 测试还检查显式配置保护、第三顶点改变对边迁移率、独立高精度弱第三顶点列、全部矩阵列的步长差分、两载流子守恒符号及端口与残差的一致性。相关 [test_element_box_transport.cpp](../../tests/test_element_box_transport.cpp) 共 7 个测试、179 项断言通过；原有 PhuMob 的 18 个测试、171 项断言也通过。

证据：[固定态合同](../../reference_tcad/simplemos_sentaurus2022/phumob_box_20260909/qualified_state/contract.json)、[汇总](../../reference_tcad/simplemos_sentaurus2022/phumob_box_20260909/qualified_state/fixed_summary.json)、[第三顶点贡献](../../reference_tcad/simplemos_sentaurus2022/phumob_box_20260909/qualified_state/third_vertex.csv)、[探针一致性](../../reference_tcad/simplemos_sentaurus2022/phumob_box_20260909/qualified_state/probe_summary.json)、[最终态 Jv](../../reference_tcad/simplemos_sentaurus2022/phumob_box_20260909/qualified_state/post_summary.json)。这些结果证明本轮候选的离散实现及已测导数正确，不能替代原生 G 分支或弱 SRH 的独立响应校准。

## 4. 八点自洽电流比较

两组各 16 次首次尝试合格、各八个双初始化比较通过，没有失败 DC 或恢复重载。原门槛仍为全载流子逐行 1e-6、KCL/Id 1e-8、端口重放 1e-8；双初始化电势最大差 1e-6 V、密度相对差 1e-4、Id 相对差 1e-6。

32 个状态的最大 KCL/Id 为 1.82471e-14，最大端口重放相对差为 1.06581e-14。最接近逐行门槛的是 n23、Vd=1 V、Vg=0.8 V、Vela 初始化的候选态，最大行比 **7.35496e-7**，仍通过原 1e-6 门槛。该点双初始化的最大电势差约 2.93802e-8 V、密度相对差约 1.13648e-6，也通过原门槛。全局 SRH 源的 64 个分量全部低于继承合同的 1e-10 下限，全局相对源资格本轮未激活。

下表电流误差为 `(Id_Vela/Id_Sentaurus−1)×100%`，Vela 电流选取各组 Vela 初始化的合格状态；两组使用同一原生 PhuMob 参考。

| NWell | Vg (V) | Vd (V) | legacy 误差 | 单元平均候选误差 | 候选 Id (A/μm) |
|---|---:|---:|---:|---:|---:|
| n19 | 0.8 | 0.05 | −0.11722918% | +0.00330975% | 7.83593182242e-7 |
| n19 | 1.0 | 0.05 | −0.73372453% | +0.00098502% | 8.25233540089e-6 |
| n19 | 0.8 | 1 | −0.51417267% | +0.00286514% | 6.62384693815e-6 |
| n19 | 1.0 | 1 | −1.13093731% | +0.00122127% | 4.78032482991e-5 |
| n23 | 0.8 | 0.05 | −0.04299754% | +0.00569362% | 3.81020928315e-9 |
| n23 | 1.0 | 0.05 | −0.05265768% | +0.00507750% | 3.68622738139e-7 |
| n23 | 0.8 | 1 | −0.03533183% | +0.00652926% | 2.43732902010e-8 |
| n23 | 1.0 | 1 | −0.07229945% | +0.00510983% | 2.02973280853e-6 |

八点绝对相对误差全部改善，包含低 NWell 控制和高 Vd 点。高低 NWell 配对另以 `[(Id_n23/Id_n19)_Vela/(Id_n23/Id_n19)_native−1]×100%` 定义，保留原生及 Vela 的高低电流绝对增量，不把绝对误差和配对误差混为一项。

| Vg (V) | Vd (V) | legacy 配对比值误差 | 候选配对比值误差 |
|---|---:|---:|---:|
| 0.8 | 0.05 | +0.07431876% | +0.00238380% |
| 1.0 | 0.05 | +0.68610094% | +0.00409244% |
| 0.8 | 1 | +0.48131563% | +0.00366401% |
| 1.0 | 1 | +1.07074734% | +0.00388851% |

四组配对比值误差也全部改善，见 [配对账本](../../reference_tcad/simplemos_sentaurus2022/phumob_box_20260909/qualified_state/nwell_pairs.csv)。

明细：[求解汇总](../../reference_tcad/simplemos_sentaurus2022/phumob_box_20260909/qualified_state/dc_summary.json)、[电流及双初始化](../../reference_tcad/simplemos_sentaurus2022/phumob_box_20260909/qualified_state/dc_comparison.csv)、[legacy 尝试](../../reference_tcad/simplemos_sentaurus2022/phumob_box_20260909/qualified_state/legacy/attempts.csv)、[候选尝试](../../reference_tcad/simplemos_sentaurus2022/phumob_box_20260909/qualified_state/candidate/attempts.csv)、[源门槛激活](../../reference_tcad/simplemos_sentaurus2022/phumob_box_20260909/qualified_state/source_activation.csv)。

## 5. 回归与证据保存

完整 `ctest --preset windows-ucrt64-release --parallel 4` 为 **772/788 通过，479.24 秒**。16 项失败名称与前置阶段完全一致；复核各失败子测试后，仍为 13 项历史身份/替代链校验和 3 项缺失历史 `tests/test_mos_mixed_material.cpp`。没有新增失败，也没有修改旧哈希或删除失败记录；完整测试集仍不是全绿。三个新增候选测试在全量运行中再次通过。见 [回归比较](../../reference_tcad/simplemos_sentaurus2022/phumob_box_20260909/qualified_state/ctest_comparison.json) 和 [完成复核证据](../../reference_tcad/simplemos_sentaurus2022/phumob_box_20260909/qualified_state/review_evidence.json)。

本轮旧二进制和旧源码继承前置阶段的不可变快照，候选输入、初始失败、修正后的固定态与 32 次 DC 分别保存。最终源码、二进制和文档另存于浅层忽略目录 `build-release/phumob_box_snapshot_20260909/`，不覆盖历史失败证据。主脚本为 [validate_simplemos_phumob_box_20260909.py](../../scripts/validate_simplemos_phumob_box_20260909.py)，输入修正及门槛前置见 [validate_simplemos_phumob_box_qualified_20260909.py](../../scripts/validate_simplemos_phumob_box_qualified_20260909.py)。本轮未提交或推送 Git。

## 6. 结论和下一阶段

本轮证据支持：单元顶点本构与 box 单元/边平均是这八点 PhuMob legacy 电流差异的重要来源；完整载流子模板与稳定链式导数使该方案能够进行合格的自洽求解。此前仅修复 legacy 弱导数时 Id 基本不变，本轮切换离散方式后八点误差均显著降低，两次实验的作用已分开验证。

候选继续显式启用。原生 G 截断的精确形成算法尚未确定，其空穴单元 1e-7 对比资格不能被此次 Id 改善覆盖；源下限未激活的 SRH 相对闭合也仍需独立验证。下一阶段可在保持这些记录的前提下扩展低 Vg 控制，重点检查接近门槛的少数载流子态，并继续量化原生 G 分支差异，再依次恢复 Enormal 和高场饱和。本轮没有重算完整 0–1 V 曲线，也没有放行完整原始模型。
