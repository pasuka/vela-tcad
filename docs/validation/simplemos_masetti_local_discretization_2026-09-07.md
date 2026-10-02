# SimpleMOS 简化模型的局部迁移率、SG 与 Poisson 耦合审计

2026-09-07。本轮已完成八个合格 Masetti 状态的局部审计、两份 Sentaurus 原生单元/box 导出，以及四工况共 20 次局部守恒扰动自洽验证。**相同导入状态上的端口电流只差约 +0.070%，高 NWell 自洽电流仍差 +9.81%～+13.70%；目前证据将重点指向状态形成及内部输运离散，尚未确定唯一根因。** 沟道方向响应已校准，正式算法与接受门槛没有改动。

本轮承接[简化模型对照](simplemos_reference_subset_comparison_2026-09-06.md)。八点范围是 n19/n23 × Vd=0.05/1 V × Vg=0.8/1 V，不代表整个 0–1 V 曲线；新扰动范围进一步限定为 Vg=1 V 的四点。n19/n23 背景 Boron 为 1e17/2e17 cm⁻³，节点数 1480/1482，三角形 2742/2746；Si 均为 1750 个三角形、942 个节点。材料为 Si、SiO₂ 和两个 Nitride 侧墙。

## 1. 状态差与端口公式差已经分开

`mapped` 保留原生电势、电子/空穴准费米势，以冻结的 Vela ni/热电压生成一致密度；`strict` 为上一轮通过全部门槛的 Vela 自洽状态。下列固定状态数值是同一 Vela 守恒端口泛函作用于 `mapped`，不是原生内部边电导的直接导出。

| n23 工况：Vd / Vg，V | 自洽 Id 相对原生误差 | mapped 端口泛函相对原生误差 | 沟道 SG κ：strict / mapped |
| --- | ---: | ---: | ---: |
| 0.05 / 0.8 | +13.6995% | +0.070190% | 0.979845 |
| 0.05 / 1.0 | +12.9842% | +0.070188% | 0.960219 |
| 1.0 / 0.8 | +11.8812% | +0.070472% | 0.969001 |
| 1.0 / 1.0 | +9.8124% | +0.070283% | 0.937040 |

八点 mapped 端口误差范围为 +0.070011%～+0.070472%。16 次边通量到端口求和恒等式最大相对差 4.44e-16。这个结果说明，单纯更换输出电流重建方式没有足够证据解释约 10% 的差异；**它没有排除内部迁移率/SG 离散通过自洽状态影响端口电流。**

以高漏压反例 n23、Vd=Vg=1 V 为例：

| 节点 | 位置 | strict−mapped 电势 | strict−mapped 电子准费米势 | 电子密度差，dex |
| --- | --- | ---: | ---: | ---: |
| 320 | 浅沟道，x=0.00776362 μm | −1.259281 mV | +0.615651 mV | −0.0314975 |
| 324 | 同深度相邻沟道点 | −0.303104 mV | +1.293877 mV | −0.0268281 |
| 338 | Si/氧化层界面，x≈0.00504594 μm | −1.376614 mV | +0.616631 mV | −0.0334850 |

沟道边 320–324 的电子迁移率为 729.636313 cm²/(V·s)，strict/mapped 完全相同；准费米势降从 16.317391 mV 变为 16.995617 mV，增加 4.15646%；共同几何下的 SG κ 减少 6.29600%，两项合成的该边通量反而减少 2.40123%。漏接触邻边 1091–1092 的 κ 基本不变，势降/通量增加约 10.098%，因此“最大漏端流出变化”和“形成差异的局部原因”必须区分。

κ 使用共同 Vela ni、热电压、几何与 SG 对数均值定义，不能标成实际原生边 κ。跨 n19/n23 对照按坐标核对：沟道节点相同，界面对应点坐标相差约 9.23e-10 μm；仅跨器件对应使用 1e-8 μm 距离检查，同网格原生导出的 1e-12 μm 门槛未改变。

## 2. Masetti 参数相同，单元平均仍未对齐

在虚拟机独立目录执行 `sdevice -P:Silicon`，取得当前 T-2022.03-SP2 硅材料默认参数。电子/空穴均选择 Formula 1，18 个参数与 Vela 的 Masetti 参数逐项完全相同。本轮 300 K、冻结掺杂下，μ 随空间变化，但对 Newton 状态 ψ、φn、φp 的导数为零。

2022 用户手册 Mobility Averaging（p464）说明先计算每个半导体单元各顶点的迁移率，再形成单元或单元边平均；这不等于节点绘图字段。新增原生 `eMobility/Element`、`hMobility/Element` 导出保持两工况电流与冻结值逐位一致，坐标完全相同；所有材料的 `BM_ElementVolume` 与各自 cell_id 的坐标三角形面积通过逐单元检查。

然而，以下简单公式仍不能逐单元复现原生导出：

| 单元公式候选 | n19 电子最大差 | n23 电子最大差 | n19 空穴最大差 | n23 空穴最大差 |
| --- | ---: | ---: | ---: | ---: |
| 三顶点局部 μ 的算术平均 | 14.6208% | 13.0244% | 11.1011% | 10.1247% |
| 将三顶点平均掺杂代入 μ 公式 | 48.4817% | 44.6495% | 39.8551% | 37.0316% |
| 三顶点原生节点绘图 μ 的平均 | 44.5066% | 38.7943% | 29.3617% | 26.1455% |

另外检查了局部 μ 的调和及几何平均；所有候选均未通过逐单元最大相对差≤1e-8 的等价性检查，失败数据保留。面积/编号检查不能单独认证另一字段的内部含义。**不能将上述最大差当作 Id 贡献，也不能把“μ(平均掺杂)”开关直接称为原生单元平均修复。** 原生实际边 SG 电导、单元内取值及绘图映射仍需进一步对齐。

## 3. 原生 Poisson 表达式已可重建

两份原生运行均在独立空目录添加 `BoxMeasureFromFile(GrdNumbering)`，取得程序自身生成的 `MeasureCoefficients.debug`，没有输入人为修改的 box 文件。用已独立核对的 40 个 Si/SiO₂ 界面三角形确认局部编号，系数最大误差≤5.12e-13，体积映射误差≤4.35e-21 μm²。

其余单元保留原生系数：n19/n23 分别有 128/124 个单元与未经原生处理的逐单元 cotangent/体积分配不同。原生日志记录 `CVPL_AverageBoxMethod` 和钝角处理，不能把这种差异默认为编号错误，也不能只依赖“非 Delaunay 元素数为零”推断逐单元等价。按节点累加后的原生 Si 体积与 signed Si 体积一致到约 6.31e-18 μm²；Vela 当前使用的全材料重心体积是另一种量。原生 K 与 Vela legacy K 的相对 Frobenius 差约 0.1134%/0.1137%，该全局范数不是局部电流误差界。

使用原生 ψ/n/p/掺杂、导出的系数和 Si 节点体积，重建自由行

`R_native = K_native ψ + q_native (n − p − ND + NA) V_Si`。

起初套用 Vela 的 q、ε0 得到残差/电荷范数约 8.21426e-6。改用 2022 用户手册内置 PMI 示例给出的 q=1.602192e-19 C、ε0=8.8542e-12 F/m，八点该比值为 6.64e-14～1.75e-13；**常数直接来自手册，没有按残差拟合。** 这是对原生 Poisson 表达式的数值复核，原生逐行残差本身没有导出。

再将 mapped 状态代入 Vela Poisson，把原生重建余项、介电几何、电荷体积、物理常数和密度约定分别记账，对实际 Vela 自由行残差的相对重建余项≤4.27e-9。通过固定 K 的条件求解重建 strict−mapped 电势，八点最大相对余项为 1.03e-7，低于冻结的 1e-5 门槛。

n23、Vd=Vg=1 V，界面节点 338 的条件电势响应：电荷体积项 −8.42273 mV，介电几何项 +0.69401 mV，密度约定项 +0.41027 mV，q/ε0 约定合计约 +0.00423 mV。它们构成 mapped 状态的 Poisson 不闭合；载流子自洽重分布会抵消其中相当部分，最终该点实际电势差仅 −1.37661 mV。这些是保留边界和余项的条件分解，不是逐项自洽修复效果。

## 4. 沟道守恒扰动的四工况校准通过

复用已核对的常量电子边通量注入实现，在自由节点 320、324 施加等大反向的连续性源。重新计算当前 Masetti 的完整 DD 端口伴随，包含全部耦合块；该自由边没有直接接触电流项。注入幅度在运行前冻结为各基线 Id 的 0.05%，另测半幅、正负和零控制。

四点共 20 次 DC 全部合格。每态 1814 个有效载流子行，无零尺度排除；最大行比值 8.56e-9，小于 1e-6；KCL/Id 最大 1.02e-14，小于 1e-8。原 global 检查全部通过，其 source_floor 资格语义保持不变，不将其误写成任意小净 SRH 源都达到相对 1e-6 闭合。

| Vg=1 V 工况 | 预测 ΔId / 注入电流 | 双幅度最大预测相对误差 |
| --- | ---: | ---: |
| n19，Vd=0.05 V | 0.134947 | 8.45e-11 |
| n19，Vd=1 V | 0.207788 | 4.46e-10 |
| n23，Vd=0.05 V | 0.169265 | 1.24e-11 |
| n23，Vd=1 V | 0.329813 | 3.18e-11 |

全部 8 个幅度结果通过预先冻结的预测误差≤0.1%、双幅度线性差≤0.1%、偶分量/奇分量≤1%、符号、零控制漂移及信噪比检查。最坏双幅度差约 3.32e-10；最小信号/零漂移约 2.47e7。此校准限于上述 **Vela 沟道守恒方向**，不外推到原生连续性源、任意分布源或有限幅度模型替换。

## 5. 线性筛查提示相互抵消，不能单项宣判根因

以下是统一自由行支持上的 `−λᵀ ΔR` 筛查，单位 A/μm；Dirichlet 替换行单列，保持总和一致。这些不同源尚未逐一自洽校准。

| Vg=1 V 工况 | 实际 strict−mapped Id | 电子输运残差投影 | 电荷体积投影 | 介电几何投影 |
| --- | ---: | ---: | ---: | ---: |
| n19，Vd=0.05 V | +1.39635e-7 | +5.76007e-7 | −4.36820e-7 | +1.43981e-9 |
| n19，Vd=1 V | +5.25270e-7 | +2.76182e-6 | −2.29661e-6 | +2.93836e-8 |
| n23，Vd=0.05 V | +2.69539e-8 | +2.75097e-8 | −1.20591e-9 | +3.29584e-10 |
| n23，Vd=1 V | +1.20908e-7 | +1.50626e-7 | −3.52133e-8 | +3.26562e-9 |

电子输运残差是主要正向项，电荷体积项方向相反，且高低 NWell 的抵消比例不同。全项线性预测的有限状态差余项占实际状态电流差，在上述四点分别约 −0.816%、+5.65%、+1.13%、+1.77%；这些余项保留，没有用沟道小扰动的通过结果替代有限差验证。因此，局部电势对体积敏感并不能推出单独换体积将改善当前 Id 或高低配对。

后续应首先明确原生单元内迁移率与边通量取值，构造同状态、相同单位及节点支持的保守输运差；再对该明确的分布源单独做正负、双幅度自洽响应校准，并与原生 Poisson 几何/体积方向交叉检查。候选必须同时解释高低 NWell 和高 Vd 反例后，才进入模型替换及正式 A/B。本轮不放行 M82/M83，也没有修改正式端口电流算法。

## 证据、复现与质量记录

Vela 使用隔离 Release runner、Eigen SparseLU、既有 4 次线性迭代修正（100 位残差累加；矩阵和状态仍为 double）；300 K、matched-ni=1.0750038488844236e10 cm⁻³、Boltzmann/no-BGN、Masetti 总杂质及原掺杂相关 SRH。未开启 HFS、表面迁移率、Auger、雪崩、DG。保持 max_iter=200、reltol=1e-7、abstol=1e-12、ψ/准费米截断 0.35/0.025 V、原线搜索、逐行和事后 global 门槛。电流为 A/μm、宽度 1 μm；坐标 μm，状态载流子 m⁻³，材料/配置掺杂按 unit_scaling 使用 cm⁻³。

原生为 T-2022.03-SP2、128 位 Super，远端目录 `/tmp/vela_simplemos_masetti_local_20260907`。两次补充 DC 和参数导出均完成并取回。实际 HDF5 导入成功。本轮执行 56 次局部只读探针、8 次扰动预检、20 次 DC 和 20 次全行事后探针；既有守恒扰动单位、方向及接触反馈 5 项单元测试通过。未新增 C++ 核心改动，无需全量 CTest；工作树原有 runner 的 127 行新增/1 行删除保持不变。

两项审计更正保留为证据：初版要求所有锐角单元等于未修正 cotangent 的假设失败，改用已有独立界面编号见证并逐项保存其余差异；初版物理 Poisson 账本的接触行存在大项相消，其分项不适合归因，最终使用自由行账本并单列实际 Dirichlet 行。原始文件没有覆盖。扰动 prepare 的首次哈希冻结因外部编译器路径失败，在任何仿真执行前单独记录工具链哈希并完成冻结，输入配置和门槛未变。

- [八点端口/状态账本](../../reference_tcad/simplemos_sentaurus2022/masetti_local_discretization_20260907/case_ledger.csv)、[局部状态](../../reference_tcad/simplemos_sentaurus2022/masetti_local_discretization_20260907/local_nodes.csv)、[SG 对照](../../reference_tcad/simplemos_sentaurus2022/masetti_local_discretization_20260907/local_sg.csv)。
- [18 参数核对](../../reference_tcad/simplemos_sentaurus2022/masetti_local_discretization_20260907/native_parameter_check.csv)、[单元迁移率失败候选](../../reference_tcad/simplemos_sentaurus2022/masetti_local_discretization_20260907/element_mobility_candidates.csv)、[单元支持检查](../../reference_tcad/simplemos_sentaurus2022/masetti_local_discretization_20260907/native_element_support_check.csv)。
- [原生几何](../../reference_tcad/simplemos_sentaurus2022/masetti_local_discretization_20260907/native_geometry.csv)、[物理常数复算](../../reference_tcad/simplemos_sentaurus2022/masetti_local_discretization_20260907/native_poisson_constant_replay.csv)、[最终自由行筛查账本](../../reference_tcad/simplemos_sentaurus2022/masetti_local_discretization_20260907/poisson_transport_free_row_screening.csv)及[解释边界](../../reference_tcad/simplemos_sentaurus2022/masetti_local_discretization_20260907/free_row_screening_scope.json)。
- [扰动合同](../../reference_tcad/simplemos_sentaurus2022/masetti_local_discretization_20260907/channel_response/contract.json)、[20 次 DC](../../reference_tcad/simplemos_sentaurus2022/masetti_local_discretization_20260907/channel_response/dc.csv)、[8 组校准](../../reference_tcad/simplemos_sentaurus2022/masetti_local_discretization_20260907/channel_response/calibration.csv)。
- [独立资格复核](../../reference_tcad/simplemos_sentaurus2022/masetti_local_discretization_20260907/independent_qualification.csv)、[最终封存](../../reference_tcad/simplemos_sentaurus2022/masetti_local_discretization_20260907/validation_evidence.json)。原始网格、状态、程序、参数和日志保留于忽略目录 `build-release/simplemos_masetti_local_20260907/`；旧证据未改写。
