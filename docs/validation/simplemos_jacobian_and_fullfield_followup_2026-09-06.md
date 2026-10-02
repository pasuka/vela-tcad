# SimpleMOS Jacobian 完整性审计与 0–1 V 场对照执行结果

2026-09-06。承接 [物理量与 Jacobian 复盘](simplemos_physics_and_jacobian_review_2026-09-06.md) 的后续任务。用户已授权向 Sentaurus 虚拟机上传、下载及运行 sdevice。本轮完成独立列差分审计、校准沟道方向的 Poisson 分解，以及 n19/n23 × 两种漏压 × 21 个 Vg 的场导出与比较。

**主要结果：84 对目标状态满足既定 Vela 接受门槛及两程序端口 KCL 检查；高 NWell 的平滑 Id 差仍存在。Jacobian 审计进一步定位了 Lombardi 数值偏导的步长误差，隔离步长对照通过，但弱交叉块和少数载流子场的精度资格仍不完整，因此没有进入正式算法修改或 M82/M83。**

这里的“状态通过”有明确范围：现有载流子行筛选、带 floor 的全局闭合及端口 KCL 通过，不表示所有载流子行都被检查，也不表示所有物理场均已达到同等精度。此次发现的空穴场资格缺口必须随数据一起使用。

**工况、单位与实际执行。** n19/n23 的背景 Boron 分别为 1e17/2e17 cm^-3，节点数为 1480/1482。采用原网格，Vd=0.05/1 V，Vg=0–1 V、步长 0.05 V；matched-ni、no-BGN、Boltzmann、SRH，迁移率为 PhuMob + Lombardi + 向量准费米梯度高场饱和。宽度按 1 um，电流单位 A/um。Vela 使用 UCRT64 Release 隔离程序、SparseLU；HDF5/TDR 导出已实际成功。Sentaurus 为 T-2022.03-SP2。

| 本轮工作 | 实际数量及结果 |
| --- | --- |
| 独立 Jacobian 列审计 | 15 个只读配置；每配置 18 个局部邻接节点、三种势，共 54 列；对所有残差行作三档中心差分 |
| Lombardi 步长对照 | 4 个额外只读配置；本轮共 19 次独立列审计 |
| 校准方向 Poisson 诊断 | 8 次只读残差探针；使用已封存的双幅度正负沟道守恒扰动态 |
| 原生全场扫描 | 4 条曲线，84 个目标状态，全部导出 |
| 原生精度加严 | n23 的 42 个保存态逐点重闭合，全部通过；保留原扫描 |
| Vela 全场重算 | 46 次升序尝试 + 40 次原生一致初值恢复 + 2 次收紧停止条件，共 88 次 DC 尝试，最终覆盖 84 个合格目标点 |

原生原始目录为 `/tmp/vela_simplemos_fullfield_20260906`，精度对照位于其 `precision` 子目录。两批任务均已完成，结果已取回。此次没有改动正式 `CoupledDDAssembler.cpp`、物理模型、端口电流定义或默认配置；保留了工作树已有的 runner 修改。没有提交代码或生成数据。

**严格接受与恢复过程。** Vela 原始停止条件为 reltol=1e-7、abstol=1e-12、max_iter=200，载流子行 eps=1e-6；独立读取最终态后检查 global tolerance=1e-6、source_floor=1e-10，以及 `abs(sum(Icontact))/abs(Id) <= 1e-8`。行筛选的源项和通量门槛保持原值。

升序求解在 n23 低漏压 Vg=0.1 V 和高漏压 Vg=0 V 停止，失败态没有用于后续升序续算。40 个缺失点改用各自原生一致映射态作为初值，仍求解完整 Vela 方程；其中两种漏压的 Vg=0.1 V 仍仅因 KCL 未过。独立冻结的最后两次对照将 reltol 收紧到 1e-10、abstol 收紧到 1e-14，原物理模型及所有事后门槛不变，两点均通过。

| n23、Vg=0.1 V | 原生初值恢复时 KCL/Id | 收紧 Newton 后 KCL/Id |
| --- | ---: | ---: |
| Vd=0.05 V | 1.3154e-8 | 9.6147e-12 |
| Vd=1 V | 6.6985e-8 | 5.3076e-9 |

最终 84 点 Vela 最大 KCL/Id 为 7.2666e-9。它们包含独立恢复点，不能称为四条全部成功的升序 continuation。全部状态的 global electron/hole **source-qualified 标记仍为 false**；通过的是带源项 floor 的闭合检查，不能宣称微小净 SRH 已获得百万分之一的相对闭合。

**原生精度检查。** 最初 Digits=8、ErrRef(E/H)=100 的扫描均报告原生收敛，但 n23 深关断点的 KCL/Id 最差为 0.005694，即 0.5694%。根据本地 [Sentaurus 用户手册](D:/software/ATCNIS01-202203/sdevice_ug.pdf) 第 223–225 页，另行冻结了 `ExtendedPrecision(128) Method=Super Digits=12 ErrRef(E/H)=1e-2 RhsMin=1e-20 Iterations=40` 的保存态重闭合。运行日志确认实际使用 128 bit double-double 和 Super。

42 点加严结果的最大 KCL/Id 为 **5.6557e-15**。原生 Id 最大变化为 +0.97924%，在 Vd=0.05 V、Vg=0 V；三种势的导出值最大变化只有 3.11e-15 V。这说明极小端口电流的数值精度变化可以显著超过可见节点场变化。它不等于已证明原生每个残差行完整收敛。n19 原始参考最大 KCL/Id 为 7.1331e-15，最终比较保留其原扫描，n23 则使用加严参考；两者物理模型、网格、偏压及输出定义一致。

**最新 Id–Vg 对比。** 下列数值来自新扫描和加严参考，不再引用旧 M66 曲线。误差为 `100*(Id_Vela/Id_native-1)`。

| 高 NWell 的 Vg 区间 | Vd=0.05 V：中位 / 最大误差 | Vd=1 V：中位 / 最大误差 |
| --- | ---: | ---: |
| 0–0.5 V | 9.9951% / 12.1172% | 9.0068% / 10.7936% |
| 0.55–1 V | 13.6535% / 14.3336% | 11.7866% / 12.2793% |
| 1 V 单点 | 13.4315% | 10.1108% |

低漏压最差点为 Vg=0.9 V；高漏压最差点为 Vg=0.85 V。低 NWell 在整个区间的最大误差分别为 8.9109% / 6.5217%。高 NWell 差异仍然覆盖 0.5–1 V，当前也不能把低 NWell 当作零误差基线。

![新 84 点电流、沟道势与 SRH 对照](../../reference_tcad/simplemos_sentaurus2022/fullfield_validation_20260906/native_precision/comparison/fullfield_comparison.png)

**逐 Vg 物理场差及其资格。** 在相同 Si 节点上比较，不插值，不移动电势零点。RMS 使用共同 Si 重心面积加权；浅沟道为深度不超过 0.05 um、横向 |y|<=0.125 um，n19/n23 分别 184/186 节点。下表为高 NWell 的 21 点 RMS 最小–最大包络，包含该几何区域内全部节点；它是观测到的场差，不能直接视为均已通过逐行精度检查。

| 物理量 | Vd=0.05 V | Vd=1 V |
| --- | ---: | ---: |
| 沟道 psi RMS，mV | 0.7442–1.0089 | 1.0959–1.6487 |
| 沟道 phi_n RMS，mV | 0.1429–0.2168 | 1.4097–1.5538 |
| 沟道 phi_p RMS，mV；资格不足 | 0.000363–4.3959 | 1.0818–15.3929 |
| 全 Si SRH 净率归一化 L1 差 | 2.9309%–3.1464% | 1.1834%–1.2307% |
| 全 Si 最大 psi 差，mV | 13.8075 | 13.8075 |
| 全 Si 最大 phi_n 差，mV | 0.6153 | 17.8807 |
| 全 Si 最大 phi_p 差，mV；资格不足 | 75 | 125 |

所有场、分区、峰值节点及 n/p 对数差均写入 CSV。SRH 为带符号净产生复合率；当前没有独立的总 G 和总 R 导出。以实际配置重算广义 Boltzmann SRH，与生产载流子项的 84 点核对最大相对 L1 差为 7.32e-15；单位换算和实际 all-cell 体积已独立核验。跨程序的率场比较统一使用 Si 重心面积，不冒充原生内部体积。

本轮额外按 `NewtonSolver.cpp` 的实际源项/通量筛选规则重建了每点载流子行资格，84 点的数量与生产探针完全一致。**所有 84 点的浅沟道空穴行都未进入现有行门槛检查。** 电子沟道被检查节点数随工况变化，为 5–186；也不能把整个电子场视为统一精度保证。

例如 n23、Vd=1 V、Vg=0.35 V 的节点 987：phi_p 差为 -0.125 V，原生/Vela 空穴浓度约 0.5067/0.004650 cm^-3，通量与源项尺度约 1.82e-19/6.70e-20，低于现有筛选门槛。n19 高漏压 Vg=1 V 的沟道 phi_p RMS 差达 226.63 mV，全部落在未被检查的空穴行。加严原生后这些导出场差没有消失；现有证据不能将其直接解释为真实物理模型误差。更小的 Id/KCL 误差也不能替代少数载流子场的独立精度资格。

**Jacobian 独立审计的新增定位。** 以焦点边 1091–1092、320–324 的完整三角邻接节点为列集合，逐列扰动 psi、phi_n、phi_p，比较全部三个方程块的完整残差差分。物理步长为 1e-5、1e-6、1e-7 V。按冻结规则，信号必须超过块内峰值的 1e-10，1e-6 与 1e-7 V 差分的相对变化须不超过 1e-4，Jacobian 相对误差须不超过 1e-3。

| 当前物理分支的四点审计 | 数量 |
| --- | ---: |
| 活跃“列—方程块”组合 | 509 |
| 残差差分稳定 | 476 |
| 原补项关闭时通过 | 452 |
| 向量场链式补项开启后通过 | 472 |
| 补项后仍稳定超差 | 4 |
| 活跃但差分不稳定，未认证 | 33 |

稳定组合没有检出 Jacobian 零位置上遗漏的可分辨残差信号。33 个不稳定组合分别为 phi_p→电子 17 个、phi_p→Poisson 12 个、psi→空穴 4 个；它们不能按低于有效分辨率的相对误差放行。这里覆盖的是焦点区域全部邻接列，不是整个器件每一列、所有 Vg 或每个非光滑切换点的证明。

剩余四项均为高漏压下的 psi 列，位于漏端邻域：n19 的 1088→电子行 1087、1089→空穴行 1087；n23 的 1089→电子行 1088、1090→空穴行 1088。误差为 0.2120%–0.2737%。改变向量补项步长不改变它们；移除 Lombardi 后相应误差降至约 1e-10，edge_projection 分支仍出现同类偏差。

因此另行冻结四个局部输运差分步长对照，把源码现有 `1e-6*max(1,abs(potential))` 的系数分别改为 1e-7、1e-8，仅作用于隔离副本。八个缺陷对照全部通过：相对 1e-6 V 的独立残差差分，误差不超过 6.42e-5；相对更细的 1e-7 V 差分，误差不超过 9.39e-7。物理残差与基线逐字节相同。证据支持这四项来自局部表面迁移率数值偏导的步长误差，不支持把它们解释为新的物理残差修正。

constant_field、无表面迁移率及场导数冻结开关均做了只读控制；冻结开关下补项 on/off 的结果相同。表面步长对照未解决弱交叉块的差分分辨率，也没有形成全域生产修复。正式实现仍混合手写解析项与局部数值差分，不能称为全符号解析 Jacobian。

**校准沟道方向的 Poisson/输运耦合。** 使用已校准的指定沟道连续性守恒扰动，在相同偏压下分解

`delta_psi = 电子电荷响应 + 空穴电荷响应 + 边界响应 + K^-1 delta_F_Poisson`。

这里 K 和电荷体积均为现有 Vela 算子。全幅/半幅扰动的电荷加边界响应相对实际 delta_psi 的重构误差为 2.81e-5–2.21e-4，全部低于 1e-3；全幅与两倍半幅势响应的相对差不超过 8.06e-8。最大实际势变化仅约 9.94e-9/5.68e-9 V，须与跨程序 mV 状态差区分。

以界面节点 338 为例，全幅中心响应中：

| 电势响应分量，V | 低漏压 | 高漏压 |
| --- | ---: | ---: |
| 沟道电子电荷 | -2.6415e-9 | -1.0050e-9 |
| 衬底空穴电荷 | +5.7149e-10 | +3.8973e-10 |
| 实际总 delta_psi | -2.0310e-9 | -5.7326e-10 |

这说明该校准方向中衬底空穴反馈会部分抵消沟道电子的势响应；它不能直接证明跨程序差异由某一电荷项造成。

对严格 Vela 态减去一致映射的原生态作同样分解时，忽略算子残差的相对重构误差达到 **5.7126 / 4.7520**。节点 338 的实际势差为 -1.0328/-0.8894 mV，而算子残差对应项为 -6.3701/-2.3446 mV，其余电荷响应互相抵消。加入残差项后代数重构闭合。这是 Vela 算子在映射态上的账本，不是 Sentaurus 内部算子；不能据此把某一分量直接命名为原生离散错误。

**验证与后续边界。** 两个隔离 C++ 审计程序编译成功；本轮新增数值证据回归 [test_simplemos_jacobian_fullfield_followup.py](../../tests/regression/test_simplemos_jacobian_fullfield_followup.py) 的 10 项检查全部通过，相关 Python 脚本语法检查通过，图已检查。由于正式求解器没有改动，没有为此次诊断重跑整个 CTest。

接下来应先处理本轮已明确的两个资格缺口：以避免残差消减的独立方法复核 33 个弱交叉块；为少数载流子场另行定义并冻结逐行及初始化不变性验证，保持本轮门槛和失败记录不变。在这些资格明确后，再把向量场链式导数与表面数值偏导步长策略整理为生产候选并作回归。物理差异定位仍应以校准沟道方向为依据，对局部介电耦合、源项体积和界面势形成做定义清晰的同扰动对照。当前没有候选同时改善高 NWell 绝对误差、高低配对及高漏压反例，M82/M83 继续保留。

**证据入口。** 最终比较以加严版本为准；原始扫描和所有恢复失败均保留。

- [最终 84 点电流及资格](../../reference_tcad/simplemos_sentaurus2022/fullfield_validation_20260906/native_precision/comparison/points.csv)、[物理场统计](../../reference_tcad/simplemos_sentaurus2022/fullfield_validation_20260906/native_precision/comparison/fields.csv)、[SRH](../../reference_tcad/simplemos_sentaurus2022/fullfield_validation_20260906/native_precision/comparison/srh.csv)、[载流子行可见性](../../reference_tcad/simplemos_sentaurus2022/fullfield_validation_20260906/native_precision/comparison/visibility/fields.csv)。
- [Jacobian 列审计](../../reference_tcad/simplemos_sentaurus2022/jacobian_columns_20260906/result.json)、[表面步长缺陷对照](../../reference_tcad/simplemos_sentaurus2022/jacobian_columns_20260906/surface_step/defect_controls.csv)、[Poisson 耦合检查](../../reference_tcad/simplemos_sentaurus2022/channel_poisson_coupling_20260906/checks.csv)。
- [原生精度稳定性](../../reference_tcad/simplemos_sentaurus2022/fullfield_validation_20260906/native_precision/precision_result.json)、[统一封存清单](../../reference_tcad/simplemos_sentaurus2022/fullfield_validation_20260906/followup_evidence.json)。
- 完整节点原值、原始 TDR、重启态和日志存于忽略目录 `build-release/simplemos_fullfield_20260906/`；最终节点 CSV 为其中 `native_precision/comparison_nodes.csv`。完整配置、源码及文件 SHA-256 由上述 manifest 追溯。
