# SimpleMOS M45-M56 对比进展与 Sentaurus 导出量汇总

## 结论摘要

截至 M56，SimpleMOS 官方原始 deck 的默认端口观测仍保持 M8/M46 基线：8 个器件、2 个漏压、16 条 Id-Vg 曲线、816 个直接偏压点全部通过，M46 的 16 份逐点比较 CSV 与 M8 位级一致。默认基线的最大 Vela-Sentaurus 差异仍位于 n23、Vd=0.05 V、Vg=0.05 V，为 0.109418680924 dex（Vela 高 28.653%）。

M47-M52 已把该目标点的主要异常从“求解状态差异”收缩到“Sentaurus 默认 substrate 端口电流的观测定义”：

- 电势势垒、电子浓度、SRH 和准费米差异在相邻栅压上平滑变化，不能复现 n23、Vg=0.05 V 的尖峰。
- 四端连续性账本把目标漏电流差异定位到 substrate 电子电流分配，并由 source 反向补偿；SRH 体源项不是主导项。
- Sentaurus 原生接触面 `eCurrentDensity · ContactSurfaceNormal` 积分与独立节点场边界积分逐点一致，但与默认 substrate `eCurrent` 明显不同。
- 仅切换 Sentaurus 的 `CurrentWeighting` 或 `DirectCurrent` 时，六状态求解场逐值不变；两种端口算法在六状态内最多相差 9.260e-31 A/um，并消除目标点默认端口与原生面通量差异的 99.999854401%。
- M53 把 `DirectCurrent` 扩展到完整16曲线后否证了“可直接替换默认 Id-Vg 观测”的假设：全局最大误差增至 2.126075 dex，只有3/16条曲线保持原资格门槛通过；72个 M52 共享端口分量仍逐值完全回放。
- M54 表明 M53 目标点和最差点的漏端变化分别有93.38%和99.40%来自source/drain共模，但全816状态共模分数中位数仅43.08%，因此完整矩阵属于混合端口观测变化。
- M56 证明本算例的Si/SiO2界面使用共享几何节点/边和分区域场记录：零坐标重复节点、界面电势严格连续、载流子密度/电流只支持在Silicon；gate电压通过附着在Oxide边界上的电极施加。

因此，现有证据不支持重新调节 HFS、SG、BGN、SRH、迁移率、接触模型或网格。剩余约 2.876e-4 的 Direct/Weighted 对原生 Tcl 面积分相对差，已限定为内部离散接触积分与绘图节点场后处理积分之间的定义差异，而不是自洽状态变化。

## 算例范围与比较矩阵

原始 M8/M46 矩阵固定 Lg=0.25 um，覆盖以下三因子组合：

| 因子 | 水平 |
|---|---|
| NWell | 1e17、2e17 cm^-3 |
| GOxTime | 10、15 min |
| LDD Dose | 1e14、2e14 cm^-2 |
| 漏压 | 0.05、1.0 V |
| 栅压格点 | 0.00-2.50 V，步长 0.05 V，共 51 点 |

八个器件节点为 n17-n24。M47-M52 的归因矩阵固定 Vd=0.05 V，使用 n23（NWell=2e17 cm^-3、GOxTime=15 min、LDD=1e14 cm^-2）作为目标器件，使用同 GOxTime/LDD 的 n19（NWell=1e17 cm^-3）作为低 NWell 控制，并精确比较 Vg=0.00、0.05、0.10 V 六个状态。

## 阶段进展

| 阶段 | 目的 | 关键结果 | 当前状态 |
|---|---|---|---|
| M45 | 对 M43 SG 核统一和 M44 准费米参考/增量保存后的历史诊断重新定基 | 准费米补偿切线/报告最大差 2.314e-15 dex；旧混合算子失败证据保留为历史结论 | 闭环 |
| M46 | 重跑完整生产矩阵 | 16/16 曲线通过，816 点；所有比较 CSV 与 M8 位级一致；最大误差 0.109419 dex | 闭环，生产矩阵未改变 |
| M47 | 默认 BGN-on 自洽状态归因 | 势垒、浓度、SRH、准费米量均未满足目标局部化规则；状态差异是平滑背景 | 未归因于局部状态异常 |
| M48 | 四端分配与连续性闭合 | substrate 重分配为漏电流差的 189.2%，source 为 -101.4%，KCL 差为 12.2%；SRH 仅 1.57% | 定位到 substrate/source 电子分配 |
| M49 | substrate 边界电子输运因子分解 | 99.994% 的目标端口差留在 `terminal_boundary_residual`；节点场边界积分不携带尖峰 | `boundary_observable_limited` |
| M50 | Sentaurus 原生 substrate 面电子通量导出 | 原生 Tcl 面积分与 M49 独立边界积分最大差 4.333e-34 A/um，但与默认端口值不一致 | `native_observable_mismatch` |
| M51 | `CurrentWeighting` 单变量回放 | 六状态场完全不变；目标端口变为 -2.376113673567e-19 A/um，消除默认-原生差的 99.999854401% | `weighted_third_observable` |
| M52 | `DirectCurrent` 单变量回放 | 与 `CurrentWeighting` 在舍入尺度一致；六状态场完全不变 | `direct_matches_weighted` |
| M53 | `DirectCurrent` 完整矩阵回放 | 16曲线、816点完整；全局误差由0.109419增至2.126075 dex，原门槛3/16通过 | `direct_increases_global_error` |
| M54 | 完整矩阵四端共模/反对称分解 | 816状态、9792端口分量；目标共模93.38%，全矩阵中位43.08% | `direct_shift_mixed` |
| M56 | Si/SiO2双节点、场支持和接触拓扑审计 | 8个TDR零重复坐标；20条共享边、21个界面节点；六状态电势差为0 | `shared_topology_continuous_potential_insulator` |

## Sentaurus 脚本与导出量

### 固定物理与偏压路径

M47-M53 共用 T-2022.03-SP2 的同一输入 TDR、网格和默认 BGN-on 物理：

- `EffectiveIntrinsicDensity(OldSlotboom)`；
- Silicon `Mobility(PhuMob HighFieldSaturation Enormal)`；
- Silicon `Recombination(SRH(DopingDependence))`；
- 先求 Poisson，再求 Poisson/Electron/Hole；漏极从 0 V 连续升至 0.05 V，随后栅极从 0 V 连续扫至 2.5 V；
- 默认 `Math` 为 `Extrapolate Iterations=20 ExitOnFailure`。M51 只增加 `CurrentWeighting`，M52 只增加与其互斥的 `DirectCurrent`。

四个电极为 source、drain、gate、substrate，初始均为 0 V。栅压施加在 gate 电极/接触上，不是把 SiO2 当作导电区施压。SiO2 中仍求解和导出静电势、电场，但 Electron/Hole 输运量只在 Silicon 半导体区比较。

### TDR 场量

`Plot` 块导出的关键场量如下：

| 类别 | Sentaurus 量 | 用途与当前比较 |
|---|---|---|
| 静电 | `Potential`、`ElectricField/Vector`、`SpaceCharge` | 势垒、电势分布和界面电场；M47 六状态共节点直接比较 |
| 载流子 | `eDensity`、`hDensity` | 势垒点、源势垒支撑区和全场浓度差；M47 目标点未呈现统一浓度偏移 |
| 电流 | `TotalCurrent/Vector`、`eCurrent/Vector`、`hCurrent/Vector` | Silicon 节点电流密度与源/漏/substrate 分量；M48-M50 用于连续性和边界积分 |
| 输运 | `eMobility`、`hMobility`、`eVelocity`、`hVelocity` | M49 的密度-迁移率-准费米梯度精确桥接 |
| 准费米 | `eQuasiFermi`、`hQuasiFermi`、`eGradQuasiFermi/Vector`、`hGradQuasiFermi/Vector` | M47 准费米状态和 M49 substrate 邻层驱动力比较 |
| 场方向分量 | `eEparallel`、`hEparallel`、`eENormal`、`hENormal` | 高场和界面法向驱动力诊断；未用于重新调参 |
| 复合 | `SRH` | M47 空间相关与积分，M48 电子/空穴连续性账本 |
| 掺杂 | `Doping`、`DonorConcentration`、`AcceptorConcentration` | 固定 TDR 掺杂和 NWell 控制识别 |
| 能带与 BGN | `BandGap`、`BandGapNarrowing`、`Affinity`、`ConductionBand`、`ValenceBand` | 默认 BGN-on 势垒、能带和参考坐标归因 |

栅扫的 `CurrentPlot(Time=(Range=(0 1) Intervals=50))` 导出完整 51 点端口曲线；`Plot(... Time=(0;0.02;0.04))` 在 0-2.5 V 栅扫中保存 Vg=0.00、0.05、0.10 V 三个状态 TDR。M47-M50 的共节点分析不做插值。

M51 和 M52 各对六状态的 17 个导出文件进行逐值不变性检查，共 102 个状态-字段文件：四个材料区域的静电势和电子/空穴准费米势，以及 Silicon 的电子/空穴浓度、电子/空穴电流密度和 SRH。两次检查的节点、坐标和全部数值均完全一致，最大绝对差和最大对称相对差均为 0。

### Si/SiO2界面与gate边界

M56对8个输入TDR和n19/n23六个求解状态进行了只读审计。所有TDR均以相同全局节点ID连接Si和SiO2界面，没有为同一坐标创建两个全局节点；每个器件有20条共享界面边和21个界面节点。求解结果在同一界面节点上分别保存Silicon和Oxide区域电势记录，六状态最大跨区域差为0 V。

`eDensity`、`hDensity`、`eCurrentDensity`和`hCurrentDensity`只存在于`Silicon_1`。gate接触附着于`Oxide_1`，source、drain、substrate附着于`Silicon_1`。所以SiO2绝缘表示不在氧化层求解电子/空穴漂移扩散，并不妨碍在氧化层外边界的gate电极施加静电势。默认deck不需要显式双节点、`HeteroInterface`、`Thermionic`或`Discontinuity`设置来建立本算例的标准Si/SiO2静电耦合。

### 端口量与原生接触面通量

`.plt` 端口文件提供四个电极的 `OuterVoltage`、`TotalCurrent`、`eCurrent` 和 `hCurrent`，用于 Id-Vg、四端 KCL 及电子/空穴连续性账本。

M50-M52 额外加入 `CurrentPlot Tcl`，在 `substrate` 接触域积分：

```text
sum_d(eCurrentDensity[d] * ContactSurfaceNormal[d])
Operation = Integrate Contact="substrate" IntegrationUnit=cm
```

该量保留原始法向符号，不做逐状态符号或幅值拟合。M50 初次运行中 `Contact` 与 `Region` 同时指定被 T-2022.03-SP2 拒绝；执行勘误只删除冗余 `Region=Silicon_1`，没有改变接触域、被积函数或器件状态。

## 目标点定量对比

### 自洽状态层

| 指标 | n23，Vd=0.05 V，Vg=0.05 V |
|---|---:|
| Sentaurus drain Id | 8.586250997030e-17 A/um |
| Vela drain Id | 1.104643803536e-16 A/um |
| Id 差异 | 0.109418680924 dex，Vela 高 28.653% |
| Vela 相对 Sentaurus 的准费米参考源势垒变化 | -1.291129 mV |
| 势垒 Maxwell-Boltzmann 电流代理 | 0.021690 dex，占观测对数差 19.8% |
| 势垒代理后未归因量 | 0.087729 dex |
| 积分 SRH 幅值比 | -0.009905 dex |
| SRH 空间 Pearson 相关 | 0.999817614 |
| 电子/空穴准费米参考加增量往返误差 | 6.939e-18 V |

n23 相邻栅压的电流误差为 0.017239、0.109419、0.039123 dex，而势垒代理为单调的 0.020488、0.021690、0.022607 dex。低 NWell n19 的三个电流误差为平滑的 0.023190、0.025249、0.026272 dex。该控制关系是排除“目标点局部 BGN-on 自洽状态异常”的核心证据。

### 端口观测层

| substrate 电子电流观测 | 目标值（A/um） | 与原生面通量的关系 |
|---|---:|---|
| Sentaurus 默认端口 `eCurrent` | 4.670984402844e-17 | 明显不一致，携带目标尖峰 |
| Tcl 原生法向面通量 | -2.376797223030e-19 | 与 M49 独立节点边界积分逐点一致 |
| `CurrentWeighting` | -2.376113673567e-19 | 残差 6.835494626e-23，相对 2.875927e-4 |
| `DirectCurrent` | -2.376113673558e-19 | 与 Weighted 差 9.260e-31；相对原生残差 6.835494719e-23 |

M48 的同一点四端恒等分解为：substrate 重分配 +4.654871563240e-17 A/um、source 重分配 -2.495023008915e-17 A/um、KCL 残差差 +3.003384840028e-18 A/um；积分 SRH 差仅为漏电流差的 1.57%。这些项是精确账本归因，不应被解释为把全部 0.087729 dex 对数余量线性分配给 substrate。

## NWell 趋势的当前解释

在 M46 的四组相同 GOxTime/LDD 配对中，将 NWell 从 1e17 提高到 2e17 cm^-3 后，Vd=0.05 V 和 1.0 V 的每一对曲线最大误差都增加；低 NWell 组全局最大值为 0.045587 dex，高 NWell 组达到 0.109419 dex。

M47-M52 支持的局部结论仍成立：更高 NWell 改变了与 substrate 接触关联的掺杂阱及其电子电流分配；n23 在 Vg=0.05 V 的默认 substrate `eCurrent` 尖峰不是 SRH、势垒或 BGN 自洽场尖峰，`CurrentWeighting`/`DirectCurrent` 可在六状态内消除该 substrate 观测异常且不改变求解场。

M53 证明该局部结论不能外推为“完整 drain Id-Vg 应统一采用 `DirectCurrent`”。完整矩阵中全局最大误差由 0.109419 dex 增至 2.126075 dex，n23 目标点由 0.109419 增至 1.109761 dex；13/16 条曲线未通过原资格门槛。八组 NWell-漏压配对中，`DirectCurrent` 仅有1组没有增加高 NWell 相对低 NWell 的最大误差放大。因此，默认端口算法确实参与了 substrate 分配异常，但整体 NWell 趋势不是通过全局切换 drain `DirectCurrent` 就能闭合的单一观测问题。

## 已闭环、未闭环与使用边界

已闭环：

- SG 电流核一致性；
- Vela 通用接触电流提取；
- 准费米参考/增量冻结打包；
- M46 对 M8 的 16 曲线逐点一致性；
- 目标六状态的 BGN-on 状态归因；
- substrate 原生面通量可观测性；
- `CurrentWeighting` 与 `DirectCurrent` 的端口算法归因。

仍保留的边界：

- M53 已完成 `DirectCurrent` 的完整16曲线诊断，但结果否证了其作为官方默认 Id-Vg 基线替代项；
- M54将目标/最差点定位为强共模变化，但完整矩阵仍是共模与反对称项混合，不能采用单一全局修正；
- M56只证明当前SimpleMOS的标准Si/SiO2表示，不外推到需要能带不连续、热发射、隧穿或界面陷阱的其他异质结；
- Direct/Weighted 与 Tcl 节点场面积分之间仍有约 2.876e-4 的相对离散定义差；
- 现有证据定位了 Sentaurus 端口观测差异，但不据此声明 Sentaurus 或 Vela 存在生产代码缺陷；
- 所有生产默认物理、网格和接触设置均未修改。

## 证据索引

- M45：`docs/validation/simplemos_m45_post_qf_rebaseline_2026-09-01.md`
- M46：`docs/validation/simplemos_m46_full_matrix_requalification_2026-09-01.md`
- M47：`docs/validation/simplemos_m47_default_bgn_self_consistent_attribution_2026-09-01.md`
- M48：`docs/validation/simplemos_m48_terminal_partition_continuity_closure_2026-09-01.md`
- M49：`docs/validation/simplemos_m49_substrate_electron_transport_attribution_2026-09-01.md`
- M50：`docs/validation/simplemos_m50_native_substrate_face_flux_export_2026-09-01.md`
- M51：`docs/validation/simplemos_m51_current_weighting_attribution_2026-09-01.md`
- M52：`docs/validation/simplemos_m52_direct_current_attribution_2026-09-02.md`
- M53：`docs/validation/simplemos_m53_direct_current_full_matrix_2026-09-02.md`
- M54：`docs/validation/simplemos_m54_terminal_common_mode_attribution_2026-09-02.md`
- M56：`docs/validation/simplemos_m56_si_oxide_interface_topology_2026-09-02.md`
- M46 机器报告：`reference_tcad/simplemos_sentaurus2022/full_matrix_requalification/m46_full_matrix_requalification_report.json`
- M47-M56 机器报告：相关目录位于 `reference_tcad/simplemos_sentaurus2022/{default_bgn_state_attribution,terminal_partition_continuity_closure,substrate_electron_transport_attribution,native_substrate_face_flux_export,current_weighting_attribution,direct_current_attribution,direct_current_full_matrix,terminal_common_mode_attribution,si_oxide_interface_topology}/`。
