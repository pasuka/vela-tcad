# SimpleMOS M45-M58 外部分析交接：NWell 误差增长尚无统一的冻结结构解释

## 交接目的

本文面向后续接手分析的其他大模型或研究人员。目标不是重新浏览全部历史，而是在保持已闭环主题不动的前提下回答两个问题：

1. 为什么八个匹配工况中，提高 Workbench 参数 `NWell` 后，Vela-Sentaurus Id-Vg 最大误差全部增大？
2. 什么最小判别实验能够区分“掺杂补偿/网格交互”“Sentaurus 默认端口观测语义”和“尚未对齐的离散物理”三类剩余解释？

当前最强结论是：n23 的局部近完全补偿、单节点浮置同号分量和零等值线贴近网格顶点都真实存在，但这些结构特征只在低 LDD 组合系统出现，不能解释全部八个 NWell/漏压误差增长行。默认 Sentaurus 端口电流算法同时存在无法由导出节点场逐项重构的接触关联阱分配语义。现有证据还不足以在两者之间作统一因果归因。

## 必须先纠正的术语

本算例的 Workbench 参数名为 `NWell`，但原始工艺 deck 用它设置：

```text
init concentration=<NWell> field=Boron
```

因此它在本算例中表示全局 p 型 body/substrate 背景硼浓度，不是传统意义上的掩膜 n-well 注入。下文保留 `NWell` 作为参数名，但“提高 NWell”应理解为背景活性硼由 `1e17` 增至 `2e17 cm^-3`。

工艺脚本锚点：`build-release/reference_tcad/simplemos_sentaurus2022/source/simplemos_nominal_fps.cmd`。

## 冻结比较基线

### 工况矩阵

固定 `Lg=0.25 um`，八个器件覆盖三因素组合：

| 器件 | NWell/Boron (cm^-3) | GOxTime (min) | LDD Dose (cm^-2) |
|---|---:|---:|---:|
| n17 | 1e17 | 10 | 1e14 |
| n18 | 1e17 | 10 | 2e14 |
| n19 | 1e17 | 15 | 1e14 |
| n20 | 1e17 | 15 | 2e14 |
| n21 | 2e17 | 10 | 1e14 |
| n22 | 2e17 | 10 | 2e14 |
| n23 | 2e17 | 15 | 1e14 |
| n24 | 2e17 | 15 | 2e14 |

匹配 NWell 对为 n17/n21、n18/n22、n19/n23、n20/n24。每个器件计算 `Vd=0.05, 1.0 V` 两条曲线；`Vg=0.00...2.50 V`，步长 `0.05 V`，每条 51 点，总计 16 条曲线和 816 个直接偏压点。

### 默认 Sentaurus 物理与求解路径

- `EffectiveIntrinsicDensity(OldSlotboom)`；
- Silicon `Mobility(PhuMob HighFieldSaturation Enormal)`；
- Silicon `Recombination(SRH(DopingDependence))`；
- 先 Poisson，再 Poisson/Electron/Hole；先升漏压，再扫栅压；
- 四电极为 source、drain、gate、substrate；gate 接触附着在 Oxide 外边界；
- 默认 `Math` 为 `Extrapolate Iterations=20 ExitOnFailure`。

M46 已确认 M45 后的生产实现没有改变 M8 基线：16/16 条曲线通过，816 点完整，16 份逐点比较 CSV 与 M8 位级一致。默认基线最大差位于 n23、`Vd=0.05 V, Vg=0.05 V`：

| 量 | 数值 |
|---|---:|
| Sentaurus drain Id | `8.586250997030e-17 A/um` |
| Vela drain Id | `1.104643803536e-16 A/um` |
| 对数差 | `0.109418680924 dex` |
| 相对差 | Vela 高 `28.653%` |

## M45-M58 已完成工作

| 阶段 | 判别问题 | 结果/分类 | 对后续分析的约束 |
|---|---|---|---|
| M45 | SG 核统一和准费米打包修正后是否需重新定基 | 历史诊断重新定基，往返误差处于舍入尺度 | SG 核和准费米打包闭环 |
| M46 | 完整默认矩阵是否回归 | 16/16 曲线、816 点通过，与 M8 位级一致 | 生产基线不可随意更改 |
| M47 | n23 目标点是否是默认 BGN-on 自洽状态尖峰 | 势垒、浓度、SRH、准费米随相邻栅压平滑 | 不支持局部 BGN/SRH 状态异常 |
| M48 | 漏电流差由哪个端口/连续性项承担 | substrate 电子分配为目标漏差的 189.2%，source 反向补偿；SRH 仅 1.57% | 主问题转为端口分配语义 |
| M49 | substrate 节点场输运因子能否携带尖峰 | 99.994% 留在端口边界观测残差 | `boundary_observable_limited` |
| M50 | 原生 substrate 面通量是否等于默认端口 eCurrent | Tcl 面积分与独立节点边界积分一致，但与默认端口不一致 | `native_observable_mismatch` |
| M51 | `CurrentWeighting` 是否改变状态或端口值 | 六状态字段完全不变，目标 substrate 异常几乎消失 | `weighted_third_observable` |
| M52 | `DirectCurrent` 是否复现 Weighted | 两者舍入尺度一致，字段完全不变 | `direct_matches_weighted` |
| M53 | Direct 能否作为完整 Id-Vg 默认替代 | 最大误差增至 2.126075 dex，仅 3/16 曲线过原门槛 | `direct_increases_global_error` |
| M54 | Default→Direct 是共模还是反对称变化 | 目标点共模 93.38%，全矩阵中位仅 43.08% | 完整矩阵是混合端口变化 |
| M55 | NWell 是否改变默认接触关联 `DopingWell` | source/drain 阱面积配对变化最高 12.56%；SRH 总电流项抵消 | `nwell_changes_well_topology_surface_redistribution` |
| M56 | Si/SiO2 是否需要显式双节点设置 | 共享几何节点和边，界面势连续，载流子只在 Si | 标准界面拓扑闭环 |
| M57 | 导出节点场能否重构完整阱边界通量 | 分区无歧义、内部界面抵消，但 Direct 接触锚点仅 144/288 通过 | `boundary_quadrature_limited` |
| M58 | 冻结掺杂补偿、零轮廓和网格拓扑能否解释 NWell 趋势 | 误差 8/8 增长，但补偿仅 3/4、拓扑/近顶点变化仅 2/4且限低 LDD | `no_systematic_frozen_tdr_relation` |

## 当前最重要的证据

### 1. 自洽状态不能解释 n23 的 Vg 局部尖峰

n23 在 `Vg=0.00, 0.05, 0.10 V` 的 Id 对数误差分别为 `0.017239, 0.109419, 0.039123 dex`，但势垒代理为平滑单调的 `0.020488, 0.021690, 0.022607 dex`。同一 GOxTime/LDD 的低背景控制 n19 三点误差为 `0.023190, 0.025249, 0.026272 dex`。

目标点的势垒代理只解释观测对数差约 19.8%，积分 SRH 幅值变化为 `-0.009905 dex`，SRH 空间 Pearson 相关为 `0.999817614`。所以不能用 BGN、势垒或 SRH 的自洽场尖峰解释目标误差。

### 2. 默认 substrate 端口电子电流是第三种观测，不等于原生面通量

目标点 substrate 电子电流：

| 观测 | 数值 (A/um) |
|---|---:|
| 默认端口 `eCurrent` | `+4.670984402844e-17` |
| Tcl `eCurrentDensity·ContactSurfaceNormal` 面积分 | `-2.376797223030e-19` |
| `CurrentWeighting` | `-2.376113673567e-19` |
| `DirectCurrent` | `-2.376113673558e-19` |

Weighted 与 Direct 在六状态内最多相差 `9.260e-31 A/um`。两者相对 Tcl 面积分仍有约 `2.876e-4` 的小型离散定义差，但这不是原始 0.109419 dex 漏端误差本身。

重要反例是 M53：把 Direct 扩展到完整矩阵后，最大跨 TCAD 误差由 `0.109419` 增至 `2.126075 dex`，13/16 曲线未过原门槛。因此不能通过全局切换 Direct 修正默认 Id-Vg。

### 3. NWell 会改变默认接触关联阱，但外部节点场不能逐项复刻其内部积分

M55 中，n19→n23 的 source/drain 关联阱面积由 `5.824327e-10` 降至 `5.108676e-10 cm^2`，缩小 `12.287%`；substrate 阱面积增加 `1.629%`。全部匹配对的 source、drain、substrate 阱支持都发生离散变化。

阱内 SRH 电荷项对电子和空穴符号相反，在总电流中严格抵消；默认/Direct 总差由阱表面相对物理接触面的重分配承担。M57 能重构 `DopingWells` 节点标签、阱边界和内部界面抵消，但 source/drain 物理接触面节点场线积分无法普遍通过 Direct 锚点。因此：

- “NWell 改变默认关联阱拓扑”有证据；
- “某个外部可计算阱面通量项精确等于默认端口差”尚无证据；
- 不能把 M57 未闭合解释成新的自洽物理差异。

### 4. n23 的补偿/网格特征真实，但不是完整矩阵统一机制

M58 从八个冻结工艺 TDR 读取 `NetActive, BActive, AsActive, PActive`。全部节点精确满足：

```text
NetActive = AsActive + PActive - BActive
```

公共坐标上，高/低 NWell 的 BActive 中位变化为 `0.30103 dex`，施主总量 P95 变化约 `0.0011-0.0012 dex`，说明主要干预确实是背景硼翻倍。

n19/n23 低 LDD 对：

| 指标 | n19 | n23 |
|---|---:|---:|
| 全局节点数 | 1480 | 1482 |
| Silicon 节点数 | 942 | 942 |
| 栅边最小 `abs(NetActive)` (cm^-3) | `1.772434e15` | `1.596681e15` |
| 最大补偿指标 | 68.337 | 99.906 |
| 浮置 n 型同号分量 | 0 | 2 |
| 单节点浮置分量 | 0 | 2 |
| 零轮廓到顶点最短距离 (um) | `4.552107e-5` | `2.233663e-5` |

n23 两个最临界节点位于 `x=0.005045938 um, y=±0.125 um`：

```text
AsActive = 7.965630e16 cm^-3
PActive  = 9.012476e14 cm^-3
BActive  = 7.896086e16 cm^-3
NetActive= +1.596681e15 cm^-3
compensation metric = 99.906
```

然而 n20/n24 高 LDD 对是决定性反例：两个漏压下误差都增加，但最大补偿度从 71.65 降至 26.03，没有新增浮置同号分量、拓扑变化或更近的零轮廓。两个低 LDD 配对的四个漏压行全部出现拓扑/轮廓变化，两个高 LDD 配对的四行全部没有，表明存在明显的 `NWell × LDD` 交互，而不是单一 NWell 主效应。

## 已闭环主题：后续不要重复

除非出现新的独立反证，外部分析不应建议重新执行以下工作：

- 全局 HFS/RefDens/速度饱和参数扫描；
- SG 电流核或通用 Vela 接触电流提取排查；
- 准费米参考/增量冻结打包排查；
- 用关闭 BGN、关闭 SRH 或重新拟合迁移率作为默认修正；
- 用 `DirectCurrent` 全局替换默认 drain Id；
- 把 SiO2 当导电区，或为当前标准 Si/SiO2 静电界面强行添加显式双节点/热发射设置；
- 用 M57 的导出节点电流密度线积分冒充 Sentaurus 内部默认 `DopingWell` 面积分；
- 仅凭 M58 描述性相关系数宣布求解器缺陷或直接修改生产网格。

## 尚未解决的问题

### A. 八组误差增长缺少统一解释

事实是 8/8 个匹配 NWell/漏压行误差增长；但冻结结构指标没有同样的 8/8 共变。可能是多个机制叠加，也可能是当前指标没有测到真正相关的离散量。尤其需要解释为什么高 LDD 的 n20/n24 在没有补偿/拓扑增强时仍然恶化。

### B. Sentaurus 默认端口算法内部语义仍不完全可观测

默认 `.plt eCurrent`、Direct/Weighted、原生接触面通量是不同观测。现有导出量只能证明默认算法与接触关联阱变化有关，尚不能获得其内部每个阱面、单元或权重贡献。需要手册级语义、可导出的原生离散贡献，或一个能严格区分候选算法的最小 deck。

### C. 结构变化与物理变化耦合

提高背景硼同时改变：净掺杂场、补偿位置、工艺扩散结果、网格点位、`DopingWell` 标签及自洽载流子状态。只比较原始 n19/n23 无法区分“字段变化”和“网格/标签变化”。下一实验必须保持一个因素不变，不能再做多因素同时变更。

### D. 当前误差指标是每条曲线最大值

M58 用 M46 每条曲线最大 Id 对数误差与静态 TDR 指标配对。最大值发生的 Vg 可因器件改变，这种汇总可能掩盖偏压局部机制。下一分析应同时保留固定 `Vg`、最大值位置和完整 51 点误差形状，避免只相关八个最大值。

### E. 极弱电流下的相对/对数指标可能放大端口语义差

目标电流约 `1e-16 A/um`。必须同时报告绝对电流差、对数差、四端 KCL、电子/空穴分量和偏压邻点连续性。单独优化 dex 误差可能把极小的绝对端口重分配误判为主要器件物理误差。

## 建议其他模型重点评估的候选方向

这些是待评估方向，不是已授权的生产修改。

1. **固定网格的背景硼字段替换/插值控制**：在一套共同 Silicon 网格上构造低/高背景硼场，保持 As/P、几何、接触和网格完全一致，再分别回放 Sentaurus 与 Vela。它能区分掺杂物理主效应和网格/标签副作用。必须先定义电荷守恒、插值误差和 TDR 合法性门槛。
2. **固定掺杂的局部网格双运行**：只在 n19/n23 临界补偿轮廓附近进行成对局部细化，同时以 n20/n24 作为高 LDD 负控制；物理、偏压路径和端口算法保持不变。若只改善低 LDD，支持补偿-网格交互；若四组一致改善，说明现有 M58 指标可能过窄。
3. **完整 51 点的分层响应建模**：按 `NWell × LDD × GOxTime × Vd × Vg` 对绝对误差和 dex 误差分别建模，检查 NWell 效应是否集中在特定 Vg/电流区间，而不是继续相关每条曲线最大值。该方向可先只读执行。
4. **最小 Sentaurus DopingWell/端口语义 deck**：设计一维或极小二维多阱器件，使物理接触面通量、内部阱界面通量、SRH 体积分和默认/Direct/Weighted 端口量都可解析闭合，用来辨识默认算法。不要在完整 SimpleMOS 上继续猜测不可见权重。
5. **离散算子共同支持比较**：若 Sentaurus 能导出原生 edge/box 电流或单元面贡献，再在完全相同几何支持上比较 Vela SG 通量；没有共同支持数据前，不应重新打开通用 SG 实现。

优先级建议：先做方向 3 的只读完整偏压分层；它成本最低，并能决定方向 1/2 应瞄准哪些 Vg。若目标仍强烈局限于低 LDD 和极弱电流，再冻结一个严格的方向 2 合同。方向 4 与前两者并行价值高，但取决于 Sentaurus 手册或原生输出能力。

## 给其他大模型的直接分析请求

可将以下文本连同本报告交给其他模型：

> 请把仓库中的 M46-M58 机器报告视为事实来源，不要重新建议 HFS、SG、BGN、SRH、准费米打包、通用接触提取或全局 DirectCurrent 调参。请重点完成四项工作：第一，解释为什么背景 Boron 从 1e17 增至 2e17 cm^-3 时八个 NWell/漏压配对的最大 Id-Vg 误差都增加，但 M58 的补偿/零轮廓/连通拓扑只在低 LDD 系统变化；第二，根据 Sentaurus Device 手册阐明默认端口 eCurrent、CurrentWeighting、DirectCurrent、ContactSurfaceNormal 积分和 DopingWell 之间可能的离散语义关系，并明确哪些是文档事实、哪些是推断；第三，参考有源码的 TCAD 实现，判断背景掺杂、近完全补偿节点和网格符号翻转可能怎样影响控制体、接触盆地或电流分配；第四，给出最多三个单变量判别实验，每个实验必须列出固定量、干变量、控制器件、预期正反结果和停止门槛。请同时解释 n20/n24 这个高 LDD 反例，不要只拟合 n19/n23。

要求外部模型在引用 Sentaurus 手册、论文或开源代码时给出可核验章节、文件路径或链接，并把“官方文档事实”“源代码事实”“从本项目证据作出的推断”分栏陈述。

## 证据入口

先读：

- `docs/validation/simplemos_m45_m52_comparison_summary_2026-09-02.md`
- `docs/validation/simplemos_m58_compensation_contour_mesh_audit_2026-09-02.md`
- `reference_tcad/simplemos_sentaurus2022/compensation_contour_mesh_audit/m58_compensation_contour_mesh_audit_report.json`
- `reference_tcad/simplemos_sentaurus2022/simplemos_m58_compensation_contour_mesh_audit_evidence.json`

核心机器账本：

- M46：`reference_tcad/simplemos_sentaurus2022/full_matrix_requalification/`
- M47：`reference_tcad/simplemos_sentaurus2022/default_bgn_state_attribution/`
- M48：`reference_tcad/simplemos_sentaurus2022/terminal_partition_continuity_closure/`
- M50：`reference_tcad/simplemos_sentaurus2022/native_substrate_face_flux_export/`
- M51/M52：`reference_tcad/simplemos_sentaurus2022/{current_weighting_attribution,direct_current_attribution}/`
- M53/M54：`reference_tcad/simplemos_sentaurus2022/{direct_current_full_matrix,terminal_common_mode_attribution}/`
- M55：`reference_tcad/simplemos_sentaurus2022/doping_well_terminal_attribution/`
- M56：`reference_tcad/simplemos_sentaurus2022/si_oxide_interface_topology/`
- M57：`reference_tcad/simplemos_sentaurus2022/doping_well_boundary_flux_closure/`
- M58：`reference_tcad/simplemos_sentaurus2022/compensation_contour_mesh_audit/`

## 本地复核命令

在工作树根目录执行：

```powershell
D:\msys64\ucrt64\bin\python.exe scripts\run_simplemos_m58_compensation_contour_mesh_audit.py --verify
D:\msys64\ucrt64\bin\python.exe tests\regression\test_simplemos_m58_compensation_contour_mesh_audit.py
ctest --test-dir build-release --output-on-failure
```

M58 复核结果为 `M58 verified: no_systematic_frozen_tdr_relation`；专用回归 7/7 通过。提交前的完整 CTest 为 772/772 通过。

## 当前结论的使用边界

- M58 是冻结输入结构归因，不是因果干预；
- 描述性相关不能证明网格或求解器缺陷；
- `NWell` 在此算例中的物理含义不能外推到其他 MOS 工艺 deck；
- 默认端口观测差异不能等同于默认求解状态错误；
- 在没有单变量控制和共同离散支持前，不应修改生产物理、网格或端口默认值。
