# SimpleMOS 高 NWell：文献、开源实现与 Sentaurus 手册驱动的后续 debug 计划

日期：2026-09-05。工作树：`D:/code-repo/vela-tcad/.worktrees/simplemos-sdevice-validation`。

本文是经过资料核对的执行建议，**不是已冻结的 M79 合同，也不是新实验结果**。本轮完成联网检索、本地手册阅读和源码核查；没有执行新器件求解、修改生产参数或改动 M78 冻结证据。后续各阶段应在执行前分别冻结合同。

建议主线：**复用现有完整 DD 端口伴随，先用已知 M74 干预验证电流响应预测，再解释真正的 Vela–Sentaurus 绝对误差；必要时用 Sentaurus 零频 IFM 获取独立响应证据，最后只运行一个有物理依据的候选 A/B。**

## 1. 当前证据带来的约束

- `NWell` 是背景 Boron 由 `1e17` 增到 `2e17 cm^-3`，不是传统 n-well 掩膜注入。
- M60 已闭环深关断端口 burst；当前目标为 `Vg=0.55–0.90 V` 的平滑误差。深关断亚 fA 点不参与后续相对误差主判据。
- M65/M66 matched-ni/no-BGN 是诊断参考支路，M8/M46 BGN-on 是生产支路。M60 原物理参考电流不能直接代入 no-BGN 误差。
- M74 的配对差缩小来自低 NWell 电流增加更多，16 个诊断端点绝对误差均恶化；M77 组合也未证明精度改善。
- M78 的 99.96% 是**电子体积干预的配对电势响应**在 gate 阶段的份额；91.57% 是 **M73 电子预测与 M74 实际电势响应的差**中由缺少自洽反馈解释的份额。二者都不是原始 Id–Vg 差异的闭合率。

需要收紧 M78 的外推边界：界面节点只承担约 12% 的硅区 `|delta psi|`，足以否定其合同中“80% 电势响应局限于混合节点”的预设，但**不能排除界面源项对 Id 具有很高的影响权重**。一个局部电荷扰动可以产生分布式电势响应；电势范数、源项大小和端口电流敏感性是三个不同量。[R1–R3]

M78 标出的 Si/Nitride 支撑也值得保留，但总 forcing 中的大份额不能自动解释高/低 NWell 差分。后续必须分开记录静态共同项、gate 增量、配对差分，以及正负贡献的抵消；材料三重点不得重复计数。

## 2. 已核对的资料与结论

### 2.1 Sentaurus Device User Guide

本地原件：[sdevice_ug.pdf](D:/software/ATCNIS01-202203/sdevice_ug.pdf)，封面为 **T-2022.03, March 2022**，1739 页。本算例执行版本为 T-2022.03-SP2；本文不把同主版手册当成补丁版私有实现的证明。1178、817 页已渲染核对，提取文本与关键公式一致。

| 页码 / 章节 | 手册明确说明的内容 | 对本算例的影响 |
|---|---|---|
| 61–62，Specifying Doping Species | 已存在的 net/total doping 数据集可优先于重新求和；active 与 chemical 掺杂有不同读取规则 | 只读核对 Poisson 实际用的净离化掺杂与进口 donor/acceptor；不能只比较名义 NWell 数值 |
| 229，式 37 | Poisson 由移动载流子与离化掺杂等电荷共同决定 | 电子-only 体积 A/B 是诊断干预，不等于完整物理离散修正 |
| 304，式 151–152 | 本征浓度由带边 DOS、带隙与温度决定；BGN 另进入有效本征浓度 | 复用 M65 的一致性结果，只审计约定；不重开 Nc/Nv 或 ni 拟合 |
| 1177–1178，式 1265–1266、表 187–188 | 实际非线性方程按单元贡献装配，介电常数等参数可逐单元处理 | 应审计各材料的源项体积与介电通量，再累加到共享电势节点 |
| 1183–1185 | AverageBox/NaturalBox/MixAverageBox 的定义及钝角、非 Delaunay 情形不同 | barycentric、signed circumcentric、截断体积不能混称“Sentaurus 体积” |
| 1187–1188 | `MeasureCoefficients.debug` 支持单元控制体积与边系数的读写；`GrdNumbering` 指定网格编号 | 复用 M34 直接几何证据。必须分别验证 Measure/Coefficients 的局部槽位排列 |
| 781–784 | IFM 可计算固定电压下电流响应；`ACCoupled` 可使用单频 0 Hz；零频下不能提取电容，普通结构的电压 Green 函数可能不适定 | 可以尝试零频**电流**响应观察器，不把氧化层隔离栅的电压 Green 函数或零频电容当参考 |
| 808–810 | `DeterministicVariation(DopingVariation(...))` 支持空间定义；默认也包含 doping 对 mobility/BGN 的作用；可同时使用 `-Mobility -BandgapNarrowing` 去掉这两条影响 | 必须声明“总掺杂响应”还是“电荷通路响应”，不能把混合参数响应当纯 Poisson 响应 |
| 814，式 803；817，表 142 | 端口响应是 Green 函数与分布扰动的积分；提供 `CurPotReACGreenFunction`、`CurECReACGreenFunction` 等 | 可与 Vela 完整耦合伴随做独立对照；需先校准单位、符号、体积权重和观察节点 |

手册中的这些 IFM 名称已核实，但尚未在当前 SP2 deck 上执行。因此其解析支持、模型限制及加载工作点身份是未来 pilot 的验收项。零频 IFM 是小信号线性分析，不能直接代表有限幅度的体积替换。

### 2.2 文献

| 来源 | 与当前问题直接相关的内容 | 使用边界 |
|---|---|---|
| [R1：Ghione & Filicori, 1993](https://amsacta.unibo.it/id/eprint/2087/)，IEEE TCAD 12(3), 425–438，DOI 10.1109/43.215004 | 把分布扰动通过伴随/Green 函数映射到器件端口响应 | 原文多数载流子器件例子提供方法依据，不证明本 SimpleMOS 物理模型等价 |
| [R2：Bonani et al., 1998](https://iris.polito.it/handle/11583/1397422)，IEEE TED 45(1), 261–269，DOI 10.1109/16.658840 | 多维双载流子器件的 Green 函数与线性扰动方法；也是手册第 23 章引用的依据 | 使用其线性响应思想，不引入噪声模型作为待调参物理 |
| [R3：Giles & Pierce, 2000，§2](https://people.maths.ox.ac.uk/~gilesm/files/ftc00.pdf)，An Introduction to the Adjoint Approach to Design | 离散伴随的转置系统、目标函数导数和参数直接项；正问题/伴随点积恒等式 | 属于通用数值方法，不能替代器件干预验证 |
| [R4：TU Wien，Box Discretization §7.2](https://www.iue.tuwien.ac.at/phd/triebl/node30.html) | box 法的守恒来源，以及右端源项积分的不同离散方式 | 有限体积守恒不意味着任意两种控制体积给出相同电流；不据此重写 SG |

本轮已阅读 R1 的机构存档摘要、R2 的机构页面及可检索摘要、R3 的离散伴随推导、R4 的离散章节。没有把检索摘要写成未实际完成的全文审阅。

### 2.3 开源代码

DEVSIM 已联网核对官方仓库，同时阅读本地 checkout `58a9a87083db00c6cadc0b4011c801db2cec5844`：

- [simple_physics.py](https://github.com/devsim/devsim/blob/58a9a87083db00c6cadc0b4011c801db2cec5844/python_packages/simple_physics.py)：`CreatePE` 按区域定义电荷；`CreateOxidePotentialOnly` 单独建立绝缘区电势方程；`CreateSiliconOxideInterface` 施加连续电势。
- [TriangleNodeVolume.cc](https://github.com/devsim/devsim/blob/main/src/GeomModels/TriangleNodeVolume.cc)：二维单元节点份额由 `0.25 * ElementEdgeCouple * EdgeLength` 形成，不能假定是三角形面积的 1/3。
- [本地 Equation.cc](D:/code-repo/devsim/src/Equation/Equation.cc:831)：node source 与该 region 的 `NodeVolume` 相乘，导数也使用相同体积。其含义与[官方模型文档](https://devsim.net/models.html)一致。

Genius 已联网核对 Cogenda 官方开源仓库，同时阅读本地 checkout `543da8452d5dfd33e6f8c457f962f6f670f0fce7`：

- [ddm1_semiconductor.cc](https://github.com/cogenda/Genius-TCAD-Open/blob/master/src/solver/ddm1/ddm1_semiconductor.cc)：区域内电荷项为 `e*(doping+p-n)*volume`；Jacobian 使用同一表达式的自动微分。
- [本地 ddm1_boundary_is_interface.cc](D:/code-repo/Genius-TCAD-Open/src/solver/ddm1/ddm1_boundary_is_interface.cc:78)：绝缘区 Poisson 行向半导体侧累加，另施加电势相等约束。

这两个实现共同支持审计“区域贡献—共享界面条件—一致 Jacobian”这一结构。它们各自的数据结构和边界实现不等于 Sentaurus 的私有实现；**不能从逻辑双节点推导本算例必须添加物理双节点**。M34 已提供当前 SimpleMOS 共享 Si/SiO2 电势节点的直接证据。

本轮不建议新建整套第三求解器 MOS 对比。若需要独立验证，可在后面的几何分支采用已有开源代码构造少量二维局部装配算例。

## 3. 与已有工作相比，真正新增什么

仓库已有以下能力，后续应复用：

- [NewtonSolver.cpp](D:/code-repo/vela-tcad/.worktrees/simplemos-sdevice-validation/src/solver/NewtonSolver.cpp:3258)：`makeArclengthContactCurrentFunctional` 从未替换边界行的连续性残差构造端口电流与导数。
- 同文件的 `evaluateTerminalCurrentAdjoint`：已求解 `J^T lambda = dI/dx`。
- [vela_example_runner.cpp](D:/code-repo/vela-tcad/.worktrees/simplemos-sdevice-validation/src/tools/vela_example_runner.cpp:3827)：已有 `terminal_current_adjoint_probe`。
- `newton_carrier_term_probe` 与 [M15 脚本](D:/code-repo/vela-tcad/.worktrees/simplemos-sdevice-validation/scripts/run_simplemos_m15_operator_adjoint.py)：已有求解方程分项残差输出。

M14/M15 使用早期生产支路和不同器件/偏压矩阵，早于 M43/M44 的回放一致性修正；不能直接复用旧数值权重。新增内容是把**现有算法重新资格化到当前 matched-ni/no-BGN 状态**，并用 M74 的真实单轴自洽响应作为答案已知的校准，不再重演通用 SG、HFS 或接触提取排查。

## 4. 分阶段执行计划

### P0：冻结比较口径与输入身份

新求解：0。

输入保留三套独立身份：生产 M8/M46–M60；no-BGN matched-ni M65/M66；诊断干预 M74/M78。

主工况为 n19/n23、n20/n24，两种 Vd，共 8 个器件/漏压工况；诊断 Vg 沿用 M74：n19/n23 为 0.90/0.80 V，n20/n24 为 0.85/0.80 V，分别对应 Vd=0.05/1 V。通过 pilot 后加入 n17/n21、n18/n22，补足 16 工况、8 配对。

分别冻结：state/curve/config/materials/mesh/doping、Sentaurus exports、M34 几何证据、当前 probe 源码和二进制哈希。对 48+48 阶段状态沿用 M78 门禁，不自动重算失败输入。

只读核对 NetDoping、ionized/active/chemical 字段、ni、温度、有效栅接触条件和单位。已有审计通过的项目仅做身份核验；发现新字段不等价才进入相应解释，不能将它变成参数扫描。

主要指标固定为每工况 `e=log10(|Id_Vela|/|Id_Sentaurus|)`、高/低各自的绝对误差、配对 `e_high-e_low`。不再单独用配对闭合率评判修复。

### M79：完整 DD 端口响应观察器资格化

目的：验证观察器能够预测**已知 M74 电子-only 干预**的实际 Id 变化，而不是只解释 delta psi。

初期读取 8 pilot 工况的 M65/M74 gate 状态，16 次伴随回放；扩展后总计 32 次 gate 伴随回放。必要时增加 drain 状态，最多另 32 次。equilibrium 保留静电身份检查，不对近零电流计算 log-current 伴随。无新 Sentaurus 或 Vela 自洽求解。

设完整耦合状态 `u=(psi, phin, phip)`，实际求解方程为 `F(u,a)=0`，端口目标为 `g(u,a)=Id`。使用：

```text
J = dF/du
J^T lambda = dg/du
dg/da = (partial g/partial a) - lambda^T(partial F/partial a)
delta log10|Id| ~= delta Id / (Id * ln(10))
```

这里 J 必须包括 Poisson–电子–空穴的所有耦合块、边界条件及生产材料/源项导数，不能继续用 M73 的单独介电矩阵逆代替。电流目标与实际基线 CSV 的提取路径必须数值一致；若不一致，停止该观察器，不另选“更方便”的端口算法。

M74 只改电子 Poisson 体积。直接残差源项在固定状态下计算为 `q*n*(V_Si_bary-V_legacy)`，正负号以当前 residual 定义为准。源项应分别在基线与候选状态取值，生成前向/反向端点预测；目标直接项也应实际检查，不能凭假设遗漏。

先用方向有限差分核对 `J*v` 和 `dg/du*v`，再对照已存的 M74 实际 ΔId。必须输出每工况电流预测误差、每配对误差以及 n19/n23 和 n20/n24 在 Vd=1 V 的反例表现。

建议在合同中预冻结的门槛：

- 完整状态、电流身份和单位通过；电流回放相对误差不超过 `1e-8`。
- 伴随相对残差不超过 `1e-10`；有信号方向的有限差分相对误差不超过 `1e-3`，并报告步长减半结果。
- 正问题/伴随点积相对闭合误差不超过 `1e-8`。
- 有限 M74 干预的 log-current 预测误差不超过 `max(0.1*|实际增量|, 1e-4 dex)`。该项是本计划的判别阈值，不是文献保证。

若线性预测因有限幅度不合格：不能宣布“模型根因已找到”，也不能按相关性放行。优先对比两端线性化；仍不合格时单独冻结小幅度/中点校准合同，最多先对一组低/高工况做新增自洽验证，记录新增求解数。未经该步骤，不把失败线性模型用于跨求解器定量归因。

### M80：当前绝对 Id 误差的分项与空间账本

目的：从“解释 M74 响应”转向“解释 M65 剩余绝对误差”。只在 M79 资格化后执行。

读取已有 M65/M69 的 Sentaurus 节点状态与对应 Vela 状态。先确认映射后的 Boltzmann 关系、接触条件、节点身份；绝不把导出的 n/p 和重新计算的 n/p 混成同一个完整状态。

对每个工况分开记账：

```text
实际误差 = g_V(u_V) - g_S(u_S)
         = [g_V(u_V) - g_V(u_S)] + [g_V(u_S) - g_S(u_S)]

第一括号的局部线性预测 ~= -lambda^T F_V(u_S)
第二括号 = 同一映射状态下的端口目标差，单独保留
```

第二个等号是代数分解；伴随预测只是局部近似。`F_V(u_S)` 是 **Vela 对映射 Sentaurus 状态的反应**，不是 Sentaurus 内部残差，也不自动等于某个 Vela 模型错误。必须与实际总误差并列报告余项。

上述局部预测在映射状态 `u_S` 处组装 `J_V(u_S)` 与电流导数，并求对应伴随；不能不加说明地代入另一个工作点的伴随。基态残差不为零时保留 `F_V(u_V)` 修正，同时比较两端线性化以识别有限状态差的非线性余项。

将残差按 Poisson 介电、电子电荷、空穴电荷、掺杂电荷、电子连续性、空穴连续性及必要的边界项拆开。Poisson 的原始分项可能很大而互相抵消，应先减去 Vela 基态对应分项，并同时保存有符号和绝对值贡献。

空间按材料单元贡献与物理区域两套维度统计：Si/SiO2、Si/Nitride、其它界面和体区；浅源、沟道、浅漏、深体。共享节点上的不同材料项保留所属单元标签，三重点不重复计数。报告 `-lambda_i*delta F_i` 对电流的贡献；不再以节点数、总电荷量或 `|delta psi|` 份额做根因排序。

必须分别回答：

1. 哪一项解释高 NWell 的绝对误差？
2. 哪一项解释高减低的增量？
3. 哪些大项仅为共同偏移或抵消项？
4. n20/n24 高 LDD 控制是否保持结论？

候选方向至少应在 6/8 配对保持符号一致，并在分区边界移动及两个 Vd 下稳定；同时有独立手册/代码/原生数据依据。仅有漂亮的伴随分项和不小的剩余项时，结论仍是“观察器定位”，不能进入生产修改。

### M81：补一个独立 Sentaurus 响应证据

优先方案是 **0 Hz 电流 IFM pilot**。先做 n19/n23、Vd=0.05 V、既定 Vg 两个工作点，成功后最多扩到 8 个 pilot 工况。

- 从现有工作点保存文件加载。若工具要求 DC 再闭合，单独计数并验证 Id/状态重闭合误差；不得称作完全没有 DC 求解。
- 使用 `ACCoupled` 的 `ObservationNode` 观察 drain 电路节点，单频 0 Hz，输出电流 Green 函数。具体 mixed-mode 电路与端口偏置在执行合同中审阅，不能改变器件的 DC 边界条件。
- 输出候选包括 `CurPotReACGreenFunction`、`CurECReACGreenFunction`、`CurHCReACGreenFunction`。先用一个归一化、空间可解析的小扰动校准符号与单位，再比较同一扰动下的 δId，而不是直接逐点比较两个求解器的原始 Green 数值。
- 可考虑 `DeterministicVariation(DopingVariation(...))`。若目标是去掉 mobility/BGN 的掺杂影响，成对指定 `-Mobility -BandgapNarrowing`；这不等价于全局关闭对应物理。仍需检查当前 SRH 等掺杂依赖是否进入实际响应，不能未经验证把它称作纯 Poisson RHS 注入。
- 若目的确实是纯 RHS 响应，优先从 `CurPot...` 与已明确单位的 Poisson 源项做积分；体积是否已包含必须以校准实验判定，禁止重复乘体积。
- 零频率不用于提取 Cgg；如后续需要栅电容，应另开小正频率/准静态合同。

建议资格门槛：pilot 的 DC 电流重闭合不超过 `1e-5 dex`；相同小扰动的 IFM 响应与显式有限差分在两档幅度下相差不超过 5%，且符号一致。幅度应由工作点信号与数值噪声决定并事先冻结。有限差分需要新增 DC 求解，只有 IFM pilot 成功后才批准该最小矩阵，不能隐藏在“只读回放”计数中。

若 IFM 在当前 SP2 构建、模型组合或工作点加载上不可用，明确记录原因，转为 Vela 内部经过有限差分资格化的伴随证据；不杜撰 Sentaurus 输出，也不重开全局收敛调参。

几何分支只在 M80 指向相关源项时执行：先复用 M34 的 n23 `MeasureCoefficients.debug`；仅对缺少数据的控制器件补最小导出，按实际启动记录 Poisson 初始化次数。审计 region-local 系数、signed/native 节点份额、边界截断、Si/Nitride 及三重点。**M34 的 n23 界面几何已经闭合，不重做其推导，也不重复 M35 的深关断 factorial。** M74/M77 的 barycentric 负结果不能直接否定尚未等价的原生离散候选。

手册的 `BoxMeasureFromFile`/`BoxCoefficientsFromFile` 具有双向行为：指定文件不存在时计算并输出，文件已存在时则读取。补导出必须使用隔离的新目录并记录文件的初始不存在状态，避免把遗留文件意外变成求解输入。

### M82：一个因果明确的候选 A/B

仅当 M79/M80 给出稳定方向且有 M81 或原生几何证据支撑时，冻结一个候选。

候选可能是具体的材料源项/系数映射修正或已证实的输入语义修正；不是新的电子/空穴/掺杂组合扫描。若数学一致性要求多个装配项同时改动，应作为有明确定义的离散方案，并用局部装配测试证明，不能称作电子单轴。

先运行 8 个控制工况，再扩展到 16；所有其它物理、网格、偏压路径、端口函数和收敛设置保持原值。每次新三阶段工作流均记录真实求解数。

建议筛选条件：高 NWell 端点绝对误差中位至少减少 20%，高/低任一控制点不新增超过 `0.001 dex` 的退化，至少 6/8 配对增长改善。20% 与 0.001 dex 是待冻结的工程筛选值，不是“达到物理正确”的充分条件。

出现 M74 型“低 NWell 变差导致配对靠近”立即判为不通过。若候选的线性预测与自洽结果不符，返回观察器有效性/非线性误差检查，不再添加第二个补偿参数。

### M83：完整曲线及生产支路回归

仅对通过候选执行 16×51 点曲线。no-BGN 与原 BGN-on 支路分别使用同物理参考；没有完整同物理参考的支路不能标记为已验收。

报告逐点绝对/相对误差、平滑窗口最大/P95误差、每个高/低 NWell 工况、配对增长、强反型端点、数值收敛和 KCL。深关断使用 M60 资格化参考与单独的绝对电流误差门槛。

改进必须同时满足绝对精度与配对差要求，不能只报告一个全局中位数。必须覆盖强反型、两个氧化时间、两个 LDD 剂量；既有 M78 和生产结果不覆盖、不替换。核心装配测试、相关回归及完整候选矩阵通过后，才讨论默认值或生产基线更新。

## 5. 优先级和停止条件

| 次序 | 工作 | 初始新增器件求解 | 继续条件 |
|---|---|---:|---|
| 1 | P0 输入与口径门禁 | 0 | 所有冻结身份通过 |
| 2 | M79 现有完整 DD 伴随在 M74 上校准 | 0；仅线性回放 | 能预测真实 ΔId，含高 Vd 反例 |
| 3 | M80 当前跨求解器绝对误差账本 | 0 | 有稳定、可独立检验的候选方向 |
| 4 | M81 原生 IFM/缺失几何资料 pilot | 最小 2 工作点；DC 重闭合、FD另计 | 单位、线性幅度、参考身份通过 |
| 5 | M82 单个候选 | 8 控制工况起步 | 绝对误差与配对差同时改善 |
| 6 | M83 完整曲线及生产回归 | 16×51 点/合格支路 | 全矩阵、核心测试、保留证据通过 |

优先交付物是 **M79 冻结合同与电流响应资格报告**，而不是新默认参数。若 M79 无法用现有 M74 状态验证电流预测，先停止根因外推。

本计划不重启：HFS/速度饱和扫描、SG 内核排查、通用接触电流提取、准费米打包、BGN/Nc/Nv 拟合、全局电荷体积组合，以及仅凭手册添加 Si/SiO2 物理双节点。

## 6. 来源与复现记录

- 手册正文与两页渲染仅留在 ignored `build-release/m79_research/`，不把整本商业手册复制到版本控制目录。
- 文献链接、手册 SHA-256、本地开源源码版本/哈希及本仓库复用入口记录在 `reference_tcad/simplemos_sentaurus2022/debug_research/post_m78_sources_2026-09-05.json`。
- 当前计划没有新求解结论，也没有冻结未来验收门槛；P0/M79 执行前必须把以上建议转为独立合同。
