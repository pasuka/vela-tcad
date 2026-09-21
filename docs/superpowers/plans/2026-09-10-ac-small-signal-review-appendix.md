# Vela AC 小信号审核、源码与资料附录 v4

日期：2026-09-10；代码基线：`d76434089dcd8564069be4b822260d7aa9c72662`。
配套 [实施计划](2026-08-31-ac-small-signal-simulation-development-plan.md)、
[规范合同](../specs/2026-08-31-ac-small-signal-contract.md)、
[差异账本](2026-09-10-ac-known-difference-ledger.json)。

保留原审核条号便于检索。§1.1/1.2 是历史处置，不覆盖 v4 的 §1.3 或规范合同。
下文报告是证据而非任务队列；本次没有运行仿真，拟新增功能仍需阶段放行。

### 1.1 外部审核意见处置记录

| 历史审核编号 | 处置 | 当时决定（现以 v4 合同为准） |
| --- | --- | --- |
| B1 | 采纳，改变 P4 核心架构 | 端电流和线性化改为未约束连续性残差/J/M 接触行提取；SG kernel 降为 oracle |
| B2 | 采纳，改变 P1 核心架构 | 电极电荷改为未约束 Poisson 残差/J 接触行反力；不独立重推生产 `D dot n` |
| B3 | 采纳并前移到 P0 | ACCompute 目标偏压成为强制 DC 落点；P5 只比较精确共同目标点 |
| B4 | 采纳，设为 P3 前置硬任务 | `LinearSolver` 必须新增 `numericFactorizationCount()`；AC 使用专用 solver 和不可变块矩阵 |
| B5 | 采纳 | 冻结 QF reference/scales/QP；FD 只局部扰动 BC，不经过 Newton 重求解；QP 默认 fail-closed |
| N1 | 采纳 | 储存项实现归属 `CoupledDDAssembler`，独立头文件只承载接口/结果类型 |
| N2 | 采纳 | 增加 essential、thermionic/Schottky、insulating pin、sheet charge 等边界行分类表 |
| N3 | 采纳 | P3 前新增 P2.5 纯介电 Poisson AC |
| N4 | 采纳 | stacked/interleaved 块布局由 fill ratio 决策，不在方案阶段冻结 |
| N5 | 采纳 | 单位激励只是导数归一化；不暴露为具有非线性含义的物理旋钮 |
| N6 | 采纳 | P3 增加 `1e3--1e12 Hz` 重掺杂 MOS 调理扫描 |
| N7 | 采纳 | 测试矩阵增加 avalanche-on 近击穿二极管 AC |
| N8 | 采纳 | 使用只读 `onAcceptedPoint` 观察者，避免侵入式重构 `DCSweep.cpp` |

这张表记录的是方案修订决定，不是功能已经实现的声明。

### 1.2 v3 复核修订记录（优先于 v2 表中的简写）

| 编号 | 修订 | 实施影响 |
| --- | --- | --- |
| R1 | Essential 反力与 Robin 边界通量分开提取 | WP-R 增加 boundary stamp 分解，未资格化的 Robin AC 拒绝执行 |
| R2 | 明确接触储存与 Poisson 电荷的局部抵消恒等式 | P2/P4 增加非零接触体积的动态守恒测试 |
| R3 | 区分电极电荷导数 `C_Q` 与低频总端口 `C_AC` | 禁止对导通器件普遍要求二者相等或 `C_Q` 列和为零 |
| R4 | 区分物理 Jacobian 与 Newton 稳定化/近似 Jacobian | P0.5 增加导数资格审计，冻结坐标不等于冻结物理响应 |
| R5 | 数值分解以矩阵缓存失效为计数边界 | 新矩阵一次；同矩阵后续 RHS 零次；M=0 跨频可继续复用 |
| R6 | P0 仅冻结落点合同，P4 才验收 DC 调度实现 | 消除 P0 必须依赖尚未实现的观察者/裁剪逻辑的循环 |
| R7 | 补充可计算误差指标、FD 条件和低频展开 | 分离严格结构门、暂定跨工具门和性能目标 |
| R8 | 对齐当前 reference_tcad 和单位版本约定 | 不再新增已退出的 examples 工程目录；区分 depth_m/depth_um |
| R9 | 当前已有 `runDCSweepStepControl` 的规则目标裁剪 | P4 扩展目标清单并复用步进器，不从零写第二套延续循环 |

## 1.3 2026-09-10 审核处置（v4）

| 编号 | 处置与限定 | 落点 |
| --- | --- | --- |
| A1 | 采纳；几何总面积本身有效，问题是材料电荷/储存支持。DC-only Poisson 修正不直接视为 AC 合格 | 合同4.3.4，WP-R，AC-DIFF-001 |
| A2 | 采纳 DC 前置门与 A/C 分门；不接受“DC 误差必然是导数误差下界”的泛化。SimpleMOS 已有 M79--M81 资料，不写成仍停在 M12--M18 | P5，AC-DIFF-002 |
| A3 | 采纳；P0 从 TDR/区域/掺杂/方程确认 poly，金属替代另列合同，不猜耗尽百分比 | P0/P5，AC-DIFF-003 |
| A4 | 采纳有效 Math/导数合同；不能仅凭 Newton 开关推定 AC/IFM 默认。M81 报告也提示此风险，版本手册证据仍须归档 | P0/P5，AC-DIFF-004 |
| A5 | 采纳；P3 三选一，不改全局 Types.h；bytes/time/原复方程误差共同决定 | 合同4.5，WP3b/P3 |
| B1 | 采纳；独立 AC 工厂关闭构造期 floor，bitwise residual replay，pattern builds<=2。makeArclengthAssembler 当前不会自动关闭 floor | 合同4.6，WP-R |
| B2 | 采纳；scaled RHS=+1/V0，physical 同时处理缩放前后两处 carrier pin | 合同4.5/4.6，P3 |
| B3 | 采纳；明确四个 step-control 调用点和 arclength 的不同范围，验收实际分派 | P4 路由表 |
| B4 | 采纳；纯介电沿用 3N solved gauge，M=0；另需电势参考和连通性 | P2.5 |
| B5 | 采纳；使用正交 sweep.ac，保留 iv/bv_reverse，名称来自实际 mesh | §5.4，P4 |
| B6 | 采纳基线卫生；分开记录代码基线、未提交 v3 输入与文档版本；本轮文档独立提交，不篡改代码 hash | 计划文首，Git 文档提交 |
| C | 拆为合同/计划/附录，并附预登记账本；统一 P0.5=WP-R 与 AC-L 映射；修正接触节点求和下标 | 三份 Markdown 与账本 |
| D/E | 新增非对称材料体积、快照复现、模式计数、poly/DC/路由测试，依赖前移 | WP-R/P3/P5、§8 |


### 2.3 Vela 当前代码基线

| 能力 | 当前证据 | 结论 |
| --- | --- | --- |
| DC 偏压扫描与延续 | `include/vela/simulation/DCSweep.h`、`src/simulation/DCSweep.cpp` | 可复用 |
| 耦合 DD 残差和 Jacobian | `CoupledDDAssembler::residual/assembleJacobian` | 复用装配路径；并非所有配置的 Newton J 都是完整物理导数 |
| 未知量布局 | `psi, phin, phip`，共 `3*N` | AC 储存导数必须尊重现有缩放和参考坐标 |
| 实数稀疏直接求解 | `LinearSolver` 支持模式和相同矩阵数值分解复用 | 可复用；P3 与局部 complex SparseLU 比较 |
| 多端接触和稳态电流 | `ContactCurrent` | 可复用，但需增加线性化 |
| 区域体电荷 C--V | `TerminalCharge` 和 `cv_quasistatic` | 仅为 legacy/prototype，不能直接认领 AC |
| 单变量 DC 外接负载 | `CoupledLoadLine` | 只能参考 bordered-system 设计，不是一般 MNA |
| Fermi/OldSlotboom/迁移率/SRH | 已有相关模型 | 可复用；须审计 AC 所需导数 |
| 局部复数后端 | Eigen SparseLU 支持复数，当前未有独立 AC wrapper | 新建局部 `ComplexLinearSolver` 可行，无须改 `Types.h` |
| 动态储存/质量矩阵 | 核心中无对应算子 | 缺失 |
| 频扫、Y 矩阵、ACExtract | 无执行路径 | 缺失 |
| SDE 高斯植入 | 前端 fail-closed，内部 mesh builder 仅常数掺杂 | 非 AC 核心；短期用 TDR 绕过 |

`ContactCurrent::computeFromResidual` 已经证明接触连续性残差行求和可产生端电流；
`NewtonSolver::makeArclengthContactCurrentFunctional` 也已经通过
`J.transpose()*residualWeights` 构造过同一函数的状态导数。AC 应把这个先例升级为
受支持的统一接触反力接口，而不是为 mobility/Fermi/BGN/high-field/avalanche 再写
第二份 SG 导数。现有方法传入空边界条件，只是 ideal/ohmic 回归锚，不能自动推广为
Robin AC 实现。P1 前必须提供物理行视图及可分离的自然边界 stamp（见 4.3）。

2026-09-10 当前主线源码复核锚（行号随提交变化，以符号为准）：

| 事实 | 当前符号/位置 | 对方案的约束 |
| --- | --- | --- |
| DC 残差行端电流 | `ContactCurrent::computeFromResidual` | B1 有仓库内先例 |
| 残差函数状态导数 | `NewtonSolver::makeArclengthContactCurrentFunctional` | 可用 `J^T*w`，无需第二套 SG 导数 |
| per-node QF reference | `CoupledDDAssembler::setQuasiFermiReferenceFields` | AC FD 必须冻结 reference field |
| contact-basin repartition | `NewtonSolver` warm-start/repartition 路径 | FD 不得经过正常 Newton solve |
| thermionic natural flux | `CoupledDDBoundaryConditions::thermionic` 和 assembler residual/J | physical view 不能传空 bcs |
| numeric factor cache | `LinearSolver::valuesMatchFactorization` | 同矩阵值必须 bitwise 不变 |
| 可见对象计数器 | 只有 `patternAnalysisCount()` | P3 必须新增 numeric count getter |
| 已有全局性能计数 | `LinearSolver.cpp` 的 `linear.factorize_calls/cache_hits` | 并非完全不可测；新增对象计数用于隔离 AC/DC 和精确验收 |
| Newton 数值对角修正 | `CoupledDDAssembler::assembleJacobian` 的 `carrierDiagonalFloor_` | 不能当作物理 dR/dx 消费 |
| 不完整源项导数路径 | 同文件 legacy avalanche 注释明确省略 driving-field/mobility 导数 | AC 必须按配置资格化，不能仅凭“已有 Jacobian”认领 |
| 当前 Schottky 支持边界 | `docs/config_schema.md` contacts；residual 的 `bcs.thermionic` | Newton thermionic_robin 原型存在，但不代表自然端口 AC 已验证 |
| 当前基准目录 | `tests/regression/README.md` | 新器件基准进入 `reference_tcad/`，小型单测进入 tests |
| 多材料节点体积 | `BoxGeometryBuilder.cpp` 累加所有相邻 Tri3；DD Poisson 使用 `vol[i]` | 共享 Si/SiO2 节点可能把氧化层面积计入半导体电荷，守恒自洽不等于物理正确 |
| 单一 Jacobian 模式缓存 | `fixedJacobianBoundarySignature_`；pattern/add 都跳过约束行 | 独立 AC 装配器、批量提取 physical J，再生成 solved J；每快照模式构建至多 2 次 |
| 独立装配器先例 | `NewtonSolver::makeArclengthAssembler()` | 可复用工厂配置路径，但当前传入 DC floor，AC 工厂必须明确禁用并恢复终态 reference |
| 当前 DC 目标步进 | `DCSweepStepControl.h` 与 `DCSweep.cpp::runDCSweepStepControl` | 已有规则 nominal target/stop 裁剪、retry 和 recorder；缺的是 AC 专用目标/快照合同 |

优先读取 `docs/README.md`、`docs/architecture.md`、`docs/config_schema.md` 和当前源码；
架构文档中更粗粒度/旧的 Schottky 状态不覆盖当前配置文档的受限支持说明。旧故障
记忆只作为风险线索，本次没有重现“3 个预存失败”等历史数量，不能写成当前测试状态。

现有 `cv_quasistatic` 在相邻 DC 点计算：

```text
C = (Q[k] - Q[k-1]) / (V[k] - V[k-1])
```

其中 Q 是用户指定区域中的 `q*(p-n+Nd-Na)` 体积分，还可按接触距离裁剪。它既不是
严格的金属电极高斯通量电荷，也不保证所有端电荷闭合。开发 AC 时不得悄悄改变
现有输出语义；应新增明确的电极电荷方法，并保留 legacy 模式兼容性测试。

## 2.5 独立工作树证据：只读复核，不等同主线功能

当前主线代码基线保持 d764340。以下是本机可读的其他工作树报告，本文未运行它们的
历史命令、未切换 checkout、未合并其代码。引用这些报告的方法/证据范围，不直接复制
其阈值或“通过”状态到 AC。P0/各阶段应将需要的中性证据固定到可移植 manifest；
其他模型无法读取本机路径时，先标记 evidence_unavailable 而不是假定通过。

### LDMOS 材料体积

- 工作树 `.worktrees/templates-ldmos-phase-a`，读取时 HEAD
  `d4285815f648be3f0b20719cbbe4d33e6a1b38b8`。
- [2026-09-02 G3 material-local qualification](../../../.worktrees/templates-ldmos-phase-a/docs/validation/templates_ldmos_g3_material_local_charge_volume_qualification_2026-09-02.md)，
  文件 SHA-256 `C1431FEC8BB455618C9D8EE78FEE0B4F8EA82C8863630B5F64F96C1AA939831E`。
- 报告中固定两点 gm/Sentaurus 比由 0.756924 变为 0.984547；这是该 G3/profile 的
  报告证据，不能称为主线 AC 实测。候选仅改变 Poisson 移动/掺杂体积，保持 continuity
  sources、固定/界面电荷、储存与 SG couples 不变，并拒绝未资格化 mixed-Voronoi 组合。
  因而不能只移植这个 DC 补丁便宣称合同4.3的动态抵消成立。
- 当前主线 BoxGeometryBuilder 最近相关提交 f5a3a6d，仍累加所有相邻 Tri3 面积；
  该事实已在本轮源码核对，不等于各 AC 配置误差已定量归因。

### SimpleMOS 方法与后续进展

- 工作树 `.worktrees/simplemos-sdevice-validation`，读取时 HEAD
  `c4e2a4f1228fd33a6abfdb30030666c90689ae7d`。
- [参考 README](../../../.worktrees/simplemos-sdevice-validation/reference_tcad/simplemos_sentaurus2022/README.md)，
  SHA-256 `8EDCA21542789A10266F0C0F6FD327212AE00F50D2FD90CCDA830D0738AF7E5C`：
  复用不可变 TDR/输入哈希、独立 equilibrium/drain/gate 状态链、精确网格禁止插值、
  受控物理阶梯与冻结比较器。其 SimpleMOS 原始物理/阈值不等于 AC 示例。
- [M79--M81 进展报告](../../../.worktrees/simplemos-sdevice-validation/docs/validation/simplemos_m79_m81_execution_progress_2026-09-05.md)
  描述 M79c 数值导数校准、M80b 原始状态与重建密度分开核对、配对误差账本；
  其中“M81 待上传”已被下一份执行报告取代，不采用该过期状态。
- [M81 原生 IFM 与固定电荷校准执行报告](../../../.worktrees/simplemos-sdevice-validation/docs/validation/simplemos_m81_native_ifm_execution_2026-09-05.md)，
  SHA-256 `27B960375321D6B5ACBF1CAE5183DB6D7A7D6C4CEE26810B15BD6F717B107935`：
  两点原生响应方向的独立校准已完成，但报告仍声明高 NWell 绝对 Id--Vg 差异未解决。
  不能把“观察器/IFM 校准通过”改写为“器件跨工具已对齐”。
- 该报告引用版本手册对 IFM 完整导数默认的说明；本轮未重新取得该手册页，故 AC P0
  仍必须归档自己对应版本/模式的 effective derivative 证据，而不是从 Derivatives
  是否出现推断省略了哪些项。也不把该零频 IFM 的支持推广为本 MVP 的 0 Hz C 输出。

以上 SHA 标识实际读取的文件内容，不推断其他工作树整体干净或其文件与 HEAD 完全一致。


## 12. 供其他大模型重点审核的问题

请同时读取规范合同和阶段计划，区分主线事实、独立工作树报告、推导与尚未执行的门。
审核结论使用 接受/有条件接受/阻塞，逐项提供章节或 gate_id、证据、反例、建议与测试。
没有源码/VM 访问时明确限制，不把拟新增接口或报告中的历史命令当成现有功能/执行授权。

1. 材料支持策略是否使 Poisson 移动电荷、S/M 和连续性源项在所有资格材料上闭合？
   是否误把 DC-only 的 Poisson 体积补丁当作完整 AC 支持？
2. 快照复现是否真正捕获 DC 终态 packed/reference/config/bcs，关闭 floor 是否仅影响 J？
   bitwise 门的同构建范围和旧 checkpoint 重闭合条件是否足够明确？
3. physical J 完整模式、两处 insulating pin 和每 snapshot<=2 次模式构建是否可测？
4. scaled BC 的 1/V0、电子/空穴行权重、Poisson 单位及接触储存抵消有无反例？
5. poly gate 的 DD 区域/边界资格是否完整，金属替代会在哪些评分项破坏原始等价？
6. Sentaurus 有效 AC/IFM 导数默认是否有版本证据？是否把 Newton 开关误当 AC 语义？
7. 三路后端比较是否公平衡量原复方程误差、factor bytes、ordering 和多 RHS 计数？
8. DC 前置门与 A/C 独立阈值能否隔离材料/迁移率差异？是否仍有选择性遗漏失败点？
9. 四个 DCSweep 调用点及独立 arclength 路径的支持边界是否清楚，实际分派有无测试？
10. 低频 C_AC 展开、C_Q Gauss 微分、Robin 自然通量与完整端口 gauge/KCL 是否成立？
11. 独立工作树来源和数字是否被错误宣称为主线已验证？证据过期或缺文件时如何阻塞？
12. stage_acceptance 的依赖失效、partial 输出、模型拒绝和 AC-L 映射是否存在漏洞？

推荐回复格式：

```text
总体结论：接受 / 有条件接受 / 阻塞
阻塞项：[gate_id/章节] 证据、问题、建议修订
非阻塞改进：[章节] 建议和收益
新增测试：fixture、预期关系、误差定义及容差依据
依赖调整：原顺序 -> 新顺序，原因
证据限制：未访问的源码/VM/报告与仅推导的结论
```

## 15. 参考资料

公开资料复核日期为 2026-09-09；网页/main 分支会变化，P0/P5 必须固定版本/commit/hash。
4.3 和 4.8 是针对 Vela 离散合同的推导，不冒充下述文献直接给出的公式或已运行结果。

- [Sentaurus Device Training：MixedMode and Small-Signal AC](https://ghzphy.github.io/Sentaurus_Training/sd/sd_3.html)：
  公开教学镜像仅作语义导航；T-2022.03-SP2 的 VM 原始 deck、对应版本手册和实际
  输出才是 P0 版本合同来源。本次未重新获取 VM 运行证据。
- [DEVSIM 官方 solver documentation](https://devsim.net/solver.html)：DC 工作点上的 AC
  与电路激励；用于流程比较，不证明与 Vela 的状态变量/残差单位相同。
- [DEVSIM small-signal diode 源码](https://github.com/devsim/devsim/blob/main/examples/diode/ssac_diode.py)：
  `solve(type="ac")`、单位 AC 源和 `V1.I` 虚部提取；对照时必须反转源/器件电流方向。
- [Genius-TCAD-Open 官方仓库](https://github.com/cogenda/Genius-TCAD-Open)：README 声明
  small-signal AC 与 mixed device/circuit 支持；作为架构参考，不把 README 当作模型
  数值等价证据。若借鉴实现，固定具体源码版本并审核 GPL 许可边界，不直接复制。
- S. E. Laux, “Techniques for Small-Signal Analysis of Semiconductor Devices,”
  IEEE TCAD, 1985, DOI `10.1109/TCAD.1985.1270145`：
  [IBM 作者机构摘要](https://research.ibm.com/publications/techniques-for-small-signal-analysis-of-semiconductor-devices--1)。
  比较瞬态 Fourier、增量电荷分区和正弦稳态；本次查阅摘要，不声称重验全文推导。
- C.-K. Lin 等，“Frequency Domain Analysis of the Distribution Function by Small
  Signal Solution of the Boltzmann and Poisson Equations,” SISPAD 1999：
  [会议论文原文](https://in4.iue.tuwien.ac.at/pdfs/sispad1999/00799254.pdf)。
  支持频域线性化方法背景；其 Boltzmann 模型范围不等同于本计划 DD MVP。

当前仓库复核入口：

- [CoupledDDAssembler.cpp](../../../src/equation/CoupledDDAssembler.cpp)：`residualImpl`、
  Poisson/连续性源项、thermionic stamps、`assembleJacobian` 和 carrier diagonal floor。
- [ContactCurrent.cpp](../../../src/post/ContactCurrent.cpp)：`computeFromResidual` 的单位/符号回归锚。
- [UnitScaling.cpp](../../../src/core/UnitScaling.cpp)：`continuitySourceIntegralFactor` 与输入单位版本。
- [LinearSolver.cpp](../../../src/solver/LinearSolver.cpp) 和
  [LinearSolver.h](../../../include/vela/solver/LinearSolver.h)：数值缓存、全局 profiler 和缺失的对象计数。
- [DCSweep.cpp](../../../src/simulation/DCSweep.cpp) 和
  [DCSweepStepControl.h](../../../include/vela/simulation/DCSweepStepControl.h)：已有规则目标裁剪与失败/成功 recorder。
- [配置合同](../../config_schema.md)、[回归规范](../../../tests/regression/README.md)、
  [PN2D 专项资格边界](../../validation/pn2d_bv_validation.md)。

- [Eigen SparseLU 官方模板文档](https://libeigen.gitlab.io/eigen/docs-nightly/classEigen_1_1SparseLU.html)：
  原生 real/complex 支持、compressed column-major 输入、ordering 和分解接口；
  2026-09-10 查阅。P3 固定本机实际 Eigen 版本，文档不替代三路性能实测。
