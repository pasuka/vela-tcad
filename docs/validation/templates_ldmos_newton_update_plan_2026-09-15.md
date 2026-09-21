# LDMOS 单点 Newton 更新映射与阻尼改进方案

日期：2026-09-15。代码基线：工作树 `codex/templates-ldmos-phase-a` @ `41abef4`
（已提交外推 A–E）；性能对照仍使用独立冻结的 R7 程序与输入，实施前记录两者
源码与程序哈希。本文记录用户已确认的设计选择及拟议阶段合同，**已获准实施；F0 完成，F1 和阻尼/刷新对照未晋级**；
本轮结果见[执行记录](templates_ldmos_newton_update_execution_2026-09-15.md)。
本文保留阶段合同与设计依据，冻结 R7 仍为生产默认。

前置研究与结果：[外推算法研究](vela_extrapolation_solver_research_2026-09-15.md)、
[外推 A–E 执行记录](templates_ldmos_extrapolation_execution_2026-09-15.md)、
[R7 重复性与 Newton 定位](templates_ldmos_r7_stability_newton_2026-09-14.md)、
[准费米表示小更新丢失诊断](templates_ldmos_stage4_newton_work_gap_2026-09-10.md)。

## 授权状态

| 事项 | 状态 |
|---|---|
| 优先级：先研究单点 Newton 的密度更新映射与阻尼，暂停扩展外推方案 | 用户确认 |
| F1 设计默认：投影下限 ε=0.01；密度映射覆盖全部硅节点自由载流子未知量 | 用户确认 |
| F1 同时运行 V1/V2 两个启用范围变体 | 用户确认 |
| F2 停止规则 R-A/R-B 在 F0 日志之后决定；R-B 属接受规则变化，届时需另行明确同意 | 用户确认该口径，未预先批准 R-B |
| F0 及后续阶段的实施 | 用户已明确授权，先 F0、F1，再比较 NLEQ_ERR 与按需 Jacobian 刷新 |
| 本机 Sentaurus VM 文件传输、构建与仿真 | 用户已明确授权 |
| 生产默认变更 | 未授权；需完成 62 点联合验收与同 VM 配对 |

最新授权覆盖构建、插桩重跑与本机 VM 实验；不自动变更生产默认或批准 F2 R-B。V1/V2“同跑”指纳入同一对照矩阵，计时仍串行执行。
实施字段与冻结参数见下文及执行记录；剩余阶段不得跳过其晋级条件。

## 1. 证据基础

以下数字来自冻结 R7 Linux Release UMFPACK 双栅压曲线的逐次迭代记录
（`reference_staging/templates_ldmos_joint_20260914/linux_environment/copied_vm/results/r7_full_vg{4,8}`），
只读复算，未重跑求解器。

### 1.1 相位分解

分类按 floor → stall → damped → linear → quadratic 的优先顺序，仅使用已记录的
`scaled_l2_before/after`、`alpha`、`initial_alpha`、`limiter_component`、
`max_abs_direction_by_block`：

| 相位（定义） | Vg=4 V | Vg=8 V |
|---|---:|---:|
| floor：迭代起始 merit < 1e-9 但门限未全部通过 | 23 | 24 |
| stall：初始全局 α < 1e-2 | 22 | 13 |
| damped：未归入 floor/stall，且接受 α < 1（初始 α ≥ 1e-2） | 90 | 78 |
| ……其中受全局限幅 / 由 QF 分量限幅 / 由 ψ 限幅 | 66 / 65 / 1 | 58 / 58 / 0 |
| linear：α=1 且 merit 比 > 0.1 | 107 | 85 |
| quadratic：α=1 且 merit 比 ≤ 0.1 | 124 | 121 |
| 合计更新 / 线搜索候选 | 366 / 496 | 321 / 446 |

边界说明：`quadratic`/`linear` 是描述性标签，不证明二次收敛，也不是新算法
更新数的严格下界；stall 全部由 QF 分量限幅，因此"受限幅 88/71 次"包含 stall。
linear 更新中最大 QF 方向落在 [0.25, 4]·$V_t$ 的比例为 103/107 与 79/85，
但该统计使用固定 300 K 的 $V_t$，且历史没有记录方向最大节点，需 F0 补节点级证据。

代表轨迹：Vg4 5.333 V 第 1–4 次 α=9.6e-5/2.0e-3/2.0e-2/9.0e-2，merit 比
0.9998/0.992/0.996/0.887；第 8–12 次 α=1 而 ψ、T 方向已 ≤1e-7，fn/fp 最大方向
0.030/0.046 → 0.024/0.026 → 0.019/0.024 → 0.011/0.020 → 0.0033/0.012 V。
Vg4 30.667 V 的五次 floor 更新 merit 约为 1.8e-12、各块方向约 1e-14；
现有逐次日志没有块/行门限快照，不能由终态通过反推此前各次阻塞的是哪项门限。

### 1.2 R8 密度模式的失败机理

生产 R8（`copied_r8_initguard` 快照）为：密度迭代预算 60、`density_update_maximum_bias_V=1.0`、
`density_update_requires_prediction=true`，且仅在 merit > 1e-6 时启用。
Windows 代表点对照 `production/r8_density_full_controls`（预算 60，无偏压/预测限制）：

| 点 | 基线更新 | R8-60 更新 | 结果 |
|---|---:|---:|---|
| Vg4 0.375 / 5.333 / 16 / 30.667 V | 14 / 15 / 15 / 10 | 5 / 7 / 48 / 60 | 通过 / 停滞拒绝 / 未收敛 / 未收敛 |
| Vg8 0.375 / 5.333 / 16 / 30.667 V | 13 / 15 / 12 / 7 | 5 / 8 / 11 / 14 | 通过 / 停滞拒绝 / 通过 / 通过 |

Vg4 16 V 轨迹：第 14 次起 α=5.60e-3 逐次降至第 43 次 1.35e-4，空穴最大方向由
5.1 V 增至 210 V，merit 仅由 0.555 降到 0.412，终态空穴 QF 与基线差 4.5 V。
机理在源码中确定：正性约束 `-0.99/relative` 在
[ElectrothermalSimulation.cpp](../../src/simulation/ElectrothermalSimulation.cpp)
的密度模式中取全局最小 α（约第 442–450 行），单个线性化目标 ≤0 的少子节点
把全部未知量的步长压到很小。若该正性限制控制 α 且没有后续回溯，对应节点
的线性密度目标为旧值的 0.01；实际限制节点及每次密度变化仍需 F0 节点日志确认。

### 1.3 其他基线数据

- 冻结 R7 deck：`initial_step_V=0.1`、`minimum_step_V=1e-4`、`maximum_step_V=4/3`、
  `max_newton=60`、**`growth_newton=12`**（代码默认 8）、`predictor=linear`；
  点输入 `defer_recentered_candidate_jacobian`、`skip_equilibrium_poisson_transport`、
  `reuse_ialmob_local_preparation`、`residual_ialmob_values_only` 均开启，
  `electrothermal_linear_solver=umfpack`，`diagnostic_stagnation_window=5`，
  电压修正上限默认 0.2 V，温度 30 K。
- R7 漏压阶段（Vg4 / Vg8）：装配 799 / 702 次（其中仅残差 393 / 341）、分解 366 / 321、
  装配 395.2 / 313.8 s、分解 139.5 / 111.2 s、回代 16.4 / 12.6 s。
  每次更新约 1.11/1.12 次全 Jacobian 装配加 1.07/1.06 次仅残差装配。
  内部端到端墙钟 626.5 / 516.0 s **包含初始化**，不能直接扣除仅漏压阶段的
  三项计时后归因输出开销。原生日志的装配计时也不直接提供与 Vela 对应的调用计数。
- 线搜索接受判据为任意下降（约第 472 行）；Bank–Rose 型比值
  $(1-\|F_{k+1}\|/\|F_k\|)/\alpha$ 在 R7 5.333 V 首步为 1.635（Vg4）与 12.0（Vg8），
  说明充分下降判据不能识别被限幅压小的步。
- A–E 的固定目标序列结果（140/139/137/134 等）不是同条件三方对照：切线仅在每档
  9.333 V 一点替换初值。可用结论是"已测外推改动无总成本收益"，不能推出初值无关。

行号以当前工作树版本为准，实施时以函数与语义定位。

## 2. 可证伪假设

| 编号 | 假设 | 证伪方式 |
|---|---|---|
| H1 | 部分 linear 更新来自非简并、近定 ψ/T 下的过量密度修正：`relative=Δn/n` 从大于 −1 的一侧接近 −1，QF 原始修正幅值接近局部 Vt，密度约除 e | F0 同节点逐步核对密度、ψ/T 贡献及统计偏离；80%覆盖率作为待冻结的诊断目标，未达到即弱化；`relative≤−1` 属非正线性目标，另行分类 |
| H2 | stall 与 QF 限幅 damped 来自单节点退化 QF 分量决定全局 α | F1 局部投影后两类计数应接近消失 |
| H3 | floor 更新由某一门限在可计算浮点下限附近摆动造成 | F0 逐次门限日志定位阻塞门限并与下限比较 |
| H4 | F1 节省的分解/装配成本足以覆盖投影、反算、回溯和恢复成本 | 与同输入基线比较完整成本；R7 无密度反算，不能要求新增反算次数≤0，须报告其新增成本并检验净收益 |

## 3. 阶段总览

| 阶段 | 固定条件 | 单项变化 | 主要判据 | 退出/停止 |
|---|---|---|---|---|
| F0 | 冻结 R7 原程序/输入保留；旧日志只读分析 | 独立插桩候选增加 opt-in 字段，获授权后重跑以采集旧日志缺失数据 | H1/H3 证据及插桩前后数值轨迹一致性 | 字段完整、来源可核验、轨迹一致后完成 |
| F1 | R7 物理、门限、线搜索、UMFPACK、`growth_newton=12` | 密度更新改逐节点投影；取消近收敛/偏压/预测前提 | 相位计数、失败、试算、反算、装配、分解；状态一致性 | 代表点出现基线无的失败且在范围内不可解 → 停止记录 |
| F2 | F0 日志 | 无（决策） | 阻塞门限归因；两种规则语义 | 用户确认后才实施 |
| F3 | F1 冻结候选 | 线搜索接受判据 | 试算、damped 计数、终态原门限 | 无总成本收益 → 不晋级 |
| F4 | F1（+F3）冻结候选 | 外层增长策略 | 含失败的完整推进总成本 | 同上 |
| F5 | 任意冻结候选 | 重装配/输出整理实现 | 装配次数与时间、有限精度等价 | 等价性不成立 → 回退 |

## 4. F0 诊断补全

opt-in 开关（建议 `diagnostic_iteration_trace: true`），默认关闭且不改任何数值路径，
在 `history[]` 每项追加：

| 字段（建议名） | 内容 | 用途 |
|---|---|---|
| `gate_blocks_weighted_l2[3]`、`gate_row_max_ratio`、`gate_row_violations`、`gate_qualified_rows` | 迭代起始的三块与逐行门限值；现有终止判据在 merit ≥ 1e-9 时短路，不求值 | H3、F2 |
| `max_direction_node[4]`、`max_direction_node_T_K[4]` | 各块最大方向节点及其温度 | 局部 $V_t$ |
| `max_direction_qf_step_Vt[2]`、`max_direction_relative_density[2]`、`max_direction_density_m3[2]` | 该节点原始 QF 步（局部 $V_t$，带符号）、线性化 $\Delta n/n$、当前密度 | H1 |
| `density_limiter{node,carrier,relative}` | 密度模式决定 `densityAlpha` 的节点；现有 `limiter_node` 只是 QF/ψ/T 限幅节点 | H2 |
| `projection{count,max_violation,per_carrier[2],charge_ratio}` | F1：按试算及最终接受候选分别记录投影节点与违反量；电子/空穴修正分别按 qA 加权，保存原始电荷量（C/m）。比值分母须使用同单位的未归一化 Poisson 残差及冻结的零分母处理，不能除以无量纲块门限 | H4，电荷一致性不预设可忽略 |
| `density_fallback_to_qf` | 8 次试算后回退 QF 方向 | 成本转移 |

门限字段同时记录上限、通过状态和候选接受后的值；节点日志保留带符号方向，
分开记录 ψ/T 与载流子对密度线性增量的贡献。逐次门限求值为 O(N)，新增时间
单独测量，不预设可忽略。F0 日志用于诊断，生产计时应关闭日志或采用匹配开关。

分析脚本（建议 `scripts/analyze_templates_ldmos_newton_phases.py`）：相位分类、
damped 拆分（限幅/未限幅，QF/ψ/T）、局部 $V_t$ 带、H1 节点级判定、floor 段阻塞门限
分类、按点/区间/尝试结果汇总，输出机器可读 `summary.json`；结论只从汇总读取。

G0：先完成旧日志分析；获授权后用独立插桩版本重放两档，报告 H1 覆盖率与
H3 阻塞门限分布，并核验非计时数值轨迹与基线一致。假设未获支持用于调整预期，
不把缺失字段填零，也不将只读旧日志统计认作新增节点证据。

## 5. F1 逐节点密度投影候选

### 5.1 设计默认（已确认）

- 点输入开关 `diagnostic_density_projection: "off"|"v1"|"v2"`，缺省 off 保持现行行为；
  本轮冻结投影比例 0.01，不另引入可调比例以缩小实验维度。
- 密度目标 $n_i^{new}=\max\big(n_i(1+\alpha\,\mathrm{rel}_i),\ \varepsilon n_i,\ n_{floor}\big)$，
  $n_{floor}$ 为极小绝对下限，保证 `electrothermalDensityQf` 的 `ratio>0` 与
  `inverseFermiDiracHalf` 定义域；其 SI 数值、与旧密度的关系及是否影响零步恒等性
  冻结为 `min(n_old,1e-250 m^-3)`，零步恒等，低于该下限的旧状态不被抬高；
  正密度除以有效态密度若仍下溢，则拒绝候选，不静默修正。非有限 `rel_i` 是无效数值，不能按普通负密度静默投影；应记录
  原因并拒绝密度候选/回退。局部投影只处理有限的越界线性目标。
- 全局 α 仅由 ψ 0.2 V 与 T 30 K 限幅决定；Dirichlet/接触约束行、有限空穴接触的
  自由空穴行、温度链式项、Fermi/BGN 精确反算沿用 R8 现有实现。
- 密度映射覆盖全部硅节点自由载流子未知量。
- 线搜索保持现行 merit 判据、8 次试算后回退 QF 方向、停滞监视 window 5，使 F1 与 F3 隔离。
- 两个启用范围变体同时运行：
  - **V1**：全部迭代、全部偏压、不要求预测。
  - **V2**：同 V1，但所有活动节点 $|\Delta\varphi_{raw}|<10^{-2}V_t$ 时退回原 QF 方向。
    此处 Δφ 指自由载流子的原始 QF 方向，Vt 使用对应节点温度。该阈值是候选
    策略，不是全物理范围内的误差保证：定 ψ/T 的 Boltzmann 极限下，0.01 Vt
    对应的密度归一化二阶差约 5e-5，而相对于 QF 更新本身的差约 0.5%。Fermi/BGN
    及 ψ/T 联动需另做实际映射误差检查；尾部反算舍入可能影响收敛，由 V1/V2 对照判断。

### 5.2 实施要点

- 投影只改载流子候选；ψ、T 候选仍为 $x+\alpha d$。电荷一致性由 `projection.charge_ratio`
  度量。若其持续 >1 且伴随 ψ 块残差回升，备选 F1b：投影后对 ψ 做一次 Poisson 块
  重校正，单独计时。
- 投影后 QF 增量超过 1e-3 触发参考重定心与重装配，按现有路径处理并计入装配成本。
- Catch2 用例：Boltzmann 极限单节点线性行 V1 一步命中 $n^*$；$n^*\le0$ 时落到
  $\varepsilon n$ 且其余未知量 α 不受影响；Fermi–Dirac 与 300/515 K 反算一致；
  约束行不被映射。Python 回归：开关缺省时输出与 R8 逐值相同。

### 5.3 四级验证（逐级 go/no-go）

| 级 | 输入与基线 | 对照组 | 通过条件 |
|---|---|---|---|
| G1a 代表点 | `r8_predictor_inputs` 的 8 个 linear_baseline 点（Vg4/8 × 0.375/5.333/16/30.667 V），基线 101 更新/182 试算；R8-60 见 1.2 | 基线、R8-60、F1-V1、F1-V2 | 8 点全过原门限；状态差 ≤1e-8 V / 1e-7 K；linear+stall+QF 限幅 damped 下降；无新增失败，回退/试算/反算逐项报告并验证净成本收益；投影计数随迭代趋零 |
| G1b 固定目标序列 | A 的 R7 已接受序列 0–9.333 V，基线 140/137 更新（`sweep.step_policy="fixed_targets"`） | 基线 vs F1 | 更新、装配、分解下降；14/13 个状态在容差内一致 |
| G1c 前 8 精确点自适应 | B 的自适应基线 152/137 更新、15/13 次尝试（含失败） | 基线 vs F1 | 含失败与恢复的总成本下降；接受状态一致 |
| G1d 完整 62 点 | R7 366/321 更新、38/36 尝试；同 VM 串行配对 | R7 vs F1 | 原电学/热学/2% 局部场/30 meV 全过；墙钟 ≤1.5 倍原生；报告更新、装配、分解、CPU |

墙钟只作门限。R7 三组重复墙钟 CV 为 26%/30%，单组配对不能认领加速；认领需
≥2 组交错重复且区间分离。主指标为计数与内部阶段计时。G1a 使用 Windows UCRT64
Release 同平台配对；G1c/G1d 使用冻结 Linux GCC11 Release UMFPACK 环境。

### 5.4 风险与对策

| 风险 | 对策 |
|---|---|
| 投影候选被 merit 拒绝，8 次后回退 QF，成本反升 | `density_fallback_to_qf` 进入 G1a 判据；普遍出现则转 F3 讨论判据，不放宽 |
| 退化少子节点每步固定除 100，空穴 QF 仍漂移 | 记录投影节点集合是否随迭代收缩；不收缩即 H2 失败，停止 |
| 尾部舍入（V1）拖长 floor 段 | V2 对照 |
| 极小密度反算越界 | $n_{floor}$ 与极端值用例 |
| 路径变化导致外层失败模式不同 | G1b 固定序列先隔离单点效应，再进入 G1c |

## 6. F2 floor 段：先观测，后决策

用 F0 日志把 23/24 次 floor 更新按阻塞门限分类，并与各门限的可计算下限比较
（ψ 块：残差约 4e-21 C/m 经 `poissonScale` 归一后逐行约 1.5e-9，L2 约 6.6e-9，
上限 5e-8；电子/空穴块与逐行比值同法估计）。

两种规则供用户选择，本方案不预设：

- **R-A**：停止条件不变（全部门限通过）；仅当方向已低于按块浮点分辨率
  （ψ：eps·max|ψ|；QF：增量表示的 ULP；T：eps·max T）时跳过分解、只重评门限。
  当前方向通常需分解后才可知，不能用上一方向直接认定本次也小；跳过分解需有
  状态/算子未变等可验证依据。门限仍未通过时只能记录停滞失败，不能作为合格点。
  分辨率还须包含零值保护与 ref/inc 重定心误差；节省有限，不改接受规则。
- **R-B**：方向低于分辨率且连续两次门限值不变时以"floor 停止"接受。
  **即使容差数值不变，这也是接受规则变化**，需用户届时明确同意，并在报告中标注该类点。

## 7. F3 线搜索判据（F1 之后独立）

- 对照组：Bank–Rose 型判据 $(1-\|F_{k+1}\|/\|F_k\|)/\alpha\ge\sigma$，σ 取
  0.1/0.25/0.5 三档作控制，不引用为固定值（Bank & Rose, *Global approximate
  Newton methods*, Numer. Math. 37, 1981, doi:10.1007/BF01398257；原文参数与近似
  Newton 方向误差相关）；Deuflhard 自然单调性检验（每次试算多一次回代，复用已分解
  $J_k$；线性求解约占总时间 3%；Deuflhard, *Newton Methods for Nonlinear Problems*,
  Springer 2004）。
- 判据：剩余 damped 计数、试算数、终态原门限；在 F1 真实候选路径上评估下降性。
- 已知限制：两类判据都不解决限幅压步，只能减少无效回溯与预测阻尼因子。

## 8. F4 外层步长（F1 稳定后）

- 基线为冻结 R7 的 `growth_newton=12`、`maximum_step_V=4/3`、增长 ×1.5 / >20 次减半
  （[ElectrothermalSweep.cpp](../../src/simulation/ElectrothermalSweep.cpp) 接受分支）。
- 单步变便宜后"≤12 即增长"几乎恒真。先在 0–4 V 复现 R8 记录的 1.333→2.472 V
  外推比 9.94 场景，比较：现策略、B 的 `actual_step`、"≤N1 增长 / >N2 减半"
  （N1、N2 由 F1 后的更新分布确定）。先固定目标序列，再自适应；不与 F1/F3 混合比较。

## 9. F5 实现开销（可独立进行）

- 参考系重装配：同一物理状态在 (ref, inc) 与重定心后 (ref′, inc′) 下装配，逐元素比较
  残差与 Jacobian（预期相对差约 1e-15）；成立则比较"`defer_recentered_candidate_jacobian`
  现行路径"与"候选全 Jacobian 装配、跳过重定心重装配"的装配时间；数值路径差异按
  G1a 容差验证。
- 输出整理：扣除初始化及漏压两阶段已记录的装配/分解/回代后，Vg4 内部
  626.5 s 中约 67.8 s（10.8%）尚未归因；其中不全是 JSON。分别计时 `output.json`
  写入、`predict` 读取上一输出、ledger/curves 保存与门限求值，再决定是否改为内存
  传递或紧凑格式。

## 10. 证据、冻结与报告口径

- 每阶段独立证据目录（沿用 `reference_staging/templates_ldmos_<阶段>_<日期>/`）、
  冻结程序 SHA256、源码归档与哈希清单、`summary.json`；失败尝试全部保留。
- 复用现有入口：代表点对照沿用 `production/run_r8_density_full_controls.py` 的
  同种子/同期望状态口径；固定序列沿用
  [run_templates_ldmos_predictor_study.py](../../scripts/run_templates_ldmos_predictor_study.py)；
  完整曲线与配对沿用 `linux_environment/benchmark_r8_initguard.py` 及 62 点验收/局部场评分脚本。
- 每阶段执行记录文档；结论只从机器可读汇总读取；相位统计随附完整分类规则与局部 $V_t$ 口径。
- 报告区分：环境失败与产品失败；同目标序列的算法对照与自适应完整路径对照；计数收益与计时收益。
- 本方案执行前不修改物理模型、容差数值、门限模式；获准实验仍沿用原门限，
  唯 F2 R-B 若日后另获明确同意，才可独立改变接受规则并重新验收。不从 Windows
  单点计数推断 Linux 曲线；G1d 通过是晋级必要条件，不自动授权替换 R7 默认。

## 11. 参考实现与文献

- Genius `src/solver/ddm1/ddm1.cc`：`potential_damping` 对 ψ 用
  $f=\ln(1+\Delta V/V_{ut})/(\Delta V/V_{ut})$；`check_positive_density` /
  `projection_positive_density_check` 对负密度候选逐节点取 0.01·|old|
  （[源码](https://raw.githubusercontent.com/cogenda/Genius-TCAD-Open/master/src/solver/ddm1/ddm1.cc)）。
- DEVSIM：未知量 Potential/Electrons/Holes，`variable_update ∈ {default, log_damp, positive}`，
  收敛用更新范数 `absolute_error`/`relative_error`，`maximum_divergence` 检测发散
  （[手册](https://devsim.net/CommandReference.html)、[示例](https://raw.githubusercontent.com/devsim/devsim/main/python_packages/simple_physics.py)）。
- Bank & Rose 的[出版信息](https://link.springer.com/article/10.1007/BF01398257)与
  [作者全文](https://ccom.ucsd.edu/~reb/m270b/newton.pdf)：充分下降判据见式（3.1）。
- Bank, Rose, Fichtner, *Numerical methods for semiconductor device simulation*,
  IEEE Trans. Electron Devices 30(9), 1983：变量选择与阻尼的经典讨论。
- 这些资料支持实验设计，不证明本算例的收益。
