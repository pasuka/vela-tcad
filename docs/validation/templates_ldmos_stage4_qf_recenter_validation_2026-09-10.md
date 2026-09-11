# LDMOS Stage 4：Newton 内部重定心、步长与预测器验证（2026-09-10）

## 范围与结论口径

本次执行顺序为：单次 Newton 内部按需重定心；保持其他条件验证最大步长 0.2 V；恢复最大步长 0.1 V 后验证 secant 初值预测。仅验证 Vg=8 V、Vd 从 0 V 开始的前 8 个精确参考点，未运行完整 0–40 V 曲线。

参考点为 `0, 1.33333333333333, 2.66666666666667, 4, 5.33333333333333, 6.66666666666667, 8, 9.33333333333333 V`。从原合格零压种子开始，排除栅压预偏置；统计包括失败尝试、重闭合及密度恢复。Newton 更新不包含迭代 0 的检查，也不将仅改变参考坐标的重定心记为更新；重定心引起的额外分解和线搜索开销均计入时间和性能计数器。

保持 D5 物理配置：300 K、Fermi–Dirac、OldSlotboom BGN、SRH/Auger、constant_field vector-QF HFS、contact_node_cell 场回退；没有引入表面迁移率、量子、雪崩或热模型。网格为 10,241 节点、19,782 三角形，其中 5,723 个 Si 节点；使用 templates_ldmos_external_averagebox（16,237 条记录）。输入几何单位 μm、载流子密度 m⁻³，电流比较单位 A/μm。

原门限全部保留：Poisson/electron/hole 块分别为 `5e-8/1e-11/3e-10`，carrier-row 比率与端口 KCL 比率均为 `1e-8`，Newton 最大更新预算 160。保留外层重闭合、密度恢复和步长回退，没有将不合格中间态作为曲线结果。

## 实现

新增默认关闭的 `solver.quasi_fermi_recenter_on_stall`。在经典耦合 Newton 中，Poisson 块已合格且某个尚未解决的自由载流子行出现非零校正 `delta` 被浮点加法吞掉（`x + delta == x`）时，重分配 QF 的参考值和增量。

- 使用 TwoSum 保存 reference + increment 的低位，保持电势与物理态在浮点误差范围内等价；并非粗略清零全部增量。
- 清空依赖坐标的 continuity 缓存，使 active-branch fingerprint 失效，重算残差和闭合诊断，再重算 Newton 方向。
- 每次求解最多允许 8 次重定心，不重置初始残差归一化或已接受更新预算，不放宽任何验收条件。
- 触发检查在本轮线搜索之后执行，可避免后续多轮长尾，但本轮已完成的线搜索仍产生开销。
- 记录 `newton.qf_recenters`、`newton.qf_recenter_discarded_line_search_trials`；开启 local-update diagnostics 时额外观测密度、残差及 Jacobian 差异。

涉及 [NewtonSolver.cpp](../../src/solver/NewtonSolver.cpp)、[CoupledDDAssembler.cpp](../../src/equation/CoupledDDAssembler.cpp)、对应头文件、[数值测试](../../tests/test_newton_solver.cpp)及[配置文档](../config_schema.md)。本次没有修改既有物理模板默认值，也没有将实验外层预测器加入生产配置。

## Release 与证据身份

使用 `windows-ucrt64-release`，GCC/UCRT64，`CMAKE_BUILD_TYPE=Release`、`-O3 -DNDEBUG`；gprof 与 LTO 关闭。本次用于算法对照的冻结可执行文件不是 Debug，也不是带 `-pg` 的计时构建。实际线性后端为 Eigen SparseLU/COLAMD。

- 本次冻结 runner SHA256：`be486a03d52184b2944e435567c42ecfa643f7fee36a27449aa2bb4f810485f2`。
- 三项热点优化后、内部重定心前的比较基线：`1da0409e4af23ba2bdf87d2eb23c9f59d132819f592b6020a7c7b9ed60615744`。
- manifest 记录基于提交 `d8ff26804a59735ff3f446846e446e0a62add914` 的源码快照和工作区差异；不能仅凭提交号代表带有本地改动的二进制。

原始证据位于未跟踪目录 `reference_staging/templates_ldmos_qf_recenter_validation_20260910/`。`binary/manifest.json`、各 curve 的 `plan.json`、`fixed/ledger.json`、子进程 profile 和 `audit_summary.json` 提供二进制、配置、输入及结果哈希。`audit_results.py` 复核原门限、点位、配置/状态哈希和计数器，结果为 `final_audit.json`。

## 冻结父状态探针

四个相同父状态探针对比原 direct + reclose 总更新数和内部重定心求解：

| Vd / V | 原更新数 | 重定心后 | 重定心次数 |
|---|---:|---:|---:|
| 4 | 9 | 6 | 1 |
| 6.43333333333333 | 53 | 8 | 1 |
| 7.56666666666667 | 37 | 8 | 1 |
| 9.33333333333333 | 8 | 6 | 1 |

四点均通过原块、逐行及 KCL 门限；端电流相对原结果差异最大约 `1.90e-13`。重定心瞬间密度最大相对变化约 `1.43e-14`，Jacobian 相对范数变化最大约 `8.08e-13`。残差无穷范数的绝对变化最大约 `1.65e-10`，因此不声称残差逐位相同，必须重新检查最终门限。探针打开了额外诊断，不用于报告生产计时收益。证据：`point_probes/results.json`、`point_audit.json`。

## 完整前 8 点对照

| 策略 | 精确点数 | Newton 更新 | 物理推进 | 子进程 | 回退 | child wall / s | child CPU / s |
|---|---:|---:|---:|---:|---:|---:|---:|
| 三项热点优化后的原基线 | 8 | 1284 | 107 | 214 | 0 | 1024.84 | 896.64 |
| 内部重定心，最大步长 0.1 V | 8 | 822 | 107 | 110 | 0 | 618.43 | 529.28 |
| 内部重定心 + 外层 secant，最大步长 0.1 V | 8 | 599 | 110 | 117 | 3 | 517.49 | 425.67 |

单独重定心相对原基线：更新减少 35.98%，child wall 减少 39.66%。叠加外层 secant 相对单独重定心：更新减少 27.13%，child wall 减少 16.32%；相对原基线分别减少 53.35% 和 49.51%。

单独重定心：非零参考点电流误差中位数 1.594856%，P95 1.624118%，低压差分电阻误差 1.032119%；原工程/最终可用门限全部通过。外层总时间 648.32 s。

叠加外层 secant：非零参考点电流误差中位数 1.594856%，P95 1.624118%，低压差分电阻误差 1.032119%；原工程/最终可用门限全部通过。外层总时间 572.59 s。

时间为子进程 wall/CPU 累计，不含外层 Python 的启动、CSV 初值外推和评分等工作；外层总时间另见 progress.json。旧基线是之前独立运行的证据，并非本次交替 AB 配对测量；本次两个完整候选也为顺序单次运行，期间有轻量状态读取及一次证据审计。时间降幅是实测观察值，没有重复测量的置信区间。更新数、物理推进次数和分解计数更直接支撑算法工作量结论。

单独重定心保持 107 次物理推进、0 回退；仍在 1.43333333333333 V 使用了一次密度恢复，说明本策略不能完全替代原恢复流程。原有 214 个子进程降至 110 个；107 次内部重定心丢弃了合计 126 个已执行的线搜索试探。数值分解从 1,612 次降至 1,148 次。

## 更大步长：负向性能实验

在重定心开启的条件下，仅将外层最大步长从 0.1 改为 0.2 V（QF cap 仍与实际步长匹配）。到 2.90166666666667 V 时累计 38 次成功推进、75 个子进程、693 次更新、11 次回退；child wall 为 646.71 s、CPU 为 531.91 s，只完成了 0、1.333333、2.666667 V 三个精确点。

该时间已超过 0.1 V 策略跑完全部 8 点所需的 618.43 s，因此通过 STOP 文件在物理尝试边界行政停止。结果为部分曲线负向证据，不能称为完整 8 点耗时，也不说明继续运行必然不能完成。未生成完整曲线评分，不推荐采用 0.2 V。证据：`curve_step02/budget_stop.json`、`STOP`、`fixed/ledger.json`。

## 预测器的对照设计与局限

内置 `none/secant` 配对使用相同物理配置、0.1 V 最大步长、0.0025 V 初始步长与同一内部恢复设置；除了 mode 与输出路径，配置一致。两者都在 0.01043125 V 的恢复过程抛出 `equilibriumCarrierState: failed to bracket charge-neutral potential`，最后成功电压为 0.005875 V。完成的 attempts 只有 24 次更新；日志各记录 138 条正序号 Newton 更新，包含抛异常的求解。异常导致最终 profile 缺失，不能把缺失计数记为零，也不能只用完成的 attempts 统计总成本。

进一步发现，[精确 bias_points 分支](../../src/simulation/DCSweep.cpp)仅在整个精确区间完成后调用 acceptPredictorHistory，内部小步只更新 localPreviousSolution。首个参考区间的 predictorPreviousSolution 因此为空；[预测器](../../include/vela/simulation/DCSweepPredictor.h)直接返回当前态。已完成 attempts 的 secant predictor hash 均等于 initial hash。这组失败结果没有实际检验到 secant 外推能力，而且其固定 0.1 V QF cap 与已合格外层流程的随步长 cap 不同，不能用来否定预测器。

另一次最初的 native 设置误将名义 step 写为 0.0025 V，产生过密目标网格，已停止并排除；保留在 `native_first8/setup_abort.json` 作为实验设置错误记录。正式配对 `native_first8_v2` 使用名义 step=4/3 V，不混入前一组数据。

为单独验证预测器，使用 `run_outer_secant.py` 在原已合格外层 0.1 V 流程中对相邻两个成功态的 psi/phin/phip 做 secant 外推，比例上限为 2；首次推进及发生回退后的重试不外推，恢复步骤保持原逻辑。CSV 的 QF reference/increment 与预测后的物理 QF 保持一致，每次外推记录前后状态及预测初值哈希。预测初值可以不满足方程，最终接受态仍须通过全部原门限。

外层 secant 完整通过 8 点验收，共 109 次预测初值、599 次更新、3 次物理回退，数值分解 794 次。它减少了大部分正常步的更新数，但仍引入少数长尾/恢复，收益需按整条限定曲线累计。结果支持在当前限定场景继续使用内部重定心加外层 secant 的候选流程；内置精确点扫压预测器尚不能直接替代此实验适配器。证据为 `curve_outer_secant/`、`run_outer_secant.py` 与 `final_audit.json`。

3 次回退分别出现在 2.666667→2.766667、5.333333→5.433333、8→8.1 V，预测比例均触及 2 的上限；direct 分别耗费 29、30、27 次更新。它们都紧随为命中精确参考点而缩短的最后一步。该相关性提示下一轮可单独检验在上一步被参考点截短后禁用外推或调整比例保护，但本次没有验证此改动，不能认定为已证实的根因。预测器运行只触发 6 次内部重定心，正常步工作量进一步下降。

完整 8 点相对原基线的最大物理场差异：单独重定心 phin 约 8.77e-14 V、载流子密度相对差异约 3.40e-12；叠加 secant 的 phin/phip 最大差异约 7.58e-11 V、载流子密度最大相对差异约 2.94e-9（分母下限 1 m⁻³）。最大 KCL 比率分别为 1.80e-13 与 9.85e-14，均远低于 1e-8 原门限。

## 验证与边界

新增 Catch2 数值测试覆盖缩放及不缩放两种坐标下的密度、通量、残差与 Jacobian 等价性，以及 1e-20 V 低位保存、原本被舍入吞掉的校正重新可表达，共 40 条断言通过。配置测试检查默认关闭和 JSON 开启。

Release 全 CTest 首次并行运行 745/746 通过；失败的 dd 聚合测试报 VTK 标量缺失，单独串行重跑通过，符合聚合测试与独立注册测试并发写同名临时 VTK 的冲突特征。随后将完整 Release CTest 以 `--parallel 1` 重跑，**746/746 通过**。串行回归在全部计时曲线结束后执行，避免与求解计时争用 CPU。日志与 JUnit：`build-release/qf-recenter-ctest.log`、`build-release/build-release/qf-recenter-dd-recheck.xml`、`build-release/qf-recenter-serial-ctest.log`、`build-release/build-release/qf-recenter-serial-ctest.xml`。

本次只支持上述网格、D5 配置、Eigen SparseLU 后端和前 8 点的结论，不外推到量子求解、其他线性后端或 40 V 全曲线。新增开关保留默认关闭；内置精确点扫压的预测器历史更新与恢复异常为独立后续修复事项，本次仅定位和记录，未混入重定心实现。
