# LDMOS 前 8 个参考点：连续性项复用与性能验证（2026-09-10）

本轮按用户要求将曲线范围限定为 Vg=8 V、Vd=0–9.33333333333333 V，
包含 0 V 在内的前 8 个精确参考点。已完成连续性项缓存实现、Release 743 项测试、
固定检查点交错对照和前 8 点完整配对审计。两版的 108 个接受状态逐字节一致；
本次配对子进程总墙钟减少 **5.42%**，CPU 减少 **5.48%**，均保持 0 回退。

## 范围与程序身份

- 工作树 `D:/code-repo/vela-tcad/.worktrees/templates-ldmos-phase-a`，
  基础 HEAD `d8ff26804a59735ff3f446846e446e0a62add914`，本报告对应其上的未提交改动。
- 基线程序为上一轮已经完成全曲线的边缓存版，SHA256
  `06301e6d44359dffa8889262524401b499c8a619cdb6b2bcb462d8c2f0734b6e`。
- 本轮候选 SHA256
  `e0b4193efc70cb7ce8072240afbc6eda840dbc24b0df21b0757503fb8f1e854f`。
- MSYS2 UCRT64 GCC 16.2.0、C++20、Release `-O3 -DNDEBUG`，LTO/gprof 关闭。
  实际后端 Eigen SparseLU/COLAMD，L2 行列均衡；配置及生成的构建文件确认
  HDF5、SPQR、UMFPACK 可用，本算例未采用 SPQR/UMFPACK。
- 网格仍为 10241 节点、19782 三角形、30022 边、5723 硅节点；
  external average-box 输运系数与 barycentric 节点体积。
  几何 μm、密度 m⁻³、电势 V、电流 A/μm。
- 300 K、Fermi–Dirac、OldSlotboom BGN、SRH/Auger、constant_field /
  transport_cell_vector QF 高场迁移率及接触电场回退；表面迁移率、雪崩、量子、热关闭。
- 沿用合格 0 V 种子；初始子步 0.0025 V、最大子步 0.1 V、增长因子 1.35，
  保留原严格同偏压再求解，外推关闭。原块门限 psi/electron/hole 为
  5e-8/1e-11/3e-10，局部行和 KCL 比值门限为 1e-8。

8 个参考电压为 `0, 1.33333333333333, 2.66666666666667, 4,
5.33333333333333, 6.66666666666667, 8, 9.33333333333333 V`。
运行器逐一到达这些点，在第 8 点验收后停止；不要求 40 V 端点、
26.6667 V 平移仿真或双栅压电流比，不将区间结果表述为全曲线资格。

## 实现与失效边界

`NewtonSolver::solveClassicalWithFrozenElectronQuantumPotential` 在完成装配器、
QF 参考系、冻结量子势及边界配置后创建一个局部 `ContinuityTermCache`。
它共享行权重计算与局部载流子行验收请求的连续性方程项，消除同状态重复求值。

- 缓存只保留一个完整状态，以原始 double 位模式比较；亚 ULP 增量和有符号零
  不被合并。切换试算状态会重新求值，恢复旧状态同样按当前缓存内容检查。
- 装配器及其物理上下文在一次 classical Newton solve 内固定；下一次求解或
  量子外迭代重新创建缓存。边界以值保存。未来若在此作用域中修改参考系、量子势、
  网格或模型，必须先失效/重建缓存，不能继续使用现有生命周期假设。
- 复用的是未加权方程项；行权重仍使用当前配置独立计算。
  全局连续性闭合使用无边界替换的不同语义，仍独立求值。
  可观测雪崩、反馈替换诊断及其他诊断 API 均未接入此缓存。
- 重算前先失效，防止异常后留下旧命中；非有限输入状态不进入有效缓存。
  原密度表达式、舍入顺序、残差、Jacobian、线搜索、步长、门限均未修改。
- 新增 `newton.continuity_terms_cache_hits/misses` 计数器，便于核对工作量。

本轮没有修改 Poisson 再校正求解器生命周期或 SparseLU 数值分解策略。

## 排序实验与固定检查点对照

先以独立 `build-release-ldmos-amd` 构建测试已有的 AMD 排序选项。
同一低压检查点 COLAMD 预热运行约 4.150 s 完成，AMD 达到 180 s 超时后被终止，
因此未进入后续 AMD 曲线实验。该结果是运行时间筛选失败，不构成 AMD 数值错误结论。
本轮保留 COLAMD；UMFPACK 在当前代码中只有独立辅助入口，未接入主 Newton 后端。

缓存候选的每组短测先对两版各预热一次，再执行 ABBA/BAAB（每版 4 次）。
同一检查点固定父状态和全部原配置，只替换程序和输出路径；关闭内部计时的实验
对两版同时关闭。总计 80 次缓存对照运行，全部输出状态与原检查点逐字节一致，
返回码、收敛标志、Newton 更新数保持一致。直接求解的非零返回码对应原流程中
需要严格再求解的步骤，比较的是其原始拒绝终态，不将其计为已合格参考点。

| 检查点/阶段 | 内部计时 | 基线墙钟中位数 | 候选墙钟中位数 | 降低 |
| --- | --- | ---: | ---: | ---: |
| 0.01043125 V / direct | 开 | 4.171 s | 3.939 s | 5.57% |
| 4 V / direct | 开 | 7.720 s | 7.444 s | 3.58% |
| 9.33333333333333 V / direct | 开 | 8.330 s | 7.104 s | 14.71% |
| 4 V / reclose | 开 | 1.913 s | 2.077 s | −8.61% |
| 9.33333333333333 V / reclose | 开 | 1.943 s | 1.822 s | 6.21% |
| 0.01043125 V / direct | 关 | 4.120 s | 4.025 s | 2.31% |
| 4 V / direct | 关 | 7.535 s | 7.128 s | 5.41% |
| 9.33333333333333 V / direct | 关 | 6.621 s | 6.178 s | 6.70% |

局部连续性诊断调用在低压 direct 从 12 次降至 6 次，4 V direct 从 20 次降至 9 次，
两次 reclose 均从 4 次降至 2 次；未改动阶段和很短的求解仍有计时波动。
不能把表中的单个降幅直接外推为区间加速，也不能跨不同时段的计时开关实验
简单相减来估算 profiler 开销。

## 前 8 点配对验收

两条曲线独立保存父状态和推进账本，以逐子进程 AB/BA 顺序串行求解；
每对子步骤即时核对目标、父状态 SHA256、返回码及 Newton 更新数。
正式运行双方均开启内部计时，以便核对消除的计算；未同时运行构建或其他性能实验。
每条控制器墙钟包含等待另一条曲线的时间，不作为单条曲线速度指标。
性能比较采用各自子进程墙钟/CPU 总和；配对整体历时单独记录。

| 指标 | 基线 | 候选 |
| --- | ---: | ---: |
| 精确参考点 | 8 | 8 |
| 物理推进 / 唯一接受状态 | 107 / 108 | 107 / 108 |
| 子进程 / Newton 更新 / 回退 | 214 / 1284 / 0 | 214 / 1284 / 0 |
| 子进程墙钟合计 | 1306.019 s（21 min 46 s） | 1235.173 s（20 min 35 s） |
| 子进程 CPU 合计 | 1106.625 s | 1045.969 s |
| 内部 dc_sweep.total | 1290.158 s | 1219.333 s |
| 连续性诊断调用 | 3288 | 1525 |
| 连续性诊断耗时 | 132.657 s | 60.836 s |
| 残差 / Jacobian / 数值分解次数 | 4528 / 1390 / 1612 | 4528 / 1390 / 1612 |

每条曲线节省子进程墙钟 **70.846 s**。两条曲线及原审计、性能汇总的配对协调器
历时为 2598.940 s（43 min 18.940 s）；该值不是单条曲线耗时。
两个子控制器约 2598 s 的计时均包含等待对方，不能拿来计算候选加速。

缓存命中 1763 次、未命中 1525 次，精确解释消除的 1763 次诊断调用（减少 53.62%）；
诊断阶段时间减少 54.14%。行权重调用仍为 1606 次，局部行验收调用仍为 1682 次，
原检查均保留。线搜索试算 4312 次、接受更新 1284 次、未接受阶段 106 次，
实际符号分析 429 次、分析缓存命中 1183 次，以及边物理 FermiHalf 调用
529147968 次在两版间均相同。0 回退不表示每次 direct 子进程均一次收敛，
原严格再求解的全部成本均已计入。

未修改阶段仍有波动：数值分解 392.540→397.748 s，Jacobian 291.375→290.518 s，
残差 187.375→184.408 s。阶段计时包含嵌套，不能把诊断、行权重和局部验收行相加。
本次完整区间的节省与减少的重复诊断成本一致；结论限定于本机、该输入和该测量协议，
不声称每个很短的子进程都更快，也不外推 0–40 V 加速比例。

### 区间物理判据

7 个非零电流参考点参与相对误差统计，全部 8 点参与完整性与 KCL 检查。
以下适用于该区间的原工程/最终门限均通过，两版指标相同：

| 指标 | 实测 | 工程门限 | 最终门限 |
| --- | ---: | ---: | ---: |
| 电流相对误差中位数 | 1.594856% | 15% | 5% |
| 电流相对误差 P95 | 1.624118% | 25% | 12% |
| 低压电阻误差 | 1.032119% | 20% | 10% |
| 最大归一化 KCL | 9.721110e-12% | 1% | 0.1% |

最大电流误差 1.627051%；9.33333333333333 V 区间末端电流为
Vela `3.12466896071033e-4 A/μm`、Sentaurus `3.07561730866064e-4 A/μm`。
这是区间末端读数，不是 40 V 端点验收。零偏压原始端电流及绝对 KCL 均为零。

## 验证与证据

- 新增数值测试覆盖 Fermi/vector-HFS、标量与节点参考场、不同接触边界、
  亚 ULP 状态变化和恢复、异常后恢复；连续性项对照独立残差及未缓存诊断。
  单独运行通过，512 项断言。
- Release 构建通过，完整 CTest **743/743 通过**，无失败，总耗时 106.89 s。
  已有参考系平移、精度、密度梯度和数值回归包含在本次测试内。
- 两条曲线的原完整性审计通过；独立复核重算电流误差和 KCL，核对冻结文件、
  428 个子步骤的程序/配置/父状态/性能文件身份、108 个逐字节一致的接受状态、
  全部子步骤工作量及上述未修改阶段计数。428 个执行票据与起止时间核对通过，
  确认两版子进程未重叠运行。
- 收尾复核确认当前 121 份源码/构建输入与受测冻结快照逐项哈希一致，当前构建程序
  与冻结候选程序 SHA256 一致，JUnit 确认 743 项通过、0 跳过。使用仓库原有
  Git 换行配置的 `git diff --check` 通过；MSYS2 Git 默认配置产生的 CRLF 显示差异
  未作为源码变化处理，也未为此改写受测文件。
- 本地实验根目录为 `reference_staging/templates_ldmos_first8_perf_20260910/`，
  其中二进制、源码快照和仿真输出均为忽略文件，不提交生成物。

证据入口（依赖保留本机工作树及历史检查点工具）：

- [基线冻结清单](../../reference_staging/templates_ldmos_first8_perf_20260910/baseline_binary/manifest.json)
- [候选冻结清单及源码快照](../../reference_staging/templates_ldmos_first8_perf_20260910/candidate_binary/manifest.json)
- [AMD 超时记录](../../reference_staging/templates_ldmos_first8_perf_20260910/ordering_replay/timeout.json)
- [计时开启检查点](../../reference_staging/templates_ldmos_first8_perf_20260910/cache_replay/summary.json)
- [计时开启直接求解补测](../../reference_staging/templates_ldmos_first8_perf_20260910/cache_direct_replay/summary.json)
- [计时关闭检查点](../../reference_staging/templates_ldmos_first8_perf_20260910/cache_replay_profiling_off/summary.json)
- [前 8 点运行适配器](../../reference_staging/templates_ldmos_first8_perf_20260910/run_first8.py)
- [配对协调器](../../reference_staging/templates_ldmos_first8_perf_20260910/run_pair.py)
- [独立完成比较](../../reference_staging/templates_ldmos_first8_perf_20260910/compare_first8.py)
- [前 8 点比较与完整性结果](../../reference_staging/templates_ldmos_first8_perf_20260910/first8_pair/comparison.json)
- [候选区间评分](../../reference_staging/templates_ldmos_first8_perf_20260910/first8_pair/candidate/score/summary.json)
- [候选原完整性审计](../../reference_staging/templates_ldmos_first8_perf_20260910/first8_pair/candidate/audit_summary.json)
- [配对协调器计时](../../reference_staging/templates_ldmos_first8_perf_20260910/first8_pair/timing.json)
- [源码与测试复核脚本](../../reference_staging/templates_ldmos_first8_perf_20260910/verify_sources.py)
- [源码、程序和测试最终复核](../../reference_staging/templates_ldmos_first8_perf_20260910/completion_verification.json)
- [CTest 日志](../../reference_staging/templates_ldmos_first8_perf_20260910/ctest.log)
- [CTest JUnit](../../build-release/first8-continuity-cache-ctest.xml)
