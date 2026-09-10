# LDMOS Stage 4 交接：边缓存优化与 Vg=8 全曲线完成（2026-09-10）

本轮任务已完成：边通量导数的重复计算优化已提交，Release 742 项测试通过，
优化版 Vg=8、Vd=0–40 V 全曲线通过原单栅压工程及最终判据。
当前没有本次待续跑的仿真任务。下一对话先读取本文件和当前 AGENTS.md，
再根据用户的新请求选择工作；下方建议不是自动执行队列。

## 1. 工作目录、提交和运行状态

- **工作树**：`D:/code-repo/vela-tcad/.worktrees/templates-ldmos-phase-a`
- **分支**：`codex/templates-ldmos-phase-a`。新聊天即使默认目录为仓库根目录，
  也应先确认分支和工作树，再在以上目录执行本任务命令。
- **代码及验证报告提交**：`b6552d0914a69f425c9afabb0f52980dd738209b`
  （`perf(dd): reuse edge mobility, geometry and endpoint states`）。本交接文档作为后续独立文档提交。
- 基础提交为 `d4285815f648be3f0b20719cbbe4d33e6a1b38b8`，是此前合并 main 的提交；
  前序精度修复和推进策略在 `8fc0fdd`，前序调查/全曲线记录在 `64c923d`。
- 整理时本地 main 为 `4938999e8eae2a4318e5ac412575efe5a6ee9332`，尚有 1 个提交
  不在本分支。本轮用户要求本地提交及交接，未要求再次合并 main 或推送远端。
  若后续合并/改变源码，需重新选择相应验证，不能直接沿用旧程序的验收结论。
- 全曲线控制器 PID 18976 已退出，本次 `candidate.exe` 无活跃进程。
  完成时 `ldmos-d5` 监测已由工具暂停；交接时其本地 automation.toml 已不存在，
  不应依旧 PID/旧监测配置重启任务。后续如需监测，应查询当时实际任务配置。

## 2. 已完成的实现及精度约束

修改文件为 [CoupledDDAssembler.h](../../include/vela/equation/CoupledDDAssembler.h)、
[CoupledDDAssembler.cpp](../../src/equation/CoupledDDAssembler.cpp) 和
[test_newton_solver.cpp](../../tests/test_newton_solver.cpp)。

1. **体区迁移率复用**：vector-QF 高场模式中，即便全局开启接触电场回退，
   无实际接触回退且无表面迁移率的体区边也能复用基态迁移率；给定迁移率时
   跳过重复驱动场求值。接触边及其第三节点 psi 扰动仍重新计算电场/迁移率，
   原解析高场反馈保留。
2. **几何与当前梯度复用**：初始化时缓存输运接触掩码、面积和 stencil 几何敏感度；
   每次 Jacobian 装配重新计算电子/空穴实际单元 QF 梯度，并在基态字段及解析反馈间共享。
   实际梯度依赖当前状态，不能跨 Newton 迭代冻结。
3. **端点状态缓存**：按当前边、载流子、端点区分，键包含 psi、局部 QF 和相对 QF
   的 double 位模式；每端最多 8 项，换边/重新装配即失效，满容量正常重新计算。
   电子 `electronDensityAt` 与 `Nc*FermiHalf(eta)` 各自保留原表达式和舍入顺序，
   不合并两条密度路径。原 QF 参考值/增量、亚 ULP 精度、通量稳定公式及差分步长不变。
4. 新测试 `Fermi vector HFS Jacobian follows changing states on bulk and contact edges`
   覆盖体区/接触回退，改变并恢复状态，核对 Jacobian 与残差差分及缓存生命周期。

SparseLU 生命周期、外推、物理模型与门限均未在本轮修改。
此前 26.6667 V 附近的精度停滞背景见旧精度调查；本轮最新证据是优化版已通过
该点的参考系平移并完整达到 40 V，不能再将旧停止状态当成当前状态。

## 3. 验证结果：正确性成立，整体加速尚未确认

### 短区间配对：Vg=8、Vd=0–1.33333333333333 V

两版均 45 个子步骤、230 次 Newton 更新、0 回退；24 个接受状态逐字节一致。

| 指标 | 优化前 | 优化后 |
| --- | ---: | ---: |
| 子进程墙钟 | 400.783 s | 410.621 s |
| 子进程 CPU | 275.328 s | 266.141 s |
| Jacobian 阶段 | 108.074 s | 90.067 s |
| 边物理阶段 | 71.446 s | 49.811 s |
| 边阶段 FermiHalf 调用 | 173,689,980 | 95,551,170 |

局部边阶段减少 30.28%，Fermi 调用减少 44.99%；整体墙钟增加 2.45%，CPU 减少 3.34%。
ABBA/BAAB 同检查点短测每版 4 次，状态完全相同，边阶段中位数减少 28.72%，
整次墙钟中位数仍未改善。未修改阶段也存在明显耗时波动；不能将局部收益外推为整段加速。

### 优化版全曲线：Vg=8、Vd=0–40 V

- 2026-09-10 09:47:32–11:46:04（Asia/Shanghai），控制器 **7111.840 s，约 1 h 58 min 32 s**。
- 子进程墙钟 **6908.023 s**，子进程 CPU **5193.203 s**；这三种时间口径不要混用。
- **31/31 精确参考点、429 次推进、431 个唯一合格状态、859 个子进程、5209 次 Newton 更新、0 回退**。
- 对 Sentaurus 电流误差：中位数 **1.330933%**、P95 **1.612869%**、最大 **1.627051%**；
  低压电阻误差 **1.032119%**，40 V 端点误差 **1.099059%**。
- 40 V 电流：Vela `4.2363315296624625e-4 A/μm`；Sentaurus `4.19027790569749e-4 A/μm`。
- 31 点最大归一化 KCL **1.372946e-10%**，分母为最大绝对端电流。
- **26.6667 V 平移 28 V 的检查合格**：最大电势/QF 差 7.11e-15 V，密度相对差 2.14e-13 以下，
  电流相对差 3.64e-14；四项原门限全部通过。
- 单栅压工程/最终五项判据全部通过。未合并旧 Vg=4 数据，未证明本版本双栅压 D5/电流比通过。
- 内部累计热点：数值分解 2106.529 s（30.96%），Jacobian 1499.598 s（22.04%），
  残差 1006.510 s（14.79%），边物理 865.031 s（12.71%，已包含在 Jacobian 内）。阶段不可相加。
- 全曲线没有优化前配对，因此尚无整体加速比例；也没有声称与优化前全曲线逐状态相同。

## 4. 程序身份、物理配置与证据入口

- 冻结优化版程序：`reference_staging/templates_ldmos_edge_cache_20260910/candidate.exe`
- SHA256：`06301e6d44359dffa8889262524401b499c8a619cdb6b2bcb462d8c2f0734b6e`
- 全曲线目录：`reference_staging/templates_ldmos_edge_cache_full_20260910/vg8_full/`
- 运行脚本：`reference_staging/templates_ldmos_edge_cache_full_20260910/run_full.py`
  （完整 31 点循环和参考系平移；不能用短区间脚本只调用 advance(40) 来替代）。
- 独立完成复核：`reference_staging/templates_ldmos_edge_cache_full_20260910/verify_completed.py`。
  20 项冻结程序/源码/输入哈希、431 个状态、859 个子步骤及性能文件均复核通过，曲线指标/KCL重算一致。
- 局部配对目录：`reference_staging/templates_ldmos_edge_cache_20260910/`；首轮 gprof 目录：
  `reference_staging/templates_ldmos_gprof_20260909/`。

这些 `reference_staging` 目录、二进制及仿真输出均为本机忽略文件，未加入提交；
报告已提交不意味着全新 clone 会带有原始证据。保留该工作树及相关历史脚本依赖。
运行器拒绝覆盖既有输出；若用户要求新实验，应创建独立目录并冻结对应程序和输入。

物理配置：Vg=8 V、300 K、Fermi–Dirac、OldSlotboom BGN、SRH/Auger、
constant_field / transport_cell_vector QF 高场迁移率及接触电场回退，表面迁移率、雪崩、量子、热关闭。
网格 10241 节点、19782 三角形、30022 边、5723 硅节点，external average-box 输运系数及
barycentric 节点体积。几何 μm、密度 m⁻³、电势 V、报告电流 A/μm。
实际后端 Eigen SparseLU/COLAMD + L2 行列均衡，内部计时开启、gprof 关闭、外推关闭。
原块门限 psi/electron/hole=5e-8/1e-11/3e-10，局部行及 KCL 比值门限=1e-8。

建议阅读顺序：

1. [优化版全曲线验收报告](templates_ldmos_stage4_edge_cache_full_2026-09-10.md)
2. [缓存实现与低压配对报告](templates_ldmos_stage4_edge_cache_2026-09-10.md)
3. [边导数重复计算与 SparseLU 核对](templates_ldmos_stage4_edge_derivative_analysis_2026-09-10.md)
4. [首轮 gprof 报告](templates_ldmos_stage4_gprof_0_1p3333_2026-09-09.md)

## 5. 后续可选工作及验收边界

已完成用户本轮要求，没有需要自动续跑的剩余任务。用户如继续要求性能优化，可按以下线索选择独立变更：

1. **Poisson 再校正求解器生命周期**：`src/solver/NewtonSolver.cpp` 的
   `recorrectPoissonStepForClippedQuasiFermi` 仍每次局部创建 `LinearSolver`。
   主 Newton 循环的 `LinearSolver::analyzePatternIfNeeded` 已缓存模式；固定网格
   不等于所有子系统共用一次分析，边界/结构变化还需失效处理。前序短区间主系统
   251 次求解只分析 44 次，Poisson 则 37 次求解分析 37 次。
   可研究让 Poisson 缓存持续于一次 Newton solve，同时保持与主矩阵的缓存独立。
   符号分析远小于数值分解成本，不应期待仅此改动解决主要热点。
2. **数值分解成本与性能测量稳定性**：先固定程序、输入、计时开关、工作量和运行顺序，
   同时观察墙钟/CPU/未修改阶段。替换后端、排序或跨偏压保持进程属于新的实验，
   需独立验证，不能与本轮缓存结果混合归因。
3. **残差与连续性诊断同状态复用**：后续可调查缓存收益，但状态、行权重、接触及模型
   改变时的失效条件需证明。外推、模型修改和门限变化应按用户另行指定的范围处理。
4. **Vg=4/双栅压验收**：只有用户要求扩展验证范围时再跑本版本，不能直接采用历史曲线替代。

未来涉及求解器源码时遵循当前 AGENTS.md，新增必要数值测试并运行相关回归；
不能为通过而放宽原门限、改变物理默认参数或合并不同参考系密度表达式。

## 6. 新对话操作提示

Windows 默认工具链为 `D:/msys64/ucrt64`；从上述工作树执行：

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git -c core.fsmonitor=false status --short
git branch --show-current
```

若源码发生变化，需要构建/验证时使用匹配的 preset：

```powershell
cmake --preset windows-ucrt64-release
cmake --build --preset windows-ucrt64-release --parallel 2
ctest --preset windows-ucrt64-release --output-on-failure
```

本轮已完成的 742 项测试日志在 `reference_staging/templates_ldmos_edge_cache_20260910/ctest.log`，
JUnit 在 `build-release/edge-cache-20260910-ctest.xml`；整理提交时再次核对源码与受测程序清单一致，
没有仅为文档重复运行求解器测试。运行 `verify_completed.py` 还会检查冻结源码哈希，
之后若主动修改源码，哈希失败须区分源码漂移和原始仿真结果问题。

PowerShell 的 Get-CimInstance 曾拒绝访问，进程活性可用 Get-Process 路径/启动时间/CPU与日志共同判断。
Git 索引在主仓库 `.git/worktrees/` 下，写入可能需要正常的沙箱权限提升；本地提交已获用户授权。
