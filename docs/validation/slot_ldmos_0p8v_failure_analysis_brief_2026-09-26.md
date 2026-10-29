# SLOT 0.8 V 外压续算中止：供外部模型分析的故障说明

日期：2026-09-26（北京时间）。本次仅整理既有证据、读取代码与计算已有 CSV 的统计；
未重启求解器、未修改物理配置、数值门槛或生产代码。定时监测已暂停。

## 1. 需要分析什么

请分别分析两个问题，不要把它们合并成一个“收敛失败”：

1. **数值问题**：固定静电势的载流子密度 Newton 修正，第 4 次为什么找不到
   满足冻结权重下降条件的候选？应如何在保持最终资格门槛的前提下恢复？
2. **控制流问题**：一个候选生成失败为什么中止整个扫描，而没有进入后续
   QF 修正或有界缩步？怎样区分已知数值拒绝与真实环境、输入和证据故障？

第二项的异常传播路径已由代码和日志确认。第一项只有候选失败与停滞证据，
尚无第 4 步的方向导数、线性系统后向误差及逐次回溯 merit，不能认定根因。

## 2. 运行身份与物理背景

- Worktree：`D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next`。
- 分支：`codex/slot-ldmos-bvds-next`；HEAD：`3e34336f4ecb6589294372599797fdb8ede13792`。
  工作区有未提交修改，HEAD 本身不能完整代表本次可执行程序或工具。
- 主驱动：`scripts/run_slot_ldmos_bv_campaign.py`；运行目录：`build/s26full2`。
- 可执行文件：`build-release/vela_example_runner.exe`（Windows MSYS2 UCRT64 Release）。
- 冻结 runner SHA256：`c4ab5e76bafc85915cba90a775763fb2068646b8fd1897be82863a3fb33fcc6f`。
- 网格：13881 节点；24272 个 Si 三角形；36681 条 Si 输运边。
- 温度 300 K；source/gate/SLOT/substrate 物理偏压为 0 V；本段参考原点为 0 V。
- 从本地平衡态启动，目标是补齐新配置的完整 0→BV 路径。
  全程外电阻 `R=1e12 Ω·µm`，`Vout=Vin+R·Id`，电流为 `A/µm`。
- 低压使用 Python 标量负载线搜索；每个内部试探点先做独立器件 Newton。
  **本次失败不是高压耦合外电阻 QR 求解器失败**；后者仅在既定高压阶段使用。
- 当前受控物理：IALMob，完整场导数；自洽 van Overstraeten 雪崩，
  `ialmob_element_density_grad_qf`、`local_ad`；SRH、Auger 密度增强；恢复掺杂的
  PolySilicon 与 ohmic 接触。仍使用原生导出的输运、迁移率和雪崩几何。
- 几何：`slot_ldmos_external_averagebox` 输运边系数；电荷和复合体积
  `material_local_circumcentric`；介电耦合 `per_cell`。
- `scaling.mode=unit_scaling`：坐标 µm、浓度 cm⁻³、电压 V。
  若字段名保留 `_m3` 等历史后缀，不能据名称再乘一次单位换算。
  报告中的原始方程残差及 flux/source 数值保持输出口径，不能直接当作端电流 A/µm。
- 构建证据确认 `VELA_HAS_UMFPACK=1`、`VELA_HAS_SPQR=1`、HDF5 状态支持开启。
  本次 invocation 中 `VELA_LINEAR_SOLVER` 未设置；按当前代码，块探针的默认
  `LinearSolver` 选择 UMFPACK。失败步骤没有单独打印后端，以上是由编译宏与选择逻辑推导，
  不应描述为从失败步骤日志直接读到了后端名。
- OMP/OpenBLAS/MKL 线程数均为 1；只运行一个重负载求解器。

背景：[当前扫描说明](slot_ldmos_full_resistor_campaign_2026-09-26.md)。
此前独立高压支路已有 BV 判据结果，但不能替代本轮从 0 V 开始的连续路径。

## 3. 最后合格点与失败坐标

截至异常退出，顶层 ledger 有 **11 个接受状态**。

| 项目 | 数值 |
| --- | ---: |
| 最后合格外压 | 0.6999991911127921 V |
| 最后合格内压 | 0.6995407268647013 V |
| 最后合格 Id | 4.584642480908619e-16 A/µm |
| 最后合格行比最大值 | 2.3838442615555347e-8 |
| 有效行数量 | 电子 12634、空穴 12634 |
| 本次提议外压 | 0.799999191112792 V |
| 失败的内部器件试探偏压 | 0.7995407268647012 V |
| 外压提议跨度 | 0.1 V |

最后合格点经独立器件、有效行、全局闭合、状态绑定和负载线确认复核通过。
失败内压只是外电阻迭代的试探坐标，尚未完成负载线闭合，不能计入接受曲线。
失败 `iv.csv` 中电流字段为 0，是未收敛行的输出占位，**不是已计算的零漏电流**。

本次目录简称：

```text
RUN  = build/s26full2
TRY  = RUN/a0011/scalar/attempts/attempt_000003
DEV  = TRY/d
CORR = DEV/local_corrector
FAIL = CORR/step_004
```

## 4. 完整器件 Newton 先发生了什么

`DEV/iv.csv`、`DEV/attempts.csv` 和 `DEV/iv_newton_failure_diagnostics.json` 表明：

- 初始整体残差范数：`9.4095297465671912`。
- 执行 35 次 Newton 更新；诊断失败迭代号 36（不要混淆两种计数）。
- 末态整体缩放残差范数：`1.5731599281743357e-10`。
- 失败分类：`carrier_row_convergence_line_search_rejected`。
- 有 504 个载流子行超门槛；最大行比 `1.6556191428002671`。
- `global_continuity_closure_satisfied=1`，但局部行不通过，不能接受。
  CSV 中另外两个 `global_*_continuity_closure_qualified=0` 对应此点采用
  `absolute` 判据而非相对源量资格；不能用这些字段互相替代。

| 原始分块残差 L2（按输出口径） | 数值 |
| --- | ---: |
| Poisson | 1.8714703794213239e-10 |
| 电子连续性 | 2.7818975492451092e-15 |
| 空穴连续性 | 6.6003910869268932e-12 |

这些分块原始范数与整体缩放范数的尺度不同，不能横向作相同单位比较。
完整 Newton 的最后一次线搜索尝试了 13 个步长：

| 阻尼 | 试探范数 | 充分下降目标 | 结果 |
| ---: | ---: | ---: | --- |
| 1 | 2.8530114187091783e-10 | 1.5730026121815183e-10 | 拒绝 |
| 0.000244140625 | 1.5731628856281013e-10 | 1.5731598897671108e-10 | 拒绝 |

这 13 次为完整 Newton 的记录，**不是**下面密度块第 4 步的 20 次回溯记录。

## 5. 密度修正阶段：改变量、行比与停滞

驱动配置 `carrier_corrector_passes=8`、`carrier_corrector_mode=density_then_qf`。
原定第 1–4 次为 `carrier_density`，第 5–8 次为 `carrier_only`（QF 坐标）。
每个成功候选先做受守护的 Poly 平衡恢复，再审计所有有效行；只有审计通过才做
完整器件复验，最后再检查负载线。块探针 stdout 的 `converged=true` 只表示探针完成。

| 密度步 | 进程结果 | 电子超限行 / 最大比 | 空穴超限行 / 最大比 | 后续完整器件复验 |
| --- | --- | --- | --- | --- |
| 1 | exit 0 | 9 / 1.2261051911 | 452 / 1.6547995094 | 未触发，局部审计失败 |
| 2 | exit 0 | 9 / 1.2261057199 | 452 / 1.6547963529 | 未触发，局部审计失败 |
| 3 | exit 0 | 9 / 1.2261226589 | 453 / 1.6546951794 | 未触发，局部审计失败 |
| 4 | exit 1 | 无候选输出 | 无候选输出 | 未执行 |

表中统计来自 `step_00N/local_audit/summary.json`。
`source_qualified_rows=12306` 是附加统计，并非此处筛掉其余有效行：
`rows_above_eps` 和 `ratio_max` 在全部 12634 个有效行上计算。

| 密度步 | 未裁剪 step 范数 | 接受 step 范数 | Δφn 最大绝对值 (V) | Δφp 最大绝对值 (V) |
| --- | ---: | ---: | ---: | ---: |
| 1 | 15017.526626422004 | 4.927199381279499 | 4.4878965e-4 | 0.0595264293 |
| 2 | 2390874.2212181636 | 5.0403542078393855 | 1.7392691e-6 | 0.0595264293 |
| 3 | 89270.38088094101 | 6.06402158742517 | 5.5714892e-5 | 0.0595264293 |

step 范数是代码变量空间的欧氏范数；未裁剪值为从密度解转换后的 QF 线性方向，
不是物理密度变化的范数。三步 `Δψ=0`。空穴残差 L2 仅从
`6.6003911e-12` 降到 `6.5585693e-12`（约 0.634%）；同时 Poisson 残差由
`1.8714704e-10` 增至 `1.9153802e-10`。固定 ψ 不等于 Poisson 残差不变，
因为密度改变了电荷。

350 个 Poly 节点的恢复改变量第 1 步仅 `3.27e-23 V`，随后为零。
据当前证据，不能把显著扰动归因于这三次 Poly 投影。

第 3 步后两个最大行比位置：

| 载流子 / 节点 | (x,y) µm | 局部密度 cm⁻³ | residual | scale | ratio |
| --- | --- | ---: | ---: | ---: | ---: |
| 电子 / 649 | (1.280514294, 3.920546875) | 8.766910906e-9 | -1.806594320e-15 | 1.473420548e-15 | 1.226122659 |
| 空穴 / 6841 | (1.167105253, 2.564843750) | 1.372806942e-6 | -5.955480594e-15 | 3.599140596e-15 | 1.654695179 |

两处均为极低密度的少数载流子：649 的净掺杂为受主占优；6841 为施主占优。
对应局部 impact 项绝对值约 `2.32e-27` 和 `1.24e-26`，远小于此两行的输运/复合项。
这支持优先检查局部少数载流子修正，但不能据两个节点排除其他区域或完整雪崩耦合影响。
这些是**第 3 步后已记录的最差行**，不一定就是第 4 步拒绝的决定性节点。

## 6. 第 4 步的准确错误与数学判据

`FAIL/config.log`：开始 21:54:20.986；失败 21:54:38.853；耗时约 17.868 s。
`FAIL/stderr.log` 原文：

```text
vela_example_runner failed: carrier_density has no residual-decreasing candidate
```

第 4 步输入为 `step_003/poly_equilibrium_candidate.h5`。
没有生成 `step.csv`、`trial_carrier_density.h5`，stdout 为空。
不是执行 300 s 超时，也没有 DLL 缺失、非正基态密度或线性求解抛错的信息。

当前 `NewtonSolver::evaluateBlockStep` 实现（`src/solver/NewtonSolver.cpp:5066–5144`）：

1. 装配当前状态的原始残差 r、Jacobian J；抽取 2N×2N 载流子子块，冻结 ψ。
2. 对 Boltzmann 密度用列变换 D：电子 `-Vt/(potentialScale*n)`，
   空穴 `+Vt/(potentialScale*p)`；求 `(Jcc D) δd = -rc`，再形成 QF 线性方向。
3. 固定当前状态每行权重：

   `wi = 1 / max(flux_abs_sum_i, |recomb_i|, |impact_i|, 1e-30)`，
   非有效行权重为零；`merit(r)=||W rc||₂`。

4. 尝试 20 次 `alpha=2^(-k), k=0…19`（最小约 `1.90735e-6`）。
   将密度相对变化下限裁到 `-0.9`，再通过 `log1p` 映射回 QF；
   再应用配置的 QF 更新限制。
5. 只有 `merit(trial) <= merit(base)*(1-1e-4*alpha)` 且有限才接受；否则继续。
   20 次全不满足则抛出上述 `runtime_error`。

**注意**：配置虽为 `componentwise_poisson_recorrect`，本块探针调用的是
`applyConfiguredStepCaps`，没有调用完整 Newton 的 Poisson 重新修正函数。
实测三步 Δψ=0 与此一致。不要假设块探针自动执行完整 Newton 的全部恢复机制。

此处 merit 是冻结权重 L2；局部验收是重新计算尺度后的最大行比。
因此“探针接受一步”不保证最大行比下降，更不保证完整器件合格。
源核使用完整导数的配置本身，也不能证明这次状态的 Jacobian 方向差分或线性解精度正确。

## 7. 为什么单次候选失败变成全任务中止

已确认调用链：

```text
bv_campaign.backend_for: 8 passes, density_then_qf
  → loadline.SolverBackend: 将 helper.execute_config 替换为 self._execute
  → noavalanche.correct_local_candidate: 第 4 步调用 execute_config
  → C++ evaluateBlockStep: 20 次候选均被拒绝，throw runtime_error
  → runner main catch: stderr + exit 1；写 CSV 的位置在 evaluateBlockStep 之后
  → SolverBackend._execute: code != 0 且 output_csv 不存在
  → raise RuntimeError("runtime_environment_or_invalid_input: ...")
  → scalar Campaign.run: 保存 failure 后重新抛出普通异常
  → bv_campaign 外层捕获: environment_or_evidence_failure，结束 campaign
```

- `scripts/run_slot_ldmos_loadline_continuation.py:537–546`：
  分类条件只检查非零退出码与缺失输出，未识别这个已知数值拒绝。
- `scripts/run_slot_ldmos_noavalanche_continuation.py:419–433`：
  原本有收到非零退出码后 break 的处理，但包装器先抛异常，未返回到该分支。
  因而 `CORR/summary.json` 只包含前三步，不包含第 4 步的退出记录。
- 同 helper 的 `density_then_qf` 是按迭代序号切换；即使简单返回非零退出码，
  现有 break 仍不会自动跳到 QF 阶段。因此单改“分类”为数值失败，不等于已实现 QF 回退。
- 内部恢复和顶层缩步只捕获 `NumericalFailure`；本次为普通 `RuntimeError`，
  所以本目标没有走到数值失败缩步处理。之前 0.7 V 目标确实出现过正常缩步，
  不能说整个工具没有缩步功能。

所以“环境/输入失败”是控制器的标签，**不是已经证实的环境根因**。
终止时两个进程均退出；未观察到用户中断或电脑重启证据。

## 8. 关键门槛与不应更改的条件

| 配置 | 值 |
| --- | --- |
| 完整 Newton abstol / reltol / stall_residual_floor | 1e-9 / 1e-8 / 1e-8 |
| 完整 Newton 最大迭代 | 40 |
| carrier_row eps / scale_floor / min_source_scale | 1e-3 / 1e-30 / 0 |
| min_carrier_density（输入 TCAD 数值） | 1e-200 |
| residual_scales ψ / φn / φp | 1.189625 / 8476.074054 / 3518.300497 |
| QF majority / minority 更新限制 | 0.1 / 0.05 V |
| 全局闭合 tolerance / source_floor / absolute_tolerance | 0.01 / 1e-11 / 1e-12 |
| 全局稳定化 / 恢复 | replace_free_row_sum / final_linear_qf_projection |
| 密度探针 carrier_diagonal_floor | disabled |
| 负载线容差 | 1e-6 V |
| 单次子进程超时 / 总调用预算 | 300 s / 28800 s |

不能通过关闭有效行检查、忽略少数载流子、放宽门槛、把失败电流占位当证据或
修改物理参数“解决”当前错误。候选恢复后的接受必须仍通过完整器件、全部有效行、
全局闭合、状态绑定和负载线确认。

## 9. 请外部模型回答的具体问题

1. 在固定 ψ、Boltzmann 坐标变换与冻结权重条件下，什么条件保证足够小步长下降？
   逐分量裁剪、90% 耗尽限制、极低密度列缩放、浮点抵消或 Jacobian 不一致，
   各会怎样破坏该条件？请区分理论条件与这里已有证据。
2. `raw_step_norm` 从 1.50e4 到 2.39e6，随后回到 8.93e4，实际最大空穴更新却
   连续三次约 59.526 mV，是否提示耗尽限制主导？还需要哪个最小诊断才能证实？
   `Vt ln(10)` 与此量级相近仅为线索，不证明第 4 步根因。
3. 应先切到 QF 块、退回合格父态缩小偏压，还是做 Poisson/载流子交替修正？
   请给出有界状态机、候选回退点、停止条件和完整复验方式，避免无限循环。
4. 如何给 runner 输出结构化的“候选未找到”结果，使其可恢复，
   同时真实 DLL/输入/文件损坏/证据不完整仍必须停止？不要把所有 exit 1 一律降级。
5. 哪些最小回归能验证：已知数值失败可缩步；环境错误仍停止；失败步记录不丢失；
   从 QF/密度恢复得到的候选不绕过完整器件资格？
6. 建议诊断优先级：第 4 步 `(Jcc D)δd+rc` 的归一化后向误差、
   `W Jcc d` 与方向差分、裁剪前后下降性、回溯 merit/有限性、最差节点贡献。
   哪一项能最少重算地排除最大范围的假设？

尚未执行这些实验。尤其没有第 4 步的 J/δx 捕获、条件数、裁剪计数、逐次 merit，
不能声称已证明病态、舍入底限、非下降 Newton 方向或物理导数错误。

## 10. 可供复查的文件与附件

所有相对路径相对于本 worktree。重点证据：

- `build/s26full2/checkpoint.json`、`a0011/scalar/checkpoint.json`：顶层/标量异常与合格 ledger。
- `DEV/iv.csv`、`attempts.csv`、`iv_newton_failure_diagnostics.json`、`iv_newton_history.csv`。
- `CORR/summary.json`；前三步 stdout、局部审计 summary/terms；第 4 步 config/log/stderr。
- `FAIL/config.json` 引用的 HDF5 输入即第 3 步 Poly 候选。
- 相关源代码：`NewtonSolver.cpp`、`LinearSolver.cpp`、`vela_example_runner.cpp`，
  三个 `run_slot_ldmos_*` 驱动和 `audit_slot_ldmos_local_continuity.py`。

轻量分析包：`build/s26_failure_brief_20260926/slot_0p8v_failure_review.zip`。
包含本报告、关键原始日志/配置/诊断、裁剪后的热点行、源代码相关摘录、
输入/运行身份 SHA256 清单和结构化统计。未复制全套网格、可执行文件或 HDF5，
因此**它是离线分析包，不是可以脱离本工作区直接运行的复现包**。
原始状态和完整大 CSV 保持在原位置，未改写。

最小重放对象是 `FAIL/config.json`；其输出仍指向原证据目录。
未来若授权重放，必须先复制配置并把所有输出/诊断路径改到新目录，保留原输入状态和
资源哈希，并仅启动一个 Release 进程；不得直接覆盖本次失败现场。
