# Templates/LDMOS Stage-4 D5 交接（2026-09-09）

本文件应用户“在新聊天继续当前算例验证对比工作”的请求编写。现场核验时间：2026-09-09 08:13–08:18，Asia/Shanghai。它是此刻的工作交接；新聊天仍应先读项目当前参考并核验现场，不能把过时日志里的命令直接当作任务队列。

## 1. 当前结论与接续重点

**D5 整体验收尚未通过。当前首要工作是修复二维梯度 HFS 连续性 Jacobian 的迁移率反馈导数遗漏，再做同初态局部性能／电流回归，最后完成 Vg=8 与 Vg=4 两条完整曲线和两栅压电流比验收。** D5 合格前不进入 IALMob。

- 旧 `edge_projection` 算法：Vg=4 的 0→40 V、31/31 参考点已完成并通过单栅压 final 门限；Vg=8 完整曲线已数值收敛，但电流偏大，精度不通过。
- 新 `transport_cell_vector` 算法：Vg=8、Vd=4 V 的电流误差从 12.40768% 降至 1.52738%，已有原门限下自洽解。完整曲线仍未验证。
- 2026-09-08 晚局部 HFS 开／关和 198 个 JVP 方向证明：vector HFS 的连续性 Jacobian 没有完整反映准费米梯度对迁移率的依赖。该源码修复**尚未实施**。
- **昨晚全曲线队列目前没有进程继续执行。** 最后接受 Vd=10.533333333333326 V；原 JSON 仍显示 `running`，不能据此报告“仍在运行”。停止原因未查明，不能推断为数值不收敛、关机或人工终止。
- 当前用户请求仅为整理交接。本轮未重启仿真、改写旧状态文件、修改求解器或提交 Git。

## 2. 工作树、构建与未提交内容

```text
Worktree: D:\code-repo\vela-tcad\.worktrees\templates-ldmos-phase-a
Branch:   codex/templates-ldmos-phase-a
HEAD:     959c1c91513c62caeb700c2f71d8cef194879419
```

最早交接中的 HEAD `345fda3` 已过时，不要 reset 回去。当前工作区已有未提交的 Newton 迭代感知步长控制修改，必须保留：

```text
docs/config_schema.md
include/vela/simulation/DCSweep.h
include/vela/simulation/DCSweepStepControl.h
src/simulation/DCSweep.cpp
tests/test_dc_sweep.cpp
```

上述五文件合计 422 行增加、11 行删除。不是本轮 HFS Jacobian 修复；当前 CoupledDDAssembler 尚未为该问题修改。另有以下未跟踪报告，亦需保留：

```text
docs/validation/templates_ldmos_stage4_matched_step_validation_2026-09-07.md
docs/validation/templates_ldmos_stage4_two_gate_validation_2026-09-07.md
docs/validation/templates_ldmos_stage4_vg4_step_growth_ab_2026-09-08.md
docs/validation/templates_ldmos_stage4_vg8_rerun_2026-09-08.md
docs/validation/templates_ldmos_stage4_current_error_localization_2026-09-08.md
docs/validation/templates_ldmos_stage4_hfs_on_off_numerics_2026-09-08.md
```

此交接文件也为新增未提交文件。输出和试验脚本主要存放在 ignored `reference_staging/`，Git diff 不会包含这些证据；不要清理该目录。

Windows 工具链固定为 MSYS2 UCRT64。性能实验使用 Release：

```powershell
Set-Location 'D:\code-repo\vela-tcad\.worktrees\templates-ldmos-phase-a'
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git branch --show-current
git rev-parse HEAD
git -c core.fsmonitor=false status --short
# 实施源码修复后再构建；仅阅读本交接不需要构建。
cmake --preset windows-ucrt64-release
cmake --build --preset windows-ucrt64-release
```

当前已验证 runner：`build-release/vela_example_runner.exe`，SHA-256 为
`9638411ec0d38434f05afc765d573060adde0d406c6818f3428c76777a066391`。
实际后端为 Eigen SparseLU/COLAMD、L2 行列平衡。Python 使用 `D:\msys64\ucrt64\bin\python.exe`。

修复后必须使用新输出目录、记录新 runner SHA 和源码 diff；不要绕过旧实验的哈希断言，或将新旧二进制混入同一曲线后当作同一性能实验。旧 runner 如需 A/B 使用，应在重建前保存到新的 ignored 实验目录并验证 SHA。

## 3. 物理与验收约束

- 2D 网格 10,241 节点、19,782 三角形、30,022 边；5,723 Si 节点，5,674 个自由 Si 节点。`th_lat` 是热边界，不作为电学端口。
- 外部 AverageBox 输运耦合 16,237 条（其中 3,548 零耦合），barycentric 体积，material-local Poisson 电荷、legacy node-local 接触重建。
- 300 K，Fermi–Dirac、OldSlotboom BGN、SRH/Auger。保持 predictor、IALMob、雪崩、量子修正和热耦合关闭。
- `constant_field` 表示低场常数迁移率加 HFS，**并非关闭高场迁移率**；真正关闭 HFS 的局部对照使用 `constant`。9 月 8 日白天的 Vg=4/Vg=8 完整扫描均为 `constant_field + edge_projection`。
- HFS 当前候选为 `constant_field + transport_cell_vector`，`jacobian_field_derivatives=true`；接触电场回退保留原配置。
- 原全局残差上限：Poisson `5e-8`，电子 `1e-11`，空穴 `3e-10`；局部 carrier-row `eps_row=1e-8`、违规数为零；端口 KCL 相对比值 `<=1e-8`。**不可放宽门限解决失败。** 零偏 KCL 仍用既有评分器的单独规则，不自行改变分母。
- 电压 V，几何 μm，密度 m⁻³，评分电流 A/μm。CSV 的 `current_total` 和 `current_total_A_per_um` 不是同一单位；用后者比较。
- 标量 `residual_norm` 可受行缩放影响，不可拿它代替 `block_psi/block_phin/block_phip`。失败 CSV 中的零电流占位不能评分。
- D5 对齐 31 个精确参考点，不用曲线插值冒充精确求解；单栅压和两栅压比都通过才关闭 D5。

参考输入路径可从任何本轮 `control.json` 复用，主要为：

```text
reference_staging/templates_ldmos_sentaurus2022/phase01_original_20260826_02/stage1_v4/vela_exact_topology/mesh.json
reference_staging/templates_ldmos_sentaurus2022/phase01_original_20260826_02/stage1_v4/vela_exact_topology/doping.csv
reference_staging/templates_ldmos_averagebox_full_mesh_20260831/transport_couples.csv
reference_tcad/templates_ldmos_sentaurus2022/contracts/materials.json
```

## 4. 已完成的曲线与电流误差定位

| 试验 | 完整性／结论 | 关键结果 |
| --- | --- | --- |
| 9/8 Vg=4 旧算法固定步长增长 | 0→40 V，31/31，通过 | 中位误差 2.130296%，P95 9.312277%，40 V 误差 2.031173% |
| 9/8 Vg=4 迭代感知增长 | 0→40 V，31/31，通过 | 未证明加速：固定／感知为 429／456 个接受推进步、4761／4837 次更新；子进程墙钟约 99.079／98.632 min |
| 9/8 Vg=8 旧算法独立重跑 | 0→40 V，31/31，精度失败 | 429 步、4985 次更新；子进程墙钟 142.819 min；中位误差 17.280895%，P95 18.643299%，40 V 16.206049% |
| 两栅压 D5 | 失败 | 40 V 电流比误差 13.892691%，超过 final 8% |
| vector HFS Vg=8、Vd=4 V | 局部通过 | Id=2.3062466308295835e-4 A/μm；SDevice 2.27155151227272e-4；误差 1.527375% |

旧完整结果证据：

- [Vg=4 A/B 报告](templates_ldmos_stage4_vg4_step_growth_ab_2026-09-08.md)，[评分](../../reference_staging/templates_ldmos_d5_vg4_step_growth_20260908/fixed/score/summary.json)。
- [Vg=8 独立重跑报告](templates_ldmos_stage4_vg8_rerun_2026-09-08.md)，[两栅压评分](../../reference_staging/templates_ldmos_d5_vg8_rerun_20260908/d5_score/summary.json)。
- [电流误差定位报告](templates_ldmos_stage4_current_error_localization_2026-09-08.md)。

误差定位中的重要边界：用 SDevice 内部状态在 Vela 重放漏端 SG 电流，4 V／40 V 误差仅 0.02178%／0.02472%；主要差异在内部状态及输运算子，不能优先归因于端口积分单位。切换 vector HFS 后，4 V 自洽结果改善，但旧 40 V 状态直接切换两轮各 160 次失败，因此没有合格的新 40 V 电流。

SDevice 外来状态中的 4,065 个非输运节点 QF 占位曾制造大残差；只可在诊断中规范化这些非活跃行，不能改动 Si 物理状态冒充收敛。固定状态残差探针也不是自洽解。详细证据均在上述报告及其 staging 链接中。

## 5. 已确认的 Jacobian 问题与局部对照

完整报告：[HFS 开／关数值对照](templates_ldmos_stage4_hfs_on_off_numerics_2026-09-08.md)。
完整证据目录：`reference_staging/templates_ldmos_hfs_on_off_20260908/`。

两模型在 Vg=8、Vd=4 V 各自建立合格初态，再独立推进同一电压步；总计十个控制、198 个 JVP 方向、两个最多 12 次更新的线性诊断前缀均已完成。

| 控制 | HFS 关 | vector HFS 开 |
| --- | ---: | ---: |
| 4→4.0025 V，QF 上限 0.0025 V | 5 次通过 | 79 次通过 |
| 4→4.1 V，QF 上限 0.1 V | 8+1=9 次通过 | 129+1=130 次通过 |
| 4→4.1 V，QF 上限 1 V | 8 次通过 | 145 次失败 |
| 4→4.1 V，40 次预算后同偏压重启 | 9 次通过 | 40+40+22=102 次通过 |
| 4→4.5 V，QF 上限 0.5 V | 46 次失败 | 231 次失败 |

提前重启减少 21.54% 更新、同模型电流及 KCL 一致；并发墙钟未下降，不能宣称端到端加速。确认实际 `max_iter` 与 `min_newton_max_iter` 都为 40，曾有覆盖疑虑但已排除。0.5 V 和 1 V 上限不能直接推广。

JVP：节点 2951 的 phin→电子块无下限归一化误差，HFS 关约 `1.74e-10`，vector HFS 开 **87.0266%**，同一 HFS 状态采用“冻结导数＋冻结迁移率残差”后约 `3.58e-10`。节点 2949、3432 对应 live 误差约 12.71%、29.73%，对步长缩小保持稳定。节点 5569 有额外差分幅度敏感性，不能与稳定大偏差混为一谈。

源码入口（行号可能随修复移动，以符号搜索为准）：

- [AssemblerUtils.h](../../include/vela/equation/AssemblerUtils.h)：`transportCellVectorEdgeGradientMagnitudes`（约 2966 行），相邻输运三角形梯度矢量面积加权后取模。
- [CoupledDDAssembler.cpp](../../src/equation/CoupledDDAssembler.cpp)：`edgeElectronTransportFlux`（约 4309 行）及空穴对应路径；vector 分支读取预计算 `electronVectorMobilityFields[e]`／hole 字段，端点 QF 扰动没有更新该字段。
- 同文件 `analyticQfMobilityFeedback`（约 4573 行）排除 vector 模式；连续性通量列没有补全二维梯度的邻域 QF 依赖。接触回退扩展的 psi 列与此不同。
- 同文件雪崩辅助路径有 `transportVectorMobilityField` 扰动重算，不能把它误认为当前关闭雪崩时的连续性 Jacobian 已完整。
- [test_newton_solver.cpp](../../tests/test_newton_solver.cpp) 已有带 vector 配置的 **avalanche** 差分测试；仍需增加直接针对连续性、非零二维 QF 梯度、第三顶点依赖的测试，避免测试被 Poisson／零电流尺度掩盖。

JVP 数据解释：原 CSV `*_relative_error` 分母含下限 1，微小扰动下可能只是绝对小量。分析脚本重算 `||Jv-FD||/max(||Jv||,||FD||)`，应看分块误差及绝对尺度。测试滞后导数时必须同时 `jacobian_field_derivatives=false` 和 `freeze_transport_mobility=true`。

共有尾段问题仍存在：HFS 关的 4→4.1 V 直接步电子块停在 `2.55e-10`，同偏压重建一次更新后 `7.64e-13`；HFS 开由 `1.67e-10` 到 `4.19e-12`。其根因尚未严格拆分为参考重定位、重打包或浮点累加。原始线性系统范数后向误差小，但部分微小行的分量后向误差大，不能宣称线性求解所有局部行完美。

## 6. 可复用的合格 checkpoint

以下路径均相对于工作树根。旧状态只能用于诊断／初始化，新二进制需重新核验同偏压残差和电流。

| 用途 | 路径 | SHA-256 |
| --- | --- | --- |
| vector HFS 原始合格 4 V | `reference_staging/templates_ldmos_current_error_20260908/vector_reclose_vd4_retry/state.csv` | `9d2b18cde163eb949ecfb53e1f5f86d71bfca644f0d37552362e7282688412f8` |
| HFS 开局部对照 4 V | `reference_staging/templates_ldmos_hfs_on_off_20260908/initialize_on/stage_0/state.csv` | `27d90ee9be95a323d431501d6a1e6dbc88d41479aba05a68e2612e44b9b3e8c1` |
| HFS 关局部对照 4 V | `reference_staging/templates_ldmos_hfs_on_off_20260908/initialize_off/stage_1/state.csv` | `4e83f6cd233929505ef65a1ca441ce0bc7386c94ebf70a31991a6336cfd33dda` |
| vector HFS 最后接受 10.533333 V | `reference_staging/templates_ldmos_vector_hfs_validation_20260908/vg8_resume4/fixed/child_00153_reclose/state.csv` | `1eb4bcb9e38a8cd7fe74b3ac9b7f117ff26f3396cf16679ebfa5de4addc91857` |

零漏压独立起点：Vg=8 用 `reference_staging/templates_ldmos_vg8_20260907/continuation_r2/initial/state.csv`；Vg=4 用 `reference_staging/templates_ldmos_shortstep_20260907/zero_closed_psi/state.csv`。按目标模型重新核验，不能只依赖历史文件名。

## 7. 9 月 9 日现场运行状态

旧队列目录：`reference_staging/templates_ldmos_vector_hfs_validation_20260908/`。

- `queue_progress.json` 留存 queue PID 6712，首任务 controller PID 14204；`child_00154_direct/status.json` 留存 solver PID 19720。本轮逐一 `Get-Process` 均未找到；按名称检查 python／vela_example_runner 也未找到。
- ledger 最后更新 `2026-09-08T22:39:15.195127+08:00`，接受电压 10.533333333333326 V，77 个接受推进步、0 rollback、5 个精确参考点，累计已记账 10,898 次 Newton 更新。
- 已完成子进程墙钟累计 6267.174 s（104.453 min），CPU 5763.219 s（96.054 min）。这是从 4 V 起的旧候选部分段，**不含未完成 child_00154 的完整成本，也不是 0→40 V 总耗时**。
- `child_00153_reclose` 正常完成；Id=`3.2393894245239208e-4` A/μm，psi/e/h 残差分别 `1.74199e-9`／`5.45521e-12`／`8.61261e-25`，局部违规 0、最大比值 `1.03997e-10`；其状态 SHA 已本轮重算。
- `child_00154_direct` 只有部分 `control.log`，curve/attempt/iteration CSV 等为零字节，没有合格输出，不能跳过并算成接受步。
- `vg8_full/`、`vg4_full/` 目录为空，说明此队列的两个独立全程任务尚未开始。
- `postprocess_status.json` 的 `waiting_for_queue` 同样不是活进程证据；`result_summary.json` 仍是 9/8 20:57 的早期 4.005875 V 快照，不能用它替代最新 ledger。
- 本轮不改写上述原文件。新聊天再读实时进程和时间；本交接快照也可能随用户操作过时。

脚本注意：`run_queue.py`、`run_candidate.py` 和局部 `run_controls.py` 对已存在 plan/output 有拒绝覆盖断言，**不是可直接重执行的续跑入口**。`run_candidate.py` 还要求 `--start` 属于精确参考点，10.533333 V 不满足，不能照抄命令硬续跑。修复后若从任意 checkpoint 续跑，需要在新目录增加明确的接续适配并保留父状态和二进制血缘，或从合格 4 V 精确点建立新实验；完整评分最终仍需独立 0→40 V 曲线。

## 8. 新聊天建议执行顺序

1. 核验本工作树／HEAD／未提交修改，阅读 `AGENTS.md`、[当前文档索引](../README.md)、[配置](../config_schema.md)、[架构](../architecture.md)，再读第 5 节报告。确认旧进程是否仍不存在，不凭旧 JSON 的 `running` 重复启动任务。
2. 在原工作树实施 vector HFS 连续性 Jacobian 修复：包括端点和相邻输运单元顶点的 phin/phip 依赖、对应稀疏模板；保持接触电场回退分支正确，避免迁移率反馈重复计入。保留原 edge projection 和 HFS-off 行为。
3. 增加有判别力的 Catch2 测试：非零且非沿边的二维梯度、第三顶点 QF 扰动、电子／空穴连续性分块、接触回退开／关、live 与匹配 frozen 差分。先使测试复现旧问题，再验证修复；运行相关 Newton／mobility 回归，按实际影响扩展 CTest。不要只运行雪崩差分测试便宣布通过。
4. 用新的 Release runner／新输出目录复跑本轮同初态 JVP 及 4→4.0025、4→4.1 V 控制。核验大导数误差消失，比较包括失败重试的更新数、原门限、各端口电流和 KCL。低偏压也需单独回归；高偏压局部改善不能证明全程稳定。
5. 导数修复验证后，再评估提前同偏压重收敛、QF 表示精度及大步长。不要把 40 次预算重启直接视为最终自适应策略；不要推广已失败的 0.5 V 大步或 1 V QF 上限。
6. 完成修复版 Vg=8、Vd=0→40 V 独立曲线，再跑 Vg=4 回归；最终分别评分 31 点、零偏 KCL、Ron、40 V 端点和两栅压电流比。若任一步失败，保存失败状态，报告最后合格电压，不计零值占位、不放宽门限。
7. 全程继续使用已验证的电势坐标变换策略（旧驱动在物理 26.666667 V 处整体换参考，偏移 28 V）；新驱动需要重新验证状态／端口等价，不能把坐标平移当作改变器件偏压。
8. 记录结论、可重复命令、编译器／后端／runner SHA、初态 SHA、全程子进程成本与失败成本。若用户随后要求提交，再整理本地提交；不要把未完成结果写成通过。

## 9. Sentaurus 和脚本使用注意

用户已明确授权向本机 `sentaurus` 虚拟机上传／下载文件和执行有关仿真，无需再次询问相同范围授权。当前优先修复本地 Jacobian，不需要为了交接启动远端任务。需要远端时先读当前 SSH 工作流，复用已有输入包和 manifest，避免覆盖参考证据。

Vg=8 既有 SDevice Extrapolate ON／OFF 自适应参考分别为 49／63 个接受电压步、216／656 次更新；相应 drain 段日志时间 97.23／213.34 s、进程总计 320.40／418.44 s。参见 [SDevice 对照报告](templates_ldmos_stage4_two_gate_validation_2026-09-07.md)。不要将全程序、漏极扫描段、CPU 和墙钟混用，也不能把其整条曲线均值直接当作 Vela 4→4.1 V 的同序列对照。

局部脚本 `run_controls.py` 通过 `run_ab → controls → investigate → 历史 run_controls` 导入助手。新脚本若直接 `import run_controls` 会有同名模块循环导入风险；使用不同模块名或明确路径加载，并核验实际导入路径。原入口作为脚本直接执行没有该问题，但输出目录已存在，应复制／适配到新实验而非覆盖重跑。

Windows 排查经验：

- `git` fsmonitor 有时提示 daemon terminated；用 `git -c core.fsmonitor=false ...`。
- `rg` 目录路径不要写 `docs/validation/templates_ldmos*`；用 `rg ... docs/validation -g 'templates_ldmos*.md'`。
- `Get-CimInstance Win32_Process` 先前受限；用 ledger 中 PID 配合 `Get-Process`，并检查进程路径／启动时间和日志增长，防 PID 复用。
- 新启动后台进程使用隐藏窗口，保留 stdout、PID、原始状态与可审计输出；不要依赖聊天工具 session ID 跨新聊天恢复。

## 10. 可直接粘贴到新聊天的接续指令

> 请进入 worktree `D:\code-repo\vela-tcad\.worktrees\templates-ldmos-phase-a`，阅读 `docs/validation/templates_ldmos_stage4_handoff_2026-09-09.md`，核验分支 `codex/templates-ldmos-phase-a`、HEAD `959c1c91513c62caeb700c2f71d8cef194879419` 和现有未提交修改。请按第 8 节继续：先修复已定位的二维梯度 HFS 连续性 Jacobian 迁移率反馈导数遗漏，完成数值单测和同初态 Release 局部 A/B，再推进 Vg=8 与 Vg=4 的完整 D5 曲线及两栅压比验收。保持原残差／局部／KCL 门限，predictor 和 IALMob 关闭，保留旧证据。交接时昨晚队列进程已不存在，最后合格点为 Vg=8、Vd=10.533333333333326 V；先核验实时状态，不要根据残留 running 标志直接重跑或覆盖旧目录。必要的本机 Sentaurus VM 文件上传、下载和仿真已获授权。
