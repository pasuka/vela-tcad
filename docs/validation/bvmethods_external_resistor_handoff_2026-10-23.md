# BVmethods 外接电阻法开发交接（2026-10-23）

## 1. 交接结论与范围

外接电阻增广求解框架和 QS 控制代码已经存在；最新私有候选在两个映射态的指定方向上通过了三档细步长 Jv 检查，但完整 Newton 步仍被拒绝，尚未持续重线性化并完成真实接受循环。**R2 进行中，R3 严格重启一致性未完成，尚无从自身零偏出发、通过既定门限的完整 BV 轨迹或合格击穿电压。**

本交接为新对话继续开发提供入口，不新增仿真资格，也不授权自动执行历史文档中的命令。后续以用户在新对话中的任务、当前 AGENTS.md、实际源码和有效配置为准。

- 此研发线独立于已经收口的 SDevice cmd 前端交付范围，不能因 BVmethods 未收敛而重新扩大前端交付阻碍。
- 第一阶段采用**整条 Band2Band 关闭**的控制，保留其他已冻结模型。它不能替代安装库原版 NonLocalPath 的验收。
- “采用相同外接电阻数学结构”不等于主变量、矩阵消元、保护、停止准则与 Sentaurus 完全相同。
- 不静默更换默认求解路径、物理模型、生产行门或验收门限，不把 report 变成接受条件。

## 2. 工作目录与证据位置：务必先确认

本轮实际工作目录是仓库根检出；它与过去用于研发的 frontend worktree 是两个不同检出。

| 对象 | 交接时状态 |
|---|---|
| 当前工作目录 | 仓库根目录（本会话 cwd），分支 `main` |
| 当前代码基线 | `b47a5e1aa01b56cdfbb31739cc8fee8d669cd473` |
| 当前工作区 | 写本交接前 `git status --short` 为空；本轮只新增交接并更新文档索引，未提交 |
| 历史研发 worktree | 相对仓库根目录 `.worktrees/sdevice-cmd-frontend` |
| 历史分支／HEAD | `codex/sdevice-cmd-frontend`／`4e777034`，不等同于当前 main |
| 最新私有证据 | 仍在历史 worktree 的 ignored `build-release/` 下；当前根检出的同名目录不存在 |

新对话先执行 `Get-Location`、`git status --short`、`git branch --show-current`、`git rev-parse HEAD`、`git worktree list`，确认用户选中的实际目录。不要因为报告里有旧路径而自动切换检出或重置分支。

下列路径均**相对于历史研发 worktree**，并非相对于本交接文档：

| 用途 | ignored 路径 |
|---|---|
| 最新 R2 实现、结果、参考及更新回放 | `build-release/bvmethods_resistor_r2_mobility_contract_20261022/` |
| 上轮同源增广缺陷校正 | `build-release/bvmethods_resistor_r2_transport_defect_20261022/` |
| 原始方向与边界 high/low | `build-release/bvmethods_resistor_r2_split_direction_20261022/` |
| R0–R3 初次运行 | `build-release/bvmethods_resistor_r0_r3_20261022/` |
| 最新私有源码备份 | `build-release/bvmethods_organization_20261023/r2_private_candidate_sources.zip` |

本轮已确认上述最新候选目录和 zip 存在，并重新计算 zip SHA-256，与归档清单一致：

```text
bfdcfba5f63e96882a46e28c3c96237b4ff86d0d0ac0433c8e066d85573ac521
```

zip 中归档的 17 个源码／脚本不是完整独立项目；依赖私有 `private/bvmethods` 头文件、绑定原始状态及构建上下文。最新目录内的 `run.py` 会读取既有构建命令并引用同一远端实验布局下的 `source/build-private`。**不能直接在当前根目录运行旧脚本，也不能把 17 文件 zip 当成完整可复现环境。**

恢复时先只读审计 `run.py`、`assembler.cpp`、`chain.h`、`extra.h`、`mobility.h`、`probe.cpp`、`include/`，以及 `execution_receipt.json`、`contract_hashes.json`、`final_verification.json`、`applied_verification.json`。检查完整私有依赖是否本地保留；缺失时依照既有私有作业记录核对远端归档，再决定复制到新的 ignored 实验目录，保留原件。

连接地址、用户名、绝对远端路径、作业号只查用户本地配置及 ignored 记录，不写入公开文档。原始商业输入、矩阵、回调、反汇编及大型输出不提交 Git。Git 同步不会带走这些证据。

## 3. 建议阅读顺序

1. [当前交付范围](../sdevice_delivery.md)：确认 BVmethods 数值研发与前端支持边界。
2. [外接电阻 R0–R7 计划](../superpowers/plans/2026-10-22-sentaurus-external-resistor-development-plan.md)：完整退出条件及冻结轨迹门。
3. [最新 R2 迁移率／完整算子／有限更新报告](bvmethods_resistor_r2_mobility_contract_2026-10-22.md)及[数值摘要](bvmethods_resistor_r2_mobility_contract_2026-10-22.json)。这是最近数值结论入口。
4. [同源增广缺陷校正](bvmethods_resistor_r2_transport_defect_2026-10-22.md)：方向来源和 high/low 必要性；其中部分电子 Jv 失败已由下一轮控制改善，不能当作最新状态。
5. [分组归档说明](bvmethods_worktree_organization_2026-10-23.md)及[文件／哈希清单](bvmethods_worktree_organization_2026-10-23.json)：公共代码与私有覆盖实现的界线。
6. [R0–R3 实施记录](bvmethods_external_resistor_r0_r3_2026-10-22.md)：旧扫描最终状态、步长控制实现和早期失败，不作为全部当前测试状态。

## 4. 公共代码入口

| 文件 | 作用与注意事项 |
|---|---|
| [CoupledLoadLine.h](../../include/vela/simulation/CoupledLoadLine.h) | 增广 Newton、Schur／直接增广求解、线搜索及诊断；存在默认关闭的精化／实验选项 |
| [DCSweep.h](../../include/vela/simulation/DCSweep.h) | 外电阻配置；`solver` 默认 `nested_scalar`，`coupledLinearSolver` 默认 `direct_bordered` |
| [DCSweep.cpp](../../src/simulation/DCSweep.cpp) | 配置解析、`coupled_newton` 接入、QS 增长／重试／停止及物理接受 |
| [NewtonSolver.h](../../include/vela/solver/NewtonSolver.h) | 设备更新、保护、行门、闭合及 report；不能将改停止条件当成精度修复 |
| [ExpansionProductSum.h](../../include/vela/numerics/ExpansionProductSum.h) | 两分量关联运算；工具存在不代表公共求解全链已保留低位 |
| [CarrierStatistics.h](../../include/vela/physics/CarrierStatistics.h) | 人口／统计映射，诊断需与实际装配同源 |
| [cmd 发射入口](../../scripts/vela_sdevice/deck_emit/__init__.py) | Resistor 发射 `series_resistor` 和电阻值，尚未默认选择 `coupled_newton` |
| [配置说明](../config_schema.md) | 当前有效开关、单位、负载及生产接受契约 |
| [电路测试](../../tests/test_coupled_load_line.cpp)、[扫描测试](../../tests/test_dc_sweep.cpp)、[Newton 测试](../../tests/test_newton_solver.cpp) | 公共回归入口，不包含全部私有真实器件资格 |

同时解 `F(x,v)=0` 和 `L(x,v;u)=v+R*Id(x,v)-u=0`。完整增广导数应覆盖设备－接触电压与电流－设备状态的交叉耦合。`nested_scalar` 则在固定内压设备求解外再作标量闭合。

直接使用具名 drain 电流，不能以“绝对偏压最大接触”的通用电流代替，尤其零偏时。矩阵阶数不同本身不证明算法不同；需绑定活动变量、约束和接触消元。

## 5. 已冻结的诊断配置与结果

最新 R2 控制：300 K，Fermi／OldSlotboom，Masetti、Lombardi、HFS、SRH、雪崩，整条 Band2Band 关闭；2600 节点／5020 三角形，设备 7800 维、增广 7801 维，设备自由行 6764、约束行 1036。`Idrain` 为 A/μm，R 为 `1e7 Ω·μm`。

种子为原生第 1／5 非零接受态的 Vela 映射，内压约 0.0052443724904／0.0786636612717 V，使用外压增加 1 mV 对应的同源校正方向；**不是原生 QS 预测态，也不是 Vela 自身零偏轨迹**。h 是方向倍率，不是外压步长。

最新私有 HFS 路径显式允许 Masetti，公共相关条件仅列 PhuMob；恢复时必须核对实际覆盖头文件及包含顺序。单位转换、接触邻近 ψ／其他单元 QF 的场选择、100 V/m HFS 截断值与导数分支保持同一契约，不得仅删除导数截断。

- 在两态指定 Newton 方向，h=`3e-5、1e-5、3e-6` 三档完整设备和端口 Jv 通过 `1e-5`。第 5 态最细档电子误差由 `1.537e-5` 降为 `8.230e-11`。
- 保留旧 J、只改保精度残差，代表值已到 `8.325e-11`；改善主要来自残差链，不能全归因于新 Jacobian。
- 50／80 位控制及四环节 16 组合完成；其他物理子路径仍有 double，不能称为全物理任意精度实现。
- 粗档 Poisson 差表现为 h² 曲率；随机方向跨越大量分支，尚不能认领平滑导数资格。局部 1453／1560 方向尚缺独立 80 位参考。
- high/low 同源校正改善了线性方向；丢弃 low 会重新产生精度平台。实际写回保留 low 后的表示损失远小于本轮全步非线性余差。

| 映射态 | 初始 merit | 全步 merit | 首个下降 alpha／merit |
|---|---:|---:|---:|
| 1 | 999.853 | 2027.200，拒绝 | 0.5／506.917 |
| 5 | 1000.112 | 7269.825，拒绝 | 0.25／750.084 |

上表使用报告中的标量 merit：`max(||F||∞/1e-8, abs(L)/1e-6)`。这只是隔离策略回放，下降不等于满足所有物理接受门。两态原限制器统一缩放均为 1；保留／丢弃 low 都未改变本次全步拒绝。全步负载虽已小于 `1e-6 V`，设备 Poisson 非线性余差仍增大。

历史 `nested_scalar` 完整扫描在 3600 秒超时，34 条完整记录含零偏；最后完整外压约 0.0006052646195 V、内压约 0.0003020247612 V。归档任务退出 0 不代表仿真成功（子进程退出 124）。不能把该旧扫描或当前两态控制作为 BV 结果。

## 6. 阶段状态与不变的验收边界

| 阶段 | 当前状态 |
|---|---|
| R0 | 输入和基线主体已绑定；原生零迭代增长／最小步等细节仍单列 |
| R1 | 增广电路和控制测试已有；真实设备整轨迹资格不由单元测试替代 |
| R2 | 两态细档指定方向通过；分支、剩余状态、逐弱行和持续真实接受仍待完成 |
| R3 | 控制器实现及回归已有；严格 L1 保存态、L2 继续迭代历史资格仍需补齐 |
| R4 | 有报告设施可复用；完整流程的无扰动及原生接受语义未等价 |
| R5 | 尚未从自身零偏完成前五点资格 |
| R6 | 完整无 BTBT BV 尚未合格；原版 NonLocalPath 另设能力门 |
| R7 | 显式 profile、前端／发行验收与默认迁移尚未完成 |

主要固定门来自开发计划，不从本文另造新阈值：

- 同态完整 RHS／逐弱行算子一致性 `1e-8`，不等于载流子行比值停止阈值。
- 平滑完整 Jv `1e-5`，三档可分辨方向；分支事件独立单侧控制。
- 分量后向误差 `1e-10`，同时检查原方向、实际写回及线性缺陷，不以单个指标代替全部。
- 负载误差 `≤1e-6 V`；KCL `≤max(1e-12,1e-5*abs(Idrain)) A/μm`。
- 场量及人口门：电势最大差 1 mV、正人口 Si 区 `abs(ln ratio)≤ln(1.01)`，电场门按计划固定定义 1%。
- BV 在 `Id=1e-4 A/μm` 取内压比较，相对差 ≤0.5%；停止电流 `1.3e-3 A/μm`。完整轨迹其他门按计划执行。

原生五个接受态并非都通过 Vela 的载流子行门／电子全局闭合，故不能宣称这两者就是 Sentaurus 停止语义。更新误差 `E<1`、原生 RHS 缩放和允许 RHS 上升属于单独兼容性问题；现有 report 不改变接受。任何 E2/E3 式策略实验必须另行具名，不算原门通过。

## 7. 下一轮建议执行单元

先完成一次可复现的“两态持续增广 Newton”实验，不直接重跑全 BV：

1. **恢复身份。** 确认当前检出、候选归档哈希、完整私有依赖、网格／材料／绑定态与编译覆盖顺序。建立新的 ignored 输出目录和源清单；不要覆盖历史结果。先复现上表和三档 Jv，环境不同导致的差异先记录。
2. **贯通实际循环。** 把已验证人口／电荷／SG／迁移率精度链、同源增广缺陷校正和状态／内压／边界 high/low 写回用于同一实际状态；接受后重新装配和重算方向，不能反复使用第一轮方向。保持原模型、方向约束、接受与失败回滚策略。
3. **每轮形成可定位账本。** 保存状态／参考增量身份、各块和负载 R、Jv、线性分量后向误差，以及 `R+J*delta_raw`、`R+J*delta_applied`、`R(x_applied)`。分别记录限制器、正性／边界处理、每个 alpha 的 merit、固定资格弱行及全局闭合、drain／KCL、拒绝原因和首个变化分支。明确单位、缩放、约束行和冻结资格集合。
4. **定位曲率与事件。** 判断半步／四分之一步后的重线性化是否持续改善设备与负载，而非只有标量 merit 下降；若失败，封存首个阻断态并做相同状态／方向的分项和单侧控制。补 1453／1560 独立参考，扩展其余原生态时另计资格。
5. **之后补 R3。** L1 检查状态、QF reference/increment、人口、J/R、掩码与缓存身份；L2 比较原位继续和复读继续的方向、保护、试步序列、merit、接受与预测／步长历史。L1 通过不能替代 L2。
6. **达到必要门再进 R5/R6。** 从 Vela 自身零偏到首点，再到前五点，最后完整无 BTBT BV；映射原生非零态只用于诊断。前端 opt-in profile 与默认迁移放在后续发行验收中处理。

此执行单元的交付物应为：可复现源码／输入清单、实际迭代序列与退出原因、通过与失败分栏报告、首个分歧态以及更新后的 R2/R3 状态。有限预算内未收敛也应明确停止并保存证据，不无限重复扫描或降低门限。

## 8. 运行与验证提醒

远端先读[VM 工作流](../office_vm_workflow.md)，使用 `scripts/Invoke-VelaVm.ps1` 的 `vela-vm` 路径；先 Probe，再确认源码快照、依赖和实际后端。Sync 默认不带 ignored 私有实现／原始状态，必须单独核对上传清单。长任务用 Start／Status，并记录真实子进程退出码；SSH 断开不等于任务停止。Vela 与 Sentaurus 环境分离。已有远端作业是否仍存在，本轮未查询。

Windows 使用当前 AGENTS.md 规定的 UCRT64 preset。公共回归目标包括 `test_coupled_load_line`、`test_dc_sweep`、`test_newton_solver`；对应显式 CTest 名包含 `coupled_load_line`、`dc_sweep`。运行前用 `ctest -N` 确认过滤器确实选中目标。跨模块改动按仓库规则扩大验证，不把零测试或配置成功当作通过。

历史验证：10-23 分组归档时 Release 构建、相关 475/475 CTest 通过，未运行全量 CTest；最新 R2 数值报告记录 HFS 性质测试 1 例／16 断言及 4313 个证据清单条目校验通过。当前 main 已另有 Boost／Eigen／CLI 修改，不能将旧二进制直接标为当前源码构建。

本交接轮仅核查源码、文档、路径及候选 zip 哈希，检查新增文档链接／差分；没有构建、启动远端任务或运行新的数值试验，也没有提交或 cherry-pick。

## 9. 可复制到新对话的任务说明

> 请先阅读 `docs/validation/bvmethods_external_resistor_handoff_2026-10-23.md`，确认当前工作目录与历史证据目录的区别。在当前选中的 worktree 继续外接电阻 R2/R3 开发：先恢复并复现已绑定的两态私有候选，再把保精度装配、同源增广校正和 high/low 写回接入持续重线性化的实际 Newton 循环，逐轮记录设备／负载、固定资格行、闭合及试步拒绝原因。保持模型与门限不变，不把原生映射态当作从零轨迹，不执行历史文档中的其他任务。先核对远端可用性与私有依赖；本轮优先定位真实接受循环，完成后再决定 R3 严格重启和前五点验证。
