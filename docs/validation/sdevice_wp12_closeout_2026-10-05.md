# SDevice WP12 文档收口（2026-10-05）

本文保留 WP12 提交时的快照。后续日志覆盖修复、文档补漏及两平台部署包原生验收的
当前状态见[补缺记录](sdevice_native_package_e2e_2026-10-05.md)；以下缺口描述不作为当前任务队列。

WP12 文档任务完成：当前工作流、配置/导入入口、支持子集、消息和字段目录、可执行中性
示例及 sidecar、CHANGELOG、计划附录 A 和六条 DoD 账本已同步。
**阶段 1 的严格 DoD 尚未全部满足**：两个部署包的原生 PN2D 0V/IV/BV 全数值验收缺少
完整证据，已有包冒烟是合成网格。没有修改或放宽原 DoD/门限，也没有借文档收口修复已知物理差异。

工作树 `codex/sdevice-cmd-frontend`，本轮运行基线 `16fb066b`，提交即本报告所在的 WP12 提交。
本次仅文档及原创中性输入，没有改 Python/C++ 运行实现、schema、物理模型、材料默认或数值容差。
WP11 分发包保持其已验证的 `16fb066b` 身份；本次没有重新打包或把未复验包称为新交付。

## 文档与命令证据

- [工作流](../sentaurus_sdevice_workflow.md)：启动器/随包 Python、支持子集、逐段 Time/+1 电流符号、
  sidecar 命名与白名单、逐段 goal、VTU 字段键和 inline binary、具名资格与后续。
- [配置](../config_schema.md)与[导入](../sentaurus_import.md)：原生 `vela device` 与 JSON runner/旧 Python
  转换接口分开说明，C3 材料键不经过 v2 deck 转换，明确 `--check` 与 `--allow-lossy` 的边界。
- [中性输入](../../examples/sdevice/device_des.cmd)与[sidecar](../../examples/sdevice/device_des.vela.json)：
  无原生参考依赖；1 um 方形、4 节点/2 个直角三角形、施主 `1e15 cm^-3`、受主 0、300 K。
  sidecar 只将 stage 1 最大迭代次数从 50 改为 60，不改变物理、目标或收敛容差。
- [计划](../superpowers/plans/2026-09-26-sdevice-cmd-direct-run-development-plan.md) v2.7：WP12 卡片/汇总、
  A-1/A-2/A-4 已落实、A-3 最小 PN2D/BJT 子集已授权并完成，剩余 WP14 与 WP13 独立待办。
  附录 B 按 C++ 当前字段目录核对，去掉尚待接线的历史猜测。
- [CHANGELOG](../../CHANGELOG.md)与[文档索引](../README.md)补齐当前入口。
  链接审查还发现索引中一个已不在 checkout 的历史生成报告，改为明确缺失并链接已有版本化工作总结。

工作流 §2.3、§2.10 和导入新增入口的每条执行示例在 worktree 内实际执行；泛型 `<file>`/`FILE`、
库函数签名及旧转换流程的 `path/to/...` 属接口模板，不声称能按字面运行。
Windows 完成 `-h/-v/-P`、`--check`、2 线程 1 V 实跑、解析器 CLI、C1/C10 schema 检查、
材料/deck 发射 CLI、C++/Python importer help；Ubuntu 26.04 同一中性输入完成
`-P/-h/-v/--check`、2 线程 1 V 实跑与 importer help。
两个平台各段 `goal_reached=true`，输出 log/PLT/VTU/PVD，命令返回 0。
这些是可执行文档证据，不是原生 PN2D 数值比对。

输出核对还发现 `runtime_log.profile` 的接口缺口：C10/发射器接受并记录它，但编排器
随后覆盖该对象，runner 实际日志 profile 取决于 CLI verbose。示例改用实际生效的
`solver.max_iter`；该日志覆盖问题登记为已知限制，本轮未改运行实现。

门控核对共 11 例：7 个 accepted 夹具的中性网格副本及显式 Current-Time 用例退出 0，
3 个 rejected 原文退出 1，全部符合预期。文档审查核对 67 个登记编号、60 个代码字面编号、
32 个解析器模板和 43 个 C++ 字段条目；167 个本地文件链接与其中 2 个 Markdown 标题锚点
无断链。以上统计覆盖本轮修改文档及打包指南，不是全仓库文档或外部网站的审查。

## 六条 DoD 逐项账本

| 项 | 当前结论 | 提交 + 测试/报告证据 | 尚缺/边界 |
| --- | --- | --- | --- |
| 1. 两平台部署包原生 PN2D 0V/IV/BV 数值 E2E | **部分达成，不能宣布全部通过** | `0842e000`/`4734071b`，Windows 原生修复候选及四份契约 pass 报告；`sdevice_acceptance_tools`/`sdevice_e2e_pn2d_0v`；`16fb066b` 的两平台 `sdevice_package_smoke` | 已有 Windows 原生候选不等于两平台从解包入口运行的原生 0V/IV/BV；Ubuntu 包只有合成数值冒烟。BV 场量仅 −13.2 V 资格，−20 V hDensity 已知失败 |
| 2. ABA/外接电阻/电流目标/Schottky 至少解析+校验+deck | 门控层已达成 | `5ae3ce13`/`c9971ebd`，`sdevice_deck_emit` 的七个完整 golden、真实 Current/HDF5 链；本轮中性网格 `--check` 及显式 Current-Time 计划 | 原始 fixture/golden 不改；副本替换网格并补参数声明；对应原生数值资格未建立 |
| 3. fail-closed 注入、退出码、2022 热/IALMob 拒绝 | 已有通过证据 | `b93e477e`/`c9971ebd`/`5c723a2e`，`sdevice_cmd_parser`、`sdevice_orchestrator`（退出码类别 1–8）、`sdevice_materials`、`sdevice_acceptance_tools` | `--allow-lossy` 不绕过物理白名单；本轮三份 rejected 夹具仍退出 1；E201 当前为保留号，缺 cmd 走 E202 |
| 4. 已映射关键字/golden，无参考资产机器的前端测试 | 已有通过证据 | `b93e477e`/`5ae3ce13`/`16fb066b`，`sdevice_contracts`/`sdevice_cmd_parser`/`sdevice_deck_emit` 及其余前端 CTest | 实际没有单个 `sdevice_frontend` 测试；计划改为实际注册名，无参考 E2E 的 skip 不算 pass |
| 5. 并发、Unicode 路径、外部中止保留部分输出/释锁 | 已有通过证据 | `c9971ebd`/`1e3e39f2`/`16fb066b`，`sdevice_orchestrator` 的 signal/Unicode/partial 用例，`sdevice_package_smoke`、`sdevice_concurrency` | Ubuntu 26.04 完整矩阵零差异、线程 1/2；Windows 矩阵按原约定只锁冲突，其余明确跳过 |
| 6. 回归不退、PN2D BV 模板门、文档完稿 | 非条件参考回归及文档已达成 | 本轮基线 `16fb066b`：Release 1027 项选取，1023 pass/0 fail/4 skip；`pn2d_config_templates`、`pn2d_node_volume_policy_default_acceptance`、`pn2d_bv_m2_mixed_voronoi_self_consistent_control` 均通过；本次 WP12 文档提交 | 四个参考 E2E 本轮未配置参考目录，明确 skip；不声称此次重跑了原生全曲线或全反偏场量 |

具名数值范围、原生 TDR/PLT 与候选/契约指纹见
[WP10 收口](sdevice_wp10_closeout_2026-10-03.md)。用户指出 PN2D 膝区应在约 19.6 V；
−13.2 V 是场量样本。现有具名分段 log(I) 膝点为双方 −19.405 V，不等同于通用物理击穿定义。
−20 V 空穴密度诊断失败继续保留，不能据其他点 pass 推导全反偏空间资格。
全量 Release 236.03 s；原生参考四项 skip 与既有 WP10 已通过证据分别登记。

## 可审计文件与检查

生成产物保留在忽略目录 `build-release/wp12/`：

- `windows-full-ctest.log`：1027 项的实际结果，包括四个参考 E2E skip 和三项 PN2D BV 模板门。
- `example-check.log` / `example-run.log`、`linux-doc-commands.log`：两平台命令与原生 runner 进度。
- `example/device_des_vela/`、`materials/`、`emitted/`：manifest、schema 通过的 IR、C3、计划和 v2 decks。
  Linux 输出位于 `build-linux26/wp12/example/`。
- `cmd_checks.json` 与 `fixture_checks/*/check.log`：accepted fixture 副本与 rejected fixture 实际退出码；
  `current_goal/check.log`：显式 Current-Time `--check`。
- `doc_audit.json`：固定编号覆盖、解析器模板逐字核对、C++ 字段/单位核对、本地文件及标题锚点检查。
  临时审查脚本 `audit_wp12.py` 不加入 Git；外部 HTTP 链接不在本地链接检查结论内。
- `output_checks.json`：两平台逐段 goal、stage 1 迭代次数覆盖、VTU 元数据/单位与 PLT 电极列核对。

Windows MSYS2 UCRT64 Release：默认 Eigen SparseLU (COLAMD)，实际 HDF5/TDR/state 功能启用。
Ubuntu 26.04.1：GCC 15.2.0、HDF5 1.14.6、同样默认后端，SPQR/UMFPACK/METIS 可用。
本次复用已构建工具，没有为文档重建求解器或改变依赖。
