# SDevice WP4 stage/deck emitter validation, 2026-10-02

工作树：`.worktrees/sdevice-cmd-frontend`，分支 `codex/sdevice-cmd-frontend`。
已先合并本地 `main` 的 `c37f9049`，合并提交为 `cf10b3f1`，无冲突。
WP4 已完成。用户已同意按当前 WP3 默认库修正栅接触验收基准，修正后的验证全部通过。

## 实现与接口

`scripts/vela_sdevice/deck_emit` 提供 `emit_stage_plan`、`rebase_stage`、`write_stage_plan`、
`StagePlanResult`、`EmitOptions` 和编号错误 `DeckEmitError`；详细调用和规则见
[workflow §2.7](../sentaurus_sdevice_workflow.md#27-stagedeck-发射器c4c10wp4)。
输入为 accepted C1、C2 与 WP3 报告，可选 C10 sidecar；输出 C4、各段 v2 deck 和映射决策。
生产接口只用 Python 标准库。已安装 `jsonschema` 时，测试额外用 draft 2020-12 校验器交叉验证。

七个 accepted 夹具对应 14 个 runner deck 和七个 stage_plan golden（21 个 JSON）：
PN2D 平衡/IV/BV、ABA Poisson、外接电阻、Schottky、BJT 多段扫描。
它们位于 `tests/fixtures/sdevice/<case>/<cmd stem>/expected_stage*.json`。
测试给缺 Parameter 的原始 cmd 副本补声明再解析，未改原始 cmd 与 IR golden。

发射覆盖 HDF5 stage 串接、Save/Load、电压步长与 Time 必经点、全局 plot 序号、
实测电压 rebase、v2 单位、WP3 寿命与 Auger 合成、接触、固定电荷、Math 决策与 sidecar 白名单。
手册 T-2022.03 Table 348 的步长默认值为 0.1/0.001/1/2/2。
C4 保留 Decrement 除数，runner 使用倒数；公开 schema 补齐 v2 SRH reference_doping_cm3，
并按 runner 接受范围允许 growth_factor=1。

## 授权的最小 C++ 修复

v2 deck 加载 C3 时，通用 deck 键转换曾把 C3 原生 `_cm3` 字段转成 `_m3`，随后被 C3 白名单拒绝。
用户已明确允许 WP4 加入最小修复和回归测试。`MaterialDatabase` 现在先识别带版本的 C3，
按材料自身 schema/单位读取；未带版本的材料覆盖文件继续沿用 deck 的键转换。
未更改物理模型、内置材料默认值或容差。

Catch2 用例 `C3 material schema keeps its units under a v2 deck` 验证 v1/v2 下相同的
物理密度、迁移率、Nc/Nv 与带隙，且错误的 `_m3` 键仍拒绝：12 个断言通过。

## 验证

Windows MSYS2 UCRT64 Debug，`windows-ucrt64-debug` configure/build preset。
配置检测启用了 HDF5/TDR、HighFive 状态文件与 UMFPACK；其他可用直接后端包括
SPQR、STRUMPACK、MUMPS、SuperLU。冒烟使用 runner 默认后端选择，未强制特定直接后端。

- Debug 全目标构建成功；修正 C++ 测试夹具 provenance 后，`test_scaling` 重建并通过。
- 已有 SDevice CTest：8/8 通过。
- WP4：22/22 项通过，包括修正后的栅接触基准与旧输入的预期偏移。
- 真 runner 冒烟通过四段：平衡、0→0.001 V、电流目标 1e-7 A/µm、Load 恢复平衡快照。
  每段退出 0，CSV/VTU/HDF5 均存在；电流终点与目标的绝对偏差小于 5e-10 A/µm。
  网格为 1 µm 正方形、4 节点、2 三角形、Si、均匀施主 1e15 cm⁻³、受主 0，300 K，
  单位按 v2；DD 段使用 Newton，Load 使用 frozen_state。材料由 WP3 当前默认 models.par 发射。
- 基准确认前，全 Debug CTest（排除 WP4 模块）：1017/1017 通过，572.55 s。
  命令：`ctest --preset windows-ucrt64-debug -E '^sdevice_deck_emit$' --output-on-failure -j 2`。
  最长单项 HFS 雅可比一致性测试 367.49 s，通过。
- 基准确认后，`ctest --preset windows-ucrt64-debug -R 'sdevice|C3 material schema|ascii_sources' --output-on-failure -j 2`：
  11/11 通过，60.62 s，包括真实 runner 冒烟。至此全部 1018 个已注册 CTest 均有通过记录；
  基准修正仅涉及测试与文档，没有重复前一轮已通过的无关数值检查。
- Python compileall 与 `git diff --check` 通过。

## 已确认的基准修正与剩余边界

计划 WP4 原验收写 `WorkFunction=4.6228018936 eV → work_function_eV≈0`。
该参考值来自遗留内置 Si；当前 WP3 默认库已经采用另一个参数集。
按已冻结公式 `χ+Eg/2+(kT/2)ln(Nc/Nv)`，当前库在 300 K 给出
`Φ_S,i=4.611003610676079 eV`。故旧输入应映射到 `0.011798282923921377 eV`。
2026-10-02 用户明确回复“我同意修改”。零偏移基准已改为当前库值，并保留旧输入预期偏移的断言；
Gold 的报告功函数也减去当前参考值。实现继续遵循冻结公式和 WP3 报告，未更改物理参数或容差。
WP4a 的历史内置材料验证与 WP13 的内置默认对齐保持原范围。

WP5 尚未开始：本次未交付 `vela_sdevice.pyz` 或完整 `vela device` 编排。
电流目标中显式 Time 中间点需要实测初始电流，目前 E407；现有控制器要求目标不低于
所选支路初始电流。一个 stage 的多个 Save 别名、非扫描接触上的串联电阻、
各材料互异但需要全局求解器参数的配置均 fail closed。
本次未访问或重跑 Sentaurus VM（WP14）。
