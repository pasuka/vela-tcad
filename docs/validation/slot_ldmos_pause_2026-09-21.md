# SLOT LDMOS 用户暂停记录

2026-09-21 22:59 +08:00，用户明确要求暂停。已停止本任务的 Python
续算驱动及其子进程，确认本地投影 runner 不再运行。没有继续提交算例。

## 最后合格结果

- 新物理配置的无雪崩本地起始续算已到 10 V；保持既有数值门槛。
- 自洽雪崩界面驱动力投影为默认关闭的独立开关
  `ialmob_boundary_field_projection`，没有调整物理系数或验收阈值。
- 投影配置的 2/5/10 V 检查点全部合格，相对原始 Sentaurus 电流差
  分别为 −0.0384703%、−0.0377308%、+0.7831518%。10 V 原配置差为
  +4.3576893%；投影结果与原生 SG 控制差为 +0.1044350%。
- IALMob 43 用例、37457 断言通过；完整 Release 构建及 888/888
  CTest 通过。实际投影 10 V 雪崩 Jacobian 分块差分相对差 7.466e-9；
  极弱列精度限制见主报告，不把全块范数解释为每列的相对精度。
- 仅接触电场回退的 10 V 对照合格，基本不改变原配置电流。
- 同时投影和接触回退、从旧状态启动的第一次纠正未合格。
  从已合格投影状态启动的第二次纠正在本次暂停时被终止，不能算合格。

## 暂停位置及恢复要求

工作区：`D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next`。
证据根目录：`build/slot_energy_reference_20260921/`，以下简称 W。

- 合格投影 10 V：
  `W/av_boundary_projection/polish/local_corrector/validation/final_state.csv`。
- Release runner：`W/vela_av_projection_runner_20260921.exe`；运行时 PATH
  仍需 UCRT64 DLL 目录，线性后端 `VELA_LINEAR_SOLVER=umfpack`。
- `W/projection_av_10_to_20v/` 续算被暂停。10.25 V 第一候选未合格；
  高压链最后可恢复的已验收点仍为上述 10 V。恢复时不要把中断文件
  当作收敛检查点，也不要直接重跑会覆盖现有目录的脚本。
- `W/av_boundary_both_from_projection/` 是被终止的第二次组合控制。
- `W/native_av_high/` 的 15/20 V 原生 GradQF/SG decks 已上传到虚拟机
  `/home/tcad/sentaurus_runs/vela_oracle_2022/slot_energy_fields_20260921/`。
  `run_av_high.sh` 因 Windows CRLF 在 shell 解析阶段失败，**仿真没有启动**。
  暂停检查确认该脚本及两个 deck 无活动进程。恢复后需先将脚本转为
  LF、重新上传，再启动；本次暂停时没有执行这一步。
- 原生 Tcl R2 顶点数据有效；R3 元素电离系数请求被原生拒绝。
  元素密度 raw 值归一化尚未完全核实，不能按 cm^-3 直接使用。
- 原始 AvalDensGradQF 元素密度构造仍未严格还原；投影高压链属于
  SG 源项对照。高压、外电阻和 BVDS 验证仍未完成。

详见 [误差预算与验证报告](slot_ldmos_energy_reference_and_source_volume_2026-09-21.md)。
工作区保留所有代码和证据；未提交、未清理其他本地改动。
