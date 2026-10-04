# SimpleMOS M62 活性物种拓扑干预

## 结论

M62 分类为 `e3_stopped_no_identity_preserving_writer`。T-2022.03-SP2 的 `tdx` Tcl接口在手册和运行时清单中均能看到逐值写回能力，n23 TDR也同时含 `BActive`、`PActive`、`AsActive` 与 `NetActive`。但是无编辑 `TdrFileSave` round-trip 未通过冻结的字段恒等门禁。

网格、单元、接触、`doping.csv` 和掺杂元数据保持字节一致；132个导出字段中有 `36` 个无关机械字段改变。代表例 `Displacement_region0.csv` 的非零值尺度因子为 `10000`，最大对称相对差为 `0.9999`。因此不能证明局部活性物种写回不会同时改变其他状态。

按照冻结停止规则，没有创建变异TDR、没有执行SDevice E3求解，也没有只编辑派生 `NetActive`。这项结果关闭的是“当前工具链能否合法执行E3”，不是浮置单节点拓扑的因果效应本身。burst的决定性判别仍由M60提供；平滑NWell成分仍由M63独立归因。

机器报告：`reference_tcad/simplemos_sentaurus2022/active_species_topology_intervention/m62_active_species_topology_intervention_report.json`。
