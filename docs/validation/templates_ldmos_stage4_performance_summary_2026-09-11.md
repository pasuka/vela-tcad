# LDMOS Stage 4 性能优化汇总（2026-09-11）

本次整理包含求解器优化、数值回归测试和 13 份阶段报告。最新结果是 Vg=8 V、Vd=0–40 V 的 31 个精确参考点全部通过原门限；完整结果见[联动全曲线报告](templates_ldmos_stage4_linked_full_curve_2026-09-11.md)。本文件记录本次提交范围，早期报告保留当时的结论与证据。

## 纳入源码的改动

- 连续性诊断缓存：每次 Newton 求解拥有独立缓存，按状态逐位匹配，在逐行缩放与逐行收敛判据之间复用；重定心时清空，异常后不保留有效命中。
- 三项热点优化：复用 SRH/Auger 载流子状态；残差及连续性诊断复用既有边迁移率、材料与几何静态缓存；按模型需求构造雪崩及单元 Jacobian 散布表。量子耦合导致 SRH 与输运密度不同时仍分别计算。
- 按需 QF 重定心：新增默认关闭的 `quasi_fermi_recenter_on_stall`，在受保护条件下重新分配参考值与增量，保留低位并刷新缓存；每次求解最多重试 8 次，不放宽原判据。配置细节见[配置文档](../config_schema.md)。
- 恢复阶段计数修正：递归密度恢复后的 `NewtonResult.iters` 和 `solution.iters` 累计恢复前已接受的更新，与合并 trace 一致。
- 数值测试覆盖源项可加性和导数、状态变化下的独立通量对照、缓存上下文、参考坐标变换不变性，以及一层/两层恢复的计数一致性。

带保护预测器、0.1/0.2 V 联动步长、direct 内部密度恢复组合以及跨参考系预测保护，仍由本地实验适配器驱动；没有改动通用扫描器默认策略或物理模型默认值。

## 验证与结果

使用 Windows UCRT64 Release（`-O3 -DNDEBUG`，gprof 关闭），运行后端为 Eigen SparseLU/COLAMD。冻结 runner SHA256：`1e7464bb5df420cf02b37527e19098ede0b0300bc64100b0a9137837f79b1331`。

整理提交时，329 个当前源码/头文件/测试及构建文件与冻结清单逐字节一致，runner 哈希一致。既有验证记录为：

- Release CTest 串行 **747/747 通过**，0 失败、0 跳过，日志总时间 265.77 s。
- carrier-row recovery 专项 **5 个测试、28 条断言通过**。
- 同程序前 8 点配对对照及完整 31 点曲线全部通过原数值门限；全曲线步长回退 0 次。

本次提交整理仅新增文档入口与汇总，没有再次修改已验证求解器，也没有重复长时仿真。测试日志位于 `build-release/low-bias-ctest.log`、`build-release/low-bias-ctest.xml` 和 `build-release/low-bias-recovery-tests.log`。

| 完整曲线指标 | 历史边缓存版 | 本次联动方案 |
|---|---:|---:|
| 精确参考点 | 31/31 | 31/31 |
| Newton 更新 | 5209 | 1089 |
| 子进程数量 | 859 | 253 |
| 数值分解次数 | 6479 | 1245 |
| 子进程墙钟 / s | 6908.02 | 826.08 |
| 控制器总时间 / s | 7111.84 | 930.75 |

历史整体对照的子进程墙钟减少 88.04%；这包含程序与推进策略的累计变化，不能解释为单项优化的同日配对加速。最新曲线相对 Sentaurus 的电流误差中位数 1.330933%、P95 1.612869%、40 V 端点 1.099059%。两端时间统计与运行环境不同，具体比较口径见全曲线报告。

## 证据入口与复现范围

- 热点：[gprof](templates_ldmos_stage4_first8_gprof_2026-09-10.md)、[连续性缓存](templates_ldmos_stage4_first8_continuity_cache_2026-09-10.md)、[优化途径](templates_ldmos_stage4_hotspot_routes_2026-09-10.md)、[三项实施](templates_ldmos_stage4_hotspot_execution_2026-09-10.md)、[SparseLU 专项](templates_ldmos_stage4_sparselu_matrix_comparison_2026-09-10.md)。
- 更新差异：[Sentaurus 前 8 点时间](templates_ldmos_stage4_sentaurus_first8_timing_2026-09-10.md)、[Newton 工作量差异](templates_ldmos_stage4_newton_work_gap_2026-09-10.md)、[按需重定心](templates_ldmos_stage4_qf_recenter_validation_2026-09-10.md)。
- 推进与恢复：[预测器/0.2 V](templates_ldmos_stage4_guarded_predictor_step02_validation_2026-09-11.md)、[联动步长](templates_ldmos_stage4_predictor_step_link_validation_2026-09-11.md)、[低压恢复](templates_ldmos_stage4_low_bias_recovery_optimization_2026-09-11.md)、[前 8 点配对](templates_ldmos_stage4_linked_recovery_first8_comparison_2026-09-11.md)、[完整曲线](templates_ldmos_stage4_linked_full_curve_2026-09-11.md)。

原始仿真输出、冻结二进制、实验适配器和审计脚本保留于忽略目录 `reference_staging/`；构建与测试输出保留于 `build-release/`，均不纳入 Git 提交。全曲线入口为 `reference_staging/templates_ldmos_linked_full_20260911/run_full.py`，命令及证据清单见完整曲线报告。该适配器依赖本地历史辅助脚本、参考数据和合格种子；仅检出本提交不能直接重跑这些实验或显示报告中的本地图像，需同时保留相应本地证据目录。源码与单元测试按仓库常规构建流程使用。
