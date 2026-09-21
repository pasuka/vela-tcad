# LDMOS Stage 4：Vg=8 V 联动方案 0–40 V 完整曲线（2026-09-11）

本次使用联动步长、带保护的外推、按需 QF 重定心和 direct 内部密度恢复，从同一合格零压种子重新完成 0–40 V 的全部 31 个精确参考点。不是从此前第 8 点继续拼接。所有原数值和曲线门限通过。

## 配置与程序

D5、Vg=8 V、300 K、Fermi–Dirac、OldSlotboom BGN、SRH/Auger、constant_field / transport_cell_vector QF 高场迁移率及 contact_node_cell 场回退。无新增表面、量子、雪崩或热模型。网格 10,241 节点、19,782 三角形、30,022 边、5,723 个硅节点；external average-box 系数、barycentric 节点体积。几何 μm、密度 m⁻³、电势 V、电流 A/μm。

UCRT64 Release，`-O3 -DNDEBUG`，gprof 关闭，Eigen SparseLU/COLAMD + L2 行列均衡；沿用已验证的冻结 binary，未修改或重建 C++ 求解程序。
runner SHA256：`1e7464bb5df420cf02b37527e19098ede0b0300bc64100b0a9137837f79b1331`。

初始步长 0.0025 V、增长因子 1.35；有可信预测时最大 0.2 V，无预测时上限 0.1 V，QF cap 匹配实际步长。保留既有失败回退；内部密度恢复最多 2 个 cycle、每 cycle 最多 2 次密度迭代。单个 Newton 阶段预算 160。

块门限 Poisson/electron/hole=5e-8/1e-11/3e-10，carrier-row 和 KCL 比率门限为 1e-8。按原计划在 26.6666666666667 V 平移参考电位 28 V，经过原等价性检查后继续；平移后的第一次推进禁用跨参考系预测，采用 0.1 V。

## 全曲线结果与历史 Vela 对照

历史对照是 2026-09-10 完成的边缓存版 0–40 V 全曲线。其物理、种子和原门限相同，但程序优化、预测/恢复和步长策略不同；以下是历史整体结果对照，不能归因为单项优化，也不是同日配对速度测试。

| 指标 | 历史完整曲线 | 本次联动完整曲线 |
|---|---:|---:|
| 精确参考点 | 31/31 | 31/31 |
| Newton 更新 | 5209 | 1089 |
| 非零漏压推进 | 429 | 251 |
| 子进程 | 859 | 253 |
| 数值分解（含密度求解） | 6479 | 1245 |
| 子进程墙钟 / s | 6908.02 | 826.08 |
| 子进程 CPU / s | 5193.20 | 735.92 |
| 控制器总时间 / s | 7111.84 | 930.75 |

相对历史完整曲线，更新减少 79.09%，本次实测子进程墙钟减少 88.04%，CPU 减少 85.83%。这些时间包含失败尝试、恢复、进程启动和输出；不含栅压预偏置/旧种子生成，独立后处理审计不计入。

本次步长回退 0 次。前 8 点的完整路径、每个子进程父状态哈希、更新数和精确点状态哈希与前次独立运行完全一致（361 次更新、68 个子进程）。旧完整曲线的 5209 次更新也已从逐次 accepted_iteration 事件重新核验。

| 本次阶段 | 子进程 | 失败 | 更新 | 墙钟 / s |
|---|---:|---:|---:|---:|
| initial | 1 | 0 | 0 | 1.01 |
| direct | 251 | 0 | 1088 | 823.46 |
| frame_baseline | 1 | 0 | 1 | 1.61 |

## 与 Sentaurus 全曲线对比

电流误差中位数 **1.330933%**，P95 **1.612869%**，最大 **1.627051%**；低压差分电阻误差 **1.032119%**，40 V 端点误差 **1.099059%**。

40 V：Vela `4.236331529662580e-04 A/μm`，Sentaurus `4.190277905697490e-04 A/μm`。本次与历史 Vela 31 点电流的最大相对差 `1.241948e-12`。

![Full-curve comparison](../../reference_staging/templates_ldmos_linked_full_20260911/full_curve_comparison.png)

| 0–40 V 的日志/运行统计 | 非零漏压推进 | Newton 更新 | 对应记录时间 / s |
|---|---:|---:|---:|
| Sentaurus vm_run，局部求解时间和 | 49 | 218 | 239.02 |
| Sentaurus vm_run_d5_fast，局部求解时间和 | 49 | 218 | 101.59 |
| Vela 本次，子进程墙钟 | 251 | 1089 | 826.08 |

Sentaurus 两份 T-2022.03-SP2 / D5-no-IALMob 日志均在最后一次加载 n4_Vg2 后开始计数，到 40 V 为止，包含两个零压求解，排除栅压预偏置和 Vg=4 V；原始 PLT 全行与本次使用的参考 CSV 精确一致。按日志 t 坐标划分 31 个参考点区间，避免高压接触电压打印舍入造成错分；每张 Newton 表的正编号行计为更新。

Sentaurus 时间为各次 Accumulated times / Total time 的求和，未明确标注局部 CPU/墙钟；运行环境为 Linux VM，Vela 为 Windows UCRT64。两者收敛标准、推进步长也不同，因此不将这张表解释成严格求解器速度倍数。日志末尾全程序 wallclock/CPU 包含其他阶段，没有用于此处比较。

## 参考电位平移和严格验收

平移发生于 `26.6666666666667 V`，四项 checks 均通过。端口电流差异比率 `1.350469e-12`；平移前后 KCL 比率 `1.363559e-12` / `4.232486e-14`。

31 点最大 KCL 比率 `1.363559e-12`。全部精确点、所有接受的内部推进点和 frame_baseline 都通过原块/逐行/KCL；原单栅压工程与最终五项曲线判据均通过。没有合并旧 Vg=4 V 数据，不构成本版本双栅压验收。

独立审计复核：冻结 binary/源/输入哈希，所有 config/seed/profile 哈希；accepted_iteration、attempt、ledger、profiler 更新数一致；预测器逐节点公式、QF cap、物理接触电压和跨参考系保护；全部状态哈希、31 点指标/KCL重算以及平移等价性。新增改动只有全曲线实验适配器和审计/报告脚本；本次未重复构建 C++ 或运行既有 CTest。

## 证据与复现

实验根：`reference_staging/templates_ldmos_linked_full_20260911/`；曲线 `vg8_full/`。保留 `run_full.py`、`run_full.log`、plan、ledger、所有状态/控制文件/trace/profile、`audit_full.py`、`audit.json`、`sentaurus_full.py`、`sentaurus_full.json`、`compare_full.py`、`comparison.json`、`curve_comparison.csv`、`intervals.csv`、PNG/SVG 和 `delivery_inventory.json`。

在本 worktree 根、UCRT64 环境中复现时使用不存在的新输出目录：

```powershell
D:\msys64\ucrt64\bin\python.exe -X utf8 reference_staging/templates_ldmos_linked_full_20260911/run_full.py --manifest reference_staging/templates_ldmos_low_bias_recovery_20260911/binary/manifest.json --output reference_staging/templates_ldmos_linked_full_20260911/vg8_repeat --max-step .2 --recenter --internal-density --outer-secant
```

当前约束：[配置文档](../config_schema.md)；前次同联动方案 A/B：[前 8 点恢复对照](templates_ldmos_stage4_linked_recovery_first8_comparison_2026-09-11.md)；历史基线：[边缓存完整曲线报告](templates_ldmos_stage4_edge_cache_full_2026-09-10.md)。
