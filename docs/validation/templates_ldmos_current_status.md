# Templates/LDMOS 当前验证状态

更新：2026-09-11。当前工作分支 `codex/templates-ldmos-phase-a`。
本页区分历史资格和最新版本验证，阶段报告保留原始结果。
每项曲线资格对应其报告中的冻结程序和输入清单，不随工作区重编译自动转移。

| 范围 | 状态 | 证据 |
|---|---|---|
| G3 Id–Vg | 本日专项 Release：原梯度、向量梯度+重定心、再修正 Auger 单位三组均 31 点、六门限通过 | [后续执行](templates_ldmos_followup_execution_2026-09-11.md) |
| D5 Vg=8 V，0–40 V | 单位修正配置从零压完整 31 点通过，1116 次更新，0 回退 | [本次执行](templates_ldmos_worker_and_ialmob_execution_2026-09-11.md) |
| D5 Vg=4 V及双栅压电流比 | 单位修正配置从零压完整 31 点通过，1139 次更新，双栅压原工程/最终门限通过 | [本次执行](templates_ldmos_worker_and_ialmob_execution_2026-09-11.md) |
| 正确 Auger 单位配置 | 已完成 62 点完整推进联合验收，晋级运行入口默认；历史系数另行显式回放 | [本次执行](templates_ldmos_worker_and_ialmob_execution_2026-09-11.md) |
| 常驻 DC 进程 | 显式 `--worker` 复用输入准备，保留逐请求种子与原门限；验证范围见报告 | [本次执行](templates_ldmos_worker_and_ialmob_execution_2026-09-11.md) |
| IALMob / D4 | 全局耦合及双栅压 0–40 V 各 31 点完成，原工程/最终联合门限通过；Vg4/Vg8 最大电流误差 1.895%/1.802% | [完整接入及 D4 验收](templates_ldmos_ialmob_global_coupling_2026-09-11.md) |
| 自热 | 后续独立阶段，当前 D4/D5 均为 300 K 等温 | 同上 |

## 两个配置合同的适用范围

`contracts/physics_contract.json` 保留 G3 历史资格使用的 `edge_projection`。
最新 D5 使用显式的
[`linked_d5_auger_units_config.json`](../../reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_auger_units_config.json)，
选择 `transport_cell_vector`、接触电场回退、按需 QF 重定心及内部密度恢复。
该独立配置不修改通用默认值。向量梯度在 G3 的资格来自本次独立 31 点复核。
数值块、逐行和 KCL 门限分别沿用各算例原值；G3 六门限并不等同于 D5 门限。

2026-09-11 参数审计新增发现：历史生成器把 SI Auger 系数直接写入按 cm⁶/s
读取的 `unit_scaling` 字段，历史曲线的显式系数比目标小 10¹² 倍。
生成器已修正并通过回归；冻结 `linked_d5_config.json` 保留旧值以复现历史对照，
不得把它视为 Auger 参数完全对齐的最终推荐配置。
另提供 [`linked_d5_auger_units_config.json`](../../reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_auger_units_config.json)
及对应 inputs 清单；该配置已通过从零压完整联合验收并成为运行入口默认。
H/N0 增强仅完成固定状态源项评估，尚未实现自洽项，资格范围见后续执行报告。

## 独立运行入口

[`run_templates_ldmos_linked_d5.py`](../../scripts/run_templates_ldmos_linked_d5.py)
不导入历史实验 Python 脚本。配置和哈希清单纳入版本管理；外部网格、掺杂、
AverageBox 系数、合格零压种子及 Sentaurus 参考曲线仍需单独保留。
[`linked_d5_auger_units_inputs.json`](../../reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_auger_units_inputs.json)
列出各文件相对工作区的路径及 SHA256；缺失或变化会拒绝运行。

从 worktree 根运行，例如：

```powershell
D:\msys64\ucrt64\bin\python.exe -X utf8 scripts/run_templates_ldmos_linked_d5.py --manifest reference_staging/templates_ldmos_low_bias_recovery_20260911/binary/manifest.json --output reference_staging/unique_vg4_first8 --gate 4 --points 8 --preflight
```

`--preflight` 只读检查，不创建输出；删除该选项执行。输出目录必须不存在。
默认 bundle 使用正确 Auger 单位。历史回放另加
`--bundle reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_inputs.json`，
不能把两组参数的运行时间混为同一基准。
`--worker` 需要包含 `--dc-worker` 的新版冻结 Release；历史 runner 不支持该选项。
该模式只复用输入准备，每个请求仍重新加载其显式状态；未启用时沿用逐请求子进程。
`--points 2` 用于迁移对照，`8` 为前八点，`31` 为全曲线。
31 点评分才包含 40 V 端点门限，短曲线不能认领完整 D5 通过。
默认可信预测最大步长 0.2 V、无预测最大 0.1 V；显式 `--max-step` 的变化
属于新性能实验。本次 0.4 V 前八点虽减少推进，却增加 46.26% 更新，已不予晋级。
高压参考系平移保留电势、密度、电流及 KCL 等价性检查。
`STOP` 文件仅在子步骤边界请求停止，不接受未通过的状态。

运行记录包含输入和程序哈希、每个求解请求配置/父状态/trace、原门限检查、
精确点状态哈希、全部 Newton 更新和子进程时间。Windows 使用 UCRT64
Release、Eigen SparseLU/COLAMD；此入口不将 Sentaurus 日志时间解释成
相同环境下的墙钟基准。冻结二进制清单需随外部证据保存。

## D4 独立配置入口

同一受版本管理运行器新增显式 `--physics-profile D4`，要求独立
`vela.templates_ldmos.linked_d4_inputs.v1` 清单。配置/几何准备入口为
[`prepare_templates_ldmos_ialmob_transport.py`](../../scripts/prepare_templates_ldmos_ialmob_transport.py)。
D5 默认及其合同保持不变。D4 零压初值必须重新闭合，不能继承 D5 解的物理资格。
当前验证进度和角点处理限制见[全局接入报告](templates_ldmos_ialmob_global_coupling_2026-09-11.md)。

D4 完整 62 点已通过原联合门限；仓库内的
[`linked_d4_inputs.json`](../../reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d4_inputs.json)
为其独立输入合同。清单默认参考系偏移分别为 Vg4：28 V、Vg8：4 V，
均须通过原等价性门限。`--resume-from` 支持经审核的只读前缀续跑到新目录；
必须保留源证据目录。具体失败尝试、修正和资格边界见完整报告。

## 后续任务顺序

1. **定位并减少 IALMob 结果记录开销。** 前八点 `dc.record_point` 占
   `dc_sweep.total` 的 40.79%；先分开测量端电流、连续性诊断及几何/状态准备，
   再评估单个已接受状态内的准备结果复用。现有计时不能证明全部记录开销
   都来自几何构造。以同一 Release、独占运行的前八点对照筛选，收益成立后
   重跑 D4 双栅压完整曲线，保持原物理和数值门限。
2. **定位高压逐行闭合的参考系敏感性。** 使用相同合格保存状态对比参考系
   平移后的块残差、逐行比和线搜索，查明导致迭代长尾的运算路径。
   Vg8 的 4 V 偏移只作为本算例已验证策略，不推广为全局默认。
3. **收敛剩余物理差异。** 分别评估 Auger H/N0 密度增强、高场驱动场及
   局部 IALMob 剩余偏差；先固定状态隔离因素，再实现和验证自洽影响。
   当前约 1.8% 的 D4 电流偏差尚未唯一归因。角点规则及其他晶向/几何
   的推广需要独立证据。涉及共享输运路径的修改需复核 G3/D5。
4. **另立自热阶段。** 在等温电学基线保留完整证据的前提下，明确热边界、
   材料热参数和原生配置，再制定相应验收；当前通过不包含自热资格。

常驻 worker 和更大步长继续保留为显式实验选项。已有对照中 worker 请求
墙钟增加 9.09%，0.4 V 步长增加 46.26% Newton 更新，均不作为已实现的
性能收益。下一阶段优先处理已定位的记录开销及数值闭合问题。
当前等温电学后续分析可使用现有资料；进入自热或扩展几何时再确认缺失输入。

## 本地提交前检查

最终工作区 Release 构建无待编译项；完整 Release CTest **767/767 通过，
99.49 s**，包含最终角点用例及运行器回归。18 份相关文档的 68 个本地链接
检查通过。日志保存在
`reference_staging/templates_ldmos_followup_20260911/ialmob_global_r4/commit_release_ctest.log`。
此次整理没有重跑长曲线；上述数值资格继续对应各报告中的冻结证据。
