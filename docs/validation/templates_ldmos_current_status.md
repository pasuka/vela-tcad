# Templates/LDMOS 当前验证状态

更新：2026-09-12。当前工作分支 `codex/templates-ldmos-phase-a`。
本页区分历史资格和最新版本验证，阶段报告保留原始结果。
每项曲线资格对应其报告中的冻结程序和输入清单，不随工作区重编译自动转移。

| 范围 | 状态 | 证据 |
|---|---|---|
| G3 Id–Vg | 高精度Fermi/BGN、来源一致材料及Auger H/N0配置：31点、原六门限通过 | [本轮实现](templates_ldmos_auger_and_stage_chain_2026-09-12.md) |
| D5 Vg=8 V，0–40 V | 新H/N0配置从零压完整31点通过，1096次更新，0回退 | [本轮实现](templates_ldmos_auger_and_stage_chain_2026-09-12.md) |
| D5 Vg=4 V及双栅压电流比 | 新H/N0配置完整31点、1126次更新；双栅压原工程/最终门限通过 | [本轮实现](templates_ldmos_auger_and_stage_chain_2026-09-12.md) |
| 正确 Auger 单位配置 | 已完成 62 点完整推进联合验收，晋级运行入口默认；历史系数另行显式回放 | [本次执行](templates_ldmos_worker_and_ialmob_execution_2026-09-11.md) |
| 常驻 DC 进程 | 显式 `--worker` 复用输入准备，保留逐请求种子与原门限；验证范围见报告 | [本次执行](templates_ldmos_worker_and_ialmob_execution_2026-09-11.md) |
| IALMob / D4 | 新H/N0配置双栅压0–40 V各31点完成，原工程/最终联合门限通过；Vg4/Vg8最大电流误差1.8918%/1.8015%，1119/1073次更新 | [本轮实现](templates_ldmos_auger_and_stage_chain_2026-09-12.md) |
| IALMob 记录性能 | 共享几何及 referenced-QF 状态，前八点651.584→456.333 s，减少29.97%；双栅压完整62点复核通过，状态哈希全部不变 | [性能及差异定位](templates_ldmos_ialmob_performance_followup_2026-09-11.md) |
| 高精度 Fermi/BGN 及材料一致性 | 函数、解析导数、逆映射及材料审计已实现；Release原768项及新增回归通过，G3/D5/D4共155个精确点复核通过 | [统计精度实现](templates_ldmos_fermi_accuracy_2026-09-12.md) |
| Auger H/N0 | 自洽源项及导数已实现，G3/D4/D5共155点完整复核通过；独立新冻结资格，未覆盖旧证据 | [本轮执行](templates_ldmos_auger_and_stage_chain_2026-09-12.md) |
| D3/D2/D1 | 原生脚本和曲线重核后，依次通过原工程/最终等效曲线门限；限原Solve未启用量子方程且等温接触消融不可分辨的路径 | [本轮执行](templates_ldmos_auger_and_stage_chain_2026-09-12.md) |
| 自热 / D0 | 独立四方程冻结Release完成300 K双栅压62点回放、代表点及0–40 V完整62点；原电学及批准热学门限全部通过；最大电流误差1.7744%/1.5825%，峰温误差1.5051/2.1893 K | [电热实现与本轮证据](templates_ldmos_d0_electrothermal_2026-09-12.md) |

## 配置合同与冻结版本的适用范围

`contracts/physics_contract.json` 保留 G3 历史资格使用的 `edge_projection`。
此前完成资格的 D5 使用显式的
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
H/N0增强已新增显式自洽项及导数，默认关闭。独立
`linked_d4_auger_density_inputs.json`、`linked_d5_auger_density_inputs.json`
已使用新Release完成原门限复核；使用本轮冻结清单，不能自动继承旧程序资格。

2026-09-12高精度统计实现另提供
[`linked_d5_fermi_accurate_inputs.json`](../../reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_fermi_accurate_inputs.json)
及独立D4清单，使用来源审计的300 K材料和新模型重新资格的零压种子。
本轮新版G3/D5/D4完整资格见统计精度报告。运行时显式传入
该清单及报告中的新冻结程序，不能把旧默认bundle、旧种子或旧二进制的
资格自动转移到新模型。

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
为此前统计版本的独立输入合同。新版采用
[`linked_d4_fermi_accurate_inputs.json`](../../reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d4_fermi_accurate_inputs.json)，
亦完成62点原门限复核。清单默认参考系偏移分别为 Vg4：28 V、Vg8：4 V，
均须通过原等价性门限。`--resume-from` 支持经审核的只读前缀续跑到新目录；
必须保留源证据目录。具体失败尝试、修正和资格边界见完整报告。

## 本轮执行与后续顺序

1. **IALMob 记录优化已完成短曲线筛选。** 几何构造408→68次、状态构造
   772→567次，记录时间261.130→45.247 s。原348次Newton及八个状态哈希
   不变。完整双栅压已在独立 `_full_r2` 目录完成，62个状态哈希全部与
   原资格一致，原联合门限通过；并行曲线墙钟不作为独占性能比较。
2. **高压长尾的完整步拒绝原因已核实。** 同一26.666667→26.766667 V
   重放：16 V参考系32次更新中观察到19个完整步被电子块1e-11门限拒绝；
   4 V参考系8次更新，无这类拒绝。观察结果与历史状态哈希相同。
   继续保留本算例4 V策略，不放宽门限，不把浮点步长尾数的路径变化当收益。
3. **统计近似与BGN联动已实现，完整资格复核通过。** 原生D4状态的Vela漏电流
   回放误差仅0.0121%，但全域方程不合格。独立积分和逆函数对照已定位
   Fermi近似及其BGN修正的密度偏差。高精度函数、导数、逆映射及独立
   材料配置已完成，G3原六门限通过；D5/D4使用新模型重新资格的零压种子
   从零完整推进，原工程/最终联合门限均通过，共155个精确参考点。
   固定原生状态密度误差中位数降至约10 ppm，完整曲线电流差异仅略降，
   剩余局部输运/离散差异仍需定位，见[本轮实现](templates_ldmos_fermi_accuracy_2026-09-12.md)。
   H/N0已接入自洽项；新固定状态对照仍不支持其为电流差主因。
   电子残差平方和97.7%集中在距Si/氧化层界面10 nm内；强行将密度逆映射
   为准费米势反而放大残差，不能作为修正。角点及其他晶向的推广仍需独立证据。
4. **独立四方程 D0 完整联合验收已通过。** 用户指定先D3/D2/D1
   后原D0；批准0.1%热平衡及1 K/5%温升标准。独立热传导、热边界及固定原生
   热源回放已完成；Vg4/Vg8峰温误差0.0801/1.3822 K，场RMS误差0.1052/0.4048 K，
   均通过该限定范围热门限。后续已新增逐节点温度物性、与离散电流一致的
   默认晶格源及四方程Jacobian，使用独立实验入口验证；生产DC默认仍为等温。
   新冻结Release完成300 K双栅压62点回放，电流最大相对变化8.8263e-14；
   代表偏压及完整D0双栅压62点全部通过原电学和批准热学门限。两曲线各
   38次非零推进、0拒绝，Newton更新369/366次；全部78次求解重审通过。
   最终证据及适用范围见[电热报告](templates_ldmos_d0_electrothermal_2026-09-12.md)。

常驻 worker 和更大步长继续保留为显式实验选项。已有对照中 worker 请求
墙钟增加 9.09%，0.4 V 步长增加 46.26% Newton 更新，均不作为已实现的
性能收益。完整的新证据见[本轮报告](templates_ldmos_ialmob_performance_followup_2026-09-11.md)。
现有资料足以继续局部差异定位和电热耦合；自热范围及热门限已获确认。

## 已提交基线检查（5ea9b85）

最终工作区 Release 构建无待编译项；完整 Release CTest **767/767 通过，
99.49 s**，包含最终角点用例及运行器回归。18 份相关文档的 68 个本地链接
检查通过。日志保存在
`reference_staging/templates_ldmos_followup_20260911/ialmob_global_r4/commit_release_ctest.log`。
此次整理没有重跑长曲线；上述数值资格继续对应各报告中的冻结证据。

2026-09-12 未提交工作区另通过完整 Release CTest **768/768，183.92 s**。
该检查包含几何/状态复用及高压 profiling 观察，日志为
`ialmob_perf_r5/resumed_release_ctest.log`；不能与上面已提交基线混为一次检查。

本轮高精度统计实现另通过完整Release **768/768，101.17 s**，日志为
`templates_ldmos_fermi_accuracy_20260912/r1/release_ctest_r2.log`。
随后增加材料审计CTest注册，当时为769项；新增材料及原运行器定向2/2通过，
内部13个Python用例通过。此次测试注册变化没有替换曲线使用的冻结程序。

本轮Auger实现首先通过完整Release CTest 773/773。加入独立热传导后全量
780项中779项首次通过，唯一失败为新热合同误入冻结三文件合同目录。
将其移到独立`thermal/`目录后相关两项2/2通过，未放宽测试或物理门限。
本轮最终曲线及热学结果分别对应各自冻结程序，详见新执行报告。


## D0 温度相关迁移率后续（2026-09-12）

独立 IALMob 低场及 alpha=0 高场温度项、局部温度偏导已实现。
Vg4/Vg8、Vd40 V 原生热态对照覆盖全部10515个Si单元和5723个Si节点；
低场最大误差≤6.92 ppm，高场节点最大误差≤4.09 ppm。
该局部迁移率阶段当时未新增D0自洽曲线资格，也没有闭合D4剩余界面残差。
后续温度相关统计/复合、热源及电热耦合已在独立四方程入口实现并完成
完整D0资格，见下节；生产主电学入口仍默认300 K。
见[温度迁移率实现与证据](templates_ldmos_d0_mobility_temperature_2026-09-12.md)。

本轮最终完整 Release CTest 786/786 通过（253.81 s）。

## D0 四方程完整验收（2026-09-12）

新增逐节点温度物性及解析偏导、低/高场迁移率温度链式导数、与实际离散
电流一致的默认晶格热源及完整四方程Jacobian。独立入口使用UCRT64 Release、
Eigen SparseLU/COLAMD，原D0不启用Thermodynamic/Peltier/RecGenHeat。
两栅压0–40 V共62个精确点全部通过；最大温升场RMS误差为0.9552/1.5001 K。
完整Release CTest 799项中798项首次通过，唯一新增Python中文合同解码问题
改为显式UTF-8后重跑通过；后续输入/门限定向回归8个用例通过。
求解、冻结源码、原生62个温度场及最终评分证据见
[完整电热报告](templates_ldmos_d0_electrothermal_2026-09-12.md)。

后续可优先研究四方程Newton更新限幅成本，继续定位局部带边/界面输运差异，
并补原生Auger温度场逐节点对照。将独立电热入口整合到生产DC运行器及推广
到其他网格/边界时需要独立验证；本算例资格不会自动覆盖这些扩展。

## 本次提交前检查及日报素材

2026-09-12 最终工作区 Release 构建无待编译项；完整 CTest **799/799 通过，
777.26 s**，包含最终输入保护和四方程停止原因回归。日志为
`reference_staging/templates_ldmos_d0_electrothermal_20260912/commit_release_ctest.log`。
测试时间受宿主负载影响，不作为求解器性能对比。本次未重新运行长曲线，
数值资格继续对应报告中的冻结程序。

[工作日报素材及三组配图](templates_ldmos_daily_report_2026-09-12.md)提供完整
电流/峰值温升对比、自热电流下降和高偏压温度场，附 PNG、PDF 和复现入口。
