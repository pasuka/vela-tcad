# SimpleMOS PhuMob 低栅压与剩余源项验证

日期：2026-09-09。分支：`codex/simplemos-sdevice-validation`。前置结果见 [PhuMob 单元平均八点验证](simplemos_phumob_box_candidate_validation_2026-09-09.md)。

**已补齐 Vg=0/0.2 V 的八个低栅压点及两种初始化，并完成原生单元迁移率、导出精度和 SRH 积分账本。** 与已有 Vg=0.8/1 V 八点合并后，候选双初始化合格 15/16，legacy 合格 16/16。n23、Vd=1 V、Vg=0 V 的候选 Id 差仍为 **−0.80813524%**，绝对差约 **−3.00069e-18 A/μm**，不能由高栅压误差小而放行整个曲线。

SRH 体积差的电荷等价值约为 +3.00624e-18 A/μm，与上述差的绝对量级接近；同体积速率差仅约 −6.15312e-22 A/μm。这是优先级很高的源项候选，**尚不是经同源扰动校准的漏端电流因果归因**。本轮未修改生产 C++、模型默认、物理常数、SRH 体积或接受条件。

## 1. 范围与输入资格

同一 SimpleMOS Si/SiO2 网格：n19 为 1480 节点、n23 为 1482 节点；每点独立检查 907 个自由 Si 节点、1814 个载流子行。Windows UCRT64 Release，UMFPACK；沿用上轮已冻结的 runner 和静态库。Sentaurus 使用 T-2022.03-SP2，原生独立目录为 `/tmp/vela_simplemos_phumob_lowvg_20260909`，16 次 DC（8 个 PhuMob、8 个同偏置 Masetti 控制）已结束、取回并全部通过版本、偏置和 KCL 检查。16 份 TDR 场导出完成。

两端 plain PhuMob、OldSlotboom、独立匹配的基础 ni、300 K、Boltzmann、掺杂相关 SRH；Enormal/HFS 仍关闭。Vela 两组只切换 `legacy` 与显式 `element_box_phumob`；配套输运、Poisson 电荷及介电几何保持原组合。输入采用 `unit_scaling`；电流 A/μm、电势 V，场表密度 m⁻³。原生迁移率比较使用 cm²/(V·s)，SRH 速率使用 cm⁻³·s⁻¹；转换在脚本中显式执行。

沿用原两条初始化路径和最多一次同偏置重载。逐行 1e-6、KCL/Id 1e-8、端口一致性 1e-8、双初始化最大势差 1e-6 V、密度相对差 1e-4、Id 相对差 1e-6 均未改变。全局源的相对条件继续保留原 1e-10 下限；低于下限不算独立的相对源闭合通过。

## 2. 低栅压电流与初始化

误差定义为 `(Id_Vela/Id_Sentaurus−1)×100%`。表中 Vela 电流选自合格的 Vela 初始化；两条路径的最终资格单独列出，未通过的比较不能按正式双初始化合格点使用。

| NWell | Vg (V) | Vd (V) | 原生 Id (A/μm) | legacy 误差 | 候选误差 | legacy 双初始化 | 候选双初始化 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| n19 | 0.0 | 0.05 | 2.8061726898e-15 | +0.06096772% | -0.03145671% | True | True |
| n19 | 0.2 | 0.05 | 3.4462397566e-13 | +0.06770445% | +0.00380212% | True | True |
| n19 | 0.0 | 1.0 | 4.3707534017e-14 | +0.10428257% | +0.00466686% | True | True |
| n19 | 0.2 | 1.0 | 5.6404330146e-12 | +0.08720880% | +0.00531431% | True | True |
| n23 | 0.0 | 0.05 | 6.9636622735e-17 | +0.09905011% | -0.01819672% | True | True |
| n23 | 0.2 | 0.05 | 3.1976593339e-15 | +0.04334751% | +0.00358452% | True | True |
| n23 | 0.0 | 1.0 | 3.7131055650e-16 | -0.75009859% | -0.80813524% | True | True |
| n23 | 0.2 | 1.0 | 1.4119943360e-14 | +0.02847033% | -0.01632131% | True | False |

新增低栅压实际尝试 35 次，其中 4 次未通过，全部保留。首次或重载失败如下：

| 组 | NWell | Vg | Vd | 初始化 | attempt | 迭代 | 最大行比 | 失败原因 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| legacy | n19 | 0.0 | 1.0 | native | 0 | 6 | 4.233393e-04 | carrier_row_convergence_line_search_rejected |
| legacy | n19 | 0.2 | 1.0 | native | 0 | 200 | 4.801126e-06 | carrier_row_convergence |
| candidate | n23 | 0.2 | 1.0 | native | 0 | 29 | 2.842266e-05 | carrier_row_convergence_line_search_rejected |
| candidate | n23 | 0.2 | 1.0 | native | 1 | 200 | 5.450973e-06 | carrier_row_convergence |

合并 16 点的最大绝对相对误差：legacy 1.13093731%，候选 0.80813524%；15/16 点的误差幅度减小。资格仍以逐点标记为准。完整值与高低 NWell 配对见 [16 点账本](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/sixteen_points.csv)、[配对](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/nwell_pairs.csv)、[尝试记录](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/all_attempts.csv)。这些是四个 Vg 的控制点，未重算 51 点完整曲线。

## 3. Jacobian、少数载流子与场差

八个合格 Vela 初始化候选状态恢复 SRH 后，各用 ψ/φn/φp 方向和四个幅度检查 Jv。较小三幅度的 168 个正式分块全部通过，最大相对差 6.54302096e-07，门槛 1e-4；弱电子/空穴交叉源块仍单列，不用全残差差分代替独立弱源校准。该导数检查独立于另一初始化是否通过。

候选选取态的最大载流子行比为 5.45097322e-06。所有高低点按载流子分别保存最差五行、局部密度、是否为少数载流子、原始残差、绝对边通量及 SRH 项；源下限激活 0/128。见 [行排序](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/minority_row_ranking.csv)、[源下限与未截底比值](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/source_support.csv)、[Jv](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/post_summary.json)。

n23 候选在自由 Si 节点上的最大原生场差如下；密度列为无量纲最大相对差，不是百分数：

| Vg | Vd | max abs(Δψ) (V) | max abs(Δφn) (V) | max abs(Δφp) (V) | max 相对 Δn | max 相对 Δp |
| --- | --- | --- | --- | --- | --- | --- |
| 0.0 | 0.05 | 1.612715e-06 | 1.205674e-07 | 3.847242e-05 | 6.369351e-05 | 1.487822e-03 |
| 0.2 | 0.05 | 1.681911e-06 | 1.246028e-07 | 3.851312e-05 | 6.668962e-05 | 1.489394e-03 |
| 0.0 | 1.0 | 2.759832e-06 | 3.827840e-03 | 4.547826e-03 | 1.375441e-01 | 1.923662e-01 |
| 0.2 | 1.0 | 2.759854e-06 | 3.827681e-03 | 4.263030e-03 | 1.375388e-01 | 1.792989e-01 |

完整场最大值、RMS 和最差节点见 [场差账本](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/native_field_comparison.csv)。

## 4. G 截断与原生端口重放的边界

56000 个单元/载流子比较及 59360 个边系数比较完成。Masetti 单元公式与原生最大相对差 4.44e-16，八个 Masetti 原生电流控制均在原 1e-8 条件内重闭合。PhuMob 未截断单元最大差 2.22e-16；电子最大差约 1.38265e-9；空穴最大差 7.35772248e-06，2960 个含截断顶点的空穴单元仍未通过 1e-7。

使用上轮高 Vg 冻结的两个有效 G 下限作为诊断，新增低 Vg 单元最大差 6.66133815e-16。本轮拟合参数数为零。该结果验证已定位的差异形式能预测新状态，仍不能确定原生内部极小值搜索或截断算法，也未将反推参数写入生产。

低 Vg 的 16 个原生准费米势 SG 端口重放均未通过 1e-6，最大相对差 0.99819378。失败同时出现于 Masetti 与 PhuMob，不能将原生直接输出 Id 的资格与导出势场重放资格混用。

对每个端口先按节点合并有符号准费米势敏感度，再估计 binary64 最近舍入的半 ULP 影响。估计范围为参考电流的 0.0063424819～202.95247 倍，16 个重放偏差均在该范围内。正负一个 ULP 的直接重放与线性预测最大差 1.11881274e-14。这支持导出舍入足以解释量级；它不是未知原生高精度值的恢复或因果证明。上述失败不能调宽门槛后计为通过，固定态端口替换量也不能直接当作 Id 归因。

证据：[单元分组](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/native_groups.csv)、[校准汇总](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/calibration_summary.json)、[端口重放](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/native_ports.csv)、[导出精度](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/export_precision.csv)。

## 5. SRH 速率和体积的独立账本

使用原生 runtime 单元顶点 box 测度形成 Si 源体积；Vela 保留生产 all-cell SRH 体积。生产行积分通过同一边探针的物理粒子通量/内部通量比例换回物理单位，再除以当前体积取得 SRH 速率。比较只覆盖同一批自由 Si 节点，分解采用恒等式 `R_V V_V − R_S V_S = (R_V−R_S)V_S + R_V(V_V−V_S)`；没有改变源体积或重算候选。

下表后两列为积分源差乘 q 后的电荷等价值，**不是已校准的漏端 ΔId**：

| NWell | Vg | Vd | 积分源相对差 | 同体积速率差 (A/μm) | 体积差 (A/μm) |
| --- | --- | --- | --- | --- | --- |
| n19 | 0.0 | 0.05 | -4.9568423% | +9.481269e-23 | +9.859051e-19 |
| n19 | 0.2 | 0.05 | -4.8086651% | +8.508198e-23 | +9.803714e-19 |
| n19 | 0.0 | 1.0 | -0.2896708% | -1.678816e-22 | +3.130050e-19 |
| n19 | 0.2 | 1.0 | -0.3043575% | -1.788749e-22 | +3.343884e-19 |
| n23 | 0.0 | 0.05 | -0.0818165% | +8.279517e-23 | +1.382381e-20 |
| n23 | 0.2 | 0.05 | -0.0443078% | +8.554560e-23 | +7.652348e-21 |
| n23 | 0.0 | 1.0 | -1.5571476% | -6.153120e-22 | +3.006236e-18 |
| n23 | 0.2 | 1.0 | -1.5434777% | -6.262711e-22 | +2.994923e-18 |

n23、Vd=1、Vg=0 的体积贡献集中节点排序为 795, 792, 1114, 794, 791, 1115, 1203, 1188。源账本闭合到浮点舍入量级；生成源积分相对差约 −1.55715%，与漏端相对差并不相等。应先校准这些局部源到端口的响应，再考虑体积候选自洽修改。见 [逐节点源账本](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/srh_fields/srh_nodes.csv)、[源汇总](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/srh_fields/srh_summary.csv)。

初次 SRH 辅助审计的字段名大小写检查失败，原输入和脚本保留；修正版按实际 `srhRecombination` 字段名写入独立目录，未改变数据或单位。修正版开始前的输入路径错误也保留说明，见 [辅助审计失败记录](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/srh_failure_ledger.json)。

## 6. 验证与后续

新增脚本完成语法检查；既有 BGN 模型恢复 5 项和完整曲线驱动 4 项回归共 9 项通过。本轮没有修改 C++ 或重建二进制，未重复全量 CTest；前置阶段 772/788、16 项历史身份/缺失文件失败仍是上一轮记录，不能当成本轮新运行。输入、失败、原生归档、导出、数值账本和 QA 分别冻结，见 [复核汇总](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/review_summary.json)、[脚本 QA](../../reference_tcad/simplemos_sentaurus2022/phumob_lowvg_20260909/qa_summary.json)。未提交或推送。

下一阶段优先在已合格状态上，以 n23 高 Vd 体积贡献热点及低 NWell 控制做局部 SRH 源/体积的双幅度正负响应校准；原生端应直接记录守恒电流或足够精度的源响应，避免依赖已失去小梯度分辨率的导出准费米势重放。只有同扰动响应、收敛、双初始化及端口守恒通过后，才开展有限体积替换。G 截断算法资格、弱源资格及初始化失败轨迹继续保留。Enormal、高场饱和和完整 0–1 V 曲线尚未恢复或放行。
