# LDMOS IALMob 局部对照闭合 — 2026-09-11

## 本轮结果

在固定 D4、Vg=8 V/Vd=40 V 状态上，已解释前轮的百分比级差异，
将对照扩大到全部 **10515 个 Si 单元、5723 个 Si 节点**。
未调整 IALMob 参数、原曲线验收门限或生产求解器配置。

| 对照范围 | 电子绝对相对误差 P95 / 最大值 | 空穴绝对相对误差 P95 / 最大值 |
| --- | ---: | ---: |
| 原生单元 → 原生节点，低场输出 | 1.14e-14 / 1.00e-12 | 1.64e-14 / 1.18e-12 |
| 原生单元 → 原生节点，高场输出 | 7.99e-15 / 1.09e-12 | 1.02e-14 / 1.18e-12 |
| 全部单元，独立 IALMob 低场重建 | 1.21e-9 / 3.93e-6 | 2.65e-6 / 5.07e-6 |
| 全部节点，独立低场重建后回映射 | 2.51e-9 / 3.88e-6 | 2.66e-6 / 4.95e-6 |
| 全部单元，加入原生默认 HFS 策略重建 | 1.44e-9 / 1.02e-5 | 2.59e-6 / 3.97e-6 |

表中全部为无量纲相对值，不是百分数。低场最大偏差约 **0.00051%**，
高场最大偏差约 **0.00103%**。没有删去最大误差点，也没有把 P95 当作全体通过。
这里的“闭合”表示原百分比级差异已经得到解释、全场重建达到上述实测精度；
**不表示逐位一致、耦合 Jacobian 已通过或 D4 自洽曲线已验收**。
少量 ppm 级偏差的来源仍未确定，不能未经证据归因于手册常数舍入。

## 1. 修正 Measure 的局部索引

前轮直接把 debug 文件的 `Measure` 三个局部槽对应为导出三角形
`[v0,v1,v2]`，实际本算例 T-2022.03 数据需要映射为 **[v0,v2,v1]**。
这是本数据格式/网格的已验证映射，不能推成任意 TDR 或三维单元的通用规则。

独立于迁移率模型，用 **12956 个非对称、无钝角邻接顶点的三角形**核对解析
Voronoi 体积：映射后三个分量的误差除以单元面积，最大 9.56e-12，
P95 9.71e-14；不映射时最大及 P95 约为 0.25。
原生网格与状态输出的单元连通性也已检查。

映射后，前轮 4054 个体区单元的最大迁移率误差从电子 **2.53185%**、
空穴 **0.97988%** 降至 **3.79e-6 / 4.46e-6**（相对值）。
原钝角单元及其邻居的百分比级长尾消失。

因此修正[前轮报告](templates_ldmos_ialmob_difference_location_2026-09-11.md)中
“Box Measure 不能直接替代迁移率平均权重”的解释：当时的失败由诊断索引
映射错误造成；正确映射的 Measure 已能重建当前样本的原生单元迁移率。
前轮报告保留为历史证据，以本页为最新结论。

已修复 `audit_templates_ldmos_averagebox_node4492.py` 的目标顶点体积取值。
该修复影响诊断的 `target_local_measure_um2` 及汇总体积；不改边耦合索引、
已冻结的 transport_couples.csv、生产离散或已有 D5 曲线状态。
过去报告引用的该诊断体积不能未经重放继续作为正确值使用。

## 2. 明确迁移率 Plot 的节点映射

本例的原生节点迁移率满足：

`mu_node(i) = sum(mu_element(T) / area(T)) / sum(1 / area(T))`，
其中 T 遍历包含 i 的 Si 单元。

这是按单元面积倒数加权的算术平均。低场、高场两套独立原生输出的全部节点
均吻合至约 1e-12；普通面积加权、等权和调和平均均不能替代它。
这条规则用于迁移率输出对照，不应不经验证套用到所有 Plot 字段。

新增受版本管理的
[`audit_templates_ldmos_ialmob_support.py`](../../scripts/audit_templates_ldmos_ialmob_support.py)，
负责单元到节点的对照，并可选重放体区局部计算核输出的 Measure 加权。
输出目录必须是新目录；JSON/CSV 保存所有点及输入 SHA256。

## 3. 界面法向场要按单元内离散距离梯度求值

使用原网格计算每个 Si 顶点到 Si/SiO2 线段的最近距离 d，在三角形内
计算线性插值的 `grad(d)`。本次重建的有效低场输入为：

`F_normal(T) = abs(E(T) dot grad(d)(T))`。

**不额外把离散的 grad(d) 归一化**。它一般不严格满足单位长度；强行归一化
会破坏当前原生离散行为。距离衰减仍使用各顶点自己的距离，Nd/Na/n/p
分别使用各顶点值；各顶点晶向使用原生 NearestInterfaceOrientation 对应参数组。
局部求值后按正确 Measure 加权，再按上一节规则回映射节点。

| 法向处理，同一 10515 单元 | 电子 P95 / 最大相对误差 | 空穴 P95 / 最大相对误差 |
| --- | ---: | ---: |
| 几何节点法向 + 原生晶向/单元电场 | 1.34e-3 / 1.31e-1 | 1.97e-3 / 1.65e-1 |
| 归一化的单元距离梯度 | 7.07e-5 / 7.51e-2 | 1.40e-4 / 8.99e-2 |
| 未归一化的单元距离梯度 | 1.21e-9 / 3.93e-6 | 2.65e-6 / 5.07e-6 |

保留了绝对值、只取正向、只取负向三种诊断假设。只有绝对值处理能同时覆盖
两种场方向；没有为了匹配删除负 Enormal 区域。
这次平直界面样本扩至 **642 个节点**。原先只有两个正 Enormal 节点
4669/4670 的约 9% / 13% 偏差，经单元求值和正确回映射，降至约
3.70e-6 / 1.86e-9（最大相对值）。

原生 InterfaceOrientation 输出包含 619 个 {100}、143 个 {110} 界面顶点；
Si 节点的 NearestInterfaceOrientation 为 4968 个 {100}、755 个 {110}。
仅按最近线段推断参数组会在拐角附近发生误配，本轮最终重建使用原生标签。
后续主求解器接入仍需独立实现并核对晶向标签，不能在自洽求解时依赖参考结果。

## 4. 高场组合的默认行为

依据本地 T-2022.03 手册第 448、451–453 页，并以原生高场输出交叉验证：

- 原始 HFS 驱动为 GradQuasiFermi；触及有效电极顶点的单元默认回退到电场。
- `RefDens=1e12` 的式 376 是 QF 梯度模长与界面平行电场模长的**标量插值**，
  不能先混合向量再取模长。
- `ParallelToInterfaceInBoundaryLayer(PartialLayer)` 默认开启：与界面共边的
  单元只取切向 QF 分量。忽略该行为时，电子最大相对误差达 11.25%。
- 本网格接触标记的边不能当作普通外边界进行上述投影。特别是 `th_lat`：
  本次不把它当作有效电极触发电场回退，也不把它的接触边当作普通外边界。
  错将其视作普通外边界时，单元 676 等位置残留最高 0.735% 电子误差。

加入上述处理后，全部单元高场重建达到首页统计。Caughey–Thomas 仍使用实际
参数文件的 beta0=1.109/1.213、alpha=0、vsat0=1.07e7/8.37e6 cm/s。
HFS 在每个单元顶点低场迁移率上计算，再进行 Measure 加权。
仅看低场单元平均值、在平均后统一施加 HFS 不等同于这个过程。

其中 `th_lat` 边的区别是本例固定状态对照支持的诊断结论，后续接入前需要将
网格接触与有效电极分类写入明确合同，避免依据热接触名称硬编码。
本轮没有开启自热，也没有据此修改生产电场回退配置。

## 环境、原生控制及验证

延续同一网格、参数和 Vg8/Vd40 保存状态：10241 节点、19782 三角形，300 K，
Fermi、OldSlotboom、SRH/Auger、IALMob(AutoOrientation FullPhuMob
PhononCombination=1)。本地独立核使用 UCRT64 **Release，-O3 -DNDEBUG**，
gprof 关闭，程序 SHA256 延续前轮冻结清单。本轮局部求值不调用 SparseLU。
状态转入计算核时用 SI；本报告比较的原生迁移率单位为 cm²/(V·s)。

在用户已授权的本机 Sentaurus VM 隔离目录新增两组 Load/Plot 作业：
`fields_low`、`fields_high`，分别关闭/开启 HFS，输出单元 ElectricField
和节点晶向。两组成功，均 **0 Newton**；状态和原有迁移率与前轮逐项一致。
从电势计算的 P1 电场与原生单元电场也已独立检查。

第一版请求 eEparallel 的 `/Element` 输出被 Sentaurus 拒绝，退出码 5；
保留失败日志，按手册的字段位置改为节点输出后重跑成功。
不将失败尝试计为完成的对照，也不使用失败 TDR。

本轮代码检查包含诊断回归 5/5、Release 的 ascii_sources 和
ialmob_probe_protocol；未修改生产 C++，没有重跑整套曲线或将前次
758/758 CTest 归为本轮新增结果。

## 复现及证据

证据根目录：`reference_staging/templates_ldmos_followup_20260911/ialmob_closure_r2/`。
原始网格、TDR 和生成输出保持 Git 忽略。

- `bundle_v2/`、`raw/`、`fields_low_export/`、`fields_high_export/`：成功原生作业；
  `raw_failed/` 保留首次失败；原始输入沿用前轮的 bundle。
- `support_low/`、`support_high/`：受版本管理的新审计入口重放结果。
- `measure_geometry_qualification.json`：12956 单元的独立几何核对。
- `raw_distance_gradient_replay/`：最终低场局部输入、计算核输出及全体单元/节点 CSV。
- `hfs_contact_edge_exclusion_summary.json`、`e_contact_edge_exclusion.csv`、
  `h_contact_edge_exclusion.csv`：最终高场重建；其它 hfs 文件保留诊断对照。
- `closure_manifest.json`：状态不变性、0 Newton 检查及最终证据/工具哈希。

单独复现高场节点输出映射（从 worktree 根，输出目录必须不存在）：

```powershell
D:\msys64\ucrt64\bin\python.exe -X utf8 scripts/audit_templates_ldmos_ialmob_support.py --mesh reference_staging/templates_ldmos_sentaurus2022/phase01_original_20260826_02/stage1_v4/vela_exact_topology/mesh.json --element-export reference_staging/templates_ldmos_followup_20260911/ialmob_difference_location/native/baseline_element_export --node-export reference_staging/templates_ldmos_followup_20260911/ialmob_difference_location/native/baseline_node_export --output reference_staging/unique_ialmob_support
```

低场全局重建实验脚本为证据父目录中的
`probe_ialmob_raw_distance_gradient_r2.py`，高场为
`probe_ialmob_hfs_contact_edges_r2.py`。这些脚本依赖本次原生输出、外部网格
及上游诊断文件，并非独立发布的自洽运行入口；应连同哈希清单保存。

## 下一阶段边界

下一阶段可以依据本次已定位规则设计单元顶点 IALMob 接口和耦合 Jacobian。
必须包含独立晶向/距离构造、Nd/Na 与两类载流子导数、单元场依赖、HFS
投影/回退的分支验证及局部到全局装配；之后再执行 D4 原曲线门限验收。
本轮没有实现该耦合接口，也没有启动 D4 自洽曲线。无需用户补充现有参数或手册。
