# SimpleMOS Enormal 完整曲线续算

日期：2026-09-13。**本轮完成：剩余 102 个状态全部补齐，408/408 状态、204/204 全场双初始化及 204/204 密度/SRH 重构通过。原有 306 态文件逐字节不变，四次中断现场保留。** 已完成 Id–Vg 与逐栅压物理场汇总和图形检查。HFS 和原始全部工况不在本轮续算范围内。

## 续算起点与冻结条件

9 月 13 日检查时，原求解器、两条驱动路径和后处理协调器已停止。原生 204 个目标已齐；Vela 完成 306/408 个状态，完成账本无数值失败记录。最后结果约在 9 月 12 日 22:57，现有日志不能确定停止原因。

| 工况 | 连续初始化完成数 | 原生相容独立初始化完成数 |
|---|---:|---:|
| n19，Vd=0.05 V | 51/51 | 51/51 |
| n19，Vd=1 V | 51/51 | 51/51 |
| n23，Vd=0.05 V | 24/51 | 31/51 |
| n23，Vd=1 V | 24/51 | 23/51 |

四个没有最终状态码和 `result.json` 的目标为：n23 两个 Vd 的连续路径 Vg=0.32 V，以及独立路径低 Vd 的 Vg=0.62 V、高 Vd 的 Vg=0.46 V。它们属于中断记录，未计为合格或数值失败。原目录逐文件复制到 `build-release/enormal_curves_20260912/resume_20260913/interrupted`，校验哈希一致后，才允许按原输入种子重算。

使用原来的显式 PhuMob + Enormal、300 K、OldSlotboom、掺杂 SRH、原 n19/n23 网格、硅侧输运和已验证 Poisson/介电几何组合。Vela 为 UCRT64 Release、Eigen SparseLU/COLAMD、四次线性修正及统一 split 状态；电流单位 A/μm、等效宽度 1 μm。源代码、隔离二进制、物理参数和全部接受门槛保持不变。

续算入口为 [resume_simplemos_enormal_curves_20260913.py](../../scripts/resume_simplemos_enormal_curves_20260913.py)。保留原驱动的完成状态缓存及最多一次同偏置数值失败重载规则。306 个旧完成态的全部文件另行冻结，两个续算分支结束时分别检查其逐字节不变性。

现场归档与封存已完成，见 [续算前置证据](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/resume_20260913/preflight_evidence.json)、[中断目录清单](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/resume_20260913/interrupted.json)和 [306 态身份](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/resume_20260913/completed_evidence.json)。9 月 13 日 09:54 启动两条驱动和后处理协调器；进程记录在 `build-release/enormal_curves_20260912/resume_active_run.json`。这只表示运行已启动，不表示 102 态已完成。

## 本轮输出与完成条件

目标仍为 204 个偏置 × 两条初始化路径，共 408 个状态。全部状态、原生对应点及全场双初始化通过后，才生成最终 Id–Vg 对比、逐 Vg 物理场和 SRH 重构汇总，以及独立源代码/二进制归档。

物理场比较使用全部硅节点；密度误差使用 log10(Vela/native)，电势及准费米势使用 Vela−native。SRH 公共正硅体积积分用于比较速率和状态差异，不等同于端口 Id 差或原生自身的源项体积。

11:02 完成全部续算：102 个剩余目标补齐，合计 408/408 个完成态合格；完成账本数值失败和重载均为零。两条路径分别完成了续算后的身份复核，原有 306 个状态的全部冻结文件逐字节不变，见 [连续路径证据](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/resume_20260913/continuation_resume_evidence.json)及 [独立路径证据](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/resume_20260913/native_resume_evidence.json)。四次行政中断仍单独保留。204 个工作点均有配对电流，最大初始化相对差 7.677e-10；此时全场配对复核和原生物理场后处理仍在进行。

## 完整 Id–Vg 对比

204/204 个原生目标、408/408 个 Vela 状态及 204/204 个完整双初始化对均通过原数值条件。误差定义为 100×(Id_Vela/Id_Sentaurus−1)，下表为延续路径结果，单位 %。

| 器件 | Vd/V | 点数 | 最小误差/% | 最大误差/% | 最大绝对误差处 Vg/V |
|---|---:|---:|---:|---:|---:|
| n19 | 0.05 | 51/51 | -0.036697355 | +0.004437068 | 0 |
| n19 | 1 | 51/51 | +0.001306597 | +0.005448560 | 0.16 |
| n23 | 0.05 | 51/51 | -0.020969929 | +0.005789136 | 0 |
| n23 | 1 | 51/51 | -0.870875909 | +0.006605729 | 0 |

本轮完整曲线没有在控制点之间出现更大的 n23 电流误差。n23 高 Vd 深关断仍是最差点；相较 PhuMob 阶段的 -0.808135237%，恢复 Enormal 后为 -0.870875909%。全部实测点的绝对百分比误差小于 1%，但本轮没有新增或放宽误差门槛，`comparison_qualified` 仍只表示原生资格、Vela 数值资格与双初始化条件通过。不能据此宣称全部局部物理场、HFS 或原始 0–2.5 V 工况已经验收。

逐点结果见 [comparison.csv](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/comparison.csv)，汇总及身份见 [summary.json](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/summary.json)和 [comparison_evidence.json](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/comparison_evidence.json)。

![完整 Id–Vg 与电流误差](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/curves.png)

## 原门槛下的最终数值检查

| 项目 | 本轮最大值 | 原门槛 |
|---|---:|---:|
| 自由载流子逐行残差比值 | 4.134775e-7 | 1e-6 |
| 端口残差积分/电流提取相对差 | 4.440892e-16 | 1e-8 |
| 双初始化 Id 相对差 | 7.676834e-10 | 1e-6 |
| 双初始化 ψ 最大差/V | 1.509315e-11 | 1e-6 |
| 双初始化 φn 最大差/V | 1.086362e-10 | 1e-6 |
| 双初始化 φp 最大差/V | 3.482722e-9 | 1e-6 |
| 双初始化密度最大相对差 | 1.347177e-7 | 1e-4 |

以上是两种 Vela 初始化之间的差异，不是 Vela 对原生的物理场误差。全部原生偏置、KCL、Vela 独立全局闭合、源项及端口资格沿用原判据。低于原 source floor 的项不因此获得新的逐项相对误差资格。

## 全部 0–1 V 的原生物理场差异

下表在全部硅节点、每条 51 点曲线上取各指标的独立最大值；各列极值不一定来自同一 Vg 或节点。电势和准费米势为绝对差，密度为 |log10(Vela/native)|；SRH 为公共正硅体积权重的归一化 L1 差。

| 工况 | ψ/μV | φn/μV | φp/μV | n/dex | p/dex | SRH L1 |
|---|---:|---:|---:|---:|---:|---:|
| n19，Vd=0.05 V | 2.5696 | 0.14386 | 38.337 | 4.3340e-5 | 6.4436e-4 | 1.3286e-5 |
| n19，Vd=1 V | 3.4898 | 2757.34 | 4420.66 | 0.046365 | 0.074274 | 4.0966e-6 |
| n23，Vd=0.05 V | 2.6334 | 0.14186 | 38.722 | 4.4956e-5 | 6.5082e-4 | 1.5900e-5 |
| n23，Vd=1 V | 3.5556 | 3827.84 | 4845.89 | 0.064263 | 0.081416 | 4.5158e-6 |

低 Vd 的场差很小；高 Vd 下局部准费米场仍有毫伏差，不能将小于 1% 的 Id 误差理解为所有局部物理量已对齐。n23 的两个具体例子为：

- Vd=1 V、Vg=0 V、节点 792：φn 差 +3.82784 mV，原生/Vela 电子浓度为 1137.016/980.626 cm^-3，约 -13.75%；该处 SRH 约为 -4.2154e16 cm^-3 s^-1。
- Vd=1 V、Vg=0.08 V、节点 1056：φp 差 +4.84589 mV，原生/Vela 空穴浓度为 340.942/411.242 cm^-3，约 +20.62%。这些是局部低载流子密度位置，本轮没有新增同扰动因果校准来将它们直接归因为 Id 偏差来源。

204/204 个状态的密度重构检查通过，在导出精度内密度重构误差为零；稳定 SRH 公式与实际生产积分项的最大归一化 L1 差为 4.082451e-16。原生对比的 SRH 公共体积归一化 L1 最大为 1.590002e-5；公共体积净源电荷等效差最大绝对值为 6.868743e-22 A/μm，位于 n23、Vd=1 V、Vg=0.4 V。该电荷等效差不是端口电流差，也不是两套求解器各自原始源项体积之差。

![原生物理场误差曲线](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/fields/field_error_curves.png)

完整数据含逐 Vg 加权 RMS、均值、最大值和节点，见 [fields.csv](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/fields/fields.csv)、[SRH](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/fields/srh.csv)、[热点节点](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/fields/nodes.csv)及 [各曲线极值位置](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/fields/curve_maxima.csv)。

## 完成与保留限制

11:07 已完成 [续算总证据](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/resume_20260913/completion_evidence.json)和 [完成计数](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/resume_20260913/completion_summary.json)。完成尝试为 408、数值失败与数值重载为零；四次行政中断单独保留，不声称整个计算从未中断。源代码和隔离二进制已独立归档，SHA256 为 `dfdb096a9fe75c0be11c9381a46b345deb1ef7f9f9e5a06c05bb468dc35a6f46`。电流及物理场 PNG 已实际查看，标签和曲线显示正常，并保留对应 PDF。

本轮未修改 C++ 求解器、模型参数、默认开关或接受门槛；验证采用实际 102 态续算、204 对全场检查、204 态密度/SRH 重构及全部冻结身份核验，没有为本轮脚本和报告变更重复构建求解器。原生 PhuMob 内部截断搜索及数学截断下的弱空穴单元迁移率资格缺口仍保留。本轮未启动 HFS、原始 n17–n24 的 0–2.5 V 全矩阵，也未提交或推送分支。
