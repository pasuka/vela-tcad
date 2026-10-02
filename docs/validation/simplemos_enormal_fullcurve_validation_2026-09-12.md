# SimpleMOS Enormal 十六点控制与完整曲线

**后续完成记录：2026-09-13 已补齐全部剩余状态，408/408 状态、204/204 全场双初始化及密度/SRH 重构通过。最终曲线、物理场数据和中断保留证据见 [9 月 13 日续算完成报告](simplemos_enormal_curve_resume_2026-09-13.md)。下文保留此前执行及中断时的记录。**

日期：2026-09-12；状态检查更新于 2026-09-13。十六点控制通过，原生完整曲线已取回并通过资格检查；Vela 完整曲线目前中断，未完成。此页不将部分结果记为完整曲线完成。

2026-09-13 09:34（北京时间）重新扫描得到 306/408 个完成且合格的状态，完成账本无失败记录；最后结果约在 2026-09-12 22:57。原求解器、两个驱动程序及后处理协调器均已不在运行。另有四个无 `result.json` 的中断目录，不能计为合格或求解失败；当前日志不能确定停止原因。剩余 102 个状态全部属于 n23：低 Vd 的延续/独立路径分别完成 24/51、31/51，高 Vd 分别完成 24/51、23/51。n19 两个 Vd 的两条路径均为 51/51。后处理文件中残留的 `waiting_for_two_curve_arms` 是旧状态，不代表协调器仍然运行。

本次另完成 [122 个已有配对点的全场复核](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/partial_dual/status_20260913/evidence.json)：122/122 通过原双初始化门槛。最大 Id 相对差 7.677e-10，ψ/φn/φp 最大差分别为 1.510e-11/1.087e-10/2.559e-9 V，密度相对差 9.897e-8。该快照同时核对了已有输入身份；不替代剩余目标和最终全曲线对照。本次为状态检查，没有重启仿真。

前置为 [Enormal 生产候选与同扰动](simplemos_enormal_production_candidate_2026-09-12.md)。使用同一冻结的显式 `phumob_lombardi / element_box_phumob / element_distance_gradient` 程序，原 n19/n23 网格，300 K、OldSlotboom、掺杂 SRH；HFS 关闭。Vela 为 UCRT64 Release、Eigen SparseLU/COLAMD、四次线性修正和 split 高低位状态。原生为 T-2022.03-SP2、ExtendedPrecision(128)、Super。电流 A/μm、等效宽度 1 μm。

## 已完成的原生曲线

n19/n23 × Vd=0.05/1 V，Vg=0–1 V、步长 0.02 V。四次原生仿真均退出 0；204/204 目标点通过偏置和 KCL 检查，全部 TDR 已导出。与既有 Vg=0.8/1.0 V 八点重合的 Id 完全一致。未改变物理参数或加入高场模型。

输入 SHA256：`d45d5faefb5974a6579d9ebcaddb166d7432a1ff47bf8ae16c7b6cbc971ef215`。
结果 SHA256：`01b7f3a4df9c1bd39e493b94e4308f99001a6997aa52a68ee6f15c5401d8424d`。

证据：[原生点](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/native_points.csv)、[原生资格](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/native_evidence.json)、[重合点](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/native_overlap.csv)、[字段导出](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/export_evidence.json)。

另以两批合格原生完整曲线比较恢复 Enormal 的有限物理影响：四个工况的网格逐字节一致，去除 Enormal 后的其余 Physics 与 Math 一致。n23 低/高 Vd 的原生 Id 在 Vg=0 时分别下降 11.795%/7.175%，在 Vg=1 V 时分别下降 30.146%/29.181%；n19 全曲线变化范围为 -30.610% 至 -8.961%。这是 **原生 Enormal 对原生 PhuMob 的模型效应**，不是 Vela 对原生的误差，也不是小扰动导数。数值、条件身份与图见 [逐点模型效应](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/native_restoration_effect/points.csv)、[证据](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/native_restoration_effect/evidence.json) 和 [图](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/native_restoration_effect/model_effect.png)。

## 已完成的十六点双初始化

新增 Vg=0/0.2 V 八点，分别从同偏置合格 PhuMob 保存态和原生势场相容初始化；16/16 状态首次通过，无失败、无重载。与此前 Vg=0.8/1.0 V 八点合并，同一候选共 **32/32 状态、16/16 双初始化通过**。

新增八点最大逐行比值 5.104e-7，低于原 1e-6 门槛。两初始化 Id 在导出精度内一致；最大电子/空穴准费米势差分别为 3.795e-11/6.083e-9 V，最大密度相对差 2.353e-7，均通过原门槛。电势的导出差为零不代表任意精度下严格相等。

新增控制点电流误差为 100×(Id_Vela/Id_Sentaurus−1)，单位 %：

| 器件 | Vd (V) | Vg=0 V | Vg=0.2 V |
|---|---:|---:|---:|
| n19 | 0.05 | -0.036697355 | +0.003831977 |
| n19 | 1.0 | +0.004726399 | +0.005442785 |
| n23 | 0.05 | -0.020969928 | +0.003566582 |
| n23 | 1.0 | -0.870875909 | -0.021325359 |

n23 高 Vd 深关断误差比前一阶段 PhuMob 的 -0.808135237% 略大。数值资格通过不代表每次恢复物理模型都会单调改善相对 Id；不因该差异重设门槛、调物理常数或采用拟合截断值。

证据：[逐次账本](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/low_control_attempts.csv)、[低栅压对照](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/low_control_comparison.csv)、[十六点汇总](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/control_summary.json)、[资格](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/control_evidence.json)。

## 尚未完成的完整曲线

两条路径分别为从合格 0.8 V 状态连续向下/向上扫压，以及每个偏置独立使用原生相容初始化。目标为 204 点 × 两条路径、408 个状态。每次保留全部原始尝试；仅允许同偏置一次重载，失败态不得传递到下一个偏置。全行、端口、全局源闭合和双初始化条件保持不变。

本轮 `vela_contract.json` 明确未新增绝对电流误差门槛。因此汇总中的 `comparison_qualified` 表示对应原生点及两条 Vela 路径通过数值和双初始化条件，可以用于对比；它本身不表示 Id 误差小于 1%，也不等于原始全部物理模型和工况已验收。原生 Id 误差需单独报告实测最大值与覆盖范围。

源文件、所有当时使用的 Python 辅助脚本及隔离二进制另有 [补充身份冻结](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/supplemental_source_evidence.json)。完成后将统一生成电流曲线、逐 Vg 的原生物理场差、独立归档和最终资格；此时尚不声明 204/204 双初始化通过。

中途首个重合点 n19、Vd=0.05 V、Vg=0.42 V 已通过原双初始化门槛：Id 相对差 4.219e-15，ψ/φn/φp 最大差分别为 2e-15/6e-17/1.566e-14 V，密度相对差 6.056e-13。这个 [独立快照](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/partial_dual/first_overlap/evidence.json)仅覆盖该已完成点，不能替代最终 204 点资格。

第二个 [中途快照](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/partial_dual/mid_overlap/evidence.json)覆盖 31 个已完成配对点，31/31 通过原全部双初始化门槛。最大 Id 相对差 1.092e-12，ψ/φn/φp 最大差分别为 4.286e-14/4.992e-11/2.559e-9 V，密度相对差 9.897e-8。结果为两种 Vela 初始化之间的差，不能解释为 Vela 对原生物理场的误差；尚未完成的目标不在这个快照中。

高场饱和当前只做了 [公式、接口和参数预检](simplemos_hfs_restore_preflight_2026-09-12.md)，没有 HFS 器件仿真或生产启用。实际恢复仍以本阶段通过为前置，随后才进入原始 n17–n24、0–2.5 V 的全部工况。
