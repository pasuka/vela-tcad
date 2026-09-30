# HFS 十六点双初始化与低栅压验证

日期：2026-09-14。分支 `codex/simplemos-sdevice-validation`。本页冻结十六点阶段结果；完整曲线执行情况另行汇总。

## 结论与范围

沿用 [HFS 八点生产候选](simplemos_hfs_production_candidate_2026-09-14.md) 的同一冻结程序，补齐 n19/n23 × Vd=0.05/1 V × Vg=0/0.2 V。合并原 Vg=0.8/1 V 八点后，**32/32 状态、16/16 全场双初始化通过**，32 次首次尝试均合格，无失败、无重载。未修改 C++、材料参数、模型默认或接受阈值。

这项通过是收敛、全行闭合、端口与初始化一致性的资格。**n23、Vd=1 V、Vg=0 的 Id 误差为 −1.050499105%，已超过 1%**，因此不能宣称全范围电流达到 1% 精度目标。两种初始化得到的该点 Id 相对差为 8.882e-16，该差异具有自洽状态的可复现性。

证据根目录：`reference_tcad/simplemos_sentaurus2022/hfs_curves_20260914`；执行目录：`build-release/hfs_curves_20260914`。

## 固定模型、几何与数值契约

原始 n19/n23 网格，300 K，Boltzmann、OldSlotboom、PhuMob→Enormal→逐顶点 Canali HFS→box 平均，SRH(DopingDependence)。原生默认 PartialLayer 与任意接触顶点的电场替代保持；低于 1 V/cm 的驱动力不饱和。沿用数学 PhuMob 截断底、现代 Vela 常数、signed Si 三项 Poisson 电荷体积、cell-material 介电装配及生成的 Delaunay transfer box 几何；SRH 体积不变。

显式 split 状态契约、100 位宽残差和迁移率/端口求值、double 解析 Jacobian、Eigen SparseLU/COLAMD 与四次线性修正保持。冻结程序：`build-release/hfs_candidate_20260914/v1/bin/vela_example_runner.exe`，SHA256 `d9848e77096647530f1c3f0622542a2ce6e9f8be6dae8b8507d9ff27cb36dbe3`。独立源码/程序归档沿用前一报告，不以当前目录时间戳代替身份检查。

原接受门槛：载流子全行比值 1e-6，KCL/Id 与端口一致性 1e-8；独立全局闭合容差 1e-6、既定源项 floor 1e-10；双初始化电势/准费米势最大差 1e-6 V、密度相对差 1e-4、Id 相对差 1e-6。每失败点最多一次同偏压重载，失败结果必须保留。本批未用重载。

## 原生参考与身份

从已合格 HFS Vg=0.8 V 原生检查点重闭合，降至零，再以 0.02 V 间隔升至 1 V，保存每个目标的状态和物理场。T-2022.03-SP2，ExtendedPrecision(128)、Super、Digits=12、RhsMin=1e-20、Iterations=40。四条原生曲线 **204/204 目标合格**；四次进程退出均为零，目标偏压与 KCL 均通过。它们是高精度参考，不是原 Applications Library 默认求解参数的身份复现。

下载归档 SHA256：`e9a076253731eb6b417759bb69bcbf6dcb9105f3d01eb2107ad8a4062cb3e7fb`。与前八点的八个原生 Id 重叠差均为零。字段通过本地 HDF5/TDR 导入器导出，没有插值。虚拟机目录 `/tmp/vela_simplemos_hfs_curves_20260914` 的本批任务均已完成。

## 十六点数值与电流结果

| 指标 | 最大值 |
|---|---:|
| 全行比值 | 2.211846194e-7 |
| 端口相对差 | 4.440892099e-16 |
| 双初始化 Id 相对差 | 6.596740931e-10 |
| 双初始化电势/准费米势差 | 2.37097259e-9 V |
| 双初始化密度相对差 | 9.171331739e-8 |

电流单位 A/μm，误差定义 `100*(Id_Vela/Id_Sentaurus−1)`：

| 器件 | Vd (V) | Vg=0 | Vg=0.2 | Vg=0.8 | Vg=1 |
|---|---:|---:|---:|---:|---:|
| n19 | 0.05 | −0.053293189% | +0.003705655% | +0.003453364% | +0.001073446% |
| n19 | 1 | +0.004052368% | +0.005421344% | +0.002790587% | +0.001087522% |
| n23 | 0.05 | −0.025931291% | +0.003439335% | +0.005774106% | +0.005197263% |
| n23 | 1 | **−1.050499105%** | −0.036689876% | +0.006626328% | +0.005098975% |

n23、高 Vd、Vg=0 的 Vela Id 为 2.8295265668e-16 A/μm。相对上一阶段 Enormal，原生 HFS 自洽 Id 变化 −17.034193%，Vela 变化 −17.184528%，原电流误差由 −0.870875909% 增至 −1.050499105%。这不是一个额外拟合因子或已定位的离散根因。早期单路径差值证据保留其当时的 dual-pending 标记，最终资格以十六点证据为准。

## 低栅压关键场差

8/8 低栅压密度及 SRH 重建通过：密度重建最大相对差为零，SRH 生产积分重建 L1 最大相对差 3.5833e-16。场指标比较全部 Si 节点，使用共同正的 barycentric Si 体积加权；这些 SRH 电荷等价量不等同于端口电流归因，也不是两个求解器各自源积分的直接差。

| n23 | 最大 ψ 差 (μV) | 最大电子 QF 差 (mV) | 最大空穴 QF 差 (mV) | SRH 加权 L1 相对差 |
|---|---:|---:|---:|---:|
| Vd=.05, Vg=0 | 1.612715 | 0.000137809 | 0.0399641 | 1.16773e-5 |
| Vd=.05, Vg=.2 | 1.681911 | 0.000141318 | 0.0400020 | 1.20863e-5 |
| Vd=1, Vg=0 | 2.759832 | 4.190201 | 4.776897 | 3.70393e-6 |
| Vd=1, Vg=.2 | 2.759854 | 4.190128 | 4.216703 | 3.76364e-6 |

高 Vd 最大电子 QF 差在节点 792，最大空穴 QF 差在节点 1056；Vg=0 最大电子/空穴密度对数差为 0.0703507/0.0802577 dex。局部场差仍待归因，不能因双初始化一致而当作跨求解器物理场已完全一致。

## 原始矩阵范围与后续边界

原始矩阵为 n17–n24 八器件 × 两个 Vd，0–2.5 V、0.05 V 步长，共 816 点。当前四条 0–1 V、0.02 V 曲线计划仅有 44 个点与原始矩阵的器件/偏压严格重合；不可把插值计为新增通过点。

独立核对八种网格后，薄氧化层 n17/n18/n21/n22 有 958 个自由 Si 节点、1916 个载流子行；原四曲线辅助脚本固定 907/1814。因此新增独立、按实际网格冻结行数的适配器，未改旧冻结脚本。在 24 个既有状态上，其完整独立验收记录与旧计算逐项一致，八种网格行集合检查通过。这只完成验收工具适配，不代表其余器件的自洽计算通过。

下一步按已通过的十六点状态资格，计算四条完整 HFS 曲线的两条初始化路径，共 408 个状态，确定深关断超过 1% 的偏压范围并做完整场汇总。原始全部工况仍待此阶段结果；不把前一阶段的空穴单元截断缺口或 1e-6 V Jv 跨分支失败重新标为通过。

## 可复核入口

- [十六点比较](../../reference_tcad/simplemos_sentaurus2022/hfs_curves_20260914/sixteen_comparison.csv)、[汇总](../../reference_tcad/simplemos_sentaurus2022/hfs_curves_20260914/sixteen_summary.json)、[证据](../../reference_tcad/simplemos_sentaurus2022/hfs_curves_20260914/sixteen_evidence.json)。
- [原生资格](../../reference_tcad/simplemos_sentaurus2022/hfs_curves_20260914/native_evidence.json)、[重叠身份](../../reference_tcad/simplemos_sentaurus2022/hfs_curves_20260914/native_overlap.csv)。
- [低栅压字段](../../reference_tcad/simplemos_sentaurus2022/hfs_curves_20260914/low_fields/fields.csv)、[SRH](../../reference_tcad/simplemos_sentaurus2022/hfs_curves_20260914/low_fields/srh.csv)、[字段证据](../../reference_tcad/simplemos_sentaurus2022/hfs_curves_20260914/low_fields/evidence.json)。
- [原始矩阵范围](../../reference_tcad/simplemos_sentaurus2022/hfs_curves_20260914/original_scope_audit.json)、[网格行集合](../../reference_tcad/simplemos_sentaurus2022/hfs_curves_20260914/original_mesh_scope.json)、[适配器核验](../../reference_tcad/simplemos_sentaurus2022/hfs_curves_20260914/matrix_support_evidence.json)。

验证：新脚本语法检查和实际执行通过；参考工具单元回归 163/163 通过。此前 HFS 程序的 209 个 C++ 测试用例为冻结前置证据，本轮未改 C++，未重复构建。未提交或推送。
