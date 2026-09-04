# Genius NPN BJT：Sentaurus 2022 与 Vela 对比算例

## 1. 算例定位

本算例将 Genius TCAD 自带的二维 NPN BJT（`examples/BJT/step1.inp`、
`step2.inp`）转换为 Sentaurus SDE/SDevice 工程，以 SDevice
T-2022.03-SP2 的结果作为基准，再把同一器件结构、节点掺杂和物理参数导入
Vela。算例用于验证三类能力：器件建模转换是否正确、端口直流特性是否一致，
以及空间场、守恒通量和复合模型是否闭合。

这是一套可复现的参考算例，但不是“所有指标在所有网格和所有后处理方法下均已
通过”的校准声明。当前状态必须按下表中的验收层级理解。

## 2. 固定输入

| 项目 | 设置 |
|---|---|
| 器件 | 二维硅 NPN BJT，6 μm × 2 μm，300 K |
| 接触 | 顶部基极与发射极，底部集电极 |
| 偏置 | 发射极 0 V，基极 0.70 V，集电极 0–3 V，步长 0.1 V |
| 粗网格 | 5611 节点、10940 三角形 |
| 局部加密网格 | 15561 节点、30780 三角形 |
| 正式物理模型 | Fermi 统计、OldSlotboom 带隙变窄、掺杂相关迁移率、掺杂相关 SRH、经典 `np-ni_eff²` Auger |
| 电流单位 | A/μm；Genius 的 100 μm 拉伸宽度仅作为元数据，不隐含乘入结果 |

`M0` 是用于排查数值流程的简化物理模型，不承担跨软件精度验收；`M1` 是完成
Sentaurus 参数对齐后的正式对比模型。后续复用本算例时，应以 M1 及其唯一正式
接受态为准。

## 3. 当前验收结果

状态快照日期：2026-09-04。

| 验收层级 | 当前结果 | 主要证据 |
|---|---|---|
| SDE 结构与掺杂 | 通过 | 尺寸、接触、解析掺杂和粗/细网格同源性审计均通过 |
| SDevice 基准曲线 | 通过 | 两套模型各 31 个精确偏置点，三端电流直接导出，KCL 通过 |
| M1 端口曲线 | 通过 | 0.5–3 V 内 Ic、Ib、Ie、β 最大误差分别为 0.00344、0.00210、0.00341、0.00135 decade，均低于 0.05 decade |
| 粗网格空间状态 | 通过 | 3 V 电势、有效电子浓度和有效空穴浓度门槛全部通过 |
| SRH、Auger 与守恒截面 | 通过 | SRH 积分比约 0.996；Auger 积分比 0.9878；0、1、2、3 V 的守恒截面均通过 |
| 粗网格默认节点电流恢复 | 部分通过 | 电子电流通过；弱空穴电流 P95 为 0.828 decade，高于 0.5 decade 门槛 |
| 局部加密网格 3 V | 通过 | 空穴电流 P95 降至 0.270 decade，端口、空间状态、输运和复合单点门槛全部通过 |
| “先单元、后节点”诊断恢复 | BJT 内通过，暂不设为默认 | 粗网格 3 V 空穴 P95 降至 0.0985 decade；跨 MOS/PN 回归结果不一致，尚不足以替换全局默认值 |

因此，本算例可以作为以下两种参考基准使用：

1. **正式端口基准**：使用粗网格 M1 的完整 31 点曲线，状态为通过。
2. **高保真空间基准**：使用局部加密网格检查 3 V 空间场、节点电流和复合率，状态为通过。

粗网格加默认节点恢复的“全量严格验收”仍为未通过，唯一失败项是对总电流贡献
很小的基区空穴电流尾部。这个失败不影响端口电流、三端 KCL、主要空间状态、
SRH/Auger 或守恒截面已经获得的通过结论。

## 4. 关键结论

- 参数对齐后，Vela 与 SDevice 的正式 M1 端口结果已经达到百分之一量级以内；
  3 V 的 Ic 和 Ib 比值分别为 0.99215 和 0.99517。
- 有效区域内的电势、电子浓度和空穴浓度均通过预先登记的空间门槛。全域空穴
  浓度的大误差主要位于 SDevice 参考浓度约 1–100 cm⁻³ 的极低浓度区，不作为
  有效载流子区域的否决项。
- SRH 的早期约 22% 差异来自后处理遗漏 Fermi–BGN 有效本征浓度修正，并非生产
  SRH 模型本身的差异；修复后 SRH 参数与局部公式已闭合到约 0.1% 量级。
- 局部网格加密显著改善弱空穴电流和低浓度尾部，说明粗网格分辨率是差异来源
  之一，不能把残差完全归因于 SDevice 的节点矢量构造语义。
- “先单元重构、再投影到节点”对当前 BJT 和多数 MOS 工况有明显改善，但在 PN
  二极管上没有一致收益，因此仍是默认关闭的绘图诊断选项，不改变有限体积残差、
  守恒边通量或端口电流。

## 5. 推荐使用方式

| 使用目的 | 推荐输入/结果 |
|---|---|
| 快速端口回归 | `vela/configs/m1_collector_sweep.json` 与 `reference_curves/bjt_m1_output.csv` |
| 建模转换审计 | `source/bjt_sde.cmd`、`contracts/geometry_contract.json`、`reports/wp0_wp2_validation.json` |
| 3 V 空间状态对比 | 粗网格正式接受态及 `comparison/spatial_comparison_summary.json` |
| 弱电流与网格敏感性复核 | `source/bjt_sde_local_refined.cmd` 与 `reports/mesh_sensitivity_validation.md` |
| 复合和守恒验证 | `reports/conservative_flux_srh_alignment_report.md` |
| 节点电流恢复研究 | `reports/cell_first_recovery_ab.md` 与 `reports/cell_first_cross_device_validation.md` |

不要使用 M0 评价 Vela 与 SDevice 的物理精度，也不要使用极低载流子浓度节点的
无掩膜最大对数误差代替端口、守恒通量和有效区域空间门槛。

## 6. 可复现流程

主要入口按执行顺序如下：

```text
SDE/SDevice 基准： source/bjt_sde.cmd -> source/bjt_m1_des.cmd
Vela 正式扫描：   vela/configs/m1_model_relaxation.json -> vela/configs/m1_collector_sweep.json
最终比较：        scripts/compare_genius_bjt_sentaurus_vela.py
守恒与 SRH：      scripts/audit_genius_bjt_conservative_sections.py
局部网格验证：    scripts/run_genius_bjt_mesh_sensitivity.py
电流恢复 A/B：    scripts/run_genius_bjt_cell_first_recovery_ab.py
跨器件回归：      scripts/run_cell_first_cross_device_validation.py
```

源代码库仅保存文本输入、归一化曲线、验收报告和精选图片。Sentaurus TDR/PLT、
Vela 状态、VTK 与运行日志保存在忽略目录
`build-release/reference_tcad/genius_bjt_sentaurus2022/`，不纳入版本控制。

## 7. 结果与图片索引

- 端口与总体验收：`comparison/comparison_summary.md`
- 3 V 空间状态：`comparison/spatial_comparison_summary.md`
- 局部网格敏感性：`reports/mesh_sensitivity_validation.md`
- 守恒通量和 SRH：`reports/conservative_flux_srh_alignment_report.md`
- 当前 BJT 电流恢复 A/B：`reports/cell_first_recovery_ab.md`
- MOS/PN 跨器件回归：`reports/cell_first_cross_device_validation.md`
- 器件网格、空间场和曲线图片：`figures/`，其中粗网格综合图片位于
  `figures/coarse_m1_comparison/`

机器可读的算例范围、状态、入口文件和 SHA-256 指纹见
`genius_bjt_sentaurus2022_reference.json`。
