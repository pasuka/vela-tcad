# SimpleMOS 常数约定与局部 SRH 体积的有限自洽验证

日期：2026-09-08。四个 Vg=0.8 V 控制点完成了常数组合的有限自洽对照，最大 Id 相对差降至 **1.12624132e-10**；经固定三步延续，两条初始化路径全部通过原严格门槛。SRH 两节点联合有限修改的 **8/8 端点、4/4 双初始化**通过，局部准费米势明显改善，但全场仍有毫伏级差异。

三个 SRH 候选合计 24 个端点中最终 22 个合格、10/12 双初始化合格；两个单节点分支仍因空穴行未闭合而失败。保留全部初次失败记录，没有降低门槛或修改生产常数、SRH 默认策略。本轮不扩大到 16 工况或完整 0–1 V 曲线。

## 1. 工况、独立性与接受条件

延续[上一轮独立小响应校准](simplemos_constants_srh_response_validation_2026-09-07.md)。n19/n23 对应 NWell Boron 1e17/2e17 cm⁻³，Vg=0.8 V，Vd=0.05/1 V。网格分别 1480/1482 节点、2742/2746 单元；各有 942 个 Si 节点、907 个自由 Si 节点，检查全部 1814 个自由载流子行。几何坐标 μm、状态密度 m⁻³、宽度 1 μm，Id 单位 A/μm。

模型为 300 K、Boltzmann/no-BGN、ni=1.0750038488844236e10 cm⁻³、总杂质 Masetti、掺杂 SRH（τn=1e-5 s、τp=3e-6 s、Nref=1e16 cm⁻³、γ=1）；HFS、表面迁移率、Auger、雪崩、DG 关闭。沿用已验证的自动 Delaunay box、element-box 输运/迁移率、分材料介电系数及三项独立 signed Si Poisson 电荷体积。

Windows MSYS2 UCRT64、C++20 Release，DC 实际使用 Eigen SparseLU/COLAMD；HDF5/TDR、UMFPACK、SPQR 编译可用，本轮未切换 DC 后端。独立切向使用 SciPy SuperLU/COLAMD，加四次 long-double 残差修正。常数程序沿用全部 47 个 core 编译单元统一常数头的隔离构建；SRH 程序只替换其中装配器对象，继续链接其余同约定对象，不与旧生产 core 常数混用。

所有 DC 保持 max_iterations=200、relative_tolerance=1e-7、absolute_tolerance=1e-12、stall floor=1e-9、ψ 更新上限 0.35 V、准费米势上限 0.025 V、contact_basin 初始化、原标量线搜索、四次线性修正。全部载流子行要求 |R|/scale≤1e-6，不排除少数载流子；KCL/Id≤1e-8。全局源闭合继续使用 1e-6 比值及原 1e-10 源下限，低于下限不能宣称相应相对精度。双初始化要求势差≤1e-6 V、密度相对差≤1e-4、Id 相对差≤1e-6，且**两端各自严格合格**。

## 2. 常数有限对照：剩余电流差得到解释

诊断运行同时采用已有 Sentaurus 手册/原生回放给出的 q=1.602192e-19 C、kb=1.380662e-23 J/K、ε₀=8.8542e-12 F/m，未拟合数值。现代 SI 默认仍为 1.602176634e-19、1.380649e-23、8.8541878128e-12。此处只改变常数约定，原 SRH 体积保持不变。

| NWell | Vd/V | 现代 SI 基态 Id 误差 | 手册常数有限实算 Id 相对差（无量纲） |
| --- | --- | --- | --- |
| n19 | 0.05 | +0.00348848% | -8.09796674e-13 |
| n19 | 1.0 | +0.00298960% | -3.04201109e-14 |
| n23 | 0.05 | +0.00588005% | -8.75688411e-12 |
| n23 | 1.0 | +0.00675240% | -1.12624132e-10 |

最后一列为 (Id_Vela/Id_native−1)，不是百分数。两个 Vd 下高 NWell 绝对差、低 NWell 控制及高低配对差均改善；说明这四点在已修正几何/介电/迁移率之后的剩余 Id 差主要来自常数约定。四点 ψ 最大绝对差为 3.05e-15 V，接近双精度和导出精度下限；这不表示器件物理具有该绝对精度，也不能外推到未经计算的偏置。

直接运行 12 次中 10 次合格，原生来源初值的 n19/Vd=0.05、n23/Vd=1 两分支失败。之后从各自原始原生来源初值做冻结 α=0、0.5、1 三步常数延续，12/12 合格，四个终点与直接合格的 Vela 来源路径形成 4/4 严格双初始化。原两次失败及原始状态保留。[直接 DC](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/constants/dc.csv)、[比较](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/constants/comparison.csv)、[延续双初始化](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/constants/continuation/dual.csv)。

原始 `constants/dual.csv` 的 qualified 字段只表示状态距离检查，不能替代两端 DC 合格。最终采用延续路径的严格结果；未回写旧汇总。

## 3. SRH 候选：先独立小扰动，再有限替换

在上述合格手册常数状态上，用几何直接给出的比值定义 V_i(α)=V_original_i×[1+α(V_signedSi_i/V_original_i−1)]。n23 节点 792/1057 对应 n19 的空间匹配节点 791/1056；位置分别约 (0.1707979, 0.3)、(0.005045938, 0.1127039) μm，最大跨网格匹配距离 9.23e-10 μm，小于冻结的 1e-8 μm。

第一节点体积比为 1.17958462195；第二节点 n19/n23 分别为 0.217276367556/0.217274908657。执行第一节点、第二节点和联合三支路，α=±0.001、±0.0005 后才到 α=1。只改变 SRH 源积分及对应六项状态导数；三项 Poisson 体积、输运几何、固定状态通量均保持原定义。[空间与体积表](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh/spatial_scope.csv)。

12/12 预检查通过：零扰动残差逐位一致，独立源公式与固定状态残差差分最大相对差 3.88591e-09（门槛 1e-5）；完整线性切向残差最大 6.03425e-11（1e-8）；弱交叉块与预期体积缩放最大相对差 2.22e-16（1e-12）。同时通过固定输运不变性、强块 Jv 与端口一致性。

52/52 小扰动 DC、24/24 响应检查通过。主响应预先冻结为全部 907 个自由 Si 节点的 (ψ,φn,φp) 向量，独立完整 J 解的预测最大相对差 2.20988e-05，双幅度比例差最大 2.38655e-05，均低于 1e-3；最大偶/奇比 8.08863e-05，最小响应/零漂移比 956.879。未以不敏感的 Id 代替少数载流子场校准。[预检查](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh/preflight.csv)、[响应](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh/response.csv)。

## 4. 有限 SRH 结果与资格

最初 24 次有限 DC 的低 Vd 12 次通过，高 Vd 12 次失败。冻结一次重启，逐支使用原失败状态继续同一有限参数，数值配置与门槛不变，仅增加只读更新诊断。12 次重启中 10 次通过，合并初次合格状态后为 22/24 严格端点、10/12 严格双初始化；全部 24 个选定端点的装配器、端口提取器和探针电流一致，但端口一致不抵销 DC 失败。

联合支路两种初始化的 8/8 端点和 4/4 双初始化均通过：

| NWell | Vd/V | 联合后 Id 相对差（无量纲） | 两初值最大载流子行比值 | 严格双初始化 |
| --- | --- | --- | --- | --- |
| n19 | 0.05 | -3.553824e-13 | 1.365604e-11 | True |
| n19 | 1.0 | +7.238654e-14 | 7.480836e-07 | True |
| n23 | 0.05 | -7.799428e-12 | 2.243897e-11 | True |
| n23 | 1.0 | -5.273681e-11 | 4.287463e-08 | True |

仅就已修改位置，高 Vd 的有符号准费米势差如下。基线为手册常数、原 SRH 体积；联合结果取选定的合格 Vela 来源支路，另一初值已通过上述严格检查。

| NWell | 本地节点 | 场 | 基线差/mV | 联合后差/mV |
| --- | --- | --- | --- | --- |
| n19 | 791 | phin | +0.0038905 | -0.0013185 |
| n19 | 791 | phip | -3.9269620 | -0.0416246 |
| n19 | 1056 | phip | +2.5569416 | +0.9406770 |
| n23 | 792 | phin | +3.8013279 | -0.1016942 |
| n23 | 792 | phip | -0.8760795 | -0.0132953 |
| n23 | 1057 | phip | +3.9013811 | +1.0601269 |

n23 节点 792 电子密度相对差从约 −13.674% 降至 +0.3941%，1057 空穴密度相对差从 +16.289% 降至 +4.1860%。联合后这两点局部 SRH 体积率与原生的相对差分别约 +2.30e-10、−4.53e-9；体积率一致不意味着积分源与全场已全局一致。n23 节点 1000/1009 的 φp 差仍约 19.239/20.802 μV，变化很小。

全部 907 个自由 Si 节点的最大准费米势绝对差仍为：

| NWell | Vd/V | 场 | 基线最大差/mV | 联合后最大差/mV | 联合后最差节点 |
| --- | --- | --- | --- | --- | --- |
| n19 | 0.05 | phin | 9.476354e-06 | 9.477949e-06 | 1162 |
| n19 | 0.05 | phip | 0.1200972 | 0.120097 | 1086 |
| n19 | 1.0 | phin | 2.751718 | 2.751719 | 793 |
| n19 | 1.0 | phip | 3.926962 | 3.853526 | 794 |
| n23 | 0.05 | phin | 2.241089e-05 | 2.248962e-05 | 1163 |
| n23 | 0.05 | phip | 0.1201939 | 0.1201939 | 1087 |
| n23 | 1.0 | phin | 3.801328 | 2.892679 | 794 |
| n23 | 1.0 | phip | 3.901381 | 1.697433 | 1056 |

因此，局部两体积能解释被选热点的大部分偏差，尚未消除邻近节点或其他 SRH 支撑范围的差异。联合状态高 Vd 的 n23 最大电子/空穴密度相对差仍约 11.84%/6.79%；只看电势和端口电流会遗漏这些少数载流子误差。[最终选定状态](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh/restart/selected_states.csv)、[严格双初始化](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh/restart/dual.csv)、[局部场](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh/restart/local_fields.csv)、[全场误差](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh/restart/fields.csv)。

初次 `srh/dual.csv` 同样只是状态距离指标；最终严格判定使用 `srh/restart/dual.csv`，失败态在场表中带 qualified_state=False，未混入上面的合格联合结果。

## 5. 原生固定状态 SRH 行回放与剩余失败

另用已导出的原生电势、准费米势、密度、SRH 体积率、box 几何及原生单元加权 Masetti 独立重建局部 SG 流出项，分别使用准费米势和密度两种公式，共保留 32 个行记录。下述源归一化指标统一使用 max(|净流出|, |原生 SRH×signed Si 体积|) 为分母。这是原生固定导出状态回放，**不是原生同源正负自洽响应校准**，未启动新 sdevice。

n23/Vd=1 的节点 792 电子/空穴行用 signed Si 体积后，源归一化残差约 2.40e-15/8.58e-16；沿用原 Vela SRH 体积则约 0.152244。节点 1057 的空穴行分别约 2.42e-12 与 3.60246，低 NWell 对应少数载流子行也支持同一几何解释。节点 1057 的多数电子行有严重大通量消减：绝对边通量和约 4.68e12 particles/(cm·s)，而原生积分源约 −0.112 particles/(cm·s)，准费米势导出回放净通量约 −4.83。该行无法在源尺度上反推体积；其异常/负推断体积保留在表中，不当作物理体积结论。[32 行回放](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh/native_rows/row_replay.csv)。

最终仍失败的是 Vd=1、仅改第二节点、原生来源初始化的两个分支：

| NWell | 空穴节点 | 最终行比值（门槛 1e-6） | 原始=截断更新/V | 最后 Poisson 块残差 |
| --- | --- | --- | --- | --- |
| n19 | 1056 | 1.61955773e-06 | -5.470904e-10 | 8.823324e-11 |
| n23 | 1057 | 4.05484928e-06 | -1.089119e-09 | 8.794431e-11 |

两点最后均为 carrier_row_convergence_line_search_rejected，13 次尝试未接受；原始更新没有触及 0.025 V 上限。对应行 raw J·dx+R 约 −6.31e-41/+1.10e-39（求解器缩放单位），最优但被拒候选仍稍微降低该空穴行残差。diagnostic 的 applied_step_V 在 best_rejected_candidate=1 时描述候选更新，并未提交到最终状态。这表明更新截断不是这里的限制；下一步应独立核对接近 Poisson 求值精度下限时的全局 merit/线搜索决定。单行改善及本次双精度线性诊断仍不足以证明任意完整步可接受。[剩余失败](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/srh/restart/remaining_violations.csv)、[原始/截断/候选更新](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/remaining_update_trace.csv)。

## 6. 完成范围与后续边界

本轮有限常数对照、SRH 分项/联合对照、双初始化、原生固定状态源行审计已执行。九项诊断输入拒绝测试通过（网格数、缺项、接触节点、重复、非正比值、多余记录、幅度越界、缺少幅度、旧/新接口冲突），各阶段哈希复核通过。未更改生产默认，未重跑生产全量 CTest；上一轮生产测试状况仍见[几何生产报告](simplemos_generated_box_mobility_validation_2026-09-07.md)，本报告不将其重新计为今日测试。

下一阶段先对上述两个失败态做固定原始方向的 trial residual/merit 独立精度扫描，区分 Poisson 消减与少数载流子真实下降，不改接受阈值；并在合格联合状态上审计 n23 的 794/1056、n19 的 793/794 及低 Vd 1086/1087 所在邻域的 SRH 支撑。只有邻域源定义与同扰动校准支持后，再冻结更完整的 SRH 体积候选；原生同源自洽响应资格尚未取得。常数兼容策略如进入生产，应做显式配置并保留 SI 默认。Vg=1、16 工况和全曲线仍需另行逐步验证。

脚本：[常数有限对照](../../scripts/validate_simplemos_constants_finite_20260908.py)、[常数延续](../../scripts/continue_simplemos_constants_finite_20260908.py)、[SRH 隔离构建](../../scripts/build_simplemos_srh_finite_20260908.py)、[SRH 校准与有限对照](../../scripts/validate_simplemos_srh_finite_20260908.py)、[一次重启](../../scripts/restart_simplemos_srh_finite_20260908.py)、[原生行审计](../../scripts/audit_simplemos_native_srh_rows_20260908.py)、[输入校验](../../scripts/test_simplemos_srh_finite_inputs_20260908.py)、[报告生成](../../scripts/report_simplemos_constants_srh_finite_20260908.py)。[最终汇总](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/final_summary.json)与[最终证据](../../reference_tcad/simplemos_sentaurus2022/constants_srh_finite_20260908/final_evidence.json)关联全部冻结合同、原失败记录及本报告。
