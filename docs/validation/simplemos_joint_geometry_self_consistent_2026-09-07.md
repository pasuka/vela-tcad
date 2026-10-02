# SimpleMOS 硅侧输运几何与 Poisson 体积联合自洽验证

2026-09-07。**联合修改在四个 Vg=1 V 工作点均通过冻结的比较门槛。高 NWell 的 Id 相对误差由 +12.9842%/+9.8124% 降至 +0.02581%/+0.07311%；低 NWell 的误差绝对值也减小。** 新增 76 次自洽计算全部合格，24 个小扰动幅度结果、12 组双初始化对照均通过。

本轮还在前置检查中发现了原基线接近平衡空穴电流的 Jacobian 求值消减问题。已在独立验证程序中采用稳定解析等价式，保持物理残差及全部收敛门槛不变后，再完成几何对照。原失败记录保留，正式求解器源文件未改写。

本轮承接[原生单元迁移率、边通量与分布源校准](simplemos_native_mobility_flux_calibration_2026-09-07.md)。结论限于 Masetti 简化模型、n19/n23 × Vd=0.05/1 V、Vg=1 V 四点；不能据此宣称整个 0–1 V 曲线或全部物理场已验证通过。M82/M83 尚未放行。

## 1. 完整替换结果

以下为实际 α=1 自洽重算结果，从原合格 Vela 状态出发；不是将小扰动导数线性外推至 α=1。误差统一定义为 `100×(Id_Vela/Id_Sentaurus−1)`。

| 工况，Vg=1 V | 基线 | 仅硅侧输运几何 T | 仅 Poisson 体积 P | 联合 T+P |
| --- | ---: | ---: | ---: | ---: |
| n19，Vd=0.05 V | +2.22309% | −6.77419% | +8.95850% | **−0.76256%** |
| n19，Vd=1 V | +1.43805% | −6.30758% | +7.30357% | **−1.01375%** |
| n23，Vd=0.05 V | +12.98419% | −0.30986% | +13.43616% | **+0.02581%** |
| n23，Vd=1 V | +9.81243% | −2.25670% | +12.50076% | **+0.07311%** |

联合修改后的四点 Id 分别为 6.43592767161e-6、3.80090469196e-5、2.08772160107e-7、1.24198948570e-6 A/μm。对应原生值为 6.48538281417e-6、3.83983089058e-5、2.08718282997e-7、1.24108214100e-6 A/μm。

单独使用 T 会改善高 NWell，却使低 NWell 控制变差；单独使用 P 会使四点 Id 误差增大。因此不能仅凭一个高 NWell 点，或者只看 Poisson 场更接近原生，就决定启用单项修改。

高低配对定义为 `log10(Id_n23/Id_n19)−log10(Id_native_n23/Id_native_n19)`：

| Vd | 基线配对差 | 联合修改配对差 |
| --- | ---: | ---: |
| 0.05 V | +0.0434687 dex | +0.00343655 dex |
| 1 V | +0.0344506 dex | +0.00474251 dex |

联合方案同时满足运行前冻结的三个比较条件：高 NWell 绝对相对误差降低、高低配对差减小、低 NWell 误差不增大。两个单项方案均未满足全部条件。这里“绝对相对误差”指上述带符号误差的绝对值。

证据：[完整替换和双初始化电流](../../reference_tcad/simplemos_sentaurus2022/joint_stable_20260907/replacement_comparison.csv)、[高低配对](../../reference_tcad/simplemos_sentaurus2022/joint_stable_20260907/NWell_pairing.csv)、[四点比较门槛结果](../../reference_tcad/simplemos_sentaurus2022/joint_stable_20260907/promotion.csv)。

## 2. 本轮实际修改的离散量

采用上一轮独立导出并校准的原生 box 系数和 Si 节点体积，没有重新拟合参数，也没有新增原生仿真。

- **T：输运几何。** 每条边的电子、空穴通量均使用原生相邻 Si 单元系数之和。端口提取器同步使用同一几何。Masetti 边迁移率保持 Vela 原定义，没有在本轮替换成原生 box 加权单元迁移率。
- **P：Poisson 电荷体积。** 将电子、空穴和净掺杂三个电荷项的体积同时改为原生 Si 节点体积；其对应 Jacobian 项同步使用该体积。没有只修改电子项。
- **T+P：联合。** 同时执行上述完整修改。Poisson 介电矩阵 K、连续性/SRH 源体积、材料参数、边界和物理常数均保持原基线定义。

系数按 `baseline×[1+α(native/baseline−1)]` 插值。α=0 为零控制，α=1 为完整替换；原生有效系数缺失或无效时不会凭空补边。完整系数计算在每次 Newton 残差、Jacobian 和端口求值中生效，不是固定状态电流后处理，也不是上一轮的常量连续性源注入。

原有 `transportEdgeCoupling` 等开关的存在并不保证提供原生已处理的 box 体积。因此本轮在忽略目录中构建带有外部系数比值输入的隔离 runner，直接使用已核对的数据；没有把实验开关变成全局默认。

证据：[原始实验合同](../../reference_tcad/simplemos_sentaurus2022/joint_geometry_20260907/contract.json)、[稳定数值版本合同](../../reference_tcad/simplemos_sentaurus2022/joint_stable_20260907/contract.json)、[实际几何比值](../../reference_tcad/simplemos_sentaurus2022/joint_stable_20260907/geometry.csv)、[隔离几何实现](../../scripts/diagnostics/simplemos_native_geometry.hpp)。

## 3. Jacobian 前置缺口及独立处理

首版固定状态预检确认了基线残差、边通量逐位不变，也确认了 T/P 残差差分及端口一致性。但在所有四工况、四个几何变体上，空穴连续性行对 ψ 的 Jv 块与中心差分严重不符：两个较小步长的相对误差约为 1，共保留 32 个失败结果，因而没有启动该版本计划中的完整替换 DC。

该问题位于原 Masetti 分支的解析导数求值。其 ψ 导数由含大密度的漂移、扩散项相减；接近平衡时真实导数很小，double 舍入会主导结果。它是本次换几何前就存在的数值问题，不能归因于新几何。在 n19、低 Vd 的基线 ψ→空穴方向中，原 Jv 范数约 5.76e-8，而残差差分约 2.91e-13，单位沿用该探针的缩放残差约定。

独立数值版本保留同一物理 SG 通量 F，改写为稳定解析等价式。令 η=(ψ₁−ψ₀)/Vt：

\[
\partial_{\psi_0}F_n=\frac{F_n}{V_t}\left[1+\frac{B'(-\eta)}{B(-\eta)}\right],\qquad
\partial_{\psi_1}F_n=-\frac{F_n}{V_t}\frac{B'(-\eta)}{B(-\eta)},
\]

\[
\partial_{\psi_0}F_p=-\frac{F_p}{V_t}\left[1+\frac{B'(\eta)}{B(\eta)}\right],\qquad
\partial_{\psi_1}F_p=\frac{F_p}{V_t}\frac{B'(\eta)}{B(\eta)}.
\]

F 使用原有稳定 SG 求值和同一准费米参考/增量约定。此修改限定于当前 equal-ni、Boltzmann/no-BGN、无状态相关迁移率的工况；没有改变准费米导数、SRH、物理残差或线搜索。本轮未认证指数截断激活状态、变 ni、Fermi 或场相关迁移率下的适用性。

独立 Catch2 测试以 100 位精度的**密度形式 SG 通量中心差分**作为参照，覆盖正负势、零及 1e-18/1e-10/0.03 V 准费米差，32 项断言通过。原几何性质测试另有 3 个测试、9 项断言通过。

随后重新执行原门槛预检：对三个输入块使用全自由节点正弦方向，步长 1e-5、5e-6、2.5e-6 V，分输出块检验，不使用 floor=1 的归一化。七个可能变化的块在两个较小步长上最大相对差 5.56e-9，原失败的 ψ→空穴块最大 5.49e-9，均小于冻结门槛 1e-4。电子—空穴、空穴—电子 SRH 交叉块在同状态几何变体间的解析值逐位相同，其弱差分噪声保留，未据此重复宣称任意弱列已全部通过。

28 个基线/完整替换最终状态也执行了 Jv 检查，七块最大相对差为 6.21e-9。稳定数值版本的零控制 Id 相对原合格基线最大漂移约 1.37e-12，固定状态物理残差与边通量逐位相同。因此该数值修正没有造成后文约 10% 的电流变化，但解决了本次前置 Jacobian 资格缺口。

证据：[首版保留失败](../../reference_tcad/simplemos_sentaurus2022/joint_geometry_20260907/preflight_jvp.csv)、[数值版本补充合同](../../reference_tcad/simplemos_sentaurus2022/joint_stable_20260907/numerical_addendum.json)、[稳定导数实现](../../scripts/diagnostics/simplemos_stable_sg_psi_derivative.hpp)、[100 位独立测试](../../tests/diagnostics/test_simplemos_stable_sg_derivative.cpp)、[稳定版预检](../../reference_tcad/simplemos_sentaurus2022/joint_stable_20260907/preflight_jvp.csv)、[最终状态 Jv](../../reference_tcad/simplemos_sentaurus2022/joint_stable_20260907/replacement_jvp.csv)。

## 4. 小扰动、完整替换及初始化资格

每工况测试共同零控制；对 T、P、T+P 各测试 α=±0.001、±0.0005，并在 α=1 时分别从原合格 Vela 状态、原生电势/准费米势配合 Vela 一致密度的状态出发。合计每工况 19 次、四工况 76 次 DC。

固定状态残差差分先与独立构造的保守边通量差及 Poisson 电荷体积差核对。预检的载流子误差界明确包含 `64×double epsilon×基线绝对边流量和`，用于限定残差相减的舍入误差；它不是非线性接受门槛，也没有覆盖上一轮保留的不可分辨小源失败记录。

小扰动采用原完整 DD 权重的 `−λᵀδR_free+δId_direct` 预测，接触替换行剔除、直接电流项单列。通过新数值版本的实际正负重算验证这些预测，而非假定旧权重自然适用于修改后矩阵。

| 项目 | 本轮结果 | 冻结门槛 |
| --- | ---: | ---: |
| 合格自洽状态 | 76/76 | 每态全部资格通过 |
| 有效载流子行独立复核 | 137,864 行 | 每态 1814 行，无零尺度排除 |
| 最大载流子行比值 | 9.42203e-7 | ≤1e-6 |
| 最大 KCL/Id | 1.993e-14 | ≤1e-8 |
| 合格正负幅度结果 | 24/24 | 下列响应检查全部通过 |
| 最大预测相对误差 | 1.087e-8 | ≤0.001 |
| 最大双幅度线性差 | 6.510e-9 | ≤0.001 |
| 最大偶分量/奇分量 | 8.597e-5 | ≤0.01 |
| 最小信号/零控制漂移 | 1.0109e7 | ≥100 |
| 双初始化合格配对 | 12/12 | 电流、势及密度均通过 |
| 最大双初始化势差 | 1.727e-8 V | ≤1e-6 V |
| 最大双初始化密度相对差 | 6.680e-7 | ≤1e-4 |
| 最大双初始化 Id 相对差 | 2.662e-10 | ≤1e-6 |

最大行比值出现在 n19、Vd=1 V 的联合方案原生初始化终态，已接近 1e-6 门槛；仍按原门槛记录通过，没有扩大容差或排除该行。所有事后 global 检查保持原 source_floor 语义并通过。

证据：[76 次原始 DC 摘要](../../reference_tcad/simplemos_sentaurus2022/joint_stable_20260907/dc.csv)、[独立逐行资格](../../reference_tcad/simplemos_sentaurus2022/joint_stable_20260907/independent_qualification.csv)、[24 个响应结果](../../reference_tcad/simplemos_sentaurus2022/joint_stable_20260907/response_calibration.csv)、[双初始化对照](../../reference_tcad/simplemos_sentaurus2022/joint_stable_20260907/initialization_invariance.csv)。

## 5. 电流接近不等于全部物理场接近

下表比较高 NWell 基线与联合方案相对原生电势/准费米势的差。RMS 按自由 Si 节点的 Si box 体积加权，max 为相同 907 节点支持上的最大绝对差；单位均为 mV。原生准费米势与 Vela 物理参考+增量在同一约定下比较，未对差值另作任意常数平移。

| Vd | 物理量 | 基线 RMS → 联合 RMS | 基线 max → 联合 max |
| --- | --- | ---: | ---: |
| 0.05 V | ψ | 1.35251 → 0.006516 | 13.80745 → 0.73105 |
| 0.05 V | φn | 0.041380 → 0.002053 | 0.43880 → 0.01546 |
| 0.05 V | φp | 0.034782 → 0.034818 | 0.43837 → 0.43818 |
| 1 V | ψ | 1.65063 → 0.074946 | 13.80758 → 1.73311 |
| 1 V | φn | 1.34070 → 0.772323 | 28.42217 → 25.74067 |
| 1 V | φp | 0.974237 → 0.779912 | 22.96270 → 20.44617 |

高漏压下，界面节点 338 的 ψ 差由 −1.37661 mV 降至 +0.06018 mV，φn 差由 +0.61663 mV 降至 +0.04444 mV；沟道节点 320、324 的结果也已记录。相比之下，全自由 Si 支持上的少数载流子准费米势最大差仍有约 20–26 mV，不能用很小的 Id 误差替代全场一致性验证。低漏压 φp 的 RMS 也没有改善。

本轮的证据支持：在这四个合格简化模型状态上，硅侧输运几何与 Poisson 电荷体积的联合作用解释并消除了大部分高 NWell 电流差，同时明显改善主要电势差。它不证明这是所有偏置、所有模型的唯一根因；本轮未替换的单元迁移率平均、介电 K、SRH 源体积及约定差仍可能影响剩余场差。

证据：[四方案物理场对照及局部节点](../../reference_tcad/simplemos_sentaurus2022/joint_stable_20260907/field_comparison.csv)。本轮没有新增生成复合项的逐场等价性认证。

## 6. 环境、复现与下一阶段

沿用原网格和模型：n19/n23 背景 Boron=1e17/2e17 cm⁻³；节点 1480/1482、三角形 2742/2746；Si 各 1750 个三角形、942 个节点，另含 SiO₂ 和两个 Nitride 侧墙。300 K、matched-ni=1.0750038488844236e10 cm⁻³、Boltzmann/no-BGN、Masetti 总杂质与掺杂相关 SRH。未启用 HFS、表面迁移率、Auger、雪崩或 DG。

使用 Windows UCRT64 C++20、隔离 Release runner、Eigen SparseLU、既有 4 次线性迭代修正（100 位残差累加；矩阵和状态为 double）。max_iter=200、reltol=1e-7、abstol=1e-12，ψ/准费米更新上限 0.35/0.025 V，contact_basin 及原标量线搜索保持不变。逐行 eps=1e-6，不按小源或少数载流子排除有效行。求解阶段 global 仍为 off，事后探针以 enforce、tolerance=1e-6、source_floor=1e-10 检查，不将其解释为任意微小净 SRH 源都达到相对 1e-6 闭合。

电流单位 A/μm、宽度 1 μm；网格坐标 μm，Vela 状态载流子 m⁻³，unit_scaling 下配置浓度按 cm⁻³ 解释。原生参照为已封存的 Sentaurus T-2022.03-SP2 结果及成功导入的 TDR/box 数据。本轮未向虚拟机新增上传或启动 sdevice。

新增内容为隔离诊断实现、Catch2 性质测试、验证脚本与报告。正式核心源码未改；原工作树 runner 的 127 行新增/1 行删除保持不变。两个隔离程序编译成功，4 个 Catch2 测试共 41 项断言通过。稳定版本首次编译的头文件接口错误已修正，仅重编失败单元后链接；失败编译日志保留。未运行全量 CTest。

下一阶段应先把同一联合方案扩展到 Vg=0.8 V 的四个既有控制点，形成完整八点对照，并检查门槛接近状态及剩余准费米场差。通过后，再核对生产实现中输运、Poisson、端口和诊断支路的一致性，扩展 16 工况及完整 0–1 V 曲线。当前没有将隔离方案设为默认，也没有以四点通过替代 M82/M83 的全部要求。

- [原隔离程序构建](../../scripts/build_simplemos_joint_geometry_20260907.py)、[稳定数值版本构建](../../scripts/build_simplemos_joint_stable_20260907.py)、[76 次执行入口](../../scripts/validate_simplemos_joint_stable_20260907.py)。
- [固定状态独立预检](../../scripts/check_simplemos_joint_geometry_20260907.py)、[独立结果分析与封存](../../scripts/analyze_simplemos_joint_stable_20260907.py)。
- [几何性质测试](../../tests/diagnostics/test_simplemos_native_geometry.cpp)、[分析指标封存](../../reference_tcad/simplemos_sentaurus2022/joint_stable_20260907/analysis_evidence.json)、[最终输入和输出哈希](../../reference_tcad/simplemos_sentaurus2022/joint_stable_20260907/validation_evidence.json)。

原始数据、源副本、可执行程序及日志位于忽略目录 `build-release/simplemos_joint_geometry_20260907/` 和 `build-release/simplemos_joint_stable_20260907/`。合同、失败与成功结果位于对应 `reference_tcad/simplemos_sentaurus2022/` 子目录。所有已冻结证据按只追加方式保存，未覆盖上一轮或首版失败数据。
