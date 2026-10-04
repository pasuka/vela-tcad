# SimpleMOS 与已验证参考算例的差异及简化模型对照

2026-09-06 冻结并运行原生任务，2026-09-07 完成下载核验、Vela 重算及汇总。**模型简化对照已完成，主差异仍在：高 NWell 四点的电流误差从 +10.11%～+14.06% 变为 +9.81%～+13.70%，只下降 0.30～0.45 个百分点。** 原模型与开启线性修正的 Masetti 模型均为 8/8 点严格合格；关闭线性修正的 Masetti 模型为 0/8。本轮审计覆盖材料、物理、算法和参数，原始证据和失败记录保留。

**已有验证按器件、网格、偏置、指标和接受门槛成立。不能将当前差异笼统解释为“打开了从未验证的模型”，也不能将 BJT 或其他 MOS 的通过结论推广到当前整个组合。** [参考目录说明](../../reference_tcad/README.md)只对普通 CSV fixture 承诺符号、趋势、有限性和数量级；较强的结论需要下列专门合同与报告支持。

## 材料与物理模型

| 项目 | 当前 SimpleMOS 定位支路 | Genius BJT M1 | TransportModels MOS |
| --- | --- | --- | --- |
| 器件/网格 | Lg=0.25 μm；n19/n23 为 1480/1482 节点、2742/2746 三角形 | 6×2 μm 硅 NPN；5611 节点，局部加密 15561 节点 | 50 nm MOS；独立网格及掺杂 |
| 实际材料区域 | Si、SiO₂、两个 Nitride 侧墙区域；没有 PolySilicon 区域 | 纯 Si | Si、氧化层、多晶硅及氮化物材料合同 |
| Si 介电/带参数 | εr=11.7，Eg=1.12 eV，亲和势=4.0727403846 eV，300 K | 相同基础数值 | 相同基础数值 |
| Si 本征浓度 | matched-ni=1.0750038488844236e10 cm⁻³ | 1.4638914958767616e10 cm⁻³，结合 Fermi/BGN 使用 | 1.4638914958767616e10 cm⁻³，结合 Fermi/BGN 使用 |
| 统计/带隙变窄 | Boltzmann、no-BGN | Fermi–Dirac、OldSlotboom、Fermi 修正 | Fermi–Dirac、OldSlotboom、Fermi 修正 |
| 迁移率 | PhuMob＋Lombardi Enormal＋高场饱和，向量准费米梯度驱动 | Masetti 掺杂迁移率；无表面/高场退化 | Masetti＋Lombardi＋高场饱和，向量准费米梯度驱动 |
| 复合 | 掺杂相关 SRH；τe,max=1e-5 s，τh,max=3e-6 s | 同样的 300 K SRH 寿命参数，另开经典 `np−ni_eff²` Auger | SRH τe,max=3e-8 s，τh,max=3e-6 s |
| 其他功能 | Auger、雪崩、DG 量子势均已关闭 | Auger 开；无雪崩、DG | DD 与电子 DG 是两个独立分支 |

当前“NWell”变量实际施加背景 Boron：n19/n23 为 1e17/2e17 cm⁻³。n19/n23 的 Si 均为 1750 个三角形；氧化层和氮化物部分网格略有不同。不能把两器件的节点号直接当作全材料一一对应，也不能把这种高低配对称为同一网格上的纯参数导数。

当前材料文件保留 PolySilicon 和量子参数条目，**条目存在不代表模型或区域启用**。SiO₂/Nitride 的 εr 分别为 3.9/7.5。当前 Si 的 Nc/Nv=2.8e19/1.04e19 cm⁻³，BJT 为约 2.8567e19/3.1046e19 cm⁻³；当前 Boltzmann/no-BGN 支路采用固定 ni，因此本轮没有为模仿 BJT 而改动 ni/DOS 或开启 Fermi/BGN。

SRH 的共同 300 K 设置为 τmin=0、Nref=1e16 cm⁻³、γ=1、总杂质浓度；BJT 的温度依赖开、当前支路关，但在冻结的 300 K 数值相同。`unit_scaling` 下旧 JSON 字段虽然部分仍以 `_m3` 结尾，物理参数值按 cm⁻³解释；导出的状态文件 n/p 则为 m⁻³，不能混用。

## 算法及“已验证”的边界

| 项目 | 当前严格定位流程 | BJT M1 / 其他参考证据 | 本轮处理 |
| --- | --- | --- | --- |
| 求解与初值 | 原生一致初态的独立耦合 Newton 重闭合，contact_basin | BJT 使用 Gummel/Newton 交接、contact_majority 与连续延续；TransportModels 存在独立重闭合与连续扫描的不同合同 | 保留当前初值及坐标定义，明确不认证完整扫描 |
| 步长/范数 | max_iter=200，reltol=1e-7，abstol=1e-12；ψ 截断 0.35 V，准费米截断 0.025 V；标量范数线搜索 | BJT reltol=1e-8、abstol=1e-9、ψ=0.2 V、准费米=0.1 V；TransportModels 连续数值合同用 block_filter | 不复制其他器件的截断或放宽范数门槛 |
| 连续性行 | 所有 1814 个有效载流子行；eps=1e-6，无尺度地板、无小载流子/小源/小通量排除 | BJT eps=1e-3，附带密度和源/通量筛选；空间载流子比较又使用 ≥1e10 cm⁻³ 的有效区域，属于另一层指标 | 所有旧门槛和失败记录不变 |
| 行缩放 | flux_fraction=1e-3，max_weight=1e12 | BJT 为 0、1e18；TransportModels 连续合同为 1e-3、1e12 | 固定，不将缩放改善和物理简化混为一谈 |
| 事后资格 | global tol=1e-6、source_floor=1e-10；KCL/Id≤1e-8 | BJT global tol=1e-6、source_floor=1e-18；TransportModels 深关断 Id/|KCL|≥10 且 global tol=0.1 | 检查尺度不同，不能只比较同名 tolerance 数值 |
| 迁移率 Jacobian | 隔离程序启用完整向量链式修正及局部差分步长；另有固定 4 次线性迭代修正 | BJT Masetti 只依赖冻结掺杂/温度，没有载流子、界面场、高场驱动的状态链式项 | 简化支路关闭实验性向量/输运差分控制，另设线性修正开/关对照 |
| 电流/空间离散 | 守恒端口电流；原网格 Poisson/SG 及体积策略 | BJT cell-first 仅为重建显示诊断；PN2D 的 mixed-Voronoi/SG-GSS-Laux 是指定非钝角网格上的原子配置 | 不改端口算法，不套用 BJT 显示重建或 PN2D 专用默认 |

固定温度和掺杂下，简化后的 Masetti 迁移率仍随空间掺杂变化，**不是全器件常数**；但对 Newton 的 ψ、φn、φp 状态导数为零。此时关闭高场/界面/载流子散射链式修正有物理依据，不会因为人为丢掉一个仍非零的迁移率偏导而制造新的近似 Jacobian。本轮没有声称剩余全部 Jacobian 已获得任意分支的符号解析认证。

已有较强证据如下：

- [BJT 当前验收](../../reference_tcad/genius_bjt_sentaurus2022/CASE_SUMMARY.md)：M1 在 VCE=0.5–3 V 的 Ic/Ib/Ie/β 最大误差分别为 0.00344/0.00210/0.00341/0.00135 dex；SRH/Auger 与守恒截面通过。粗网格弱空穴节点电流恢复仍失败，局部加密 3 V 通过。简化 M0 不承担精度验收，所以“改成常迁移率就一定已有验证”不成立。
- [TransportModels DD 21 点](transportmodels_dd_contact_basin_v1_2026-08-24.md)：深关断/过渡/导通最大电流误差 2.409%/7.610%/2.930%；前 5 点连续，后 16 点独立重闭合。[空间复核](transportmodels_three_regime_spatial_2026-08-24.md)另有弱电流与场差异，不能将端口通过当作全场通过。Lombardi/HFS 已有这里的 MOS 证据，缺的是当前 PhuMob 组合、网格和严格弱行覆盖。
- [PN2D BV 当前合同](pn2d_bv_validation.md)：指定 M0/M2 非钝角网格的雪崩开启 RMSE 约 0.0018–0.0019 dex，正向最大误差 0.4066%。不支持将 mixed-Voronoi、源项映射等拆开移植为当前 MOS 默认。
- [Schottky 验收](../../reference_tcad/schottky_charon_sentaurus2018/README.md)：697 节点纯硅、热发射 Robin 接触，24 点最大误差 0.478652 dex，在其 0.5 dex 合同内通过；属于接触/延续功能验证，不是当前 MOS 的百分之一精度基线。其准费米截断必须关闭的经验也不能直接覆盖当前 0.025 V 设置。
- [SingleDevice 当前验收](../../reference_tcad/singledevice_sentaurus2018/singledevice_validation_20260817.json)：3584 节点/6972 三角形，原 deck 开启电子 DG、OldSlotboom、DopingDep/HFS/Enormal 和 SRH。Vd=0.1/1.1 V 的两条 Save/Load 自洽支路通过，Ion 相对误差约 0.654%/0.260%，低电流另用混合绝对/对数规则；三点空间场明确仅作诊断，KCL 规则为绝对残差≤1e-14 A/μm 或相对值≤1%。这也提供 HFS/Enormal 的 MOS 端口证据，但不是当前全行资格。
- [BVmethods](../../reference_tcad/bvmethods_sentaurus2018/bvmethods_nontransient_validation_20260817.md)的 IIC、外电阻和电压转电流已有对应范围的通过记录，NMOS continuation 尚 pending、瞬态不在该验收内；当前 0.8/1 V 直流对照未启用这些功能。[nmos2d 普通 fixture](../../reference_tcad/nmos2d_sentaurus2018/nmos2d_sentaurus2018_reference.json)允许 9 个数量级、无需趋势一致，不能把它作为精密 MOS 基线；PN2D 其他粗化/畸变网格目录也不能替代正式 M0/M2 合同。

## 本次实际对照定义

先进行一组有证据支持的复合简化，而非一次关闭所有数值稳定措施。两侧同时将 `PhuMob HighFieldSaturation Enormal` 换成 `DopingDependence`（Vela `masetti`，总杂质浓度）。保留网格、掺杂、接触、300 K、matched-ni/Boltzmann/no-BGN 和原 SRH。它能检验原迁移率组合是否为差异的必要原因；不能单凭这一组区分 PhuMob、Lombardi 和 HFS 各自贡献。

工作点为 n19/n23 × Vd=0.05/1 V × Vg=0.8/1 V，共 8 点。Sentaurus 原模型和简化模型均重新计算，共 8 份 deck、16 个 DC 状态，同用 `ExtendedPrecision(128) Method=Super Digits=12 ErrRef(e/h)=1e-2 RhsMin=1e-20 Iterations=40`，精确加载对应偏置保存态后重新耦合求解。

Vela 共三组、24 次独立自洽计算：原模型＋现有向量链式修正＋4 次线性修正；Masetti＋实验性向量/输运差分控制关闭＋线性修正关闭；同一 Masetti 初态及配置＋4 次线性修正。两组 Masetti 只差线性修正开关。保留原生 ψ/电子和空穴准费米势，以 Vela 的冻结 ni/热电压重算初始密度，随后完整自洽；不是把原生密度直接代入后只算一次电流。

Sentaurus T-2022.03-SP2 运行目录为 `/tmp/vela_simplemos_reference_subset_20260906`。Vela 复用隔离 Release 程序、Eigen SparseLU，4 次修正使用 100 位残差累加而矩阵/状态仍为 double。TDR 通过现有 HDF5 导入程序导出。电流为 A/μm、宽度 1 μm。本轮无生产源代码或默认值修改。

## 完成结果

16 个原生状态均退出成功、偏置精确匹配，KCL/Id 最大约 3.66e-15，输入文件返回哈希一致。Vela 的原模型 8 点和 Masetti＋修正 8 点全部满足 Newton 成功、1814 行 eps=1e-6、零尺度排除数为 0、原 global 检查及 KCL/Id≤1e-8。它们是同模型、同偏置的独立自洽对比。

下表为有符号误差 `100×(Id_Vela/Id_Sentaurus−1)`；所有列均来自上述合格状态：

| 器件 | Vd / Vg，V | 原模型误差 | Masetti 简化误差 | 误差下降，百分点 |
| --- | --- | ---: | ---: | ---: |
| n19 | 0.05 / 0.8 | +7.7159% | +7.5587% | 0.1571 |
| n19 | 0.05 / 1.0 | +2.5068% | +2.2231% | 0.2837 |
| n19 | 1.0 / 0.8 | +4.1404% | +4.0900% | 0.0503 |
| n19 | 1.0 / 1.0 | +2.0012% | +1.4381% | 0.5632 |
| n23 | 0.05 / 0.8 | +14.0642% | +13.6995% | 0.3647 |
| n23 | 0.05 / 1.0 | +13.4315% | +12.9842% | 0.4474 |
| n23 | 1.0 / 0.8 | +12.2626% | +11.8812% | 0.3814 |
| n23 | 1.0 / 1.0 | +10.1108% | +9.8124% | 0.2984 |

简化不是无效开关：例如 n23、Vd=0.05 V、Vg=1 V 的原生 Id 从 1.50105654e-7 增至 2.08718283e-7 A/μm；Vela 从 1.70267163e-7 增至 2.35818667e-7 A/μm。两侧电流都明显变化，但相对误差只由 13.4315% 降到 12.9842%。因此在这 8 个点上，PhuMob/Lombardi/HFS 组合不是产生大部分偏差的必要条件；此实验不能量化各子模型的独立贡献，也不能外推到整个 0–1 V。

高低 NWell 配对误差定义为 `log10(Id23/Id19)_Vela−log10(Id23/Id19)_Sentaurus`，每种模型分别使用自己的原生基线：

| Vd / Vg，V | 原模型，dex | 简化模型，dex | 变化 |
| --- | ---: | ---: | --- |
| 0.05 / 0.8 | +0.02486968 | +0.02411296 | 改善 |
| 0.05 / 1.0 | +0.04398136 | +0.04346868 | 改善 |
| 1.0 / 0.8 | +0.03261596 | +0.03134788 | 改善 |
| 1.0 / 1.0 | +0.03322474 | +0.03445059 | 恶化 |

所以“每个器件的绝对误差略降”并不等于“高低配对一致改善”；高漏压反例仍须保留，本次不放行 M82/M83。

**数值开关对照。** Masetti 两组初态哈希相同，配置除输出位置外相同；实验性向量/输运差分控制都关闭。关闭线性修正的 8 态均以 `carrier_row_convergence_line_search_rejected` 退出：仍有 6～237 个超差行，全部为空穴行，最大行比值为 2.48e-4～1.21。它们虽全部满足事后 global/KCL，仍然不合格。开启固定 4 次修正后，8 态经 3～5 次接受迭代全部通过，端口电流相对关闭修正的结果最大变化仅 1.37e-12。

这一对照进一步区分了两件事：线性修正对严格少数载流子闭合很重要，但不解释当前约 10% 的端口差；反过来，端口电流看起来稳定也不足以接受少数载流子状态。此处没有完成多初值不变性或完整扫描回归，4 次修正仍是隔离诊断功能。

**附加导出语义审计。** 简化模型的原生节点 e/hMobility 在每器件四个偏置上逐值完全相同，符合固定掺杂/温度的模型设置。但“将节点总掺杂代入 C++ 默认 Masetti 公式”与“原生节点绘图迁移率”不是已经校准的等价量：n19/n23 电子最大差分别为 48.17%/45.16%，空穴为 39.11%/36.74%。例如 n23 节点 982 的电子局部公式值为 131.8192、节点导出值为 240.3709 cm²/(V·s)。这些只作为字段语义诊断保留；原生内部边/单元计算及到节点的投影尚未对齐，不能据此宣称 Masetti 解析式错误或把该百分比当作 Id 误差贡献。

## 后续定位与证据

本次简化模型可作为后续局部离散审计的更简单基线：先在已合格的 Masetti 状态上校准原生边/单元迁移率与节点导出的关系，随后对齐 SG 电导、Poisson 电荷体积及 Si/绝缘体界面势形成，沿既有校准的守恒扰动方向验证电流响应。这样可在没有场依赖迁移率链式项的条件下检查剩余耦合差异。此次没有更改电流算法、材料参数、阈值或继续启动这些新增扰动仿真。

完成的检查包括 16 个原生状态返回验证及 HDF5 字段导出、24 次 Vela 自洽和 24 次全行/旧 global 事后检查；既有 Masetti Catch2 3 项测试、10 条断言通过。新增脚本语法、报告链接、冻结输入及最终状态/结果哈希另由封存脚本核验。无新增 C++ 核心改动，因此未运行全量 CTest。工作树已有 `vela_example_runner.cpp` 的 127 行新增/1 行删除保持不变。

- [原生合同](../../reference_tcad/simplemos_sentaurus2022/reference_subset_20260906/native_contract.json)、[Vela 合同](../../reference_tcad/simplemos_sentaurus2022/reference_subset_20260906/vela_contract.json)、[结果摘要](../../reference_tcad/simplemos_sentaurus2022/reference_subset_20260906/summary.json)。
- [原生电流](../../reference_tcad/simplemos_sentaurus2022/reference_subset_20260906/native_points.csv)、[24 态结果及失败记录](../../reference_tcad/simplemos_sentaurus2022/reference_subset_20260906/vela_runs.csv)。
- [模型差异](../../reference_tcad/simplemos_sentaurus2022/reference_subset_20260906/model_comparison.csv)、[NWell 配对](../../reference_tcad/simplemos_sentaurus2022/reference_subset_20260906/nwell_pairs.csv)、[线性修正开关对照](../../reference_tcad/simplemos_sentaurus2022/reference_subset_20260906/numerical_comparison.csv)。
- [节点迁移率诊断](../../reference_tcad/simplemos_sentaurus2022/reference_subset_20260906/mobility_nodal_diagnostic.csv)、[语义边界及公式参数](../../reference_tcad/simplemos_sentaurus2022/reference_subset_20260906/mobility_nodal_diagnostic.json)。
- [最终封存](../../reference_tcad/simplemos_sentaurus2022/reference_subset_20260906/validation_evidence.json)。原始 TDR/SAV、程序、状态和日志保留于忽略目录 `build-release/simplemos_reference_subset_20260906/` 及前轮隔离程序目录。
