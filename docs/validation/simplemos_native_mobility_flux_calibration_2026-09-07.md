# SimpleMOS 原生单元迁移率、守恒边通量与分布源校准

2026-09-07。**本轮已确认原生单元迁移率的 box 加权定义，复现原生守恒 SG 通量，并完成组合输运、几何、迁移率三个电子分布源的独立自洽响应校准。** 原生补充导出最终 4/4 成功；Vela 新增 60/60 个自洽状态、24/24 个双幅度响应结果通过原有门槛。独立审计重新计算了 108,840 个载流子行。

主要定位进展是：高响应的 Si/SiO₂ 界面边上，当前 Vela 输运几何含两侧材料的贡献，原生载流子通量使用硅侧贡献；这些边的硅侧几何约为 Vela 值的 20.064%。对所定义的分布源，几何差异主导降低 Id 的方向，迁移率取值差异方向相反、影响较小。这已成为有校准依据的局部离散差异，但尚未证明它单独解释全部高 NWell 电流误差。

本轮承接[上一轮局部离散审计](simplemos_masetti_local_discretization_2026-09-07.md)。范围为 n19/n23 × Vd=0.05/1 V、Vg=1 V 四工作点；没有新增整个 0–1 V 扫描，也没有修改正式电流算法或放行 M82/M83。

## 1. 工况、单位与资格保持不变

n19/n23 的背景 Boron 为 1e17/2e17 cm⁻³；网格节点 1480/1482、三角形 2742/2746，Si 均为 1750 个三角形、942 个节点。材料为 Si、SiO₂ 和两个 Nitride 侧墙。每态核对 907 个自由 Si 节点，共 1814 个电子/空穴行。

物理模型为 300 K、Boltzmann/no-BGN、matched-ni=1.0750038488844236e10 cm⁻³、Masetti 总杂质浓度及掺杂相关 SRH；未启用 HFS、表面迁移率、Auger、雪崩、DG。上一轮已确认的 18 个 Masetti 参数保持不变。此模型下迁移率随冻结掺杂在空间变化，对当前 Newton 的 ψ、φn、φp 状态导数为零。

Vela 使用隔离 Release runner、Eigen SparseLU 和既有 4 次线性迭代修正，线性残差以 100 位精度累加，矩阵和状态仍为 double。保留 max_iter=200、reltol=1e-7、abstol=1e-12、ψ/准费米更新上限 0.35/0.025 V、contact_basin 参考及原标量线搜索。逐行 eps_row=1e-6，不按零尺度、小源或少数载流子排除有效行。求解配置的 global 检查仍为 off，**事后探针**按原有 enforce、tolerance=1e-6、source_floor=1e-10 检查；没有把它改写成任意微小净 SRH 源都达到相对 1e-6 闭合。

端口电流为 A/μm、宽度 1 μm；网格坐标 μm，Vela 状态载流子 m⁻³，unit_scaling 下材料/配置浓度按 cm⁻³ 解释。原生为 Sentaurus T-2022.03-SP2、128 位 Super；实际 TDR/HDF5 读入和原生状态重闭合成功。补充任务使用独立目录 `/tmp/vela_simplemos_masetti_runtime_20260907`，输入、结果和失败记录均已取回。

## 2. 原生单元迁移率定义已确认

对半导体单元 T，原生电子、空穴迁移率均满足：

\[
\mu_T=\frac{\sum_{i\in T}M_{T,i}\,\mu_{\rm Masetti}(N_{D,i}+N_{A,i})}
                  {\sum_{i\in T}M_{T,i}}.
\]

其中 M 是**原生处理后的 element-vertex box measure**。分母必须使用该单元实际导出分量之和；不能换成坐标三角形面积，也不能换成未修正的几何分量。原生对部分单元的处理会使两种分母不同。上一轮算术平均、μ(平均掺杂)等失败候选与原始数据均保留。

| 核对项 | 数量或范围 | 最大差 |
| --- | --- | ---: |
| box 加权公式对原生单元 μ | 2 网格 × 2 载流子 × 1750 单元 | 相对 4.44e-16 |
| 本轮运行时 μ 对既有 TDR 单元 μ | 4 工况，电子/空穴 | 逐值完全相同 |
| 运行时坐标映射回原网格 | 含材料界面重复顶点 | 1.11e-16 μm |
| 运行时 element-edge 系数对原生 debug 导出 | 四工况硅单元 | 绝对 5.68e-14 |
| 运行时累计 Si 节点体积对原生 debug 导出 | 四工况 | 3.16e-18 μm² |

运行时拓扑含 1562/1564 个内部顶点，通过坐标与 Si 单元顶点集合映射回原始网格；没有假设内部序号等于 TDR node_id。四次成功补充运行的 Id 与冻结原生值完全一致，KCL/Id 最大 1.55e-15。

证据：[单元恒等式](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/box_mobility_identity.csv)、[运行时检查](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/geometry_runtime/checks.csv)、[拓扑映射](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/geometry_runtime/cell_mapping.csv)。

## 3. 原生保守 SG 边通量校准通过

令 g(T,e) 为原生单元对边的 box 系数，定义

\[
W_{n,e}=\sum_{T\in Si,\,T\supset e}g(T,e)\mu_{n,T},\qquad
W_{p,e}=\sum_{T\in Si,\,T\supset e}g(T,e)\mu_{p,T}.
\]

W 单位为 cm²/(V·s)，g 无量纲。取原生 ψ、n、p，η=(ψ₁−ψ₀)/Vt、B(x)=x/(exp(x)−1)：

\[
F_{n,e}=W_{n,e}V_t[n_0B(-\eta)-n_1B(\eta)],\qquad
F_{p,e}=W_{p,e}V_t[p_0B(\eta)-p_1B(-\eta)].
\]

这里 F 为每单位厚度的粒子线通量，单位 1/(cm·s)。端口约定 s(e,C)=1(node₀∈C)−1(node₁∈C)，则

\[
I_C=-q_{native}\,10^{-4}\sum_e s(e,C)(F_{n,e}-F_{p,e})\quad[\mathrm{A/\mu m}].
\]

本轮使用原生常数 q=1.602192e-19 C、kB=1.380662e-23 J/K；没有调整参数来拟合电流。另用准费米势等价式，以 60 位 Decimal 计算小势差下的通量，避免漂移/扩散大项相减。两种形式独立复现四工况端口电流。

| Vg=1 V 工况 | 密度形式 Id 相对原生差 | 准费米形式 Id 相对原生差 |
| --- | ---: | ---: |
| n19，Vd=0.05 V | +9.73e-10 | −1.86e-10 |
| n19，Vd=1 V | +3.06e-10 | −4.46e-10 |
| n23，Vd=0.05 V | −9.98e-9 | −1.82e-9 |
| n23，Vd=1 V | −6.49e-9 | −8.44e-9 |

八个端口复算全部通过运行前冻结的相对 1e-6 门槛，复算 KCL/Id 均小于 1e-8。将边通量散度与原生 SRH、Si 节点体积一起重建自由 Si 连续性表达式，以漏端粒子通量归一化，准费米形式的电子残差 L1 最大为 4.06e-8，空穴为 4.14e-24。后一个数值与电子端口使用同一归一化量，**不能解释成空穴自身相对精度达到 24 位**；密度形式的少数载流子残差受消减影响更明显，原值保留。这是原生表达式复算，不替代原非线性逐行接受检查。

本轮没有将绘图 `eCurrentDensity/Element` 直接当作保守边通量，也没有将 Tcl `ReadFlux` 当作载流子通量。所校准的是上述原生单元系数、迁移率和 SG 组合，通过端口求和与自由行连续性获得支持。

证据：[端口复算](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/native_sg_port_checks.csv)、[连续性复算](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/native_continuity_replay.csv)、[公式及单位实现](../../scripts/audit_simplemos_masetti_box_mobility_20260907.py)。

## 4. 三个定义明确的分布源均通过自洽校准

在每个合格 Vela 基态上，用相同 Vela 状态、常数和 SG 公式，构造

`δF_combined = F_Vela × (W_native / (g_Vela μ_Vela) − 1)`。

先统一 cm/m 单位。这个差分只比较输运系数；原生热电压、密度约定和 Poisson 修正没有混进该源。δF 在每次 Newton 求解中保持为常量，通过 +node₀/−node₁ 注入连续性行，同时按相同方向更新接触电流。它的状态导数为零。使用现有完整 DD 端口伴随计算

`dId/dα = −λᵀ δR_free + δId_direct`，

剔除被 Dirichlet 替换的接触残差行，并明确保留接触电流直接项。

另外冻结两条可加方向：几何源将 g_Vela 改为原生 Si 的 g、保持 Vela μ；迁移率源在该 Si 几何下比较 Vela μ 与原生系数加权单元 μ。两源逐边相加等于组合源，最大浮点相对余项 2.54e-16；端口线性预测之和的相对余项不超过 4.45e-16。这是一种明确的有序分解，迁移率项并非“在原 Vela 全材料几何上单独换 μ”的结果。

每条方向都测试 α=0、±0.001、±0.0005，四工况各五次自洽重算。**没有运行 α=1 的系数替换。** 下表列的是 α=0 处的归一化导数 `100×(dId/dα)/Id_Vela`，不是实际误差改善百分比。

| Vg=1 V 工况 | 几何方向 | 迁移率方向 | 组合方向 |
| --- | ---: | ---: | ---: |
| n19，Vd=0.05 V | −8.77215%/α | +1.01084%/α | −7.76131%/α |
| n19，Vd=1 V | −7.60753%/α | +1.42244%/α | −6.18509%/α |
| n23，Vd=0.05 V | −11.76173%/α | +0.117025%/α | −11.64470%/α |
| n23，Vd=1 V | −10.98130%/α | +0.189705%/α | −10.79160%/α |

例如 n23、低 Vd 的组合源在 α=0.001 时，奇分量约为基线 Id 的 −0.0116447%，不能据此写成“电流误差已降低 11.64%”。组合响应的接触直接项占总导数约 0.62%～1.16%，其余主要来自自洽状态反馈。

| 验证组 | 合格 DC / 幅度结果 | 最大预测相对误差 | 最大双幅度线性差 |
| --- | --- | ---: | ---: |
| 组合输运源 | 20/20；8/8 | 6.76e-11 | 5.91e-11 |
| 几何源 | 20/20；8/8 | 5.78e-11 | 6.70e-11 |
| 迁移率源 | 20/20；8/8 | 2.28e-9 | 2.83e-9 |

60 态最大载流子行比值为 8.56e-9，小于 1e-6；最大 KCL/Id 为 1.02e-14，小于 1e-8；无零尺度排除。所有事后 global 检查通过。24 个幅度结果最大偶/奇分量为 1.29e-6，最小信号/零漂移约 3.70e6。预测误差≤0.1%、双幅度差≤0.1%、偶/奇≤1%、符号、零漂移≤1e-5 dex、信噪比≥100 的原门槛全部保持。

证据：[组合合同](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/distributed_response/contract.json)、[几何合同](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/factor_resolved_response/geometry/contract.json)、[迁移率合同](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/factor_resolved_response/mobility/contract.json)、[独立资格复核](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/independent_qualification.csv)、[独立响应复核](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/independent_calibration.csv)、[分解导数](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/factor_summary.csv)。

## 5. 主要响应来自界面硅侧几何

组合源的最大响应边位于 x≈0.00504594 μm 的 Si/SiO₂ 界面。n23、Vd=0.05 V 的最大项为边 1014（342–343），占总导数约 30.96%；Vd=1 V 的最大项为边 1007（338–342），占约 34.64%。这些是完整分布方向下的伴随贡献排序，尚未逐边独立做有限差分。

以 n23 的边 1014 为例，Vela 几何约 0.01315389；原生硅侧约 0.002639179、氧化层侧约 0.01051471，两者之和复现 Vela 几何。硅侧占 20.064%，其余约 79.936% 来自氧化层侧。边 1007 也呈现相同几何比例。合并迁移率后，n23 两条对应主导边的 `W_native/(g_Vela μ_Vela)` 分别约为 0.200894、0.201063。

这是可直接核对的材料支持差异，并非将单元 μ 的百分差当作 Id 贡献。当前 [CoupledDDAssembler](../../src/equation/CoupledDDAssembler.cpp) 通过 `transportEdgeCoupling` 选择材料侧输运几何，否则使用完整网格边几何；本轮冻结配置沿用后者。**该开关的存在不等于其余迁移率、Poisson、源体积及端口支路已经通过组合资格。** 本轮只在隔离 runner 加入了常量分布源，未启用生产开关。

低 NWell 控制说明仍需保留耦合解释：其原始 Id 误差约 +2.22%/+1.44%，而组合方向导数已有 −7.76%/α、−6.19%/α；高 NWell 原误差为 +12.9842%/+9.8124%。这些原始误差以原生 Id 为分母，导数以 Vela Id 为分母。上一轮的 Poisson 电荷体积项与输运项存在明显抵消，不能将本轮小扰动外推成单项完整替换一定改善高低配对。

证据：[逐边贡献](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/distributed_edge_contributions.csv)、[主导边各材料几何及迁移率](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/dominant_edge_material_support.csv)、[材料支持复核脚本](../../scripts/explain_simplemos_native_transport_support_20260907.py)、[支持证据封存](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/material_support_evidence.json)。

## 6. 失败记录、分辨率定义与工作边界

原生导出前两版各四次尝试以 exit 5 结束：第一版请求了运行时不支持的 edge eMobility，第二版使用了运行时不识别的 Potential 绘图别名。第三版限定为受支持的单元迁移率、拓扑、系数、体积分量导出，4/4 重闭合成功。前八次失败没有作为有效物理状态纳入对比；脚本、日志和输入均保留，见[原生尝试账本](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/native_attempts.csv)。

第一版独立几何源保留了 `g_native/g_Vela−1` 中接近浮点舍入的量。单位源预检在四工况分别有 522、525、475、479 行无法通过既定残差相减误差界，因而没有进入 DC。这是诊断源在大流量背景下不可分辨的问题，原失败记录和合同没有改写。

后续新建并冻结了有明确分辨率的方向：当 `|g_native/g_Vela−1|≤1e-10` 时仅将这条**诊断几何源**的比值设为 1，其余不变；迁移率源取组合源减该几何源，保持可加。预检误差界、Newton 接受条件和响应门槛没有放宽。被滤除量对单位 α 响应的逐边绝对和上界，归一化 Id 后最大 1.08e-14；实际总预测变化最大约 1.53e-16。它远小于本轮信号，但属于源定义变化，不能省略说明或声称第一版通过。见[保留失败](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/retained_preflight_failures.csv)及[可加性与分辨率影响](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/factor_additivity_and_resolution.csv)。

当前可以据此安排下一阶段：以已校准的原生硅侧系数及 box 加权 μ 构造完整守恒输运候选，核对电子/空穴与端口的一致使用；同时对原生 Poisson 系数、Si 电荷/源体积方向做相应独立校准，再在原四工况冻结门槛下进行单项和联合自洽 A/B。只有高 NWell 绝对误差、高低配对及高 Vd 控制均改善，才扩展工作点与完整曲线。本轮未执行这些完整模型替换，也没有完成空穴或 Poisson 分布源响应校准。

## 复现与封存

本轮隔离 C++ runner 编译成功；5 项 SG 数学性质测试通过，覆盖平衡零通量、方向反转、密度/准费米形式等价、共同势平移不变性和电子/空穴电流符号。三个方向共 24 次源预检、60 次 DC、60 次全行事后探针；另保留第一版几何源的 8 次失败批次预检。独立脚本从原始端口状态与逐行 CSV 重算资格、正负响应、直接项、自由行投影和分解可加性。

新增的是验证脚本、测试及证据；正式核心求解代码未改。工作树先前已有的 runner 修改仍为 127 行新增/1 行删除，未覆盖。没有因这些诊断任务运行全量 CTest。

从本工作树根目录使用 `D:/msys64/ucrt64/bin/python.exe`：准备脚本先冻结合同，执行脚本按合同运行，分析/封存脚本核对输入哈希。输出目录为只追加证据，已有文件不会被静默覆盖；重做实验应另建运行目录。生成数据位于忽略目录 `build-release/simplemos_masetti_runtime_20260907/`，报告与表格位于 `reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/`。

- [运行时准备](../../scripts/prepare_simplemos_runtime_geometry_only_20260907.py)、[运行时审计](../../scripts/audit_simplemos_runtime_geometry_20260907.py)。
- [隔离分布源实现](../../scripts/build_simplemos_distributed_flux_20260907.py)、[组合校准](../../scripts/validate_simplemos_distributed_transport_20260907.py)、[分解校准](../../scripts/validate_simplemos_transport_factors_resolved_20260907.py)。
- [SG 性质测试](../../tests/regression/test_simplemos_native_sg_calibration.py)、[独立复核与封存脚本](../../scripts/seal_simplemos_transport_calibration_20260907.py)。
- [独立检查及指标](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/independent_check.json)、[最终文件哈希封存](../../reference_tcad/simplemos_sentaurus2022/masetti_runtime_calibration_20260907/validation_evidence.json)。
