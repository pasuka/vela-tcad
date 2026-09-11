# LDMOS IALMob 局部差异定位 — 2026-09-11

## 结论及资格边界

本轮定位出一个重要的比较口径问题：**Sentaurus 输出的节点迁移率不能直接
当作该节点状态代入 IALMob 公式的结果。** 关闭高场饱和、改用原生单元输出后，
4054 个体区三角形的迁移率与各顶点计算值的解析 Voronoi 加权平均高度一致。
这一结论有同状态原生对照支持，没有拟合模型参数。

它支持当前独立计算核在该体区样本内的正确性，但不等于整个 IALMob 已通过：
钝角单元及邻居仍有残差，节点输出回映射尚未重构，界面项仍未闭合。
本轮未修改生产求解器、物理参数、门限或默认几何策略，也未开展耦合 Jacobian
和 D4 自洽曲线验收。前次实现及测试见
[计算核执行记录](templates_ldmos_worker_and_ialmob_execution_2026-09-11.md)。

## 配置、运行及控制检查

- 原生状态：D4-classical 的 `n4_des.tdr`，Vg=8 V、Vd=40 V，源极/衬底 0 V。
- 网格：10241 节点、19782 三角形；5723 个 Si 节点、10515 个 Si 三角形。
  使用原始 TDR 和完全相同的全局节点/单元编号，坐标逐项核对。
- 300 K，Fermi、OldSlotboom、SRH/Auger，IALMob(AutoOrientation)，
  原生日志确认 FullPhuMob、PhononCombination=1；使用实际 Siliconc100.par。
- 原生工具：本机 VM 的 Sentaurus T-2022.03-SP2。经用户授权上传至隔离目录
  `/home/tcad/sentaurus_runs/vela_oracle_2022/ialmob_location_20260911_r1`。
- 本地计算核：UCRT64 Release，`-O3 -DNDEBUG`，gprof 关闭。
  `ialmob_probe.exe` SHA256 与前次冻结清单一致：
  `90f127bc15a847823e829909845d18200a2cb10aa279c769e23d553a5bea9245`。
  本轮局部求值不调用线性求解器；既有曲线环境为 Eigen SparseLU/COLAMD。
- 几何以 μm 表示；TDR 密度 cm⁻³、场 V/cm、迁移率 cm²/(V·s)。
  计算核输入转换为 SI、输出乘 10⁴ 转回 cm²/(V·s)。

| 原生诊断 | 相对原始 Physics 的变化 | 输出 |
| --- | --- | --- |
| baseline_node | 无 | 节点迁移率 |
| baseline_element | 无 | `eMobility/Element hMobility/Element` |
| low_field_node | 仅关闭 HighFieldSaturation | 节点迁移率 |
| low_field_element | 仅关闭 HighFieldSaturation | 单元迁移率 |
| low_field_geometry | 同 low_field_element，增加缺失时导出几何系数的 Math 选项 | 单元迁移率及 MeasureCoefficients.debug |

五组 Solve 均仅含 Load、Plot，**每组 0 次 Newton 更新**。
前四组接触电压完全一致；势最大绝对变化小于 1e-10 V，密度最大相对变化小于
1e-10，掺杂和温度不变。第五组与 low_field_element 的状态字段及 10515 个
电子/空穴单元迁移率逐项完全相同。

原始 D4 与 baseline_node 的电子/空穴节点迁移率最大相对变化仅
8.01e-12 / 9.79e-11，因此可以排除保存状态重载造成此前百分之几十长尾。
这些 Load/Plot 作业包含许可证等待，约 200 s 的原生墙钟不是曲线性能指标。

## 差异来源一：驱动场与节点输出位置

原 deck 未指定驱动方式，日志确认 HighFieldSaturation 使用 GradQuasiFermi。
输出的 Eparallel 不能直接代替它；RefDens=1e12 还引入低密度下向界面平行场
的过渡。逐点使用输出 QF 梯度也未必等于实际单元求值位置的驱动。

沿用前次体区、高密度 ≥1e16 cm⁻³、无量纲 QF 驱动 <0.01 的节点样本：

| 节点比较方式 | 电子绝对相对误差 P95，253 点 | 空穴绝对相对误差 P95，311 点 |
| --- | ---: | ---: |
| 前次逐点公式 + 输出 QF 高场估算 | 24.9336% | 5.2251% |
| 本次原生关闭 HFS，逐点低场公式 | 19.3748% | 4.8103% |

关闭 HFS 后长尾仍存在，不能把问题全部归因于高场公式。
长尾集中在源/漏附近掺杂快速变化区域：电子邻接掺杂跨度 ≤0.001 decade
的 112 个样本，P95 为 0.02736%；跨度 >1 decade 的 26 个样本为 52.7463%。
这是空间关联证据，不能单独作为因果证明。

## 差异来源二：体区单元迁移率平均

选择三个顶点均距 Si/SiO2 界面至少 0.2 μm 的 4054 个 Si 三角形。
局部计算以独立 FullPhuMob 体区极限求值（法向场 0、界面距离 1 m），原生对照
只关闭 HFS；有限界面距离效应没有被人为从原生值中扣除。
该样本不同于上面的 253/311 个节点，不能直接跨表计算误差改善倍数。

对顶点 i，采用解析二维 Voronoi 子体积权重
`w_i = Σ(j≠i) |x_i-x_j|² cot(theta_k) / 8`，k 为第三顶点，
并计算 `mu_element = Σ(w_i mu_i) / Σ(w_i)`。权重完全来自原网格，无拟合参数。

下表误差均为 **无量纲相对值**；P95 对绝对误差取分位，最大值未排除异常：

| 方法，同一 4054 单元 | 电子 P95 | 空穴 P95 | 电子最大绝对误差 | 空穴最大绝对误差 |
| --- | ---: | ---: | ---: | ---: |
| 三顶点等权平均 | 1.189943e-2 | 8.812896e-3 | 1.545087e-1 | 1.223403e-1 |
| 解析 Voronoi 加权 | 7.425160e-12 | 9.832180e-12 | 2.531851e-2 | 9.798789e-3 |

手册第 464 页说明每个单元顶点单独求值、默认按 Element 平均。
**具体解析权重的吻合是本轮实验推断，不是手册文字对权重公式的明示。**
近零 P95 不能覆盖最大残差：绝对相对误差 >1e-4 的电子 6 单元、空穴 7 单元，
均为钝角三角形或其共边邻居。该 1e-4 仅用于定位，不是新增验收门限。
若计数阈值降为 1e-6，分别为 12/97 个，完整记录保留，未删点评分。

| 代表单元 | 质心 (x,y)，μm | 电子相对偏差 | 空穴相对偏差 |
| --- | --- | ---: | ---: |
| 10232 | (-9.97702, 10.25523) | -2.53185% | -0.97988% |
| 10238 | (-9.97539, 10.25521) | -0.74669% | -0.21574% |
| 8680 | (-9.47043, 0.05729) | -0.32538% | -0.24509% |
| 3971 | (-8.01276, 4.15365) | +0.23069% | +0.07895% |

剩余位置支持继续检查钝角几何处理，但尚不足以认定全部误差由某一种几何修正引起。
其中只有一部分位于负对边 cotangent 和判定的非 Delaunay 边附近，不能统称为
“非 Delaunay 残差”。

## 差异来源三：Box 几何权重、节点映射及界面项尚未闭合

第五组新建空子目录，启用 `BoxMeasureFromFile(GrdNumbering)` 和
`BoxCoefficientsFromFile(GrdNumbering)`，目录内事先不存在 debug 文件，
按手册第 1187–1188 页导出当前几何系数。

**复核修正：未显式写 AverageBoxMethod 不代表它关闭。** 本次所有相关原生日志
均显示 `CVPL_AverageBoxMethod = TRUE`；新导出的文件与旧 AverageBox 证据
SHA256 完全一致：`05ec1a43936e0ff71bea7cdfd956d80c7d692b228187ccae8625ec53ce58d49e`。
因此排除“新旧 Box 选项不同”这一假设。

直接把导出的 Measure 当作迁移率平均权重，体区电子/空穴 P95 反而为
0.09908% / 0.03759%，最大误差为 11.9561% / 6.5788%。
这说明 Box 装配体积不能未经验证就替代迁移率的全部平均权重。
从原生单元迁移率回算原生节点输出时，面积、等权、Voronoi 及若干常见平均方式
均未精确闭合；没有用拟合权重掩盖这一差异。

原先两个严格筛选的平直 {100}、正 Enormal 节点 4669/4670，在 HFS-off
原生对照中，局部公式电子偏差 -8.9432% / -9.0289%，空穴
-13.0692% / -13.1708%。这两个重掺杂样本仍受输出位置问题影响，
不能归因于某个界面散射系数，也不代表整个沟道的资格。
界面法向场符号、单元内法向/晶向选择以及节点回映射仍需独立对照。

## 下一步建议

1. 在上述钝角单元及共边邻域上核对迁移率平均的实际截断/重分配规则；
   同时保留普通单元的解析权重对照，避免误改已吻合部分。
2. 增补原生单元法向场、QF 梯度和 InterfaceOrientation 诊断，明确界面
   求值位置与节点输出映射，再扩大沟道局部样本；现阶段不调整散射参数。
3. 局部低场和高场对照闭合后，接入主求解器的 Nd/Na、n/p、法向场及晶向依赖，
   补齐耦合 Jacobian 的数值验证，再开展 D4 曲线原门限验收。

本轮已有网格、状态、参数及本地 T-2022.03 手册足够继续，不需要重复提供资料。

## 证据和本轮检查

证据根目录（本机生成物，不纳入 Git）：
`reference_staging/templates_ldmos_followup_20260911/ialmob_difference_location/`。

- `native/plan.json`、`bundle/`、`raw/`：四组原生输入、SHA256、TDR 和日志；
  `native/controls_audit.json` 核对状态、偏压、输出支持及 0 Newton。
- `native/geometry/`、`native/geometry_export/`：第五组日志、debug 系数和导出字段。
- `geometry_weights_summary.json`、`native_cells_summary.json`、
  `default_measures_summary.json`：同样本平均方式及 HFS-off 对照。
- `e_remaining_cells_v2.csv`、`h_remaining_cells_v2.csv`：残差、节点编号、
  质心、钝角及非 Delaunay 邻域标记；`doping_strata.json` 为节点掺杂分层。
- `surface_node_control.json`、`final_summary.json`、`evidence_hashes.json`：
  界面两点对照、最终范围和证据哈希；`tail_locations.png` 已修复标签遮挡并检查。

实验脚本位于上述证据根目录的父目录，包括
`audit_ialmob_native_controls.py`、`analyze_native_ialmob_cells.py`、
`probe_ialmob_geometry_weights.py`、`probe_ialmob_native_measures.py --default`、
`finalize_ialmob_location.py`。它们依赖已冻结的外部证据，并非独立发布的运行入口；
部分脚本拒绝覆盖 CSV，重放须使用新的证据目录。

本轮实际完成五组原生诊断、状态/输出位置检查、Release 探针哈希校验、
逐单元对照、图像检查和文档链接/差异检查。
未改动生产 C++，没有重复全套构建；前次 758/758 Release CTest 是前次实现
的结果，不作为本轮新增曲线或 IALMob 集成资格。
