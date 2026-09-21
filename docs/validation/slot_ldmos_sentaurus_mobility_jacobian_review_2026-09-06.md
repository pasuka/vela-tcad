# SLOT-LDMOS 验证进展与迁移率 Jacobian 完整性审查

日期：2026-09-06。范围：当前 `codex/slot-ldmos-bvds-next` 工作树的只读分析与诊断探针；本轮未修改求解器、模型或验收标准。

## 结论

当前已取得按既定交接规则验收、可重启的 1.25 V 无雪崩检查点，但尚未完成同物理模型的 Sentaurus/Vela BVDS 对照。原 Sentaurus SLOT 参考启用 IALMob(FullPhuMob)，当前 Vela 检查点使用 Masetti 高场迁移率，必须先区分模型差异和求解器差异。

完整 Newton Jacobian 必须包含迁移率对所有实际参与计算的未知量的依赖。当前代码已包含高场迁移率导数及邻接三角形节点耦合，`jacobian_field_derivatives=true` 不能作为完整性证明。本次确认了表面迁移率与单元向量 GradQF 组合下的一处端点电势导数漏项，并发现边迁移率调用没有传入实际载流子密度。这些问题影响后续表面迁移率/IALMob 对照，不能直接用来解释当前 IALMob 关闭的 1.25 V 尾部。

## 已有验证进展

| 项目 | 当前证据 | 判断 |
|---|---|---|
| Sentaurus 二维电流与串联电阻 | 1 V 直接偏压及 1e12 ohm*um 校准完成 | 电流用 A/um；不能再乘一次宽度 |
| Sentaurus 高压参考 | 无雪崩、IIC、自洽雪崩及外电阻曲线已完成 | 在 1e-7 A/um 准则下，内部 BVDS 为 38.5209012 V |
| Vela 严格 1 V | 10 次 Newton；Id=5.1772757291342995e-16 A/um | 相对 Sentaurus 直接 1 V 的 5.90522510576319e-16，低约 12.33%，约 -0.05714 decade；模型不完全相同 |
| Vela 严格 1.25 V | 2 次真实 Newton；Id=6.5889954852964324e-16 A/um | `final_global_continuity_restoration`；已验证重启 |
| 1.25 V 全局连续性 | 电子/空穴失配 2.4601233527902148e-19 / 2.3598226144186885e-19，使用现有内部单位 | 通过原 1e-12 绝对准则；不是 A/um 数值 |
| 1.25 V 重启 | 电流及电子/空穴密度保持一致 | 参考电势重组舍入修复已验证 |
| Vela 1.5 V、2 V 及后续 BVDS | 本批尚无符合新合同的证据 | 不能引用旧高压分支作为当前 Vela BVDS |

严格验收沿用既定局部行资格筛选；1.25 V 的 qualified rows 为零，仍有 1186 条处于 warn 模式的 QF 范围诊断。因此，通过全局闭合不等于所有局部方程都通过独立绝对误差检查。

Sentaurus 无雪崩高压曲线在内部 0.9015892251 V 与 1.3539269051 V 之间跨过 1.25 V。本轮没有把插值值称作精确的 Sentaurus 1.25 V 检查点。已有 0.90625 V 与 1 V 的 Vela 电流也存在下降，应在统一模型和验收下复核低压曲线连续性，不能仅凭单点通过就宣称整条曲线已验证。

## 对照原始日志的差异

### 原始日志来源不能混用

1. SLOT 参考日志：
   `D:/code-repo/vela-tcad/build-release/reference_tcad/slot_ldmos_sentaurus2022/run01/sentaurus_bvds_result/logs_curves/`。
   本次读取 `unit_resistor_1v_direct.log_des.log`、`avalanche_off_60v.log_des.log`、对应 deck 及归一化曲线。
2. 普通 LDMOS 原始日志：`D:/gpt-tmp/LDMOS/output_des.log`。
   该日志是 **R-2020.09、三维四面体算例**，33153 顶点、176839 四面体，使用 ILS(21)、64-bit double、Lombardi。
   它可以帮助发现模型/导数问题，但不是二维 SLOT 电流/BVDS 的同器件基准。

| 方面 | SLOT Sentaurus 日志/输入 | 当前 Vela 1.25 V | 对差异的解释力度 |
|---|---|---|---|
| 低场与表面迁移率 | Silicon 下 IALMob(FullPhuMob PhononCombination=1)；高场为 GradQF 驱动的 Caughey–Thomas 饱和 | `masetti_field`，无 IALMob | **已确认的物理模型差异**；单点电流误差不能全部归于数值算法 |
| 参数 | 日志列出自定义 IALMob mumin、C、delta、alpha_sr 等覆盖值，例如电子 alpha_sr=5e-22 cm3 | 当前 Masetti 默认/导入参数集 | 必须逐项比对有效参数及单位，仅模型名称相近不够 |
| 数值精度 | 日志第 218 行：80-bit long double | `Real=double`，稀疏矩阵与状态为 double；部分诊断/归约采用 long double/Boost | 与已发现的微小 QF 更新丢失一致；并不证明应全面更换精度 |
| 离散几何 | AverageBoxMethod；27275 三角形，2 个 non-Delaunay 元素 | 同三角形数；Delaunay 聚合回退、mixed Voronoi；运行日志显示 882 次边聚合修复、1 次回退 | 仍需逐项比较边耦合和控制体积，不能直接把 obtuse 数量当作 non-Delaunay 数量 |
| 线性解法 | Pardiso；普通 LDMOS 日志则为 ILS | 当前 Newton 使用 LinearSolver，默认 Eigen SparseLU/COLAMD；UMFPACK 已检测到但不等于此路径采用它 | 后端不同已确认，尚未证明是主要误差源；历史单次运行未单独保存后端名称 |
| 停止判据 | Digits=5、RHSmin=1e-10、ErrRef(e/h)=1e10 等 | 自定义缩放残差、stall floor、局部资格筛选和双区间全局闭合 | 数值不能直接一一比较；需共同的物理残差与电流判据 |

SLOT 日志顶点数为 14069，Vela 导入节点数为 13881，三角形数相同。这里仅记录计数差异；未核实其内部节点展开/接口表示映射，不能据此断言导入漏节点。

输入 deck 中 `#Derivatives` 被注释。日志明确列出 `With avalanche derivatives` 与 `No incomplete Newton`，但这些行不能单独证明每条迁移率求导链均已打开。不能把“雪崩导数”解释成“完整迁移率导数”，也不能仅凭注释行判断对应版本的默认行为。

## 当前 Jacobian 的覆盖与缺口

### 已有高场迁移率求导

- `include/vela/physics/MobilityModel.h:88`：`jacobianFieldDerivatives` 默认 true；最终检查点配置也显式为 true。
- `src/equation/CoupledDDAssembler.cpp:4202`：有场依赖时，边端点输运通量采用局部中心差分，步长为 `1e-6*max(1,abs(value))`。
- 同文件 `4265` 附近：对 transport-cell-vector GradQF 的完整邻接节点 stencil 求迁移率响应并散入 Jacobian；对应稀疏结构在 `3039` 附近保留。

因此，当前 `analytic` 选项实际包含局部数值微分，不能称为处处解析求导。是否足够准确，要看具体状态、步长与活跃分支；尚无本轮真实 1.25 V 状态的完整逐列导数验收。

### 已确认：vector GradQF 与表面迁移率组合漏端点法向场响应

`4208` 附近在 vector GradQF 模式令 `fixedMobility=mun`；`4030` 附近只要 fixedMobility 非负就不重新计算迁移率。因此，端点 psi 扰动没有更新表面法向场对迁移率的作用。补充表面 psi stencil 的分支在 `4315` 附近又跳过了边的两个端点。空穴存在对称路径（`4350` 之后及 `4448` 附近）。

本次使用已有六节点 Si/SiO2 测试拓扑，在每格 **10 nm** 的几何下做独立中心差分。保持默认 Lombardi 参数、无复合、无雪崩，比较两个载流子方程对六个 psi 自由度的矩阵块：

| 配置 | 中心差分步长 | 块相对 Frobenius 误差 |
|---|---:|---:|
| Lombardi + edge_projection | 1e-7 V | 1.1850e-10 |
| Masetti 高场、无表面模型 + transport_cell_vector | 1e-7 V | 9.7649e-11 |
| Lombardi + transport_cell_vector | 1e-7 V | **4.2998823e-3（0.4300%）** |
| 同上，步长扫描 | 1e-5 至 1e-8 V | 稳定在约 4.29988e-3 |
| 同上，仅诊断性放大声子散射项 100 倍 | 1e-7 V | 7.5881e-2；不是实际 SLOT 参数或物理结果 |

默认参数下最大差异位于矩阵行 9、列 2：装配值 4.855707226352e18，中心差分值 4.83310989776e18。这是小网格内部残差导数值，不是端电流。

探针使用 SI 输入、300 K、Silicon 掺杂 1e23 m^-3；改变为无表面的 Masetti 高场控制后，vector 分支恢复差分一致性。这进一步把该反例定位到表面响应，而不是普遍否定 vector GradQF 的实现。

原有 Lombardi 测试网格为每格 1 um，相对 10 nm 的界面作用长度，其表面贡献容易被压低；整块相对范数测试也可能掩盖局部漏项。既有 `[mobility]` 两个测试本次仍通过（4 条断言），但不足以排除此问题。

### 已确认：边迁移率未接入实际 n、p

`include/vela/equation/AssemblerUtils.h:1256` 附近调用 electronMobility/holeMobility 时传入的是 `n=0,p=0`。而 `src/physics/MobilityModel.cpp:364` 之后的 Lombardi 公式支持由同类/异类载流子密度控制的粗糙度指数。

当 alpha 或相关载流子系数非零时，问题同时涉及 **残差中的物理模型和 Jacobian**：不能只加上非零的 dmu/dn、dmu/dp，却仍让残差使用零密度。两者必须使用同一个真实载流子状态和一致的边/单元密度映射。默认 Masetti 对自由载流子密度无此依赖，默认 Lombardi alpha=0，因此该缺口不直接解释当前 Masetti 无表面分支。

### 条件性缺口：旧节点雪崩近似分支

`src/equation/CoupledDDAssembler.cpp:5102` 附近明确省略了场和迁移率导数，仅保留局部密度导数。它不是所有雪崩路径的描述：triangle-GSS/local-AD 等其他路径另有完整局部源项导数。本次 1.25 V 关闭雪崩，因此此处不是其失配来源；以后启用雪崩必须核对实际走到的 source/Jacobian 分支。

## “完整偏导数”应覆盖什么

设一条边的离散通量为 `F=mu*S`，`S` 为剔除迁移率因子的 SG 通量。对于邻域任意未知量 u：

\[
\frac{\partial F}{\partial u}
=\mu\frac{\partial S}{\partial u}
+S\left(
\frac{\partial\mu}{\partial E_d}\frac{\partial E_d}{\partial u}
+\frac{\partial\mu}{\partial E_\perp}\frac{\partial E_\perp}{\partial u}
+\frac{\partial\mu}{\partial n}\frac{\partial n}{\partial u}
+\frac{\partial\mu}{\partial p}\frac{\partial p}{\partial u}
\right).
\]

- 电场或 GradQF 驱动力的重建、归一化、界面投影及低密度平滑都要按实际实现求导；邻接三角形的第三个节点也在依赖范围内。
- 对 Boltzmann 状态，在固定 ni/T 的平滑区间，`dn/dpsi=n/Vt`、`dn/dphin=-n/Vt`，`dp/dpsi=-p/Vt`、`dp/dphip=p/Vt`；Fermi、限幅或其他统计必须使用其实际导数。
- 固定掺杂、固定网格、固定温度不额外引入 Newton 未知量导数；“完整”是对当前实际求解变量的完整导数。
- `abs/max/clamp` 交界处不一定存在唯一经典导数，应一致处理活跃分支，并通过单侧或避开切换点的差分核对。
- 开启雪崩时，还需沿 `G=alpha_n|Jn|/q+alpha_p|Jp|/q` 等实际源项定义继续传播 `dalpha` 和 `dJ`，不能只在输运通量补齐后就称整个 Jacobian 完整。

## 建议的验证顺序

1. 优先补充能激活表面作用、隔离端点和第三节点响应的回归；修复 vector GradQF 表面端点漏项。此项是纯 Jacobian 一致性问题，可保持当前残差模型。
2. 在启用载流子依赖的表面模型前，先明确与 Sentaurus 相同的密度映射/参数/单位，再使残差与其导数同时使用真实 n、p。这会改变现有简化模型的有效行为，应单独评估。
3. 对当前 1.25 V 固定状态，分别检查无表面输运、复合和以后启用的雪崩源项，使用多步长逐列/局部基方向核对；避免只看近零 Jv 的相对误差。
4. 建立同物理的 Masetti 无表面对照，或在两侧完整对齐 IALMob。现有 12.33% 的 1 V 差异目前不能归因到某一个参数或导数。
5. 导数与模型核对后再按现有标准验证 1.5 V、2 V；保留检查点的重启、电流与局部/全局残差证据，不通过放宽阈值消除差异。

以上是分析建议，本轮没有执行模型/求解器修复或新的升压计算。

## 本轮证据

- 小网格探针源码及 JSON：`build/slot_ldmos_analysis_20260906/mobility_probe_central.cpp`、`mobility_probe_central.json`；使用 UCRT64 C++20 链接本工作树已有 Debug `libvela_core.a`，不调用线性求解器。
- 无表面控制：同目录 `mobility_probe_bulk_control.cpp`、`mobility_probe_bulk_control.json`；源码及链接库哈希存于 `source_evidence_hashes.json`。
- 1 um 原测试尺度及 10 nm 前向差分探索也保存在同一目录；上文定量结论采用独立中心差分结果。
- 已有测试：`build/test_newton_solver.exe "[mobility]"`，2 用例/4 断言通过；未重复跑全部测试。
- 严格 1 V 原始 CSV：`build/slot_ldmos_noavalanche_1v_post_globalization/outputs/stages/01_delaunay_noavalanche_adaptive/attempt_005_v_1/iv.csv`。
- 严格 1.25 V 配置、原始 CSV、运行日志和重启对照：`build/slot_ldmos_e3_20260905/final_strict/`、`final_restart/`、`restart_validation.json`。
- 前次验证完整结果：[1.25 V 报告](slot_ldmos_e3_strict_1p25_2026-09-05.md)。
