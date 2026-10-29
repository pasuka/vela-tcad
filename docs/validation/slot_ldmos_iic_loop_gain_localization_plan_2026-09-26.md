# SLOT LDMOS：35/37/38 V 碰撞电离自洽源项偏差——原因分析与源头定位执行计划（2026-09-26）

本文是对 [IIC 源项差异预算](slot_ldmos_iic_source_budget_2026-09-25.md) 的
只读分析和后续执行计划。分析基于 worktree `slot-ldmos-bvds-next`
（HEAD `3e34336f`，工作区含未提交文档）、`build/s25iic/`、`build/s25to37/`、
`build/s25to38/` 既有证据，以及两段只读 Python 核对（第 1.2 节）。
本文没有修改源码、配置或任何既有证据目录，也没有启动求解。

计划面向由大模型代理执行：每一步给出目的、输入、工具、输出、判据和禁止事项。
执行时先读 [AGENTS.md](../../AGENTS.md)、[SLOT 状态交接](slot_ldmos_status_and_handoff_2026-09-24.md)
和 [配置模式](../config_schema.md)。

**2026-09-26 审阅修订**：按审阅意见修正五处方法表述——(1) 由 `M` 反推的 `K`
改称“等效标量反馈指标 `K_eff`”，它复现约 9% 差异只是代数一致性检查；(2) 不再以
IIC 总积分接近排除基础源分布的放大贡献，改由伴随权重评估（`ΔS_src`），并把
“差异只能在响应侧”降为待检验假设（源值匹配 ≠ 导数匹配，增加反向冻结回放）；
(3) S4.2 的 α 缩放灵敏度改为 `d ln S_SC/d ln s = M`、`d ln M/d ln s = K·M`；(4) S0 的
α 比较改用实际电子/空穴驱动场，一次 Newton 步改称局部误差估计，新增 S0.6
固定电场假设检查；(5) 伴随核对式补上负号与单位系数 `c`，归因同时报告有符号与
绝对贡献及抵消指数，Python 复制算子改用按行操作数舍入界验收。

**2026-09-26 第二轮审阅修订**（审阅方已用 Wolfram 做符号/数值核对）：(a) §1.1 补充标量
模型的精确关系与分布式线性模型 `z = b + Lz`，注明非正规算子的响应放大不能仅凭特征值
判断，真实环路谱分析推到伴随定位之后（S4.4）；(b) §3 改用精确拆分“状态差 + 冻结源核差”，
一阶式写为 `wᵀ[F_V(x_V) − F_V(x_S)]`，仅在 `wᵀF_V(x_V)` 可忽略时才简化为 `−wᵀF_V(x_S)`；
(c) S1 增加五项工具自身校验（两个线性系统的数值残差、`gᵀδx = −wᵀr`、总源方向差分、
Jacobian–残差一致性、正负/绝对贡献），分组定义冻结；(d) S4.2 改为先用已有伴随解计算
离散灵敏度 `d ln S/d ln s = 1 + wᵀQ/S`（不要求单模态），再用少量 ±α 重收敛核对；(e) S0.1
改称“剩余收敛误差估计”，增加同方向不同步幅的线性性检查；(f) 验收改为同时看源差、
端电流、势场、局部残差与配置一致性，不以 `K_eff` 接近原生作为根因成立的单一条件。

**2026-09-26 第三轮细节修订**：(i) S4.2 的 α 灵敏度 `σ = 1 + wᵀQ/S` 要求 `J`、`w`、`Q`、`S` 来自
同一合格基态：在 Vela 自洽态 `x_V` 上另做一次伴随解（S0.1 同时捕获 `x_V` 上的系统），不复用
`x_S` 上的 `w`；(ii) 线性校验改用归一化后向误差 `η = ‖Jδx + r‖/(‖J‖‖δx‖ + ‖r‖)`（伴随同形），
内积恒等式用混合绝对/相对判据；(iii) S4.4 注明 `J_SC⁻¹J_IIC = (I − L)⁻¹` 的成立条件，同时检查
最大模特征值与最接近 1 的特征值，范数比较仅作诊断（正规矩阵二范数也是 `1/min|1 − λ_i|`）；
(iv) ±1% 对照差不单独称为非线性余项（含重收敛误差、差分截断、分支切换），差异明显时才补
±0.5% 检查步长依赖。

**执行顺序（冻结）**：

| 序 | 工作 | 对应阶段 | 目的 |
| ---: | --- | --- | --- |
| 1 | 修正计划中 K、α 灵敏度、残差符号、“误差上界”表述 | 本文 | 避免把待验证假设写成结论 |
| 2 | 35 V 校验源梯度、线性捕获与伴随恒等式 | S0，S1 第 1～2 步及五项校验 | 先验证诊断工具本身 |
| 3 | 35 V 方程块/界面/空间窗口归因 | S1 第 3～4 步 | 判断 H1 是否占主导 |
| 4 | H1 单因素冻结控制 + 独立重收敛 | S2.1，S3（35 V） | 形成可重复的修正证据 |
| 5 | 扩展到 38 V，核对 α 灵敏度与非线性余项 | S1/S2/S3（38 V），S4.2，S4.4 | 检验近击穿放大机制 |

## 1. 现象重述：数据与“约 1% 等效反馈差被 M 放大”的图像相容（待检验，非已证明）

### 1.1 用等效标量反馈指标重述既有数字

若器件近似满足“固定电场、单一主模态、线性反馈”，则 `S_SC = S_IIC/(1 − K)`。
本文只把 `K_eff = 1 − 1/M`（`M = S_SC/S_IIC`）定义为**等效标量反馈指标**：它是
既有积分结果的另一种写法，不是已证明的环路主特征值；本批小电流（~1e-14 A/µm）
也不足以单独证明局部空间电荷影响可忽略，固定场假设由 S0.6 定量检查。

标量模型中 `dM/dK = M²`（M 的**绝对**变化），`dM/M = M·dK`（相对变化的一阶放大系数
是 `M`）；以原生为参考的精确关系为 `(M_V − M_S)/M_S = M_SΔK/(1 − M_SΔK)`。

| Vds / V | M_Vela | M_原生 | K_eff,Vela | K_eff,原生 | δK_eff | δK_eff/K_eff | 一阶 M_S·δK_eff | 精确 M_SδK/(1−M_SδK) | 实测 δM/M |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 35 | 2.209472 | 2.190339 | 0.547403 | 0.543449 | +0.003954 | +0.728% | +0.86595% | +0.87352% | +0.874% |
| 37 | 3.582107 | 3.500900 | 0.720835 | 0.714359 | +0.006476 | +0.907% | +2.26702% | +2.31961% | +2.320% |
| 38 | 8.945177 | 8.201328 | 0.888208 | 0.878068 | +0.010140 | +1.155% | +8.31565% | +9.06986% | +9.070% |

精确列与实测列相等是因为 `K_eff` 本就由 `M` 反推，属代数必然，不能独立验证反馈模型。

分布式线性模型应写为 `z = b + Lz`，`S = cᵀ(I − L)⁻¹ b`：总量放大同时依赖算子 `L`、
源分布 `b` 与观测权重 `c`；非正规算子还可能表现出仅凭特征值难以判断的预解式放大
（伪谱，参见 Trefethen 关于预解算子与伪谱的工作）。因此当前统一使用 `K_eff` 作描述性
指标，真实环路谱分析（主特征值与预解式范数）放在伴随定位之后（S4.4）。

关键推论（均为待检验的工作假设，不是结论）：

- 由 `δM/M ≈ M·δK_eff` 复现 +0.87% / +2.32% / +9.07% 是**代数一致性检查**：它说明
  三点数据与“约 0.7%～1.2% 的等效反馈差被 `M` 放大”的图像相容，但不独立证明
  存在同量级的真实环路特征值差。标量近似是否成立由 S0.6、S4.2 和 S4.4 检验。
- **基础源分布的放大贡献不能凭总积分接近排除**。令 `R = (I − L)⁻¹`，固定观测权重下
  一阶扰动为 `δS = cᵀRδb + cᵀR·δL·Rb`：第一项是源分布差，第二项是反馈算子差；
  即使 `1ᵀδb = 0`，第一项也不一定为零。审阅方的 2×2 数学反例（`L = diag(0.99, 0)`，
  `b = (0.01, 0.99)`，`δb = (0.001, −0.001)`，`c = (1, 1)`）：基础源总量变化 0、源分布
  L1 扰动 0.2%，而同一反馈算子下总输出变化 +4.9749%（`M = 1.99`、`K_eff = 0.49749`，
  主特征值 0.99）。它不是当前器件模型，但说明 IIC 总差 ≤0.06%、节点 L1 约 0.38%
  不构成反馈后贡献的上限，也说明 `K_eff` 与主特征值可以相差很远。需要 S1 中用
  伴随权重显式评估（第 3 节 `ΔS_src`，对应 `cᵀRδb`）。
- **同状态源值匹配到 ppm 不等于源项导数匹配到 ppm**。既有冻结回放只约束了原生
  自洽态与 IIC 态上的源值；源函数在 Vela 态邻域是否一致，由 S1 的反向冻结回放和
  线性化有效性门检验。因此“差异在响应侧（Poisson 场与输运离散）”是当前**优先
  假设**，不是排除其他项后的结论。
- `δK_eff/K_eff` 随偏压增长（0.73% → 0.91% → 1.15%），需要单独解释（第 2 节 H4）。
- 因此问题改写为：**在 Vela 与原生的完整自洽算子差中，哪些行、哪些项经近奇异
  系统放大后构成了 35 V 的 +0.887% 与 38 V 的 +9.005%？** 答案由第 3 节的伴随归因
  给出，而不是由标量模型预先指定。

### 1.2 本次新增的只读核对（可由 S0 脚本复现）

只读读取 `build/s25iic/local_iic_density/`（Vela IIC 合格态的 observable 源项）、
`build/s25iic/native/fields_{35,38}V.npz`（原生 IIC Plot 场）、
`build/slot_element_ialmob_20260921/local_control_r2/native_geometry.json`
（原生 Si 顶点体积）和 `comparison.json` 的单位换算，得到 **无反馈** 状态下
Vela/原生逐节点源项比：

| 节点 | (x, y) / µm | 说明 | 35 V IIC 比 | 38 V IIC 比 | 38 V 电子比 / 空穴比 |
| ---: | --- | --- | ---: | ---: | --- |
| 8237 | (0.0472, 2.1618) | SC 差异最大节点 | 1.00346 | 1.00368 | 1.0022 / 1.0043 |
| 8102 | (0.1181, 2.1252) | SC 差异第 2 | 1.00200 | 1.00224 | 1.0018 / 1.0035 |
| 8166 | (0.0472, 2.1984) | SC 差异第 4 | 1.00333 | 1.00348 | 1.0019 / 1.0041 |
| 8247 | (0.0330, 2.1618) | 热点最浅节点 | 1.00380 | 1.00403 | 1.0023 / 1.0045 |
| 10490 | (0.0049, 3.5953) | Si/SiO₂ 界面节点，同算子电场最大差 | **1.48826** | **1.48672** | 1.131 / 1.488 |
| 10492 | (0.0048, 3.6054) | Si/SiO₂ 界面节点，电势最大差 | **1.30988** | **1.30832** | 0.834 / 1.323 |
| 6626 | (0.0472, 3.5908) | 10490 下方约 42 nm | 1.01362 | 1.01350 | 0.9992 / 1.0221 |

而同一批热点节点在 **自洽** 状态下的比为 1.108～1.114（`hotspots_38V_total.csv`），
几乎均匀，略高于全局 `M_Vela/M_原生 = 1.0907`。即：热点处无反馈源项只差
+0.2%～+0.4%，自洽 +11% 主要来自反馈后的状态响应（其构成待 S1 分解）；而器件
表面 y≈3.6 µm 的界面节点存在一个 **与偏压无关、无反馈也存在** 的局部源项差异
（+31%～+49%，以空穴为主），积分权重小。它由驱动场、密度、迁移率和电流分配
共同决定，在 S0.3 分解之前只能称为“局部差异”，不能直接认定为电场缺陷。

器件几何（读取 `build/slot_ldmos_noavalanche_1p25_resume02/mesh.json`；x 为深度，
Si 表面在 x≈0.0049 µm，y 为横向）：

| 结构 | x / µm | y / µm |
| --- | --- | --- |
| drain 接触 | 0.0049 | 0.00–0.35 |
| Nitride_2.1（SLOT 漏侧相邻介质） | −0.204～−0.085 | 0.45–2.15 |
| SLOT 电极 | −0.413～−0.075 | 2.15–3.60 |
| Nitride_1.1 | −0.203～−0.011 | 3.497–3.590 |
| PolySilicon_1 栅 | −0.207～−0.007 | 3.600–4.250 |
| source 接触 | 0.0049 | 4.353–4.690 |
| substrate 接触 | 8.0 | 0–4.69 |

- 自洽差异热点 (0.03–0.12, 2.09–2.20) µm 位于 **SLOT 场板漏侧边缘 (y=2.15)
  下方 28～113 nm**；其正上方的 Si/SiO₂ 界面节点为 8552、8551、8470、8469、
  8402、8401、8311、8310（x=0.00495，y=2.015–2.272）。
- 电势/电场最大差节点 10490/10492 位于 **SLOT 栅侧边缘、多晶硅栅漏侧边缘与
  Nitride_1.1 的三重材料拐角** 正下方的 Si 表面。
- 两个差异位置都是场板边缘、多材料拐角下方的表面区。这是支持 H1（界面 Poisson
  离散差异）的空间线索，与同状态源核已匹配到 ppm 的事实相容，但本身不是证明。

### 1.3 既有证据已排除或降级的项

- 同状态源核**值**（α、极值谐均密度、最佳顶点 IALMob 迁移率、线性 QF 梯度、原生
  截断源体积、接触分支）：原生自洽态回放差 ≤2.2 ppm，原生 IIC 态 ≤2.13 ppm。
  源核**导数**与 Vela 态邻域的一致性尚未检验（S1 反向冻结回放）。
- 积分口径、接触节点遗漏、P1/box 权重：≤0.05%。
- IIC 基础源项分布：总差 ≤0.06%、节点 L1 ≤0.38%——只约束无反馈源值，其反馈后
  贡献待 S1 伴随权重评估，此处**不**作排除。
- 偏压/参考系错位：若为等效漏压偏移，需 δV≈36 mV 才给出 δK_eff=0.004
  （dK_eff/dV≈0.11 /V），而三点电势 RMS 差恒为 ~5 mV 且最大差在界面节点而非接触，
  不成立；但接触行残差仍在 S1 归因中单独列出。
- IALMob 边迁移率投影：`IalTransport.cpp` 以元素 box 系数占比加权
  （`c.weight/=coefficientSum`），`μ_edge·Σc_e ≡ Σ c_e μ_e`，与原生逐元素装配在
  代数上等价，不是候选。

## 2. 碰撞电离项偏差的可能原因（按优先级）

### H1（首选）Poisson 界面电荷体积策略为 `global`，造成近表面势场/电场系统差

- 机制：`build/s25to38/to38_o35/current/base.json` 没有 `discretization` 块，
  `ConfigParsing::parsePoissonChargeVolumePolicy` 因此取 `global`：Si/绝缘体
  界面节点的掺杂与载流子电荷按 **Si+介质** 的 mixed-Voronoi 体积积分，而原生
  按材料分区的元素体积装配。掺杂电荷误差与偏压无关，恰对应三点恒定的
  电势 RMS ~5 mV、最大 ~56 mV（在界面节点 10492）。
- 直接证据：[2 V 误差预算](slot_ldmos_energy_reference_and_source_volume_2026-09-21.md)
  的“同时改 Poisson 电荷体积”控制把 Si 势差 P95/最大值从 7.956/44.293 mV 降到
  **1.146/16.107 mV**，多晶硅从 52.6/72.6 mV 降到 12.4/23.2 mV；当时因 2 V 电流
  不变（−0.041006% → −0.040437%）保留了 `global`，文档已提醒“后续高压比较应
  明确所用 Poisson 几何，而非继承低压电流接近的结论”。
- 量级估计（非证明）：α 对驱动场的对数灵敏度 `b/F` 在 F≈3e5 V/cm 时约为
  电子 4.1、空穴 6.8。δK_eff/K_eff≈1% 只需热点处 δF/F≈0.15%～0.25%（约 400～750 V/cm），
  与 1.2 节热点 IIC 差 +0.2%～+0.4% 的量级相容。是否真由驱动场差造成，由 S0.2/S0.3
  的实际驱动场分解和 S1 伴随归因判断。
- 预测：切换 `discretization.poisson_charge_volume_policy =
  material_local_mixed_voronoi` 后，(a) 冻结原生态上 Vela 的界面节点 Poisson 行
  残差下降一个量级以上；(b) 自洽 35 V 的 δK_eff/K_eff 由 0.73% 降到 ≤0.25%，ψ RMS 差
  降到 ≤2 mV；(c) 10490/10492 的 IIC 比从 1.49/1.31 明显趋近 1。
- 判别：S2.1 冻结变体 + S3 自洽确认。

### H2 Poisson 边平均介电常数与原生分材料 AverageBox 耦合不同

- 机制：`AssemblerUtils.h::edgeEpsilon` 用相邻单元 `eps_r` 算术平均乘总耦合，
  界面边（Si|SiO₂、Si|Nitride、Poly|SiO₂）的正确装配是 `Σ ε_cell·cot_cell`。两侧
  三角形形状不同（介质网格粗）时误差可达数十百分点，作用于沿界面的横向耦合，
  正是场板边缘的电场集中区。2 V 控制在修 H1 后仍剩 P95 1.1 mV/最大 16 mV。
- 现状：无配置开关；Templates 的 `poisson_couple_profile:
  templates_ldmos_region_averagebox` 只接受 barycentric + Templates 载流子 profile。
- 判别：S2.2 先用 Python Poisson 复制算子做 2×2 变体（电荷体积 × 介电耦合）在
  原生态上评估界面残差；只有 Python 证明介电项主导时才提议新增诊断策略。

### H3 载流子输运边耦合（`assembled_delaunay_fallback` vs 原生 AverageBox）

- 机制：883 个严格钝角单元、日志 882 次边聚合修复。产生载流子在热点附近的
  通量分配、沿表面的横向输运受耦合系数影响；IIC 中电子/空穴积分差异号相反
  （−0.20% / +0.44%）说明输运分布确有差异。2 V 原生耦合控制
  （`carrier_transport_couple_profile: slot_ldmos_external_averagebox`，
  `build/slot_element_ialmob_20260921/budget/native_transport.csv`，36,681 条边）
  电流几乎不变，但 2 V 不代表击穿附近的反馈响应。
- 限制：`DCSweep.cpp`/`vela_example_runner.cpp` 在 `impact_ionization.model != none`
  时拒绝该 profile，需小幅解除守卫（仅探针/IIC 诊断，附测试）或用 Python 回放。
- 判别：S1 归因中输运通量行的权重；S2.3 变体。

### H4 Vela 38 V 自洽态收敛度不足（解释 δK_eff/K_eff 随偏压增长）

- 事实：38 V 最大行比 8.59e-4 / 7.94e-4（门槛 1e-3），35 V 只有 2.4e-5 / 1.4e-5。
  近奇异算子把局部残差放大约 `M` 倍，8.9×8.6e-4≈0.8%，与 1.15%−0.73% 的增量
  同量级；不能解释 35 V 的 0.73%。
- 判别：S0.1 从已接受态做一次原始 Newton 步，给出源项收敛误差的**局部估计**
  （非严格上界），并用第二步的收缩比作佐证。

### H5 原生参考自身精度

- 事实：原生 Digits=5、ErrRef(e/h)=1e10 cm⁻³ 远大于雪崩区载流子密度、
  RHSMin=1e-10；但重复求解电流变化 9.3e-11，端口 KCL 不闭合 1.15e-5，
  提示影响 ≪1%。
- 判别：S4.1 收紧判据重算 38 V 原生点（虚拟机），预期变化 <0.05%。

### H6 接触内建势/边界模型差（低概率）

- 判别：S0.5 比较接触节点 ψ、φn、φp 差与 S1 边界行归因权重。

## 3. 定位方法核心：冻结交叉残差 + 一阶伴随归因

跨求解器总源差先作精确拆分：

```text
S_V(x_V) − S_S(x_S) = [S_V(x_V) − S_V(x_S)]  +  [S_V(x_S) − S_S(x_S)]
                      └── 状态差（伴随估计对象）──┘   └── 冻结源核差（38 V 已测 +0.000157 pp）──┘
```

伴随只估计第一项。记 Vela 自洽态 `x_V`、原生导入态 `x_S`，`F_V` 为 Vela 残差，
`J_V` 为 `x_S` 处的 Vela Jacobian，`g = ∂S/∂x`，`J_Vᵀ w = g`，则一阶有

```text
S_V(x_V) − S_V(x_S) ≈ gᵀ(x_V − x_S) ≈ wᵀ[F_V(x_V) − F_V(x_S)]
```

只有在 `wᵀF_V(x_V)`（本地终态剩余残差的加权贡献，由 S0.1 用同一 `w` 直接计算）确认
足够小时，才简化为

```text
ΔS_lin = −wᵀ r,   r = F_V(x_S),   a_i = −w_i r_i,   Σ_i a_i = ΔS_lin
```

`r = F_V(x_S)` 包含 Vela 与原生的算子差（源项、输运、Poisson、复合、边界各项，不预设
哪一项为零），但也包含状态转换误差（密度重算、导出精度）与原生自身的未收敛残差，
不能无条件称为纯算子差；噪声底与原生 KCL 不闭合用于给出这两部分的量级。`a_i` 是
第 i 行（节点×方程）对状态差的一阶贡献，可按方程块、节点类别、空间窗口和残差项
类型（通量/复合/雪崩/规范/边界）汇总。近奇异的 `J_V` 正是要捕捉的放大机制；该方法
不依赖“单一主模态”或固定场假设，但依赖线性化有效性与伴随近似误差的控制（伴随
残差加权的依据及非线性余项讨论见 Pierce & Giles 的伴随误差修正工作），因此必须先过
下列校验。全部量都可由现有工具获得，无需改 C++：

- `simulation_type: newton_step_probe`（`src/tools/vela_example_runner.cpp`，
  `NewtonSolver::evaluateStep`）在给定状态上装配 J、r，解 `J δx = −r`，输出
  逐节点 `psi/phin/phip_residual`、`delta_*`、`trial_*`。探针配置中把
  `quasi_fermi_update_limit_V`、`quasi_fermi_update_limit_minority_V` 和 `max_update`
  都设为 `0`（`clampQuasiFermiStep` 在两个上限 ≤0 时直接返回，`max_update>0`
  才生效），并核对输出 `raw_step_norm == step_norm`（否则 δx 被截断并触发
  Poisson 重校正）。
- 环境变量 `VELA_LINEAR_CAPTURE_DIR=<新目录>` 让 `LinearSolver` 把该次线性
  系统写为 `input_0.bin`（`src/solver/LinearCapture.cpp`，只捕获 rows≥30000 且
  rows%3==0 的系统，每进程最多 4 个，拒绝覆盖）。格式：`"VELALU02"`、u64 index×2、
  u64 rows/cols/nnz、int32 outer[cols+1]、int32 inner[nnz]、f64 values[nnz]、
  u64 len + f64 b[len]（= −r）、u64 len + f64 x[len]（= 原始 δx）。Python 用
  `scipy.sparse.csc_matrix((values, inner, outer), shape=(rows, cols))` 读取。
- **符号与单位校验（先于一切归因）**：同一状态、同一配置下分别以
  `coupling_mode: self_consistent` 和 `postprocess_only` 捕获 `J_SC, r_SC` 与
  `J_IIC, r_IIC`；两者只差雪崩源项及其 `local_ad` 导数。源项以 `−S_i`（内部单位）
  同时进入电子行与空穴行，因此

  ```text
  S = −c · 1ᵀ (r_SC − r_IIC)[电子自由 Si 行]
  ```

  `c` 为内部单位到 A/µm 的换算系数：用 observable 探针的源积分（已与原生匹配到
  2 ppm）反算 `c`，要求 `c>0`（得到负值说明符号约定与此处相反，必须先核对装配
  代码再继续），并与 `comparison.json → avalanche_source.term_to_A_per_um` 交叉核对；
  不一致时以反算值为准并报告差异。用空穴行重复，要求得到同一 `S` 与 `c`。
  随后 `g = −c · 1ᵀ (J_SC − J_IIC)[电子自由 Si 行, :]`（单位 A/µm 每内部势单位，
  待求量为捕获系统的内部变量，`potential_scale` 已包含在 J 中），空穴行重算一致；
  于是 `wᵀr` 直接是 A/µm。`g` 的独立差分校验见 S1 第 2 步。
- 伴随解：`w = spsolve(J_SC.T.tocsc(), g)`（或 `splu`）；41,643 阶稀疏 LU 秒级。
- **工具自身的五项校验（全部在捕获系统的内部变量空间完成，只在最后报告时换算为 A/µm 与 V）**：
  1. 两个线性系统的**归一化后向误差**（不用绝对范数阈值）：

     ```text
     η_primal  = ‖Jδx + r‖ / (‖J‖·‖δx‖ + ‖r‖)        （δx 为捕获文件中的 x）
     η_adjoint = ‖Jᵀw − g‖ / (‖Jᵀ‖·‖w‖ + ‖g‖)
     ```

     ‖·‖ 用 2-范数，`‖J‖` 用稀疏可算的 `‖J‖_1` 或 `‖J‖_∞` 作上界代替（报告所用范数）。
     直接 LU 的期望量级为 O(ε_mach·增长因子)；验收 `η ≤ 1e-12`，否则说明矩阵读取或近奇异
     求解有问题，先修正。同时报告 `cond` 估计（`scipy.sparse.linalg.onenormest` 于 `J⁻¹`）。
  2. 伴随恒等式 `gᵀδx = −wᵀr`，采用混合绝对/相对判据，避免接近零或强抵消时误判：

     ```text
     |gᵀδx + wᵀr| ≤ 1e-8 · (|gᵀδx| + |wᵀr|)/2 + 16·ε_mach · (Σ_i|g_i δx_i| + Σ_i|w_i r_i|)
     ```

     并报告两个内积各自的抵消比 `Σ_i|g_iδx_i| / |gᵀδx|`、`Σ_i|w_i r_i| / |wᵀr|`（后者即 `C`）。
  3. 总源方向差分：取 2～3 个方向 `d`（含 `δx/‖δx‖` 和热点节点单位向量），用 observable
     探针计算 `[S(x+εd) − S(x−εd)]/2ε` 与 `gᵀd` 比较（目标相对误差 ≤1e-3，两个步长确认
     差分稳定）。
  4. Jacobian–残差一致性：同一组方向上用 `newton_residual_probe` 计算
     `[F(x+εd) − F(x−εd)]/2ε` 与 `Jd` 比较（目标 ≤1e-4）。这同时区分两类问题：代码
     Jacobian 与自身残差导数不一致只影响求解路径与收敛质量（若残差严格收敛到同一根，
     不改变该离散方程的根），而源函数在邻域的真实导数差才会改变反馈响应（由反向冻结
     回放检验）。
  5. 正负贡献与抵消：报告 `Σ_{a_i>0} a_i`、`Σ_{a_i<0} a_i`、`Σ|a_i|` 与 `C`（下条）。
  方程分组（方程块、节点类别、空间窗口）在 35 V 首次归因时写成分组文件并冻结，38 V 与
  所有变体沿用同一分组，否则逐行归因不具可比性。
- **归因报告格式**：对任意行组 G 同时给出有符号份额 `Σ_{i∈G} a_i / ΔS_lin`、绝对份额
  `Σ_{i∈G}|a_i| / Σ_i|a_i|`，以及该组噪声底 `Σ_{i∈G} a_noise,i`；全局报告抵消指数
  `C = Σ_i|a_i| / |ΔS_lin|`。`C` 大时有符号份额易被抵消误导，第 4 节的阈值同时要求
  两种份额。
- **线性化有效性门**：`ΔS_lin`（含 `wᵀF_V(x_V)` 修正项）与实测状态差 `S_V(x_V) − S_V(x_S)`
  （35 V ≈+0.887%、38 V ≈+9.005%，均已扣除 ppm 级冻结源核差）比较，35 V 偏差 ≤10%、
  38 V ≤30% 视为一阶归因可用；否则以 35 V 为主归因偏压、38 V 作趋势确认。另用
  `trial_*` 列构造试探态（以及 0.5×δx 的半步态）并以 `newton_carrier_term_probe` 计算
  非线性源项，估计非线性余项。一阶预测与实测相符只证明线性化可用，不证明单模态模型。
- **噪声底**：把原生态按导出精度扰动（密度重算 4.8 ppm、QF ±1 ULP）得到 `r_noise`
  与 `a_noise`，要求 `|Σ a_noise| ≪ |ΔS_lin|`，并按组报告。原生 KCL 不闭合 ~1e-5 为
  原生自身残差的量级参考。
- **基础源分布的放大贡献 `ΔS_src`**（对应 1.1 节的 `cᵀRδb` 项）：以 IIC 逐节点源差
  `δs0_i = S_i^{Vela,IIC} − S_i^{原生,IIC}`（A/µm，共同 12,306 行）作为等效注入扰动，
  电子行与空穴行同时取 `−δs0_i/c`，则
  `ΔS_src = −Σ_i (w_{φn,i} + w_{φp,i}) δs0_i / c`。报告 `ΔS_src/ΔS_lin` 及其有符号/
  绝对份额。`δs0` 本身是两边输运/静电算子差的结果，`ΔS_src` 是“分布效应”的量级
  检查，不与 `a_i` 相加成因果预算。
- **反向冻结回放**（检验源函数在 Vela 态邻域是否一致，回应“源值匹配 ≠ 导数匹配”）：
  用 `scripts/make_sentaurus_tdr_field_substitution.py`/`perturb_sentaurus_tdr_node_field.py`
  把 Vela 35 V（再 38 V）自洽态 ψ/φn/φp 写入原生 TDR，原生仅 Load/Plot/CurrentPlot
  （无 Newton），把原生核在 Vela 态上的逐节点源与 Vela 自身装配源比较。仍在 ppm
  级 → 源函数（含导数）差可排除；出现百分级差 → 源核仍是候选，S1 归因需把源行
  差单列。属虚拟机任务，可与 S1 并行。

## 4. 执行计划

统一环境：worktree 根目录，PowerShell：

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
$env:VELA_LINEAR_SOLVER = "umfpack"
# Release 运行器：build-release/vela_example_runner.exe
# 先核对 SHA-256 == fd364d29b0c9327aa5a8fba2f502a5dba42b961e4d8ec0b589cae92083c2ee0b（s25iic 记录）；不同则记录新哈希并说明是否重编译
```

所有输出写入新目录 `build/s26loop/<step>/`，每个目录附 `manifest.json`
（运行器/配置/状态/网格/材料/几何 SHA-256、命令行、退出码、wallclock）。
两套电位约定不要混用：生产扫描/已接受状态用 `potential_origin_V: 35`（物理电压 =
求解器电压 + 35 V）；`build/s25iic/frozen_*` 的冻结探针配置用 `potential_origin_V: 0`
且接触偏压直接为物理值（drain = 35/37/38，其余 0）。报告统一用物理 Vds。

### S0 零代码快速检查（成本低，全部只读或探针）

S0.1 **剩余收敛误差估计（H4）**。对 Vela 已接受的 35/37/38 V 自洽态
（`build/s25to37|s25to38/*/qualification_audit.json` 指向的 `final_state.h5`）
和 38 V IIC 态（`build/s25iic/local_iic_density/`），以生产配置（只改
`simulation_type`、`output_csv`，并把三个步长上限置 0）运行 `newton_step_probe`；
在 Vela 自洽态 `x_V` 上同时设置 `VELA_LINEAR_CAPTURE_DIR`，并以 `postprocess_only`
再跑一次，得到 `J_SC(x_V)`、`J_IIC(x_V)`、`r_SC(x_V)`、`r_IIC(x_V)`（供 S4.2 的同基态
灵敏度与 S1 的另一端点线性化使用）。用 `trial_psi/phin/phip` 构造试探态（参照
`build/s25iic/frozen_native.py` 用 `scripts/state_archive.py` 写 HDF5），再用
`newton_carrier_term_probe` 计算源积分；沿同一方向另取 0.5×δx、0.25×δx 的部分步，
检查源变化是否随步幅线性（线性预测是否稳定）；从整步试探态再做第二步，记录两步
残差块范数的收缩比。一次原始 Newton 步来自局部一阶展开，给出的是**剩余收敛误差
估计**；成为严格误差界还需逆算子与非线性余项的控制，此处不作该声称。S1 得到伴随
向量 `w` 后，再报告 `wᵀF_V(x_V)`（同一权重下终态残差对总源的加权贡献），作为进入
`−wᵀF_V(x_S)` 简化式的前提。输出 `ΔS_conv/S`、部分步线性性、两步残差块范数、
`carrier_row_convergence.csv` 中热点节点行比。判据：`|ΔS_conv/S| < 0.1%`、部分步
线性、且第二步残差再收缩 ≥10× → H4 不作为主因（不表述为“排除”）；若 38 V 达
0.3%～0.8% 则 H4 解释增量部分，并在 S3 中对 38 V 追加更多 Newton 更新（不改门槛）
后重新取值。

S0.2 **实际驱动场、α 比与静电场对照（H1/H2 直接物理量）**。生产配置的雪崩与 HFS
驱动力是 `quasi_fermi_gradient`（`cell_gradient`，接触处 `contact_electric_field_fallback`，
界面 `ialmob_boundary_field_projection`），因此 α 的比较必须用两边状态的**实际电子/
空穴驱动场**：对每个 Si 单元用同一线性三角形梯度算子分别求 |∇φn|、|∇φp|
（Vela 自洽态 vs 原生导出 QF，同一物理参考），按生产规则处理接触/界面单元，代入
同一 van Overstraeten 参数得 `α_n(F_n)`、`α_p(F_p)` 的 Vela/原生比值图。由 ψ 恢复的
电场 |∇ψ| 及其 α 只作为 Poisson 差异诊断另列，不能替代生产驱动场。参考对齐、
Si 节点集合与梯度算子复用 `build/s25to38/compare.py` 与 `prepare_native.py`。列出：
热点单元、热点正上方界面节点列（8552…8310）、y∈[3.5,3.7] 界面区、全 Si 的分布；
沿穿过热点的若干电场线用两边场分别计算电子/空穴电离积分并给出相对差（诊断量，
非原生 IIC 路径积分）。注意 `fields_*V.npz` 中的 `eIonIntegral/hIonIntegral`
（0.1287/0.0984）三个偏压完全相同，是早期 Plot 的静态字段，禁止当作 35～38 V
的路径积分。判据：热点及其上方 100 nm 内 `α_p` 比 ≥1.003 且电离积分相对差与
`δK_eff/K_eff` 同号同量级 → 驱动场差是主因候选；α 比 ≈1.000 → 转向输运（H3）。

S0.3 **IIC 热点与界面节点的因子分解**。用 Vela IIC 态和原生 IIC 态（`frozen_iic_h5_*V`
或 `native/fields_*V.npz` 加导出的 n/p/QF），在热点单元以及 10490/10492/6626 相邻
单元上分别比较 α、`D=2n_min n_max/(n_min+n_max)`、元素迁移率、|∇φ_qf| 及源分配份额
的比值，交叉核对 S0.2。判据：确定 1.2 节热点 +0.2%～+0.4% 与界面节点 +31%～+49%
的局部差异分别由哪些因子承担；在此之前不将界面节点差异称为电场缺陷。

S0.4 **网格与界面审计**。用 `scripts/audit_slot_ldmos_delaunay_mesh.py` 思路统计
热点窗口（x<0.3, y∈[1.9,2.4]）与 y∈[3.4,3.8] 表面区内的钝角单元、聚合修复边、
界面单元及其两侧三角形形状比（cot 比），为 H2/H3 提供先验。

S0.5 **接触/边界差（H6）**。接触节点 ψ、φn、φp 两边差的统计；判据 ≤2 mV。

S0.6 **固定电场假设的定量检查**。同一引擎内比较自洽态与 IIC 态：Vela 38 V 自洽 vs
Vela 38 V IIC 的 ψ 差（全 Si RMS/最大，热点窗口、y∈[3.5,3.7] 界面区分列），原生同样
比较（`frozen_h5_38V` vs `frozen_iic_h5_38V`）；同时给出热点区 `q(n+p)` 与 `q|N_D−N_A|`
之比。判据：热点区 |δψ| ≤0.1 mV 且电荷比 ≤1e-3 → 固定场近似成立，1.1 节的标量模型
可作解释框架；否则伴随分析仍有效，但 `K_eff` 只作描述性指标使用。

### S1 冻结交叉残差与伴随归因（成本中，Python + 探针，不改 C++）

输入：`build/s25iic/frozen_h5_{35,37,38}V/native_state.h5`（原生自洽态）与
`frozen_iic_h5_*V/`（原生 IIC 态）及各自 `config.json`。

1. 以 `frozen_h5_38V/config.json` 为模板生成两份探针配置（删去 `carrier_term_probe` 键）：
   `simulation_type: newton_step_probe`，`coupling_mode` 分别为 `self_consistent`
   与 `postprocess_only`，`quasi_fermi_update_limit_V`、`quasi_fermi_update_limit_minority_V`、
   `max_update` 均置 `0`，`output_csv` 指向新目录；
   分别在独立进程、独立 `VELA_LINEAR_CAPTURE_DIR` 下运行，得 `J_SC, r_SC, δx_SC`
   与 `J_IIC, r_IIC`。运行前后核对状态文件哈希不变。
2. 先做第 3 节的符号与单位校验（反算 `c>0`，电子行/空穴行一致，与 `term_to_A_per_um`
   交叉核对），再计算 `g`、`w`、`a_i`，并完成第 3 节的五项工具校验（两个线性系统残差、
   `gᵀδx = −wᵀr`、总源方向差分、Jacobian–残差一致性、正负/绝对贡献）。`g` 的独立差分
   校验：选 3 个热点节点（如 8237、8102、8402），分别对其 ψ/φn/φp 做 ±1e-6 V 中心差分，
   用 observable `newton_carrier_term_probe`（`coupling_mode: postprocess_only`，
   `solved_equation_terms: false`）求总源积分（A/µm）的差商，与 `g` 的对应分量比较
   （目标相对误差 ≤1e-3；差商自变量需换成捕获系统的内部势单位，即除以
   `potential_scale`）。同时用 `newton_carrier_term_probe`（`solved_equation_terms: true`，
   以实际表头为准）和 `newton_poisson_term_probe` 把 `r_i` 拆成通量/复合/雪崩/规范/
   边界与 Poisson 各项，得到项级归因 `a_i^{term} = −w_i r_i^{term}`（先核对各项之和与
   探针残差逐行一致）。同时计算 `ΔS_src`（第 3 节），并用同一 `w` 对 Vela 已接受态的
   终态残差计算 `wᵀF_V(x_V)`（接 S0.1），确认其相对 `ΔS_lin` 可忽略后才使用简化式。
3. 汇总维度（首次归因时写成分组文件并冻结）：方程块（ψ/φn/φp）；节点类别（Si 内部、
   Si/绝缘体界面、接触、多晶硅、介质规范行）；空间窗口（热点 y∈[2.0,2.3]；SLOT 栅侧
   边缘 y∈[3.5,3.7]；栅/源结 y>3.8；漏侧 y<0.5；其余）；钝角单元邻接与否。每组同时报告
   有符号份额、绝对份额与噪声底，全局报告抵消指数 `C`。输出逐行 CSV（node、block、x、
   y、class、`r_i`、`w_i`、`a_i`、项拆分）、`summary.json`、按 |a_i| 排序的前 50 行及累计占比。
4. 有效性：`ΔS_lin` vs 实测；试探态非线性一步源项（用 `trial_*` 列写状态后跑
   observable 探针，电子/空穴发起部分可分开读出）；噪声底。再在原生 IIC 态
   （`frozen_iic_h5_*V`）上重复：此时步长用 `J_IIC`，`g` 仍由同一态上的
   `J_SC − J_IIC` 得到，预测 `ΔS_IIC,lin = −gᵀ J_IIC⁻¹ r_IIC`；M≈1，一阶预测应复现
   IIC 差（38 V 总 −0.059%，试探态读出的电子 −0.20%/空穴 +0.44%），并比较
   IIC 态与自洽态的输运通量行残差模式是否一致（算子差应与雪崩状态无关）。
5. **顺序约束**：先在 35 V 完成第 1～4 步并通过有效性门（符号/单位、`g` 差分、
   `ΔS_lin` 偏差 ≤10%、噪声底），再在 38 V 重复；38 V 与 35 V 的归因分布对比用于解释
   `δK_eff/K_eff` 增长。

判据（决策树；“份额”指第 3 节定义的有符号份额，并要求绝对份额同向支持；每条
结论附 `C` 与该组噪声底）：

- Poisson 行（尤其界面节点）有符号份额 ≥50% 且绝对份额 ≥30% → 进入 S2.1/S2.2，
  H1/H2 为主。
- 输运通量行有符号份额 ≥50% 且绝对份额 ≥30%，集中在热点/钝角单元 → 进入 S2.3，
  H3 为主。
- `ΔS_src` 有符号份额 ≥30% → 基础源分布效应不可忽略，S2 各变体需同时报告其对 IIC
  逐节点源差的影响。
- 雪崩源行（`r_SC − r_IIC` 部分）或反向冻结回放出现百分级差 → 源核重新列为候选，
  先排查源核再进入 S2。
- `C > 5` 或各组有符号份额相互抵消 → 不按份额下结论，改按 S2 各变体的 `ΔS_lin`
  实际下降量排序。
- 分布弥散、无单一族达到阈值 → 记录为多因素离散差；S2 三个变体逐个测 `ΔS_lin`
  下降量。
- `ΔS_lin` 与实测严重不符且 35 V 也不符 → 复核 `c`、`g` 的符号/单位、边界行处理、
  状态转换，不得进入 S3。

### S2 单因素冻结变体（成本低～中）

S2.1 **H1：Poisson 电荷体积策略**。在 S1 的探针配置中加入
`"discretization": {"poisson_charge_volume_policy": "material_local_mixed_voronoi"}`
（要求 `node_volume_policy: mixed_voronoi`，本配置已满足；参考 2 V 控制
`build/slot_energy_reference_20260921/budget/material_charge/`）。重跑 S1 第 1～3 步，
报告：界面节点（8552…8310、10490、10492 及全部 Si/绝缘体界面节点）ψ 行残差
变体前后之比；`ΔS_lin` 变体前后（有符号与绝对份额均列）；`w` 的变化。判据：界面
ψ 行残差下降 ≥5× 且 `ΔS_lin` 下降 ≥50% → H1 确认进入 S3。

S2.2 **H2：介电耦合**。改造 `scripts/audit_templates_ldmos_g3_interface_charge_volume.py`
的 Python Poisson 复制算子为 SLOT 版本。算子等价性验收在 Vela 自洽态上逐行进行：
`|r_replica,i − r_Vela,i| ≤ 64·ε_mach·Σ_k|term_{i,k}|`（操作数舍入界，各项取自
`newton_poisson_term_probe` 输出），不得用接近零的残差本身作相对分母。通过后在
原生态上评估 2×2 变体：电荷体积 {global, material_local_mixed_voronoi} × 介电耦合
{edgeEpsilon 算术平均, 逐单元 ε·cot}。报告各变体界面节点残差 RMS/最大值（及各自
的操作数舍入界）。若逐单元 ε 使残差再下降 ≥3×，提议新增默认关闭的诊断策略
`discretization.poisson_permittivity_coupling: per_cell`（C++ 附加改动，需两材料
解析夹具测试，默认路径逐位不变），否则不改代码。

S2.3 **H3：输运耦合**。优先 Python：用 Vela 自洽/原生态、`native_transport.csv`
与 Vela 网格耦合，回放电子/空穴 SG 通量散度（形式见 `CoupledDDAssembler.cpp`
残差边循环、`sg_edge_flux_probe` 输出可作校验），比较两套耦合在热点/钝角单元
的通量散度差与 S1 输运行残差的相关性。若需实际探针，可在 `DCSweep.cpp`
与 `vela_example_runner.cpp` 中把 `slot_ldmos_external_averagebox` 的雪崩拒绝
放宽到 `ialmob_element_density_grad_qf`（附测试，默认不变），重跑 S1 第 1～3 步。

### S3 自洽确认（成本中～高，逐点 300 s 预算；前提：S1 在 35 V 的有效性门已通过，且 S2 至少一个变体达到判据）

1. 配置 = `base.json` + 获胜变体（首选 S2.1 的 Poisson 策略），不改任何门槛
   （`abstol=1e-9`、`reltol=1e-8`、`eps_row=0.001`、全局闭合 0.01、
   `stall_residual_floor=1e-8`）。
2. 35 V：以 Vela 既有 35 V IIC 态和自洽态为初值（不是原生态），分别以
   `postprocess_only` 与 `self_consistent` 做固定偏压重平衡（`newton_solve_from_state`
   或 `scripts/run_slot_ldmos_loadline_continuation.py` 的单点模式，沿用
   `build/s25iic/run_local_iic_density.py` 的验收流程：完整 Newton、每类
   12,634 有效行、全局闭合、独立审计）。计算 `M_Vela'`、`δK_eff/K_eff'`、ψ RMS/最大差、
   10490/10492 IIC 比、热点 SC 比、IIC 电子/空穴积分差。
3. 判据（五个维度同时看，不以 `K_eff` 接近原生作为根因成立的单一条件）：
   - 源差：`δK_eff/K_eff` 由 0.73% 降到 ≤0.25%；IIC 与 SC 的总/电子/空穴源积分差与节点 L1
     均不劣化，10490/10492 的 IIC 比与热点 SC 比趋近 1；
   - 端电流：IIC 电流差 |Δ| ≤0.15%（不劣于 −0.076%），SC 电流差由 +0.240% 下降；
   - 势场：Si ψ RMS 差 ≤2 mV，界面节点（含 10492）最大差明显下降，实际电子/空穴驱动场
     在热点的比趋近 1（S0.2 同一算子）；
   - 局部残差：完整 Newton、每类 12,634 有效行、全局闭合全部通过，`wᵀF_V(x_V')` 仍可忽略；
   - 配置一致性：变体配置与 `base.json` 的字典差仅为预期键（字典相等性检查）。
   五项均满足 → 35 V 根因定位成立；再做 38 V（必要时 37 V），期望 `δM/M` 由 +9.07% 降到
   ≤+3%。若 38 V 残余仍明显大于 `M·δK_eff(35 V)` 的标量重述，结合 S0.1 追加 Newton 更新
   后复测，把收敛贡献、非线性余项与算子贡献分开报告。
4. 若 S2.1 不足而 S2.2/S2.3 有效，按同样流程叠加；每次只加一个因素，记录交互项。
5. 所有新点必须同时保留失败/超时证据；不得覆盖 `s25iic`、`s25to37`、`s25to38`。

### S4 原生侧对照与标量近似检验（S4.1/S4.3 可与 S1～S3 并行；S4.2 依赖 S0.1 在 `x_V` 的捕获与一次新伴随解，S4.4 依赖 S1 的捕获）

S4.1 **原生精度控制（H5）**：从原生 38 V 保存态重算，仅改 `Math` 中
`Digits=7`、`ErrRef(electron)=ErrRef(hole)=1e6`、`RHSMin=1e-14`（其余逐字节相同，
新输出前缀），比较电流与源积分变化；预期 <0.05%。若不收敛，仅放宽
`Iterations`，不改物理。
S4.2 **α 缩放灵敏度（检验标量近似；先用伴随解一次线性计算，再少量重收敛核对）**。
若缩放 `s` 同时作用于基础源与反馈，标量近似下 `S_SC(s) = s·S_IIC/(1 − sK)`，两种观测量
必须分开：

```text
d ln S_SC / d ln s |_{s=1} = 1 + K·M = M（总源灵敏度）,
 d ln[S_SC(s)/(s·S_IIC)] / d ln s |_{s=1} = K·M（倍增因子灵敏度）
```

更适合非线性离散方程的是直接的离散灵敏度：若 `F(x, s) = F_0(x) − s·Q(x) = 0`（`F_0 = r_IIC`
无反馈残差，`Q = r_IIC − r_SC` 为源向量，内部单位），则在 `s = 1`

```text
J · dx/d ln s = Q,      d ln S / d ln s = 1 + gᵀJ⁻¹Q / S = 1 + wᵀQ / S
```

这是方程可微、`J` 可逆且缩放确实只按上式进入时成立的局部灵敏度，**不要求单模态标量
假设**（标准离散伴随方法，参见 Giles & Pierce）。但它由隐函数定理在**方程的根**处导出：
`J`、`w`、`Q`、`S` 必须来自**同一合格基态**。S1 在原生导入态 `x_S` 上的 `w`（那里
`F_V(x_S) ≠ 0`）不能拿来预测 Vela 合格自洽态 `x_V` 的 α 灵敏度。执行：

1. 在 Vela 合格自洽态 `x_V` 上计算：用 S0.1 已捕获的 `J_SC(x_V)`、`J_IIC(x_V)`、`r_SC(x_V)`、
   `r_IIC(x_V)`，以同一 `c` 得 `g_V = −c·1ᵀ(J_SC − J_IIC)(x_V)[电子行, :]`，**额外做一次伴随解**
   `J_SC(x_V)ᵀ w_V = g_V`（非“零额外求解”），`Q_V = r_IIC(x_V) − r_SC(x_V)`，`S = S_V(x_V)`，
   则 `σ = 1 + w_VᵀQ_V/S`（35 V 与 38 V），与 1.1 节的 `M` 比较；差 ≤10% 则标量近似可作局部
   解释框架，否则 `K_eff` 只作描述性指标。前提：用 observable 探针确认
   `impact_ionization.A_scale = s` 在固定状态上使 `ialmob_element_density_grad_qf` 源向量逐节点
   精确变为 `s·Q`（两类载流子、高/低场分支均成立），否则 `F_0 − sQ` 的形式不成立，需改用
   实际的 `∂F/∂s`。`w_V` 同时可用于 S1 的“另一端点线性化”：`w_Vᵀ[F_V(x_V) − F_V(x_S)]` 与
   `x_S` 端的 `wᵀ[F_V(x_V) − F_V(x_S)]` 两个估计包夹非线性余项。
2. 少量重收敛核对：Vela 以 `A_scale = 1 ± 0.01` 在 35 V（再 38 V）固定偏压重平衡，分别取
   SC 与 IIC，中心差分得 `d ln S_SC/d ln s` 与 `d ln M/d ln s`，与 `σ`、`M`、`K_eff·M` 比较。
   ±1% 对照与 `σ` 的差**不能单独称为非线性余项**：它还包含重收敛误差（以 S0.1 的
   `ΔS_conv/S` 为参考）、差分截断误差与可能的分支切换（如 α 高/低场分支、选点、线搜索
   路径）。仅当差异明显（如 >3%）时，再补 `±0.5%` 检查步长依赖；不预先增加其他重算。
3. 原生侧（可选）：van Overstraeten `a` 系数乘 `s`，其余不改，同样测两种灵敏度，与 Vela
   对照。
S4.3 若 S2.1 确认 H1，可在原生侧用 `Plot` 导出场板边缘界面节点电荷/势对照，
补充空间证据；不要求。
S4.4 **环路谱诊断（伴随定位之后，可选，纯 Python）**。在**同一状态、同一变量缩放**下，用
捕获的矩阵定义线性化环路算子 `L = J_IIC⁻¹ Q_x`，`Q_x = J_IIC − J_SC`（源项 Jacobian，
`J_SC = J_IIC − Q_x`），此时 `J_SC⁻¹ J_IIC = (I − L)⁻¹` 成立。以稀疏 LU 为 `LinearOperator`，
用 `scipy.sparse.linalg.eigs` 同时求：(a) `L` 的最大模特征值 `λ_max`；(b) **最接近 1 的特征值**
`λ_near1`（对 `(I − L)⁻¹ = J_SC⁻¹ J_IIC` 取最大模特征值 `μ`，`λ_near1 = 1 − 1/μ`）；两者均与
`K_eff` 对照。预解式范数 `‖(I − L)⁻¹‖₂`（`svds` 或幂迭代于 `AᵀA`）一般**不能**用 `1/(1 − λ_max)`
代表：即使对正规矩阵，二范数对应的也是 `1/min_i|1 − λ_i|`。因此报告 `‖(I − L)⁻¹‖₂` 与
`1/min_i|1 − λ_i|` 的比值仅作诊断：比值 ≫ 1 提示非正规放大可能存在，但不能单独据此证明。
两边引擎都在原生态（或 Vela 态）上用 Vela 算子评估（原生算子不可直接得到），因此这是
Vela 环路的谱，只用于判断“单模态/正规”近似是否合理，不用于跨引擎归因。

### S5 报告、策略与记忆

- 新报告 `docs/validation/slot_ldmos_iic_loop_gain_localization_<日期>.md`：
  沿用本仓库格式（结论先行、配置/单位/后端、门槛未改、证据路径、SHA-256），
  并在 `docs/README.md` 与 `slot_ldmos_bvds_validation.md` 增加链接。
- 若 H1 成立：提议把 `discretization.poisson_charge_volume_policy:
  material_local_mixed_voronoi` 写入 **显式 SLOT 配置 profile**，不改全局默认；
  按 AGENTS.md 要求重新资格化 0→38 V 链（或至少 2/10/35/37/38 V 锚点 + 外电阻
  起点），并在 `tests/regression` 增加该策略的配置契约测试。
- 更新 `slot_ldmos_status_and_handoff` 的“Poisson/输运几何”一行。

## 5. 约束与禁止事项

- 不修改生产求解器默认、物理参数、任何验收门槛；不用调参使 M 逼近原生。
- 原生保存态只用于冻结诊断，不作 Vela 求解初值。
- 探针中关闭 QF/整体步长上限仅用于取原始 Newton 步，不得写入器件配置。
- C++ 改动仅限附加诊断（策略开关、守卫放宽），默认路径逐位不变，附 Catch2
  测试，Release 构建，运行相关 CTest（如 `pn2d_config_templates`、
  `ialmob_element.*protocol`、`^Element avalanche`）并报告；未运行的检查明确写出。
- 每个数值结论标明状态来源（Vela 自洽/IIC、原生自洽/IIC）、偏压、单位
  （源项 A/µm、密度 m⁻³ 或 cm⁻³、坐标 µm）和后端（UMFPACK/SparseLU）。
- 不把 `ΔS_lin` 当作精确因果分摊；报告线性化有效性与噪声底；任何份额结论同时给出
  有符号份额、绝对份额与抵消指数 `C`。
- 术语：`K_eff` 只称“等效标量反馈指标”，不称主特征值或环路增益；在 S1 归因与 S3
  确认完成前，不得写“已证明约 1% 的环路增益差导致约 9% 的源项差”，只能写“与等效
  标量反馈差被放大的图像相容”；界面节点源比在 S0.3 分解前只称“局部差异”；一次
  Newton 步的结果只称“剩余收敛误差估计”；`r = F_V(x_S)` 称“交叉残差”，不称纯算子差。
- 验收同时看源差、端电流、势场、局部残差与配置一致性（S3 第 3 步），不以 `K_eff` 接近
  原生作为根因成立的单一条件。
- 符号、单位换算与变量缩放（`potential_scale`、`residual_scales`）在 `J`、`r`、`g`、`δx`、`w`
  中保持同一约定；方程分组一旦冻结不再改动。
- 不将本计划任何中间结果解释为 BVDS 偏差；BVDS 仍以 `|Id|=1e-7 A/µm` 判据比较。

## 6. 交付物与报告格式

每阶段一个 `build/s26loop/<step>/summary.json` 与一节 Markdown，包含：

| 项 | 内容 |
| --- | --- |
| 状态 | 偏压、来源、哈希 |
| 主结果 | S0.1 `ΔS_conv/S`、部分步线性性与两步收缩比；S0.2 实际驱动场 α 比与电离积分差；S0.6 固定场检查；S1 `c`、五项工具校验结果、`wᵀF_V(x_V)`、`ΔS_lin`/实测、噪声底、`ΔS_src`、`C`、前 50 行及各组有符号/绝对份额；S2 各变体残差与 `ΔS_lin` 下降；S3 五维度验收（源差、端电流、势场、局部残差、配置一致性）、`M`、`δK_eff/K_eff`；S4.2 `σ = 1 + w_VᵀQ_V/S`（同基态）与重收敛差分；S4.4 `λ_max`、`λ_near1`、`‖(I−L)⁻¹‖₂` 与 `1/min|1−λ_i|` 的比（诊断） |
| 判据结论 | 通过/未通过及依据 |
| 未完成/限制 | 明确列出 |

最终报告以第 1.1 节表格的形式给出修正前后的 `M`、`K_eff`、`δK_eff/K_eff`，并给出剩余
未归因预算（含噪声底与线性化偏差）。

## 附录 A：关键路径

- 分析对象：`docs/validation/slot_ldmos_iic_source_budget_2026-09-25.md`。
- 既有脚本：`build/s25iic/{integral_scope,frozen_native,frozen_iic_native,compare_frozen,compare_iic,run_local_iic_density}.py`；
  场/势对照与原生提取：`build/s25to38/{compare,prepare_native,audit}.py`、`build/s25to37/{compare,verify}.py`。
- 原生冻结态：`build/s25iic/frozen_h5_{35,37,38}V/native_state.h5`、`frozen_iic_h5_*V/`。
- 原生场：`build/s25iic/native/fields_{35,37,38}V.npz`、`native/manifest.json`。
- 原生 Si 顶点体积：`build/slot_element_ialmob_20260921/local_control_r2/native_geometry.json`。
- 原生输运耦合：`build/slot_element_ialmob_20260921/budget/native_transport.csv`。
- 生产配置：`build/s25to38/to38_o35/current/base.json`；单位换算
  `comparison.json → avalanche_source.term_to_A_per_um`。
- 探针实现：`src/tools/vela_example_runner.cpp`（`runNewtonStepProbe`、
  `runNewtonCarrierTermProbe`、`runNewtonPoissonTermProbe`、`runSgEdgeFluxProbe`、
  `runNewtonFeedbackSubstitutionProbe`）；`src/solver/NewtonSolver.cpp::evaluateStep`；
  `src/solver/LinearCapture.cpp`。
- 策略解析：`src/simulation/ConfigParsing.cpp`（`parsePoissonChargeVolumePolicy`、
  `applyCarrierTransportCoupleProfile`、`applyPoissonCoupleProfile`）；
  `include/vela/equation/AssemblerUtils.h::edgeEpsilon`；`src/equation/IalTransport.cpp`。
- 方法先例：`scripts/diagnose_pn2d_bv_exact_carrier_term_states.py`、
  `scripts/diagnose_pn2d_bv_loop_gain_sensitivity.py`、
  `scripts/audit_templates_ldmos_g3_interface_charge_volume.py`。

## 附录 B：线性捕获读取示例

```python
import numpy as np, scipy.sparse as sp
def read_capture(path):
    with open(path, 'rb') as f:
        assert f.read(8) == b'VELALU02'
        u64 = lambda: int(np.frombuffer(f.read(8), '<u8')[0])
        _idx0, _idx1 = u64(), u64()
        rows, cols, nnz = u64(), u64(), u64()
        outer = np.frombuffer(f.read(4 * (cols + 1)), '<i4')
        inner = np.frombuffer(f.read(4 * nnz), '<i4')
        vals = np.frombuffer(f.read(8 * nnz), '<f8')
        nb = u64(); b = np.frombuffer(f.read(8 * nb), '<f8')
        nx = u64(); x = np.frombuffer(f.read(8 * nx), '<f8')
    return sp.csc_matrix((vals, inner, outer), shape=(rows, cols)), b, x
# J, minus_r, delta_x = read_capture('build/s26loop/s1_38V_sc/capture/input_0.bin')
# 校验：np.allclose(J @ delta_x, minus_r) 且 delta_x 与 newton_step_probe CSV 的 delta_* 列一致（注意 potential_scale）
```
