# Templates/LDMOS G3 最大 gm 差异：后续调查方向与任务计划（WP3）

日期：2026-09-01
状态：提案（未执行）
修订：R1（2026-09-01）——合并外部评审四项修正：T1 状态集合合同化、
T2 守恒约束与正则化稳健性防护、ElementEdgeCurrent 表述限定、T4 变体状态再生协议
修订：R2（2026-09-01）——T4 增设 `2–10x` 不确定带、T6 增设保护门、
AreaFactor/33.5x 归一化措辞收敛为“排除面积/宽度归因 + 支持性解释”
修订：R3（2026-09-01）——修订后批准整改：2.2 改为“所需修正”表述（撤回
`delta_F` 直接测得声明）、T4 增设跨引擎物理资格门、2.3 饱和灵敏度加
上/下风限定、T3/T5 增设热点迁移保护门、T1 增设共线性防护与逐边灵敏度表、
T2 重定位为新工具开发（可行性/筛选）、停止分类更名
修订：R4（2026-09-01）——执行合同细化（批准前置）：T7 拆分为 T7a/T7b
（ledger 证据包不再依赖 T6）、T2 图定义冻结（节点/边域、跨域边、规范
定向）、T4 资格门量化为显式不等式并预注册机器可读 summary schema
范围：`Vd=0.1 V`、`Vg=1.0→1.1666667 V` 最大 gm 段，24.3076% 自洽 gm 差异的根因收敛与停止分类
分支：`codex/templates-ldmos-phase-a`，基线提交 `b61f2ba`

本文可独立阅读：第一、二节复述必要输入事实，第三节起为新分析与计划。

---

## 一、输入事实摘要（来自已冻结报告）

1. 31 点 G3 曲线全部收敛，除最大 gm 门（24.308% > 20%）外全部门通过；
   两引擎最大 gm 位于同一精确段 `Vg=1.0→1.1666667 V`。
2. 冻结算子闭合：Vela 算子作用于完整 Sentaurus 状态（SSS），两端点端口电流
   比均为 `1.000202`，段 gm 误差仅 `+0.0202%`；自洽（VVV）gm 为 `0.756924x`。
   差异归因于电子准费米势 `phin` 自洽反馈（`-0.108 dex`；psi `+0.0014`，phip `0`）。
3. 七个热点节点 `3721, 3974, 3973, 4091, 3727, 4021, 3771`（六个直接
   Si/介质界面节点 + 4021 为 Si 侧一环）：
   - Sentaurus 自身七节点电子 RHS 闭合到 `3.5e-15 / 4.4e-15 A`（相对端口电流 ~1e-9）；
   - Vela 在相同 SSS 状态上的七节点物理化残差为 `1.24e-9 / 6.98e-9 A/um`
     （相对端口电流 `9.9e-4 / 2.5e-3`）；
   - 残差相对入射边绝对通量和仅 `5.75e-5–2.56e-3`，即大通量相消后的结构化小量；
   - 七行符号成组相反，随段增长 `2.82–3.16x`（端口电流增长 `2.25x`）；
   - 七节点复合（SRH/Auger）项 `<=4.52e-23`，输运通量完全主导。
4. Jacobian 证据：node 3721 `phin` ±1 µV 单模扰动，Sentaurus iteration-0
   NewtonPlot 与 Vela 中心差分列余弦 `-0.999999831`，最小二乘标量 `-33.7806`
   （与未扰动基线 `33.46–33.59` 漂移 0.96%），去标量后相对 L2 仅 `0.058%`。
   注意：**只探测过这一个 phin 列，psi 列从未探测**。
5. `AreaFactor=1→2` A/B：端口电流精确 2x，但完整 `eContinuityRhs` 逐字节不变
   ⇒ 排除面积/宽度对 `~33.5x` 的归因；结合手册与其他 A/B，支持（但未直接
   证明）NewtonPlot RHS 为实现相关内部归一化的解释。
6. 已排除（各专项报告）：源项体积与 master/slave 表示（三级消融 A=B=C 逐位相同）；
   10 种 generalized-Einstein/Fermi 边平均（最优 0.98x，该界面簇 Fermi 修正幅度本身
   仅 0.3–3.6%）；氧化层 couple 污染（当前 profile 已是 Si-only external AverageBox
   couple；node 4492 精确系数重放已消化该类）；Poisson AverageBox edge couple；
   Measure 归一化；ErrRef；1/Vt；线性预条件；状态精度；Vela 解析 Jacobian；
   2D 厚度/单位；`Math{ElementEdgeCurrent}`（限本次 VSV A/B 观测：未改变可观测
   RHS、终态与端口电流；不作普遍断言）。
7. G3 物理段：Fermi 统计、OldSlotboom BGN、SRH/Auger、Silicon 中
   `Mobility(HighFieldSaturation)`（HFS **开启**，内部 GradQF 驱动、接触单元
   ElectricField fallback）。只读 A/B 表明：**关闭 HFS 使绝对电流改变约一个数量级**，
   但不改变冻结态 log 斜率。
8. 约束：T-2022.03-SP2 无 assembly 级逐边电流真值可导出（`eCurrentDensity`
   仅顶点后处理场）；收敛 equation-balance 只约束行和。

---

## 二、证据链统一判读

### 2.1 行缩放不可能是根因（数学定性）

若 `F_vela(x) = D(x) · F_sent(x)`（D 为正对角/行缩放），则两算子零点集**必然相同**。
自洽固定点不同这一事实本身即证明差异不是行缩放。`~33.5x` 的归属：AreaFactor
A/B 已排除面积/宽度归因，结合手册与其他 A/B 支持内部归一化解释（非直接证明）；
但无论其确切归属如何，正行缩放不改零点集，故该标量降级为非根因，不再投入归属工作。

### 2.2 根因对象：`x_sent` 不是 Vela 的零点（所需修正视角）

不得声称“差异算子 `delta_F = F_vela − F_sent` 已被直接测得”：Sentaurus
NewtonPlot 行缩放未知，两引擎输出不在共同归一化空间，不能直接相减成函数。
归一化无关的成立事实只有一条——零在任何正行缩放下都是零，而：

```
x_sent 不是 Vela 的零点；
r_required = F_vela(x_sent) = 已测得的 SSS 七节点残差（Vela 单位）
B · delta_Phi_required = −r_required
```

`delta_Phi_required` 是“使 Sentaurus 状态成为 Vela 零点所需的一组可行边修正”，
**不是**实际测得的 Sentaurus/Vela 逐边通量差。因此 T2 的定位是可行性分析与
候选筛选，不能单独宣称完成根因分解；根因确认必须由 T3/T4/T5 的正交控制
实验承担。由事实 3，`r_required`（1.24e-9 / 6.98e-9 A/um）是纯输运项、
界面局域、偏压相干增长的结构化小量。

### 2.3 关键机制假设（可证伪，统一解释全部观测）

**通量差集中在 Bernoulli 饱和的界面法向边**（反型层节点与体侧节点密度比达
多个数量级）。注意 SG 饱和的方向性：通常只压制**下风端**的 QF 灵敏度，
上风端灵敏度仍可保持 `O(flux/Vt)`；“单模列探不到”只对特定的扰动节点 ×
边方向组合成立，须由 T1 逐边解析灵敏度表逐边判定，而非默认成立。
该假设同时解释四组表面矛盾的观测：

| 观测 | 饱和边机制下的解释 |
| --- | --- |
| 单模 Jacobian 列去标量后匹配 0.058% | 若 node 3721 处于相关饱和边的**下风端**，则该边通量对其 phin 的灵敏度被指数压制，±1 µV 扰动探不到该差异（上风端仍为 `O(flux/Vt)`，上/下风身份须 T1 逐边核实）；0.058% 匹配**不自动构成排除** |
| SSS 端口电流/gm 闭合到 0.02% | 界面法向饱和边不携带净沟道电流，drain-cut 积分不受影响 |
| 七节点行残差非零且结构化 | 饱和边通量全额进入界面节点行平衡 |
| 自洽 phin 偏移 mV 级、gm 掉 24% | 行失衡驱动 phin 沿低刚度的沟道方向漂移，Boltzmann/Fermi 因子放大为密度/电流差 |

### 2.4 候选空间盘点

结合第一节事实 6 的排除清单，**尚未测试**的逐边通量因子族只剩：

- **H1：逐边 HFS GradQF 驱动力离散约定。**
  G3 中 HFS 是 ~10x 量级因子（事实 7）。Vela 有
  `high_field_gradient_discretization = transport_cell_vector / edge_projection`
  两种现成实现；Sentaurus 的 element 级驱动力约定未知。栅边缘/沟道界面处二维
  梯度方向剧变，是三种离散分歧最大的位置；逐边 µ 差几个百分点即可产生观测到的
  0.006%–0.26% 行失衡，而 drain-cut 加权聚合（已验证匹配 0.016%）对逐边分布差
  不敏感。
- **H2：Fermi 简并密度支撑 + OldSlotboom 带边差在 SG 参量中的联合放置。**
  已测的 10 个变体只动了 generalized-Einstein 因子 g；带边差（VariableNi 类
  `log(ni1/ni0)` drift 项在 Fermi 下的对应物）与简并因子在 Bernoulli 参量中的
  放置从未变过。沟道横向掺杂梯度与强反型简并恰好都落在热点簇。
- **H3：未枚举的联合约定**（由 H1/H2 敲除结果收敛后再定义）。

---

## 三、对报告第七节六个问题的直接回答

1. **多节点独立 QF 基扰动如何重建局部 Jacobian？**
   7 节点 phin 基 ±ε 各一次 iteration-0 NewtonPlot（约 14–28 次 VM 运行）可重建
   7×7 块；无单位比较量 = 行归一化矩阵 + 列内比值 + 跨列标量恒定性检验。
   **必须补 psi 列**（现只探过一个 phin 列）。但按 2.3 的饱和边机制，该实验对
   根因判别力低于物理敲除（第四节 D2），建议降为条件任务（D5）。
2. **行缩放等价为何固定点仍不同？**
   正行缩放保零点集（2.1），等价关系必然在某处破缺。三类候选：
   (a) 未探测方向（psi 列、其他节点行）；(b) 状态相关缩放；
   (c) **饱和边通量差，且扰动节点恰处于下风端使 phin 灵敏度被指数压制**
   （候选机制，见 2.3；上/下风身份须 T1 逐边核实）。
   不需要"未覆盖常数项"这种额外假设——(c) 在探测过的方向上自然表现为近似常数项。
3. **Sentaurus 界面顶点 box 语义？**
   按 element-wise 装配（手册式 1266），载流子方程只累计半导体单元贡献——
   Vela 已通过 Si-only external couple 对齐这一层。剩余未知是逐边通量核
   （驱动力放置/密度支撑/简并修正），公开文档不可得，只能靠敲除实验反推
   （D2）或 Synopsys 渠道（D5）。
4. **ElementEdgeCurrent 不影响 NewtonPlot 意味着什么？**
   仅能作观测范围内的表述：在本次 VSV A/B 中，激活它未改变可观测 RHS、终态
   与端口电流。现有证据**不足以**普遍断言它只用于后处理、完全不参与任何 SG
   装配路径。结论不变——对本问题判别力为零、不再投入——但所有引用必须保留
   该观测边界，不得升级为对 Sentaurus 内部实现的普遍声明。
5. **能否构造零点不变量绕开逐边电流不可导出？**
   能，且是本计划核心：Sentaurus 收敛行和恒为零 ⇒ 对每个节点 i、每个同合同
   状态 `sum_e(delta_phi_e) = −r_i^vela`。以 Vela 逐边分解为基，在节点—边关联
   矩阵约束下做加权最小范数逆 + 假设族单参数拟合（第四节 D1/T2，含正则化
   稳健性防护）。注意先例边界：PN2D Direction-2 脚本
   （`scripts/diagnose_pn2d_bv_flux_reformulation_feasibility.py`）只实现了
   逐节点 `s_i = −residual_i/flux_i` 一致性检验，是**方法论**先例；关联矩阵
   逆、共享边反对称、多正则化与置换检验均为 T2 新开发内容。
6. **24.31% 能否直接停止分类？**
   现在**不能**：通量核假设空间未穷尽（H1/H2 未测）、无敲除因子分解、
   逆问题一致性未做。停止判据见第六节。

---

## 四、调查方向（优先级排序）

### D1（P0，纯 Vela 本地，最便宜）：残差标度律 + 逐边归因逆问题

在与 G3 合同**同物理**的状态序列上批量计算 SSS 七节点/界面带残差。状态集合
必须满足 `Vd=0.1 V`、G3 物理段完全一致：本地现有 `0–1 V` 七个同合同状态
加 `1.1667 V` 端点共八个，基本满足入口数量；**不得混入 21 个代表状态中的
IdVd/BV/全物理态**（物理合同不同，残差跨合同不可比）。如需加密序列，在 VM
上以同一 G3 物理补采（一次运行多点导出，成本约等于一次敲除运行）：

- **热点边解析灵敏度表（前置产出）**：每条热点入射边对两端点 `phin` 的
  解析灵敏度、Bernoulli 参量 `eta` 与上/下风分类——先判定单模扰动对哪些
  通量族确实“探不到”，再决定 2.3 机制假设与单模 Jacobian 证据的相互约束；

- 回归 `log ||r||` 对 `log n_interface`、`log Id`、`Vg`、局部 `|grad phin|`，
  得到判别性标度指数（∝密度 / ∝电流 / ∝驱动力 三类假设的指纹）；
- **共线性防护**：同一 Vg 路径上 `n_interface`、`Id`、`|grad phin|` 高度共线。
  必须报告设计矩阵条件数与 VIF，并做留一交叉验证；共线性超限（如 VIF>10）
  时，T1 只作描述性证据，不得单独锁定 H 族——因果判别由 T3/T4 正交控制承担；
- 对每状态解逐边修正 `{delta_phi_e}`，**必须按守恒结构建模而非自由变量**：
  以节点—边关联矩阵 `B` 表达（每条内部边以 ±1 同时进入两端点行，强制共享边
  反对称性），求解 `B · delta_phi = −r` 的加权最小范数问题；
- **图定义（冻结，保证符号指标可复现）**：节点域 = 评估集内的自由硅
  电子 continuity 行（排除 Dirichlet/接触行与氧化层节点）；边域 = 与节点域
  相接的硅输运边（与 Si-only couple 合同一致）；规范定向 = 每条边从小节点
  ID 指向大节点 ID，`B(tail)=+1`、`B(head)=−1`；**跨域边**（仅一个端点在
  节点域内）只进入域内端点行，域外端点不设方程；所有符号翻转比例、反对称
  检验均在该规范定向下计算与报告；
- **假结构防护**：该逆问题高度欠定，解依赖正则化选择。归因结构必须在至少
  三种正则化下稳定（unweighted L2、按 `|phi_e|` 加权、按 couple/边长加权），
  并通过边属性标签置换检验给出显著性，否则不得作为定向证据；
- 检验修正与边属性的相关结构：`|eta|` 饱和度、界面法向/平行取向、钝角单元
  归属、掺杂梯度、HFS 驱动力大小。

**双重产出**：若修正呈单属性结构（且通过稳健性与置换检验）⇒ **提名** H 族
并给出 2.3 假设的直接检验（筛选证据，根因确认仍需 T3/T4/T5 正交控制）；
若逐节点不一致（变号、跨数量级、且在全部测试正则化下均无属性结构）⇒
构成停止分类证据之一（PN2D Direction-2 标准）。

### D2（P0，需 VM，判别力最高）：物理敲除矩阵

三个隔离 Sentaurus VM 变体（仅物理段单项开关，网格/参数/其余物理不动，
deck 与 TDR SHA-256 留痕）。**状态再生协议（防混淆的核心）**：Fermi/BGN 类
敲除会同时改变接触中性条件、载流子映射和自洽状态本身，因此每个变体必须
**各自完整生成新的 Sentaurus 自洽收敛态**（2 个端点偏压），Vela 再以与该变体
完全相同的物理**对该变体自身的状态**做 SSS 重放；严禁把基线状态与变体物理
交叉配对，否则会把状态变化误判为离散变化：

| 敲除 | 变体 | 判别对象 |
| --- | --- | --- |
| K1 | HFS 关（常低场迁移率） | H1 全族 |
| K2 | Fermi → Boltzmann（HFS 保持） | H2 简并支撑分支 |
| K3 | OldSlotboom BGN 关 | H2 带边差放置分支 |

**跨引擎物理资格门（每变体先决条件，预注册阈值）**：K2/K3 会改变
`ni/Nc/Nv` 与载流子重构、接触中性/QF 边界、OldSlotboom 能带分配，并切换
Vela 实际进入的 SG 公式分支；“相同物理”声明不足以保证可比。每个变体
进入判定前必须通过：

- (a) 节点 `n/p` 重构（Vela 由该变体 `psi/phin/phip` 重构 vs Sentaurus
  导出场）：自由硅节点上 `|log10(n_vela/n_sent)|` 与 `|log10(p_vela/p_sent)|`
  的 P95 `≤ 1e-3 dex`、max `≤ 1e-2 dex`；
- (b) 接触行：接触节点 `max|phin − V_c| ≤ 1e-9 V`、`max|phip − V_c| ≤ 1e-9 V`，
  接触中性残差 `|n − p − N_net| / max(n, p, |N_net|) ≤ 1e-6`；
- (c) 固定状态端口电流：两端点均 `|I_vela/I_sent − 1| ≤ 1e-3`
  （G3 基线 SSS 为 `2.02e-4`，留约 5x 裕度）；
- (d) 网格/边集合：顶点数、坐标哈希与输运边集合精确相等（零容差）。

任一资格门失败 ⇒ 该 K 分支记为**“合同不等价”**，不得使用 `>10x/<2x` 做
离散归因或洗脱。

预注册判据：判定量为**变体内部**引擎间归一化残差 `||r_7|| / sum_e|phi_e|`
（在该变体自身收敛态上评估）；跨变体只比较此归一化量、不比较绝对残差。
**分母定义（冻结）**：`phi_e` = Vela SG 电子**粒子**线通量经既有 particle
scale 物理化为 `A/um` 的逐边通量（与 `r` 同单位、同一次装配导出）；边集合 =
评估节点集的全部入射硅输运边（固定七节点即 40 边集合；变体 top-7 用其
自身入射集）；零分母保护：`D = max(sum_e|phi_e|, 1e-6 × |I_port^variant|)`，
若 `sum_e|phi_e|` 低于该下限则指标记无效、分支自动落入不确定带。
相对 G3 基线下降 `>10x` ⇒ 该物理项的**离散**是载体；下降 `<2x` ⇒ 整族洗脱；
下降落在 `2–10x` ⇒ 判定为**不确定**：既不得归因、也不得洗脱，该族保持
开放，交由 T2 归因结构与 T5 定向变体裁决（必要时增补同合同状态再评）。
为防热点随状态转移，固定七节点与变体自身重算的 top-7 热点两套集合并列报告。
**机器可读产出（预注册 schema）**：每变体每偏压输出 `analysis/summary.json`，
字段至少含 `{variant, bias_V, gate_np_p95_dex, gate_np_max_dex,
gate_contact_qf_max_V, gate_neutrality_max_rel, gate_port_rel, gate_mesh_equal,
r7_l2, sum_abs_phi, denominator_floor_hit, normalized_residual,
ratio_vs_baseline, hotspot_set, verdict}`；判定与复核只读该 summary，
不读中间产物。
方法论直接移植自 PN2D no-impact 控制实验。三次运行即可对联合约定空间做
因子分解，远快于在 Vela 侧枚举变体。

### D3（P1，可立即并行）：H1 快速 A/B

Vela 现成开关 `high_field_gradient_discretization`（transport_cell_vector /
edge_projection）在 2 个端点 SSS 上只读重放。验收与 T4/T5 统一为“0.5x 主门 +
迁移保护门”：固定七节点通过 0.5x（L2、max、逐点，两端点）的同时，候选
自身 top-7、完整界面带 L2/max、全硅电子残差 L2/max 均不得恶化，热点集合
转移必须显式报告——防止候选把残差推到相邻节点换取表面改善（node 4492
审计的误差转移即为先例）。成本最低、先验最强，若直接命中可短路 D2。

### D4（P1）：晋级链

命中变体 → 定 psi carrier reclose 复现两端点 gm（廉价预测器，符合冻结边界
"固定状态实验通过预注册改善门"）→ 门通过才允许 31 点曲线重跑。改善门为
"gm + 保护门"联合门：段 gm 误差 ≤20%，且两端点 log 电流误差、最坏 KCL
误差、reclose 收敛性（迭代数/收敛标志）均不得劣于修改前基线——防止变体
以恶化电流幅值或守恒性为代价换取 gm。

### D5（P2，条件执行）：补全 Jacobian 重建 + Synopsys 渠道

仅当 D1+D2 均不定向时执行：7 节点 phin 基 + psi 列 NewtonPlot 重建（问题 1）；
向 Synopsys 提交界面顶点载流子装配核与 NewtonPlot RHS 归一化的版本专属询问
（问题 3/4）。不阻塞主线。

---

## 五、任务计划（WP3）

| # | 任务 | 输入 | 方法 | 验收门（预注册） | 成本 | 依赖 |
| --- | --- | --- | --- | --- | --- | --- |
| T1 | 残差标度律 + 热点边灵敏度表 | G3 同合同 IdVg 态（本地 0–1 V 七态 + 1.1667 V 端点；加密则 VM 同物理补采） | `newton_carrier_term_probe` 批量 SSS 重放 + 回归 + 逐边解析灵敏度/上下风分类 | 标度指数置信区间可区分 ∝n / ∝I / ∝grad(phin)，并报告 VIF/条件数与留一稳定性；共线性超限则降级为描述性证据、不得锁定 H 族 | 小；若补采 +1 次 VM | 无 |
| T2 | 逐边所需修正逆问题（可行性/筛选） | T1 逐边表 | 关联矩阵约束逆（共享边反对称）+ ≥3 种正则化稳健性 + 置换检验 + 假设族单参拟合 | 单属性族解释 ≥90% 残差能量（两端点同时）且结构对全部正则化稳定 ⇒ 提名候选族（非根因确认）；逐节点不一致 ⇒ 记入停止证据 | 中（新工具开发：合成图/守恒/正则化稳定性测试为验收前置） | T1 |
| T3 | H1 快速 A/B | 2 端点 SSS | 切换现有 HFS 梯度离散实现，只读重放 | 0.5x 主门（固定七节点 L2/max/逐点，两端点）+ 迁移保护门：候选自身 top-7、完整界面带 L2/max、全硅电子残差 L2/max 均不得恶化，热点转移显式报告 | 小 | 无（并行） |
| T4 | 敲除矩阵 K1–K3 | VM + 2 端点 | 每变体独立生成 Sentaurus 收敛态，Vela 以相同物理重放**该变体自身状态**；先过跨引擎资格门（n/p 重构、接触 QF/中性、端口电流、网格/边集合），再判定变体内归一化残差 | 资格门失败 ⇒ 记“合同不等价”、不得归因/洗脱；通过后：>10x 归因 / <2x 洗脱 / 2–10x 不确定——变体内，双热点集合 | 3 次 VM 运行 + 导入 | 无 |
| T5 | 族内定向变体 | T2/T3/T4 提名族 | 只读通量核变体重放 | 与 T3 相同：0.5x 主门 + 迁移保护门（候选 top-7、界面带、全硅残差不恶化，转移显式报告） | 中 | T2–T4 |
| T6 | reclose gm 门 | T5 命中变体 | 定 psi carrier reclose 两端点 | 段 gm 误差 ≤20%；保护门：两端点 log 电流误差、最坏 KCL、收敛性（迭代数/标志）均不得劣于基线 | 小 | T5 |
| T7a | 31 点复审 | T6 通过的候选 | 全曲线重跑 | 31 点全部门通过（含 max-gm ≤20%） | 中 | T6 |
| T7b | ledger 证据包 | T2/T4 产出 + 第六节判据 | 停止分类材料汇编与批准申请 | 第六节三条判据同时满足 | 中 | T2、T4（**不依赖 T6**：无候选通过 T5/T6 时亦可触发） |

**建议启动顺序**：T1 + T3 立即并行（纯本地）；T4 三个 VM deck 同步准备提交；
T2 在 T1 逐边表落地后一天内出结论。最快路径 = T3 直接命中；最稳路径 = T4 因子分解。

---

## 六、停止分类判据（ledger draft → approved）

同时满足以下三条方可申请批准为 known difference：

1. **敲除洗脱**：T4 三项敲除均 <2x（HFS/Fermi/BGN 的离散均非载体；
   任一项落入 `2–10x` 不确定带即不满足本条，须先行裁决）；
2. **逆问题不可逆**：T2 所需逐边修正呈逐节点不一致结构
   （参照 PN2D Direction-2 定量标准：变号比例 >15%、幅值跨 ≥2 个数量级、
   且在全部测试正则化下均无单一边属性可解释 ≥50% 残差能量）；
3. **影响定量封装**：仅 max-gm 单门失败 24.31%；中位 0.070 dex、P95 0.100 dex、
   31/31 收敛、最坏 KCL 1.1e-9；分类为"受支持诊断范围内未解析的界面装配
   差异"（H3 仍是开放集合，不得表述为已证明的引擎内在差异），
   IALMob 解锁决定移交批准人。

---

## 七、冻结边界符合性声明

本计划全程不触碰以下冻结项：

- 不调整迁移率、阈值、flatband、bulk mobility 参数吸收差异；
- 不继承 PN2D `element_edge_sg_gss_laux` 原子 profile；
- 不把 NewtonPlot 绝对 RHS 用于跨引擎验收（仅用其无单位形状/比值）；
- 不启用 IALMob 或 predictor；
- 31 点曲线仅在 T6 预注册改善门通过后重跑；
- ledger 在第六节判据满足或修复落地前保持 draft。

所有新实验均为只读重放或隔离 VM 诊断运行，不改生产 C++ 默认。

---

## 八、关键工件索引

- 31 点曲线：`docs/validation/templates_ldmos_g3_averagebox_curve_qualification_2026-08-31.md`
- gm 冻结状态：`docs/validation/templates_ldmos_g3_gm_segment_fixed_state_audit_2026-08-31.md`
- Sentaurus 方程余额：`docs/validation/templates_ldmos_g3_sentaurus_equation_balance_2026-09-01.md`
- 单模 QF 扰动：`docs/validation/templates_ldmos_g3_qf_single_mode_perturbation_2026-09-01.md`
- 七节点三级消融：`docs/validation/templates_ldmos_g3_seven_node_interface_ablation_2026-08-31.md`
- Fermi 边平均 A/B：`docs/validation/templates_ldmos_g3_fermi_edge_average_ab_2026-09-01.md`
- 界面主从节点审计：`docs/validation/templates_ldmos_g3_interface_pair_box_audit_2026-08-31.md`
- node 4492 审计：`docs/validation/templates_ldmos_g3_averagebox_node4492_audit_2026-08-31.md`
- RHS 量纲链：`docs/validation/templates_ldmos_g3_rhs_dimension_chain_2026-09-01.md`
- 诊断合同：`reference_tcad/templates_ldmos_sentaurus2022/contracts/diagnostics/templates_ldmos_external_averagebox_profile.json`
- 差异账本：`reference_tcad/templates_ldmos_sentaurus2022/known_difference_ledger.json`
- PN2D 方法论先例（一致性检验）：`scripts/diagnose_pn2d_bv_flux_reformulation_feasibility.py`（主仓库）
