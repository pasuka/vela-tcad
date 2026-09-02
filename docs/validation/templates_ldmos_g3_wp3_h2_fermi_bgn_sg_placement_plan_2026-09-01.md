# Templates/LDMOS G3 WP3 H2：Fermi/BGN 联合 SG 参量放置固定状态控制计划

日期：2026-09-01
状态：草案（未执行，待审查）
依据提交：`dbcece1`（T1–T5 诊断工件冻结）
上位规范：`templates_ldmos_g3_gm_wp3_investigation_plan_2026-09-01.md` R4
范围：`Vd=0.1 V`，`Vg=1.0 V` 与 `1.1666667 V` 两个最大 gm 端点的固定状态电子连续性算子控制

本文只制定 H2 的预注册实验合同，不执行实验，不修改生产默认值，也不授权
T6 或 31 点曲线。所有量化门槛以上位 R4 计划为唯一规范来源；本文复述值仅为
执行便利，发生冲突时以上位计划为准。

---

## 一、决策背景

T1–T5 已把当前 WP3 状态收敛到以下边界：

1. T1 的八态标度回归受共线性限制，只能作为描述性证据；热点边的 HFS 驱动力
   与 `|eta|` 饱和度有较高但不足以锁定根因的解释率。
2. T2 逐边逆问题没有形成跨正则化稳定、满足守恒并可通过置换检验的唯一候选；
   它不能单独确认根因，也尚不足以触发停止分类。
3. T3/T5 的 HFS 驱动力与幅值变体未通过两端点 `0.5x` 主门和热点迁移保护门；
   H1 现有枚举未命中，不进入 T6。
4. T4 的 K2/K3 变体未通过跨引擎载流子重构、接触或固定状态端口电流资格门，
   因而是“合同不等价”，不能用于洗脱 H2。
5. 既有 Fermi edge-average A/B 只改变 generalized-Einstein 因子 `g` 的平均方式；
   10 个变体均未命中，但它没有改变 Fermi 修正与 OldSlotboom 带边项在端点
   简并坐标和 Bernoulli 漂移参量之间的联合放置。

因此，H2 仍是开放假设。不过，下一轮不得再次做无界参数扫描；应使用有限、
可证伪的二因子表示矩阵，并在解释残差前先证明候选的代数不变量与合同资格。

## 二、目标、假设与明确排除项

### 2.1 目标

在完全相同的两个 Sentaurus 完整收敛状态上，回答一个有限问题：

> Vela 当前由 Fermi + OldSlotboom 有效 `ni/Nc` 显式构造的端点统计坐标和
> SG 带边漂移参量，是否与同一冻结 `n/phin` 状态所要求的热力学闭合表示存在
> 结构化差异；该差异是否足以关闭七节点电子连续性残差，并且不会把热点转移
> 到界面带或全硅区域？

### 2.2 本轮不测试的内容

- 不重新枚举 `g` 的算术、几何、调和或端点平均；既有 10 变体结论保持冻结。
- 不关闭全部 OldSlotboom BGN，不关闭 Fermi 统计，不重跑 K2/K3 VM 状态。
- 不改变 `n/p`、`psi/phin/phip` 冻结输入，不调整接触中性条件。
- 不改变 HFS、迁移率、AverageBox couple、网格、复合、源项或端口积分。
- 不引入连续权重、经验比例、阈值平移或迁移率标定。
- 不把任一固定状态改善直接解释为 Sentaurus 内部实现真值。

## 三、冻结输入合同

H2 必须复用 T3/T5 已审计的两套完整 Sentaurus 状态及同一 Vela 基线：

| 项目 | 冻结值 |
| --- | --- |
| 偏压 | `Vd=0.1 V`；`Vg=1.0 V, 1.1666667 V` |
| 状态 | 各端点完整 Sentaurus `psi/phin/phip/n/p`，不得使用 VSV 混合态 |
| 统计/BGN | Fermi–Dirac + OldSlotboom；合同 revision 4 |
| BGN 参数 | `Nref=1e17 cm^-3`、`coefficient=0.009 eV`、`smoothing=0.5`、`offset=0` |
| Fermi 修正 | 生产基线开启；按 physics contract revision 4，`dEg0=-0.01595 eV` 由生产 BGN field 处理一次，`materials[].intrinsic_carrier_density_cm3` 不含该项 |
| 输运几何 | `legacy_node_local + external AverageBox`，Si-only 输运边 |
| HFS | G3 生产基线 `edge_projection`，保留既有接触 ElectricField fallback |
| 禁用项 | IALMob、predictor |
| 主集合 | 固定七节点 `3721,3974,3973,4091,3727,4021,3771` 及其 40 条入射硅边 |
| 保护集合 | 候选自身 top-7、完整栅沟道界面带、全部自由硅电子 continuity 行 |

执行前生成 `protocol.json`，记录状态文件、网格、材料合同、配置、基线边表和
代码提交的 SHA-256。两个端点的节点数、坐标哈希、Si-only 边集合及 40 边集合
必须精确相等；任何不一致均停止实验。

## 四、数学定义与候选矩阵

### 4.1 两套端点统计坐标与两套 Bernoulli 参量

本轮不直接打开或关闭任一 BGN/Fermi 分量。先比较两种对**同一物理状态**的表示：

- `X = explicit_band`：生产公式由 Fermi + OldSlotboom 有效 `ni/Nc` 显式构造；
- `N = density_inverted`：由冻结载流子密度反演 Fermi 积分得到，不猜测 BGN
  在带边中的分配。

对电子端点 `i` 定义：

```text
z_i                = ln(ni_eff_i / Nc_i)
eta_i^X            = (psi_i - phin_i) / Vt + z_i
eta_i^N            = inverse_Fermi_half(n_i / Nc_i)
g^S                = secant_generalized_Einstein(n_0, n_1, eta_0^S, eta_1^S), S in {X, N}
D_band             = (psi_1 - psi_0) + Vt (z_1 - z_0)
A_band^S           = D_band / (Vt g^S)
A_closure^S        = ln(n_1 / n_0) + (phin_1 - phin_0) / (Vt g^S)
```

`A_closure` 不是拟合式；它来自 SG 左右端密度比恒等式，并保证平坦物理准费米势
时 `A=ln(n_1/n_0)`。生产合同完全一致时，应有 `eta^X=eta^N` 且
`A_band=A_closure`。两者差值正是本计划要测量、再按 OldSlotboom/Fermi 分量
做只读归因的“联合放置闭合差”。

每个候选的电子边通量使用同一对称左右端形式计算：

```text
Phi_SA = C_edge * mu_edge * Vt * g^S
         * [B(-A^S_A) n_0 - B(A^S_A) n_1]
```

其中 `C_edge`、`mu_edge`、`n_0/n_1` 全部来自同一次冻结基线装配。实现必须采用
稳定 Bernoulli/高精度参考，不能依赖把 `phin_0 == phin_1` 直接短路为零来通过
平衡态测试。

### 4.2 预注册 2×2 矩阵

| 变体 | `eta/g` 来源 | Bernoulli 参量来源 | 解释角色 |
| --- | --- | --- | --- |
| H2-P00 | X：显式带边 | band：显式 `psi+ni/Nc` | 生产公式重构控制；必须逐边复现当前算子 |
| H2-P01 | X：显式带边 | closure：密度/QF 恒等闭合 | 只替换 Bernoulli 参量表示；主候选 |
| H2-P10 | N：密度反演 | band：显式 `psi+ni/Nc` | 只替换端点统计坐标/`g`；诊断因子 |
| H2-P11 | N：密度反演 | closure：密度/QF 恒等闭合 | 完全闭合表示；主候选 |

矩阵只允许这四格。不得在结果出来后添加 `lambda` 混合、半分裂、按边选择或
其他连续调参。OldSlotboom 基础项、Fermi correction 和 `dEg0` 的逐分量贡献只
作为 `A_closure-A_band`、`eta^N-eta^X` 的标签/回归解释量输出，不作为第五个
候选，也不得在本轮按经验重新分配。

### 4.3 二因子效应量

除逐候选门槛外，固定输出：

```text
E_eta         = r(P10) - r(P00)
E_argument    = r(P01) - r(P00)
E_interaction = r(P11) - r(P10) - r(P01) + r(P00)
```

效应量只用于说明“统计坐标、Bernoulli 参量表示、交互项”哪一项改变残差，
不取代资格门，
也不代表已经测得 Sentaurus/Vela 的逐边通量差。

## 五、资格门：先判合同，再看残差

### 5.1 输入与生产重构门

1. `H2-P00` 必须在两个端点逐边复现当前生产电子 SG 通量：相对 L2
   `<=1e-10`，最大绝对差 `<=1e-18 A/um`。
2. 用 P00 边差经冻结关联矩阵重构的电子 continuity 行必须以相同阈值复现
   现有基线残差；端口切面电流相对差 `<=1e-10`。
3. 显式 `eta^X` 与密度反演 `eta^N` 对冻结 Sentaurus `n/p` 的重构门沿用 T4：
   自由硅节点 P95 `<=1e-3 dex`、max `<=1e-2 dex`；同时报告
   `eta^N-eta^X` 与 `A_closure-A_band` 的 P50/P95/max 和空间热点。
4. P00 任一门失败即停止，不得评价 P01/P10/P11。

### 5.2 每个候选的离散不变量

在合成边和器件边上分别验证：

- 交换端点后通量严格反对称：相对误差 `<=1e-12`；
- 平坦准费米势的左右端通量自然相消，不能依赖专用早退分支；
- Boltzmann 极限退化到已资格化的 variable-`ni` SG 公式，误差 `<=1e-10`；
- `z_0 == z_1` 时四格精确退化到同一常 `ni/Nc` Fermi 公式；
- `g`、Bernoulli 参量和通量均有限，`g>0`，不得静默钳位；
- 解析/自动微分 Jacobian 与中心差分一致，默认相对门 `<=1e-6`。

任一不变量失败，候选记为 `algebraically_invalid`，不进入器件评分。

### 5.3 候选合同资格

四格均使用同一冻结 `n/p/psi/phin/phip`，不重新解释接触或载流子状态。进入
残差判定前仍须共同通过：

- 显式与密度反演映射满足 P95 `<=1e-3 dex`、max `<=1e-2 dex`；
- 接触 `phin/phip` 与电压偏差 `<=1e-9 V`；中性残差相对值 `<=1e-6`；
- P00 固定状态两端点端口电流相对 Sentaurus `<=1e-3`；
- 候选端口电流/KCL 不超过第六节保护门。

若共同映射门失败，整个 H2 固定状态实验记为 `contract_not_equivalent`，不能用
矩阵结果做因果归因。若某一格单独违反反对称、自然平衡或退化极限，则该格记为
`algebraically_invalid`；即使七节点残差改善，也禁止提名生产候选、禁止进入 T6，
亦不得用于批准 ledger。

## 六、固定状态残差计算与验收门

### 6.1 只读回放

对每个合格候选，仅替换电子 SG 边通量核：

```text
r_candidate = r_P00 + B * (Phi_candidate - Phi_P00)
```

`B` 使用 R4/T2 已冻结的节点—边规范定向：小节点 ID 到大节点 ID，
`B(tail)=+1`、`B(head)=-1`；接触/Dirichlet 行不进入自由节点评分。空穴通量、
Poisson、SRH/Auger、HFS、mobility、couple 和状态字段均保持 P00。

### 6.2 主门与迁移保护门

候选必须在两个端点同时满足：

1. 固定七节点 `r7_l2`、`r7_max` 与七个逐节点绝对残差均 `<=0.5x P00`；
2. 候选自身 top-7 的 L2/max 不得恶化；
3. 完整栅沟道界面带 L2/max 不得恶化；
4. 全部自由硅电子 continuity 行 L2/max 不得恶化；
5. drain/source 端口电流、全器件 KCL 与基线相比不得超过既有冻结门；
6. 热点集合、热点迁移距离和新增 top-20 节点必须显式报告。

这里的“不恶化”按 R4/T3/T5 的冻结实现读取，不另行放宽。任何只改善固定七节点、
却把残差转移至 node 4492 类邻域的候选均判失败。

## 七、执行顺序

### H2-0：协议与 schema 冻结

- 写入 `protocol.json`、候选枚举、公式版本、输入哈希和所有门槛；
- 扩展机器可读 `analysis/summary.json` schema；
- 先提交方案/schema，再运行数值实验，避免看结果后改门槛。

### H2-1：生产公式重构与字段审计

- 导出两端点每条 Si 输运边的 `ni_eff/Nc/eta_X/eta_N/g_X/g_N/`
  `A_band/A_closure`，并附 OldSlotboom、Fermi correction、`dEg0` 分量标签；
- 完成 P00 逐边、行、端口三层重构；
- 完成显式/密度反演映射资格表与闭合差空间热点表；
- P00 或输入合同失败则停止。

### H2-2：公式级测试

- 合成均匀边、variable-`ni` 边、强简并边、Bernoulli 饱和边；
- 端点交换、平衡态、Boltzmann 极限、常 `ni/Nc`、Jacobian–FD；
- 只允许通过的候选进入 H2-3。

### H2-3：两端点 2×2 固定状态回放

- 一次性运行四格，不增加其他变体；
- 输出固定七节点、候选 top-7、界面带、全硅及端口/KCL 指标；
- 计算二因子效应量与跨端点方向一致性；
- 判定程序只读取 `analysis/summary.json`。

### H2-4：结论与分支

- **P01 或 P11 通过全部资格、0.5x 与迁移保护门**：提名“显式 Fermi/BGN
  带边表示与密度/QF 恒等闭合存在算子级差异”；再依据预注册分量标签判断差异
  是否与 OldSlotboom/Fermi correction 空间梯度一致。此时仍不能声称已证明
  Sentaurus 的具体带边分配，也不得直接启动 T6。
- **仅 P10 改善但自然平衡或映射资格失败**：记录为 `algebraically_invalid` 或
  `contract_not_equivalent`，不提名、不洗脱 H2，仅说明 `g` 对残差有代数敏感性。
- **所有合格候选均失败**：洗脱本计划中有限的 H2 放置矩阵；禁止追加连续调参，
  回到 H3 开放集合与上位计划停止判据。
- **仅一个端点通过或热点迁移**：结论为 `inconclusive`，不得进入 T6。

### H2-5：条件性后续（不属于本计划授权）

只有 H2-4 提名候选后，另行审批：

1. 默认关闭的 LDMOS 专用 C++ 诊断 profile；
2. 电子公式与解析 Jacobian 的单元/FD 测试；
3. 空穴侧对称语义资格，避免只修电子而破坏统一物理合同；
4. 生产探针固定状态复现；
5. T6 两点自洽 reclose，随后才可能执行 T7a 31 点复审。

## 八、机器可读产出

每候选、每偏压至少输出：

```text
variant
bias_V
input_hashes_equal
edge_set_equal
production_edge_reconstruction_rel_l2
production_edge_reconstruction_max_A_per_um
production_row_reconstruction_rel_l2
gate_np_p95_dex
gate_np_max_dex
gate_contact_qf_max_V
gate_neutrality_max_rel
gate_port_rel
gate_antisymmetry_rel
gate_boltzmann_rel
gate_constant_band_rel
gate_jacobian_fd_rel
r7_l2
r7_max
r7_ratio_vs_p00
r7_per_node_ratio
candidate_top7_l2
candidate_top7_max
interface_band_l2
interface_band_max
all_silicon_l2
all_silicon_max
hotspot_set
hotspot_migration
effect_eta_l2
effect_argument_l2
effect_interaction_l2
qualification_verdict
residual_verdict
overall_verdict
```

建议目录：

```text
reference_staging/templates_ldmos_phase_a/
  wp3_h2_fermi_bgn_sg_placement/
    protocol.json
    analysis/summary.json
    analysis/edge_metrics.csv
    analysis/node_metrics.csv
    analysis/factorial_effects.csv
```

原始大文件保持 ignored，不提交仓库；仓库只提交脚本、schema、紧凑 summary、
报告与测试。

## 九、测试、成本与停止纪律

### 9.1 自动测试

- 公式反对称、平衡态自然相消、Boltzmann 极限、常带边退化；
- P00 对现有生产通量的逐边精确重构；
- 2×2 因子效应恒等式；
- 图定向、跨域边和热点集合的确定性；
- schema 缺字段、阈值边界和 `contract_not_equivalent` 禁止晋级测试。

公式/Jacobian/schema 测试作为 CI 硬门；完整器件固定状态输出为计划性人工证据门。

### 9.2 成本预算

本计划主要复用现有两端点本地冻结状态和边表，不需要新的 Sentaurus VM 运行。
预计新增一个离线分析脚本、一组公式/schema 测试和一份报告；执行成本低于 T4，
明显低于任何 31 点自洽曲线。若输入字段不足，只允许补充本地只读 probe；需要
新 VM 状态时必须暂停并重新审批输入合同。

### 9.3 停止纪律

- 不因单个节点改善而放宽 0.5x 或迁移保护门；
- 不因 P10/P11 的未资格或代数无效变体命中而修改生产物理；
- 不把 P01 命中表述为已证明 Sentaurus 内部采用该公式；
- 不在 H2 完成前批准 known-difference ledger；
- 不在 H2-5 独立审批前修改 C++ 默认值、启动 T6 或 31 点曲线。

## 十、审查问题

执行前只需审核以下四点：

1. “显式带边表示 × 密度/QF 恒等闭合表示”是否是无需猜测 Sentaurus BGN
   分配的最小合同等价 H2 控制？
2. 2×2 矩阵是否覆盖 H2 的最小联合表示空间，且没有与既有 10 个 `g` 平均变体重复？
3. 自然平衡、映射和迁移保护门是否足以防止把不自洽的反事实误升格为物理结论？
4. P01/P11 若通过，是否仍应先形成有参数来源的分量归因与默认关闭诊断 profile，
   并完成空穴侧对称资格，
   再进入 T6？
