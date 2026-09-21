# Vela AC 小信号规范合同 v4

修订日期：2026-09-10。代码复核基线：`d76434089dcd8564069be4b822260d7aa9c72662`。
本文规定候选开发合同，不声明数值门已经执行。

阅读配套 [阶段实施计划](../plans/2026-08-31-ac-small-signal-simulation-development-plan.md) 和
[审核与源码证据附录](../plans/2026-09-10-ac-small-signal-review-appendix.md)。
为延续审核引用，本文件保留原 §4.1--4.9 和 §11.1 编号；新增材料支持细则为 §4.3.4。
P0.5 与 WP-R 是同一阶段；各 AC-L 映射在实施计划 §3.3。
材料支持、poly 栅和有效参考导数等 pending 决策按阶段门处理，不能用本文件替代审批。

## 4. 必须先冻结的数学和数据契约

### 4.1 端口方向和矩阵索引

建议冻结如下约定：

- `Y[row_terminal, column_terminal] = dI_row / dV_column`；
- 端电流流入器件为正；
- phasor 采用 `Re(x_hat*exp(+j*omega*t))`，因此电容电流为 `+j*omega*C*V_hat`；
- 每列激励幅值默认为 `1 V` 复幅值，其他端口 AC 电压为零；
- 输出始终带显式端口顺序；
- 用户可指定参考端用于 reduced matrix，但生产结果保留 full matrix；
- `A = Re(Y)`；
- `C = Im(Y)/omega`，仅允许 `frequency_Hz > 0`；
- `omega=0` 使用单独的准静态接口，不在 AC 中除零；
- 不假定 `Yij=Yji` 或 `Cij=Cji`，除非 fixture 本身是平衡无源互易系统。

### 4.2 二维单位

二维器件的规范输出使用每米器件宽度：

- 电流：`A_per_m`；
- 电荷：`C_per_m`；
- 导纳/电导：`S_per_m`；
- 电容：`F_per_m`。

可同时输出 `*_per_um` 便利列，换算必须为规范每米值乘 `1e-6`。如果用户提供有限的
器件宽度，则另输出总量并在 manifest 中记录归一化方式。配置遵循当前单位版本：
legacy 使用 `depth_m`，`format_version: 2` 使用 `depth_um`，解析后统一为米。
禁止在同一列中混用总量和每宽度量；P6 电路端必须用乘过实际宽度的 A/S/F 总量。

### 4.3 电极反力、接触载流子储存和自然通量端口

以下是依据当前残差符号的离散推导，仍须由 P0.5/P2/P4 的独立测试验证。先将行恢复
为物理每米单位：Poisson 为 `C/m`，连续性为粒子数流率每米，储存为粒子数每米。
记这些行是 `r_psi, r_n, r_p, s_n, s_p`，不是已经覆盖的 Dirichlet 行。

#### 4.3.1 Essential 接触与符号

当前 Poisson 边贡献为 `G*(psi_i-psi_j)`，源项为负的节点净电荷，因此：

```text
r_psi,i = F_D,i - rho_i
rho_i = q*(p_i-n_i+Nd_i-Na_i)*V_i + fixed_charge_i
Q_k = sum_(i in essential Poisson contact k) r_psi,i
Qx_k = w_psi,k^T * J_psi,phys       # 行向量；列梯度是 J_psi,phys^T*w
```

`F_D` 是装配器的内部边电位移通量和，`V_i` 为该载流子/掺杂项已批准的物理控制面积，
不默认等于跨材料几何总面积；多材料推广和门槛见 4.3.4。对包围半导体的外法向，
金属电极电荷等于负的器件外向 D 通量。故上述反力符号为正，不再留一个任意待选的
`s_Q=+/-1`。缩放模式下 Poisson 行乘 `epsilon_ref*V0` 恢复 `C/m`，legacy SI 因子为 1。
固定体电荷和 sheet charge 均须来自同一 Poisson 装配源项并且只计一次。

连续性行共享同号净复合源 `U_i`，从当前 SG flux/残差和 DC 电流权重得：

```text
r_n,i = F_n,i + U_i;       r_p,i = F_p,i + U_i
s_n,i = V_i*n_i;           s_p,i = V_i*p_i
i_carrier,k = sum_(i in contact k)[-q*(r_n,i + ds_n,i/dt) + q*(r_p,i + ds_p,i/dt)]
i_total,k = i_carrier,k + dQ_k/dt
```

DC 时退化为现有 `computeFromResidual` 的 `-q*r_n+q*r_p`。代码命名的
`holeCurrent` 并非可直接与 `electronCurrent` 相加的物理分量；必须以其
`totalCurrent=electronCurrent-holeCurrent` 和上述行权重为准。

若使用数值归一化行，电流恢复因子是
`continuityResidualScale()*currentDensityLineIntegralFactor`，然后施加电子 `-q`、
空穴 `+q`；不得对已经恢复的物理行再次乘因子。储存内部原始单位使用
`vol*n*continuitySourceIntegralFactor()`（p 同理），再除连续性行尺度。
该 source factor 与 current line factor 的乘积必须等于 `chargeAreaFactor`。
时间统一为秒，不凭空引入 `t0` 或漏乘 `2*pi`。

对静态掺杂/固定电荷和成对 G/R 源，接触局部恒等式为：

```text
i_carrier,k = sum_(i in contact k)(-q*F_n,i + q*F_p,i) + d rho_contact,k/dt
Q_k = F_D,contact,k - rho_contact,k
i_total,k = sum_(i in contact k)(-q*F_n,i + q*F_p,i) + d F_D,contact,k/dt
```

因此“physical M 接触反力 + Poisson 反力时间导数”没有重复计数：接触控制体电荷
项恰好抵消。若只保留 DC 反力导数再加 `j*omega*Qx`，却丢掉 physical M 接触行，
反而会遗漏这项抵消。测试必须使用非零接触控制体积，不能靠令接触体积为零掩盖错误。
纯绝缘电极没有载流子储存/流入；不把绝缘区 phin/phip pin 当作物理电流。

对于 essential 载流子接触，定义物理权重 `w_c=(-q,+q)` 后，2.2 的输出算子为：

```text
H = w_c^T*J_cont,phys;        D = w_c^T*Ju_cont,phys
K = w_c^T*M_cont,phys + Qx;   E = w_c^T*Mu_cont,phys + Qu
```

接触行集合按方程分别确定；不能把只有 psi Dirichlet 的电极自动归为三方程 Dirichlet。

#### 4.3.2 Robin/thermionic 不能用完整残差反力代替边界电流

自然边界满足 `r_bulk+r_boundary+ds/dt=0`，其完整 solved/physical 残差响应趋于零，
但外部端电流一般不为零。该端口应从同一装配器返回的 boundary stamp 计算：

```text
i_carrier,k = -w_c^T * r_boundary,k
i_carrier,k_hat = -w_c^T*(J_boundary,k*x_hat + B_boundary,k*u_hat)
```

Poisson 若仍是 essential，则继续用其电荷反力。保留所有自然项的 physical view
仍然必要，但必须另提供按接触标识的 `boundary residual/J/du`，不能在后处理重写
thermionic 导数，也不能传空 bcs 后对全部自然接触盲目求和。自然通量的显式电压依赖
会使 RHS 出现在连续性行；“仅 BC 行非零”的结论仅限 4.5 的理想 Dirichlet 子集。

P2 必须验证自然行的 M 策略；生产 MVP 先支持 ideal ohmic/metal gate。Robin AC
只有在电压导数、表面通量符号、DC/AC KCL 的专门资格门全部通过后才开放，否则
fail-closed。内部异质界面不是外部端口，内部 sheet/interface flux 不能再算成端电流。

#### 4.3.3 Gauss 闭合与几何范围

对电边界完整、其他外边界零通量的 fixture：

```text
sum_k Q_k + Q_domain = 0
```

若存在非零已知外向自然 D 通量 `Phi_N`，则检查
`sum Q_k + Q_domain - Phi_N = 0`，或把其对应电极纳入完整端口集。内部界面两侧
通量抵消、界面电荷按装配归属计一次。重叠接触节点必须有唯一权属或在解析时拒绝。

Q/I/S/M 复用 WP-R 批准的材料支持测度和 edge coupling；不独立复制几何，
也不把 mixed-Voronoi 变为全局默认。PN2D BV 的 SG/GSS-Laux 原子 profile 和非钝角
资格边界仅适用于原模板。一般钝角几何一致性测试不得冒充该 profile 的资格扩展。

#### 4.3.4 材料支持体积：WP-R 前置决策

几何节点总面积 `V_geom,i=sum_region V_i,region` 本身不是错误；错误风险来自把它
无条件乘以共享 Si/SiO2 节点的半导体载流子/掺杂密度。几何划分（barycentric 或
mixed-Voronoi）与材料支持（legacy_global 或 material_local）是两个独立策略字段。

WP-R 交付 `material_volume_decision.json`，明确选择、范围、几何组合和 DC 重验条件：

- `legacy_global` 保留既有离散用于回归/诊断，已知差异登记为 open；不能仅凭 Q/M
  自洽或 Gauss/KCL 通过就认领 Sentaurus nMOS 等价。
- `material_local` 是待批准候选，按材料/方程支持计算载流子电荷和储存。概念上
  `rho_mobile=q*sum_r V_i,r*(p_i,r-n_i,r)`，
  `s_n=sum_r V_i,r*n_i,r`、`s_p=sum_r V_i,r*p_i,r`。
  当前每节点一套状态不能自动支持任意异质结；先限于已明确映射的 Si/绝缘体及
  经资格化的栅材料组合，不能从此公式推断已有 region-local 状态自由度。
- 移动电荷与 S/M 必须满足 `d rho_mobile=-q*ds_n+q*ds_p`；连续性 G/R 的材料
  支持、掺杂归属、域内净电荷诊断和 Jacobian 同时审计。固定氧化层电荷使用自己的
  材料面积，sheet charge 使用界面长度，不能把所有源项一律换成硅区体积。
- LDMOS 独立工作树已验证的候选只改变 Poisson 移动/掺杂电荷，未改变连续性/储存，
  并限制与 mixed-Voronoi 组合。它是 DC 证据，不是可直接移植的 AC 守恒合同。

在 P1/P2 数值门冻结前完成策略决策；若需生产修改另立 WP-R 子工作包并重新收敛 DC。
相同 snapshot 的残差复现只比较同一已选物理策略，不能拿新体积策略去要求复现旧
策略终态。任何策略变化使旧 Q/M/AC 门失效，重新验收 WP-R、P1/P2、P2.5 和 P5。
不对称 Si/氧化层纵向网格 MOSCAP 的测度/DC 资格在 WP-R、C_Q 在 P1、S/M 在 P2、
AC Cgg 过渡区/Cgs/Cgd 和频率趋势在 P5 分阶段验收。合并报告两策略差异和 Gauss/
网格趋势，不要求 WP-R 先依赖未实现的 P1/P2；不通过拟合阈值、迁移率或宽度吸收差异。

### 4.4 动态储存算子

储存算子归属 `CoupledDDAssembler`，因为它必须共享 node volume、CarrierStatistics、
BGN、frozen QP、DDScalingSpec 和 QF reference。可用独立头文件声明结果类型，但
生产实现不得在独立类中复制这些状态。装配器应返回与稳态残差行单位和缩放一致的
`S(x,u)` 及：

```text
M = dS/dx
Mu = dS/du
```

原则：

- Poisson 内部行没有载流子时间储存项；
- 电子和空穴连续性物理行包含各自控制体储存；
- essential Dirichlet 行在 solved system 中是代数约束，对应 M 行为零，但未约束
  physical row 的 M 必须保留供接触反力提取；
- thermionic/Schottky 是自然通量边界，不替换连续性物理行，储存项必须保留；
- 绝缘材料中为消除零行而生成的 phin/phip pin 是代数 gauge 行，M 为零；
- fixed/interface sheet charge 进入 Poisson 反力，但静态 sheet charge 不产生 M；
- `psi/phin/phip` 使用现有缩放和 quasi-Fermi reference 坐标；
- n、p 对三个未知量的导数必须复用 CarrierStatistics 的一致导数；
- 电子/空穴储存均取正的粒子数，按 4.3 的单位恢复和动态恒等式验收；
- inactive 模型不得改变 M；未来动态陷阱/热/量子方程通过独立 block 扩展。

边界行冻结表：

| 行类型 | solved J 行 | solved M 行 | physical reaction J/M |
| --- | --- | --- | --- |
| Ohmic/metal-gate essential Dirichlet | 替换为代数约束 | 0 | 跳过替换，保留物理行供端口提取 |
| Thermionic/Schottky Robin | 保留自然通量 | 保留连续性储存 | 完整行用于平衡；端电流另取该接触 boundary stamp |
| 绝缘区内部 phin/phip pin | 代数 gauge | 0 | 无物理贡献的 carrier 行为零，pin 不进入反力 |
| Fixed/interface sheet charge | 进入 Poisson 行 | 0 | 进入 Poisson 端口反力和 Gauss closure |
| 未接触内部 DD 行 | 正常物理行 | 正常储存 | 不直接进入端口权重 |

### 4.5 频域系统与三种等价表示

物理系统 `(J+j*omega*M)z=b` 可使用以下三种局部表示；`Ndof=3*Nnode` 是完整 DD
自由度数，不要把 complex Ndof 系统误写成仅有 Nnode 个未知量。

```text
real stacked:
[ J       -omega*M ] [Re(z)] = [Re(b)]
[ omega*M  J       ] [Im(z)]   [Im(b)]

real interleaved: 对上述 2*Ndof 系统施加对应的行/列排列
complex:         (J+j*omega*M) z = b，Ndof 个复数未知量
```

P3/WP3b 必须比较三者，再选择生产后端。局部
`Eigen::SparseLU<Eigen::SparseMatrix<std::complex<double>>>`（参见 [Eigen 官方说明](https://libeigen.gitlab.io/eigen/docs-nightly/classEigen_1_1SparseLU.html)）与独立 wrapper 不要求
修改全仓库 `Types.h`。人工实块构造继续作为符号/交叉验证 oracle，不能将实数表示
既定为 MVP、只把 complex 当作 AC-L6 后的可选优化。

比较固定 DC snapshot、物理算子、矩阵模式、ordering 策略和等化策略。记录 build/
analyze/factor/multi-RHS time、nnz(L/U)、fill ratio、因子及矩阵实际字节/峰值内存、
回代误差；复数元素字节数不同，不能仅以 fill ratio 跨表示决定优劣。三者均在
原复方程上计算同一 backward error。生产选择与计数合同由 `backend_decision.json`
冻结，未选路径可只保留小型 oracle，不要求长期维护三套生产求解器。

理想 Dirichlet 激励且端口电压只进入 BC 时，`Mu=0`，单位 phasor RHS 是实向量，
仅激励接触的相应 essential 行非零。BC 为 `x-(u/V0-ref/V0)`，故
`dR/du=-1/V0`、`b=+1/V0`；legacy SI 的 V0=1。不是一律 stamp +1。
若 BC 电压映射不是单位斜率，使用映射链式导数；Robin 的直接电压项按 4.3.2 处理。


### 4.6 AC 工作点冻结快照

每次 AC 计算必须从已接受的 DC 点创建不可变 `ACOperatingPointSnapshot`，至少包含：

- packed `x` 和物理 `psi/phin/phip/n/p`；
- per-node electron/hole quasi-Fermi reference fields；
- `DDScalingSpec`、单位系统和装配行恢复物理单位所需 scale；
- Newton 残差诊断使用的固定 block scales；
- frozen electron quantum correction field 及其来源；
- essential 和 natural boundary conditions；
- mesh/material/doping/physics/discretization identity；
- DC bias 和收敛摘要。

J、M、接触反力、解析 `dR/du` 和 finite-difference oracle 必须复用同一个 snapshot。
FD oracle 只在快照内对 BC 数值做 `+/-delta V` 局部算子评估，不运行 Newton，不
重新建立 contact-basin partition，不重新分配 QF reference，也不重新选择 residual
scale。此限制针对固定 x 的偏导 oracle（J/M/B/Qx）；P1 的总 `dQ/dV` 和 P5 的 DC
电导验证可另做 `+/-delta V` 收敛工作点求解，比较物理量，并明确标记为另一类测试。

冻结对象是数值坐标、单位换算、诊断尺度和明确声明的外部固定参数，不是所有物理
量的取值。对扰动状态重新评估迁移率、高场、Fermi/BGN、G/R 等已启用模型，按它们
的真实状态依赖求导；若复用缓存，必须按扰动状态失效，不能强制沿用 DC mobility。
不可微分支记录活动区和离分支边界距离；跨分支的中心 FD 不作为光滑导数硬门。

物理 `J=dR/dx` 与 Newton 求解辅助矩阵分开声明。`carrierDiagonalFloor_`、pseudo-time、
阻尼或其他只为迭代稳定加入的矩阵项不得进入 AC 物理 J。现有近似 Jacobian 分支
（例如省略场/迁移率导数的 avalanche 路径）必须修齐同一装配器导数、选择已审计的
一致路径，或对该 AC 配置 fail-closed；不得建立第二份端电流导数作为补丁。
物理导数可使用局部一致 FD，不要求一律手写解析式，但必须声明方法、步长及误差。

如线性求解需要行/列等化，它属于等价变换：若固定 reference 下 `dx_phys=T*dx`、行尺度为 W，则
`J_scaled=W*J_phys*T`、`M_scaled=W*M_phys*T`、`B_scaled=W*B_phys`，端口导数与解
也要对应转换。Newton 收敛诊断 norm 的 block scale 不一定是 W，二者不得混同。
验收同时报告原方程与求解变换后的残差；不能以等化掩盖物理 J/M 不一致。

AC 使用独立装配器实例及独立线性求解实例，不借用 DC Newton 的可变缓存。
沿 `makeArclengthAssembler()` 的配置工厂路径建立 AC 工厂，但不能直接照搬其当前
`cfg_.carrierDiagonalFloor`：AC 构造时 floor 关闭；其他物理配置保持同一 snapshot，
复制终态逐节点 reference 和必要固定场，不重新执行 contact-basin repartition。
数据依赖需有明确所有权/生存期，不能让装配器引用临时 mesh/material/config。

前置 `snapshot_residual_replay` 门：在同一进程/构建/序列顺序下，保存 DC 终态
原始 solved residual 向量，AC 装配器同 x/bcs 的 solved residual 必须逐字节一致。
它比较原始向量而非日志里的加权 norm；只关闭 Jacobian floor 不应改变 residual。
失败时停止 AC，报告首个不同节点/方程、配置与状态哈希。跨平台/重序列化的回放另有
误差诊断，不冒充 bitwise 门；使用旧 checkpoint 没有原始 residual 时先做可追溯的
DC 重闭合并捕获快照。

模式纪律：当前 pattern 构造和 add 都忽略约束行，且只有一个 boundary signature
缓存。AC 每快照一次批量装配完整 physical J，然后由明确的 BC/gauge stamping
生成 solved J；或至多使用两个固定视图实例。禁止多端/多频循环切换同一实例的视图。
独立计数 `ac.jacobian_pattern_builds<=2`（不含 DC、FD 诊断和换快照），重复频率/
RHS 增量为 0。若复用全局 profiler 验收，必须隔离作用域，不能混入 DC 计数。
绝缘 carrier pin 在当前 residual 的缩放前后各写一次，physical view 必须同时绕过；
solved view 保留单位 gauge 行。

Density-gradient quantum potential 当前在 `CoupledDDAssembler` 中是外层迭代给定的
冻结修正，不是 AC 自洽动态未知量。因此：

- P3 默认对 `electronQuantumPotentialEnabled()==true` fail-closed；
- 后续若允许显式 `quantum_response="frozen"`，必须在 manifest 和每行输出中标记
  `frozen_qp_approximation=true`；
- frozen-QP 结果不能用于认领自洽量子 AC 等价；
- `quantum_response="self_consistent"` 在增加 QP 自洽线性化方程前保持未支持；
  自洽线性化的 QP 可以是代数方程，不应预设必须有量子时间储存项。

### 4.7 输出 schema

建议至少生成：

1. `ac_matrix.csv`，长表，一行一个矩阵元素；
2. `ac_bias_summary.csv`，每个 DC/频率点一行的求解和守恒摘要；
3. `ac_manifest.json`，记录版本、输入哈希、物理合同、端口顺序、单位和门槛；
4. 可选的复数状态场输出，默认关闭。

`ac_matrix.csv` 最小列：

```text
analysis_id
operating_point_id
dc_sweep_contact
dc_bias_V
ac_target_index
is_forced_target
landing_error_V
frequency_Hz
row_terminal
column_terminal
y_real_S_per_m
y_imag_S_per_m
a_S_per_m
c_F_per_m
y_real_S_per_um
y_imag_S_per_um
c_F_per_um
solve_residual_norm
terminal_kcl_residual_S_per_m
qf_reference_snapshot_id
frozen_qp_approximation
```

`ac_manifest.json` 必须记录：

- Vela 提交、构建类型、线性后端和平台；
- 输入配置、网格/TDR、材料和参数文件哈希；
- 端口排序、参考端、激励幅值和正方向；
- 二维归一化；
- DC 工作点状态来源和是否插值；
- 频率列表与 ACCompute 选点规则；
- AC 目标 bias 列表、每个目标是否强制落点及落点误差；
- 物理模型和离散 profile；
- QF reference/scaling snapshot identity 和 `frozen_qp_approximation`；
- 每频率分解次数、右端数、矩阵维度和非零元；
- KCL、charge closure、row/column sum 等质量指标。

另须保存 `schema_version`、`analysis_id`、`operating_point_id`、完整 DC 电压向量、
输入单位版本、全端口/激励子集、boundary/derivative qualification、矩阵/快照哈希、
实际后端和 ordering、成功/失败分解计数、重试次数、误差指标定义版本。
若声明可重建工作点，须保存检查点及其路径/哈希，包含 packed state、逐节点 reference
和必要固定场；只有 snapshot 哈希不能重建状态。不保存检查点的运行显式标记不可重放。
QP 在 MVP 一律拒绝；未来近似字段不意味着已有可用配置。
每个成功 `(operating_point_id, frequency_Hz)` 必须有完整 `N*N` 元素；失败列不填零、
不输出伪完整矩阵。先临时写入、验证完整性再发布，partial manifest 列出缺失目标及原因。

### 4.8 电极准静态 C 与低频 AC 的不同合同

固定 DC 点，记 `B=Ju`，`B_s=Mu`，单位端口列 `e_l`。对已约束且可逆的 J：

```text
J*x0 = -B*e_l
J*x1 = -(M*x0 + B_s*e_l)
x_hat = x0 + j*omega*x1 + O(omega^2)

G0[:,l] = H*x0 + D*e_l
C_AC,0[:,l] = H*x1 + K*x0 + E*e_l
C_Q[:,l] = Qx*x0 + Qu*e_l
```

P1 输出应命名 `C_Q`/`electrode_quasistatic_capacitance`，不可直接写入 AC 的 `c_F_per_m`。
导通器件一般有 `H*x1`，且 K 含接触储存，因此 `C_AC,0 != C_Q` 并不是错误。
纯介电 fixture 二者相等；阻挡栅的相应行在低频极限亦可使用 `C_Q` 对照。一般 PN/
导通 MOS 必须用上述两次实线性求解的低频展开作 oracle，并用频率减半序列验证趋近，
不能用任意固定的“低频 1 Hz”替代松弛时间判断。J 奇异/接近分岔时不使用这个展开。

完整封闭器件的 Gauss 微分恒等式是
`sum_k C_Q[k,l] + dQ_domain/dV_l = 0`，不是所有 `C_Q` 列和都为零。
完整总 Y 的列和为零对应 KCL，行和为零对应全端共同电位平移。后者仅在所有电参考
一起平移、无漏掉的电端口/隐式外接参考时成立。只对无源平衡 fixture 检查互易和
Hermitian 电导部非负；偏置晶体管/近击穿器件不能强加无源性、对称性或电容元素全正。
初版直接求全 N 列；禁止通过强制补行/补列或裁剪元素修复守恒。Reduced 是选定参考端
接地后删对应行列的视图；仅在完整端口、KCL/gauge 已独立验证时才可无歧义重建 full。

### 4.9 验收误差指标（必须可计算）

所有 floor 在看最终 Vela 对比曲线前写入阈值文件，携带单位和适用 fixture：

- 矩阵全局误差：`max_abs(X-Xref)/max(max_abs(Xref), X_floor)`；Y、A、C 分别评分，
  不用大的电导掩盖小电容，也不用整体曲线峰值掩盖某个工作点。
- 显著元素：`abs(Xref_ij)>=max(1e-4*max_abs(Xref), X_floor)`；其他元素使用绝对误差。
- Gauss：`abs(sum Q+Q_domain-Phi_N)/max(sum abs(Q)+abs(Q_domain)+abs(Phi_N), Q_floor)`；
  记录原始 C/m 残差，DC 内部 Poisson 残差必须足够小，不能只看总和偶然抵消。
- 第 l 列 KCL：`abs(sum_k Y_kl)/max(sum_k abs(Y_kl), Y_floor)`；gauge 用行和同式。
- 线性 normwise backward error：`eta=norm_inf(A*z-b)/max(norm_inf(A)*norm_inf(z)+norm_inf(b), tiny)`；
  A 是实块矩阵，复原系统另验。记录 componentwise backward error（逐行绝对值分母）
  与可得的条件估计；小 backward error 不证明病态问题的前向精度。
- FD：冻结快照，至少三档物理步长，检查截断区趋近和舍入平台；光滑小 fixture 的
  显著导数元素门 `1e-6`、全局门 `1e-8`。一般大状态不能普遍要求 FD 达到 `1e-12`。
  `1e-12` 仅用于仿射 Dirichlet BC 的增量式/零基准 fixture：无大数相减，尺度和
  电压换算完全固定，同时检查 BC 行支撑集；正常 DC 基线中心差分使用通用 FD 门。

未指定 floor、工作点、单位、归一化公式的“相对残差通过”不算验收。解析精度门只
用于其离散能精确表示的 fixture；一般器件另报告网格误差，不混同舍入误差。

### 11.1 v4 统一 checklist 合同

本节在拆分后的规范合同中是唯一的跨阶段纪律定义。各阶段引用 gate_id，不重复
定义相互冲突的例外。required gate 只有 evidence-backed pass 才放行；pending/fail/
blocked 均不算通过。not_applicable 必须说明与 qualified_scope 的对应关系，由审核批准。
诊断值可以保留，但不进入正式通过统计；入账本不自动豁免验收。变更材料、状态、
单位、导数、布局或门槛按 dependency 使旧证据失效并重验。

| gate_id | required_by | 最小证据 |
| --- | --- | --- |
| material_support_decision | WP-R/P1/P2/P5 | 几何与材料支持分离、联合电荷/储存恒等式、决策和重验范围 |
| gate_material_contract | P0/P5 | TDR 接触/区域映射、poly DD 资格或明确的金属替代范围 |
| reference_derivative_contract | P0/P5 | 原始/有效 Math、版本默认、AC 专用语义与受控变体 |
| snapshot_residual_replay | WP-R/P3/P4 | 同构建原始 residual bitwise 复现、state/config/reference 哈希 |
| physical_derivatives | WP-R/P3/P4 | 无 Newton-only floor，live residual FD 与配置资格矩阵 |
| storage_charge_identity | P2/P4 | 移动电荷与 S/M 一致、接触体积抵消、单位恢复 |
| boundary_row_policy | WP-R/P2/P4 | essential/Robin/gauge/sheet 分类及完整导数支撑集 |
| ac_pattern_cache | WP-R/P3 | 独立 AC 作用域，每 snapshot pattern builds<=2，RHS/频率不重建 |
| backend_equivalence_and_reuse | P3 | 三路比较、选中后端 count 语义、同矩阵多 RHS 复用 |
| forced_target_and_observer | P4 | 实际调用路径、目标列表/命中、只读发布、拒绝/初始/恢复事件 |
| dc_alignment | P5 | 每目标端电流/psi/log n,p 合格且覆盖率完整 |
| separate_ac_metrics | P5 | A/C 分区门及原始/受控合同，完整 Y/KCL/gauge/低频门 |
| complete_output_compatibility | P4/P7 | full matrix/单位/partial 状态、原始数据、legacy 行为 |
| scope_and_model_qualification | 各阶段 | 不支持配置 fail-closed，QP 拒绝，非默认体积/后端范围显式 |

`stage_acceptance.json` 最小结构（示意，不是执行结果）：

```json
{
  "schema_version": 1,
  "plan_version": "v4",
  "stage_id": "P0.5",
  "work_package_id": "WP-R",
  "source_commit": null,
  "input_hashes": {},
  "qualified_scope": [],
  "required_gates": ["material_support_decision", "snapshot_residual_replay", "ac_pattern_cache"],
  "gates": {
    "snapshot_residual_replay": {
      "status": "pending",
      "evidence": [],
      "metrics": {},
      "reason": "not executed"
    }
  },
  "tests_selected": 0,
  "decision": "pending"
}
```

执行版 schema 要求 required_gates 每项都有记录，证据带哈希/单位/阈值版本，缺项即
未完成；上例是待填草稿。阶段 ID 与 WP ID 分开保存，`P0.5 ≡ WP-R`，P2.5 是独立
介电资格门。AC-L 映射以实施计划 §3.3 为准；本规范不单独宣布任何级别通过。
