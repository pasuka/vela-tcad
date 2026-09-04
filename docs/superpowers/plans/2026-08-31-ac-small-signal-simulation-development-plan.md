# Vela TCAD AC 小信号仿真功能开发方案

日期：2026-08-31

状态：外部审核修订稿 v2；5 个阻塞项已纳入冻结决策，仍不授权直接实施

Vela 基线：`270f5258a3c6354b2b3564af98ae4ab9434de2ca`

Sentaurus 版本：T-2022.03-SP2

Sentaurus 参考工程：
`/atctools/Synopsys/tcad/T-2022.03/tcad/T-2022.03-SP2/Applications_Library/GettingStarted/sdevice/AC`

## 1. 执行摘要

本方案拟为 Vela 增加二维半导体器件的频域小信号 AC 能力。第一主验收对象是
Sentaurus `GettingStarted/sdevice/AC` 中的 nMOS 示例。该示例在直流栅压工作点上
对 Poisson--电子--空穴漂移扩散系统做 1 MHz 正弦小信号线性化，提取四端口
导纳矩阵：

```text
Y(omega) = A + j * omega * C
```

其中 `A` 是同相电导/跨导矩阵，`C` 是正交电容矩阵。示例同时覆盖：

- 显式 MixedMode 器件和理想电压源系统：`AC_des.cmd`；
- 隐式自动端口系统：`AC1_des.cmd`；
- DC 栅压延续和 ACCompute 子采样；
- 多端口 `Node`、电压源 `Exclude`、单频和一般频扫语义；
- `ACExtract` 全端口矩阵输出；
- 显式和隐式系统的总栅电容一致性检查。

Vela 已经具备 DC 工作点、耦合漂移扩散残差与解析 Jacobian、多端接触、稳态
端电流、稀疏实数线性求解和大部分相关物理模型，因此不需要另写一套半导体 DC
求解器。真正缺失的是：

1. 与连续性方程一致的动态储存算子及其 Jacobian；
2. 物理电极电荷、位移电流和端口电流导数；
3. `J + j*omega*M` 频域线性系统；
4. 多端口激励、完整 Y/A/C 矩阵、频率扫描和稳定输出契约；
5. 隐式端口系统和后续的一般 MixedMode MNA；
6. 解析、开源 TCAD 和 Sentaurus 三层验证链。

推荐总体依赖修订为：

```text
P0a Sentaurus oracle、符号/单位/强制落点/schema 冻结
P0b 公开解析 fixture（与 P0a 并行）
  -> WP-R 统一的未约束残差/J/M 接触行访问接口
  -> [P1 残差行电极电荷与准静态 C
      || P2 装配器内动态储存 M
      || WP3a 人工矩阵实块 builder]
  -> [P1 + WP3a] P2.5 纯介电 Poisson AC
  -> [P2 + P2.5 + WP3b/3c] P3 单器件隐式 AC 核心 J+jwM
  -> P4 全 Y/A/C、频扫、配置与输出
  -> P5 DEVSIM/Sentaurus 交叉验证
  -> P6 显式 MixedMode MNA
  -> P7 性能、文档、前端与长期回归
```

方案的核心决策是：

1. **不把当前 `cv_quasistatic` 当成 AC。** 它目前是区域体电荷的相邻偏压差分，
   不包含频率、相位、位移电流或完整多端口响应。
2. **先做单器件隐式系统，再做一般 MixedMode。** `AC1_des.cmd` 是 AC 数值核心的
   第一验收对象，`AC_des.cmd` 是一般电路耦合的下一阶段验收对象。
3. **MVP 使用实数 2x2 块系统，但不预先冻结内存排列。** 暂不要求把全仓库
   数值类型改成复数；stacked 和 interleaved 由 P3 的 fill ratio 实测决定。
4. **储存算子与未来瞬态共享设计，但不先实现完整瞬态积分器。** AC 只要求
   状态储存项及其一致线性化。
5. **端口语义先于曲线拟合。** 先冻结端口顺序、参考端、正方向、二维宽度归一化、
   A/C/Y 定义和输出 schema，再比较 Sentaurus 数值。
6. **所有跨工具数值门槛先标为暂定。** P0 导出官方 oracle 后必须生成并审核
   `threshold_freeze.json`，不得为通过最终曲线临时改门槛。
7. **端电流和端电荷的生产主路径统一为装配一致的反力提取。** 在保留自然边界项、
   但不覆盖 essential 行的物理残差/J/M 上，对接触节点对应方程行加权求和；独立
   SG 电流和显式 `D dot n` 只作交叉验证 oracle。
8. **AC 目标偏压必须强制落点。** 自适应 DC 步长应裁剪到下一个待计算 AC 目标，
   中间接受点可以存在但不替代目标点，也不允许靠插值评分。
9. **AC 线性化快照冻结数值坐标。** QF reference field、DD scaling、残差诊断尺度
   和 frozen quantum correction 在 J/M/dRdu/FD oracle 之间不得重分配或重算。

### 1.1 外部审核意见处置记录

| 编号 | 处置 | v2 决定 |
| --- | --- | --- |
| B1 | 采纳，改变 P4 核心架构 | 端电流和线性化改为未约束连续性残差/J/M 接触行提取；SG kernel 降为 oracle |
| B2 | 采纳，改变 P1 核心架构 | 电极电荷改为未约束 Poisson 残差/J 接触行反力；不独立重推生产 `D dot n` |
| B3 | 采纳并前移到 P0 | ACCompute 目标偏压成为强制 DC 落点；P5 只比较精确共同目标点 |
| B4 | 采纳，设为 P3 前置硬任务 | `LinearSolver` 必须新增 `numericFactorizationCount()`；AC 使用专用 solver 和不可变块矩阵 |
| B5 | 采纳 | 冻结 QF reference/scales/QP；FD 只局部扰动 BC，不经过 Newton 重求解；QP 默认 fail-closed |
| N1 | 采纳 | 储存项实现归属 `CoupledDDAssembler`，独立头文件只承载接口/结果类型 |
| N2 | 采纳 | 增加 essential、thermionic/Schottky、insulating pin、sheet charge 等边界行分类表 |
| N3 | 采纳 | P3 前新增 P2.5 纯介电 Poisson AC |
| N4 | 采纳 | stacked/interleaved 块布局由 fill ratio 决策，不在方案阶段冻结 |
| N5 | 采纳 | 单位激励只是导数归一化；不暴露为具有非线性含义的物理旋钮 |
| N6 | 采纳 | P3 增加 `1e3--1e12 Hz` 重掺杂 MOS 调理扫描 |
| N7 | 采纳 | 测试矩阵增加 avalanche-on 近击穿二极管 AC |
| N8 | 采纳 | 使用只读 `onAcceptedPoint` 观察者，避免侵入式重构 `DCSweep.cpp` |

这张表记录的是方案修订决定，不是功能已经实现的声明。

## 2. 背景与问题定义

### 2.1 Sentaurus AC 示例解决什么问题

该工程构造二维四端 nMOS。`nMOS_dvs.cmd` 给出约 50 nm 栅长、薄栅氧、p 型体区、
重掺杂多晶硅栅和高斯源漏/延伸区，并通过半结构镜像生成 source、drain、gate、
substrate 四端器件。

`AC_des.cmd` 的工作点路径为：

1. 求解 Poisson；
2. 求解 `Poisson Electron Hole`；
3. 栅压降至 `-3 V`；
4. 栅压从 `-3 V` 扫至 `+3 V`；
5. 在 `ACCompute(Time=(Range=(0 1) Intervals=40))` 选定的位置做 AC；
6. 频率固定为 `1e6 Hz`；
7. 端口为 s/d/g/b，外接理想电压源在 AC 中被 `Exclude`；
8. 输出全端口 A/C 元素。

主动 DC 物理包括 Fermi、OldSlotboom、掺杂依赖迁移率、准费米梯度高场饱和、
Enormal 和掺杂/温度相关 SRH。Plot 列表中的 avalanche、BTBT 或 quantum potential
字段不表示这些模型已经激活。

`AC1_des.cmd` 使用同一器件和偏压流程，但通过 `ImplicitACSystem` 自动建立端口节点。
两个 visual 脚本分别绘制 `c(g,g)` 和 `c(N_gate,N_gate)`，再检查结果重合。

### 2.2 小信号 AC 与准静态 C--V 的区别

设离散 DC 残差、储存量和端口电流为：

```text
R(x, u) = 0
S(x, u)
I(x, u)
Q(x, u)
```

`x` 是器件未知量，`u` 是端口电压。围绕 DC 工作点做一阶频域线性化：

```text
(Jx + j*omega*Mx) * x_hat = -(Ju + j*omega*Mu) * u_hat

i_hat = (Ix + j*omega*Qx) * x_hat
      + (Iu + j*omega*Qu) * u_hat
```

其中：

- `Jx = dR/dx`，Vela 已有主要实现；
- `Mx = dS/dx`，当前缺失；
- `Ix/Iu` 不另写平行 SG 导数，而由接触连续性反力对未约束 J/M 的行提取得到；
- `Qx/Qu` 由接触 Poisson 反力对同一未约束 J 的行提取得到，当前 `TerminalCharge`
  不能直接替代；
- 对每个被激励端口求解一次，可得到 Y 的一列；
- 同一频率的所有端口列应复用同一次数值分解。

准静态 `dQ/dV` 只描述 `omega -> 0` 的储存响应，不能给出一般频率下的载流子滞后、
导电跨导或相位。因此 P1 的全端口准静态 C 是 AC 的资格门和低频 oracle，不是最终
AC 求解器。

### 2.3 Vela 当前代码基线

| 能力 | 当前证据 | 结论 |
| --- | --- | --- |
| DC 偏压扫描与延续 | `include/vela/simulation/DCSweep.h`、`src/simulation/DCSweep.cpp` | 可复用 |
| 耦合 DD 残差和解析 Jacobian | `CoupledDDAssembler::residual/assembleJacobian` | 可复用 |
| 未知量布局 | `psi, phin, phip`，共 `3*N` | AC 储存导数必须尊重现有缩放和参考坐标 |
| 实数稀疏直接求解 | `LinearSolver` 支持模式和相同矩阵数值分解复用 | 可复用；MVP 用实块系统 |
| 多端接触和稳态电流 | `ContactCurrent` | 可复用，但需增加线性化 |
| 区域体电荷 C--V | `TerminalCharge` 和 `cv_quasistatic` | 仅为 legacy/prototype，不能直接认领 AC |
| 单变量 DC 外接负载 | `CoupledLoadLine` | 只能参考 bordered-system 设计，不是一般 MNA |
| Fermi/OldSlotboom/迁移率/SRH | 已有相关模型 | 可复用；须审计 AC 所需导数 |
| 复数稀疏类型 | `Types.h` 只有 `double` | 缺失；MVP 不要求新增全局复数类型 |
| 动态储存/质量矩阵 | 核心中无对应算子 | 缺失 |
| 频扫、Y 矩阵、ACExtract | 无执行路径 | 缺失 |
| SDE 高斯植入 | 前端 fail-closed，内部 mesh builder 仅常数掺杂 | 非 AC 核心；短期用 TDR 绕过 |

`ContactCurrent::computeFromResidual` 已经证明接触连续性残差行求和可产生端电流；
`NewtonSolver::makeArclengthContactCurrentFunctional` 也已经通过
`J.transpose()*residualWeights` 构造过同一函数的状态导数。AC 应把这个先例升级为
受支持的统一接触反力接口，而不是为 mobility/Fermi/BGN/high-field/avalanche 再写
第二份 SG 导数。现有方法传入空边界条件，P1 前必须扩展为“保留 thermionic 等自然
边界，但跳过 Dirichlet essential 行覆盖”的明确装配模式。

v2 审核时已复核的代码事实：

| 事实 | 当前符号/位置 | 对方案的约束 |
| --- | --- | --- |
| DC 残差行端电流 | `ContactCurrent::computeFromResidual` | B1 有仓库内先例 |
| 残差函数状态导数 | `NewtonSolver::makeArclengthContactCurrentFunctional` | 可用 `J^T*w`，无需第二套 SG 导数 |
| per-node QF reference | `CoupledDDAssembler::setQuasiFermiReferenceFields` | AC FD 必须冻结 reference field |
| contact-basin repartition | `NewtonSolver` warm-start/repartition 路径 | FD 不得经过正常 Newton solve |
| thermionic natural flux | `CoupledDDBoundaryConditions::thermionic` 和 assembler residual/J | physical view 不能传空 bcs |
| numeric factor cache | `LinearSolver::valuesMatchFactorization` | 同矩阵值必须 bitwise 不变 |
| 可见对象计数器 | 只有 `patternAnalysisCount()` | P3 必须新增 numeric count getter |

现有 `cv_quasistatic` 在相邻 DC 点计算：

```text
C = (Q[k] - Q[k-1]) / (V[k] - V[k-1])
```

其中 Q 是用户指定区域中的 `q*(p-n+Nd-Na)` 体积分，还可按接触距离裁剪。它既不是
严格的金属电极高斯通量电荷，也不保证所有端电荷闭合。开发 AC 时不得悄悄改变
现有输出语义；应新增明确的电极电荷方法，并保留 legacy 模式兼容性测试。

### 2.4 开源实现和文献基线

- Sentaurus Training 定义 `Y=A+j*omega*C`，并说明显式/隐式 AC 系统、Node、
  Exclude 和 ACCompute；
- DEVSIM 支持 DC 工作点上的 small-signal AC，并在二极管示例中由端口电流虚部
  除以 `omega` 得到电容；
- Genius-TCAD-Open 提供 DDMAC、ACSWEEP 和器件--电路混合模式；
- Laux 1985 比较瞬态 FFT、增量电荷分区和正弦稳态方法；
- SISPAD 1999 给出将时间导数写为 `j*omega` 项并在 DC 点解复线性系统的路径。

这些来源用于约束数值架构，不用于假定与 Sentaurus 私有离散完全相同。

## 3. 目标、非目标和成功等级

### 3.1 总目标

建立一套可重复、可审计的二维器件 AC 流程，能够：

1. 从 Vela 已收敛的 DC 工作点构造频域小信号系统；
2. 对一个或多个正频率求解端口响应；
3. 输出完整多端口复导纳 Y、同相矩阵 A 和电容矩阵 C；
4. 支持固定工作点和 DC 偏压扫上的 ACCompute 子采样；
5. 支持第一阶段隐式端口系统；
6. 后续支持理想电压源、电阻和电容的显式 MixedMode；
7. 在解析 fixture、DEVSIM 和 Sentaurus 示例上形成分层回归；
8. 使动态储存算子未来可被 BDF/TRBDF 等瞬态功能复用。

### 3.2 非目标

第一版不包括：

- 完整瞬态时间积分器；
- RF 噪声、S 参数、谐波平衡或大信号周期稳态；
- 电感、互感、受控源、传输线和一般 SPICE 网表；
- 3D AC；
- 热 AC、量子势动态方程或陷阱动态；
- 自动复制 Sentaurus SDE 的完整几何/掺杂语言；
- 通过经验缩放电荷、电流或频率来拟合 Sentaurus；
- 把偏置晶体管强制要求为互易矩阵；
- 把 Plot 中存在但未激活的物理量纳入等价范围。

### 3.3 成功等级

| 等级 | 定义 |
| --- | --- |
| AC-L0 | Sentaurus 官方工程、输入、输出、哈希、端口和单位已封存，阈值已冻结 |
| AC-L1 | 物理电极电荷和全端口准静态 C 通过解析与守恒测试 |
| AC-L2 | 动态储存算子及其 Jacobian 通过有限差分和缩放审计 |
| AC-L3 | 单器件隐式 AC 可解单频和频扫，并输出完整 Y/A/C |
| AC-L4 | DEVSIM、Sentaurus `AC1_des.cmd` 和低频准静态极限通过 |
| AC-L5 | 显式 MixedMode V/R/C 系统通过，复现 `AC_des.cmd` |
| AC-L6 | 性能、文档、CI、错误诊断和长期回归全部固化 |

报告和发布说明必须声明达到的最高等级，不得仅写“已支持 AC”。

## 4. 必须先冻结的数学和数据契约

### 4.1 端口方向和矩阵索引

建议冻结如下约定：

- `Y[row_terminal, column_terminal] = dI_row / dV_column`；
- 端电流流入器件为正；
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

可同时输出 `*_per_um` 便利列，换算必须为规范每米值乘 `1e-6`。如果用户提供有限
`depth_m`，则另输出总量并在 manifest 中记录归一化方式。禁止在同一列中混用
总量和每宽度量。

### 4.3 物理电极电荷

生产 AC 不使用区域归属电荷，也不另写一套独立边界几何积分作为端电荷。冻结定义为
接触 essential Poisson 方程的离散反力：

```text
Q_k = s_Q * w_psi,k^T * R_physical_unconstrained(x, natural_bcs)

dQ_k/dx = s_Q * J_physical_unconstrained^T * w_psi,k
```

`w_psi,k` 在接触 k 的 Poisson 行为 1，其余行为 0；`s_Q` 恢复装配器的物理单位和
最终正号。这里的 `unconstrained` 是“保留体电荷、FVM 边耦合、fixed/interface
sheet charge 和自然边界，跳过 essential Dirichlet 行替换”，不是简单传空 bcs。

这个反力与显式 `integral(D dot n)` 在连续层面等价，但生产实现直接复用
`CoupledDDAssembler` 的 node volume、mixed-Voronoi/非钝角策略、边 coupling 和
缩放。独立 `D dot n` 计算仅作为规则网格解析 fixture 的 oracle，不得成为生产主
路径。

具体 `s_Q` 正号由“电流流入器件为正”约定、平行板电容器和全域 Gauss 闭合冻结。
对半导体 ohmic 端口，总端电流包含连续性反力和 Poisson 反力的时间导数；对绝缘
栅，连续性分量应为零而位移分量非零。

v2 的候选离散端口公式为：

```text
i_carrier,k_hat = s_I * w_cont,k^T
                  * ((J_phys + j*omega*M_phys) * x_hat
                     + (Ju_phys + j*omega*Mu_phys) * u_hat)

i_displacement,k_hat = j*omega * s_Q * w_psi,k^T
                       * (Jpsi_phys * x_hat + Jpsi_u_phys * u_hat)

i_total,k_hat = i_carrier,k_hat + i_displacement,k_hat
```

这里 `w_cont,k` 内含电子/空穴符号组合。P0.5 必须从现有 residual 单位推导
`s_I/s_Q`，P2.5/P4 必须用解析平板和 DC `computeFromResidual` 锁定符号。还必须由
外部审核或 manufactured transient fixture 排除“接触控制体 M 反力”和 Poisson
位移项的双计数；在这一点关闭前，公式是候选冻结式而不是完成声明。

区域体电荷可保留为诊断量，用于检查：

```text
sum_k Q_electrode,k + Q_domain = 0
```

但不能用任意“接触半径”分区替代电极电荷。

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
- 具体电子/空穴符号以 Vela 残差定义和 manufactured transient identity 冻结，
  不允许仅凭教科书符号手填；
- inactive 模型不得改变 M；未来动态陷阱/热/量子方程通过独立 block 扩展。

边界行冻结表：

| 行类型 | solved J 行 | solved M 行 | physical reaction J/M |
| --- | --- | --- | --- |
| Ohmic/metal-gate essential Dirichlet | 替换为代数约束 | 0 | 跳过替换，保留物理行供端口提取 |
| Thermionic/Schottky Robin | 保留自然通量 | 保留连续性储存 | 与 solved physical row 相同 |
| 绝缘区内部 phin/phip pin | 代数 gauge | 0 | 不认作端口电流 |
| Fixed/interface sheet charge | 进入 Poisson 行 | 0 | 进入 Poisson 端口反力和 Gauss closure |
| 未接触内部 DD 行 | 正常物理行 | 正常储存 | 不直接进入端口权重 |

### 4.5 实数块系统

MVP 将：

```text
(J + j*omega*M) z = b
```

在数学上可展开为：

```text
[ J       -omega*M ] [Re(z)] = [Re(b)]
[ omega*M  J       ] [Im(z)]   [Im(b)]
```

这只是数学排列，不冻结内存布局。P3 必须对至少一个重掺杂 MOS fixture 比较：

- stacked `[Re(all); Im(all)]`；
- interleaved `[Re(x0), Im(x0), Re(x1), Im(x1), ...]`。

以 symbolic fill、numeric fill、factor time、solve time、实现复杂度和残差共同冻结
布局。两种布局必须代表同一个块系统，不能以改变方程缩放换取 fill 优势。

隐式理想 Dirichlet 激励下 `Mu=0`，单位 phasor RHS 是实向量，且只在被激励接触的
essential BC 行非零。实现必须提供独立单元测试确认块符号，不能依赖最终 Cgg 曲线
间接发现错误。

是否在 AC-L6 后引入 `SparseMatrix<std::complex<double>>`，应由性能数据决定，不在
MVP 中预先扩散复数类型。

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
scale。

Density-gradient quantum potential 当前在 `CoupledDDAssembler` 中是外层迭代给定的
冻结修正，不是 AC 自洽动态未知量。因此：

- P3 默认对 `electronQuantumPotentialEnabled()==true` fail-closed；
- 后续若允许显式 `quantum_response="frozen"`，必须在 manifest 和每行输出中标记
  `frozen_qp_approximation=true`；
- frozen-QP 结果不能用于认领自洽量子 AC 等价；
- `quantum_response="self_consistent"` 在增加 QP 动态/线性方程前保持未支持。

### 4.7 输出 schema

建议至少生成：

1. `ac_matrix.csv`，长表，一行一个矩阵元素；
2. `ac_bias_summary.csv`，每个 DC/频率点一行的求解和守恒摘要；
3. `ac_manifest.json`，记录版本、输入哈希、物理合同、端口顺序、单位和门槛；
4. 可选的复数状态场输出，默认关闭。

`ac_matrix.csv` 最小列：

```text
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

## 5. 目标软件架构

### 5.1 分层

```text
物理/离散层
  CoupledDDAssembler physical/solved row views
  CoupledDDAssembler storage S/M
  ContactReactionFunctional (Poisson + continuity row weights)

频域数值层
  ACBlockSystemBuilder
  ACSmallSignalSolver

工作流层
  ACConfig / ACResult / ACWriter
  ACSweepRunner

电路层（后续）
  CircuitMNA
  MixedModeACSolver
```

规则：

- 物理算子不得依赖 JSON、CSV 或 `DCSweep.cpp`；
- 数值层接受已收敛工作点、装配器、边界条件和端口描述；
- 工作流层负责 DC 扫描、ACCompute 选点、频率调度和输出；
- 电路层只通过稳定的器件端口线性化接口耦合；
- 每层必须能被小型 Catch2 fixture 独立调用。

### 5.2 建议新增文件

| 文件 | 责任 |
| --- | --- |
| `include/vela/equation/CoupledDDStorage.h` | `S/M/Mu` 结果类型和行分类合同 |
| `include/vela/post/ContactReactionFunctional.h` | Poisson/连续性接触行权重、物理单位和反力接口 |
| `src/post/ContactReactionFunctional.cpp` | 从装配器 physical residual/J/M 提取 Q/I 及导数 |
| `include/vela/solver/ACBlockSystem.h` | 实 2x2 频域块矩阵构造 |
| `src/solver/ACBlockSystem.cpp` | 频率块装配和质量检查 |
| `include/vela/simulation/ACSmallSignal.h` | 配置、结果和单器件 solver API |
| `src/simulation/ACSmallSignal.cpp` | 多端激励、频扫和 Y/A/C 提取 |
| `include/vela/io/ACWriter.h` | AC CSV/JSON 输出契约 |
| `src/io/ACWriter.cpp` | long-form matrix 和 manifest writer |
| `include/vela/circuit/CircuitMNA.h` | P6 的 V/R/C MNA 数据结构 |
| `src/circuit/CircuitMNA.cpp` | P6 的电路 stamping |

建议新增测试：

```text
tests/test_contact_reaction_functional.cpp
tests/test_coupled_dd_storage.cpp
tests/test_ac_block_system.cpp
tests/test_ac_small_signal.cpp
tests/test_ac_output.cpp
tests/test_ac_mixed_mode.cpp        # P6
tests/regression/test_sentaurus_ac.py
tests/regression/test_devsim_ac.py
```

### 5.3 预计修改的现有文件

- `CMakeLists.txt`：加入新源码和测试；
- `include/vela/simulation/DCSweep.h`：仅增加共享工作点/选点接口或 AC 配置挂接；
- `src/simulation/DCSweep.cpp`：只增加只读 `onAcceptedPoint` 观察者和目标偏压强制落点，
  不侵入式重构约 8400 行的求解循环；
- `include/vela/solver/LinearSolver.h`、`src/solver/LinearSolver.cpp`：P3 前必须增加对象级
  `numericFactorizationCount()`；如需显式 `factorize/solveFactored` 或矩阵 RHS API，
  另以测试和性能证据决定；
- `include/vela/equation/CoupledDDAssembler.h`、`src/equation/CoupledDDAssembler.cpp`：
  增加 physical/solved row view、储存 S/M 和边界行分类；
- `include/vela/post/ContactCurrent.h`、`src/post/ContactCurrent.cpp`：使既有
  `computeFromResidual` 委托统一 `ContactReactionFunctional`，SG edge 实现继续作为
  独立 DC oracle；
- `include/vela/simulation/CurveSweep.h`、`src/simulation/CurveSweep.cpp`：增加
  `ac_small_signal` 类型，或由独立 analysis type 路由；
- `src/tools/vela_example_runner.cpp`：接入 AC 工作流；
- `examples/`：加入无授权限制的小型 MOSCAP/PN AC 示例。

不建议把所有新功能继续堆入已经很大的 `DCSweep.cpp`。

### 5.4 建议配置草案

```json
{
  "sweep": {
    "mode": "ac_small_signal",
    "contact": "gate",
    "start": -3.0,
    "stop": 3.0,
    "initial_step": 0.02,
    "ac": {
      "system": "implicit",
      "nodes": ["source", "drain", "gate", "body"],
      "reference_node": "body",
      "frequencies": {
        "start_Hz": 1000000.0,
        "stop_Hz": 1000000.0,
        "points": 1,
        "spacing": "decade"
      },
      "compute": {
        "mode": "normalized_intervals",
        "intervals": 40,
        "landing": "force"
      },
      "output": {
        "matrix_csv": "ac_matrix.csv",
        "summary_csv": "ac_bias_summary.csv",
        "manifest_json": "ac_manifest.json"
      }
    }
  }
}
```

配置规则：

- `frequencies_Hz` 显式数组和范围对象二选一；
- 所有频率必须有限且大于零；
- `nodes` 不得重复，必须解析为真实电接触；
- `reference_node` 必须在 nodes 中；
- AC 求导内部固定使用单位 phasor，不提供会暗示大信号含义的 excitation amplitude
  物理旋钮；有限差分步长只存在于 diagnostics/test 配置；
- P3 仅接受 `system=implicit`，其他值 fail-closed；
- P6 才接受显式 circuit/netlist；
- `compute` 生成强制 AC 目标偏压；DC 自适应步必须裁剪到下一个待命中目标。中间
  接受点仍可用于延续，但不得替代目标或参与 P5 评分；
- `normalized_intervals` 在开始 DC 扫描前解析为不可变物理 bias 列表，并写入
  manifest；还应支持显式 `bias_points_V`，两者二选一；
- 输出默认 full matrix，reduced matrix 只是派生视图；
- legacy `cv_quasistatic` 的默认和列名在迁移期保持不变。

## 6. 分阶段实施内容和验收标准

## P0：Sentaurus oracle、合同和阈值冻结

### 目标

在写 AC 数值代码前，消除端口、符号、单位、记录数和参考数据不确定性。

### 实施内容

1. 在 VM 原样运行 `nMOS_dvs.cmd -> AC_des.cmd -> AC1_des.cmd`；
2. 记录 Sentaurus 版本、主机、输入文件 SHA-256、材料参数文件和日志；
3. 导出显式和隐式 AC 的全部 bias/frequency/A/C 元素，不只导出 Cgg 图片；
4. 确认 `Intervals=40` 实际产生的 AC 记录数和端点行为；
5. 确认 `.plt` 中端口行列顺序、正方向、A/C 单位和二维宽度归一化；
6. 确认 `Y=A+j*omega*C` 的数值关系；
7. 比较显式与隐式系统，并形成逐元素差异报告；
8. 把 Sentaurus normalized ACCompute 位置映射为精确物理 bias 目标，确认端点、目标数、
   落点精度和中间 DC 点行为；
9. 建立 `reference_staging/sentaurus_ac_t2022_03_sp2/`，不提交受许可约束的原始
   工件；
10. 创建 `ac_contract.schema.json`、`threshold_freeze.json` 和读取器测试；
11. 建立一个不依赖 Sentaurus 的公开解析 fixture 数据集供 CI 使用；该任务与 VM
    复跑并行，不等待受许可 oracle；
12. 冻结 Vela 的目标落点规则：自适应步裁剪到待命中 AC bias，精确接受后才触发
    `onAcceptedPoint`，严禁用插值结果触发 AC。

### 必须交付

```text
reference_manifest.json
input_hashes.json
sentaurus_ac_full_matrix.csv
sentaurus_ac1_full_matrix.csv
sentaurus_explicit_implicit_compare.json
ac_contract.schema.json
threshold_freeze.json
```

### 验收标准

- 同一官方工程可从干净工作目录重复运行；
- 所有输入和输出都有哈希，原始文件只读保留；
- 每条矩阵记录能唯一定位 bias、frequency、row、column；
- 显式/隐式 Sentaurus 的全矩阵全局归一差异暂定不大于 `1e-6`；若官方自身高于
  此值，必须先解释并按实测冻结，不能扩大 Vela 容差掩盖；
- `A + j*omega*C` 重建值与导出的复 Y（若可导出）在输出精度内闭合；
- schema 对缺端口、重复端口、非正频率、混合单位和缺哈希 fail-closed；
- 用初始步约 `0.03 V`、目标点不在自然步长网格上的 fixture，证明所有 AC 目标均
  被步长裁剪精确命中，且中间 DC 接受点不改变目标顺序；
- `onAcceptedPoint` 是只读观察者；回调不得改写接受状态、下一步预测器或 retry
  决策；
- P0 未通过不得进入最终跨工具等价声明。

## P0.5：统一 physical row view 和接触反力接口

### 目标

先冻结 B1/B2 的共同基础，使 P1 电荷和 P2 储存可以安全并行。

### 实施内容

1. 为 `CoupledDDAssembler` 明确区分：
   - `solved` view：包含 essential 行替换，供 Newton/AC 内部状态求解；
   - `physical` view：保留体项、边通量、自然边界和静态 sheet charge，跳过
     essential 行替换，供端口反力提取；
2. physical view 必须接受完整 bcs，并只移除 essential row replacement，不能用
   `CoupledDDBoundaryConditions{}` 代替；
3. 定义 `ContactReactionFunctional` 的 Poisson、electron continuity、hole continuity
   权重和物理单位恢复因子；
4. 增加接收完整 bcs 和 `ResidualView::Physical` 的
   `ContactCurrent::computeFromResidual` overload 并委托统一接口；保留旧 overload
   供无自然边界的兼容调用，但不得在 AC 中使用；
5. 暴露 `value`、`stateDerivative` 和后续 `storageDerivative`，不暴露第二套 SG
   production derivative；
6. 为 ohmic、metal gate、thermionic/Schottky、insulating pin 和 interface sheet
   charge 建立行分类测试。

### 验收标准

- 现有 ohmic fixture 中，统一反力端电流与 `computeFromResidual` 逐接触一致到机器
  舍入；
- 统一反力端电流与独立 SG `computeDetailed` 在适用 DC fixture 上通过现有 KCL
  门；
- natural thermionic 项在 physical view 中保留，essential replacement 不存在；
- solved view 的现有 residual/Jacobian 测试字节语义不变；
- 接触权重、block offset、行类型和物理 scale 都有显式 API，不由调用者手写 `N`
  和 `2*N`；
- P0.5 未冻结前，P1/P2 不进入生产实现。

## P1：物理电极电荷与全端口准静态 C

### 目标

建立不依赖区域人为分区的端电荷定义，为位移电流和低频极限提供 oracle。

### 实施内容

1. 通过 P0.5 的 Poisson contact reaction functional 计算 Q，不另写生产边界积分；
2. 从 physical Poisson residual/J 行提取 Q 和 Qx，并恢复物理单位；
3. 支持半导体 ohmic、金属栅/绝缘层边界和二维宽度归一化；
4. 输出电极反力电荷、由同一 physical Poisson 装配源项求和得到的域内净电荷和
   Gauss closure；
5. 使用冻结 snapshot 的中心有限差分作为 Qx oracle；
6. 独立规则网格 `D dot n` 只作为平行板 oracle；
7. 由每个端口单位扰动生成 full quasi-static C matrix；
8. 保留当前 `TerminalCharge` 和 legacy `cv_quasistatic`，新增明确的
   `charge_method=poisson_reaction`；
9. 为 full/reduced C 矩阵增加输出和守恒诊断。

### 测试

- 一维/二维平行板电容器；
- 有氧化层 MOSCAP；
- 两端 PN 结；
- 多端对称结构；
- 不同网格细化和不同 `depth_m`；
- 电位整体平移的 gauge invariance；
- legacy C--V 列名和数值不变测试。

### 验收标准

- 均匀平行板 fixture 的 Q 和 C 与解析值相对误差不大于 `1e-8`，或达到该离散
  fixture 的机器精度预期；
- 每个 fixture 的 `|sum(Q_electrode)+Q_domain|` 归一残差不大于 `1e-10`；
- Poisson reaction J 行提取的 Qx 与中心差分在显著元素上相对误差不大于 `1e-6`，全局归一
  误差不大于 `1e-8`；
- C 矩阵共同模式/电荷守恒行列和归一残差不大于 `1e-8`；
- 逐次减半扰动时，中心差分在截断误差区呈二阶收敛；
- legacy `cv_quasistatic` 测试全部通过；
- 无接触半径、区域归属参数或独立复制的 FVM geometry 参与生产电极电荷定义；
- mixed-Voronoi/钝角网格 fixture 的 Q/closure 使用与 Poisson 残差相同的 node volume
  和 coupling，不能通过降低 closure 门槛绕过。

## P2：动态储存算子和一致 Jacobian

### 目标

构造 `S/M/Mu`，使 DD 频域系统在现有 `psi/phin/phip` 坐标、缩放和边界条件下
数学闭合。

### 实施内容

1. 在 `CoupledDDAssembler` 上定义 storage `S/M/Mu` 公共接口，独立头文件只放结果
   类型和行分类；
2. 使用与连续性源项完全相同的 node volume 策略，为电子/空穴控制体装配储存量；
3. 使用 CarrierStatistics 一致导数建立对 `psi/phin/phip` 的 M block；
4. 正确处理 quasi-Fermi reference field 和 `DDScalingSpec`；
5. 按 4.4 的边界表区分 solved/physical M：Dirichlet solved 行和 insulating pin 行
   清零，thermionic/Schottky physical 行保留储存；
6. 提供 finite-difference storage Jacobian；
7. 输出每个 block 的范数、非零元和非有限值诊断；
8. 建立 MB、Fermi、OldSlotboom on/off、混合材料和量子势 frozen correction 测试；
9. 明确 SRH、迁移率和雪崩不直接产生瞬态储存项，但会通过 J 影响 AC。
10. 所有 FD 在同一冻结 `ACOperatingPointSnapshot` 上执行，不允许调用 Newton
    重新分配 contact-basin reference 或重算尺度。

### 验收标准

- `M` 尺寸与 J 完全一致，模式稳定且无非有限值；
- 约束行的 M 行严格为零；
- 解析 M 与中心差分在显著元素上相对误差不大于 `1e-6`，全局归一误差不大于
  `1e-8`；
- scaling on/off 转回物理单位后矩阵和响应一致，归一差异不大于 `1e-9`；
- quasi-Fermi reference 整体变换不改变物理 M 和最终响应，归一差异不大于
  `1e-9`；
- 电子/空穴储存符号通过 manufactured conservation test，而不是只靠源码审阅；
- thermionic/Schottky fixture 的连续性 M 行保留，essential 和 insulating pin 行为零；
- frozen reference FD M 与解析 M 的归一差异满足同一 `1e-6/1e-8` 门；
- 现有 full test suite 无回归。

P1 和 P2 在 P0.5 接触行接口冻结后可以并行；二者不得各自复制 node volume、边界
分类或 scale 恢复逻辑。

## P2.5：纯介电 Poisson AC 资格门

### 目标

在没有 DD 储存、迁移率、复合或接触导电电流干扰的情况下，端到端验证 Poisson
反力电荷、隐式 Dirichlet 激励、实块系统 M=0 和位移电流提取。

### 实施内容

1. 建立带两个金属电极的均匀氧化层平板 fixture；
2. 求解 DC Poisson 工作点；
3. 由 physical Poisson reaction/J 得到 Q/Qx；
4. 以 M=0 通过 AC block 路径求状态扰动；
5. 输出 `Y=j*omega*C`；
6. 对 stacked/interleaved 两种候选布局采集 fill ratio 和解残差，但暂不据单一小
   fixture 冻结最终布局。

### 验收标准

- `A` 在双精度舍入内为零；
- `C` 与解析平行板值相对误差不大于 `1e-8`；
- `Y_imag=omega*C` 的归一误差不大于 `1e-11`；
- 两端电流等大反向，KCL 归一残差不大于 `1e-10`；
- 从 `1e3` 到 `1e12 Hz`，C 与频率无关、Y_imag 与频率成正比；
- P2.5 未通过不得进入含载流子的 P3。

## P3：单器件隐式 AC 数值核心

### 目标

在给定收敛 DC 工作点上构造和求解实数 2x2 频域块系统，并支持多端 RHS。

### 实施内容

1. 先给 `LinearSolver` 增加对象级 `numericFactorizationCount()` 和精确单元测试；
2. AC solver 持有专用 `LinearSolver`，不得与 3N DC Newton solver 共享实例；
3. 实现 `ACBlockSystemBuilder`；
4. 组装 `J`、`M` 和实块矩阵；
5. 对 stacked/interleaved 采集 fill 和时间，以代表性 MOS 而非平板 fixture 冻结布局；
6. 为 implicit contact excitation 建立解析 `dR/du`；隐式 Dirichlet 下 `Mu=0`，
   单位 RHS 为实且仅在被激励 essential BC 行非零；
7. 以同一冻结 snapshot 上的局部边界有限差分作为 RHS oracle；不得通过正常
   Newton `+/-delta V` 重求解；
8. 每个频率只装配一次块矩阵并数值分解一次，对每个端口顺序求解 RHS；
9. 多 RHS 循环持有同一个 `const SparseMatrixd`，不得 coeffRef、compress、copy-back
   或以任何方式触碰矩阵值；
10. 返回复状态扰动、线性残差、条件/失败诊断和性能计数；
11. 支持单频、显式频率数组、linear 和 decade 频扫；
12. QP 默认 fail-closed；非正频率、未知端口、重复端口、未收敛工作点和奇异系统
    同样 fail-closed。

### 测试

- 人工 1x1 和 2x2 `J+j*omega*M`；
- M=0 的纯电阻极限；
- J=0、适当约束下的纯电容极限；
- 多 RHS 与逐列独立求解一致；
- 稀疏模式相同而 omega 改变时只复用 symbolic analysis，不误复用 numeric factor；
- 同一 omega 多端 RHS 只做一次 numeric factorization；
- 边界解析 RHS 与中心差分 RHS 对照。
- 冻结 per-node QF reference、DD scaling、residual diagnostic scales 和 frozen QP
  field 后的 FD RHS 对照；
- 重掺杂 MOS 在 `1e3--1e12 Hz` 的块矩阵调理扫描；
- `AC + density-gradient QP` 默认配置 fail-closed 并给出明确错误。

### 验收标准

- well-scaled 人工系统的相对线性残差不大于 `1e-11`；
- 实块解与 Eigen dense complex reference 的归一差异不大于 `1e-11`；
- 每频率专用 solver 的 `numericFactorizationCount()` 增量严格等于 1；
- 每频率 RHS 数等于实际激励端口数；
- 改变频率会重新数值分解但不重复 symbolic analysis（稀疏模式不变时）；
- 冻结 snapshot 上解析边界 RHS 与有限差分归一差异不大于 `1e-12`；
- 重掺杂 MOS 全频段相对线性残差不大于 `1e-10`，无未登记 pivot/奇异告警；
- 布局决策报告同时给出 stacked/interleaved 的 fill ratio、factor/solve time 和残差；
- 所有失败都包含 bias、frequency、矩阵维度、零行/列、非有限值和后端信息。

## P4：端口线性化、完整 Y/A/C、配置和输出

### 目标

把复状态响应变成物理端口响应，并接入可重复的 DC 偏压扫描工作流。

### 实施内容

1. 使用 P0.5 `ContactReactionFunctional`，不实现第二套生产 SG 导数 kernel；
2. 从 physical electron/hole continuity residual、J 和 M 接触行提取端口载流子反力
   及其小信号响应；
3. 从 physical Poisson residual/J 接触行提取 Q/Qx，并形成位移项
   `j*omega*dQ`；
4. 按冻结符号组合电子、空穴和位移响应，得到总端口复电流；
5. 继续保留 `ContactCurrent::computeDetailed` 的 SG edge 路径，作为 DC 值的独立
   交叉验证，不作为 AC production derivative；
6. 逐端单位 phasor 激励生成 full Y，再派生 A 和 C；
7. 计算 terminal KCL、Y row/column sum、charge closure 和共同模式检查；
   对 closed fixture，接触 reaction weights 的并集必须覆盖全部物理边界反力，使 KCL
   成为装配/激励 stamping 验证，而不是独立 SG 导数能否偶然闭合的测试；
8. 实现 `ACSmallSignalSolver` 和 `ACSweepRunner`；
9. 在 DCSweep 只增加只读 `onAcceptedPoint` 观察者，并把 P0 冻结的 AC 目标注册为
   forced landing points；自适应步先裁剪到目标，精确接受后才触发 AC；
10. 中间 DC 接受点继续服务延续，但不运行 AC、不参与 P5 评分；
11. 实现配置解析、long-form CSV、summary 和 manifest；
12. 支持保存失败前已完成的 bias/frequency 数据，并明确 partial 状态；
13. 可选输出复 `psi/phin/phip` 场，默认关闭以控制体积。

### 验收标准

- 统一残差行 DC 端电流与既有 `computeFromResidual` 在 ohmic fixture 上逐接触一致；
- 统一残差行端电流与 SG `computeDetailed` 在适用 fixture 上通过现有 DC KCL 门；
- 接触反力 J/M 行提取与冻结 snapshot 中心差分在显著元素上相对误差不大于
  `1e-6`；
- 每个电压列的 terminal KCL 归一残差不大于 `1e-8`；
- 共同模式电压激励的端口响应归一残差不大于 `1e-8`；
- full matrix 与指定参考端生成的 reduced matrix 可相互一致转换；
- CSV 中 `Y_real=A`、`Y_imag=omega*C` 在 17 位有效数字写回后闭合；
- 每米与每微米列换算精确到输出舍入误差；
- bias/frequency/terminal 排序确定，重复运行文件字节级稳定，时间戳字段除外；
- ACCompute 不改变 DC 方程、接受判据、预测器或 retry 语义，但会
  通过步长裁剪**受控改变接受点序列**以强制命中目标；测试必须证明所有目标精确
  落点且观察者只读；
- 第一版 `system=explicit` 必须明确报“尚未支持”，不能静默按 implicit 运行。

## P5：解析、DEVSIM 和 Sentaurus 分层验证

### 目标

证明 Vela 的 AC 不仅能运行，而且在物理、数值和端口语义上可与独立实现比较。

### 验证层次

#### L-A：解析和 manufactured fixtures

- 平行板电容；
- 线性 RC 网络等价块；
- 平衡 MOSCAP 低频/高频极限；
- 低注入 PN 结 small-signal conductance/capacitance；
- avalanche-on 近击穿二极管：雪崩只进入 J，不进入 M；
- 对称无源结构的互易性；
- 偏置晶体管只检查 KCL/gauge，不强制互易。

#### L-B：DEVSIM 开源交叉验证

- 使用公开可提交的 PN 二极管和 MOSCAP；
- 固定网格、材料参数、边界、DC 点和频率；
- 比较 DC 状态、端电流、Y、C 和频率趋势；
- 保存 DEVSIM 版本、脚本和输出，不依赖联网运行 CI。

#### L-C：Sentaurus `AC1_des.cmd`

- 使用同一 TDR、端口、偏压和 1 MHz；
- 使用 P0 冻结的强制目标 bias；若任一工具未在目标上收敛，该点记为失败而不是
  插值；
- 先比较 DC 工作点，再比较 full A/C；
- 按显著元素、全局矩阵、Cgg 曲线和守恒分别评分；
- 不用最终 Cgg 一条曲线替代全矩阵比较；
- 不在陡峭区跨不同 bias 点插值后评分。

### 暂定跨工具门槛

P0 后必须冻结；以下仅作为 v2 审核起点：

- 只在双方共同的精确 bias/frequency 点评分；
- “共同精确点”来自强制落点合同，不是两个自适应接受点集合的偶然交集；
- 对 `|X_ref| >= 1e-4 * max|X_ref|` 的显著矩阵元素，复 Y、A、C 相对误差不大于
  `5%`；
- 全矩阵 `max|X-Xref|/max|Xref|` 不大于 `2%`；
- Cgg 曲线显著区中位相对误差不大于 `3%`、P95 不大于 `5%`；
- 接近零的元素用全局归一绝对误差，不报告失真的巨大相对误差；
- 频率趋势、符号、拐点顺序、KCL 和低/高频极限是独立硬门；
- 任何拟接受差异必须进入 `known_difference_ledger.json`，注明范围、证据和到期条件。

### 验收标准

- L-A 所有解析硬门通过；
- DEVSIM 至少一个二极管和一个 MOSCAP 通过冻结门槛；
- Sentaurus 隐式工程的全部端口矩阵和 Cgg 曲线通过冻结门槛；
- 低频 AC C 与 P1 的准静态全端口 C 在选定低频 fixture 上归一差异不大于 `1%`；
- 提高/降低频率时结果连续，无非物理符号翻转或数值尖峰；
- avalanche-on fixture 的 M 与 avalanche-off 完全一致，而 Re(Y) 按冻结 J 路径
  响应；对应边级 avalanche Jacobian 审计必须通过；
- 所有比较报告同时展示参考值、Vela 值、差异、阈值和 pass/fail，禁止只输出结论。

## P6：显式 MixedMode V/R/C

### 目标

在稳定器件端口线性化之上增加一般电路节点，使 `AC_des.cmd` 成为验收对象。

### 最小范围

- ground；
- 独立 DC/AC 理想电压源；
- 电阻；
- 电容；
- 一个或多个 Vela device instance；
- `Node` 选择和理想源 `Exclude` 语义；
- blocked device/circuit assembly；
- full node Y 和选定端口输出。

第一版不含电感、受控源、二极管 compact model、SPICE parser 或时域电路。

### 实施内容

1. 定义 `CircuitMNA` 节点、支路未知量和 element stamping；
2. 为 V/R/C 建立 DC、J 和 M stamps；
3. 将器件端口电压扰动和端电流响应接入 MNA；
4. 明确理想电压源为什么必须从 AC unknown/equation 提取集合中 Exclude；
5. 实现显式 system 配置和 fail-closed parser；
6. 先验证解析 RC，再验证单器件显式/隐式等价；
7. 最后复现 Sentaurus `AC_des.cmd`。

### 验收标准

- R、C、串并联 RC 和多节点网络与解析 Y 的相对误差不大于 `1e-10`；
- MNA 每个非 ground 节点 KCL 归一残差不大于 `1e-10`；
- Vela 显式系统与 Vela implicit 系统全矩阵归一差异不大于 `1e-9`；
- Sentaurus `AC_des.cmd` 和 `AC1_des.cmd` 在同一冻结门槛内分别通过；
- 错误 Node、浮空电路、重复源、未 Exclude 的理想源给出可诊断错误；
- 多器件实例的端口命名无碰撞，输出包含 instance 和 node 标识。

## P7：性能、文档、CI 和长期扩展

### 实施内容

1. 为 J/M/block build、symbolic analysis、numeric factor、RHS solve、terminal extraction
   增加分阶段计时；
2. 记录矩阵维度、非零元、fill ratio、内存、频率数和端口数；
3. 比较实块与可选 complex backend spike，只有有证据时才决定迁移；
4. 将小型解析/DEVSIM fixture 放入常规 CI；
5. 将 Sentaurus 对照放入人工或夜间回归，不把受许可文件提交仓库；
6. 编写用户配置、输出 schema、端口符号、错误诊断和示例文档；
7. 增加 Python API 仅在 C++ API 稳定后进行；
8. 评估 SDE 高斯植入支持，但不与 AC-L6 绑定；
9. 记录未来 transient、noise、S 参数和热 AC 的扩展接口，不提前实现。

### 性能验收

- 同一频率 N 个端口不允许做 N 次 numeric factorization；
- AC 单频单工作点峰值内存暂定不超过同状态 DC Jacobian 求解峰值的 `3x`；
- 实块矩阵维度为 `2*Ndof`，非零元增长须有报告，禁止无界 dense 转换；
- 1 个频率、4 个端口的 nMOS AC wall time 暂定不超过一次同点 DC Newton 完整求解的
  `2x`；该值在 P3 性能实测后冻结；
- 100 个频率时，运行时间应主要随 numeric factorization 线性增长，端口 RHS 不得
  主导；
- 所有性能指标都在固定 Release 构建、固定主机、固定网格上测量。

### AC-L6 验收

- Debug 和 Release 全测试通过；
- ASan/UBSan 可用平台无新错误；
- 示例配置从干净构建可运行；
- schema、配置错误和部分输出恢复有测试；
- 文档明确支持边界和非目标；
- Sentaurus/DEVSIM 基准报告可由单一命令重建；
- 性能和已知差异账本已签署。

## 7. 建议工作包和提交拆分

每个工作包应保持可独立审查，禁止在一个提交中同时更改端电荷定义、M 符号、
线性后端和跨工具容差。

| WP | 内容 | 建议提交边界 |
| --- | --- | --- |
| WP0a | Sentaurus oracle、schema、目标 bias、读取器、比较报告 | 与 WP0b 并行；不改求解器 |
| WP0b | 公开解析 fixture | 与 VM 复跑并行，可提交 CI |
| WP-R | physical/solved row view、边界分类、接触反力权重 | P1/P2 共同前置 |
| WP1a | Poisson 接触反力 Q 和 Gauss closure | 不含独立生产 D-flux kernel |
| WP1b | Qx、full quasi-static C 和 FD oracle | 不含动态 AC |
| WP2a | 装配器内 `S(x)` 和 finite-difference M | 可与 P1 并行 |
| WP2b | 解析 M、边界行和 snapshot 审计 | 不含频域 solve |
| WP2.5 | 纯介电 Poisson AC | P1 + 最小 AC block；P3 资格门 |
| WP3a | 人工矩阵实块 builder | 不接器件 |
| WP3b | solver factorization 计数、布局 spike、immutable multi-RHS | AC 专用 solver |
| WP3c | frozen snapshot implicit RHS | 给定工作点 API |
| WP4a | continuity/Poisson 反力端口响应 | SG 路径仅作 oracle |
| WP4b | Y/A/C、频扫和输出 | 接入工作流 |
| WP5a | DEVSIM 基准 | 公开 fixture |
| WP5b | Sentaurus AC1 基准 | 受许可工件留在 staging |
| WP6a | V/R/C MNA | 解析电路测试 |
| WP6b | device--circuit coupling | Vela explicit/implicit 对照 |
| WP6c | Sentaurus AC 显式验收 | 最后进入 |
| WP7 | 性能、CI、文档、可选 complex spike | 不改变已冻结物理合同 |

每个提交至少包含：

- 对应测试；
- 变更前失败、变更后通过的证据；
- 新接口的单位和符号注释；
- 不改变的兼容行为说明；
- 相关 manifest/schema 版本变更；
- 若涉及数值门槛，附批准记录。

## 8. 测试矩阵和命令

### 8.1 测试矩阵

| 维度 | 必须覆盖 |
| --- | --- |
| 统计 | Maxwell--Boltzmann、Fermi |
| BGN | off、OldSlotboom |
| 材料 | 纯半导体、半导体/氧化层/金属栅 |
| 端口 | 2 端、3 端、4 端 |
| 频率 | 单频、linear、decade、重复频率拒绝/去重策略 |
| 工作点 | 平衡、低偏置、导通晶体管 |
| 边界行 | Ohmic/metal-gate Dirichlet、thermionic/Schottky、insulating pin、sheet charge |
| 缩放 | on/off、quasi-Fermi reference 变化 |
| QP | off；on 默认 fail-closed；显式 frozen approximation 仅在后续单独测试 |
| 雪崩 | off、avalanche-on 近击穿二极管（只改变 J，不改变 M） |
| 调理 | 重掺杂 MOS，`1e3--1e12 Hz` |
| 线性后端 | SparseLU；可选 UMFPACK/QR 仅作兼容和诊断 |
| 输出 | per-m、per-um、finite depth、full/reduced matrix |
| 失败 | 非正频率、未知端口、奇异矩阵、非有限值、未收敛 DC |

### 8.2 外部审核指定的回归锚

| Fixture | 预期关系 | 所属阶段 |
| --- | --- | --- |
| 纯介电氧化层平板 AC | `Y=j*omega*C_analytic`，`A=0` | P2.5 |
| frozen reference FD RHS vs analytic stamping | 归一差不大于 `1e-12` | P3 |
| continuity reaction DC current vs `computeFromResidual` | 逐接触机器精度一致 | P0.5/P4 |
| 同频多 RHS | `numericFactorizationCount()` 增量严格为 1 | P3 |
| 重掺杂 MOS `1e3--1e12 Hz` | 相对解残差不大于 `1e-10`，无未登记 pivot 告警 | P3 |
| AC + density-gradient QP 默认配置 | fail-closed，明确说明 frozen QP 非自洽 | P3 |
| thermionic/Schottky 和 insulating fixture | 自然通量行保留 M；essential/gauge 行 M=0 | P2 |
| avalanche-on 近击穿二极管 | M 与 avalanche-off 相同；Re(Y) 消费已审计 J | P5 |
| 自适应 DC 强制 AC 落点 | 每个目标精确命中；中间点不触发 AC；无插值 | P0/P4 |

### 8.3 Windows UCRT64 基本命令

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
Set-Location "D:\code-repo\vela-tcad"
cmake --preset windows-ucrt64-debug
cmake --build --preset windows-ucrt64-debug
ctest --preset windows-ucrt64-debug
```

AC 定向测试建议：

```powershell
ctest --test-dir build --output-on-failure -R "terminal_electrode_charge|dynamic_storage|ac_"
```

Release 性能：

```powershell
cmake --preset windows-ucrt64-release
cmake --build --preset windows-ucrt64-release
ctest --preset windows-ucrt64-release -R "ac_"
```

## 9. 数据、可复现性和报告契约

### 9.1 文件分层

```text
reference_staging/sentaurus_ac_t2022_03_sp2/
  original/          # VM 原始/不可再生输入和输出，只读、ignored
  normalized/        # 统一 long-form CSV
  manifests/         # 哈希、版本、合同、阈值

examples/ac_*/       # 可公开、可提交的小型 Vela fixture
tests/fixtures/ac/   # 小型解析/DEVSIM 固化数据
build*/ac_runs/      # 可再生运行缓存
```

不得把 Sentaurus 受许可文件提交仓库。不能再生的 oracle 不应只留在 `build/`。

### 9.2 每次运行必须记录

- 代码提交和 dirty 状态；
- 编译器、构建类型、Eigen/SuiteSparse 后端；
- 配置、网格、材料、参数哈希；
- 物理模型和离散 profile；
- DC 收敛指标；
- 频率、端口顺序、激励和参考端；
- J/M/block 维度、非零元、分解次数；
- Y/A/C、KCL、charge closure、row/column sums；
- wall time 和峰值内存；
- 失败分类和 partial output 状态。

### 9.3 报告必须回答

1. 比较的是哪个工作点、频率、矩阵元素和单位？
2. DC 状态是否先闭合？
3. 差异来自 J、M、端电流、端电荷还是边界激励？
4. 数值是否通过 KCL、Gauss、gauge 和线性残差硬门？
5. 参考和 Vela 是否使用相同网格、模型和参数？
6. 哪些差异已解释，哪些仍是 blocker？
7. 达到哪个 AC-L 等级？

## 10. 风险、阻塞条件和停止规则

### 10.1 主要风险

| 风险 | 影响 | 缓解 |
| --- | --- | --- |
| 接触反力权重或物理 scale 错误 | Q/I 符号和单位全错 | 统一 residual/J/M 行提取 + 解析 fixture + 既有 DC current 对照 |
| solved/physical 行视图混淆 | essential 行掩盖真实端口反力 | P0.5 双 view API 和边界分类测试 |
| 连续性储存符号错误 | 频率响应相位错误 | assembler-owned S/M + manufactured conservation + FD M |
| 缩放/reference 未进入 M | scaling on/off 不一致 | P2 强制 invariance 测试 |
| 用区域体电荷或独立几何替代电极反力 | 多端 C 不守恒 | physical Poisson residual/J 接触行提取 |
| QF reference 或 scale 在 FD 中重分配 | 产生假性 dR/du 不匹配 | immutable AC snapshot；禁止 FD 经过 Newton |
| frozen QP 被误认为自洽 AC | 量子响应声明错误 | 默认 fail-closed；显式 approximation 标记 |
| 每端重复分解 | 多端/频扫性能失控 | profiler 硬门 |
| DC/AC 共享 solver 冲刷 pattern cache | 重复 symbolic/numeric 工作 | AC 专用 `LinearSolver` 实例 |
| 多 RHS 间触碰块矩阵 | bitwise cache 静默失效 | `const` block matrix + 对象级 factorization count |
| 自适应 DC 不命中 AC 目标 | P5 无共同点可评分 | P0 冻结目标，步长裁剪强制落点 |
| 重掺杂行缩放导致 AC 块病态 | 高频/低频 pivot 或大残差 | `1e3--1e12 Hz` MOS 调理扫描，先登记不擅自改缩放 |
| Sentaurus 端口/单位解释错误 | 错误 oracle | P0 全矩阵和 explicit/implicit 自检 |
| 过早实现 MNA | 核心错误难定位 | AC1 先于 AC explicit |
| 为拟合扩大容差 | 无法证明等价 | threshold freeze + difference ledger |
| SDE 前端缺口被误认为 AC blocker | 范围膨胀 | 先使用生成 TDR |

### 10.2 阶段阻塞条件

遇到以下任一情况应停止进入下一阶段：

- Sentaurus explicit/implicit 结果自身不能解释；
- 端口方向或二维单位尚未冻结；
- P1 Gauss closure 或 P2 M finite-difference gate 失败；
- physical/solved row view 或自然边界分类尚未通过；
- scaling/reference invariance 失败；
- 同频率发生多次不必要 numeric factorization；
- KCL 只能通过删掉某个端口或经验缩放满足；
- 跨工具比较依赖不同 DC 工作点或未经批准的插值；
- AC FD oracle 经过 Newton、重新分配 reference 或重新选择 scale；
- 阈值在看过最终结果后被临时放宽；
- 生产代码只能依靠受许可 Sentaurus 文件才能通过常规 CI。

### 10.3 回退原则

- P1 失败：停在端电荷定义，不进入频域求解；
- P2 失败：保留准静态 full C，不宣称 AC；
- P3/P4 失败：保留人工矩阵 solver，不接工作流；
- P5 Sentaurus 未达门槛但解析/DEVSIM 通过：登记差异，状态为 AC-L3，不认领
  Sentaurus 等价；
- P6 失败：保留隐式 AC，显式 MixedMode 继续标记未支持；
- complex backend 无收益：继续使用实块系统，不因“架构美观”迁移。

## 11. 已冻结决策和剩余批准项

### 11.1 v2 已冻结，不得由实现者改写

1. **电极电荷**：生产定义是未做 essential 行替换的 physical Poisson residual/J
   接触行反力；显式 `D dot n` 仅作解析 fixture oracle。
2. **端电流**：生产定义是 physical electron/hole continuity residual/J/M 接触行
   反力；SG `computeDetailed` 仅作独立 DC oracle。
3. **离散一致性**：Q/I/S/M 必须共享 `CoupledDDAssembler` 的 node volume、edge
   coupling、mixed-Voronoi 策略、CarrierStatistics、BGN、QP、scale 和 QF reference。
4. **边界视图**：solved view 执行 essential replacement；physical view 保留自然
   边界并跳过 essential replacement。空 bcs 不等于 physical view。
5. **动态储存**：实现在 `CoupledDDAssembler` 内；独立文件只定义结果类型/接口。
6. **工作点快照**：J/M/dRdu/FD 固定 QF reference fields、scaling、residual diagnostic
   scales 和 frozen QP field；FD 不经过 Newton。
7. **QP**：默认 fail-closed；显式 frozen-QP 只能作为带 manifest 标记的近似，不能
   认领自洽量子 AC。
8. **ACCompute**：目标 bias 是强制落点；自适应步裁剪命中目标，中间接受点不评分，
   禁止插值。
9. **DCSweep 接入**：只读 `onAcceptedPoint` 观察者加 forced landing，不做侵入式
   operating-point 循环重构。
10. **线性求解**：AC 使用专用 `LinearSolver`；P3 前新增
    `numericFactorizationCount()`；同频多 RHS 使用不可变同一矩阵对象。
11. **块系统**：采用实 2x2 数学系统，但 stacked/interleaved 由 P3 fill/time 数据
    冻结。
12. **激励幅值**：内部单位 phasor 仅用于导数归一化，不作为用户大信号物理旋钮。
13. **2D 单位**：规范值使用 per-meter，per-um 和 finite depth 是显式派生量。
14. **矩阵输出**：full matrix 始终保存，reference 只生成 reduced view；manifest
    必须保存端口顺序。
15. **P6 范围**：V/R/C、独立电压源、ground、device instance、Node/Exclude 是第一
    版最小充分 MixedMode。

### 11.2 仍需在对应阶段批准

1. `s_Q`、电子/空穴 reaction 的最终正号和物理 scale 数值；由 P0.5/P1 fixture
   冻结，不允许凭源码直觉决定；
2. legacy `cv_quasistatic` 的弃用周期和新 `poisson_reaction` 方法何时成为默认；
3. P0 后跨工具显著元素 floor 和最终误差门槛；
4. 首个公开 DEVSIM fixture 的网格、材料和参数；
5. stacked 或 interleaved 的最终块布局；
6. P6 的具体配置/网表格式；
7. AC-L6 后是否值得增加 native complex backend；
8. Sentaurus oracle 的存储、访问权限和长期备份责任人。

## 12. 供其他大模型重点审核的问题

请审核者以 v2 已冻结决策为前提，不要重复建议平行 SG/D-flux 生产 kernel。应逐项
给出“接受 / 修改 / 阻塞”及证据。重点问题如下：

1. physical/solved 双 view 是否足以覆盖 essential、thermionic/Schottky、insulating
   pin、metal gate 和 static sheet charge，是否还缺其他行类型？
2. continuity reaction 使用 `J+j*omega*M` 接触行、Poisson reaction 使用
   `j*omega*Qx` 时，是否存在接触控制体储存或位移电流双计数？请给出离散推导。
3. physical Poisson reaction 的符号和 scale 应如何从现有 residual 单位严格推导？
4. `S/M` 对 electron/hole continuity 的符号是否与现有 residual 定义一致？
5. frozen snapshot 还需要封存哪些会影响装配的可变缓存或 branch choice？
6. QF reference field 变化的 gauge invariance 测试是否还应覆盖 contact-basin 边界
   节点重新归属？
7. P1/P2 finite-difference 步长和 `1e-6/1e-8` 门是否足以区分截断与舍入误差？
8. KCL、Gauss closure、row/column sum、互易和 passivity 在哪些 fixture/偏置下才是
   合法硬门？
9. P2.5 纯介电路径是否还需要显式测试多介质界面和 interface sheet charge？
10. stacked/interleaved 的评测器件和 freeze 指标是否足够，是否需要 ordering 交叉项？
11. heavy-doping `1e3--1e12 Hz` 扫描应使用何种无量纲条件指标和 pivot 告警标准？
12. forced landing 与 retry/step growth/预测器交互还需要哪些状态机测试？
13. long-form 输出与 manifest 是否足以重建 full/reduced matrix 和 frozen snapshot？
14. P0 暂定 Sentaurus self-check 与 P5 跨工具门槛是否合理？
15. P1/P2 并行后，哪一个共同接口变更必须触发双方重新评审？
16. 是否存在本方案仍遗漏的守恒律、规范自由度、接触反力或 2D 单位风险？

审核回复建议使用：

```text
总体结论：接受 / 有条件接受 / 阻塞

阻塞项：
- [章节] 问题、证据、建议修订

非阻塞改进：
- [章节] 建议和收益

建议新增测试：
- fixture、预期关系、容差依据

建议调整依赖：
- 原顺序 -> 建议顺序，原因
```

## 13. 实施批准后的首轮任务

批准本方案并完成外部审核后，只并行启动 P0a/P0b，不直接编码 WP-R 或 P1--P6：

1. 在 VM 复跑并封存官方 AC 工程；
2. 导出 explicit/implicit 全矩阵；
3. 写最小解析器和 normalized CSV；
4. 导出并冻结 ACCompute 物理目标 bias 和强制落点语义；
5. 并行建立公开平行板/Poisson 解析 fixture；
6. 冻结端口、单位、符号和记录数；
7. 生成 `threshold_freeze.json` 草案；
8. 由第二个模型或人工独立复核 P0；
9. P0 通过后先单独批准 WP-R；WP-R 通过后才可并行批准 WP1/WP2。

## 14. 新任务启动指令

可将下面内容直接交给后续实现代理：

```text
请阅读：
docs/superpowers/plans/2026-08-31-ac-small-signal-simulation-development-plan.md

当前只执行 P0a/P0b：Sentaurus oracle、强制目标 bias、合同/阈值冻结和公开解析
fixture。不要实现 physical row view、AC solver、动态储存矩阵、MixedMode 或修改
cv_quasistatic 语义。

开始前：
1. 记录 git commit 和 dirty 状态；
2. 复核 Sentaurus VM 路径、版本和输入哈希；
3. 原始工件只读保存到 ignored reference_staging；
4. 导出 explicit/implicit 的完整 A/C/Y，而不只是 Cgg 图；
5. 记录 normalized ACCompute 到物理 bias 的精确映射和实际落点；
6. 公开解析 fixture 与 VM 复跑并行，但不得混淆来源；
7. 所有结论必须带端口、偏压、频率、单位和来源；
8. 先提交 P0 产物和审计报告，等待审批后再进入 WP-R。
```

## 15. 参考资料

- Sentaurus Device Training, MixedMode and Small-Signal AC Simulation：
  <https://ghzphy.github.io/Sentaurus_Training/sd/sd_3.html>
- DEVSIM solver documentation：<https://devsim.net/solver.html>
- DEVSIM small-signal diode example：
  <https://github.com/devsim/devsim/blob/main/examples/diode/ssac_diode.py>
- Genius-TCAD-Open：<https://github.com/cogenda/Genius-TCAD-Open>
- S. E. Laux, “Techniques for Small-Signal Analysis of Semiconductor Devices,”
  IEEE TCAD, 1985, DOI `10.1109/TCAD.1985.1270145`：
  <https://research.ibm.com/publications/techniques-for-small-signal-analysis-of-semiconductor-devices--1>
- SISPAD 1999 frequency-domain small-signal paper：
  <https://in4.iue.tuwien.ac.at/pdfs/sispad1999/00799254.pdf>
