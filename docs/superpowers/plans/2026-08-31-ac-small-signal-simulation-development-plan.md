# Vela TCAD AC 小信号仿真功能开发方案

创建日期：2026-08-31；本次复核：2026-09-10

状态：外部审核修订稿 v4；按 2026-09-10 审核修订。设计修订不等于物理策略已批准或功能已实现

代码复核基线：`d76434089dcd8564069be4b822260d7aa9c72662`。

文档基线：本轮开始时该计划已有未提交的 v3 修改，v4 在其上增量修订；不能再将
工作区称为干净。代码基线与文档修订版本分别标识，文档提交不表示求解器已经改变。

证据边界：本次复核了当前源码、配置/回归文档与公开一手资料；未重新运行 VM、
Sentaurus 或数值回归。下文“验收标准”均为未来必须执行的门，不是已经通过的结果。

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
  -> P0 合同：材料/多晶栅/有效 AC 导数设置、DC 与 A/C 分离阈值
  -> P0.5 = WP-R：材料体积决策、独立装配器快照复现、物理导数和接触接口
  -> [P1 残差行电极电荷与准静态 C
      || P2 装配器内动态储存 M
      || WP3a 人工矩阵与等价表示 builder]
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
3. **P3 比较三种表示再选生产路径。** real stacked、real interleaved 与局部
   complex SparseLU 使用同一工作点/方程/验收门；不改全局 `Types.h`。
4. **储存算子与未来瞬态共享设计，但不先实现完整瞬态积分器。** AC 只要求
   状态储存项及其一致线性化。
5. **端口语义先于曲线拟合。** 先冻结端口顺序、参考端、正方向、二维宽度归一化、
   A/C/Y 定义和输出 schema，再比较 Sentaurus 数值。
6. **所有跨工具数值门槛先标为暂定。** P0 导出官方 oracle 后必须生成并审核
   `threshold_freeze.json`，不得为通过最终曲线临时改门槛。
7. **端口响应统一复用同一装配器。** Essential 接触提取未覆盖的物理行反力；
   Robin 接触提取同一装配器的边界通量 stamp，不能对已满足平衡的完整自然边界
   残差求和充当端电流。独立 SG 电流和显式 `D dot n` 只作交叉验证 oracle。
8. **AC 目标偏压必须强制落点。** 自适应 DC 步长应裁剪到下一个待计算 AC 目标，
   中间接受点可以存在但不替代目标点，也不允许靠插值评分。
9. **AC 线性化快照冻结数值坐标。** QF reference field、DD scaling、残差诊断尺度
   和 frozen quantum correction 在 J/M/dRdu/FD oracle 之间不得重分配或重算。

### 1.1 文档导航与修订范围

| 文档 | 唯一职责 |
| --- | --- |
| 本文 | 背景/目标、P0--P7、工作包、测试和阶段验收 |
| [规范合同](../specs/2026-08-31-ac-small-signal-contract.md) | §4.1--4.9 的数学/数据合同与 §11.1 checklist |
| [审核附录](2026-09-10-ac-small-signal-review-appendix.md) | 审核处置、当前源码锚点、独立工作树来源、审核问题/文献 |
| [预登记差异账本](2026-09-10-ac-known-difference-ledger.json) | 体积/DC/多晶栅/参考导数的 open 项，不是已执行的数值结果 |

本次改变：WP-R 前移材料支持决策和快照残差复现；P3 三路后端比较；P4 正交 sweep.ac
及路由边界；P5 逐目标 DC 门和独立 A/C 门。原规范条号保留用于外部审核引用，
本文提及“合同 §4.x/§11.1”均指规范合同，不是本文件中的重复定义。
审核附录区分历史决定与 v4；当前未决项见 §11.2。


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

i_hat = (H + j*omega*K) * x_hat + (D + j*omega*E) * u_hat
```

其中：

- `Jx = dR/dx`，Vela 已有主要实现；
- `Mx = dS/dx`，当前缺失；
- `H/D` 是端口稳态导电响应的状态/电压导数，来自同一装配器；
- `K/E` 包含接触连续性储存反力与 Poisson 反力导数，不能一概写成电极 `Qx/Qu`；
- 当前 `TerminalCharge` 不能直接替代 Poisson 反力；具体公式见 4.3；
- 对每个被激励端口求解一次，可得到 Y 的一列；
- 同一频率的所有端口列应复用同一次数值分解。

准静态电极 `dQ/dV` 只描述电荷随 DC 状态的变化。导通器件的低频总电流还包含
载流子响应滞后引起的导电项；即使 `omega -> 0`，它也不必等于 `Im(Y)/omega`。
P1 的 `C_Q` 是位移电荷资格门；只有适用的阻挡端口/纯介电 fixture 可直接作低频
电容 oracle。一般器件使用 4.8 的完整低频展开。

### 2.3 当前能力与实施缺口

代码仍是 DC 基线，尚无 AC 执行路径。可复用 residual/J、端电流、参考坐标、现有
实数求解缓存和 DC 步进器；需开发统一接触接口、S/M、局部频域求解与工作流。

高风险前置缺口：共享 Si/oxide 节点的材料支持体积、完整物理导数、装配模式缓存/
快照复现、多晶栅建模和参考有效 AC 导数。源码表及 SimpleMOS/LDMOS 方法来源见
[附录 §2.3](2026-09-10-ac-small-signal-review-appendix.md)。独立工作树候选不是当前主线已实现功能。


### 2.4 开源实现和文献基线

- Sentaurus Training 定义 `Y=A+j*omega*C`，并说明显式/隐式 AC 系统、Node、
  Exclude 和 ACCompute；
- DEVSIM 支持 DC 工作点上的 small-signal AC；其二极管示例从电压源支路电流虚部
  提取电容，带负号是源电流与器件电流方向相反，不能直接照搬到 Vela 端口；
- Genius-TCAD-Open 提供 DDMAC、ACSWEEP 和器件--电路混合模式；
- Laux 1985 比较瞬态 FFT、增量电荷分区和正弦稳态方法；
- Lin 等 SISPAD 1999 讨论 Boltzmann--Poisson 小信号频域求解，是“DC 点线性化 +
  `j*omega`”方法参考，不是 Vela DD 储存符号或 THz 物理有效性的直接证明。

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
| AC-L0 | P0：参考工件、材料/多晶栅/有效导数合同、目标点及阈值分区冻结 |
| AC-L1 | P0.5=WP-R + P1：材料策略、快照复现/物理导数资格、电极 Q/C_Q 通过 |
| AC-L2 | P2 + P2.5：S/M 审计及纯介电端到端 AC 通过；P1 为共同前置 |
| AC-L3 | P3 数值核心和 P4 工作流均通过：隐式 AC 单频/频扫、完整 Y/A/C |
| AC-L4 | P5：全目标 DC 前置门、DEVSIM/原始 AC1 分合同 A/C 门、适用低频极限通过 |
| AC-L5 | 显式 MixedMode V/R/C 系统通过，复现 `AC_des.cmd` |
| AC-L6 | 性能、文档、CI、错误诊断和长期回归全部固化 |

报告和发布说明必须声明达到的最高等级，不得仅写“已支持 AC”。

## 4. 数学与数据合同入口

[规范合同 §4.1--4.9](../specs/2026-08-31-ac-small-signal-contract.md) 是唯一公式/单位/误差定义来源：
端口与二维单位、接触反力及储存抵消、材料支持、S/M、三种等价表示、冻结快照、
输出 schema、C_Q/C_AC 低频区别和可计算误差指标。

统一采用其中 §11.1 gate_id checklist；实现者不在不同阶段另写一套符号/单位或例外。


## 5. 目标软件架构

### 5.1 分层

```text
物理/离散层
  CoupledDDAssembler physical/solved row views
  CoupledDDAssembler storage S/M
  ContactReactionFunctional (essential row weights + natural boundary stamps)

频域数值层
  ACFrequencySystemBuilder (real block / local complex candidates)
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
| `include/vela/post/ContactReactionFunctional.h` | 按方程的接触行权重、自然边界 stamp、单位和反力接口 |
| `src/post/ContactReactionFunctional.cpp` | Essential 提取 residual/J/M；Robin 提取装配器 boundary stamp |
| `include/vela/simulation/ACOperatingPointSnapshot.h` | WP-R 冻结快照和导数资格元数据；名字为拟新增 |
| `include/vela/solver/ACBlockSystem.h` | 实块等价表示与 oracle，不预定为生产后端 |
| `src/solver/ACBlockSystem.cpp` | 实块装配和质量检查 |
| `include/vela/solver/ComplexLinearSolver.h` / `src/solver/ComplexLinearSolver.cpp` | P3 局部 complex SparseLU 候选及同等计数合同 |
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
  不侵入式重构当前约 8680 行的求解流程；复用现有 `runDCSweepStepControl`，
  并按需扩展 `include/vela/simulation/DCSweepStepControl.h` 的目标接口；
- `include/vela/solver/LinearSolver.h`、`src/solver/LinearSolver.cpp`：P3 前必须增加对象级
  `numericFactorizationCount()`；如需显式 `factorize/solveFactored` 或矩阵 RHS API，
  另以测试和性能证据决定；
- `include/vela/equation/CoupledDDAssembler.h`、`src/equation/CoupledDDAssembler.cpp`：
  增加 physical/solved row view、储存 S/M 和边界行分类；
- `include/vela/post/ContactCurrent.h`、`src/post/ContactCurrent.cpp`：使既有
  `computeFromResidual` 委托统一 `ContactReactionFunctional`，SG edge 实现继续作为
  独立 DC oracle；
- `include/vela/simulation/CurveSweep.h`、`src/simulation/CurveSweep.cpp`：保持 DC mode
  枚举，新增正交 `sweep.ac` 配置和显式支持组合校验；不新增 AC sweep mode；
- `src/tools/vela_example_runner.cpp`：接入 AC 工作流；
- `reference_tcad/ac_small_signal/`：加入可公开的 MOSCAP/PN 基准；
- `configs/schema/vela-simulation.schema.json`、`docs/config_schema.md`：实现阶段同步
  JSON schema、运行时解析和文档，不只给 runner 增加私有字段。

不建议把所有新功能继续堆入已经很大的 `DCSweep.cpp`。

### 5.4 建议配置草案：AC 是正交分析块

以下为拟新增字段的局部示意，不是当前可运行 deck；mesh/material 和前置偏置路径省略。
接触名必须解析到实际网格名称，`substrate` 仅是官方器件映射后的示例，不能自动改成 body。

```json
{
  "simulation_type": "dc_sweep",
  "solver": { "method": "newton" },
  "sweep": {
    "mode": "iv",
    "contact": "gate",
    "start": -3.0,
    "stop": 3.0,
    "step": 0.15,
    "initial_step": 0.02,
    "ac": {
      "enabled": true,
      "system": "implicit",
      "nodes": ["source", "drain", "gate", "substrate"],
      "reference_node": "substrate",
      "frequencies_Hz": [1000000.0],
      "compute": { "mode": "normalized_intervals", "intervals": 40, "landing": "force" },
      "output": {
        "matrix_csv": "ac_matrix.csv",
        "summary_csv": "ac_bias_summary.csv",
        "manifest_json": "ac_manifest.json"
      }
    }
  }
}
```

- `sweep.mode` 保持 `iv` 或 `bv_reverse`；AC 作为 `sweep.ac` 附加，不复制两种
  DC 模式的雪崩/停止/延续语义。单工作点先用 P3 snapshot API；P4 可用已支持的
  显式 bias_points 单点路径，但不能绕过 DC 接受和快照门。
- P4 初版支持电压驱动 Newton、普通规则目标循环及显式 bias_points 外层接受点；
  排除 arclength、外接电阻/电流控制、Gummel、legacy CV 混用。后续逐条资格化，
  schema 与运行时必须一致拒绝不支持组合。受限 bv_reverse 仍可覆盖近击穿前
  稳定支路的 avalanche-on fixture，不代表支持跨 snapback 的 AC。
- 端口名唯一且覆盖全部电极，reference 在 nodes 中；只对已资格化的接触/材料/
  导数配置开放。reduced 为 full 的派生视图，不因选择 reference 丢掉端口。
- 显式 frequencies_Hz 与范围对象二选一，有限且正，重复报错。范围 points 指总点数，
  start=stop 仅允许 points=1；linear/decade 仅定义间隔，不隐式改变点数。
- intervals=N 在 Vela 中是含两端 N+1 个目标；在扫描前用索引生成物理列表。
  `ac.compute.bias_points_V` 是 AC 选点请求，与 DC 已有 `sweep.bias_points` 区分；
  显式 DC 列表模式下 AC 目标必须来自该列表，否则报错，不静默插点重排 DC 路径。
- 普通循环合并原 DC 规则目标与 AC 目标，保留方向和 target_id；步长裁剪只改变
  接受点序列，不改变方程或收敛标准。频率可升序调度，bias 不统一升序重排。
- P3 只验收给定工作点 API；JSON/输出/调度在 P4 交付。QP 与未支持模型 fail-closed；
  单位 phasor 是导数归一化，有限差分步长只在 diagnostics 中。


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
12. 冻结 Vela 的目标落点状态机合同及测试输入：未来由自适应步裁剪命中目标、接受后
    才触发 AC。P0 不实现 DCSweep 回调/步长裁剪，其端到端验收属于 P4。
13. 从 `nMOS_dvs.cmd`、预处理 deck 和 TDR 核对 gate 接触所属区域/材料、栅掺杂、
    厚度、是否求电子/空穴以及 WorkFunction。交付 `gate_material_contract.json`：
    明确 retained_semiconductor_poly 或 metal_equivalent，不凭 quantum 模块的 poly
    名称识别推断 DD 栅支持。金属替代须另列模型差异并禁止认领原始多晶栅等价。
14. 封存原始及有效 Math 段、全局/Device scope、Derivatives/AvalDerivatives 等设置、
    版本默认、AC 专用设置和有效模型导数；交付 `reference_derivative_contract.json`。
    不能假定 Newton 开关完全决定 AC/IFM 的导数默认。对缺少有效设置证据的分支标
    unknown；若需受控开关 A/B 则单独冻结输入，保留原始 oracle，不覆盖原 deck。
15. 在 P0 预登记多材料体积、多晶栅替代、DC 状态未对齐及有效导数差异账本；数值
    影响未量化填 pending，不虚填预计百分比作为豁免。

### 必须交付

```text
reference_manifest.json
input_hashes.json
sentaurus_ac_full_matrix.csv
sentaurus_ac1_full_matrix.csv
sentaurus_explicit_implicit_compare.json
ac_contract.schema.json
threshold_freeze.json
gate_material_contract.json
reference_derivative_contract.json
known_difference_ledger.json
```

### 验收标准

- 同一官方工程可从干净工作目录重复运行；
- 所有输入和输出都有哈希，原始文件只读保留；
- 每条矩阵记录能唯一定位 bias、frequency、row、column；
- 显式/隐式 Sentaurus 的全矩阵全局归一差异暂定不大于 `1e-6`；若官方自身高于
  此值，必须先解释并按实测冻结，不能扩大 Vela 容差掩盖；
- `A + j*omega*C` 重建值与导出的复 Y（若可导出）在输出精度内闭合；
- schema 对缺端口、重复端口、非正频率、混合单位和缺哈希 fail-closed；
- 交付初始步约 `0.03 V`、目标不在自然步长网格上的调度测试规格和目标清单；
  本阶段只验证清单/索引/容差，实际裁剪和观察者测试在 P4 执行；
- P0 未通过不得进入最终跨工具等价声明。

阈值文件分区管理：解析/结构硬门、Sentaurus 自身一致性门、跨工具物理差异门、
性能目标分别记录来源。P0 可依据参考自身的重跑噪声和输出精度冻结 floor，但不能
用尚未出现的最终 Vela 结果调门；未知值保持 pending 并阻止相关等价声明。

## P0.5（即 WP-R）：材料策略、快照资格和接触接口

### 目标

先冻结材料支持体积和 B1/B2 的共同接口，使 P1 电荷和 P2 储存可以安全并行。
本阶段明确登记独立工作树候选与当前主线的差异，不把 DC-only 材料局部修正直接
当成 AC 实现。材料策略未决时可以做接口/解析诊断，不能冻结生产 P1/P2 数值门。

### 实施内容

1. 先交付材料体积决策，按合同 4.3.4 审计 Poisson/S/M/连续性源项及固定电荷支持；
   如选择改动，另批准物理子工作包并重收敛 DC，保持 legacy 默认和专项回归边界。
2. 为 `CoupledDDAssembler` 明确区分：
   - `solved` view：包含 essential 行替换，供 Newton/AC 内部状态求解；
   - `physical` view：保留体项、边通量、自然边界和静态 sheet charge，跳过
     essential 行替换，供 essential 端口反力提取；
   - `boundary stamps`：同一装配过程按外部接触分组的自然通量 residual/J/du；
   - 内部 insulating pin 只保留在 solved 代数行，physical 的无载流子贡献行为零；
3. physical view 接受完整 bcs，跳过 essential replacement 和内部 gauge pin，保留
   真实体项与自然边界；不能用 `CoupledDDBoundaryConditions{}` 代替；
4. 定义 `ContactReactionFunctional` 的 Poisson、electron continuity、hole continuity
   权重和物理单位恢复因子；
5. 增加接收完整 bcs 和 `ResidualView::Physical` 的
   `ContactCurrent::computeFromResidual` overload 并委托统一接口；保留旧 overload
   供无自然边界的兼容调用，但不得在 AC 中使用；
6. 暴露 `value`、`stateDerivative` 和后续 `storageDerivative`，不暴露第二套 SG
   production derivative；
7. 为 ohmic、metal gate、thermionic/Schottky、insulating pin 和 interface sheet
   charge 建立行分类测试。
8. 在 P1/P2 之前定义不可变 snapshot、物理单位恢复和状态/方程坐标变换；不能把
   两阶段都依赖的类型留到 WP3c 才设计；
9. 逐模型/配置登记 `physical_derivative`、`newton_only_approximation` 或
   `unsupported_ac`，审计对角 floor、冻结 mobility、源项和非局部分支；
10. 现有模式构建提前屏蔽 essential 行，physical J 需要独立模式或联合模式，
    不能仅关闭最后一次行覆盖；测试保留接触行 off-diagonal 与完整导数支撑集。
11. 建立独立 AC 装配工厂：参考 makeArclengthAssembler 配置路径，floor 关闭构造，
    精确恢复终态 state/reference/bcs；先过 solved residual bitwise replay 再做任何 AC。
12. physical J 一次装配、派生 solved J 或固定两视图缓存；按 snapshot 隔离模式计数，
    不在频率/端口循环切换 signature。完整处理缩放前后两处 insulating pin。

### 验收标准

- `snapshot_residual_replay`：同构建、同 snapshot 的 AC solved residual 与捕获的
  DC 原始终态 residual bitwise 一致；失败前不允许计算 Y；
- `ac.jacobian_pattern_builds<=2`/snapshot；重复频率/RHS 增量为 0，DC 计数不混入；
- 不对称 Si/oxide MOSCAP 先交付面积/电荷支持审计及策略决策，选中 DC 路径可复现；
  本阶段用 manufactured 支持测度 oracle 验证粒子数/电荷恒等式。生产 C_Q 敏感度
  在 P1 验收、生产 S/M 在 P2、AC Cgg/Cgs/Cgd 在 P5，避免 WP-R 反向依赖 P1/P2；
- 现有 ohmic fixture 中，统一反力端电流与 `computeFromResidual` 逐接触一致到机器
  舍入；
- 统一反力端电流与独立 SG `computeDetailed` 在适用 DC fixture 上通过现有 KCL
  门；
- natural thermionic 项在 physical view 中保留，essential replacement 不存在；
- `bulk+boundary=physical` 逐行重建一致，Robin 收敛全残差为零但 boundary stamp
  允许非零端电流；不会把“零残差”误验成“零电流”；
- 物理 J directional FD 通过 4.9；构造会触发 Newton floor 的 fixture，证明 floor
  不进入 AC J，且去除迭代正则后重新检查 DC 物理残差；
- 冻结坐标不冻结物理 mobility；live residual 方向导数必须与申报的 AC J 匹配；
- 原有 Newton solved view 的默认行为/兼容测试不变，新 AC 物理导数模式单独验收；
- 接触权重、block offset、行类型和物理 scale 都有显式 API，不由调用者手写 `N`
  和 `2*N`；
- P0.5 未冻结前，P1/P2 不进入生产实现。

## P1：物理电极电荷与全端口准静态 C

### 目标

建立不依赖区域人为分区的端电荷定义，为位移响应提供 oracle；低频使用范围见 4.8。

### 实施内容

1. 通过 P0.5 的 Poisson contact reaction functional 计算 Q，不另写生产边界积分；
2. 从 physical Poisson residual/J 行提取 Q 和 Qx，并恢复物理单位；
3. 支持半导体 ohmic、金属栅/绝缘层边界和二维宽度归一化；
4. 输出电极反力电荷、由同一 physical Poisson 装配源项求和得到的域内净电荷和
   Gauss closure；
5. 使用冻结 snapshot 的中心有限差分作为 Qx oracle；
6. 独立规则网格 `D dot n` 只作为平行板 oracle；
7. 对 `J*x0=-B*e_l` 求 DC 状态灵敏度，生成 `C_Q=Qx*x0+Qu*e_l`，并用受控 DC
   重求解差分检查总导数；这种 DC 灵敏度测试与固定 x 的 Qx/BC FD 必须分别命名；
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
- `C_Q` 共同模式行和、含 `dQ_domain/dV` 的 Gauss 微分恒等式归一残差不大于
  `1e-8`；仅纯介电等适用 fixture 才要求 `C_Q` 列和为零；
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
8. 建立 MB、Fermi、OldSlotboom on/off 和混合材料测试；frozen QP storage 的单元
   诊断如保留须标记非生产 AC，不能据此开放 QP AC 配置；
9. 明确 SRH、迁移率和雪崩不直接产生瞬态储存项，但会通过 J 影响 AC。
10. 所有 FD 在同一冻结 `ACOperatingPointSnapshot` 上执行，不允许调用 Newton
    重新分配 contact-basin reference 或重算尺度。

### 验收标准

- `M` 尺寸与 J 完全一致，模式稳定且无非有限值；
- 约束行的 M 行严格为零；
- 解析 M 与中心差分在显著元素上相对误差不大于 `1e-6`，全局归一误差不大于
  `1e-8`；
- scaling on/off 按 4.6 同时恢复行和列坐标后算子作用一致，归一差异不大于 `1e-9`；
  P2 不依赖尚未实现的 AC solve，最终响应不变性在 P3/P4 再验；
- quasi-Fermi reference 整体和逐节点重分配后，重新打包同一物理状态的 M 作用不变，
  归一差异不大于 `1e-9`；最终 AC 响应在 P3/P4 再验；
- 电子/空穴储存符号通过 manufactured conservation test，而不是只靠源码审阅；
- 在任意光滑 manufactured 状态/方向上验证
  `-q*ds_n+q*ds_p=d rho_mobile` 及 4.3 的接触抵消，覆盖非零接触体积和单位版本；
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
4. 默认沿用完整 3N solved 系统：纯介电网格的 phin/phip 已由单位 gauge 行固定，
   gauge M=0 即可；无需为本阶段先开发 Poisson-only N 自由度路径。检查至少一个
   电势参考/Dirichlet、介电图连通和 gauge 完整，不能把全 Neumann Poisson 称为非奇异；
5. 输出 `Y=j*omega*C`；
6. 实块最小路径先完成端到端资格；P3 再加入 complex 三路对比，不据这个小 fixture
   冻结最终后端，保留它作为选中生产路径的回归。

### 验收标准

- `A` 在双精度舍入内为零；
- `C` 与解析平行板值相对误差不大于 `1e-8`；
- `Y_imag=omega*C` 的归一误差不大于 `1e-11`；
- 两端电流等大反向，KCL 归一残差不大于 `1e-10`；
- 从 `1e3` 到 `1e12 Hz`，C 与频率无关、Y_imag 与频率成正比；
- 增加双介质串联平板和固定 sheet charge fixture；静态 sheet 改变 Q 的偏置项，
  在固定线性介电常数下不产生额外 M 或 AC 电容；
- P2.5 未通过不得进入含载流子的 P3。

## P3：单器件隐式 AC 核心与三路后端决策

### 目标和前置

在 WP-R 已复现的 DC snapshot 上求解 `J+j*omega*M`。P1/P2/P2.5 通过后进入。
本阶段交付 C++ 单点 API、frequency/RHS 调度和 backend_decision，不提前开放 JSON。

### 实施内容

1. 给现有 LinearSolver 增加对象级数值分解计数，并在局部 ComplexLinearSolver 候选
   提供同一语义：attempt 含失败，success/failure/cache_hit 独立；清缓存不归零
   生命周期计数。修正现有头文件“每次 solve 都分解”的过期注释。
2. AC 装配器和求解器均为独立实例；三种表示消费完全相同的物理 J/M/B 和快照。
   real stacked、real interleaved、complex Ndof 均比较，不扩散全局复数类型。
3. 使用合同 4.5 的 scaled BC RHS：单位电压对应 `+1/V0`，不一律置 1；
   frozen-reference 仿射增量 FD 和一般物理工作点 FD 分别验收。
4. 按合同 4.5 固定三路 benchmark：均匀介电、PN、代表性重掺杂 MOS，以及小型 dense
   complex oracle。固定 Release/主机/ordering，测 factor bytes/峰值内存和 time，
   不用不同标量大小的 fill ratio 单独选型。所有正确性门通过后才按成本选择生产路径。
5. 每 snapshot 批量装配 J/M/端口算子，每频率一次建立数值矩阵，多 RHS 使用同一
   immutable compressed 矩阵/分解。新 bias 或物理/缩放/布局变更使相应缓存失效。
6. 三种路径都报告 count、backend、ordering、矩阵哈希、backward error；非受支持的
   环境变量诊断后端不能绕过计数门后自称支持 AC。

### 验收标准

- 首先满足 `snapshot_residual_replay` 和 `ac.jacobian_pattern_builds<=2`。
- well-scaled 人工系统及三种表示在原复方程上的 backward error `<=1e-11`；
  与 dense complex oracle 的全局误差 `<=1e-11`。
- 新建求解器首次新矩阵所有 RHS 合计 1 次分解；后续 RHS 0 次；改变矩阵值 1 次，
  bitwise 相同缓存命中 0 次。M=0 可跨频复用；失败/显式重试单列而非隐藏。
- 模式不变时 symbolic count 不增加；同矩阵仅改 RHS 不分解，改一个值必须重分解。
  实/复选中后端均接受上述精确计数测试，未选路径不作生产性能承诺。
- 仿射 BC 增量 FD 与 stamping 归一差 `<=1e-12`，覆盖 V0!=1、非零逐节点 reference；
  普通 DC 基线 FD 按合同 4.9 的 `1e-6/1e-8` 与步长扫描门。
- 重掺杂 MOS `f in [1e3,1e12] Hz` 的原复方程 backward error `<=1e-10`；
  记录 componentwise error、可得条件/pivot 诊断和后端实际字节数。三路测试使用相同
  等化策略，不通过新增假导电/储存项改善调理；缺失诊断写 unavailable。
- 纯介电 3N solved carrier gauge 为单位 J、零 M；选中后端重跑 P2.5。
- scaling/reference 变换后物理响应不变；QP、未收敛 DC、奇异系统、未知端口、
  非有限/非正频率和未审计模型给出具体错误。
- 交付 `backend_decision.json`，列出三路正确性/性能数据及选择理由；
  高频扫描只作数值压力测试，不是 DD/电静态近似在 THz 的物理有效性声明。


## P4：端口线性化、完整 Y/A/C、配置和输出

### 目标

把复状态响应变成物理端口响应，并接入可重复的 DC 偏压扫描工作流。

### 实施内容

1. 使用 P0.5 `ContactReactionFunctional`，不实现第二套生产 SG 导数 kernel；
2. Essential 载流子接触从 physical residual/J/M 提取反力；自然接触必须按 4.3.2
   使用 boundary stamp，尚未资格化则拒绝该配置；
3. 从 physical Poisson residual/J 接触行提取 Q/Qx，并形成位移项
   `j*omega*dQ`；
4. 按冻结符号组合电子、空穴和位移响应，得到总端口复电流；
5. 继续保留 `ContactCurrent::computeDetailed` 的 SG edge 路径，作为 DC 值的独立
   交叉验证，不作为 AC production derivative；
6. 逐端单位 phasor 激励生成 full Y，再派生 A 和 C；
7. 计算 terminal KCL、Y row/column sum、charge closure 和共同模式检查；
   对 closed fixture，端口提取必须覆盖全部边界反力/自然通量。离散守恒提供结构
   保证，但仍需检查内部解残差、完整物理导数、源项配对和激励 stamp，不能仅凭
   “同一装配器”宣称 KCL 已自动通过；
8. 实现 `ACSmallSignalSolver` 和 `ACSweepRunner`；
9. 在 DCSweep 只增加只读 `onAcceptedPoint` 观察者，并把 P0 冻结的 AC 目标注册为
   forced landing points；自适应步先裁剪到目标，精确接受后才触发 AC；
10. 中间 DC 接受点继续服务延续，但不运行 AC、不参与 P5 评分；
11. 实现配置解析、long-form CSV、summary 和 manifest；
12. 支持保存失败前已完成的 bias/frequency 数据，并明确 partial 状态；
13. 可选输出复 `psi/phin/phip` 场，默认关闭以控制体积。

当前步进器已对规则 `nominalTarget` 和终点做裁剪；AC 扩展范围按下表固定。P4 优先扩展该接口以合并原 DC 目标与不可变 AC 目标（保留两类 ID），
复用 `DCSweepStepControlState` 的 adaptiveStep、retry 和 stopRequested；recorder
还会接收失败事件，不能直接把它当成 AC 的 onAcceptedPoint。AC 快照应在接受状态
更新完成后只读发布；初始 DC 点在循环外接受的路径也必须发布一次。
现有内部终点容差为 `1e-12`，并以 `nominalTarget += step` 浮点累加；AC 目标应由
索引生成且独立判定完成，不改变 AC 关闭时的 legacy 语义。

| 当前调用点（HEAD 行号，仅定位） | 路径 | P4 初版 AC 处理 |
| --- | --- | --- |
| DCSweep.cpp:7300 | 外边界控制的内部电压延续 | 不挂 AC；external resistor/current-control 组合明确拒绝 |
| DCSweep.cpp:7731 | 耦合外电阻的 outer voltage 延续 | 不挂 AC；P6/后续资格化再支持 |
| DCSweep.cpp:8038 | 显式 bias_points 间的内部延续 | 内部点不触发 AC；在外层已接受请求点发布快照 |
| DCSweep.cpp:8430 | 普通 iv/bv_reverse 电压扫描 | 合并 DC/AC 目标，在最终接受后触发，覆盖循环外初始点 |

arclength 是另一条路径，不由上述四处调用自动覆盖；初版拒绝该组合。测试必须记录
`dc_execution_path`、target_id、observer_count，证明实际进入哪条路径，不能只测试
helper 而遗漏 runner 分派。已有 breakdown/stopRequested 提前停止时剩余 AC 目标
记录 incomplete，未执行目标不进入评分。

### 验收标准

- 统一残差行 DC 端电流与既有 `computeFromResidual` 在 ohmic fixture 上逐接触一致；
- 统一残差行端电流与 SG `computeDetailed` 在适用 fixture 上通过现有 DC KCL 门；
- 接触反力 J/M 行提取与冻结 snapshot 中心差分在显著元素上相对误差不大于
  `1e-6`；
- 每个电压列的 terminal KCL 归一残差不大于 `1e-8`；
- 共同模式电压激励的端口响应归一残差不大于 `1e-8`；
- full/reduced 转换按 4.8 条件验证；保存原始全矩阵，禁止用重建算法掩盖 KCL 误差；
- CSV 中 `Y_real=A`、`Y_imag=omega*C` 在 17 位有效数字写回后闭合；
- 每米与每微米列换算精确到输出舍入误差；
- bias/frequency/terminal 排序确定，固定平台/后端重复运行数值表字节级稳定；
  时间戳、耗时和内存等非确定性遥测分离，不要求跨平台 bitwise 相等；
- ACCompute 不改变 DC 方程、接受判据、预测器或 retry 语义，但会
  通过步长裁剪**受控改变接受点序列**以强制命中目标；测试必须证明所有目标精确
  落点且观察者只读；
- Forced landing 专项：覆盖升/降扫、起终点、相邻目标小于通常 min_step、目标处
  Newton 失败后退步重试、重复请求、恢复运行和多段 bias ramp；目标失败不得标为
  已完成/已评分。剩余距离过小时允许专门的终点裁剪，但不放宽 DC 收敛判据；
- 落点以不可变 target_id 和目标 BC 赋值为准，不能把邻近状态重命名成目标。
  `landing_tolerance_V` 建议 `32*epsilon*max(1 V,abs(Vtarget))`，按序列化精度在 P0
  冻结；完整 DC 电压向量也必须满足合同。观察者只读、拒绝点不触发、每目标恰好一次；
- 观察者失败默认令 AC 分析 fail-closed，保留已接受 DC 状态和 partial 产物；继续
  DC/跳过 AC 需显式策略并标记失败，不能将 AC 错误伪装成 DC 重试来改变结果；
- 相同强制目标清单下，观察者开/关的 DC 接受状态一致；AC 关闭且无目标时保持 legacy
  自适应流程。比较时不能要求“有强制目标”和“无强制目标”的 DC 路径完全相同；
- 第一版 `system=explicit` 必须明确报“尚未支持”，不能静默按 implicit 运行。

## P5：分层验证，DC 对齐先于 A/C 评分

### 目标和复用依据

解析/DEVSIM/原始 Sentaurus AC1 是三类不同证据。先通过每个强制落点的 DC 前置门，
再分别评分 A 与 C。复用 SimpleMOS 的不可变 TDR、精确偏压网格、独立偏置状态链、
模型消融、状态身份/残差回放和误差账本方法，不新建一套含义不同的比较器。
具体来源/commit 及适用范围见审核附录；SimpleMOS 不同于本 AC 器件，其阈值与
“已通过”状态不能直接继承。0.1 dex 的历史电流误差不是 AC 导数误差的数学下界。

### L-A / L-B：解析和开源资格

- 纯介电平板、双介质/sheet charge、线性 RC 等价系统。
- MOSCAP 低/高频极限明确少数载流子供给/G-R/工作点条件；PN 和导通 MOS 使用合同
  4.8 的完整低频展开，不强制所有端口 `C_Q=C_AC`。
- avalanche-on 近击穿稳定电压支路：同快照 M 不受源开关直接影响；重收敛的不同
  状态 M 可以不同。通过 live residual Jv、KCL 和频率趋势，不能只看 Re(Y) 改变。
- 无源平衡 fixture 才要求互易/无源；导通或近击穿器件不强制电容元素非负。
- DEVSIM 至少一个 PN 和一个 MOSCAP，固定代码版本、脚本、网格/材料/参数/
  单位/端口。公开固化 fixture 用于离线 CI；带许可运行独立管理。

### L-C0：每个强制落点的 DC 前置门

先检查输入映射/物理合同和本工具收敛，再按相同物理电压、共同电势参考和节点/
区域映射对比，不插值、不经验平移 Vth/psi。交付 `dc_alignment_gate.json`：

| 指标 | 定义及要求 |
| --- | --- |
| 端电流 | 有信号时分别检查符号和 `abs(log10(abs(Iv)/abs(Is)))`；近零使用预冻结 A/m floor 与绝对误差，不能把符号吞进 log |
| 电势 | 同物理参考下 psi 的 max/P95/RMS（V）；不通过拟合常数偏移消除差异 |
| 载流子 | 各已匹配半导体区域的 log10(n)、log10(p) 差（dex），冻结密度 floor、活动掩码、覆盖率 |
| 状态身份 | 记录原始导出 psi/QF/n/p 与 Vela 重建密度的区别，审计 ni/Vt/参考电位，不混同物理差异与坐标差异 |
| 残差/路径 | DC 物理残差合格、snapshot replay 通过、target 命中、状态链与模型 fingerprint 一致 |
| DC 灵敏度诊断 | 在预定相邻点/独立小扰动上检查 gm/gds 和导数信号/数值噪声；DC 电流吻合不单独证明 A 合格 |

每项容差以单位/区域/端口写入 `threshold_freeze.json.dc`，P0 建合同、P5 评分前批准；
未知容差保持 pending，不把 SimpleMOS 的宽松历史曲线门直接当作本器件 DC 门。

若任一强制点 DC 门失败，该点 A/C 为 `blocked_dc_alignment`，不得标为通过或纳入
正式误差统计；仍可计算并保存明确标记的诊断 AC。报告所有目标、DC 通过数、AC
可评分数及缺失原因；不能只挑通过点集合认领 AC-L4。DC 前置门也不保证 AC 等价。

### L-C1：材料、栅和导数语义资格

- 材料支持体积必须是 WP-R 冻结的策略；后续变更触发 DC/P1/P2/AC 全链重验。
- Gate 从 P0 的材料合同确认：若保留多晶半导体栅，核对 DD 材料参数、掺杂、电极
  位置、载流子未知量与 poly depletion MOSCAP/解析极限；若替换为 metal_gate，
  另建 metal-equivalent fixture 和差异项，不认领原始 poly AC1 等价。
- 完整记录 Fermi、OldSlotboom、迁移率/高场 driving force/Enormal、SRH、温度和
  接触功函数。模型同名不证明参数、离散和导数语义相同。
- P0 的原始/effective Math 与 AC derivative 合同必须明确。只凭 Derivatives/
  AvalDerivatives 文本出现或缺席不能推断有效 AC J。若参考有效导数与完整物理
  Vela J 不同，作为受控 derivative-matched 变体单列比较；原始 oracle 不覆盖，
  Vela 不静默丢导数“追平”。无法确认时原始等价门 pending。

### L-C2：A/C 独立阈值与最终验收

P0/P5 按合同 4.9 的带单位 floor 定义独立键，不再用一个全矩阵门覆盖所有模型/区间：

| 阈值命名 | 冻结依据 |
| --- | --- |
| `A.significant/global` | DC 与模型/导数对齐、gm/gds 信号和参考数值噪声 |
| `C.cox_dominated` | 已资格化的积累/强反型条件、几何/介电/栅模型；可比 A 更严，但非普适 |
| `C.transition` / `C.off_diagonal` | 体积/多晶耗尽/电荷分配敏感度，独立 floor 和网格趋势 |
| `Y.complex` | 与独立 A/C 门联合验收，不能用大的电导淹没小的正交响应 |
| `Cgg.curve` | 分区间的中位/P95 和完整目标覆盖率 |

v3 的显著元素 5%、全矩阵 2%、Cgg 中位 3%/P95 5% 仅保留为历史起点，不是已批准
的统一新门；最终值在查看目标 Vela AC 曲线前按对应合同冻结，变更有版本和理由。

验收：L-A 硬门通过、DEVSIM PN/MOSCAP 合格；原始 AC1 所有目标 DC 前置门通过，
原始材料/导数合同闭合，全 A/C/Y 及 Cgg 分区指标达标；适用低频与合同 4.8 展开
差异目标 1% 并验证频率趋近。未关闭的 legacy volume、metal substitution 或未知
reference derivative 差异不能凭“已入账本”豁免 AC-L4。报告显示参考/Vela 值、误差、
阈值、资格状态和数据覆盖，保留失败点及部分产物。


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
4. 区分 DUT Y 提取与受源驱动网络分析：前者排除 DC 钳位偏置源，后者理想 AC
   电压源本来就是合法 MNA 支路，不能把所有未 Exclude 电压源都视作配置错误；
5. 实现显式 system 配置和 fail-closed parser；
6. 先验证解析 RC，再验证单器件显式/隐式等价；
7. 最后复现 Sentaurus `AC_des.cmd`。
8. 电路含电阻负载时先求一致的 device+circuit DC 工作点（可复用器件 DC 内核），
   不能把无负载器件偏压直接当成网络内部电压；将 per-m Y 乘实例物理宽度后 stamping；
9. Node 端口选择与内部节点消元分别定义。内部块非奇异时可用 Schur complement，
   含理想源/浮空约束时按完整 MNA 求解并诊断秩，不对奇异 Y 子块直接求逆。

### 验收标准

- R、C、串并联 RC 和多节点网络与解析 Y 的相对误差不大于 `1e-10`；
- MNA 每个非 ground 节点 KCL 归一残差不大于 `1e-10`；
- Vela 显式系统与 Vela implicit 系统全矩阵归一差异不大于 `1e-9`；
- Sentaurus `AC_des.cmd` 和 `AC1_des.cmd` 在同一冻结门槛内分别通过；
- 错误 Node、未约束的浮空网络、冲突理想源给出可诊断错误；合法 AC 源正常求解，
  DUT Y 提取中的偏置源钳位冲突有专门错误与测试；
- 多器件实例的端口命名无碰撞，输出包含 instance 和 node 标识。

## P7：性能、文档、CI 和长期扩展

### 实施内容

1. 为 J/M/block build、symbolic analysis、numeric factor、RHS solve、terminal extraction
   增加分阶段计时；
2. 记录矩阵维度、非零元、fill ratio、内存、频率数和端口数；
3. 维护 P3 已选后端与三路基准证据；后续更换表示/后端须重新验收，不重复后置选型；
4. 将小型解析/DEVSIM fixture 放入常规 CI；
5. 将 Sentaurus 对照放入人工或夜间回归，不把受许可文件提交仓库；
6. 编写用户配置、输出 schema、端口符号、错误诊断和示例文档；
7. 增加 Python API 仅在 C++ API 稳定后进行；
8. 评估 SDE 高斯植入支持，但不与 AC-L6 绑定；
9. 记录未来 transient、noise、S 参数和热 AC 的扩展接口，不提前实现。

### 性能验收

- 同一不可变矩阵不得随每个新增端口重复 numeric factorization，计数规则按 P3；
- AC 单频单工作点峰值内存暂定不超过同状态 DC Jacobian 求解峰值的 `3x`；
- real 路径为 `2*Ndof`、complex 为 `Ndof`，报告实际存储字节和非零元；dense 仅限小型 oracle；
- 1 个频率、4 个端口的 nMOS AC wall time 相对固定初始猜测/接受条件下同点 DC
  Newton 求解暂定目标为 `2x`；同时报告 DC 迭代数，不能用已收敛零迭代状态作分母。
  性能比例在 P3 实测后冻结，不作为物理正确性的替代门；
- 100 个频率时分别报告装配、分解、RHS 和端口提取耗时；非零 M 的新矩阵预期
  主要增加分解成本，M=0 允许跨频复用，不能强制要求其耗时随频率线性增长；
- 所有性能指标都在固定 Release 构建、固定主机、固定网格上测量。

### AC-L6 验收

- Debug 和 Release 全测试通过；
- 实现开始先记录现有失败基线；任何预存失败单列，不能认作新回归或静默豁免。
  若尚未达到全通过，只报告无新增回归和明确限制，不认领完整 AC-L6；
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
| WP-R (=P0.5) | 材料体积决策、独立装配器/replay、视图/导数资格 | P1/P2 共同前置；策略变化先重收敛 DC |
| WP1a | Poisson 接触反力 Q 和 Gauss closure | 不含独立生产 D-flux kernel |
| WP1b | Qx、full quasi-static C 和 FD oracle | 不含动态 AC |
| WP2a | 装配器内 `S(x)` 和 finite-difference M | 可与 P1 并行 |
| WP2b | 解析 M、边界行和 snapshot 审计 | 不含频域 solve |
| WP2.5 | 纯介电 Poisson AC | P1 + 最小 AC block；P3 资格门 |
| WP3a | 人工矩阵实块 builder | 不接器件 |
| WP3b | stacked/interleaved/complex 三路比较及统一计数 | P3 冻结生产后端，不后置到 AC-L6 |
| WP3c | 复用 WP-R snapshot 的 implicit RHS | 给定工作点 API；不后置设计共享快照 |
| WP4a | essential 反力/已资格化自然通量端口响应 | SG 路径仅作 oracle；含动态抵消和低频展开 |
| WP4b | Y/A/C、频扫、forced landing、输出/schema | 接入工作流并验收状态机 |
| WP5a | DEVSIM 基准 | 公开 fixture |
| WP5b | Sentaurus AC1：逐点 DC 门 -> 材料/导数合同 -> 分离 A/C 门 | 复用 SimpleMOS 方法；保留全部目标覆盖率 |
| WP6a | V/R/C MNA | 解析电路测试 |
| WP6b | device--circuit coupling | Vela explicit/implicit 对照 |
| WP6c | Sentaurus AC 显式验收 | 最后进入 |
| WP7 | 已选后端性能、CI、文档 | 不改变已冻结物理合同 |

每个提交至少包含：

- 对应测试；
- 变更前失败、变更后通过的证据；
- 新接口的单位和符号注释；
- 不改变的兼容行为说明；
- 相关 manifest/schema 版本变更；
- 若涉及数值门槛，附批准记录。

每阶段另交付 `stage_acceptance.json`：阶段 ID、输入/提交哈希、已实现接口、实际选择
的测试数、逐门结果、未支持模型、回归变化和审阅决定。WP-R 交付导数资格矩阵与
离散合同；P1/P2 交付 Q/M 数值审计；P3 交付计数/布局报告；P4 交付状态机和输出
合同报告；P5 交付完整覆盖率及差异账本。没有实际执行证据的门保持 pending。

本计划不将“纯介电一天”等估计当作交付承诺。WP-R 的物理导数补齐、P0 的许可/VM
可用性和 P6 的一致 DC 电路工作点是工期主要不确定项；实现团队在 WP-R 审计完成
后再按缺口估时，不能只按文件数量估算 P1/P4。

## 8. 测试矩阵和命令

### 8.1 测试矩阵

| 维度 | 必须覆盖 |
| --- | --- |
| 统计 | Maxwell--Boltzmann、Fermi |
| BGN | off、OldSlotboom |
| 材料 | 纯半导体、不对称 Si/oxide MOSCAP、已资格化多晶半导体栅与金属替代对照 |
| 端口 | 2 端、3 端、4 端 |
| 频率 | 单频、linear、decade、重复频率拒绝/去重策略 |
| 工作点 | 平衡、低偏置、导通晶体管；每个跨工具目标先过 DC 对齐门 |
| 边界行 | Ohmic/metal-gate Dirichlet、thermionic/Schottky、insulating pin、sheet charge |
| 缩放 | on/off、quasi-Fermi reference 变化 |
| QP | off；on 默认 fail-closed；显式 frozen approximation 仅在后续单独测试 |
| 雪崩 | off/on；同状态下无直接 M 项，重收敛状态的 M 可不同 |
| 调理 | 重掺杂 MOS，`1e3--1e12 Hz` |
| 线性表示/后端 | real stacked、real interleaved、local complex SparseLU；选择后端统一计数验收 |
| 输出 | per-m、per-um、finite depth、full/reduced matrix |
| 失败 | 非正频率、未知端口、奇异矩阵、非有限值、未收敛 DC |

### 8.2 外部审核指定的回归锚

| Fixture | 预期关系 | 所属阶段 |
| --- | --- | --- |
| 纯介电氧化层平板 AC | `Y=j*omega*C_analytic`，`A=0`；完整 3N solved gauge J=1/M=0 | P2.5 |
| 不对称 Si/oxide MOSCAP | 分离几何/材料策略；报告 Cgg 过渡区及 Cgs/Cgd 敏感度，登记差异 | WP-R/P1/P5 |
| 快照残差复现 | AC solved residual 与捕获的 DC 终态向量 bitwise 相同 | WP-R |
| 模式重建计数 | 每 snapshot AC pattern builds<=2，后续 RHS/频率为0 | WP-R/P3 |
| 多晶栅 MOSCAP | 半导体模型验证耗尽/串联电容极限；金属替代单独账本、不等价声明 | P0/P5 |
| DC 前置门 | 端电流符号/dex、psi、log n/p 和全部目标覆盖满足冻结合同 | P5 |
| V0!=1 的 BC 激励 | 指定方程行 RHS=1/V0；非零 reference 不改变该导数 | P3 |
| DC 路由覆盖 | 普通循环/显式 bias_points 各实际触发；受限控制/arclength 明确拒绝 | P4 |
| frozen reference FD RHS vs analytic stamping | 仿射增量 fixture `1e-12`；一般 FD 按 4.9 | P3 |
| continuity reaction DC current vs `computeFromResidual` | 逐接触机器精度一致 | P0.5/P4 |
| 同频多 RHS / M=0 跨频 | 首次新矩阵 1 次；缓存命中 0 次；失败单计 | P3 |
| 重掺杂 MOS `1e3--1e12 Hz` | backward error 不大于 `1e-10`；条件诊断如实报告 | P3 |
| AC + density-gradient QP 默认配置 | fail-closed，明确说明 frozen QP 非自洽 | P3 |
| thermionic/Schottky 和 insulating fixture | 自然通量行保留 M；essential/gauge 行 M=0 | P2 |
| avalanche-on 近击穿二极管 | 同快照 M 相同；全物理 J directional FD 与 Y 守恒通过 | P5 |
| 自适应 DC 强制 AC 落点 | P0 冻结清单；P4 验收命中、失败重试、只读观察者、恢复 | P0/P4 |
| 非零接触体积 manufactured 状态 | 接触 M 与 Q 中局部电荷导数抵消 | P2/P4 |
| Robin 自然接触 | 完整残差为零但边界电流非零；boundary stamp/J/du 一致 | P0.5/开放该功能前 |
| 触发 carrier diagonal floor 的 DC 状态 | AC J 不含 Newton floor；live residual FD 匹配 | P0.5/P3 |
| 导通器件低频展开 | `C_AC,0=H*x1+K*x0+E`；不强制等于 `C_Q` | P4/P5 |
| 电极电荷微分 Gauss | `sum C_Q[:,l]+dQ_domain/dV_l=0` | P1 |
| 缺列/失败的输出 | manifest 为 partial，失败列不填零，不参与完整矩阵评分 | P4 |

### 8.3 Windows UCRT64 基本命令

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --preset windows-ucrt64-debug
cmake --build --preset windows-ucrt64-debug
ctest --preset windows-ucrt64-debug
```

从当前 checkout/worktree 根目录运行，不跳回固定主目录。文档修改本身不需要构建
求解器。下列 AC 测试名/正则是拟新增的注册约定，当前尚不存在；先检查 `ctest -N`
确实选中预期测试，零测试不算通过：

```powershell
ctest --test-dir build -N -R "contact_reaction|coupled_dd_storage|ac_"
ctest --test-dir build --output-on-failure -R "contact_reaction|coupled_dd_storage|ac_"
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

reference_tcad/ac_small_signal/  # 拟新增的公开器件基准，遵循当前回归目录规范
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
| solved/physical 行视图混淆 | essential 行掩盖真实端口反力 | P0.5 独立 AC 装配器、视图与缓存门 |
| Si/oxide 材料支持体积混用 | 反型电荷、S/M、Cgg/端口分配与 A 偏差 | WP-R 前置策略决策，Poisson/S/M 联合资格；AC-DIFF-001 |
| DC 未对齐仍评分 A/C | AC 错误归因/选择性统计 | 每目标 DC 前置门、A/C 分合同、完整覆盖率；AC-DIFF-002 |
| 多晶栅被金属替代 | 栅耗尽与 Cgg/过渡区不等价 | P0 Gate 材料合同、poly fixture；AC-DIFF-003 |
| 参考有效导数未知 | 合法模型差异被误判为代码缺陷 | 原始/effective Math 与受控变体；AC-DIFF-004 |
| solved/physical 反复切换同一签名缓存 | 接触行被漏装或每频重复 rebuild | snapshot 独立作用域，每快照模式重建<=2 |
| 对 Robin 全残差提取端电流 | 收敛后错误输出零电流 | 同装配器 boundary stamp 分解和导数资格门 |
| Newton J 冒充物理 J | 遗漏物理响应或引入假电导 | 剥离迭代正则，按模型分支检验 live residual FD |
| 把 C_Q 当作所有端口低频 AC C | 错误验收/错误守恒约束 | 4.8 低频展开和 Gauss 微分测试 |
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

统一由规范合同 §11.1 的 required checklist 判定；下列风险出现时对应门为 fail/blocked：

- Sentaurus explicit/implicit 结果自身不能解释；
- 端口方向或二维单位尚未冻结；
- 材料支持体积尚未决策，或 P1 Gauss closure/P2 M 与移动电荷支持不一致；
- physical/solved row view 或自然边界分类尚未通过；
- scaling/reference invariance 失败；
- 同频率发生多次不必要 numeric factorization；
- KCL 只能通过删掉某个端口或经验缩放满足；
- 跨工具目标未通过 DC 前置门，或依赖不同物理电压/未经批准的插值；
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
- 三路后端均达正确性门后按 P3 实测选择；某候选失败不以放宽门槛掩盖。

## 11. 已冻结决策和剩余批准项

### 11.1 统一阶段纪律

见 [规范合同 §11.1](../specs/2026-08-31-ac-small-signal-contract.md)。
阶段产物 stage_acceptance.json 引用 required_gates 和证据哈希；pending/fail/blocked
不放行。入差异账本不是容差豁免，修改材料支持或导数合同会使旧数值证据失效。



### 11.2 仍需在对应阶段批准

1. WP-R 的 legacy_global/material_local 选择、适用材料/几何组合、是否实施物理子包；
2. P0 的 poly 栅映射与有效 reference AC derivative 语义，不能以 unknown 进入原始等价门；
3. 每点 DC 前置门、A/C 分区 floor 和最终阈值；参考变体与原始 AC1 的声明范围；
4. Robin AC 开放阶段、导数补齐范围、首个公开 DEVSIM fixture；legacy CV 默认保留；
5. P3 三路实测后的生产后端；P6 电路配置格式及额外 DC 路径资格；
6. Sentaurus 工件存储/访问/备份责任与中性证据移植，独立工作树结果需固定来源；
7. 独立数值测试证据与阶段放行。本文修订不替代这些批准。


## 12. 外部审核

请将本文、[规范合同](../specs/2026-08-31-ac-small-signal-contract.md)、[审核附录](2026-09-10-ac-small-signal-review-appendix.md) 和
[差异账本](2026-09-10-ac-known-difference-ledger.json) 一起交给其他模型；附录 §12 提供审核问题和回复格式。



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

以下模板仅供用户批准首轮实施后使用；阅读本计划本身不触发执行：

```text
请阅读：
docs/superpowers/plans/2026-08-31-ac-small-signal-simulation-development-plan.md
docs/superpowers/specs/2026-08-31-ac-small-signal-contract.md
docs/superpowers/plans/2026-09-10-ac-small-signal-review-appendix.md

当前只执行 P0a/P0b：Sentaurus oracle、强制目标 bias、合同/阈值冻结和公开解析
fixture。只冻结强制落点合同，不实现 DCSweep 观察者/裁剪。不要实现 physical row
view、AC solver、动态储存矩阵、MixedMode 或修改 cv_quasistatic 语义。

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

## 15. 证据入口

公开资料和当前源码/独立工作树来源集中在 [审核附录 §2.3、§2.5、§15](2026-09-10-ac-small-signal-review-appendix.md)。
本轮只修订文档并核对来源；未构建求解器、未运行 VM 或宣称数值门通过。
