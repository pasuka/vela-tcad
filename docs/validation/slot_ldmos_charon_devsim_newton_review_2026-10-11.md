# SLOT：本地 Charon / DEVSIM Newton 实现核查

2026-10-11。仅检索用户指定的两个本地仓库及其手册，没有编译、运行或修改这两个工程，也没有据此改变 Vela 物理、默认更新规则或验收。

## 结论及当前问题

找到了与当前问题直接相关的处理方式：DEVSIM 的**浓度正性更新与电势对数阻尼分开处理**；Charon 的**变量缩放、带停滞保护的更新误差检查**；以及 Charon 显式提供的**低密度电离驱动力抑制**。最后一项改变模型评估，应与纯数值策略分开。

当前 SLOT 已经把行比改为报告，但仍保留全局闭合及边界保护。最新无限幅 Bank–Rose 路径面临的是巨大 QF 方向、试探人口出域以及极小统一阻尼下的闭合停滞。不能再把全部困难归结为 `eps_row=0.001`。同一冻结系统的首次最大空穴 QF 方向约为 −1.01e8 V；这不是合适的完整有限步。

本轮同时完成的有效域回退复验见 [Bank–Rose 人口有效域恢复](slot_ldmos_bank_rose_population_recovery_2026-10-11.md)。本文件的外部代码结论不以该复验最终是否成功为前提。

## 来源身份与核查范围

| 本地工程 | HEAD | 主要内容 |
|---|---|---|
| `D:/code-repo/devsim` | `58a9a87083db00c6cadc0b4011c801db2cec5844` | `src/math/Newton.cc`、`src/Equation/Equation.cc`、标准漂移扩散 Python 模型；本地手册标记 Release 2.10.0 |
| `D:/code-repo/tcad-charon` | `7cc38745625a6011ae3584ed111ec7ee74fb890e` | NOX 配置与适配层、缩放、迁移率/电离 evaluator、用户手册第15章 |

这是本地检出版本，不代表上游最新版；实际读取文件逐一计算 SHA-256，见 `build/s11external_newton_review/manifest.json`。HEAD 从本地 Git 引用核对。Charon 的完整 Newton/线搜索实现依赖外部 Trilinos NOX；本次两个仓库中没有查到 NOX 核心实现，因此不宣称核对了其全部算法细节。

## 1. DEVSIM：分别限制不同变量的更新

标准模型的 `CreateECE/CreateHCE` 以 `Electrons/Holes` 为未知量，配置 `variable_update="positive"`；电势 `Potential` 使用 `log_damp`。

来源：[simple_physics.py](/D:/code-repo/devsim/python_packages/simple_physics.py:269)、[Equation.cc](/D:/code-repo/devsim/src/Equation/Equation.cc:151)。

### 浓度正性

`PositiveSolutionUpdate` 首先计算 `n_trial = n_old + delta_n`。父态必须为正；若试探非正，改为 `n_new = 0.001*n_old`，并把实际增量改为 `n_new-n_old`。极端下溢时还有 `0.5*n_old`、保持父值的后备分支。空穴同理。

这里的系数是 **0.001，不是此前 Sentaurus 观察到的0.01**；这一函数没有给正浓度增长设置 `1e5` 比例上限。该更新避免非正人口，但不能据这段代码宣称任意溢出/NaN 均有保护，也不能保证最终残差下降。

### 电势对数阻尼

源码在 `|delta_psi| > 0.0259` 时采用

\[
\delta\psi_{app}=\operatorname{sign}(\delta\psi)\,0.0259\,
\ln\left(1+\frac{|\delta\psi|}{0.0259}\right),
\]

其余使用原位移。这里 `0.0259` 是实现中的固定常数，不能描述为自动随温度变化的热电压。它是对电势的更新规则，不是给所有 QF 分量统一乘一个小系数。

### 对 SLOT 的意义与限制

这种分变量处理不会仅因一个低密度节点的巨量 QF 方向，就把所有电子输运和电势更新按相同比例压小。这是值得检验的数值策略，但标准 DEVSIM 示例并不是当前 SLOT 完整 Fermi/PhuMob/电离系统的等价模型，不能直接搬用其 SG 公式、门槛或预测成功率。

此前已经做过[同线性系统人口有限步 C/C′ 对照](slot_ldmos_pure_population_update_comparison_2026-10-10.md)：空穴坏行由 A 的4194降为3414，但仍有新增坏行，未通过当时完整资格。那一批使用另一冻结配置、6 mV目标和 `[0.01,1e5]` 投影，不是这次最新5→10 mV系统。因此不能把“人口更新”作为从未试过的新方法，更不能说已有代码证明它能解决问题。

有辨别力的新对照应明确差异：在当前相同父态、相同线性系统上，比较已存在的 QF 更新与 **DEVSIM 式逐节点正性人口更新 + 电势单独阻尼**；保留 Fermi 相对逆映射和边界掩码，记录最终实际人口，而不是只记录投影目标。若另做浓度坐标线性求解，必须同时核对物理 Jacobian、正则化、闭合替换和约束的坐标变换。

建议候选状态的构造顺序是：从同一原始方向按完整链式关系得到人口线性增量；分别确定阻尼后的 ψ 和正性处理后的 n/p；再以这个**最终 ψ 与最终 n/p**共同恢复 Fermi QF。不要先恢复 QF、随后独立覆盖 ψ，使实际人口又偏离正性投影目标。若需要 Poisson 修复，应重做协调更新或重新求解，不能把投影前后不同人口混为一个状态。此处是基于已有 C/C′ 证据提出的 Vela 实验设计，并非 DEVSIM 源码中已有的 Fermi-QF 适配实现。

`0.001` 以及是否设置正增长上限都应作为明确的对照参数；仍需检查非有限值、源核有效域、邻行和真实残差，不建议因为 DEVSIM 这一个函数没有上限就无条件删除所有保护。

设 QF 变量为 x、人口变量为 y、x=G(y)，则物理 Jacobian 应满足 `J_y = J_x * dG/dy`。对同一个稳定化系统做精确链式变换，在精确算术下不会凭空产生另一条物理线性方向；可能的收益来自条件数、有限步更新及一致的可行域处理，而不是“换变量就一定修好”。Fermi、本构及电势耦合必须纳入变换，不能只对载流子列乘一个任意常数。

## 2. DEVSIM：更新误差、失败恢复与稳定求差

`DefaultUpdate` 在应用上述规则后计算每个方程的

\[
e_{abs}=\max_i|\delta x_{app,i}|,\qquad
e_{rel}=\max_i\frac{|\delta x_{app,i}|}{|x_{new,i}|+e_{min}}.
\]

`e_min` 默认 `1e-10`，含义随变量及其单位而定；不是 Sentaurus 的 `ErrRef=1e10 cm^-3`。区域/器件再组合这些方程误差，不能把全器件指标直接写成一个统一的节点最大值。

DC Newton 循环每步重新装配；对所有器件要求**绝对误差和相对误差同时满足**。检查的是更新误差，不能把它等同于 Vela 的整体 RHS、逐行通量比或分载流子守恒。检查的 DC 循环没有通用 Armijo/Bank–Rose 试探线搜索。失败后 `RestoreSolutions()` 恢复父态；`ramp.py` 对识别出的收敛失败恢复偏压、将步长减半，其他异常照常抛出。这个 helper 在成功时不自动增步。

来源：[Equation.cc](/D:/code-repo/devsim/src/Equation/Equation.cc:237)、[Region.cc](/D:/code-repo/devsim/src/Geometry/Region.cc:1305)、[Device.cc](/D:/code-repo/devsim/src/Geometry/Device.cc:338)、[Newton.cc](/D:/code-repo/devsim/src/math/Newton.cc:673)、[ramp.py](/D:/code-repo/devsim/python_packages/ramp.py:55)。手册 PDF第54页及60页介绍 DC Newton 参数与 `maximum_error/maximum_divergence`。

此外，[simple_dd.py](/D:/code-repo/devsim/python_packages/simple_dd.py:13) 使用 Bernoulli 恒等式以及 `kahan3` 组合 SG 电流项；Poisson 电荷也使用补偿求和。这支持继续把“近零通量求差”与“大人口更新”分开诊断，但不能据此认定 Vela 缺少对应稳定处理，更不能未经 Fermi 推导直接替换当前有效带边 SG。

## 3. Charon：缩放、NOX 与小步停滞保护

### 原生变量和尺度

标准漂移扩散方程以电势、电子浓度、空穴浓度为自由度。缩放代码定义

`psi/V0, n/C0, p/C0`，其中 `V0 = kB*T0/q`；同时定义 `R0=D0*C0/X0^2`、`J0=q*D0*C0/X0` 等源项/电流尺度。代码的 kB 使用 eV/K，故实现为 `kb*T0`。

来源：[EquationSet_DriftDiffusion](/D:/code-repo/tcad-charon/src/equation_sets/Charon_EquationSet_DriftDiffusion_impl.hpp:344)、[Scaling_Parameters.cpp](/D:/code-repo/tcad-charon/src/Charon_Scaling_Parameters.cpp:11)。Vela 已有自身单位缩放，本项不是“Vela 尚未归一化”的证据；应比较缩放后各块、正则化和停止指标的一致性。

### 收敛和停滞区分

Charon 通过 NOX 配置 Newton，可选择 Full Step 或 Backtrack。不同 solver pack 的停止树不同：例如 pack5/6 把 NormF 和 NormWRMS 放在 AND 成功分支；还有 NormF 为主的配置。树外层的 OR 包含 MaxIters 停止项，**达到迭代上限不是成功**，主程序会检查 NOX 返回的 Converged 状态。

用户手册第154页给出

\[
\|\delta x\|_{WRMS}=\sqrt{\frac1N\sum_i
\left(\frac{x_i^k-x_i^{k-1}}
{\epsilon_{rel}|x_i^{k-1}|+\epsilon_{abs}}\right)^2}.
\]

同时列出两个可选保护：线性求解达到的误差需足够小；线搜索实际步长比例需大于指定下限。后者避免把**阻尼压出来的小更新**误认为已接近方程根。手册 alpha=0、beta=1 分别禁用这两项；不能将示例 solver pack 数值直接抄作 SLOT 推荐值，具体 NOX 版本实现还需另核。

来源：[createSolverPack5.inp](/D:/code-repo/tcad-charon/scripts/charonInterpreter/parseGenerator/parseInputs/SolverBlock/createSolverPack5.inp:44)、[createLineSearchMethod.inp](/D:/code-repo/tcad-charon/scripts/charonInterpreter/parseGenerator/parseInputs/SolverBlock/createLineSearchMethod.inp:1)、[Charon_Main.cpp](/D:/code-repo/tcad-charon/src/Charon_Main.cpp:1068)，以及 `D:/code-repo/tcad-charon/Charon_UserManual.pdf` 第153–154页；第154页公式已渲染核对。

对当前 SLOT，最实用的借鉴是把原始方向误差、实际有限更新、阻尼比例、真实 RHS 与全局闭合并列记录。现有诊断已区分前两者；后续可据连续多步趋势识别停滞，不能简单改为“实际更新小即成功”。这不要求重新把所有弱行设为阻断条件。

### 续算与正性代码的边界

Charon 暴露 LOCA Natural/Arc Length 续算以及步长 Aggressiveness 配置；本次只确认接口，未取得其外部 LOCA 核心算法。当前低漏压困难没有已证实的分支折点，不应只因存在弧长选项就认定必须改为弧长续算。

另发现 `Charon_RythmosObserver_ClipSolution.hpp` 可以把负值置零，但调用位置是瞬态**已完成时间步观察器**，不能说这是稳态 Newton 默认正性算法。Vela 当前源核需要严格正人口，直接照搬“置零”反而会触发其有效域检查。QuasiFermiLevels evaluator 对派生 QF 的浓度下限也不能当作全部 Newton 方程的正性保证。

## 4. Charon：低密度驱动力处理必须单独标作模型变体

### 电离源

`Avalanche_vanOverstraeten` 明确提供

\[
F_{n,eff}=\frac{n}{n+n_0}F_n,\qquad
F_{p,eff}=\frac{p}{p+p_0}F_p.
\]

参数 `eDrForceRefDens/hDrForceRefDens` 默认均为0；注释明确针对低密度区强产生/复合导致的收敛困难。正参考密度会抑制低密度驱动力，改变电离源值及其导数；**不是仅冻结 Jacobian，也不是向电场 |E| 插值回退**。

来源：[Avalanche_vanOverstraeten_impl.hpp](/D:/code-repo/tcad-charon/src/evaluators/Charon_Avalanche_vanOverstraeten_impl.hpp:31)，手册第94–95页有同类驱动力阻尼说明。当前 SLOT 曾在无电离配置仍出现弱行困难，因此它是独立物理消融候选，不能作为唯一解释或直接启用的修复。

### 高场迁移率

Philips–Thomas 的边上 GradQuasiFermi 路径用

\[
\frac{\nabla n}{n}\bigg|_{edge}
\approx \frac{2(n_j-n_i)}{h(n_i+n_j)}
\]

构造 QF 梯度；若两端人口为正，这个**密度梯度比值项**的绝对值小于 `2/h`。它与 `log(n_j/n_i)/h` 在极端人口比下不同；整个驱动力还包含有效电场，不能说整体被 `2/h` 限住。这里的有界性是从源码公式推导的，不是 SLOT 实验结论。

同 evaluator 的积分点路径对平方范数加 `1e-20` 后开方，注释说明避免零点 Jacobian NaN；边上路径虽有相同初值注释，但后续会覆盖，不能笼统宣称所有边都有这个下限。默认驱动力还是 ElectricField，非本例指定的 QF 梯度。

来源：[Mobility_PhilipsThomas_impl.hpp](/D:/code-repo/tcad-charon/src/evaluators/Charon_Mobility_PhilipsThomas_impl.hpp:697)。这些都是需要重新验证残差、导数、Fermi 兼容性和原生映射的模型/离散变化，不作为保持原方程不动点的数值改进。

## 5. 建议的执行优先级

| 优先级 | 独立对照 | 判断依据 |
|---|---|---|
| 1 | 当前同父态的分变量人口正性更新与电势阻尼；与无限幅 Bank–Rose 保持同物理、后端、约束和最终闭合保护 | 能否避免极小统一步，同时改善真实 RHS、电子/空穴闭合；记录所有投影与新增坏行，不以仅减少弱行数量判断 |
| 2 | 补充停滞分类与更新坐标审计 | 区分原始方向小、投影后小、阻尼后小；核对线性后向误差及连续多步实际改善，未闭合不宣称成功 |
| 3 | 对关键弱行冻结态比较 Charon 式 HFS 梯度及电离参考密度 | 仅作模型变体，分别观察源/输运值与导数、端口分量及热点，不能与纯更新对照混在同一臂 |

可先用已有捕获系统验证链式转换和一次有限试探，发现方向或邻行明显恶化即停止该方案；再选择少量重线性化。无需先安装/移植整个 NOX，也没有证据要求先做全链路 float80。本次仅提出这些对照，未执行新的策略或模型试算。

## 证据文件

- `build/s11external_newton_review/manifest.json`：源码、PDF、HEAD与哈希。
- `devsim.txt`、`charon.txt`：带 PDF 页码的本地手册提取文本。
- `charon_newton154.png`：停止条件公式的可视核对。
- 本次没有证明 Charon 或 DEVSIM 在此 SLOT 网格及完全相同材料/源核下能收敛；方法存在与对本算例有效是两个不同结论。
