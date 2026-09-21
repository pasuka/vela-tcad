# 外推与延续算法研究：Vela 电热 Newton 的可借鉴方案

检索日期：2026-09-15。研究基于本地提交 `718ecdb`，只提出验证方案，
未引入新算法、第三方依赖或改变生产默认配置。引用为作者论文、项目官方文档和源码；
上游 `main/master/dev` 链接可能变化，实施时应再冻结实际版本。

## 结论与当前证据

建议顺序：受保护的历史外推及步长协调 → 切线预测 → 局部低阶拟合 →
停滞阶段 Anderson/非线性 GMRES。线性方程历史投影留待迭代线性后端阶段。
这些是结合本算例的工程判断，不是文献对 Vela 或 Sentaurus 加速效果的保证。

当前 R7 已启用全场线性外推：电势、物理电子/空穴准费米势、晶格温度。
它不是只预测端电流。生产后端为 UMFPACK 直接求解；减少外层 Newton 次数，
可以减少 Jacobian 数值分解，而给直接求解器一个更好的初值本身没有同样作用。

| 双栅压完整曲线，Vg=4 / 8 V | R7 | 可选 R8 |
|---|---:|---:|
| 漏压 Newton 更新，含失败尝试、不含初始化 | 366 / 321 | 354 / 314 |
| 线搜索候选 | 496 / 446 | 448 / 401 |
| 非零推进失败尝试 | 1 / 0 | 1 / 1 |

R8 总更新减少 2.77%，线搜索候选减少 9.87%；没有证明稳定墙钟加速。
R7 仍为默认。R8 在 Vg=8 V 的 0–4 V 区间仍为 86 次更新，局部改进被
恢复开销抵消。细节见[重复性与 Newton 定位](templates_ldmos_r7_stability_newton_2026-09-14.md)。

本地实现位置：

- [`ElectrothermalSweep.cpp`](../../src/simulation/ElectrothermalSweep.cpp)：
  `predict`、历史点选择、步长增长和回退。
- [`ElectrothermalSimulation.cpp`](../../src/simulation/ElectrothermalSimulation.cpp)：
  单点 Newton、线搜索、状态表示与密度更新。
- [`ElectrothermalIterationControl.h`](../../include/vela/solver/ElectrothermalIterationControl.h)：
  直接求解后端、结构复用与停滞控制。

## 必须区分的三层问题

设器件残差为 `F(x,Vd)=0`，`x=(psi,fn,fp,T)`；准费米势在程序中采用参考值与
增量分离存储。以下公式里的状态差必须先统一物理参考。

| 层级 | 被预测或加速的对象 | 代表方法 | 对当前任务的意义 |
|---|---|---|---|
| 扫描延续 | 下一个电压点的非线性初值 | 割线、切线、局部多项式 | 最直接减少漏压 Newton 更新 |
| 单个电压点 | 非线性迭代序列 `x(k)` | Anderson、非线性 GMRES | 针对固定点/分块迭代或停滞恢复 |
| 单次 Newton 内部 | `J dx=-F` 的线性解 | 历史子空间投影、Krylov 复用 | 需迭代线性后端及预条件；不是 UMFPACK 的初值开关 |

预测值只作初值，最终仍需原四方程、逐行、KCL、热平衡和参考结果验收，
不能用预测值或曲线插值代替精确参考点的求解。

## 1. 受保护的割线预测和步长协调：优先验证

两点割线预测为

`x_pred = x_n + (V_target-V_n)/(V_n-V_previous) * (x_n-x_previous)`。

LOCA 将预测、非线性校正和步长控制分离，提供首步预测器与割线预测器。
其自适应步长实现根据校正迭代数平滑调整，并裁剪到上下限；该实现本身仍
使用迭代数，不能将其描述为已经具备预测残差质量门限。
参见 [LOCA 割线源码](https://raw.githubusercontent.com/trilinos/Trilinos/master/packages/nox/src-loca/src/LOCA_MultiPredictor_Secant.C)
及[自适应步长源码](https://raw.githubusercontent.com/trilinos/Trilinos/master/packages/nox/src-loca/src/LOCA_StepSize_Adaptive.C)。

本地问题：当前历史选择只要求已接受、正间距、外推步长比≤3；近点不满足时可
改用更老的点。`predict` 检查有限值、温度范围和参考电势兼容，但未比较候选
的残差质量。计划步长在命中精确输出点后也没有按实际截短步长重置。

R8 的一次具体轨迹为：1.21875 → 1.333333 V，然后尝试 2.472396 V。
若用最近历史点，步长比约9.94，故改用0.7125 V历史点，比例约1.835。
Vg=8 V 该尝试10次更新后停滞，随后减步。它证明该处值得控制变量比较；
尚不能单凭一次轨迹证明“老历史点”是唯一原因。

建议独立实验：

1. 记录计划步长、实际步长、输出点截短、历史跨度及历史点年龄，避免只看更新数。
2. 在相同目标电压、相同边界投影、相同参考及缩放下，比较割线初值与恒定状态初值
   的残差；结合各方程块和逐行指标拒绝变差预测。残差检查自身也计入耗时。
3. 对输出点截短和历史跨度突变增加增长保护。先分别测试历史策略和步长策略，
   再与 R8 密度更新联动，避免多个变化掩盖原因。
4. 预测失败回退到合格状态或缩小步长；不提高 Newton 验收容差。

## 2. 切线/灵敏度预测：可利用现有 Jacobian

对 `F(x(V),V)=0` 求导，得到

`J(x_n,V_n) t = -F_V(x_n,V_n)`，`x_pred = x_n + deltaV * t`。

LOCA 的 `Tangent::compute` 先计算参数导数和 Jacobian，再通过线性求解得到
预测方向；这是可以借鉴的 C++ 实现。
见[切线预测源码](https://raw.githubusercontent.com/trilinos/Trilinos/master/packages/nox/src-loca/src/LOCA_MultiPredictor_Tangent.C)。

对 Vela 的潜在收益是减少对不均匀历史点的依赖。若当前合格状态的 Jacobian
恰好已数值分解，多一个右端求解可能较便宜；但 Newton 最后保存的分解未必对应
最终状态，特别是参考重定心或末次更新之后。必须核对矩阵与状态，额外装配/分解
成本不能省略计时。结构分析可复用，不等于数值因子可以跨状态沿用。

实现前应在固定状态核验 `F_V` 和切线预测误差随小步长的阶次。参数偏导需采用
真实漏极边界的定义，包括中性接触及有限复合边界的电压依赖；温度链式项已在
完整 `J` 中时不能重复加入。有限差分参数导数可作诊断对照，不直接认领解析实现。

## 3. 局部低阶多项式和多候选预测：处理不均匀历史

BifurcationKit 支持割线、切线、局部多项式和多候选预测。其多项式方案对历史
弧长归一化后做低阶最小二乘拟合；多候选方案将候选与校正残差下降检查结合。
见[官方预测/校正文档](https://bifurcationkit.github.io/BifurcationKitDocs.jl/dev/Predictors/)
及[延续主循环源码](https://github.com/bifurcationkit/BifurcationKit.jl/blob/master/src/Continuation.jl)。

可借鉴归一化和局部拟合思想，在 Vela 的电压参数上构造低阶拟合，这属于改编，
不等同于原项目的弧长算法。采用少量近期合格状态、QR/秩检测、变量缩放与降阶
回退；不直接使用跨全曲线高阶插值。先比较常值/线性/二次候选的准备成本和
真实残差，再决定是否值得多启动一次 Newton 校正。

若未来研究折返点或回滞，可另考虑伪弧长延续。1989年的半导体延续研究已将
预测—校正用于 CMOS 的 V–I 曲线及触发/保持极限点，说明该方法与器件问题直接相关；
当前已能完成的单调漏压扫描还没有必须改为伪弧长的证据。
见 Bell Labs 原始研究条目
[Continuation Methods in Semiconductor Device Simulation](https://www.nokia.com/bell-labs/publications-and-media/publications/continuation-methods-in-semiconductor-device-simulation/)。

## 4. Anderson / 非线性 GMRES：作为停滞恢复候选

Walker 与 Ni 的论文研究了固定点 `x=g(x)` 的历史组合加速，并在特定线性、
无截断、非停滞条件下建立与 GMRES 的联系，不能把该结论直接扩展为任意
非线性 Newton 的等价关系。
见 H. F. Walker, P. Ni, *Anderson acceleration for fixed-point iterations*,
SIAM J. Numerical Analysis 49 (2011), 1715–1735，
[作者全文](https://users.wpi.edu/~walker/Papers/Walker-Ni,SINUM,V49,1715-1735.pdf)。

对 `f_k=g(x_k)-x_k`，一种写法为
`gamma = argmin ||f_k - DeltaF gamma||`，
`x_next = g(x_k)-DeltaG gamma`。短历史 QR 更新、病态列删除、阻尼及延迟启动
可参考 Walker 的[算法与实现说明](https://users.wpi.edu/~walker/Papers/anderson_accn_algs_imps.pdf)。

可直接阅读的开源实现：

- **SUNDIALS KINSOL**：`KINFP`、`KINPicardAA`、`AndersonAcc`。
  [C 源码](https://raw.githubusercontent.com/LLNL/sundials/main/src/kinsol/kinsol.c)；
  [使用文档](https://sundials.readthedocs.io/en/latest/kinsol/Usage/index.html)。
  Anderson 用于固定点/Picard 策略，不能把 `KIN_LINESEARCH` 说成同一个 AA 策略。
- **PETSc SNESNGMRES**：有候选选择、线搜索与历史重启接口；外层算法不使用
  用户提供的 Jacobian，可与非线性预条件器组合。
  [官方说明](https://petsc.org/main/manualpages/SNES/SNESNGMRES/)，
  [C 源码](https://petsc.org/main/src/snes/impls/ngmres/snesngmres.c.html)。
  组合方法的论文为 Brune、Knepley、Smith、Tu，
  *Composing scalable nonlinear algebraic solvers*, SIAM Review 57(4), 535–565 (2015)，
  [DOI](https://doi.org/10.1137/130936725)。

对 Vela，先定义合适的固定点/分块校正映射，再试加速；不直接混合不同电压的
Newton 残差。偏压、物理模型、状态坐标或约束改变时清理历史，或证明转换后一致。
混合后必须重算物理残差、检查温度与密度可行性，并保留原 Newton 回退。

AA 不保证改善已进入二次收敛的 Newton。Evans、Pollock、Rebholz、Xiao
专门讨论了线性收敛与二次收敛的区别，支持将此项优先放在缓慢固定点或
停滞阶段，而非默认施加到每一步 Newton。
见[论文及作者摘要](https://arxiv.org/abs/1810.08455)。

## 5. 线性矩阵历史投影：迭代后端的后续路线

P. F. Fischer, *Projection techniques for iterative solution of Ax=b with
successive right-hand sides*, Computer Methods in Applied Mechanics and
Engineering 163(1–4), 193–204 (1998)，研究从历史解构造后续线性问题初值。
PETSc 的[官方实现说明](https://petsc.org/main/manualpages/KSP/KSPGUESSFISCHER/)
提供论文信息和使用接口；[fischer.c](https://petsc.org/release/src/ksp/ksp/guess/impls/fischer/fischer.c.html)
实现历史空间、投影、更新和重置。

它能为迭代线性求解提供非零初值，尤其适合相关右端序列。PETSc 的模型2要求
正定能量范数，不能直接套用本算例一般非对称的耦合 Jacobian。
Vela 若进入该路线，须先建立迭代后端和预条件，再处理变化的 Jacobian、
参考坐标以及旧历史失效；基于当前矩阵重算投影作用或清空历史。
当前 UMFPACK 基线不因启用这一类“初值外推”自动减少分解次数。

## 建议的验证矩阵与停止条件

以下为下一阶段可执行设计，不代表本轮已实施。

| 顺序 | 固定条件 | 单项变化 | 主要判据 |
|---|---|---|---|
| A | R7物理、坐标、线性后端、既定低压目标序列 | 历史/候选残差保护 | 更新、失败、装配、候选检查成本 |
| B | A所用物理和原验收标准 | 自适应增长与精确输出点截短策略 | 完整推进成本、不得只挑成功尝试 |
| C | 相同目标序列及后端 | 割线与切线预测对照 | 节省的Newton分解能否抵消切线准备 |
| D | 独立候选配置 | 与R8密度更新联动、局部拟合 | 是否消除约1.33–2.47 V的额外恢复 |
| E | 明确的单点停滞样本、原门限 | AA/NGMRES恢复 | 真实残差、正性、失败回退和总成本 |

从两个栅压的0–4 V问题区间及原前8个精确点开始，保留全部尝试记录。
区分同目标序列的算法对照与自适应完整路径对照。局部有效后，再独立完成
双栅压0–40 V完整62点及同VM重复配对计时；原电学/热学、局部密度2%、
带边30 meV及墙钟≤原生1.5倍口径不变。涉及共享输运/物性时再复核G3/D5/D4。

必须统计初始化、失败、候选筛选、Jacobian装配/分解、线搜索、输出和父运行器
成本；保存/恢复需维持预测历史，不能只验证连续运行。已有R7的资格不能自动
转移到新预测器。若更新减少但总时间增加，保留诊断结果，不晋级默认配置。

## 本轮提交检查

`718ecdb` 保存了生产集成、R7/R8与既有报告/图件、环境研究；原始大型
仿真证据仍留在忽略目录，版本化配置记录输入及程序来源。
Windows UCRT64 Release全目标构建成功；Release CTest **817/817**通过，
用时139.22 s；`git diff --check`通过。此文是检索与设计记录，没有新增
仿真耗时或算法收益的实测认领，也不能据公开实现推断Sentaurus内部源码相同。
