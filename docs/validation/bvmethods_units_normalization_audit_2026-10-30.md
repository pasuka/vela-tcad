# BVMethods 单位转换与归一化往返审计（2026-10-30）

**可以直接归一化。当前问题是已经归一化的状态又被转成物理量再归一化，
以及某些物理场经过不必要的单位恒等往返；这会损失双精度最后一位。**
本审计轮次未修改这些单位转换或状态归档实现。后续修复与验证另见
[归一化修复记录](bvmethods_normalization_fix_2026-10-30.md)。
[此前SRH稳定导数修复及复验](bvmethods_srh_stable_derivative_2026-10-30.md)独立保留。

## 1. 为什么乘回V再除回V0不能保证恢复

物理输入电势可一次转换为 `u=psi/V0`，这一步有明确单位含义。
但若已拥有求解坐标u，不应为了保存检查点先算 `psi=fl(u*V0)`，随后再用
`fl(psi/V0)` 重建u。浮点乘除各自舍入，代数恒等不保证位模式恒等：

\[
\operatorname{fl}\!\left(\operatorname{fl}(uV_0)/V_0\right)\ne u.
\]

本例 `V0=0.025851999786435535 V`。真实预测态的节点3：
`21.200072448909626 → 21.20007244890963`，相差1 ULP。
归一化电势最大差 `3.552713678800501e-15`，折合电势约 `9.18e-17 V`。
这很小，但对近零SG电流及由其方向参与的AV线性化，不能仅凭绝对电压小就假定位模式等价。
尚未证明这些往返是当前全部收敛失败的唯一原因。

当前BVMethods人口模式的坐标契约：

| 数据 | 本例含义 |
|---|---|
| 求解向量psi块 | `psi/V0`，归一化电势 |
| 求解向量n/p块 | 物理内部浓度；本例为cm^-3，并非n/C0、p/C0 |
| 外接电阻内压lambda | V |
| C0／D0及线性行列缩放 | 方程行尺度及线性求解预处理，独立于上述状态存储 |

因此不能对整个状态向量再次统一除C0/V0；必须按明确的变量和单位契约处理。

## 2. 已确认的转换路径

### A. 人口状态pack/unpack存在多余电势往返

`src/equation/CoupledDDAssembler.cpp:1173` 的人口 `pack(CoupledDDState)` 分支先做
`state.psi * V0`，再交给 `packPopulation`；后者在1233行附近又除V0。
而 `unpack(x)` 的 `state.psi` 直接来自x的psi块，本来就是归一化坐标。
因此该组合的 `pack(unpack(x))` 不能保证逐位等于x，可直接复制归一化psi避免这一往返。
不能把物理电势输入接口 `packPopulation(psi_V,...)` 的必要单次转换混同为错误。

### B. DDSolution／归档的人口电势尚无原始归一化坐标保障

`src/solver/NewtonSolver.cpp:3530` 附近的 `makeSolution` 将psi乘V0输出为V；
人口模式的 `packReferencedSolution`（42行附近）又用物理psi调用 `packPopulation`，
并提前返回，未走QF模式已有的packed坐标恢复分支。

`src/io/DDSolutionState.cpp` 的人口归档已经额外保存n/p内部数值；同单位恢复为
`exact_internal` 时可逐位保留n/p。这是已有的有效保护。
但psi目前仍只以物理V存储；人口归档还明确拒绝旧QF的 `packed_psi` 组合，
所以不能把旧QF packed schema直接当作人口坐标使用。

现有 `Population continuation restart...` 测试明确使用 `1e-10` 等相对容差比较归一化psi，
并不承诺其位模式相同；人口归档位模式测试中的psi为零，也不能暴露非零电势往返。
HDF5双精度编码本身可逐位保存数据，问题发生在其上层的状态乘除重建。

### C. 热路径中的派生准费米势也有往返

`CoupledDDAssembler::unpack` 在约1216行计算物理准费米势后除V0，放入归一化的
`state.phin/phip`。`src/equation/BvMethodsCoupledRuntime.cpp:108` 附近的
`makeState` 又将它们乘V0，供SG／GradQF／AV使用。
这条 `physical_QF → normalized_QF → physical_QF` 链会重复舍入派生场。
它没有直接覆盖主未知量n/p，但可能影响微小QF差和电流方向，属于需做算子隔离的热路径问题。
本轮未独立导出全部物理QF原值，因此不捏造本算例QF改变节点数。

### D. 本例浓度进入SG节点统计时有可省略的SI往返

`include/vela/equation/BvMethodsTri3Operator.h:66` 使用
`units.internalConcentrationToM3(n) / 1e6`，DOS／ni初始化也使用同类链。
本例内部浓度已是cm^-3，换算因子为1e6，因此该表达式为 `(n*1e6)/1e6`，
总体单位因子为1，却仍可能改变最后一位。SRH直接使用原内部人口，SG节点统计经此转换，
可能在同次组装中看到不同的最后一位人口数值。

用已保存的真实packed数据按当前算式作binary64回放：

| 状态 | psi往返改变槽位 | n经SI改变槽位 | p经SI改变槽位 |
|---|---:|---:|---:|
| 第五点预测态 | 250 | 79 | 100 |
| 第五点失败态 | 246 | 102 | 92 |

这些改变均为1 ULP，浓度最大相对差约 `2.12e-16`。
这是单位恒等路径的表示误差，不能把它误报为几个数量级的浓度误差。
可先计算合并单位因子，再统一用于人口、DOS、ni及对应切线；本单位制中因子应精确为1。
需要同时复核值与导数，不能只改一个局部乘除次序就宣称整个SG契约已闭合。

### E. Poisson自由行可进一步采用原归一化坐标构造差分

`BvMethodsCoupledRuntime::populationBase` 当前先用物理psi差，汇总后再除
`epsRef*V0`；对应Jacobian含 `g*V0/(epsRef*V0)`。
在解析式上扩散部分可直接用归一化psi差，减少先乘V0再相减的舍入；电荷部分的V0仍须保留。
这是潜在的算术一致性改进，不是本轮已证明的量纲错误，修改前必须保持完整残差／Jacobian契约并验证Jv。
边界源、介电材料支撑及精确消元的RHS修正不能随之删除。

### F. 同一二进制的算子与方向控制

额外用同一个 `packed_probe` 可执行文件，对两个真实检查点分别输入物理重建态、原始packed态，
并重复原始packed态作为控制。结果：

- 两个检查点的原始RHS、bootstrap及完整Jacobian均逐位相同。
- 相同packed态的重复原始方向也逐位相同。
- 物理重建与原始态的方向不逐位相同。预测态各块以其无穷范数为尺度的差约
  `1.58e-8`（bootstrap）／`1.34e-7`（full）；失败态完整方向的psi／n块最大约 `2.78e-6`。

所以本例没有证明“往返一定改变RHS或Jacobian”。核心线性预处理的列尺度依赖当前状态
（例如psi尺度含 `|x_psi|+VT/V0`），最后一位坐标差可进入后续缩放。
本轮测到了最终数值方向的敏感性；尚未用完全固定的行列缩放独立闭合其全部来源，
也未证明它是冷第五点失败的唯一原因。不能把数值方向差直接解释为同等大小的物理电势误差。
同程序控制数据见ignored的 `unit_pair_analysis.json`，相同状态重复方向及矩阵均保存校验摘要。

## 3. 已检查的保护和范围

- `makePopulationTrial` 的主变量更新直接在packed向量中做，不要求每次更新都从物理场重新pack。
- 标准耦合延续在 `DCSweep.cpp:8655` 保存 `pendingCoupled.state`，下一点直接使用该packed状态。
  因此不能说当前完整扫描在每个接受点都必然重新往返psi。
- 一次启动、切换经DDSolution的求解接口、重启以及外部诊断重建仍会触及上述装包路径。
- 残差／Jacobian两侧的行列缩放及回代是合法预处理，不能与不必要的状态重建混为一谈；
  同态比较必须使用相同变量、行尺度、列尺度和约束消元口径。
- `UnitScalingConfig.mode` 与 `physicalUnitSystem` 是两个成员。手动给mode赋值不会自动改物理单位，
  构造函数传入mode则会选择对应单位制；需读取实际unitSystem，不能只靠mode字符串推断cm^-3。
- `DDSolution.n/p` 的老注释仍写m^-3，而实际归档代码按内部浓度处理并在SI字段写出前转换，
  是需要澄清的接口文档；不能据该注释对TCAD人口重复乘1e6。

本轮额外运行已有 `test_scaling`（31例／120断言）与 `test_state_archive`（13例／372断言），
全部通过。通过这些容差及归档测试不等于所有归一化坐标逐位不变；上面的零电势／容差覆盖边界需要补强。
本次未完成所有材料和跨工具绝对连续性行单位的全面审计，不能据此宣称全部物理归一化已经闭合。

## 4. 建议修复顺序

1. 给人口状态的psi建立与n/p相同的权威坐标契约，记录V0及参考／单位身份；
   同尺度同参考恢复时直接复制归一化psi，尺度或参考改变时明确转换。
2. 消除人口 `pack(unpack(x))` 的往返，并增加非零、易失1ULP输入的严格位模式回归。
3. 在同一packed状态中分别控制热路径QF往返与浓度SI恒等往返，量化SG电流、RHS、Jacobian与Jv变化。
4. 最后评估Poisson归一化差分重写，并复验重启、首五目标和完整BV；不以放宽容差代替算术闭合。

证据来自当前源码、已归档BVMethods原始检查点及本地unit/archive测试。
数值回放和日志位于ignored的 `build-release/bvmethods_srh_stable_20261030/`，
其中 `unit_roundtrip_counts.json` 保存逐项计数及具体例子。
