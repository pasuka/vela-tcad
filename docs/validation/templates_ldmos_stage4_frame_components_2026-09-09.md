# LDMOS Vg=8：参考系切换后电子残差的分项因果对照

日期：2026-09-09，Asia/Shanghai。承接
[前一轮参考系定位](templates_ldmos_stage4_frame_diagnosis_2026-09-09.md)。

## 结果

在已经补偿外部参考平移低位、修正 Newton QF 重打包顺序的基础上，剩余电子
残差的主要障碍已定位到 **vector HFS 场强计算先重建绝对 QF 电势，再求空间差分**。
两项独立干预支持同一结论：固定为原参考系的场强，或直接用参考值与增量的
相对差计算每个三角形梯度，都使平移态的电子残差通过原 `1e-11` 门限。

相对梯度诊断版本从补偿平移态出发，0 次 Newton 更新即通过原全局与局部门限；
保存、重读仍通过。从上一轮已经停滞的状态出发，1 次更新后电子残差降到
`3.83047e-12`，也通过。电势、密度、端电流和 KCL 的原换参考等价性门限同时通过。

这次是独立诊断副本中的单点因果验证，**没有修改正式源码、控制器或原证据，
没有续跑 26.6667–40 V，也不是 D5 整体验收通过声明**。

## 配置与隔离

- 同一 10,241 节点、19,782 三角形、30,022 边网格；几何输入 μm、密度输入
  m⁻³、电势 V、端电流 A/μm；`unit_scaling`。
- Vg=8 V、物理 Vd=26.6666666666667 V；frame 0 与整体减 28 V 的 frame 28。
  frame 28 的 source/substrate=-28 V、gate=-20 V、drain≈-1.33333 V。
- AverageBox 输运，barycentric 体积、material-local Poisson、legacy node-local
  接触；300 K、Fermi–Dirac、OldSlotboom BGN、SRH/Auger、constant_field 和
  transport_cell_vector HFS；原接触场强回退策略保留。量子、雪崩和 predictor 关闭。
- UCRT64 GCC、Release `-O3 -DNDEBUG`、Eigen SparseLU/COLAMD、L2 行列平衡。
- 全局 psi/e/h 门限 `5e-8 / 1e-11 / 3e-10`，局部 `eps_row=1e-8`、违规零；
  KCL ratio≤1e-8。单点上限40次，失败对照在13次后因线搜索拒绝结束，没有用尽预算。
- 原合格父状态 SHA256：
  `f861982783bbdc730d5eb742aa6406430a85488244d15f72574d8c64d6a2dec2`。
- 新诊断二进制 SHA256：
  `a93745d56562d553b22fd58bf1da2b27158dead045710ddef9c15d9549e85846`。

所有产物位于[本轮诊断目录](../../reference_staging/templates_ldmos_frame_components_20260909/)。
复用上一轮 `NewtonSolver.repack.obj`，只在新的 `CoupledDDAssembler.components.cpp`
副本加入运行时可选干预。无干预时，两个参考系的完整 residual CSV 与上一轮
repack 诊断逐字节相同；正式 assembler、NewtonSolver、静态库、runner 和原计划
全部冻结输入的哈希检查通过。

## 1. 分离密度、场强、SG 求值和残差累加

以下均为**同一补偿平移输入、同一接触投影、同一重定位策略**的原始残差。
表中的“通过”只描述电子块是否≤1e-11；完整闭合另列于后。

| 单项干预 | psi 残差 | 电子残差 | 电子块 |
| --- | ---: | ---: | --- |
| 无干预 control | 4.925889e-8 | 1.460242e-10 | 未通过 |
| 密度重建中的 `scaled_psi*V0-reference` 用 long double 中间量 | 4.085550e-8 | 1.463662e-10 | 未通过 |
| 固定 n/p 为原参考系重建值 | 4.978718e-9 | 1.462987e-10 | 未通过 |
| 固定 vector 场强为原参考系值 | 4.925889e-8 | 9.679528e-12 | 通过 |
| 每个三角形用 reference/increment 相对差计算 vector 场强 | 4.925889e-8 | 8.801839e-12 | 通过 |
| 仅电子 SG drop 与通量末端求值用 long double | 4.925889e-8 | 1.460283e-10 | 未通过 |
| 仅残差累加器用 long double | 4.925889e-8 | 1.460246e-10 | 未通过 |
| 同时固定原参考系 n/p 和场强 | 4.978718e-9 | 1.055642e-11 | 未通过 |
| 密度、相对 vector、SG、累加四项组合 | 4.085550e-8 | 8.995769e-12 | 通过 |

固定密度显著恢复 Poisson 残差，却没有消除电子块问题；固定场强则在密度误差
原样存在时让电子块通过。两者的组合并非简单相加，且残差接近门限时可能跨界，
不能声称所有舍入项都已消失。固定输入仅用于归因，没有用于实际闭合求解。

密度干预仅改变上述中间算术，不是把 Fermi–Dirac 及整个密度函数改成高精度；
累加干预仅改变残差存储及累加精度，通量本身仍先按对应分支求值。
这限定了两项阴性对照的结论范围。

## 2. 场强误差进入连续性残差的具体路径

正式 `CoupledDDAssembler::residualImpl()` 先通过
`referencedValueRelativeTo(reference, increment, 0.0)` 生成 double 的
`phinPhysical`，再传入 `transportCellVectorEdgeGradientMagnitudes()`。
`cellScalarGradient()` 随后计算节点间电势差。
平移后参考值约为 −28 V，小量在绝对电势合并时已经舍入；后续差分不能找回它们。
中间相加曾使用 long double 也不足以保证结果，因为物化值最终仍是 double。

诊断保持梯度几何、相邻输运单元面积加权、取模、HFS 模型及接触回退逻辑不变。
只在每个三角形选本地锚点 a，直接构造：

```text
qf_relative(i) = (reference(i) - reference(a))
              + (increment(i) - increment(a))
```

参考差和增量差在 long double 中计算，最后转成供梯度使用的 Real。
因此避免先生成约 −28 V 的绝对 QF 值。残差、局部项诊断、SG 边诊断和 Jacobian
基准场强的对应调用一起切换，便于闭合检查使用一致的基准场强。

残差扰动最大的边提供了具体证据：边 **11992**，节点 **3974→3984**。
相对梯度干预使电子驱动场从 `5094.907243560144` 变为 `5094.907243665755`
（当前内部电场单位 V/cm），仅差 `1.05611e-7`；边迁移率差
`−7.35622e-9`（内部 cm²/(V·s)）。大边通量的微小变化在相邻节点上以相反符号
进入残差，正对应 control 中节点3984的 `+5.37853e-11` 和节点3974的
`−5.35419e-11`。完整边记录使用内部通量单位，不能直接当作 A/μm 电流。

两个参考系间电子残差向量的差异范数从 `1.48840e-10` 降至 `1.13492e-11`
（两个参考系均使用相对梯度）。它没有严格变成零，说明其他舍入项仍存在。

## 3. 原门限闭合、重启及等价性

| 求解对照 | 更新数 | 最终电子残差 | 完整闭合 |
| --- | ---: | ---: | --- |
| 相对 vector，从原参考系合格态 | 0 | 9.224689e-12 | 通过 |
| 相对 vector，从补偿平移态 | 0 | 8.801839e-12 | 通过 |
| 上一行保存后重新读取 | 0 | 8.801839e-12 | 通过 |
| 相对 vector，从上一轮停滞态恢复 | 1 | 3.830474e-12 | 通过 |
| 仅 SG long double，从补偿平移态 | 13 | 3.891929e-11 | 线搜索拒绝 |
| 四项组合，从补偿平移态 | 0 | 8.995769e-12 | 通过 |

相对 vector 的补偿平移态 psi=`4.925889e-8`、hole=`1.62024e-23`，局部违规0、
最大局部ratio=`2.88926e-12`。全局和局部条件均由原求解入口实测，不是只根据
residual probe 的返回状态判断。

按原换参考门限将输出反平移，与原父状态比较：

- 最大电势绝对差 `7.10543e-15 V`，门限1e-7 V。
- 最大密度相对差 `2.13296e-13`，门限1e-6，沿用原 `max(|old|,1)` 分母。
- 最大端电流变化／原最大端电流 `3.63754e-14`，门限1e-8。
- 新状态 KCL ratio `7.66060e-13`，门限1e-8。

保存重读、停滞态恢复及四项组合的等价性也都通过。
这里复核了原 guard 的指标，没有把结果写回正式 frame ledger。
端电流仍由原 `ContactCurrent` 实现计算；正式修复仍需统一并验证后处理路径，
不能以本次 guard 对照替代跨模块回归。

## 4. Jacobian 与线搜索交叉验证

新增18个 JVP 方向：两个参考系，各检查 phin3984、phin4535、psi4535，幅度
1e-7、1e-11、1e-13 V。按 `||Jv-FD||/max(||Jv||,||FD||)` 重算，不加1的分母下限。

在 frame28、phin3984→电子块、1e-13 V 方向，上一轮误差 `7.60432e-3`，
相对 vector 后为 `1.36824e-5`，改善约556倍；1e-7 V 时为 `1.29320e-7`。
phin4535→电子块在1e-11和1e-13 V时分别为 `5.19319e-11`、`1.23247e-8`。
这与场强物化导致小 QF 更新响应失真的证据一致。

psi4535→Poisson 块的小步长误差未被这个干预解决：frame28、1e-13 V
仍为 `5.16435e-3`。密度/Poisson 的舍入敏感性确实存在，但本次被测单点不需要
消除它就能满足原全局门限。初始 psi 已接近5e-8门限，后续偏压的裕量仍须验证。

源码核对表明，启用绝对块门限时，block_filter 要求至少一个**当前仍未达标**的
残差块满足下降条件，并检查各块包络；已经达标的其他块下降不能单独放行。
因此原停滞不能用“总残差或 Poisson 下降”解释为应该接受。
诊断没有修改线搜索、门限或该规则。

## 后续实现边界

单点因果链现在有通过原门限的反事实支持。适合进入正式修复的最小方向是：
保持外部参考平移低位、保持重打包增量、让 vector QF 梯度直接消费参考值/增量。
应以共同的相对梯度算术覆盖残差、Jacobian 场强/反馈、局部诊断和后处理，避免
多个绝对电势重建分支再次产生差异。

本诊断只替换了8处 vector 基准场强调用；Jacobian feedback 内部梯度及其他
潜在重建路径还没有作为正式公共接口重构。本次18个方向和单次更新通过不足以
替代完整 Jacobian、参考系不变性、重启和高场数值回归。

正式落地后，再以原门限验证26.6667 V换参考与后续连续点，最后重跑0–40 V的
Vg=8和双曲线D5评分。Vg=4原完整通过证据和Vg=8原失败证据均保留。

## 复现及证据

从本工作树根运行，诊断脚本的首次运行拒绝覆盖已有 case 目录；重复实验应使用
新的输出目录。`components.py build probes solves validate` 会编译隔离副本，运行
24组残差对照、6组闭合和18个 JVP 方向。两个只读归约脚本可重复运行。

- [干预脚本](../../reference_staging/templates_ldmos_frame_components_20260909/components.py)
- [局部源码 diff](../../reference_staging/templates_ldmos_frame_components_20260909/source.diff)
- [构建命令和源码哈希](../../reference_staging/templates_ldmos_frame_components_20260909/build_manifest.json)
- [残差、闭合、JVP与冻结检查](../../reference_staging/templates_ldmos_frame_components_20260909/analysis.json)
- [状态电流等价性及边证据](../../reference_staging/templates_ldmos_frame_components_20260909/verification.json)
- [归约脚本](../../reference_staging/templates_ldmos_frame_components_20260909/analyze_components.py)
- [等价性核验脚本](../../reference_staging/templates_ldmos_frame_components_20260909/verify_components.py)
