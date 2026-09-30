# SimpleMOS packed 检查点一致性修复与剩余停滞定位

日期：2026-09-10；分支 `codex/simplemos-sdevice-validation`。

## 结论

检查点保存/加载一致性已修复：原失败态的 4446 个内部坐标及 4446 个残差分量逐项完全复现。两个零迭代控制的原始残差块一致，32 次首次 DC 的物理状态和求解状态也与修复前完全一致。

但原一次重载协议下，本版为 **26/32 状态、10/16 双初始化通过**，此前为 31/32 和 15/16。准确保存坐标移除了旧重载的舍入扰动，六个原生初始化失败态均不能借重载继续前进。这不是物理电流发生大幅改变，也不应将旧版本的合格点拼接进本版验收。

独立重放六个失败态的 78 个线搜索候选：空穴块均下降，但原加权总平方范数均增加。五个状态的电势候选增量全部在 double 加法中丢失；另一个控制还出现电子残差回升。本轮解决了检查点一致性，**尚未解决小于坐标舍入间隔的 Newton 电势更新**。完整 PhuMob 曲线及 Enormal/HFS 仍未放行。

## 配置与门槛

沿用 n19/n23 原匹配网格（1480/1482 节点），原生 box 修改、Si 输运几何、signed Si 三项 Poisson 电荷体积、逐单元介电系数。300 K、Boltzmann、OldSlotboom/匹配 ni、plain PhuMob `element_box_phumob`、掺杂相关 SRH；Enormal/HFS 关闭，SRH 体积和物理常数不变。

实际 UCRT64 Release、Eigen SparseLU、四次线性修正，显式 `poisson_residual_precision: binary128`、`stable_merit_comparison: true`、`exact_dirichlet_updates: true`。电流 A/μm、电势 V、密度 m⁻³；packed 势坐标为电势除以 potentialScale 的无量纲值。三项选项默认仍保持原路径。

门槛保持全部自由 Si 载流子逐行 1e-6、KCL/Id 及端口一致性 1e-8，双初始化势差 1e-6 V、密度相对差 1e-4、Id 相对差 1e-6。仍只有失败后一次同偏压重载，未追加重载或放宽门槛。[冻结输入与协议](../../reference_tcad/simplemos_sentaurus2022/phumob_packed_restart_20260910/cohort_v2/validation_contract.json)。

## 根因隔离与实现

前置输出型诊断将保存坐标按自由/接触节点及 ψ/电子/空穴分块往返。最终失败态中自由 ψ 有 185 个坐标变化，电子 90 个、空穴 73 个，均为 1 ULP；接触坐标没有变化。

仅将自由 ψ 往返就使 Poisson 范数由 6.06299775284e-11 变为 1.09733661698e-10，残差差范数为 8.72661046867e-11；电子、空穴各自往返的残差差范数仅 3.75416e-26、9.98157e-26。因此排除了接触更新作为此次保存态跳变的来源。

`DDSolution` 新增可选原始 packed 三块坐标及缩放值。binary128 Poisson 求解结果导出四个可选 CSV 列：`packed_psi`、`packed_electron_qf_increment`、`packed_hole_qf_increment`、`packed_potential_scale_V`。只在坐标、物理场、准费米参考与增量自洽时写入；不完整或损坏的 CSV 元数据拒绝加载；内存物理场已修改时不输出失效提示。保留旧 CSV 和次正规数归零行为。

重启缩放一致时精确保留自由 ψ；对应准费米参考相同才保留载流子坐标，不同参考仍走原稳定转换。显式接触继续由边界条件设置。此次不改变残差公式、Jacobian、线搜索权重、状态 double 存储或生产默认。详见 [配置说明](../config_schema.md)。

两个零迭代对照中 packed 文件均复现保存前原始残差块；去掉新列的 legacy 对照仍出现旧跳变。合格控制重新初始化时自动块归一化系数不同，因此应比较原始残差块，不能把两个不同归一化下的总数直接视为表示误差。[往返对照](../../reference_tcad/simplemos_sentaurus2022/phumob_packed_restart_20260910/roundtrip_v2.csv)、[4446 坐标及残差检查](../../reference_tcad/simplemos_sentaurus2022/phumob_packed_restart_20260910/vector_summary.json)。

## 十六点与电流结果

32 次首次 DC 的物理 CSV 字段和 iterations/failure_reason/final_residual/exit_code 全部与修复前一致；四个新增元数据列另存，不参与该物理字段身份比较。[首次求解身份表](../../reference_tcad/simplemos_sentaurus2022/phumob_packed_restart_20260910/first_attempt_identity.csv)。

本版共 38 次尝试，其中 12 次失败保留；最终 26/32 状态、10/16 双初始化通过。六个失败全部来自原生初始化，16 个 Vela 初始化全部通过。双初始化 Id 最大相对差 6.01740879e-14，但不能代替行及场资格。[逐点比较](../../reference_tcad/simplemos_sentaurus2022/phumob_packed_restart_20260910/cohort_v2/comparison.csv)、[全部尝试](../../reference_tcad/simplemos_sentaurus2022/phumob_packed_restart_20260910/cohort_v2/attempts.csv)。

| NWell | Vd / V | Vg / V | 最差空穴节点 | 行残差比 | 双初始化最大 φp 差 / μV | 密度相对差 |
|---|---:|---:|---:|---:|---:|---:|
| n19 | 0.05 | 0.0 | 801 | 0.00035296665 | 4.76282 | 0.00018421716 |
| n19 | 0.05 | 0.2 | 801 | 0.00035297731 | 4.763 | 0.00018422391 |
| n23 | 0.05 | 0.2 | 802 | 0.00031106492 | 4.23057 | 0.00016363243 |
| n23 | 1.0 | 0.0 | 800 | 0.023412702 | 928.982 | 0.035296644 |
| n23 | 1.0 | 0.2 | 983 | 0.00071616201 | 18.8031 | 0.00072707149 |
| n19 | 0.05 | 0.8 | 801 | 0.00012041875 | 1.62486 | 6.2850335e-05 |

n23 合格 Vela 初始化路径的电流误差 `100*(Id_Vela/Id_Sentaurus-1)`：

| Vg / V | Vd=0.05 V | Vd=1 V | 双初始化：低 / 高 Vd |
|---|---:|---:|---|
| 0 | -0.01819672% | -0.80813524% | 通过 / 未通过 |
| 0.2 | +0.00358452% | -0.01632131% | 未通过 / 未通过 |
| 0.8 | +0.00569362% | +0.00652926% | 通过 / 通过 |
| 1 | +0.00507750% | +0.00510983% | 通过 / 通过 |

## 剩余停滞的独立证据

只增加输出的 Newton 副本重放六个失败重载配置，状态 CSV 和求解状态 6/6 与生产一致；每个状态记录 13 个候选，总计 78 个。90 位和 140 位 Decimal 对实际残差 double 值独立计算原块权重与缩放下的平方范数差，符号全部一致；没有下降候选被拒绝。原始与截断更新逐项相同，更新上限没有削去本次方向。[重放身份](../../reference_tcad/simplemos_sentaurus2022/phumob_packed_restart_20260910/stalls/identity.csv)、[全部候选贡献](../../reference_tcad/simplemos_sentaurus2022/phumob_packed_restart_20260910/stalls/trial_merit.csv)、[最差行更新](../../reference_tcad/simplemos_sentaurus2022/phumob_packed_restart_20260910/stalls/summary.csv)。

五个低栅压失败态的最大无量纲 ψ 更新为 1.76415e-15～3.54842e-15，而所有候选的 ψ 坐标变化均为零；准费米势可以更新并减少空穴残差，却引起 Poisson 电荷变化，缺少相应的可表示 ψ 更新来补偿。它们的 65 个候选均为 Poisson 块增加，抵消载流子块下降。

n19、Vd=0.05、Vg=0.8 的最大无量纲 ψ 更新为 4.70546e-15；首个候选有 132 个 ψ 坐标改变，小阻尼候选则为零。13 个候选中 9 个 Poisson 块增加；其余 4 个虽 Poisson 下降，但电子块增加超过该改善。空穴块 13/13 下降。该反例不能简单归结为所有状态的 Poisson 残差都上升。

因此下一步应隔离验证能保留微小 ψ 更新的坐标表示，并同步检查 Poisson 电荷、输运势差、Jacobian 与 CSV 往返一致性；Vg=0.8 控制还需保留电子通量精度核验。不能通过改变全局归一化、降低 Poisson 权重或新增重载次数把此次失败改记为通过。

## 验证和边界

- Release 全目标构建完成。两项新增 Catch2 测试共 18 个断言通过：CSV 精确往返/失效元数据，以及自由 ψ 零迭代重启。
- 完整 CTest **781/797**；16 项失败名称与前一轮完全一致，为历史证据身份或缺失文件检查，非全绿。未修改旧基线以消除失败。
- 16 个合格 Vela 状态的正式 Jv 分块 **336/336** 通过，最大相对差 6.54302096e-07，门槛 1e-4。该检查不覆盖弱交叉块的独立相对资格。[Jv 汇总](../../reference_tcad/simplemos_sentaurus2022/phumob_packed_restart_20260910/cohort_v2/jvp_summary.json)。
- 64 个全局 SRH 分量仍未激活原 1e-10 源下限，无新增相对源资格。前一轮 20/20 局部源等结果属于前一版实现，本轮未重跑该批校准或新增 Sentaurus 仿真。
- 初次构建的新测试曾误用显式构造类型的赋值，修正后重建通过；首版往返验证误把不同自动归一化总范数作为同态门槛，失败日志和脚本保留，v2 改为逐原始残差块精确比较。没有改动正式求解接受门槛。

原证据及修改前源码/二进制存档均保留，完整产物由 [完成清单](../../reference_tcad/simplemos_sentaurus2022/phumob_packed_restart_20260910/completion_evidence.json) 索引。本轮未提交或推送；下一阶段处理上述更新表示限制，达到统一 16/16 后再开展完整曲线和后续物理模型恢复。
