# SLOT IdVd：Fermi 多晶硅偏压预测与 Poisson 自洽验证（2026-10-05）

**结论：0.1 V 同初值两臂均通过，组合路径在 0.2/0.3 V 继续通过原完整资格。**
本批获得 3 个不同门压的合格验证点，另有 0.1 V 对照态；不重复计点。
全部求解已结束，两个串行批次累计 449.92 s（约 7.50 分钟），低于 1,500 s 上限；
冻结输入与父态哈希复核均未改变。未取得非零漏压的 IdVd 曲线。

## 任务与边界

本批承接[组合实现及冷启动复验](slot_ldmos_idvd_combination_validation_2026-10-05.md)，
验证此前建议的 Fermi 兼容多晶硅初值与 Poisson 自洽建立。
原来 Vg=0.1/0.05/0.025/0.0125 V 从冷父态直接启动均失败。本批采用同一个
完整合格冷父态，同一冻结 runner、物理模型、网格及器件验收门槛。
只检验 Vd=0 下的门压建立；不是 IdVd 漏压全路径或 BVDS 验证。

配置为 300 K、Fermi/Bennett-Wilson、Si PhuMob/Lombardi/HFS 和 Si-only 雪崩，
Poly 常数迁移率；13,881 节点、41,155 边、27,275 Tri3 单元。
元素盒输运、mesh_default 几何；不宣称与原生外部 AverageBox 离散等价。
电压为物理 V（origin=0），密度文件单位 m⁻³，电流 A/µm。
Windows UCRT64 Release，沿用 HDF5/SuiteSparse 构建及 UMFPACK 后端。

原完整门槛不变：abstol=1e-9、reltol=1e-8、eps_row=0.001、全局闭合 tolerance=0.01，
同时检查全部有效载流子行、QF 边界和完整 Newton 收敛原因。
全局闭合在本批零漏压条件下走原有 absolute 分支，不能称为 relative 资格通过。

## 预测方法

新增显式初值工具 [slot_fermi_poly_predictor.py](../../scripts/slot_fermi_poly_predictor.py)。
对连接同一 ohmic gate、与其他输运材料不共享节点的孤立 Poly 分量，执行

\[
 \psi' = \psi+\Delta V_g,\qquad
 \varphi_n'=\varphi_n+\Delta V_g,\qquad
 \varphi_p'=\varphi_p+\Delta V_g.
\]

材料、掺杂、温度保持不变，Fermi 占据依赖的势差不变，保存的 n/p 不作指数重置。
参考量同步平移，QF 增量保留。非 Poly 节点保持原值。
这是初值预测，不是局部平移后的 Poisson 方程自动满足；与周围氧化层/Si 的电场
变化仍必须求解。工具检查 mesh/material/doping 哈希、父态角色、接触连通性、
源作用范围和 Fermi 配置，并将输出明确标为 `qualified=false`。
目前只允许每次不超过 0.1 V 的有界诊断，未接入默认生产续算。

实际预测改变 350 个 Poly 节点，两种密度逐值不变；最大 `psi-QF` 舍入变化为
5.55112e-17 V。两个回归测试通过，覆盖正/负偏压变换的统计势差、未选节点与密度
保持，以及非有限输入/未经支持的打包坐标拒绝；测试不代替器件资格。

## 同初值路径对照：Vg=0.1 V

两臂使用同一个预测态，目标接触偏压及最终门槛相同。

| 路径 | 更新次数 | 墙钟 | 最终超标行 | 最大行比 | 原完整资格 |
|---|---:|---:|---:|---:|---|
| A：预测→完整 Newton | 4 | 139.39 s | 0 | 1.57461e-9 | 通过，abstol |
| B：预测→Poisson→完整 Newton | Poisson 4；完整 Newton 0 | 85.79+17.25=103.04 s | 0 | 2.98489e-9 | 通过，initial_abstol |

B 的 Poisson 使用 `preserve_initial`，固定预测态 QF，只修复电势。
Poisson 输出本身不计完整合格点；后续完整 Newton 独立执行原始终态检查，
零更新是因为输入已满足原全部门槛，不是绕过载流子行检查。
两臂 QF 边界违规均为 0。与旧未预测 .1 V 尝试相比，这支持“Poly 初值不一致是
本点启动障碍”的判断；不能据此断言以后所有非收敛都具有相同原因。

两臂最终场量高度一致：Si/Poly 最大电势差均为 1.11022e-16 V，
电子/空穴密度最大相对差不超过 7.27e-15；Poly QF 逐值一致。
漏电流均处于数值零水平，绝对差 3.15e-53 A/µm，**不作相对电流精度评价**。
组合路径本次耗时较低，但只有一次样本，不能推广成普遍性能结论。
同初值两条算法路径对照也不等同于不同物理初值的分支独立性证明。

Poisson 修复并非空操作：相对预测态，Si 最大电势改变 61.03 mV，
Poly 最大改变 86.98 mV；载流子密度随 Fermi 统计自洽更新。
这与“在固定电势下重置 QF 并直接保留旧密度”的操作不同。

## 相邻点及成本

在本批 1,500 s 总求解预算的剩余部分，从 B 的完整合格 .1 V 状态串行补验 .2/.3 V，
每点使用同一预测→Poisson→完整资格路径。Poisson 单次上限 240 s，
相邻点完整 Newton 单次上限 180 s，失败即停；预测态和 Poisson 态不计点。
结果如下：

| Vg / V | Poisson 墙钟/s | 完整验收墙钟/s | 载流子超标行 | 最大行比 | 完整验收 |
|---|---:|---:|---:|---:|---|
| 0.2 | 83.97 | 16.77 | 0 | 4.31118e-9 | initial_abstol，0 更新 |
| 0.3 | 85.90 | 19.55 | 0 | 4.85098e-9 | initial_abstol，0 更新 |

两点 QF 边界违规均为 0，全局闭合 absolute 分支通过，漏电流仍为数值零。
相邻点批次记录 `two_adjacent_points_completed`，并非超时或失败退出。
完整链为冷平衡→.1→.2→.3 V，只有各点完整 Newton 验收后的状态才用于下一次预测。

## 证据与剩余限制

- [同初值对照与资源哈希](../../build/s05idvd_init/summary.json)
- [两路径场量对照](../../build/s05idvd_init/same_seed_comparison.json)
- [场量变化](../../build/s05idvd_init/field_changes.json)
- [相邻点记录](../../build/s05idvd_init/adjacent/summary.json)
- [可复算分析](../../build/s05idvd_init/analyze.py)
- [预测工具回归](../../tests/regression/test_slot_fermi_poly_predictor.py)

冻结 runner SHA-256：`d137894b21332b59f20afedba08e99daf334d6f0a931a517a54b38401b920c35`；
冷父态 SHA-256：`88ace85993926bb84af75469c3d0ef4c8aaaa2a5808b4ce1287559b7d2680424`。
本批没有修改 C++ 求解器、模型参数或旧证据，没有重建求解器或恢复 BVDS。
前批暴露的高补偿 BGN 差异、原生几何兼容性仍未消除。
即使门压建立通过，也仍需门压全路径及 Vd>0 的 IdVd 定量复验。

下一步可把本规则作为显式初始化选项接入有界 IdVd 驱动，先建立 Vg=5 V 全门压路径，
再从合格态逐级建立非零漏压并与原生同偏压点比较。进入非零漏压后，不能假定
固定 QF 的 Poisson 初始化仍能替代完整输运迭代；全部原资格检查继续执行。
