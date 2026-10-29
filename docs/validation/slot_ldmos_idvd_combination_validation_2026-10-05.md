# SLOT IdVd 组合实现与冷启动复验（2026-10-05）

**结论：显式组合与材料分流已实现，针对性回归通过；冷平衡态完整合格，
门压建立四次有界尝试均失败，尚未得到非零门压合格点或 IdVd 曲线。**
最终批次耗时 2,273.63 s（37.89 分钟），因门压严格验证重试耗尽停止，
不是达到预算上限、进程异常或完成原生曲线对比。所有冻结资源哈希复核不变。

## 范围与资格边界

本批执行“补齐组合能力及材料分流，再做定量复验”。前置证据见
[IdVd 原生基线与能力预检](slot_ldmos_idvd_preflight_2026-10-05.md)。
原生主脚本为 `IdVd_des.cmd`：300 K、Vg=5 V、Vd=0–60 V，
不是关态 BVDS。本批没有恢复 BVDS 扫描，也没有修改旧结果或验收门槛。

本批首先实现显式的 **元素盒 Fermi/PhuMob/Lombardi/HFS/雪崩组合**。
它不是对原生离散格式已经等价的声明：这里使用 `mesh_default`、
`transport_edge_geometry=element_box`、cell-material Poisson 系数；
此前 IDL 所用原生外部 AverageBox 边耦合系数未用于这条元素装配路径。
外部几何加 PhuMob 雪崩的配置保护继续保留，不以读入但被覆盖的外部系数宣称等价。

独立输入：

- [IdVd 控制配置](../../configs/regression/slot_ldmos_idvd_fermi_control.json)
- [IdVd 材料配置](../../configs/regression/slot_ldmos_idvd_fermi_materials.json)
- [有界复验驱动](../../scripts/run_slot_ldmos_idvd.py)

网格 13,881 节点、41,155 边、27,275 三角形，沿用已归档 SLOT 网格和节点掺杂。
输入 `unit_scaling` 使用 µm、cm⁻³、cm²/(V·s)，输出电流 A/µm。
HDF5 状态的载流子字段使用 m⁻³。原点为 0 V，不使用 BVDS 的 36 V 偏移。
构建为 Windows UCRT64 Release，HDF5 与 SuiteSparse 后端已编入；
性能记录确认实际调用 UMFPACK。数值复验冻结 runner SHA-256 为
`d137894b21332b59f20afedba08e99daf334d6f0a931a517a54b38401b920c35`。
源码、配置、程序及测试记录见 [实现清单](../../build/s05idvd_impl/implementation_manifest.json)
与 [测试汇总](../../build/s05idvd_impl/test_results.json)。工作树含未提交改动，不能仅用 HEAD 标识程序。

随后仅将新元素源诊断的 `flux_sign` 从未适用的边通量零值标签改为
`not_applicable`，不改 residual、Jacobian 或验收逻辑。最终构建 SHA-256 为
`f0e575001b8c2556123fe1439eea4ba1bb1cece6b5087fdf3d3533c209ba18b2`；
该标签改动后的两个新源测试（142 断言）及 BV probe 五项（202 断言）再次通过。
正在进行的数值复验始终使用前述冻结副本，没有中途更换程序。

## 完成的组合能力

1. **Fermi 输运 Jacobian**：元素 PhuMob 密度链由 Boltzmann 的单位对数响应改为
   Fermi 的 `dn/deta / n`、`dp/deta / p`；保留相邻元素第三顶点及异种载流子列。
   端点广义 SG 偏导仍使用既有中心差分，固定迁移率后再加完整迁移率链，避免漏列或重复计入。
   本实现不能称为“整个 Jacobian 都由符号推导精确计算”。
2. **材料分流**：新增 `mobility.constant_materials`，IdVd 指定 `PolySilicon`。
   多晶硅直接取材料常数 μn=1417、μp=470.5 cm²/(V·s)，不经过 PhuMob/Lombardi/HFS。
   Si 保留完整迁移率模型；SRH/Auger 沿用 Si-only 配置。
3. **显式雪崩离散**：新增 `phumob_element_density_grad_qf`，仅 Si 有源。
   元素顶点积分为
   `M_k [alpha_n mu_n n_k |grad(phi_n)| + alpha_p mu_p p_k |grad(phi_p)|]`，
   包括单位转换、Lombardi 和输运 HFS 场；VO 雪崩使用配置中的 GradQF 驱动场。
   local AD 组合 Fermi 密度响应及解析构成偏导。量子、RefDens 插值、额外源修正、
   非 300 K 等未经本批验证的组合仍拒绝。
4. **诊断一致性**：全局源积分、BV process probe 和活动分支记录使用相同局部源计算。
   修复首次冷启动暴露的诊断接口遗漏，避免诊断回落到缺少完整 QF 向量的旧边接口。

Fermi+Bennett 材料参数已独立固化，但 **材料映射仍未全部通过原生定量验证**。
前置十节点 BGN 检查中九点差值小于 8e-9 eV，高补偿节点 11404 仍有
6.9216e-4 eV 差异；本批没有拟合或掩盖该差异。因此通用 Sentaurus 导入器的
Bennett+Fermi 保护没有解除。

## 回归检查

| 检查 | 结果 |
|---|---|
| `test_element_box_transport` | 17 项、5,223 断言通过 |
| `test_bv_process_probe` | 5 项、202 断言通过 |
| `test_impact_ionization` | 61 项、2,702 断言通过 |
| `test_box_geometry [carrier_transport_profile]` | 4 项、34 断言通过 |
| IdVd 配置约束 | 3 项通过 |
| IdVd 驱动拒绝态隔离/缩步父态绑定 | 1 项通过，mock 控制流，不是仿真证据 |

新增数值检查覆盖全列 Jacobian 中心差分的两个扰动尺度、电子/空穴成对生成、
氧化层/多晶硅无雪崩、Poly 常数迁移率，以及诊断源积分与实际连续方程源积分相等。
另已执行载流子统计 13 项/900 断言、PhuMob 18 项/171 断言、复合 43 项/234 断言。
这些测试不替代 SLOT 全网格曲线资格。

构建遇到本机 GCC 标准头搜索路径重复问题；仅在生成的 `build-release/build.ninja`
移除多余的 UCRT64 根 `-isystem` 项后完成构建。未提交构建产物或改系统工具链。
一次 Python 驱动测试因沙箱默认临时目录权限失败，改用工作树内临时目录后通过。

## 冷启动与成本诊断

所有证据写入新目录 [build/s05idvd_impl](../../build/s05idvd_impl/)。
驱动从全零偏压冷态开始，仅完整器件、有效载流子行、全局闭合及 QF 边界同时通过的
Newton 状态才能成为下一点父态。eps_row=0.001、closure=0.01、
abstol=1e-9、reltol=1e-8 均不变。未使用原生状态热启动或 Boltzmann poly 重置。

- `cold_run`：51.39 s 因诊断路径缺少完整 QF 而退出，0 个合格态；这是已修复的接口错误。
- `cold_run2`：修复后达到 300 s 子进程上限，0 个合格态；此结果只能记为超时，
  不据此判定器件方程无解、门槛过严或线性求解失败。
- `cold_diagnostic/iter0`：从上述冷初始化的 Poisson 种子重放，0 次 Newton 更新，20.79 s。
- `cold_diagnostic/iter1`：相同种子，最多 1 次更新，58.99 s；残差归一值 1→0.579354，
  ψ 块 147.712261→85.577709，更新被接受。因诊断的迭代上限退出，非合格点。
  Jacobian 17.19 s、线性分析/分解/求解合计 0.526 s；连续方程诊断 22.42 s、
  residual 10.51 s、活动分支记录 7.12 s。这些 profiler 项可能嵌套，不能直接全部相加。

`cold_run3` 的冷平衡 **通过**：7 次 Newton 更新，338.36 s，原停止原因为
`abstol`，归一残差 5.75132e-10，载流子超标行 0，最大行比 1.87267e-9，
QF 边界违规 0。全局闭合为 **absolute/absolute** 通过，不是 relative 资格通过。
输出电流 6.61070e-53 A/µm，处于零偏压零电流水平，不能据此判断开态输运精度。
ψ 块原始范数 8.49541e-8 与归一停止残差口径不同，不应直接拿它与 abstol 比较。

## 冷平衡场量对照

已逐项核对两侧全部节点 ID/坐标以及三角形连接/材料一致；状态密度由 m⁻³ 转为 cm⁻³
后与原生同节点对照。以下是逐节点相对差的绝对值分位数，未以统计平均掩盖最大差异。

| 区域/场量 | 95% 分位 | 最大值 | 最大差异节点 |
|---|---:|---:|---:|
| Si 电子密度 | 0.13369% | 2.70259% | 3604 |
| Si 空穴密度 | 0.01795% | 2.71279% | 11404 |
| Poly 电子密度 | 0.01589% | 0.13638% | 12551 |
| Poly 空穴密度 | 0.01588% | 0.13627% | 12551 |

证据与可复算脚本：[cold_density_comparison.json](../../build/s05idvd_impl/cold_density_comparison.json)、
[compare_cold.py](../../build/s05idvd_impl/compare_cold.py)。
11404 同时是前述 BGN 公式差异节点，但其密度差的完整归因尚未做受控分离；
3604 的最大电子差也不能直接归给某个迁移率模型。零偏压对照不替代非零漏压曲线。

## 门压建立中的阻塞

所有门压尝试的漏压均为 0 V，父态均为上述完整合格冷平衡态；失败状态没有成为续算父态。

| 门压 | 墙钟/s | 接受 Newton 更新 | 最后超标行数 | 最大行比 | 结果 |
|---|---:|---:|---:|---:|---|
| 0.1 V | 190.91 | 1 | 1,168 | 约 1 | 全局闭合线搜索拒绝 |
| 0.05 V | 599.34 | 10 | 3,825 | 1.247596 | 载流子行线搜索拒绝 |
| 0.025 V | 494.40 | 9 | 73 | 1.000004 | 载流子行线搜索拒绝 |
| 0.0125 V | 650.45 | 8 | 10,006 | 1.000015 | 载流子行线搜索拒绝 |

0.05 V 末态的归一残差已到 4.45836e-13，但局部行检查没有通过；这说明不能仅用
总残差接受该点。超标行分为 158 条 Si 电子行、3,603 条 Si 空穴行和 64 条 Poly
电子行。最大绝对残差分别约 2.53e-24、5.12e-29、1.57e-14（装配器内部残差单位）。
多晶硅多数载流子相对接触 QF 仍有约 0.470 mV 偏移；0.025 V 该偏移约 0.131 mV。
Si 弱行中有部分达到 scale 下限，但 Poly 超标电子行没有达到下限。
因此本批不能把全部失败归因为数值底限，更不能据此判定 eps_row=0.001 过严。
细分证据见 [gate005_row_groups.json](../../build/s05idvd_impl/gate005_row_groups.json)。

最后一级 0.0125 V 仍未通过，且超标行数不随跨度单调下降。该末态 ψ/电子/空穴
原始块范数为 2.38769e-11 / 3.56563e-15 / 1.79061e-14，全局闭合仍按绝对分支通过，
不能用这些指标覆盖局部行失败。驱动结果为 `completed=true`、`goal_reached=false`、
`stop_reason=gate_strict_failure_after_bounded_retries`、`resources_unchanged=true`。
完整记录见 [summary.json](../../build/s05idvd_impl/cold_run3/summary.json)。
本批只有 1 个合格冷平衡记录；非零门压及漏压曲线合格点均为 0。
未继续启动求解，也不给出开态电流相对误差或原生 IdVd 定量验证通过结论。

下一步有针对性的对照是 **Fermi 兼容的孤立多晶硅偏压初值和 Poisson 自洽建立**：
只把保持统计关系的状态变换作为预测初值，再通过同一完整器件资格；不得套用
Boltzmann 的载流子指数重置或直接接受预测态。该方案尚未在本批验证。
同时保留高补偿节点的 BGN 差异及外部几何兼容性为独立待验证项；
不能将冷平衡通过替代材料映射、离散几何和开态输运的全部验证。
