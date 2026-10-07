# BVmethods 从 0 V 的纯 QS 扫描诊断，2026-10-06

原固定栅压、固定电阻、从 0 V 的全程纯 QS 路径**仍未修复**。
本轮完成 13 组原生对照；增大外电压步长、二阶预测、分组求解和电流控制都没有越过
约 `1.409265 mA/um` 的原支路终点。同一外电压 13 kV 下取得了两个不同的稳态解，
支路获取方式确实影响结果。该证据不证明动态稳定性，也不是原支路不存在后续解的严格证明。

## 当前契约与环境

当前仓库 `reference_tcad/bvmethods_sentaurus2018/source/external_resistor_sdevice.cmd`
和 `models.par` 为本轮基准：300 K；漏极固定 `R=1e7 ohm*um`；源、栅、衬底外电压均为 0 V，
栅 Barrier 为 −0.55 eV；`OldSlotboom`、掺杂/高场/Enormal 迁移率、
`SRH(DopingDep)`、`Band2Band(E2)`、`Avalanche(Eparallel)` 和 Fermi 保持原声明。
原停止值 `1.443e-3 A/um`，漏极外电压 Goal 为 100000 V。
Math 中 `Transient=BE` 只是离散选项；没有 Solve Transient 的控制仍是稳态扫描。

原生环境为现有 Sentaurus VM（RHEL 7.9），`T-2022.03-SP2`，默认 Blocked/Super 求解，
具名 monolithic 对照为 Method=Super。统一原网格 2600 节点、5020 单元，原生报告非 Delaunay 单元为零。
网格 SHA-256：`7cd3968f5dd0e3cd63ef62ed0fcced465dd33bea465167f0be8139c16025902e`。
本轮没有更换原生环境，也没有声称完成 Ubuntu 或 Vela 数值重验。

历史 `v1/source/pp5_des.cmd` 不是上述模型的逐字原件：它声明
`Band2Band(Model=NonLocalPath)`、`Avalanche(GradQuasiFermi)`，停止值为 `1.3e-3 A/um`。
旧最终采样电流不能替代旧声明的停止阈值。README 已纠正此前“全部保留原模型和目标”的过宽表述；
本轮没有通过改回历史模型或降低当前停止值来消除失败。
历史 cmd SHA-256：`30f025ed061f0c113cc4e1caa7c00f56003a6ac9ad01a4199e014b08602108a5`。

## 对照结果

| 控制目录 | 变化/用途 | 最后接受外电压 V | 内电压 V | 结果 |
| --- | --- | ---: | ---: | --- |
| monolithic_jump | 先到 13 kV，再尝试 14.45 kV；Super、无外推 | 14099.811214 | 7.160258 | MinStep；85 个状态负载/KCL检查通过，停止值未达到 |
| wide_jump | 先到 1 kV，再尝试 14.5 kV；Super | 14099.811219 | 7.160256 | MinStep |
| wide_blocked | 同上，Blocked/Super | 14099.811224 | 7.160254 | MinStep |
| quadratic | 二阶外推，80 位扩展精度 | 14099.811222 | 7.160255 | MinStep |
| hole_predictor | 仅空穴二阶外推 | 14099.811213 | 7.160258 | MinStep |
| gummel_single | EP80，带 BreakOnFailure 的分组求解 | 14099.811219 | 7.160256 | MinStep |
| current_trace | 从零起源的 13 kV 保存态切换电流控制 | 14099.811225 | 7.160251 | MinStep；电流 `0.00140926509745 A/um`，不作为原电压路径资格 |
| frozen_reseed | 从同一纯 QS 13 kV 态重载，冻结雪崩导数 | 13000.002766 | 7.134220 | MinStep；仅重启诊断 |
| avaldens_gradqf | 替换雪崩电流近似的单独数值对照 | 1579.598828 | 6.258297 | MinStep；不属于原离散近似 |
| segmented_jump | 分段放大步长并线搜索 | — | — | 有界诊断人工终止，PLT 不完整；不作为完成轨迹 |
| gummel_jump | 首个多子对象 QS 写法 | — | — | 原生语法拒绝，无数值资格 |
| gummel_sequence | 修正为 Plugin 容器但未 BreakOnFailure | 14500 | 2.941402 | 无效“成功”：巨大电流抵消，被独立校验拒绝 |
| reverse_branch | v4 高电流完整保存态反向稳态扫到 13 kV | 13000 | 6.637744 | 11 点负载/KCL/稳态检查通过；含 v4 种子，不是从零获取 |

`gummel_sequence` 打印 Good Bye 和停止标记，但漏极总电流为 `8.689e129 A/um`，
电子/空穴电流约 `±8.501e146 A/um`，同时存在无法求解的线性系统。
不能用结束字符串、进程 exit=0 或停止标记单独认定通过。MinStep 控制中的部分进程也返回 0。
所有原始日志、完整和不完整 PLT 均保留；没有修补中止 PLT 来冒充完整输出。

## 13 kV 下的两个稳态解

| 获取方向 | 内电压 V | 漏极 A/um | 源极 A/um | 衬底 A/um |
| --- | ---: | ---: | ---: | ---: |
| 从零纯 QS 正向 | 7.1342022476712 | 0.00129928657977523 | −3.12302891714165e-6 | −0.00129616355085809 |
| 从 v4 保存态反向 | 6.63774406432838 | 0.00129933622559357 | −0.000273914332499853 | −0.00102542189309371 |

两者的负载残差分别为 `2.91e-11 V` 和 `2.73e-11 V`，终端 KCL 相对残差均小于 `5.4e-15`，
位移电流及电子/空穴分量闭合通过既定门限。内电压相差 `0.496458 V`，源/衬底电流分配显著不同。
因此改变数值初值可能选中不同稳态根；这些结果支持多稳态/支路切换诊断，
仍不能唯一归因于寄生双极、特定离散缺陷或负载线折点。

用户在 SLOT LDMOS BVDS 中观察到“外步长增大、内步长缩短”，与
`ΔVouter = ΔVinner + R ΔI` 一致：电流陡增时，外电压变化主要落在电阻上。
本轮已测试跨越较大外电压间隔，仍回到同一终点；电流控制也在相近电流失败。
所以该 NMOS 问题不能仅归结为外电压步长不足；本轮证据不能替代 SLOT LDMOS 自身网格和日志的验证。

## 校验工具、归档与剩余工作

新增 `scripts/validate_wp14_zero_qs.py` 与注册回归 `wp14_zero_qs`。
它检查源模型/电阻/参数/停止值、从零起点、全程固定其他电极偏置、单调外电压、
所有接受点的负载关系/KCL/位移电流/电流分量闭合、原停止值和原生日志状态。
沿用负载 `1e-6 V`、KCL `1e-5`、位移比例 `1e-6`、分量闭合 `1e-8`；
低电流 KCL 使用明确 `1e-12 A/um` 绝对下限，未放宽任何高电流门限。
辅助栅扫描、重载、电流切换和瞬态种子不能取得原固定栅压纯电压 QS 的通过资格。

Windows UCRT64 Release 的四项相关注册 CTest **4/4 通过**：
`sdevice_reference_artifacts`、`wp14_remaining_sources`、`wp14_external_recovery`、`wp14_zero_qs`。
新增纯 QS 回归含 9 个用例，覆盖原生假成功、坏中间态、辅助栅偏置及失败分类。
首次 sources 检查因默认 TEMP 写权限失败，改用当前工作树 TEMP 后复验通过。
本轮未改 C++ 求解器、未运行全仓 CTest，也未扩大 v4 已有资格。

VM 归档：`/home/tcad/vela_reference_artifacts/bvmethods_zero_qs_controls/v1`。
13 组控制和元数据共 **149 文件**（含清单），下载镜像成员集合、大小和 SHA-256 逐项一致。
tar SHA-256：`c0952f7c0d5e654538ef79c72b9c6a56349405e3920a92eea0d93d34d9d47497`；
inventory SHA-256：`36e177e70b21b519b08ae02f8ec6976f53cfb82eb08a2a98b6ea963a7ce4328c`。
结构化结果见[诊断摘要](bvmethods_zero_qs_diagnosis_2026-10-06.json)。原生原始输出留在 VM 与忽略的
`build-release/qs_zero_20261006/`，当前默认 v4 manifest 未修改。

## 已授权辅助栅偏置纯 QS 候选：执行失败

用户明确授权后，已完成候选：从 0 V 到漏极外电压 13 kV，声明栅压 0→0.2→0 V，
随后继续原漏极扫描。另完成同一路径的大步长对照，只把两段栅扫描的 InitialStep/MaxStep
从 0.05 改为 1，并设置局部 `Extrapolate(MinStep=1)`；失败时仍允许缩步。
模型、材料、网格、原电阻和停止值与本报告基准一致，原生环境及默认 Blocked/Super 后端相同。
两个源文件分别在
`reference_tcad/bvmethods_sentaurus2018/zero_qs_controls/gate_qs_seed_v1/` 与 `gate_qs_jump_v1/`。

| 控制 | 三阶段接受状态数 | 实际最高栅外电压 V | 最后漏极外电压 V | 最后内电压 V | 最后漏极电流 A/um |
| --- | --- | ---: | ---: | ---: | ---: |
| 原候选 gate_qs_seed | 71 / 87 / 45 | 0.0966954759163124 | 14099.8112218879 | 7.16025526464873 | 0.00140926509666232 |
| 大步长 gate_qs_jump | 71 / 75 / 45 | 0.0966954759244808 | 14099.8112218879 | 7.16025526464873 | 0.00140926509666232 |

**两组均失败**：升栅未到 0.2 V，最终漏极电流也未达到原 `1.443e-3 A/um`。
每组升栅段和最终漏极段各报一次 MinStep。原生在升栅失败后继续执行降栅、Save 及后续漏极扫描，
最终打印 Good Bye 且进程退出码为 0；输出文件或 Save 的存在不能单独证明各阶段目标完成。
两组降栅恢复到 0 V 后，内电压均回到约 `7.1342022476712 V`，
源极电流约 `−3.12302891714e-6 A/um`，未取得 v4 高源极电流分支。
大栅步长对照的失败支持该路径无法仅靠这项步长修改完成；它不是没有其他纯 QS 初始化方式的证明。

新增独立 `scripts/validate_wp14_gate_qs.py`：检查源/参数/网格、四段目标次序、
无重载和瞬态、三个输出阶段衔接、升栅最高值和恢复值、保持 13 kV 的栅扫描，
以及全部接受点的物理闭合和最终原停止条件。两组 **394 个状态**的物理闭合均通过，
最大负载残差 `1.465e-9 V`，最大高电流 KCL 相对残差 `7.42e-15`，
位移电流及电流分量闭合通过既定门限；阶段目标/日志成功/停止值检查仍明确失败。
原固定栅压校验器保持原范围，辅助路径没有获得替代资格。

Windows UCRT64 Release 五项相关注册 CTest **5/5 通过**（上述四项加 `wp14_gate_qs`）。
补充“升栅失败后仍继续并达到后续停止标记”的回归后，`wp14_gate_qs` **1/1 注册测试、8 个用例**通过；
三个状态校验器的最终 Python 回归 **21/21**通过。没有改 C++、重跑全仓或降低数值门限。

两组完整原生数据单独归档为
`/home/tcad/vela_reference_artifacts/bvmethods_zero_qs_controls/v2`：**36 文件**（含清单），
下载镜像的完整成员集合、大小及 SHA-256 全部一致。
tar SHA-256：`8059f2d3957ae937ef294b621c2ffb3fba8167f93ede8c08d055d0f2fc640f7b`；
inventory SHA-256：`bfd581366d2c9573e4fde648dac34eb8b6513cc4dbe59cc25fe375bb14c153dc`。
完整检查项与指纹见[辅助栅对照摘要](bvmethods_auxiliary_gate_qs_2026-10-06.json)。
当前 v4 manifest 未修改；原固定栅压从零纯 QS 和本次辅助栅路径均保留未通过状态。
