# 合并版本的 HFS Codespaces 基线

日期：2026-09-26。工作分支：`codex/simplemos-sdevice-validation`。

本轮按用户要求，将合并后的版本放到既有 Codespace
`vela-tcad-compute-69r6rj7pvvjc54g9` 重新完成 HFS 基线。大型状态、节点场和逐边数据保留在云端；验收、双初始化、原生场比较和汇总在云端执行，仅回传摘要、失败记录及证据哈希。

## 版本和范围

- HEAD：`06d0acb425b052f0df424341ab86a5994ed45bf5`，包括本地 main `04e07713a270a259064393bcb3b183294aa19f72` 的合并。
- 从合并后的工作区快照出发，包含此前已验证但尚未提交的 SimpleMOS 生产修改；不是只构建 HEAD。最初的 206 个源码/CMake 文件与上一轮合并验证快照一致。实际预检随后发现 HDF5 适配缺口，本轮增加下述状态 I/O 修复，并在新的快照中重编译、验证；没有把最初二进制与修复后结果混记。
- 首先运行 n19/n23 × Vd=0.05/1 V × Vg=0/0.2/0.8/1 V，共 16 点、32 个独立初始化状态。全部通过后，才进入四条 Vg=0–1 V、步长 0.02 V 的完整曲线，共 204 点、408 状态。
- 原生 204 点直接复用经历史哈希校验的 Sentaurus HFS 曲线。旧 Windows HFS 已完成的状态不计入这次合并版本结果。
- 本轮不包括 SRH 体积修复或原始 n17–n24、0–2.5 V 的 816 点矩阵。

## 配置、几何和单位

300 K、Boltzmann、OldSlotboom、PhuMob + Enormal + HighFieldSaturation、掺杂依赖 SRH。显式保留 element-box PhuMob、element-distance-gradient Enormal、element-vertex-partial-layer HFS 及其场导数。

几何保留逐材料介电系数、Delaunay box 转移、硅侧输运与独立 signed Si Poisson 电荷体积；**SRH 仍为 all_cell 体积**。网格 n19 为 1480 节点，n23 为 1482 节点；各有 942 个 Si 节点、907 个自由 Si 节点，每态核验 1814 个载流子行。

`unit_scaling` 保持原设置。状态 CSV 密度是 m⁻³；原生场密度及此配置下 all_row 的旧 `_m3` 浓度列为内部 cm⁻³；SRH 率为 cm⁻³ s⁻¹；端口电流为 A/μm。物理场比较使用共同的正 barycentric Si 体积，生产 SRH 重建则使用原 all_cell 体积，二者不可混用。

云端 GCC 16 Release/Ninja 编译成功。HDF5/TDR、HDF5 状态、UMFPACK、SPQR、METIS 被实际检测到。由于合并后默认线性后端会优先选择已安装的 UMFPACK，本轮显式设置 `VELA_LINEAR_SOLVER=sparselu`，保持原 HFS 的 Eigen SparseLU/COLAMD 后端；每进程线性/BLAS 线程为 1，初始并发 4。

最初二进制 SHA256（仅用于失败的格式预检，不是最终 HFS 求解器）：
`2a137261153e036941e7c1ee649b9644b57ddf99c0327555faf42dfdd1eccd2c`。

## 验收与续算协议

- 自由载流子行 `|R| / max(sum|edge flux|, |SRH|, |impact|) ≤ 1e-6`。
- 源闭合容差 `1e-6`，源激活底限 `1e-10`；KCL/|Id| 和端口提取相对差各 `≤1e-8`。
- 双初始化三势最大差 `≤1e-6 V`，载流子密度相对差 `≤1e-4`，Id 相对差 `≤1e-6`。
- 每个目标最多允许一次同偏置保存态重载；首次失败始终保留。只有合格状态可以传播到下一偏置。
- Id 相对 Sentaurus 的百分比单独汇总；数值资格通过不等于所有点的电流差都小于 1%。
- 全曲线数值资格通过后，再进行 Decimal100 密度/SRH 重构，以及电势、两种准费米势、密度和 SRH 的原生场比较。密度重构门槛 `1e-12`、生产 SRH 积分 L1 重构门槛 `1e-7`。

## 工作流预检与已保留的失败

新增可移植执行器的 9 项失败路径/资格测试在 Windows 和 Codespaces 都通过。一个历史合格态的 6 项验收结果与旧执行器逐项一致；8 组历史双初始化的场差计算也一致。

第一批 32 个启动全部在读材料时退出，尚未进入 Newton：新版 `MaterialDatabase` 拒绝历史自定义标签 `vela.transportmodels.sentaurus2022.materials.v1`。这属于输入格式兼容问题，不是 32 个数值不收敛态。完整记录保存在第一批云端目录和本地 `build/hfs_merged_20260926/attempt1_summary/`。

第二批在新的 `replay_v2` 目录中，将输入材料副本转为当前支持的无 schema 旧格式：只移除根对象的 `schema`，其余 JSON 内容完全相等，四种材料的每个参数、数值、单位和继承规则不变。每个副本均有迁移说明及前后哈希。没有修改生产读取器，也没有覆盖历史材料文件。

第二批也在求解前停止：生产重启接口只接受 `.h5`。进一步审计发现，原 HDF5 DD 适配器还会遗漏 `packedState`、`packedLow`、势标度及 split 网格标识，单纯转换文件后缀不能保持此前资格。

为此，本轮补齐 `StateArchive`/`DDSolutionState` 与 Python codec 的可选三块坐标字段；低位余量逐位保存，包括 subnormal 尾数。完整性、物理投影、元数据、网格和潜势原点检查保持严格。迁移工具单独解析 CSV 中重复的标度/标识列，生产入口仍只使用 HDF5。四个历史诊断入口补齐 HDF5 独立网格/原点上下文。

执行器的种子、DC 输出、同偏置重载和后处理探针全部使用 HDF5。供分析使用的 `state.csv` 是校验 HDF5 独立网格身份后导出的诊断表，不作运行时回退。

验证结果：

- 220/220 个初始状态完成迁移，全部数值字段在读取后逐位相等。
- C++ 状态合同 10 用例、333 断言通过，包含旧 CSV 导入、三块余量和失效写入保留旧检查点。
- HDF5 修复后完整 Codespaces CTest 为 967/983 通过，无新增失败；16 项均属于此前已复现的历史证据失败。
- 补齐诊断上下文后，`state_archive_contract`、`state_archive_python_interop`、`state_archive_dc_worker` 三项复验全部通过。
- 真实四工况的普通/packed 共 8 个种子全部通过 C++/Python 双向逐位比较，生产端口提取最大相对差 `6.4393e-15`，通过 `1e-8` 原门槛。

最终批次的 **32/32 控制状态、16/16 双初始化已全部通过**，无数值失败、无重载；已进入完整曲线。2026-09-26 13:29:30（北京时间）的云端快照为 **7/408 曲线状态通过**。这只是启动后进度，完整基线尚未完成。

最终求解器 SHA256 为
`1a612009119e7187716e9bbfd883638c7f798c11fc825fe4e97a45d1ee23cd3a`；
输入清单 SHA256 为
`cc606900b3d72a96155d306c9c71968857b0565e343d391c377bbc51fba5ab67`。
驱动继续完成 408 态曲线，再做物理场分析。n23、Vd=1 V、Vg=0 的控制 Id 差为 **−1.0504991046%**，与旧 HFS 基线一致；本轮没有处理已知 SRH 体积差。

合并版本 n23 控制点的 `(Id_Vela / Id_Sentaurus − 1) × 100%`：

| Vg (V) | Vd=0.05 V | Vd=1 V |
|---|---:|---:|
| 0 | −0.025931291% | −1.050499105% |
| 0.2 | +0.003439335% | −0.036689876% |
| 0.8 | +0.005774106% | +0.006626328% |
| 1.0 | +0.005197263% | +0.005098975% |

16 组双初始化最大差：ψ 为 `2.2e-16 V`，电子/空穴准费米势为 `2.3710e-9 / 2.0817e-9 V`，密度相对差 `9.1713e-8`，Id 相对差 `6.5967e-10`。全部通过冻结门槛。这些是离散控制点，不能代替完整曲线验收。

## 证据与操作位置

- 持久根目录：`/workspaces/simplemos-hfs-merged-20260926`。
- 当前批次：`/workspaces/simplemos-hfs-merged-20260926/replay_v4`；前面三批的失败/预检记录保留。
- 冻结身份：当前批次 `run_seal.json`，输入 `inputs/manifest.json`、`inputs/source_hashes.json`。
- 执行日志：当前批次 `run.log`；结束码：`run.exit`；完成后仅打包 `summary.tgz`。
- 大文件：`controls/`、`curves/`、`field_nodes/` 留在云端。
- 本地准备与精简证据：`build/hfs_merged_20260926/`（忽略目录）。
- 已取回的控制摘要：`build/hfs_merged_20260926/control_summary/summary/`；压缩包 SHA256 为 `967ae2c5b56306c616fd41b9b0254000bbeeda1a00c151d5de26d867b4e91f2f`，下载后已验证。
- 本地 `build/hfs_merged_20260926/watch_and_fetch.ps1` 正保持实时日志连接；任务结束后只下载并校验最终 `summary.tgz`，写入 `retrieval_status.json`。该连接不自动重启中断任务，也不修改 Codespace 的机型或超时设置。Codespace 停止或本地连接中断仍需要依据已保存的检查点恢复，不能把后台运行当作完成保证。

实现：[输入准备脚本](../../scripts/prepare_simplemos_hfs_cloud_20260926.py)；[云端执行及比较脚本](../../scripts/simplemos_hfs_cloud_20260926.py)；[校核测试](../../tests/regression/test_simplemos_hfs_cloud.py)。

本轮已针对新增 I/O 修复重跑完整 CTest，随后对诊断上下文修改进行相关复验。此前合并测试的历史证据失败继续保留，不因此宣称全套测试全绿。最终物理基线是否完成，必须读取当前批次的 `summary/control_summary.json`、`summary/progress.json`；不能用旧 HFS 结果替代。
