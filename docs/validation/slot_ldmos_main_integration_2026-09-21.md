# SLOT-LDMOS 合并 main 与 Templates/LDMOS 进展对照

日期：2026-09-21。范围是合并本地 `main`、核对模板 worktree 的已实现能力，
并验证合并兼容性；本次没有将模板的物理参数或通过资格移植到 SLOT BVDS。

## 代码来源与保留措施

- 当前 worktree：`.worktrees/slot-ldmos-bvds-next`，分支
  `codex/slot-ldmos-bvds-next`。
- 合并来源：本地 `main` 的 `60f28b44dfe4a993ac54ca4d455043853dccec74`，
  相对本轮开始的 SLOT HEAD 新增 81 个提交；没有访问或推送远端。
- 查阅对象：`.worktrees/templates-ldmos-phase-a`，HEAD
  `c91b567ed36c95d554f8ca64d42ae52730e2b93b`。
  其 `src/`、`include/` 与上述 main **没有差异**。模板提交已迁移到 main，
  提交 ID 因 cherry-pick 不同，不需要再次逐个移植。
- 合并前的 42 个本地变更保存为 `9d53a6f0`，恢复分支为
  `codex/slot-ldmos-bvds-next-pre-main-20260921`；生成算例输出保持原位。
- `e9444a9` 将旧 `IALMobility.h/.cpp` 改名为 `SlotIALMobility.h/.cpp`，
  避免 Windows 与 main 新增的 `IalMobility.h/.cpp` 发生大小写路径冲突。
  旧 SLOT 参数单独保存为 `slotIalMobility`；既有 `ialmob_field` deck 继续
  选择旧端点实现，低场入口明确命名为 `slot_ialmob`。main 的 `ialmob`
  保持其显式元素几何及参数合同，二者不会隐式互换。

## 已实现或改进的功能

| 功能 | Templates/LDMOS 已完成的工作 | 对当前 SLOT 算例的意义和边界 |
|---|---|---|
| IALMob 与完整 Jacobian | 独立 IALMob 公式、晶向/界面距离、顶点迁移率及 Measure 平均；每个三角形的 9 个电势列通过载流子密度、法向场及高场驱动链式反馈，包含跨载流子项；电热接口另有温度导数 | 已有比 SLOT 端点平均更完整的空间对齐实现。SLOT 原有密度/表面偏导仍保留。迁移需要为 SLOT 网格准备并验证元素几何与原始参数，不能只改模型名字 |
| 高场与接触场导数 | 修复向量 GradQF 的非端点列，先平均向量再求模；接触电场回退包含相邻顶点反馈；普通无表面模型的高场反馈使用显式链式导数 | 与本轮 SLOT 端点修复存在重叠，合并消除了重复导数叠加，保留表面场和 n/p 变化。是否缩小 SLOT 电流差异仍需同模型对照 |
| 网格与控制体积 | 独立 AverageBox 输运系数、材料局部 Poisson 电荷体积；电热接口也接入来源一致的材料加权 Poisson 边系数及 Measure 体积 | 可用于量化 SLOT 的混合 Voronoi/assembled-Delaunay 与原生离散差异；模板网格权重不能直接用于 SLOT |
| 统计与产生复合 | 高精度 Fermi 函数、解析导数、逆映射及 BGN 联动；材料来源核对；Auger 单位、H/N0 增强及生成项开关统一到源项和导数 | SLOT 仍需按其原始 Sentaurus deck 核对参数和开关。模板曾发现 Auger 系数单位相差 10^12，不能据此断言 SLOT 同样存在该错误 |
| 混合接触与低电流精度 | p+/n+ 共用金属接触采用节点局部极性；准费米参考系重定心、亚 ULP 增量保存、严格绝对块和局部行门限下的推进改进 | 适合 SLOT 的微小漏电流诊断；原接触物理、局部行和全局闭合门限需要保留，不能用更宽停止条件解释为精度提升 |
| 连续扫描与恢复 | 预测、检查点、真实暂停恢复、普通 DC worker、准备缓存与主 Newton 跨偏压符号分析复用 | 可以减少重复准备成本。每个新点仍重新数值分解，结构、输入或后端变化会使缓存失效；不保证解决 SLOT 在 2.85 V 附近的推进困难 |
| 稀疏后端与性能 | 默认优先 UMFPACK，未编译时用 SparseLU；另有 STRUMPACK、MUMPS、SuperLU_MT；显式选择失败不会静默换后端 | 已随 main 合入。本次实际构建支持五种后端；旧 SLOT 资格来自 SparseLU/COLAMD，合并对照显式固定旧后端以便比较 |
| 自热 | 独立四方程生产入口、温度相关硅物性和 IALMob、有限空穴交换边界、热收支审计及 R11 显式高场求值 | 普通 `dc_sweep` 仍是等温；SLOT BVDS 的自热或击穿后支路没有因此通过验证 |
| 可复现算例 | D5 中性网格、掺杂、权重、材料、双栅压零压种子和参考曲线已入库，统一脚本可运行 | 新环境无需 Sentaurus VM 即可检查 D5。D5 关闭 IALMob、雪崩、量子和自热，不能作为 SLOT 完整物理合同 |

原生元素 IALMob 的普通 `CoupledDDAssembler` 入口当前明确限制为 **300 K、
量子和雪崩关闭**（见 `refreshIalMobility`）。因此，本次虽然接入新模型源码，
仍不能宣称“严格原生元素 IALMob + 自洽雪崩 + 外电阻 BVDS”已实现并完成验证。

## 模板验证进展与 SLOT 现状

模板最终 R11 等温复核：G3 为 Vd=0.1 V、Vg=0–5 V 的 31 点；
D5/D4 各在 Vg=4/8 V 下由零漏压连续推进至 40 V，共 155 个精确点，
全部通过各自原门限。D5 最大电流误差为 1.592417%/1.625083%，
D4（IALMob）为 1.891761%/1.801547%；四条 Id–Vd 曲线均无外层失败回退。
这些是导通态曲线，网格为 10,241 节点、19,782 三角形，电流单位 A/µm。

SLOT 为关断态反向漏压 DC，其他端子为 0 V，13,881 节点、27,275 三角形；
当前配置采用 Boltzmann 统计，不能把模板的 Fermi 精度提升直接归因于 SLOT 的差异。
截至本次检查，原冻结 Release 证据显示：

| SLOT 资格 | 实际完成状态 |
|---|---|
| IALMob 无雪崩连续低压 | `build/ir2` 完成 0→2 V，52 次尝试；2 V 漏电流 8.8425679161e-16 A/µm，最大局部行比 3.25694e-7 |
| 无雪崩独立检查点 | `in2/in3/in5/in10` 的 2/3/5/10 V 通过；10 V 为 2.8667550308e-15 A/µm |
| 自洽雪崩独立检查点 | `ia2r/ia3/ia5/ia10` 的 2/3/5/10 V 通过；10 V 为 2.9714412803e-15 A/µm，最大局部行比 5.73045e-4 |
| 2→10 V 连续推进 | `ih11` 最后合格点 2.85 V，任务曾中断；独立 10 V 点不能替代连续扫描资格 |
| 外电阻 | `irx2r/summary.json` 为 passed，外部 2 V、R=1e12 Ω·µm 的负载线复验通过；不是击穿后分支验证 |

模板性能资格也有范围：五后端统一单线程完整曲线批次完成第二轮后由用户暂停，
20/30 条曲线、620 点、10 组双栅压联合验收通过，第三轮未启动。
最新默认值收尾版本另有全量 843 项测试和双栅压前八点资格，
不是本次合并后二进制的自动继承资格。

## 本次合并验证

使用 MSYS2 UCRT64 GNU 16.2.0、Release；CMake 确认 HDF5/TDR、SPQR、
UMFPACK、STRUMPACK、MUMPS、SuperLU_MT 与 METIS 可用，Eigen 排序为 COLAMD。
Sentaurus 重启转换及 SLOT IALMob 工具 7 项 Python 回归通过；SLOT 和入库 D5
相关 CTest 7/7 通过，D5 九项输入依赖哈希预检通过。

新 Release 显式使用 SparseLU，读取旧 `in2` 的合格 2 V 状态、保持原物理和
门限，零次 Newton 更新通过完整初始验收。漏电流
`8.8425574366608805e-16 A/µm` 与冻结结果逐值相同；全局电子/空穴闭合比
分别为 `1.01955e-12` / `1.02075e-12`。另行生成连续性方程项并独立审计，
每载流子 12,306 个有效行、1,575 个无输运行；电子/空穴最大局部比为
`3.96911e-12` / `7.56483e-6`，均小于原 `eps_row=0.001`，无超限行。
该检查证明已存状态的兼容性，不认领重新完成低压连续扫描或高压推进。
证据：`build/main_merge_20260921/slot_2v_final/` 的 `iv.csv` 与
`local_audit/summary.json`。最终运行器 SHA-256 为
`c8a6f05f82b75dfc23f195412a66d97a3e6b6562d93ba8157ec58a3946f19b0e`。

全量首轮 874/875 通过；失败的停止判据回归发现自动合并保留了一个旧初始
残差下限分支，绕过 main 新增的接触准费米势检查。移除该重复分支，最终
全局闭合修正也执行相同的绝对块及接触门槛，未放宽任何容差。受影响目标
重新构建后，**完整 Release CTest 875/875 通过，138.29 s**。
最终构建检查为 `ninja: no work to do`。构建与两轮测试记录分别保留于
`build/merge_main_20260921_build.log`、`build/merge_main_20260921_rebuild.log`、
`build/merge_main_20260921_ctest.log` 和 `build/merge_main_20260921_ctest_final.log`。

## 后续采用顺序

建议先以合并后的 Release 保留旧 SLOT 模型完成同偏压兼容性检查，再独立准备
SLOT 的元素 IALMob 几何及参数，完成无雪崩迁移率/完整 Jacobian/2 V 状态对照。
随后分离 AverageBox、材料局部体积、统计及产生复合的影响，形成误差预算；
元素 IALMob 与雪崩的耦合及其全部导数需要单独实现、审计和回归，再推进连续
高压与外电阻击穿支路。本次不自动启动这项物理迁移。

## 依据

- [模板当前状态](templates_ldmos_current_status.md)及[性能收尾](templates_ldmos_performance_closeout_2026-09-20.md)
- [R11 完整等温复核](templates_ldmos_r11_isothermal_2026-09-17.md)
- [生产电热复现入口](templates_ldmos_production_reproduction.md)
- [原生 Poisson 与准备复用](templates_ldmos_native_poisson_and_preparation_2026-09-13.md)
- [Auger 生成语义对照](templates_ldmos_generation_alignment_2026-09-13.md)
- [SLOT IALMob 历史对照](slot_ldmos_ialmob_alignment_2026-09-07.md)
- [入库 D5 输入](../../reference_tcad/templates_ldmos_sentaurus2022/validation_d5/README.md)

本地运行证据位于忽略的 `build/`；上述历史结果均保留其原模型、程序、网格和
门限边界。本轮重新读取了 SLOT 各 `qualified_progress_snapshot.json`、
低压 summary 及外电阻 summary，没有重跑原 Sentaurus 对照。
