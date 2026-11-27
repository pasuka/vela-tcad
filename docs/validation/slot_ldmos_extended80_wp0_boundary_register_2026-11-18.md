# SLOT LDMOS 端到端 80 位验证边界登记

日期为 2026-11-18。本文是 [ExtendedPrecision 开发计划](../superpowers/plans/2026-11-18-sd-sim-extendedprecision80-development-plan.md) 的 WP0 交付物。表中的验证方法用于后续工作包。实际运行结果单独记录。

## 首批范围

电势＋浓度 SLOT LDMOS，13,881 节点，0 V Poisson 到首次耦合，再推进到 0.10 V 外电源并保存复读。基线沿用既有 SLOT 配置，包括 IALMob、碰撞电离、增广外接电阻和当前接受口径。外电源电压、器件漏端电压、接触电流和负载线残差分别记录。冻结系统以本批选定矩阵的实际维数和映射为准。

## 精度边界与验证

| 边界 | 入口与消费方 | 验证 |
|---|---|---|
| 脚本、IR 和构建分派 | [Math 解析](../../scripts/vela_sdevice/cmd_parser/sections.py) → [阶段配置](../../scripts/vela_sdevice/deck_emit/__init__.py) → [runner 选择](../../scripts/vela_sdevice/orchestrator/__init__.py) → [C++ 配置检查](../../src/tools/vela_example_runner.cpp)。现有检查围绕编译的中间算术模式。 | 先测既有解析及分派，再用真实 CLI 核验每个阶段的能力和错误退出。检查 `VELA_BIN`、`VELA_EXAMPLE_RUNNER`、`VELA_EXAMPLE_RUNNER_EXTENDED80` 的覆写及错配。 |
| 数值首次解析和输入量化 | [CMD 数值](../../scripts/vela_sdevice/cmd_parser/syntax.py)、[材料派生计算](../../scripts/vela_sdevice/materials/__init__.py)、[网格和掺杂导入](../../scripts/vela_sdevice/mesh_import/__init__.py) → C++ 配置、[物理常数](../../include/vela/core/PhysicalConstants.h)与[单位换算](../../src/core/UnitScalingSystem.cpp)。 | 文本使用在 binary64 合并、在 binary80 可区分的一对值。对外部文件分别记录文件字段的原生精度与升精度后的消费类型。检查实际参与方程的温度、DOS、尺度和接触平衡人口。 |
| 冷态和权威状态 | [Poisson 装配](../../src/equation/PoissonAssembler.cpp) → [DDSolution](../../include/vela/solver/GummelSolver.h) → [耦合状态](../../include/vela/equation/CoupledDDAssembler.h) → [增广状态](../../include/vela/simulation/PseudoArclength.h)。 | 同一次冷启动中对比 Poisson 末态、耦合入口和首次装配实际读取值，含归一化电势、权威人口、参考势及接触未知量。 |
| 本构、输运和导数 | [SLOT 浓度运行路径](../../src/equation/LegacyPopulationRuntime.cpp) → [SG](../../src/discretization/ScharfetterGummel.cpp)、[IALMob](../../src/equation/IalTransport.cpp)、[雪崩](../../src/equation/IalElementAvalanche.cpp)与[局部 AD](../../include/vela/equation/Tri3LocalForwardAD.h)。 | 同态回放 SG 边差、单元 7095 的带符号梯度与偏导、残差和 Jacobian 消费身份。物理模型复现与精度传递分别验收。 |
| 全局装配、线性接口与回退 | [约束消元](../../src/solver/ConstraintElimination.cpp)、[增广系统及选择后端](../../include/vela/simulation/CoupledLoadLine.h) → [线性求解](../../src/solver/LinearSolver.cpp) → [缺陷与修正](../../include/vela/solver/LinearRefinement.h)。 | 冻结两套系统的矩阵、RHS、全局自由度映射和父态。核对行列缩放、消元、方向返回、恢复及失败回退。覆写 `VELA_LINEAR_SOLVER` 时检查实际后端类型和诊断。 |
| 写回、停止与接触 | [浓度试步](../../include/vela/solver/PopulationState.h)、[停止与坐标转换](../../include/vela/simulation/PopulationCoordinateNewton.h)、[电流](../../src/post/ContactCurrent.cpp) → [偏压扫掠](../../src/simulation/DCSweep.cpp)。 | 对比同一方向的 80 位副本与降精度副本，记录实际写回、保正、E 判据、负载线、尝试次数和接受理由。核对失败分支及父态资格。 |
| 独立进程的持久化和交接 | [状态写入](../../src/io/DDSolutionState.cpp) → [HDF5 字段契约](../../src/io/StateArchive.cpp) → [Save/Load、阶段重绑定](../../scripts/vela_sdevice/deck_emit/__init__.py)。当前阶段重绑定会消费 [从输出解析的电压和电流](../../scripts/vela_sdevice/orchestrator/__init__.py)。 | 同一父态对比分段和连续运行。验证无损归档的有效位、尺度和端口实参。测试显式升精度的旧存档，以及精度不匹配和缺字段时的错误。 |

## 基线比较口径

- 首次重放以已有 VM 验证程序、冻结用例和本次生成的输入哈希清单绑定，记录程序的实际哈希。后续新增功能产生的可执行文件独立记录。
- 固定所选后端、线程、物理和报告模式。逐点比较偏压、Newton 次数与停止理由，比较 Poisson 保存态、逐点权威字段、末态和端口记录。源文件路径、日志时间和 HDF5 容器填充不用于数值等价比较。
- 固定系统证据保存 MatrixMarket、RHS 和全局索引、父态、配置、端口读数与对应哈希。记录重放的墙钟时间和峰值内存可用性；在同资源环境测量才用于性能比较。
- VM 原始数据按独立目录保留。运行记录区分未提交、运行中、成功、明确失败和超时。
