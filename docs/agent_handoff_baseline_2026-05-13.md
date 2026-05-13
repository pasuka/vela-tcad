# Vela TCAD Baseline Handoff (2026-05-13)

## 1. 目标
本文件用于给后续 agents 提供统一的事实基线，覆盖：
1. 当前代码架构。
2. 已实现功能清单。
3. 测试用例清单。
4. 后续详细开发方案（分阶段、分工、交付物、验收标准）。

## 2. 架构清单

### 2.1 构建与入口
1. 核心静态库：vela_core，汇总全部求解与物理模块。
证据：[CMakeLists.txt](CMakeLists.txt#L30)
2. 命令行入口：vela_example_runner。
证据：[CMakeLists.txt](CMakeLists.txt#L104), [src/tools/vela_example_runner.cpp](src/tools/vela_example_runner.cpp#L64)
3. 可选 Python 绑定：VELA_ENABLE_PYTHON + pybind11 模块 _core。
证据：[CMakeLists.txt](CMakeLists.txt#L25), [CMakeLists.txt](CMakeLists.txt#L78), [bindings/pyvela.cpp](bindings/pyvela.cpp#L120)

### 2.2 源码分层
1. core：尺度系统与常量，数值尺度转换。
目录：[src/core](src/core)
2. mesh：网格拓扑、边构建、box 几何。
目录：[src/mesh](src/mesh)
3. material：材料参数数据库与覆盖加载。
目录：[src/material](src/material)
4. physics：掺杂、载流子统计、迁移率、复合、BGN 接口。
目录：[src/physics](src/physics)
5. discretization：Bernoulli 与 Scharfetter-Gummel 通量基础。
目录：[src/discretization](src/discretization)
6. equation：Poisson、DD、CoupledDD 装配。
目录：[src/equation](src/equation)
7. solver：线性求解、Gummel、Newton。
目录：[src/solver](src/solver)
8. numerics：残差范数与回溯线搜索。
目录：[src/numerics](src/numerics)
9. simulation：Poisson 任务流与 DC 扫描任务流。
目录：[src/simulation](src/simulation)
10. io/post/tools：Mesh/VTK/CSV、接触电流后处理、CLI 工具。
目录：[src/io](src/io), [src/post](src/post), [src/tools](src/tools)

### 2.3 典型算法路径
1. Poisson-only：配置解析 -> 网格/材料/掺杂 -> Poisson 装配 -> 线性解 -> VTK。
证据：[src/simulation/PoissonSimulation.cpp](src/simulation/PoissonSimulation.cpp#L58), [src/equation/PoissonAssembler.cpp](src/equation/PoissonAssembler.cpp#L95)
2. DD + Gummel：每个偏置点执行 Poisson + 电子连续 + 空穴连续，检查三变量收敛。
证据：[src/solver/GummelSolver.cpp](src/solver/GummelSolver.cpp#L222), [src/solver/GummelSolver.cpp](src/solver/GummelSolver.cpp#L287)
3. Coupled Newton：组装残差/Jacobian，线性步 + 回溯线搜索，双阈值停机。
证据：[src/solver/NewtonSolver.cpp](src/solver/NewtonSolver.cpp#L219), [src/solver/NewtonSolver.cpp](src/solver/NewtonSolver.cpp#L238), [src/solver/NewtonSolver.cpp](src/solver/NewtonSolver.cpp#L278)
4. DC Sweep：自适应步长、失败回退、重试与 CSV/VTK 输出。
证据：[src/simulation/DCSweep.cpp](src/simulation/DCSweep.cpp#L236), [src/simulation/DCSweep.cpp](src/simulation/DCSweep.cpp#L256), [src/simulation/DCSweep.cpp](src/simulation/DCSweep.cpp#L274)

## 3. 功能清单（当前实现状态）

| 功能 | 状态 | 关键证据 |
|---|---|---|
| Poisson 装配与求解 | implemented | [src/equation/PoissonAssembler.cpp](src/equation/PoissonAssembler.cpp#L95), [src/solver/LinearSolver.cpp](src/solver/LinearSolver.cpp#L7) |
| 固定电荷/界面片电荷 | implemented | [src/simulation/PoissonSimulation.cpp](src/simulation/PoissonSimulation.cpp#L107), [src/simulation/PoissonSimulation.cpp](src/simulation/PoissonSimulation.cpp#L133) |
| 漂移扩散装配（Poisson+载流子、n/p 连续） | implemented | [src/equation/DDAssembler.cpp](src/equation/DDAssembler.cpp#L63), [src/equation/DDAssembler.cpp](src/equation/DDAssembler.cpp#L120), [src/equation/DDAssembler.cpp](src/equation/DDAssembler.cpp#L196) |
| Bernoulli 函数 | implemented | [src/discretization/Bernoulli.cpp](src/discretization/Bernoulli.cpp#L6) |
| Scharfetter-Gummel 通量 | implemented | [src/discretization/ScharfetterGummel.cpp](src/discretization/ScharfetterGummel.cpp#L6) |
| Gummel 非线性迭代 | implemented | [src/solver/GummelSolver.cpp](src/solver/GummelSolver.cpp#L306) |
| Newton 全耦合求解 | partially implemented | 内核已实现 [src/solver/NewtonSolver.cpp](src/solver/NewtonSolver.cpp#L149)，CLI 尚未作为独立 simulation_type 暴露 [src/tools/vela_example_runner.cpp](src/tools/vela_example_runner.cpp#L64) |
| 自适应 DC 扫描 | implemented | [src/simulation/DCSweep.cpp](src/simulation/DCSweep.cpp#L129), [src/simulation/DCSweep.cpp](src/simulation/DCSweep.cpp#L236) |
| 材料数据库与外部覆盖 | implemented | [src/material/MaterialDatabase.cpp](src/material/MaterialDatabase.cpp#L62), [src/material/MaterialDatabase.cpp](src/material/MaterialDatabase.cpp#L98) |
| 区域掺杂映射 | implemented | [src/physics/DopingModel.cpp](src/physics/DopingModel.cpp#L39) |
| 迁移率模型（常数/CT） | implemented | [src/physics/MobilityModel.cpp](src/physics/MobilityModel.cpp#L69) |
| 复合模型（SRH/Auger） | implemented | [src/physics/RecombinationModel.cpp](src/physics/RecombinationModel.cpp#L36), [src/physics/RecombinationModel.cpp](src/physics/RecombinationModel.cpp#L111) |
| Mesh JSON 读取与校验 | implemented | [src/io/MeshReader.cpp](src/io/MeshReader.cpp#L90) |
| VTK/CSV 输出 | implemented | [src/io/VTKWriter.cpp](src/io/VTKWriter.cpp#L11), [src/io/CSVWriter.cpp](src/io/CSVWriter.cpp#L14) |
| 接触电流后处理 | implemented | [src/post/ContactCurrent.cpp](src/post/ContactCurrent.cpp#L22) |
| Python API（load_mesh/run_poisson/run_dc_sweep/write_vtk） | implemented | [bindings/pyvela.cpp](bindings/pyvela.cpp#L120), [bindings/pyvela.cpp](bindings/pyvela.cpp#L151), [bindings/pyvela.cpp](bindings/pyvela.cpp#L153) |
| Bandgap Narrowing 实体模型 | not found | 目前为零增量接口 [src/physics/BandgapNarrowing.cpp](src/physics/BandgapNarrowing.cpp#L5) |

## 4. 测试用例清单

### 4.1 C++ 单元测试目标（Catch2）
注册于 [CMakeLists.txt](CMakeLists.txt#L117) 到 [CMakeLists.txt](CMakeLists.txt#L159)。

1. test_scaling: [tests/test_scaling.cpp](tests/test_scaling.cpp)
2. test_mesh: [tests/test_mesh.cpp](tests/test_mesh.cpp)
3. test_poisson: [tests/test_poisson.cpp](tests/test_poisson.cpp)
4. test_box_geometry: [tests/test_box_geometry.cpp](tests/test_box_geometry.cpp)
5. test_bernoulli: [tests/test_bernoulli.cpp](tests/test_bernoulli.cpp)
6. test_sg_flux: [tests/test_sg_flux.cpp](tests/test_sg_flux.cpp)
7. test_dd_gummel: [tests/test_dd_gummel.cpp](tests/test_dd_gummel.cpp)
8. test_newton_solver: [tests/test_newton_solver.cpp](tests/test_newton_solver.cpp)
9. test_mobility: [tests/test_mobility.cpp](tests/test_mobility.cpp)
10. test_recombination: [tests/test_recombination.cpp](tests/test_recombination.cpp)
11. test_dc_sweep: [tests/test_dc_sweep.cpp](tests/test_dc_sweep.cpp)

### 4.2 关键测试主题与代表用例
1. Bernoulli 稳定性。
代表：[tests/test_bernoulli.cpp](tests/test_bernoulli.cpp#L9), [tests/test_bernoulli.cpp](tests/test_bernoulli.cpp#L55)
2. SG 通量守恒与大偏置有限性。
代表：[tests/test_sg_flux.cpp](tests/test_sg_flux.cpp#L84), [tests/test_sg_flux.cpp](tests/test_sg_flux.cpp#L111)
3. MeshReader 校验与错误分支。
代表：[tests/test_mesh.cpp](tests/test_mesh.cpp#L200), [tests/test_mesh.cpp](tests/test_mesh.cpp#L393)
4. Poisson 装配、边界、固定电荷/片电荷、材料覆盖。
代表：[tests/test_poisson.cpp](tests/test_poisson.cpp#L246), [tests/test_poisson.cpp](tests/test_poisson.cpp#L474), [tests/test_poisson.cpp](tests/test_poisson.cpp#L568)
5. Gummel 收敛健壮性与输出有效性。
代表：[tests/test_dd_gummel.cpp](tests/test_dd_gummel.cpp#L81), [tests/test_dd_gummel.cpp](tests/test_dd_gummel.cpp#L120), [tests/test_dd_gummel.cpp](tests/test_dd_gummel.cpp#L175)
6. Newton 收敛、Jacobian 对齐、数值稳定。
代表：[tests/test_newton_solver.cpp](tests/test_newton_solver.cpp#L94), [tests/test_newton_solver.cpp](tests/test_newton_solver.cpp#L220)
7. DC Sweep 自适应步长、失败重试、终点无过冲。
代表：[tests/test_dc_sweep.cpp](tests/test_dc_sweep.cpp#L150), [tests/test_dc_sweep.cpp](tests/test_dc_sweep.cpp#L217), [tests/test_dc_sweep.cpp](tests/test_dc_sweep.cpp#L249)
8. 迁移率与复合模型选择及耦合运行。
代表：[tests/test_mobility.cpp](tests/test_mobility.cpp#L50), [tests/test_mobility.cpp](tests/test_mobility.cpp#L86), [tests/test_recombination.cpp](tests/test_recombination.cpp#L11)

### 4.3 Python 与回归测试
1. Python API：CTest 条目 python_api。
证据：[CMakeLists.txt](CMakeLists.txt#L91), [tests/python/test_python_api.py](tests/python/test_python_api.py#L69)
2. 工程回归：CTest 条目 regression。
证据：[CMakeLists.txt](CMakeLists.txt#L174), [scripts/run_regression.py](scripts/run_regression.py#L1), [tests/regression/README.md](tests/regression/README.md)
3. ASCII 约束测试。
证据：[CMakeLists.txt](CMakeLists.txt#L113)

## 5. 给其它 Agents 的分析输入包

### 5.1 建议直接提供的上下文
1. 本文档：[docs/agent_handoff_baseline_2026-05-13.md](docs/agent_handoff_baseline_2026-05-13.md)
2. 构建与测试定义：[CMakeLists.txt](CMakeLists.txt)
3. 求解器核心：[src/solver/GummelSolver.cpp](src/solver/GummelSolver.cpp), [src/solver/NewtonSolver.cpp](src/solver/NewtonSolver.cpp), [src/numerics/LineSearch.cpp](src/numerics/LineSearch.cpp)
4. 关键测试：[tests/test_dd_gummel.cpp](tests/test_dd_gummel.cpp), [tests/test_newton_solver.cpp](tests/test_newton_solver.cpp), [tests/test_dc_sweep.cpp](tests/test_dc_sweep.cpp)

### 5.2 可直接复制给子 agent 的任务模板
1. 架构风险审计任务：
请基于 docs/agent_handoff_baseline_2026-05-13.md，对 solver/equation/numerics 三层做接口一致性与可维护性风险审计，输出高/中/低风险项、证据位置、最小改动建议。
2. 算法稳定性任务：
请聚焦 Gummel/Newton/LineSearch，分析收敛判据、步长策略、失败回退机制在高偏置和高掺杂下的潜在失效模式，并给出可执行的实验矩阵。
3. 测试缺口任务：
请对 tests 目录做覆盖映射，输出未覆盖行为清单，按优先级给出新增测试计划（文件名、测试名、断言点）。

## 6. 详细开发方案（建议）

### 阶段 A：基线冻结与可观测性增强（1 周）
1. 目标：固化当前行为，避免后续改动引入不透明回归。
2. 任务：
- 为 Newton history 增加更明确统计出口（残差、阻尼、步长）。
- 为 DCSweep 失败路径补充更细粒度 reason 字段（线性失败/线搜索失败/重试耗尽）。
3. 交付物：
- 新增或扩展测试，覆盖失败路径与诊断输出。
- 一份基线报告（残差曲线与 IV 结果摘要）。
4. 验收：
- 相关测试全通过。
- 对同一输入可重复得到一致收敛结论。

### 阶段 B：求解器稳定性与物理一致性强化（2-3 周）
1. 目标：降低高偏置/高掺杂场景不收敛风险。
2. 任务：
- 引入可配置的收敛判据组合（相对 + 绝对 + 分量门限）。
- 线搜索增加策略开关（纯 Armijo、宽松下降、固定阻尼兜底）。
- 梳理 SG helper 与 DD 装配中的公式一致性，补一致性测试。
3. 交付物：
- 配置项文档与默认参数建议。
- 稳定性回归测试组。
4. 验收：
- 相比基线，在困难算例上收敛率提升或迭代次数降低。
- 无新增 NaN/Inf 或负载流子异常。

### 阶段 C：功能边界补全（2 周）
1. 目标：补齐已识别的能力缺口。
2. 任务：
- 明确 Newton 在 CLI 的能力边界（文档先行，必要时接入 simulation_type）。
- 实装 Bandgap Narrowing 物理模型并校准测试。
3. 交付物：
- 功能文档更新。
- 新增 BGN 单元与集成测试。
4. 验收：
- BGN 相关测试通过且结果趋势符合预期。
- 用户入口能力与文档一致，无歧义。

### 阶段 D：多前端一致性与回归自动化（1-2 周）
1. 目标：确保 C++ CLI 与 Python API 行为一致。
2. 任务：
- 对 Poisson/DC sweep 关键输出建立跨前端一致性检查。
- 扩展 regression 数据集，覆盖正向/反向扫描与失败重试场景。
3. 交付物：
- 一致性测试脚本。
- 增强版 regression 清单。
4. 验收：
- 前端输出在容差内一致。
- regression 对关键路径全绿。

## 7. 协作分工建议（面向多 Agents）
1. Vela Architecture Analyst：负责阶段 A/B 的结构和算法审计报告。
2. Explore：负责阶段 B/C 的跨目录证据检索与差异核对。
3. 默认编码 agent：负责阶段实现、补测试、修复编译与回归。

## 8. 当前已知高优先级风险
1. Newton 内核已实现但入口边界不直观，可能造成使用误解。
证据：[src/tools/vela_example_runner.cpp](src/tools/vela_example_runner.cpp#L64)
2. Bandgap Narrowing 仍为零实现，限制高掺杂物理真实性。
证据：[src/physics/BandgapNarrowing.cpp](src/physics/BandgapNarrowing.cpp#L5)
3. 线搜索行为有实现但缺少针对回溯过程的独立细粒度测试。
证据：[src/numerics/LineSearch.cpp](src/numerics/LineSearch.cpp#L11), [tests/test_newton_solver.cpp](tests/test_newton_solver.cpp)
