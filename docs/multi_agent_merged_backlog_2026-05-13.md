# Multi-Agent Merged Implementation Backlog (2026-05-13)

## 1. 输入来源
1. 架构/算法风险审计（Vela Architecture Analyst）。
2. 测试覆盖缺口审计（Explore）。
3. 基线文档：[docs/agent_handoff_baseline_2026-05-13.md](agent_handoff_baseline_2026-05-13.md)

## 2. 合并优先级规则
1. 先修复“功能不可达/易误用”的阻断项。
2. 再补“收敛可靠性 + 关键测试缺口”。
3. 最后做性能和可观测性增强。

## 3. 合并后的实施 Backlog

| ID | 优先级 | 类型 | 事项 | 关键改动位置 | 测试任务 | 依赖 |
|---|---|---|---|---|---|---|
| B01 | P0 | 功能缺口 | 暴露 Newton 入口（CLI + 可选 DCSweep 分派） | [src/tools/vela_example_runner.cpp](../src/tools/vela_example_runner.cpp#L64), [src/simulation/DCSweep.cpp](../src/simulation/DCSweep.cpp#L173), [include/vela/simulation/DCSweep.h](../include/vela/simulation/DCSweep.h#L1) | 新增 Newton 路径集成测试，确保配置可达 | 无 |
| B02 | P0 | 稳定性 | Gummel 收敛判据增加 abstol 保底，避免三重相对误差锁死 | [include/vela/solver/GummelSolver.h](../include/vela/solver/GummelSolver.h#L30), [src/solver/GummelSolver.cpp](../src/solver/GummelSolver.cpp#L287) | 新增边界收敛测试（高掺杂/强阻尼/maxIter耗尽） | 无 |
| B03 | P0 | 测试缺口 | LineSearch 回溯失败路径专测（alpha 衰减/acceptFunction 拒绝/minDamping 边界） | [src/numerics/LineSearch.cpp](../src/numerics/LineSearch.cpp#L11), [include/vela/numerics/LineSearch.h](../include/vela/numerics/LineSearch.h#L9) | 新文件 tests/test_line_search_backtrack_failure.cpp | 无 |
| B04 | P0 | 测试缺口 | Newton 线搜索失败回退行为专测（返回 acceptedX） | [src/solver/NewtonSolver.cpp](../src/solver/NewtonSolver.cpp#L238), [src/solver/NewtonSolver.cpp](../src/solver/NewtonSolver.cpp#L248) | 扩展 [tests/test_newton_solver.cpp](../tests/test_newton_solver.cpp#L261) | B03 |
| B05 | P1 | 稳定性 | Newton 残差按分块归一化/加权，减轻混合量纲停机偏差 | [src/solver/NewtonSolver.cpp](../src/solver/NewtonSolver.cpp#L219), [src/numerics/ResidualNorm.cpp](../src/numerics/ResidualNorm.cpp#L16), [include/vela/solver/NewtonSolver.h](../include/vela/solver/NewtonSolver.h#L18) | 新增分块残差停机行为测试 | B04 |
| B06 | P1 | 可维护性 | SG 通量实现收敛到单一代码路径（DDAssembler 与 CoupledDDAssembler 一致） | [src/equation/DDAssembler.cpp](../src/equation/DDAssembler.cpp#L143), [src/equation/CoupledDDAssembler.cpp](../src/equation/CoupledDDAssembler.cpp#L173) | 新增一致性对比测试 | B02 |
| B07 | P1 | 稳定性/性能 | Newton warm-start 选项化（避免无条件清零 phin/phip） | [src/solver/NewtonSolver.cpp](../src/solver/NewtonSolver.cpp#L128), [src/solver/NewtonSolver.cpp](../src/solver/NewtonSolver.cpp#L194), [include/vela/solver/NewtonSolver.h](../include/vela/solver/NewtonSolver.h#L18) | 扩展 Newton 扫描迭代数对比测试 | B01 |
| B08 | P1 | 功能缺口 | 温度参数化（300K硬编码改为配置透传） | [src/solver/GummelSolver.cpp](../src/solver/GummelSolver.cpp#L21), [src/solver/NewtonSolver.cpp](../src/solver/NewtonSolver.cpp#L16), [include/vela/solver/GummelSolver.h](../include/vela/solver/GummelSolver.h#L30), [include/vela/solver/NewtonSolver.h](../include/vela/solver/NewtonSolver.h#L18) | 新增温度敏感性测试 | B02 |
| B09 | P1 | 测试缺口 | DCSweep 步长边界专测（growth 后失败回缩、minStep 触发中止） | [src/simulation/DCSweep.cpp](../src/simulation/DCSweep.cpp#L256), [src/simulation/DCSweep.cpp](../src/simulation/DCSweep.cpp#L266) | 扩展 [tests/test_dc_sweep.cpp](../tests/test_dc_sweep.cpp#L217) | 无 |
| B10 | P1 | 测试缺口 | 高掺杂 Gummel 稳定性专测（1e24 级） | [src/solver/GummelSolver.cpp](../src/solver/GummelSolver.cpp#L222), [tests/test_dd_gummel.cpp](../tests/test_dd_gummel.cpp#L175) | 新文件 tests/test_gummel_high_doping.cpp | B02 |
| B11 | P2 | 性能 | DDAssembler 缓存几何量（edgeCells/nodeVolumes/edgeCouplings） | [src/equation/DDAssembler.cpp](../src/equation/DDAssembler.cpp#L69), [src/equation/CoupledDDAssembler.cpp](../src/equation/CoupledDDAssembler.cpp#L54) | 回归测试确保数值等价 | B06 |
| B12 | P2 | 性能 | LinearSolver 稀疏结构复用（analyzePattern 缓存） | [src/solver/LinearSolver.cpp](../src/solver/LinearSolver.cpp#L9) | 压测/性能回归（不改物理结果） | B11 |
| B13 | P2 | 功能缺口 | Bandgap Narrowing 实体模型实现 | [src/physics/BandgapNarrowing.cpp](../src/physics/BandgapNarrowing.cpp#L5), [src/physics/RecombinationModel.cpp](../src/physics/RecombinationModel.cpp#L36) | 新增 BGN 单测 + 回归样例 | B08 |
| B14 | P2 | 可观测性 | LineSearch/Newton 增强诊断历史（可选） | [src/numerics/LineSearch.cpp](../src/numerics/LineSearch.cpp#L11), [src/solver/NewtonSolver.cpp](../src/solver/NewtonSolver.cpp#L269) | 输出历史结构测试 | B03 |

## 4. 串行执行建议（里程碑）

### M0（开工前，半天，基线冻结）
1. 重新配置 Debug 构建并运行当前测试，记录已知失败、耗时和主要样例输出。
2. 验收：形成一份本地基线记录，至少包含构建命令、CTest 结果、`examples/pn_diode/simulation.json` 与 `examples/moscap/simulation.json` 的当前运行状态。

### M1（本周，P0 阻断项）
1. B01, B02, B03, B04。
2. 验收：Newton 可被配置调用；关键失败路径有测试；P0 相关测试全绿；新增配置字段在示例或文档中可被发现。

### M2（下周，稳定性主线）
1. B05, B06, B07, B08, B09, B10。
2. 验收：困难用例收敛率上升或失败原因更明确；新增边界测试全绿；温度参数从配置透传到 Gummel/Newton，且默认行为兼容 300K。

### M3（后续，性能与物理扩展）
1. B11, B12, B13, B14。
2. 验收：性能提升可量化；新增物理模型与文档一致；诊断历史不会改变默认求解结果。

### M4（发布收口，1 周）
1. 梳理 README、docs/examples.md、Python API 文档与示例配置，确认新入口和默认参数一致。
2. 运行全量 CTest、回归脚本和代表样例，输出 release note 草稿。
3. 验收：用户从文档可以复现 Poisson、Gummel DC sweep、Newton DC sweep、Python API 四条路径。

## 5. 建议 agent 分工
1. Vela Architecture Analyst：B01/B02/B05/B06/B07/B08 的设计评审与风险复核。
2. Explore：B03/B04/B09/B10/B14 的测试缺口回归检查。
3. 默认编码 agent：按 backlog 实施代码、补测试、跑构建与回归。

## 6. 每项统一 DoD
1. 代码改动有对应测试。
2. 不引入新编译告警/错误。
3. CTest 全量通过（至少相关分组通过）。
4. 文档与配置字段同步更新。

## 7. 后续开发计划（执行版）

### 7.1 推荐 PR 切片
1. PR-0：基线与测试脚手架。只补测试辅助函数、CTest 标签或文档化基线，不改求解行为。
2. PR-1：Newton 入口可达（B01）。最小化暴露 CLI/DCSweep 分派，保持 Gummel 作为默认路径。
3. PR-2：收敛判据与失败路径测试（B02/B03/B04）。先补 LineSearch 独立测试，再接 Newton 失败回退测试，最后修改 Gummel abstol。
4. PR-3：稳定性主线（B05/B07/B09/B10）。聚焦 Newton 残差、warm-start 和 DCSweep 边界，避免同时改 SG 公式。
5. PR-4：SG 通量去重与一致性（B06）。将 DDAssembler/CoupledDDAssembler 的 SG 计算收敛到同一 helper 或同一调用约定。
6. PR-5：温度参数化（B08）。统一 solver options、simulation config、Python 入口的默认值和透传路径。
7. PR-6：性能优化（B11/B12）。先建立性能基线，再引入缓存；每一步都保留数值等价测试。
8. PR-7：BGN 与诊断历史（B13/B14）。BGN 先做物理模型与单测，诊断历史保持 opt-in。

### 7.2 每个 PR 的进入条件
1. 明确本 PR 只处理哪些 backlog ID，非目标问题只记录不顺手修。
2. 先跑对应现有测试，确认不是在未知红灯上开发。
3. 新增配置项时先确定默认值、单位、兼容行为和文档落点。
4. 改 solver/equation/numerics 时必须列出至少一个最小复现实例或单元测试入口。

### 7.3 测试矩阵
| 变更范围 | 必跑测试 | 追加验证 |
|---|---|---|
| LineSearch/Newton | `ctest --test-dir build --output-on-failure -R "newton|line_search"` | 检查失败回退不会覆盖最后 accepted state |
| Gummel/DDAssembler | `ctest --test-dir build --output-on-failure -R "dd_gummel|sg_flux"` | 高掺杂、强阻尼、maxIter 耗尽场景 |
| DCSweep | `ctest --test-dir build --output-on-failure -R dc_sweep` | growth 后失败回缩、minStep 中止、CSV 行数 |
| 温度/BGN/物理模型 | `ctest --test-dir build --output-on-failure -R "mobility|recombination|dd_gummel"` | 300K 默认兼容，温度变化趋势符合预期 |
| 性能缓存 | `ctest --test-dir build --output-on-failure` | 固定样例输出数值等价，记录耗时对比 |
| Python API | `ctest --test-dir build --output-on-failure -R python_api` | 仅在启用 `VELA_ENABLE_PYTHON=ON` 后执行 |

### 7.4 风险门禁
1. 数值容差门禁：任何收敛判据调整都要同时记录绝对残差、相对残差和迭代次数，避免只靠“测试通过”判断质量。
2. 默认行为门禁：新增 solver option 必须保持旧配置可运行；默认温度为 300K；默认 DCSweep 路径仍为 Gummel，除非配置显式选择 Newton。
3. 物理一致性门禁：SG 通量、温度参数和 BGN 改动不能只看单元值，还要覆盖一个端到端 PN 或 MOSCAP 样例。
4. 性能门禁：缓存类优化不能改变矩阵装配顺序之外的数值语义；若浮点差异不可避免，需要给出容差说明和基线对比。
5. 可观测性门禁：诊断历史必须可选，默认输出保持简洁，避免破坏 CSV/VTK 文件格式兼容性。

### 7.5 交付物清单
1. 代码：每个 backlog ID 对应独立 commit 或清晰 PR 分段。
2. 测试：新增测试文件或测试用例名称写入 PR 描述，说明覆盖的失败模式。
3. 文档：配置 schema、示例 JSON、README/docs/examples.md 与 Python API 说明同步。
4. 结果：至少保留一次全量 CTest 摘要和代表样例运行摘要。
5. 决策记录：对 solver 默认值、容差、BGN 模型公式来源和性能缓存策略做短 ADR 或 PR 备注。

### 7.6 推荐执行顺序
1. 先完成 M0，确认构建工具链和当前测试基线。
2. 按 B03 -> B04 -> B02 -> B01 处理 P0：先锁住失败路径，再改收敛判据和入口行为。
3. 按 B09 -> B10 -> B05 -> B07 -> B06 -> B08 处理 P1：先补边界测试，再做稳定性和接口扩展。
4. 按 B11 -> B12 -> B14 -> B13 处理 P2：先做可验证的性能优化，再做诊断与物理模型扩展。
5. 最后进入 M4，做文档、示例和发布验收。

### 7.7 每周检查点
1. 周初：确认本周目标 backlog ID、依赖是否满足、是否需要拆 PR。
2. 周中：同步测试新增情况和第一个可运行样例结果。
3. 周末：汇总 CTest、样例、性能或收敛数据，并更新剩余风险。
