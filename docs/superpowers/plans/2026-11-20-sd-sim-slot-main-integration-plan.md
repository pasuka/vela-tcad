# SLOT 与 main 的 sd-sim 统一合并方案

## 基线与目标

分析基线为本地 `main@3d56db35` 和 `codex/slot-ldmos-bvds-next@55e0fbd9`，共同祖先为 `a535991c`。从共同祖先起两边共有 57 个路径都被修改。按 Git 实际合并算法只读预检，有 26 个冲突文件，涉及配置、发射器、装配、Newton、DC 扫描和测试。合并对象为已提交版本，本地 main 工作树上的开发中改动由现有工作树继续管理。

面向用户的器件仿真统一从 `sd-sim <device.cmd>` 启动。由 `File` 段解析相对路径，经同一套 PAR 物理参数和 DF-ISE GRD/DAT 网格、掺杂导入，生成可检查的阶段 JSON，执行 Solve 阶段，最终发布结果、日志和运行清单。C++ runner 是该流程的内部执行器。器件类别不决定代码分派，已解析的物理模型、离散配置、方程和精度要求决定装配与执行策略。同一个模型或数值算法在同一套实现中分别实例化为 `double` 和 `long double`。

## 两边实现的取舍

| 层 | 当前差异 | 合并方向 |
|---|---|---|
| CMD 与输出 | main 有共同的 deck 发射、`Extrapolate` 映射、日志级别、结果元数据与 Windows Workbench/SDS 集成。SLOT 增加原始数值字面量、掺杂判定、宽状态身份及基于参考 JSON 的 IALMob 绑定。 | 使用 main 的入口、阶段编排与输出元数据框架。把 SLOT 的 PAR 映射、DF-ISE 导入和精度检查纳入同一流程。当前参考 JSON 用于过渡期配对，合格 CMD 最终直接生成配置。 |
| 模型与离散 | main 已共享 Tri3 梯度、电流重建、配对 SG，并增加 Masetti/Fermi/雪崩的可配置能力。SLOT 增加 IALMob/HFS、直接浓度偏导和稳定的边驱动，以及 80 位几何、本构和装配。 | 从 main 的 Tri3/SG 公共接口出发，把 SLOT 的稳定算式和标量能力合入相应算子。统一残差、电流及 Jacobian 对操作数的消费，保留各自已验证的物理模型、积分方式和几何权重。 |
| Newton 与接受 | main 的浓度路径由 `PopulationCoordinateNewton` 至 `CoupledLoadLineNewton`。SLOT 80 位由 `ExtendedInitialStage.inl` 调用另一套固定偏压与外电阻校正器。 | 共同定义正性写回、导数恢复、停止和接受事件，之后把两种精度接到同一个标量化 Newton 驱动。固定偏压和外电阻通过方程与回调选择；物理残差、保护和报告模式由配置驱动。 |
| 扫描与预测 | main 的普通 DC 步长实现仍在 `DCSweep.cpp`，包含最新 `Extrapolate` 接入。SLOT 已将自适应步长抽为 `DCSweepStepControl.cpp`，供普通与宽精度调用。 | 以 SLOT 已标量化的控制器作为唯一实现，与 main 对照移入仍缺的状态和事件语义。预测器负责从已接受历史产生初猜，校正器提供更新次数和结果，步长控制器负责增减、落点、重试和回滚。 |
| 精度和存档 | SLOT 具有受支持配置的 80 位输入、SparseLU、实际写回和权威状态字复读。main 有既有输出合同、普通状态档和多种线性后端。 | 同一阶段合同按能力选择标量和后端。扩展精度从权威输入传到装配、写回、电流和接受档，显示文件按各自格式导出。同一发布器按精度选择存档编码，输出语义和失败处理一致。 |

## 实施顺序

1. **建立合并基线。** 从已提交的 SLOT 头建立临时集成分支，引入固定的 main 头并保留双方祖先。按配置和发射器、物理和装配、Newton 和步长、测试四组处理 26 个冲突。记录两边当前的算例清单、配置、接受点、失败重试和停止原因，作为趋势对照。合并阶段先保留两条已验证的求解路径，确认功能可重复后再逐层归并。
2. **先贯通统一入口。** 将 CMD/PAR/GRD/DAT 读入、阶段 JSON、清单和输出作为一个生产流程验证。沿用统一的 TDR/DF-ISE 掺杂来源判定，记录物种数据或净／总掺杂重建、参数作用域、接触来源和几何身份。用实际 SLOT CMD 的 Physics/Math/Solve 驱动 IALMob、接口几何及材料选择。验证 `Extrapolate`、`Save/Load`、`ExtendedPrecision(80)`、`File` 输出命名与各阶段覆盖。用同一入口验证 Workbench/SDS 启动。成功时发布曲线、场量以及适用的检查点，失败时保留日志与诊断，并以实际退出码和状态角色标识。
3. **归并物理和离散。** 逐项接上公共 Tri3 几何、SG 边驱动、接触、IALMob/HFS、Avalanche、SRH/Auger 和端口积分算子。对共享算子使用一个标量化实现，使生产残差和 Jacobian 对电势／浓度使用同一组物理量、有效分支和单位；直接以电势和浓度为未知量求偏导，准费米势仅作为模型所需的派生量。保持原算例的统计模型、材料参数和离散配置。PN2D 的 `element_edge_sg_gss_laux` 验证配置作为整体保留。优先用冻结态比较逐行 J/RHS、接触电流和方向差分，定位差异后再更新生产调用。
4. **统一步长与非线性执行。** 优先合并已共享的 `DCSweepStepControl` 与 main 的合法外推和回滚语义，再归并宽／普通 Newton 校正器。以同一个接受或拒绝事件驱动 CSV、曲线、正式存档和父态提升。逐分量保正、全局保正、保正后导数省略、RHS 增长保护及负载线报告由原配置与同一规则实现控制。线性后端按精度和能力选择，宽方向持续用于缺陷评估、试步、实际写回和下一轮装配。完成调度、校正及归档的配对后移除重复生产实现。
5. **全量验收后合回。** 完整构建集成版本并运行下面的全量 CTest。保存 CTest 清单、各配置实际启用的后端、JUnit 和失败日志。对每个失败定位并复验，最终重跑全量。完成 SLOT 完整运行及与基线的审查后，再将集成分支合入目标分支。

## 验收标准

- **全部 CTest 算例。** Windows UCRT64 Debug 完整构建后运行 `ctest --preset windows-ucrt64-debug --output-on-failure --no-tests=error`。在参考环境按 main CI 配置启用 PN2D、BVMethods、U-2022.12、TransportModels DD/DG 与 SingleDevice 参考注册项，并满足 HDF5 和相关后端的实际构建要求，运行 `ctest --test-dir <build-directory> --output-on-failure --no-tests=error`。对独立构建才注册的测试（例如 Python API）使用相应构建配置补全套件。每份 CTest 结果均覆盖该配置的注册清单。检查重复或缺失注册及跳过原因，必需算例全部执行并通过。开发过程中的筛选测试用于定位故障。
- **真实入口。** CTest 包含以真实 `sd-sim` 读取 CMD/PAR/GRD/DAT 的端到端用例，检查阶段 JSON、运行清单、日志、接受点、结果文件和失败退出码。由 CMD 中的精度请求触发真实分派。请求无法由模型、构建或后端满足时，启动前给出明确诊断。
- **数值行为。** 同一标量、状态、单位和模型设置下，共用路径的 Poisson 与连续性装配、电流、线性方向、接受／拒绝和步长事件与冻结基线相符。已保证逐位一致的字段维持逐位标准，其余采用既有算例容差。任何超差都定位到具体模型、求解或输出变化，确认后再更新基线。
- **SLOT 与参考算例。** 将 0 V Poisson、首次耦合、0→0.10 V 外电阻、Save/Load 同目标续算及拒绝点不可发布或成为父态纳入 CTest。用同一生产入口另做 SLOT 0→1 V 和完整 BVDS，核对接受点、不收敛重试、停止原因、更新次数、器件漏压、分载流子电流、负载线和终态物理场。0 V 和非零偏压的权威 80 位档都需验证无损复读。保留节点 10561 的少子浓度偏差作为现有基线限制，合并前后单独对照。
- **完成判定。** SLOT 与其他支持的输入共享一个生产入口、一个阶段与输出合同，以及各模型、离散、步长和 Newton 规则的公共实现。允许不同算例选择不同物理或离散配置，允许按精度选择兼容的线性后端。临时参考配置和旧校正器的删除以对应的通过证据为依据。

## 源码入口

- [CMD 阶段编排](../../../scripts/vela_sdevice/orchestrator/__init__.py)和[SLOT IALMob 参数绑定](../../../scripts/vela_sdevice/orchestrator/slot_ialmob.py)
- [普通浓度 Newton](../../../include/vela/simulation/PopulationCoordinateNewton.h)和[宽精度生产阶段](../../../src/tools/ExtendedInitialStage.inl)
- [DC 步长控制](../../../src/simulation/DCSweepStepControl.cpp)和[CTest 注册](../../../CMakeLists.txt)
- [PN2D BV 验证](../../validation/pn2d_bv_validation.md)和[SLOT 80 位及完整扫压现有证据](../../validation/slot_ldmos_extended80_wp1_wp4_2026-11-18.md)
