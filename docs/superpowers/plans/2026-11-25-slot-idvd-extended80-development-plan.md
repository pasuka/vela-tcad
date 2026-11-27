# SLOT LDMOS IdVd 端到端 80 位开发计划

2026-11-25。本文在 [sd-sim 80 位总计划](2026-11-18-sd-sim-extendedprecision80-development-plan.md) 的基础上，定义 SLOT IdVd 的实施范围和验收。

## 背景与目标

原始 `IdVd_des.cmd` 的 `Math` 段已能请求 `ExtendedPrecision`。当前运行器只为 Boltzmann＋IALMob 的外接电阻 BVDS 和固定端电压 IdVg 声明了端到端 80 位能力。IdVd 的 Fermi＋PhuMob/Lombardi＋雪崩阶段生成后，缺少对应能力声明和实际消费链，执行前报 `E202`。已有 VM 短扫使用单独核验的 JSON 和 binary64 主状态，可作为对照基线。

最终交付入口为 `sd-sim`。使用原始 SLOT Workbench CMD 展开的等温输入，在 300 K、13,881 节点网格上，从 \(V_G=V_D=0\) V 冷启动，保持 \(V_D=0\) V 栅预扫并保存 \(V_G=5\) V 父态，再按原脚本 `Save/Load` 进行漏压扫描。未知变量为电势和电子、空穴浓度，使用 Unit 缩放和长双精度 SparseLU。以 **0–2 V 短扫及复读**作为首轮整链验收，以原始 **0–60 V** 扫描作为完整算例验收。

所称 80 位涵盖权威输入及派生量、主状态、物理模型、残差、Jacobian、约束消元、线性方向、写回、收敛判断、预测与扫压控制，以及存档和跨进程复读。Linux VM 上以 `numeric_limits<long double>` 的基数 2、64 位有效尾数、最大指数 16384 核验。原始文件中已经量化的数据保留其来源精度，整数拓扑和仅用于展示的输出沿用现有表示。

## 实施步骤

| 工作包 | 实施内容 | 完成证据 |
|---|---|---|
| **WP0 基线和物理契约** | 绑定原始 CMD、展开后的 PAR、网格、掺杂和运行程序的身份。登记 Fermi、Bennett–Wilson 的 Fermi 修正、PhuMob＋Lombardi、高场饱和、硅区 SRH(DopingDep)／Auger、`Avalanche(GradQuasiFermi)` 及接触电场回退。核对修正的 PolySi 带隙、DOS 和参考势，记录单位和参数来源。归档已核验 Vela 0 V 父态、6 mV 开关对照及 Sentaurus 对应状态和迭代日志。 | 实际配置可追溯到输入。0 V 和 6 mV 的状态、接触边界、坐标和报告单位均已绑定，供后续同态比较。 |
| **WP1 CMD 和材料接口** | 复用 CMD 解析、材料导入、阶段生成和 `Save/Load`。补齐 IdVd 的 Bennett–Wilson 默认选择及 Fermi 修正，核对 PolySi、PhuMob/Lombardi、HFS、SRH、Auger 和雪崩参数从 PAR 到实际消费的单位。保留网格、掺杂、材料、温度、接触与步长的权威数值。为 IdVd 增加有明确物理范围的能力绑定，等 WP2–WP4 接通后发布。 | 原始 CMD 生成冷启动、栅预扫、复读、漏压扫描及电流文件切换。诊断能指出未支持的阶段、模型或后端。精度请求在各阶段均有记录。 |
| **WP2 Fermi 状态与冷启动** | 复用宽 Poisson、`ExtendedPopulationState`、接触约束和装配骨架。扩展原先面向 Boltzmann 的材料／接触初始化接口。实现宽精度的 \(F_{1/2}\)、逆函数及导数，统一实际使用的 Bennett–Wilson、有效带边、DOS、热电压和欧姆接触目标。Poisson 末态直接交接给浓度耦合求解。 | 稀薄区、过渡区和简并区的正反换算及导数通过独立高精度参考。Vela 自生零偏压 Poisson 和首次耦合均通过，状态低位在阶段交接中保留。 |
| **WP3 输运、源项与直接偏导** | 将有效带边 SG、PhuMob、Lombardi、HFS、SRH、Auger、雪崩及电流归约接入同一宽状态。先在内部坐标构造边电势差，再进行单位换算。残差与 Jacobian 使用同一套操作数和分支，由 \(\psi,n,p\) 直接播种 AD 或计算解析偏导。保留单元散射、接触电场回退及接触电流定义。 | 用冻结的 0 V、6 mV 同态做分项与方向差分，记录自由连续性行和接触行。重点核对单元 21242 与节点 10848–10850，并复查漏端节点 10237。对近零电流模分别核对分支、一侧方向和可解析的电流变化。 |
| **WP4 Newton、扫压、CLI 和存档** | 复用 `ExtendedFixedBiasNewton`、`ExtendedSparseLu`、宽状态档与 `Save/Load` 交接。接通宽矩阵和 RHS、约束消元、精化、线性质量检查、写回、保正、电流、更新误差和实际扫压值。复用宽扫压和外推，满足电势线性、正浓度对数外推及拒绝后的现有回退规则。 | 固定矩阵和低位哨兵验证线性方向及写回。权威 HDF5 复读后的状态有效位一致，同配置下下一次装配、端口电流和接受决定可复现。真实 CLI 通过能力检查并进入求解。 |
| **WP5 VM 验收与回归** | 先做新父态的 6 mV 重放，再完成 0–2 V 外推开、关配对及续算，最后用原始步长设置执行 0–60 V。运行受影响的专项回归和完整 CTest，并复跑现有 80 位 BVDS、IdVg。 | 发布阶段精度／后端清单、验收日志、接受态及结果文件。在共同漏压范围报告工作点、回退、Newton 更新和物理量差异。 |

WP0 先固定证据。WP1 的前端准备和 WP2–WP3 的物理实现可按共享输入、状态及装配接口推进，WP4 复用现有 80 位流程完成接入，WP5 验收。实现时将模型和标量选择接入共享算子，避免新增独立的 IdVd Newton 求解器。

## 数值策略与试验开关

- 对 6 mV 分叉，先对照同一父态的原始 SG、稳定边驱动和近零电子电流斜率处理。现有试验中，组合策略用了 7 次更新，仅斜率保护用了 15 次。两者是低压验收锚点，迭代次数可随精度和停止规则变化。
- 近零电子电流斜率策略若拟纳入正式路径，应独立配置并记录触发单元、作用的偏导和状态近零尺度。雪崩残差继续按当前状态求值。用可测电流下的方向差分和高偏压复验检查其适用边界。
- `VELA_IDVD_PARENT_AUDIT` 保留在历史试验记录中。新路径验收通过后，从后续运行脚本和求解入口移除这项依赖。日志开关配对须复现同态矩阵、RHS 和接受轨迹。
- 命令脚本的 `Digits=5`、电子与空穴 `ErrRef=1e10`、25 次 Newton 上限、外推和回退语义应能核查。记录实际使用的停止分支。`RHSmin` 的等价性按独立的原生范数映射证据验收，转换数值前明确单位和归一化。

## 验收标准

1. **输入与能力**。原始等温 CMD 经 `sd-sim` 用同一输入包完成分派，每个阶段实际使用 80 位状态和同类型 SparseLU。日志记录各阶段输入、程序版本、精度、模型和线性后端。缺少受支持后端、阶段精度不匹配及热输运配置有明确的失败原因。
2. **精度传递**。用 \(1\) 与 \(1+2^{-60}\) 或相应物理尺度的对照验证权威文本输入及派生量、工作状态、矩阵／RHS、求解方向、实际写回和存档复读。冻结同一个 \(A,b\)，在原始坐标重算 \(Ad-b\) 和分量后向误差，对照方向降精度产生的差异。通过现有验收条件并报告实际裕量。
3. **物理与导数**。使用同一网格和同一权威状态，比较电势、电子／空穴浓度和准费米势、接触电流、连续性分项、Jacobian 作用与方向差分。参考 Sentaurus 时注明取样偏压、缩放和同态关系，分别报告精度误差、模型差异及停止判据差异。
4. **连续运行和续算**。从修正材料的 0 V 冷态完成零漏压栅预扫，保存并复读 \(V_G=5\) V、\(V_D=0\) V 父态。首个 6 mV 和 2 V 终点均须接受，重启配对从同一接受档、同一目标继续。统计时排除零偏压初始点，分别记录接受和失败试探的 Newton 更新、回退次数、接受理由和耗时。
5. **完整扫描和回归**。在 VM 上按原始 CMD 的绝对步长运行，首步 6 mV、最大 1.2 V、最小 \(6\times10^{-11}\) V、增步系数 1.2。在原判据下接受 60 V 并正常退出，复核最后工作点、电流和权威存档。全部 CTest 通过，并单列已启用的 SLOT 80 位 BVDS、IdVg、IdVd 运行证据。未达到的漏压记录为未通过，保留最后接受点及失败原因。

首轮交付先报告 **80 位精度链通过情况**和**器件扫压通过情况**。高精度格式本身只回答精度链问题。完整端点能否接受，以 VM 的原条件验证结果为准。

## 代码入口及证据

- [现有 80 位 SLOT 阶段](../../../src/tools/ExtendedInitialStage.inl)、[宽状态装配接口](../../../include/vela/equation/ExtendedPopulationDevice.h)、[固定端压 Newton](../../../include/vela/solver/ExtendedFixedBiasNewton.h)、[宽 SparseLU](../../../include/vela/solver/ExtendedSparseLu.h)及[宽状态存档](../../../include/vela/io/ExtendedPopulationArchive.h)。
- [Fermi 统计](../../../src/physics/CarrierStatistics.cpp)、[PhuMob 雪崩 SG 电流](../../../include/vela/equation/PhuMobSgSourceCurrent.inl)及[单元雪崩装配](../../../include/vela/equation/PhuMobElementAvalancheAssembly.inl)。
- [CLI 能力检查](../../../scripts/vela_sdevice/orchestrator/__init__.py)、[阶段生成](../../../scripts/vela_sdevice/deck_emit/__init__.py)及[现有 IALMob 绑定](../../../scripts/vela_sdevice/orchestrator/slot_ialmob.py)。
- [配置说明](../../config_schema.md)、[IdVd 预检](../../validation/slot_ldmos_idvd_preflight_2026-10-05.md)、[VM 工作流程](../../office_vm_workflow.md)及[回归指南](../../../tests/regression/README.md)。
