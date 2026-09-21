# Templates/LDMOS WP0、阶段 0/1、阶段 1.5 与 WP1.75 执行报告（2026-08-26）

## 结论

WP0、阶段 0、阶段 1、阶段 1.5 和 WP1.75 的技术工作已完成。原始 Sentaurus oracle、代表状态、结构/掺杂/接触导入、确定性复跑、exact-mesh Poisson 成本下界、状态重启/回放资格和版本化配置合同均通过相应技术门。

预算已于 2026-08-26T12:41:25Z 完成双签，当前机器汇总状态为 `pass`、最高正式等级仍为 `L1`。阶段 1.5 与 WP1.75 已关闭进入经典 DD 前的资格门；它们不声称经典 DD 曲线或任何新增物理已经通过。

## 执行边界与版本

- 工作分支：`codex/templates-ldmos-phase-a`
- 原始 oracle run id：`phase01_original_20260826_02`
- Sentaurus：T-2022.03-SP2；已现场封存 sprocess、sdevice、svisual 三项版本横幅。
- VM 原始源归档 SHA-256：`f20cb16396ce5e1f938d8e93a6a0f4eb6e2a3e2e2a976d0cec3cd3b904e6288f`
- 正式原始运行只物化 Workbench token，不改变 Applications Library 物理 deck。
- TDR、PLT、日志与压缩包全部保存在 ignored 的 `reference_staging/`；仓库只跟踪中性工具、合同、测试和本报告。

## WP0 交付

已交付并验证以下基础设施：

- `validation_summary`、`known_difference_ledger`、`threshold_freeze`、`budget_freeze` 和 phase-0/1 manifest 的 schema、校验器与 Markdown 渲染器；
- 隔离式 VM runner、源文件/产物 SHA-256 清单、PLT 中性 CSV 归一化和 fail-closed oracle 审计；
- 仅增加非 loadable `Plot` 的代表状态派生 deck、状态捕获、状态合并和端口/代表点分类；
- Sentaurus TDR 导入对缩写活性掺杂字段、显式 `cm -> um` 坐标合同和精确接触边对的支持；
- 阶段 1 结构审计、确定性复跑、exact-mesh Poisson 成本探针和预算草案生成。

## 阶段 0：原始 oracle

原始链路四个 stage 均以退出码 0 完成：SProcess 约 99.4 分钟，IdVg 约 5.9 分钟，IdVd 约 15.2 分钟，BV 约 15.6 分钟。最终 `n1_fps.tdr` SHA-256 为 `07af6353c47c03a740583b873bac58e733a377a8bb9fa7e08baa3d30bc2c67f3`。

| 曲线 | 点数 | 偏置范围 | 最大绝对漏极电流（A/um） | 审计 |
| --- | ---: | ---: | ---: | --- |
| IdVg | 31 | 0–5 V | 5.3543439e-6 | 有限、偏置非降 |
| IdVd, Vg=4 V | 31 | 0–40 V | 1.2817568e-4 | 有限、偏置非降 |
| IdVd, Vg=8 V | 32 | 0–40 V | 2.4928770e-4 | 有限、偏置非降；保留初始重复点 |
| BV | 91 | 0–50.333494 V | 1.2737045e-8 | 有限、保留 continuation 原始行序 |

BV 在 `Iadapt=6.5e-13 A/um` 后从电压控制转入电流 continuation，存在 30 次电压回退。这是合法路径，不得排序，也不得在陡峭段跨点线性插值评分。

最终代表状态集合为 21 个：IdVg 5 个、IdVd 12 个、BV 4 个。BV 状态覆盖：

| 角色 | 状态 | Vd（V） | Id（A/um） |
| --- | --- | ---: | ---: |
| pre-Iadapt / near-Iadapt | `state_bv_path_0000_des.tdr` | 50.3283893 | 6.2577957e-13 |
| avalanche growth | `state_bv_path_0001_des.tdr` | 50.3118115 | 9.1122936e-11 |
| criterion pre | `state_bv_path_0002_des.tdr` | 50.3108940 | 8.5865643e-9 |
| criterion post | derived final device state | 50.3113603 | 1.1264309e-8 |

第一次派生捕获得到的 56 个 TDR 中，49 个文件名以 `state_` 开头、7 个为输入/最终结构；其中 32 个 BV `state_bv_path` 实际只覆盖 0–1 V 预段。该集合被审计识别为无效 BV 代表集并保留作诊断证据，没有进入最终验收。修正版使用原始 PLT 的 continuation `time` 精确选点；派生运行在越过 1e-8 A/um 后比原始运行稍早结束，故门槛后状态使用该派生运行的最终 TDR，并记录其实际 PLT 电压/电流，而非冒用未命中的原目标值。

阶段 0 最终八项 oracle gate 全部为 `true`。

## 阶段 1：结构、网格和离散合同

结构审计在相同输入上重复执行并得到逐字节一致结果。

- 10,241 个全局顶点、19,782 个三角形、3 个体区域；坐标范围 x=[-10.32, 10] um、y=[0, 11] um。
- Sentaurus 原始坐标被显式声明为 cm，统一乘以 1e4 导出 um；最大坐标误差为 0。
- Silicon_1、Oxide_1.1、Oxide_1.2 的三角形集合和面积逐区域完全一致，最大面积相对误差为 0。
- 五个接触的精确边集合全部闭合：gate 452 边、drain 24 边、source 13 边、th_lat 101 边、substrate 9 边；`th_lat` 明确为 thermal-only，不写入电学偏置。
- NetActive 的 p95/max 对数误差均为 0 dex，符号不一致为 0，所有导入掺杂有限且非负。
- 最终结构中不存在 PolySilicon 体区域，gate 边由 Oxide_1.1 所有；后续 PolySi 功函数必须作为接触/界面语义验证，不能假设存在 PolySi bulk region。

网格不是非钝角/Delaunay 网格：有 1,016 个钝角三角形、13 条非 Delaunay 内边和 1,041 个负的原始 half-cotangent 权重；最小角 0.0630404°、最大角 174.451474°。因此阶段 1 草案将 `legacy_cell_reconstructed`、barycentric control volumes 和精确导入接触边冻结为候选公共离散合同，但仍标记 `physics_use_authorized=false`。PN2D 的 `element_edge_sg_gss_laux` 原子捆绑不得据此推断为 LDMOS 或全局默认。

## 成本试跑与已批准预算

成本探针是真正的线性 Poisson exact-mesh 结构下界，不调用 Newton，也不是材料或物理验收：

- 10,241 nodes、30,022 unique edges；Poisson 结构 nnz 上界 70,285；耦合 DD 未知量估计 30,723。
- 单偏置 1.608393 s；五点合计 4.561418 s；峰值 working set 26.95 MiB。
- 已批准入口预算：方案 A base 0.55 h、worst 1.37 h；方案 B base 4.56 h、worst 11.40 h；组合 base 5.66 h、worst 14.14 h；并发度 1，存储预算 5/15/40 GiB。

此前 stage1_v3 曾误用 `dc_sweep + poisson_only`，但该路径仍先进入耦合平衡态并出现 `nonfinite_residual`。它仅作为阶段 1.5 的先验输入，不构成求解器缺陷结论；正式成本数据来自 stage1_v4 的线性 Poisson 探针。

## 阶段 1.5：重启、冻结回放和同偏压重闭合

资格运行使用 10,241 节点 exact mesh、真实节点掺杂和三个代表状态：0 V 平衡态、Vg=0/Vd=0.1 V 的 IdVg 低漏压点，以及 Vg=4/Vd=0.1 V 的可分辨 IdVd 预偏置点。状态 CSV 使用完整 double 精度（17 位有效数字）。oxide-owned gate 在本阶段按 `metal_gate`、临时 `flatband_voltage=0` 处理；这只用于求解器资格，PolySi 功函数映射仍由阶段 2 的 `G-contact/poly` 控制负责。

| 资格门 | 结果 | 门槛 |
| --- | ---: | ---: |
| 五个经典持久化场 binary64 round-trip 最大差 | 0 | 0 |
| `frozen_state` 最大变化 | 0 | 0 |
| 0 V Vela→Vela `psi` 最大差 | 3.3014052 µV | ≤10 µV |
| IdVg settled→repeat `psi` 最大差 | 0.0020366 µV | ≤10 µV |
| IdVd Vela→Vela `psi` 最大差 | 0.0226950 µV | ≤10 µV |
| IdVd 漏极电流相对差 | 8.6705685e-9 | ≤1e-3 |
| IdVd 归一化 KCL 不平衡 | 1.9775148e-9 → 7.3097225e-11 | 不恶化 |

IdVg 首次 imported→settled 的 `psi` 变化为 80.8376 µV；首次求解由相对步长条件提前结束，后续 residual-floor 收敛后再次重闭合仅变化 2.04 nV。因此资格门使用 settled→repeat，不把 Sentaurus→Vela 的边界、物理和离散差异误判为序列化误差。该 IdVg 点的电流约 1.58e-19 A/um，低于预注册的 1e-15 A/um 绝对分辨率地板，其相对电流和 KCL 只报告、不评分。

重掺杂 exact mesh 的冷 Gummel 初始化仍可在 Silicon/Oxide 共享节点产生空穴密度溢出；现在会报告具体节点、载流子和缩放量，而非下游非有限残差。使用仓库既有的 Fermi、OldSlotboom、准费米更新限幅、block-filter line search 和 continuity row scaling 后，从 Sentaurus 保存状态进入 Vela 固定点的三组 Newton 运行全部收敛。本阶段没有修改核心 Newton 算法。

## WP1.75：版本化材料、物理和离散合同

已新增并通过三份 `additionalProperties:false` 的版本化 schema 和 golden 合同：

- `materials.json`：浓度固定为 `cm^-3`、迁移率固定为 `cm^2/(V*s)`，能量、温度和三类热参数也显式声明单位；C++ `MaterialDatabase` 会拒绝未知字段、错误单位和越界值，并在 legacy SI 与 TCAD internal 两套内部单位下正确换算。
- `physics_contract.json`：分开记录 Fermi、OldSlotboom、SRH/Auger、bulk/high-field mobility、电子/空穴量子、Okuto 和 thermal 的启用状态与实现状态；待开发功能不能伪装成已实现运行时键。
- `discretization_contract.json`：冻结 phase-A classical 的 SG edge flux、barycentric control volume、cell-reconstructed field/source 和精确接触边积分为一个原子 profile；明确禁止在 phase A 启用 avalanche，也不从 PN2D profile 推断全局默认。

校验器覆盖单位、范围、未知键、canonical JSON round-trip 和显式版本迁移。legacy material 迁移必须由调用者声明源单位为 `legacy_si` 或 `tcad_internal`，不允许猜测；并包含针对历史上迁移率单位误写导致约 3000 倍电流坍缩风险的数量级回归。

## 最终 gate 与后续动作

| Gate | 状态 |
| --- | --- |
| official oracle and representative states | pass |
| exact topology structure | pass |
| exact-mesh cost probe | pass |
| budget double approval | pass |
| exact-mesh restart qualification | pass |
| WP1.75 versioned contracts | pass |

双签记录：benchmark owner 为 `Ted Chin (explicit Codex task authorization)`；independent reviewer 为 `OpenAI Codex evidence reviewer (non-human)`。复核确认 oracle manifest 哈希一致、探针运行全部成功、三档 scenario 等于 scope 预算求和，且当前约 2.11 GiB staging 小于最低 5 GiB 存储预算。

WP0、阶段 0/1、阶段 1.5 和 WP1.75 现以 L1 关闭。方案 A 的下一执行项是 WP2 / 阶段 2–3 的经典低压 DD、`G-contact/poly` 控制和单因素曲线校核；hRecVelocity/IALMob/hQP、热和 Okuto 仍需各自证据门与范围约束。

## 可复现证据位置

本机 ignored 根目录：`reference_staging/templates_ldmos_sentaurus2022/phase01_original_20260826_02/`。

- 正式 oracle 审计：`reports/oracle_audit.json`
- 最终验收摘要：`reports/validation_summary.json`
- 最终代表状态：`representative_states_final/state_inventory.json`
- 阶段 1 v4：`stage1_v4/reports/structure_audit.json`
- 成本与预算：`stage1_v4/cost_probe/cost_probe.json`、`stage1_v4/contracts/budget_freeze.json`
- 阶段 1.5：`stage1_v4/qualification/qualification_summary.json`
- WP1.75：`stage1_v4/contracts/wp175/qualification_report.json`

这些路径仅用于本地复核；其中的专有或大型产物不进入 Git。

后续执行结果见
`docs/validation/templates_ldmos_phase23_execution_2026-08-27.md`：Sentaurus 单因素链已完成，
但 Vela G4 经典平衡态同偏压 reclose 触发停止门，当前仍只认领 L1。
