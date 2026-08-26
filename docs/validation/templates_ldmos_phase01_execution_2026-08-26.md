# Templates/LDMOS WP0、阶段 0/1 执行报告（2026-08-26）

## 结论

WP0、阶段 0 和阶段 1 的技术工作已完成。原始 Sentaurus oracle、代表状态、结构/掺杂/接触导入、确定性复跑和 exact-mesh Poisson 成本下界均通过相应技术门。

当前机器汇总状态为 `unresolved`、最高正式等级为 `L0`，唯一未关闭项是 `budget_freeze.approval.status=draft`：预算仍需 benchmark owner 与 independent reviewer 双签。因此本报告不授权进入阶段 1.5，也不声称经典 DD 或任何新增物理已经通过。

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

## 成本试跑与预算草案

成本探针是真正的线性 Poisson exact-mesh 结构下界，不调用 Newton，也不是材料或物理验收：

- 10,241 nodes、30,022 unique edges；Poisson 结构 nnz 上界 70,285；耦合 DD 未知量估计 30,723。
- 单偏置 1.608393 s；五点合计 4.561418 s；峰值 working set 26.95 MiB。
- unsigned 草案：方案 A base 0.55 h、worst 1.37 h；方案 B base 4.56 h、worst 11.40 h；组合 base 5.66 h、worst 14.14 h；并发度 1，存储草案 5/15/40 GiB。

此前 stage1_v3 曾误用 `dc_sweep + poisson_only`，但该路径仍先进入耦合平衡态并出现 `nonfinite_residual`。它仅作为阶段 1.5 的先验输入，不构成求解器缺陷结论；正式成本数据来自 stage1_v4 的线性 Poisson 探针。

## 最终 gate 与后续动作

| Gate | 状态 |
| --- | --- |
| official oracle and representative states | pass |
| exact topology structure | pass |
| exact-mesh cost probe | pass |
| budget double approval | unresolved |

建议批准 WP0、阶段 0/1 的技术产物，同时维持阶段 1.5 禁入，直至预算由 benchmark owner 和 independent reviewer 双签。双签后应重新生成 `validation_summary`，正式把最高等级从 L0 提升到 L1；阶段 1.5、经典 DD、hRecVelocity/IALMob/hQP、热和 Okuto 均需另行授权与资格门。

## 可复现证据位置

本机 ignored 根目录：`reference_staging/templates_ldmos_sentaurus2022/phase01_original_20260826_02/`。

- 正式 oracle 审计：`reports/oracle_audit.json`
- 最终验收摘要：`reports/validation_summary.json`
- 最终代表状态：`representative_states_final/state_inventory.json`
- 阶段 1 v4：`stage1_v4/reports/structure_audit.json`
- 成本与预算：`stage1_v4/cost_probe/cost_probe.json`、`stage1_v4/contracts/budget_freeze.json`

这些路径仅用于本地复核；其中的专有或大型产物不进入 Git。
