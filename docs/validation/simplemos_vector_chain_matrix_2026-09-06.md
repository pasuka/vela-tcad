# 向量驱动场链式导数四工作点扩展验证：通过

2026-09-06。接续已封存的 [n23 低 Vd 隔离验证](simplemos_vector_chain_validation_2026-09-05.md)，将同一候选扩展至 n19/n23 × 两个 Vd。新增 30 个自洽状态全部严格合格；连同复用的 10 个状态，共核对 40 个状态。四工作点、两档幅度的补项后伴随预测全部通过原 0.1% 门槛，最大相对差为 **5.5628e-8**。新增实际方向 Jv 的 18 项检查全部通过，电流加权导数缺陷至少缩减 **99.992964%**。

这使浅沟道固定电荷方向的 Jacobian/伴随校准覆盖了高低 NWell 和高 Vd。自洽电流最大相对变化仅 **1.1380e-9**，高 NWell 的 Vela–Sentaurus 电流偏差仍为 **+14.3336% / +12.2626%**。正式源码、端口电流定义和默认模型本轮未修改；M82/M83 未放行。

## 工况、网格及计算方法

| 器件 | Vd (V) | Vg (V) | 网格节点 | 浅沟道自由节点 | 本轮计算 |
|---|---:|---:|---:|---:|---|
| n19 | 0.05 | 0.9 | 1480 | 184 | 新增开/关各 5 个状态 |
| n19 | 1.0 | 0.8 | 1480 | 184 | 新增开/关各 5 个状态 |
| n23 | 0.05 | 0.9 | 1482 | 186 | 复用已封存开/关各 5 个状态 |
| n23 | 1.0 | 0.8 | 1482 | 186 | 新增开/关各 5 个状态 |

网格沿用 `build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics/vela/{n19,n23}/mesh.json`，材料沿用 M65 的 `matched_materials.json`，初态为 `build-release/simplemos_strict_flux_precision/<case>/state.csv`。完整路径与哈希见本轮 `freeze.json`。低/高 Vd 分别固定上述不同 Vg，不能把本轮四点结论当作完整 Id–Vg 扫描。

沿用 Boltzmann、无 BGN、`equal_ni_flux_evaluation=compensated_log_expm1`、SRH，以及 `phumob_field_lombardi / total_impurity / quasi_fermi_gradient / transport_cell_vector`。温度及其他物理设置保持原配置。复用已封存的 `chain_runner.exe`，以环境变量 `VELA_VALIDATE_VECTOR_CHAIN_ENABLE` 控制数值链式导数补项。新增进程显式使用 `VELA_LINEAR_SOLVER=sparselu`；伴随仍用既有校准过的 SparseLU 迭代精化。UMFPACK/SPQR 虽在构建中可用，本轮未选用。复用已导出的网格和状态，无新增 TDR 导入或 Sentaurus 运行。

扰动为固定体电荷数密度 0、±1e13、±5e12 cm^-3，按 `Q_i=q·amplitude·1e6·area_i` 转为 C/m，只加入 Poisson 固定电荷项；电流均为 A/um。浅沟道沿用深度不超过 0.05 um、横向绝对坐标不超过 0.125 um 的 Si 节点，排除 Poisson Dirichlet 节点。n19 采用几何 signed-Si 面积，n23 采用已校准的原生 signed-Si 面积；总选定面积分别为 1.25840665237e-14 / 1.25842653017e-14 m²。n19 尚无对应原生局部有限差分资格。

配置、源项、可执行文件和既有判据在计算前冻结。三个新工作点各 6 次只读预检：关补项的五个电荷源，加开补项的零电荷对照。全部证明载流子残差、掺杂与本征浓度不变；Poisson 源的相对 L2 误差不超过 1.122e-16，选定范围外无源项，基态残差开/关逐条一致。

## 收敛与独立接受检查

保持 `max_iter=200`、`reltol=1e-7`、`abstol=1e-12`、ψ 阻尼 0.35、准费米更新上限 0.025 V，以及原有 carrier-row 的源/通量筛选和下限。载流子行 `eps_row=1e-6` 强制接受；每个输出状态另用原电荷诊断程序进行只读检查，确认载流子行零违规、全局配置满足及 `KCL/Id≤1e-8`。本轮新增状态的最大 KCL/Id 为 1.1639e-10。

状态通过的是既有严格组合判据，收敛原因仍包括 `reltol` 和 `stall_residual_floor`，不表示所有状态都达到了绝对残差 1e-12。全局独立审计仍采用 `tolerance=1e-6, source_floor=1e-10`；所有状态的电子/空穴净源强度未获得 `global qualified` 资格，不能宣称净 SRH 达到百万分之一相对闭合。

响应同时检查正负方向、双幅度线性、偶次分量、零控制漂移及信号/漂移比，全部沿用上一轮门槛。没有根据计算结果调整阈值。

## 伴随预测与实际自洽响应

下表响应为全幅中心差 `(Id(+1e13)-Id(-1e13))/2`。误差定义统一为 `abs(FD/prediction-1)`，百分数列已乘 100；最后一列取两档幅度中的最大相对误差，未乘 100。

| 工况 | 原预测误差 (%) | 补项后预测 (A/um) | 补项后实际响应 (A/um) | 补项后最大相对误差 |
|---|---:|---:|---:|---:|
| n19 / Vd=0.05 | 1.357347 | 5.91149318997e-10 | 5.91149286113e-10 | 5.5628e-8 |
| n19 / Vd=1 | 4.155719 | 1.04376185283e-9 | 1.04376182379e-9 | 2.8939e-8 |
| n23 / Vd=0.05（复用） | 0.340785 | 9.31760017220e-12 | 9.31760043514e-12 | 2.8220e-8 |
| n23 / Vd=1 | 0.568835 | 5.07844893199e-12 | 5.07844913542e-12 | 4.0058e-8 |

原预测在低 Vd 偏低、高 Vd 偏高；补项在两种方向都恢复真实响应。原算法 8/8 幅度未达预测门槛，补项后 8/8 通过。原算法自洽状态本身全部严格合格，预测失败不能解释为这些状态未收敛。

三次新伴随的相对残差不超过 6.013e-17，报告基态电流与相应严格初态一致。补项开启后，全幅正/负扰动的迭代数变化为 n19 低 Vd：8/8→3/3；n19 高 Vd：17/17→5/5；n23 高 Vd：26/8→7/5。全部工况详见 `iteration_ledger.csv`，不据此声称通用收敛速度提升。

## 沿实际状态方向的 Jv 复核

每个新工作点从原算法实际正负解构造全幅、半幅方向。准费米变化先以十进制合并 reference 与 increment，再转为双精度；在同一严格基态比较补项开/关的 Jv 与残差中心差分。测试完整方向及 ψ、φn、φp 分量，h=1、0.5、0.25。冻结来源状态及方向文件后运行，共 12 次只读调用、639648 条行记录。

全部残差差分和端点残差差分在开/关之间逐条相同。用原伴随权重投影完整方向的导数缺陷，两档幅度、三档步长、三个新增工作点的 18 项检查全部达到原 99% 缩减门槛；最差仍缩减 99.992964%。复核同时检查完整行数、唯一键和有限数值。此投影验证不自动赋予其他空间源、连续性源或微小空穴信号相同资格。

## 高低 NWell 配对与真实误差

以原已审计 Sentaurus 参考和本轮零控制电流计算；本轮没有重新运行商业程序。

| 工况 | 补项后零控制 Id (A/um) | 相对 Sentaurus 偏差 (%) |
|---|---:|---:|
| n19 / Vd=0.05 | 2.05202952590e-6 | +4.838749 |
| n19 / Vd=1 | 3.02197531931e-6 | +4.140359 |
| n23 / Vd=0.05 | 1.78761432612e-8 | +14.333599 |
| n23 / Vd=1 | 9.22887930098e-9 | +12.262584 |

高减低的电流误差仍为 0.03765204645 dex（低 Vd）和 0.03261596583 dex（高 Vd），开/关无有意义变化。全幅真实 `ΔId/Id` 的高/低 NWell 比为 1.80932545 / 1.59320495；这只是当前两种网格、既定源项下的局部敏感度，不是绝对误差的归因量。

这轮定位并消除的是伴随预测中的导数缺项。物理残差未变，绝对 Id 偏差与高低配对差没有改善。n23 低 Vd 与已封存 Sentaurus 原生有限差分的归一化响应差仍约 -0.709%；其他三点没有新增原生局部有限差分，不能把本地校准推广为跨程序一致性。

## 后续边界

下一步可使用这四个已校准状态及权重，继续建立沟道/界面输运与迁移率、Poisson 离散源项的误差贡献账本；明确每项状态、单位、空间支持，再选一个可由原生同扰动复核的候选。连续性源必须另做独立校准，不能沿用固定 Poisson 电荷的资格。

正式纳入 Jacobian 补项前，仍需覆盖其他迁移率/驱动分支及完整稀疏邻接模板。本轮未修改正式电流算法、未启用 cell-first 默认，也未放行 M82 的物理候选 A/B 或 M83 的完整 16×51 点回归。

## 产物与验证

便携账本位于 `reference_tcad/simplemos_sentaurus2022/vector_chain_matrix/`；输入、原始状态、日志和 Jv 逐行输出位于忽略目录 `build-release/simplemos_vector_chain_matrix_20260906/`。本轮新增 30 次非线性重闭合、30 次独立接受检查、18 次源项预检、3 次完整伴随、12 次 Jv 调用；无新 Sentaurus 任务。

- `contract.json / freeze.json`：冻结输入及门槛。
- `dc_ledger.csv / calibration.csv / current_identity.csv`：严格状态、两档响应、开/关电流一致性。
- `jvp_checks.csv / observer_checks.csv`：实际方向导数与伴随资格。
- `physical_error_ledger.csv / paired_error_ledger.csv`：绝对误差及高低 NWell 配对。
- `result.json / review.json / evidence.json`：结果、独立逐项审阅及哈希封存。

执行的 Python 回归共 11 项全部通过：新增响应判定 6 项，加既有 M79 资格判定 5 项。只增诊断编排与报告，没有新编译求解器或宣称完成其他生产分支回归。

从当前 worktree 根目录复核：

```powershell
& D:/msys64/ucrt64/bin/python.exe -m unittest tests.regression.test_simplemos_m79_terminal_response_qualification tests.regression.test_simplemos_vector_chain_matrix
& D:/msys64/ucrt64/bin/python.exe scripts/validate_simplemos_vector_chain_matrix.py verify
& D:/msys64/ucrt64/bin/python.exe -c "import sys; sys.path.insert(0, 'scripts'); import validate_simplemos_vector_chain_matrix as m; m.a.verify(m.OUT/'review.json')"
```

`prepare/preflight/run/analyze/jvp` 已执行，封存目录禁止覆盖重跑；后续复现实验应使用新的输出目录及冻结契约。
