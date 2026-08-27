# Templates/LDMOS WP1.5 重新资格与阶段 3 入口审计（2026-08-27）

## 结论

WP1.5 的 0 V exact-mesh 平衡态与状态重启资格已经关闭，但有偏压固定点仍未关闭：

- G4 状态经 `poisson_only` bootstrap 后，以近冻结 QF 的 coupled Newton 关闭；
- 17 位状态 repeat 在 10241 个节点上满足 `max |dpsi|=0 V`、
  `max |dphin|=2.12e-30 V`、`max |dphip|=0 V`，电子/空穴密度逐节点不变；
- 0→0.1 V 漏压 ramp 可首次相对收敛，但同偏压 repeat 仍以
  `line_search_non_decrease` 失败，并出现 `0.826 V` 接触多数载流子 QF 跃迁；
- `contact_basin` 表示可把该跃迁降到约 `0.217 V`，但严格 repeat 仍失败；关闭
  初态相关的 continuity row scaling 也未关闭该门；
- 因此不运行 31 点 G3 Id-Vg，不认领阶段 3，也不开始 IALMob 实现。

当前最高正式等级仍为 `L1`。阶段 2 的 0 V 子门通过，WP1.5 继续以“有偏压
restart/continuation 鲁棒性”状态开启。

## 根因与修复

### 0 V：block filter 与数值地板

原始 `block_filter` 从同一 G4 oracle 启动时，在 9 次更新后以
`line_search_non_decrease` 失败。Poisson 块为 `2.274e5`，电子/空穴块已约为
`7.37e-13 / 9.44e-13`，接受步长多次降到 `1/2048--1/1024`。改为标准 merit
line search 后，Poisson 块可快速下降，但无保护的 coupled 更新会沿耗尽区病态 QF
方向产生约 `0.65 V` 的接触跃迁，仍不是合格平衡态。

生产路径最终采用：

1. `solver.method=poisson_only`，以接触 basin QF 关闭电荷/Poisson；
2. 0 V coupled cleanup 将 `quasi_fermi_update_limit_V` 固定为 `1e-12 V`；
3. 同偏压 repeat 与迭代 0 安全门验证序列化后的固定点。

| deck | 结果 | 迭代 | 终止原因 |
| --- | --- | ---: | --- |
| `g_contact_polysi_poisson_eq` | 通过 | 5 | `poisson_only_contact_basin_converged` |
| `g_contact_polysi_eq` | 通过 | 4 | `stall_residual_floor` |
| `g_contact_polysi_eq_repeat` | 通过 | 0 | `initial_stall_residual_floor` |

Newton floor 的初始、line-search 和 max-iteration 三条接受路径现在都要求载流子行、
全局连续性和接触多数载流子 QF 安全门通过。Poisson-only 在严格下降被舍入噪声阻塞时，
也只按相同的 Poisson/载流子块上限接受数值地板。

### MetalGate 扫描映射缺陷

首次 0.1 V repeat 曾出现“收敛”，但逐节点比较发现 453 个栅节点从
`0.5609636105202798 V` 被错误投影到 `0 V`。原因是 `DCSweep` 在被扫描接触为
`metal_gate` 时，用原始扫描电压覆盖了已经包含 flatband 的有效 Dirichlet 势。

现已统一为：

`psi_gate = applied_bias - flatband_voltage - work_function_eV`

修复后新增小型 gate-sweep 回归测试；先前的假通过作废。该缺陷不影响漏极等欧姆接触
扫描，也不改变非扫描金属栅的既有映射。

### 0.1 V：初态相关行缩放与 QF 病态仍开放

从合格 0 V 状态做 G3 高场 0→0.1 V ramp，首次 0.1 V 点在 4 次 Newton 更新后按
相对残差收敛，原始块残差为：

- Poisson：`1.895e-5`；
- electron：`4.96e-11`；
- hole：`2.82e-12`。

同偏压重新加载后，continuity row scaling 由新初态重新计算；求解 14 次后在
Poisson `3.70e-7` 处 line-search 非下降，并测得 `0.826 V` 接触多数载流子 QF
跃迁。`contact_basin` 路线的首次 0.1 V 点同样可相对收敛，但 repeat 在 Poisson
`2.89e-6`、接触 QF 跃迁 `0.217 V` 处失败。禁用 continuity row scaling 的受控
消融也没有通过 repeat。

这说明首程基于大偏压跳变初态建立的相对范数会接受一个无法在新行权重下重闭合的
状态；仅调整 QF reference 或移除行缩放不足以关闭问题。完整 G3 曲线必须继续停止。

## 代码与证据

- `NewtonSolver`：安全初始 floor、三条 floor 的统一接触 QF 门、Poisson-only 数值
  地板接受，以及真实接触 QF 失败诊断；
- `DCSweep`：扫描 metal gate 时保持 flatband/work-function 有效势；
- `prepare_templates_ldmos_phase23.py`：正式 Poisson bootstrap、0 V 近冻结 QF、
  0.1 V repeat deck 和合格 seed 依赖；
- `prepare_templates_ldmos_wp15_diagnostics.py`：block filter、merit、QF reference、
  QF 限幅、row-scaling、repeat 与 `max_iter=0` guard 消融矩阵；
- 大体量日志与状态位于 ignored `reference_staging/.../stage1_v4/`。

## 停止决定

阶段 3 入口保持关闭。下一步仍属于 WP1.5：使 continuity row scaling/收敛范数在
continuation 与同偏压 restart 间具备固定点不变性，并阻止 relative convergence
接受接触多数载流子 QF 不安全状态。只有 0.1 V strict repeat 通过后，才恢复 G3
精确点 Id-Vg；IALMob 不得用于掩盖该求解器问题。

## 回归验证

- `test_newton_solver`：98 cases / 1301 assertions；
- `test_dc_sweep`：101 cases / 3516 assertions；
- Templates/LDMOS Python regression：13 tests；
- 三组测试全部通过。
