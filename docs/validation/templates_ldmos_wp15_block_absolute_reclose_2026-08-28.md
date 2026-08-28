# Templates/LDMOS WP1.5 分块绝对收敛资格（2026-08-28）

## 结论

WP1.5 的本轮阻断已关闭。Vela 新增的分块绝对收敛、最佳接受迭代保存和
contact-majority QF 分支保护，在 Templates/LDMOS exact mesh、
`Vg=0 V / Vd=0.0231559774221138 V` 上通过资格验证。全过程保持 revision-4
物理、网格、SG 离散、无 predictor 和无 IALMob 不变。

正式 reclose 在第 5 次 Newton 迭代以 `block_abstol` 接受，得到
`2.1587697850e-15 A/um` 的漏极总电流。Sentaurus PLT 真值为
`2.0831669600e-15 A/um`，相对误差 `3.6292%`、幅值误差 `0.01548 dex`，通过
阶段 3 的 `0.15 dex` 门槛。

## 实现合同

新增的 `solver.block_absolute_convergence` 为显式、默认关闭的权威合同。
`mode=enforce` 时，标量 `reltol`、`abstol` 和 stall floor 不能绕过三个 raw L2
block ceiling。本次 exact-mesh 资格值为：

| block | ceiling | 接受态 |
| --- | ---: | ---: |
| Poisson `psi` | `5.0e-8` | `4.756945579e-8` |
| electron continuity `phin` | `2.0e-9` | `1.628469303e-9` |
| hole continuity `phip` | `3.0e-10` | `2.982312431e-10` |

新增 `contact_majority_qf_branch_drop_limit_V` 覆盖所有正式收敛路径和最佳迭代
资格。`contact_majority_qf_branch_guard_contacts` 可显式限定接触作用域。本次只保护
`drain`，门槛为 `5e-11 V`：初态漏极多数载流子 QF drop 为
`8.21284e-12 V`，接受态为 `1.17961e-16 V`，均通过。

全器件未限定作用域的最大 drop 初态约为 `1.90 V`，来自其他接触的既有状态；它
不能作为漏极 fA 分支的资格门。实现同时排除了两个不同接触之间的 boundary edge，
避免把 contact-to-contact 边误报为 contact-to-interior 分支。

## 最佳迭代保存

Newton 现在持续保存“残差最低且通过 QF guard”的接受态及其 iteration、block residual
和 QF drop。失败尝试启用 `rejected_state_directory` 时，DCSweep 除 parent、initial、
final 外还写出 `*_best.csv`，并在 `newton_attempts.csv` 中登记相应元数据。

故意保持旧 strict `reltol=1e-10` 的 exact-mesh 控制仍在第 35 次尝试因
`line_search_non_decrease` 失败；新链路保存了 iteration 34、残差
`4.352647783e-8` 的 `attempt_1_bias_0p023156_best.csv`。该控制的 best 与 final
SHA-256 相同，说明该路径的最后接受态恰好也是最低残差态；保存语义另由
DCSweep 回归覆盖。

## exact-mesh 守恒与局部残差

接受态固定回放结果：

- 原 checkpoint 的 electron continuity L2 `1.055080877e-4` 降至
  `1.628469303e-9`，改善约 `64790` 倍；
- 漏极第一圈 electron residual L2 为 `4.198419476e-12`，不再主导全局范数；
- SG cut 与 ContactCurrent 相对差 `4.93e-15`；
- ContactCurrent 与曲线电流相对差 `1.83e-16`；
- free-node SG divergence 与 continuity flux 最大绝对差 `1.03e-25`；
- stable SG 与 long-double 端口电流相对差 `4.02e-15`。

漂移/扩散拆分仍有约 `2.54e17` 的端口相消条件数，因此只保留作诊断，正式评分继续
使用 stable SG、long-double 闭合和 Sentaurus PLT。

## 测试与判定

- `test_newton_solver.exe`：101 cases、1328 assertions 全部通过；
- `test_dc_sweep.exe`：101 cases、3528 assertions 全部通过；
- `tests.regression.test_templates_ldmos_phase23`：21 tests 全部通过；
- exact-mesh block reclose：通过，`block_abstol`，5 iterations；
- exact-mesh strict failure control：按预期失败并保存 best state。

由此可重新进入 G3 曲线校核，但仍不得启用 IALMob 或 predictor。下一步应把相同
block/QF 合同应用到无 predictor 的精确 G3 偏压序列，检查各点是否保持单调、端口
守恒和 Sentaurus 对比门槛。

## 可复现产物

- 审计器：`scripts/audit_templates_ldmos_g3_continuity_sg_contact.py`
- exact-mesh 汇总（ignored）：
  `reference_staging/ldmos_g3_align_20260827/continuity_sg_contact_wp15_v4/summary.json`
- 正式配置（ignored）：同目录
  `vela_wp15_block_reclose/wp15_block_reclose.json`
- 正式曲线/状态（ignored）：同目录
  `vela_wp15_block_reclose/curve.csv` 与 `state.csv`
- strict best-state 控制（ignored）：同目录
  `vela_strict_reclose/rejected_states/attempt_1_bias_0p023156_best.csv`
