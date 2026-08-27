# Templates/LDMOS G3 Sentaurus 重闭合与偏压路径探针（2026-08-27）

## 结论

T-2022.03-SP2 在 Templates/LDMOS `G3-no-IALMob` 的 exact mesh 上可稳定完成
`Vg=0 V, Vd=0.1 V` 的同进程与 Save/Load 重闭合。状态重载不是该工程当前的
阻塞项；Vela 的有偏压失败不能归因于 Sentaurus 状态文件精度或该固定点本身不存在。

- Sentaurus 原始漏压 ramp 实际接受 15 个点（含 0 V），不是 `0 -> 0.1 V` 单步；
- 最终 `0.097722431245 -> 0.1 V` 步只需 2 次 Newton 更新；
- 同进程重闭合和 Save/Load 重闭合均只需 1 次更新；
- Save/Load 前后逐节点势和准费米势差远低于 WP1.5 的 `10 uV` 门；
- Vela 复刻同一接受点序列时只到达 `0.023155977422 V`，在下一点
  `0.031641657941 V` 经 100 次迭代失败，因此路径粒度只是放大因素，不是完整根因。

WP1.5 仍保持开启，阶段 3 入口继续关闭。

## Sentaurus 运行合同

独立派生目录：

`/root/sentaurus_runs/vela_oracle_2022/templates_ldmos/g3_vd0p1_reclose_20260827_01`

派生 deck 保持 G3 的 Electrode、Physics、材料参数和初始求解不变，只执行漏压
ramp，并追加 `CNormPrint`、`NewtonPlot(Error Residual Update)`、Save/Load 和两次
同偏压 Coupled。原始 G3 目录没有修改。

有效 Math/Solve 语义为：

| 项目 | T-2022.03-SP2 有效值 |
| --- | ---: |
| `Digits` | 5 |
| `Iterations` | 25 |
| `NotDamped` | 100 |
| `ErrRef(Poisson)` | `0.025852 V` |
| `ErrRef(Electron/Hole)` | `1e10 cm^-3` |
| 非线性求解器 | Bank/Rose |
| 线性求解器 | blocked decomposition |
| drain ramp | `InitialStep=0.01`, `Increment=1.35`, `MinStep=1e-5`, `MaxStep=0.2` |

原始 G3 接受的 drain 电压为：

`0, 0.001, 0.002144666667, 0.003668599556, 0.005697462208,
0.008398554687, 0.011994609139, 0.016782156301, 0.023155977422,
0.031641657941, 0.042938927273, 0.057979358509, 0.077722431245,
0.097722431245, 0.1 V`。

## Newton 与端口证据

最终 ramp 步的 RHS 序列为 `60.5 -> 30.6 -> 0.0604`，最终 error 为
`1.9221e-2`。同进程重闭合从 `|Rhs|=0.0604` 出发，一次更新后 error 为
`2.1827e-8`；Save/Load 重闭合从 `|Rhs|=0.985` 出发，一次更新后 error 为
`2.1816e-8`。

| 状态 | drain eCurrent (A) | drain hCurrent (A) | drain conduction current (A) |
| --- | ---: | ---: | ---: |
| ramp 到 0.1 V | `3.29654212850941e-15` | `1.47808651812696e-19` | `3.29668993716122e-15` |
| 同进程重闭合 | `3.29664613907816e-15` | `1.47808643324888e-19` | `3.29679394772149e-15` |
| Save/Load 重闭合 | `3.29664599425547e-15` | `1.47808642050069e-19` | `3.29679380289752e-15` |

同进程与 Save/Load 重闭合的电子电流相对差约 `4.39e-8`。ramp 首次接受值与
重闭合值相差约 `3.16e-5`，仍明显小于当前器件曲线门槛，但说明只看低精度终端
日志不足以验证状态等价。

## 逐节点状态比较

Vela `sentaurus_import` 从常规 Plot TDR 导出了 Potential、e/h QuasiFermiPotential
和密度；NewtonPlot TDR 另外包含 Poisson/e/h continuity RHS、Error 和 Newton
Update。常规 Plot 的逐节点 max-norm 如下：

| 比较 | max `|dpsi|` (V) | max `|dphin|` (V) | max `|dphip|` (V) | max rel `dn` | max rel `dp` |
| --- | ---: | ---: | ---: | ---: | ---: |
| ramp -> 同进程重闭合 | `5.05e-15` | `1.04e-12` | `1.10e-13` | `4.04e-11` | `4.12e-12` |
| 同进程 -> Save/Load 重闭合 | `1.11e-16` | `7.77e-16` | `1.94e-16` | `2.03e-15` | `2.38e-15` |

准费米势相对误差在接近零的节点没有判别意义，故资格以绝对 max-norm 为准。
Save/Load 结果相对同进程结果处于 double 舍入量级。

## Vela 同路径对照

生成器新增独立 `g3_drain_prebias_sentaurus_path` 诊断 deck；生产
`g3_drain_prebias` 仍保持原配置。Vela 成功接受到 `0.023155977422 V`，随后在
`0.031641657941 V` 失败：

- 终止：`max_iterations`，100 次更新；
- 最终 Poisson 块：`1.61997e-3`；
- 最终 electron continuity 块：`6.99440e-4`；
- 最终 hole continuity 块：`4.09691e-9`；
- 失败前 line search 多数只接受 `1/1024--1/256` 的步长。

在 `0.023155977422 V`，Vela drain current 为 `2.6161e-11 A/um`，而 Sentaurus
同点约为 `2.0832e-15 A`。两边 2-D 电流归一化合同仍需在最终定量引用前复核，但
Vela 的电子 drift/diffusion 分量各约 `2.676e-2 A/um`、靠高位相消得到净电流，
已经表明深截止区存在严重的通量闭合/状态差异，不能用 continuation 步长修复。

## 后续决定

1. 把 Sentaurus 的常数参考势 `4.6337 V`、`ErrRef` 和 Bank/Rose 归一化纳入
   Vela 固定状态 residual/Jacobian 审计，而不是继续盲调 line-search 参数。
2. 在已接受的 `0.023155977422 V` Vela 状态上回放 SG edge flux、节点 continuity
   RHS 和端口积分，并与同偏压 Sentaurus TDR 对齐，优先解释 drift/diffusion 的
   高位相消误差。
3. 以 `0.023155977422 -> 0.031641657941 V` 为最小失败转移，导出 Vela 第 0/1 次
   Newton 状态，与本次 Sentaurus NewtonPlot 的 residual/update 字段做行级比较。
4. 上述固定状态差异关闭前，不放宽 restart 门、不启动 IALMob，也不运行完整 G3
   Id-Vg。

## 证据封存

本地 ignored 目录：

`reference_staging/.../phase23_t2022_contract_v2/g3_vd0p1_reclose_probe_20260827_01/`

证据包 `g3_vd0p1_reclose_evidence.tgz` 的 SHA-256 为
`486a542cf5ddcb064770d2a587411e41e51f0824421575509b62f8e028e3b735`。包内包含
派生 deck、完整 log、timing、三组 PLT/Plot TDR、Save 状态和两次重闭合前后
NewtonPlot TDR。专有/大体量运行产物不提交仓库。
