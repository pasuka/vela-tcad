# Templates/LDMOS G3 WP1.5 无 predictor 精确序列执行记录

日期：2026-08-28
状态：**按冻结合同失败关闭；阶段 3 曲线评分未授权**

## 结论

WP1.5 的 block-aware 绝对收敛合同已接入 G3 deck 生成器，并在无 predictor、
Sentaurus 原始偏压点上执行。两条预注册入口均在进入 31 点 Id-Vg 之前失败：

1. 从 Vela 映射平衡态沿 15 点 Sentaurus drain 轨迹运行时，`Vd=0` 和
   `0.001 V` 通过，在 `0.00214466666666667 V` 被 Poisson 绝对块门拒绝；
2. 从 Sentaurus `Vg=0 V / Vd=0.1 V` 精确状态做同偏压 G3 reclose 时，载流子块
   和 drain-majority QF 分支已合格，但 Poisson 块停在冻结门槛的约 12 倍。

因此没有生成或评分不合格种子派生的 Id-Vg 曲线，也没有计算 Vth、gm、SS 或强反型
误差。当前阻断项是有限偏压 Poisson/Jacobian 数值地板，不是 IALMob、高场迁移率
单位、predictor 或接触 QF 分支。

## 本轮执行合同

- 物理层：G3-no-IALMob，即 Fermi、OldSlotboom、SRH(DopingDep)、Auger 和
  HighFieldSaturation；无 IALMob、量子势、自热或 Okuto；
- 偏压：drain 使用 15 个 Sentaurus 已接受点；Id-Vg 保留 G3 oracle 的 31 个原始
  CurrentPlot 点 `0...5 V`，不插值、不跳点；
- `solver.block_absolute_convergence`：`psi <= 5e-8`、`phin <= 2e-9`、
  `phip <= 3e-10`；
- `reltol=1e-10`、`abstol=1e-14`、`stall_residual_floor=1e-12`；
- drain-majority QF 分支保护 `5e-11 V` 仅用于深关态种子和 drain 预偏置。
  它不用于导通区 Id-Vg，因为有限电流需要物理的准费米势梯度；
- sweep 中不存在 predictor 配置；所有失败尝试均保存 initial/final/best state 和
  逐迭代 CSV。

## 结果

### 路线 A：映射平衡态到 Vd=0.1 V

| 偏压 | 状态 | 最佳迭代 | psi 块 | electron 块 | hole 块 | drain QF drop |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 0 V | 通过 | 0 | `5.90675e-10` | `6.26877e-24` | `5.30492e-19` | `2.74e-29 V` |
| 0.001 V | 通过 | 10 | `6.94010e-10` | `8.06846e-12` | `3.06001e-12` | `4.12e-18 V` |
| 0.0021446667 V | 失败 | 27 | `5.42318e-7` | `8.39616e-11` | `5.57870e-11` | `7.94e-17 V` |

失败点的 carrier 块分别只占冻结上限的 `4.20%` 和 `18.60%`，QF 分支也远低于
保护限；只有 psi 块为上限的 `10.846x`。第 28 次迭代的 line search 因
`line_search_non_decrease` 拒绝，失败分类为
`block_absolute_convergence_line_search_rejected`。

在唯一可比较的非零接受点 `Vd=0.001 V`：

- Vela：`6.134079418203e-16 A/um`；
- Sentaurus：`1.325935383432e-16 A/um`；
- Vela/Sentaurus：`4.6262x`，绝对 log 误差 `0.66523 dex`。

该点尚低于可靠曲线分辨区，数值只作诊断，不作为阶段 3 电流硬门。

### 路线 B：Sentaurus Vg=0/Vd=0.1 精确种子 reclose

最佳接受迭代为 8：

| 指标 | 观测值 | 冻结上限 | 上限占比 |
| --- | ---: | ---: | ---: |
| psi 块 | `6.0248457e-7` | `5e-8` | `12.0497x` |
| electron 块 | `1.2098255e-10` | `2e-9` | `6.05%` |
| hole 块 | `2.0684965e-10` | `3e-10` | `68.95%` |
| drain-majority QF drop | `5.5511e-17 V` | `5e-11 V` | `1.11e-4%` |

第 9 次迭代同样因 `line_search_non_decrease` 结束。两条独立入口的最大 Poisson
残差都落在 node `4601`（`x=-9.9210837 um, y=0.9453125 um`，重掺杂接触邻域）：

- 路线 A：`2.5077043e-7`；
- 路线 B：`2.5452969e-7`。

该空间重合把后续工作收敛到接触邻域 Poisson/Jacobian 缩放和 line-search 地板。

## 生成器与证据

`scripts/prepare_templates_ldmos_phase23.py` 现在：

- 对所有 G3 deck 物化冻结的 block 绝对合同；
- 对深关态 seed/drain 预偏置物化 drain-majority QF 保护；
- 将生产 Id-Vg 种子绑定到 Sentaurus 精确状态的两次同偏压 reclose，而不是已知失败的
  直接 `0 -> 0.1 V` 跳变；
- 为所有 G3 运行默认启用 Newton history、attempt ledger 和 rejected best-state；
- 在 manifest 中显式记录 no-predictor 与 QF 保护的适用域。

未提交的可再生成运行证据位于：

`reference_staging/templates_ldmos_sentaurus2022/phase01_original_20260826_02/`
`stage1_v4/phase23_t2022_contract_v5_wp15/`

其中关键文件为：

- `g3_drain_prebias_sentaurus_path.csv`；
- `g3_drain_prebias_sentaurus_path_newton_attempts.csv`；
- `g3_drain_prebias_sentaurus_path_newton_iterations.csv`；
- `g3_drain_prebias_sentaurus_path_rejected_states/attempt_3_bias_0p002145_best.csv`；
- `g3_idvg_seed.csv`；
- `g3_idvg_seed_newton_attempts.csv`；
- `g3_idvg_seed_newton_iterations.csv`；
- `g3_idvg_seed_rejected_states/attempt_1_bias_0p000000_best.csv`。

## 判定与下一步

本轮完成了“将合同应用到完整精确序列”的资格执行，但触发了预注册 stop rule；阶段 3
仍停在 G3 种子资格门，不能进入 IALMob 或开启 predictor。下一开发项应继续属于 WP1.5：

1. 对 node 4601 邻域做 Poisson 行/列尺度、Jacobian 条件数和线性解增量审计；
2. 以两份保存的 best-state 做同状态残差/Jacobian 回放，确认是线性系统缩放还是
   nonlinear merit 地板；
3. 修复后原样重跑两条入口；只有 psi 块达到冻结 `5e-8` 门，才启动 31 点 Id-Vg；
4. 不放宽冻结阈值，不把失败 best-state 当作已收敛种子，不改变 G3 物理。

## 回归

- `python -m unittest tests.regression.test_templates_ldmos_phase23`
- 21 tests passed。
- `python -m unittest discover -s tests/regression -p "test_templates_ldmos*.py"`
- 44 tests passed；同时补齐 revision-4 BGN 合同中 `sentaurus_dEg0_eV` 和
  `sentaurus_qf_reference_mapping` 的 schema 声明。
