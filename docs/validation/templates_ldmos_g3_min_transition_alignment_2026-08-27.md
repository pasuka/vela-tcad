# Templates/LDMOS G3 最小失败转移算子对齐（2026-08-27）

## 结论

本轮完成了 `Vd=0.0231559774221138 V` 固定状态上的 Sentaurus/Vela SG edge
flux、continuity RHS 与端口积分对齐，并提取了
`0.0231559774221138 -> 0.0316416579413075 V` 的 Sentaurus Newton 第 0/1 次
状态。阶段 3 仍保持阻断，但原因已经从泛化的“Newton 不稳”收敛为三个可独立验证的
问题：

1. G3 生成 deck 的 mobility 语义和单位均不符合 oracle。Sentaurus G3 只启用
   `Mobility(HighFieldSaturation)`，Vela 却生成 `masetti_field`；同时生成器把 TCAD
   内部的 `cm2/(V s)`、`cm^-3`、`cm/s` 数值又转换为 SI，而 `unit_scaling`
   解析器实际要求直接保留 TCAD 内部数值。
2. Sentaurus 的历史外推将目标点第 0 次 electron continuity RHS 的 L2 范数从
   `19.012 A` 降至 `2.2244e-8 A`，影响约 `8.55e8` 倍；但 Vela 现有 secant
   predictor 单因素控制在 `0.0021447 V` 已失败，不能直接作为等价修复。
3. 把 Sentaurus 已收敛状态载入 Vela 后，carrier residual 很小，但 Vela Poisson
   residual L2 达 `2.9303e5`。这与 mobility 无关，说明仍有电势参考、材料/BGN、
   介质或接触边界合同没有闭合。

因此，当前证据不支持启动 IALMob，也不支持通过放宽 Newton 门槛绕过失败。应先修正
G3 mobility 合同/单位，再分别资格验证 predictor 和固定状态 Poisson。

## 对齐合同与数据

- Sentaurus oracle：T-2022.03-SP2，exact mesh，G3-no-IALMob 派生 deck；原始模板
  未修改。
- 固定状态：Sentaurus 正常 continuation 在 `0.0231559774221138 V` 接受的 TDR，
  与 Vela 无 predictor 路径最后接受的同偏压 restart CSV。
- Vela 算子：生产 `sg_edge_flux_probe`、`edge_mobility_probe`、
  `newton_residual_probe`、`newton_carrier_term_probe` 和
  `newton_step_probe`，未另写替代公式。
- SG 空间对齐：将 Sentaurus 节点电流密度投影到 Vela primal edge，再按 Vela
  dual-couple 长度积分。这是空间诊断，不是两个离散体系的逐边恒等式。
- 端口对齐：Vela 为 SG contact-cut 的导电电流、默认 `1 um` 深度；Sentaurus TDR
  的 `ContactCurrentFlux` 可能包含与 PLT 导电电流不同的分量，故当前只作诊断项。
- 所有大型导出、派生 deck、日志和 CSV 均位于 ignored 的
  `reference_staging/ldmos_g3_align_20260827/`，不进入版本库。

## Sentaurus 最小转移与 predictor 证据

正常历史路径在 `0.023155977422 V` 重闭合误差约为 `9.98e-8`，随后到
`0.031641657941 V` 的 Newton 过程为：

| 路径 | 求解器显示初始 RHS | 第一次更新后 RHS | 结果 |
| --- | ---: | ---: | --- |
| 保留 continuation 历史/Extrapolate | `498` | `340` | 两次更新后收敛 |
| Save/Load 后直接跳点、无 predictor | `1.17e12` | `8.11e11` | 单步探针按预期终止 |

原生 NewtonPlot 的未缩放 continuity 数据也给出相同结论：

| 状态 | history electron RHS L2 | no-predictor electron RHS L2 | 放大倍数 |
| --- | ---: | ---: | ---: |
| iteration 0 | `2.2244e-8 A` | `19.012 A` | `8.55e8` |
| iteration 1 | `9.0616e-9 A` | `12.179 A` | `1.34e9` |

两条路径的 iteration-0 状态并不相同：硅区 `psi` 最大差恰为本次漏压步长
`8.48568 mV`，电子/空穴密度最大差分别为 `0.1422/0.1425 dex`；到
iteration 1，`psi` 最大差增至 `17.62 mV`，电子/空穴密度最大差为
`0.5668/0.6039 dex`。所以 Save/Load 直跳不是正常 continuation 的等价复现。

Vela secant predictor 控制只改变
`sweep.continuation.predictor={mode: secant, fields: [psi,phin,phip]}`。它接受
`0` 和 `0.001 V`，在 `0.0021446667 V` 第 12 次更新发生
`line_search_non_decrease`；无 predictor 基线则接受到 `0.0231559774 V` 后才在
下一点失败。结论是：Vela predictor 接口已经存在，但当前实现/阻尼策略尚未通过
此 exact-mesh 路径的资格门，不能直接写入生产 deck。

## Mobility 与 SG edge flux

当前生成器把电子低场 mobility `1417` 写为 `0.1417`，把饱和速度 `1.07e7`
写为 `1.07e5`。但 `unit_scaling` 的 `mobilityToInternal()` 是恒等映射，现有单测也
以 `1417`、`9.68e16` 和 `1.07e7` 作为内部输入。固定 Sentaurus 状态回放显示：

| 指标 | 当前 `masetti_field` | HFS-only 可表示代理控制 |
| --- | ---: | ---: |
| 全部可比边电子 mobility 中位误差 | `4.0000 dex` | `0.0000 dex` |
| 主载流边电子 mobility 中位误差 | `4.0000 dex` | `0.0000 dex` |
| 主载流边 SG 线电流中位误差 | `3.8843 dex` | `0.1193 dex` |
| 主载流边 SG 线电流 p95 误差 | `4.2702 dex` | `0.9027 dex` |

代理控制使用现有 `caughey_thomas_field`，令 `mu_min` 等于材料低场 mobility，
从而消除 doping limiter，仅保留 high-field limiter。它不是最终生产语义，但足以证明
当前约四个数量级的局部 SG 差异主要来自生成 deck 的 mobility 合同，而不是 SG
公式本身。

## Continuity RHS 与固定状态闭合

Sentaurus history iteration 0 的原生 RHS 极值位于漏极侧沟道/结区：electron RHS
最大 `4.0041e-9 A`，坐标 `(-9.9659, 3.6416) um`；无 predictor 时最大值变为
`9.3879 A`，坐标 `(-9.9764, 10.4844) um`。Vela 基线在 100 次更新后的电子残差
热点也位于漏极顶部附近，说明它更接近 Sentaurus 的无 predictor 分支，而不是正常
history 分支。

在各自已接受的 `0.023155977422 V` 状态上，Vela residual L2 为：Poisson
`1.5902e-6`、electron `1.1651e-12`、hole `4.8315e-11`。把 Sentaurus 同点状态
输入相同 Vela 算子后，electron/hole 仍仅为 `8.2006e-13/1.3618e-13`，但 Poisson
跃升至 `2.9303e5`。HFS-only 控制不改变这个 Poisson 数值，因此下一次 fixed-state
审计应拆解电荷、介电、BGN/统计、电势零点和接触边界项，而不是继续调整 mobility。

硅区状态差也支持这一判断：中位 `psi` 差为 `6.93 mV`、电子密度差仅
`0.00155 dex`，但空穴密度中位差达 `1.878 dex`；最大范数受重掺杂接触/极低载流子
节点支配，只保留作定位，不作为曲线验收门。

## 端口积分

固定 Sentaurus 状态经当前 Vela mobility 回放的 drain SG cut 为
`-9.6284e-18 A/um`，Sentaurus TDR `ContactCurrentFlux` 为
`-3.0517e-15 A`，幅值低 `2.501 dex`。HFS-only 代理则得到
`-9.8292e-14 A/um`，反向超出约 `32.2` 倍。另一方面 Sentaurus PLT 在同偏压的
电子加空穴导电电流约为 `2.0832e-15 A`，与 TDR 标量也不完全相同。

所以本轮将逐边空间趋势作为 mobility/SG 的主要证据，端口数值保持诊断等级；在把
Sentaurus displacement/conduction 定义和 Vela `1 um` 深度合同写成同一积分前，
不得用端口倍率反推材料参数。

Vela 自身的一致性已复核：最后接受状态的 SG drain cut 为
`2.61607e-17 A/um`，与 DCSweep 的 `current_total_A_per_um` 相同。之前报告中的
`2.6161e-11 A/um` 实为未乘 `1e-6 m` 深度的内部电流列，现已纠正。

## 决策与后续实施顺序

1. 修正 WP1.75/生成器的单位合同，并为 G3 增加明确的 HFS-only mobility 语义；
   `masetti_field` 只能用于 deck 明确包含 `DopingDependence` 的分支。
2. 增加生成 deck 单测：在 `unit_scaling` 下必须写入 `1417`、`9.68e16`、
   `1.07e7` 等 TCAD 内部值，禁止再次转换为 SI。
3. 在修正 mobility 后重跑无 predictor 的精确偏压序列；再单独调查 Vela secant
   在第三个偏压点的 line-search 失败。两项不得合并改动。
4. 用 `newton_carrier_term_probe`/Poisson 分项探针完成 Sentaurus 固定状态的
   电荷与边界审计；Poisson 闭合前不启动 IALMob。
5. 最后建立同语义的 Sentaurus 导电端口积分，关闭 TDR/PLT/SG cut 三者的定义差异。

## 可复现交付物

- 分析器：`scripts/analyze_templates_ldmos_g3_min_transition.py`
- 回归测试：`tests/regression/test_analyze_templates_ldmos_g3_min_transition.py`
- ignored 汇总：`reference_staging/ldmos_g3_align_20260827/analysis/summary.json`
- Sentaurus 派生 deck 与原始输出：
  `reference_staging/templates_ldmos_sentaurus2022/phase01_original_20260826_02/`
  `stage1_v4/phase23_t2022_contract_v2/g3_min_transition_alignment_20260827_01/`

分析器的 7 项单元测试全部通过；本轮只新增验证工具和报告，并修正文档中的单位误读，
没有修改生产 solver、物理模型或生成器。
