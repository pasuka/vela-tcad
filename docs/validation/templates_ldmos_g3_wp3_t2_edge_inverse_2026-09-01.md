# Templates/LDMOS G3 WP3 T2 逐边所需修正逆问题（2026-09-01）

## 判定

T2 已完成，但**没有提名候选族，也不构成停止证据**。

在两个冻结端点，最强的单属性单参模型为 HFS 驱动力与 Bernoulli
`|eta|`：

| 属性族 | 1.0 V 残差能量解释 | 置换 p | 1.1666667 V 残差能量解释 | 置换 p |
| --- | ---: | ---: | ---: | ---: |
| HFS 驱动力大小 | 80.8544% | 0.003 | 68.4140% | 0.013 |
| Bernoulli `|eta|` | 62.7502% | 0.019 | 73.5976% | 0.009 |
| 钝角单元归属 | 40.7555% | 0.083 | 39.4729% | 0.109 |
| 掺杂梯度 | 16.2185% | 0.325 | 24.1429% | 0.212 |
| 界面平行取向 | 5.3991% | 0.575 | 10.2768% | 0.420 |
| 界面法向取向 | 4.4083% | 0.606 | 8.4789% | 0.470 |

HFS 与 `|eta|` 的置换结果显示它们并非随机标签下的普通水平，但没有任何单一
属性在两个端点同时达到冻结的 `>=90%` 门。最优属性还从 1.0 V 的 HFS
驱动力切换为 1.1666667 V 的 `|eta|`，因此不能提名 H1/H2，更不能宣称根因。

三种加权最小范数解的边级结构也不稳定。1.0 V 的三组两两余弦为
`0.8902, 0.6576, 0.3779`，1.1666667 V 为
`0.9153, 0.6781, 0.4670`；冻结稳定门要求全部 `>=0.9`。两端点的 top-5
Jaccard 最低均为 `0.25`。这与图的高度欠定性质一致，禁止把某一正则下的
热点边当作定向证据。

停止标准也未通过：两个端点共 14 个逐节点
`required/sum(incident |phi|)` 样本以合并中位数符号为参照，反号比例为 `42.8571%`
（通过严格 `>15%` 子门），但幅值只跨 `1.64952` 个数量级（未达到 `>=2`）；
并且最强单属性在较弱端点仍可解释 `68.4140%` 残差能量，不满足“无属性
`>=50%`”条件。三项必须同时成立，所以 T2 不得作为停止分类证据。

## R4 图合同

- 节点域严格为冻结七个自由 Silicon 电子 continuity 行：
  `3721, 3974, 3973, 4091, 3727, 4021, 3771`；所有行均验证
  `electron_gauge=electron_boundary=0` 且存在电子输运通量。
- 每态 T1 的 40 条拓扑入射边中，24 条满足 `couple>0` 且电子迁移率大于零，
  构成 Silicon 输运边域。八态边集合逐条相同，共分析 192 条 edge-state 记录。
- 规范方向严格为小节点 ID 指向大节点 ID，`B(tail)=+1`、
  `B(head)=-1`。图为 `7x24`、秩 7，其中 1 条域内共享边以反对称符号同时进入
  两行，23 条跨域边只进入域内端点。
- 规范定向逐边重构 T1 节点电子输运项，八态最大相对 L2 误差低于
  `1e-10`；三种逆解的 `B delta=-r` 最大相对闭合误差为
  `3.385e-16`。

T1 `state_manifest.json` 中的 mesh、状态、carrier/SG probe config 均重新做
SHA-256 核验；carrier/SG CSV 必须位于同一 ignored T1 包内。任一 schema、
八偏压、冻结物理、节点域、哈希、边集合、行自由度或矩阵秩不符，工具都会
fail closed。

## 三种正则化

求解统一写成：

```text
min sum_e (delta_e / scale_e)^2
subject to B delta = -r
```

- `unweighted_l2`：`scale_e=1`；
- `abs_phi_weighted`：`scale_e=|phi_e|`；
- `couple_over_length_weighted`：`scale_e=couple_e/length_e`。

后两者先按态内最大值无量纲化；精确零或极小尺度采用最大值的 `1e-9` 显式
floor，避免把零基线边偷偷变成硬约束。每个解均要求节点域满行秩和
相对守恒闭合 `<=1e-8`，实际远低于该限值。

边属性的单参模型为
`delta_phi_e = alpha * phi_e * normalized_attribute_e`，用
`B delta_phi` 对 `-r` 的残差能量解释率评分。每个属性、每个状态执行 999 次
固定种子的边标签置换；常数通量比例没有作为候选属性，因为
`r` 本身由同一 T1 SG 输运表重构，纳入它会产生 `alpha=-1`、100% 解释的
代数恒等式而非物理发现。

## 属性派生与数据边界

- `|eta|`、物理电子边通量、`couple/length`、端点坐标和 HFS 驱动力直接来自
  T1 `sg_edges.csv`。
- 掺杂梯度来自同一 T1 `carrier_terms.csv` 的两端净掺杂；carrier 表的
  `cm^-3` 数值按既有诊断合同换算到 `m^-3` 后除以边长。
- 界面法向由 T1 manifest 哈希绑定 mesh 中每个界面节点的 Silicon 相邻单元
  质心方向派生，边方向分解为法向/平行绝对投影。
- 钝角标签由同一 mesh 的 Silicon 三角形边长平方判据派生，并按每条边相邻
  Silicon 单元的钝角比例记录。

这些 mesh 派生字段不是 Sentaurus assembly 级逐边真值。T2 输出只是“使
Sentaurus 状态成为 Vela 零点所需的一组可行修正”，不是实际测得的
Sentaurus/Vela 逐边通量差。

## 复现

从本 worktree 根目录运行：

```powershell
D:\msys64\ucrt64\bin\python.exe scripts\audit_templates_ldmos_g3_edge_inverse.py `
  --t1-root reference_staging\templates_ldmos_g3_wp3_t1_residual_scaling_20260901 `
  --output-dir reference_staging\templates_ldmos_g3_wp3_t2_edge_inverse_20260901 `
  --permutations 999
```

ignored 输出：

- `reference_staging/templates_ldmos_g3_wp3_t2_edge_inverse_20260901/summary.json`；
- `edge_corrections.csv`；
- `nodal_required_scales.csv`；
- `single_parameter_fits.csv`；
- `report.md`。

## 验证

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest `
  tests.regression.test_audit_templates_ldmos_g3_edge_inverse -v
D:\msys64\ucrt64\bin\python.exe -m py_compile `
  scripts\audit_templates_ldmos_g3_edge_inverse.py `
  tests\regression\test_audit_templates_ldmos_g3_edge_inverse.py
```

6 个合成/合同测试覆盖：规范方向与跨域边、共享边反对称、守恒闭合、三正则
稳定夹具、强属性置换检出、未知 T1 schema，以及非规范边/秩亏矩阵的
fail-closed 拒绝；与 T1 聚焦测试合计 9 项通过。
本轮未运行 reclose、31 点曲线或 VM，也未修改生产默认、IALMob、predictor 或
冻结计划。
