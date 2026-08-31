# Templates/LDMOS G3 box 系数与状态反馈审计（2026-08-31）

> 2026-08-31 更正：本报告的 `0.403147 dex / 2.530x` phin 反馈只适用于旧的
> dominant-contact-mean + mesh-default 自洽状态。使用已资格化的
> `legacy_node_local + external AverageBox` 31 点状态重跑后，中位反馈仅为
> `0.005615 dex / 1.0130x`；旧平台归因不再适用于修正合同。contact HFS 在同合同
> 冻结回放中仍是必要算子项。新证据见
> `templates_ldmos_g3_averagebox_curve_qualification_2026-08-31.md`。

## 结论

在 IALMob 与 predictor 持续关闭、接触 HFS 与 WP1.5 Jacobian 修复保持不变的
条件下，中段 `2.61--2.73x` 电流差异已分解为两个层次：

- 同一 Sentaurus 状态上的生产 Vela 算子地板在 `Vg=1/3--2/3 V` 仅为
  `1.046--1.049x`；
- 自洽状态反馈额外放大中位 `0.403147 dex`，即 `2.53015x`，四个审计点均由
  电子准费米势 `phin` 单变量主导。

现有 barycentric、负局部 cotangent truncation、mixed-Voronoi 节点体积不能
解释该反馈。一次短期 signed-cotangent 聚合原型也未改变主残差热点，且同偏压
reclose 电流只改变约 `0.12%`；该原型已撤回，没有进入生产配置。因而本轮不再
运行 31 点 coefficient 曲线，阶段 3 P95 门仍失败，known-difference ledger
继续保持 draft。

## 范围与冻结控制

- exact SProcess topology，`Vd=0.1 V`；
- G3 no-IALMob，内部 GradQF HFS + 已资格化的接触 ElectricField 回退；
- predictor、impact ionization 和 PN2D 私有 `element_edge_sg_gss_laux` 捆绑均未启用；
- 不调整 bulk mobility、不平移门压、不放宽任何 Newton/KCL 门槛；
- 只比较 exact CurrentPlot 偏压，不做陡峭区跨点插值评分。

## 1. 固定状态 coefficient A/B

在 `Vg=1/6, 1/2, 5/6 V` 的同一 Sentaurus G3 状态上，比较四个显式组合：

1. barycentric node volume + 生产正 barycentric fallback；
2. barycentric node volume + negative local cotangent truncation；
3. mixed-Voronoi node volume + 生产 fallback；
4. mixed-Voronoi node volume + truncation。

| Vg (V) | 生产算子 Vela/Sentaurus | truncated | mixed volume | mixed + truncated |
| ---: | ---: | ---: | ---: | ---: |
| 0.166667 | 1.79505 | 1.79 | 1.79505 | 1.79 |
| 0.500000 | 1.04621 | 1.05 | 1.04621 | 1.05 |
| 0.833333 | 1.04570 | 1.05 | 1.04570 | 1.05 |

30,022 条输运边中，truncation 改变 1,042 条边，mixed node volume 对固定 SG
edge flux 严格不变。端口积分没有出现可解释 `2.6x` 的变化，因此现有 box
开关被排除为主因。

Sentaurus 日志确认 `CVPL_AverageBoxMethod = TRUE`。网格包含 1,021 个钝角
单元，但只有 14 个 non-Delaunay 单元，其中硅区 7 个；最大
`CoeffIntersection` 异常位于 element 10508/8438、node 4492。该节点确实进入
continuity 热点清单，但 truncation 没有消除最大热点，说明不能把一个局部
non-Delaunay 修复外推为全曲线根因。后续若继续比较 AverageBox，必须取得
LDMOS 本算例的逐边 `MeasureCoefficients.debug` 作为直接 oracle，而不是继承
PN2D 模板策略或从后处理电流场反推。

## 2. 八组合状态反馈分解

在 `Vg=1/6, 1/3, 1/2, 2/3 V`，分别由 Sentaurus (`S`) 或 Vela (`V`)
提供 `psi/phin/phip`，穷举八种组合。所有组合使用同一个生产 SG/HFS 算子。
下表列出自洽 Vela 状态 `VVV`、固定 Sentaurus 状态 `SSS`，以及只替换电子
QF 的 `VSV`。

| Vg (V) | VVV/Sentaurus | SSS/Sentaurus | VSV/Sentaurus | 反馈放大 (dex) | 主导状态族 |
| ---: | ---: | ---: | ---: | ---: | :---: |
| 0.166667 | 2.73853 | 1.79505 | 1.79865 | 0.183442 | phin |
| 0.333333 | 2.73074 | 1.04900 | 1.05258 | 0.415504 | phin |
| 0.500000 | 2.68708 | 1.04621 | 1.04978 | 0.409664 | phin |
| 0.666667 | 2.60632 | 1.04568 | 1.04925 | 0.396630 | phin |

`psi` 单独替换只恢复约 `0.00145 dex`，`phip` 的 Id 敏感性近零；同时替换两种
QF 与只替换 `phin` 的结果等价。`phin` 的全场中位 offset 为零，去中心 P95
只有 `0.103--0.110 mV`，但最大局部差约 `11.3 mV`。这说明当前平台不是统一
参考能级偏移，而是小范围空间 QF 差异经亚阈值输运指数放大。

## 3. continuity 与空间定位

把 Sentaurus 的 `SSS` 状态代入 Vela continuity 算子后，电子残差几乎完全由
边 flux 构成，SRH/Auger 项的 L1 约 `6.3e-12`，不足以解释差异。以
`Vg=0.5 V` 为例：

- `SSS` electron residual：L1 `3.94960e-5`、L2 `2.42852e-6`；
- 已收敛 `VVV` electron residual：L2 `3.78e-10`；
- 最大行位于 node 3747/10233（`x≈-9.90 um, y≈3.88 um` 的 Si/SiO2
  界面），其次包括 node 4492、4538；
- 最大 flux-feedback 边包括 source 接触边 4570--4571，以及
  `x≈-9.979 um, y≈3.36--3.42 um` 的源侧/沟道侧边；其典型
  `couple/length` 为 `0.18--0.99`，并不是此前 drain 侧 `18--83` 的异常边族。

因此开放问题已从“drain 端口积分或全局 box 权重”收窄为：源侧 Si/SiO2
界面附近电子 continuity 的 QF/边界/界面离散语义，以及 Sentaurus AverageBox
在该局部区域的真实逐边系数。

## 4. 同偏压 reclose 停止门

从完全相同的 `Vg=0.5 V` Sentaurus 冻结状态启动单偏压 Newton reclose：

| profile | 收敛迭代 | drain Id (A/um) |
| :--- | ---: | ---: |
| production legacy | 6 | 1.73647e-10 |
| truncated negative cotangent | 6 | 1.73453e-10 |
| 短期 signed 聚合诊断原型 | 5 | 1.73434e-10 |

三者均回到约为 Sentaurus `2.68x` 的同一 Vela 分支。signed 原型虽将 SSS
continuity L1 从 `3.95e-5` 降到 `2.52e-5`，却不改变 node 3747 主热点；其
reclose 电流改善仅约 `0.12%`。根据单因素停止规则，该原型被撤回，31 点
coefficient 曲线不具备启动条件。

## 决策与下一开发项

1. 不保留新的全局 coefficient profile；现有生产离散默认不变。
2. Stage-3 L2 仍因 P95 `0.432753 dex > 0.20 dex` 失败；ledger 不审批。
3. 下一开发项继续属于 WP1.5/离散资格：导出 LDMOS 本算例 AverageBox 的逐边
   系数 oracle；对 node 3747/10233/4492/4538 的电子 continuity 行做边界类型、
   界面邻接、系数、密度平均和第一步 Newton 增量对齐。
4. 在上述局部算子证据关闭前，继续关闭 IALMob 与 predictor，也不把该差异
   登记为网格收敛的引擎固有地板。

## 可复现工件

- 状态八组合与 continuity/edge 审计：
  `scripts/audit_templates_ldmos_g3_state_feedback.py`；
- 同偏压 reclose：`scripts/run_templates_ldmos_g3_same_bias_reclose.py`；
- fixed coefficient A/B：
  `scripts/audit_templates_ldmos_g4_fixed_state_coefficients.py`；
- ignored 大型输出：
  `reference_staging/templates_ldmos_g3_state_feedback_20260831/`、
  `reference_staging/templates_ldmos_g3_contact_hfs_coefficient_ab_20260831/`、
  `reference_staging/templates_ldmos_g3_reclose_coefficient_ab_20260831/`。
