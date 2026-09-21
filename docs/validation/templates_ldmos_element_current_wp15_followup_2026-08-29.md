# Templates/LDMOS element-current 与 WP1.5 后续执行记录

日期：2026-08-29

状态：WP1.5 低电流 KCL 门通过；element-current 资格失败；阶段 3 仍因曲线
P95 门失败而保持 fail。

## 范围与不变量

本轮保持 T-2022.03-SP2 exact 10,241-node topology、G3-no-IALMob 物理、
`Vd=0.1 V`、31 个精确 CurrentPlot 门压点不变。IALMob、predictor、量子、
自热和雪崩均关闭；没有修改 flatband、迁移率参数或曲线阈值，也没有施加
后处理门压平移。

## Element/cell-aware carrier-current profile

新增显式 opt-in 的
`mobility.carrier_current_discretization = "element_qf_gradient"`：

- 对每个输运 Tri3 单元重构 P1 electron/hole quasi-Fermi gradient；
- 以单元平均载流子密度和 mobility 计算矢量 particle current；
- 用 P1 test-function gradient 组装守恒 continuity residual；
- terminal current 在精确 contact boundary edge 的相邻输运单元上积分；
- Jacobian 包含完整 cell stencil；SG edge profile 仍为默认；
- 当前实验 profile 禁止 surface mobility 和 coupled impact ionization。

合成测试证明仿射场精确性、单元守恒、terminal boundary integration 和
Jacobian/finite-difference 一致性。但在三个 Sentaurus 固定状态上的生产算子
回放没有改善误差：

| Vg (V) | Sentaurus terminal (A/um) | SG error (dex) | element error (dex) | element / Sentaurus |
| ---: | ---: | ---: | ---: | ---: |
| 0.166667 | 6.77043e-14 | 1.38913 | 1.39141 | 24.6269 |
| 0.500000 | 6.46225e-11 | 1.04443 | 1.05700 | 11.4026 |
| 0.833333 | 9.88947e-8 | 1.04392 | 1.05653 | 11.3902 |

几何实现与使用 Vela bulk mobility 的独立手工积分一致。早期仅抽取少数样本后
写下的 `58--528 cm2/(V s)` 区间不完整，完整 G3 nodal 范围为
`22.66--1014.32 cm2/(V s)`，接触一环 element 范围为
`20.78--1015.71 cm2/(V s)`。deck 原文、日志声明、G4 平衡态全域
`1417 cm2/(V s)` 和全场 Masetti 拟合共同排除了 DopingDependence；低值来自
HighFieldSaturation 的接触支撑语义。因此不能把早期使用 Sentaurus mobility
的 `1.72x` Python 重构当作 Vela 自洽 profile 的结果，但原因也不能再表述为
缺少 low-field mobility spatial support。

按停止规则，本 profile 没有进入自洽 G3 曲线。它保留为实验 profile，后续若要
继续资格，必须先实现有来源的“硅体 GradQF + 接触边界 ElectricField 回退”
HFS 支撑并做 G3/G4 单因素回放；不得全局切换 ElectricField、启用 Masetti 或
标定 bulk mobility 数值取得表面拟合。

固定状态生成物位于 ignored 路径：

`reference_staging/templates_ldmos_g3_shift_kcl_20260829/audit_element_profile/`。

## WP1.5 低电流 continuity/line-search 修复

严格 reclose 将 electron block ceiling 从 `2e-9` 收紧到 `1e-11` 后，定位到
node 5564：相邻 node 5569 的绝对 electron QF 都打印为 `0.1 V`，但一 ULP
差异经细长边 coupling 放大为 `1.16819e-10` continuity residual。逐节点补偿
求和和 global closure 均不能关闭该误差。

有效的单因素数值合同是：

1. `quasi_fermi_reference = contact_basin`，以接触盆地参考值加小增量保存 QF；
2. `line_search_mode = block_filter`；
3. 对启用 `block_absolute_convergence` 的 filter，只要求当前仍违反 ceiling 的块
   取得充分下降；已经合格的块使用其绝对 ceiling 作为 envelope 下界；
4. 保持 `psi <= 5e-8`、electron `<= 1e-11`、hole `<= 3e-10`，不放宽硬门。

严格 `Vg=1/6 V` 同偏压 reclose 在一次 Newton 后达到：

- psi block：`9.29470e-10`；
- electron block：`1.71518e-13`；
- hole block：`3.80127e-13`；
- four-terminal KCL residual：约 `1.18e-21 A/um`；
- relative KCL error：约 `9.3e-9`。

对比控制中，`contact_basin + merit` 停在 electron block `5.34e-11`；把 ceiling
临时放宽到 `6e-11` 虽可收敛，但 KCL 仍为 `1.10%`，没有通过冻结门。因此最终
采用的是 block-aware globalization 修复，而不是容差放宽。

## 31 点 G3 exact-point 结果

新的合同已物化到 `prepare_templates_ldmos_phase23.py` 的全部 G3 decks。
Sentaurus seed reclose、Save/Load repeat 和无 predictor 31 点 Id-Vg 均完成：

- 31/31 点收敛；总 Newton iteration 290，单点最大 21；
- 最大 electron block residual：`9.40516e-12`；
- 最大 psi block residual：`1.10649e-9`；
- 最大 relative KCL error：`0.00103787`（`0.103787%`），通过 `1%` 门；
- `Vg=1/6 V` relative KCL error：`5.23801e-8`；
- resolved points 从 30 增至 31。

曲线结果保持物理差异不变：median log error `0.00666064 dex`，P95
`0.287023 dex`，固定电流 Vth error `29.6262 mV`，强反型端点相对误差
`1.45137%`。除 P95 `0.20 dex` 门外，其余 L2 指标均通过。因此 WP1.5 数值
阻塞已关闭，但阶段 3 总状态仍为 fail；剩余问题首先是已独立复现的接触 HFS
支撑语义，其次才是 G4 自洽状态下仍存在的 coefficient/current-support 差异；
不能再归因于低电流 KCL 或 line-search 地板。

本轮 ignored 资格目录：

`reference_staging/templates_ldmos_sentaurus2022/phase01_original_20260826_02/stage1_v4/phase23_t2022_contract_v6_wp15_block_filter/`。

## 验证

- `test_newton_solver`：102 test cases、1340 assertions passed；
- `test_sg_flux`：28 test cases、227 assertions passed；
- `test_mobility`：30 test cases、142 assertions passed；
- `python -m unittest tests.regression.test_templates_ldmos_phase23`：21 tests passed；
- Release 全量 CTest：723/723 passed；
- exact-point qualification：31/31 converged，KCL gate pass，P95 gate fail。
