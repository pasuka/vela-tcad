# Templates/LDMOS G3 亚阈值偏移与低电流 KCL 审计

日期：2026-08-29

状态：根因已定位；阶段 3 仍未通过 L2 硬门；形成候选已知离散差异，等待审批。

## 范围与不变量

本审计承接 31 点 G3 `Id-Vg` 资格运行，保持以下合同不变：

- T-2022.03-SP2 Templates/LDMOS exact 10,241-node topology；
- `Vd=0.1 V`，七个精确公共门偏压点 `0, 1/6, ..., 1 V`；
- Fermi-Dirac、OldSlotboom BGN、SRH/Auger、G3 HighFieldSaturation；
- IALMob 关闭，predictor/continuation 关闭；
- SG edge current、barycentric node volume、正 barycentric 负-cotangent 回退；
- 不引入阈值平移、迁移率标定、量子、自热或雪崩项。

基线 31 点运行已经全部收敛，但 Sentaurus G3 曲线的 resolved P95
log 误差为 `0.287219 dex`（门限 `0.20 dex`），`1e-8 A/um` 固定电流
阈值偏移为 `29.6261 mV`，最低 resolved 点 `Vg=1/6 V` 的四端 KCL
相对误差为 `6.16280%`（门限 `1%`）。

## Sentaurus 精确状态真值

在 Sentaurus VM 上新增只读诊断 deck，复用封存的 `n1_fps.tdr` 与
`sdevice.par`，重新运行 G3-no-IALMob 并在七个公共偏压保存 TDR。每个 TDR
导入 ElectrostaticPotential、electron/hole QF、density、mobility、current
density 等字段。生成物封存在 ignored 路径：

`reference_staging/templates_ldmos_g3_shift_kcl_20260829/`。

88 个门响应 Si/oxide 界面节点的 Vela-Sentaurus 差异为：

| Vg (V) | psi 中位差 (mV) | psi P95 绝对差 (mV) | electron density 中位差 (dex) |
| ---: | ---: | ---: | ---: |
| 0.0000 | 2.78 | 11.26 | 0.047 |
| 0.1667 | 1.87 | 12.23 | 0.032 |
| 0.3333 | -0.41 | 13.32 | -0.006 |
| 0.5000 | -2.77 | 14.02 | -0.046 |
| 0.6667 | -4.56 | 14.02 | -0.074 |
| 0.8333 | -6.01 | 12.28 | -0.085 |
| 1.0000 | -8.99 | 13.44 | -0.106 |

由界面 electron-density 斜率换算的等效门压偏移中位数为 `-11.47 mV`，
P95 绝对值为 `96.40 mV`，且偏差随门压改变符号。因此终端曲线的
`29.6 mV` 不是一个可用固定 psi、flatband 或能级参考修正关闭的均匀
电势偏移。

## 低电流 KCL

`Vg=1/6 V` 时：

- drain current：`1.19460e-13 A/um`；
- four-terminal KCL residual：`-7.36205e-15 A/um`；
- native、compensated 和 long-double 四端求和的最大差仅
  `6.56e-28 A/um`；
- electron continuity residual signed sum 为 `-1.66697e-10`，absolute
  sum 为 `2.28025e-9`；
- 原资格点以 electron block residual `2.86022e-10` 通过现有
  `2e-9` 绝对门。

将 electron block ceiling 单因素收紧到 `1e-11` 后，Newton 只接受两个
更小步长；electron residual 停在 `2.85437e-10`，随后发生
`block_absolute_convergence_line_search_rejected`。因此 6.16% KCL 不是端口
浮点求和误差，也不能由单纯收紧停止门关闭；它是当前 exact-mesh 低电流
Newton 地板的可观测后果。KCL 硬门仍失败，不能豁免。

## 同状态电流算子回放

生产 `sg_edge_flux_probe` 在三个代表 Sentaurus 状态上回放。首先用 Vela
自身状态验证 drain-cut：三个偏压的 SG cut 与生产端口电流最大相对差仅
`3.47e-15`，证明积分方向、单位和接触边集合正确。Sentaurus 状态在 drain
边的 SG 重构 electron density 与 TDR supplied density 最大差小于
`4.94e-4 dex`，排除 density mapping。

| Vg (V) | Sentaurus terminal | native Jn 边界积分 | cell QF-gradient 重构 | Vela edge-SG 同状态 | SG / terminal |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.1667 | 6.7704e-14 | -5.4415e-14 | 2.7133e-13 | 1.6586e-12 | 24.50 |
| 0.5000 | 6.4622e-11 | -4.9230e-11 | 1.1154e-10 | 7.1584e-10 | 11.08 |
| 0.8333 | 9.8895e-8 | -7.5337e-8 | 1.7047e-7 | 1.0942e-6 | 11.06 |

电流单位均为 `A/um`；native Jn 的符号来自导入接触边定向，比较使用幅值。
在后两个偏压，native Jn 边界积分约为 terminal 的 `0.76x`，简化 cell
QF-gradient 重构约为 `1.72x`，均保持原生量级；edge-SG 则约为 `11.07x`。
从 cell gradient 到 edge-SG 额外放大约 `6.4x`。

Sentaurus 状态在 edge-SG 下的主导 drain 边具有 `couple/length=18--83`，
而 Vela 自洽状态的主导边约为 `6--9`；两者 drain 邻域 electron mobility
均约为 `0.1417 m2/(V s)`。关闭负局部 cotangent 的正 barycentric 回退后，
中点误差只从 `1.04443` 变为 `1.04398 dex`，说明 1,042 个负局部
cotangent 回退不是主因，差异集中在细长、强正耦合边上的 edge-only
current support。

## 结论与处置

1. `29.6261 mV` 是终端曲线对 exact-mesh 电流离散差异的等效表示，不是
   需要修改 flatband、BGN、Fermi 或 mobility 参数的物理阈值偏移。
2. G3 HighFieldSaturation 关闭控制没有改善亚阈值误差；IALMob 和 predictor
   继续保持关闭。
3. 当前 SG+barycentric profile 在该极端非正交网格上可完成 31 点运行并在
   强反型收敛，但亚阈值 P95 误差构成候选引擎固有离散地板。
4. 阶段 3 仍为 **fail**：在 independent reviewer 批准 known-difference
   ledger 前，不修改 `0.20 dex` P95 和 `1%` KCL 硬门。

后续开发应拆成两个独立任务：

- 为一般 Tri3/非正交网格开发并资格验证 element/cell-aware carrier current
  support，使同状态 current replay 接近 Sentaurus native Jn；这属于新的
  discretization profile，不得静默替换现有 SG profile；
- 继续 WP1.5 的低电流 electron-continuity line-search/残差地板工作，并增加
  可选 terminal-KCL acceptance 诊断；不得用放宽 KCL 门或 predictor 掩盖。

## 实现与验证

- 新增 `audit_templates_ldmos_g3_idvg_shift_kcl.py`，自动完成状态比较、KCL
  精度审计、SG 自状态闭环、Sentaurus 同状态回放、native/current-gradient
  重构和单因素几何控制；
- 暴露默认兼容的
  `mesh_geometry.fallback_negative_cotangent` 资格开关；`false` 表示把负局部
  cotangent 贡献置零，并非保留有符号负 coupling；
- 新增 Python 合成边界积分/单元梯度测试，并扩展 box-geometry 配置解析测试。

已执行：

- `python -m unittest tests.regression.test_templates_ldmos_g3_shift_kcl_audit tests.regression.test_templates_ldmos_phase23 tests.regression.test_templates_ldmos_phase01`：47 tests passed；
- `build/test_box_geometry.exe`：10 test cases、65 assertions passed；
- UCRT64 Debug 全量构建：96 targets built。
