# Templates/LDMOS WP2、阶段 2–3 执行报告（2026-08-27）

## 结论

本轮完成了 WP2 所需的 T-2022.03 参数合同、PolySi gate 控制、G4 经典平衡态
oracle、G0--G4 Id-Vg 单因素链和 Vela 生产算子固定状态回放。后续 WP1.5 重新资格
通过 Poisson bootstrap + 近冻结 QF coupled merit 关闭了 0 V 同偏压 reclose；同时
修复了 gate sweep 丢失 PolySi flatband 的边界映射缺陷。但 G3 的 0.1 V 状态仍不能
通过严格 repeat，因此阶段 3 的 Vela 自洽曲线未启动。详细证据见
`templates_ldmos_wp15_requalification_2026-08-27.md`。

- 当前最高正式等级：`L1`；
- WP2：部分完成，`stop_rule_triggered=true`；
- 阶段 2：0 V 平衡态/重启子门通过，偏压入口未关闭；
- 阶段 3：Sentaurus 消融完成，Vela 曲线未运行、未评分；
- 不认领 `L2`、`L3`、A+ 或方案 B。

## 已执行范围

### T-2022.03 显式合同

通过虚拟机中的 `sdevice -P` 查询并审计 T-2022.03-SP2 Silicon 参数，材料与物理
合同升级为 revision 2：

- 300 K `Eg=1.1241592307692307 eV`、`Nc=2.856665228e19 cm^-3`、
  `Nv=3.104641908e19 cm^-3`；
- OldSlotboom `Ebgn=0.009 eV`、`Nref=1e17 cm^-3`、`C=0.5`；
- `dEg0(OldSlotboom)=-0.01595 eV` 按既有 split-ni 合同冻结进材料 `ni`，不在
  BGN 中重复计数；
- DopingDependence Formula 1 完整映射为 Vela `masetti`，不再用简化
  Caughey-Thomas 冒充；
- SRH 掺杂与温度依赖已显式启用，并写入电子/空穴 `tau_min/tau_max/Nref/gamma`；
- HighFieldDependence 的 `vsat0/beta0` 显式写入；
- Sentaurus Auger 的 300 K `A+B+C` 冻结为 Vela 标量系数。由于 Vela 尚不支持
  `H/N0` 密度增强，合同将重组复合标为 `partial`。

对应 schema 现在要求上述带单位字段且继续拒绝未知键。

### PolySi gate 控制

从封存等势 gate 的 453 个节点得到：

- Sentaurus gate 电势：`+0.5609636105202798 V`；
- 节点间 spread：`0 V`；
- Vela `psi_gate = bias - flatband_voltage`，故
  `flatband_voltage=-0.5609636105202798 V`。

Sentaurus `Material="PolySi"(N)` 与显式
`Barrier=-0.5609636105202798 V` 的 31 个 Id-Vg 点逐点完全一致，电流 log 差、
强反型相对差和三个诊断 Vth 交点差均为 0。这关闭了 gate 功函数映射，而不是用
端电流拟合 flatband。

### G4 经典平衡态 oracle

从 `G4-no-highfield` 只派生输出/路径控制 deck：Physics 段不变，在原始 Poisson 和
Poisson/Electron/Hole 初始耦合后停止，并保存非 loadable TDR。该运行成功，生成
10241 节点的 exact-topology 经典平衡态；此状态已转换为 17 位有效数字 Vela restart
CSV，可用于阶段 2 空间评分。

## Sentaurus Id-Vg 单因素结果

六条曲线均为 31 个相同 CurrentPlot 点，不对曲线分数做插值。单条墙钟为
`287.65--403.91 s`。

| 单因素 | 父/子层 | 中位 log 增量 | 强反型相对增量（中位 / 最大） |
| --- | --- | ---: | ---: |
| hQP | G0 / G1 | `6.31e-14 dex` | `1.55e-11% / 3.42e-11%` |
| eQP | G1 / G2 | `5.06e-14 dex` | `2.22e-11% / 4.07e-11%` |
| IALMob | G2 / G3 | `0.2761 dex` | `83.69% / 84.34%` |
| HighFieldSaturation | G3 / G4 | `0.00613 dex` | `1.261% / 1.274%` |
| PolySi/barrier | G0 / control | `0 dex` | `0% / 0%` |

原 deck 虽在 Physics 中声明 e/h QuantumPotential，但 Solve 序列只耦合
Poisson/Electron/Hole，没有量子势方程；消融结果证明两项在当前官方执行路径中低于
数值分辨。hQP 强反型判据明显低于 2%，但 `threshold_freeze.json` 仍未指定并批准
固定电流 Vth 的具体电流水平，因此本轮只报告 `1e-12/1e-10/1e-8 A/um` 的诊断
交点，不能正式签署 `hqp_secondary`。

## Vela 阶段 2 首次结果与后续重新资格

完整 revision-2 合同从 G4 经典状态启动时：

- 0 V mapped equilibrium 在第 9 次更新后失败；
- 原因：`line_search_non_decrease`；
- 最终势块残差：`2.2744539146463174e5`；
- 电子/空穴连续性残差分别约 `7.37e-13`、`9.44e-13`；
- 失败集中在势块，最高残差节点位于约 `2.9e20--4.0e20 cm^-3` 的重掺杂接触邻域；
- 由于首个 fixed point 未关闭，同偏压 repeat、0→0.1 V drain ramp 和完整 G3 Id-Vg
  均按停止规则未运行。

作为独立证据，生产探针在相同 G4 状态上成功完成：

- `edge_mobility_probe`：30022 edges，T-2022.03 Masetti 参数；
- `sg_edge_flux_probe`：30022 edges；
- `newton_carrier_term_probe`：10241 rows；
- 固定状态公式重建相对 Sentaurus 密度的中位误差为电子 `0.1203 dex`、空穴
  `0.1269 dex`，P95 分别为 `0.4330 dex`、`0.1290 dex`。

后续 WP1.5 已关闭 0 V 子门：Poisson bootstrap 后的 coupled/repeat 在 10241 个节点
上除 `2.12e-30 V` 的电子 QF 舍入差外完全一致。执行 gate 扫描 repeat 时还发现并
修复了 DCSweep 覆盖 metal-gate flatband 的缺陷。修复后，0→0.1 V 首程可收敛，
但严格 repeat 在 `3.70e-7` Poisson 块和 `0.826 V` 接触多数载流子 QF 跃迁处失败；
`contact_basin` 与关闭 continuity row scaling 的消融均未关闭。因此阶段 2 只通过
0 V 平衡态/重启子门，偏压入口仍未关闭。

补充 Sentaurus 重闭合探针确认 T-2022.03-SP2 在相同 G3 exact mesh 上的同进程和
Save/Load 重闭合均只需一次 Newton 更新，Save/Load 前后势/QF 差处于 double 舍入
量级。Vela 复刻 Sentaurus 实际 15 点漏压路径后仍在
`0.023155977422 -> 0.031641657941 V` 转移失败，因此 restart 文件精度和单纯粗步长
均已排除。详见 `templates_ldmos_sentaurus_g3_reclose_probe_2026-08-27.md`。

## 后续开发项

1. 继续 WP1.5 的有偏压 exact-mesh 资格，使 continuity row scaling 与收敛范数在
   continuation/restart 间保持固定点不变，并阻止接触 QF 不安全的 relative convergence。
2. 对 split-ni、constant-reference potential、OldSlotboom/Fermi 和 TDR potential/QF
   参考零点做逐节点公式审计，关闭约 `0.12 dex` 的固定状态密度中位差。
3. 0.1 V strict repeat 通过后才运行 G3 exact-point Id-Vg。IALMob 是当前最大的已证实功能缺口，
   应以 G2-G3 的约 84% 强反型增量为 WP 优先级依据；不得先调 bulk mobility 拟合。
4. 在 hQP 决策前补齐并双签固定电流 Vth 电流水平。若继续沿用官方 Solve 序列，应把
   “QP 已声明但量子方程未耦合”写入 physics contract 和差异账本。

## 可复现产物

专有 TDR/PLT、运行日志和大体量 CSV 保存在 ignored
`reference_staging/.../stage1_v4/phase23_t2022_contract_v2/`。仓库只提交生成器、
合同、schema、回归测试和本报告，不提交生成仿真输出。
