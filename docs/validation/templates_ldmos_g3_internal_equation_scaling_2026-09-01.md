# Templates/LDMOS G3 Sentaurus 内部方程缩放审计

## 结论

本轮关闭了 `33.7806012x` 标量的三个候选归因：Electron `ErrRef`、逐行
AverageBox `Measure` 归一化，以及线性求解器预条件。简单的 `1/Vt` 解释也不
成立。Sentaurus NewtonPlot `eContinuityRhs` 的绝对幅值仍只能归类为未公开的
统一方程归一化或统一 assembly coefficient；它不能作为跨引擎绝对电流门，
但空间形状、偏压相对变化和受控扰动响应仍可用于诊断。ledger 继续保持
`draft`。

最强的直接证据是：在完全相同的高端点 VSV `±1 uV` 状态上，仅把
`ErrRef(Electron)` 从 `1e10` 改为 `1e8 cm^-3`，iteration 0/1 的完整
`eContinuityRhs` 和 iteration 1 的 `NewtonStepEDensityUpdate` 均逐字节不变。
只有 `eDensityError` 改变，而且 5674 个有效节点的候选/基线比例以最大
`4.44e-16` 相对误差精确服从
`(|n|+1e10)/(|n|+1e8)`。

## 手册合同

T-2022.03-SP2 *Sentaurus Device User Guide* 第 192--194 页把 `|Rhs|` 定义为
残差范数，把 `ErrRef` 放在相对更新误差
`Delta x / [epsilon_R (|x| + x_ref)]` 的分母中；Electron `ErrRef` 的单位是
`cm^-3`。第 208--210 页说明 CNorm/NewtonPlot 暴露的是内部且
implementation-dependent 的 RHS、error 和 update 数据。手册没有给出普通
drift-diffusion 方程的公开逐行缩放公式。本轮全文检索到的 `RhsScale` 仅属于
MDFT/SDFT 模式，不适用于该 G3 经典 DD deck。

这给出可证伪预测：若 `ErrRef` 是 `33.78x` 的来源，修改它必须改变
`eContinuityRhs`；若它只属于更新误差，RHS 和 Newton 增量应保持不变，而
`eDensityError` 应按上式变化。实验明确支持后者。

## ErrRef 单因素 A/B

| 合同项 | 基线 | 候选 |
| --- | ---: | ---: |
| 状态 | high-endpoint VSV，node 3721 `±1 uV` | 同左 |
| Physics、网格、Math/Solve 其余字段 | 冻结 | 同左 |
| `ErrRef(Electron)` | `1e10 cm^-3` | `1e8 cm^-3` |
| production default | 不修改 | 不修改 |

两个符号分支结果一致：

| 观测 | minus | plus |
| --- | ---: | ---: |
| iteration-0 `eContinuityRhs` CSV | 同哈希 | 同哈希 |
| iteration-1 `eContinuityRhs` CSV | 同哈希 | 同哈希 |
| iteration-1 electron Newton update | 同哈希 | 同哈希 |
| iteration-1 `eDensityError` CSV | 不同 | 不同 |
| error 比例公式最大相对误差 | `4.44e-16` | `4.44e-16` |
| CNorm 最大节点 | `3867 -> 4227` | `3867 -> 4227` |
| iteration-0 `|Rhs|` | `7.28e5 -> 7.28e5` | `1.26e6 -> 1.26e6` |
| iteration-1 `|Rhs|` | `4.72e9 -> 4.72e9` | `4.72e9 -> 4.72e9` |

因此 `ErrRef` 是更新误差/停止判据合同，不是 continuity RHS 的方程幅值合同。

## 控制盒 Measure

从冻结的 `MeasureCoefficients.debug` 按 `(element, region)` 重建五节点硅区
AverageBox Measure，并与前轮 node-3721 中央差分的逐行 Sentaurus/Vela
比例对齐：

| node | silicon Measure (`um^2`) | `abs(Sentaurus/Vela)` |
| ---: | ---: | ---: |
| 3714 | `5.97170e-5` | `33.8308` |
| 3720 | `4.08821e-6` | `33.7635` |
| 3721 | `1.06908e-5` | `33.7835` |
| 3722 | `2.92776e-5` | `33.7090` |
| 3723 | `7.42921e-6` | `33.8761` |

Measure 跨 `14.6071x`，而绝对行比例只跨 `1.00496x`；二者 log-space Pearson
相关系数为 `-0.07625`。若 NewtonPlot 在这些行上逐行乘或除 Measure，比例也
应随 Measure 显著变化。观测到的近统一列标量排除了这一简单语义。该结论不
否定 Measure 参与物理有限体积装配，只否定它作为 `33.78x` 的逐行输出
归一化因子。

## 行范数、预条件与输出顺序

四份日志都给出同一事件顺序：

1. 打印 iteration-0 `|Rhs|`；
2. 写出 `newton_0_0.tdr`；
3. 打印 CNorm update error；
4. 执行第一轮线性求解并打印 iteration 1（含 `#inner/#iterative`）。

日志同时明确写出 `Without diagonal preconditioning`。因此线性求解器的
预条件不可能倒因果地改变 iteration-0 NewtonPlot RHS。本轮没有为一个已被
事件顺序排除的因素再运行 Method/ILS A/B。

此外，minus/plus 的导出 electron RHS L2 分别为 `2.65426e-5` 和
`2.76339e-5`，而日志全局 `|Rhs|` 为 `7.28e5` 和 `1.26e6`，数值相差十个
数量级以上。日志范数不是导出 electron field 的直接 L2；NewtonPlot 输出前
确实还存在未公开的变量/方程表示转换。这个事实支持“内部归一化”，但现有
公开资料仍不能推导其精确 `33.7806x` 系数。

## 热电压候选与 Vela 量纲链

300 K 下 `1/Vt = 38.6817 V^-1`，比观测标量高 `12.6704%`；
`33.7806012 * Vt = 0.873296`，不是单位值，因此简单漏乘/多乘一个 `Vt`
也被排除。Vela 源码明确以 `V0=kT/q`、`D0=mu0*V0` 构造缩放，并将 carrier
residual 除以 `C0*D0`；SG 工件再以逐边完全一致的物理 particle-flux 比例
恢复 `4.416428645843634e-5 A/um` 每 scaled residual。故不能为了追随
NewtonPlot 的内部数字再人为调整 Vela 的物理量纲链。

## 冻结解释边界

- 已排除：二维宽度/AreaFactor、Electron ErrRef、逐行 silicon Measure、
  线性求解器预条件、简单 `1/Vt`。
- 仍开放：Sentaurus 未公开的统一非线性方程归一化，或尚无 assembly-level
  导出可验证的统一系数。
- 可继续使用：NewtonPlot RHS 的局部支持、方向、归一化后的空间形状、偏压
  相对增长，以及同一 Sentaurus 合同内的 A/B。
- 禁止使用：把 NewtonPlot `eContinuityRhs` 标签 `A` 当作端口 `A` 或 Vela
  `A/um`，以其绝对幅值建立跨引擎验收门。
- 不启用 IALMob，不调整阈值/迁移率参数，不修改 production default；31 点
  曲线仍不因本轮内部标量而重跑。

若后续要继续追查绝对标量，唯一有辨识力的证据是受支持的 assembly-level
element-edge current/row-scale 导出或 Synopsys 对该版本内部格式的明确说明；
重复改变收敛参数、线性 Method 或器件 AreaFactor 不再具有诊断价值。

## 可复跑工件

- 审计脚本：`scripts/audit_templates_ldmos_g3_internal_scaling.py`；
- 回归测试：`tests/regression/test_audit_templates_ldmos_g3_internal_scaling.py`；
- ignored oracle：`reference_staging/templates_ldmos_g3_internal_scaling_20260901/`；
- 汇总：上述目录的 `analysis/summary.json`；
- Measure/行比例表：上述目录的 `analysis/measure_row_scale.csv`。
