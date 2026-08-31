# Templates/LDMOS G3 Si/SiO2 界面主从节点最小验证

日期：2026-08-31  
结论：**局部假设成立，但整体硬门失败；不启动同偏压 reclose，不修改生产求解器。**

## 验证问题与边界

本轮只验证一个最小、可证伪的问题：Vela 把几何上共享的 Si/SiO2 节点作为
一个全局节点，并在 carrier continuity 中使用相邻所有材料单元累计得到的
`edge.couple`；如果把共享节点概念性拆为：

- semiconductor master：`psi/phin/phip`；
- oxide slave：仅 `psi`；
- 约束：`psi_slave - psi_master = 0`；
- carrier edge coupling：只保留半导体区域的局部 box contribution；

是否能在完全相同的 Sentaurus `Vg=0.5 V` 冻结状态上，将电子 continuity
残差至少改善 2 倍。

这不是 Sentaurus 内部实现的逆向声明，也没有真的扩展 Vela 自由度。probe
只在离线重放中建立概念性重合节点对，保持 Poisson 不变，并把每条 SG 电子
通量按 `semiconductor_couple / all_material_couple` 重标定。

## 输入与资格自检

- exact SProcess topology：10,241 个全局节点；
- 冻结状态：G3 no-IALMob、`Vd=0.1 V`、`Vg=0.5 V`、SSS；
- conceptual Si/oxide pairs：762；概念性 region-local 总节点数 11,003；
- Si/oxide 邻接输运边：760；
- 自由电子 continuity 行：5,674，其中共享界面行 759；
- 全局 `edge.couple` Python/C++ 重现最大相对误差：`4.37e-16`；
- 自由行逐边 SG 通量重装配最大绝对误差：`2.51e-22`，最大相对误差
  `4.75e-8`。

因此结果不是 couple 公式、长度单位、边方向或接触 Dirichlet 行混入造成的。
Dirichlet/纯氧化层载流子行不参与残差评分，且保持原边界残差不变。

## 固定状态结果

| 指标 | baseline | region-local candidate | candidate / baseline | 2x 硬门 |
| --- | ---: | ---: | ---: | :---: |
| electron residual L1 | `3.94960e-5` | `2.27109e-5` | `0.575018` | 报告项 |
| electron residual L2 | `2.42852e-6` | `1.57776e-6` | `0.649679` | 失败 |
| maximum absolute row | `6.04064e-7` | `5.40895e-7` | `0.895426` | 失败 |

冻结热点逐点结果：

| node | 界面共享 | baseline residual | candidate residual | 绝对值比 |
| ---: | :---: | ---: | ---: | ---: |
| 3747 | 是 | `6.04064e-7` | `-1.23418e-11` | `2.04e-5` |
| 10233 | 是 | `-5.77592e-7` | `9.07185e-8` | `0.157062` |
| 4538 | 是 | `-5.10737e-7` | `2.74112e-10` | `5.37e-4` |
| 4492 | 否 | `-5.40895e-7` | `-5.40895e-7` | `1.0` |

三个共享界面热点全部超过 2 倍改善，证明当前全局 carrier couple 确实把氧化层
一侧 box contribution 带入了这些界面行，而且该语义对局部残差有实质影响。
但是 node 4492 是仅属于硅区的非界面节点，候选算子对它严格无影响；它成为新
的最大残差行。其余非界面残差也使整体 L2 只能降至 `0.649679x`。

## 停止决定

冻结硬门要求：L2 比值不大于 0.5、最大行比值不大于 0.5、每个自由界面热点
比值不大于 0.5。第三项通过，前两项失败，所以：

1. 不运行成本更高的同偏压自洽 reclose 和 31 点曲线；
2. 不把 region-local carrier couple 作为生产 profile，也不改变 Poisson/自由度；
3. 不把现有 P95 差异登记为已批准的引擎固有地板；
4. 保留结论：界面主从/重合节点是一个已证实的局部缺口，但不是当前
   `2.61--2.73x` 自洽电流平台的充分根因。

## 下一最小审计

下一步应保持 IALMob 与 predictor 关闭，针对 node 4492 及其一环硅区单元，
对齐 Sentaurus AverageBox 的 region-local node volume、逐边
`CoeffIntersection/couple`、电子密度平均和第一步 Newton 增量。若无法取得
`MeasureCoefficients.debug`，只能把候选规则保持为实验性 A/B，不得宣称与
Sentaurus 等价。界面主从自由度的生产实现应等该非界面硬门关闭后再评估，
避免为局部改善引入新的全局方程结构。

## 可复现工件

- probe：`scripts/audit_templates_ldmos_interface_pair_box.py`；
- 回归测试：`tests/regression/test_audit_templates_ldmos_interface_pair_box.py`；
- ignored 输出：
  `reference_staging/templates_ldmos_g3_interface_pair_box_20260831/summary.json`、
  `interface_pairs.csv`、`interface_edges.csv`、`node_residuals.csv`。
