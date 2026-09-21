# Templates/LDMOS G3 node 4492 AverageBox 审计（2026-08-31）

## 结论

在 IALMob 与 predictor 继续关闭、`Vd=0.1 V`、`Vg=0.5 V` 的同一 G3
Sentaurus 状态上，已取得 T-2022.03-SP2 直接生成的
`MeasureCoefficients.debug`，并完成非界面硅节点 4492 的逐单元、逐边、密度平均
和第一次 Newton 更新对齐。

直接 AverageBox 系数回放把 node 4492 的电子 continuity 残差从
`-5.40895e-7` 降到 `-1.31409e-10`，绝对值改善约 4,116 倍；但一环残差 L2
只从 `6.93948e-7` 降到 `6.17998e-7`，即 `0.890554x`，未通过冻结的 `0.5x`
整体门。误差主要转移到 node 10235。因此，**4492 本点的根因已定位为
AverageBox 逐边系数语义，但局部替换不是可批准的全局修复**；不启用生产
AverageBox profile，不运行 31 点曲线，known-difference ledger 继续保持 draft。

## 直接 oracle 与局部编号资格

- exact SProcess topology：10,241 个节点；node 4492 仅属于 Silicon_1；
- 一环节点：4487、4488、4491、4492、4493、4529、4530、10235；
- 相邻硅单元：8433、8436、8437、8438、8514、10507、10508；
- Sentaurus 版本：T-2022.03-SP2；
- `Math`：`AverageBoxMethod`、`BoxMeasureFromFile(GrdNumbering)`、
  `BoxCoefficientsFromFile(GrdNumbering)`；
- 直接文件 SHA-256：
  `05ec1a43936e0ff71bea7cdfd956d80c7d692b228187ccae8625ec53ce58d49e`。

T-2022.03-SP2 用户手册的式 (1266) 明确：离散按 element-vertex 装配，二维系数
为 `kappa_ij = d_ij/l_ij`。debug 文件的 `Measure` 三值按局部顶点排列，而三角形
`Coefficients` 遵循 TDR 边序 `[v2-v0, v1-v2, v0-v1]`。本轮首先错误地按“对顶
边”解释，得到约 18 倍恶化；用普通直角/锐角单元的解析 cotangent 系数资格后
纠正边序。纠正后的 Vela couple 重装配与 C++ 输出最大绝对误差为 0，且边序已
加入回归测试，错误的 18 倍结果不作为证据保留。

## 逐边系数

下表为 node 4492 的七条 incident transport edge；比值为
`Sentaurus AverageBox couple / Vela production couple`。

| edge | Vela couple (m) | AverageBox couple (m) | 比值 |
| --- | ---: | ---: | ---: |
| 4487--4492 | `1.61133e-8` | `1.57156e-8` | `0.975318` |
| 4488--4492 | `3.85550e-9` | `3.85550e-9` | `1.000000` |
| 4491--4492 | `1.51140e-8` | `1.51140e-8` | `1.000000` |
| 4492--4493 | `0` | `0` | 同为零 |
| 4492--4529 | `3.02087e-9` | `1.51044e-9` | `0.500000` |
| 4492--4530 | `3.90914e-10` | `0` | `0` |
| 4492--10235 | `1.11043e-7` | `5.74517e-8` | `0.517385` |

主要闭合来自 4492--10235 的约半权重以及 4492--4529/4530 的截断，并非统一
比例或端口积分误差。node 4492 的 AverageBox measure 为
`8.62181057e-5 um^2`。

## 密度平均与固定状态残差

新增的只读 SG probe 输出“单位 couple 响应”，包括生产 couple 为零的边；它不
修改装配或求解默认。由粒子通量、迁移率和 QF 梯度反推的有效 SG 密度，与端点
均值的中位相对误差为：

| 端点密度平均 | 中位相对误差 |
| --- | ---: |
| arithmetic | `0.287869` |
| geometric | `0.0677988` |
| harmonic | `0.114666` |
| logarithmic | `0.140194` |

这只用于分类现有 SG 权重，不能据此宣称 Sentaurus 使用 geometric mean。保持
当前 SG 密度语义、仅替换直接 AverageBox couple 时：

| 指标 | baseline | AverageBox replay | 比值 | `0.5x` 门 |
| --- | ---: | ---: | ---: | :---: |
| node 4492 `|electron residual|` | `5.40895e-7` | `1.31409e-10` | `0.000242948` | 通过 |
| 一环 electron residual L2 | `6.93948e-7` | `6.17998e-7` | `0.890554` | 失败 |
| 一环 maximum absolute row | `5.40895e-7` | `4.54606e-7` | `0.840461` | 失败 |

一环最大候选残差转移到 node 10235：baseline `7.28178e-8`，候选
`-4.54606e-7`。因此不能把仅围绕 4492 的系数替换写入生产求解器。

## 第一次 Newton 更新对齐

为避免把不可重启的 `Plot(-Loadable)` TDR 当作 Save，本轮在独立远端目录中：

1. 按原 G3 合同重新扫到 `Vd=0.1 V, Vg=0.5 V`，生成包含 SLP 的原生 `.sav`；
2. 在第二个 SDevice 进程中加载该 Save，启用 AverageBox，只执行一次
   `Poisson Electron Hole` Newton；
3. 导出 `NewtonStep*Update`、Poisson/e/h RHS 和更新后 QF 场。

Sentaurus 该状态的 drain electron current 为 `6.462e-11 A`。node 4492 对齐为：

| 量 | Sentaurus AverageBox | Vela production |
| --- | ---: | ---: |
| `delta psi` (V) | `-3.63629e-17` | `1.70682e-3` |
| `delta phin` (V) | 双精度下 `0` | `1.96514e-11` |
| electron density update (`cm^-3`) | `3.80938e4` | 不同未知量，不直接比较 |
| post electron density (`cm^-3`) | `2.24326e19` | 状态来自同一 Sentaurus oracle |

Sentaurus 电子密度相对更新仅 `1.70e-15`，一环更新后的 `psi/phin/phip` 最大差均
处于 `0--1.84e-16 V` 量级，说明其状态对自身 AverageBox 算子已经闭合。Vela
同状态的毫伏级 Poisson 更新和非零 continuity 更新不是状态文件精度噪声。

## 决策与下一步

1. node 4492 的“状态/密度错误”假设关闭；直接证据支持 Box 系数语义差。
2. 局部一环实现门失败，不增加生产配置，不启动同偏压 reclose 或 31 点曲线。
3. 后续全网格只读实验已经完成并通过：全自由电子行 L2 和 maximum 分别降到
   baseline 的 `0.0173075x` 和 `0.0433710x`，没有实质热点转移。完整结论见
   `templates_ldmos_g3_averagebox_full_mesh_audit_2026-08-31.md`。下一步才实现可选
   的全网格 external AverageBox 诊断 profile 和真正的 Newton A/B。
4. 不继承 PN2D 私有 `element_edge_sg_gss_laux` 捆绑；IALMob 与 predictor 继续
   关闭；ledger 仍为 draft。

## 可复现工件

- oracle 运行器：`scripts/run_templates_ldmos_averagebox_probe.py`；
- node 4492 审计：`scripts/audit_templates_ldmos_averagebox_node4492.py`；
- 回归测试：
  `tests/regression/test_run_templates_ldmos_averagebox_probe.py`、
  `tests/regression/test_audit_templates_ldmos_averagebox_node4492.py`；
- ignored 大型输出：
  `reference_staging/templates_ldmos_averagebox_node4492_20260831/`、
  `reference_staging/templates_ldmos_averagebox_node4492_newton1_capture_20260831/`、
  `reference_staging/templates_ldmos_averagebox_node4492_audit_20260831/`。
