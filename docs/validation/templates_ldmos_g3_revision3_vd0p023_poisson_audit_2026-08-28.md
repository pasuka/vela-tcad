# Templates/LDMOS G3 revision 3：0.023156 V checkpoint 与 Poisson 分项审计

## 结论

revision 3 已按 Sentaurus 接受的无 predictor 精确偏压前缀保存
`Vd=0.0231559774221138 V` checkpoint，9/9 点全部由 Vela 标记为收敛。
该状态可作为诊断 checkpoint，但尚不能作为 G3 曲线验收状态：终点电子 continuity
残差 L2 为 `1.0300e-4`，有 334 个 QF bounds violation，漏极电流为
`4.5484e-9 A/um`，相对 Sentaurus 同点固定状态回放电流
`-9.8292e-14 A/um` 仍相差约 4.66 dex。

Poisson 五项审计逐节点闭合，排除了审计公式和 Poisson 边通量漏项。主要差异位于
Sentaurus 状态进入 Vela 后的载流子状态重构：Sentaurus CSV 中保存的 `n/p` 与 Vela
根据同一 `psi/phin/phip`、当前 Fermi/OldSlotboom 合同重构的 `n/p` 中位误差均约
`0.108 dex`，最大分别为 `1.071 dex` 和 `1.211 dex`。因此当前证据把首要问题从
“Poisson 离散/缩放”收窄为“Fermi + BGN 的电势—准费米势—载流子密度映射合同”，
并同时保留重掺杂接触边界势映射为次要待审项。

## 1. 范围与不变量

- 物理层：G3 classical bulk + high-field，revision 3 mobility 合同；
- 偏压：`Vg=0 V`，`Vd=0.0231559774221138 V`；
- continuation：严格使用 Sentaurus 已接受的前 9 个偏压点；
- predictor：关闭；
- 本轮未引入 IALMob、QP、热、Okuto，也未调整 Newton 门槛；
- Sentaurus 与 Vela 状态均通过同一个 Vela 生产装配器回放。

## 2. checkpoint 资格记录

精确偏压序列为：

`0, 0.001, 0.00214466666666667, 0.00366859955555556,`
`0.00569746220829630, 0.00839855468664514, 0.0119946091394869,`
`0.0167821563010369, 0.0231559774221138 V`。

终点记录：

| 指标 | revision 3 checkpoint |
|---|---:|
| 扫描点 | 9/9 |
| Newton 迭代 | 30 |
| Poisson residual L2 | `2.1557e-8` |
| electron continuity residual L2 | `1.0300e-4` |
| hole continuity residual L2 | `3.3527e-11` |
| QF bounds violations | 334 |
| drain current | `4.5484e-9 A/um` |
| electron drift / diffusion | `+274.273966788 / -274.273966784 A/um` |

最后一行显示漏极电流来自约 `274 A/um` 两个大项的高位相消，故 solver 的
`converged=true` 不能提升为曲线验收通过。

## 3. Sentaurus/Vela 状态差

仅统计 5723 个 silicon 节点：

| 状态量 | 中位绝对差 | 最大绝对差 |
|---|---:|---:|
| `psi` | `6.878 mV` | `1.1456 V` |
| `phin` | `9.77e-12 V` | `2.5870 V` |
| `phip` | `7.957 mV` | `7.7212 V` |
| `n` | `0.00156 dex` | `17.8199 dex` |
| `p` | `0.29889 dex` | `19.3210 dex` |

中位 `phin/n` 已较好对齐，但 `phip/p`、少量极端节点和接触/介质 gauge 节点仍未
关闭。最大值不得用于代表主体 silicon 区的整体偏差。

## 4. Poisson 分项闭合

探针按生产公式输出：介电边通量、电子电荷、空穴电荷、净掺杂电荷、固定界面
电荷，以及 Dirichlet 接触行替换。每行满足：

`production = dielectric + electron + hole + doping + fixed + boundary replacement`。

| 回放状态 | Poisson residual 最大值 | L2 | 分项闭合误差最大值 |
|---|---:|---:|---:|
| Sentaurus checkpoint | `7.9810e4` | `2.9303e5` | `2.91e-11` |
| Vela revision 3 checkpoint | `8.7242e-9` | `2.1557e-8` | `2.80e-11` |

闭合误差相对主贡献项约为机器精度量级。Vela 自洽状态的非接触 Poisson 行也已
闭合，因此没有证据支持在此阶段修改 Poisson 边通量、box 体积、符号或缩放。

## 5. 代表节点审计

### 5.1 p+ 重掺杂非接触节点 4670

| 量 | Sentaurus 输入 | Vela 从同一势/QF 重构 |
|---|---:|---:|
| 净掺杂 | `-2.91681e26 m^-3` | 同一输入 |
| 空穴密度 | `2.91661e26 m^-3` | `1.24907e26 m^-3` |
| 介电项 |  | `-6.0167` |
| 空穴电荷项 |  | `-5.9779e4` |
| 掺杂电荷项 |  | `+1.39596e5` |
| Poisson residual |  | `+7.98104e4` |

Sentaurus 保存的空穴密度与净掺杂近似电中性；Vela 重构值仅为其约 42.8%，直接
留下主导 Poisson 残差。这是非接触行，不能由 Dirichlet 行替换解释。

### 5.2 n+ 接触节点 4573

Sentaurus 输入电子密度为 `5.16606e26 m^-3`，Vela 重构为
`2.93213e26 m^-3`（约 56.8%）。该节点是接触行，物理 Poisson 方程会被边界残差
替换；它同时暴露载流子重构差异和接触平带/边界势差异，不能单独用于判断内部
Poisson 离散。

### 5.3 全网格载流子重构误差

| 载流子 | 中位绝对误差 | 最大绝对误差 |
|---|---:|---:|
| electron | `0.10868 dex` | `1.07138 dex` |
| hole | `0.10848 dex` | `1.21122 dex` |

同一探针对 Vela 自身 checkpoint 的输入密度与重构密度误差为严格 `0 dex`，证明
CSV round-trip 与探针读取没有制造该差异。

## 6. 门控结论与后续顺序

1. **checkpoint 保存：通过。** 文件完整，精确偏压前缀和 no-predictor 合同成立。
2. **Poisson 分项闭合：通过。** 生产残差与五项加接触行替换逐节点闭合。
3. **Sentaurus 状态在 Vela 中 reclose：失败。** 主残差由载流子重构差异驱动。
4. **G3 曲线验收：继续阻断。** continuity 残差、QF violation 和端口高位相消未关闭。
5. **下一单因素任务：** 固定该 Sentaurus checkpoint，分别审计/回放 Fermi 统计和
   OldSlotboom BGN 的 `psi/QF -> n,p` 映射参数与参考能级；优先在重掺杂 n+/p+
   节点做解析/固定状态单元测试。仅在载流子重构误差关闭后，重新评估接触边界势。

在上述映射合同关闭前，不进入 IALMob，不启用 predictor，也不修改 Poisson
离散或网格 profile。

## 7. 产物

- deck 生成器：`scripts/prepare_templates_ldmos_phase23.py`
- 分项探针：`newton_poisson_term_probe`
- 状态/算子审计器：`scripts/analyze_templates_ldmos_g3_min_transition.py`
- ignored checkpoint：
  `reference_staging/templates_ldmos_sentaurus2022/phase01_original_20260826_02/`
  `stage1_v4/phase23_t2022_contract_v3/`
  `g3_drain_prebias_vd0p023_checkpoint_state.csv`
- ignored 审计汇总：
  `reference_staging/templates_ldmos_sentaurus2022/phase01_original_20260826_02/`
  `stage1_v4/reports/g3_revision3_vd0p023_poisson_audit_20260828/summary.json`

