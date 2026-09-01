# Templates/LDMOS G3 界面逐边状态分解与 Poisson A/B（2026-09-01）

## 结论

在 `Vd=0.1 V`、`Vg=1.0/1.1666667 V`、IALMob 与 predictor 关闭的冻结合同下，
已完成七个栅沟道界面热点节点全部 40 条入射边的生产 SG 核分解，以及默认关闭的
region-local AverageBox Poisson 自洽 A/B。

逐边分解以 `1.34e-14` 以内的相对误差精确重构生产电子通量。在同偏压
`VSV -> VVV` 闭合中，保持 VSV 的准费米驱动力会分别保留 `92.29%`、`80.87%`
的初始七节点残差；保持旧 couple 的残差仅为初始值的约 `1e-11`，单独保持旧
mobility、广义 Einstein、Bernoulli 或右端密度均只留下不超过 `0.91%`。这说明
本轮残差关闭由 `phin` 状态更新承载，并非 couple 或 HFS 幅值随偏压漂移。

Poisson 候选将 30,022 条网格边改为按材料分区加权的 AverageBox couple，其中
760 条为异质材料邻接边。两个端点均在 4 次 Newton 内通过原冻结门，但
`gm/Sentaurus` 从 `0.7569242` 变为 `0.7559671`，略微变差。因而 P1 级
region-local Poisson edge couple 不是剩余 `24.31%` 最大 gm 差的主因，不进入
生产默认；known-difference ledger 继续保持 draft。

## 逐边 SG 核合同

生产 Fermi–Dirac SG 通量按以下乘积导出并逐边记录：

```text
Jn = (couple/length) * mobility * Vt * generalized_Einstein
     * right_density * Bernoulli(argument) * expm1(qf_argument)
```

诊断新增输出只暴露生产公式的既有中间量，不改变 residual 或 Jacobian：

- HFS 驱动力与 edge mobility；
- 两端 Fermi `eta`；
- drift potential、广义 Einstein 因子；
- Bernoulli argument 与 `B(+x)/B(-x)`；
- quasi-Fermi argument。

40 条唯一入射边中有 13 条 oxide/interface 零载流子边，保持零通量并单独计数，
不把其非载流子 Fermi 输入混入统计。

## 偏压增长分解

VSV 七节点 residual L2 从 `3.15921017e-3` 增至 `9.42720525e-3`，增长
`2.98403865x`。将高偏压的单个因子替换为低偏压值所得 L2/高偏压基线为：

| 保持低偏压因子 | L2 / 高偏压基线 |
| --- | ---: |
| geometry | `1.000000` |
| mobility | `1.002964` |
| generalized Einstein | `0.990863` |
| Bernoulli | `1.168116` |
| right density | `0.591374` |
| QF drive | `0.498119` |

该表是非线性反事实，不能把各项作线性贡献相加。它能排除 geometry、mobility
和 Einstein 的偏压增长是 `2.984x` 的直接来源；density 与 QF drive 随状态共同
变化，Bernoulli 单独冻结会改变大项相消并使残差反而增大。

七行高偏压 cancellation condition 约为 `12--30`。残差是相邻大通量的差，
不是单条异常边：例如 node 3721 的主导边贡献为 `-0.02669`、`+0.01586`、
`+0.01481` 和 `-0.00914`。

## 同偏压 VSV 到 VVV 分解

| Vg (V) | VSV residual L2 | VVV residual L2 | 保持旧 QF drive / VSV | 其余单因子最大 / VSV |
| ---: | ---: | ---: | ---: | ---: |
| 1.000000 | `3.15921e-3` | `4.64302e-14` | `0.922875` | `0.003452` |
| 1.166667 | `9.42721e-3` | `6.25827e-14` | `0.808728` | `0.009084` |

该结果只证明 Vela 固定点关闭路径主要由 `phin` edge drive 携带；它不等同于
“Sentaurus 的 SG 公式已知不同”。真正的跨引擎因果仍需直接比较广义 Einstein/SG
边平均语义或取得 equation-balance 证据。

## Region-local Poisson P1 A/B

候选 profile 对每条边构造等效 couple：

```text
couple_eff = sum_cell(eps_r(cell) * AverageBox_local_couple(cell))
             / arithmetic_mean_edge_cells(eps_r)
```

Vela 现有 Poisson 装配仍以 `edgeAvg(eps_r) * couple_eff / length` 使用该值，因此
乘积严格等价于分区求和。carrier transport 独立保持已资格化 external
AverageBox couple；节点电荷体积仍为 barycentric。本轮是 P1 couples-only 诊断，
不是 P2/P3 全一致 Poisson 合同。

| 指标 | Vg=1.0 V | Vg=1.166667 V |
| --- | ---: | ---: |
| baseline Id (A/um) | `1.021797549e-6` | `2.213015335e-6` |
| candidate Id (A/um) | `1.020859307e-6` | `2.210570765e-6` |
| candidate / baseline | `0.999082` | `0.998895` |
| Newton iterations | `4` | `4` |
| max `|delta psi|` (V) | `0.048318` | `0.056976` |
| max `|delta phin|` (V) | `4.61e-5` | `8.64e-5` |

| gm | A/(um V) | Sentaurus 比值 |
| --- | ---: | ---: |
| Sentaurus | `9.442565684e-6` | `1.000000` |
| baseline Vela | `7.147306711e-6` | `0.756924` |
| candidate Vela | `7.138268744e-6` | `0.755967` |

候选 region-local interface 一环的 analytic-versus-symmetric-FD JVP 最大相对误差为
`9.08e-8` 和 `4.76e-8`。因此负结果不是候选 Jacobian 与其 primal residual
不一致造成的。

## 决策

1. 七节点 continuity source Measure、严格约束 double-node 表示、carrier couple、
   HFS/mobility 幅值与 P1 region-local Poisson couple 均已完成单因素审计。
2. `templates_ldmos_region_averagebox` 只保留为显式、默认关闭的诊断 profile；不改变
   全局 C++ 默认，也不继承 PN2D 私有原子捆绑。
3. 下一最小靶点是 Fermi–Dirac generalized-Einstein/SG 的边平均语义：在相同
   `psi/phin/n` 上比较 chemical-potential derivative 的端点、算术、对数或积分平均，
   先做七节点固定状态算子 A/B，再决定是否运行自洽端点。
4. 在该比较完成前，不启用 IALMob，不调整阈值、迁移率或 flatband，不批准
   known-difference ledger，也不重跑 31 点曲线。

## 可复现工件

- 逐边脚本：`scripts/audit_templates_ldmos_g3_interface_edge_state.py`；
- Poisson A/B：`scripts/audit_templates_ldmos_g3_side_local_poisson_ab.py`；
- 默认关闭 profile 实现：`src/simulation/ConfigParsing.cpp`、
  `src/mesh/DeviceMesh.cpp`、`src/tools/vela_example_runner.cpp` 与
  `src/simulation/DCSweep.cpp`；
- C++ 隔离测试：`tests/test_box_geometry.cpp`、`tests/test_newton_solver.cpp`；
- 单元测试：
  `tests/regression/test_audit_templates_ldmos_g3_interface_edge_state.py`、
  `tests/regression/test_audit_templates_ldmos_g3_side_local_poisson_ab.py`；
- ignored 运行证据：
  `reference_staging/templates_ldmos_g3_interface_edge_state_20260901/`、
  `reference_staging/templates_ldmos_g3_side_local_poisson_ab_20260901/`。
