# Templates/LDMOS G3 Fermi–Dirac SG 边平均固定状态 A/B（2026-09-01）

## 结论

在 `Vd=0.1 V`、`Vg=1.0/1.1666667 V`、IALMob 与 predictor 关闭的冻结
VSV 合同下，已对七个栅沟道界面热点节点的 24 条有效输运边比较生产 secant、
经典 Einstein、eta 中点、端点算术/几何/调和、密度加权、对数密度中点、eta
积分近似以及端点上下界共 10 种替代语义。

没有候选通过预注册门槛。最优的 `classical` 仅将两个端点的 residual L2 降至
生产值的 `0.9831x` 和 `0.9728x`，最大残差降至 `0.9859x` 和 `0.9787x`，但
最差节点反而增至 `1.0209x` 和 `1.0340x`。距离“两个端点的 L2、maximum 和
每个选定节点均不超过生产值 `0.5x`”的门槛很远。

因此，已枚举的 generalized-Einstein/SG 边平均选项不是七节点残差的充分修正；
按停止规则不实现 C++ profile、不进行自洽端点 reclose，也不启动 31 点曲线。
该负结果不等同于证明 Sentaurus 与 Vela 的完整 SG 离散相同；下一步需取得
Sentaurus equation-balance 或直接 element/edge current 证据，继续审计 QF drive
的离散定义。known-difference ledger 维持 draft。

## 冻结合同与门槛

- 状态：`VSV = (Vela psi, Sentaurus phin, Vela phip)`；
- carrier couple：已资格化的 external AverageBox；
- 温度：`300 K`；IALMob、predictor 均关闭；
- 节点：`3721, 3974, 3973, 4091, 3727, 4021, 3771`；
- 每个偏压含 40 条唯一入射边，其中 24 条为有效电子输运边；
- 通过条件：两个端点的 L2、maximum 与七个逐点残差都必须 `<=0.5x` 生产值。

本轮是只读固定状态回放，不修改生产 residual、Jacobian 或默认配置。生产 secant
因子的离线重构最大相对误差为 `2.11e-14`，生产通量重构最大相对误差为
`1.30e-16`，排除了诊断公式自身的符号或缩放错误。

## 候选定义

生产广义 Einstein 因子为 Bessemoulin–Chatard 型 secant：

```text
g_secant = (eta1 - eta0) / log(n1 / n0)
```

替代项只改变同一冻结边上的 `g`，同时一致地重算 SG prefactor、Bernoulli
argument 与 quasi-Fermi argument。候选包括：

- `classical = 1`；
- `eta_midpoint` 与 `log_density_midpoint`；
- 两端局部 `F_{1/2}/F'_{1/2}` 的算术、几何、调和及密度加权平均；
- eta 区间的 Simpson 积分近似；
- 两端局部因子的 minimum/maximum bounds。

## 结果

生产七节点 residual 为：

| Vg (V) | L2 | maximum abs |
| ---: | ---: | ---: |
| 1.000000 | `3.15921017e-3` | `1.82887663e-3` |
| 1.166667 | `9.42720525e-3` | `5.15991936e-3` |

候选相对生产结果如下，排序依据为两个端点中较差的 L2 比值：

| 候选 | L2 ratio @ 1.0/1.1667 V | worst maximum ratio | worst node ratio | 通过 |
| --- | ---: | ---: | ---: | --- |
| classical | `0.9831 / 0.9728` | `0.9859` | `1.0340` | 否 |
| endpoint minimum | `0.9919 / 0.9871` | `0.9914` | `1.0049` | 否 |
| log-density midpoint | `0.9990 / 0.9987` | `0.9988` | `1.0000` | 否 |
| eta midpoint | `0.9991 / 0.9988` | `0.9988` | `1.0000` | 否 |
| eta integral (Simpson) | `1.0000 / 1.0000` | `1.0000` | `1.0001` | 否 |
| endpoint harmonic | `1.0019 / 1.0024` | `1.0035` | `1.0063` | 否 |
| endpoint geometric | `1.0019 / 1.0025` | `1.0036` | `1.0064` | 否 |
| endpoint arithmetic | `1.0019 / 1.0026` | `1.0037` | `1.0066` | 否 |
| density weighted | `1.0063 / 1.0099` | `1.0123` | `1.0228` | 否 |
| endpoint maximum | `1.0119 / 1.0180` | `1.0210` | `1.0387` | 否 |

`classical/production` 的边因子范围在低端点为 `0.9768--0.9974`，高端点为
`0.9650--0.9957`。也就是说，该界面簇的 Fermi 修正幅度只有约 `0.3%--3.6%`；
不同合理平均的影响更小，无法关闭当前残差。

## 决策与范围限制

1. 排除“上述 10 种 generalized-Einstein 因子平均中的任一种可单独关闭七节点
   VSV 残差”的假设；不把结果泛化为所有 Sentaurus SG 语义均已排除。
2. 不新增 generalized-Einstein C++ 配置项，不触发端点 reclose，不重跑 31 点。
3. 保持 external AverageBox、`legacy_node_local` 和已资格化接触 HFS 合同不变；
   IALMob 与 predictor 继续关闭。
4. 下一最小证据是 Sentaurus 两端点的 electron equation-balance 或直接
   element/edge current 导出，用同一七节点簇对齐逐边 QF drive、密度支撑与残差。
5. 在直接跨引擎算子证据取得前，差异账本不批准为固有引擎地板。

## 可复现工件

- 审计脚本：`scripts/audit_templates_ldmos_g3_fermi_edge_average_ab.py`；
- 单元测试：`tests/regression/test_audit_templates_ldmos_g3_fermi_edge_average_ab.py`；
- ignored 运行证据：
  `reference_staging/templates_ldmos_g3_fermi_edge_average_ab_20260901/`；
- 输入逐边状态：
  `reference_staging/templates_ldmos_g3_interface_edge_state_20260901/`。
