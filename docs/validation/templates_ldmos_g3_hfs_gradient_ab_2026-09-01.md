# Templates/LDMOS G3 HFS 梯度离散 A/B（WP3 T3，2026-09-01）

## 结论

在冻结提交 `01f20ac` 上，对 `Vd=0.1 V`、`Vg=1.0/1.1666667 V` 两个完整
Sentaurus SSS 状态运行生产 `newton_carrier_term_probe`。唯一变化是
`high_field_gradient_discretization` 从基线 `edge_projection` 切换为
`transport_cell_vector`。

候选显著压低固定七节点的整体残差：低、高端点 L2 分别为基线的
`0.2099999x` 和 `0.0686689x`，max 分别为 `0.2301665x` 和 `0.0571528x`。
候选自身 top-7、完整 Si/介质界面一环带和全硅节点的 L2/max 均未恶化，
所以迁移保护门通过。

但是预注册主门还要求七节点逐点均 `<=0.5x`。低端点 node 3721、4091、3771
分别为 `1.22176x`、`0.792687x`、`1.99388x`；高端点 node 3721 为
`0.539315x`。因此 **T3 总门失败**：该现成开关不能直接进入 T6 reclose，
但它对 H1 给出强筛选支持，应与 T2/T4 结果共同用于 T5 的族内定向变体。

## 冻结合同与输入

- 基线：`edge_projection`；候选：`transport_cell_vector`；
- 状态：两个端点均为完整 Sentaurus `psi/phin/phip` 的 SSS 固定状态；
- 物理：Fermi、OldSlotboom BGN、SRH/Auger、GradQF HFS、接触单元
  ElectricField fallback、Si-only external AverageBox、barycentric node volume；
- IALMob、predictor、impact ionization 均关闭；
- 未执行自洽 reclose、31 点曲线或任何生产默认修改。

输入配置与状态：

- `Vg=1.0 V`：
  `reference_staging/templates_ldmos_g3_averagebox_state_feedback_vg1_20260831/vg_1p000000/SSS/carrier_probe.json`
  与同目录 `hybrid_state.csv`；
- `Vg=1.1666667 V`：
  `reference_staging/templates_ldmos_g3_averagebox_state_feedback_vg1p166667_20260831/vg_1p166667/SSS/carrier_probe.json`
  与同目录 `hybrid_state.csv`。

两个重新运行的 `edge_projection` CSV 与上述既有基线逐字节相同：

| 端点 | SHA-256 |
| --- | --- |
| `Vg=1.0 V` | `a2082f0b724b57218e657ae75fdb612291bc1160ca411e7ad306fdbe87a78826` |
| `Vg=1.1666667 V` | `1c49fc2e1d23f3be4f2cf359f164c660a5ca587624eef0790a5d314103636c1d` |

## 主门

固定节点为 `3721, 3974, 3973, 4091, 3727, 4021, 3771`。

| 端点 | 七节点 L2 比 | 七节点 max 比 | 逐点超限 | 主门 |
| --- | ---: | ---: | --- | --- |
| `Vg=1.0 V` | `0.2099999x` | `0.2301665x` | 3721 `1.22176x`; 4091 `0.792687x`; 3771 `1.99388x` | fail |
| `Vg=1.1666667 V` | `0.0686689x` | `0.0571528x` | 3721 `0.539315x` | fail |

逐点符号在若干位置翻转。这正是整体 L2/max 改善不能替代逐点门的原因：现成
`transport_cell_vector` 大幅消除了主导七节点模态，但仍在少数节点留下重分配。

## 热点迁移保护

“完整界面带”定义为所有 Si/非 Si 直接界面节点及其 Si 侧一环，共 1474 节点；
“全硅”定义为至少属于一个 Silicon 三角形的全部 5723 个节点，包括生产探针的
边界替换行。top-7 从每个变体自身的全硅残差重算。

| 端点 | 保护集合 | L2 比 | max 比 | 结果 |
| --- | --- | ---: | ---: | --- |
| `1.0 V` | 候选自身 top-7 | `0.999999516x` | `0.999999849x` | 不恶化 |
| `1.0 V` | 完整界面带 | `0.951533529x` | `0.999999849x` | 不恶化 |
| `1.0 V` | 全硅 | `0.947405509x` | `0.999999849x` | 不恶化 |
| `1.1666667 V` | 候选自身 top-7 | `0.999998807x` | `0.999999629x` | 不恶化 |
| `1.1666667 V` | 完整界面带 | `0.842892852x` | `0.999999629x` | 不恶化 |
| `1.1666667 V` | 全硅 | `0.840494760x` | `0.999999629x` | 不恶化 |

低端点 top-7 集合完全不变：`4571, 4541, 739, 738, 732, 699, 4569`。高端点
基线 top-7 中的 `4676,4677` 被候选的 `699,4569` 替换，交集为五个节点；
候选集合上的 L2/max 仍略低于基线，因此该集合变化已报告但不构成恶化。

## 判读边界

本 A/B 排除了“现成 `transport_cell_vector` 开关完整满足冻结 T3 门”。它没有
洗脱 H1：固定七节点整体改善幅度很大且没有把全局残差推高，说明 HFS GradQF
离散约定与目标界面模态高度耦合。由于少数逐点门失败，现阶段只能把 H1 作为
T5 候选族，不能声明根因确认或启动 T6。

## 复现

```powershell
D:\msys64\ucrt64\bin\python.exe scripts\audit_templates_ldmos_g3_hfs_gradient_ab.py `
  --runner build-release\vela_example_runner.exe `
  --endpoint vg_1p000000=reference_staging\templates_ldmos_g3_averagebox_state_feedback_vg1_20260831\vg_1p000000\SSS\carrier_probe.json `
  --endpoint vg_1p166667=reference_staging\templates_ldmos_g3_averagebox_state_feedback_vg1p166667_20260831\vg_1p166667\SSS\carrier_probe.json `
  --output-dir reference_staging\templates_ldmos_g3_hfs_gradient_ab_20260901
```

- 驱动：`scripts/audit_templates_ldmos_g3_hfs_gradient_ab.py`；
- 回归：`tests/regression/test_audit_templates_ldmos_g3_hfs_gradient_ab.py`；
- ignored 原始输出与机器可读分析：
  `reference_staging/templates_ldmos_g3_hfs_gradient_ab_20260901/`，其中
  `analysis/summary.json`、`fixed_seven.csv`、`metric_ratios.csv` 为裁决入口。
