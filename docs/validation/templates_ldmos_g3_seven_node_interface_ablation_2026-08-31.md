# Templates/LDMOS G3 七节点三级界面消融（2026-08-31）

## 结论

针对最大 gm 误差段的七个重点节点 3721、3974、3973、4091、3727、4021、
3771，在 `Vd=0.1 V`、`Vg=1.0/1.1666667 V` 的同一 VSV 冻结状态上完成
三级消融：

- A：当前共享节点、Si-only external AverageBox 输运 couple、全局重心源项体积；
- B：保持 A 的边通量不变，仅把连续性源项体积改为 Si 侧 AverageBox Measure；
- C：在 B 上引入概念性的 Si-master/oxide-slave 电势对，并精确施加
  `psi_si=psi_oxide` 后静态消元。

七节点电子残差 L2 在低、高端点分别为 `3.15921017e-3` 和
`9.42720525e-3`，增长 `2.98403865x`。A、B、C 三层的残差逐位相同。
因此，这组冻结证据排除了“连续性源项体积”和“严格连续约束下仅改变单双节点
表示方式”是七节点误差的直接来源。

六个节点是直接 Si/介质界面节点；4021 只属于其 Si 侧一环。Sentaurus 两端点
区域导出中，六对节点可观察到的 `psi_si-psi_oxide` 均为严格零。该事实支持 C
层的精确连续约束，但不能证明 Sentaurus 内部没有保留 region-local vertex，也不能
排除两侧 Poisson box、介电通量或 Jacobian 自洽装配不同。

## 三级定义与结果

| 层级 | 单一变化 | Vg=1.0 V residual L2 | Vg=1.166667 V residual L2 | 高/低增长 | 结论 |
| --- | --- | ---: | ---: | ---: | --- |
| A | 当前共享节点 + Si-only AverageBox couple + 全局源项体积 | `3.15921e-3` | `9.42721e-3` | `2.98404x` | 基线 |
| B | 源项体积改为 Si 侧 AverageBox Measure | `3.15921e-3` | `9.42721e-3` | `2.98404x` | 相对 A 为 `1.0x` |
| C | 精确约束双节点并静态消元 | `3.15921e-3` | `9.42721e-3` | `2.98404x` | 相对 B 为 `1.0x` |

B 层虽使各节点的体积比例覆盖 `0.4342--1.0120`，但七节点上的 SRH/impact
电子源项只有约 `1e-24--4e-23`，而边通量残差为 `1e-3` 量级，故在双精度结果中
没有可见变化。它不是“体积刚好相同”，而是“被缩放的物理项在该状态下可忽略”。

C 层采用约束系统

```text
R_n(psi_si, phin, phip) = 0
R_interface = psi_si - psi_oxide = 0
```

精确消去 `psi_oxide` 后，电子连续性主行与 B 完全相同。这是代数恒等关系，不是
新拟合结果。只有允许界面电势跳跃，或令两侧 Poisson couple、控制体及介电通量
以不同方式参与自洽 Jacobian，显式双节点才可能给出 B 之外的新响应。

## 逐节点结果

| node | 分类 | A residual @ 1.0 V | A residual @ 1.166667 V | 绝对值增长 | Si Measure / global volume |
| ---: | --- | ---: | ---: | ---: | ---: |
| 3721 | 直接界面 | `-1.82888e-3` | `-5.15992e-3` | `2.82136x` | `0.43423` |
| 3974 | 直接界面 | `1.26470e-3` | `3.89619e-3` | `3.08073x` | `0.48376` |
| 3973 | 直接界面 | `1.16617e-3` | `3.64393e-3` | `3.12471x` | `0.84756` |
| 4091 | 直接界面 | `9.88618e-4` | `3.12451e-3` | `3.16048x` | `0.77638` |
| 3727 | 直接界面 | `9.86427e-4` | `3.05270e-3` | `3.09471x` | `0.50520` |
| 4021 | Si 侧一环 | `-9.39821e-4` | `-2.71774e-3` | `2.89176x` | `1.01201` |
| 3771 | 直接界面 | `-9.18011e-4` | `-2.70579e-3` | `2.94745x` | `0.90122` |

七行的符号成组相反，且每行随栅压增长约 `2.82--3.16x`，更符合保守输运行
在局部边通量失衡下形成的空间模态，而不是单节点源项或一个缺失约束造成的孤立
异常。

## 因果边界与后续靶点

本轮可以否定：

1. 将界面节点源项体积从全局重心体积换成 Si 局部 Measure 能关闭误差；
2. 在电势严格连续、carrier 只在 Si master 上装配时，单节点与受约束双节点的
   纯表示差异能关闭误差；
3. Sentaurus 已导出的两侧界面电势存在可见跳跃。

本轮不能否定：

1. Si/oxide 两侧 Poisson couple、介电通量、控制体积或界面电荷行的装配差异通过
   自洽反馈改变 Vela `psi`；
2. 七节点各条 Si 入射边上的 `psi/phin` 差、Bernoulli 因子、迁移率或 HFS 驱动力
   组合仍有差异；
3. Sentaurus 内部 region-local vertex/Jacobian 结构虽导出相同电势，但其线性化响应
   与共享节点实现不等价。

下一阶段应先对七节点的全部 Si 入射边做电子通量状态分解，再做一个保持
AverageBox、contact 与 mobility 合同不变的自洽 side-local Poisson/interface
装配试验。IALMob 和生产默认继续不变。

## 可复现工件

- 脚本：`scripts/audit_templates_ldmos_g3_seven_node_interface_ablation.py`；
- 回归：`tests/regression/test_audit_templates_ldmos_g3_seven_node_interface_ablation.py`；
- ignored 大型输出：
  `reference_staging/templates_ldmos_g3_seven_node_interface_ablation_20260831/`；
- AverageBox oracle SHA-256：
  `05ec1a43936e0ff71bea7cdfd956d80c7d692b228187ccae8625ec53ce58d49e`。

该任务未修改 C++、未增加生产自由度、未更改默认模型，也未运行或声称完成全自洽
显式 double-node 求解。
