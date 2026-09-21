# Templates/LDMOS G3：Fermi、OldSlotboom 与参考能级审计

## 结论

固定 Sentaurus `Vg=0 V, Vd=0.0231559774221138 V` checkpoint 的单因素审计已关闭
`psi/QF -> n,p` 映射合同。根因不是 Fermi 公式或 QF 符号，而是两个独立的单位/基线问题：

1. phase-2/3 生成器把 `Nref=1e17 cm^-3` 乘以 `1e6` 后写入旧式
   `unit_scaling` deck，生产解析器因而实际收到 `1e23 cm^-3`，OldSlotboom 项几乎消失；
2. 材料 `ni=1.4638914958767616e10 cm^-3` 已包含 `dEg0=-15.95 meV`，同时生产
   BGN 场又承担该能量项，造成半带隙对应的约 `7.975 mV` 参考偏差。

revision 4 将 BGN/SRH 的参考浓度按 TCAD 内部 `cm^-3` 写入，并把材料基线改为
`ni=1.07531272317754e10 cm^-3`，从 `ni` 中去掉 `dEg0`，同时保持 Sentaurus 与 Vela
物理 QF 恒等映射和 0 V 欧姆接触 QF 不变。修正后，生产组合
`Fermi + OldSlotboom + Fermi correction` 是六个单因素分支中的明确最优项。

## 1. 范围与不变量

- 固定状态：Sentaurus G3，`Vg=0 V`，`Vd=0.0231559774221138 V`；
- 网格、掺杂、接触、温度、Poisson 离散和 mobility 均固定；
- predictor 关闭；
- 未启用 IALMob、QP、热或 Okuto；
- 六个物理分支只改变 Fermi/Boltzmann、OldSlotboom 和 Fermi correction；
- 所有分支均通过同一个 Vela 生产 `newton_poisson_term_probe` 回放。

## 2. 单因素矩阵

重掺杂多数载流子统计阈值为 `|Nnet| >= 1e25 m^-3`。

| 分支 | 多数载流子中位误差 (dex) | QF 中位误差 (mV) | BGN 中位/最大误差 (eV) | 自由 Poisson L2 |
|---|---:|---:|---:|---:|
| Fermi + OldSlotboom + correction | `5.18e-4` | `0.0969` | `1.16e-5 / 2.96e-4` | `1.245e3` |
| Fermi + OldSlotboom，无 correction | `0.178` | `29.10` | `1.55e-5 / 0.190` | `1.093e5` |
| Boltzmann + OldSlotboom | `0.485` | `28.90` | `1.55e-5 / 0.190` | `1.123e6` |
| Boltzmann，无 BGN | `0.663` | `39.48` | `2.62e-3 / 0.351` | `3.906e5` |
| Fermi + correction only | `0.500` | `68.50` | `2.62e-3 / 0.161` | `3.758e5` |
| Fermi，无 BGN | `0.817` | `97.48` | `2.62e-3 / 0.351` | `4.103e5` |

三个代表节点的生产 BGN 绝对误差分别为 `0.1461 meV`、`0.1045 meV` 和
`0.04975 meV`；因此 OldSlotboom 参数、Fermi correction 和 Sentaurus 导出的 BGN
场已在亚 meV 尺度闭合。

## 3. 参考能级判定

旧材料基线下，电子和空穴的有符号 QF 映射误差在全部 silicon 节点分别近似
`-7.99 mV` 和 `+7.98 mV`，正好对应 `dEg0/2`。把 Sentaurus QF 临时平移
`+7.975/-7.975 mV` 可在固定状态上关闭该误差，但这一控制实验不能成为生产方案：
它会把 0 V 欧姆接触的物理 QF 从零移开。实际自洽试跑因此在 100 次 Newton 后失败，
hole continuity L2 为 `20.7353`。

保持 QF 不变、改用去除 `dEg0` 的材料基线后，原始 Sentaurus QF 的中位映射误差为
`0.0969 mV`；反向再施加半 `dEg0` 平移会把误差恶化到 `8.07 mV`。这构成了对
“材料基线承担一次、生产 BGN 再承担一次”的双计数诊断，也确定了 revision 4 的
参考能级合同。

## 4. 自洽资格与 checkpoint

revision 4 的 exact-mesh 结果：

| 运行 | 结果 | Newton 迭代/点 | 关键指标 |
|---|---|---:|---|
| PolySi Poisson bootstrap | 通过 | 1 点 | `converged=true` |
| 0 V coupled reclose | 通过 | 1 | electron continuity `5.89e-24` |
| 0 V independent repeat | 通过 | 0 | initial residual floor 闭合 |
| 无 predictor 精确偏压前缀 | 9/9 通过 | 终点 30 | 到达 `0.0231559774221138 V` |

终点仍不是曲线验收状态：漏极电流为 `4.65890e-9 A/um`，而 Sentaurus 固定状态
回放约为 `-9.8292e-14 A/um`；电子 continuity L2 为 `1.05508e-4`，QF bounds
violation 为 195。相较 revision 3 的 334 个 violation 有改善，但约 4.68 dex 的电流
差仍存在。因此本轮关闭的是 Fermi/BGN/参考能级映射，不把剩余差异归因于 IALMob，
下一诊断应转回有限偏压 continuity/SG 大项相消和接触邻域分支。

## 5. 交付与证据

- 固定状态审计器：`scripts/audit_templates_ldmos_fermi_bgn_mapping.py`；
- 生产诊断扩展：Poisson probe 同时输出输入/重构密度、`ni/Nc/Nv`、所需 QF 和映射误差；
- 合同：`physics_contract.json` revision 4、`materials.json` revision 3；
- ignored 证据：
  `reference_staging/.../reports/g3_fermi_bgn_mapping_audit_contract_v4_20260828/summary.json`；
- 自洽输出：
  `reference_staging/.../phase23_t2022_contract_v4/`。
