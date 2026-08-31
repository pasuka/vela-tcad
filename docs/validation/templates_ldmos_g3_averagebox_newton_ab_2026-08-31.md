# Templates/LDMOS G3 AverageBox Newton、接触 BC 与 reclose 审计（2026-08-31）

## 结论

LDMOS 专用、显式选择且默认关闭的
`templates_ldmos_external_averagebox` carrier-transport profile 已完成实现、
Jacobian 资格、一次 Newton A/B 和 `Vd=0.1 V, Vg=0.5 V` 同偏压 reclose。

直接 AverageBox 输运 couple 将 Sentaurus 固定状态上的电子 continuity residual
L2 降至基线的 `0.0173075x`，maximum 降至 `0.0433710x`，且不改变初始
Poisson residual。一次 Newton 的两个分支均产生有效状态并显著降低 residual。

自洽 reclose 同时暴露了更大的接触合同缺陷：Templates/LDMOS 的 source metal
物理上短接 p+ body pickup 与 n+ source，而 Vela 默认
`dominant_signed_contact_mean` 会把与整个接触平均符号相反的 p+ 节点重构为 n+
边界。改为本模板显式 `legacy_node_local` 后，`1.08185 V` 的最大 psi 差、约
`40.84 mV` 的 hole-QF 台阶及主要电流误差同步关闭。

最终的 `legacy_node_local + external AverageBox` 组合在 3 次 Newton 后收敛，
Drain 电流为 `6.65923e-11 A/um`，相对 Sentaurus 的
`6.46225e-11 A/um` 为 `1.03048x`、`0.0130407 dex`，单偏压阶段 3
门通过。该结果授权完成本轮诊断 profile 与 LDMOS 接触合同修订；31 点曲线仍需
由独立的 no-predictor 曲线任务验证，不由本报告宣称通过。

## 1. 数据对账

正式 31 点 G3 资格工件在 `Vg=0.5 V` 记录：

- Vela production：`1.736470037555399e-10 A/um`；
- Sentaurus：`6.46224808521937e-11 A/um`；
- 比值：`2.687099x`；
- log 误差：`0.429284 dex`。

因此同偏压 A/B 中作为 legacy reclose 输入的 `1.73647e-10 A/um` 与生产曲线
完全一致，不存在不同基线或 warm-start 多解对账问题。本轮 AverageBox-only 的
`0.292972 dex` 是候选结果，不是旧 31 点曲线的最大误差。

## 2. profile 实现合同

配置入口为：

```json
"mesh_geometry": {
  "node_volume_policy": "barycentric",
  "carrier_transport_couple_profile":
    "templates_ldmos_external_averagebox",
  "external_averagebox_couples_file": "transport_couples.csv",
  "external_averagebox_expected_edges": 16237
}
```

实现只覆盖 carrier SG couple；Poisson couple、节点控制体积和 source volume
保持基线。CSV 必须严格为 `node0,node1,couple_m`，包含全部 16,237 条边及零
couple 记录。非法节点、重复边、缺边、负值和非有限数立即拒绝。

profile 只允许 exact-mesh classical G3 诊断：IALMob、impact ionization、量子势
和 predictor 均关闭；节点体积必须为 barycentric。它不改变全局默认，也不继承
PN2D 私有 `element_edge_sg_gss_laux` 捆绑。

## 3. fixed-state 与一次 Newton

| 指标 | mesh default | external AverageBox | candidate / baseline |
| --- | ---: | ---: | ---: |
| electron residual L2 | `2.42852e-6` | `4.20316e-8` | `0.0173075` |
| electron residual maximum | `6.04064e-7` | `2.61989e-8` | `0.0433710` |
| initial psi residual difference | - | `0` | unchanged |

冻结热点 3747、10233、4492、4538 的候选/基线绝对残差比分别为
`2.04e-5`、`5.36e-5`、`2.43e-4`、`5.37e-4`，全部通过 `0.5x` 门。

在原接触合同下，一次 Newton 的 combined residual：

| 分支 | initial | post-step | ratio |
| --- | ---: | ---: | ---: |
| mesh default | `1248.8704` | `2.96583` | `0.0023748` |
| external AverageBox | `1248.8704` | `2.88174` | `0.0023075` |

在 node-local 接触合同下，两个分支也均通过一次 Newton residual-decrease 门。

## 4. `1.08185 V` 家族定位

最大差异节点为 source contact node 4663：

- 坐标：`x=-9.9787315 um, y=0.5 um`；
- 邻接：Silicon_1 与 Oxide_1.2；
- 局部净掺杂：`-2.9179583e20 cm^-3`，真实 p+；
- Sentaurus：`psi=-0.587587 V`、`p=2.9179583e26 m^-3`；
- dominant-contact Vela：`psi=0.494263 V`、`n=6.8391490e25 m^-3`；
- `delta psi=+1.0818503857 V`。

source 接触的 14 个节点从约 `-2.92e20 cm^-3` p+ 连续跨越到
`+5.17e20 cm^-3` n+。默认 `dominant_signed_contact_mean` 将前七个真实 p+
节点当成与接触平均极性相反的 outlier，并以 n 型均值替换；这对 PN2D tie-node
保护有用，但不适用于 LDMOS source short。

显式 `legacy_node_local` 后 node 4663 的 `delta psi` 为 `-17.0 uV`，多数载流子
密度与 Sentaurus 对齐；全场最大 `delta psi` 降为 `29.23 mV`。因此本家族是
已关闭的真实 Ohmic BC 缺陷，不应静默掩码。后续统计应把 Ohmic contact、
oxide-only 与 Silicon transport interior 分层报告。

## 5. 两因素 reclose 矩阵

| contact reconstruction | transport couple | Id (A/um) | Vela/Sentaurus | error (dex) |
| --- | --- | ---: | ---: | ---: |
| dominant contact mean | mesh default | `1.73647e-10` | `2.68710` | `0.429284` |
| dominant contact mean | AverageBox | `1.26869e-10` | `1.96323` | `0.292972` |
| node local | mesh default | `9.19737e-11` | `1.42325` | `0.153280` |
| node local | AverageBox | `6.65923e-11` | `1.03048` | `0.013041` |

log 空间贡献近乎正交：

- AverageBox 在两种接触合同下分别关闭 `0.136312`、`0.140239 dex`；
- node-local 接触在两种 couple 下分别关闭 `0.276004`、`0.279931 dex`。

这排除了“Sentaurus seed 与 sweep warm-start 收敛到不同解”的解释，也证明
AverageBox 的 fixed-state 改善不是只能改变 residual 排名的修饰性结果。

最终候选：

- 3 次 Newton，`block_abstol`；
- combined residual `1.69464e-8`；
- terminal KCL relative error `1.09615e-11`；
- 相对 Sentaurus state 的 max `delta psi/phin/phip`：
  `29.23/5.247/0.231 mV`。

## 6. phip/SRH 证据调和

既有 phip/SRH 审计发现 dominant-contact 结果具有约 `20--40 mV` 的空间一致
hole-QF/psi 台阶。node-local 单因素使 max phip 差从 `40.84 mV` 降为
`0.231 mV`，同时把电流误差从 `0.292972 dex` 降为 `0.013041 dex`。
因此该工件不是孤立的比较约定：它记录了 p+ source pickup 被错误多数型重构后，
hole-QF 锚点丢失向 p-body 传播的自洽反馈。

FD、OldSlotboom BGN、net/total doping 仍需作为接触公式单元测试覆盖，但在本轮
约 `0.28 dex` 主差异中不再是首要开放根因；局部 node doping 已直接重现
Sentaurus source 接触多数载流子状态。

## 7. 决策

1. LDMOS phase-2/3 deck 显式使用 `legacy_node_local`；不修改 Newton 全局默认，
   避免影响 PN2D tie-node 合同。
2. external AverageBox profile 保持模板私有、显式选择、默认关闭。
3. known-difference ledger 继续为 draft；本轮关闭的是可修复合同缺口，不登记为
   引擎固有地板。
4. 不在本任务启动 31 点曲线。下一资格任务使用显式 node-local + AverageBox、
   IALMob/predictor 关闭，运行精确 31 点并重新计算 P95、Vth、gm 和 KCL。
5. Poisson AverageBox P1/P2/P3 暂不开发。当前组合已达 `0.013 dex`，应先由 31 点
   判断剩余误差是否具有系统性，再决定是否需要全一致 Poisson/Measure 合同。

## 8. 可复现工件

- profile/门控脚本：`scripts/run_templates_ldmos_averagebox_newton_ab.py`；
- oracle 生成：`scripts/audit_templates_ldmos_averagebox_full_mesh.py`；
- 配置合同：
  `reference_tcad/templates_ldmos_sentaurus2022/contracts/diagnostics/templates_ldmos_external_averagebox_profile.json`；
- ignored 运行目录：
  `reference_staging/templates_ldmos_averagebox_newton_ab_release_20260831/`、
  `reference_staging/templates_ldmos_averagebox_newton_ab_node_local_release_20260831/`、
  `reference_staging/templates_ldmos_node_local_mesh_default_reclose_release_20260831/`。
