# Templates/LDMOS G4、HFS coefficient 与 phip/SRH 审计

日期：2026-08-29

状态：冻结顺序中的三项审计已完成；接触 HFS 支撑语义已定位，现有
coefficient 开关被排除，phip/SRH 被独立分类；阶段 3 P95 门仍失败，
known-difference ledger 保持 draft。

## 范围与纪律

本轮继续使用 T-2022.03-SP2 Templates/LDMOS exact 10,241-node topology、
`Vd=0.1 V`、Fermi、OldSlotboom、SRH/Auger，并保持 IALMob 与 predictor
关闭。G4 只相对 G3 删除 `Mobility(HighFieldSaturation)`。Sentaurus 新任务
全部在隔离目录运行，未修改封存 oracle。

没有启用 Masetti、标定 bulk mobility、施加门压平移，也没有把 PN2D 的
`element_edge_sg_gss_laux` 私有原子捆绑移植到 LDMOS。

## G4 三点固定状态回放

Sentaurus G4 新保存七个 exact-bias TDR；固定状态回放选取
`Vg=1/6, 1/2, 5/6 V`。基线使用常数 mobility、barycentric node volume 与
正 barycentric 负-cotangent fallback。

| Vg (V) | Sentaurus terminal (A/um) | Vela SG / Sentaurus | 绝对误差 (dex) |
| ---: | ---: | ---: | ---: |
| 0.166667 | 7.95617e-14 | -1.76170 | 0.24593 |
| 0.500000 | 7.81197e-11 | 1.06229 | 0.02624 |
| 0.833333 | 1.19060e-7 | 1.05924 | 0.02499 |

最低点处电流符号受低电流边相消影响，只作诊断；中高两点已从 G3 的
`11.08/11.06x` 坍缩到 `1.062/1.059x`。三点 log-ratio 中位离 G4 自洽曲线
平台比离 G3 固定状态回放近 `0.758 dex`，支持 **assembly-real** 分支：
Sentaurus HFS 接触语义确实进入了产生状态/端口电流的离散装配，不是纯
plot-time TDR 重评估。

## LDMOS 专用 coefficient A/B

在完全相同的三个 G4 状态上运行四个显式组合：

1. barycentric node volume + 正 barycentric 负-cotangent fallback；
2. barycentric node volume + truncated Voronoi（负局部 cotangent 置零）；
3. mixed-Voronoi node volume + 正 barycentric fallback；
4. mixed-Voronoi node volume + truncated Voronoi。

逐边输出包含 `couple`、`couple/length`、平均掺杂、两端密度、mobility 和
particle line flux。30,022 条输运边中，truncation 改变 1,042 条边，最大
couple 变化 `7.32761e-8 m`；但端口回放变化不超过约 `0.17%`，中高偏压约
`0.07%`。mixed-Voronoi node volume 对固定状态 SG 的逐边 flux 严格不变，
符合“节点体积不进入固定 edge flux”的正交控制预期。

因此当前已有的 barycentric/truncated 开关不能解释 G4 自洽曲线中段
`~1.9x` 平台。该平台仍需在实际 box coefficient、状态反馈和边/单元装配层
继续审计，不能登记为已知引擎地板。

## G3 element HFS 驱动力

Sentaurus G3 重新导出了三个代表状态的 element `eMobility`、
`eGradQuasiFermi` 与 `ElectricField`。T-2022.03 不支持
`eEparallel/Element`，首次隔离运行在写状态时明确报
`undefined Element-Scalar eEparallel`；修正任务保留 nodal eEparallel，并
成功生成 7 个状态。

使用合同冻结的 Caughey--Thomas 参数（`mu0=1417 cm2/(V s)`、
`vsat=1.07e7 cm/s`、`beta=1.109`）拟合：

| 区域 | 样本 | 保存 element mobility | GradQF 预测中位误差 | ElectricField 预测中位误差 |
| --- | ---: | ---: | ---: | ---: |
| 全硅（代表点） | 10,515 | 中位约 1417 | 约 0 dex | 0.046--0.116 dex |
| drain 节点一环 | 48 | 20.78--1015.71 | 1.295 dex | 约 1e-16 dex |
| drain 边界单元 | 24 | 20.78--1015.71 | 1.299 dex | 约 1e-16 dex |

三个偏压的接触结果稳定。正确分类为：硅体内部使用 GradQF，drain 接触边界
在零/未定义 GradQF 支撑下使用 ElectricField 型回退。不能据此把全局 HFS
驱动力改成 ElectricField。deck 与日志仍明确声明无 DopingDependence，G4
平衡态全域 mobility 为 1417；Masetti 归因继续被排除。

## phip 参考与 SRH 敏感性

对七个公共偏压的 5,723 个硅节点比较：

- `Vela - Sentaurus phip` 的跨偏压中位为 `+20.163 mV`，最大中位
  `20.245 mV`；
- 全硅去中心残差 P95 约 `15.6--16.3 mV`，但这些尾部集中在近零空穴区域；
- 使用 `max(pVela,pSentaurus) >= 1e16 m-3` 的非验收诊断掩码后，1,208 个
  hole-populated 节点的去中心 P95 仅 `0.217--0.269 mV`；
- hole-populated 节点的原始 p 中位误差为 `0.00004--0.00084 dex`，而把
  phip 统一减去约 20 mV 并按 Vela 映射重构 p 后，误差恶化为
  `-0.340 至 -0.327 dex`；
- 全硅原始 p 的 `+0.220--0.233 dex` 中位差因此属于近零密度节点统计污染，
  不能用于物理标定。

这证明约 +20 mV 是 hole-populated 区域中的 QF 参考/推导约定差，而不是应
加入 flatband、BGN 或 density mapping 的物理位移。P1 phip 梯度比较对统一
offset 不敏感，单独保存在审计 JSON 中。

SRH 仍对该参考约定敏感。原生 TDR `srhRecombination` 按 element-lumped 面积
积分，并与 Vela frozen-state `srh_balance` 对比：

| Vg (V) | Sentaurus SRH net | Vela 原状态 | 仅平移 phip | 平移 phip + 映射一致 p |
| ---: | ---: | ---: | ---: | ---: |
| 0.0 | -2.73011e-16 | -1.83036e-16 | -2.42819e-16 | -2.77529e-16 |
| 0.5 | -2.79981e-16 | -1.86392e-16 | -2.48855e-16 | -2.84702e-16 |
| 1.0 | -2.71234e-16 | -1.77085e-16 | -2.38922e-16 | -2.73066e-16 |

单位均为 `A/um`。映射一致的固定状态平移使 SRH 积分接近 Sentaurus，但使有效
空穴密度变差，说明需要修正的是 generalized SRH 内 QF-split 与能带参考合同，
不能修改载流子密度映射。该固定状态控制对 drain Id 的最大相对影响约
`1.84e-5`，对当前 Id-Vg 次要，但会影响后续体区复合/产生校核。

## 结论与下一开发项

1. 优先实现 LDMOS 接触边界 HFS 回退：内部保持 GradQF，仅在接触边界的
   未定义/零 GradQF 支撑上使用合同化 ElectricField 回退；增加固定状态
   G3/G4 单元测试与三点回放门。
2. 通过 HFS 回放后重跑无 predictor 的 31 点 G3；再审计 G4 自洽 `~1.9x`
   平台的状态反馈和 box coefficient。现有 cotangent fallback 与 node-volume
   开关已排除，不应重复。
3. 独立开发 generalized SRH QF-reference 合同测试；先以固定状态 SRH
   integral 闭合，不把约 20 mV 写成 flatband 或 density 修正。
4. 完成上述资格前不启动 IALMob，不启用 predictor，不批准 known-difference
   ledger。

## 工件与验证

可重复脚本：

- `scripts/audit_templates_ldmos_g4_fixed_state_coefficients.py`；
- `scripts/audit_templates_ldmos_g3_element_hfs.py`；
- `scripts/audit_templates_ldmos_phip_srh.py`。

大型生成物位于 ignored `reference_staging/`：

- `templates_ldmos_g4_fixed_state_coefficient_audit_20260829/`；
- `templates_ldmos_g3_element_hfs_audit_20260829/`；
- `templates_ldmos_phip_srh_audit_20260829/`。

已执行 4 个新增 Python regression tests，全部通过。Sentaurus G4 状态任务
通过；G3 element 任务在删除不受支持的 element eEparallel Plot 项后通过，
失败任务作为导出能力证据保留。
