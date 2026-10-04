# SimpleMOS M59 端口观测 burst 与平滑 NWell 增幅分解

## 结论

M59 分类为 `port_observable_burst_association_plus_smooth_nwell_growth`。它只读取冻结的 M46、M54 和 M58 账本，没有执行新的 Sentaurus 或 Vela 求解，也没有改写历史账本。

在 `|substrate(default-Direct)eCurrent / Sentaurus Id|` 阈值从0.005扫描到0.03时，八个匹配NWell/漏压行的burst-free每曲线最大误差增幅始终位于 `0.0211` 到 `0.0260` dex。这证明M58使用的每曲线最大值混合了平滑NWell增幅与低电流端口观测burst，但平滑成分的BGN、静电或迁移率归因仍留给M63。

在主阈值0.02下，低漏压标志点的零截距斜率为 `0.370`、Pearson为 `0.996`；高漏压斜率为 `-0.107`。因此0.36量级的正关联只用于Vd=0.05 V，不能外推到全部15个混合漏压标志点。

`default-Direct` 是M54直接记录的端口观测差，不是独立导出的节点连续性残差。把它解释为p阱内部残差和是等待M60收敛控制实验检验的box恒等式假说。M59不声明Sentaurus求解器缺陷，也不授权生产参考修改。

M58中仅n21/n23含单节点浮置n型分量；它们的默认四端KCL中位幅值相对所有无浮置器件的最大中位值至少高 `334.6` 倍。这支持浮置阱调制默认KCL观测，但不是内部算法实现的独立证明。

## 下一门槛

M60只允许冻结的Sentaurus收敛控制干预，检验收紧收敛后端口观测burst是否下降，同时非burst曲线点保持稳定。M60之前不得把端口观测差称为已证实的Newton残差。

机器报告：`reference_tcad/simplemos_sentaurus2022/port_observable_decomposition/m59_port_observable_decomposition_report.json`。
