# SimpleMOS M66 matched-ni 完整曲线复核

M66 分类为 `matched_ni_full_curve_material_but_not_complete`。16 条 Vela matched-ni 曲线均从 M65 冻结诊断点续算，形成每条 51 点的 0--2.5 V 曲线；未启动新的 Sentaurus、平衡态或漏压态求解。组合曲线在检查点拼接处逐位复用 M65 电流行，最大身份误差为 `0.000e+00 dex`；Newton 续算起点的再次闭合变化单独记录，最大为 `1.854e-05 dex`，不被隐去或误作拼接身份误差。

在 0.55--0.90 V 平滑窗口，8 个低/高 NWell 配对的中位绝对误差增幅为 `0.026335 dex`，高 NWell 工况逐点最大绝对误差为 `0.064469 dex`（相对误差 `16.00%`）。这确认 M65 的基础 ni 差异在完整曲线上是实质成分，但没有闭合剩余平滑差异。

机器报告：`reference_tcad/simplemos_sentaurus2022/matched_ni_full_curve/m66_matched_ni_full_curve_report.json`。
