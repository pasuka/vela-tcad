# SimpleMOS M69 阶段残差定位

M69 分类为 `gate_stage_dominant`。在平衡、仅漏压、诊断栅压三个状态上复用同一 noBGN 物理路径，势垒配对增量的中位绝对份额分别为：平衡 `22.10%`、漏压 `2.70%`、栅压 `75.20%`。主导阶段为 `gate`。

Sentaurus 最终电流相对 M65 的最大重放误差 `2.666e-10 dex`。该结论定位状态差异形成阶段，不把势垒代理等同于完整漂移扩散电流闭合。机器报告：`reference_tcad/simplemos_sentaurus2022/stage_residual_localization/m69_stage_residual_localization_report.json`。
