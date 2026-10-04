# SimpleMOS M61 路径敏感性门禁

M61/E2判定为 `not_required_m60_decisive`，因此没有创建或执行步长计划deck。

M61原本只在M60/E1结果不确定或无实质作用时执行。M60已将15个burst标志点全部消除，n23目标点观测差削减98.8585%，默认Id与Vela差降到0.0243983 dex，已对收敛残差机制给出决定性阳性判别。M63独立闭合平滑NWell成分为阈值样水平平移，M62则按写回门禁安全停止；两者都不重新触发M61条件。

机器报告：`reference_tcad/simplemos_sentaurus2022/path_sensitivity_gate/m61_path_sensitivity_gate_report.json`。
