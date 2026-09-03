# SimpleMOS M64 BGN 平滑 NWell 干预

## 结论

M64 分类为 `bgn_independent_threshold_shift_dominant`。在 M60 收紧收敛策略、默认端口观测和完整偏压路径不变的条件下，Sentaurus 与 Vela 同时从 OldSlotboom BGN 切换到显式 no-BGN；其余物理、网格、接触和生产默认值均未改变。

8 个 NWell 配对的误差增幅中位闭合率为 `0.000000`，其中 `0/8` 个配对达到 80% 闭合。BGN-on 配对增幅中位数为 `0.038263` dex，BGN-off 后为 `0.111220` dex，差分中的差分为 `-0.073282` dex。

8/8 配对中，启用 BGN 的响应都与 no-BGN 的 NWell 增幅反向；它对 no-BGN 增幅的中位抑制比例为 `66.708390%`。因此，BGN 不是 M63 平滑增幅的来源，而是在抵消一个更大的 BGN-independent 阈值样失配。

BGN-on 的拟合水平平移范围为 `2.709`--`5.673` mV；BGN-off 后为 `13.171`--`19.912` mV。全部平滑窗口曲线保持单调，差分恒等式最大数值残差为 `4.441e-16` dex。

这一干预归因的是两个求解器对“启用 BGN 后的自洽响应”之差，不等于 OldSlotboom delta-Eg 公式错误；M39 的公式级闭环与 M8/M46/M60 参考资格均保持不变。

机器报告：`reference_tcad/simplemos_sentaurus2022/bgn_smooth_nwell_intervention/m64_bgn_smooth_nwell_intervention_report.json`。
