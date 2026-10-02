# SimpleMOS M63 平滑 NWell 误差归因

## 结论

M63 分类为 `threshold_like_horizontal_shift_dominant`。分析只读取M60收紧收敛后的16条默认曲线、M46/M8的Vela点值以及既有M13/M40账本，没有新增求解或修改生产参数。

在冻结的Vg=0.55--0.90 V窗口内，16条曲线的水平平移模型误差能量解释率中位数为 `0.997962`，8个NWell配对增幅的解释率中位数为 `0.973200`。拟合的Sentaurus-minus-Vela等效栅压平移范围为 `2.709` 至 `5.673` mV。

M13在Vg=0.8 V的n17/n21固定状态迁移率链有 `0/2` 个Vd工况与端口误差增长同号，因此现有固定状态证据不支持“迁移率误差随NWell增大”作为主导解释。

该结论只闭合到曲线级的阈值样/电静力学平移；它不能单独证明BGN根因。M40仅作为深关断BGN/SRH交互背景保留，BGN、载流子统计与自洽电势对这一平移的各自贡献仍需独立干预才能区分。

机器报告：`reference_tcad/simplemos_sentaurus2022/smooth_nwell_attribution/m63_smooth_nwell_attribution_report.json`。
