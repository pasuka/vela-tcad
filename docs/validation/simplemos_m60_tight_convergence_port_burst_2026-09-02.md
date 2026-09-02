# SimpleMOS M60 收紧收敛端口 burst 判别

## 结论

M60 分类为 `tight_convergence_suppresses_port_bursts`。32个T-2022.03-SP2 deck使用冻结的同一收敛策略包，覆盖默认/Direct两种观测、16条曲线和每种算法816个精确偏压点。物理、网格、接触和偏压路径未改变；M8/M46继续作为生产基线。

目标n23、Vd=0.05 V、Vg=0.05 V的default-minus-Direct substrate电子电流由 `4.694745539580e-17` 变为 `-5.358921274820e-19` A/um，差值削减 `98.858528%`。收紧后默认Id与Vela的差为 `0.024398` dex，观测差/Id为 `-5.131607e-03`。

完整矩阵的收紧后burst标志点为 `0`，M59原标志点为 `15`；原非burst点相对M46默认曲线的最大变化为 `2.689137e-03` dex。默认与Direct配对快照字段不变性为 `True`。

该结果只检验冻结的收敛策略包，不分离Digits与ErrRef各自贡献。次级亚fA参考资格为 `True`，且无论结果如何都不替换M8/M46。

机器报告：`reference_tcad/simplemos_sentaurus2022/tight_convergence_port_burst/m60_tight_convergence_port_burst_report.json`。
