# SimpleMOS M51 CurrentWeighting 归因

## 结论

M51 分类为 `weighted_third_observable`。唯一干预是在 M50 的全局 `Math` 中加入 `CurrentWeighting`；六状态物理、网格、偏压路径和诊断面通量保持冻结。

目标点 n23、Vd=0.05 V、Vg=0.05 V：默认 substrate `eCurrent` 为 `4.670984402844e-17` A/um，`CurrentWeighting` 后为 `-2.376113673567e-19` A/um，变化 `-4.694745539580e-17` A/um；同次原生法向面通量为 `-2.376797223030e-19` A/um。

相对于原生法向面通量，`CurrentWeighting` 消除了默认端口差异的 `99.999854401%`；剩余绝对差为 `6.835494626104e-23` A/um，相对原生量为 `2.875926713e-04`。它非常接近原生量，但未达到冻结合同的严格闭合阈值，因此保持 `weighted_third_observable`，不改写为 `weighted_matches_native_face`。

## 六状态

| 器件 | Vg (V) | 默认 substrate eCurrent | CurrentWeighting eCurrent | 原生法向面通量 |
|---|---:|---:|---:|---:|
| n23 | 0.00 | -6.851949722781e-20 | -2.376074738380e-19 | -2.376758665068e-19 |
| n23 | 0.05 | 4.670984402844e-17 | -2.376113673567e-19 | -2.376797223030e-19 |
| n23 | 0.10 | 9.946215855645e-18 | -2.376158298818e-19 | -2.376841412840e-19 |
| n19 | 0.00 | -4.674613063737e-19 | -4.308385770749e-19 | -4.309957860788e-19 |
| n19 | 0.05 | 1.655476854780e-17 | -4.308536874293e-19 | -4.310107320227e-19 |
| n19 | 0.10 | 4.208355494245e-18 | -4.308705775908e-19 | -4.310274365636e-19 |

## 状态不变性

共比较 `102` 个状态-字段文件；失败值数量为 `0`。最大绝对差为 `0.000000e+00`，最大对称相对差为 `0.000000e+00`。

## 边界

- 手册定义默认端口电流为关联掺杂阱表面通量与阱体生成率积分之和；`CurrentWeighting` 是减小数值误差的加权域积分。
- 本任务没有执行 `DirectCurrent`，没有修改 HFS、SG、准费米、BGN、SRH、迁移率、接触模型、网格或生产默认值。
- M51 只归因 Sentaurus 端口观测算法，不授权修改 Vela 或 Sentaurus 生产配置。

机器报告：`reference_tcad/simplemos_sentaurus2022/current_weighting_attribution/m51_current_weighting_attribution_report.json`。
