# SimpleMOS M53 DirectCurrent 完整矩阵重资格

## 结论

M53 分类为 `direct_increases_global_error`。在冻结 M8/M46 TDR、物理、网格、偏压路径和 Vela 曲线后，唯一干预是在 16 个 Sentaurus 原始 deck 的全局 `Math` 中加入 `DirectCurrent`。

完整执行 8 个器件、2 个漏压、16 条 Id-Vg 曲线和 816 个精确栅压点。默认端口观测下的全局最大 Vela-Sentaurus 误差为 `0.109418680924` dex；`DirectCurrent` 下为 `2.12607506889` dex，变化 `+2.01665638796` dex。按原 M8 数值/趋势门槛，`3` 条曲线通过、`13` 条失败。

最差 DirectCurrent 点为 `n24_vd_0p05`、Vg=`0.00` V，误差 `2.12607506889` dex。原 M46 目标 n23、Vd=0.05 V、Vg=0.05 V 的默认、DirectCurrent 和 Vela 漏电流分别为 `8.586250997030e-17`、`1.422273555882e-15`、`1.104643803536e-16` A/um。

## M52 锚点回放

六个共享状态的四端电子、空穴和总电流共 `72` 行全部回放：`True`。最大绝对差为 `0.000000000000e+00` A/um，最大对称相对差为 `0.000000000000e+00`。

## NWell 配对

| 低/高 NWell 器件 | Vd (V) | 默认高减低最大误差 (dex) | Direct 高减低最大误差 (dex) | 未增加放大 |
|---|---:|---:|---:|---|
| n17/n21 | 0.05 | 0.032021 | 0.892247 | False |
| n17/n21 | 1.00 | 0.025863 | 1.266539 | False |
| n18/n22 | 0.05 | 0.022949 | 1.545543 | False |
| n18/n22 | 1.00 | 0.026044 | 1.691050 | False |
| n19/n23 | 0.05 | 0.072343 | 1.135157 | False |
| n19/n23 | 1.00 | 0.023063 | -0.222823 | True |
| n20/n24 | 0.05 | 0.013016 | 1.814475 | False |
| n20/n24 | 1.00 | 0.023341 | 1.102445 | False |

`DirectCurrent` 在 `1/8` 个 NWell-漏压配对中没有增加高 NWell 相对低 NWell 的最大误差放大。

## 边界

- M53 检验的是 Sentaurus 端口观测算法，不修改任何求解状态物理或生产默认值。
- M8/M46 官方原始 deck 默认结果仍是主要兼容性基线；`DirectCurrent` 结果不会静默替代它。
- 假设被否证也是合同允许的完成结果；未针对运行结果修改阈值。
- 未重新排查 HFS、SG、准费米、BGN、SRH、迁移率、接触模型或网格。

机器报告：`reference_tcad/simplemos_sentaurus2022/direct_current_full_matrix/m53_direct_current_full_matrix_report.json`。
