# SimpleMOS M56 Si/SiO2 界面拓扑审计

## 结论

M56 分类为 `shared_topology_continuous_potential_insulator`。八个输入TDR均未发现坐标完全相同但全局ID不同的几何重复节点；Si/SiO2 界面由两区域共享同一组全局节点ID和边来表示。求解场仍按区域分别导出，因此同一个界面全局节点在 Silicon 和 Oxide 电势文件中各有一条区域记录。

六个求解状态的所有共享界面节点均同时具有 Silicon/Oxide 电势值，最大跨区域差为 `0.000000000000e+00` V。`eDensity`、`hDensity`、`eCurrentDensity`、`hCurrentDensity` 只在 `Silicon_1` 上有支持；氧化层中的准费米绘图记录不能解释为氧化层启用了载流子输运方程。

gate 接触在八个TDR中均附着于 `Oxide_1`，source、drain、substrate 均附着于 `Silicon_1`。因此“SiO2绝缘”意味着没有氧化层电子/空穴漂移扩散输运，不意味着不能在氧化层外边界的 gate 电极施加静电势。默认 deck 不含显式 `HeteroInterface`、`Thermionic`、`Discontinuity` 或双节点开关；本算例的区域场支持由TDR材料/区域拓扑和接触边界自动确定。

## 八个输入TDR

| 器件 | 全局节点 | Si/Ox共享边 | 界面节点 | 重复坐标组 | gate区域 |
|---|---:|---:|---:|---:|---|
| n17 | 1542 | 20 | 21 | 0 | Oxide_1 |
| n18 | 1542 | 20 | 21 | 0 | Oxide_1 |
| n19 | 1480 | 20 | 21 | 0 | Oxide_1 |
| n20 | 1480 | 20 | 21 | 0 | Oxide_1 |
| n21 | 1540 | 20 | 21 | 0 | Oxide_1 |
| n22 | 1540 | 20 | 21 | 0 | Oxide_1 |
| n23 | 1482 | 20 | 21 | 0 | Oxide_1 |
| n24 | 1482 | 20 | 21 | 0 | Oxide_1 |

## 六个求解状态

| 器件 | Vg (V) | 界面节点 | 最大电势差 (V) | 载流子密度/电流支持区域 |
|---|---:|---:|---:|---|
| n19 | 0.00 | 21 | 0.000e+00 | Silicon_1 |
| n19 | 0.05 | 21 | 0.000e+00 | Silicon_1 |
| n19 | 0.10 | 21 | 0.000e+00 | Silicon_1 |
| n23 | 0.00 | 21 | 0.000e+00 | Silicon_1 |
| n23 | 0.05 | 21 | 0.000e+00 | Silicon_1 |
| n23 | 0.10 | 21 | 0.000e+00 | Silicon_1 |

## 边界

- 结论只针对冻结的 SimpleMOS TDR/deck，不外推到需要能带不连续、热发射或隧穿模型的其他异质结。
- M56 没有新增求解、修改网格、复制节点或添加界面物理。
- 准费米量在绝缘区域可作为绘图/派生记录存在，但载流子密度和电流密度的实际支持仍限定于 Silicon。

机器报告：`reference_tcad/simplemos_sentaurus2022/si_oxide_interface_topology/m56_si_oxide_interface_topology_report.json`。
