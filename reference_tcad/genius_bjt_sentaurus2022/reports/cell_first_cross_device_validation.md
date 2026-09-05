# 跨器件先单元后节点电流恢复验证

## 结论

共完成 14 个工况、28 个载流子 A/B。MOSFET 中 17/18 组 P95 改善，且全部无实质退化；PN 中 0/10 组 P95 改善、8/10 组矢量 RMSE 改善。六节点微夹具存在点态 P95 与全局矢量 RMSE 的权衡；候选有明显跨 MOSFET 收益，但尚不具备替换全局默认值的证据，应保持默认关闭。

本轮只改变电流矢量的后处理恢复顺序，不改变已接受状态、有限体积残差、守恒边通量或端口电流。比较区域为仅与该工况半导体区域相邻的节点，活跃节点门槛为各工况 SDevice 电流峰值的 `1e-6`。

## A/B 结果

| 器件 | 工况 | 载流子 | 直接恢复 P95 [dec] | 先单元后节点 P95 [dec] | P95 变化 [dec] | 矢量 RMSE 变化 |
|---|---|---|---:|---:|---:|---:|
| pn2d_sentaurus2018_minimal6 | mirror_m12V | electron | 0.399797 | 0.588597 | +0.1888 | -0.0642742 |
| pn2d_sentaurus2018_minimal6 | mirror_m12V | hole | 0.443305 | 0.712331 | +0.269026 | -0.0688039 |
| pn2d_sentaurus2018_minimal6 | mirror_m19V | electron | 0.363169 | 0.540847 | +0.177679 | -0.0613893 |
| pn2d_sentaurus2018_minimal6 | mirror_m19V | hole | 0.392604 | 0.640759 | +0.248155 | -0.0643568 |
| pn2d_sentaurus2018_minimal6 | sketch_m12V | electron | 0.399797 | 0.588597 | +0.1888 | -0.0642742 |
| pn2d_sentaurus2018_minimal6 | sketch_m12V | hole | 0.443305 | 0.712331 | +0.269026 | -0.0688039 |
| pn2d_sentaurus2018_minimal6 | sketch_m19V | electron | 0.363169 | 0.540847 | +0.177679 | -0.0613893 |
| pn2d_sentaurus2018_minimal6 | sketch_m19V | hole | 0.392604 | 0.640759 | +0.248155 | -0.0643568 |
| pn2d_sentaurus2022_same_mesh | reverse_m20V | electron | 1.43699 | 1.46635 | +0.0293666 | +0.0304225 |
| pn2d_sentaurus2022_same_mesh | reverse_m20V | hole | 1.11722 | 1.13173 | +0.0145107 | +0.117288 |
| singledevice_sentaurus2018 | lin_0000 | electron | 0.623032 | 0.497371 | -0.125661 | -0.123445 |
| singledevice_sentaurus2018 | lin_0000 | hole | 0.847629 | 0.372835 | -0.474794 | -0.0067243 |
| singledevice_sentaurus2018 | lin_0010 | electron | 0.423732 | 0.202109 | -0.221623 | -0.0521258 |
| singledevice_sentaurus2018 | lin_0010 | hole | 0.374677 | 0.293141 | -0.0815363 | -0.00673603 |
| singledevice_sentaurus2018 | lin_0020 | electron | 0.434385 | 0.220923 | -0.213461 | -0.0379576 |
| singledevice_sentaurus2018 | lin_0020 | hole | 0.478817 | 0.30171 | -0.177107 | -0.0423764 |
| singledevice_sentaurus2018 | sat_0000 | electron | 0.608706 | 0.456911 | -0.151796 | -0.108074 |
| singledevice_sentaurus2018 | sat_0000 | hole | 1.02756 | 0.412369 | -0.615188 | -0.000218407 |
| singledevice_sentaurus2018 | sat_0010 | electron | 0.561538 | 0.261389 | -0.300149 | -0.0495148 |
| singledevice_sentaurus2018 | sat_0010 | hole | 0.474237 | 0.293362 | -0.180875 | -0.00335572 |
| singledevice_sentaurus2018 | sat_0020 | electron | 0.495574 | 0.259363 | -0.236211 | -0.0222272 |
| singledevice_sentaurus2018 | sat_0020 | hole | 0.423811 | 0.296307 | -0.127504 | -0.00326335 |
| transportmodels_sentaurus2022 | dd_deep_off | electron | 2.70082 | 2.72082 | +0.0199999 | -0.0104792 |
| transportmodels_sentaurus2022 | dd_deep_off | hole | 1.28632 | 0.42562 | -0.860698 | +0.00214967 |
| transportmodels_sentaurus2022 | dd_on | electron | 0.30195 | 0.231869 | -0.0700812 | -0.004084 |
| transportmodels_sentaurus2022 | dd_on | hole | 0.303263 | 0.206554 | -0.0967086 | -0.00610585 |
| transportmodels_sentaurus2022 | dd_threshold | electron | 0.351062 | 0.271145 | -0.0799164 | -0.0120078 |
| transportmodels_sentaurus2022 | dd_threshold | hole | 0.417262 | 0.251534 | -0.165728 | -0.00682126 |

## 证据边界

- `transportmodels_sentaurus2022` 覆盖普通漂移扩散 NMOS 的关断、阈值和导通状态。
- `singledevice_sentaurus2018` 覆盖含电子量子势状态的 NMOS，分别检查低漏压线性区和高漏压饱和区。
- `pn2d_sentaurus2018_minimal6` 是 6 节点 PN 二极管算法微夹具；`sketch/mirror` 两种三角剖分用于检查方向敏感性，不代表生产网格精度。
- `pn2d_sentaurus2022_same_mesh` 使用显式 DF-ISE 21 节点、24 三角形网格；SDevice 状态与电流矢量逐节点映射，无插值。
- 绝对误差同时包含 Vela 与 SDevice 输运模型差异；本报告主要用同一冻结状态下 direct/cell-first 的相对变化判断恢复方法的可迁移性。
