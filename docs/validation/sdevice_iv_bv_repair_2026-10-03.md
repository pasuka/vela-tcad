# SDevice IV/BV 前端修复，2026-10-03

此文保留修复阶段快照；后续原契约收口、−13.2 V 的有限场量资格、击穿区判据和
−20 V 空穴密度限制见 [WP10 收口报告](sdevice_wp10_closeout_2026-10-03.md)。

本修复基于 `5c723a2e`，响应首次 WP10 验收中的 IV 密度/接触失败及 BV 提前停扫。
源命令、网格、300 K 材料参数、扫描目标和已有验收门限保持原值。
本阶段结果表及资格边界见本文末尾；[首次 WP10 记录](sdevice_acceptance_2026-10-03.md)保留为历史对照。

## 根因与修复

1. **普通 SRH 被误激活掺杂寿命。** WP3 报告携带 Scharfetter 参数，C++ 在收到
   `srh_doping_dependence` 对象时缺省启用它。前端现在把参数声明与模型激活分开：
   普通 `SRH` 明确写 `enabled=false`；显式 `SRH(DopingDependence)` 仍启用原参数。
   这与原生 BV 日志的 `Without doping dependent lifetimes` 一致。
2. **欧姆接触的少数电子准费米电势被固定到 0。** 旧前端继承了 Newton 的
   minority-electron relaxation 默认值。p 侧加正偏压时，这使接触电子密度虚增。
   用本例 `ni=1.6551776920e10 cm^-3`、`NA=1e17 cm^-3` 与 300 K 计算，
   `2 VT ln(NA/ni)=0.8073159325 V` 正是错误电子密度超过空穴密度的交叉点，
   与旧 IV 的停点 0.8073157685 V 一致；校验随之发现 `phin=0` 不等于接触偏压。
   前端现在关闭该放松，固定两种准费米电势到欧姆接触电压。
3. **已有状态的 Coupled 延续被逐点 Gummel 初始化破坏。** 仅修物理映射仍会在高正偏
   出现小步重试；保留已接受状态的直接 Newton 对照能到达 10 V 和 −20 V。
   后续 Coupled 段改用 warm-start Newton 与 `contact_majority` 准费米参考坐标。
   增加现有局部连续性行检查、行缩放、源/接触守恒检查和 0.1 V 更新上限，防止小的
   归一化块残差掩盖错误的小电流。迭代上限仍取源 Iterations，不放宽验收。
4. **高场/雪崩继承了不同的物理默认值。** 原生日志明确打印高场迁移率和雪崩均使用
   Gradient Quasi Fermi；旧 deck 却继承电场驱动和 `carrier_density` 雪崩源。
   前端现在显式映射 DD 默认驱动力和 `current_density` 源，保留显式驱动力的优先级。

## BV 显式数值配置

PN2D BV 随附 [数值 sidecar](../../reference_tcad/pn2d_sentaurus2018/source/pn2d_bv_sdevice.vela.json)，
显式选择当前 [PN2D BV 验证参考](pn2d_bv_validation.md) 的 `element_edge_sg_gss_laux` 原子配置：
SG/GSS-Laux 电流支撑、element-vertex-box 源映射、Bernoulli 中点密度、mixed-Voronoi
节点体积及 non-obtuse 运行时资格同时应用。使用现有 local AD 源导数，保留接触元素电场回退。
该配置保留 van Overstraeten 系数、材料参数和完整扫描；C++ 全局默认和 IV 几何配置保持原值。

C10 增加显式 profile 选择及完整 `legacy_cell_reconstructed` 回退，拒绝拆分几何覆盖和
没有 GradQF avalanche 的配置。C2 的本例最大角为 `90.00000000000001°`、钝角数 0；
旧前端要求严格锐角会拒绝合法的直角网格，现按当前 non-obtuse 契约接纳 90°，仅容许
1e−8° 浮点误差，并同时要求 C2 非钝角资格和 C++ 运行时资格。

验收驱动复制并记录 sidecar SHA-256，复用完成运行时也核验它，旧无 sidecar 的 BV 结果不能冒充新配置。
这只为所选 PN2D BV 案例配置该组合，不为任意设备或任意网格建立同等资格。

## 验证方法

使用 UCRT64 Release 的真实 `vela.exe` → zipapp → `vela_example_runner.exe` 链路，
HDF5/TDR 实际启用；普通 DD Newton 使用 Eigen SparseLU，单线程。
同一原生 v2 参考、1943 节点/3680 三角形、2×0.5 µm Silicon PN 结；密度 cm^-3、
电场 V/cm、电流 A/µm、300 K。每次最终候选从平衡态开始，保留源 Plot 的 0.05 V 必经点。
原生源码/参数/网格的六项核验先通过，继而核对目标、曲线、电流守恒及偏置匹配的场量。

测试包含真实 p 型欧姆接触在 0.2 V 的两种准费米电势和少数载流子密度、普通/掺杂依赖
SRH 开关、默认/显式驱动力、延续状态、profile 原子性及回退、非钝角拒绝、sidecar 篡改与缺失。
相关 Release CTest 和三个 PN2D profile 契约测试共 **11/11** 通过。
原生产物启用的 `sdevice_e2e_pn2d_0v` 与 `sdevice_e2e_gates` **2/2** 通过；
最后复核 `sdevice_deck_emit` 与 `sdevice_acceptance_tools` **2/2** 通过。
后者覆盖 27 个 deck 测试（含两个真实 runner 测试）及 11 个验收工具测试。

## 最终结果

**原 IV 密度失败及 IV/BV 提前停扫均已修复。IV 正式验收通过；BV 完整扫描、曲线和
端电流守恒通过，但偏置雪崩场量阈值仍未获资格，WP10 整体不能标为通过。**

| 项目 | IV | BV |
| --- | --- | --- |
| 完整扫描 | 0 → 10 V，206 个接受点 | 0 → −20 V，811 个接受点 |
| 各阶段 runner / CLI 返回码 | 均为 0 | 均为 0 |
| 曲线验收窗口 | 0.2…0.3 V，3 个参考点 | −20…−0.05 V，409 个参考点 |
| 最大电流对数误差 | 0.000204181 dex | 0.00140364 dex |
| 最大电流相对误差 | 0.0470033% | 0.322678% |
| 窗口内最大端电流 KCL 相对误差 | 0.0216649% | 0.337588% |
| 电流符号 / 趋势 / 平衡场量 | 通过 | 通过 |
| 偏置场量比较 | 0.2 V：通过 | −13.2 V：诊断数值通过，阈值未获资格 |
| 验收工具状态 / 返回码 | `pass` / 0 | `fail` / 1（仅 `unqualified_threshold`） |

接受点计数包含扫描初始点；曲线门限保持 0.4 dex，KCL 门限保持 2%。IV 的
10 V 仅作为完整扫描到达性证据，电流一致性资格仍限于契约的 0.2…0.3 V 窗口。

| 偏置场量指标 | IV 0.2 V | BV −13.2 V（诊断） |
| --- | --- | --- |
| 电势 RMS | 6.24638e−6 V | 7.28655e−6 V |
| 电子密度 log10 RMS / p95 | 0.000150666 / 0.000223383 dex | 0.0000534149 / 0.000150510 dex |
| 空穴密度 log10 RMS / p95 | 0.000152043 / 0.000223387 dex | 0.0000534145 / 0.000150510 dex |
| 矢量电场 global relative L2 | 0.103805% | 0.502622% |
| 活跃密度节点数 | 各 1926 | 各 442 |

场量比较使用对应原生 TDR，活跃密度门槛为 1e5 cm^-3。BV 的上述诊断数值虽低于
现有静态非雪崩门限，契约明确禁止把该资格继承到偏置雪崩状态；工具忠实保留
`unqualified_threshold`。本次未改 `configs/sdevice_acceptance_contract.json`，未放宽门限。

## 本地产物与复现

最终报告位于 worktree 的忽略目录，生成场量和扫描输出未入 Git：

- IV：`build/iv_bv_repair_acceptance/pn2d_iv_ih52_kve/acceptance_report.json`
  （SHA-256 `d0adbfdf63ed35c9984b2f974f754c500713e478e463f8acf5d51ea0b9fe3ee2`）。
- BV：`build/iv_bv_repair_final/pn2d_bv_8hyd8uui/acceptance_report.json`
  （SHA-256 `2234fbf855ef121cfbadd5d29cfae4bd4e1f6eb96310d8b9e9236b0088951b70`）。
- 完整扫描 `stage01.csv`：IV SHA-256
  `bcd94f248e745c9c165823a9b978d99dec9d3cd21fe026431374f043d6cc75ff`；BV
  `34694e3846d40897604d9a9cf9023d8fb7b2a2178ae99cbd3d4c513288505589`。
- 原生 v2 cmd：IV `4b61ebc1c9feaf21307397e2ef732b092a90eb366082016d1a03befc36eb4b2d`；
  BV `e562e5407e60350f64b5fd6c820b4d1f2e7594fb6718cff1abca63fcbe0d6286`。
- 共同参数 `b4b3ebfdefba530f756f3855d43d7d587720689771d8badc747b61439ed42742`；
  网格 `920b55f739477e37154aa7d33af5db08bded38377d675172a6f77c7027e8d287`。
- BV sidecar `fcf449b2071751426f63a36f40b981ccf41c87953bdb53dfbef817599d26a017`。
- runner `0a4e78d01863a6e1343b695d32e6a03c23a1dd951f4a788c6d160c9ead18b6f0`。
  原生参考版本为 `T-2022.03-SP2 (0.7745337)`；报告还保留参考清单、TDR、PLT 和契约指纹。

报告内版本提交为构建时的 `5c723a2e`，运行实际包含本次未提交的前端修复；此处报告与
CSV 指纹对应最终证据，不把此前部分扫描、超时或中断的诊断运行计作通过。

在 worktree 根目录配置并构建 Release 后运行：

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --preset windows-ucrt64-release
cmake --build --preset windows-ucrt64-release
python scripts/sdevice_acceptance.py --case pn2d_iv --reference-dir reference_staging --output-dir build/iv_bv_verify --launcher build-release/vela.exe --tdr-importer build-release/sentaurus_import.exe
python scripts/sdevice_acceptance.py --case pn2d_bv --reference-dir reference_staging --output-dir build/iv_bv_verify --launcher build-release/vela.exe --tdr-importer build-release/sentaurus_import.exe
```

BV 保留全部 401 个源 Plot 偏置，且每步执行严格连续性检查，单线程完整运行耗时较长。
驱动默认壁钟上限由 1800 s 改为 7200 s；超时仍失败，求解器的源 Iterations 上限未增加。
最后一条命令因上述场量资格限制预期返回 1，应同时检查 `goals_reached`、`curve` 和
`biased_fields`，不能只用 CLI 返回码判断求解是否完成。
