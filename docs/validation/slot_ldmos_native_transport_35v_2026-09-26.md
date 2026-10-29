# SLOT 完整输运几何、Jacobian 与 35 V 重收敛（2026-09-26）

## 结果及范围

完整原生 Si 输运几何的受控配置已完成装配复制、Jacobian 检查及 **35 V
IIC/自洽雪崩独立重收敛**。两点均从此前本地合格态启动，各执行 4 次完整
Newton 更新；电子/空穴各 12634 个有效行和全局闭合均通过，未用原生状态
初始化，未使用额外候选修正或 Poisson 预热。

| 35 V 指标 | 原外心体积配置 | 本轮完整原生输运几何 |
| --- | ---: | ---: |
| IIC 电流差 | −0.000489% | −0.000713% |
| SC 电流差 | +0.020205% | **−0.001768%** |
| IIC 总雪崩源差 | +0.004561% | **−0.001627%** |
| SC 总雪崩源差 | +0.063135% | **−0.002976%** |
| IIC 总源节点 L1 相对差 | 0.028277% | **0.001903%** |
| SC 总源节点 L1 相对差 | 0.093213% | **0.003173%** |
| Si 势差 RMS | 0.039319 mV | 0.039319 mV |

IIC 电流的极小差异略增，不能声称每个指标都单调改善；局部源分布和 SC
指标有明显改善。IIC/SC 电流分别为 `5.28308262732e-15` 和
`7.85667248749e-15 A/µm`。源积分分别为 `1.17501955433e-15` 和
`2.57365640168e-15 A/µm`，源比较采用共同有效 Si 12306 节点。

本轮资格仅为新输运配置的两个 35 V 固定偏压点。先前 35→38 V 连续路径
属于上一配置，不能直接沿用；新配置的连续高压、外电阻及 BVDS 仍待验证。

## 显式配置与实现

本轮只切换完整 Si 输运几何：

```json
"mesh_geometry": {
  "edge_coupling_policy": "assembled_delaunay_fallback",
  "node_volume_policy": "mixed_voronoi",
  "require_non_obtuse": false,
  "carrier_transport_couple_profile": "slot_ldmos_external_averagebox",
  "external_averagebox_couples_file": "<absolute path>/native_transport.csv",
  "external_averagebox_expected_edges": 36681
}
```

同时把 `solver.mobility.ialmob.geometry_file` 指向匹配的
`local_control_r2/native_geometry.json`，包含原生元素边系数及顶点迁移率平均
测度。36681 条 CSV 边恰好覆盖全部 Si 元素边；全网格共有 41155 条边，
其他边保留既有耦合。只换 CSV 或只换迁移率测度均不等价于本轮配置。

保持外心电荷/复合体积、`per_cell` 介电耦合、IALMob 参数和原测得雪崩几何
不变。迁移率顶点平均、输运边权重、雪崩迁移率选点、雪崩源分配是不同配置
对象，不作隐式联动或全局默认替换。

此前外部输运 profile 的入口统一拒绝雪崩。本轮将探针与 DCSweep 的检查
合并为 `validateCarrierTransportCouplePhysics`，仅为 SLOT profile 放行：
IALMob、`van_overstraeten`、`ialmob_element_density_grad_qf`、`local_ad`。
源核原有细项约束继续检查；量子势和其他雪崩源仍拒绝。Templates profile
原有限制及 `mesh_default` 路径保持。未修改源公式或导数算法。

## 装配与 Jacobian 验证

1. 新二进制在原配置、原生 35 V IIC 冻结态上的 SG 边探针输出与上一轮
   **逐字节相同**。
2. 开启完整原生输运几何后，C++ 实际装配与上一轮 Python/元素核离线控制
   在 24 个目标节点一致；误差按原边通量绝对和归一化，最大 `3.12e-15`。
   原生状态只用于这项冻结诊断，不授予本地收敛资格。
3. `Observed IALMob avalanche assembles every live potential column in both
   unit systems` 通过：包含任意元素输运权重、live IALMob、源 Jacobian、
   全部 18 个势变量列、Si/oxide 界面及 SI/TCAD 两种尺度。
4. 在实际原生 35 V SC 状态上，检查节点 **10484、10494、8237** 的
   psi/phin/phip 三类扰动，以 `1e-6/5e-7 V` 两个步长构成 18 个方向。
   每次比较完整器件装配的 Jv 与残差中心差分，不冻结迁移率或雪崩源。

实际器件方向检查的电子、空穴方程块最大归一化误差分别约
`1.35e-7 / 1.36e-7`，低于 `3e-5`；psi 扰动对应的 Poisson 块最大
`1.03e-9`。未采用工具输出默认的 `max(1,norm)` 分母掩盖小电流误差，
而是按解析/差分块范数的较大者归一化。

少数 QF→Poisson 交叉块很弱，有限差分受抵消限制、相对误差可达 1；
**不声称这些近零项逐项达到同样相对精度**。相应完整列的物理导数由上述
可分辨解析夹具覆盖；实际器件检查也不是对全部 41643 个列逐列差分。
方向探针关闭算法性载流子对角地板以检查物理 Jacobian；重收敛保留原求解器
地板配置。不能把地板当成物理残差导数。

## 35 V 完整器件资格

| 模式 | 完整 Newton 更新 | 最大有效载流子行比 | 行门槛 | 全局闭合 |
| --- | ---: | ---: | ---: | --- |
| IIC | 4 | 2.00455e-5 | 0.001 | 通过 |
| SC | 4 | 1.63756e-5 | 0.001 | 通过 |

本轮 Newton 上限及预算资格下限为 40。`abstol=1e-9`、`reltol=1e-8`、
stall `1e-8`、载流子行门槛 `0.001`、全局闭合门槛 `0.01` 均未放宽。
IIC 求解方程的电子/空穴雪崩反馈项确认为零，源比较由独立 observable 探针
计算；SC 则使用完整自洽源及其导数。

近表面 IIC 空穴源对比：

| 节点 | 修改前源差 | 本轮源差 |
| --- | ---: | ---: |
| 10487 | −3.317717% | −0.001430% |
| 10490 | −3.473858% | −0.001483% |
| 10494 | +11.456612% | −0.002866% |
| 6687 | +2.267707% | −0.001737% |

10494 的物理空穴 QF 差从 −18.1286 mV 降到 **−0.026489 mV**，密度比值
从 0.496477 恢复为 **1.000003989**。这与上一轮冻结输运残差的定位一致，
支持近表面输运权重/迁移率平均测度是该局部状态和源差的主要原因。

SC 电子/空穴源积分差分别为 −0.002851%/−0.003288%，节点 L1 相对差
为 0.003106%/0.003352%，改善不限于电子/空穴相互抵消。
积分倍增比 M 为 2.190309，原生 2.190339；仍只作积分指标，不称主特征值。
没有对剩余约 30 ppm SC 源差给出严格收敛误差上界，也未重新计算原生参考。

## 构建、回归及证据

UCRT64 Release 构建通过，检测到 UMFPACK；沿用 HDF5 状态格式，实际读写及
种子哈希检查通过。网格、300 K、A/µm 单位与上一轮一致。原点为 35 V，
本轮求解器内部漏端坐标为 0 V，对应物理内部 Vds=35 V。串行单重负载执行。

6 项相关 CTest 最终通过：外部几何及默认拒绝契约、新 SLOT 物理范围守卫、
完整 IALMob 源/装配 Jacobian、PN2D 模板契约及源探针协议。
首轮新增守卫测试仍运行构建期间旧夹具对象中的未识别迁移率名；将夹具明确
设为有效的 `masetti_lombardi` 并加入识别断言后重编，单项复验通过。
旧失败日志保留，未重复运行已通过的其他测试或全量回归。

新 runner SHA256：
`c4ab5e76bafc85915cba90a775763fb2068646b8fd1897be82863a3fb33fcc6f`。
证据目录 `build/s26transport/`：

- `frozen_acceptance.json`：默认逐字节兼容、装配复制和几何资源哈希。
- `jvp/output.csv`、`jvp_normalized.json`：真实器件方向与块归一化结果。
- `reclosure.json`、`iic/`、`sc/`：本地种子、完整更新历史、状态及独立资格审计。
- `comparison.json`：共同支撑上的原生电流、势、密度和源对比。
- `verification.json`：完整 Si 边覆盖、配置只改输运、未放宽门槛、资源哈希检查。
- `build_release.log`、`ctest.log`、`ctest_retry.log`：构建与测试记录。
- `frozen.py`、`reclose.py`、`compare.py`、`verify.py`：复核入口；重新执行需
  使用新输出目录，不覆盖已有证据。

下一步应以当前新输运配置做 38 V 独立重收敛和原生场/源对照，通过后再补齐
35→38 V 连续路径；外电阻和 BVDS 资格随后重新建立。

关联：[近表面定位](slot_ldmos_near_surface_hole_transport_2026-09-26.md)、
[上一配置连续路径](slot_ldmos_circumcentric_continuation_38v_2026-09-26.md)、
[配置契约](../config_schema.md)。
