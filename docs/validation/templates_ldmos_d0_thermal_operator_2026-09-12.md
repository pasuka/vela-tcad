# D0独立热传导实现及原生热源对照

本轮已实现独立稳态热传导组装及诊断入口，尚未接入电热自洽求解。
下述通过结果使用固定原生热源，不能认领完整D0电流曲线通过。
用户已批准[热验收合同](../../reference_tcad/templates_ldmos_sentaurus2022/thermal/d0_thermal_acceptance.json)：
全域热平衡≤0.1%；峰值温升及温升场RMS误差分别≤max(1 K,原生对应温升量的5%)；
原电学门限不变。零功率不计算相对误差，另查恒温零源。

## 实现和验证范围

[`LatticeHeatAssembler`](../../include/vela/equation/LatticeHeatAssembler.h)
采用Tri3单元、显式区域导热率、单元恒定体热源及Robin表面热阻边界。
残差单位W/m、Jacobian单位W/(m K)，均按单位器件宽度；坐标到米的
换算必须显式提供。导热率取单元平均温度，并包含dk/dT导数。
热边界与电学端口分开，其他外边界绝热；拒绝内部边、重复热边和无材料参数的区域。

[`lattice_heat_probe`](../../src/tools/lattice_heat_probe.cpp)仅求解或回放指定热源，
不改变电学求解器默认设置。冷热材料、界面热流、Robin边界、全域热平衡、
SI/um等价、导热率Jacobian、恒温零源及均匀体热源的一维解析收敛性均有测试。
[`热评分脚本`](../../scripts/analyze_templates_ldmos_thermal.py)要求完整一致的节点
映射，拒绝非有限量和无效体积，不用峰值正确代替全场检查。

## 原生诊断与热源口径

在已授权本机Sentaurus虚拟机独立目录
`/home/tcad/sentaurus_runs/vela_oracle_2022/d0_heat_diagnostics_20260912_r1`
执行6个T-2022.03-SP2诊断：Vg=4/8 V、Vd=40 V各自进行只加载输出、
原模型重闭合、关闭hRecVelocity后重闭合。原文件未覆盖，脚本和参数哈希
见`reference_staging/templates_ldmos_auger_followup_20260912/thermal_native_r1/plan.json`。

原D0只求解Temperature，不启用Thermodynamic、Peltier、RecGenHeat。
用户指南p.248方程71及p.253--255用于区分默认热方程和可选源项。
实际输出显示：该模型的eJEHeat/hJEHeat/ColHeat数据均为零，不能用它们
拼出当前热源；仅Load/Plot也未重算TotalHeat。执行一次原模型Coupled后
才得到可用TotalHeat，Vg4温度改变仅约1e-12 K。

第一次按零分项拼接热源的试验回到300 K、未通过热门限，完整保留。
后续使用重闭合后的TotalHeat，仅向Si单元沉积：理想绝缘氧化层没有电学
热源，界面共享节点在氧化层绘图数据中的非零值不能作为新增热源。
采用原生节点值的单元平均，因此这是绘图热源转移实验，未证明与原生
求解器内部热源离散完全相同。Vg4转移源积分4731.76 W/m，日志总热约4724 W/m。
不能用分项绘图之和替代原生内部方程；用户指南p.250也明确区分二者精度。

## 固定原生热源的热传导对照

同一10241节点/Si-SiO2网格；坐标um→m；Si
k(T)=100/(-0.0393+0.00155T+1.82e-6T²) W/(m K)，氧化层1.4 W/(m K)。
th_lat=300 K，表面热阻5e-7 m²K/W，101条边；体积权重使用全域Tri3面积/3。
采用UCRT64 Release、Eigen SparseLU。此处温度回放不包含电学Newton。

| Vg / Vd | 原生峰温K | Vela峰温K | 峰值误差K | 温升场RMS误差K | 热平衡相对误差 | 热Newton更新 |
|---|---:|---:|---:|---:|---:|---:|
| 4 / 40 V | 401.112357 | 401.032237 | 0.080120 | 0.105174 | 2.50e-15 | 2 |
| 8 / 40 V | 513.568615 | 512.186451 | 1.382164 | 0.404779 | 1.74e-15 | 3 |

两点通过已批准热学门限，实际网格的恒温零源检查也通过。
Vg8峰值门限为10.6784 K，不是固定1 K。机器结果分别在
`thermal_source_probe_r2_totalheat/closed`和`thermal_source_probe_vg8_r3/closed`。
这里的热平衡比较转移后源与Vela边界散热，不是声称原生绘图积分与其
内部热源达到同样误差。

## 热态接触复合端点对照

关闭hRecVelocity与原热态重闭合控制相比：

| 栅压 | 漏电流相对变化 | 峰温变化K |
|---|---:|---:|
| 4 V | 3.419995e-7（0.00003420%） | 9.02916e-5 |
| 8 V | 3.046043e-7（0.00003046%） | 1.45242e-4 |

两个40 V端点变化很小，不能外推为完整热态消融曲线或通用接触模型验收。
数据见`thermal_native_r1/hot_hrec_summary.json`。

## 剩余工作

接入逐节点温度相关统计、BGN/材料、SRH/Auger、低场/高场IALMob，构造
与实际电流离散一致的热源及电热耦合Jacobian，再从零压开展D0曲线。
先检查冻结300 K严格复现等温资格；按低/中/高偏压检查热场、电学块、
逐行闭合、KCL和功率，最后扩展双栅压0--40 V。不得把固定原生热源测试
中的合格温度场替代上述自洽实现。

局部温度相关低/高场 IALMob 已在[后续实现](templates_ldmos_d0_mobility_temperature_2026-09-12.md)
中完成双栅压热态核对；上述剩余工作中的主方程接入仍未完成。
