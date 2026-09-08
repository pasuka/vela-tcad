# 局部电化学势、SG 状态电导与守恒连续性扰动验证

2026-09-06。承接[校准输运分解](simplemos_calibrated_transport_decomposition_2026-09-06.md)，完成高 NWell n23 的两个工作点、两条指定边的局部场核对，以及 26 次 Vela DC 重算。

**沟道边 320–324 的规定连续性扰动通过双幅度正负响应校准，最大预测相对误差 8.58e-8。漏接触相邻边 1091–1092 的直接端口项与状态反馈近乎抵消，净响应在原幅度及补充幅度下均未通过完整校准。** 因此，上一轮漏边贡献排名不能直接作为误差根因证据；沟道及邻近界面仍是下一轮应检查的区域。尚未证明界面势差的形成原因，未得到可进入正式修改的物理候选。

本轮没有新增 Sentaurus 运行、生产算法修改或默认设置修改，未放行 M82/M83。新增隔离程序与脚本用于诊断；原有 `src/tools/vela_example_runner.cpp` 工作区改动保留。

## 工况与方法边界

| 项目 | 本轮设置 |
| --- | --- |
| 网格 | 原 n23 Si/SiO2 网格，1482 节点；沿用原掺杂、材料及拓扑 |
| 两工作点 | Vd=0.05 V、Vg=0.9 V；Vd=1 V、Vg=0.8 V |
| 基态 | 已封存 `simplemos_strict_flux_precision` 状态，保留准费米 reference/increment |
| 模型 | M65 matched-ni、Boltzmann、无 BGN、SRH；`phumob_field_lombardi / total_impurity / quasi_fermi_gradient / transport_cell_vector`；稳定 `compensated_log_expm1` 通量 |
| 权重与导数 | 沿用已校准的完整 DD 漏端伴随及隔离的向量迁移率链式导数修正 |
| 实际线性后端 | 显式选择 SparseLU；编译可用 SuiteSparse，但本次未选用其求解后端 |
| 单位 | 电势 V、电流 A/um、边长度 m；位置标签 um；物理粒子线通量单位为粒子/(m·s) |

局部场使用已有 Sentaurus 导出，未重新运行或导入 TDR。对比原生电势、电子准费米势和浓度；Vela 状态准费米势用 reference+increment 重建。局部范围含焦点边的一跳 Si 邻域及最近 Si/SiO2 节点。

为避免把不同热电压与原始浓度混用，SG 状态代理分别按各程序热电压及电势/准费米势重建 Boltzmann 浓度，并使用相同 Vela 边几何和已匹配 ni。Decimal 60 位用于差分及对数平均，不能恢复原生导出中已丢失的精度。

令 `eta=(psi1-psi0)/Vt`，`L=n0 B(-eta)`、`R=n1 B(eta)`，定义不含迁移率的状态系数

`kappa = Vt * (couple/length) * (L-R)/log(L/R)`。

它乘迁移率及对数失衡给出相应 SG 粒子线通量。**原生侧 kappa 是原生状态代入共同 Vela 几何的代理，未导出 Sentaurus 原生边电导，也未验证其内部离散算子相同。** 因此下表不能称为两程序真实电导误差。

## 局部场核对结果

电子准费米势差及梯度均按 node1−node0 定向。

| Vd / 边 | 原生准费米势差 V | Vela 准费米势差 V | 原生 / Vela 梯度 V/m | kappa 原生态/Vela态 |
| --- | ---: | ---: | ---: | ---: |
| 0.05 / 漏边 2392，1091–1092 | -1.89678245e-10 | -2.17302769e-10 | -1.46123 / -1.67405 | 1.000004218 |
| 1 / 漏边 2392，1091–1092 | -9.983692e-11 | -1.12281603e-10 | -0.769119 / -0.864990 | 1.000004218 |
| 0.05 / 沟道边 965，320–324 | 0.00330608424 | 0.00347208860 | 150899 / 158476 | 1.024895214 |
| 1 / 沟道边 965，320–324 | 0.00897475693 | 0.00945704233 | 409634 / 431646 | 1.029839851 |

漏边两端的原生−Vela 静电势偏差约为 -0.109 / +0.118 uV；状态系数差仅约 4.2 ppm，而准费米势差明显不同。这与不同电流在低阻漏侧流出的表现相容，不能由此证明漏端无误差。其原生准费米势差只有约 1e-10 V：端点双精度 ULP 和占势差的比例分别为 7.32e-8、2.22e-6。这只是数值分辨率尺度，不是已测得的原生求解误差。

沟道边附近的原生−Vela 静电势差如下：

| 节点 | Vd=0.05 V，mV | Vd=1 V，mV |
| --- | ---: | ---: |
| 320 | +0.968084 | +0.826684 |
| 324 | +0.137115 | +0.0293903 |
| 最近界面节点 338 | +1.032831 | +0.889380 |

节点 338 到沟道边中心约 0.01129 um；漏边最近界面节点 1041 距离约 0.175 um，不能把该节点当成漏接触处的局部界面。沟道静电势差、状态系数差及准费米失衡同时存在；本轮只核对终态空间关系，尚未通过原生 Poisson 项或受控路径确定界面势形成因果，也未把它归因于迁移率。

## 守恒扰动与独立校准

在一条既有 Si 边上添加与状态无关的电子粒子线通量 `s`：Dirichlet 行替换前，node0 电子连续性残差加 `s`，node1 减 `s`。Poisson、空穴、掺杂和原有物理模型保持原值。其状态导数为零，因此不新增 Jacobian 项。

同一通量加入载流子项诊断与端口电流提取。接触边必须计入显式端口项；自由节点偶极则无此直接项。若 `lambda` 为完整 DD 伴随，接触残差行的有效权重设为零，预测为

`delta_Id = direct - (lambda0_free-lambda1_free) * s_scaled`。

1091 为漏接触、1092 为自由节点，正向注入的直接常规电流项为 `-I_inj`；320、324 均为自由节点。这里只改变指定恒定线通量，未把它解释为真实的迁移率或复合模型变化，也不能自动推广到任意连续性源。

原生内部线通量到粒子/(m·s) 的换算为 `1e4 * 1e-6 = 0.01`；相应 `I_inj[A/um]=q*1e-6*0.01*s_native`。归一化残差到 A/um 的系数为约 1.903556886e-5，另由原有边通量核对。

固定状态预检共 12 次只读调用：两个偏压下分别检查零源、漏边源、沟道源的残差与端口函数。零源残差与封存基态逐项一致；四个源的残差增量最大相对误差 3.42e-14，Poisson/空穴/掺杂/ni 未变。端口直接项符号和单位通过，六个端口函数值与独立 ContactCurrent 提取结果一致。此项验证没有用伴随预测去代替真实源和端口实现检查。

## 接受条件和调用数量

继承 max_iter=200、reltol=1e-7、abstol=1e-12、psi damping=0.35、准费米单步限制 0.025 V；逐行载流子接受 eps_row=1e-6，原有源/通量资格下限未放宽。DC 配置逐行门槛为 enforce；全局闭合在 DC 内仍为 off，随后独立只读接受检查使用 enforce、tolerance=1e-6、source_floor=1e-10，并要求 KCL/|Id|<=1e-8。

原生 Newton 收敛、独立逐行和全局接受检查必须同时通过。全局电子/空穴净源的 `qualified` 仍为 false，`satisfied` 表示按冻结的源强度门槛接受，**不表示净 SRH 源已达到百万分之一相对闭合**。

响应额外要求：预测相对误差<=0.001，双幅度差<=0.001，偶对称偏移/奇响应<=0.01，奇响应/零源漂移>=100，两符号方向正确，零源 Id 漂移<=1e-5 dex。

| 批次 | 新增 DC | 原有接受条件通过 | 附加状态门槛通过 | 响应幅度通过 |
| --- | ---: | ---: | ---: | ---: |
| 原幅度：2 零源 + 2偏压×2边×4符号幅度 | 18 | 18 | 未增加 | 沟道 4/4；漏边 0/4 |
| 漏边独立补充幅度：2偏压×4符号幅度 | 8 | 8 | 6/8 | 0/4 |

26 次 DC 均有独立只读接受检查；另有上述 12 次固定状态预检。不存在新增 Sentaurus 求解或新增伴随求解。响应通过数按 full/half 计数，每个幅度均需要正负两次 DC，不能将 4/4 理解为仅四次求解。

## 沟道边响应：通过

执行前冻结 full/half 注入幅度为基态 Id 的 0.05% / 0.025%。

| Vd | full 预测 A/um | full 实测奇响应 A/um | full / half 预测相对误差 | 双幅度相对差 | 响应/注入电流 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0.05 | 1.250865214e-12 | 1.250865321e-12 | 8.58e-8 / 8.40e-8 | 1.81e-9 | 0.139948 |
| 1 | 1.505911301e-12 | 1.505911339e-12 | 2.53e-8 / 1.62e-8 | 9.14e-9 | 0.326348 |

两偏压所有符号、漂移、偶分量、收敛条件均通过。这补齐了 **Vela 指定沟道边连续性扰动** 的独立同扰动校准。它没有校准 Sentaurus 连续性源，也没有证明有限幅度的真实物理候选会沿相同方向响应。

## 漏边响应：直接项与反馈抵消，净响应未获资格

原 full/half 注入幅度为 Id 的 5% / 2.5%。full 预测分项为：

| Vd | 直接项 A/um | 状态反馈 A/um | 净预测 A/um |
| --- | ---: | ---: | ---: |
| 0.05 | -8.93807163045e-10 | +8.93807161844e-10 | -1.201846591e-18 |
| 1 | -4.61443965103e-10 | +4.61443965026e-10 | -7.679348352e-20 |

净增益仅约 -1.34e-9 / -1.66e-10。原幅度 full 实测为 -1.20177e-18 / -7.68178e-20 A/um，但奇响应/零源漂移仅 3.92 / 0.0715，偶分量比约 1.95 / 1002，均不满足门槛。只看奇响应与预测相近会误判。只看 `-lambda·delta_F` 又会遗漏几乎等大的直接端口项。

保留原始失败后，另行冻结补充控制：选取使预测净响应达到至少 `1e-7 Id` 的最小十进制幂注入，即低 Vd full/half=100/50 Id，高 Vd=1000/500 Id。幅度从冻结预测计算，未按新增 DC 拟合。另加所有节点 psi/电子及空穴准费米变化<=0.001 V 的限制，其余原门槛不变。

| Vd | full / half 预测相对误差 | 双幅度相对差 | full / half 偶分量比 | 结论 |
| --- | ---: | ---: | ---: | --- |
| 0.05 | 0.1895% / 0.0495% | 0.1398% | 0.000519 / 0.001835 | full 误差及双幅度差未通过 |
| 1 | 0.8517% / 5.1318% | 4.2439% | 0.04778 / 0.08379 | 误差、双幅度、偶分量未通过；full 状态变化超限 |

高 Vd 正负 full 状态最大变化约 1.259 / 1.999 mV，分别发生在节点 947 / 958 的空穴准费米势。两次 Newton 及原有接受检查仍通过，失败的是新增 1 mV 状态限制。半幅度最大变化约 0.215 / 0.573 mV。不得将这两种资格混为一谈。

补充控制改善了信号/零漂比，但未使净响应获得完整资格。强抵消后的净量对数值闭合非常敏感；现有 KCL 接受额度也不足以单独保证约 1e-15 A/um 净量的千分之一精度。反馈单项可由实测净响应减去直接项复核，不能因此放行抵消后的净响应。本轮停止放大幅度，没有放宽响应或状态门槛。

## 结论、剩余问题与下一步

本轮已排除把漏边空间排名直接当作成因的推断：必须同时计算直接端口项和自洽反馈。数据支持其主要体现电流流出位置的解释，但漏边净响应未完全校准，不能宣布漏侧因素已被彻底排除。

后续应优先围绕已校准的沟道边及节点 338 附近，将静电势形成分解到定义明确的局部 Poisson 电荷/介电项，并与沟道输运状态系数变化作对照；候选仍须覆盖高低 NWell 配对与高 Vd 反例。当前两点高 NWell 结果不替代 n19 控制及完整曲线。

若继续追究漏边微小净响应，需要先建立针对直接项/反馈抵消的误差预算，验证更严格闭合或更高精度能否降低零源漂移和偶分量，不能继续单纯加大注入。原生连续性源的同扰动校准、真实原生边电导及界面势形成原因均仍待验证。当前没有足够证据修改正式电流算法。

## 可复核证据

- [冻结输入与原幅度约定](../../reference_tcad/simplemos_sentaurus2022/local_conservative_flux/contract.json)、[输入哈希](../../reference_tcad/simplemos_sentaurus2022/local_conservative_flux/freeze.json)。
- [局部边场](../../reference_tcad/simplemos_sentaurus2022/local_conservative_flux/local_fields.csv)、[局部节点场](../../reference_tcad/simplemos_sentaurus2022/local_conservative_flux/local_nodes.csv)。
- [原幅度 DC](../../reference_tcad/simplemos_sentaurus2022/local_conservative_flux/dc.csv)、[校准明细](../../reference_tcad/simplemos_sentaurus2022/local_conservative_flux/calibration.csv)。
- [补充幅度约定](../../reference_tcad/simplemos_sentaurus2022/local_conservative_flux/contact_resolution/contract.json)、[补充校准](../../reference_tcad/simplemos_sentaurus2022/local_conservative_flux/contact_resolution/calibration.csv)。
- [原始状态及接受检查复核](../../reference_tcad/simplemos_sentaurus2022/local_conservative_flux/review.json)、[完整证据清单](../../reference_tcad/simplemos_sentaurus2022/local_conservative_flux/evidence.json)。
- 脚本：`scripts/validate_simplemos_local_conservative_flux.py`、`scripts/resolve_simplemos_contact_flux_signal.py`、`scripts/audit_simplemos_local_electrochemical_fields.py`、`scripts/review_simplemos_local_conservative_flux.py`。隔离源、可执行文件、编译命令及逐点原始输出位于 `build-release/simplemos_local_conservative_flux_20260906/`。

验证记录：隔离 C++20/O2 程序编译成功；`tests/regression/test_simplemos_local_conservative_flux.py` 的 5 个单位/守恒/方向测试通过；12 次固定态预检、26 次 DC 和各自接受探针已完成。复核脚本从原始 JSON 重新检查 26 个状态、六个端口提取一致性及 12 个幅度的通过/失败原因；未新增生产代码，未运行全量 CTest。

在本工作树使用 `D:/msys64/ucrt64/bin/python.exe scripts/validate_simplemos_local_conservative_flux.py verify` 检查封存证据。已冻结目录供核验，不应直接重跑 prepare 或改写历史失败记录。
