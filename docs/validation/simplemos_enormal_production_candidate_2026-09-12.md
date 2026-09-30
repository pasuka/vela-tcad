# SimpleMOS Enormal 生产候选与同扰动验证

日期：2026-09-12。Enormal 显式候选及新增低栅压控制已通过十六点、32 状态及 16/16 双初始化；四工况双幅度原生及 Vela 同扰动校准全部通过。原生完整曲线 204/204 点合格，Vela 完整曲线正在计算。本记录承接已完成的 [PhuMob 四条完整曲线](simplemos_phumob_fullcurve_validation_2026-09-12.md) 和 [Enormal 原生形成方式校准](simplemos_enormal_native_calibration_2026-09-12.md)，不以 PhuMob 曲线代替 Enormal 验收。

## 候选实现和范围

显式选择 `phumob_lombardi`、`element_box_phumob` 和 `surface.discretization=element_distance_gradient`。每个硅单元内使用电势梯度与未单位化的界面距离梯度点积，结合顶点距离和实时 PhuMob 迁移率形成 Lombardi 倒迁移率，再按原生 box 权重形成单元及边迁移率。

同一输入作用于装配、split 高低位状态求值、端口和诊断；Jacobian 包括全部单元邻接电势列、PhuMob 的两种载流子链式偏导和守恒散布。该候选限定 300 K、alpha=0、默认界面选择；高场饱和仍关闭。普通路径默认设置不变，q、kb、ε0、SRH 体积和原接受门槛不变。

采用原 n19/n23 网格，Vd=0.05/1 V、Vg=0.8/1 V；坐标 μm，浓度 m⁻³，电流 A/μm，等效宽度 1 μm。Vela 为 UCRT64 Release、Eigen SparseLU/COLAMD、四次线性修正和显式 split 状态；原生为 T-2022.03-SP2、ExtendedPrecision(128)、Super。双初始化分别使用相同偏置的合格 PhuMob 状态和从原生势场重建的相容密度状态。原生原始密度另存，仅用于固定公式重放。

## 已完成的本地检查

- 一致构建后的八个相关测试程序全部通过：207 个用例、4561 项断言。新增三项单元测试中的 960 项断言覆盖全部 15 列、多步长、弱非局部交叉项、少数载流子、低于 1 ULP 的响应、高低位分拆不变性及守恒端口。此前编译过程中使用不一致对象的预备测试失败日志仍保留；最终结果来自后续一致构建。
- 八个真实工作点的固定状态 C++ 单元迁移率重放全部通过，最大相对差 3.331e-15，对照对象是独立计算的同一数学截断公式。该检查不消除既有的原生空穴截断差异。
- 首个工况在节点 1000 的 ψ、φn、φp 三列，使用 1e-6/1e-20 V 两档扰动，18 个真实归一化分块检查全部通过，最大相对差 2.496e-10。
- 扩展的八态逐行列检查全部通过：144 组、213264 个行比较，784 个非零项，最大相对差 5.027e-9，无弱项豁免。该范围是每个真实矩阵的节点 1000 三列，不代表全部真实矩阵列已审计；小型受控网格另有全列测试。证据见 [逐行列检查](../../reference_tcad/simplemos_sentaurus2022/enormal_candidate_20260912/v1/real_column_evidence.json)。

候选源代码、二进制和测试已独立归档：[归档身份](../../reference_tcad/simplemos_sentaurus2022/enormal_candidate_20260912/v1/archive_evidence.json)。固定状态与初步导数证据分别见 [fixed_evidence](../../reference_tcad/simplemos_sentaurus2022/enormal_candidate_20260912/v1/fixed_evidence.json) 和 [smoke_evidence](../../reference_tcad/simplemos_sentaurus2022/enormal_candidate_20260912/v1/smoke_evidence.json)。原 PhuMob 合格程序仍保留在上一阶段归档，不受本候选构建覆盖。

## 原生小扰动参数方向的纠正

第一批误将原生 `a_ac/a_sr` 作为直接倒迁移率乘子；20 个原生点正常完成，但正负两档电流完全不变。T-2022.03 手册第 995 页公式 1074 表明它们是应力增强权重，当前无应力模型时不起该作用。原分析器因此出现零斜率除法错误；该批次明确记为响应资格失败，保留 [失败账本](../../reference_tcad/simplemos_sentaurus2022/enormal_response_20260912/failure_summary.json) 和原日志。

修正方向为共同倒散射因子 f=1+δ：将原生 B、C、delta、eta 同时除以 f，两种载流子都相同。这样声学与粗糙散射的倒迁移率都乘以 f；不改变原生应力参数。Vela `surface.acoustic_factor/roughness_factor` 是直接诊断乘子，不与原生同名应力权重等同。

修正批次先在 n19、Vd=0.05 V、Vg=0.8 V 做 δ=0、±0.001、±0.0005 五点，再根据先导结果决定另外三个工况。检查保留零扰动身份、非零信号、双幅度斜率、偶/奇响应和相对零漂移门槛；零信号不可记为线性度通过。输入合同见 [material_scale_v2](../../reference_tcad/simplemos_sentaurus2022/enormal_response_20260912/material_scale_v2/native_contract.json)。

最终四工况的 20 个原生状态、4 个零幅度身份检查、8 组幅度响应检查全部通过。先导零扰动电流、势场、密度和迁移率与基准完全一致；原生归一化 Id 响应为 -0.2742860029/-0.2742859857。四工况最大双幅度斜率相对差 6.298e-8，最大偶/奇比 2.868e-4。零漂移记为“未观测到”，不解释为无限导出精度；跨求解器核对另按原生 15 位有效数字计入噪声下限。证据见 [先导](../../reference_tcad/simplemos_sentaurus2022/enormal_response_20260912/material_scale_v2/pilot_evidence.json) 与 [其余工况](../../reference_tcad/simplemos_sentaurus2022/enormal_response_20260912/material_scale_v2/rest_evidence.json)。

## 自洽试算

首个工况在本工况原生响应和逐行列检查通过后，双初始化试算通过。两种初始化各用两次 Newton 更新；最大行比 1.386e-9、KCL/Id 小于 2.192e-16、端口提取相对差小于 2.221e-16。两初始化 Id 相对差 1.910e-14，ψ、φn、φp 最大差分别为 5.88e-15、4.439e-12、6.312e-13 V，密度相对差 1.719e-10；相对原生电流差 +0.003387987%。证据见 [先导双初始化](../../reference_tcad/simplemos_sentaurus2022/enormal_candidate_20260912/v1/pilot_dual_evidence.json)。

最终八点批次 **16/16 状态、8/8 双初始化通过**，16 次首次尝试，无失败、无重载。最大逐行比值 3.432835e-8、最大端口提取相对差 2.221e-16；最大双初始化 Id 相对差 2.045e-11，ψ/φn/φp 最大差分别为 1.719e-12/3.790e-11/7.467e-11 V，密度相对差 2.889e-9。这些场差是两种 Vela 初始化之间的差异，不是对原生场的误差。

误差定义为 100×(Id_Vela/Id_Sentaurus−1)，单位 %：

| 器件 | Vd (V) | Vg=0.8 V | Vg=1.0 V |
|---|---:|---:|---:|
| n19 | 0.05 | +0.003387987 | +0.001045113 |
| n19 | 1.0 | +0.002964473 | +0.001306597 |
| n23 | 0.05 | +0.005757409 | +0.005136330 |
| n23 | 1.0 | +0.006605273 | +0.005178238 |

完整数值与证据见 [八点对照](../../reference_tcad/simplemos_sentaurus2022/enormal_candidate_20260912/v1/comparison.csv)、[汇总](../../reference_tcad/simplemos_sentaurus2022/enormal_candidate_20260912/v1/summary.json) 和 [资格证据](../../reference_tcad/simplemos_sentaurus2022/enormal_candidate_20260912/v1/comparison_evidence.json)。跨求解器扰动使用同一已冻结程序，最终 Vela 20/20 状态首次通过，无失败、无重载；8/8 幅度对照通过。最大原始电流响应相对差 6.398122e-5（0.006398122%），最大归一化响应相对差 3.156430e-6，低于原 1e-3 门槛。没有用归一化结果掩盖原始响应，也没有使用拟合截断值。见 [同扰动汇总](../../reference_tcad/simplemos_sentaurus2022/enormal_same_response_20260912/summary.json) 和 [证据](../../reference_tcad/simplemos_sentaurus2022/enormal_same_response_20260912/comparison_evidence.json)。

## 合格八点的原生物理场对照

全部硅区（含接触）的势场使用保存状态的高低位重构，密度和 SRH 使用同一状态重新求值。八态密度重构在 double 导出精度内一致；对实际生产积分源的 SRH 重构最大 L1 相对差 3.194e-16。首版误用排除接触的验收掩码，触发集合断言；纠正为原生硅区全集后通过，未删去接触或热点。历史 `all_row` 的 `*_m3` 列在当前 unit-scaling 下实际为 cm⁻³；状态 CSV 的密度列才是 m⁻³。

n23 的全硅区最大绝对场差如下；这些是 Vela 对 Sentaurus 的差异，区别于上面的双初始化指标：

| Vd (V) | Vg (V) | ψ (μV) | φn (μV) | φp (μV) | 电子密度最大对数差 (dex) | 空穴密度最大对数差 (dex) |
|---|---:|---:|---:|---:|---:|---:|
| 0.05 | 0.8 | 2.483 | 0.1414 | 38.69 | 0.00004248 | 0.0006504 |
| 0.05 | 1.0 | 2.630 | 0.1322 | 38.72 | 0.00004485 | 0.0006508 |
| 1.0 | 0.8 | 3.372 | 3827 | 1334 | 0.06425 | 0.02242 |
| 1.0 | 1.0 | 3.556 | 3827 | 1226 | 0.06425 | 0.02063 |

高 Vd 的电子准费米势最大差位于节点 792，空穴最大差位于 1057/795；电流吻合不代表所有局部场同样吻合。按共同正的重心硅面积权重比较 SRH 速率，n23 四点 L1 相对差为 4.490e-6 至 1.530e-5；相应电荷等价积分差为约 -6.257e-22 至 +1.006e-22 A/μm。这是共同体积下的速率差，不能当成各求解器自身源积分、漏电流差或因果响应预测。完整场表、节点值、重构与证据见 [字段汇总](../../reference_tcad/simplemos_sentaurus2022/enormal_fields_20260912/fields.csv)、[SRH](../../reference_tcad/simplemos_sentaurus2022/enormal_fields_20260912/srh.csv) 和 [证据](../../reference_tcad/simplemos_sentaurus2022/enormal_fields_20260912/evidence.json)。

例如 n23、Vd=1 V、Vg=0.8 V，节点 792 的原生/Vela 电子密度为 1137.19/980.80 cm⁻³，空穴密度为 22561.19/22014.81 cm⁻³；SRH 为 -4.215424822e16/-4.215424647e16 cm⁻³s⁻¹。节点 1057 的原生/Vela 空穴为 7.6233/8.0271 cm⁻³，电子约 2.352e13 cm⁻³。局部较大的准费米势和相对密度差，与很小的速率差可以并存；在低载流子密度区，SRH 生成率接近饱和形式，对小密度变化并不同比例响应。这一公式解释不能替代局部源或输运方向的自洽因果校准。

## 后续放行条件

导数、原生及 Vela 同扰动、八点双初始化已经完成。已按授权启动四条原生 Enormal 0–1 V、步长 0.02 V 曲线，输入合同见 [曲线合同](../../reference_tcad/simplemos_sentaurus2022/enormal_curves_20260912/native_contract.json)。原生曲线已通过 204/204 点资格，新增 Vg=0/0.2 V 八点也全部通过，形成十六点双初始化；Vela 四条完整曲线已经启动。进展见 [十六点与完整曲线](simplemos_enormal_fullcurve_validation_2026-09-12.md)。逐行、端口、全局源闭合及双初始化门槛保持原值，保留全部失败尝试；之后恢复高场饱和及原始 n17–n24 全部工况。原生内部 PhuMob 截断算法仍有独立未完成项，不将诊断拟合值写入生产。
