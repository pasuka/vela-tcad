# SLOT LDMOS：Sentaurus 与当前 Vela 仿真脚本对照

核对日期：2026-09-21。Vela 源码基线：`76817440`，工作分支：
`codex/slot-ldmos-bvds-next`。

当前 Vela 已采用原始参数的元素 IALMob，并实现迁移率及自洽雪崩源对局部
势变量的完整链式导数。300 K 的主要迁移率、SRH 和雪崩电离系数已对齐，
但整个算例还不能视为与 Sentaurus 等价：Auger 密度增强未启用，材料与
多晶硅/接触处理有明确差异，AverageBox 和雪崩源项离散也未完全对齐。
外电阻功能已实现，当前 SLOT 的本地击穿支路仍未通过验收。

本文整理输入、实际启用行为和已有验证证据；不将参数文件中所有模型、
Plot 输出名称或历史实验配置视为当前启用的物理模型。

**对照基准与可追溯性**

原 Workbench 路径 `D:/gpt-tmp/slot-LDMOS/V3E/_NLDMOS/_40V/_Copy`
当前不可用。本次使用虚拟机中已执行并冻结的原始物理参考 deck、原始
`pp2_des.par` 和对应运行日志；不声称该参考 deck 与最早 Workbench 输入
逐字节相同。以下比较限定在 SDevice 器件仿真阶段，SProcess 的工艺结果
作为网格/掺杂输入，不涉及工艺步骤的重新实现。

| 编号 | 核对对象 | 用途 |
| --- | --- | --- |
| S1 | [冻结参考 deck](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/build/slot_deck_comparison_20260921/original_reference/bvds_external_resistor_final.cmd) | 原始物理、外电阻、原始扫描与停止条件 |
| S2 | [原始 pp2_des.par](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/build/slot_deck_comparison_20260921/original_reference/pp2_des.par) | 模型参数；SHA256 `8f3c94fbca5ac9cd462dadee00549e032466a72bfb02fa092e4d3bac9d2e2919` |
| S3 | [冻结参考实际运行日志](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/build/slot_deck_comparison_20260921/original_reference/bvds_external_resistor_final.log_des.log) | 接触、激活模型、精度、AverageBox 等实际行为 |
| S4 | [9 月 21 日原生外电阻对照 deck](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/build/slot_element_ialmob_20260921/oracle/slot_external_bv.cmd) | 保留原物理，显式 Derivatives、分段续算和元素迁移率输出 |
| S5 | [近期原生运行日志](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/build/slot_element_ialmob_20260921/oracle/completed/slot_external_bv_console.log) | 最近外电阻对照的执行证据 |
| V1 | [当前合格 10 V 雪崩点完整配置](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/build/slot_element_ialmob_20260921/av_10v_umfpack/outputs/stages/02_delaunay_avalanche_adaptive/attempt_001_v_10/local_corrector/validation_003/config.json) | 当前元素 IALMob 物理和耦合求解基线 |
| V2 | [当前 10 V 外电阻尝试配置](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/build/slot_element_ialmob_20260921/external_10v_umfpack/config.json) | 独立的外电阻功能验证；不能当作已合格支路 |
| V3 | [材料覆盖文件](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/build/low_voltage_20260906/intrinsic_aligned_materials.json) | ni 校准及 PolySilicon/Nitride 定义 |
| E1 | [已有执行与验收记录](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/docs/validation/slot_ldmos_element_ialmob_execution_2026-09-21.md) | 检查点、Jacobian、预算及失败边界 |
| E2 | [本次证据哈希清单](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/build/slot_deck_comparison_20260921/evidence_manifest.json) | 输入快照及 IALMob 参数逐项一致性结果 |

S1/S2/S3 取自虚拟机目录
`/home/tcad/sentaurus_runs/vela_oracle_2022/slot_ldmos_bvds_t2022_20260819_run01/payload_3c32b254/`。
上述 `build/` 文件是本地忽略的证据文件，不随本报告纳入版本管理。
“当前 Vela”以 V1 为准；无雪崩扫描与外电阻验证分别使用独立配置，不能
用最早的 `ialmob_field` 配置代表本轮实现。

**功能与启用状态总表**

| 功能 | 原始 Sentaurus | 当前 Vela | 对照结论 |
| --- | --- | --- | --- |
| 仿真类型 | 二维、等温、稳态 Poisson + Electron + Hole；Quasistationary 扫漏极 | `dc_sweep`、`bv_reverse`；耦合漂移扩散，检查点及续算 | 类型相同；目前本地没有完整连续击穿扫描 |
| 几何与掺杂 | `n1_fps.tdr` | 导入后的 `mesh.json`、逐节点 `doping.csv` | 13,881 节点、27,275 三角形；坐标和拓扑已核对，材料/掺杂适用范围仍有差异 |
| 温度 | `Temperature=300` | 默认 300 K | 当前等温温度一致 |
| 统计 | Boltzmann，未启用 Fermi | Boltzmann | 一致；Fermi 仅做过单因素灵敏度实验 |
| 带隙窄化 | `EffectiveIntrinsicDensity(OldSlotboom)` | `bandgap_narrowing: old_slotboom` | 模型主参数一致；不等于全部材料能带参数一致 |
| 低场/界面迁移率 | `Enormal(IALMob)`；日志为 FullPhuMob、PhononCombination=1 | `mobility.model: ialmob`，元素核，unnamed isotropic 参数 | 已对齐；不再是旧 `ialmob_field` 或 Enhanced Lombardi 替代模型 |
| 高场迁移率 | `HighFieldSaturation` | 开启 HFS，准费米梯度驱动 | 300 K 参数一致，见下表 |
| 低密度场插值 | 原 deck 的 RefDens 行被注释 | IALMob `reference_density_m3=0`；雪崩 RefDens=0、插值 none | 当前禁用状态匹配；不能继承其他模板的非零默认值 |
| SRH | `SRH(DopingDep)` | SRH + 掺杂依赖寿命 | 300 K 寿命参数一致；温度依赖开关不同 |
| Auger | 普通 `Auger`，无 WithGeneration；S2 含 H/N0 密度增强 | `auger_with_generation:false`；仅常数 Cn/Cp，未配置密度增强 | 生成开关已修正，密度增强仍未对齐 |
| 雪崩 | `Avalanche(GradQuasiFermi)`，van Overstraeten；`AvalDensGradQF`、`ElementVolumeAvalanche` | 自洽 van Overstraeten，元素 SG/GSS-Laux 电流，顶点 box 源映射 | 电离系数一致；产生率的离散构造不完全相同 |
| Jacobian | 原 deck `#Derivatives`；实际日志明确 With avalanche derivatives | 迁移率场导数开启，元素九列导数；雪崩 `source_jacobian: local_ad` | 本地导数已实现和审计；原生导数状态应以日志结合 deck 判断 |
| 接触与多晶硅 | 五个接触运行日志均为 Ohmic；原生 PolySilicon 有掺杂及载流子 | source/drain/substrate 为 ohmic；gate/SLOT 为 metal_gate；PolySilicon 仅静电 | 明确物理差异，尚未由同偏压验证关闭 |
| 控制体积/边耦合 | 原生日志 `CVPL_AverageBoxMethod=TRUE` | assembled-Delaunay 回退边耦合 + mixed-Voronoi 节点体积 | 保留本地兼容几何；不能声称已整体替换成原生 AverageBox |
| 外电阻 | 漏极 R=1e12，随外电源升压 | 已实现 `series_resistor`、`nested_scalar`，R=1e12 Ω·µm | 参数可对应；本地 SLOT 非零步进仍未合格 |
| 热、量子、BTBT、陷阱 | 当前 Physics/Solve 未启用 | 当前配置未启用 | 这些不是本次原始物理对齐所缺失的激活项 |

S3 全局 Physics 摘要中的 “Without SRH/Auger” 不能单独代表 Silicon：
后续 `Material="Silicon"` 覆盖明确启用了 SRH/Auger/IALMob/雪崩。
参数文件列出模型、日志列出参数改写，也不单独证明该模型已经参与求解。

**材料、电极和边界参数**

| 项目 | Sentaurus 参考 | Vela 当前值/行为 | 判断 |
| --- | --- | --- | --- |
| 源、衬底电压 | 0 V | 0 V | 相同 |
| gate、SLOT 电压 | 0 V，Ohmic 接触 | 0 V，metal_gate，`flatband_voltage=0` | 电压数字相同，边界物理不等价 |
| 漏极 | 初始 0 V，外电阻升压 | V1 的 `bias_points:[10]` 实际覆盖接触初始 0 V；V2 扫外电源 | 需要区分初值、实际点和内外电压 |
| Si 相对介电常数 | 11.7 | 11.7 | 相同 |
| Si 带隙，300 K | 1.12416 eV | 默认 1.12 eV | 相差 −0.00416 eV |
| Si 电子亲和势 | 4.0727 eV | 默认 4.05 eV | 相差 −0.0227 eV |
| Si 本征浓度 | 由原生材料模型得到的校准目标 | 显式 `ni=1.0749877484352535e10 cm^-3` | 已单独校准 ni；不能再直接用带隙差推算 ni 或漏电差 |
| SiO2 相对介电常数/带隙 | 3.9 / 9 eV | 3.9 / 9 eV | 相同 |
| SiO2 电子亲和势 | 0.9 eV | 默认 0.95 eV | 相差 +0.05 eV；对当前解的影响未单独量化 |
| PolySilicon | 原生保存态有非零掺杂与载流子 | eps=11.7、Eg=1.12 eV、Chi=4.05 eV；没有 ni/迁移率/DOS，载流子行不参与输运；当前掺杂置零 | 结构保留不等于半导体多晶硅物理保留 |
| Nitride | 结构绝缘区 | eps=7.5、Eg=5 eV、Chi=1.9 eV，静电材料 | 本文未给出原生该材料全部默认参数的逐项一致性结论 |
| Ohmic 重构 | 原生接触实现 | `contact_boundary_reconstruction: legacy_node_local` | 尚无完整等价证明 |

材料默认值见 [MaterialDatabase.cpp](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/src/material/MaterialDatabase.cpp:253)，
原生 Si 带隙/亲和势见 [S2](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/build/slot_deck_comparison_20260921/original_reference/pp2_des.par:1249)，
Ohmic 类型见 [S3](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/build/slot_deck_comparison_20260921/original_reference/bvds_external_resistor_final.log_des.log:235)。

原生 2 V 元素回放的 PolySilicon 区域共 350 个节点、549 个三角形：
电子密度约 `8.02e2～4.81e19 cm^-3`、空穴约 `4.51～4.65e18 cm^-3`，
施主约 `4.20e13～3.34e19 cm^-3`、受主约 `3.01e12～1.00e19 cm^-3`。
对应的本地逐节点掺杂文件在这 350 个节点均为 Nd=Na=0。这是本次发现的
实际输入差异，而不只是材料名称不同；其对器件解的贡献尚未做独立控制。
原生证据见 [PolySilicon 电子密度](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/build/slot_element_ialmob_20260921/native_elements/fields/eDensity_region3.csv)
及同目录 DonorConcentration/AcceptorConcentration/hDensity 区域文件。

**IALMob 参数完整对照**

S2 中生效的是 unnamed isotropic `IALMob` 参数，不能替换为命名的
`IALMob "110"/"111"` 或模板 AutoOrientation 参数。当前 Vela 的 100、110、111
三组均填入相同的 unnamed 参数，因此此配置没有引入晶向差异。

下表为原生电子/空穴参数；当前 Vela 两类载流子、三组参数逐项一致。
`electron_parameters_cm`/`hole_parameters_cm` 使用原生 cm 制参数约定，
不同经验系数的量纲不同，不能把整张表都按 cm²/(V·s) 解释。

| 参数 | 电子 | 空穴 | Vela 对照 |
| --- | ---: | ---: | --- |
| mumax，cm²/(V·s) | 1417 | 486 | 相同 |
| mumin，cm²/(V·s) | 56 | 46 | 相同 |
| alpha | 0.68 | 0.719 | 相同 |
| theta | 2.285 | 2.247 | 相同 |
| n_ref，cm⁻³ | 9.68e16 | 2.23e17 | 相同 |
| B | 9e5 | 8.1e5 | 相同 |
| C | 4.84e4 | 8.36e3 | 相同 |
| delta | 7.15e14 | 2.58e14 | 相同 |
| eta | 1e50 | 1e50 | 相同 |
| N1 | 1 | 1 | 相同 |
| N2 | 1 | 1 | 相同 |
| lambda | 0 | 0 | 相同 |
| lambda_sr | 0 | 0 | 相同 |
| A | 2 | 2 | 相同 |
| alpha_sr | 5e-22 | 0 | 相同 |
| k | 1 | 1 | 相同 |
| nu | 0 | 0 | 相同 |
| D1_inv | 300 | 150 | 相同 |
| D2_inv | 300 | 150 | 相同 |
| nu0_inv | 1.5 | 1.5 | 相同 |
| nu1_inv | 2 | 2 | 相同 |
| nu2_inv | 0.5 | 0.5 | 相同 |
| alpha1_inv | 0 | 0 | 相同 |
| alpha2_inv | 0 | 0 | 相同 |
| l_crit，cm | 1e-6 | 1e-6 | 相同，等于 10 nm |
| l_crit_c，cm | 1e-6 | 1e-6 | 相同，等于 10 nm |
| mass，m/m0 | 1 | 1.258 | 由 me_over_m0/mh_over_m0 映射 |
| other_mass，m/m0 | 1.258 | 1 | 对方载流子质量映射 |
| EnormMinimum | 0 | 0 | 当前核按零下限路径匹配，非独立 JSON 参数 |

FullPhuMob 的共同筛选参数 `nref_D=4e20`、`nref_A=7.2e20 cm^-3`、
`cRef_D=0.21`、`cRef_A=0.5` 由当前核内默认值提供，与 S2 相符。
参数映射文件为 [parameters.json](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/build/slot_element_ialmob_20260921/local_control_r2/parameters.json)，
原始块从 [S2 第 1819 行](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/build/slot_deck_comparison_20260921/original_reference/pp2_des.par:1819) 开始。

当前 `effective_electrodes=[source,drain,substrate]`，
`contact_electric_field_fallback=false`，局部几何文件为
`local_control_r2/local_geometry.json`。相同保存态下全部 24,272 个 Si 单元
迁移率回放的最大相对差：电子 `2.98469e-5`、空穴 `5.63143e-6`，这里为
无量纲比例。这证明当前状态的模型回放高度一致，不能替代自洽 I–V 验证。

| 高场饱和参数 | 原生电子 / 空穴 | Vela 当前 | 结论 |
| --- | --- | --- | --- |
| 模型 | Caughey–Thomas | 对应高场饱和核 | 当前配置对应 |
| beta0 | 1.109 / 1.213 | 1.109 / 1.213 | 300 K 一致 |
| vsat0 | 1.07e7 / 8.37e6 cm/s | 1.07e5 / 8.37e4 m/s | 单位换算后一致 |
| betaexp | 0.66 / 0.17 | 当前只验收 300 K | 不据此认领变温一致性 |
| vsatexp | 0.87 / 0.52 | 当前只验收 300 K | 不据此认领变温一致性 |
| 驱动力 | GradQuasiFermi 路径 | `quasi_fermi_gradient`、`transport_cell_vector` | 元素场实现已回放，非所有离散量逐点完全相等 |
| RefDens 插值 | 未启用 | IALMob `reference_density_m3=0` | 一致 |

**OldSlotboom、SRH 与 Auger**

OldSlotboom 主参数：`Ebgn=0.009 eV`、`Nref=1e17 cm^-3`、根号内常数
`0.5`、OldSlotboom 的 `dEg0=0`，当前 Vela 使用对应值。

| SRH 参数 | Sentaurus 电子 / 空穴 | Vela 电子 / 空穴 | 判断 |
| --- | --- | --- | --- |
| tau_min | 0 / 0 s | 0 / 0 s | 相同 |
| tau_max | 1e-5 / 3e-6 s | 1e-5 / 3e-6 s | 相同 |
| Nref | 1e16 / 1e16 cm⁻³ | 1e16 / 1e16，unit_scaling 输入 | 相同 |
| gamma | 1 / 1 | 1 / 1 | 相同 |
| 掺杂依据 | DopingDep | `total_impurity`，Nd+Na | 当前目标对应 |
| 陷阱能级 | Etrap=0 eV | 中间能级 SRH | 当前对应 |
| 温度依赖开关 | 未启用 TempDep；运行日志明确关闭 | `temperature_dependence:true` | 开关不同；300 K 的乘数均为 1 |
| Talpha | S2 为 −1.5 / −1.5，但此开关未激活 | −1.5 / −1.5，参考温度 300 K | 不能推断离开 300 K 后仍等价 |

原生 Auger 系数不是单纯常数。S2 给出的形式为：

`Cn(n,T) = [An + Bn(T/300) + Cn(T/300)^2] × [1 + Hn exp(−n/N0n)]`，
空穴系数用对应的 p 和参数。普通 Auger 不包含负产生项。

| Auger 参数 | 原生电子 / 空穴 | 当前 Vela | 判断 |
| --- | --- | --- | --- |
| A，cm⁶/s | 6.7e-32 / 7.2e-32 | 未以 A/B/C 多项式配置 | 当前固定 300 K |
| B，cm⁶/s | 2.45e-31 / 4.5e-33 | 同上 | 当前固定 300 K |
| C，cm⁶/s | −2.2e-32 / 2.63e-32 | 同上 | 当前固定 300 K |
| 300 K 基础系数 A+B+C | 2.9e-31 / 1.028e-31 cm⁶/s | `auger_cn_m6_per_s=2.9e-31`、`auger_cp_m6_per_s=1.028e-31`，unit_scaling | 基础系数相同，字段后缀不能代替单位规则 |
| H | 3.46667 / 8.25688 | 未配置 `auger_density_dependence`，默认关闭 | **明确缺少密度增强** |
| N0 | 1e18 / 1e18 cm⁻³ | 未配置 | 同上 |
| WithGeneration | 关闭 | `auger_with_generation:false` | 已对齐；旧配置默认 true 不代表当前配置 |

Vela 已有 H/N0 模型及系数对密度的解析导数，见
[配置契约](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/docs/config_schema.md:1115) 与
[默认配置](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/include/vela/physics/RecombinationModel.h:60)。
所以这是当前 SLOT 配置的遗漏，而不是需要从零实现的功能。对应的候选配置为：

```json
"auger_density_dependence": {
  "enabled": true,
  "electron": {"enhancement": 3.46667, "reference_density_m3": 1e18},
  "hole": {"enhancement": 8.25688, "reference_density_m3": 1e18}
}
```

以上 N0 数字仅适用于当前 `unit_scaling` 的 cm⁻³ 输入语义；本次仅记录候选，
没有修改 V1。形式上低密度系数增强倍数分别趋向 4.46667、9.25688，但普通
Auger 在产生侧关闭，因此不能把这些倍数当成漏电倍数，也不能直接认定它
解释了剩余约 2.75% 的电流差异。此前五项预算没有覆盖这项 H/N0 开关。

**雪崩参数、源项与完整 Jacobian**

van Overstraeten 参数按原生 cm 制列出；当前 Vela 内部 SI 系数已作单位换算，
没有用系数缩放拟合现有 I–V。

| 参数 | 电子 | 空穴 | 当前 Vela |
| --- | ---: | ---: | --- |
| a_low，cm⁻¹ | 7.03e5 | 1.582e6 | 相同 |
| a_high，cm⁻¹ | 7.03e5 | 6.71e5 | 相同 |
| b_low，V/cm | 1.231e6 | 2.036e6 | 相同 |
| b_high，V/cm | 1.231e6 | 1.693e6 | 相同 |
| 分段场 E0，V/cm | 4e5 | 4e5 | 相同 |
| 声子能量，eV | 0.063 | 0.063 | 相同 |

系数来源：[原始参数文件](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/build/slot_deck_comparison_20260921/original_reference/pp2_des.par:2766)
和 [Vela 模型](D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next/include/vela/physics/ImpactIonizationModel.h)。

| 雪崩/导数细节 | Sentaurus | Vela 当前 |
| --- | --- | --- |
| 耦合 | Poisson、电子、空穴自洽 | `coupling_mode:self_consistent` |
| 驱动场 | `GradQuasiFermi` | 准费米梯度；零 RefDens、无插值、无额外载流子截断 |
| 用于产生率的电流/密度 | `AvalDensGradQF` | `current_approximation:element_edge_sg_gss_laux`，Bernoulli 中点密度 |
| 源积分 | `ElementVolumeAvalanche` | `source_mapping_mode:element_vertex_box_measure`，当前本地元素分配几何 |
| 迁移率导数 | S1 中 `Derivatives` 被注释；S4 显式开启 | `jacobian_field_derivatives:true` |
| 雪崩导数 | S3/S5 实际运行均报告 With avalanche derivatives | `source_jacobian:local_ad` |
| 偏导范围 | 不能仅凭一行 deck 推断内部全部实现 | 元素三个顶点各 ψ、φn、φp，共九列；包含 n/p、法向场、高场及迁移率反馈 |

在当前等温、固定掺杂配置中，“完整”是相对于求解未知量 ψ/φn/φp 的链式
导数；不表示已经验证温度、量子势或所有统计模型的扩展。当前元素雪崩实现
验收范围为 300 K、Boltzmann、无量子势及上述 canonical 源项路径。
`source_jacobian_mode` 不是当前解析器读取的键，不能拿历史文件中该键
作为激活导数的证据。

已有实际 10 V 状态的雪崩源有限差分审计相对矩阵差为 `5.2416e-10`，
ψ/φn/φp 分块分别约 `5.42e-10 / 4.96e-10 / 8.27e-10`。同状态导数通过
并不代表原生与本地的源项离散相同。原生关闭 AvalDensGradQF 的独立实验
曾使 10 V 电流变化约 `+0.6780%`，仍不能把这个变化量直接等同于本地误差。

**网格、AverageBox 与源项体积**

Si 区域为 12,410 个节点、24,272 个三角形。原生日志记录整个网格有 2 个
non-Delaunay 单元；non-Delaunay 与普通钝角三角形不是同一个分类。
V1 设置 `require_non_obtuse:false`，没有通过严格非钝角资格。

V1 保留 `assembled_delaunay_fallback` 边耦合和 `mixed_voronoi` 节点体积。
元素迁移率使用原生正权重分配并缩放到旧边耦合；原生为零、本地为正的
4 条边保留本地回退。这属于拆分对照，不能把“导入了原生元素几何”写成
“全部方程均已使用原生 AverageBox”。V1 没有启用材料局部 Poisson 电荷体积；
该选项和原生边耦合仅做过独立诊断对照。

此前 2 V 自洽单因素变化：原生元素 Measure `+0.000515%`、原生输运边耦合
`−0.000757%`、材料局部 Poisson 体积 `−0.014959%`、关闭 Auger 产生
`−0.000811%`、Fermi 灵敏度 `−0.001188%`。各项分别从同一基线出发，
不能直接相加，也没有解释剩余约 2.75%。连续性源项材料体积和接触处理
仍需独立预算。

**求解器、扫描与停止条件**

| 项目 | 冻结原始物理参考 S1 | 近期 Sentaurus 对照 S4 | 当前 Vela |
| --- | --- | --- | --- |
| 初始化 | Poisson 最多 1000 次、LineSearchDamping=1e-8，然后耦合 DD | Load 已有 10 V 保存态，再耦合 DD | 导入保存态；部分点需 Poisson 预热及载流子修正，再完整耦合验收 |
| Newton 上限 | 每步 25 | 每步 25 | V1 15；连续控制等配置可为 30，不能混为一个固定值 |
| 全局化 | Extrapolate、Notdamped=100 | 同左 | V1 merit 线搜索；无雪崩连续控制使用 block_filter |
| 误差设置 | RelErrControl，Digits=5，RHSmin=1e-10，ErrRef(n,p)=1e10 | 同左 | reltol=1e-8，abstol=1e-9；另有局部行及全局闭合门槛 |
| 数值精度 | ExtendedPrecision；日志明确 80-bit long double | 同左 | UCRT64 Release，主求解 Real=double；多精度诊断不是整个求解器扩展精度 |
| 线性后端 | ParDiSo | ParDiSo | 无雪崩与 2/5 V 雪崩使用 SparseLU；10 V 雪崩成功控制使用 UMFPACK |
| 线程 | 5 | 5 | 未在这些 JSON 中承诺相同线程数；Release 是构建优化设置 |
| 电压目标 | 外电源 0→200000 V，漏极串 R=1e12 | 从 10 V 保存态出发，外电源目标 20/40/60/100/1000/10000/200000 V | 独立器件内电压点 2/3/5/10 V；连续无雪崩已到 3 V，连续雪崩 10→10.1 V |
| 步长 | Initialstep=1e-6，MaxStep=0.05，Minstep=1e-7 | 每段 Initialstep=1e-4，MaxStep=0.025，Minstep=1e-9 | 以 V 为单位并由具体驱动流程决定；外电阻内层最大试步 1e-4 V |
| 放大/缩小步长 | Increment=1.6，Decrement=2 | 1.4 / 2 | 通用 sweep 为 growth=1.2、shrink=0.5；分段驱动另有控制 |
| 电流停止阈值 | AbsVal=1e-7 | AbsVal=2e-7 | 当前验证配置未设置等价的 1e-7 A/µm BV 验收阈值；高阈值 breakdown 占位项不能代替它 |
| 实际完成范围 | 冻结参考有外电阻击穿支路 | 394 点，完成到外电源 200000 V | 外电阻已知起点负载线一致，非零步进未合格 |

Sentaurus Quasistationary 的 Step 是归一化扫描进度，不是直接以伏特表示。
对 S1 的 0→200000 V 目标，`Initialstep=1e-6` 对应外电源初始增量约 0.2 V；
不能与 Vela 的 `initial_step=1e-4 V` 原样比较。外电源电压也不等于器件内
漏极电压，二维单位下负载线为：

`V外 = V内 + R(Ω·µm) × Id(A/µm)`。

原生 AreaFactor=1，本报告电流按 A/µm 对比；不要再次额外乘除 1 µm 宽度。
S4 最终内电压约 `38.528094 V`、电流约 `1.9996147e-7 A/µm`，并非器件承受
200000 V。按统一 `1e-7 A/µm` 标准插值得到原生 BV 约 `38.5209015 V`；
这个数值尚没有合格的本地 BV 与之对应。

V1 的 `bias_points:[10]` 优先定义实际求解点，残留的 `start/stop=0.105`
不能当作当前扫描范围；V1 虽序列化了 R=1e12，`external_circuit.enabled=false`，
因此 10 V 检查点不是外电阻扫描。V2 才启用外电阻，其目标外电压是
`10.003028145555882`、`10.003128145555882 V`，对应已知内电压 10 V 的起点
和一个小增量。外层负载线残差门槛 `1e-6 V`、内电压容差 `1e-8 V`，
最多 30 次括区间/迭代；内部非零试点未通过，因此不认领支路完成。

Vela 独立验收还要求局部 `eps_row=0.001`，以及全局连续性相对门槛 0.01、
绝对门槛 1e-12、source_floor=1e-11。后两项是内部连续性量的门槛，不应
直接标为端口 A/µm。其归一化方式不同于 Sentaurus Digits/RHSmin，不能
仅按数字大小断言哪一套标准更严格。当前还启用了载流子对角稳定项、
准费米更新限制和最终闭合恢复，需与物理模型差异分开理解。

**输出功能及未启用物理**

| 内容 | Sentaurus | 当前 Vela/对照工具 |
| --- | --- | --- |
| I–V | `.plt` 端口曲线 | `iv.csv` 与外电路/边界评估 CSV |
| 空间状态 | `.tdr`，密度、势、掺杂、场、电流等 | 状态 CSV、网格/逐节点掺杂；支持 VTK，但 V1 `write_vtk:false` |
| 迁移率 | S1 `eMobility hMobility`；S4 改为 `/Element` | 元素回放/审计输出 |
| 求解诊断 | Newton 表、RHS、误差、线性求解日志 | attempts/iterations CSV、局部行审计、全局闭合 JSONL、Jacobian 审计 |
| 温度、隧穿、陷阱等 Plot 名称 | 原 deck 列出，不等于已启用方程或机制 | 当前不作为必须补齐的激活模型 |

原始 Thermode/Thermodynamic 已注释，Solve 没有温度或载流子能量方程；
没有启用量子势、BTBT 或陷阱物理。`Band2Band`、`eDirectTunnel`、
`eTrappedCharge`、`eTemperature` 出现在 Plot 中，不能据此把这些模型
列入本次原始算例的“已启用功能”。原始 par 中有 SurfaceRecombination
参数也不表示表面 SRH 已被激活。

**已有验证范围与对差异的判断**

| 器件内电压 | 无雪崩 Vela Id，A/µm | 相对原生 | 自洽雪崩相对原生 |
| --- | ---: | ---: | ---: |
| 2 V | 9.0118343756e-16 | −2.74887% | −2.729706% |
| 3 V | 1.1889774634e-15 | −3.01762% | 本表无对应合格点 |
| 5 V | 1.6839859459e-15 | −5.13559% | −4.837721% |
| 10 V | 2.9214614951e-15 | −2.14127% | −0.751459% |

这些是同偏压独立检查点，不是完整连续 0→10 V 扫描的证据。连续无雪崩
2→3 V 已完成，最终电流与独立 3 V 点相对差约 `1.86e-10`；连续雪崩
10→10.1 V 已完成。最近代码验收为 Release CTest 877/877，通过记录见 E1；
本次文档核对没有重新运行数值仿真或 CTest。

对后续严格对齐，依据本次证据的建议优先级如下。这是差异分析建议，
不是已执行或已证实能够消除电流差的修正：

1. 单独启用现有 Auger H/N0 模型，保持普通 Auger 开关及原验收门槛，核对
   系数、产生复合量及其 Jacobian。该项是可直接定位的配置遗漏。
2. 单独对齐 Si/SiO2 能带参数并保留 ni 校准对照，检查电势参考；另行验证
   多晶硅掺杂/载流子和 gate/SLOT 接触，避免同时改变多项而无法归因。
3. 继续分离 AverageBox 边耦合、材料局部体积及连续性源积分，随后对照
   AvalDensGradQF 与本地 SG/GSS-Laux 的空间产生率。
4. 将数值精度、停止准则、初始化和外电路连续步进单独验收，最后使用
   统一内电压和 `1e-7 A/µm` 标准比较 BV；不能把失败缩步点当作击穿。

本次最关键的更新是：迁移率“参数未对齐”已不是主要未完成项；当前仍有
确定的配置和边界物理差异，且还没有证据把剩余电流差唯一归因于某一项。
