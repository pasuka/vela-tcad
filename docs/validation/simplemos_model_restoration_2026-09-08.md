# SimpleMOS 原始物理差异与逐步恢复验证

本轮完成原始指令差异审计、OldSlotboom 单项开关对照，以及独立的本征浓度约定对齐对照。主要发现是：恢复 OldSlotboom 时原生基础 ni 同时变化，单改 Vela 的 BGN 名称会造成约 7%–33% 的电流差；对齐 ni 后 16 个状态通过原有自洽、守恒和双初始化门槛，n23 在 Vg=0.8/1 V 的差异降至约 0.0050%–0.0065%。BGN 分支的电势方向 Jv 仍未通过，所以尚未进入 PhuMob 恢复或完整模型放行。

## 原始脚本与当前模型

原始是 T-2022.03-SP2 Applications Library 的 `GettingStarted/swb/SimpleMOS/sdevice_des.cmd`，不是后期关闭 BGN 的 PhuMob/Lombardi/HFS 诊断版本。来源和原始 SHA256 见 [M8 原始物理合同](../../reference_tcad/simplemos_sentaurus2022/simplemos_m8_original_physics_contract_v1.json)；本地按 n23 展开的[原始指令](../../build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics/sentaurus_bundle/n23/n23_vd_0p05_des.cmd)可直接核对。

| 项目 | 原始 Sentaurus | 当前简化对照 | 后续恢复内容 |
|---|---|---|---|
| 有效本征浓度 / BGN | `EffectiveIntrinsicDensity(OldSlotboom)` | 显式 `NoBandGapNarrowing`；Vela `model=none` | 单独开启 OldSlotboom，检查 ΔEg、有效 ni、Boltzmann 状态关系与电流 |
| 基础迁移率 | `PhuMob` | 原生 `Mobility(DopingDependence)`；Vela `masetti` | 恢复 PhuMob；两者都是非恒定迁移率，Masetti 依赖掺杂，PhuMob 还需核对载流子依赖 |
| 法向场 / 表面迁移率 | `Enormal` | 关闭 | 核对实际 Lombardi 公式、界面法向场及单元权重后恢复 |
| 高场饱和 | `HighFieldSaturation` | 关闭 | 恢复高场模型，核对准费米梯度驱动力与全邻接导数 |
| 复合 | `SRH(DopingDependence)` | 保留相同掺杂相关 SRH | 不属于被关闭的项目；SRH 体积保持当前策略 |
| 统计 | Boltzmann；原始无显式 Fermi | Boltzmann | 无需“恢复 Fermi” |
| Auger、雪崩、量子修正 | 原始未启用 | 未启用 | 不列为本次简化导致的差异 |

原始 SRH 的电子/空穴 τmax 分别为 1e-5/3e-6 s，τmin=0、参考掺杂 1e16 cm⁻³、指数 1；当前仍保留。温度 300 K，n19/n23 对应 NWell=1e17/2e17 cm⁻³；这些条件相同。已有 Si、SiO₂、氮化物区域及原生网格/掺杂未因简化而更换。

此外有三类差异，不能与“关闭模型”混在一起：

1. **材料与常数约定**：当前 Vela 使用无 BGN 对照匹配的 Si 基础 ni=1.0750038488844236e10 cm⁻³，材料拷贝来自 M65；开启 BGN 时本轮先保持该值，独立核对原生有效 ni。q、kb、ε₀ 仍用当前 Vela SI 数值，未按手册舍入值替换。
2. **离散与求解设置**：Vela 当前显式启用 `element_box` 输运/迁移率平均、Delaunay box 转移、逐单元介电系数和 Poisson 的独立 signed Si 电荷体积。SRH 体积不变。原始脚本没有指定这些内部离散细节，不能据此声称内部实现完全相同。当前强化的原生 Math 和 Vela 严格逐行闭合继续保留。
3. **验证范围与路径**：原始 M8 是 n17–n24 × Vd=0.05/1 V，Vg=0–2.5 V、步长 0.05 V；当前完整简化扫描只覆盖 n19/n23 × 两个 Vd，Vg=0–1 V、步长 0.02 V。原始从平衡态扫偏压；本阶段从已保存的目标偏置独立重闭合，尚不等于从零启动的完整扫描回归。

## 第一阶段的冻结设计

覆盖 n19/n23 × Vd=0.05/1 V × Vg=0/0.2/0.8/1 V，共 16 个偏置点，每点 BGN 开/关两组：32 个原生目标状态、64 个 Vela 初次 DC 尝试。Vela 分别从该点的合格无 BGN Vela 保存态、以及对应原生电势/准费米势重建的相容状态启动。失败最多在同偏置、相同设置下重载一次；不把失败态传播到其他目标。

原生使用显式 `EffectiveIntrinsicDensity(BandGapNarrowing(OldSlotboom))`；仅删除该指令不能构成无 BGN 对照。当前两组均为 Masetti 加掺杂相关 SRH，其他配置保持一致。

原生 Math 为 `ExtendedPrecision(128) Method=Super RelErrControl Digits=12 ErrRef(Electron/Hole)=1e-2 RhsMin=1e-20 Iterations=40 ExitOnFailure CNormPrint`。Vela 使用当前 Release 生产程序、Eigen SparseLU/COLAMD、4 次高精度线性修正；最大 200 次 Newton 迭代。网格坐标 μm，状态浓度 m⁻³，电流 A/μm，等效宽度 1 μm。

独立检查全部 1814 个自由 Si 载流子行，逐行相对比值 ≤1e-6，KCL/Id ≤1e-8；原有全局源闭合、端口提取、双初始化电势 ≤1e-6 V、密度 ≤1e-4、电流 ≤1e-6 的门槛不变。不新增电流对参考的绝对验收阈值。另在 Vg=0.8 V 对两组模型进行三分量、三幅度 Jv 差分诊断，去掉探针内置 `max(1,FD)` 归一化地板后报告实际分块误差；弱交叉源块不能仅靠该差分证明完整性。

- [原生合同](../../reference_tcad/simplemos_sentaurus2022/bgn_restore_20260908/native_contract.json)
- [Vela 合同](../../reference_tcad/simplemos_sentaurus2022/bgn_restore_20260908/vela_contract.json)
- [原生运行脚本](../../scripts/simplemos_bgn_restore_native_20260908.py)
- [Vela 运行及资格脚本](../../scripts/simplemos_bgn_restore_vela_20260908.py)

## 实际执行结果

所有计算已结束，结果已取回。原生 T-2022.03-SP2 的 8 份指令、32 个目标状态全部正常退出并满足原生电压/KCL 门槛；回传网格、初态、指令均与冻结输入逐文件核对。原生场身份检查覆盖 32×942 个 Si 节点，原生 Boltzmann 关系残差为浮点舍入量级。

### R1：只切换 BGN 名称

原生日志明确选择 `Bennett/Wilson without bandgap narrowing` 或 `OldSlotboom with bandgap narrowing (no Fermi)`。OldSlotboom 的总杂质浓度公式与导出 ΔEg 最大差 2.77556e-17 eV；问题出在基础本征浓度约定：

| 参数 | 显式无 BGN | OldSlotboom 原生导出反推 | R1 Vela |
|---|---:|---:|---:|
| 基础 ni，cm⁻³ | 1.075003848884424e10 | 1.463880541255919e10 | 两组均保持 1.075003848884424e10 |
| 相对无 BGN 基础 ni 的变化 | 0 | +36.17444652% | 0 |

反推采用 `ni_base = ni_eff × exp(−ΔEg/(2VT_native))`；VT 仅从导出电势、准费米势和密度关系检查，约 0.0258519952664849 V，没有写回 Vela。16 个 OldSlotboom 状态内、节点间、偏置间的基础 ni 一致。因此上述 36.17% 是模型约定切换，不能用 q/kb 小数舍入差解释；Vela 既有原始材料 ni=1.463891495876762e10 与该原生 OldSlotboom 值的约 7.5e-6 相对差，才属于另一项很小的常数约定差。

R1 共 66 次 Vela DC 尝试，62 次合格。16 个 BGN 开启目标均通过双初始化，但基础 ni 语义不一致，不能计为同模型对比通过。16 个无 BGN 目标中 14 个通过本轮重载双路径检查；n23、Vg=0 V、两个 Vd 的 Vela 保存态路径各两次未通过，原生初始化路径通过。两点首次失败的最差行比分别为 1.32577e-4、6.01683e-4；一次重载后仍为 6.51938e-5、2.18402e-5。失败是本次额外重载试验结果，不抹去此前完整曲线的合格保存态。

使用错误基础 ni 的 BGN 组相对原生电流差为 −33.2174%～−6.7596%；这是有意保持其他参数不变的诊断结果，不是可放行的物理模型比较。详见 [语义账本](../../reference_tcad/simplemos_sentaurus2022/bgn_restore_20260908/semantics.csv)、[R1 电流与双初始化结果](../../reference_tcad/simplemos_sentaurus2022/bgn_restore_20260908/comparison.csv)。

### R1b：独立对齐 OldSlotboom 基础 ni

随后冻结独立配置，仅将材料副本的 Si.ni 改为 1.4638805412559193e10 cm⁻³。选择依据是 16 个原生场反推值的一致性，不使用 Id 拟合；q、kb、ε₀、BGN 公式、Masetti、SRH、几何和求解设置均不变。没有改写 M65 材料、原始失败结果或全局默认。

32 次初次 DC、1 次同偏置重载，共 33 次尝试，最终 16/16 个目标通过全部状态、逐行、守恒、端口及双初始化门槛。唯一首次失败是 n19、Vd=0.05 V、Vg=0 V 的原生初始化：内部收敛且逐行通过，但独立 KCL/Id=1.55350e-8 超过 1e-8；重载一次后通过。两个接近门槛的合格态是 n19、Vd=1 V、Vg=0.8/1 V 的 Vela 初始化，最大行比为 9.24546e-7/9.87986e-7，不能称为宽裕裕量。

误差定义为 `100 × (Id_Vela/Id_Sentaurus − 1)`，如下表。全部是本轮 OldSlotboom + Masetti + SRH、基础 ni 对齐后的合格点：

| 器件 | Vg，V | Vd=0.05 V 电流差，% | Vd=1 V 电流差，% |
|---|---:|---:|---:|
| n19 | 0 | −0.03658419 | +0.00455261 |
| n19 | 0.2 | +0.00376020 | +0.00531155 |
| n19 | 0.8 | +0.00328417 | +0.00283074 |
| n19 | 1 | +0.00091411 | +0.00117916 |
| n23 | 0 | −0.02132499 | −0.88268614 |
| n23 | 0.2 | +0.00352636 | −0.02084582 |
| n23 | 0.8 | +0.00569486 | +0.00652913 |
| n23 | 1 | +0.00504174 | +0.00504624 |

双初始化最大 ψ、电子准费米势、空穴准费米势差分别为 2.2e-15、3.820694e-11、5.0838023e-10 V；密度最大相对差 1.96650e-8，Id 最大相对差 4.07061e-9。本阶段证明本征浓度约定解释了新增 BGN 大部分电流差，没有解释 n23 深关断的剩余差异，也没有新增 Id 验收阈值。

来源：[R1b 合同](../../reference_tcad/simplemos_sentaurus2022/bgn_restore_20260908/intrinsic_control/contract.json)、[R1b 比较表](../../reference_tcad/simplemos_sentaurus2022/bgn_restore_20260908/intrinsic_control/comparison.csv)、[全部尝试](../../reference_tcad/simplemos_sentaurus2022/bgn_restore_20260908/intrinsic_control/attempts.csv)、[R1b 脚本](../../scripts/simplemos_bgn_intrinsic_control_20260908.py)。

### Jacobian 资格仍有缺口

R1 在 8 个 Vg=0.8 V 状态进行 72 个方向、216 个分块检查。无 BGN 的两个较小幅度共 56 个非交叉分块全部满足 1e-4 门槛，最大真实相对差 2.61881e-7。BGN 组对应的 56 项有 12 项不通过。随后在 R1b 合格状态上用四个幅度 1e-4、3e-5、1e-5、3e-6 V 独立重复，48 个方向、144 个分块检查中，较小三个幅度的 84 个非交叉分块有 18 项不通过。

R1b 在 1e-5 V 幅度的关键结果如下。数值是分块误差范数除以 `max(||Jv||,||FD||)`，不是除以带 1 地板的量：

| 器件 | Vd，V | ψ→电子行，相对差 | ψ→空穴行，相对差 |
|---|---:|---:|---:|
| n19 | 0.05 | 9.74343e-6 | 约 1 |
| n19 | 1 | 6.97116e-7 | 约 1 |
| n23 | 0.05 | 1.14180e-3 | 约 1 |
| n23 | 1 | 2.15389e-4 | 约 1 |

空穴块解析变化范数约 1.93e-11～1.29e-10，而残差差分约 1.76e-16～1.82e-15；减小步长后误差保持平台。这里应同时看绝对量和相对量，不能将“约 100%”解释为电流也有 100% 误差。电子/空穴交叉源块虽有小差分误差，仍未替代独立高精度弱块完整性审计。

源码中稳定 SG ψ 偏导替换目前限定为 `!bgnEnabled_ && niI == niJ`，BGN 会走可变 ni 分支。平台误差和该分支条件提供了复核方向，但本轮没有通过独立高精度计算判定是公式遗漏还是求值消减，也没有修改 Jacobian。**自洽与电流比较通过，不等于 BGN Jacobian 完整性已通过。**

来源：[R1 Jv 账本](../../reference_tcad/simplemos_sentaurus2022/bgn_restore_20260908/jvp.csv)、[R1b Jv 合同](../../reference_tcad/simplemos_sentaurus2022/bgn_restore_20260908/intrinsic_control/jvp/contract.json)、[R1b Jv 分块结果](../../reference_tcad/simplemos_sentaurus2022/bgn_restore_20260908/intrinsic_control/jvp/blocks.csv)、[复核脚本](../../scripts/audit_simplemos_bgn_matched_jvp_20260908.py)。

## 分阶段恢复顺序

| 阶段 | 独立工作 | 进入下一阶段所需证据 |
|---|---|---|
| R1 | 已完成 BGN 开关与独立 ni 对齐；下一步复核可变 ni SG 的稳定 ψ 偏导 | 保持此次全部失败记录，独立高精度公式/分块 Jv 通过后，再复查重载和双初始化 |
| R2 | 在相同几何下恢复 PhuMob 基础迁移率 | 原生单元到边形成方式、电子/空穴载流子依赖与交叉 Jacobian 校准通过 |
| R3 | 加入 Enormal/Lombardi | 法向场与表面距离/单元权重明确，小幅度正负响应及全部邻接导数通过 |
| R4 | 加入高场饱和 | 准费米梯度驱动力、场重建及交叉块导数通过，覆盖高 Vd 反例 |
| R5 | 完整原始物理与工况回归 | 上述阶段合格后再扩展 n17–n24 和原始 0–2.5 V 曲线，同时保留无 BGN / Masetti 控制 |

当前 `element_box` 迁移率生产分支明确只接受 constant/Masetti，见 [MobilityModel.cpp](../../src/physics/MobilityModel.cpp) 和 [配置说明](../config_schema.md)。因此 R2–R4 不能仅把字符串替换成 `phumob_field_lombardi`：那会触发拒绝；切回旧平均方式又会同时改变离散，失去单项恢复的可解释性。当前未修改该保护、生产默认或接受门槛。

## 本轮验证与保存

累计完成 32 个原生 DC 目标和 99 次 Vela DC 尝试，保留其中 5 次失败记录；另完成 120 个 Jv 方向、360 个分块检查。所有后台任务均已结束。新增实验控制的 5 项单元检查全部通过，5 个 Python 源/测试文件语法检查通过；32 个原生状态的坐标映射误差均为 0，16 个无 BGN 原生重算点相对上一轮原生 Id 的漂移为 0，R1b 材料副本逐字段核对仅 Si.ni 改变。

本轮未修改 C++、构建配置、原有材料或接受门槛，未重跑完整 CTest。原有完整测试中的历史证据身份失败记录仍保留，不能用本轮诊断检查代替全套回归。新增脚本、报告保留在工作区，生成网格、状态、TDR、CSV 与冻结证据保留本地，本轮未提交、推送或合并分支。
