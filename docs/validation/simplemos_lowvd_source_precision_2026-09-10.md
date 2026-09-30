# SimpleMOS 低 Vd 带源 Poisson 精度定位与隔离验证

日期：2026-09-10。工作分支：`codex/simplemos-sdevice-validation`。

本轮完成实际带源拒绝轨迹定位、独立高精度参考、两种 Poisson 求值候选及双幅度低 Vd 自洽对照。求和精度候选为 2/6 DC 通过；从保存坐标高精度重算物理量的候选为 10/10 通过。再加入精确接触行更新，扩展四个控制工况，共 20/20 DC、8/8 全场响应、16/16 原生端口响应通过。未将不同实现的合格状态拼接成统一 16/16，也未开始完整 PhuMob 曲线、Enormal 或高场饱和。

## 配置与冻结条件

复用 n19/n23 的既有二维网格（分别 1480/1482 节点），Vg=0 V、300 K，Boltzmann、OldSlotboom、plain PhuMob `element_box_phumob`、掺杂相关 SRH。低 Vd 隔离使用 0.05 V，联合控制另覆盖 Vd=1 V。保持已验证联合几何、单元介电系数及原 SRH 体积。电流单位 A/μm。

八个局部节点的等量电子/空穴源、正负幅度 ±0.001/±0.0005 与前轮相同；零源基态来自已冻结的第二次零漂移重启。源独立于状态，不新增 Jacobian 项。该源试验不等于有限 SRH 体积替换。

实际 Windows UCRT64 Release、Eigen SparseLU、四次线性修正，稳定范数差和重启坐标修复仍为诊断叠加实现。保留 200 次迭代、绝对/相对残差 1e-12/1e-7、既有残差平台 1e-9 退出支路、载流子逐行 1e-6、KCL/Id 1e-8、行缩放、步长上限和回溯网格。平台接受仍须通过原逐行资格，不能用绝对小残差覆盖载流子失败。未改变物理常数或接受门槛。

## 实际拒绝原因

四个正负全幅失败轨迹均逐位复现最终状态、退出状态与重复残差，共 56 个快照。对每组基态及实际 trial 0/5/12 共 16 个快照做全部方程行的 Decimal 60/100 位检查，明确在两个载流子方程中保留同一独立源项。精度一致性最大 2.1793e-41，舍入项重建误差最大 3.06888e-16，均按通量/源项规模归一化。

原始步与截断步四组完全相同，独立完整 J 线性缺陷/残差范数最大 4.56536e-17。四个实际全步均已将载流子行比降至约 2.98e-11～1.16e-10，但 Poisson 块使总范数增加而拒绝；不能归因为载流子不下降或更新截断。

n19/n23 正扰动全步最大电势更新约 9.03e-17/8.99e-17 V。最小回溯时自由节点电势更新全部被状态加法舍入吸收，剩下 39 个发生极小变化的电势节点均为 Dirichlet 接触。只提高范数比较精度不足以修复这类候选。

保持生产双精度内核操作数的高精度求值，四组全步仍不下降；从实际保存坐标高精度形成 ψ/n/p 的诊断使两个 n23 全步下降，但 n19 仍不下降。无舍入的数学试探状态四组均下降。后者不是可直接接受的生产状态，仅用于区分求值误差与更新表示误差。

## 两种独立求值候选

- `kernel`：用 binary128 重算 Poisson 边差、乘积与求和，ψ/n/p 仍使用现有双精度值。
- `packed`：额外从实际双精度保存坐标、准费米参考和增量，以 binary128 重算物理电势及 Boltzmann 浓度，再形成 Poisson 残差。

两者保持原双精度状态存储及原 Jacobian；迁移率、连续性通量、端口提取、源项和体积不变。采用 Boost `cpp_bin_float_quad` 的软件 binary128，仅作为有保护的隔离候选，不是全系统高精度求解。其 47,392 个 Poisson 节点检查通过独立 100 位参考，最大按项规模误差约 2.31e-33，小于冻结 1e-25。关闭选择器时，n19/n23 正扰动重放均逐位保持原轨迹。

| NWell | 源幅度 | 求值方式 | DC 资格 | 最大载流子行比 | KCL/Id | 迭代 |
| --- | --- | --- | --- | --- | --- | --- |
| n19 | 0.0 | kernel | 通过 | 4.68864141e-11 | 2.44626e-15 | 13 |
| n19 | 0.0 | packed | 通过 | 4.68882091e-11 | 3.85119e-15 | 4 |
| n19 | 0.001 | kernel | 失败 | 0.00017719478 | 2.8461e-15 | 2 |
| n19 | 0.001 | packed | 通过 | 3.95697549e-11 | 1.03804e-15 | 15 |
| n19 | -0.001 | kernel | 失败 | 0.000177976704 | 1.18084e-15 | 0 |
| n19 | -0.001 | packed | 通过 | 4.68882091e-11 | 1.00838e-15 | 97 |
| n23 | 0.0 | kernel | 通过 | 8.72241181e-11 | 5.75367e-16 | 0 |
| n23 | 0.0 | packed | 通过 | 9.47459304e-11 | 3.6735e-15 | 171 |
| n23 | 0.001 | kernel | 失败 | 6.80300176e-05 | 5.75367e-16 | 0 |
| n23 | 0.001 | packed | 通过 | 8.88300353e-11 | 1.41629e-15 | 200 |
| n23 | -0.001 | kernel | 失败 | 6.26660553e-05 | 1.32777e-16 | 39 |
| n23 | -0.001 | packed | 通过 | 1.00252972e-10 | 1.28351e-15 | 181 |
| n19 | 0.0005 | packed | 通过 | 1.32857525e-11 | 1.75643e-15 | 30 |
| n19 | -0.0005 | packed | 通过 | 2.04145824e-11 | 9.45771e-16 | 8 |
| n23 | 0.0005 | packed | 通过 | 1.77625581e-10 | 9.73698e-16 | 89 |
| n23 | -0.0005 | packed | 通过 | 5.65680142e-11 | 6.19626e-16 | 200 |

所有失败和原接受预算保留。DC 通过仍不自动等于双幅度全场响应、初始化不变性或生产实现资格；后续必须依据完整对照和相应门槛决定是否移植。

Poisson 对 ψ/φn/φp 三个方向的三档步长检查，两个低 Vd 基态共 18/18 通过；最大相对差 1.39453e-9，小于 1e-6。该结论限于此次 Poisson 块，不能代替任意载流子交叉块审计。重新导出完整 J、固定源、输运通量的预检通过；源项没有新增导数。

## 接触约束与四点联合验证

四个实际失败基态各有 246 个严格单位矩阵边界行，其中 49 个线性求解更新带有约 1e-45 的约束舍入偏差。只把这些单位行的步长设为其残差的相反数，即精确求解原边界方程，自由行步长逐位不变；约束行线性缺陷变为零。没有修改边界偏置、自由状态、接受门槛或物理模型。

三个隔离对照全部通过：n19 负全幅从 97 次迭代降为 10 次；n23 零源从 171 次降为 3 次，正全幅从 200 次降为 3 次。两项 n23 电流逐位相同，n19 相对变化约 1.4e-16。原长迭代轨迹仍保留，不能以此覆盖初始失败。

随后固定该数值组合，在 n19/n23 × Vd=0.05/1 V 上重算完整 J、源插入、正负双幅度及零源。继续复用已经运行并验证的原生同源实验，未新增原生仿真。

| 工况 | DC 通过 | 最大载流子行比 | 最大 KCL/Id | 迭代范围 |
| --- | --- | --- | --- | --- |
| m65_n19_vd_0p050000_endpoint | 5/5 | 4.68882e-11 | 3.85119e-15 | 4–14 |
| m65_n19_vd_1p000000_endpoint | 5/5 | 3.68894e-10 | 4.67663e-15 | 5–20 |
| m65_n23_vd_0p050000_endpoint | 5/5 | 1.77626e-10 | 3.6735e-15 | 3–38 |
| m65_n23_vd_1p000000_endpoint | 5/5 | 3.65461e-10 | 8.96903e-15 | 2–7 |

联合全场按原预测相对误差 1e-3、双幅度 1e-3、偶/奇比 1e-2、信号/漂移 ≥100 判断；端口按同扰动原生相对误差 1e-3 判断。原生漏端/衬底响应最大相对差 9.9080264e-06。注意：通用源分析器的 `summary.json` 保留旧字段 `native_response_qualified=false`；本轮独立原生端口结果应读取 `native_ports.csv` 与 `combined_summary.json` 的明确端口计数，不能把通用字段解释为这 16 项端口检查失败。

当前完成数值候选与源响应资格，尚未把这些诊断叠加实现移入正式求解器；下一阶段是生产移植、统一实现下的 16 个控制点与双初始化复验，通过后再扩展四条 PhuMob 完整曲线及恢复 Enormal/HFS。此次合格源响应不等于有限 SRH 体积替换已获验证。

## 证据与范围

- [输出身份](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/identity.csv)、[分块实际范数差](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/precision/actual_deltas.csv)、[高精度参考](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/precision/merit.csv)。
- [独立 C++ 内核验证](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/poisson_candidate/kernel_check/summary.json)、[全幅及零源 DC](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/poisson_candidate/dc.csv)、[半幅 DC](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/poisson_candidate/half_amplitude/dc.csv)。
- [Poisson 方向导数](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/poisson_candidate/response/jvp/checks.csv)、[低 Vd 全场响应](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/poisson_candidate/response/response.csv)、[接触行隔离](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/poisson_candidate/contact_step/dc.csv)。
- [四工况联合 DC](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/poisson_candidate/contact_step/full_source/dc.csv)、[联合全场响应](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/poisson_candidate/contact_step/full_source/response.csv)、[原生端口响应](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/poisson_candidate/contact_step/full_source/native_ports.csv)。
- [阶段汇总](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/stage_summary.json)。

本轮新增诊断脚本、诊断 C++ 头文件和报告；未修改生产 C++，未新增 Sentaurus 仿真，未提交或推送。联合零源 n23、Vg=0 的电流误差为 -0.01819672%（Vd=0.05 V）和 -0.80813524%（Vd=1 V）；后者 Vela 电流为 3.68309865058046e-16 A/μm，原生为 3.71310556502558e-16 A/μm。数值修复没有消除既有深关断电流差。[当前零源电流对比](../../reference_tcad/simplemos_sentaurus2022/phumob_source_precision_20260910/zero_current_comparison.csv)。完整 PhuMob 16 点、双初始化及四条曲线仍待统一实现后验证。
