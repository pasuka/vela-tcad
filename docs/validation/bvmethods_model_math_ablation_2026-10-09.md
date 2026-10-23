# BVMethods 物理模型与数学字段隔离对照（2026-10-09）

## 结果与范围

已完成 11 组 Sentaurus 原生从零纯 QS 单因素对照，以及本机、T470p、C560 的 Vela
低偏置交叉检查。最有用的新线索是高场迁移率的 Newton 线性化：Vela 保留全部物理
residual，仅设置 `mobility.jacobian_field_derivatives=false`，三平台均接受 10 µV 首点；
开启导数的匹配基准均被原有载流子逐行门限拒绝。这是后续数值候选，尚未检验指定外目标负载或完整 QS。

Sentaurus 关闭高场饱和或带隙收窄后均达到原电流停止值，独立整条轨迹检查通过；
关闭 BTBT 虽达到停止值，但一个接受点未通过负载门限。关闭预测器没有消除失败，
关闭雪崩/迁移率导数反而更早失败。由此不能将原生高压失败归因于唯一一个模型，
也不能仅凭简化模型通过证明 Vela 与 Sentaurus 的全部差异已经定位。

本次按用户请求执行诊断消融，覆盖此前计划中的“物理模型不变”诊断范围限制；
原始模型的验收契约、固定门限、生产默认值和已资格化 v4 的范围保持原定义。
未修改 C++、canonical `.cmd`、网格、材料参数或参考基线，未提交、推送或换包。
机器可读结果见[JSON](bvmethods_model_math_ablation_2026-10-09.json)。

## 原生共同条件与单因素设置

- Sentaurus VM：RHEL 7.9，SDevice T-2022.03-SP2，版本 0.7745337；每进程单 assembly/solver 线程。
  原始矩阵两进程并行，后加迁移率导数案例短时成为第三进程；本轮不比较耗时性能。
- 原始 `bvmethods_nmos_msh.tdr`：2600 节点、5020 三角形；300 K。
  每组保持相同网格、`models.par`、四个接触、Plot、Solve、零偏初始化和 QS 扫描。
- 漏极串联电阻 `1e7 Ω·µm`，外目标 `100000 V`，原停止 `|Id|≥1.443e-3 A/µm`。
  曲线提取在首个向上越过 `|Id|=1e-4 A/µm` 的两点间，对 **InnerVoltage** 线性插值。
- 原生模型基准含 Fermi、OldSlotboom、掺杂/Enormal/GradQuasiFermi 高场迁移率、SRH(DopingDep)、
  Avalanche(Eparallel) 和 Band2Band(E2)。Math 基准含 Extrapolate、AvalDerivatives。
- Physics 对照仅删除或替换表中一个设置；Math 对照保留完整 Physics。
  未删掉整个 Recombination，也未同时更改材料和数学容差。

| 案例 | 相对原脚本的唯一变动 | 100 µA/µm 内电压 (V) | 相对基准 (mV) | 扫描及独立整轨迹检查 |
|---|---|---:|---:|---|
| baseline | 无 | 6.318404 | 0 | MinStep；已有接受点检查通过，未到停止值 |
| no_btbt | 删除 `Band2Band(E2)` | 6.514858 | +196.454 | 到停止值；一处负载检查失败 |
| no_enormal | 删除 `Enormal` | 6.317734 | −0.670 | MinStep；负载检查失败 |
| no_hfs | 删除 `HighFieldsaturation(GradQuasiFermi)` | 6.282514 | −35.890 | 到停止值；整轨迹检查通过 |
| no_bgn | `OldSlotboom` → `NoBandGapNarrowing` | 6.253173 | −65.230 | 到停止值；整轨迹检查通过 |
| boltzmann | 删除 `Fermi` | 6.318960 | +0.556 | MinStep；已有接受点检查通过 |
| constant_srh | `SRH(DopingDep)` → `SRH` | 6.316321 | −2.082 | MinStep；负载检查失败 |
| no_avalanche | 删除 `Avalanche(Eparallel)` | 9.915732 | +3597.329 | 到停止值；整轨迹检查通过；此值不是雪崩 BV |
| no_extrapolate | 删除 `Extrapolate` | 6.320317 | +1.913 | MinStep；已有接受点检查通过 |
| no_aval_derivatives | `AvalDerivatives` → `-AvalDerivatives` | 未到提取电流 | — | MinStep；已有接受点检查通过 |
| no_mobility_derivatives | Math 添加 `-Derivatives` | 未到提取电流 | — | MinStep；负载检查失败 |

全部案例已结束，没有超时或缺失输出。Math 开关实际生效：原生日志分别报告
`Without avalanche derivatives`，以及 `Without exact derivatives / With avalanche derivatives`。
这里的 `-Derivatives` 是 SDevice 的数学开关，覆盖范围与 Vela 的场导数开关不同。

独立原生门限沿用既有校验：`|Vouter−Vinner−R·Id|≤1e-6 V`，
`|ΣIterminal|≤max(1e-12,1e-5·|Id|) A/µm`，位移电流和电流分量闭合、
有限性、零偏起点、固定其他接触与外电压单调性均逐点检查；原停止值/成功退出另行检查。
`Good Bye` 不能掩盖 MinStep 失败，达到停止值也不能代替物理闭合检查。

关键失败边界：

| 案例 | 最后外电压 (V) | 最后内电压 (V) | 最后电流 (A/µm) | 最大负载残差 (V) |
|---|---:|---:|---:|---:|
| baseline | 14099.811222 | 7.160255 | 1.4092651e-3 | 1.2396e-9 |
| no_extrapolate | 14099.811222 | 7.160255 | 1.4092651e-3 | 1.8317e-9 |
| no_aval_derivatives | 240.364786 | 5.758013 | 2.3460677e-5 | 1.0673e-9 |
| no_mobility_derivatives | 3.410818 | 3.072959 | 3.3785924e-8 | 1.6640e-6 |

`no_btbt` 的失败接受点为第 54 个零基行：外电压 2142.78514598 V、内电压
6.58382181 V、电流 2.13620054488e-4 A/µm；负载残差 7.79295e-4 V。
其 100 µA/µm 插值的两个括点独立检查均通过，故该阈值处的曲线敏感性仍可使用，
但整条扫描不能宣布通过。Enormal/SRH 两组最大负载残差为 7.34472e-5/9.23519e-4 V。

小于数 mV 的提取变化受原始自适应 QS 采样与线性插值影响；JSON 保留两个括点及检查，
本轮未做网格/步长误差收敛研究，因此不将这些数值当成高精度模型效应。
关闭 Avalanche 仍保留 BTBT(E2)，9.915732 V 只是相同电流阈值，不能称作雪崩击穿电压。

## Vela 匹配首点检查

各案例先重新获得零偏态，再直接求内电压 `10 µV`；此筛查移除外电阻控制器，
没有预偏置、瞬态获取或修改载流子门限。三平台均使用当前 Release runner、
SparseLU/COLAMD 和单线程；T470p 为 Windows/UCRT64，C560 为 Ubuntu 26.04.1 LTS（现场核对 `/etc/os-release`）。
网格和材料 JSON 在两节点的序列化字节不同，但解析后内容逐项相等，掺杂 CSV 哈希相同。

沿用已记录 Q2/Q3 候选：contact_basin、L2 行列平衡、10 次线性精化、精确 Dirichlet 更新、
transport/avalanche 内部未知量列差分参数 1e-9/1e-11、componentwise_poisson_recorrect。
Newton 20 次，无首点回退；`eps_row=0.001`，全局连续性闭合门限仍为 0.01。
这些是诊断候选配置，尚未作为发行默认。

| 平台/实际控制 | 高场物理模型 | 场迁移率 Jacobian 导数 | 10 µV 结果 | 最大载流子行比 |
|---|---|---|---|---:|
| 本机 baseline | 开 | 开 | carrier_row_convergence 拒绝 | 1.98566027 |
| 本机 no_btbt | 开，BTBT 关 | 开 | carrier_row_convergence 拒绝 | 1.98566069 |
| 本机 no_hfs | 关，Enormal/Masetti 保留 | 开关已无高场作用 | 接受 | 7.30685e-5 |
| 本机 frozen_hfs_jacobian | 开 | 关 | 接受 | 9.98823e-4 |
| T470p baseline | 开 | 开 | carrier_row_convergence 拒绝 | 1.98566145 |
| T470p frozen_hfs_jacobian | 开 | 关 | 接受 | 9.98823e-4 |
| C560 baseline | 开 | 开 | carrier_row_convergence 拒绝 | 1.58484659 |
| C560 frozen_hfs_jacobian | 开 | 关 | 接受 | 8.37067e-4 |

本表的全局电子/空穴连续性均在原门限内。独立终端 KCL 的证据需补充限定：
拒绝 sweep CSV 的电流可能为占位零，不能据此宣称独立通过；后续
[实际保存态四端审计](bvmethods_mobility_loadline_audit_2026-10-09.md)已重新提取两节点
完整导数拒绝态及冻结导数接受态的实际电流，确认 KCL 在原门限内。
拒绝态通过这些检查仍不代表接受。冻结场导数的 Windows 结果接近载流子行门限，
不能从一次首点通过推断整段扫描稳定，也不能证明完整 Jacobian 公式必然错误。
本轮优先线索是场迁移率线性化、弱行与更新限幅之间的配合。

**证据标签更正：** 首次远程 worker 将目录名写反：原目录 `frozen_hfs_jacobian`
实际写入 `jacobian_field_derivatives=true`，原目录 `baseline` 实际为 `false`。
计算配置与输出未被改写，原归档完整保留；本报告及 JSON 根据实际 `config.json`
进行分类。已更正后续 worker，并新增“依据实际 selector 分类”的契约测试。
该标签错误曾导致进度更新中的反向判断；以上表格为核对后的最终结果。

## 后续定位与仍未验证项

1. 优先对开启/冻结场导数的实际 Newton 态做同方向输运导数和线性缺陷对照，
   分开检查电子/空穴低密度弱行、GradQuasiFermi 分支、差分尺度与更新限幅。
   不直接把关闭导数设为全局默认。
2. 将保持全部物理 residual 的数值候选接回外电阻首点，独立检查指定外目标负载、
   载流子行、全局连续性、KCL。通过后才能进入 Q4 完整从零纯 QS 轨迹。
3. 原生高压 MinStep 单独调查：当前关闭预测器没有修复，关闭导数更早失败；
   高场/BGN/BTBT 单因素敏感性可用于之后匹配场量，但不能据此替换原模型。
4. 尚未做完整消融矩阵的 Vela BV 曲线或空间场量匹配、历史 2018/当前 2022 版本归因、
   细化电流采样、网格收敛或两平台全仓回归。本轮结果不能直接解释全部历史 BV 偏差。

## 工具、复核与归档

新增工具[脚本](../../scripts/bvmethods_native_model_math_ablation.py)会拒绝源 token 漂移、
校验输入身份、保留失败日志、记录阈值括点并执行独立轨迹检查。
新增[契约测试](../../tests/regression/test_bvmethods_native_model_math_ablation.py) 5 项通过，
既有零 QS 契约 9 项通过，共 14 项；Release runner 目标构建确认无需重编译。
未修改求解器实现，本轮未重跑完整 CTest。三个节点的运行器身份、远程 24 项源文件及
下载归档成员均核验；所有本轮计算已结束。

```powershell
$env:Path = "D:/msys64/ucrt64/bin;D:/msys64/usr/bin;$env:Path"
python tests/regression/test_bvmethods_native_model_math_ablation.py
python tests/regression/test_wp14_zero_qs.py
python scripts/bvmethods_native_model_math_ablation.py summarize `
  --root build-release/bvmethods_native_ablation_20261009 --plot
```

本地完整证据位于 `build-release/bvmethods_native_ablation_20261009/`。
`results.json`、各 `trajectory.csv`、Physics 总览/近阈值放大及 Math PNG/SVG 为生成证据，
不加入 Git。图中展示原生接受点；部分点的独立审计失败状态见本报告和 JSON。
VM 持久归档为 `<REMOTE_HOME>/vela_reference_artifacts/bvmethods_model_math_ablation/v1/`，
包含 11 组原生输入/输出、原始两节点和本机 Vela 对照包、执行计划及 helper/test 快照。
归档哈希和逐文件核验结果见 JSON 的 `transfer_receipts.sealed_archive_receipt.json`。

同目录另存 `neutral_inputs.tar.gz`：7 个精确中性输入及清单，上传后校验 SHA-256；
对应哈希见 JSON 的 `transfer_receipts.neutral_inputs_receipt.json`。
