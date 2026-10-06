# SDevice WP10 原契约收口，2026-10-03

本次基于 `0842e000`，补齐 PN2D BV 的独立偏置场量资格。原契约要求的偏置为
**−13.2 V**；它是场量比较点，**不是击穿电压**。用户明确选择按原契约收口，保留
额外 −20 V 空穴密度差异为已知限制。本报告不建立全反偏范围的空间场量资格。

主验收保持 PN2D 0V、IV 0.2…0.3 V 曲线及 0.2 V 场量、BV −20…−0.05 V 曲线及
−13.2 V 场量、参考六检查、C6 符号/逐段 Time 和门控案例。此前结果见
[首次 WP10 快照](sdevice_acceptance_2026-10-03.md)与[IV/BV 修复快照](sdevice_iv_bv_repair_2026-10-03.md)。

## 独立门限与资格范围

[契约](../../configs/sdevice_acceptance_contract.json)新增 `biased_spatial_qualification.pn2d_bv`，
资格 ID 为 `wp10-pn2d-bv-v2-spatial-minus13p2V-300K`。静态非雪崩门限不自动继承。

| 量 | 新独立门限 | 工程含义 |
| --- | --- | --- |
| 电势 RMS / p95 | 0.0005 / 0.001 V | RMS 0.5 mV，约为 300 K 热电压的 1.93% |
| 电子、空穴 log10 密度 RMS / p95 | 0.005 / 0.01 dex | 倍率约 1.0116 / 1.0233 |
| 矢量电场 global / junction relative L2 | 各 0.02 | 全局及结区各 2%，均必须通过 |
| 活跃密度支撑 | 原生密度 ≥1e5 cm^-3，每种至少 100 节点 | 防止空支撑或极少节点被计为通过 |
| 坐标配对 | 1e−10 µm，节点 ID 双射 | 同网格直接配对，无场量插值 |

这些是针对该案例的工程回归预算，不是商业级击穿预测精度声明。门限在额外偏置场量比较前
固定于 `2026-10-03T11:52:34.735599Z`；初始契约 SHA-256 为
`c603ffc221a619c1f3a762fd72b6a1e3b0a02e48eb7caf90e735c6f6cf8f04df`。
冻结记录在 `build/wp10_closeout/qualification_frozen.json`。此前已看过 −13.2 V，
因此不宣称这是完全盲测；其余初始六个点用于追加检验。

初始七点审查在 −20 V 失败，完整失败报告保留。随后用户明确选择保持原 WP10 要求的
−13.2 V 资格范围，其他点作为诊断。**没有改变已冻结的任何数值门限**，也没有删除
−20 V 失败。接着补充 −19.5/−19.6/−19.7 V 诊断与已有击穿膝点判据。

资格同时绑定原生 v2 cmd、参数、网格、BV sidecar 和参考清单 SHA-256，以及扫描段
`solver + mesh_geometry + contacts` 的规范 JSON 指纹
`90b87da7fa7f57750050a8b3beb914e805fce52d112f24e2ab0350616e6287e9`。
适用的是 300 K、1943 节点/3680 三角形、2×0.5 µm Si PN 结、同节点掺杂、
非钝角网格、显式完整 `element_edge_sg_gss_laux` 配置和自洽 VanOverstraeten。
输入、温度、物理模型、接触或数值组合漂移后不能沿用资格；缺资格仍返回
`unqualified_threshold`。这不会改变 C++ 全局默认或其他设备的资格。
`.gitattributes` 保留契约已记录的字节，并将 canonical BV sidecar 固定为 LF，
避免新 checkout 的换行转换改变已绑定的字节指纹。

## 多偏置场量证据

单位为 V、cm^-3、V/cm；密度误差为活跃节点上的 log10 比值；电场以几何节点面积
加权比较矢量，全局与结区分别检查。使用真实 HDF5/TDR importer 与已有原生 v2 TDR，
逐点核验原生 TDR 和候选 VTU 指纹、Bias_V、Converged、CmdSHA256、单位和节点配对。

| 偏置 (V) | 电势 RMS (V) | eDensity log RMS (dex) | hDensity log RMS (dex) | 电场 global / junction L2 | 结果 / 用途 |
| --- | --- | --- | --- | --- | --- |
| -0.05 | 6.005773e-06 | 6.325426e-05 | 6.325405e-05 | 0.0949433% / 0.0848761% | 诊断通过 |
| -1 | 5.469358e-06 | 5.304079e-05 | 5.30409e-05 | 0.287987% / 0.0584501% | 诊断通过 |
| -5 | 5.054145e-06 | 5.943143e-05 | 5.942932e-05 | 0.898269% / 0.0268889% | 诊断通过 |
| -10 | 6.246319e-06 | 4.331366e-05 | 4.331364e-05 | 0.594885% / 0.0182326% | 诊断通过 |
| -13.2 | 7.286545e-06 | 5.341493e-05 | 5.34145e-05 | 0.502622% / 0.0156126% | 通过 / 必过资格点 |
| -15 | 8.016653e-06 | 5.048204e-05 | 5.048204e-05 | 0.464813% / 0.014543% | 诊断通过 |
| -19.5 | 9.645062e-06 | 6.613418e-05 | 6.612821e-05 | 0.396231% / 0.0125779% | 诊断通过 |
| -19.6 | 9.685703e-06 | 6.71481e-05 | 6.713963e-05 | 0.395011% / 0.0125424% | 诊断通过 |
| -19.7 | 9.728508e-06 | 6.821529e-05 | 7.046285e-05 | 0.3938% / 0.0125071% | 诊断通过 |
| -20 | 9.869574e-06 | 6.865485e-05 | 0.01548827 | 0.390215% / 0.0124031% | **hDensity 诊断失败** |

只有 −13.2 V 是原契约的必过场量点。报告每个点都保留 `acceptance_required`、
独立 `status` 和原始测量；总体 `biased_fields.status` 使用必过点，
`diagnostic_status` 单独保留附加点失败。完整扫描、曲线、KCL、平衡场量仍全部是必过项。
因此不能把主验收的 `pass` 解释为表中所有诊断点通过。

## 击穿区与判据

用户指出 PN2D 的击穿区域应在约 |V|=19.6 V，而旧 13.2 V 是算法问题造成的提前表现。
本次没有把 −13.2 V 场量资格命名为击穿资格，也没有把扫描停止当作击穿判据。

复用当前 [PN2D 原子配置契约](contracts/pn2d_node_volume_policy_atomic_default_acceptance_v2.json)
及 `analyze_pn2d_avalanche_on_bv_parity.continuous_breakpoint`：在固定
−18、−18.5、−19、−19.25、−19.5、−19.7、−19.8、−19.85、−19.9、−19.95、−20 V
上，以连续两段 log10(|I|) 最小二乘拟合、0.001 V 网格搜索折点。
原契约的 V_break 参考误差上限保持 0.10 V；两条曲线须在完全相同的精确偏置取样，
不插值代替缺点。同样输出已有 1 dex/V 持续斜率交叉判据的有/无交叉结果。

原生拟合膝点为 **-19.405 V**，Vela 为 **-19.405 V**，
参考误差 **0.000 V**，通过 0.10 V 门限。斜率交叉结果为双方
`shared_no_crossing`，不能宣称该斜率判据已经触发击穿。−19.6 V 的原生电流为
−6.1500467e−16 A/µm，Vela 为 −6.1519811e−16 A/µm；该点四项场量诊断通过。

这里的拟合折点是具名曲线膝点指标，数值取决于上述窗口及定义，不能直接等同于某个
指定电流阈值、积分等于 1 的击穿判据或通用物理临界电压。不能为了对齐约 19.6 V
而修改拟合定义或硬编码答案。该结果确认当前 Vela 与同条件原生参考的膝区一致，
并保留近膝区空间场量尚未获得整体资格的限制。

## −20 V 空穴密度限制与严格复核

−20 V 空穴密度在 905 个活跃节点上 RMS 为 **0.01548827 dex**，p95 为
**0.01058107 dex**，超过 0.005 / 0.01 dex 门限。其余三项在该点通过。
示例节点 260，坐标 (0.890625, 0.015625) µm：原生空穴密度 134093.50 cm^-3，
Vela 106604.75 cm^-3；空穴准费米电势差约 −5.908 mV。节点最大误差仅作为诊断。

在原 −20 V HDF5 状态上做独立严格 Newton 复核，保持材料、网格、物理模型和
80 次源迭代上限：局部行 eps 1e−8、源门槛 1e−22、活跃密度门槛 1e5 cm^-3，
全局连续性容差 1e−6，源 floor 1e−22。4 次 Newton 更新后收敛，局部行最大比值
2.78451e−12，全局空穴连续性比值 6.66180e−15，而 hDensity RMS 仍为
0.01548827 dex。该复核是精度诊断，不替代主候选，不修改正式数值配置。

结果排除了原接受状态的收敛不足作为主要解释；具体空间算子/映射差异尚未定位，
没有把它校准掉或计为通过。该限制写入 `known_difference_ledger` 并随验收 JSON 返回。
本次仅资格化电势、两种载流子密度和矢量电场四项；不建立准费米、雪崩系数、
雪崩产生率或其他物理量的空间资格，也不推广到任意网格、材料或其他温度。

## 验证与复现

**原契约的 0V、IV、BV 和门控案例复审均通过，返回码均为 0。**

| 案例 | 复审结果 | 本地最终 JSON 报告 |
| --- | --- | --- |
| pn2d_0v | `pass` | `build/wp10_closeout_knee/pn2d_0v_u0sgfg_u/acceptance_report.json` |
| pn2d_iv | `pass` | `build/wp10_closeout_knee/pn2d_iv_2lti7st3/acceptance_report.json` |
| pn2d_bv | `pass` | `build/wp10_closeout_knee/pn2d_bv_1txlcbl2/acceptance_report.json` |
| gates | `pass` | `build/wp10_closeout_knee/gates_56amfbhm/acceptance_report.json` |

最终契约 SHA-256：`88c81240e90182d9d8a07bd07825454eea19439de02dc381357b59db7fd67e3d`。四份报告均绑定这一契约。

报告 SHA-256：

- pn2d_0v：`a016ef2fa8805b1ab06c5f951506879ee9629de42e3af2ff30741b812710302a`。
- pn2d_iv：`9e4dcae40c0d420743241289a3cabcaf12af0f9c313493f81810c8fde750094f`。
- pn2d_bv：`7f96270ae594c3de5608e65c375d2a68462a04cdffce6e9700282aa433031387`。
- gates：`7b6fe728e988055859ff796221030982dfdb16f1cecd206dfeb55a265c64b5c9`。

初始七点严格审查失败报告保留在
`build/wp10_closeout/pn2d_bv_jod_jhb2/acceptance_report.json`，SHA-256
`549f9e8c2ce594e126d3a66dd5486010b34b540544c7b230dfa2d74f81f73c43`。
−20 V 严格 Newton 精度复核位于 `build/wp10_closeout/refine_20/`。
数值数据来自既有修复候选；其 run_manifest 中记录构建时提交 `5c723a2e`，实际前端修复
已提交为 `0842e000`。本次复审不改这些历史运行记录，CSV/VTU 指纹随新报告记录。

测试包括输入/参考/配置漂移拒绝、静态资格不得继承、部分指纹和重复偏置拒绝、
缺失必过点、附加点失败必须留在报告、原生 TDR 篡改拒绝，以及故意注入
2 mV 电势、10% 载流子密度、5% 矢量电场误差的负对照。结区门控单独保持约束。
击穿膝点回归拒绝提前的拟合折点及缺失/重复的精确偏置。

Release `sdevice_contracts`、`sdevice_reference_artifacts`、`sdevice_deck_emit`、
`sdevice_acceptance_tools` **4/4** 通过；最终验收工具复核通过，含 **17** 个单元测试。
设置实际参考目录后的 `sdevice_e2e_pn2d_0v`、`sdevice_e2e_gates` **2/2** 通过，
并复审实际完整 IV/BV 候选。未改 C++、求解器、原生 cmd/par、扫描目标或原有门限。

普通 DD 使用 Eigen SparseLU，单线程，UCRT64 Release；HDF5/TDR 功能实际启用。
已有完整 IV/BV 运行保持 206/811 个接受点，分别到达 10/−20 V；本次逐项核验候选
cmd/par/grid/sidecar 后重新导出原生场量及独立计算，没有重新推进整个长扫描。

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts/sdevice_acceptance.py --case pn2d_bv --reference-dir reference_staging --output-dir build/wp10_reaudit --launcher build-release/vela.exe --tdr-importer build-release/sentaurus_import.exe --candidate-dir build/iv_bv_repair_final/pn2d_bv_8hyd8uui/candidate
ctest --preset windows-ucrt64-release -R '^sdevice_acceptance_tools$' --output-on-failure
```

不指定 `--candidate-dir` 会从平衡态重新运行，完整 BV 单线程扫描较慢，默认壁钟上限
7200 s。原生产物、候选场量/状态和精度复核输出仍在忽略目录中，不提交原件。

2026-10-06 后续：−20 V 空穴密度的向量 GradQF/局部 AD 修复及完整扫线核验在
[独立记录](sdevice_followup_2026-10-06.md)中跟踪。本报告保留原收口范围、失败证据和当时结论；
当前绑定资格以验收契约为准。
