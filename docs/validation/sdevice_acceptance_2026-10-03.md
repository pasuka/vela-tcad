# SDevice WP10 一致性验收，2026-10-03

此文保留首次 WP10 验收快照；后续 IV 密度与 IV/BV 扫描修复及当前资格边界见
[IV/BV 修复报告](sdevice_iv_bv_repair_2026-10-03.md)。以下失败数值属于修复前状态。

**WP10 工具、参考重跑、C6 固化和输入门控已经交付；整体数值验收未通过。**
PN2D 0V 通过；IV 比较窗口通过，但全扫描和 0.2 V 电子密度失败；BV 全扫描失败。
不将缺失偏置、部分扫描、参考缺失或未获资格的阈值计作通过。
未修改 C++、材料参数、物理默认、扫描目标或既有容差。

本记录覆盖 worktree `sdevice-cmd-frontend`，基于 `c9971ebd`，此前 WP4/WP6.1/WP5
提交保持原样。执行期间发现并修复一个必要的 WP5 Python 打包问题：标准库 zipapp
自动入口忽略 `main()` 的返回值，原生启动器因此在打印 E105 等拒绝错误后仍退出 0。
入口改用 `raise SystemExit(main())`；新增隔离 zipapp 的 1/2/3 返回码及原生启动器
拒绝返回码回归。此前 WP5 的成功路径证据仍成立，失败返回码以本记录为准。

## 参考资格与归档

用户明确授权在 VM 上最小重跑 PN2D 0V/IV/BV 和 BJT m0 两段扫描、复用网格并写入 v2。
原生运行目录 `/home/tcad/vela_sdevice_wp10_20261003`；正式归档
`/home/tcad/vela_reference_artifacts/{pn2d_sentaurus2018,genius_bjt_sentaurus2022}/v2`。
版本 **T-2022.03-SP2 (0.7745337)**。四个原始扫描均退出 0，PN2D IV 到 10 V、BV 到 −20 V，
BJT base 到 0.70 V、collector 到 3 V。原有 v1 未改变。

只增加显式 `Parameter="models.par"`；0V/IV 另外补充 `Plot{ElectricField/Vector}`，
以满足矢量场验收。两次输出补跑同样退出 0，原始运行及替换前的产物单独保留。
首次批量归档因两份 `console.log` 重名中止；改用案例前缀后复用成功结果继续归档，
未因此重复求解。`wp10_inputs.json` 记录首次输入，`wp10_*output_augmentation.json`
记录最终仅增加输出的输入指纹。

网格和同版本 `sdevice -P` 原生默认参数均复用 v1。参数 SHA-256：
`b4b3ebfdefba530f756f3855d43d7d587720689771d8badc747b61439ed42742`。
PN2D 同一硅网格 **1943 节点、3680 三角形**；BJT m0 网格 **5611 节点**。
坐标 µm、浓度 cm⁻³、电流 A/µm、电场 V/cm，300 K。
本地 UCRT64 Release，HDF5/TDR 和 HDF5 状态存储实际启用；UMFPACK/SPQR 编入，
本次普通 DD Newton 路径使用 Eigen SparseLU，不把可选后端编入等同于使用。

PN2D **626** 文件、BJT **21** 文件先按 VM 逐文件 SHA-256/大小及传输包 SHA-256 核验，
再用 WP9 工具完成六项核对；最终均通过，无 partial/warn。
原生 SDevice 日志回显、producer cmd/par 哈希、PLT 目标、电极数据及 TDR 可读性均核验。
新增 [PN2D v2 清单](../../reference_tcad/pn2d_sentaurus2018/artifact_manifest_v2.json) 和
[BJT v2 清单](../../reference_tcad/genius_bjt_sentaurus2022/artifact_manifest_v2.json)。
它们的 `deck_matches_repository=false` 如实表示字面源码新增了声明/输出；验收驱动另行检查
允许的机械变换，确保所有其余源 token 不变。未用这个标记绕过物理等价检查。
仅清单及小型证据进入 Git，原生 TDR/PLT/日志保存在 VM 和忽略目录。

v1 盘点也保留失败：PN2D、BJT、Schottky 的边界 TDR 被导入器拒绝；不把这些旧归档当作
数值验收通过。BVmethods 旧归档与仓库 deck 物理不一致，仍仅用于输入门控。
其 v1 六项核验也未通过：ABA PLT 的全局 `MaxField` 数据组不被当前电极型 PLT 读取器接受，
`n1_half_bnd.tdr` 的边界元素编码同样不受导入器支持。未将这类旧参考失败计入通过。

## C6 符号和 Time

全局 `CURRENT_SIGN_FACTOR=+1.0` 已固化。原生 `pn2d_iv.plt` SHA-256：
`a5f712e3ba0bfeea221e0de681ce380b3ad4ab14eb42c63b1324eac431709b25`。
比较同一 **Anode**，不借用 Cathode，不引入逐案例负号。

| 首三个正向参考偏置 V | 原生 Anode A/µm | Vela 原始 Anode A/µm | 因子 |
| --- | ---: | ---: | ---: |
| 0.01 | 5.80543692864973e−19 | 5.3887968315032304e−18 | +1 |
| 0.02284 | 1.62155570626151e−18 | 1.3977546608561714e−17 | +1 |
| 0.03942928 | 3.74720553099733e−18 | 2.9009156046579526e−17 | +1 |

Vela 两端点间采用已接受原始 CSV 的 log-current 插值核对符号，仅证明方向；这些小偏置的
幅值差异不被声明为通过。三个点进入 PLT writer 回归。原始 IV CSV SHA-256：
`46ef5f1333c64836265ebf0a351bcfb6e7600d8c0dd065a650302a4600278752`。

BJT base 文件 36 行、collector 文件 31 行，两个文件首末 Time 均为 **0、1**，逐段重置，
不累加；现有 WP5 拼接语义无需改变。两个 PLT SHA-256 分别为
`dd05b1888149055604614241312f51e8ba3da32028f5cf7b40d131d5aa391c26`、
`0aff50fa989d04249b606e9b811a1fd655a7ea545ab4dcefe2e8293bfae0ed9b`。
结论已回填 [workflow §3.3](../sentaurus_sdevice_workflow.md#33-plt-与-vtu)，并有重置/累加拒绝回归。

## 数值结果与未通过项

契约 [sdevice_acceptance_contract.json](../../configs/sdevice_acceptance_contract.json) 的曲线门限
引用原 PN2D IV **0.2～0.3 V / 0.4 dex**；守恒引用现行 E2E 的 **2% / 1e−20 A/µm**。
相对误差同时报告，并由 0.4 dex 导出上界 `10^0.4−1`，不另造或放宽门限。
空间门限引用已验证静态、300 K、无雪崩的 PN2D variants：电势 RMS 0.005 V/p95 0.01 V，
有效参考密度 ≥1e5 cm⁻³，密度 log RMS 0.05 dex/p95 0.1 dex，电场面积加权矢量 L2 10%。
参考/候选节点 ID 和坐标必须双射，字段单位、分量、有限性及 VTU 命令/偏置元数据均检查。
结区统计采用 x=1±0.1 µm，作为定位诊断，不以未经资格的区域门限替代全域门限。

| 案例 | 曲线 | 场量 | 扫描完成 | 最终状态 |
| --- | --- | --- | --- | --- |
| PN2D 0V | 最大电流 6.89e−25 A/µm，通过 | 电势 RMS 1.00e−5 V；密度 RMS 1.82e−4 dex；矢量场 L2 0.1009%，通过 | 是 | 通过 |
| PN2D IV | 0.2～0.3 V：3 点，0.373 dex、最大相对误差 57.64%、同号且同趋势；最坏 KCL 1.526%，通过 | 平衡态通过；0.2 V 电子密度 RMS **0.2764 dex**、p95 **0.6990 dex**，失败；电势/空穴/电场通过 | **否，止于 0.8073157685 V** | 失败 |
| PN2D BV | 未覆盖 −20～−0.05 V 比较窗口 | 平衡控制通过；未获得 −13.2 V 候选状态 | **否，止于 −3.3090559726 V** | 失败 |

IV 最后一次尝试的接触一致性验证失败：Anode node 2 的 `phin=0` 与偏置不匹配，
CSV 标记 `validation_failed`。BV 记录 `line_search_non_decrease`，最小步附近 Newton
继续失败，最终 runner 退出 1；两个 CLI manifest 都是 E502/退出 **5**，`goal_reached=false`。
没有截短 deck、放宽验收、用平衡态替代带偏状态，或为了通过而改变物理默认。
BV −13.2 V 历史 0.70x 偏差仅登记账本，本次无法计算该偏置的新比例；它不能豁免目标未到达。
既有无雪崩空间资格也不自动推广到雪崩带偏状态，契约明确保留该资格缺口。

Schottky 和支持的 BVmethods（ABA Poisson、ABA coupled、external resistor、voltage-to-current）
共五个 `vela device --check` 全部退出 0，无求解。2022 transportmodels `pp6_des.cmd`
退出 1，正确记录 **E105 / E109**（Thermode / Continuation），未因缺文件改判退出 3。
不把不支持的 Continuation/Transient 当作可通过的 BVmethods。

## 验证与复现

- 新验收工具单元测试 10 项通过；覆盖参考缺失 skip、六检查/partial 拒绝、源码变换、
  符号、插值同号、曲线覆盖、KCL、节点/单位/分量、ASCII/inline binary 和逐段 Time。
- 增加 PLT 三点符号回归和 zipapp/原生拒绝返回码回归。
- Release 相关 CTest 15 项：11 通过，4 参考用例因未设置环境变量明确 skip。
- Release 全量 **1025** 项：**1020 通过、1 失败、4 skip**，240.87 s。
  唯一失败为既有 IALMob 测试读取临时 JSON 失败，独立重跑 1/1 通过；未修改该测试或其容差。
- 最终新增验收/PLT writer 单元合计 **26/26** 通过；包括三点符号回归。
  最终 Release CTest `sdevice_acceptance_tools` / `sdevice_plt_writer` / `sdevice_orchestrator`
  **3/3** 通过，包含 zipapp 与原生返回码、真 runner 和状态链检查。
- 显式设置参考环境变量后，原网格 `sdevice_e2e_pn2d_0v` **1/1** 通过；
  `sdevice_e2e_gates` **1/1** 通过。第一次这两项复核中，重打包与门控启动并行，
  第一个命令读到未写完的 pyz；构建完成后的门控重试通过。没有将该次失败隐藏为通过。
  IV/BV 的完成运行经驱动严格核对输入哈希后复审，均明确返回失败，与上述数值结果一致。

从当前 worktree 根使用 Release preset。WP9 `fetch` 不支持计划旧示例中的 `--version`，
改用显式 v2 清单：

```powershell
python scripts/reference_artifacts.py fetch --case pn2d_sentaurus2018 --manifest reference_tcad/pn2d_sentaurus2018/artifact_manifest_v2.json --dest reference_staging
python scripts/reference_artifacts.py fetch --case genius_bjt_sentaurus2022 --manifest reference_tcad/genius_bjt_sentaurus2022/artifact_manifest_v2.json --dest reference_staging
$env:VELA_REFERENCE_ARTIFACTS_DIR = (Resolve-Path reference_staging).Path
ctest --preset windows-ucrt64-release -R '^sdevice_e2e_' --output-on-failure
```

门控案例还需按原 v1 清单 fetch Schottky/BVmethods。只有本地原生产物齐备并通过六检查和
源码变换检查后才启动数值候选；环境变量缺失返回 77，目录或清单损坏返回失败。
单案例运行、JSON 报告及已有完成运行的严格哈希复审：

```powershell
python scripts/sdevice_acceptance.py --case pn2d_iv --reference-dir reference_staging --output-dir build-release/sdevice_acceptance --launcher build-release/vela.exe --tdr-importer build-release/sentaurus_import.exe
# 可加 --candidate-dir <已完成运行的cmd目录>，cmd/par/grid必须逐字节哈希相同。
```

本次详细 JSON、六项核验、所有候选 CSV/VTU/状态/失败诊断位于忽略目录
`build/wp10_acceptance_final/`、`build/wp10_probe_iv_vector/`、`build/wp10_probe_bv/`、
`build/wp10_verify/`、`build/wp10_native/`；CTest 日志为 `build/wp10_ctest_*.txt`。
后续必须定位 IV minority-density/contact 验证与 BV Newton 失败，并补足雪崩空间门限的适用证据，
才能把 WP10 更新为全部验收通过；这些问题不以物理参数校准或容差放宽代替修复。
