# SimpleMOS PhuMob 原生局部成对源校准

日期：2026-09-09。分支 `codex/simplemos-sdevice-validation`。接续 [局部源与重启验证](simplemos_phumob_local_source_and_restart_validation_2026-09-09.md)。

**已完成用户明确授权的上传、20 次 Sentaurus DC 及结果取回。原生 DC 资格为 20/20，同源双幅度响应资格为 8/8，Vela–Sentaurus 端口响应对比为 8/8。** 两套求解器的原始响应最大相对差为 0.001002%。这校准的是八热点冻结成对源的局部漏端导数；没有执行有限 SRH 体积替换。

## 配置与输入身份

n19/n23 × Vd=0.05/1 V，Vg=0；网格分别 1480/1482 节点，300 K、Boltzmann、OldSlotboom、plain PhuMob、掺杂相关 SRH。沿用此前四个合格保存态。Vela 使用显式 element-box PhuMob 候选和原 all-cell SRH 体积，实际线性后端为 SparseLU；本轮只复用其已冻结的 20 次 DC 与独立端口校准，没有重新计算 Vela。

原生为 T-2022.03-SP2、Super、ExtendedPrecision(128)、Digits=12、ErrRef(n/p)=1e-2、RhsMin=1e-20、Iterations=40、ExitOnFailure。正常结束与原生求解日志、偏置及 KCL 共同构成原生资格；这不等同于 Vela 的逐行 1e-6 条件。四个零控制先通过后才启动 16 次正负计算。20 点最大 KCL/Id 为 2.656e-15，最大偏置误差 0.000e+00 V。

PMI 附加与状态无关的成对复合源，保留原 SRH。每节点物理积分源 `S_i=R_Vela(base)×(V_native−V_Vela)`，单位 particles/(m·s)；原生速率 `S_i/(V_native[m²]×1e6)`，单位 cm⁻³s⁻¹。幅度 α=0、±0.001、±0.0005。n23 八节点为 795、792、1114、794、791、1115、1203、1188，n19 按坐标匹配。全部 160 个运行节点检查通过，包含坐标、速率、幅度、源积分和 PMI 的 128 位精度标志。每个返回网格、状态和指令文件均与上传前冻结哈希核对。

上传目录 `/tmp/vela_simplemos_phumob_pair_source_20260909`。输入包 SHA256：`1c5ed4b063431dd101104f4a6ebb967f4ad0f37ff4899ae2aba0af87bd7875c6`。[本轮传输记录](../../reference_tcad/simplemos_sentaurus2022/phumob_local_source_20260909/native/transfer_completed.json)、[冻结原生合同](../../reference_tcad/simplemos_sentaurus2022/phumob_local_source_20260909/native/contract.json)、[源项身份](../../reference_tcad/simplemos_sentaurus2022/phumob_local_source_20260909/native/comparison/source_identity.csv)。

## 原生执行与编译兼容性

VM 的系统 GCC 为 4.8.5，原 cmi 命令使用 `-std=c++14` 而被编译器拒绝。保留该失败日志，沿用全部源码及其他编译参数，以 GCC 4.8 支持的 `-std=c++1y` 编译对象，再由 cmi 链接共享库。没有安装系统工具链或修改物理计算公式。零控制电流与原参考的导出值完全一致，全部任务成功加载 PMI 并记录精度标志 2（128 位）。

[原编译失败](../../build-release/phumob_local_source_20260909/native/zero_raw/compile_initial_gcc48.log)、[兼容编译日志](../../build-release/phumob_local_source_20260909/native/zero_raw/compile_cxx1y.log)、[链接日志](../../build-release/phumob_local_source_20260909/native/zero_raw/link_cxx1y.log)、[四个零控制](../../reference_tcad/simplemos_sentaurus2022/phumob_local_source_20260909/native/zero_points.csv)、[20个原生结果](../../reference_tcad/simplemos_sentaurus2022/phumob_local_source_20260909/native/response_points.csv)。复现时不能直接使用仍保留在冻结输入包中的原 `run_zero.sh` 编译步骤；应使用本轮传输记录保存的兼容编译命令。`run_response.sh` 不包含编译，按冻结脚本执行。

## 同扰动响应

`dId/dα=[Id(+α)−Id(−α)]/(2α)`。端口单位 A/μm。下表给出全幅值；两幅度全部明细见证据。

| NWell | Vd (V) | Sentaurus dId/dα | Vela dId/dα | Vela/native−1 | 双幅度对比资格 |
| --- | --- | --- | --- | --- | --- |
| n19 | 0.05 | +9.938418702e-19 | +9.938323750e-19 | -0.000955% | 2/2 |
| n19 | 1.0 | +5.127422492e-19 | +5.127373204e-19 | -0.000961% | 2/2 |
| n23 | 0.05 | +1.577348930e-20 | +1.577333824e-20 | -0.000958% | 2/2 |
| n23 | 1.0 | +3.066585853e-18 | +3.066556510e-18 | -0.000957% | 2/2 |

固定门槛为原生双幅度差 ≤1e-3、偶/奇比 ≤0.01、信号/零漂移 ≥100、原始跨求解器响应差 ≤1e-3。原生双幅度最大差 9.231e-08，最大偶/奇比 1.969e-07。电流导出为 15 位有效数字；零漂移噪声取导出分辨率下限，最小信号/噪声 2.564e+06，最大舍入误差/信号 1.950e-07。不将导出文本一致解释为无限精度。

另存按手册印刷 q=1.602192e-19 C 与 Vela q=1.602176634e-19 C 归一化的诊断结果，最大差 0.000043%。手册值可能仅为四舍五入，不能据此确认原生内部常数，也不替代原始响应门槛；本轮没有更改任何常数。

[双幅度响应对比](../../reference_tcad/simplemos_sentaurus2022/phumob_local_source_20260909/native/comparison/responses.csv)、[四端导数](../../reference_tcad/simplemos_sentaurus2022/phumob_local_source_20260909/native/comparison/terminal_derivatives.csv)、[汇总](../../reference_tcad/simplemos_sentaurus2022/phumob_local_source_20260909/native/comparison/summary.json)。

## 结论边界与后续

原生同源实验现已执行，其资格以上表为准。此前 Vela 漏端切向/伴随/实际响应为 8/8，但全场响应仍为 6/8；n23 低 Vd 的全场信号/零漂移 68.64/34.32 未过 100。原生端口通过不能覆盖该缺口。

n23、Vg=0、Vd=1 V 的现有 PhuMob 候选 Id 误差仍为 −0.80813524%，绝对差约 −3.00069e-18 A/μm。局部源正向导数只能说明这一小扰动的影响，不能将 α=1 线性外推当作有限自洽重算结果。Vg=0.2、Vd=1 的失败初始化也仍未解决：隔离重启修复后的最差电子 1089 行比为 1.045902376e-6，超过 1e-6。

下一步优先在独立重启精度轴上补齐低 Vd 全场漂移，并定位电子 1089 失败态的实际线搜索残差精度与 Poisson 小更新。只有这些资格补齐后，才考虑状态相关 SRH 体积插值、双初始化有限替换及 Enormal/HFS 恢复。

本轮未改生产 C++、默认、SRH 体积或接受门槛，未新跑全量 CTest，未提交或推送。仅原生仿真、分析脚本和证据文档发生本轮变化；之前所有失败与历史报告保留。
