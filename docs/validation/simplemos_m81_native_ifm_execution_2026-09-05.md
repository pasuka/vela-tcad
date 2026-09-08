# M81：Sentaurus 原生零频 IFM 与固定电荷校准

2026-09-05，工作分支 `codex/simplemos-sdevice-validation`。

**结论：M81 两工作点的原生电流响应校准通过。** n19/n23、Vd=0.05 V、Vg=0.9 V 的零频 IFM、空间 Green 输出及两档固定电荷中心差分均完成。四组 IFM–FD 比较最大误差为 **1.645207%**，符号全部一致，满足事先冻结的 5% 门槛。尚未确定可进入 M82 的因果明确修正，高 NWell 的 Vela–Sentaurus Id–Vg 差异仍未解决。

用户本轮明确授权上传两工况的网格、状态和指令至 Sentaurus 虚拟机并执行 M81；此前上传阻断已解除。旧 `m81_execution_status.json` 和 `simplemos_m79_m81_execution_progress_2026-09-05.md` 是保留的历史记录；当前状态以本报告及 `sentaurus_ifm_pilot/m81_fd_calibration_result.json` 为准。

## 执行范围与参考身份

- 虚拟机 SSH 别名：`sentaurus`；原生运行版本：`T-2022.03-SP2`。
- 远端隔离目录：`/tmp/vela_simplemos_m81_ifm_20260905/`。首轮、输出修订、固定电荷校准分别位于 `bundle`、`output_amendment_bundle`、`fixed_charge_bundle`，未覆盖旧运行。
- n19 是低 NWell，n23 是高 NWell；其 NWell 为 p 型 Boron 体区，分别为 1e17、2e17 cm^-3。
- 使用 M69 已保存工作点及原网格；每次上传输入均通过本地/回收副本 SHA256 核对。物理保持 M65/M69 的 no-BGN 诊断分支与原 mobility/SRH 组合，不能与生产 BGN-on 曲线混用。
- `ImplicitACSystem`，source/substrate=0 V、drain=0.05 V、gate=0.9 V；保留原端口定义和状态 Load。2D 默认厚度为 1 um。
- 本轮没有修改 Vela 核心物理、生产默认值或材料参数，没有新增 Id–Vg 扫描。

## 原生 pilot 与输出修订

首轮两个任务返回码均为 0，输出了命名的 `dIuniform_charge(N_drain)`，但没有输出空间 Green TDR，故首轮记录为 `failed_missing_green_output`，没有按成功处理。

按本地 Sentaurus User Guide pp787–788、815，冻结输出修订，仅新增未命名的 `Noise(DiffusionNoise(LatticeTemperature))` 分组以启用 NoisePlot。修订后生成两个 `m81_green_N_drain_0000_acgf_des.tdr`。两轮的 DC 电流及确定性 dI 数值完全一致；DC `.plt` 文件哈希也分别一致。

| 检查项 | n19 | n23 |
|---|---:|---:|
| 重闭合 Id（A/um） | 1.95731960115792e-6 | 1.56350743583972e-8 |
| 对既定参考的最大重闭合误差（dex） | 3.115527e-8 | 3.506672e-10 |
| +1e13 cm^-3 donor variation 原生 dI（A，默认 1 um 厚度） | 1.01255967603074e-9 | 1.09657877411713e-11 |
| dI/Id | 5.1731954e-4 | 7.0135821e-4 |
| 硅区 Green 节点数 | 942 | 942 |
| DC、输出和原生响应资格 | 通过 | 通过 |

两工况均导出了有限、非零且节点映射完整的 `CurPotReACGreenFunction`（s^-1）、`CurECReACGreenFunction`（无量纲）、`CurHCReACGreenFunction`（无量纲）。CurPot 覆盖四个材料区域；载流子 Green 函数在硅区。

另以原始 M69 TDR 对重闭合 TDR 做只读逐节点核对：电势/准费米势最大差为 **1.9573201e-9 V**，载流子密度最大差为 **3.2828693e-8 dex**。直接对 M69 的工作点电流，n19/n23 重闭合差分别为 3.1155161e-8、4.1798726e-11 dex。这是补充状态身份审计，没有事后修改原先的 Id 放行门槛。

ACExtract 的 `i(,N_drain)` 与器件 `drain TotalCurrent` 符号相反。确定性 `dIuniform_charge(N_drain)` 的符号另由下面的显式 FD 验证，未机械沿用电路电流符号。

## 两档显式固定电荷校准

pilot 成功后，先冻结 `simplemos_m81_fixed_charge_fd_contract_v1.json`，再上传和运行。硅体区加入 `Traps(FixedCharge Conc=...)`；使用 +1e13、-1e13、+5e12、-5e12 cm^-3，每个器件另加 Conc=0 对照。保留相同 Load、ImplicitACSystem 和 DC 边界，移除只用于 IFM 的噪声与确定性扰动配置；不更改掺杂剖面或 mobility、BGN、SRH 参数。语法、体浓度单位和电荷符号依据手册 pp543–545。

比较公式：`central_delta = [Id(+a) - Id(-a)] / 2`，预测值为 `native_dI × a/1e13`。没有拟合正负号、归一化倍数或物理参数。

| 器件 | 幅度（cm^-3） | 原生 IFM 预测（A/um） | DC 中心差分（A/um） | 相对误差 | 资格 |
|---|---:|---:|---:|---:|---|
| n19 | 1e13 | 1.0125596760e-9 | 1.0254640577e-9 | 1.274432% | 通过 |
| n19 | 5e12 | 5.0627983802e-10 | 5.1273203401e-10 | 1.274433% | 通过 |
| n23 | 1e13 | 1.0965787741e-11 | 1.1146185950e-11 | 1.645100% | 通过 |
| n23 | 5e12 | 5.4828938706e-12 | 5.5730988278e-12 | 1.645207% | 通过 |

正扰动相对零对照均提高 Id，负扰动均降低 Id。两档导数相对变化：n19 为 1.0053861e-8，n23 为 1.0501620e-6；最大偶次非线性贡献为中心差分信号的 3.7069e-4。最小信号/零对照漂移比超过 4.0e5。全部 10 个 DC 任务正常退出。

**限制：1.27%–1.65% 的偏差基本不随扰动幅度变化。** 不能将其解释成已通过减小幅度消除的非线性误差，也不能宣称 DopingVariation 与 FixedCharge 的离散源项逐节点完全等价。本轮仅按事先冻结的 5% 标准资格化该均匀电荷方向。手册 p196 明确说明 IFM 默认采用完整导数，不能未经证据把剩余差异归咎于忘记显式添加 `Derivatives`。

## 复用 M34 的原生 Green 积分校准

另外冻结了只读积分合同，复用 M34 已有 n23 `MeasureCoefficients.debug`，未重新运行几何导出。M34 网格哈希与当前 n23 输入完全一致，节点坐标与单元顺序核对通过。采用 M34 已资格化的局部排列 `[0,2,1]` 累加硅区 signed native measure，942 个硅节点的面积和为 0.9949540618291826 um^2。

使用 `delta_I = q × delta_N × sum(Gpot × native_Si_measure) × 1e-12`，其中 q=1.602176634e-19 C，delta_N=1e13 cm^-3，1e-12 为 um^3→cm^3 换算（2D 厚度 1 um）。体积仅乘一次，未拟合符号或系数。

- 原生空间积分：1.0965682572568285e-11 A/um。
- Sentaurus 原生 dI：1.09657877411713e-11 A/um。
- 相对误差：**9.5906109e-6，即 0.000959061%**，符号一致。

这补充资格化了 **n23、均匀硅区电荷方向** 的原生 Green 单位、符号和体积加权。没有将此结论推广到 n19 的原生几何、任意局部源形状、连续性方程源项或高 Vd。

## 计数及后续放行

本轮执行 2 次首轮显式 DC 重闭合 + 2 次首轮 AC/IFM；2 次输出修订显式 DC 重闭合 + 2 次修订 AC/IFM；8 次扰动 DC + 2 次零对照 DC。合计 **14 次显式 DC、4 个 AC/IFM 点、0 次新偏压扫描**。ACCoupled 内部的工作点预解包含在 AC/IFM 任务内，没有称作无 DC 工作。M34 积分和状态核对均为只读，0 新求解。

M81 两工况校准完成；M82/M83 尚不放行。已有 M80 的最大单项是 electron_transport，但聚合 Poisson 项在 6/8 配对更大；当前 M81 校准的是均匀电荷方向，尚不能证明某个传输或源项离散修改可同时改善绝对误差和高低 NWell 配对增量。

进入 M82 前，应先把 M80 的具体候选残差转换成定义明确的方程源项和空间支持，使用已校准的原生权重进行同扰动对照。若涉及电荷源，可先复用 n23 原生体积；若涉及载流子输运，需要另外校准连续性方程扰动。仍须覆盖既定高 Vd 反例，形成因果明确的单一候选后再冻结 A/B。此处没有额外启动矩阵扩展或生产支路修改。

## 文件与复现

主要便携账本位于 `reference_tcad/simplemos_sentaurus2022/sentaurus_ifm_pilot/`：

- `m81_initial_native_result.json`：保留首轮缺失 Green 的失败记录。
- `m81_native_pilot_result.json`：补输出后的 pilot 放行、字段统计和原始文件哈希。
- `m81_state_closure_audit.json`：源状态与重闭合状态的逐字段差异。
- `m81_fd_case_ledger.csv`、`m81_fd_calibration_ledger.csv`：10 次 DC 执行和 4 组中心差分。
- `m81_fd_calibration_result.json`、`m81_fd_execution_evidence.json`：资格结果及冻结证据。
- `m81_native_green_integral_result.json`：n23 原生体积积分结果。

本地大体积原始结果和导出字段保留在忽略目录 `build-release/m81_sentaurus_ifm_pilot/`，没有提交仿真二进制输出。

在 UCRT64 环境运行以下只读检查：

```powershell
python scripts/run_simplemos_m81_sentaurus_ifm_pilot.py --verify-inputs
python scripts/analyze_simplemos_m81_fixed_charge_fd.py --verify
python scripts/freeze_simplemos_m81_native_execution.py --verify
python -m unittest discover -s tests/regression -p "test_simplemos_m81*.py"
```

准备和分析脚本拒绝覆盖已有合同/结果；复验使用 `--verify`，需要重新仿真时应建立新修订目录。此次新增的 5 项回归覆盖中心差分非线性消除、符号误用、对照漂移、零信号及重复数据列；没有改动求解器行为。
