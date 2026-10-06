# SDevice 栅接触参考能级核对（WP4a，2026-09-29）

WP13 更新（2026-10-05）：本文保留当日旧内置材料的数值证据。当前检查器改用无 BGN 300 K 默认快照，
`Φ_S,i=4.611003610676079 eV`，平带预期也由新 `ni` 计算；冻结公式与 1 mV 门限不变。
见[材料对齐与复验记录](builtin_material_alignment_2026-10-05.md)。

## 结论

Sentaurus 绝缘体上接触（手册 Eq. 104）与 Vela `metal_gate` 接触使用**同一参考能级**
（本征半导体的 Fermi 能级：两者的欧姆接触都写成 `ψ = V + Vt·asinh(N/2nᵢ)`），因此映射是精确的：

| SDevice 写法 | Vela `contacts[]` 字段 | 公式 |
| --- | --- | --- |
| `Barrier = Φ_MS`（默认 0） | `type: "metal_gate", work_function_eV: Φ_MS` | `ψ_gate = V − Φ_MS`（Eq. 104）↔ `ψ_gate = bias − work_function_eV` |
| `WorkFunction = Φ_M` | `work_function_eV: Φ_M − Φ_S,i` | `Φ_S,i = χ + E_g/2 + (kT/2)·ln(N_c/N_v)`，取 Vela 材料文件中该参考半导体的 χ、E_g、N_c、N_v |
| `Material = "Gold"` | 同上，`Φ_M` 取 `.par` 的 `Bandgap{WorkFunction}`（WP3 报告 `electrodes[]`） | 同上 |

平带电压（无氧化层电荷、均匀掺杂、背面欧姆接触 0 V）：`V_FB = Φ_MS + Vt·asinh(N_net/(2nᵢ))`。

注意 Vela 字段名 `work_function_eV` 的语义是**相对本征参考的功函数差**（即 Φ_MS），不是金属绝对功函数；
`flatband_voltage` 与之等价（互斥）。C4 发射规则据此写入 [`docs/sentaurus_sdevice_workflow.md`](../sentaurus_sdevice_workflow.md) 第 2.5 节。

## 夹具与方法

- 脚本：[`scripts/check_sdevice_gate_contact_mapping.py`](../../scripts/check_sdevice_gate_contact_mapping.py)；
  回归：[`tests/regression/test_sdevice_gate_contact_mapping.py`](../../tests/regression/test_sdevice_gate_contact_mapping.py)
  （CTest `sdevice_gate_contact_mapping`；无 runner 时夹具用例跳过，公式用例始终运行）。
- 几何：二列条带（宽 0.05 µm），Si 厚 1 µm（24 层，界面处加密）+ SiO₂ 10 nm（2 层）；
  `substrate` 欧姆接触在 y=0，`gate` 为 `metal_gate` 接触在 y=1.01 µm。单位 µm，`scaling.mode=unit_scaling`。
- 材料：Vela 内置 Si/SiO₂（nᵢ=1e10 cm⁻³，χ=4.05 eV，E_g=1.12 eV，N_c=2.8e19，N_v=1.04e19 cm⁻³），T=300 K。
- 掺杂：Si 均匀 N_A=1e17 cm⁻³ → `Vt·asinh(−N_A/2nᵢ) = −0.416685 V`。
- 求解：`solver.method=newton`（reltol 1e-10，abstol 1e-14，无复合），`sweep.bias_points` 显式偏置；
  每点写 `.vtu`（`vtk_format=vtu`，`vtk_fields=[ElectrostaticPotential]`，WP6），读取界面节点与体节点电势。
- 提取：`f(V_g) = ψ_surface − ψ_bulk` 单调，先在预测值 ±20 mV 取 5 点做线性内插求根，再以 ±2 mV、±0.2 mV 窄窗各细化一次。

## 结果（`build/gate_contact_check/gate_contact_check.json`）

| Φ_MS / V | 预测 V_FB / V | 测得 V_FB / V | 偏差 |
| ---: | ---: | ---: | ---: |
| 0.00 | −0.4166850053 | −0.4166850053 | −2.2e-16 V |
| +0.55 | +0.1333149947 | +0.1333149947 | −2.2e-16 V |
| −0.30 | −0.7166850053 | −0.7166850053 | −2.2e-16 V |

三个 Barrier 值的偏差均在浮点精度（≪ 1 mV 验收），三次细化估计一致到 1e-14 V。
偏差恰为机器精度并非巧合：平带时 Vela 解在半导体内严格均匀，`ψ_gate = V_g − Φ_MS = ψ_bi` 精确成立。

`WorkFunction` 形式示例（Vela Si 参数，`Φ_S,i = 4.6228018936 eV`）：Φ_M=4.05 eV → Φ_MS=−0.5728 V；Φ_M=5.1 eV → Φ_MS=+0.4772 V。

## 限制

- `WorkFunction`/`Material=` 形式依赖“本征参考半导体”的定义；Vela 侧按 `χ + E_g/2 + (kT/2)ln(N_c/N_v)`
  自洽实现（与 Vela 的 `ψ=0 ↔ n=p=nᵢ` 约定一致）。Sentaurus 是否在 Φ_S,i 中包含 `(kT/2)ln(N_c/N_v)` 项
  （≈ 12.8 meV）尚无本地参考运行可核对，待 A-3 重跑时用一个 `WorkFunction=` 栅极的 MOS 算例确认；
  `Barrier` 形式不受此影响。
- 夹具只覆盖无界面电荷、均匀掺杂、等温 300 K；氧化层固定电荷（`Traps(FixedCharge)`）与多晶硅栅
  （`Material="Silicon"(N=...)`，手册以 `(kT/q)ln(N/nᵢ)` 近似）不在本核对范围。

## 后续原生补验（2026-10-06）

本报告保留当时旧内置材料夹具的结果。WP14 已用当前参数完成三个 Barrier 和两个 WorkFunction 控制，
确认 T-2022.03-SP2 的栅参考含 DOS 项且使用 `Eg(0)`；当前库的原生零偏移输入是
`4.633743995291464 eV`。物理本征功函数 `4.611003610676079 eV` 保留为非零偏移控制。
五组 Vela/原生平带差均小于 1 mV；当前转换、哈希与适用范围见
[后续记录](sdevice_followup_2026-10-06.md)和 workflow §2.5，不沿用本报告的旧材料数值。
