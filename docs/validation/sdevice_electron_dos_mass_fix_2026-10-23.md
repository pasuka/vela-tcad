# SDevice 电子 DOS 有效质量修复（2026-10-23）

已修复材料发射器 `electron_dos_mass`／`evaluate_semiconductor`：
启用 OldSlotboom 时，其常数带隙偏移同时进入电子 DOS 质量的两个参考带隙：

```text
mt = a * (Eg(0) + dEg0) / (Eg(T) + dEg0)
me/m0 = [(6 mt)^2 ml]^(1/3) + mm
```

ni 仍通过含该偏移的有效带隙求指数项，不重复计入掺杂相关 BGN。
未启用 OldSlotboom 时偏移为零。公开 `bandgap_eV` 保持原有 Varshni 语义，
材料报告新增两个 DOS 有效带隙字段，provenance 明确记载偏移用于质量和 ni。
本次保留 CODATA 2018 常数，未修改全局 C++ 常数或显式历史材料 JSON。

依据：[新原生参数导出审计](bvmethods_material_export_2026-10-23.md)。
机读验证：[修复摘要](sdevice_electron_dos_mass_fix_2026-10-23.json)。

## 具体常数与剩余差异

以下为本轮对照绑定的数值；原生常数控制不是声称由 `-P` 直接导出。
Nc/Nv 的基础因子使用 kB、h、m₀，ni 的指数和热电压另使用 q。

| 常数 | Vela 当前值 | 原生 T-2022.03 对照值 | 单位 |
|---|---:|---:|---|
| kB | 1.380649e-23 | 1.380662e-23 | J/K |
| q | 1.602176634e-19 | 1.602192e-19 | C |
| h | 6.62607015e-34 | 6.626176e-34 | J·s |
| m₀ | 9.1093837015e-31 | 9.109534e-31 | kg |

Poisson 对照还涉及 ε₀：Vela `8.8541878128e-12`、原生控制 `8.8542e-12 F/m`；
它不参与本次 DOS 质量和 ni 计算。

300 K 新发射的 Si 材料为 Nc=`2.85827787833203e19`、
Nv=`3.104657033455896e19`、基础 ni=`1.463899277199692e10 cm⁻³`。
相对新原生导出/推导值，Nc/Nv 各剩余约 `+9.05132e-6`、ni 约 `+1.27988e-5`。
独立换用原生常数后的计算误差约 `1.96e-10`；生产中没有启用该常数控制。

## 检查与影响范围

- 材料模块 35/35 单元测试通过。新增独立原生 Nc/Nv 对照、250/300/350 K 的
  模型启用条件与质量温度关系、非正有效带隙拒绝检查。无 BGN 的内置快照和
  历史显式材料文件一致性检查通过。
- 历史 PN2D ni 字面量未修改；修复后偏差约 `5.32e-6`，其旧保护门从 `5e-4`
  **收紧至 `1e-5`**。新原生 DOS 回归独立于该历史字面量，容差为 `1e-9`。
- 7 项定向 Release CTest 最终通过：`sdevice_materials`、`sdevice_deck_emit`、
  `sdevice_orchestrator`、`sdevice_gate_contact_mapping`、`sdevice_acceptance_tools`、
  `wp14_remaining_sources`、`pn2d_config_templates`。首轮 deck 测试的旧期望暴露了
  OldSlotboom DOS 修正带来的约 `7.23 μV` 栅参考变化；按独立原生 DOS 比值更新期望，
  保持原容差，定向重测通过。
- 当前生产发射器重新读取 BVMethods 用户参数，发射新材料并执行独立
  Poisson→零偏联立冒烟检查。2600 节点、5020 三角形、300 K、1e7 Ω·μm，
  整条 Band2Band 关闭，采用已明确的几何控制；Windows UCRT64 Release、UMFPACK、
  单线程。Poisson 7 次更新、联立 0 次更新，最大电势差 `0.151654 mV`；
  电子/空穴最大相对差约 `0.0104701% / 0.0104261%`，通过本轮场量门。

该冒烟没有重新验证完整 J/RHS、有限 RHS 更新、BV 扫描或原版 NonLocalPath；
此前尚未通过的资格仍未通过。本次没有改 C++ 求解器，因此未重新编译核心或运行完整 CTest。
已保存的材料 JSON 和已有发行包不会自动更新；修复作用于当前源码之后重新发射的材料。
原始日志、配置、材料、状态和验证脚本保留于 ignored
`build-release/bvmethods_material_export_20261023/`。
