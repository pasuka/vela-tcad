# SimpleMOS M45–M77 当前验证工作交接

日期：2026-09-05  
工作树：`D:\code-repo\vela-tcad\.worktrees\simplemos-sdevice-validation`  
分支：`codex/simplemos-sdevice-validation`  
handoff 基线 HEAD（不含本文提交）：`72551ec`  
已合并本地 main：`8e99e57`，合并提交 `3a46f85`

## 1. 新聊天的最短恢复路径

1. 进入上述工作树，不要切换到仓库主 checkout。
2. 先阅读本文，再阅读 M60、M65、M72、M73、M74 和 post-main M77 报告。
3. 执行本文末尾的只读复核命令。
4. 下一任务建议定义为 **M78：M74 电子 Poisson 电荷体积响应的阶段/空间归因**；必须先写合同并冻结，再分析。
5. M78 第一阶段只读取现有 M65/M69/M72/M74 状态和账本，不做新求解、不改生产默认值。

用户此前已明确授权：针对当前 SimpleMOS 验证范围，向 Sentaurus 主机上传相关 TDR/deck、执行求解并取回结果。M78 建议先不使用该授权；若后续任务改变 deck、扩大器件矩阵或涉及破坏性远端操作，应重新核对范围。

## 2. 算例与术语

Workbench 参数 `NWell` 在这个 deck 中实际控制：

```text
init concentration=<NWell> field=Boron
```

所以本文的“高 NWell”是全局 p 型 body/substrate 背景硼由 `1e17` 提高到 `2e17 cm^-3`，不是传统掩膜 n-well 注入。

固定 `Lg=0.25 um`，器件矩阵为：

| 低/高配对 | GOxTime (min) | LDD Dose (cm^-2) | 背景 Boron (cm^-3) |
|---|---:|---:|---:|
| n17 / n21 | 10 | 1e14 | 1e17 / 2e17 |
| n18 / n22 | 10 | 2e14 | 1e17 / 2e17 |
| n19 / n23 | 15 | 1e14 | 1e17 / 2e17 |
| n20 / n24 | 15 | 2e14 | 1e17 / 2e17 |

每个器件包含 `Vd=0.05 V` 与 `1.0 V` 两条 Id–Vg 曲线；每条曲线 `Vg=0.00...2.50 V`、步长 `0.05 V`，共 16 条曲线、816 个点。

## 3. 当前总体结论

原始最大差异包含两个互相独立的成分，后续不得再混为一个根因：

### 3.1 深关断端口 burst：已由 M60 闭环

- 原 M46 最大差异位于 n23、`Vd=0.05 V, Vg=0.05 V`：`0.109418681 dex`，Vela 电流高 `28.653%`。
- M48–M59 把该尖峰定位到 Sentaurus 默认端口算法的 substrate 电子电流分配，而不是 SRH、势垒、自洽浓度或 Vela SG 核。
- M60 只收紧 Sentaurus 收敛策略后，目标 `default-minus-Direct` substrate 电子电流由 `4.6947455e-17` 降到 `-5.3589213e-19 A/um`，削减 `98.8585%`。
- 原 15 个 burst 标志点全部消失；目标默认 Id 与 Vela 的差降为 `0.024398 dex`。
- 因此该 burst 是亚 fA 尺度默认端口残差可见度问题，不应再开启通用接触提取、SG、HFS 或准费米排查。

### 3.2 平滑 NWell 差异：尚未完全闭环

- 去除 burst 后，8 个 NWell/漏压配对仍有稳定的 `0.021–0.026 dex` 平滑最大误差增幅，集中在 `Vg=0.55–0.90 V`。
- M63 证明它主要表现为阈值样水平平移：曲线误差能量解释率中位 `0.997962`，配对增幅解释率中位 `0.973200`，等效栅压平移 `2.709–5.673 mV`。
- M64 证明 BGN 不是其来源；BGN 实际在抵消更大的 no-BGN 失配。
- M65 发现基础本征浓度约定是重要成分：Sentaurus no-BGN `ni=1.075003848884424e10 cm^-3`，Vela 原值 `1.463891495876762e10 cm^-3`。匹配 `ni` 后，8 对误差增长中位值由 `0.111220` 降到 `0.037029 dex`，闭合 `68.05%`。
- 剩余 `0.037029 dex` 与 gate 阶段形成的静电/准费米势垒响应有关，但现有归因还没有支持生产默认修改。

## 4. M64–M77 因果链

| 任务 | 结论 | 约束/含义 |
|---|---|---|
| M64 | `bgn_independent_threshold_shift_dominant` | BGN 抑制而非制造 NWell 平滑增幅，不再把 BGN 当主根因调参 |
| M65 | `base_intrinsic_density_convention_material_but_not_dominant` | 匹配 no-BGN `ni` 闭合配对增幅中位 68.05%，生产材料未改 |
| M66 | `matched_ni_full_curve_material_but_not_complete` | 16×51 点完整复核；高 NWell 最大误差 `0.064469 dex`（16.00%） |
| M67 | `residual_qf_barrier_material_but_not_complete` | 剩余配对增幅的势垒代理中位闭合 68.94%；电势分量大于准费米分量 |
| M68 | `dos_ratio_inactive_in_frozen_boltzmann_path` | 当前 classical no-BGN Vela 路径不响应 Nc/Nv 比值，不再用 DOS 参数拟合 |
| M69 | `gate_stage_dominant` | 剩余势垒配对增量约 75.20% 在 gate 阶段形成 |
| M70 | `transport_support_proxy_not_material` | `n*mu*|grad(phin)|` 源侧代理方向不符，不支持迁移率/输运支撑主因 |
| M71 | `gate_surface_partition_material_but_not_dominant_charge_proxy_unqualified` | 表面势代理中位闭合 36.97%；旧 P1 栅电荷代理不合格 |
| M72 | `native_charge_qualified_full_profile_material_but_not_dominant` | 原生栅反力相对 Sentaurus ContactCharge 中位误差 1.853%；完整界面代理闭合 32.89% |
| M73 | `material_partition_material_but_not_dominant` | 固定状态材料分区 Poisson 候选同号 8/8，中位闭合 73.57%，只是线性观察器 |
| M74 | `electron_poisson_volume_material_but_not_dominant` | 电子电荷体积独立自洽 A/B：8/8 配对改善，中位闭合 41.05% |
| M75 | `hole_poisson_volume_not_material` | 空穴电荷体积响应近乎为零 |
| M76 | `dopant_poisson_volume_not_material` | 掺杂电荷体积使配对增长恶化到中位 `0.051796 dex` |
| M77 | `combined_poisson_volume_not_material` | 三项组合只有 3/8 改善；电子改善被掺杂项抵消，非线性交互很小 |

这里最重要的矛盾是：M73 的固定状态线性材料分区代理很强，但真正自洽组合 M77 不成立。M74 的电子独立干预能够改善配对增长，却使高 NWell 的绝对误差略增。这提示下一步应定位 **电子体积干预在哪个阶段、哪些 Si/SiO2 界面节点产生差分响应**，而不是继续增加全局组合开关。

## 5. 最新 post-main M77 结果

合并本地 main 后重新编译，并使用独立合同重新执行全部 16 个 M77 候选工作流、48 个阶段。没有覆盖原 M77 证据，也没有新增 Sentaurus 求解。

验收：

- 16/16 候选工作流收敛；
- 8/8 配对完整；
- 基线哈希错误 0；
- 非单变量配置差异 0；
- 所有验收项通过；
- post-main `--verify` 通过。

关键指标：

| 指标 | 数值 |
|---|---:|
| 改善配对数 | 3/8 |
| 基线配对增长中位值 | `0.037028975 dex` |
| 组合候选配对增长中位值 | `0.034937626 dex` |
| 自洽配对降幅中位值 | `-0.000271306 dex` |
| 高 NWell 基线绝对误差中位/最大 | `0.058271337 / 0.064469480 dex` |
| 高 NWell 候选绝对误差中位/最大 | `0.081375077 / 0.087196723 dex` |
| 独立项降幅和中位值 | `-0.000346398 dex` |
| 非线性交互中位/最大绝对值 | `0.000080242 / 0.000116601 dex` |

合并前后以下四份核心 CSV 的 SHA-256 字节级一致：case、pair、curve-response、charge-volume-decomposition。最新 main 没有改变 M77 数值结果。

## 6. 分支、提交与合并状态

关键提交：

```text
89533eb validation: decompose SimpleMOS Poisson charge volumes
3a46f85 Merge main into codex/simplemos-sdevice-validation
72551ec validation: requalify SimpleMOS M77 after main merge
```

本地 `main` 为 `8e99e57`。在 `72551ec` 时，当前分支相对 main 为落后 0、领先 38 个提交，工作树干净。

合并时 `main` 删除了旧的 examples 体系和 `tests/test_mos_mixed_material.cpp`。冲突解决策略是保留 main 的删除，同时把 M74–M77 两个核心装配测试迁移为独立目标：

- `tests/test_poisson_charge_volume.cpp`
- CMake 目标 `test_poisson_charge_volume`

不要恢复已删除的 examples 测试文件。

本机 MSYS Git 在 `core.fsmonitor=true` 时曾把清洁工作树误报为全量 dirty。需要写 Git 索引时，优先使用：

```powershell
& 'C:\Users\qzw\.cache\codex-runtimes\codex-primary-runtime\dependencies\native\git\cmd\git.exe' -c core.fsmonitor=false status --short --branch
```

不要用 `git reset --hard` 或批量 checkout 处理该误报。

## 7. 已闭环且不得重复的方向

除非出现新的独立反证，不要重新执行或建议：

- 全局 HFS/RefDens/速度饱和参数扫描；
- SG 电流核和通用 Vela 接触电流提取排查；
- 准费米参考、冻结或增量打包排查；
- 用关闭 BGN、关闭 SRH、Nc/Nv 或迁移率拟合修正生产默认值；
- 用 Sentaurus `DirectCurrent` 全局替换默认 drain Id；
- 把 SiO2 当导电区，或为标准 Si/SiO2 静电界面强制加入双节点/热发射边界；
- 把导出节点电流密度线积分当成 Sentaurus 默认 DopingWell 内部面积分；
- 继续把深关断 burst 与 `Vg=0.55–0.90 V` 平滑差异混合分析；
- 把 M73 固定状态线性代理直接当作已验证的自洽生产修正；
- 启用 M74–M77 诊断开关作为生产默认值。

## 8. 建议下一任务：M78 合同草案

建议标题：**M78 电子 Poisson 电荷体积自洽响应的阶段与空间归因**。

### 目标

解释 M74 为什么能使 8/8 NWell 配对改善、中位闭合 41.05%，但同时使高 NWell 绝对误差略增；确定响应是否主要在 gate 阶段、Si/SiO2 混合节点以及高/低 NWell 的差分电势中形成。

### 冻结输入

- M65 matched-ni/no-BGN 基线的 16 个三阶段工作流；
- M69 阶段定义和势垒账本；
- M72 原生栅反力与完整界面观察器；
- M73 材料分区 Poisson 线性账本；
- M74 电子-only 的 48 个现有阶段状态与 16 条曲线；
- M60 收紧收敛 Sentaurus 曲线作为参考，禁止新增端口算法变量。

本机现有 `build-release/m74_electron_poisson_charge_volume`，包含 48 份 `state.csv` 和 48 份 `curve.csv`。M78 第一阶段应只读这些文件；若哈希门禁失败则停止，不要静默重算。

### 唯一比较轴

`baseline` 对 `poisson_electron_transport_node_volume=true`。保持 hole、dopant、连续性/复合体积、介电边、输运体积、网格、接触、BGN、ni、偏压路径和收敛设置不变。

### 必须输出

1. 16 个工况在 equilibrium、drain、gate 三阶段的 `delta psi/phin/phip/log10(n)/log10(p)` 范数与 NWell 配对差分；
2. 完整 Si/SiO2 界面的节点级 `delta psi`、电子 Poisson 电荷反力和栅反力响应；
3. source/channel/drain/substrate 邻域的空间分区账本，不能只给全局最大值；
4. M73 线性预测与 M74 实际自洽 `delta psi` 的节点级和配对级误差；
5. 对“配对增长改善”与“高 NWell 绝对误差恶化”分别记账；
6. 控制工况至少包含 n19/n23 与高 LDD 反例 n20/n24，并覆盖两个 Vd。

### 建议判别门槛

- 16/16 工况、48/48 阶段状态哈希身份通过；
- 无新 Sentaurus/Vela 求解；
- 任何状态差必须来自 M74 单轴，而非文件拼接或偏压点错配；
- gate 阶段份额、界面节点份额和 M73→M74 预测误差都必须按 8 个配对报告中位/P95/最大值；
- 若没有至少 6/8 配对同号且空间支撑稳定，不进入生产策略讨论；
- 若 M74 响应主要由少数混合节点承担，再单独冻结 M79 局部支撑 A/B；否则停止该路线。

M78 不应新增全局调参，也不应把 M74 的 41% 配对闭合包装成绝对误差改善。

## 9. 证据入口

优先阅读：

- `docs/validation/simplemos_m60_tight_convergence_port_burst_2026-09-02.md`
- `docs/validation/simplemos_m65_nobgn_intrinsic_density_attribution_2026-09-03.md`
- `docs/validation/simplemos_m66_matched_ni_full_curve_2026-09-03.md`
- `docs/validation/simplemos_m69_stage_residual_localization_2026-09-03.md`
- `docs/validation/simplemos_m72_native_gate_reaction_2026-09-03.md`
- `docs/validation/simplemos_m73_material_partitioned_poisson_ledger_2026-09-03.md`
- `docs/validation/simplemos_m74_electron_poisson_charge_volume_2026-09-03.md`
- `docs/validation/simplemos_m77_combined_poisson_charge_volume_post_main_2026-09-05.md`

机器报告：

- `reference_tcad/simplemos_sentaurus2022/tight_convergence_port_burst/m60_tight_convergence_port_burst_report.json`
- `reference_tcad/simplemos_sentaurus2022/nobgn_intrinsic_density_attribution/m65_nobgn_intrinsic_density_attribution_report.json`
- `reference_tcad/simplemos_sentaurus2022/native_gate_reaction/m72_native_gate_reaction_report.json`
- `reference_tcad/simplemos_sentaurus2022/material_partitioned_poisson_ledger/m73_material_partitioned_poisson_ledger_report.json`
- `reference_tcad/simplemos_sentaurus2022/electron_poisson_charge_volume/m74_electron_poisson_charge_volume_report.json`
- `reference_tcad/simplemos_sentaurus2022/combined_poisson_charge_volume_post_main/m77_combined_poisson_charge_volume_post_main_report.json`

最新冻结入口：

- `reference_tcad/simplemos_sentaurus2022/simplemos_m77_combined_poisson_charge_volume_post_main_contract_v2.json`
- `reference_tcad/simplemos_sentaurus2022/simplemos_m77_combined_poisson_charge_volume_post_main_contract_freeze_v2.json`
- `reference_tcad/simplemos_sentaurus2022/simplemos_m77_combined_poisson_charge_volume_post_main_evidence.json`

## 10. 本地复核命令

在工作树根目录运行：

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"

cmake --build build-release --parallel 4

build-release\test_poisson_charge_volume.exe
build-release\test_mobility.exe
build-release\test_newton_solver.exe "NewtonSolver: parses signed AverageBox transport node volume"

python -m unittest `
  tests.regression.test_simplemos_m73_material_partitioned_poisson_ledger `
  tests.regression.test_simplemos_m74_electron_poisson_charge_volume `
  tests.regression.test_simplemos_m75_m77_poisson_charge_volume_decomposition

python scripts\run_simplemos_m77_post_main_requalification.py --verify
```

2026-09-05 已验证：Release 全量编译成功；Poisson 电荷体积测试 2/2、9 项断言通过；mobility 28/28、138 项断言通过；Newton 专项 19 项断言通过；M73–M77 Python 回归 6/6 通过；post-main M77 verify 通过。

如仅需复核冻结结果，不要运行 `--run`。完整 M77 自洽复跑会执行 16×3 个阶段，耗时显著。

## 11. 可直接交给新聊天的请求

> 请在 `D:\code-repo\vela-tcad\.worktrees\simplemos-sdevice-validation` 的现有分支继续 SimpleMOS 验证。先完整阅读 `docs/validation/simplemos_m45_m77_current_validation_handoff_2026-09-05.md`，核对 HEAD、工作树和 post-main M77 冻结证据。不要重复 HFS、SG、接触提取、准费米打包、BGN、Nc/Nv 或全局 Poisson 电荷体积组合调参。请先制定并冻结 M78 合同，再只读比较 M65 与 M74 的 48 个三阶段状态，定位电子 Poisson 电荷体积响应在 equilibrium/drain/gate 阶段和 Si/SiO2 界面节点上的形成位置，同时分别解释配对增长改善与高 NWell 绝对误差恶化。若现有状态哈希不通过，请停止并报告，不要自动重算。
