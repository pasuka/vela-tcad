# SLOT：Auger 与材料/多晶硅的独立研制

日期：2026-09-21。源码基线 `76817440`，分支
`codex/slot-ldmos-bvds-next`，本报告对应其上的工作区改动。

本轮分开检查 Auger 密度增强、Si/SiO2 能带参数、多晶硅半导体及接触恢复。
Auger 增强在 2、5、10 V 的成对修正后没有改变输出精度内的漏电流；恢复
多晶硅后，合格 2 V 电流增加 3.407712%，相对原始 Sentaurus 的偏差由
约 −2.748866% 缩小到 +0.565172%。材料能量参考仍有约 7.93 mV 差异，
此结果不代表全部材料物理或高压/BVDS 已对齐。

## 共同条件与证据

- Windows UCRT64 C++20 Release；无雪崩、300 K、Boltzmann、OldSlotboom，
  原始各向同性元素 IALMob。未改变全局模型默认值。
- 输入网格仍为 `build/slot_ldmos_noavalanche_1p25_resume02/mesh.json`：
  13,881 节点、27,275 三角形；本地边耦合及 mixed-Voronoi 体积保持不变。
  原生导出的节点坐标和单元拓扑经过逐项检查。本轮未替换成 AverageBox。
- `unit_scaling` 使用 cm 基物理参数；网格坐标 µm，状态 CSV 中浓度为
  m⁻³，报告电流为 A/µm。虽然键名包含 `_m6_per_s` / `_m3`，该配置中的
  Auger 系数及参考密度遵循 active-unit 约定，不能再次乘入 SI 转换。
- Auger/能带首轮使用 SparseLU；多晶硅成功点及最终三组 2 V 复验使用
  UMFPACK。更换初始化、后端和全局化并非各自独立的性能实验。
- 验收沿用原配置：`abstol=1e-9`、`stall_residual_floor=1e-8`、
  `eps_row=0.001`；全局连续性闭合强制开启，容差 0.01。
  独立审计检查所有有效载流子行，包含没有 SRH/Auger 源项的多晶硅。
- 完整证据汇总：
  [independent_controls_summary.json](../../build/slot_physics_controls_20260921/independent_controls_summary.json)。
  以下 `build/` 数据为本地忽略的实验产物，不随源码提交。

原始参考参数来自
`build/slot_deck_comparison_20260921/original_reference/pp2_des.par`，SHA256
`8f3c94fbca5ac9cd462dadee00549e032466a72bfb02fa092e4d3bac9d2e2919`。
原始 Workbench 目录当前不可用；“原始”沿用已冻结并执行的参考 deck 定义，
不声称与最早 Workbench 文本逐字节相同。
背景见[脚本对照记录](slot_ldmos_sentaurus_vela_deck_comparison_2026-09-21.md)。

## A：Auger 密度增强

原生 300 K 参数电子/空穴分别为：

| 参数 | 电子 | 空穴 |
| --- | ---: | ---: |
| C(300 K)=A+B+C，cm⁶/s | 2.9e-31 | 1.028e-31 |
| H | 3.46667 | 8.25688 |
| N0，cm⁻³ | 1e18 | 1e18 |

本地已有密度增强函数及其完整 n/p 导数，因此本轮通过独立配置启用
`auger_density_dependence`，继续保持 `auger_with_generation:false`。
没有把普通 Auger 改成反向产生模型，也没有同时改变 ni、接触或材料。

首轮固定状态验证均以 0 次 Newton 更新通过。为避免把“已满足停止条件”
误当成“增强项没有作用”，又在每个偏压分别对基线和增强配置强制执行一次
载流子 Newton 修正，再执行原门槛下的完整耦合验证和全有效行审计。

| Vd | 普通 Auger，A/µm | 加入 H/N0，A/µm | 成对结果 |
| --- | ---: | ---: | --- |
| 2 V | 9.0118343756412291e-16 | 9.0118343756412291e-16 | 输出双精度值相同 |
| 5 V | 1.6839859459084896e-15 | 1.6839859459084896e-15 | 输出双精度值相同 |
| 10 V | 2.9214614951097886e-15 | 2.9214614951097886e-15 | 输出双精度值相同 |

反偏耗尽区的普通 Auger 产生支路关闭；剩余正向复合主要处于近热平衡区。
2/5/10 V 的最大正准费米分裂仅为约 1.45e-13 / 2.85e-13 / 5.57e-13 V。
2 V 两配置的节点 SRH+Auger 项最大变化约 1.61e-34（内部连续性单位）。
这说明该项不能解释当前低压百分数级电流差；并不说明 Auger 增强函数恒为零，
也不外推到高注入或击穿状态。

证据：`forced_corrector_results.json`、`auger_active_region_audit.json`，均位于
`build/slot_physics_controls_20260921/`。

## B1：Si/SiO2 能带参数

对 2 V 逐项运行 Si Eg、Si Chi、SiO2 Chi，以及三项同时匹配的对照：

| 项目 | 基线 | 原生目标 | 当前同偏压影响 |
| --- | ---: | ---: | --- |
| Si Eg，eV | 1.12 | 1.12416 | 无可见电流变化 |
| Si Chi，eV | 4.05 | 4.0727 | 无可见电流变化 |
| SiO2 Chi，eV | 0.95 | 0.9 | 无可见电流变化 |

四组均通过原门槛，全有效行审计通过，电流为
`9.0118343756415979e-16 A/µm`，0 次 Newton 更新。这里固定了已校准的
Si `ni=1.0749877484352535e10 cm^-3`，材料参考温度就是 300 K，接触为
ohmic 或指定固定 flatband 的 metal_gate。上述字段没有改变当前使用的
密度及接触方程；不能再用 Eg 差额另行推算 ni 的变化。

这些参数对照只证明当前路径的作用范围。它们没有实现不同材料之间的
统一能量参考，也没有检验 Schottky/work-function 或变温边界。
证据：`local_results.json` 和 `noav_2v/*/config.json`。

## B2：多晶硅与接触恢复

原生导出 `build/slot_element_ialmob_20260921/native_elements/` 表明：
PolySilicon 区域 3 包含 350 节点、549 单元，具有非零施主/受主及电子/空穴。
该区域使用常数迁移率，未启用 Silicon 专属的 SRH/Auger/雪崩配置。
在虚拟机另用 `sdevice -P:PolySilicon` 导出材料参数，保存在
`build/slot_physics_controls_20260921/poly_native.par`，SHA256
`7fedee3835cf715cc8c5a9db1ecf15368bbc4251f44983ddbca8de1ac33d3d55`。

| 项目 | 恢复值及来源 |
| --- | --- |
| μn / μp | 原生全部单元一致：1417 / 470.5 cm²/(V·s) |
| ni(300 K) | 原生 n/p/BGN/准费米状态拟合：1.4638805412559153e10 cm⁻³ |
| 拟合热电压 | 0.025851995266484913 V；记录原生值，不替换 Vela 物理常数 |
| 拟合残差 | 最大对数误差 3.55e-15 |
| Eg(300 K) | 原生 PolySilicon 为 1.1082092307692308 eV，不能复制 Si Eg |
| 不含半 BGN 的 Chi | 4.0727403846153845 eV |
| BGN | OldSlotboom；Ebgn=0.009 eV、Nref=1e17 cm⁻³、C=0.5 |
| 掺杂 | 逐节点恢复原生 DonorConcentration / AcceptorConcentration |
| gate / SLOT | 改为 ohmic；电压仍为 0 V |

ni 拟合采用
`log(sqrt(n*p)) = log(ni0) + (BGN + phip - phin)/(2*Vt)`。
局部材料文件还提供 Nc/Nv 占位值以满足材料定义；本轮为 Boltzmann，
不将这些 DOS 值当作原生 Fermi/DOS 对齐证据。

实现包含两项显式、可选的能力：

1. `mobility.ialmob.constant_material_mobility_cm2_per_Vs`：让 Si 保留实时
   元素 IALMob，同时孤立多晶硅使用常数迁移率。常数迁移率自身导数为零，
   SG 通量对 ψ/φn/φp 的导数仍参与完整 Jacobian。名字与 Si 共边等不支持
   的配置会报错。
2. `solver.recombination_materials:["Si"]`：一致约束 Newton/Gummel 的
   SRH/Auger 残差、Jacobian、线性化与诊断积分。默认空列表保留旧行为；
   混合半导体共享节点且源项归属不同的情况拒绝执行，不隐式近似其体积。

**求解与独立性。** 先做零掺杂半导体多晶硅，再恢复原生掺杂。
`poly_intrinsic` 在 12 轮修正后仍未通过，不能作为可比较的电流点；
`poly_doped` 从旧 Si 基线启动、使用 SparseLU 的首轮完整求解超时 180 s。
随后保持恢复后的物理配置，使用原生 2 V 完整保存态初始化、Poisson 预热、
UMFPACK 和 merit 线搜索；12 轮载流子修正、累计 7 次完整耦合 Newton 更新
后通过原门槛。这里隔离的是“多晶硅半导体/掺杂/接触”这一物理组合，
尚不能把其变化分别归因于 ni、掺杂或 ohmic 中的单一因素。

| 2 V 配置 | Id，A/µm | 相对原始 Sentaurus |
| --- | ---: | ---: |
| 原普通 Auger、静电多晶硅基线 | 9.0118343756412291e-16 | −2.748866% |
| 恢复掺杂多晶硅、Si 专属 SRH/Auger | 9.3189317483778206e-16 | +0.565172% |
| 原始 Sentaurus 参考 | 9.26655970826958e-16 | 参照 |

最终 Release 从合格状态重放，得到相同电流。此复验为 0 次更新接受已解状态，
并非从 0 V 的独立连续扫描。Poisson 残差范数 `2.2443436460188418e-9`，
接受原因 `initial_stall_residual_floor`，符合原有 `1e-8` 地板；不能表述为
达到 `1e-9` 绝对阈值。独立全有效行审计：每种载流子 12,634 行，比基线
新增 328 行；无超标行，电子/空穴最大比约 `7.996e-6 / 8.761e-12`。
全局电子/空穴闭合比约 `4.586e-16 / 3.057e-16`。

证据：`poly_controls/poly_doped_native_seed/`、`final_release_revalidation/`。
成功初始化 CSV 的 SHA256：
`2b44a3aee7253614f78bd8bb392afde4480cdd8639f649e897b878c32c5eff8a`。

**剩余材料差异。** 在原生状态比较
`psi - (Vt_local*log(n/ni0) - BGN/2 + phin)`，Si 中位数约 `5.10e-8 V`，
PolySilicon 中位数约 `7.92781756e-3 V`。这表明目前统一的局部密度势参考
尚未保留原生材料间能量偏移。本轮没有用接触电压补偿把这项差异拟合掉。
证据为 `intrinsic_energy_reference_audit.json`。剩余 +0.565% 不能全归因于
该偏移；需进一步隔离它与控制体积、接口电场及源项离散的影响。

## 导数与回归

- Release 完整构建通过；879/879 CTest 通过，用时 135.14 s。
- 新增两项 C++ 用例：材料源项配置解析；Si/孤立多晶硅在 SI 与 TCAD 两种
  单位下的电流、常数迁移率、源项掩码、Newton 全列差分及 Gummel 线性化。
  合计 72 断言通过，差分相对容差 2e-6。
- 现有 `[auger]` 复合模型测试 7 项、102 断言通过；完整套件也覆盖已有
  密度增强与耦合 Jacobian 测试。
- 配置生成脚本新增 4 项回归：独立变量、输入不可覆写、单位/温度/原生
  字段校验、多晶硅恢复时保持 Si 输入。脚本最后修改后，注册的
  `slot_ldmos_ialmob_tools` 组重新通过。

真实 2 V 网格抽查 6 个 Si 节点的全部 ψ/φn/φp 列，多晶硅配置额外抽查
3 个 poly 节点。中心差分步长 1e-6；直接源项组装避免从大矩阵相减：

| 配置 | Poisson 相对矩阵差 | 输运相对矩阵差 | SRH/Auger 相对矩阵差 |
| --- | ---: | ---: | ---: |
| Auger H/N0 | 1.282e-9 | 6.052e-9 | 1.703e-9 |
| 恢复多晶硅 | 1.328e-9 | 5.502e-9 | 8.402e-7 |

微小 ψ 输运列的单独相对差约 7.6–7.9e-5，但绝对差范数仅约 1.8–1.9e-16。
多晶硅配置源项 φn 列相对差约 2.20e-6，绝对差范数 1.66e-23；
原始逐块数字保存在 `jacobian_audits/*/blocks.csv`。抽查不能替代全域逐列
验证，真实状态微小支路也不适合仅按相对误差判断漏导数。

最终 Release runner SHA256：
`1122317ef913b0870a289b4b30bee8fa1e33fa482b25883e19ad0e306308cb99`。

## 复用方式及边界

从本工作树根目录运行，保持 UCRT64 DLL 在 PATH：

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
$env:VELA_LINEAR_SOLVER = "umfpack"
cmake --build --preset windows-ucrt64-release --parallel 4
python scripts/prepare_slot_ldmos_physics_controls.py --help
```

脚本 `--base-config` 指向合格普通 Auger 配置，`--initial-state` 指向对应
状态；`--parameters` 指向原始 pp2_des.par，生成六个 Auger/能带对照；
改用 `--poly-fields` 指向原生导出 fields 目录，生成两种多晶硅对照。
`--output` 必须是新目录。准备器检查坐标/拓扑/单位并记录源文件哈希。
最终脚本已用真实数据重新生成
`build/slot_physics_controls_20260921/poly_controls_final_bundle/`，19 个输入
依赖均记录哈希。默认 poly seed 保留输入中的 Si 状态，并只初始化多晶硅；
它不承诺仅靠该 seed 即能收敛。本次成功路径的原生完整 seed、配置及完整
续算命令保存在 `poly_controls/poly_doped_native_seed/invocation.json`。

复验已有合格点可读取 `final_release_revalidation/*/config.json`，沿用
`run_slot_ldmos_noavalanche_continuation.py` 的 `--initial-bias 2 --targets 2
--max-attempts 1 --audit-local-continuity`。所有新尝试应另建输出目录。

本轮完成了两类独立配置、必要的多材料输运/源项实现及首轮量化验收。
新多晶硅配置尚未验证 5/10 V、完整 0→2 V 扫描、自洽雪崩或外电阻支路；
旧静电多晶硅配置的高压资格不转移到该配置。Auger 对照仍保持独立，
未将其与多晶硅合并后直接宣称最终原始模型已对齐。
