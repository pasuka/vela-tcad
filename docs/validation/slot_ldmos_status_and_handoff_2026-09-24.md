# SLOT LDMOS：原始功能对照、验证状态与下一阶段交接

本日后续执行已开始：三个 SLOT 驱动完成 HDF5 适配并修复 probe 候选
输出遗漏；新增原生 2.025 V 参考。下文是执行前状态快照，最新进度及
证据见 [HDF5 续算执行记录](slot_ldmos_hdf5_continuation_2026-09-24.md)。

核对日期：2026-09-24。源码基线 `c95f8b0d`，分支
`codex/slot-ldmos-bvds-next`，工作目录
`D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next`。
本次是文档梳理，没有启动新的求解或重新运行回归。
本机检查未发现 vela/sdevice 求解进程，开发代理列表无运行子代理；
本次没有检查虚拟机进程。“进行中”指已经开展但尚未验收的工作，不表示后台在算。

## 参考范围和当前结论

“原始 Sentaurus”指已冻结并实际执行的参考，而非最早 Workbench 文件逐字节复原：

- `build/slot_deck_comparison_20260921/original_reference/bvds_external_resistor_final.cmd`
- 同目录 `pp2_des.par`，本次重算 SHA256：
  `8f3c94fbca5ac9cd462dadee00549e032466a72bfb02fa092e4d3bac9d2e2919`。
- 同目录 `bvds_external_resistor_final.log_des.log`。

当前 Vela 以 `build/slot_main_integration_20260923/verify_2p125v/config.json`
和其材料/掺杂/几何文件为准。不能继续把 9 月 21 日最早脚本对照表里的
“Auger 未启用、多晶硅仅静电、gate/SLOT 为 metal_gate”当作现状。

当前主要功能已实现，验收仍分三层：参数/局部核、完整器件检查点、连续 BVDS。
新源 `ialmob_element_density_grad_qf` 的最后合格内部偏压为 **2.125 V**；
2 V 相对原生自洽雪崩参考电流差 **−0.0410440%**。2.025～2.125 V 的五个新点
已通过本地资格，但没有这些点的新原生自洽电流对照，不能继承 2 V 的误差值。
10 V 候选已实际计算但仍不合格。严格原始模型的 BVDS/击穿后支路未完成。

旧投影 SG/GSS-Laux 配置的内部 37 V、外部 37.1 V 已通过，属于另一物理配置。
它证明旧配置下求解和外电路流程的能力，不是当前新源的高压资格。

## 原始功能参数与当前验证矩阵

表中“已完成”仅适用于所列温度、模型、网格和检查范围。

| 功能 | 原始 Sentaurus 设置 | 当前 Vela 与证据 | 状态及剩余边界 |
| --- | --- | --- | --- |
| 仿真类型 | 2D 等温稳态 Poisson/Electron/Hole，关态漏极 Quasistationary | `dc_sweep/bv_reverse`，全耦合 DD | 基础功能已完成；完整原始扫描未完成 |
| 几何与掺杂 | `n1_fps.tdr` | 13,881 节点、27,275 三角形，24,272 个 Si 三角形；坐标/拓扑核对，多晶硅掺杂恢复 | 导入与当前网格核对已完成；不包括 SProcess 工艺重实现 |
| 电极 | source/gate/SLOT/substrate=0 V；drain 经外电阻升压；日志五者均 Ohmic | 五电极均 ohmic；新源当前验收用内部偏压，外电阻关闭 | 边界恢复及混合接触保护修正已完成；整体边界数值等价性不作全范围承诺 |
| 温度/统计 | 300 K，未启用 Fermi | 300 K Boltzmann | 当前范围已对齐；变温/Fermi 不属于原始激活目标 |
| BGN | OldSlotboom，Ebgn=0.009 eV、Nref=1e17 cm⁻³、C=0.5 | `old_slotboom`，相同主参数 | 当前等温路径已对照 |
| Si 材料 | Eg=1.12416 eV、Chi=4.0727 eV；原生密度/能带 | 已恢复相同 Eg/Chi 和原生 DOS；ni=1.0749877484352535e10 cm⁻³；不重复由 Eg 改写已校准 ni | 材料参考及密度关系已验证；不等于全部空间场逐点完全一致 |
| PolySilicon | 有掺杂、载流子；常数迁移率 | μn/μp=1417/470.5 cm²/(V·s)，ni=1.4638805412559153e10 cm⁻³，Eg=1.1082092307692308 eV，Chi=4.0727403846153845 eV；参考偏移 +7.927767678211338 mV | 恢复与独立对照已完成，密度/Poisson/接触链式关系一致；温度扩展未验收 |
| SiO2/Nitride | 绝缘介质；SiO2 Chi=0.9 eV | 当前 aligned_full 材料文件没有覆盖 SiO2 Chi，仍为库默认 0.95 eV；Nitride 静电材料 | 已做 SiO2 Chi 单因素控制，当时 2 V 电流无可见变化；当前配置并非全部材料参数逐项相同，Nitride 全参数等价未证明 |
| 界面迁移率 | Enormal(IALMob)，unnamed isotropic 参数，FullPhuMob | 元素 IALMob，电子/空穴及三组晶向使用原始 unnamed 参数；密度/法向场/高场依赖实时计算 | 参数、元素回放、完整链式导数已实现并验证；不是 Enhanced Lombardi 替代 |
| 高场迁移率 | HighFieldSaturation；β0=1.109/1.213，vsat=1.07e7/8.37e6 cm/s | 准费米梯度，vsat=1.07e5/8.37e4 m/s；RefDens=0 | 300 K 核及导数已验证；参数温度指数不代表变温已验收 |
| SRH | SRH(DopingDep)，τmax=1e−5/3e−6 s，τmin=0，Nref=1e16 cm⁻³，γ=1 | 相同参数；Si 专属；总杂质 Nd+Na | 当前 300 K 对齐；本地 temperature_dependence=true、原生未启用 TempDep，300 K 乘数为 1，开关仍不完全相同 |
| Auger | 普通 Auger，无 WithGeneration；C300=2.9e−31/1.028e−31 cm⁶/s，H=3.46667/8.25688，N0=1e18 cm⁻³ | 增强已开启，含系数对密度导数，WithGeneration=false，Si 专属 | 2/5/10 V 无雪崩独立成对修正及复验完成，电流在输出精度内相同；高注入影响未验证 |
| 复合积分 | Silicon 作用范围 | `material_local_mixed_voronoi`，排除共享界面上的介质面积，残差/Jacobian/诊断一致 | 已实现、差分和 2 V 误差预算完成；无雪崩误差由约 +0.56% 降至 −0.041006% |
| Poisson/输运几何 | 日志 `CVPL_AverageBoxMethod=TRUE` | 当前本地边耦合、mixed-Voronoi 节点体积；有原生耦合及局部电荷体积单因素控制 | 进行中：已有预算，尚未整体复现 AverageBox，不把小电流差当作场量完全一致 |
| 雪崩系数 | van Overstraeten，GradQuasiFermi；电子 a=7.03e5 cm⁻¹，b=1.231e6 V/cm；空穴低/高场 a=1.582e6/6.71e5，b=2.036e6/1.693e6；E0=4e5 V/cm，ℏω=0.063 eV | 单位换算后相同，self_consistent，无 RefDens 插值 | 系数及基础核已对齐 |
| 雪崩离散 | AvalDensGradQF、ElementVolumeAvalanche | 新核 `ialmob_element_density_grad_qf`；固定网格原生选点、接触/界面分支、实时迁移率和源映射 | 进行中：24,272 Si 元素选点已识别，10/37 V 冻结源回放通过；新源整器件仅到 2.125 V；通用网格原生选点算法未推出 |
| 完整 Jacobian | deck 的 Derivatives 被注释，日志明确 With avalanche derivatives | 迁移率端点/密度/场依赖和源项三顶点 ψ/Fn/Fp 九列链式导数；local_ad | 核、小夹具全列及实际网格选定方向差分已完成；不是实际 13,881 节点全部列的穷尽差分证明 |
| 数值精度/求解 | 日志 80-bit ExtendedPrecision，ParDiSo，5 线程；Digits=5、RHSmin=1e−10、ErrRef(e/h)=1e10，Iterations=25 | Release/UMFPACK；主要 binary64，分离 QF；abstol=1e−9、reltol=1e−8，全部有效行 eps=1e−3、全局闭合 0.01 | 进行中：不能直接比较不同归一化阈值；精度表示和 10 V 停滞问题未解决 |
| 初始化/扫描 | Poisson 初始解 Iterations=1000、LineSearchDamping=1e−8，再全耦合；Extrapolate | 有本地 Poisson、预热、步长回退、候选修正及完整复验 | 旧无雪崩有不依赖原生状态的分段 0→2 V；新源没有完整冷启动连续扫描 |
| 外电阻 | R=1e12，初始 0 V，目标外电源 2e5 V | R=1e12 Ω·µm；独立器件负载线续算已实现，旧 SG 外压 37.1 V 合格 | 旧模型阶段完成；新源外电阻未验收，C++ 联立电路全局化失败不能视为已由替代流程修复 |
| BVDS 停止 | drain AbsVal=1e−7；2D area factor=1，按 A/µm 统一比较 | 尚未在严格新源模型下达到该电流并验收完整支路 | 未执行完成：统一判据 BVDS、击穿后支路及内压/外压区分 |
| 输出/续算 | I–V .plt、空间 .tdr、10/20/30/40/50/60 V Plot 触发 | I–V CSV、空间/源项/迁移率审计、HDF5 检查点 | 输出与单点 HDF5 迁移已完成；三个新 SLOT 驱动的 HDF5 适配尚未执行 |

原生归一化扫描步长为 Initialstep=1e−6、MaxStep=0.05、MinStep=1e−7，
Increment=1.6、Decrement=2；不可直接把它们理解成 Vela 的电压步长。
2e5 V 是外电源的目标上限，不是器件已验证到 200 kV，也不是击穿电压。
热、量子、BTBT、隧穿、陷阱与表面 SRH 没有在原始 Physics/Solve 中激活；
仅出现在 Plot 或参数文件里的名称不应计为待补齐物理功能。

## 已验收、已开展未完成、尚未开展的任务

已验收阶段：原始参数整理；元素 IALMob 及完整导数；Auger 独立对照；
多晶硅/材料参考恢复；Si 复合体积修正；混合接触分类修正；新源固定网格
选择与冻结回放；新源 2→2.125 V 的 0.025 V 分段续算；旧 SG 两段内部
20/37 V 及 49 个新点审计、外压 37.1 V 自动流程；main 合并与单点 HDF5 兼容性。
最新实际 Release 回归为 **912/912**，见合并报告，不是本日重跑。

已开展但未完成：

1. 新源更高偏压自洽收敛。直接 2→2.125 V 大步曾失败，小步已通过；
   3 V 有中断尝试，10 V 有完整数值失败，不能把它们归入“从未执行”。
2. 10 V Poisson/接触阻塞：ψ 范数 1.2764099707460164e−9；接触多数 QF 差
   5.35167645226901e−11 V，高于 floor 路径的 5e−11 V 守卫；全有效行/闭合通过
   不能替代完整器件资格。修正前接触平均掺杂错误已解决，剩余跨 n/p 边不能随意排除。
3. 物理/离散/数值误差预算：2 V 已形成局部预算，高压旧 SG 约 +5.13% 的
   37 V 电流差尚未形成完整归因；新源冻结误差不等于端电流误差已解释。
4. 原生筛选函数有效下限差异已有 ppm 级诊断解释，未把拟合常数装入生产核；
   固定网格观测几何不等于任意新网格的通用选择算法。

尚未开展或没有完成该配置的整体验收：三个 SLOT 驱动 HDF5 迁移；新源
单一流程冷启动扫描；2.025～2.125 V 原生同偏压新增自洽对照；新源合格
3/5/10 V 及更高压曲线；新源外电阻/BVDS/击穿后曲线；SLOT worker 缓存收益测量。
“未合格”和“未开展”应在后续记录继续分开。

## 建议下一阶段顺序及验收条件

1. **先适配 HDF5**：`scripts/run_slot_ldmos_noavalanche_continuation.py`、
   `scripts/run_slot_ldmos_loadline_continuation.py`、`scripts/slot_poly_equilibrium.py`。
   覆盖预测、carrier-only 修正、多晶硅重启、确认态数值比较、断点恢复及身份校验。
   保留 split QF；测试有效/损坏/错网格状态及中断恢复；再复验 2/2.125 V。
   不恢复生产 CSV 后备读取；容器哈希与解码数值一致性分开。
2. **整理新源低压可重复流程**：先恢复 0.025 V 路径，检查原生同偏压点；
   然后从本地平衡态重复 0→2.125 V。候选修正后仍要完整耦合、每类 12,634
   个有效行及全局闭合验收。不要把导入原生保存态称为冷启动。
3. **解决 10 V 数值阻塞**：分开检验表示精度、Newton 更新与接触守卫；
   保留全部物理偏导，使用有界实验，不凭 ULP 敏感度宣称残差下界，也不放宽门槛。
4. **新源 3/5/10 V 及高压同模型比较**：电流、势、密度、迁移率、产生率和积分
   分别比较，完善几何/接触/停止误差预算；以新结果取代旧模型资格。
5. **恢复新源外电阻与 BVDS**：先通过低压负载线，再升压到统一 1e−7 A/µm
   判据及击穿后支路，分别报告内压/外压；电路复用状态不能绕过独立器件审计。

worker 结构缓存性能评估可在流程稳定后进行；缓存数值兼容性已初检，但未测
SLOT 跨请求命中收益，不将其他 LDMOS 模板的加速比转用到当前算例。

## 可供新对话直接读取的入口

优先阅读本文件、`docs/README.md`、`docs/state_archive.md` 以及：

- [main 合并及最新复验](slot_ldmos_main_integration_2026-09-23.md)。
- [混合接触修正与小步续算](slot_ldmos_contact_guard_and_escape_controls_2026-09-23.md)。
- [新源受控装配](slot_ldmos_observed_source_integration_2026-09-23.md)：其后段“仅 2 V”已由后续 2.125 V 结果更新。
- [材料参考及体积预算](slot_ldmos_energy_reference_and_source_volume_2026-09-21.md)。
- [原始参数明细](slot_ldmos_sentaurus_vela_deck_comparison_2026-09-21.md)：仅用其原生参数/来源，旧 Vela 状态以本文件覆盖。
- [外电阻旧 SG 流程验收](slot_ldmos_loadline_review_2026-09-23.md)。

可复验 HDF5：`build/slot_main_integration_20260923/2v.h5`、`2p125v.h5`。
最后合格点对应 `verify_2p125v/config.json`；结果 `compatibility.json`。
`10v_unqualified.h5` 是诊断候选，不能当合格续算锚点。
CSV 历史链的身份清单是 `build/slot_followup_20260923/qualified_small_step_manifest.json`；
所有 `build/` 状态/几何证据不随 Git 传输，新目录/新机器必须额外保留这些输入。
当前材料文件 `_comment` 仍有“PolySilicon electrostatic-only”的旧说明，
实际已含 ni/迁移率/DOS/参考偏移；应按字段和运行行为判定，不能只读注释。

运行使用 UCRT64 Release preset、PATH 包含 `D:/msys64/ucrt64/bin`，
器件验证指定 `VELA_LINEAR_SOLVER=umfpack`。不要用旧 exe 验证新接口。
Sentaurus 原连接地址 192.168.119.131，既有 root 密钥方式；本次未复查连通性。
GCC11 路径为 `/home/tcad/sentaurus_runs/vela_oracle_2022/vela_build_20260914/toolchain/opt/rh/devtoolset-11/root/usr`。
既有授权允许上传/下载对照文件；不要恢复旧定时暂停或监测任务，除非用户再次要求。

建议新对话开场内容：

> 请继续 SLOT LDMOS 验证，使用现有工作目录
> D:/code-repo/vela-tcad/.worktrees/slot-ldmos-bvds-next，分支 codex/slot-ldmos-bvds-next。
> 先阅读 docs/validation/slot_ldmos_status_and_handoff_2026-09-24.md，核对当前 Git、配置和证据。
> 本阶段先完成三个 SLOT 驱动的 HDF5 适配及 2/2.125 V 回归，再按文档顺序推进。
> 使用 Release/UMFPACK，保留完整迁移率/雪崩 Jacobian 和原验收门槛。
> 新源仅合格至 2.125 V，10 V 未合格；旧 SG 的 37/37.1 V 资格不可转用。

这段文字是供用户选择的新对话请求模板，本次仅完成梳理，不因此启动实现或新任务。
建议在当前合并验收完成的节点换新对话，旧对话保留为审计历史；继续使用本对话也可行。
理由是划清新的状态接口阶段，减少旧模型、旧阈值说明和已被更新结论的混淆，
并非声称长对话必然失效或新对话必然更准确。
OpenAI 的[长任务实践](https://developers.openai.com/blog/run-long-horizon-tasks-with-codex)
也强调把目标、验收和状态落入可持续读取的项目文件。
