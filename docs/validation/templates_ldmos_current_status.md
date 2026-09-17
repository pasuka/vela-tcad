# Templates/LDMOS 当前验证状态

更新：2026-09-17。当前工作分支 `codex/templates-ldmos-phase-a`。

## 当前总览（2026-09-17）

最新 [R9 UCRT64 gprof 复核](templates_ldmos_gprof_r9_2026-09-17.md)：Vg8 从零
漏压起前 8 个精确点全部通过；普通 Release 166.93 s（含初始化），漏压
137 次更新，装配/分解占漏压尝试墙钟 61.35%/23.30%。完整及仅漏压插桩轨迹
与对照逐值一致；IALMob 为主要可归因采样热点，另需关注跨点几何准备。
gprof 部分调用计数归属异常，使用内部计数；本轮未修改求解策略或推荐配置。

最新[带保护的局部标量 Newton 验证](templates_ldmos_neutral_root_newton_2026-09-16.md)：
Windows 8 个代表点全部非计时输出逐值一致，平均密度评估 73→16.39 次/根，
外层更新和装配次数不变；829/829 Release 回归通过。Linux 两轮 32 次运行
全部通过原门限，非计时结果与对照及跨轮逐值一致，密度评估减少 77.50%；
墙钟仅减少 1.34%/0.46%。候选默认关闭，不转移完整曲线资格。

本轮主体 Newton 与装配对照见[新执行记录](templates_ldmos_newton_cost_2026-09-16.md)：
Windows 代表点缓存保持逐值一致；局部 QF 裁剪 3/8 点失败，不晋级。
Linux 两轮代表点及完整双栅压对照已完成；缓存的全部非计时轨迹与 R9 逐值一致，
但墙钟 525.05/520.02→527.12/504.09 s，未证明稳定加速，不纳入推荐配置。
主体非 floor 更新仍为 343/297；新编译程序的 R9 全轨迹已与原冻结 R9 核对一致。

已将通过完整双栅压、真实恢复及两轮重复验收的接触一致性方案命名为
**R9 显式推荐配置**，见[生产复现指南](templates_ldmos_production_reproduction.md#recommended-explicit-r9-contact-consistency-profile)。
导出器兼容默认仍为 R7；选择 R9 必须显式给定 profile。资格对应既有冻结
Linux 程序和输入，新编译程序及后续缓存/更新候选须单独验证。

| 原任务 | 当前完成情况 | 仍有的边界或后续事项 |
|---|---|---|
| 局部物理量差异闭合 | 有限 hRecVelocity、Auger 生成语义、温度物性及原生 Poisson 几何已实现并对照；最新 D0 双栅压 62 点通过 2% 密度 RMS、30 meV 带边门限 | 属于工程门限内一致，局部场仍有余差；不宣称商业求解器内部算法完全相同 |
| 生产 DC 电热入口集成 | 显式 `electrothermal_dc_sweep` 已支持中性初始化、栅压预偏置、四方程扫描、失败恢复、检查点及电热输出；生产/独立入口一致性和真实恢复通过 | 普通 `dc_sweep` 仍为等温；接触修正已纳入显式推荐 R9，通用默认仍关闭 |
| 整体性能与重复性 | R7 和最新接触候选均完成同 VM 配对；候选两轮数值轨迹一致，更新 366/321→344/302，墙钟中位数约原生 1.018/1.063 倍 | CPU 仍为原生 1.625/1.763 倍；非 floor 更新 343/297 未减少，时间收益幅度不稳定 |
| 最终联合验收 | R7 已完成冻结配置联合验收；最新显式接触候选另完成前8点恢复、完整62点电学/热学/局部场及两轮逐组墙钟门限 | 资格随冻结程序和输入；接触方案已晋级显式推荐 R9，新缓存/局部限幅未纳入推荐 |

G3/D5/D4 已按各阶段冻结配置完成验收及共享路径复核；本次接触候选阶段未
重新运行这些等温长曲线。D3/D2/D1 的资格限原 Solve 实际激活的等效路径。
原脚本未激活的量子方程、Thermodynamic/Peltier/RecGenHeat 是独立扩展，不列作
当前原脚本生效功能对齐的欠项。继续优化应优先针对主体 Newton 更新与装配
成本；物理余差研究、更多网格/边界推广均需独立实验。

接触方案两轮资格见[完整曲线联合验收](templates_ldmos_contact_sweep_2026-09-16.md)，
最新单组成本对照见[Newton 与装配报告](templates_ldmos_newton_cost_2026-09-16.md)。
本次本地提交前，全目标 UCRT64 Release 构建成功，完整 CTest **828/828**
通过（113.46 s）；18 个修改/新增 Python 文件语法及相关文档链接检查通过。
日志位于 `reference_staging/templates_ldmos_contact_sweep_20260916/commit_review/`。
本次构建和测试不替换完整曲线所用冻结 Linux 程序的资格或计时。
以下保留各阶段执行经过；历史段落的后续建议不自动构成当前任务队列。

## Newton 诊断与接触候选执行记录

本轮[Newton 更新映射与阻尼验证](templates_ldmos_newton_update_execution_2026-09-15.md)已完成：
F0 在 Linux Release 重放 R7 两档 74 次尝试，状态、残差和原历史逐值一致。
Windows Release 八代表点中，R7 为 8/8 通过、101 更新；V1 为 3/8、258 更新，
V2 为 4/8、170 更新，均不晋级。两个 5.333 V 样本的 NLEQ_ERR 型对照没有恢复
密度候选；保守 Jacobian 刷新更新 30→33、装配 71→78、分解 30→29，无净收益。
F0 确认 floor 仍受原逐行/块门限阻塞，不启用 F2 R-B。
后续[载流子质量矩阵与投影伪瞬态验证](templates_ldmos_pseudo_transient_2026-09-15.md)
已实现独立开关、通过四变量偏导及实际 SG 扩散检查；两个难点的 12 组对照中，
三档 PTC+V1 均不通过。1× 消除初始非正密度目标，但真实残差方向斜率为正，
分别冻结迁移率/体复合导数后仍然如此。候选不晋级，未放宽终态门限，
完整曲线候选验收未启动，R7 仍为生产默认。
随后[伪瞬态接受准则验证](templates_ldmos_pseudo_acceptance_2026-09-15.md)完成：
实际后向 Euler 缺陷 G 的充分下降判据跨过初始拒步，但 G+SER 和 G+模型误差
控制均在两档各 60 更新后未通过稳态门限，真实 F 范数仍增大。8 个定向数值
用例、27 项电热回归通过；R7 两档仍各 15 更新，状态/残差/历史逐值一致。
候选不晋级。后续 [Poisson 一致初始化对照](templates_ldmos_poisson_initialization_2026-09-15.md)
已完成：两档各 5 更新闭合 Poisson，原边界恢复和 QF/T 保持检查通过。
G+模型误差控制走到近稳态，但两档逐行比仍为 3.75e-5/1.63e-5，超过 1e-8；
加准备共 56/65 更新，均不晋级。原种子 R7 各 15 更新通过；加准备后为
24 更新通过/16 更新失败，不接为默认。29 项电热回归通过。
后续[近稳态质量项、密度与 QF 表示隔离](templates_ldmos_near_steady_isolation_2026-09-15.md)
已完成 18＋6 组：关闭质量项仍未闭合；自身状态的局部参考重表示后，
两档均只需 1 次 R7 QF 更新通过原门限。密度映射多拒绝 8 次才回退到同一
QF 步，终态逐值相同。加回原前处理和伪瞬态总更新为 57/66，未证明整段
加速；30 项电热回归通过。
[独立自动近稳态切换](templates_ldmos_near_switch_2026-09-15.md)现已实现，并完成
原种子完整代表点轨迹对照：带 Poisson 准备时总更新 32/31，均通过原门限；
不准备时仍各 60 次失败、未触发切换。当前 R7 仍各 15 次通过，状态、残差和
历史逐值保持；9 个 C++ 用例/5003 断言和31项电热回归通过。候选不晋级，
生产默认和门限不变。
[R7 高压独立局部重表示对照](templates_ldmos_r7_local_rebase_2026-09-15.md)已完成：
Windows/Linux 各四点两轮配对，共32次均通过原门限；每轮总更新分别44→44、
45→45，装配85→88、87→88。Linux Vg4/30.667 V的14次更新和5次floor长尾
完全保留，无收益、不晋级。两平台各9个数值用例及32项电热回归通过。
后续[有限接触行舍入审计](templates_ldmos_contact_row_roundoff_2026-09-16.md)已完成：
六个 Linux 原轨迹状态逐值重放，高精度求和/SG/复合重算未消除失败。
定位到有限接触中性电势目标随一个温度浮点间隔跳变、ψ 留下 8–16 个浮点
间隔偏差，经接触通量放大。五个失败固定状态恢复该代数等式后，原逐行与
三块门限均通过，fn/fp/T 不变。两平台各9个数值用例及34项电热回归通过。
该审计阶段仅完成固定状态定位。随后[接触代数一致性候选](templates_ldmos_contact_consistency_2026-09-16.md)
已实现默认关闭的一次近稳态修复，失败则原样回退；Windows/Linux 各四个高压点
两轮配对共32次全部通过原门限、重复轨迹逐值一致。Linux Vg4/30.667 V
原种子完整轨迹14→9次更新，五次接触长尾被消除；两平台单轮总更新44→40、
45→40，装配85→84、87→82。Windows总墙钟中位数33.10→30.33 s，
Linux50.26→50.95 s且两轮优劣反转，尚不认领稳定时间加速。
两平台各9个定向数值用例和36项运行器回归通过，含温度链式关系及实际重装配
后拒绝回退。原R7基线状态/残差/历史保持，生产默认及门限不变。
[跨步与完整曲线验证](templates_ldmos_contact_sweep_2026-09-16.md)已完成：
初始化与栅压预偏置保护已实现，38项运行器回归通过；双栅压前8点及候选的
真实暂停/恢复均通过，状态与R7逐值一致，漏压更新保持152/137。
两轮完整双栅压R7/候选/原生共12组串行配对完成，原电学/热学、62点局部场、
逐轮1.5倍墙钟门限均通过；两轮非计时轨迹逐值一致，前8点真实暂停/恢复与
完整前缀一致。漏压更新366/321→344/302，floor更新23/24→1/5，主体更新
343/297不变。两轮墙钟中位数539.75/550.99→520.71/511.32 s，减少3.53%/7.20%；
相对原生墙钟1.018/1.063倍，CPU仍1.625/1.763倍。中断后已补齐未入账原生
样本，原样本及额外成本保留；逐轮时间收益波动较大，不承诺稳定加速比例。
本轮显式候选取得限定配置下的62点联合资格，仍默认关闭，生产R7默认不变。

外推 A–E 已按用户要求串行完成：A六组、B四组、C四组完成，无总成本收益。
D首轮六组及反向四组复测完成，局部拟合仅令Vg8的原失败尝试从10次更新
缩短为8次，仍需减步恢复；两档计时优劣均在复测中反转，不认领稳定加速。
E实现及Windows Release822/822、Linux Release七项定向检查通过；获源码上传
授权后完成两个停滞样本的四组VM对照。NGMRES未恢复成功，更新仍为12/10，
装配34/25→41/28；原失败状态、Newton历史及残差逐值不变，不晋级默认。
本轮状态见[外推执行记录](templates_ldmos_extrapolation_execution_2026-09-15.md)，R7仍为生产默认。
经用户评审收敛的下一轮方案见[Newton更新映射与阻尼方案](templates_ldmos_newton_update_plan_2026-09-15.md)：
用户已授权并完成 F0、F1 及 NLEQ_ERR/Jacobian 的阶段对照，候选未晋级；
生产默认、物理与门限保持不变，F2 R-B 接受规则变更需届时另行明确同意。
本页区分历史资格和最新版本验证，阶段报告保留原始结果。
每项曲线资格对应其报告中的冻结程序和输入清单，不随工作区重编译自动转移。

R7双栅压62点电学、热学、2%局部场及本轮同VM墙钟验收已通过，实际
暂停/恢复和生产/独立入口数值一致性也已通过。用户恢复任务后已完成受版本
管理的输入导出入口、运行说明和回归验证，已纳入本轮本地分支整理。
用户随后要求重复性验证及进一步减少Newton更新，冻结R7已补齐三组完整配对：
数值逐值一致，六组墙钟比均≤1.5，但Vela墙钟样本CV为26.4%/30.2%，绝对时延
仍有明显波动。9月15日恢复后，R8密度坐标候选已补初始化预测保护并通过
回归；从中性300 K开始完成62点电学/热学/2%局部场及本轮同VM配对墙钟
验收。漏压更新366/321→354/314，线搜索496/446→448/401；墙钟比
1.049/1.409。R8仅各一组计时，未证实稳定墙钟加速，作为显式可选配置，
默认仍为R7。剩余重点是密度更新与步长增长、精确点附近历史预测的协调。见
[重复性与Newton定位](templates_ldmos_r7_stability_newton_2026-09-14.md)。
复现入口见[生产运行指南](templates_ldmos_production_reproduction.md)，
冻结结果见[联合验收记录](templates_ldmos_production_electrothermal_2026-09-14.md#r7完成后的暂停记录)。

此前`c64039b`提交前的Windows UCRT64 Release检查为817/817（139.22 s）。
本轮A–E实现后全目标Release构建成功，匹配CTest为822/822（129.42 s）；
E同版Linux七项定向检查及两个停滞样本对照已完成；这些检查不替代冻结Linux
程序的完整曲线和计时资格。本轮A–E均未获得可晋级的总成本收益。

| 范围 | 状态 | 证据 |
|---|---|---|
| 密度坐标Newton / R8 | 仅可信历史预测且Vd≤1 V启用，初始化关闭；双栅压完整62点电热/局部场与各一组同VM墙钟门限通过；更新354/314、线搜索448/401；Vg8多一次停滞恢复，未认领稳定墙钟加速；R7仍默认 | [本轮验证](templates_ldmos_r7_stability_newton_2026-09-14.md#r8完整曲线计时与联合验收结果) |
| 电热装配候选 / R6 | 75个迁移率/输运用例、Release814/814及G3/D5/D4复核通过；Linux Vg4到28 V已883.08 s，超过原生完整465.27 s的1.5倍，保存22精确点后停止；29次尝试数值字段与R4相同，但分解变慢抵消装配收益；Vg8未启动 | [候选验证](templates_ldmos_production_electrothermal_2026-09-14.md#r6装配候选) |
| 电热线性后端 / R7 | Linux UMFPACK双栅压从零初始化完整62点独立联合验收通过；墙钟626.65/516.07 s，原生610.06/493.71 s，分别1.027/1.045倍；真实恢复及入口一致性通过；版本化输入导出与真实零压初始化检查通过 | [生产运行指南](templates_ldmos_production_reproduction.md) |
| 电热性能基线 / R4 | Linux双栅压从中性初值开始完整62点通过原电学/热学及2%局部场门限；同环境墙钟989.36/1142.43 s，分别为原生2.34/2.26倍，性能门限未通过；G3及D5/D4本版复核通过 | [集成与性能记录](templates_ldmos_production_electrothermal_2026-09-14.md) |
| 电热生产入口 / R2 | 双栅压从中性初值开始的0–40 V完整62点通过原电学/热学及新版局部场门限；代表点检查点恢复逐项一致，Release810/810；同环境R2仅记录主动暂停的部分曲线，完整性能结果见R4 | [集成记录](templates_ldmos_production_electrothermal_2026-09-14.md) |
| D0 新增联合验收局部场门限 | 用户批准新版载流子加权相对RMS各2%，带边30 meV等其余门限不变；P3及生产R2各自完整62点全部通过。旧1%结果58/62单独保留；R4完整62点也已通过；同环境R4性能未通过 | [新合同及差异定位](templates_ldmos_joint_acceptance_2026-09-14.md) |
| 原生Poisson几何 / N1 | 新几何双栅压完整62点独立通过；最大电流误差0.02409%/0.00714%，峰温误差0.00940/0.04198 K；最大导带差约7.8/26.7 meV，局部仍有余差 | [本轮执行](templates_ldmos_native_poisson_and_preparation_2026-09-13.md) |
| 电热准备 / 停滞退出 / SparseLU / IALMob屏蔽根 | P3双栅压完整62点原门限通过，全部接受状态与N1逐项一致；Vg4完整406→358更新，Vg8保持318；累计CPU870.95/602.48 s；Release808/808通过，计时边界及串行对照见报告 | [性能验证](templates_ldmos_native_poisson_and_preparation_2026-09-13.md#p1-温度物性准备复用) |
| 共享 Auger 生成开关 / 最新等温复核 | 新显式配置关闭生成项；本轮R4b已重跑G3 31点、D5/D4各62点并通过原门限；后两项为精确状态重闭合，不替代新模型从零推进性能 | [最新执行](templates_ldmos_generation_alignment_2026-09-13.md) |
| 自热 / 历史B4 R3基准 | 有限hRecVelocity、关闭Auger生成；Vg4/8各自从零漏压完成31点，原电学/批准热学联合通过；419/325次更新，子进程累计墙钟2609.58/2107.81 s；显式全场预测，未继承旧R4资格 | [B4完整验收](templates_ldmos_generation_alignment_2026-09-13.md#r3同一冻结b4模型的全场线性预测) |
| D0 局部统计与 Poisson 几何隔离 / 历史定位 | 原生Fermi控制及40 V零迭代对照定位B4几何合同差异；原报告只认领固定状态诊断，后续N1独立实现和自洽资格见本表首行 | [历史诊断](templates_ldmos_generation_alignment_2026-09-13.md#40-v-原生状态的算子与几何隔离) |
| G3 Id–Vg | 高精度Fermi/BGN、来源一致材料及Auger H/N0配置：31点、原六门限通过 | [本轮实现](templates_ldmos_auger_and_stage_chain_2026-09-12.md) |
| D5 Vg=8 V，0–40 V | 新H/N0配置从零压完整31点通过，1096次更新，0回退 | [本轮实现](templates_ldmos_auger_and_stage_chain_2026-09-12.md) |
| D5 Vg=4 V及双栅压电流比 | 新H/N0配置完整31点、1126次更新；双栅压原工程/最终门限通过 | [本轮实现](templates_ldmos_auger_and_stage_chain_2026-09-12.md) |
| 正确 Auger 单位配置 | 已完成 62 点完整推进联合验收，晋级运行入口默认；历史系数另行显式回放 | [本次执行](templates_ldmos_worker_and_ialmob_execution_2026-09-11.md) |
| 常驻 DC 进程 | 显式 `--worker` 复用输入准备，保留逐请求种子与原门限；验证范围见报告 | [本次执行](templates_ldmos_worker_and_ialmob_execution_2026-09-11.md) |
| IALMob / D4 | 新H/N0配置双栅压0–40 V各31点完成，原工程/最终联合门限通过；Vg4/Vg8最大电流误差1.8918%/1.8015%，1119/1073次更新 | [本轮实现](templates_ldmos_auger_and_stage_chain_2026-09-12.md) |
| IALMob 记录性能 | 共享几何及 referenced-QF 状态，前八点651.584→456.333 s，减少29.97%；双栅压完整62点复核通过，状态哈希全部不变 | [性能及差异定位](templates_ldmos_ialmob_performance_followup_2026-09-11.md) |
| 高精度 Fermi/BGN 及材料一致性 | 函数、解析导数、逆映射及材料审计已实现；Release原768项及新增回归通过，G3/D5/D4共155个精确点复核通过 | [统计精度实现](templates_ldmos_fermi_accuracy_2026-09-12.md) |
| Auger H/N0 | 自洽源项及导数已实现，G3/D4/D5共155点完整复核通过；独立新冻结资格，未覆盖旧证据 | [本轮执行](templates_ldmos_auger_and_stage_chain_2026-09-12.md) |
| D3/D2/D1 | 原生脚本和曲线重核后，依次通过原工程/最终等效曲线门限；限原Solve未启用量子方程且等温接触消融不可分辨的路径 | [本轮执行](templates_ldmos_auger_and_stage_chain_2026-09-12.md) |
| 自热 / D0 历史R4 | 独立四方程冻结Release完成300 K双栅压62点回放、代表点及0–40 V完整62点；原电学及批准热学门限全部通过；最大电流误差1.7744%/1.5825%，峰温误差1.5051/2.1893 K；对应旧有符号Auger及理想空穴接触模型 | [电热实现与本轮证据](templates_ldmos_d0_electrothermal_2026-09-12.md) |
| 有限 hRecVelocity / 局部物性对照 | A→B阶段完成有限空穴边界及Jacobian、D0关闭Auger生成，双栅压各4个代表点通过，Release802项通过；后续B4完整资格见本表最新R3行，局部场仍有差异 | [A/B阶段证据](templates_ldmos_hrec_and_local_physics_2026-09-13.md) |

## 配置合同与冻结版本的适用范围

`contracts/physics_contract.json` 保留 G3 历史资格使用的 `edge_projection`。
此前完成资格的 D5 使用显式的
[`linked_d5_auger_units_config.json`](../../reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_auger_units_config.json)，
选择 `transport_cell_vector`、接触电场回退、按需 QF 重定心及内部密度恢复。
该独立配置不修改通用默认值。向量梯度在 G3 的资格来自本次独立 31 点复核。
数值块、逐行和 KCL 门限分别沿用各算例原值；G3 六门限并不等同于 D5 门限。

2026-09-11 参数审计新增发现：历史生成器把 SI Auger 系数直接写入按 cm⁶/s
读取的 `unit_scaling` 字段，历史曲线的显式系数比目标小 10¹² 倍。
生成器已修正并通过回归；冻结 `linked_d5_config.json` 保留旧值以复现历史对照，
不得把它视为 Auger 参数完全对齐的最终推荐配置。
另提供 [`linked_d5_auger_units_config.json`](../../reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_auger_units_config.json)
及对应 inputs 清单；该配置已通过从零压完整联合验收并成为运行入口默认。
H/N0增强已新增显式自洽项及导数，默认关闭。独立
`linked_d4_auger_density_inputs.json`、`linked_d5_auger_density_inputs.json`
已使用新Release完成原门限复核；使用本轮冻结清单，不能自动继承旧程序资格。

2026-09-12高精度统计实现另提供
[`linked_d5_fermi_accurate_inputs.json`](../../reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_fermi_accurate_inputs.json)
及独立D4清单，使用来源审计的300 K材料和新模型重新资格的零压种子。
本轮新版G3/D5/D4完整资格见统计精度报告。运行时显式传入
该清单及报告中的新冻结程序，不能把旧默认bundle、旧种子或旧二进制的
资格自动转移到新模型。

## 独立运行入口

[`run_templates_ldmos_linked_d5.py`](../../scripts/run_templates_ldmos_linked_d5.py)
不导入历史实验 Python 脚本。配置和哈希清单纳入版本管理；外部网格、掺杂、
AverageBox 系数、合格零压种子及 Sentaurus 参考曲线仍需单独保留。
[`linked_d5_auger_units_inputs.json`](../../reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_auger_units_inputs.json)
列出各文件相对工作区的路径及 SHA256；缺失或变化会拒绝运行。

从 worktree 根运行，例如：

```powershell
D:\msys64\ucrt64\bin\python.exe -X utf8 scripts/run_templates_ldmos_linked_d5.py --manifest reference_staging/templates_ldmos_low_bias_recovery_20260911/binary/manifest.json --output reference_staging/unique_vg4_first8 --gate 4 --points 8 --preflight
```

`--preflight` 只读检查，不创建输出；删除该选项执行。输出目录必须不存在。
默认 bundle 使用正确 Auger 单位。历史回放另加
`--bundle reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_inputs.json`，
不能把两组参数的运行时间混为同一基准。
`--worker` 需要包含 `--dc-worker` 的新版冻结 Release；历史 runner 不支持该选项。
该模式只复用输入准备，每个请求仍重新加载其显式状态；未启用时沿用逐请求子进程。
`--points 2` 用于迁移对照，`8` 为前八点，`31` 为全曲线。
31 点评分才包含 40 V 端点门限，短曲线不能认领完整 D5 通过。
默认可信预测最大步长 0.2 V、无预测最大 0.1 V；显式 `--max-step` 的变化
属于新性能实验。本次 0.4 V 前八点虽减少推进，却增加 46.26% 更新，已不予晋级。
高压参考系平移保留电势、密度、电流及 KCL 等价性检查。
`STOP` 文件仅在子步骤边界请求停止，不接受未通过的状态。

运行记录包含输入和程序哈希、每个求解请求配置/父状态/trace、原门限检查、
精确点状态哈希、全部 Newton 更新和子进程时间。Windows 使用 UCRT64
Release、Eigen SparseLU/COLAMD；此入口不将 Sentaurus 日志时间解释成
相同环境下的墙钟基准。冻结二进制清单需随外部证据保存。

## D4 独立配置入口

同一受版本管理运行器新增显式 `--physics-profile D4`，要求独立
`vela.templates_ldmos.linked_d4_inputs.v1` 清单。配置/几何准备入口为
[`prepare_templates_ldmos_ialmob_transport.py`](../../scripts/prepare_templates_ldmos_ialmob_transport.py)。
D5 默认及其合同保持不变。D4 零压初值必须重新闭合，不能继承 D5 解的物理资格。
当前验证进度和角点处理限制见[全局接入报告](templates_ldmos_ialmob_global_coupling_2026-09-11.md)。

D4 完整 62 点已通过原联合门限；仓库内的
[`linked_d4_inputs.json`](../../reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d4_inputs.json)
为此前统计版本的独立输入合同。新版采用
[`linked_d4_fermi_accurate_inputs.json`](../../reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d4_fermi_accurate_inputs.json)，
亦完成62点原门限复核。清单默认参考系偏移分别为 Vg4：28 V、Vg8：4 V，
均须通过原等价性门限。`--resume-from` 支持经审核的只读前缀续跑到新目录；
必须保留源证据目录。具体失败尝试、修正和资格边界见完整报告。

## 本轮执行与后续顺序

1. **IALMob 记录优化已完成短曲线筛选。** 几何构造408→68次、状态构造
   772→567次，记录时间261.130→45.247 s。原348次Newton及八个状态哈希
   不变。完整双栅压已在独立 `_full_r2` 目录完成，62个状态哈希全部与
   原资格一致，原联合门限通过；并行曲线墙钟不作为独占性能比较。
2. **高压长尾的完整步拒绝原因已核实。** 同一26.666667→26.766667 V
   重放：16 V参考系32次更新中观察到19个完整步被电子块1e-11门限拒绝；
   4 V参考系8次更新，无这类拒绝。观察结果与历史状态哈希相同。
   继续保留本算例4 V策略，不放宽门限，不把浮点步长尾数的路径变化当收益。
3. **统计近似与BGN联动已实现，完整资格复核通过。** 原生D4状态的Vela漏电流
   回放误差仅0.0121%，但全域方程不合格。独立积分和逆函数对照已定位
   Fermi近似及其BGN修正的密度偏差。高精度函数、导数、逆映射及独立
   材料配置已完成，G3原六门限通过；D5/D4使用新模型重新资格的零压种子
   从零完整推进，原工程/最终联合门限均通过，共155个精确参考点。
   固定原生状态密度误差中位数降至约10 ppm，完整曲线电流差异仅略降，
   剩余局部输运/离散差异仍需定位，见[本轮实现](templates_ldmos_fermi_accuracy_2026-09-12.md)。
   H/N0已接入自洽项；新固定状态对照仍不支持其为电流差主因。
   电子残差平方和97.7%集中在距Si/氧化层界面10 nm内；强行将密度逆映射
   为准费米势反而放大残差，不能作为修正。角点及其他晶向的推广仍需独立证据。
4. **独立四方程 D0 完整联合验收已通过。** 用户指定先D3/D2/D1
   后原D0；批准0.1%热平衡及1 K/5%温升标准。独立热传导、热边界及固定原生
   热源回放已完成；Vg4/Vg8峰温误差0.0801/1.3822 K，场RMS误差0.1052/0.4048 K，
   均通过该限定范围热门限。后续已新增逐节点温度物性、与离散电流一致的
   默认晶格源及四方程Jacobian，使用独立实验入口验证；生产DC默认仍为等温。
   新冻结Release完成300 K双栅压62点回放，电流最大相对变化8.8263e-14；
   代表偏压及完整D0双栅压62点全部通过原电学和批准热学门限。两曲线各
   38次非零推进、0拒绝，Newton更新369/366次；全部78次求解重审通过。
   最终证据及适用范围见[电热报告](templates_ldmos_d0_electrothermal_2026-09-12.md)。

常驻 worker 和更大步长继续保留为显式实验选项。已有对照中 worker 请求
墙钟增加 9.09%，0.4 V 步长增加 46.26% Newton 更新，均不作为已实现的
性能收益。完整的新证据见[本轮报告](templates_ldmos_ialmob_performance_followup_2026-09-11.md)。
现有资料足以继续局部差异定位和电热耦合；自热范围及热门限已获确认。

## 已提交基线检查（5ea9b85）

最终工作区 Release 构建无待编译项；完整 Release CTest **767/767 通过，
99.49 s**，包含最终角点用例及运行器回归。18 份相关文档的 68 个本地链接
检查通过。日志保存在
`reference_staging/templates_ldmos_followup_20260911/ialmob_global_r4/commit_release_ctest.log`。
此次整理没有重跑长曲线；上述数值资格继续对应各报告中的冻结证据。

2026-09-12 未提交工作区另通过完整 Release CTest **768/768，183.92 s**。
该检查包含几何/状态复用及高压 profiling 观察，日志为
`ialmob_perf_r5/resumed_release_ctest.log`；不能与上面已提交基线混为一次检查。

本轮高精度统计实现另通过完整Release **768/768，101.17 s**，日志为
`templates_ldmos_fermi_accuracy_20260912/r1/release_ctest_r2.log`。
随后增加材料审计CTest注册，当时为769项；新增材料及原运行器定向2/2通过，
内部13个Python用例通过。此次测试注册变化没有替换曲线使用的冻结程序。

本轮Auger实现首先通过完整Release CTest 773/773。加入独立热传导后全量
780项中779项首次通过，唯一失败为新热合同误入冻结三文件合同目录。
将其移到独立`thermal/`目录后相关两项2/2通过，未放宽测试或物理门限。
本轮最终曲线及热学结果分别对应各自冻结程序，详见新执行报告。


## D0 温度相关迁移率后续（2026-09-12）

独立 IALMob 低场及 alpha=0 高场温度项、局部温度偏导已实现。
Vg4/Vg8、Vd40 V 原生热态对照覆盖全部10515个Si单元和5723个Si节点；
低场最大误差≤6.92 ppm，高场节点最大误差≤4.09 ppm。
该局部迁移率阶段当时未新增D0自洽曲线资格，也没有闭合D4剩余界面残差。
后续温度相关统计/复合、热源及电热耦合已在独立四方程入口实现并完成
完整D0资格，见下节；生产主电学入口仍默认300 K。
见[温度迁移率实现与证据](templates_ldmos_d0_mobility_temperature_2026-09-12.md)。

本轮最终完整 Release CTest 786/786 通过（253.81 s）。

## D0 四方程完整验收（2026-09-12）

新增逐节点温度物性及解析偏导、低/高场迁移率温度链式导数、与实际离散
电流一致的默认晶格热源及完整四方程Jacobian。独立入口使用UCRT64 Release、
Eigen SparseLU/COLAMD，原D0不启用Thermodynamic/Peltier/RecGenHeat。
两栅压0–40 V共62个精确点全部通过；最大温升场RMS误差为0.9552/1.5001 K。
完整Release CTest 799项中798项首次通过，唯一新增Python中文合同解码问题
改为显式UTF-8后重跑通过；后续输入/门限定向回归8个用例通过。
求解、冻结源码、原生62个温度场及最终评分证据见
[完整电热报告](templates_ldmos_d0_electrothermal_2026-09-12.md)。

后续可优先研究四方程Newton更新限幅成本，继续定位局部带边/界面输运差异，
并补原生Auger温度场逐节点对照。将独立电热入口整合到生产DC运行器及推广
到其他网格/边界时需要独立验证；本算例资格不会自动覆盖这些扩展。

## 本次提交前检查及日报素材

2026-09-12 最终工作区 Release 构建无待编译项；完整 CTest **799/799 通过，
777.26 s**，包含最终输入保护和四方程停止原因回归。日志为
`reference_staging/templates_ldmos_d0_electrothermal_20260912/commit_release_ctest.log`。
测试时间受宿主负载影响，不作为求解器性能对比。本次未重新运行长曲线，
数值资格继续对应报告中的冻结程序。

[工作日报素材及三组配图](templates_ldmos_daily_report_2026-09-12.md)提供完整
电流/峰值温升对比、自热电流下降和高偏压温度场，附 PNG、PDF 和复现入口。
