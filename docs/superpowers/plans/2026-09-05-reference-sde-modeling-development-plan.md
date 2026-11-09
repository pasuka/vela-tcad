# reference_tcad SDE 建模与 GRD/DAT/MSH 导出开发计划

更新日期：2026-11-01。当前状态：P1–P5 完成；P6 建模与交付主范围完成。
main `72a5b129` 的实际 Ubuntu CI 两个 job 已通过，当前离线包及校验清单已本地归档。
现行交付与未完成扩展见 [交付说明](../../sde_delivery.md)；历史商业对照仍对应 `e1b83fd3`。
本文件替代 2026-09-05 版本的待办顺序和覆盖统计。旧实测保留在
[9 月审计](../../validation/reference_sde_audit_2026-09-05.json)，不作为当前资格。

## 目标、依据和合并基线

从原始二维 SDE 脚本离线读取几何、材料、接触、掺杂及网格约束，经 Gmsh 生成同源的
`device.grd`、`device.dat`、`device.msh`。建模入口只需要 SDE；商业 SDE/TDX、TDR 和日志
仅作内部验证。DAT 的当前必交付是掺杂场，求解场导出另行安排。

- 历史开发分支：`codex/sentaurus2022-sde-workflow`；早期合并基线 `87fc8f56`。
- 开发提交已 cherry-pick 到 main，独立合并修复收束于 `9aa5b942`；当前验收基线为 `72a5b129`。
- 当前收尾工作在 main checkout 进行；本轮更新文档并归档既有 CI 产物，不推送或发布 Release。
- 当前规范以 [文档索引](../../README.md)、[架构](../../architecture.md)、
  [配置合同](../../config_schema.md)、[回归指南](../../../tests/regression/README.md)、
  [PN2D BV 政策](../../validation/pn2d_bv_validation.md) 为准。
- 新的求解入口与发行资格依据 [SDevice 交付指南](../../sdevice_delivery.md)、
  [支持矩阵](../../sdevice_support_matrix.json)、[sd-sim CLI 计划](2026-10-23-sd-sim-cli-plan.md)。

## 最新代码带来的计划调整

| 当前事实 | 对开发计划的影响 |
| --- | --- |
| 主线 `d58dc92e` 已退役 custom NMOS、coarse7x3、Minimal6、skewed-Tri3 案例 | 删除原计划针对这些案例的里程碑；不恢复已删除的目录和合同 |
| 当前主案例源目录有 8 份 SDE 几何脚本，含新增 gate_reference | P5 已达到 8/8；文件数不代表器件族数 |
| 主线新增 `sd-sim`，保留 `vela device`，具备解析、材料、阶段发射、编排、打包 | 复用现有 SDevice 消费端，不再新增并行的求解编排器 |
| 原生导入支持 GRD + 同名 DAT，TDR 依赖检测到的 HDF5；无原生 MSH 消费路径 | 以 GRD/DAT 桥接求解端，MSH 面向建模交换和回读；暂不增加 MSH 求解导入器 |
| 自产八例 GRD/DAT 已通过主线 `mesh_import.import_grid`，P2 三例八份命令 --check 仍通过 | BJT 自主建模已验收；其数值求解和跨平台资格单列 |
| PN2D 既有包在具名输入/网格/材料/profile 上取得数值资格 | 1482/2800 的自产网格不自动继承 1943/3680 参考网格的数值资格 |
| BVmethods NonLocalPath、外电阻 R2/R3 仍有独立未收口门禁 | 建模完成不等于完整 BV 仿真完成；不把该数值研究并入本建模项目 |
| 主线要求 Eigen 5.0.1、系统 Boost/PFR；PARDISO 默认编译关闭，CI Ubuntu 26.04 | 更新构建和发行依赖清单；不以切换求解后端解决建模问题 |

仍复用分支已有 `offline-static → Device IR v2 → Gmsh → mesh-field bundle → exporters`。
原 M0–M5 是代码基础及历史证据，不是当前合并树已完成数值验收的结论。

## 2026-10-30 实测覆盖

盘点范围为 `reference_tcad/*/source/*.{cmd,scm}` 中含几何创建的主案例源文件。
变体归档、engineering 冻结副本、oracle wrapper 不混入这个分母；未来增加支持前单独登记。
初始合并基线见 [本轮审计](../../validation/reference_sde_audit_2026-10-29.json)，当前结果见
[P5 验收](../../validation/reference_sde_p5_acceptance_2026-10-30.json)。

| 主案例源脚本 | 结果 | 当前阻塞/后续资格 |
| --- | --- | --- |
| `gate_reference/source/gate_sde.cmd` | 重复合同通过，BAMG 非字节确定 | 格式回读、原生导入及 P2 gate/contact 命令检查通过 |
| `pn2d_sentaurus2018/source/pn2d_sde.cmd` | 生成成功：1482 / 2800 | 格式回读、原生导入通过；自产网格数值资格单列 |
| `schottky_charon_sentaurus2018/source/schottky_n_sde.cmd` | 生成成功：319 / 560 | 格式回读、原生导入通过；接触类型由 SDevice 材料/电极合同继续验证 |
| `genius_bjt_sentaurus2022/source/bjt_sde.cmd` | 两次生成/独立导入合同通过 | General、局部自适应、节点掺杂/积分/结位置通过；数值资格另设门禁 |
| `genius_bjt_sentaurus2022/source/bjt_sde_local_refined.cmd` | 两次生成/独立导入合同通过 | 同上；三个原始局部窗口 P95 网格尺寸均小于基础版本 |
| `singledevice_sentaurus2018/source/singledevice_sde.cmd` | 两次完整生成及半/全独立导入通过 | 圆角、占位符、区域/窗口加密与本地镜像；原 SDevice E403 边界单列 |
| `transportmodels_sentaurus2022/source/nMOS_dvs.cmd` | 两次完整生成及半/全独立导入通过 | 三类快捷加密、本地镜像及节点场；原 SDevice E105/E109 边界单列 |
| `bvmethods_sentaurus2018/source/bvmethods_nmos_sde.cmd` | 两次完整生成及半/全独立导入通过 | PolySi 别名、栅极移除、多区域接触、multibox、Erf、x 镜像；BV 数值资格另列 |

8 份原脚本均成功；未删除或简化原始命令。商业对照另有 `+z` 上下文包装，与离线原脚本建模分别记账。
本轮没有改写这些 SDE 源脚本。SimpleMOS、LDMOS 的固定工程网格/已有求解资格不等于
对应原始 SDE 已可自主建模；它们作为后续扩展范围，不擅自计入 8/8。

## 交付合同

- 三个主文件必须来自同一规范化 bundle。GRD 保存几何、区域、材料、接触；DAT 保存
  donor/acceptor、净掺杂及可保留的 active species；MSH 4.1 ASCII 保留 Physical Groups 和 NodeData。
- 长度 µm、掺杂 cm^-3。ID 映射、输入/输出哈希、Gmsh 版本及影响拓扑的参数必须可追溯。
- 明确 region-local/shared-node 掺杂与 placement 顺序；原生消费端当前拒绝 region-local
  dataset 排列（E402），扩展时必须显式合同转换，不能做无声平均。
- 检查正面积、单元方向、区域面积/覆盖、界面共形、接触边界和物种字段；格式回读逐节点比较。
- 几何生成、格式回读、原生导入、`sd-sim --check`、求解、参考数值资格分别记录状态。
  日志中的成功字样或单独退出 0 不代替完成合同；求解目标失败按现有 `--exit-on-failure` 处理。
- 未支持命令带文件/行/命令失败；失败后清理旧成功标记，防止旧 GRD/DAT/MSH 被误认为新产物。
- 保留 PN2D BV 原子 profile 和非钝角门，不推广为全局默认，不以改模型、放宽容差换取覆盖率。

## 修订后的实施顺序

### P0：合并收尾与当前基线（本轮）

保留主线案例退役和唯一 Linux workflow，合并 Gmsh 安装/probe、双方有效 CTest 注册，
以及哈希封印文件属性。旧 coarse7x3 测试迁移至现存 PN2D，步长断言遵循已有 qualified
profile 的 0.02/0.02；不更改生成器或验收阈值。保留分支完整 PN2D BV 的 Gmsh 合同，
它仍须在新合并树上重新取得数值证据，不能因合同文件存在而报已 qualified。

退出条件：main 是祖先、无冲突标记、CMake 配置成功、受影响回归通过、变更检查通过。
整合后的 main 已通过 C++ 构建及 Ubuntu 全量 CTest（4 项条件跳过）；
Windows 全量回归与新的器件数值资格不由该 CI 结果覆盖。

### P1：三个现有案例的可重复产品入口（已完成，2026-10-29）

新增 case manifest 与批量 regression：源/材料哈希、区域/接触名、尺寸和掺杂合同、成功阶段。
把本次三例直接生成和原生导入纳入自动回归；同时冻结五例的可定位拒绝边界。
验证重复运行、失败重跑、路径含空格、干净输出目录和原文件保护。
保持当前仅需 SDE 的独立 CLI；是否新增用户入口名称在统一 CLI 设计时再确定。

退出条件：三例稳定生成 GRD/DAT/MSH，独立读回与消费端身份对比通过，统计仍为 3/8。
不得因为现有 SDevice 包支持某项 Physics，就把相应 SDE 建模标成完成。

P1 实施与证据见 [使用说明](../../reference_sde_cases.md) 和
[验收摘要](../../validation/reference_sde_p1_acceptance_2026-10-29.json)。已交付版本化八例合同、
批量 CLI、三例解析几何/逐节点掺杂/完整接触/消费端身份检查、五例定位拒绝及失败恢复保护。
PN2D/Schottky 三次重复字节一致；用户已明确选择 gate_reference 按几何、接触、掺杂和网格
质量合同重复验收并保存不同网格哈希，不宣称 BAMG 字节确定性。P1 专项 18 项和相关
9 个 CTest 均获得通过证据；没有修改源脚本、物理模型或既有数值门槛。

### P2：与当前 sd-sim 的完整接口验证（已完成，2026-10-29）

复用 `scripts/vela_sdevice/mesh_import`、材料解析、deck_emit、orchestrator、sidecar/schema。
以自产 GRD/DAT 替换本地测试副本的 `File{Grid}`，记录原始与有效输入哈希，不改冻结源。
为具名 PN2D、gate、Schottky 输入分别检查 electrode/region/material/species/单位身份、
Parameter 声明和 gate/contact 边界。优先执行 `--check`；短数值控制和完整 IV/BV 另设门禁。
MSH 保持外部交换格式，不通过引入 TDR 作为强制中间文件来连接两端。

退出条件：三例自产网格可进入现有命令前端的合法检查路径；不支持的物理必须正确拒绝。
PN2D 数值比较同时保留已 qualified 网格控制组和自产网格实验组，不复用旧包资格结论。

P2 已实现 [检查工具与输入策略](../../sde_sdevice_checks.md)，验收见
[结构化记录](../../validation/reference_sde_p2_acceptance_2026-10-29.json)。三例自产网格上的
8 份具名 SDevice 输入通过真实 sd-sim --check；P2 专项 17/17、相关 CTest 8/8 通过。
修复 Interface 线记录误入体材料导致的 gate E403，保留完整接口边及两侧归属；补齐
Parameter、材料/单位、原始/有效哈希、sidecar 和阶段/实际 sweep 目标身份核验。
E301/E303/E401/E403/E106/E409 的实际拒绝边界均有回归。原始缺 Parameter 的输入仍按
E301 记录；Schottky 默认参数仅用于显式检查变体。PN2D 冻结参考控制组与自产网格实验组
分别记账，不继承数值资格；完整求解、参考重复计算和 Linux 资格均未在 P2 执行。

### P3：BJT 两份原脚本（已完成，2026-10-29）

实现有类型 body handle 与 define 中创建几何，保留执行顺序和 source map；补充实际所需
有界纯辅助函数。General 参数串和 x/y、abs、exp、幂运算使用受限 AST，禁止 Python eval。
扩展 IR/schema、迁移和能力矩阵，不能用 Gaussian 或常数替代 General 解析掺杂。
在连续掺杂求值基础上实现 MaxTransDiff：先固定阈值/变换/停止语义，再有界重网格，
记录每轮指标、节点预算和失败停止；不能以全域最细固定网格冒充自适应等价性。

退出条件：两份原脚本直接导出，三个接触、掺杂截面/积分、结位置及局部加密合同通过，累计 5/8。
旧 TDR 转换得到的 BJT 输入仅作为对照，不冒充自主建模结果。

P3 已交付 [使用说明及映射边界](../../sde_bjt_modeling.md) 与
[结构化验收](../../validation/reference_sde_p3_acceptance_2026-10-29.json)。两份原脚本均保持
原始哈希，typed body handle、有界纯函数、General 受限 AST 与 IR/schema 迁移通过回归。
MaxTransDiff 使用带最长边闭包的一致局部细分，固定 asinh 变换、源 Value、方向最小尺寸
停止语义和资源预算；逐轮统计明确保留最小尺寸停止后仍超阈值的边，不宣称 native snmesh
三角剖分等价。用户事先确认的节点误差 1e-6、积分误差 1%、局部 min_y 结位置门槛全部通过。
两次重复各自通过几何/接触/掺杂/局部尺寸与原生消费端身份合同，并逐次保留非相同哈希。
连续掺杂暴露的 DAT 补偿点精度损失已通过 17 位有效数字输出修复，没有放宽既有容差。
P2 八份 --check 仍通过；未执行 BJT IV、live SDE/TDX 或新数值资格。

### P4：SingleDevice 与 TransportModels 完整器件（已完成，2026-10-30）

补齐辅助函数、body 绑定、fillet-2d 曲边及 chord-error 合同、细氧化层/界面加密。
将 refine-box、refine-doping、refine-interface 明确转换到 IR，冻结缺省值与叠加顺序。
显式提供 @node@ 等预处理参数；已知 tdx 镜像命令解析为受限本地 bundle 操作，不执行任意 shell。
镜像包含轴上节点合并、单元方向、材料/接触重命名及全部节点场；半器件与全器件分开验收。

退出条件：两份原脚本完整器件导出及独立半/全原生导入通过，累计 7/8；未知参数仍明确失败。
用户于 2026-10-30 明确选择：保留原始 SDevice 设置，消费端按明确拒绝边界验收；不扩展材料/求解前端。

P4 实现与边界见 [使用说明](../../sde_mos_modeling.md) 和
[结构化验收](../../validation/reference_sde_p4_acceptance_2026-10-30.json)。两份源脚本未改写；
0.05 um 原始半径、1e-5 um 弦误差、边界推导面积预算、1e-6 掺杂误差和 1e-8 um 镜像/接触
容差均通过。区域局部 MaxTransDiff 与界面细分保持共享边一致；最终网格再镜像，保留所有标量场。
每次重复分别保存半/全哈希并独立导入；不宣称 BAMG 字节确定性或非钝角资格。
P2 八份命令仍正向通过；原 MOS 命令按 E403、E105/E109 明确拒绝，没有改动物理设置。

### P5：BVmethods 建模（已完成，2026-10-30）

显式 PolySi/PolySilicon 别名，保留原始材料名。实现 body contact/remove 的创建态 ownership、
保留边界及被删除区域后续引用语义；以小型内部 oracle 确认后实施。
实现 multibox 窗口、正负 grading ratio、裁剪及组合；复核现有 Erf 与 Symm 方向实现。
复用 P4 的圆角和本地镜像，确认移除 gate body 后接触与源漏对称性。

退出条件：主案例 8/8 原脚本直接建模导出。该结果不解除 NonLocalPath、外电阻 R2/R3、
纯 QS 轨迹或原生接受过程的独立限制；不得删除不支持命令或简化源脚本来报完成。

### P6：资格、发行和扩展（建模与交付主范围完成；扩展待规划）

历史本地验收完成 Windows UCRT64 与 Ubuntu 26.04 独立构建。按用户后续选择，
GitHub 仅保留 `ci.yml`，两个 job 均使用 Ubuntu 26.04；main 实际 CI 已通过，
Windows 不列为 GitHub CI 要求。确认 Gmsh runtime
版本、加载路径、许可/通知与离线依赖，现有 SDevice 包不自动包含新增建模工具和 Gmsh。
保持 `sd-sim` 已预留的 `-d`、`-n`、`--max_threads` 等状态，遵循其独立 CLI 计划。
内部 TDX 仅校验自产 GRD/DAT；官方 SDE 几何/掺杂对比按面积、界面、截面、积分和尺寸
统计验收，不要求不同网格节点编号相同。几何容差及掺杂/尺寸门槛必须实施前冻结。
先完成 8 份主案例，再按用户优先级登记 SimpleMOS/LDMOS、Polygon body、其他布尔模式、
非零 decay 等扩展；3D、ACIS、完整 Elk Scheme 不在本轮承诺范围。

旧版对 P1–P6 的 31–48 工作日估算仅为历史规划，不再作为当前待办工期。
新增扩展应按各自合同重新估算；Polygon/skewed-Tri3 的旧优先顺序不恢复为当前任务队列。

## 证据与仓库政策

TDR、`.log*`、stdout/stderr 和原始运行记录只写 ignored build/ 或用户本地目录。
公共仓库保留源、schema、小型自有 fixture、结构化审计摘要与哈希；不得改扩展名或将原始
日志完整嵌入 JSON 绕过排除规则。生成 GRD/DAT/MSH 是产品输出，默认同样现场生成而不提交。
旧 Git 对象中的历史原件本轮不重写；本计划不新提交 TDR/原始日志。
私有主机、地址、用户、安装路径和 VM job ID 不写公共计划。P5 按用户请求执行商业语义对照；
实际公共 CI 链接和产物哈希另行记录，未启动新的商业或 IV/BV 数值资格任务。

## 本轮实际验证与下一步

P4 八例完整合同回归通过：7 例各两次自产完整模型，BVmethods 保留第 32 行 PolySi 拒绝。
两份 MOS 还逐次通过半/全独立原生导入并保存各自哈希。相关 13 个 CTest 和
20 项完整案例测试通过，合计 191 项 Python 测试（1 条件跳过）；实际构建范围和证据见
[P4 验收](../../validation/reference_sde_p4_acceptance_2026-10-30.json)。未执行 T470p/Linux、
完整 C++/全量 CTest 或 MOS IV/BV；配置 VM 的 SDE 许可证不可用，不记作商业验证通过。

P5 实际验证：八份源各两次完整建模和独立原生导入通过；17 项相关 CTest 及 20 项
完整案例测试通过，共 224 项 Python 测试、3 项条件跳过。sd-sim/pyz Release 构建通过，
未执行完整 C++/全量 CTest 或新数值资格；[结构化验收](../../validation/reference_sde_p5_acceptance_2026-10-30.json) 保留实际范围和逐次哈希。

P5 将 BVmethods 提升为第八份完整建模案例。创建态 body/remove 不复活旧材料，删除区域掺杂
声明保持已知但不活动；Gate multibox 明确不活动，Channel 保留原窗口与 grading。圆角半径
0.08 um、弦误差 1e-5 um、节点掺杂 1e-6 和接触/镜像 1e-8 um 门槛未放宽。
Vela VM 的 U-2022.12 已完成受控商业对照：ratio=1 实测采用最大尺寸，Erf 端点及尺度通过
区域局部字段审计。原 BV 源由独立 +z 包装加载，源文件不改写；原始输出和私有运行记录不提交。

P6 本地结果见 [离线包与资格说明](../../sde_p6_qualification.md) 和
[结构化验收](../../validation/reference_sde_p6_acceptance_2026-10-30.json)。Windows UCRT64 与
Ubuntu 26.04 独立构建通过；两份最终离线包各八例、各两次，32 次建模与独立导入、32 次商业
网格合同对照全部通过。最终 Windows 第二次八份 GRD/DAT 的 TDX 编码/回转及存储字段审计通过。
Gaussian 积分控制和局部 red 闭包修复了实测积分/轮数失败；1% 门槛和原资源预算未放宽。
每平台相关 18 项 CTest 分阶段通过，最终 P6 专项含 19 项单测，P2 八份命令检查持续通过。

2026-10-30 的本地验收当时未运行 CI；后续 main 的
[实际 CI](https://github.com/pasuka/vela-tcad/actions/runs/41193813343) 已通过，记录在
[main 验收摘要](../../validation/reference_sde_main_acceptance_2026-11-01.json)。
1190 项全量选择中 1186 通过、4 项条件跳过、0 失败；SDE 专项 18/18 通过。
当前 Ubuntu 离线包八例各两次通过，归档文件、SHA-256、4996 项载荷哈希和版本身份已核对。
该 CI 不重新签发历史商业对照或 Windows 当前版本资格。

本轮文档与交付证据收尾完成。后续仅保留待确定优先级的 SimpleMOS/LDMOS、Polygon、
其他布尔模式和非零 decay；这些登记项尚未实现。BVmethods 的 NonLocalPath、外电阻
R2/R3、新网格 IV/BV、MOS 材料/求解前端支持仍属独立门禁，不阻塞已交付的八例建模。
