# SDevice 前端交付范围与成果审查（2026-10-22）

## 结论与边界

应恢复以“具名参考 cmd → 有效 JSON → Vela 运行及输出”为交付主线。前端已有可用实现，
2026-10-07 两平台发行包已有具名数值验收；不能把仍未完成的安装库原版 BVmethods
NonLocalPath／外电阻研究当作整个前端尚未实现，也不能把这些实验算成已发布能力。

本轮发现真实的重新交付阻碍：**当前发射器改变了 PN2D 初始化流程，尚无本轮对应的两平台
原生数值回归**。这是前端主线应优先处理的问题，优先于继续追究 BVmethods 微小 Jv 差异。
没有证明该改动错误；也不能因新 golden 已通过，就沿用旧运行时的数值资格。

本轮完成范围盘点、重新发射比较、回归、旧包哈希复核和待办排序。未修改求解器、物理默认、
接受门限或 fixture golden；未启动 VM／Sentaurus 新作业、完整 BV 或重新打包。允许使用 VM
不意味着必须重复已有本地可验证证据。本轮没有取消用户先前授权的扩展研发，下面将其单列。

机器可读的[输入、模型、哈希与比较账本](sdevice_delivery_scope_audit_2026-10-22.json)记录
10 个冻结 fixture，以及当前 `reference_tcad/*/source/*.cmd` 的全部 **43 个文件**。
后者包括变体、消融与拒绝输入，文件存在不等于承诺支持或数值验证通过。
工作树 HEAD 为 `a0ed4e073efee0214c39365d3be90a12be51c624`，存在大量未提交研发修改；
本轮结果绑定实际文件哈希，不仅绑定 HEAD。

## 1. 案例与输入版本

原范围依据为[初始计划](../superpowers/plans/2026-09-26-sdevice-cmd-direct-run-development-plan.md)、
[交接 §6.4](sdevice_cmd_frontend_handoff_2026-09-30.md)、
[fixture 索引](../../tests/fixtures/sdevice/fixtures.json)及
[现行验收契约](../../configs/sdevice_acceptance_contract.json)。本轮 10/10 fixture 的 LF
规范化哈希均匹配索引及仓库来源；账本同时保留原始字节哈希，避免混同行尾转换。

| 输入（目录/文件） | 主要物理或流程 | 原文当前解析 | 已有资格上限 |
|---|---|---|---|
| `pn2d_sentaurus2018/pn2d_0v_sdevice.cmd` | OldSlotboom、Masetti、SRH，零偏 | E301，缺 Parameter | 补声明的 v2 输入：两平台数值 E2E |
| `pn2d_sentaurus2018/pn2d_iv_sdevice.cmd` | 同上，正偏 QS | E301 | v2：曲线／具名场量 E2E；曲线窗口依原契约 |
| `pn2d_sentaurus2018/pn2d_bv_sdevice.cmd` | 加 HFS、VanOverstraeten 雪崩 | E301 | v2：到 −20 V，十个固定场量点；不是任意网格资格 |
| `bvmethods_sentaurus2018/aba_poisson_sdevice.cmd` | Fermi、Poisson-only、电离积分、ElectricField | 0 | 解析／发射／具名 `--check`；非完整数值对齐 |
| `bvmethods_sentaurus2018/external_resistor_sdevice.cmd` | Fermi＋OldSlotboom、Masetti/Enormal/HFS、SRH(DopingDep)、E2、Eparallel 雪崩、外电阻 | 0 | 解析／发射／具名 `--check`；从零纯 QS 未获资格 |
| `schottky_charon_sentaurus2018/schottky_n_des.cmd` | SRH、Schottky、功函数及复合速度、多段输出 | E301 | 补声明后 `--check`；另有具名 JSON 重启/IV 验证，不能冒充新包 cmd 数值 E2E |
| `genius_bjt_sentaurus2022/bjt_m0_des.cmd` | SRH、双 QS、Plot/Time | E301 | 补声明后发射与阶段/时间语义；原生归档不等于 Vela 数值 E2E |
| `singledevice_sentaurus2018/singledevice_sdevice.cmd` | 未展开预处理／占位符 | E101/E102 | 正确拒绝 |
| `bvmethods_sentaurus2018/full_physics_ablation_sdevice.cmd` | 含白名单外 Physics | E107 | 正确拒绝；已识别子项不代表整份可运行 |
| `transportmodels_sentaurus2022/pp6_des.cmd` | Thermode/System 等超范围结构 | E105/E109 | 正确拒绝 |

缺 Parameter 的五份原文并非本轮新增故障：既定契约要求显式材料文件，测试仅在副本补声明。
本轮七个发射控制使用合成 4 节点/2 三角形网格及 `allow_lossy=True`；这不是七个真实器件
已经重新数值通过。必须保留材料有损账本，不能把该选项解释为允许替换不支持模型。

后续包门控还覆盖 `aba_coupled_sdevice.cmd`、`voltage_to_current_sdevice.cmd`；两者在上述
43 文件清单中有精确哈希和当前解析模型，已有包证据为 `--check`，不能升级为完整 BV 数值验收。
目录名的 `2018` 不用来推断当前原生引擎版本；PN2D v2 及后续原版参考须按各自清单绑定版本。

安装库原版是另一目标：[原版冻结报告](bvmethods_original_reference_q0_q1_2026-10-15.md)记录
T-2022.03-SP2、Workbench resistor 节点 5；cmd SHA-256 为
`30f025ed061f0c113cc4e1caa7c00f56003a6ac9ad01a4199e014b08602108a5`，par 为
`618d0467ee2608a76c7d3bf1dead41d5cc4709c73ad8dcfb59e341636c2630f3`，网格为
`7cd3968f5dd0e3cd63ef62ed0fcced465dd33bea465167f0be8139c16025902e`。
它启用 **NonLocalPath 与 GradQuasiFermi 雪崩**，不是原 fixture 的 E2/Eparallel。
本轮只核对其冻结报告，未重新访问原生归档或重跑。既有原生重复轨迹成功不等于 Vela 成功。

## 2. 端到端比较：已经验证什么，缺什么

重新读取 2026-10-07 Windows 包的三个 PN2D 通过报告及其真实 cmd、par、sidecar、网格资格，
以当前 parser/material emitter/deck emitter 重新生成 JSON。三个报告哈希均匹配封存账本；
未重跑 mesh importer，因此这一部分复用了已验收网格资格，不能认领当前 importer 数值回归。

| 比较对象 | 0 V | IV | BV | 解释 |
|---|---|---|---|---|
| 当前生成 stage 计划 vs 已验收计划 | 完全一致 | 完全一致 | 完全一致 | 阶段/目标/输出计划保持 |
| 当前材料 JSON vs 已验收材料 | 完全一致 | 完全一致 | 完全一致 | 使用相同 par 与材料身份 |
| 初始化求解器 | 有变化 | 有变化 | 有变化 | `gummel_newton` → `newton`，增加 `warm_start=true` |
| Poisson 初始化预算 | 有变化 | 有变化 | 有变化 | 新增 `poisson_max_iter=100` |
| OldSlotboom Fermi 修正键 | 显式 false | 显式 false | 显式 false | 与当前 C++ 缺省 false 相同，不据此认定物理变化 |
| 运行元数据 | 写入阶段差别 | 写入阶段差别 | 写入阶段差别 | `runtime_log.append`/`vtk_field_data` 由 orchestrator `_write_plan` 注入，不属于数值参数丢失 |

完整未过滤差分在 JSON 的 `pn2d_reemission`。初始方法和预算变化来自当前发射器对
Poisson＋Coupled 合并初始化的修改；即使初值更合理，仍须用原物理验收门证明输出和运行行为。
**当前 golden 通过只证明实现与当前 golden 相符，不能独立证明与已验收 JSON 行为相同。**

对支持范围内尚只有 `--check` 的案例，不应虚构一份“已验证完整 JSON”进行比较。
若下一版要承诺这些案例可完整计算，应逐个新增 cmd→JSON 与具名 JSON 直跑的同输入对照和数值
验收，分别记录前端映射差异与后端求解能力。不能用原生 Sentaurus 归档存在替代此证据。

现行数学字段也有明确边界：`Digits`、`RelErrControl`、电子/空穴 `ErrRef` 在发射器中记录 W201，
不自动替代 Vela 收敛规则；stage `Iterations` 映射 `solver.max_iter`。report 诊断及私有接受策略
研究不代表已经复刻原生终止判据。当前外电阻默认 `nested_scalar`，`coupled_newton` 是显式选择。

## 3. 成果可交付性

| 成果层级 | 当前实例 | 可认领范围 | 未完成条件 |
|---|---|---|---|
| 已封存发行包 | Windows/UCRT64 与 Ubuntu 26.04，运行时 `0a89942e` | 历史具名 PN2D E2E、门控、打包/并发证据 | 不能代表当前 dirty 源码 |
| 已有公共入口源码 | cmd parser、材料、发射、编排、PLT/VTU、schema、launcher | 可执行转换及门控；本轮回归通过 | 新初始化及其他运行时改动须筛选、重新构建验证 |
| 当前生产路径上的未提交修改 | Poisson 初始化、Fermi/NoFermi 映射、外电阻 continuation 等 | 代码确实进入公共入口，必须纳入差分审查 | 不能仅因源文件存在就视为已发布 |
| 实验 C++ 接口 | `vela::experimental` NonlocalPath 几何/声子内核/双源/实际残差诊断 | 有自主实现及分层实验，当前文件尚未入 Git | 公开 parser 仍 E106；完整源/Jacobian/接受未获资格 |
| 私有 R2 候选 | 同源高低位装配、迁移率精度链、增广缺陷校正 | 两态具名方向和隔离有限更新证据 | 当前状态重线性化、完整迭代、原版轨迹与 R3 L2 未完成 |
| 原生参考及诊断证据 | WP14 v4、原版轨迹、矩阵/RHS/PMI 采集 | 解释差异、建立输入身份和测试基准 | 不是可分发的 Vela 功能 |

最新 [R2 记录](bvmethods_resistor_r2_mobility_contract_2026-10-22.md)中的两态细档 Jv 通过有价值，
但全步回放仍拒绝，尚未提交状态持续迭代；它不阻碍旧具名包的使用，也不能宣布原版 BV 支持。
WP14 v4 瞬态获取＋稳态续扫与原始从零纯 QS 保持不同资格。

## 4. 本轮验证及证据局限

- 10 个原始 fixture：输入身份 10/10 相符；退出码 2 个 0、5 个 E301/3、3 个拒绝/1，符合契约。
- 七个可准备输入：14 个 stage deck 与现行 golden、七份 stage plan 全部相符。
- 十个前端 unittest 模块：189 项，178 通过、11 跳过。包括 parser、contracts、mesh、materials、
  deck emitter、orchestrator、acceptance、gate mapping、PLT、concurrency。
- 首轮有 12 个错误，包含 sandbox 默认临时目录权限失败以及隔离子进程的编码读取异常；
  设置工作区 TMP/TEMP 和 UCRT64 PATH 后，同一测试代码复验通过。未改断言或门限。
- 附加 27 项全部通过：两个真实 runner 合成器件/状态链冒烟、原版参考契约、PN2D 模板契约。
  runner 是现存 2026-10-20 二进制，哈希已登记，**不是本轮当前源码新构建**。
- 两个 2026-10-07 包均在本地，SHA-256 与已发布审查记录匹配；三个 Windows PN2D 报告也匹配。
  两平台完整数值通过依据为[封存包记录](sdevice_package_closeout_2026-10-07.md)，不是本轮新运行。
- 本轮未跑完整 CTest、原生数值 E2E、Linux 并发或新包。跳过项不算通过。

归档 PN2D 数值资格：1943 节点/3680 非钝角三角形；长度 µm、人口 cm⁻³、电流 A/µm、
电场 V/cm，实际 HDF5/UMFPACK 等后端身份见封存包记录。BV 十锚点包含 −19.5/−19.6/−19.7/−20 V；
−13.2 V 是场量检查点，不是该算例的击穿结论。本轮没有生成新的物理解。

## 5. 重新排序的任务与退出条件

| 顺序 | 任务 | 退出条件 | 是否阻碍前端下一版发布 |
|---|---|---|---|
| D0，本轮完成 | 冻结输入、验收层级与源码/包身份 | 本报告及 JSON；公开支持口径可追溯 | 已完成 |
| D1，下一项 | 审查当前源码相对封存运行时的前端及 C++ 依赖改动；首先处理 PN2D 初始化差异 | 明确每项保留/隔离理由、具体源码快照及依赖；保护全部无关工作 | 是 |
| D2 | 同一 cmd/par/grid/sidecar 的当前 cmd 链与具名 JSON 直跑对照；先零偏，再 IV/BV | 当前真实 runner 下语义、阶段、终端/场量及原门全部通过；差异按前/后端登记 | 是 |
| D3 | Windows/UCRT64 与 Ubuntu 26.04 重新打包、原生 PN2D E2E、check/rejection、迁移/并发、相对 RUNPATH | 两包各自绑定源码/运行时/输入哈希，完整相关回归；skip 不代 pass | 是 |
| D4 | 用户交付入口与支持矩阵收口 | 发布对象/命令/依赖/已知限制一致；旧报告不冒充新包 | 是 |
| 扩展 S | Schottky/BJT/ABA/电流目标等由门控资格升级为完整数值承诺 | 每例具名输入和前后端对照，另定覆盖与数值预算 | 只在承诺升级时阻碍 |
| 独立 N/R | 安装库原版 NonLocalPath、R2 增广 Newton、R3 严格重启/完整纯 QS | 按原研发门逐项通过，才接回公开 cmd 能力并申请发布资格 | 不作为原 PN2D/check-only 交付的前提 |

D1 不能简单把所有近期修复回滚；有些是通用正确性修复，有些只服务诊断，有些修改公共行为。
先分类及绑定证据，再选发布集合。原门限不因本次范围审查而降低。

文档缺口：fixture 索引的 `reference_artifacts` 仍含 WP0 时“未重跑”字样，虽有后来 v2/v4
记录；它是历史元数据，不应当当前状态。工作流对 Math 的笼统描述也需明确“记录”与“实际实现”。
本轮已补充工作流/计划入口及 Math 边界，旧 fixture/golden 和历史报告保留；后续 D4 统一版本化
状态来源，避免同时维护相互冲突的完成列表。

本轮可重查产物位于 ignored `build-release/frontend_scope_audit_20261022/`；
审查脚本、首轮/复验日志、现存 runner 和关键源码哈希在伴随 JSON 中登记。
该报告是本次范围决策依据，历史报告中的“下一步”不自动成为前端发布前置任务。
