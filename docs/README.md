# Vela Documentation

- [BVmethods 外接电阻法开发交接，2026-10-23](validation/bvmethods_external_resistor_handoff_2026-10-23.md)：
  当前检出与历史私有证据位置、R0–R7 状态、两态 R2 进展及持续增广 Newton 的下一轮入口；尚无完整 BV 资格。

- [sd-sim 入口与 CLI 更新计划](superpowers/plans/2026-10-23-sd-sim-cli-plan.md)：用户确定新名称，
  保留旧入口；未实现参数具名预留，`-d` 数值导出延期。

- [BVmethods 三组源码与诊断归档，2026-10-23](validation/bvmethods_worktree_organization_2026-10-23.md)：
  交付修复、实验依赖和诊断资料分别提交；原始执行记录保留 ignored，N4/R2/R3 未取得新资格。

- [SDevice 用户交付入口与支持矩阵](sdevice_delivery.md)：D0–D4 收口；发行包哈希、两平台命令、
  PN2D 数值资格及 check/拒绝边界；[D4 两机验证](validation/sdevice_delivery_d4_2026-10-23.md)。

- [SDevice 发布差分与重新验收，2026-10-22](validation/sdevice_delivery_revalidation_2026-10-22.md)：
  D1–D3 已完成（10-23 收口）；同输入三路数值对照、两平台全套回归和新包原生 E2E 通过。
  Linux 完整并发及相对 RUNPATH、Windows 锁冲突和两包迁移/指纹通过，平台跳过单列。

- [SDevice 前端交付范围与成果审查，2026-10-22](validation/sdevice_delivery_scope_audit_2026-10-22.md)：
  冻结 10 个原始 fixture，盘点 43 个参考 cmd；区分旧包数值资格、当前源码及独立 BVmethods 研发。
  同输入重新发射发现 PN2D 初始化流程已变，下一项先审查发布差分并重新验证主线，不能沿用旧包资格。

- [BVmethods R2 迁移率精度契约、完整算子及有限更新，2026-10-22](validation/bvmethods_resistor_r2_mobility_contract_2026-10-22.md)：
  两态 Newton 三个细档的完整 Jv 满足原门，50／80 位参考稳定，16 组合析因完成。
  实际写回回放仍拒绝全步；粗档曲率、事件侧及实际迭代资格仍单列，R2 继续。

- [BVmethods R2 迁移率取样与同源增广缺陷校正，2026-10-22](validation/bvmethods_resistor_r2_transport_defect_2026-10-22.md)：
  1453／1560 细档抖动收缩到迁移率取样；两态 high/low 同源校正消除线性精度平台。
  校正方向的第五点最细档电子 Jv 仍失败，场输入低位与真实 Newton 接受待补，R2 继续。

- [BVmethods R2 原始方向、边界增量与系数低位，2026-10-22](validation/bvmethods_resistor_r2_split_direction_2026-10-22.md)：
  保精度原始方向与约束行已绑定；high/low 实际矩阵算子误差约 4–6e-17。
  第五点最细档电子响应仍失败；尚未进入增广缺陷校正及 Newton 接受，R2 继续。

- [BVmethods R2 同源电荷与 SG 实际装配，2026-10-22](validation/bvmethods_resistor_r2_actual_chain_2026-10-22.md)：
  两态电子自由行与 drain 导数改善；实际写入方向下多数块满足原门。
  packed 写回、约束行、double 矩阵舍入及第 5 点粗档 Poisson 仍失败，R2 未收口。

- [BVmethods R2 人口电荷与 SG 方向，2026-10-22](validation/bvmethods_resistor_r2_population_sg_2026-10-22.md)：
  恢复 1266 电荷微响应；匹配实际写入方向后，236／241 主要 SG 边通过。
  同源结构化 Poisson 细档参考通过，原矩阵与完整电子块仍失败，R2 保持进行中。

- [BVmethods R2 同源残差增量与端口定位，2026-10-22](validation/bvmethods_resistor_r2_increment_2026-10-22.md)：
  实际装配操作数重放恢复小 PN 微源和空穴端口响应；两态空穴自由行方向通过。
  人口／电荷内部精度、电子完整方向与第 5 点 drain 236／241 仍阻断 R2，生产接受不变。

- [BVmethods R2 原生准则对照与线性精度修复，2026-10-22](validation/bvmethods_resistor_r2_native_2026-10-22.md)：
  区分原生更新误差、全局 L1 与逐行后向误差；默认关闭的增广缺陷校正恢复两态线性精度。
  六条 SRH 分项导数通过但完整残差量化仍失败，完整 Jv 和非线性资格未收口。

- [BVmethods 外电阻 R0–R3 实施与验证，2026-10-22](validation/bvmethods_external_resistor_r0_r3_2026-10-22.md)：
  取回旧扫描最终证据，接入增广迭代数驱动步长及停止电流；T470p 控制器回归通过。
  两态分量后向误差及六条 SRH 跨块检查未过原门，R2 与严格重启 L2 未收口。

- [BVmethods 外接电阻法开发计划，2026-10-22](superpowers/plans/2026-10-22-sentaurus-external-resistor-development-plan.md)：
  复用增广耦合求解器，补齐端口／接触导数、整组更新与迭代数驱动 QS；分开验收电路算法、生产接受和原生轨迹保真。R0–R3 已开始实施，R2 严格门及 R3 重启 L2 尚未收口。

- [BVmethods VM 从零完整 BV 扫描尝试](validation/bvmethods_vm_full_bv_attempt_2026-10-21.md)：
  最终证据已取回：完整扫描超时，主表有 34 条完整记录（含零偏）；同内压控制失败。后续分析见 2026-10-22 的 R0–R3 实施报告。

- [BVmethods VM 微小电荷、空穴跨块与分支事件验证](validation/bvmethods_vm_increment_event_validation_2026-10-21.md)：
  五态两批验证及 435 项证据检查通过；保精度参考恢复微小响应，接触允许方向通过 Jv，实际低位装配和 Newton 接受仍未收口。

- [BVmethods VM 符号一致的局部分支源导数](validation/bvmethods_vm_branch_source_2026-10-21.md)：
  六次运行与 178 项证据检查通过；五态电子 QF→空穴方向恢复，电荷微增量与弱响应仍待核验，未进入 Newton 接受。

- [BVmethods VM 229／1308 源链隔离控制](validation/bvmethods_vm_qf_source_chain_2026-10-21.md)：
  六个控制完成；定位边 366／4132 的跨零差分，首态细步长候选恢复局部 Jv，完整接受未收口。

- [BVmethods VM 雪崩扰动支撑，2026-10-21](validation/bvmethods_vm_avalanche_support_2026-10-21.md)：
  私有边外电势列修复消除首态局部 ψ 空穴方向的主要偏差；五态运行正常退出。
  后续控制已取回；Poisson 体积／常数显著降低残差，完整资格仍未通过。

- [BVmethods VM 私有 C++ 候选，2026-10-21](validation/bvmethods_vm_cpp_candidate_2026-10-21.md)：
  晚间恢复后私有 C++ 五态与重复完成，端电流差降至约 0.00018%–0.00334%。
  源项隔离将局部空穴方向偏差指向雪崩链；完整 RHS／全部 Jv 未通过。

- [BVmethods vela VM 材料筛选与单元迁移率，2026-10-21](validation/bvmethods_vm_transport_2026-10-21.md)：
  已迁移并独立重复五态诊断；Si 支撑加接触电势高场候选将漏极电子通量差降至约 0.0004%–0.0039%。
  仍为冻结 SG 因子控制，Linux C++ 实际装配及导数尚未获得资格。

- [BVmethods T470p 端口输运与 Poisson 分项，2026-10-21](validation/bvmethods_port_poisson_2026-10-21.md)：
  25 次诊断、166 项证据检查；几何组合大幅降低 Poisson 残差，三条边解释约 98.8% 漏极电子通量差。
  确认界面边含氧化层支撑；离线常数／权重控制未晋升，完整 RHS 与收敛规则对齐仍待完成。

- [BVmethods T470p 状态与接触映射，2026-10-21](validation/bvmethods_state_contact_mapping_2026-10-21.md)：
  私有材料／统计／热电压契约消除毫安级投影异常；五点端电流仍差约 7.24%–7.26%，完整 RHS 与负载未通过。

- [BVmethods T470p 可比性诊断与隔离回放，2026-10-21](validation/bvmethods_comparable_replay_2026-10-21.md)：
  五组固定内压人口投影、冻结分母、增广负载及无扰动回放；投影端量不匹配，严格原生身份未通过。

- [BVmethods 原生相对残差与过程可比较性，2026-10-21](validation/bvmethods_native_comparability_2026-10-21.md)：
  补齐 60 项 RHS/实际残差绑定和 20 项人口转换；首点接触行支配初始残差。
  固定内压与增广 QS 的分母、QF 更新和试步需分别绑定，生产接受策略不变。

- [BVmethods 原生五态完整收敛规则复核，2026-10-21](validation/bvmethods_native_convergence_2026-10-21.md)：
  更正 TCAD/SI 行尺度混用：五态局部行门和电子全局闭合均失败；空穴闭合、输出负载／KCL通过。
  完整相对残差和过程性接受仍未绑定，不把输出电流与装配通量混用。

- [BVmethods 原生前五步与载流子行门，2026-10-21](validation/bvmethods_native_row_gate_2026-10-21.md)：
  15 套原生线性系统及五个实际接受态的通量／源支撑；单位更正后五步均超出 0.001 行门。
  修订该门的原生兼容性解释；保留端量身份差异，生产接受策略未替换。

- [BVmethods 电子跨块与电荷增量私有修复，2026-10-21](validation/bvmethods_crossblock_charge_fix_2026-10-21.md)：
  T470p 两批重复验证及操作数定位；修正私有零支撑采集后，电子 706／721 完整方向三档通过 1e-5，
  Poisson 电荷响应恢复。保留首批失败证据，完整空穴相消与生产资格仍未收口。

- [BVmethods 私有解析 Jacobian 与跨块诊断，2026-10-21](validation/bvmethods_private_analytic_jacobian_2026-10-21.md)：
  T470p 十二次核验；717／722 实际 SG 方向三档通过原 1e-5 门。
  电子跨块偏差仍在；Poisson 空穴电荷分项正确、合成残差丢失增量，完整 Jv／N4／纯 QS 未收口。

- [BVmethods 门限归属与实际完整 Jv，2026-10-20](validation/bvmethods_full_jv_contract_2026-10-20.md)：
  恢复 N0 的 1e-5 方向契约，保留 1e-8 RHS 门及额外导数诊断。T470p 十二次核验；
  原 SG 的 717 仍约 2.28e-5 超门，完整方向相消及电子跨块／Poisson 电荷分辨问题继续保留。

- [BVmethods 支撑边迁移率链与场增量，2026-10-20](validation/bvmethods_mobility_chain_2026-10-20.md)：
  T470p 四次诊断；低场／表面／HFS 原语及保精度几何场反馈通过独立检查。
  识别粗档截断与 double 差分舍入，完整 SG 行三档未达到额外 1e-8 诊断目标；正式方向契约为 1e-5，接受策略不变。

- [BVmethods 717 逐边操作数与解析 SG 候选，2026-10-20](validation/bvmethods_sg_analytic_candidate_2026-10-20.md)：
  T470p 八次诊断；717 固定部分偏差的 97.76% 定位到 717–718 的空穴 QF 差分列。
  高精度解析偏导 600 项通过，但完整方向未达到额外 1e-8 诊断目标，候选未采用。

- [BVmethods 隧穿独立参考与 SG／Fermi 分离，2026-10-20](validation/bvmethods_sg_fermi_reference_2026-10-20.md)：
  双机声子内核 16 例／212 断言通过；T470p 四次诊断分离 722 的反馈相消与 717 的固定迁移率
  Jacobian 偏差。采样人口链通过独立级数对照，完整 SG 导数未达到额外 1e-8 诊断目标，N4／QS 未收口。

- [BVmethods 空穴同方向分项导数核验，2026-10-20](validation/bvmethods_hole_direction_components_2026-10-20.md)：
  T470p 八次诊断确认 722／717 全步余差的 SG 占比超过 99.6%；SG 导数存在尺度敏感偏差，
  尚未达到额外 1e-8 导数目标。原轨迹逐位不变，未修改接受策略或采用候选方向。

- [BVmethods 实际写回、舍入余量与非线性阻断，2026-10-20](validation/bvmethods_writeback_constraints_2026-10-20.md)：
  T470p 四次运行；112 项身份检查、184 项精确审计及 330 项独立核验通过。
  最大空穴超界主要位于写回后的非线性余差；电子低位及固定边界控制未解除阻断，原门限不变。

- [BVmethods 退化工作集修复与可表示性核验，2026-10-20](validation/bvmethods_degenerate_working_set_2026-10-20.md)：
  T470p 三臂各两次实际装配，7 个数学控制、104 项远端及 169 项独立核验通过。
  高精度投影解除工作集停滞；精确主变量方向仍有亚 ULP 越界，实际试步未采用，原门限不变。

- [BVmethods 硬约束秩、电子块结构与投影对照，2026-10-20](validation/bvmethods_hard_rank_structure_2026-10-20.md)：
  T470p 四次运行、53 项远端和 154 项独立核验通过。活动集停于近相关阻断约束；
  电子块不满足零 RHS／无外耦合消元条件。投影下降候选仍超限，未采用，原门限不变。

- [BVmethods 独立低位 SG 与硬约束方向，2026-10-20](validation/bvmethods_low_hard_constraints_2026-10-20.md)：
  T470p 两批 8 次运行、124 项远端和 138 项独立检查通过。独立低位 SG 恢复 rows5／rows0
  的试步下降；完整低位本构尚未贯通，硬约束优化未获最优性及非线性采用资格，原门限不变。

- [BVmethods 低位入口与固定资格受限方向，2026-10-20](validation/bvmethods_tail_fixed_rows_2026-10-20.md)：
  T470p 8 次运行、148 项远端和 212 项独立检查通过。电子低位转入 reference 再次舍失；
  零方向候选恢复全步下降且违规 285→281，但仍恶化固定资格行，未采用，生产门限不变。

- [BVmethods 全步保护、写回与下降方向，2026-10-20](validation/bvmethods_full_step_protection_2026-10-20.md)：
  T470p 四态 8 次运行、236 项远端及 702 项独立检查通过。rows5／21／0 分离电子 QF 写回损失；
  rows285 的 134 个空穴分量截断反转线性下降方向，主要不利支撑为 831／714／700／655。原门限不变。

- [BVmethods 参考 Charon／DEVSIM 的终止判据实验，2026-10-20](validation/bvmethods_alternate_stopping_2026-10-20.md)：
  T470p 六态七臂共 42 次运行，三种替代终止策略均未改善收敛；完整结果／轨迹与基线一致。
  Charon 风格受 Alpha=1 条件阻断，DEVSIM 风格受空穴相对更新阻断；生产默认和原门限不变。

- [BVmethods 保精度写回与受限耦合方向，2026-10-20](validation/bvmethods_bounded_writeback_2026-10-20.md)：
  T470p 144 项分析检查、端量 6 例／174 断言及标准 runner 控制通过。单点输出可显式指定接触；
  电子 QF 写回损失已分离；受限方向降低整体 merit 但恶化固定资格弱行，尚不具备采用资格。

- [BVmethods 人口映射、实际写回与拒绝试步，2026-10-20](validation/bvmethods_update_mapping_trials_2026-10-20.md)：
  T470p 四次诊断与封存回放一致，53 项分析断言通过。rows5 分离原始人口放大与写回损失；
  rows285 的全部 13 档受保护方向已在线性预测中失去整体下降；通用端量误选有效电势最大的 gate，须明确 drain 身份。

- [BVmethods 装配散射、缓存历史与隔离回放，2026-10-19](validation/bvmethods_update_error_scatter_cache_2026-10-19.md)：
  双机 12 个 checkpoint 的实际散射逐位还原 Jacobian，162 个共同迭代数值一致；
  缓存差异明确为冷启动结构/符号分析及末态汇总访问，完整缓存身份仍不整体放行。
  补齐 172 份当前行账本、4 份最佳被拒步误差和 12 次隔离一步回放；六态仍未收敛，接受策略不变。

This directory is organized around current project behavior: build it, choose
a configuration schema, import reference fixtures, and inspect validation
evidence.

## Current References

- [BVmethods 独立 J/R 与完整结果位模式，2026-10-19](validation/bvmethods_update_error_independent_audit_2026-10-19.md)：
  双机六态 84 次运行、924 项核对通过；独立 off 侧证据覆盖 172 个迭代位置，
  完整 NewtonResult、J/R、方向与线搜索在七臂控制中一致，六态仍未收敛。
  本轮补齐限定诊断链路的无扰动证据；严格逐元素/缓存契约、行账本和生产 report 接入仍待完成。

- [BVmethods 实际人口与精确 checkpoint，2026-10-19](validation/bvmethods_update_error_population_checkpoint_2026-10-19.md)：
  新增解析 Jacobian 实际人口被动采集及实验性人口保真 sidecar，逐轮报告和重启共同区间扩展到六态。
  双机 172 轮人口对照、12 个精确 checkpoint 和 162 个 L2 事件匹配；六态仍未收敛。
  普通 SI 存档和生产接受策略不变；逐元素/缓存、独立 off 侧 J/R 与完整无扰动资格仍需补齐。

- [BVmethods 被动 Newton 采集与原位重启比较，2026-10-19](validation/bvmethods_update_error_capture_restart_2026-10-19.md)：
  新增默认关闭的有界采集入口；Newton 144 例／3232 断言及 report 8 例／59 断言通过。
  两机 21 行态连续 19 轮重启对照一致；保存人口末位、严格 L1 与完整 report 无扰动资格仍保留。

- [BVmethods report 计数、人口方向与跨机器舍入，2026-10-19](validation/bvmethods_update_error_crosshost_population_2026-10-19.md)：
  独立未知量／统计条目及分子／分母来源分开，双机 7 例／51 断言通过；六态人口方向通过 1e-8 检查。
  私有 FMA 控制解释部分跨机差异；严格 L1/L2、生产 report 接入及无扰动硬门仍未完成。

- [BVmethods 原生更新误差采集与 T470p 复核，2026-10-19](validation/bvmethods_update_error_native_collection_2026-10-19.md)：
  两组原生 65 点轨迹逐字节一致；原始更新分子＋保护后分母复现 412 轮误差，统计条目与线性未知量计数分离。
  T470p 六态同机复现及 31 断言通过；跨机器 J 摘要差异、严格 L1/L2 与完整 report 无扰动仍待处理。

- [BVmethods 更新误差定义与 L1 检查，2026-10-19](validation/bvmethods_update_error_s0_l1_2026-10-19.md)：
  阶段 0 现有资料已登记；六态重复与往返的完整 J/R、重算人口、活动掩码和资格集合一致。
  存档人口严格逐位契约及原位历史仍待补齐；仅新增独立 RMS 报告原语，未改变生产接受策略。
  近期顺序以[只报告诊断计划](superpowers/plans/2026-10-19-bvmethods-update-error-report-plan.md)为准。

- [BVmethods动态源资格与非线性邻域，2026-10-19](validation/bvmethods_dynamic_rows_2026-10-19.md)：
  资格不变的小范围及双节点控制使违规21→11→10；有限分支的0违规候选仍有退出资格后的失衡行。
  T470p五组续算均未收敛，明确保留源资格退出账本、完整残差及耦合方向检查，未改生产门限。

- [BVmethods受限耦合修正与耗尽分支，2026-10-19](validation/bvmethods_constrained_correction_2026-10-19.md)：
  候选同时改善718及全局闭合，但新增行违规，未接入生产；918存在正向有限分支零点。
  非接受恢复种子在T470p推进25轮后剩21条源主导空穴行，原QS路径及N4仍未获资格。

- [BVmethods真实阻断行与截断方向，2026-10-19](validation/bvmethods_blocking_rows_2026-10-19.md)：
  723场导数误差定位到固定差分尺度；718等逐分量截断反转全局闭合方向。
  285行中279空穴／6电子，统一缩放受918巨大更新限制；七次续算均未收敛，门限不变。

- [BVmethods参考SLOT的弱行分类，2026-10-19](validation/bvmethods_weak_row_policy_classification_2026-10-19.md)：
  区分尺度下限支配、真实不平衡与亚ULP诊断；当前285条阻断行无下限支配行，三条亚ULP电子行并非当前阻碍。
  后续优先真实阻断行及分载流子场导数／截断，不将所有弱行一概放行。

- [BVmethods MinStep落点修复与弱行方向，2026-10-19](validation/bvmethods_min_step_weak_direction_2026-10-19.md)：
  修复减半越过最小步时漏试，115例4400断言通过；T470p两条限时扫描获得非零接受态，原首目标仍未获资格。
  定位亚ULP更新及电子／空穴场导数链的不同响应；生产门限和默认重试数不变。

- [BVmethods 无 Band2Band：零偏几何、弱行与外电阻定位，2026-10-19](validation/bvmethods_no_btbt_diagnosis_2026-10-19.md)：
  T470p 15 个控制定位盒体积／钝角转移及独立场恢复差异；显式几何候选零偏 ψ/n/p 通过，原首步仍失败。
  恰好 10 µV 外步可接受，默认减半未试到 MinStep；保留更新截断、弱行舍入及全部原门限。

- [BVmethods 双方关闭整条 Band2Band 对照，2026-10-19](validation/bvmethods_no_btbt_comparison_2026-10-19.md)：
  原生从零纯 QS 单次轨迹通过；T470p 当前 Release 的零偏场量不通过严格门，首点及完整扫描均无非零接受态。
  6 条失败试探一致，保留数学映射／重试差异；原版目标、N4 和既有门限不变。

- [BVmethods N4统计舍入、保精度源装配与一致导数，2026-10-18](validation/bvmethods_n4_precision_pipeline_2026-10-18.md)：
  SG统计链分阶段31／33／35条冲突，实际36条共同门冲突保留；未调整验收门限。
  NonLocalPath实际保精度双源在匹配电子方向三档1639条非零行全部通过，158条内部响应差异消除。
  雪崩原50条在最细实际增量中通过，全自由行仍455条扩展指数范围失败；N4继续，尚无完整耦合资格。

- [BVmethods N4 SG共同门、雪崩微小响应与完整NonLocalPath方向，2026-10-18](validation/bvmethods_n4_sg_avalanche_full_jv_2026-10-18.md)：
  36条SG冲突分解为参数舍入35条及乘积舍入新增1条；原失败保留。
  雪崩原50条在独立保精度增量／细分控制中通过，1935旧导数差定位到指数次正规数舍入。
  全4156种子自主NonLocalPath三方向通过，实际C++仍160条、同投影控制仍158条失败；N4继续，门限未改。

- [BVmethods N4同态实际人口方向与匹配Jv，2026-10-18](validation/bvmethods_n4_matched_population_jv_2026-10-18.md)：
  三态实际人口精确绑定，41条目标输运Jv通过；自主细分后全2121行输运通过，35条单侧余差分离为曲率与事件。
  同态SG共同门仍有36条区间不相交，三态38/39/38条目标绝对失败保留；完整NonLocalPath Jv未获资格。
  6次Release私有求值、两平台各10项测试、3项CTest及35项证据检查通过，N4继续。

- [BVmethods N4原生非零状态与SG/RHS绑定，2026-10-18](validation/bvmethods_n4_native_sg_binding_2026-10-18.md)：
  九次原生控制完成，三态各273条SG支撑重放一致，41条目标行RHS账本误差不超过8.881e-17。
  确认边长除／乘往返产生1 ULP差；正负细步各41条实际人口方向未匹配意图，负向2156/2157量化为零。
  本批Load态不同于旧正QS RHS2，旧区间与95条失败保留；两平台各7项测试、2项CTest及46项证据检查通过，N4继续。

- [BVmethods N4 SG当前值/增量与参考态契约，2026-10-18](validation/bvmethods_n4_sg_current_contract_2026-10-18.md)：
  21配对态独立SG重放；固定参考补偿在电子方向三档失败124/116/118行，候选拒绝。
  最细档41条SG行的兼容绝对值/Jv允许区间不相交，尚需同态原生扰动确认；95条基态失败保留。
  两平台各8项回归、Windows4项CTest及25项证据检查通过；N4继续，未改门限或进入Newton。

- [BVmethods N4配对主状态实际装配与接口舍入，2026-10-18](validation/bvmethods_n4_paired_live_assembly_2026-10-18.md)：
  21态贯通当前人口/QF、迁移率、SG/SRH及Poisson；864实际方向约6.48e-15，ψ/p局部分项三档通过，95条SG绝对失败保留。
  修复私有NonLocalPath接口非最近舍入，两平台38238值一致及三态实际复核通过；double小增量分辨率仍不足。
  两平台各8例66断言、另20例149断言及64项证据检查通过；完整Jv/接触耦合和Newton未获资格，N4继续。

- [BVmethods N4主变量配对、QF／迁移率与事件侧，2026-10-18](validation/bvmethods_n4_primary_pair_mobility_2026-10-18.md)：
  私有配对人口将SRH三档意图方向120／144／178条失败降为零；SG固定系数方向全部通过，95条绝对失败仍保留。
  864同侧余差定位到绝对QF舍入，C++局部候选三档约5.98e-14；4701有限跳变及两侧控制完成。
  七态完整观察器逐值不变，两平台原语／重放、34项证据检查通过；第一项与N4继续，N5未开始。

- [BVmethods N4 SG／SRH保精度实际装配，2026-10-18](validation/bvmethods_n4_sg_srh_preserved_assembly_2026-10-18.md)：
  SRH低位贯穿实际私有速率／积分／逐行残差，电子方向三档192／230／259条实际增量失败降为零；完整基准RHS仍通过。
  95条SG绝对失败均定位为电子输运；全宽视图导数改善但绝对门仍拒绝，迁移率4701事件侧及864同侧链保留。
  主变量意图量化仍单列，21次完整装配、两平台性质测试及证据审计完成；第一项未收口，N4继续。

- [BVmethods N4保精度表示、小步弱行及五态路径链，2026-10-18](validation/bvmethods_n4_reference_increment_states_2026-10-18.md)：
  三条弱行小步实际源方向通过；五态20780个种子全部解析，45项非零数学方向通过，三类导数契约的27条远端对照通过。
  私有谱对数形式精度候选消除四条源响应失败，真实正QS完整RHS仍通过原门；默认分支不变。
  24项证据检查及两平台表示契约通过；SG／SRH全装配及五态完整Jacobian仍待验证，N4继续。

- [BVmethods N4近零驱动、SRH精度与完整路径链，2026-10-18](validation/bvmethods_n4_drive_path_chain_2026-10-18.md)：
  全部4156个Si种子的变化路径／谱／收集链及远端导数接入私有6843阶候选，含146123个非零项。
  局部电势完整响应缺项消除；SRH近零统计及最终速率输出量化分离，SG全宽精度候选95行失败保留。
  公共方向弱行仍未通过，生产配置未改；N4继续、N5与Newton接受循环未开始。

- [BVmethods N4自主主变量Jacobian及精度阻断，2026-10-18](validation/bvmethods_n4_primary_jacobian_2026-10-18.md)：
  统计／迁移率／SRH／雪崩与Poisson、局部负载组成6843阶部分矩阵；局部电势三档输运和雪崩控制通过。
  12个Catch2用例123断言、19项证据检查通过；近零源和公共电势方向仍有精度失败，NonLocalPath完整路径导数未接齐。
  完整Jacobian／Jv／自主线性方向未取得资格，尚未进入Newton循环；N4继续、N5未开始。

- [BVmethods N4实际SG增量、迁移率分支与变化路径谱原语，2026-10-18](validation/bvmethods_n4_sg_variable_path_2026-10-18.md)：
  10条粗档差异均涉及高场迁移率1 V/cm分支切换；系数不变三档的SG原语方向最大差9.532e-14。
  谱原语支持当前变化路径独立求值，三档逐行方向收敛通过；8个Catch2用例66断言、5项谱导数性质及18项证据检查通过。
  原语偏导尚未组成完整主变量Jacobian，未进入真实Newton接受循环；N4继续、N5未开始。

- [BVmethods N4谱增量、输运分项及RHS一致方向候选，2026-10-18](validation/bvmethods_n4_kernel_linear_candidate_2026-10-18.md)：
  2693固定几何保精度增量使小步长源误差降至1e-9至1e-11量级；n块10条粗档差异的变化主要来自输运。
  独立缩放／迭代改进线性控制通过1e-10；候选实际RHS方向改善三个残差块、负载及KCL，但部分弱行仍变差。
  参考矩阵方向不等于自主Jacobian，固定几何控制不等于生产内核；N4继续，完整Jv／Newton尚未取得资格。

- [BVmethods N4源精度、分块方向、原生边界及耦合，2026-10-18](validation/bvmethods_n4_ordered_closure_2026-10-18.md)：
  1317源响应定位到2693，三档可分辨方向通过；21次分块增广装配闭合通过，但分块Jv仍有超门／未分辨行。
  三个新增原生控制的同态双源通过；完整Newton步满足负载与KCL门却未共同改善载流子残差，打印线性系统7行超门保留。
  16项证据检查及六组旧封存校验通过；N4继续，完整Jv／Newton及生产晋升未通过。

- [BVmethods N4完整增广增量、载流子支撑与收集边界，2026-10-17](validation/bvmethods_augmented_increment_support_2026-10-17.md)：
  五类方向35次完整增广状态求值通过人口、源守恒及6843条残差分项闭合；基准对原生RHS最大差1.781e-9。
  ψ／n／p选定支撑与独立混合方向的可分辨响应相符，负载列全行通过；远端空穴1317两档源导数失败和Newton相消／未分辨结果保留。
  收集边界两侧各三档控制稳定，20项证据检查及五组旧封存校验通过；完整Jacobian／Jv及Newton未取得资格，N4继续、N5未开始。

- [BVmethods N4保精度边增量、同侧雪崩与正QS导数契约，2026-10-17](validation/bvmethods_edge_increment_derivative_contracts_2026-10-17.md)：
  三条电子行的SG边输入不变，微小FD差定位到源／全局累加舍入；九组保精度增量方向通过原1e-5门，最大差1.952e-8。
  七个同侧原生态的雪崩源及12组重点源导数通过；真实正QS中显式非局部导数改变71145个矩阵条目，三个较小ψ步的11条远端源行通过。
  粗档收集事件及QF重新定基失败保留；完整增广Jacobian／Jv、Newton和生产晋升仍未完成，N4继续、N5未开始。

- [BVmethods N4原生非零RHS、跨平台数学与增广方向，2026-10-17](validation/bvmethods_native_perturbed_rhs_math_2026-10-17.md)：
  14个同态原生控制、95802条独立方程通过原1e-8门，最大差2.852e-10；实际双源账本通过固定规则。
  MPFR私有标量内核的6009组位模式在Windows／Ubuntu 26.04一致，两平台各4例10197断言通过。
  三档导数及16次完整增广方向已检查；输运微差、雪崩单侧分支和Newton负载响应仍待验证，N4继续、N5未开始。

- [BVmethods N4可变化状态与增广负载，2026-10-17](validation/bvmethods_live_state_augmented_load_2026-10-17.md)：
  四类物理项共享可变化人口／QF状态，移除固定参数数学表；参考GNU通用数学下五态34215条独立方程通过原1e-8门，最大差约2.86e-10。
  两数学环境62次装配、独立人口检查及真实增广负载方向通过；Ubuntu默认数学的零偏10条弱行失败保留。
  非零试算的装配闭合不替代同态原生RHS对照；完整Jacobian／Newton及生产晋升未完成，N4继续、N5未开始。

- [BVmethods QF参考、严格人口与完整RHS映射，2026-10-17](validation/bvmethods_qf_reference_full_rhs_2026-10-17.md)：
  原始归一化QF与自主参考带状态接入私有实际装配，五态通过原128 epsilon人口入口；每态6843条独立方程的固定态映射通过原1e-8逐行门。
  80组非零人口方向、10个成功原生采集、15次装配、4个C++用例921断言及18项证据检查通过。
  SG／源仍含冻结数学控制，任意扰动的完整算子和真实增广Newton尚待验证；生产配置未晋升，N4继续、N5未开始。

- [BVmethods归一化SRH与雪崩SG支撑，2026-10-17](validation/bvmethods_srh_avalanche_alignment_2026-10-17.md)：
  自主归一化SRH及与连续性相同的SG雪崩电流接入私有实际装配；GNU/Linux冻结数学控制下，五态全部载流子自由行通过原1e-8门，最大误差约2.86e-10。
  4个原生采集、22次成功装配、5个C++用例214断言和18项证据检查通过；本机近零舍入差保留。
  严格QF人口、完整RHS、增广负载及任意扰动资格仍待补齐；候选未晋升生产，N4继续、N5未开始。

- [BVmethods源分量与严格人口精度，2026-10-17](validation/bvmethods_source_components_precision_2026-10-17.md)：
  原生SRH输入／输出、三种源数组及实际盒积分绑定同一RHS1；首点和零偏全部自由行的剩余差定位到源项。
  首点85条空穴行为SRH差超预算，另2条为雪崩弱行；零偏另有4条空穴雪崩重点行，NonLocalPath差在预算内。
  8个原生采集、五态观察装配、21520组100位统计反算及19项检查完成；SRH最终相减精确，需追查相消前操作数，N4继续、N5未开始。

- [BVmethods全自由行、SG求值与源项支撑，2026-10-17](validation/bvmethods_free_row_sg_source_2026-10-17.md)：
  保留归一化主电势并对齐Bernoulli分支／求值顺序，五态25次装配；BV三态全部载流子自由行通过原门限，首点电子在自主VM数学控制下通过。
  72项原生SG基础电流观察值逐值相同；首点87条空穴行和零偏仍失败，重点2423／2432已定位为净源项主导。
  4个C++用例4457断言、20项证据检查通过；严格人口、完整残差与增广负载保留，N4继续、N5未开始。

- [BVmethods冻结主人口与BGN顺序，2026-10-17](validation/bvmethods_population_bgn_contract_2026-10-17.md)：
  显式归一化SG主人口接口、BGN归一化／求和及平台舍入控制完成，五态35次实际装配；首点空穴2515及击穿后电子2156失败消除。
  12个C++用例10880断言、16项证据检查通过；BV三态电子仍有37／38／17条自由行失败，严格人口、完整残差与增广负载尚未取得资格，N4继续。

- [BVmethods自主SG带状态与五态弱行，2026-10-17](validation/bvmethods_autonomous_sg_band_chain_2026-10-17.md)：
  自主DOS／默认统计／Bernoulli及原版亲和势温度分量接入私有实际装配，首点电子2156由3.62e-6降至1.62e-16；五态25次控制完成。
  8个C++用例433断言、22项证据检查通过；首点空穴、击穿后人口表示及其他自由行仍失败，N4继续、N5未开始。

- [BVmethods首点近零SG与绝对输运支撑，2026-10-17](validation/bvmethods_nearzero_sg_support_2026-10-17.md)：
  24组SG操作数、两个弱行的实际单元／逐边支撑已绑定RHS1；双精度导入控制误差低于4e-16，九条非零边分别闭合。
  6个C++用例138断言、12项证据检查通过；自主有效参数生成、五态弱行、严格人口和负载耦合仍待验，N4继续。

- [BVmethods接触单元雪崩驱动力修复，2026-10-17](validation/bvmethods_contact_avalanche_force_2026-10-17.md)：
  私有源补齐已有接触电场回退；三态BV雪崩源差低于1e-15，末态2156空穴弱行由约1.37%降至1.44e-14。
  6个C++用例109断言、8项采集契约检查通过；近零弱行、严格人口及负载耦合仍待验，N4继续。

- [BVmethods自主归一化接入实际装配，2026-10-17](validation/bvmethods_normalized_geometry_actual_assembly_2026-10-17.md)：
  独立参考长度及原始TDR转换顺序复现5356个内部坐标分量、7619条边长度和4156个Si代表身份；自主规则接回私有实际残差，五态三臂完成。
  7个C++用例25断言、16项诊断检查通过；局部空穴源、近零弱行、严格人口及负载耦合仍待验，N4继续、N5未开始。

- [BVmethods代表顶点容差与原生顺序核对，2026-10-17](validation/bvmethods_representative_order_audit_2026-10-17.md)：
  识别内部归一化边方向行列式的严格比较及代表索引交换；原版固定态重绑4156个Si身份，150个并列单元全部一致，7619条边长度逐值相同。
  15项诊断检查通过；自主长度归一化及实际装配接入仍待验，N4继续、N5未开始。

- [BVmethods N4兼容统计、实际装配与逐弱行推进，2026-10-16](validation/bvmethods_nonlocalpath_n4_actual_assembly_2026-10-16.md)：
  修复冻结SRH人口契约、Poisson冻结电荷／绝对边项及端点缩放舍入；114个C++用例、7项CTest和16项诊断检查通过。
  私有统计／盒支撑／Gamma输运／局部迁移率与雪崩候选接入实际残差，71次固定态控制完成。
  三态BV整体差异约1e-8，但局部空穴源、近零弱行、严格人口身份及并列角规则仍未通过；N4继续、N5未开始。

- [BVmethods N4本构人口、RHS分项与负载行定位，2026-10-16](validation/bvmethods_nonlocalpath_n4_rhs_diagnosis_2026-10-16.md)：
  原版默认Joyce–Dixon统计与BGN链已独立绑定；定位SRH重心体积和冻结人口实际／分项契约差。
  修正私有探针载流子RHS单位，分离总雪崩源，确认一个原生电路自由度；二十次VM短控制、
  本机十次探针、十五项诊断检查及两个CTest通过。实际完整RHS、逐弱行和Vela负载耦合仍待验收。

- [BVmethods N4 C++单元源与实际残差检查，2026-10-16](validation/bvmethods_nonlocalpath_n4_cpp_closure_2026-10-16.md)：
  自主C++覆盖五态20780个Si单元状态实例，源字段门通过；本机12用例1331断言及10项CTest通过。
  T470p同源构建、五态复核及全量CTest通过（1052通过，4项E2E条件跳过）。
  实际双源入口及温度／极小事件保护已验证；本构人口、完整RHS及逐弱行仍待对齐，N4继续、N5未开始。

- [BVmethods N4自主五态全硅源与残差转移，2026-10-16](validation/bvmethods_nonlocalpath_n4_full_source_2026-10-16.md)：
  私有Python原型覆盖五态各4156个Si单元；击穿附近全源加权L1低于4.2e-8，源字段门通过。
  本机／VM各38项测试、本机10项CTest通过；C++单元源迁入和完整实际RHS仍待验收，N4继续。

- [BVmethods N4非均匀标量内核与五态，2026-10-16](validation/bvmethods_nonlocalpath_n4_nonuniform_kernel_2026-10-16.md)：
  定位带内位置插值及节点系数取样顺序；五态18条路径、608次标量对照，击穿附近选定单元双源L1低于2.3e-14。
  两机30项测试、本机10项CTest通过；零偏极小源保留舍入差，全源／完整RHS待验收，N4继续。

- [BVmethods N4递归分段与能量区间收集，2026-10-16](validation/bvmethods_nonlocalpath_n4_recursive_collection_2026-10-16.md)：
  595次标量求值、1190组两端系数及291个递归节点对照；自主均匀场340单元双源误差约2e-14，通过该态源门。
  两机24项测试及本机10项CTest通过；非均匀标量内核、其余状态和完整RHS仍待办，N4继续。

- [BVmethods N4能量窗与路径点收集，2026-10-16](validation/bvmethods_nonlocalpath_n4_energy_collection_2026-10-16.md)：
  默认无导数分支的窗口／带状态已绑定，13单元26次两端线性映射误差低于2e-15；原生系数尚未自主生成。
  单点电子收集的一阶矩反例已分离；本机10项CTest及T470p18项测试／重放通过，继续N4求积接受与系数生成。

- [BVmethods N4原生单元射线与源收集，2026-10-16](validation/bvmethods_nonlocalpath_n4_collection_2026-10-16.md)：
  原生每单元重心射线与当前顶点候选有差异；采齐均匀场340次双源增量并绑定真实RHS1。
  独立射线积分总源差约9.7e-6，节点分配未闭合；本机10项CTest及T470p重放通过，N4继续、N5未开始。

- [BVmethods NonLocalPath N4双源装配与差异，2026-10-16](validation/bvmethods_nonlocalpath_n4_2026-10-16.md)：
  私有双源、接触账本及实际残差诊断已接入，修复累计顶点舍入；N4原生门未通过。
  节点起点候选在BV附近源积分差约19%；独立盒积分缩小均匀场差异但仍未过门，继续N4收集／窗口对照。

- [BVmethods NonLocalPath N3声子辅助内核，2026-10-16](validation/bvmethods_nonlocalpath_n3_2026-10-16.md)：
  私有Si WKB／前因子、占据差及显式有限窗口完成；本机／T470p Release均通过11用例108断言及8项CTest，
  独立50位求积最大源误差约7.11e-15。后续N4实际双源／窗口对照进行中；公开前端仍E106。

- [BVmethods NonLocalPath N2自主几何，2026-10-16](validation/bvmethods_nonlocalpath_n2_2026-10-16.md)：
  路径起点／方向、穿越、端点、反射、选区和审计已实现；本机／T470p Release通过14项几何测试、
  7项CTest及五态逐路径字节一致检查。两处接触／反射角点明确保留待定；后续N3内核已完成，公开前端仍E106。

- [BVmethods NonLocalPath N1接口与单位，2026-10-16](validation/bvmethods_nonlocalpath_n1_2026-10-16.md)：
  私有配置、带状态／导数、路径、双源账本及远程稀疏支撑完成；本机与T470p Release均通过
  11项接口测试／97断言、3项E2测试及6项CTest。后续N2几何已完成；公开前端仍报E106。

- [BVmethods NonLocalPath N0契约与原生控制，2026-10-15](validation/bvmethods_nonlocalpath_n0_2026-10-15.md)：
  65状态双源守恒，真实400K确认Ppath能量作用，QS同态导数／策略对照及单位／状态工具完成；
  [v1验收预算](validation/contracts/bvmethods_nonlocalpath_n0.v1.json)已冻结。后续N1接口已完成；生产模型仍未实现。

- [NonLocalPath BTBT 开发子计划，2026-10-15](superpowers/plans/2026-10-15-nonlocalpath-btbt-development-plan.md)：
  对照2022手册、原生默认与一手文献，细化原版Q2/Q3为N0–N8；先核对双源／版本契约，
  再实现动态路径、Si声子辅助、实际源装配与非局部导数，最后进入全轨迹及两平台验收。
  N0原生控制／契约与数值预算、N1接口、N2几何及N3私有内核门已完成；N4自主全源字段门通过，C++实际装配／完整RHS仍待验收；生产模型仍未实现。

- [BVmethods 安装库原版 Q0/Q1，2026-10-15](validation/bvmethods_original_reference_q0_q1_2026-10-15.md)：
  已恢复并完成原件／三材料默认冻结、专用校验和原生从零复核，两次66点PLT字节一致；原生提取BV为6.67V。
  70项相关测试及5项CTest通过；Vela缺NonLocalPath真实能力，下一项Q2/Q3，旧负向诊断暂缓。

- [BVmethods 安装库原版对齐修订计划，2026-10-15](superpowers/plans/2026-10-15-bvmethods-installed-library-alignment-plan.md)：
  当前目标为 NonLocalPath／GradQuasiFermi、停止电流 1.3e-3；Q0审计与Q1原生参考已完成。
  NonLocalPath 前端当前报 E106，须补真实能力；旧 E2 证据保留原范围。

- [BVmethods 自主统计实际装配暂停检查点，2026-10-15](validation/bvmethods_statistics_assembly_checkpoint_2026-10-15.md)：
  用户要求暂停；基准 4,304 个统计点及接受点导数缓存逐位一致，T470p 私有装配与两端各 37 项测试完成。
  完整 RHS、第二版数学库对照分析和五态方向验证仍未完成，不认领纯 QS 资格。

- [BVmethods 自主初值、谷权重及统计值链，2026-10-15](validation/bvmethods_eta_seed_valley_2026-10-15.md)：
  独立 VM 候选解除原生初值／谷权重依赖，五态 70 个初值及最终 η、140 次更新逐位复现。
  UCRT 基准及实际统计导数／完整 RHS 仍待办；本机与 T470p 各 30 项测试通过。

- [BVmethods 积分值、导数、谷求和及 η 更新，2026-10-15](validation/bvmethods_integral_update_2026-10-15.md)：
  同源表＋独立 VM 数学库复现五态 385 次积分／导数、140 次谷求和与接受更新。
  原生初值条件下两弱行基准／方向闭合，自主完整统计与 RHS 仍待办；两端各 23 项测试通过。

- [BVmethods 权重生成、归一化与逐次舍入，2026-10-15](validation/bvmethods_weight_rounding_2026-10-15.md)：
  自主生成 240 个节点／权重及最后多项式操作数逐位一致；确认权重使用最后修正前的评估值。
  表生成补齐，严格 η 反算的弱行基准仍未闭合；本机 29 项、T470p 15 项测试通过。

- [BVmethods η 反算、积分表与非零方向，2026-10-15](validation/bvmethods_eta_inversion_2026-10-15.md)：
  五态确认多谷 240 点积分／Newton 路径，η 差异主要分解到有限阶积分和权重表示；迭代／求值项约 1.85e-18。
  自主 240 点电子弱行基准差约 −0.0464%，仍未闭合；本机 25 项、T470p 15 项测试通过。

- [BVmethods Γ 输入、缓存时序与 Bernoulli 分支，2026-10-15](validation/bvmethods_gamma_cache_bernoulli_2026-10-15.md)：
  14 个真实 Γ 端点绑定 48 次 SG 消费，未发现本首态旧缓存混用；三个 logl 末位差已定位。
  五态 480 个 Bernoulli 值／导数逐位复现；自主 η、完整 RHS 和原始纯 QS 仍未认证。

- [BVmethods SG 参数链、独立统计与非零方向，2026-10-15](validation/bvmethods_sg_argument_chain_2026-10-15.md)：
  五首态 120 个 a／反向控制逐位一致，完整首轮系统身份通过；新增参数审计与性质测试。
  独立统计方向较近但基准弱行未闭合，原生 Γ／Bernoulli 精度链续查；全 RHS 和纯 QS 未认证。

- [BVmethods 实际 SG 操作数与舍入分离，2026-10-15](validation/bvmethods_actual_sg_operands_2026-10-15.md)：
  T470p 实际装配核对人口／系数；原生操作数＋binary80 乘减控制下两弱行差约 7.15e-15／2.77e-15。
  100 位与输入量化控制分离求值路径；仍依赖原生操作数，自主 SG、全 RHS 和纯 QS 未认证。

- [BVmethods 端点物理身份与原始 RHS1 支撑绑定，2026-10-15](validation/bvmethods_endpoint_rhs1_binding_2026-10-15.md)：
  手册／实际类型确认电势、亲和能、带隙及收缩；两弱行全部 24 条元素边绑定 RHS1。
  独立参考量换算后原生弱行差约 2.69e-16／1.72e-15；Vela/native 全 RHS 和纯 QS 仍待办。

- [BVmethods 原始消费路径与逐边驱动，2026-10-15](validation/bvmethods_passive_consumer_2026-10-15.md)：
  两远端恢复；六批原生完整系统身份一致，取得实际 Bernoulli/驱动端点项并区分首态与 RHS2。
  当轮字段身份、缓存消费及弱行绑定待办已由最新记录补齐；继续 Q2/Q3，物理默认和门限未改。

- [BVmethods 输入参考、原生统计与近零驱动，2026-10-15](validation/bvmethods_reference_statistics_2026-10-15.md)：
  T470p 完成原生 logGamma/独立统计对照，新增十系统完整原生变量/RHS 身份一致，方向差约 4e-10。
  近零值仍未闭合；FermiDelta 字段不能直接等同人口导数，实际消费路径续查见最新记录。

- [BVmethods 观察器效应、拆分采集与非零驱动，2026-10-14](validation/bvmethods_observer_effect_drive_2026-10-14.md)：
  T470p 定位双 QF 组合读取效应；拆分后十个原生系统完整变量/RHS1 与原始基准一致。
  人口末位差解释为舍入路径；独立非零方向约 1e-9，近零基准仍未闭合，继续 Q2/Q3。

- [BVmethods 保精度装配回调与首态绑定，2026-10-14](validation/bvmethods_precision_callback_2026-10-14.md)：
  原生保精度人口/带边采集打通，两弱行打印 RHS 身份一致，完整 RHS 有四个末位变化。
  独立 η 重算仍未闭合输运；148 文件哈希通过，8 项性质测试通过，继续 Q2/Q3。

- [BVmethods 独立 η、带边及非零驱动，2026-10-14](validation/bvmethods_eta_band_precision_2026-10-14.md)：
  同首态补导出 DOS/带边，T470p 完成 36 个独立定义积分反算；80-bit 非零响应差约 7.66e-7/3.21e-7。
  近零值仍未闭合；扩展精度只读 PMI 采集尝试未通过，失败证据保留，完整 RHS/四臂资格仍待办。

- [BVmethods 内部驱动、人口统计与精度控制，2026-10-14](validation/bvmethods_drive_population_precision_2026-10-14.md)：
  原生两种精度共 14 个首轮系统、T470p 五次零更新探针完成；107/26 文件下载核验。
  相同导出 QF 不足以绑定人口/RHS，微小请求出现派生量化；继续审计有效统计驱动，完整 RHS 尚未对齐。

- [BVmethods RHS 时序、Load 与接触审计，2026-10-14](validation/bvmethods_rhs_load_contact_audit_2026-10-14.md)：
  T470p 上轮 293 文件下载封存通过；8 对连续 RHS 身份和 17 个实际线性系统复核通过。
  五次实际输运探针区分身份与电势映射控制，弱行余差仍在；完整 RHS 未对齐，四臂导数资格受前置条件约束。

- [BVmethods 私有动态场、完整 RHS 与弱行精度控制，2026-10-13](validation/bvmethods_live_rhs_chain_2026-10-13.md)：
  T470p 完成 47 次实际探针，动态 E2 源及解析链式项接入；同方向细步差约 6.81e-9。
  固定输入的 100 位 SG 求值不能解释两弱行余差，完整 RHS 仍未闭合。
  10 月 13 日离线时仅有摘要；10 月 14 日原始输出下载及封存已补齐，见上项。

- [BVmethods 自主节点场几何、接触方向与近零 QF 控制，2026-10-13](validation/bvmethods_field_geometry_state_2026-10-13.md)：
  T470p 自主生成盒转移及节点场规则，六区域基准场约 1e-14、接触/氧化层方向约 5.58e-11；Si 显式 E2 源约 4.95e-14。
  17 组源关闭控制分离近零 QF 导出路径差，两弱行尚未完全闭合；128-bit 初始化积分超时未认证。
  本机和 T470p 各 28 项测试通过；生产默认未改，完整 RHS／Jacobian 与从零纯 QS 仍待办。

- [BVmethods 节点场算子、显式 E2 实际残差与输运弱行，2026-10-13](validation/bvmethods_node_field_transport_2026-10-13.md)：
  T470p 完成原生一环响应的独立扰动验证和 2088 行自主场恢复；显式参数＋原生节点场的实际 BTBT 分项闭合。
  两条弱行的非零边系数由原生 Jacobian 复核至约 1e-13／1e-15；该轮尚未识别几何公式，后续自主规则见最新记录。
  本机和 T470p 各 20 项测试通过；完整 RHS／Jacobian 和从零纯 QS 仍未认证，继续 Q2/Q3。

- [BVmethods 自主局部迁移率、弱行及 BTBT 实际场，2026-10-13](validation/bvmethods_mobility_weak_btbt_2026-10-13.md)：
  T470p 自主局部 µ 及两次私有实际装配完成；HFS 接触回退、切向投影及截断闭合，输运弱行余差仍在。
  零源观察器及参数导出确认 E2 的节点场和 A2/B2 差异，原生场＋原生参数重建约 3.55e-15。
  本机 91 项、T470p 15 项辅助测试通过；自主场恢复及全系统导数未闭合，继续 Q2/Q3。

- [BVmethods 逐边输运、Poisson 与 SRH/BTBT 支撑，2026-10-13](validation/bvmethods_edge_poisson_source_2026-10-13.md)：
  T470p 完成 19 次固定态探针；原生几何＋盒平均局部 µ 的私有实际输运差降至电子约 9.76e-5、空穴约 4.82e-4。
  SRH 体积差、Poisson 边／电荷盒及常数指纹已分离；BTBT 场支撑和弱行余差仍待定位，未认领完整 RHS 或 QS 资格。
  本机 34 项辅助回归、T470p 6 项性质测试通过；生产默认未变，证据独立封存。

- [BVmethods 私有实际装配源与完整 RHS，2026-10-12](validation/bvmethods_private_assembly_source_2026-10-12.md)：
  分支／局部取样接入实际残差，源积分内部闭合至约 4e-16；同精度原生首态及源 RHS 支撑复核通过。
  空穴完整残差差降至约 0.57%，电子输运和 Poisson 仍未对齐；本机 57 项、两平台各 18 次探针和 10 项测试通过。
  生产默认未变，继续 Q2/Q3；导数、Newton 更新与纯 QS 仍待办。

- [BVmethods 单元素分支与唯一 SG／内置源，2026-10-12](validation/bvmethods_element_branch_sg_2026-10-12.md)：
  单元素控制识别 E0 ±5% 内的升余弦 alpha 过渡，原物理及七个完整网格阈值案例验证通过。
  已知代表身份下独立重建全部 SG／源；编码态约 2e-14、近零源约 1e-8，弱电流余差保留。
  本机 47 项、两平台各 6 项测试及分支／四组源复算通过；生产默认、完整 RHS 和纯 QS 资格未变。

- [BVmethods 源迁移率取样与空穴分支，2026-10-12](validation/bvmethods_source_mobility_transition_2026-10-12.md)：
  回调顺序唯一标记单元和代表顶点，采用该次局部迁移率后，两固定态的完整硬 PMI 源闭合至约 3e-16。
  空穴内置分支仍有约 0.0733% 差异；节点源反算秩亏，尚未识别原生过渡公式。
  本机 41 项、两平台各 6 项测试及各四组源复算通过；继续 Q2/Q3，完整 RHS 和纯 QS 未认证。

- [BVmethods 端点 QF、阈值分支与元素源支撑，2026-10-12](validation/bvmethods_qf_branch_element_2026-10-12.md)：
  三机验证实际 reference/increment；零／节点参考值没有消除主要 RHS 差异。
  七组原生阈值消融确认空穴分支响应；最大内角规则覆盖已有唯一主要回调，但两种迁移率位置的完整源重建均未通过。
  本机 35 项、两平台各 6 项测试通过；仍在 Q2/Q3，唯一源算子和从零高压纯 QS 未认证。

- [BVmethods 有效带边 SG、QF 精度与源积分，2026-10-12](validation/bvmethods_sg_precision_source_transfer_2026-10-12.md)：
  三机实际装配对照显示 SG 单项替换未改善完整 RHS；绝对 QF 往返暴露弱电流舍入敏感性。
  扩展精度改善近零回调，盒源到内置 RHS 闭合；0.0733% 空穴观察器差定位到低／高场分支。
  本机 29 项及两平台各 8 项测试通过，仍在 Q2/Q3，唯一元素源支撑和完整纯 QS 未认证。

- [BVmethods 同本构态实时雪崩与 SG 支撑核验，2026-10-10](validation/bvmethods_consistent_state_sg_source_2026-10-10.md)：
  含接触的人口/QF 在私有诊断态闭合，实时重算排除旧冻结雪崩大差；完整 RHS 和弱空穴行仍有差异。
  编码激励定位接触端点绝缘投影例外，主要回调匹配显著改善；全局唯一 SG 支撑未认证。
  本机 21 项、两平台新增 4 项测试及分析通过，生产物理和从零高压纯 QS 资格未变。

- [BVmethods 同人口 QF、接触与完整残差分解，2026-10-10](validation/bvmethods_qf_contact_transport_split_2026-10-10.md)：
  三个更新前状态的控制区分接触人口和 QF 映射造成的大输运差；完整 RHS、源项和弱行全部保留。
  原生 Poisson 盒、单元边系数及独立矩阵尺度可解释残差差；生产默认未变。
  冻结源控制仍有雪崩和输运余项，原生 QF 实际算子改变人口；完整物理/唯一 SG 支撑和高压纯 QS 未认证。

- [BVmethods 迁移率顺序修复、导数与完整 RHS 验证，2026-10-10](validation/bvmethods_mobility_order_fix_2026-10-10.md)：
  非 PhuMob 组合现先计算表面低场迁移率，再施加高场饱和；模型链式导数和实时残差差分通过。
  本机 Release 全量 1049 项最终通过，T470p/C560 各 95 项迁移率及 143 项 Newton 回归通过。
  两平台新鲜零偏冒烟通过；零偏行/源资格未激活，不替代非零轨迹验收。
  SG 边通量、分项和完整 RHS 在 Vela 内部闭合；原生共同自由行仍未整体对齐，继续 Q2/Q3。

- [BVmethods SG 支撑、边界投影与迁移率顺序，2026-10-10](validation/bvmethods_sg_support_followup_2026-10-10.md)：
  13 组原生固定状态观察完成；两边 SG 配合绝缘边界投影可复现有限 QF 激励回调，代表支撑仍有歧义。
  直接 HFS 回调确认原生先 Enormal 后高场，Vela Masetti 顺序相反；同输入局部空穴顺序贡献约 2.88%。
  未改生产物理或认领完整 RHS/Jacobian、从零高压纯 QS 通过，继续 Q2/Q3。

- [BVmethods 本构契约、同状态几何与源项复核，2026-10-10](validation/bvmethods_constitutive_geometry_rebaseline_2026-10-10.md)：
  已修复前端 Fermi–OldSlotboom 联动及显式 NoFermi；两平台新鲜零偏和冻结雪崩导数短程重新基线化通过。
  原生更新前状态、DOS/带边/QF/局部迁移率及源电流已捕获，Poisson 介电系数和半导体盒体积获得矩阵验证。
  DOS 质量公式、统计近似和 SG/源支撑仍有差异；完整 RHS/交叉块与原始高压纯 QS 尚未认证。

- [BVmethods 矩阵元素、RHS 与物理映射复核，2026-10-10](validation/bvmethods_matrix_physics_alignment_2026-10-10.md)：
  参考 SLOT 聊天，补齐自由人口匹配、完整 Fermi 跨块变换及 6764 个共同自由器件行的逐元素比较。
  确认前端未启用 Fermi–OldSlotboom 统计修正；只读控制显著缩小本构差，但 RHS、Poisson 几何及源导数仍未对齐。
  BVmethods 的 density_gradient 实为 SG 边通量，不等同 SLOT 旧梯度源核；未改正式模型或认领完整 QS 修复。

- [BVmethods 密度更新与源项精度调查，2026-10-10](validation/bvmethods_density_update_followup_2026-10-10.md)：
  T470p/C560 匹配状态复验；核对原生条件式 1% 行为，密度坐标候选仍失败。
  排除微小 SRH 分块相减的误读，原始空穴方向对雪崩 Jacobian 耦合显著敏感。
  冻结雪崩导数候选短程至外电压 0.01 V 通过：T470p 单次从零自适应扫描，
  C560 超时后按匹配检查点续扫；两平台保存态独立门限通过。T470p 直接首目标仍被线搜索拒绝。
  继续 Q2/Q3；完整从零纯 QS 尚未验收。

- [BVmethods 原生 Jacobian / RHS 导出与对照，2026-10-09](validation/bvmethods_sentaurus_debug_matrix_2026-10-09.md)：
  `sdevice -d` 成功导出 14 次系统，实际 6843 阶耦合矩阵回代闭合；
  完成节点、单位及弱行形状对照，发现原生非正密度更新的局部处理线索。
  同状态全 Fermi Jacobian 认证和原始完整 QS 尚未完成，继续 Q2/Q3。

- [BVmethods 耦合接受检查与弱行更新，2026-10-09](validation/bvmethods_coupled_acceptance_weak_update_2026-10-09.md)：
  真实行/连续性/KCL/负载检查及 CSV JSON 转义已补齐，接触行正则化已修复。
  完整场导数的稳态恢复候选通过首点和匹配下一点；完整高压轨迹及空穴方向导数仍待验证。
  三平台 Release 回归及两平台原始证据见记录；继续 Q2/Q3，Q4–Q7 未启动。

- [BVmethods 场迁移率导数、Newton 更新限制与外电阻复验，2026-10-09](validation/bvmethods_mobility_loadline_audit_2026-10-09.md)：
  前序快照：两平台从新鲜零偏纯 QS 取得外电压 10 µV 的严格通过点；续扫及匹配导数控制失败。
  QR 保存态 KCL 超限、耦合默认诊断缺口已在上项补修；完整 Q4 尚未启动。

- [BVmethods 物理/数学字段隔离对照，2026-10-09](validation/bvmethods_model_math_ablation_2026-10-09.md)：
  11 组原生单因素诊断完成；Vela 保留物理 residual、关闭场迁移率 Jacobian 导数，
  三平台直接内电压 10 µV 首点通过；本轮历史快照未接回外电阻，后续复验见上项，原模型高压 MinStep 保留。

- [BVmethods 自由弱行、线性精化与同方向导数审计，2026-10-08](validation/bvmethods_weak_rows_direction_audit_2026-10-08.md)：
  修复精化缩放及 Eparallel 源项 Jacobian；三平台完整 Release 回归通过，候选从零首个 10 µV 点
  仍被载流子逐行门限拒绝，全局连续性/拒绝态 KCL 通过；Q2/Q3 继续，Q4 未启动。
- [BVmethods 分支、载流子条件化及 QS 步长/失败追踪，2026-10-08](validation/bvmethods_branch_conditioning_step_trace_2026-10-08.md)：
  修复外目标自适应首步和失败状态追踪；两平台相关回归通过，六次回退仍未获得非零接受点。
  相邻边导数分支及少数载流子弱行进一步定位；数值修复仍在 Q2/Q3，Q4 未启动。
- [BVmethods 外电阻与有偏置导数对照，2026-10-08](validation/bvmethods_loadline_biased_derivatives_2026-10-08.md)：
  两平台完成标量/耦合及缩括区、QR 对照；零偏通过，未取得非零接受点。
  有偏置电子方向差异定位到雪崩源项；仍在 Q2/Q3，Q4 未启动。
- [BVmethods 零偏初始化与导数诊断，2026-10-08](validation/bvmethods_zero_bias_initialization_2026-10-08.md)：
  修复合并初始阶段丢失完整 Poisson 预算，保留暖启动 Newton 交接；两平台各两次新鲜零偏通过。
  零电流源项导数边界仍保留；后续控制器对照见上项，完整 QS 轨迹尚未验收。
- [BVmethods 从零纯 QS Q0–Q3 执行记录，2026-10-08](validation/bvmethods_zero_qs_q0_q3_2026-10-08.md)：
  首轮快照：两节点 Release 部署及相关回归通过；修复 Eparallel 映射；原生仍止于 14.1 kV，
  首轮 Vela 初始化候选均失败，既有两支路空间差异补导出；最新零偏结果见上项。
- [BVmethods 从零纯 QS 后续计划，2026-10-08](superpowers/plans/2026-10-08-bvmethods-zero-qs-followup-plan.md)：
  Q0–Q3 首轮、零偏修复、控制器/分支/弱行审计及最小实现修复已执行；仍在 Q2/Q3，Q4–Q7 未启动。
  T470p 负责 Windows，联想 C560 负责 Ubuntu 26.04，VM 负责原生对照。

- [SDevice package acceptance and commit closeout, 2026-10-07](validation/sdevice_package_closeout_2026-10-07.md):
  current passing Windows/Ubuntu 26.04 packages, full regression, ten strict BV
  anchors per platform, relative RUNPATH, concurrency and source/evidence commits.

- [BVmethods 从零纯 QS 诊断，2026-10-06](validation/bvmethods_zero_qs_diagnosis_2026-10-06.md)：
  原路径 13 组及辅助栅 2 组原生对照、13 kV 下的两个不同稳态解及独立整条轨迹校验；
  原固定栅压路径和辅助栅候选均未通过，149/36 文件归档逐项核验，当前 v4 资格保持原范围。

- [SDevice 遗留修复与 WP14 原生补验，2026-10-06](validation/sdevice_followup_2026-10-06.md)：
  Schottky 完整重启、PN2D 十点严格场量通过；原生五组栅控制确认 Eg(0)/DOS 映射；
  WP14 BVmethods v4 恢复版原生/稳态及 C7 六项通过，原始从零 QS 限制保留。

- [WP13 built-in Si/SiO2 alignment, 2026-10-05](validation/builtin_material_alignment_2026-10-05.md):
  no-BGN 300 K defaults, dependency audit and explained regression baselines;
  Windows 1026-case regression and Ubuntu 26.04 16-case affected-scope verification.

- [SDevice workflow](sentaurus_sdevice_workflow.md): current supported subset,
  CLI/gates/messages, output units and executable neutral mesh/sidecar examples.
- [SDevice log/documentation repairs and native package E2E, 2026-10-05](validation/sdevice_native_package_e2e_2026-10-05.md):
  historical package acceptance and integrity checks under the then-current contract;
  its -20 V hole-density limitation was resolved and revalidated in the 2026-10-07 closeout above.
- [SDevice WP12 documentation closeout, 2026-10-05](validation/sdevice_wp12_closeout_2026-10-05.md):
  command/link/message/field audits, Release regression and six-item DoD evidence;
  historical snapshot preceding the package E2E follow-up above.
- [Changelog](../CHANGELOG.md): SDevice frontend delivery and qualified limits.

- [SDevice WP11 Ubuntu 26.04 migration, 2026-10-04](validation/sdevice_wp11_ubuntu26_04_2026-10-04.md):
  Windows and Ubuntu 26.04 relocatable packages, retired previous Linux artifacts, offline Python/license provenance,
  both-platform frontend regression and exact Linux concurrency/thread-cap verification.
- [SDevice packaging guide](sdevice_packaging.md): offline inputs, runtime layout,
  manifests, licenses and extracted-package/concurrency checks.
- [Linux relative RUNPATH validation, 2026-10-06](validation/sdevice_linux_runpath_2026-10-06.md):
  all 66 dynamic ELF objects use relative RUNPATH; relocated package dependency,
  bundled-Python concurrency and fingerprint checks pass.

- [SDevice WP10 closeout, 2026-10-03](validation/sdevice_wp10_closeout_2026-10-03.md):
  independent, pinned BV spatial qualification at −13.2 V under the original contract;
  named breakdown-knee comparison and near-knee diagnostics, with the −20 V hole-density limitation retained.

- [SDevice IV/BV repair, 2026-10-03](validation/sdevice_iv_bv_repair_2026-10-03.md):
  corrected model activation, Ohmic contacts and Newton continuation; IV reaches 10 V and passes,
  BV reaches −20 V with passing curve/KCL gates; repair-stage snapshot preceding the scoped WP10 closeout.

- [SDevice WP10 acceptance, 2026-10-03](validation/sdevice_acceptance_2026-10-03.md):
  verified native v2 archives, frozen current sign and per-stage Time, fail-closed
  curve/field gates and zipapp exit-code repair; initial failures retained as historical evidence.

- [SDevice WP5 orchestration and WP6.1 metadata, 2026-10-02](validation/sdevice_orchestrator_validation_2026-10-02.md):
  native CLI and standalone zipapp, gates, locking, interruption and incremental outputs;
  Release 1020/1020 CTest and real-runner Windows Unicode-path smoke checks.

- [SDevice WP4 stage/deck emitter, 2026-10-02](validation/sdevice_deck_emit_validation_2026-10-02.md):
  local main merged; seven fixture goldens, v2 decks, HDF5 chaining and real-runner smoke checks.
  WP4 acceptance passed with the confirmed current WP3 gate reference baseline; CLI orchestration is now delivered by WP5.

- [Remote compute VM](office_vm_workflow.md): private SSH configuration,
  Micromamba GCC 16 builds, isolated Sentaurus environment, file transfer and
  detached job submission/status through `scripts/Invoke-VelaVm.ps1`.

- [SimpleMOS raw states/logs archive and local main migration](../reference_tcad/simplemos_sentaurus2022/raw_assets_20261001/README.md):
  complete local H5/TDR coverage, 60,954 assets including structured logs and existing bundles,
  persistent VM storage, restore tool and source-to-main commit mapping.

- [SimpleMOS historical evidence archive and recovery](../reference_tcad/simplemos_sentaurus2022/local_evidence_20260930/README.md):
  1,424 small records in Git; complete 1,946-file archive on the Sentaurus VM,
  with exact hashes and recovery of 522 detailed historical artifacts to their original paths.

- [Latest D0 Release closeout, 2026-09-25](validation/templates_ldmos_release_closeout_2026-09-25.md):
  Codespaces Release 858/858 CTest, 590 exact points, 698 identical state pairs,
  62 native local fields and six checkpoints passed. Three paired rounds reduce
  mean dual-gate wall time by 1.10%; all 6,390 archive members verified remotely,
  with the small verification report retrieved locally at the user's request.

- [Linked defaults and electrothermal cross-point linear reuse, 2026-09-25](validation/templates_ldmos_thermal_linear_reuse_2026-09-25.md):
  unified production defaults, exact-pattern reuse and failure/restart invalidation;
  UCRT64 Debug 858-test regression passed; subsequent full-curve timing is recorded
  in the latest D0 Release closeout above.

- [Shared Halley default integration, 2026-09-25](validation/templates_ldmos_halley_default_2026-09-25.md):
  C++ and isothermal/electrothermal JSON defaults unified, explicit legacy retained.
  Codespaces Release 858/858 CTest and seven full D0/G3/D5/D4 curves (217 points)
  passed. Default D0 states match the qualified explicit Halley run exactly.

- [Screening-root full D0 qualification and repeats, 2026-09-25](validation/templates_ldmos_screening_full_2026-09-25.md):
  20 full curves / 620 exact points, 10 native joint checks and 744 state pairs
  passed. Halley reduced aggregate wall time by 6.70% over three fresh
  interleaved dual-gate rounds (5.68%-7.24% per round); all six gate pairs were
  faster. Newton work is unchanged, with small assembly-count differences.
  Evidence was retrieved; default was legacy in that experiment. The subsequent
  shared-default integration above promotes Halley after full regression.

- [Screening-root candidates on Codespaces, 2026-09-24](validation/templates_ldmos_screening_candidates_2026-09-24.md):
  local main merged; Release 856/856 CTest and eight short D0 curves passed.
  Safeguarded Newton/Halley/TOMS748 retain states and Newton work; one-round
  paired wall reductions are 6.36%/7.18%/6.32%. Legacy was the default then;
  full curves and repeated end-to-end qualification are tracked in the newer
  2026-09-25 report above.

- [Codespaces build environment](codespaces.md): container dependencies and
  PowerShell/Codex commands for remote builds, tests, and result retrieval.

- [T470p electrothermal subcost measurements, 2026-09-24](validation/templates_ldmos_thermal_subcost_2026-09-24.md):
  eight full D0 curves and four joint checks passed with identical states/work;
  screening roots dominate the three measured candidates. Thermal scatter
  requires explicit clock-bias qualification; diagnostic timers default off.

- [T470p electrothermal gprof, 2026-09-24](validation/templates_ldmos_thermal_gprof_2026-09-24.md):
  latest default structure reuse, full D0 Release/gprof controls and original gates;
  assembly/IALMob and factorization dominate. Four-equation linear-object reuse
  across points remains distinct from the completed assembly cache.

- [Default electrothermal assembly structure reuse, 2026-09-24](validation/templates_ldmos_thermal_structure_2026-09-24.md):
  four-equation/thermal structure reuse; C560/T470p full regression, D0 curves,
  three paired repetitions and omitted-option checks passed. Default is on;
  measured aggregate timing improvement is small (0.70%).

- [Default DC assembly structure reuse, 2026-09-23](validation/templates_ldmos_structure_default_2026-09-23.md):
  production option, invalidation boundaries and T470p default-path regression.
- [Restart state storage](state_archive.md): current `vela.state/2` fields, units,
  mesh/source identity, atomic checkpoint writes and strict recovery semantics.
- [HDF5 production state migration plan, 2026-09-23](validation/templates_ldmos_hdf5_production_migration_plan_2026-09-23.md):
  authorized scope, state schema, restart transactions, active-fixture migration,
  electrothermal coverage and acceptance gates. See the execution record for final results.
- [HDF5 migration execution, 2026-09-23](validation/templates_ldmos_hdf5_migration_execution_2026-09-23.md):
  T470p Release validation; 850 CTest entries, G3/D4/D0 and frozen-300 K recovery
  passed. Final three-round D5 overhead is 0.4572%, below 3%; C560 Linux CTest passed.
- [HDF5 Linux CI, 2026-09-23](validation/templates_ldmos_linux_ci_2026-09-23.md):
  isolated Ubuntu 26.04/GCC 15 Debug build and 850/850 CTest passed on C560 with
  SparseLU; actual dependency detection and GitHub Actions differences are recorded.

- [LDMOS worktree cleanup, 2026-09-22](validation/templates_ldmos_worktree_cleanup_2026-09-22.md):
  duplicate artifacts and explicitly authorized historical JSON/CSV outputs removed;
  118.68 GB total logical file bytes released. Production dependencies and latest HDF5
  evidence retained; old per-point output paths may no longer be replayable.

- [TDR and Vela HDF5 schema compatibility, 2026-09-22](validation/templates_ldmos_tdr_hdf5_compatibility_2026-09-22.md):
  actual files share the HDF5 container but require different schemas; importer tests
  and bidirectional schema-rejection checks pass. No native TDR write compatibility claimed.

- [T470p HDF5 full-curve repeated pairing, 2026-09-22](validation/templates_ldmos_hdf5_full_repeats_2026-09-22.md):
  three alternating serial rounds complete: twelve full D5 curves / 372 points, six joint
  gates and all-field repeated audits pass. HDF5 adds 0.8695% total wall time versus VDS1,
  meeting the approved 3% criterion; unchanged Release program, numerical gates and defaults.

- [Public DD state formats on T470p, 2026-09-22](validation/templates_ldmos_public_state_formats_2026-09-22.md):
  optional HDF5/HighFive, FlatBuffers and structured NPY interfaces pass local/remote
  checks and twelve paired short curves / 96 exact points; all-field equality confirmed.
  One-round combined wall costs are +0.297%/+0.304%/+0.672% versus matched VDS1.
  No speed ranking or default promotion; HDF5 remains the preferred public-storage candidate.

- [T470p binary restart state and reduced seed writes, 2026-09-22](validation/templates_ldmos_binary_state_2026-09-22.md):
  eight short curves / 64 exact points pass unchanged gates and identical states
  and trajectories. VDS1 saves 12.11%/10.95% versus CSV in one paired round.
  Separate memory predictor seeds remove 60 files per curve but add 11.78%/12.52%
  wall time including transport and reconstruction audit; not promoted. Defaults unchanged.

- [T470p D5 file reuse and cost-per-volt controls, 2026-09-22](validation/templates_ldmos_d5_file_efficiency_2026-09-22.md):
  eight short curves / 64 exact points pass unchanged gates and cold final audits.
  Bounded CSV/hash reuse saves 3.00%/2.61% in one paired round with identical trajectories.
  Separate efficiency control saves six services per curve but increases combined wall time
  0.87%; cost-regression feedback never triggers. No default promotion or full-curve claim.

- [T470p D5 cost-feedback continuation and parent timing, 2026-09-22](validation/templates_ldmos_d5_cost_continuation_2026-09-22.md):
  dual-gate first-eight-point controls pass original state gates, but wall time grows
  10.57%/10.21% and updates 367/349 to 386/379. Candidate not promoted to full curves.
  Nonoverlapping parent timers locate 26.1/26.3 s outside worker calls in the controls;
  budget-per-step feedback does not minimize total continuation cost. Defaults unchanged.

- [T470p D5 assembly/Fermi reuse candidates, 2026-09-22](validation/templates_ldmos_d5_kernel_candidates_2026-09-22.md):
  immutable pattern/scatter reuse passes dual-gate first-eight-point controls;
  one paired round reduces wall time by 5.66%/5.90% with identical states and trajectories.
  Fermi endpoint reuse reduces evaluations without a clear timing gain. Combined
  62-point full curves pass original gates and exact state/trajectory checks;
  full-curve qualification is not paired timing. Defaults unchanged.

- [T470p D5 completed gprof results, 2026-09-22](validation/templates_ldmos_d5_gprof_results_2026-09-22.md):
  62 exact points pass unchanged curve/joint gates, with identical states, trajectories
  and recorded counters to ordinary Release. Cross-checks Fermi edge assembly, repeated
  pattern/scatter construction and point-service overhead; proposes progress/cost-based
  continuation without v3's fixed voltage split. No production-policy changes.

- [T470p D5 baseline gprof and adaptive continuation, 2026-09-21](validation/templates_ldmos_d5_gprof_t470p_2026-09-21.md):
  historical pause record and isolated Release/frozen DLL profiling contract;
  resumed and completed on September 22 (see results above).
  Separates v3's case-specific voltage split from a proposed
  progress/cost-driven adaptive controller; no production changes.

- [T470p D5 protected 0.8 V full curves, 2026-09-21](validation/templates_ldmos_d5_step08_full_2026-09-21.md):
  four full curves and fresh paired states pass unchanged gates. One paired round
  reduces wall time by 9.70%/29.06%, but Vg4 Newton updates increase 5.34% due to
  costly recovery points. Low-voltage and no-history safeguards are required;
  two earlier failures remain separate evidence. No production-default promotion.

- [T470p D5 targeted step/prediction/row controls, 2026-09-21](validation/templates_ldmos_d5_targeted_controls_2026-09-21.md):
  20 bounded paths and 12 node-trace controls qualify; 0.8 V steps with fixed 0.2 V QF
  cap reduce short-window work. Larger jumps can trigger internal recovery. Diagnoses
  one low-voltage prediction counterexample and corrects legacy trace-density units
  in analysis; no full-curve or production-default promotion.

- [T470p D5 Newton and prediction localization, 2026-09-21](validation/templates_ldmos_d5_newton_localization_2026-09-21.md):
  six repeated UMFPACK traces: 251 drain advances versus native 49; average updates per
  advance comparable. Sixteen fixed-state prediction controls pass original gates;
  one low-voltage prediction counterexample identified. No production default changes.

- [LDMOS performance closeout and production defaults, 2026-09-20](validation/templates_ldmos_performance_closeout_2026-09-20.md):
  cross-point main-Newton analysis reuse enabled by default; UMFPACK preferred,
  SparseLU fallback, then STRUMPACK/MUMPS/SuperLU_MT in the available-backend list.
  Portable D5 inputs and dual-gate validation entry; historical experiments remain explicit controls.

- [T470p 1/2/4-thread pilot, 2026-09-20](validation/templates_ldmos_t470p_threads_2026-09-20.md):
  39 matrix replays and 13 Vg8 first-eight-point D5 curves pass; 104 exact points,
  unchanged gates. STRUMPACK 2/1 takes 167.445 s versus 230.187 s at 1/1;
  UMFPACK/MUMPS do not benefit under tested settings. One curve round, defaults unchanged.

- [T470p identical-binary host pilot, 2026-09-20](validation/templates_ldmos_t470p_host_pilot_2026-09-20.md):
  five single-thread reuse backends, Vg8 first eight D5 points on N150 and T470p;
  all 80 points pass, remote wall times 12.79–34.48% lower in this single round.
  Full-curve repeats remain paused; host timings must not be pooled across machines.

- [Single-thread five-backend repeats, 2026-09-20](validation/templates_ldmos_single_thread_repeats_2026-09-20.md):
  STRUMPACK dual-gate qualification first, then all five reuse backends with solver/BLAS
  1+1, three fresh D5 rounds (30 curves); original gates and historical timings preserved.

- [Full-curve backend timing execution, 2026-09-20](validation/templates_ldmos_full_curve_timing_execution_2026-09-20.md):
  all 12 D5 first-round curves and six dual-gate qualifications pass; stopped for shortlist
  review. Future runs use solver/BLAS 1+1 for every backend; historical STRUMPACK 2+1
  results stay separate. Repeats and D4 pending; no production default change.

- [Full-curve backend timing plan, 2026-09-20](validation/templates_ldmos_full_curve_timing_plan_2026-09-20.md):
  revised future thread contract fixes all five solvers and BLAS to one thread.
  Original schedule retained as provenance; thread changes require a new frozen batch,
  not a resume that mixes old STRUMPACK two-thread timings with single-thread repeats.

- [Four-backend cross-point curve validation, 2026-09-19](validation/templates_ldmos_four_backend_reuse_2026-09-19.md):
  SparseLU/UMFPACK/MUMPS/SuperLU_MT, ordinary worker versus analysis reuse,
  all 16 dual-gate D5 first-eight-point curves pass with identical paired states,
  Newton updates and factorization counts. One paired round, solver/BLAS 1+1;
  default remains off, auxiliary analyses remain separate.

- [Cross-point Newton analysis reuse, 2026-09-19](validation/templates_ldmos_cross_point_analysis_2026-09-19.md):
  default-off sequential linear context; all six D5 first-eight-point controls pass.
  Incremental wall savings over ordinary worker are 8.57%/8.78%; main-only services
  add no analyses after reuse. Auxiliary solvers remain separate; one round only.

- [STRUMPACK ordering comparison, 2026-09-19](validation/templates_ldmos_strumpack_ordering_2026-09-19.md):
  METIS/AMD/MMD/AND over 36 captured systems, three analysis policies and three rounds;
  all 1,296 systems pass. METIS wins with reuse, AMD wins with per-matrix rebuilding.
  Ordering/tree timings and fill recorded; factor structure varies across some repeats.
  No production or curve qualification change.

- [STRUMPACK verified 2+1 threads: D5 first-eight-point controls, 2026-09-19](validation/templates_ldmos_strumpack_verified_threads_2026-09-19.md):
  all four dual-gate curves and runtime thread audits pass; end-to-end wall time is
  2.87%/4.60% higher than paired UMFPACK. Faster factorization is offset by analysis
  and solve costs; one paired round only, no full-curve qualification or promotion.

- [BLAS threading and STRUMPACK compression, 2026-09-19](validation/templates_ldmos_blas_compression_2026-09-19.md):
  isolated Eigen BLAS and compression candidates; runtime thread verification and
  13 configurations over 36 matrices, three rounds, all accuracy gates pass.
  Corrects the previous BLAS=1 assumption for OpenMP OpenBLAS; no production promotion.

- [Sparse backend and METIS execution, 2026-09-18](validation/templates_ldmos_sparse_backend_execution_2026-09-18.md):
  optional adapters and linear-contract tests implemented; historical/current matrix screening passes.
  All six single-thread D5 dual-gate first-eight-point controls pass; no new wall-time winner.
  STRUMPACK 2/4-thread curve controls also pass but are slower than matched UMFPACK controls;
  SuperLU_MT multithreaded replay failures retained. No winner triggered the full-curve repeat stage.
  Screening complete; default unchanged.

- [Sparse backend and METIS comparison plan, 2026-09-18](validation/templates_ldmos_sparse_backend_plan_2026-09-18.md):
  UCRT64 MUMPS/SuperLU_MT/STRUMPACK/METIS inventory and five smoke checks pass;
  staged integration and fixed-matrix/curve plan only, no new production backend qualification.

- [D5 diagnostics-disabled repeats and strategy controls, 2026-09-18](validation/templates_ldmos_umfpack_low_noise_2026-09-18.md):
  all 8 curves and original joint gates pass; four wall-time gains of 10.38%–17.93%.
  Repeated states/counters match; strict idle conditions were unavailable.
  Historical frozen-matrix controls distinguish ordering, scaling and pivot effects.

- [Isothermal full-curve backend repeats, 2026-09-18](validation/templates_ldmos_umfpack_full_repeats_2026-09-18.md):
  all 16 D5/D4 curves and original joint gates pass; repeat states are identical.
  Factor fill/ordering/work/memory quantified; D4 wall-time gains repeat, D5 timings vary.

- [Isothermal SparseLU/UMFPACK controls, 2026-09-18](validation/templates_ldmos_isothermal_umfpack_2026-09-18.md):
  optional shared backend with symbolic/numeric reuse and failure checks;
  fixed-system replay and first-eight-point controls, qualification tracked separately from R11.

- LDMOS group review, September 12–17, 2026: the generated report at
  `reports/ldmos_review_2026-09-12_17/report.md` is not present in this checkout.
  See the [versioned September 16–17 work summary](validation/templates_ldmos_work_summary_2026-09-16_17.md)
  for the corresponding scope and validation history.

- [Architecture](architecture.md): source tree map, solver paths, and supported
  implementation boundaries.
- [Config schema](config_schema.md): implementation-aligned JSON field
  reference for Poisson, DC sweeps, Newton, unit-scaling input mode, contacts,
  boundaries, solver options, and regression blocks.
- [Sentaurus import](sentaurus_import.md): HDF5/TDR import prerequisites,
  `sentaurus_import` CLI usage, and end-to-end conversion workflow.
- [Sentaurus VMware SSH workflow](sentaurus_vm_ssh_workflow.md): Host-only
  VMware networking, legacy CentOS SSH setup, and remote Sentaurus run/copy
  commands for the local Sentaurus 2018 VM.
- [SDevice cmd direct-run workflow](sentaurus_sdevice_workflow.md): frozen
  `vela device` contracts (execution IR v2, mesh qualification, stage plan,
  artifact manifest, sidecar), CLI/exit codes, output naming and the message
  catalogue; fixtures under `tests/fixtures/sdevice/`. Reference-artifact
  inventory and pending archive decisions:
  [2026-09-26 inventory](validation/sdevice_reference_artifact_inventory_2026-09-26.md).
  Implementation status, remaining work packages and handoff notes:
  [2026-09-30 handoff](validation/sdevice_cmd_frontend_handoff_2026-09-30.md).
- [Poisson unit-scaling notes](development_poisson_unit_scaling.md): developer
  notes for the scaled Poisson assembly path used by
  `scaling.mode = "unit_scaling"`.
- [PN2D BV validation](validation/pn2d_bv_validation.md): current qualified
  template policy, validation gates, limitations, and evidence map.
- [SimpleMOS engineering reference](../reference_tcad/simplemos_sentaurus2022/engineering/README.md):
  qualified configuration, run entry, external input requirements, 816-point
  evidence and a portable read-only hash/coverage audit.
- [SimpleMOS branch status](validation/simplemos_branch_status.md): explicit
  production options, diagnostic boundaries, local evidence dependencies, and
  the current n23 Id-Vg comparison.
- [Validation evidence](validation/): dated reports and machine-readable
  contracts for checked-in reference TCAD work.
- [LDMOS performance summary, 2026-09-11](validation/templates_ldmos_stage4_performance_summary_2026-09-11.md):
  implemented optimizations, Release verification, and the validated Vg=8 V
  linked-continuation experiment with its reproduction limits.
- [LDMOS current validation](validation/templates_ldmos_current_status.md):
  current qualification boundaries, configuration profiles and reproducible
  linked D4/D5 entry points.
- [IALMob explicit and generated partials, 2026-09-17](validation/templates_ldmos_symbolic_partials_2026-09-17.md):
  high-field analytic partials, branch-preserving low-field code generation,
  independent AD/finite-difference checks and frozen-seed/curve controls; implementation-stage evidence.
- [IALMob kernel costs and same-VM follow-up, 2026-09-17](validation/templates_ldmos_symbolic_vm_2026-09-17.md):
  Release kernel microbenchmarks and deadline-protected full-curve controls; qualification tracked separately from R10.
- [Explicit high-field full repeated qualification, 2026-09-17](validation/templates_ldmos_symbolic_full_2026-09-17.md):
  recommended explicit R11 profile; 12 full curves, exact repeated trajectories, real restart and 62-point physical audit; R10 rollback.
- [R11 G3/D5/D4 isothermal full sweeps, 2026-09-17](validation/templates_ldmos_r11_isothermal_2026-09-17.md):
  final-source isothermal regression: all 155 exact points and both dual-gate joint scores pass; original physics/gates, full continuation distinct from saved-point reclosure.
- [R9 UCRT64 first8 gprof, 2026-09-17](validation/templates_ldmos_gprof_r9_2026-09-17.md):
  Vg8 Release control and cold/drain-only profiles; exact trajectories, IALMob/assembly hotspots and profiler limits.
- [R9 preparation and kernel follow-up, 2026-09-17](validation/templates_ldmos_preparation_followup_2026-09-17.md):
  immutable preparation reuse, thermal high-field reuse controls and fixed-seed candidate-cost diagnosis.
- [Preparation full VM repeat study, 2026-09-17](validation/templates_ldmos_preparation_full_2026-09-17.md):
  qualified explicit R10 preparation profile; two full dual-gate rounds, native pairing and actual pause/resume.
- [LDMOS work summary, 2026-09-16–17](validation/templates_ldmos_work_summary_2026-09-16_17.md):
  daily-report material covering contact repair, root controls, hotspot profiling and R10 qualification.
- [Safeguarded neutral-root Newton, 2026-09-16](validation/templates_ldmos_neutral_root_newton_2026-09-16.md):
  original-equation scalar root acceleration, finite-budget rounding guard and isolated representative controls.
- [Newton and assembly cost controls, 2026-09-16](validation/templates_ldmos_newton_cost_2026-09-16.md):
  explicit R9 contact recommendation, exact neutral-root caching and failed local QF clipping controls.
- [LDMOS production electrothermal reproduction](validation/templates_ldmos_production_reproduction.md):
  recommended explicit R11 high-field profile, R10 rollback and compatible R7 export; Release UMFPACK execution and qualification limits.
- [Extrapolation and continuation research, 2026-09-15](validation/vela_extrapolation_solver_research_2026-09-15.md):
  primary papers and open-source algorithms, with proposed R7/R8 predictor and Newton experiments.
- [Extrapolation A–E execution, 2026-09-15](validation/templates_ldmos_extrapolation_execution_2026-09-15.md):
  guarded predictor, step-growth, tangent, local-fit and NGMRES controls; complete costs,
  frozen evidence and negative results, with R7 retained as the production default.
- [Newton update execution, 2026-09-15](validation/templates_ldmos_newton_update_execution_2026-09-15.md):
  F0 reproduces all 74 R7 attempts exactly; local density projection and NLEQ_ERR-type controls
  fail promotion, and conservative Jacobian refresh shifts cost to extra updates/assemblies.
  R7 remains the production default; projected pseudo-transient mass-matrix research is recorded.
- [Carrier pseudo-transient validation, 2026-09-15](validation/templates_ldmos_pseudo_transient_2026-09-15.md):
  storage signs, temperature/statistics partials and SG diffusion verified; twelve controls
  fail promotion, with fixed-state coefficient isolation identifying non-descent steady-residual directions.
- [Pseudo-transient acceptance validation, 2026-09-15](validation/templates_ldmos_pseudo_acceptance_2026-09-15.md):
  actual backward-Euler defect trials cross initial rejection, but both time controllers fail
  steady gates within 60 updates at each frozen point; R7 remains unchanged.
- [Poisson-consistent initialization controls, 2026-09-15](validation/templates_ldmos_poisson_initialization_2026-09-15.md):
  both preparations close with QF/temperature held; model-controlled PTC reaches near-steady
  states but fails original carrier-row gates, with all preparation costs retained. No promotion.
- [Near-steady mass, density and QF representation controls, 2026-09-15](validation/templates_ldmos_near_steady_isolation_2026-09-15.md):
  state-only local reference changes plus one QF update close both frozen failures under original gates;
  density mapping adds rejected trials. That isolation stage did not test automatic switching or full-path acceleration.
- [Automatic near-steady switch, 2026-09-15](validation/templates_ldmos_near_switch_2026-09-15.md):
  opt-in switching closes both full representative trajectories after Poisson preparation in 32/31 total
  updates; direct raw seeds still fail. Original R7 needs 15/15, so the candidate is not promoted.
- [Standalone R7 high-voltage QF rebase, 2026-09-15](validation/templates_ldmos_r7_local_rebase_2026-09-15.md):
  four points with two paired repeats on each platform pass original gates; Windows/Linux updates
  remain 44/45 per round, with extra reassembly and no reduction of the Linux five-update floor tail.
  Diagnostic only; production defaults unchanged.
- [Finite-contact row arithmetic audit, 2026-09-16](validation/templates_ldmos_contact_row_roundoff_2026-09-16.md):
  exact frozen R7 prefixes isolate the high-voltage floor tail to neutral-contact potential consistency,
  amplified by finite hole exchange; high-precision sum/SG/recombination arithmetic does not remove it.
  Existing contact projection closes the five failed frozen states under original row/block gates;
  that audit alone did not qualify an automatic Newton candidate.
- [Finite-contact algebraic consistency candidate, 2026-09-16](validation/templates_ldmos_contact_consistency_2026-09-16.md):
  an opt-in terminal repair passes 32 paired original-seed high-voltage runs on Windows/Linux;
  Linux's 14-update tail becomes 9 under unchanged gates. Per-round updates fall 44/45 to 40/40,
  but Linux representative-point wall time reverses between repeats. Subsequent full-curve
  qualification is linked below; production R7 defaults are unchanged.
- [Contact-repair sweep validation, 2026-09-16](validation/templates_ldmos_contact_sweep_2026-09-16.md):
  first8 dual-gate controls and real pause/resume pass with exact R7 states and unchanged
  152/137 drain updates; initialization is explicitly protected. All 62 points and two rounds
  of R7/candidate/native controls pass the original joint gates with exactly repeated trajectories.
  Updates fall 366/321 to 344/302; median wall time falls 3.53%/7.20%, with variable per-round
  benefit and higher CPU cost than native. The candidate remains opt-in.
- [Newton update-mapping and damping plan, 2026-09-15](validation/templates_ldmos_newton_update_plan_2026-09-15.md):
  R7 iteration-phase decomposition, the R8 global positivity-alpha failure mechanism, and the
  proposed F0–F5 stage contracts and confirmed design choices (per-node density projection first);
  implementation and local VM experiments are authorized; see the execution report for completed controls.
- [Linux/UCRT64 environment alignment study](validation/vela_linux_ucrt_environment_alignment_2026-09-14.md):
  measured package versions, Linux dependency dry-runs, source-build gaps and unqualified deployment plans.
- [LDMOS daily report, 2026-09-12](validation/templates_ldmos_daily_report_2026-09-12.md):
  physics and electrothermal qualification summary with reproducible static
  current, temperature-rise and nodal-temperature comparison figures.
- [LDMOS daily report, 2026-09-13](validation/templates_ldmos_daily_report_2026-09-13.md):
  finite contact and Auger alignment, native Poisson geometry qualification,
  and serial performance results with reproducible comparison figures.

Optional feature switches used by this repository:

- `VELA_ENABLE_HDF5` (default ON): enables Sentaurus inventory/export support
  when an HDF5 package is found by CMake.
- `VELA_ENABLE_PYTHON` (default OFF): enables the pybind11 Python module and
  `python_api` CTest target.

See `CMakePresets.json` for the shipped Windows UCRT64 preset combinations.

## External Fixture And Test References

- [Regression README](../tests/regression/README.md): engineering regression
  runner behavior, summary JSON fields, and assertion configuration.
- [Reference TCAD README](../reference_tcad/README.md): neutral CSV export
  format and comparison workflow.

## Evidence And History Policy

Dated files under `validation/` are immutable point-in-time evidence and can
contain decisions that were later superseded. Read the current validation
summary first, then follow its links to the relevant evidence. Design specs and
the small number of execution plans still referenced by retained evidence are
provenance records, not active work queues.

Treat the root `README.md`, this index, `config_schema.md`, checked-in
`reference_tcad/` fixtures, CMake targets, tests, and current source code as
the source of truth for current behavior.
