# Vela Documentation

- [SLOT main合并与BVMethods外电阻重跑，2026-11-09](validation/slot_ldmos_main_bvmethods_restart_2026-11-09.md)：保留原SLOT适配、算子默认隔离、验证结果与两小时独立参考配置重跑。

- [SDE 当前交付与验收，2026-11-01](sde_delivery.md)：main 八例建模、Ubuntu CI、离线包归档与历史商业证据边界。
- [BVMethods 正式入口四项补齐，2026-11-09](validation/bvmethods_production_closure_2026-11-09.md)：
  普通入口、当前分叉证据、原始 CMD 整链与跨阶段账本；各项终态以报告为准。

- [BVMethods 正式入口集成与回归，2026-11-08](validation/bvmethods_production_integration_2026-11-08.md)：
  Unit 坐标、保护后偏导事件、近似判停和 CMD/常数合同的正式接入；验证状态与保留限制。

- [低压失稳问题归档与剩余主线，2026-11-08](validation/bvmethods_low_failure_archive_2026-11-08.md)：
  研究暂缓、根因和门禁资格保留，回到正式默认集成与回归；不以原生同类现象豁免性能条件。

- [原生低压失败器件/电路响应，2026-11-08](validation/bvmethods_native_low_failure_response_2026-11-08.md)：
  8份冻结系统与少子聚焦归因，区分绝对浓度分量和相对人口方向；四开关恢复未放行，无生产改动。

- [原生低压失败诊断补采，2026-11-08](validation/bvmethods_native_low_failure_capture_2026-11-08.md)：
  同配置前17目标复现，8份矩阵/方向与7组实际写回；自由浓度1%保护得到确认，直接缩步机制仍不同，共同根因未闭合。

- [原生约2.51 V与Vela约2.43 V失败过程对照，2026-11-07](validation/bvmethods_native_vela_low_failure_comparison_2026-11-07.md)：
  逐轮RHS/E、保护与五帧Schur/RHS分块；原生RHS增长终止、Vela写回前绝对线性门拒绝，巨大方向来自器件块，未认定同一上游根因。

- [2.4307 V绝对线性精度定位与有界试修，2026-11-07](validation/bvmethods_defaults_absolute_precision_2026-11-07.md)：
  大空穴方向放大舍入，原坐标/Full及补偿校正均未找到通过原门的修复，保留缩步、生产不改。

- [负载线RHS修复与完整BV复验，2026-11-07](validation/bvmethods_defaults_rhsfix_fullbv_2026-11-07.md)：
  25/40组前置及两组完整BV通过；Unit377轮/4失败，消除原RHS后向拒绝，保留1次绝对门拒绝和默认放行限制。

- [全变量Unit高压回代拒绝定位，2026-11-07](validation/bvmethods_defaults_backward_audit_2026-11-07.md)：
  22拒绝与3对照逐位复现，定位电路RHS重复求值舍入；统一RHS冻结反事实25/25过原门，未改生产。

- [默认迁移完整BV复核，2026-11-07](validation/bvmethods_defaults_fullbv_2026-11-07.md)：
  原坐标/全变量Unit两组停流、BV插值及逐点质量/物理输出与终态证据复核完成，实际资格见报告，生产默认保持。

- [默认迁移0–100 V分区工作量，2026-11-07](validation/bvmethods_defaults_regional100_2026-11-07.md)：
  原坐标/全变量Unit连续扫描、同版本原生参照及必要新失败入口配对，计数/锚点/成本与质量独立核验，生产默认未切换。

- [默认迁移5–10 V缩步及线性拒绝检查，2026-11-07](validation/bvmethods_defaults_mid10_2026-11-07.md)：
  同核原坐标/全变量Unit连续0–10 V及失败入口配对，保留失败、独立质量与逐SHA终态证据，生产默认不变。

- [默认迁移0–1 V首次分叉检查，2026-11-07](validation/bvmethods_defaults_low1_2026-11-07.md)：
  原坐标/全变量Unit与最低校正四臂均到1 V，无失败；同入口2×2重放定位入口和保护路径依赖，生产默认未切换。

- [VM固定步长0–10 V Full AND/OR与原生对照，2026-11-05](validation/bvmethods_vm_fixed_full_scan_2026-11-05.md)：
  0.1/0.02 V固定网格均保留提前失败；92份同态J/RHS、24份同目标物理场及逐目标轮数，116份冻结装配和终态证据核验完成，未获得完整扫描资格。

- [兼容常数VM生产A0的N20冷0–1 V对照，2026-11-04](validation/bvmethods_vm_iteration20_low1_2026-11-04.md)：
  耦合上限与步长工作量预算40→20；AND/OR均到1V，OR工作量下降、AND增加，保留原生首次路径分叉和独立资格。

- [VM全变量Unit坐标重构验证，2026-11-04](validation/bvmethods_vm_unit_basis_2026-11-04.md)：
  全变量坐标、30组冻结态及固定/0–1 V对照；原门保持，方向门未全过且工作量增加，未切换生产。

- [Vela与原生有效缩放规则核查，2026-11-03](validation/bvmethods_native_scaling_audit_2026-11-03.md)：
  一个Jacobian校准、14帧留出与独立TDR预测；区分变量/方程/线性平衡，不改生产门。

- [VM近零电流AV导数与停止规则审计，2026-11-03](validation/bvmethods_vm_av_stopping_2026-11-03.md)：
  12行250位AV差分、15份原生更新误差重建、影子与同态/区间配对；RHS分支未冒称等价；有界对照、暖态配对和终态证据已完成，监测已停止。

- [兼容常数VM生产A0完整BV，2026-11-03](validation/bvmethods_vm_sentaurus2022_full_2026-11-03.md)：
  生产A0兼容常数AND/OR完整扫描终态均通过；保留共同区间、成本及报告型守恒比较。

- [兼容常数T470p冷0–100 V分区比较，2026-11-03](validation/bvmethods_t470p_sentaurus2022_low100_2026-11-03.md)：
  按实际外压统计接受点、目标尝试和失败；AND提前数值失败，OR到100 V通过，终态与证据已核验。

- [修复核心A1三组、入口与平台定位，2026-11-02](validation/bvmethods_repaired_a1_t1_t3_2026-11-02.md)：
  T1/A1e与T2/T3诊断结束；保留低区反例、入口恢复及投影停滞证据，原A1完整BV准入仍阻断。

- [U-2022.12常数兼容配置与T470p验证，2026-11-03](validation/bvmethods_sentaurus2022_constants_2026-11-03.md)：
  默认CODATA2018保留，显式兼容配置与材料重生成；15帧RHS明显对齐，保留首帧六空穴行导数差异；默认全回归1219项零失败。

- [原生U-2022.12参数/TDR与基本常数核查，2026-11-03](validation/bvmethods_native_parameter_constants_audit_2026-11-03.md)：
  真正原生par导出与同版本手册旧常数；无拟合重建2152节点QF至约1e-15 V，生产设置不改。

- [Vela VM原生首五点矩阵/RHS系数偏差定位，2026-11-02](validation/bvmethods_vm_static_operator_bias_2026-11-02.md)：
  Poisson与DOS/kT系数识别、独立RHS反事实和AV直接分项留出验证，生产模型不改。

- [SRH修复版S2–S5接触/端口与保护/高压复验，2026-11-02](validation/bvmethods_s2_s5_2026-11-02.md)：
  接触与统计反事实、250位SRH分项、36份Schur及保护历史复验；新AND/OR全程首步统计完成。

- [原生首五点固定方向差与根/轨迹复验S0–S1，2026-11-02](validation/bvmethods_native_roots_s0s1_2026-11-02.md)：
  五份真正原生末态各2次Full通过AND；Vela冷首五目标15接受/15尝试/0失败，限定为局部根偏移。

- [T470p SRH修复版与原生首五点逐轮矩阵/RHS对照，2026-11-02](validation/bvmethods_t470p_native_first5_srh_comparison_2026-11-02.md)：
  15帧Full装配与原始方向、6765维共同系统及边界仿射RHS分解，保留跨工具一致性限制。

- [SRH修复版VM完整BV AND/OR，2026-11-02](validation/bvmethods_vm_srh_fixed_full_2026-11-02.md)：
  AND/自动复核OR完整扫描均通过独立终态验收，保留分区成本、报告型守恒和跨版本比较限制。

- [SRH修复版T470p Full低压AND/OR对照，2026-11-02](validation/bvmethods_t470p_srh_full_low_2026-11-02.md)：
  原参数冷0–1 V、始终完整导数，两组均线性精度失败，保留共同外压范围点数/失败对照。

- [BVMethods物理算术修复与T470p复验，2026-11-02](validation/bvmethods_physics_arithmetic_fix_2026-11-02.md)：
  人口SRH、AV微场及通用复合Full导数修复；六冻结态和原门冷窗口通过，SG默认不变。

- [BVMethods人口坐标下物理模型算术与导数审计，2026-11-02](validation/bvmethods_physics_coordinate_audit_2026-11-02.md)：
  Auger共享excess/限幅、AV微场AD边界及SG派生QF消减的隔离复现，区分实际启用与人工控制。

- [BVMethods近似模型可信度与Schur窗口消融，2026-11-02](validation/bvmethods_bootstrap_trust_windows_2026-11-02.md)：
  36份合格Schur分解、30组消融与12组选择器配对窗口终态；保留中低区失败，P3未启动。

- [BVMethods保护后逐轮路径与0.51 kV失稳审计，2026-11-02](validation/bvmethods_post_protection_path_audit_2026-11-02.md)：
  同态项族／候选复核／近解续算，定位AV切线遗漏与持续局部冻结，并保留算术敏感性限制。

- [BVMethods低偏压反复保护与阶段线性化诊断，2026-11-02](validation/bvmethods_protection_ledger_2026-11-02.md)：
  逐节点保护、同态切线干预、E分解、20/40轮与高偏压反例，保留未闭合归因。

- [BVMethods与Sentaurus区域采样和BV附近步长对照，2026-11-02](validation/bvmethods_native_step_comparison_2026-11-02.md)：
  原生日志/PLT与AND、OR完整轨迹逐次对齐，定位100–1000 V密集尝试及BV前后步长变化。

- [BVMethods两组OR策略完整BV，2026-11-01](validation/bvmethods_or_full_2026-11-01.md)：
  同二进制Full/自动复核并行OR、逐点独立强制门与八小时监督。

- [BVMethods OR停止与历史AND轨迹复核，2026-11-01](validation/bvmethods_or_stopping_review_2026-11-01.md)：
  相同Full二进制与packed初态的AND/OR实测、残差/更新两分支及历史轨迹的反事实边界。

- [BVMethods始终完整Jacobian消融与完整BV资格，2026-11-01](validation/bvmethods_full_jacobian_ablation_2026-11-01.md)：
  同二进制三组P0、只报告ρ、入口逐位核验与有界P1资格。

- [BVMethods bootstrap保护后模型复核与显式历史验证，2026-11-01](validation/bvmethods_bootstrap_model_review_2026-11-01.md)：
  显式C++复核、共享探测额度、真实历史分组及VM窗口资格。

- [BVMethods 2343 V方向、保护与隔离策略审计，2026-11-01](validation/bvmethods_direction_projection_audit_2026-11-01.md)：
  实际预测态重放、合格Schur重建、保护后模型、全局化/历史控制及分项计时。

- [BVMethods bootstrap线性失败同态完整导数恢复，2026-10-31](validation/bvmethods_bootstrap_linear_restore_2026-10-31.md)：
  核心单次恢复、防循环与原门限回归；独立VM同态预检和完整BV监督。

- [BVMethods失败增广系统冻结与线性解法对照，2026-10-31](validation/bvmethods_frozen_linear_failure_2026-10-31.md)：
  同态矩阵逐行精度、16项后端控制、完整导数单步通过及bootstrap线性失败恢复缺口。

- [BVMethods SG统计反馈局部bootstrap核心控制，2026-10-31](validation/bvmethods_sg_bootstrap_core_2026-10-31.md)：
  显式局部边选择、成对冻结及恢复契约；原生逐轮和多步暖态复验。

- [BVMethods弱人口列的SG统计反馈定位，2026-10-31](validation/bvmethods_sg_statistical_tangents_2026-10-31.md)：
  统计／QF／SRH切线分解、180位复核与局部相邻边控制，恢复六个原生弱人口方向。

- [BVMethods线性化策略与阶段历史控制，2026-10-31](validation/bvmethods_linearization_policies_2026-10-31.md)：
  四类导数／事件切换／整份旧矩阵复用、76项控制及被整体误差掩盖的弱人口列。

- [BVMethods下限RHS与原生逐轮有效矩阵，2026-10-31](validation/bvmethods_floor_rhs_effective_matrix_2026-10-31.md)：
  输运主导的下限平台、迁移率切线与30项保护历史关联、同态Coupled阶段拆分控制。

- [BVMethods失败态Newton定位，2026-10-31](validation/bvmethods_newton_localization_2026-10-31.md)：
  保正相对人口导数、四策略/20与40轮控制、原生实际更新及弱行有效矩阵与RHS。

- [BVMethods完整耦合停止修复，2026-10-31](validation/bvmethods_coupled_stop_fix_2026-10-31.md)：
  C组提前接受修复、设备/电路共同残差、候选恢复诊断与三组VM重跑。

- [BVMethods 停止规则 A/B/C 三组对照，2026-10-30](validation/bvmethods_stopping_controls_full_2026-10-30.md)：
  OR／AND 与负载线 enforce／report 独立配置、核心回归和三小时 VM 完整扫描。

- [BVMethods残差范数与停止规则复核，2026-10-30](validation/bvmethods_residual_stopping_2026-10-30.md)：
  原生打印范数实测、RMS计数与原始方向契约、OR/AND停止逻辑和当前失败态逐项判断。

- [BVMethods新模式18点与Sentaurus日志对照，2026-10-30](validation/bvmethods_quasistationary_first18_2026-10-30.md)：
  逐点I–V／步长／迭代、末段导数恢复周期、局部空穴保护和停止规则差异。

- [BVMethods Quasistationary新模式VM完整扫描，2026-10-30](validation/bvmethods_vm_quasistationary_full_2026-10-30.md)：
  独立源码构建、原验收门冷Poisson起步、一小时完整扫描监督及状态记录。

- [BVMethods外接电阻Quasistationary步长接入，2026-10-30](validation/bvmethods_quasistationary_steps_2026-10-30.md)：
  外压自适应、无额外内压限幅、CMD/JSON配置及兼容回归。

- [BVMethods VM双0.5 V完整扫描，2026-10-30](validation/bvmethods_vm_inner_0p5_full_2026-10-30.md)：
  T470p对照迁移至VM，复用已验证二进制，独立冷Poisson起步、成对0.5 V上限与一小时截止。

- [BVMethods T470p双0.5 V内压控制对照，2026-10-30](validation/bvmethods_t470p_inner_0p5_full_2026-10-30.md)：
  单次Newton及预测内压上限同步设置、一小时远端扫描与独立预检；资格以实际结果为准。

- [BVMethods归一化状态及热路径修复，2026-10-30](validation/bvmethods_normalization_fix_2026-10-30.md)：
  原始psi坐标、严格重启、QF／浓度往返和Poisson归一化差分；资格以验证记录为准。

- [BVMethods单位与归一化往返审计，2026-10-30](validation/bvmethods_units_normalization_audit_2026-10-30.md)：
  检查状态装拆包、重启、热路径QF和浓度SI往返；量化真实检查点1ULP变化，并区分方程预处理与状态权威坐标。

- [BVMethods SRH稳定导数修复与VM多步复验，2026-10-30](validation/bvmethods_srh_stable_derivative_2026-10-30.md)：
  核心SRH切线修复、123例回归及同态RHS不变；tail02改善为13轮，但固定第五目标及旧轨迹冻结态仍失败，已归档packed坐标算子。

- [BVMethods VM 第二次 bootstrap 六个电子弱行，2026-10-30](validation/bvmethods_vm_second_bootstrap_rows_2026-10-30.md)：
  固定尺度 RHS／直接人口列与独立精度复核；定位 SRH 统计导数相消及 SG/GSS 弱支撑，保持完整 RHS 和原验收门。

- [BVMethods VM 边界目标与材料统计一致性，2026-10-30](validation/bvmethods_vm_boundary_consistency_2026-10-30.md)：
  一致初始化通过冷启动首五目标，tail02迭代投影后续算五点；定位浮点级边界差触发的第二轮保护分岔，独立闭合Poisson常数／行尺度隔离。

- [BVMethods VM 接触人口、参考量与 tail02 弱行对照，2026-10-30](validation/bvmethods_vm_boundary_contract_2026-10-30.md)：
  同目标浓度式／QF式边界、56套核心方向、边界仿射传播和材料参考核验；自身边界投影后tail02 34轮通过，绝对物理行尺度仍未闭合。

- [BVMethods VM 原生首五点、首次保护与 tail02 对照，2026-10-30](validation/bvmethods_vm_native_lockstep_2026-10-30.md)：
  同目标自主复算、50套核心人口保护核验、入口下限／Power控制与同态矩阵和RHS交叉求解；首五目标与完整BV资格分别记账。

- [BVMethods 内压步长数学分析与T470p冻结定位，2026-10-30](validation/bvmethods_t470p_step_math_2026-10-30.md)：
  Wolfram外推/载流子尺度核验、局部1%保护下溢定位、三种Newton上限同态对照及弱行交叉Jv。

- [BVMethods 核心稳定缩放与直接人口列线性化，2026-10-30](validation/bvmethods_relative_columns_core_2026-10-30.md)：
  核心乘积与故障分类修复、相对列全链路、真实冻结矩阵回归、局部续算与Sentaurus浓度上限核查。

- [BVMethods 弱人口行与新零列同态定位，2026-10-30](validation/bvmethods_weak_population_rows_2026-10-30.md)：
  原始非零列的缩放中间溢出、三组冻结PARDISO复算，以及tail02最差人口节点与原始Newton方向。

- [BVMethods 浓度下限与外推阈值独立验证，2026-10-30](validation/bvmethods_density_controls_2026-10-30.md)：
  四组12项VM对照、format v2单位修正、冷窗口逐位一致，以及困难态残差通过而人口更新误差仍失败的边界。

- [人口主变量限制与停止条件审查，2026-10-29](validation/potential_density_policy_audit_2026-10-29.md)：
  清理QF专用配置，统一入口校验，修正独立人口Newton的权威分块停止契约。
- [T470p BVMethods第二轮热点优化，2026-10-29](validation/bvmethods_t470p_round2_2026-10-29.md)：
  同态诊断融合与Lombardi常系数预计算，普通Release 31.92→24.97 s，数值轨迹保持一致。
- [T470p当前核心BVMethods gprof与优化，2026-10-29](validation/bvmethods_t470p_gprof_2026-10-29.md)：
  显式变量开关确认，节点统计与区域接口判断复用；普通Release短窗口配对、状态与算子一致性。

- [BVMethods自主原始量、边界人口与双机整理，2026-10-29](validation/bvmethods_autonomous_primitives_2026-10-29.md)：
  自主材料／统计／SG原型、完整RHS及Jacobian对照、边界与冷Poisson前五点；零偏弱空穴余差单列。

- [BVMethods微小驱动力与Cn表示差异，2026-10-29](validation/bvmethods_microdrive_2026-10-29.md)：
  重点边能带／统计／QF误差分解，Cn公共偏移相减的独立控制，以及当前源码VM三小时任务。

- [BVMethods 高精度统计、单元AV与弱行RHS误差账本，2026-10-28](validation/bvmethods_precision_budget_2026-10-28.md)：
  6个单元的AV方向贡献闭合、统计值／导数基准、非零边界仿射修正，以及原生off保留源项的口径更正。

- [BVMethods 原始量、JD 系数与 Poisson 基准，2026-10-28](validation/bvmethods_primitive_contract_2026-10-28.md)：
  独立能带／归一化契约、JD系数精度修正、完整Poisson残差隔离闭合，以及剩余AV／弱空穴差异。

- [BVMethods 人口 SG 算术路径与原生同态复核，2026-10-28](validation/bvmethods_population_sg_arithmetic_2026-10-28.md)：
  能带组合、统计／QF 与 Bernoulli 分项核验；显式直接人口 SG、零偏及有限态完整算子对照。

- [BVMethods 对数人口外推、无损重启与零偏 AV 同态复核，2026-10-28](validation/potential_density_predictor_archive_2026-10-28.md)：
  核心前五目标通过，JSON/CMD 与七点重启逐位一致；原生完整精度及共同 SG 电流诊断单列。

- [BVMethods 人口模式剩余问题定位，2026-10-28](validation/potential_density_remaining_diagnosis_2026-10-28.md)：
  原生对数外推、固定第五点隔离验证、材料／VT与零偏AV差异、人口单位往返重启对照。

- [BVMethods 电势／人口核心模式开发，2026-10-28](validation/potential_density_development_2026-10-28.md)：
  P0–P5 实施记录；显式 ψ/n/p、直接装配、局部保护、核心导数恢复及 JSON/CMD 路径，资格以实际测试为准。

- [BVMethods 50%人口保护多偏压续算，2026-10-28](validation/bvmethods_retained_continuation_2026-10-28.md)：
  双机共同通过0.8至1.25 V窗口，区分健康入口、持续保留比例及半步控制；全局到边界τ=0.99仍可能保留1%，未改默认或验收门。

- [BVMethods 局部人口保留1%／10%／50%对照，2026-10-28](validation/bvmethods_retained_fraction_2026-10-28.md)：
  双机50%均五轮通过同一尾部形成检查点，首步算子与原始方向保持，弱行平衡及全网格人口尾部改善；未修改默认比例。

- [BVMethods 尾部形成前候选进展与弱行导数，2026-10-28](validation/bvmethods_weak_row_jv_2026-10-28.md)：
  T470p／C560复核首次1%回退后的负输运斜率、弱行进展及近似矩阵相减精度；保持完整RHS与原验收门。

- [BVMethods 相对人口坐标完整扫描停滞定位，2026-10-27](validation/bvmethods_stable_population_failure_2026-10-27.md)：
  电子尾部与全局保护节点542、分块更新误差及候选/恢复选择链；同目标成功路径四步通过，保留物理导数核验边界。

- [BVMethods VM 自适应候选完整扫描失败定位，2026-10-27](validation/bvmethods_vm_adaptive_failure_2026-10-27.md)：
  节点498电子尾部使固定人口列尺度溢出；复用开关同态对照一致，区分中途缩步与最终不可继续接受态。

- [BVMethods 核心默认约束消元与 Jacobian 支撑对照，2026-10-27](validation/bvmethods_core_constraint_elimination_2026-10-27.md)：
  保留非零边界 RHS 和完整回代，实际自由系统 6765 维；15 套同态原装配保持，严格非零及显著支撑差异单列。

- [BVMethods 四项求值优化 VM 完整扫描完成，2026-10-27](validation/bvmethods_vm_four_switch_full_2026-10-27.md)：
  冷 Poisson、全局缩步与自动导数恢复，103 个接受点达到停止电流，正常退出；与后续新候选策略分开记账。
- [BVMethods T470p／C560 候选选择与迭代内复用，2026-10-27](validation/bvmethods_candidate_selection_2026-10-27.md)：
  按原验收尺度选取人口候选，两类起点均通过；普通 Release 配对验证共享计算收益及逐位一致性。
- [BVMethods T470p／C560 人口保护策略对照，2026-10-27](validation/bvmethods_population_protection_comparison_2026-10-27.md)：
  同态首轮算子一致；局部保护在预测态七轮通过，在接受态失败，保留原验收门及守恒报告。
- [BVMethods T470p／C560 同态求值复用优化，2026-10-27](validation/bvmethods_evaluation_reuse_2026-10-27.md)：
  冷启动前五点普通 Release 配对分别加速1.77×／2.26×；状态、轨迹、30个算子文件一致，显式开关启用。
- [BVMethods gprof 热点与诊断驱动成本，2026-10-27](validation/bvmethods_gprof_hotspots_2026-10-27.md)：
  C560 Release／gprof 前五点状态逐位一致；状态检查占55.48%，联立线性求解占1.15%，T470p复核进行中。
- [BVMethods vela VM 自动恢复策略完整扫描启动，2026-10-27](validation/bvmethods_vm_automatic_full_start_2026-10-27.md)：
  独立 Release／PARDISO 冷启动，22 项测试通过，首五个非零点接受；后台完整扫描启动快照，未宣称完成。
- [BVMethods vela VM 自动导数恢复切换验证，2026-10-27](validation/bvmethods_vm_derivative_switch_2026-10-27.md)：
  同态矩阵／RHS 保持一致，自动恢复及故障注入回退验证通过，随后五点各两次更新；仍为未安装诊断。
- [BVMethods vela VM 同检查点导数恢复验证，2026-10-26](validation/bvmethods_vm_restoration_verification_2026-10-26.md)：
  同态矩阵／RHS 一致；完整导数第 2 轮通过并复核至第 8 轮，原恢复选择持续返回近似阶段；未接入生产或完成 BV。
- [BVMethods vela VM 三小时对照结果分析，2026-10-26](validation/bvmethods_vm_predictor_analysis_2026-10-26.md)：
  两组按时限停止；仅报告版本推进更远，但导数恢复反复返回近似阶段，完整 BV 未完成。
- [BVMethods vela VM 预测内压限制双臂对照启动，2026-10-26](validation/bvmethods_vm_predictor_comparison_2026-10-26.md)：
  enforce／report_only 冷启动控制、独立三小时时限和启动时一致性核验。
- [BVMethods T470p 从零启动的自适应导数恢复试算，2026-10-26](validation/bvmethods_t470p_adaptive_bv_2026-10-26.md)：
  277 个接受点推进至外压 5.79877 V；挂起日志查询造成内存压力，完整 BV 未完成，重跑已暂停。
- [BVMethods 人口坐标、到边界步长与导数恢复试验，2026-10-26](validation/bvmethods_population_newton_trials_2026-10-26.md)：
  人口列通过物理差分核验；直接全局正值缩步未通过，分阶段恢复完整导数后续算五点；电子闭合量与方程残差之和吻合，完整 BV 尚未通过验证。
- [BVMethods 完整 RHS 下的线性化与人口更新控制，2026-10-26](validation/bvmethods_linearization_population_control_2026-10-26.md)：
  两个检查点的六组 Release 控制；近似矩阵配人口候选达到原停止门与报告 RMS＜1，守恒报告仍有超线项，尚不构成完整 BV 资格。
- [BVMethods 连续阶段 RHS 方向导数与有效线性系统，2026-10-26](validation/bvmethods_continuous_rhs_jv_2026-10-26.md)：
  17 组人口扰动支持完整导数；该步连续有效矩阵与冻结迁移率／AV 导数控制相符，不能将两类矩阵混为同一基准。
- [BVMethods 弱空穴线性化与人口保护联合核验，2026-10-26](validation/bvmethods_weak_hole_joint_2026-10-26.md)：
  固定输运尺度解释行 311 RHS 放大；原生同态阶段重启改变矩阵，差异关联受保护人口列；40 轮人口控制仍未收敛。
- [BVMethods 原生人口有限更新与弱空穴接受尺度，2026-10-26](validation/bvmethods_population_acceptance_2026-10-26.md)：
  原生四步及 1% 正值保护、RMS 停止量实测；Vela 八步人口控制仍未收敛，四态同人口矩阵差异保留。
- [BVMethods 原生阻尼配置对齐与失败归因，2026-10-26](validation/bvmethods_damping_failure_cause_2026-10-26.md)：
  CMD 默认直接完整步；空穴 RHS／AV 源分解与 T470p、C560 同检查点复现，完整 BV 仍未通过。
- [BVMethods 线搜索仅报告与完整步扫描，2026-10-26](validation/bvmethods_line_search_report_only_2026-10-26.md)：
  同 Jacobian/RHS 的首点从 15 轮降至 3 轮；自主完整步扫描在 0.216002 V 外压停止，原验收门保持。
- [BVMethods PARDISO、Poisson 与前五个联立点，2026-10-25](validation/bvmethods_pardiso_first5_2026-10-25.md)：
  显式原生 PARDISO 新运行；Vela 后端切换、七轮 Poisson、五点迭代及 15 套同态 Jacobian/RHS 对照。
- [BVMethods 历史完整扫描核对、极低人口修复与重新扫描，2026-10-25](validation/bvmethods_full_recovery_2026-10-25.md)：
  区分旧检查点路径和当前自主入口；修正极端 SG 比值/指数中间量以及 DOS 缩放统计精度，分别记录自主和重启验证。
- [BVMethods 关闭 QF 限幅与同态原始方向核验，2026-10-25](validation/bvmethods_qf_direction_2026-10-25.md)：
  对照历史测试，保留线搜索和验收门；记录弱少子方向、独立线性复算和极低浓度自动微分失效定位。
- [Windows h5py ROS3 DLL 修复，2026-10-25](validation/h5py_ros3_dll_repair_2026-10-25.md)：
  修复 h5py / HDF5 功能不匹配；原 8 项失败全部通过，完整 Release 为 1109 通过 / 0 失败 / 4 项参考数据跳过。
- [BVMethods JSON / CMD 三步一致性验证，2026-10-25](validation/bvmethods_cmd_e2e_2026-10-25.md)：
  最新 Release、实际 CMD check 与自主执行；13 个接受点 / 39 次尝试和基础场精确一致，完整 BV 仍在低偏最小步长处失败。
- [BVMethods Release 完整 BV 扫描，2026-10-25](validation/bvmethods_full_profile_2026-10-25.md)：
  当轮自主 Poisson 起步及严格残差续算达到停扫电流；BV 约 6.67049 V，采用私有续算驱动。
- [BVMethods 显式算子契约](bvmethods_operator_contract.md)：
  自主 QF 状态、统一统计／FD-BGN、SG 算术、残差／矩阵／端口入口及适用边界。
- [BVMethods 显式算子接入验证，2026-10-24](validation/bvmethods_production_profile_2026-10-24.md)：
  Release 同态检查点、有限电流对照及多步负载线回归；零偏原生方向与完整 BV 资格单列。

- [BVMethods 零偏 AV 电流方向与同态 RHS／Jacobian 核验，2026-10-24](validation/bvmethods_zero_av_direction_2026-10-24.md)：
  共同原始量复算定位 SG 舍入路径；AV 源及完整矩阵差闭合至约 1e-12，自主 QF 原臂及完整 RHS 精度边界单列。

- [BVMethods 有限电流完整算子及非 AV 差异定位，2026-10-24](validation/bvmethods_finite_transport_operator_2026-10-24.md)：
  闭合输运平均权重及 SRH 体积；原生有限检查点完整矩阵差降至约 1e-8，Poisson 常数组合控制与零偏微小电流资格单列。

- [BVMethods AV 单元材料、源残差与完整方向导数契约，2026-10-24](validation/bvmethods_av_material_contract_2026-10-24.md)：
  逐顶点迁移率独立闭合，有限激励完整 Jv 与首个有限 Newton 步完成；零偏敏感性及约 0.8 mV 的原生首步差仍未消除。

- [BVMethods 原生实际 AV 源项单侧扰动，2026-10-24](validation/bvmethods_native_av_source_response_2026-10-24.md)：
  代表点两邻边恢复及盒体积分配通过直接源核验；非直角单元区分三边候选，原生 AV 大矩阵项可重建，迁移率支撑仍待闭合。

- [BVMethods 原生 AV 导数分离与单侧残差核验，2026-10-24](validation/bvmethods_native_av_derivative_split_2026-10-24.md)：
  同检查点显式开关确认异常交叉项来自 AV；实际残差区分源支撑缺失与近零绝对电流分支抵消，原生单元权重尚待确定。

- [BVMethods 同态 Jacobian 对比，2026-10-24](validation/bvmethods_jacobian_comparison_2026-10-24.md)：
  差异集中于四个空穴行的电子／内压交叉项；方向分解及源 Jacobian 控制完成，原生 AV 分项仍待分离。

- [BVMethods 匹配基准后的 RHS 与首步更新，2026-10-24](validation/bvmethods_rhs_first_update_2026-10-24.md)：
  补齐同统计／材料／热电压的完整系统；6765 维 RHS 差主要在 Poisson，首步仍未一致，独立线性复核完成。

- [BVMethods 一致统计契约与零偏回放，2026-10-24](validation/bvmethods_statistics_contract_2026-10-24.md)：
  实验解析契约通过独立原生状态与 FD-BGN 控制；同输入源项显著降低，生产零偏及极小电流资格仍单列。

- [BVMethods 能带与费米统计直接核验，2026-10-24](validation/bvmethods_energy_statistics_2026-10-24.md)：
  新原生 Poisson／联立能带确认默认电子统计差约 0.2186%；能量参考与 ni 带隙闭合，FD-BGN 差异单列。

- [BVMethods 碰撞电离公式与代码核查，2026-10-24](validation/bvmethods_avalanche_formula_audit_2026-10-24.md)：
  全边源公式与成对装配复现；定位费米／FD-BGN 本构不一致及接触驱动场激活，低场指数底值另列。

- [BVMethods QF 更新上限依据与取消对照，2026-10-24](validation/bvmethods_qf_update_limit_control_2026-10-24.md)：
  0.1 V 是可调稳定化控制；异常导入态取消后回溯失败，自身短扫已保存轨迹相同，尚无完整 BV。

- [BVMethods 首轮空穴大方向的碰撞电离源定位，2026-10-24](validation/bvmethods_fp_avalanche_cause_2026-10-24.md)：
  十组单步控制确认电子碰撞电离空穴源主导；导入态非零电流与接触邻域电场替代共同触发，修复资格仍待验证。

- [BVMethods 共同 RHS 范数与前五轮更新对照，2026-10-24](validation/bvmethods_common_norm_updates_2026-10-24.md)：
  原生日志 RHS 的无穷范数复现；完整物理场匹配、6765 维共同归一化与 QF 增量对照，定位空穴少子方向及实际限幅差异。

- [BVMethods 停止准则与守恒诊断分离，2026-10-24](validation/bvmethods_stopping_control_2026-10-24.md)：
  增加可选浓度 RMS 或残差停止；同态严格残差对照通过，三组扫描各保留 22 条完整记录后超时，尚无 BV。

- [BVMethods 本地 Release 完整碰撞电离扫描尝试，2026-10-23](validation/bvmethods_full_impact_release_2026-10-23.md)：
  修复材料后实际运行三组严格控制；30 个接受偏压点，电子全局闭合阻断，尚无合格 BV。

- [SDevice 电子 DOS 有效质量修复，2026-10-23](validation/sdevice_electron_dos_mass_fix_2026-10-23.md)：
  OldSlotboom 常数偏移进入电子质量带隙比；保留 CODATA 常数，材料/前端回归与零偏冒烟通过。

- [BVMethods 原生材料参数新导出，2026-10-23](validation/bvmethods_material_export_2026-10-23.md)：
  `sdevice -P` 与新 TDR 复现 Nc/Nv/ni 诊断值；确认活动系数一致，差异位于 BGN 质量换算和常数。

- [BVMethods 分阶段零偏对照，2026-10-23](validation/bvmethods_staged_zero_2026-10-23.md)：
  独立 Poisson→联立对照定位几何/电荷体积差异；7801→6843 等价消元通过，场量控制通过，算子仍未一致。

- [PN2D 分阶段零偏对照，2026-10-23](validation/pn2d_staged_zero_2026-10-23.md)：
  按 Poisson→联立顺序独立求解；匹配材料/常数后终态吻合，有限更新和生产材料修复仍单独记账。

- [PN2D 零偏高场与 debug 精度跟进，2026-10-23](validation/pn2d_zero_hfs_followup_2026-10-23.md)：
  有限 RHS 和首个有限 RHS 状态的近似切线已对齐；较后少子高场、弱行闭合与写回资格仍未通过。

- [PN2D 零偏 Jacobian/RHS 对照，2026-10-23](validation/pn2d_zero_operator_2026-10-23.md)：
  当前先完成零偏；材料 DOS 质量换算及迁移率平均差异已分离，全模型一致性尚未通过。

- [PN2D 原生外电阻矩阵对照，2026-10-23](validation/pn2d_external_resistor_native_matrix_2026-10-23.md)：
  Sentaurus 调试轨迹成功；5762 维等价消元、人口/QF 变换与结区迁移率、接触/RHS、停止语义定位。

- [PN2D 外接电阻法起步验证，2026-10-23](validation/pn2d_external_resistor_start_2026-10-23.md)：
  本地 main 合并、负向外压入口修复及旧 v3 净掺杂配置证据；当前总杂质修正与重跑见原生矩阵对照。
- [BVmethods 原脚本自主建模与 P5 合同](sde_bv_modeling.md)：创建态栅极移除、multibox、Erf、x 镜像及 live Sentaurus 语义对照。
- [MOS 原脚本自主建模与 P4 合同](sde_mos_modeling.md)：圆角弦误差、快捷加密、半/全镜像、独立导入及明确消费端拒绝边界。

- [BJT 原脚本自主建模与 P3 合同](sde_bjt_modeling.md)：General 受限 AST、有界 MaxTransDiff、掺杂积分和结位置。

- [自产 SDE 网格进入 sd-sim 的 P2 检查](sde_sdevice_checks.md)：8份具名命令、接口边界分类修复、输入/材料/sidecar/阶段身份与明确拒绝门禁。

- [离线 SDE 案例回归与使用说明](reference_sde_cases.md)：八例自产网格合同与原生身份检查、重复运行及失败恢复。

- [SDE 建模与 GRD/DAT/MSH 后续计划，2026-10-30](superpowers/plans/2026-09-05-reference-sde-modeling-development-plan.md)：
  基于已合并的 main 基线；主案例 8/8 建模导出及原生导入，P1–P5 完成；P6 本地跨平台资格通过，真实 CI 按用户选择未执行。

- [P6 离线建模包与本地双平台资格](sde_p6_qualification.md)：最终包 32 次主案例运行、32 次商业网格对照、
  8 次最终 TDX 回转和两平台 P2 命令检查通过；原始 TDR/日志不提交，真实 CI 未执行。

- [PARDISO 编译开关与 HFS 测试耗时分析，2026-10-23](validation/pardiso_build_option_hfs_audit_2026-10-23.md)：
  PARDISO 默认编译关闭、GitHub CI 显式关闭；HFS 高精度 Jacobian 检查的工作量与双机 OFF／ON 验证。

- [PARDISO main 接入与双机验证，2026-10-23](validation/pardiso_cross_host_validation_2026-10-23.md)：
  T470p／C560 完成匹配源码的 Debug 构建、全套回归失败项复核、MKL 多线程契约与四组数值回归；默认后端不变。

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

- [SLOT actual SRH/Auger branches and local chain checks, 2026-10-27](validation/slot_ldmos_srh_auger_branch_chain_2026-10-27.md): live source operands and partials, fixed-branch directional checks, high-precision replay and separate native-consumer limits; production acceptance unchanged.

- [SLOT native 80-bit QF/DOS chain checks, 2026-10-26](validation/slot_ldmos_native80_qf_dos_chain_2026-10-26.md): typed cache binding, exact eta replay, frozen-state directional checks and asymmetric weak recombination responses; full native Newton chain remains separate.

- [SLOT P1 acceptance checklist, 2026-10-21](validation/slot_ldmos_p1_acceptance_checklist_2026-10-21.md): fixed shadow-diagnostics scope, verified evidence, implementation and coverage gaps, closure conditions, and separate P2/native-alignment follow-ups; no change to acceptance rules.
- [SLOT P1 closeout batch, 2026-10-21](validation/slot_ldmos_p1_closeout_batch_2026-10-21.md): independent 24-hour observer build, update-chain and rejected-trial capture, paired validation, and explicitly limited historical substitution.

- [SLOT native low-bias convergence replay, 2026-10-21](validation/slot_ldmos_native_convergence_replay_2026-10-21.md): Sentaurus iteration capture, SLOT update-formula replay, P1 shadow limits, current-definition checks and an isolated IdVd doping-mismatch attribution; production criteria unchanged.

- [SLOT P1 shadow convergence diagnostics, 2026-10-20](validation/slot_ldmos_p1_shadow_convergence_2026-10-20.md): default-off observer, single-thread invariance pairing, exact archived-state replay and explicit P1 coverage limits; legacy acceptance unchanged.

- [SLOT weak-row classification, update metrics and bounded repeats, 2026-10-20](validation/slot_ldmos_weak_rows_repeat_adjoint_2026-10-20.md): saved-state edge classification, actual-update action, three U metrics, regional conservation, independent adjoint checks and an ongoing four-hour B repeat campaign.
- [SLOT scale binding, bounded span and Fermi-state audit, 2026-10-20](validation/slot_ldmos_scale_span_fermi_audit_2026-10-20.md): matched BVDS archive/build constants, diagnostic controls, same-parent span probes, IdVd population audit, and corrected weak-row physical-current units.

- [SLOT BVDS three-arm overnight attempt, 2026-10-19](validation/slot_ldmos_bvds_three_arm_overnight_2026-10-19.md): two fresh zero-bias starts with different high-voltage step caps and one qualified 170 mV restart; common twelve-hour deadline, full qualification and explicit resistor confirmation.

- [SLOT BVDS 170 mV tail and load-line confirmation, 2026-10-19](validation/slot_ldmos_bvds_170mv_tail_loadline_2026-10-19.md): one extra qualified update per path, full readbacks and production resistor seed confirmation; prior current path screen clarified as report-only for this case.

- [SLOT BVDS structure-preserving adjacent points, 2026-10-19](validation/slot_ldmos_bvds_poly_structure_adjacent_2026-10-19.md): full qualification and zero-update readback at 160/170 mV, unchanged homogeneous electron block, and a stopped current path-reproduction screen; no long scan.

- [SLOT BVDS homogeneous polysilicon structure-preserving solve, 2026-10-19](validation/slot_ldmos_bvds_poly_structure_preserving_2026-10-19.md): exact eligibility guards, original full-system defect checks, and paired bounded updates with repeated full qualification and archive readback; no new bias scan.

- [SLOT BVDS first-direction capture and thread replay, 2026-10-19](validation/slot_ldmos_bvds_first_direction_thread_replay_2026-10-19.md): actual PARDISO CSR/RHS and affine lift identity, one/four-thread fixed-system repeats, and homogeneous polysilicon direction drift; no new bias continuation.

- [SLOT BVDS production-cache identity and step-failure association, 2026-10-19](validation/slot_ldmos_bvds_cache_identity_step_association_2026-10-19.md): same-run population/SG generation evidence, source derivative identity checks, and an independent read-only audit of repeated small steps near 0.15576 V.

- [SLOT BVDS gap/tangent probes, 2026-10-19](validation/slot_ldmos_bvds_gap_tangent_preparation_2026-10-19.md): saved-operand symmetric decomposition and bounded VM finite perturbations, with explicit state/readback and branch-evidence limits; no Newton updates or new bias scan.

- [SLOT BVDS same-state operand replay, 2026-10-19](validation/slot_ldmos_bvds_same_state_operand_probe_2026-10-19.md): frozen 25/50 mV native SG/source captures, independent arithmetic replay, population/drive attribution and verified raw contact-reaction ledgers; no Newton updates.

- [SLOT BVDS qualified 50 mV twelve-hour continuation, 2026-10-19](validation/slot_ldmos_bvds_50mv_resume12_2026-10-19.md): independent qualified-parent continuation, staged adaptive device steps, explicit production resistor confirmation and fixed global deadline.

- [SLOT BVDS 5 mV cold-parent short path, 2026-10-19](validation/slot_ldmos_bvds_5mv_cold_path_2026-10-19.md): one-hour independent zero-bias requalification, adaptive device-voltage continuation with 25/50 mV checkpoints, unchanged full qualification and separate supply-voltage ledger.

- [SLOT BVDS overnight bounded directions, 2026-10-18](validation/slot_ldmos_bvds_overnight_directions_2026-10-18.md): four independent fixed-device, augmented-circuit, direction/coordinate and operand-arithmetic groups; fixed 12-hour window and unchanged full qualification.

- [SLOT BVDS charge-constant compatibility and local residuals, 2026-10-18](validation/slot_ldmos_bvds_charge_compat_local_residuals_2026-10-18.md): isolated q/epsilon profile with fixed Vt, 25/50 mV population/QF contracts, exact residual ledgers and remaining transport-sensitive weak rows.

- [SLOT BVDS 25/50 mV same-state matrix comparison, 2026-10-18](validation/slot_ldmos_bvds_25_50mv_matrix_comparison_2026-10-18.md): matched actual device biases, population/QF state contracts, physical Jacobian/RHS blocks and augmented contact equations.

- [SLOT BVDS cross-block direction and finite-step effects, 2026-10-18](validation/slot_ldmos_bvds_crossblock_direction_effect_2026-10-18.md): latest frozen-state population/QF contracts, source-branch attribution limits, fixed-RHS derivative controls and full-physics residual/port evaluations.

- [SLOT BVDS weak-floor report policy and short paths, 2026-10-18](validation/slot_ldmos_bvds_floor_report_short_path_2026-10-18.md): explicit user-authorized weak-row reporting policy; non-floor rows, closure and boundary protections remain enforced.

- [SLOT BVDS ordered contract and short validation, 2026-10-18](validation/slot_ldmos_bvds_ordered_contract_validation_2026-10-18.md): fixed evidence contract, one-direction archived trials, constant-material SG coordinates, contact-voltage low parts and original-qualification gates.

- [SLOT BVDS boundary projection and weak-row response, 2026-10-18](validation/slot_ldmos_bvds_boundary_projection_weak_response_2026-10-18.md): default-off synchronized-boundary diagnostic, strict paired-direction stop, original qualification, and frozen PolySilicon common-mode Jacobian row-sum audit.

- [SLOT BVDS 75 mV path and bounded nonzero-supply diagnostic, 2026-10-18](validation/slot_ldmos_bvds_75mv_path_bordered_2026-10-18.md): strict path stop, poly weak-row and actual low-state audit, opt-in bordered interface regressions and contact-voltage rounding limitation.

- [SLOT BVDS independent high-precision directions and fresh cold qualification, 2026-10-18](validation/slot_ldmos_bvds_high_precision_direction_cold_2026-10-18.md): 113/160-bit frozen-system references, closure-sum rounding audit and bounded new-boundary cold qualification.

- [SLOT BVDS pivot controls and split boundary targets, 2026-10-18](validation/slot_ldmos_bvds_pivot_split_boundary_2026-10-18.md): frozen-system pivot/scaling controls, target-row defect audit, split Dirichlet high/low targets, voltage-column regression and strict repeated-read qualification.

- [SLOT BVDS frozen-system refinement and boundary synchronization, 2026-10-18](validation/slot_ldmos_bvds_fixed_refinement_boundary_sync_2026-10-18.md): unsuccessful reused-factor corrections, old/new-bias qualification controls, first-read Dirichlet low-tail clearing and idempotent repeated read.

- [SLOT BVDS regularization strength paired updates, 2026-10-18](validation/slot_ldmos_bvds_regularization_strength_2026-10-18.md): frozen-parent 1e-8/1e-10/0 controls, strict qualification failures, tail linear defects and zero-update drain low-state changes.

- [SLOT BVDS equilibrium reference and weak-row regularization, 2026-10-17](validation/slot_ldmos_bvds_equilibrium_weak_rows_2026-10-17.md): same-potential common-QF controls, 220-digit SG replay, explicit contact-constraint ledger correction, and saved-direction explanation of residual retention on all 3,878 weak rows.

- [SLOT BVDS screening off/on relinearization, 2026-10-17](validation/slot_ldmos_bvds_ial_screening_relinearization_2026-10-17.md): same-parent three-update paired control, 39,150-dimensional device/circuit solves, unchanged strict qualification failures, exact low-state and boundary audits, and verified evidence package.

- [SLOT BVDS opt-in screening compatibility, 2026-10-17](validation/slot_ldmos_bvds_ial_screening_compatibility_2026-10-17.md): default-preserving compatibility selection, moving-cutoff thermal chain rule, one-sided branch validation, actual source derivatives and explicitly reproduced baseline test failure.

- [SLOT BVDS IALMob screening and derivative-cache identity, 2026-10-17](validation/slot_ldmos_bvds_ial_screening_derivative_identity_2026-10-17.md): directly observed screening cutoff/G, independent one-step reconstruction, scattering/HFS/source replay, and normalized own-node source derivative identities.

- [SLOT BVDS low-field/HFS and carrier-source decomposition, 2026-10-17](validation/slot_ldmos_bvds_hfs_source_detail_2026-10-17.md): native vertex low-field operands, same-contract HFS reconstruction, actual carrier-resolved source caches, and resolved branch-controlled source derivatives.

- [SLOT BVDS state contract, SG/HFS and avalanche pairing, 2026-10-17](validation/slot_ldmos_bvds_sg_hfs_state_contract_2026-10-17.md): same-system cache phases, actual SG and mobility operands, explicit AvalDensGradQF semantics, source scatter and branch-aware derivative diagnostics.

- [SLOT BVDS population generation and QF cache audit, 2026-10-17](validation/slot_ldmos_bvds_population_qf_cache_2026-10-17.md): same-population versus generation-input-QF controls, runtime band/DOS/ni operands, cache conversion stages, and finite-drive source value/derivative comparisons.

- [SLOT BVDS direct source and same-input SG arithmetic, 2026-10-17](validation/slot_ldmos_bvds_direct_source_same_input_2026-10-17.md): actual SRH/Auger numerator, denominator and scatter, independently matched support/lifetime, primary-variable equilibrium consistency, and 18 frozen PolySilicon edge/carrier arithmetic controls.

- [SLOT BVDS actual SG operand back-substitution, 2026-10-17](validation/slot_ldmos_bvds_sg_operand_backsubstitution_2026-10-17.md): runtime mobility operands, 23 active-edge stencils, actual-RHS row reconstruction, source/scatter remainder limits, and native default versus DirectCurrent contact-current definitions.

- [SLOT BVDS SG diagnostic export repair, 2026-10-17](validation/slot_ldmos_bvds_sg_diagnostic_export_fix_2026-10-17.md): carrier-specific production prefactors, independent 200-digit arithmetic reference, exact split-input sidecar, zero-support availability, and unchanged fixed-state RHS/current evidence.

- [SLOT BVDS split ledger and native SG operands, 2026-10-17](validation/slot_ldmos_bvds_split_ledger_native_sg_operands_2026-10-17.md): actual RHS/contact-reaction ledger, runtime SG operands and ordered accumulation, independent decimal arithmetic, and the observer/production population-prefactor boundary.

- [SLOT BVDS zero-bias RHS and bounded follow-up, 2026-10-17](validation/slot_ldmos_bvds_zero_rhs_followup_2026-10-17.md): constitutive/contact-reaction separation, six descending coupled updates, unchanged row/closure qualification failures, explicit boundary audit and separately passing equilibrium references.

- [SLOT BVDS zero-bias bordered iteration, 2026-10-17](validation/slot_ldmos_bvds_zero_bordered_iteration_2026-10-17.md): real 39,150-dimensional PARDISO device/circuit system, affine carrier constraints, staged Poisson/Jacobian comparison and an unqualified single-update diagnostic.

- [SLOT BVDS weak flux, finite tail and equilibrium reference, 2026-10-16](validation/slot_ldmos_bvds_weak_flux_tail_equilibrium_2026-10-16.md): independent electron/hole stable-flux checks, source balance, rejected-tail writeback audit and a separately qualified zero-bias equilibrium reference.

- [SLOT BVDS near-zero Auger relinearization and manual audit, 2026-10-16](validation/slot_ldmos_bvds_auger_relinearization_manual_2026-10-16.md): residual-preserving three-arm, three-update diagnostics, unchanged strict qualification, archived low-part statistics and verified native Auger formula/branch documentation.

- [SLOT BVDS Auger equilibrium arithmetic, 2026-10-16](validation/slot_ldmos_bvds_auger_equilibrium_arithmetic_2026-10-16.md): direct-library zero-point properties, 100-digit arithmetic references, native source/derivative transition brackets and explicit limits on inferring the native internal zero convention.

- [SLOT BVDS hole-source branches and cross-carrier support, 2026-10-16](validation/slot_ldmos_bvds_compatibility_hole_source_branches_2026-10-16.md): paired hole-population directions, distinct Auger and avalanche branches, actual element support and local absolute derivative scales.

- [SLOT BVDS thermal/ni compatibility, 2026-10-16](validation/slot_ldmos_bvds_thermal_ni_compatibility_2026-10-16.md): isolated version-bound four-arm experiment, fV-scale constitutive alignment, source/matrix comparison and explicit thermal-baseline regression limitation.

- [SLOT BVDS native constitutive and source finite differences, 2026-10-16](validation/slot_ldmos_bvds_native_constitutive_source_fd_2026-10-16.md): actual DOS/BGN outputs, independent band-reference mapping reconstruction, validated native state perturbations and source-only ablation.

- [SLOT BVDS constitutive and contact-chain audit, 2026-10-16](validation/slot_ldmos_bvds_constitutive_contact_chain_audit_2026-10-16.md): corrected eliminated-contact chain transformation, thermal-scale mapping trend, fixed-population finite differences and current VM task stops.

- [SLOT BVDS same-state assembly and actual PARDISO capture, 2026-10-16](validation/slot_ldmos_bvds_same_state_pardiso_capture_2026-10-16.md): matched native populations, QF-to-density Jacobian differences, native-QF control, and same-process CSR/RHS/returned-direction verification at the actual solver interface.

- [SLOT BVDS split-state 24-hour continuation and native single iteration, 2026-10-16](validation/slot_ldmos_bvds_C_resume24_native_iteration_2026-10-16.md): independent strict-policy hot continuation, preserved original deadline after a packed-predictor failure, and compensated-doping native binary64 input/direction/writeback capture.

- [SLOT BVDS C-arm tail and load-line validation, 2026-10-16](validation/slot_ldmos_bvds_C_tail_loadline_validation_2026-10-16.md): same-target diagnostic block corrections, original full qualification, actual external-resistor control and decoded-state confirmation.

- [SLOT BVDS C-arm split-state adjacent validation, 2026-10-16](validation/slot_ldmos_bvds_C_split_adjacent_2026-10-16.md): same-parent nearby targets and sampled path reproduction under strict qualification, exact residual reload, and checkpoint QF-frame restoration.

- [SLOT BVDS C-arm production low-state validation, 2026-10-16](validation/slot_ldmos_bvds_C_production_low_state_2026-10-16.md): opt-in IALMob/SRH/Auger/avalanche integration, strict qualification and exact reload, explicit low-zeroing control, and a same-target warm-start reproduction.

- [SLOT BVDS C-arm compensated Poisson validation, 2026-10-16](validation/slot_ldmos_bvds_C_compensated_poisson_2026-10-16.md): same-direction hi/lo versus ordinary coordinates, nonlinear Poisson reduction, loss on projection into ordinary DD, and local contact-guard evidence without changing qualification.

- [SLOT BVDS C-arm paired and frozen-state validation, 2026-10-16](validation/slot_ldmos_bvds_C_paired_poisson_contact_validation_2026-10-16.md): matched hot/secant failures, direct packed-writeback quantization evidence, and source-edge current/continuity sensitivity with unchanged qualification.

- [SLOT BVDS C-arm contact-guard stop analysis, 2026-10-16](validation/slot_ldmos_bvds_C_guard_failure_2026-10-16.md): eight bounded failures at 29.233 V, independent source-edge QF reconstruction, Poisson plateau, and unchanged qualification rules.

- [SLOT BVDS disk-stop resumption and 24-hour extension, 2026-10-16](validation/slot_ldmos_bvds_disk_resume_extension_2026-10-16.md): six independently supervised warm continuations, source-specific frozen solvers and acceptance rules, preserved pending spans, and original deadlines extended by 24 hours.

- [SLOT compensated BVDS hot-start 1 V ceiling, 2026-10-15](validation/slot_ldmos_compensated_bvds_hot1_cold_2026-10-15.md): independent 24-hour zero-bias cold start, same B strategy with only the span ceiling increased.

- [SLOT compensated BVDS parallel cold starts, 2026-10-15](validation/slot_ldmos_compensated_bvds_parallel_cold_2026-10-15.md): four independently supervised 24-hour cold-start arms comparing hot/secant initialization and step ceilings under unchanged strict qualification.

- [SLOT compensated BVDS adjacent validation, 2026-10-15](validation/slot_ldmos_compensated_bvds_adjacent_2026-10-15.md): strict cold-to-100 mV qualification with one bounded step reduction, four accepted points and unchanged physics.

- [SLOT compensated cold recovery, 2026-10-15](validation/slot_ldmos_compensated_cold_recovery_2026-10-15.md): bounded recovery obtains a strictly qualified zero state; paired continuation fails, and frozen probes quantify Poisson writeback loss without changing acceptance rules.
- [SLOT compensated cold failure analysis, 2026-10-15](validation/slot_ldmos_compensated_cold_failure_analysis_2026-10-15.md): coupled cold initialization, residual plateau, weak rows, gate contact edge and the unused recovery path; no restart or tolerance change.
- [SLOT compensated-doping BVDS cold campaign, 2026-10-15](validation/slot_ldmos_compensated_bvds_cold_campaign_2026-10-15.md): independent 24-hour VM cold-start run with original active compensation, current discretization, strict BVDS qualification and detached supervision.
- [SLOT original active-doping node audit, 2026-10-15](validation/slot_ldmos_active_doping_node_audit_2026-10-15.md): original As/P/B active components match the old compensated Vela input; the recent net-split input differs despite matching net doping.
- [SLOT BVDS total-doping cold initialization, 2026-10-15](validation/slot_ldmos_bvds_total_doping_cold_2026-10-15.md): independent native control and compensated-doping zero states, preserved input fields, BGN agreement and scoped mobility comparison to the old Vela model.
- [SLOT BVDS tdx doping comparison, 2026-10-15](validation/slot_ldmos_bvds_tdx_doping_comparison_2026-10-15.md): GRD/DAT conversion independently confirms input versus cold-state donor/acceptor differences, with source-neighborhood coordinates and conversion limitations.
- [SLOT BVDS September/October history comparison, 2026-10-15](validation/slot_ldmos_bvds_history_comparison_2026-10-15.md): verified historical high-voltage states, changed doping/discretization/reference frame, and conditional Poisson/contact stall evidence.
- [SLOT IdVd 1.10 V final analysis and plan, 2026-10-15](validation/slot_ldmos_idvd_1p10_analysis_plan_2026-10-15.md): completed research chain, early-stop versus local recovery, remaining material-specific errors and scoped next steps.
- [SLOT overnight and daily work summary, 2026-10-14 to 15](validation/slot_ldmos_work_summary_2026-10-14_night_to_15.md): source-backed daily-report material separating completed diagnostics, numerical limitations and ongoing tasks.
- [SLOT BVDS Poisson writeback audit, 2026-10-15](validation/slot_ldmos_bvds_poisson_writeback_audit_2026-10-15.md): frozen linear directions, sub-ULP packed updates, residual readback and independent contact guard evidence.
- [SLOT BVDS local double/binary128 solve, 2026-10-15](validation/slot_ldmos_bvds_precision_local_solve_2026-10-15.md): packed PolySilicon restore consistency repair; paired local solves retain Poisson/contact limitations despite carrier closure passing.
- [SLOT BVDS extended-precision prefactor repair, 2026-10-15](validation/slot_ldmos_bvds_precision_prefactor_2026-10-15.md): corrected material-reference populations, unchanged double baseline, gauge/ULP limitations and frozen contact/load-line sensitivity.
- [SLOT IdVd path reproduction and BVDS Poisson/contact guard diagnostics, 2026-10-15](validation/slot_ldmos_path_reproduction_poisson_guard_2026-10-15.md): 1.030 V path comparison, bounded continuation status, matched hot/secant failures, ULP sensitivity and contact-edge threshold evidence.
- [SLOT IdVd adjacent results and BVDS stall analysis, 2026-10-15](validation/slot_ldmos_idvd_adjacent_bvds_stall_analysis_2026-10-15.md): recovered dilute populations, remaining polysilicon rows, and the contact-edge QF veto at the BVDS Poisson residual plateau.
- [SLOT BVDS overnight analysis and 24-hour continuation, 2026-10-15](validation/slot_ldmos_bvds_overnight_analysis_resume_2026-10-15.md): budget versus numerical stops, Poisson/contact-QF stall evidence, and three independently resumed strict-policy arms.
- [SLOT coordinated-source initial-state and adjacent-point validation, 2026-10-15](validation/slot_ldmos_source_coordination_adjacent_2026-10-15.md): bounded research campaign, previous-job status and contact-consistent warm-start checks; see report for completion status.
- [SLOT transport/source linearization coordination, 2026-10-15](validation/slot_ldmos_source_transport_coordination_2026-10-15.md): same-seed four-arm diagnostics, released-row balance and unchanged residual checks.
- [SLOT dilute mapping and temporary density hold, 2026-10-15](validation/slot_ldmos_inverse_hold_parallel_validation_2026-10-15.md):
  extreme-density arithmetic repair, four bounded same-seed diagnostic arms,
  and original versus actual update-direction/source-derivative attribution.

- [SLOT BVDS aligned cold start and three matched targets, 2026-10-14](validation/slot_ldmos_bvds_aligned_cold_three_points_2026-10-14.md):
  independent AverageBox and effective-doping configuration, strict fresh-cold
  qualification, and bounded nonzero-bias recovery with four-thread PARDISO.

- [SLOT BVDS zero-bias same-state comparison, 2026-10-14](validation/slot_ldmos_bvds_zero_state_comparison_2026-10-14.md):
  matched populations expose compensated-doping/BGN and contact differences;
  physical matrix comparison separates native drain auxiliary equations.

- [SLOT PhuMob hotspot and secant validation, 2026-10-14](validation/slot_ldmos_phumob_hotspot_secant_2026-10-14.md):
  exact invariant caching removes repeated G-minimum searches; fixed-state matrix,
  direction and field checks precede bounded two-point drain continuation.

- [SLOT overnight cost audit and native references, 2026-10-14](validation/slot_ldmos_overnight_cost_native_reference_2026-10-14.md):
  first-update residual guards limit IdVd spans; repeated bounded recovery suppresses
  BVDS growth; independent native reference jobs export matrices and physical fields.

- [SLOT BVDS guarded linear refinement, 2026-10-13](validation/slot_ldmos_bvds_guarded_refinement_2026-10-13.md):
  opt-in ordinary Newton refinement preserves old BVDS physics and qualification;
  fixed-matrix correction passes, while the paired failed target remains limited by nonlinear closure.

- [SLOT VM refined adjacent points, 2026-10-13](validation/slot_ldmos_vm_refined_adjacent_validation_2026-10-13.md):
  three first-attempt diagnostic continuations to 42 mV with retained-factor refinement,
  unchanged linear guards and exact saved-state residual/current readback.

- [SLOT VM refinement and shared diagnostics, 2026-10-13](validation/slot_ldmos_vm_refinement_shared_diagnostics_2026-10-13.md):
  unchanged linear guards pass with retained-factor refinement; bounded 6 to
  13.008 mV validation and bitwise paired shared-evaluation checks, with measured costs.

- [SLOT VM linear guard and diagnostic cost audit, 2026-10-13](validation/slot_ldmos_vm_linear_guard_cost_audit_2026-10-13.md):
  six rejected steps share one interface hole row; frozen-system residual
  correction passes unchanged linear guards, with dual-system timing limits.

- [SLOT Linux VM campaign and PARDISO threading, 2026-10-13](validation/slot_ldmos_linux_vm_long_campaign_2026-10-13.md):
  four-thread same-parent replay passed; detached twelve-hour IdVd and independent
  old-policy BVDS campaign launched, with remote deadline supervision.

- [SLOT Linux VM short validation, 2026-10-13](validation/slot_ldmos_linux_vm_short_validation_2026-10-13.md):
  GCC16 and Linux PARDISO verified at three adjacent drain biases with the
  diagnostic transport linearization; final readback completed in a separate
  zero-update audit. Sentaurus execution awaits its license configuration.
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

- [SLOT native derivative-switch audit, 2026-10-13](validation/slot_ldmos_native_derivative_switch_2026-10-13.md):
  same-state Sentaurus mobility/avalanche derivative controls preserve the solve RHS;
  111 negative hole diagonals and omitted terms match the Vela diagnostic closely,
  while wider source coupling and historical scan context remain distinct.

- [SLOT mobility-derivative adjacent paths, 2026-10-13](validation/slot_ldmos_mobility_adjacent_path_2026-10-13.md):
  bounded local/T470p continuation from 6 mV, different-path field and contact
  reaction checks, and field changes after the first RHS-or-update stopping point;
  diagnostic linearization only, with original residuals and guards retained.

- [SLOT projection and mobility-derivative recovery, 2026-10-13](validation/slot_ldmos_projection_mobility_recovery_2026-10-13.md):
  local/T470p bounded diagnostics distinguish repeated positivity projection,
  an unsuccessful activity-set release rule, and transport derivative feedback;
  original residuals and acceptance thresholds remain unchanged.

- [SLOT retained-potential three-step trial, 2026-10-13](validation/slot_ldmos_keep_psi_three_steps_2026-10-13.md):
  retain 126 contact-potential constraints in the 39,149-dimensional PARDISO
  system; identical initial assembly/RHS, three accepted adaptive drain points,
  and remaining low-density hole-QF path sensitivity.

- [SLOT IdVd three-day handoff, 2026-10-11 to 13](validation/slot_ldmos_idvd_handoff_2026-10-11_to_13.md):
  evidence-backed configuration and acceptance history, linear replay versus
  nonlinear field differences, current adaptive-step failures, and bounded
  diagnostic questions for an independent model review.

- [SLOT native-style steps toward 1 V, 2026-10-13](validation/slot_ldmos_native_steps_to_1v_2026-10-13.md):
  bounded local 6 mV start, iteration-aware growth and halving with Fermi QF
  extrapolation; six accepted points to 13.002 mV, numerical guards and budget stops.

- [SLOT RHS-or-update stopping policy, 2026-10-13](validation/slot_ldmos_rhs_or_update_policy_2026-10-13.md):
  current SLOT opt-in entry reports closure/boundary thresholds, retains numerical
  and saved-state guards, and verifies RHS-only and update-only acceptance.

- [SLOT single-process population loop, 2026-10-13](validation/slot_ldmos_population_memory_2026-10-13.md):
  opt-in in-memory updates and linear guards; fixed-parent 5/12.1 mV replay,
  unchanged fields and acceptance, bounded timing and output comparison.

- [SLOT acceptance and diagnostic efficiency, 2026-10-12](validation/slot_ldmos_population_acceptance_efficiency_2026-10-12.md):
  update-or-RHS acceptance with unchanged safeguards; paired frozen-state checks
  of compact matrix output and removal of repeated port/residual diagnostics.

- [SLOT population projection pair, 2026-10-12](validation/slot_ldmos_population_projection_pair_2026-10-12.md):
  fixed PARDISO four-update comparison of 0.001/0.01 positivity projection;
  0.001 qualifies, while 0.01 improves some hole populations but retains a Poisson residual gap.

- [SLOT same-direction population update discrimination, 2026-10-12](validation/slot_ldmos_same_direction_update_discrimination_2026-10-12.md):
  native and current rules diverge at first-step hole positivity projection;
  a single-factor 0.001/0.01 replay preserves the other fields but does not restore closure.

- [SLOT native system through production linear interfaces, 2026-10-12](validation/slot_ldmos_native_linear_interface_replay_2026-10-12.md):
  PARDISO and UMFPACK reproduce native float64 directions; single-thread results
  match direct library calls bitwise, with independent residual and input audits.

- [SLOT 5 mV refinement and contact-hole ledger, 2026-10-12](validation/slot_ldmos_5mv_linear_refinement_hole_ledger_2026-10-12.md):
  explicit residual correction repairs the frozen UMFPACK linear error;
  contact source accounting revises the source-hole comparison from 26.1% to 5.16 ppm.
  Low-density population projection and field differences remain unresolved.

- [SLOT exact 5 mV revalidation versus native float64, 2026-10-12](validation/slot_ldmos_5mv_native64_revalidation_2026-10-12.md):
  PARDISO passes scoped acceptance in four updates; UMFPACK is rejected by the
  third linear guard. Drain current differs by 7.59 ppm, with remaining hole-field differences.

- [SLOT 6 mV final states versus Sentaurus, 2026-10-12](validation/slot_ldmos_6mv_native_final_comparison_2026-10-12.md):
  exact-bias native comparison finds 7.59 ppm drain-current difference but
  substantial low-density hole-QF and source hole-current differences.

- [SLOT near-zero linear guard and backend pair, 2026-10-12](validation/slot_ldmos_mixed_linear_guard_pair_2026-10-12.md):
  opt-in mixed residual guard retains block and constraint checks; same-parent
  6 mV UMFPACK/PARDISO relinearization, with unchanged nonlinear acceptance.

- [SLOT T470p population update and strict elimination, 2026-10-12](validation/slot_ldmos_t470p_population_elimination_2026-10-12.md):
  paired UMFPACK/PARDISO pass scoped 10 mV acceptance; adjacent steps expose
  a finite-population inverse-mapping guard and a separate linear-error outlier.

- [SLOT population-update linear backend comparison, 2026-10-11](validation/slot_ldmos_population_linear_backends_2026-10-11.md):
  SparseLU and UMFPACK repeat four stabilization combinations against PARDISO;
  scoped acceptance and current agree while low-density hole-QF trajectories differ.

- [SLOT regularization and closure row replacement comparison, 2026-10-11](validation/slot_ldmos_stabilization_ablation_2026-10-11.md):
  four-arm population finite-step comparison with independent stabilization
  controls, unchanged physical equations, and preserved closure acceptance.

- [SLOT population finite-step verification, 2026-10-11](validation/slot_ldmos_population_finite_step_verification_2026-10-11.md):
  bounded relinearization reaches the scoped acceptance at 10 mV; a separate
  trajectory reproduces drain current but retains low-density hole-QF differences.

- [SLOT component population update and frozen Charon-inspired kernels, 2026-10-11](validation/slot_ldmos_component_population_charon_diagnostics_2026-10-11.md):
  three-step same-parent update comparison, stagnation and backward-error audit,
  and separate low-density model sensitivities; no new qualified continuation point.

- [SLOT local Charon and DEVSIM Newton review, 2026-10-11](validation/slot_ldmos_charon_devsim_newton_review_2026-10-11.md):
  source-verified population positivity, potential damping, scaling and stagnation
  checks; separates numerical updates from low-density driving-field model changes.

- [SLOT Bank-Rose population-domain recovery, 2026-10-11](validation/slot_ldmos_bank_rose_population_recovery_2026-10-11.md):
  bounded rejection of typed invalid-population trials preserves the parent and
  physics; the paired uncapped continuation audit records qualification separately.

- [SLOT uncapped Bank-Rose comparison, 2026-10-11](validation/slot_ldmos_row_report_bank_rose_2026-10-11.md):
  replacing backtracking with Bank-Rose leaves the population-domain exception unresolved;
  a larger initial K permits six tiny updates but does not establish global closure.

- [SLOT row-report continuation without QF caps, 2026-10-11](validation/slot_ldmos_row_report_nocap_2026-10-11.md):
  10 mV and the 7.5 mV retry exit during the first trial population evaluation;
  the initial 10 mV systems match the capped baseline, with no new accepted updates.

- [SLOT row-report neighboring voltages and path reproduction, 2026-10-11](validation/slot_ldmos_row_report_neighbors_2026-10-11.md):
  three new-policy points reach 20 mV with closure enforced; a different voltage path
  reproduces current and potential, while low-density carrier QF fields remain path dependent.

- [SLOT row-report stopping policy and 5 mV verification, 2026-10-11](validation/slot_ldmos_row_report_policy_2026-10-11.md):
  the user-selected case-local policy stops at update nine with original closure enforced;
  strict row diagnostics and local minority-field differences remain separately reported.

- [SLOT PARDISO/SparseLU three-update comparison, 2026-10-11](validation/slot_ldmos_mkl_relinearization_2026-10-11.md):
  an opt-in MKL backend reduces some weak-row violations while unchanged QF updates
  leave the dominant hole imbalance unqualified; initial matrices and RHS are identical.

- [SLOT oneMKL PARDISO native-system replay, 2026-10-11](validation/slot_ldmos_mkl_native_replay_2026-10-11.md):
  raw float64 native systems reproduce the saved direction closely; common population
  finite steps separate linear accuracy from remaining nonlinear qualification.

- [SLOT native-style OR stopping diagnostic, 2026-10-11](validation/slot_ldmos_native_gate_trial_2026-10-11.md):
  both eliminated backends stop at update nine on the relaxed RHS branch;
  unchanged trajectories and original zero-update audits distinguish diagnostic success from row qualification.

- [SLOT native/Vela convergence criteria comparison, 2026-10-11](validation/slot_ldmos_convergence_criteria_comparison_2026-10-11.md):
  native update-error stopping and shadow population-update metrics distinguish
  global convergence indicators from unresolved Vela carrier-row qualification.

- [SLOT eliminated-system QF cap on/off comparison, 2026-10-11](validation/slot_ldmos_eliminated_cap_toggle_2026-10-11.md):
  both uncapped precisions exit on a first-trial SG population-domain exception;
  no accepted updates, distinguishing failed trial evaluation from exhausted backtracking.

- [SLOT eliminated-system 5 mV Newton trials, 2026-10-11](validation/slot_ldmos_eliminated_newton_trial_2026-10-11.md):
  opt-in double/long-double elimination backends pass solver tests; bounded
  nonlinear trials satisfy closure but both fail original carrier-row qualification.

- [SLOT fixed-system constraint elimination comparison, 2026-10-11](validation/slot_ldmos_constraint_elimination_precision_2026-10-11.md):
  exact affine elimination agrees with the high-precision reference; double fixes
  the 10280 sign but increases other direction errors, so elimination alone is not a recovery.

- [SLOT high-precision fixed-system direction reference, 2026-10-11](validation/slot_ldmos_high_precision_direction_reference_2026-10-11.md):
  113/160-bit references and independent decimal residual checks quantify weak
  direction errors; double has the wrong sign at 10280, without implying nonlinear recovery.

- [SLOT in-memory matrix and direction capture, 2026-10-11](validation/slot_ldmos_in_memory_direction_capture_2026-10-11.md):
  unchanged production trajectories and exact matrix replay locate the 10280 sign
  split in the finite-precision base solve; both nonlinear paths remain unqualified.

- [SLOT same-state precision direction layers, 2026-10-11](validation/slot_ldmos_float80_direction_layers_2026-10-11.md):
  regularization, Woodbury, QF caps and Poisson recorrection isolated; three-neighbor
  swaps reproduce the 10280 row difference, while original qualification still fails.

- [SLOT mixed-precision three-update comparison, 2026-10-11](validation/slot_ldmos_float80_relinearization_2026-10-11.md):
  identical-parent production updates reduce some carrier-row violations with
  long double, while dominant residuals and closure remain nearly unchanged.

- [SLOT float80 trial and optional mixed-precision backend, 2026-10-11](validation/slot_ldmos_float80_trial_2026-10-11.md):
  checked long-double precision, global type migration blockers, identical-matrix
  residual comparison, and zero-update SLOT probes; production defaults unchanged.

- [SLOT Bank–Rose implementation and bounded comparison, 2026-10-10](validation/slot_ldmos_bank_rose_2026-10-10.md):
  opt-in Algorithm Global damping, persistent K and NotDamped activation;
  numerical tests and a same-parent 5 mV comparison retain original qualification.

- [SLOT undamped-window and OR stopping comparison, 2026-10-10](validation/slot_ldmos_native_style_stopping_2026-10-10.md):
  the diagnostic stops on RHS after nine updates but fails original hole-row qualification;
  distinguishes raw-direction update error, applied updates and unchanged production gates.

- [SLOT native 64-bit and 80-bit precision comparison, 2026-10-10](validation/slot_ldmos_native_precision_2026-10-10.md):
  matched 5 mV native solves both stop after four updates; current agrees closely,
  while residual and near-zero quantities retain precision sensitivity.

- [SLOT linear backends and Newton algorithm review, 2026-10-10](validation/slot_ldmos_linear_backends_newton_review_2026-10-10.md):
  six backend controls share identical matrices and RHS, with finite-step checks;
  version-bound Sentaurus and PETSc references separate linear accuracy from globalization.

- [SLOT boundary elimination, stabilization and finite steps, 2026-10-10](validation/slot_ldmos_boundary_stabilization_finite_steps_2026-10-10.md):
  compares both matrix origins in a common free QF space, audits stable closure-row
  solves, and separates finite population accuracy from actual row and conservation quality.

- [SLOT cross-block flux and source decomposition, 2026-10-10](validation/slot_ldmos_crossblock_flux_source_decomposition_2026-10-10.md):
  corrects omitted contact-column coordinate terms and separates smooth SG transport
  from near-zero avalanche-current directional sensitivity using matched-state probes.

- [SLOT zero and exact 5 mV initial matrix comparison, 2026-10-10](validation/slot_ldmos_initial_zero_5mv_matrix_comparison_2026-10-10.md):
  original matrix statistics include a subsequently corrected contact-column mapping
  error; see the erratum and cross-block decomposition above. RHS statistics are unaffected.

- [SLOT matched-system population update comparison, 2026-10-10](validation/slot_ldmos_pure_population_update_comparison_2026-10-10.md):
  three-step controls reuse the original stabilized direction, reduce hole-row failures
  without full qualification, and distinguish population updates from Poisson recorrection.

- [SLOT contact-preserving relative population inverse, 2026-10-10](validation/slot_ldmos_contact_preserving_relative_qf_2026-10-10.md):
  exact zero-change and contact preservation pass, with manufactured increment recovery;
  native population matching still fails row qualification and exposes sub-ULP QF sensitivity.

- [SLOT population inverse mapping and flux representation, 2026-10-10](validation/slot_ldmos_population_inverse_flux_representation_2026-10-10.md):
  frozen controls distinguish exact QF re-anchoring, density round-trip perturbations,
  boundary remapping and stable SG arithmetic; all qualification gates remain unchanged.

- [SLOT regional closure, relinearization and constitutive profile, 2026-10-10](validation/slot_ldmos_four_stage_relinearization_2026-10-10.md):
  three-step diagnostic replays separate transient electron closure from persistent row failures;
  bounded population directions and an independently built constants profile remain unqualified,
  with a frozen population-matching control exposing charge/transport sensitivity.

- [SLOT constitutive constants and population update coordination, 2026-10-10](validation/slot_ldmos_mapping_population_coordination_2026-10-10.md):
  native-final-state formulas explain remaining mapping and Poisson constant terms;
  bounded population and cross-block regularization trials improve some rows but
  fail full qualification, with explicit linear-backend sensitivity controls.

- [SLOT native final-state qualification and direction layers, 2026-10-10](validation/slot_ldmos_native_qualification_direction_layers_2026-10-10.md):
  zero-update qualification distinguishes residual mapping from stopping criteria;
  identical physical matrices expose regularization and clipping effects, with
  an independently reconstructed first production update and unchanged gates.

- [SLOT effective-band transport implementation and bounded validation, 2026-10-10](validation/slot_ldmos_effective_band_transport_validation_2026-10-10.md):
  opt-in Fermi edge transport, complete small-mesh derivative tests and SLOT column checks;
  frozen Si electron residual agreement improves, but both nonlinear controls stop
  on carrier-row line-search rejection with no qualified new state.

- [SLOT Fermi SG edge and constitutive mapping, 2026-10-10](validation/slot_ldmos_fermi_sg_edge_constitutive_mapping_2026-10-10.md):
  frozen edge substitutions isolate a Fermi transport-weight difference; multi-node
  thermal-voltage and base-ni checks explain spatial QF offsets without changing production defaults.

- [SLOT full residual, SG transport and contact mapping, 2026-10-10](validation/slot_ldmos_full_residual_transport_mapping_2026-10-10.md):
  frozen population/QF controls separate constitutive driving-field differences from
  transport assembly; contact sets and source/transport units are audited without new accepted states.

- [SLOT source-kernel comparison with joint stopping rules, 2026-10-10](validation/slot_ldmos_source_trend_stop_revalidation_2026-10-10.md):
  opt-in row-quality/raw-direction/line-search monitoring and longer bounded controls;
  both arms stop on original line-search failure before time limits, with no new qualified point.

- [SLOT SG source-kernel alignment, 2026-10-10](validation/slot_ldmos_sg_source_alignment_2026-10-10.md):
  opt-in effective-band SG current reconstruction with partial-layer current/force
  projection, live population derivatives, source-support checks and bounded
  self-consistent controls. Frozen source agreement is separate from device qualification.

- [SLOT native source support and residual audit, 2026-10-09](validation/slot_ldmos_native_source_residual_audit_2026-10-09.md):
  same-state source splitting verifies runtime source volumes against actual RHS;
  SG current reconstruction and interface projection differ from Vela's density-gradient source.
  Fresh native assembly passes two fixed-population psi finite differences; earlier iteration
  matrices are not interchangeable with the independently reassembled physical derivative.

- [SLOT PMI source and cross-block comparison, 2026-10-09](validation/slot_ldmos_pmi_source_crossblock_comparison_2026-10-09.md):
  matched-population diagnostics reduce QF mapping differences to sub-microvolt
  levels, but hotspot source and psi cross-block differences remain. Two Vela
  derivative columns pass finite-difference checks; native source support remains unresolved.

- [SLOT native integral BGN PMI comparison, 2026-10-09](validation/slot_ldmos_native_integral_bgn_pmi_2026-10-09.md):
  GCC 11 PMI passes scalar/derivative tests and two zero-update native observations.
  Explicit BW plus normalized-integral statistics matches all exported Si/Poly
  nodes; no new self-consistent device state or full Jacobian equivalence is claimed.

- [SLOT gate continuation and PolySilicon BGN audit, 2026-10-09](validation/slot_ldmos_gate_continuation_poly_bgn_2026-10-09.md):
  four new gate points and independent reproduction pass to 0.160859375 V.
  Native normalized quadrature aligns PolySilicon population relations; a
  NoFermi diagnostic isolates the remaining Si statistical BGN correction.

- [SLOT aligned gate path and native integral route, 2026-10-09](validation/slot_ldmos_aligned_gate_integral_route_2026-10-09.md):
  three small gate points and a separate direct cold-to-target reproduction pass
  the original gates. Native RelativeToDOSMass quadrature aligns the Si population
  integral, but BGN statistical correction and PolySilicon remain unmatched.

- [SLOT AverageBox cold state and constitutive alignment, 2026-10-09](validation/slot_ldmos_average_box_cold_constitutive_alignment_2026-10-09.md):
  a separate configuration uses the actual AverageBox Poisson path and passes
  fresh cold-state qualification. The 8.214 ppm ratio is explained by physical
  constants; native DOS matches, while Fermi approximation and BGN mapping
  differences remain. No gate/drain continuation or stabilization ablation.

- [SLOT constitutive, Poisson and source audit, 2026-10-09](validation/slot_ldmos_constitutive_poisson_source_audit_2026-10-09.md):
  the signed per-cell Poisson override bypasses AverageBox; a frozen control
  aligns potential stencils, while charge scaling and Fermi/QF mapping still
  differ. Two hotspot cross-block derivatives pass local finite differences.
  Native intermediate fields are usable despite final-output timeout; no Vela
  continuation point is accepted.

- [SLOT same-state residual and Jacobian layers, 2026-10-09](validation/slot_ldmos_same_state_matrix_layers_2026-10-09.md):
  native free populations are matched before a full Fermi chain transformation.
  Contact constitutive mismatches and hotspot cross-block differences remain;
  weak-row regularization is quantified separately. No recovery step is accepted.

- [SLOT connected population trials with neighbor guards, 2026-10-09](validation/slot_ldmos_connected_population_guard_2026-10-09.md):
  twelve same-state frozen trials compare restricted and locally recomputed
  population directions. No trial passes true neighbor guards; new electron
  failures are already visible in linear predictions, motivating explicit
  protected-row constraints rather than soft least-squares weights alone.

- [SLOT population updates with Fermi inverse mapping, 2026-10-09](validation/slot_ldmos_population_fermi_update_2026-10-09.md):
  same-state bounded population trials improve hotspot recovery and Poisson
  residuals but create thousands of new hole-row failures. Mapping and frozen
  population checks pass; no trial is accepted for continuation.

- [SLOT second-step population recovery and QF clipping, 2026-10-09](validation/slot_ldmos_second_step_recovery_2026-10-09.md):
  exact first-step replay and a bounded second update distinguish a positive
  recovery direction clipped at node 10849 from continued depletion and
  cross-block coupling at 10237. Full qualification remains unmet; comparisons
  across different states and variable coordinates do not establish a Jacobian bug.

- [SLOT native Newton matrices and iteration comparison, 2026-10-09](validation/slot_ldmos_native_newton_matrix_2026-10-09.md):
  Sentaurus `-d` exports five full sparse systems; NewtonPlot captures a
  four-update 6 mV recovery with density positivity handling. Matrix/vector
  reconstruction and observer controls pass; cross-solver variable transforms
  and Vela final qualification remain distinct, unfinished checks.

- [SLOT bounded joint least-squares directions, 2026-10-09](validation/slot_ldmos_bounded_joint_direction_2026-10-09.md):
  equilibrated constrained subproblems yield millivolt-scale coupled updates and
  reduce all three residual blocks. Original failed rows and closure remain
  unqualified: hole source reduction has not established balancing transport.
  No recovery step or continuation state is approved.

- [SLOT regularization and cross-block direction comparison, 2026-10-09](validation/slot_ldmos_regularization_block_direction_2026-10-09.md):
  four regularization strengths and conditional hole-block directions do not
  provide an admissible 6 mV first step. Removing regularization yields enormous
  directions; clipping and block damping retain thousands of failed rows.
  All trials are diagnostic, with unchanged physics and final qualification.

- [SLOT 6 mV direction, clipping and merit decomposition, 2026-10-09](validation/slot_ldmos_first_step_direction_merit_2026-10-09.md):
  zero-update replay reproduces the first accepted trial and separates carrier
  clipping from Poisson recorrection. Global merit masks hole-row deterioration;
  source growth and regularized cross-block directions require coordinated repair.
  No production physics, final threshold or qualified-point count changes.

- [SLOT zero-state population/source and 6 mV first-step audit, 2026-10-09](validation/slot_ldmos_zero_state_source_first_step_2026-10-09.md):
  frozen native fields show up to 1.35% population remapping, while density-only
  substitution changes avalanche integration by 0.278%. Two bounded first-step
  arms retain thousands of failed rows; global merit decreases while hole
  residuals increase. No qualification or production physics change is implied.

- [SLOT zero-drain native field comparison, 2026-10-09](validation/slot_ldmos_zero_drain_native_comparison_2026-10-09.md):
  exact Vg=5 V / Vd=0 native TDR and Vela state share node coordinates,
  connectivity and materials. Channel-surface potential differs by at most
  1.8 µV; localized silicon differences reach 0.346 mV and 2.71% in density.
  Near-zero currents and source outputs require separate interpretation;
  this comparison does not resolve the nonzero-drain bootstrap failure.

- [SLOT 6 mV drain bootstrap short validation, 2026-10-09](validation/slot_ldmos_drain_6mv_short_validation_2026-10-09.md):
  the latest 5 V gate / zero-drain state revalidates, but 6, 3 and 1.5 mV
  warm starts all reach their bounded soft deadlines with thousands of failed
  carrier rows, mainly holes. No new drain state is qualified; neighboring
  extrapolation and reproduction remain conditional on establishing a first point.

- [SLOT gate warning-only continuation, 2026-10-09](validation/slot_ldmos_gate_warning_resume_2026-10-09.md):
  the predictor now logs changes above 0.1 V without rejecting them; the adaptive
  policy retains its 0.25 V maximum and halves failed spans. The separately
  authorized three-hour gate-only continuation starts from a revalidated
  0.496289 V state. See the report and live campaign evidence for completion.

- [SLOT gate 0.25 V predictor comparison, 2026-10-09](validation/slot_ldmos_gate_025_comparison_2026-10-09.md):
  one 0.25 V step and a 0.1/0.1/0.05 V path both pass the original full gates
  from Vg=0.246289 V to 0.496289 V at Vd=0. Maximum potential difference is
  8.22e-15 V; elapsed costs are 157 versus 391 seconds. This local diagnostic
  did not itself change the then-production 0.1 V predictor cap or establish
  full-ramp stability; the subsequent warning-only change is documented above.

- [SLOT IdVd three-hour continuation and diagnosis, 2026-10-08](validation/slot_ldmos_idvd_three_hour_resume_2026-10-08.md):
  original gates passed at 2.4 V and reproduced through a different fine voltage
  path; SG current agreement is 4.1e-10 relative. Frozen observations distinguish
  failed electron-row inflow from source terms. The budgeted run stopped at
  3.339357 V with 81 new main-path qualified points; three reproduction points
  are counted separately. Final record and hash checks pass; no 60 V claim.

- [SLOT IdVd small-step secant path, 2026-10-08](validation/slot_ldmos_idvd_small_step_path_2026-10-08.md):
  four new qualified points through 1.194279 V without fallback; endpoint
  reproduction from an earlier qualified seed pair passes the original gates.
  Total solver wall time is 816 s; original physics and thresholds retained.

- [SLOT IdVd electron weak rows and bounded continuation, 2026-10-08](validation/slot_ldmos_idvd_electron_weak_row_2026-10-08.md):
  two new qualified drain voltages through 0.921167 V with same-target reproduction;
  frozen HFS response dominates one old failed row, while connected restoration
  is rejected for neighbor damage. Original physics and qualification retained.

- [SLOT IdVd AverageBox implementation and bounded recomputation, 2026-10-07](validation/slot_ldmos_idvd_average_box_2026-10-07.md):
  opt-in unweighted 2-D face intersections reproduce native zero edges and
  element measures; 49 C++ cases and 19 Python checks pass. Three bounded
  controls remain unqualified with carrier-row and closure failures.
- [SLOT IdVd native box geometry and update decomposition, 2026-10-07](validation/slot_ldmos_idvd_update_decomposition_geometry_2026-10-07.md):
  native coefficients and measures observed with zero updates; first-step row
  expansion reproduced; factor-one block filter is not a nonworsening guard.
  Two bounded arms remain unqualified with original final gates preserved.
- [SLOT IdVd and SimpleMOS alignment controls, 2026-10-07](validation/slot_ldmos_idvd_simplemos_alignment_2026-10-07.md):
  non-Delaunay edges block direct geometry transfer; three native ablations
  complete, while six bounded local physics/numerics arms remain unqualified.
  Frozen edge/source/port checks preserve original gates and ledger identity.
- [SLOT IdVd existing failed-row ratio guards, 2026-10-07](validation/slot_ldmos_idvd_failed_row_ratio_guard_2026-10-07.md):
  real guards detect nine worsening failed rows; five constrained inner solves
  remain unverified despite first-order descent directions. No accepted update
  or qualified point; original physics and final gates unchanged.
- [SLOT IdVd nonlinear SG neighbor guards, 2026-10-07](validation/slot_ldmos_idvd_nonlinear_sg_guard_2026-10-07.md):
  protects 4034/4036 in same-parent controls and accepts two diagnostic updates;
  still 147 failed rows and worsening 4046 ratio from faster scale reduction.
  Original gates retained; no qualified point or scan restart.
- [SLOT IdVd SG/mobility curvature split, 2026-10-07](validation/slot_ldmos_idvd_sg_mobility_curvature_2026-10-07.md):
  fixed-mobility SG factors explain 98.77%/98.15% of local flux curvature at
  4034/4036; mobility and interaction terms separately reconciled to full flux.
  Frozen probes only; no accepted update or qualified point.
- [SLOT IdVd fixed-direction derivative check, 2026-10-07](validation/slot_ldmos_idvd_fixed_direction_derivative_2026-10-07.md):
  central differences on 4034/4036 converge quadratically to Jd within
  2.57/3.03 ppm; flux curvature explains the finite-step prediction failure.
  Frozen diagnostics only, no accepted Newton update or qualified point.
- [SLOT IdVd edge subdomains and prediction margins, 2026-10-07](validation/slot_ldmos_idvd_edge_patch_margin_2026-10-07.md):
  three connected 63-DOF patches and reoptimized smaller steps remain rejected;
  isolates nonlinear prediction errors at 4034/4036 and changing row scale at
  4046. No accepted update, qualified point or production restart.
- [SLOT IdVd carrier and conservation guards, 2026-10-07](validation/slot_ldmos_idvd_carrier_conservation_guard_2026-10-07.md):
  block-balanced directions and real closure checks accept three diagnostic
  updates across two chains; scaled trial directions with nonlinear damage
  are rejected. Still 147 failed rows and no qualified point.
- [SLOT IdVd guarded coupled carrier/potential repair, 2026-10-07](validation/slot_ldmos_idvd_guarded_coupled_repair_2026-10-07.md):
  matched 63-DOF controls prevent new failed rows; four diagnostic updates
  reduce Poisson residual 42-fold but leave 147 failed rows and worsen hole
  closure. No qualified point or production restart.
- [SLOT IdVd density mapping and hotspot initialization, 2026-10-07](validation/slot_ldmos_idvd_density_mapping_hotspot_seed_2026-10-07.md):
  independent Fermi quadrature localizes the mapping discrepancy; 55-node
  population-consistent seeds increase transport but disturb Poisson and
  neighbor rows. No common fixed-state loading scale, no qualified point.
- [SLOT IdVd contact reaction, population representation and native HFS audit, 2026-10-07](validation/slot_ldmos_idvd_contact_reaction_population_hfs_2026-10-07.md):
  zero-update native observer restored and identity verified; contact source
  reactions and SG fluxes separated, saved/reconstructed populations matched.
  No Newton update or qualified point; initialization remains conditional.
- [SLOT IdVd source ledger and HFS formation audit, 2026-10-07](validation/slot_ldmos_idvd_source_ledger_hfs_audit_2026-10-07.md):
  independent edge/source conservation closes internally; native signed-Si
  source integration exposes a contact-accounting question. Native zero-update
  observer times out at startup; conditional loading remains unstarted.
- [SLOT IdVd coupled restoration, source homotopy and native RefDens controls, 2026-10-07](validation/slot_ldmos_idvd_restoration_homotopy_refdens_2026-10-07.md):
  connected 63-variable restoration remains unqualified; verified homotopy
  directions reveal a hotspot population barrier. Two native smoothing arms
  complete and the explicit-zero arm times out; no production model change.
- [SLOT IdVd low-density mechanism and native-control review, 2026-10-07](validation/slot_ldmos_idvd_low_density_mechanism_review_2026-10-07.md):
  symbolic residual-barrier analysis, primary literature/code references, and
  148 re-extracted native node records; version-specific RefDens controls proposed.
  Read-only investigation, no new solve or qualified point.
- [SLOT IdVd local nonlinear trust-region control, 2026-10-07](validation/slot_ldmos_idvd_local_trust_region_2026-10-07.md):
  fresh residuals and Jacobians enable three/two guarded diagnostic steps;
  thirty/twenty-nine failed rows remain. One closure-probe timeout, no qualified
  point or production restart; all original gates unchanged.
- [SLOT IdVd directly constrained local directions, 2026-10-07](validation/slot_ldmos_idvd_constrained_qf_direction_2026-10-07.md):
  feasible active-set directions avoid uniform clipping, but nonlinear row
  checks still fail; directional differences support finite-step curvature.
  Twelve frozen probes, no Newton restart or qualified point.
- [SLOT IdVd bounded local QF repair, 2026-10-07](validation/slot_ldmos_idvd_local_qf_repair_2026-10-07.md):
  six frozen local trials with 30/50 variables; neighboring-row guards avoid
  new failures but do not resolve the original thirty. No Newton restart or
  qualified point; unchanged physics and acceptance gates.
- [SLOT IdVd QF mobility-feedback audit, 2026-10-07](validation/slot_ldmos_idvd_qf_feedback_audit_2026-10-07.md):
  stabilized matrix/direction and finite-difference evidence, frozen single-row
  scan, and a residual-preserving diagnostic linearization; the bounded solve
  remains unqualified, with original gates and ledgers unchanged.
- [SLOT IdVd LLM technical handoff, 2026-10-06](validation/slot_ldmos_idvd_llm_handoff_2026-10-06.md):
  self-contained context, qualification boundaries, source mapping controls,
  failure trajectories and node/material/angle evidence for external analysis.
- [SLOT IdVd contact-field fallback bounded control, 2026-10-06](validation/slot_ldmos_idvd_contact_field_fallback_2026-10-06.md):
  205 contact-adjacent cells audited; source-side holes substantially restored,
  but three bounded solves produce no qualified state. Thirty local rows remain;
  original gates and ledgers unchanged, no long continuation or float80 claim.
- [SLOT IdVd spatial, model-switch and precision controls, 2026-10-06](validation/slot_ldmos_idvd_spatial_ablation_precision_2026-10-06.md):
  four matched drain points through 2.3683165 V; native avalanche/HFS/64-bit
  controls complete. Frozen probes prioritize contact avalanche driving-force
  mapping; no changed-model qualified continuation or Vela float80 claim.
- [SLOT IdVd 60 V secant campaign, 2026-10-06](validation/slot_ldmos_idvd_60v_campaign_2026-10-06.md):
  user-authorized nine-hour local run from Vd=0.310540 V, with twenty-minute
  monitoring; consult live summary and exit evidence for progress, not this launch record.
- [SLOT IdVd opt-in secant CLI integration, 2026-10-06](validation/slot_ldmos_idvd_secant_integration_2026-10-06.md):
  explicit drain predictor and validated history import; two new drain points
  qualify through 0.310540 V in 326 s with unchanged gates. Saved-evidence audit
  passes; no long scan or numerical failure-recovery claim.
- [SLOT IdVd secant neighbors and same-target reproduction, 2026-10-06](validation/slot_ldmos_idvd_secant_neighbors_2026-10-06.md):
  three adjacent drain points qualify through 0.210913 V; alternate secant
  reproduction qualifies, while anchor warm-start reproduction fails the original
  row and closure gates; diagnostic only, no production ledger update.
- [SLOT IdVd warm-start/secant bounded comparison, 2026-10-06](validation/slot_ldmos_idvd_secant_comparison_2026-10-06.md):
  secant directly qualifies Vd=0.1134313 V in 182 s; warm-start exhausts its
  1200 s arm budget after qualifying only the intermediate point. Original
  gates remain unchanged; this is not a paired-endpoint speedup measurement.

- [SLOT IdVd bounded resume, 2026-10-06](validation/slot_ldmos_idvd_resume_2026-10-06.md):
  resumes after the user pause from qualified Vg=2 V / Vd=0 V, with explicit
  parent revalidation, separate imported/new counts and the remaining budget.
  Consult the new campaign summary for current progress.

- [SLOT IdVd gate-initialization integration and bounded full-path run, 2026-10-05](validation/slot_ldmos_idvd_gate_path_integration_2026-10-05.md):
  opt-in Fermi Poly/Poisson initialization is integrated with original strict
  acceptance and retry isolation. New cold and first gate point pass; the
  campaign was paused by the user at Vg=2 V; see the October 6 resume above.

- [SLOT IdVd Fermi Poly initialization validation, 2026-10-05](validation/slot_ldmos_idvd_fermi_poly_initialization_2026-10-05.md):
  rigid Poly bias prediction passes both direct-Newton and Poisson-then-Newton
  qualification at Vg=0.1 V, followed by qualified 0.2/0.3 V points at Vd=0.
  Original gates are unchanged; the full gate/drain path and native agreement
  remain unvalidated. The predictor is explicit and not yet in the default driver.

- [SLOT IdVd combination implementation and cold validation, 2026-10-05](validation/slot_ldmos_idvd_combination_validation_2026-10-05.md):
  explicit Fermi/PhuMob/Lombardi/HFS plus Si-only avalanche and Poly constant
  mobility pass targeted regressions and cold-equilibrium gates. Four bounded
  nonzero-gate attempts fail the original gates; no qualified IdVd current curve
  is available. This element-box profile does not reuse the external AverageBox
  transport geometry, whose compatibility guard remains in force.

- [SLOT IdVd native baseline and compatibility preflight, 2026-10-05](validation/slot_ldmos_idvd_preflight_2026-10-05.md):
  native main-deck Vg=5 V / Vd=0--60 V completes with 93 drain records; Vela
  preflight identified Fermi/PhuMob/avalanche and material-scope compatibility
  blocks. See the newer implementation report above for the subsequent work;
  this preflight produced zero qualified Vela IdVd points.

- [SLOT IDL independent materials and cold-path validation, 2026-10-05](validation/slot_ldmos_idl_cold_material_validation_2026-10-05.md):
  a separate Bennett-Wilson material bundle completes cold equilibrium, Vd=0.1 V
  and Vg=0--5 V; 89 accepted states, unchanged gates, 0.01111% maximum current
  difference at 11 native matched-gate samples. Native geometry remains required.

- [SLOT IDL PolySilicon material attribution, 2026-10-05](validation/slot_ldmos_idl_poly_material_attribution_2026-10-05.md):
  coherent intrinsic density and energy-reference mapping explains the previous
  threshold-region current excess; local interventions precede the cold path above.

- [SLOT IDL Bennett-Wilson continuation, 2026-10-05](validation/slot_ldmos_idl_bennett_continuation_2026-10-05.md):
  Vd=0.1 V and the 0--5 V gate path complete under the original gates; six matched-bias
  samples expose 28--29% subthreshold current differences despite 0.135% agreement at 5 V.

- [SLOT IDL Bennett-Wilson model mapping, 2026-10-05](validation/slot_ldmos_idl_bennett_wilson_mapping_2026-10-05.md):
  explicit model, parameter/unit mapping and related regressions pass; cold zero-bias
  integration passes on the SLOT mesh. Nonzero-drain and gate results are in the continuation report above.

- [SLOT IDL reference and initialization diagnostics, 2026-10-05](validation/slot_ldmos_idl_initial_validation_2026-10-05.md):
  native 0--5 V gate reference completed; runtime Bennett/Wilson BGN differs from the
  provisional Vela none configuration. No matched-physics IDL pass is claimed.

- [SLOT bounded breakdown continuation, 2026-10-04](validation/slot_ldmos_s04_breakdown_goal_2026-10-04.md):
  38 accepted points reach 6.27688e-8 A/um under the original gates; the next point exhausts
  bounded retries. Evidence audit passes, but the 1e-7 current goal remains unmet.

- [SLOT damped adjacent points and same-supply validation, 2026-10-03](validation/slot_ldmos_s03_damped_adjacent_same_supply_2026-10-03.md):
  full recovery followed by half-step polishing passes two target points and one independent
  sample; the third target stops at the original local-patch cap, with all costs retained.

- [SLOT carrier-block damping comparison, 2026-10-03](validation/slot_ldmos_s03_carrier_damping_comparison_2026-10-03.md):
  half steps suppress the alternating correction at both qualified inputs; eight cases
  pass original gates, with actual-supply confirmation kept distinct from same-supply reproduction.

- [SLOT symmetric precision revalidation, 2026-10-03](validation/slot_ldmos_s03_symmetric_precision_revalidation_2026-10-03.md):
  two original-gate repair cycles reveal alternating current differences; own-bias controls
  reproduce the effect, tighter solves fail, and no production points are promoted.

- [SLOT same-parent workflow and binary comparison, 2026-10-03](validation/slot_ldmos_s03_same_parent_workflow_comparison_2026-10-03.md):
  both binaries reproduce the old recovery result; current recovery also passes and
  closes the fixed target, with all attempt costs and remaining path sensitivity reported.

- [SLOT native-settings and row-threshold review, 2026-10-03](validation/slot_ldmos_s03_native_settings_row_review_2026-10-03.md):
  original deck/log and archived successful branch checked; six physical inputs match,
  original row gate previously crossed the criterion, workflow regression prioritized.

- [SLOT recovered-state target correction and independent confirmation, 2026-10-03](validation/slot_ldmos_s03_recovered_target_reproduction_2026-10-03.md):
  both three-map arms pass original device and fixed-supply gates at 1666.4974 V;
  recovery costs retained, 108 regressions pass, diagnostic only and no ledger promotion.

- [SLOT saved Poisson-failure state recovery, 2026-10-03](validation/slot_ldmos_s03_poisson_failure_recovery_2026-10-03.md):
  one carrier correction and two Poisson updates pass original device and actual-supply
  gates; target supply remains 0.642 V away, diagnostic state only, no ledger promotion.

- [SLOT next continuation stopped at pre-carrier Poisson repair, 2026-10-03](validation/slot_ldmos_s03_port_continuation_stop_2026-10-03.md):
  first 10% target has a converged coupled candidate but Poisson line-search
  failure; zero new points, original gates preserved, last qualified parent retained.

- [SLOT port-corrected three-point continuation and same-supply samples, 2026-10-03](validation/slot_ldmos_s03_port_corrected_continuation_2026-10-03.md):
  three original-gate new pilot points and two independent samples pass through
  1.480e-9 A/um; full cost recorded, recovery retained, no long scan started.

- [SLOT same-supply port correction and independent reproduction, 2026-10-03](validation/slot_ldmos_s03_same_supply_reproduction_2026-10-03.md):
  both recovered seeds pass original fixed-supply and device gates after three
  bounded maps each; five coupled exits require recovery, no new bias advance.

- [SLOT carrier and Poisson recovery of independent-seed failures, 2026-10-03](validation/slot_ldmos_s03_carrier_poisson_recovery_2026-10-03.md):
  both saved failures recover through original device and actual-supply gates;
  one carrier-block correction resolves the rows, while same-inner current
  differences remain and no new bias points are added.

- [SLOT predictive-selection switch and independent-seed validation, 2026-10-03](validation/slot_ldmos_s03_switch_independent_validation_2026-10-03.md):
  same-seed on/off both pass three cases; two different-seed reproductions fail
  row and patch-budget gates, so no new bias advance; original criteria preserved.

- [SLOT predictive neighboring-row selection and three-point validation, 2026-10-03](validation/slot_ldmos_s03_predictive_neighbor_validation_2026-10-03.md):
  default-off first-cycle FD-impact expansion; three original-gate consecutive
  points accepted through 1.112e-9 A/um, unchanged production ledger and gates.

- [SLOT patch-boundary propagation and anticipatory repair, 2026-10-02](validation/slot_ldmos_s02_patch_boundary_2026-10-02.md):
  same-seed six/eight-row initial patches recover the remaining 10% input with
  original gates; boundary transfer quantified, diagnostic states only.

- [SLOT bounded neighboring-row repair, 2026-10-02](validation/slot_ldmos_s02_neighbor_row_repair_2026-10-02.md):
  default-off local QF repair recovers two of three failed inputs with unchanged
  gates and actual-supply confirmation; one remains rejected, no scan restart.

- [SLOT row-tolerance, ULP and local-repair diagnostic, 2026-10-02](validation/slot_ldmos_s02_row_tolerance_diagnostic_2026-10-02.md):
  relaxed tolerances merely accept the unchanged seed; a damped two-node repair
  passes original device gates with unchanged evaluated drain current; PETSc/Wolfram
  parameter comparison, no production promotion or scan restart.

- [SLOT adjacent-point, fixed-supply and goal continuation, 2026-10-02](validation/slot_ldmos_s02_sequence_port_goal_2026-10-02.md):
  three original-gate adjacent points pass; a same-seed repair-path control resolves
  the fixed-supply discrepancy to 0.100 microvolt; bounded goal continuation is launched.

- [SLOT six-row repair and merit controls, 2026-10-02](validation/slot_ldmos_s02_six_row_repair_2026-10-02.md):
  one carrier correction followed by Poisson repair passes original device and
  actual-supply gates in both coordinate controls; target-supply reproduction remains separate.

- [SLOT combined-path replay and actual-supply confirmation, 2026-10-02](validation/slot_ldmos_s02_repeat_loadline_2026-10-02.md):
  deterministic replay and original actual-supply confirmation pass; independent
  parent-seeded control still fails six electron rows, with no ledger promotion.

- [SLOT merit trials and poly/Poisson recovery controls, 2026-10-02](validation/slot_ldmos_s02_merit_f3_2026-10-02.md):
  rejected full trial improves Si rows; bounded combined repair passes original
  fixed-inner device checks, while target load-line confirmation remains open.

- [SLOT ULP, fixed-inner and bounded coupled-depth diagnostics, 2026-10-02](validation/slot_ldmos_s02_ulp_f1_f2_2026-10-02.md):
  controlled QF sensitivity, five unqualified fixed-inner trials, and a one-update
  coupled improvement blocked by merit on further correction; original gates retained.

- [SLOT fold current-target control and same-supply controls, 2026-10-01](validation/slot_ldmos_fold_current_target_2026-10-01.md):
  opt-in fold-safe proposal/history, 360 recorded-state decisions, and bounded
  A/B/C controls; all candidates fail local-row qualification, no new accepted point.

- [SLOT full-state Jacobian tangent predictor, 2026-10-01](validation/slot_ldmos_full_state_tangent_2026-10-01.md):
  default-off predictor, unchanged strict gates, four qualified same-parent
  controls; independent same-bias reclosure exhausts its bounded budget.

- [SLOT short continuation after validation reuse, 2026-09-30](validation/slot_ldmos_reuse_short_continuation_2026-09-30.md):
  two six-step paths and two independent controls pass unchanged gates;
  real bounded recovery succeeds, while full-run speedup remains unproven.

- [SLOT validation reuse and bounded candidate correction, 2026-09-30](validation/slot_ldmos_validation_reuse_candidate_correction_2026-09-30.md):
  unified proof adoption, content-checked process-local reuse, 145 regressions,
  and fixed-candidate controls; extra correction does not improve cost.

- [SLOT candidate-current bias and repeated validation cost, 2026-09-30](validation/slot_ldmos_candidate_bias_validation_cost_2026-09-30.md):
  read-only balance decomposition and state/proof comparison identify incomplete
  candidate closure, a repeated validation, and reusable parsing costs; no new solve.

- [SLOT conditional row repair and 10/12.5 mV short batches, 2026-09-30](validation/slot_ldmos_conditional_row_repair_batch_2026-09-30.md):
  default-off one-pass repair, 137 regressions, 12 qualified branch points and
  two independent controls; row quality improves without a throughput gain.

- [SLOT nodes 3595/3603 correction-cost check, 2026-09-30](validation/slot_ldmos_nodes3595_3603_correction_cost_2026-09-30.md):
  two one-pass carrier repairs restore row-quality margin at unchanged bias and
  gates; measured full diagnostic overhead motivates conditional integration.

- [SLOT point-6 reclosure and step-quality gate, 2026-09-30](validation/slot_ldmos_point6_reclosure_step_gate_2026-09-30.md):
  bounded reclosure removes the repeatability gap; 8/10 mV same-parent trials pass
  final gates, but the 10 mV hole-row quality prevents further step expansion.

- [SLOT matched-feedback same-parent short batch, 2026-09-30](validation/slot_ldmos_matched_feedback_batch_2026-09-30.md):
  six accepted steps per mode on T470p, 11.94% higher local supply advance per
  second, working candidate-quality backoff and a same-bias repeatability caveat.

- [SLOT matched-supply feedback implementation, 2026-09-30](validation/slot_ldmos_matched_feedback_fix_2026-09-30.md):
  opt-in frozen-model evaluation at accepted supply, conservative fallback,
  identity binding, 119 regressions and 18-point replay; no new device solve.

- [SLOT prediction feedback coordinate review, 2026-09-30](validation/slot_ldmos_prediction_feedback_review_2026-09-30.md):
  causal replay of 18 accepted steps separates actual-supply migration from
  same-supply model error; production policy and final qualification unchanged.

- [SLOT finite step expansion from 5 to 10 mV, 2026-09-30](validation/slot_ldmos_finite_step_expansion_2026-09-30.md):
  three qualified same-parent trials, six strict advances and two independent
  controls; port prediction feedback limits sustained step growth.

- [SLOT point-3 bounded same-bias reclosure, 2026-09-30](validation/slot_ldmos_point3_bounded_reclosure_2026-09-30.md):
  ordinary restart performed no update; one carrier correction followed by a full
  Newton update removed the current discrepancy under unchanged final gates.

- [SLOT current-scaled stepping and bounded growth validation, 2026-09-30](validation/slot_ldmos_current_growth_batch_2026-09-30.md):
  four same-parent step trials, 12 qualified advances with quality-triggered
  backoff, two independent controls, and a residual-related repeatability caveat.

- [SLOT daily report and figures, Sept 29–30 morning](validation/slot_ldmos_daily_report_2026-09-30.md):
  management-ready summary, actual device cross-section, 12-point continuity and
  fixed-parent step comparison; qualified voltage is not a BVDS result.

- [SLOT CLI integration and short continuation, 2026-09-30](validation/slot_ldmos_short_cli_batch_2026-09-30.md):
  12 qualified T470p advances, two independent controls, three same-parent step
  controls, 47 regressions and downloaded evidence independently revalidated.

- [SLOT higher-bias candidate conservation and repair cost, 2026-09-29](validation/slot_ldmos_candidate_conservation_high_parents_2026-09-29.md):
  three local parent advances, independent same-bias controls, complete qualification,
  and measured candidate/strict-reclosure costs.

- [SLOT frozen candidate conservation, 2026-09-29](validation/slot_ldmos_candidate_conservation_2026-09-29.md):
  contact reference/increment preservation, independent candidate closure,
  unchanged strict reclosure and bounded local verification.

- [SLOT progress/predictor repair and bounded diagnostics, 2026-09-29](validation/slot_ldmos_progress_predictor_diagnostics_2026-09-29.md):
  opt-in stagnation guard, matched predictor voltage, same-state conservation
  diagnostics, strict local step controls, and remaining candidate-current mismatch.

- [SLOT zero-device-RHS option and three-parent validation, 2026-09-29](validation/slot_ldmos_zero_rhs_multi_parent_2026-09-29.md):
  default-off candidate option, unchanged strict final gates, three additional
  T470p controls, and revalidated 35 V / latest high-voltage restart windows.

- [SLOT candidate forcing comparison on T470p, 2026-09-29](validation/slot_ldmos_candidate_forcing_comparison_2026-09-29.md):
  same-parent strict/inexact controls, unchanged final qualification, successful
  zero-device-RHS candidate at the original tolerance, and independent reclosure.

- [SLOT local 35 V matched-bias review and high-voltage continuation, 2026-09-28](validation/slot_ldmos_home_35v_to_bv_2026-09-28.md):
  matched current/field/source comparison, independently requalified origin-36
  parent window and local continuation toward the original current criterion.

- [SLOT bounded recovery and independently audited continuation through 35 V, 2026-09-28](validation/slot_ldmos_bounded_recovery_2026-09-28.md):
  unchanged qualification gates, actual-supply circuit confirmation, reference-frame
  control and 44 accepted T470p points from 31.3 to 35 V external bias.

- [SLOT hybrid prediction through 30 V on T470p, 2026-09-28](validation/slot_ldmos_hybrid_continuation_30v_2026-09-28.md):
  28 new resistor points and one frame requalification, actual rational-candidate
  use, independent audits and quantified post-frame step-recovery cost.

- [SLOT hybrid port predictor and T470p pilot, 2026-09-28](validation/slot_ldmos_hybrid_port_predictor_2026-09-28.md):
  opt-in limited endpoint baseline, causal rational-candidate selection, 87
  regressions and independently audited 26.77-to-27.2 V forward validation.

- [SLOT port step controller and offline prediction, 2026-09-27](validation/slot_ldmos_port_step_prediction_2026-09-27.md):
  opt-in independent outer, accepted-inner and trial-inner caps, qualified
  three-point slope fitting, cost/error feedback and read-only evidence replay.

- [SLOT native-goal adaptive BV campaign, 2026-09-27](validation/slot_ldmos_native_goal_campaign_2026-09-27.md):
  launched from the qualified 3 V resistor state toward 200000 V external or
  the 1e-7 A/um current criterion, with qualified-parent prediction and adaptive steps.

- [SLOT predicted resistor continuation through 3 V, 2026-09-27](validation/slot_ldmos_predictor_continuation_3v_2026-09-27.md):
  conservative adaptive continuation from the qualified 2 V state, observed
  first-step recovery, unchanged gates and full attempt-cost accounting.

- [SLOT state prediction and adaptive stepping, 2026-09-27](validation/slot_ldmos_predictor_adaptive_validation_2026-09-27.md):
  opt-in qualified-parent prediction and success-based growth, unchanged fallback
  and acceptance gates, with separate fixed-path and adaptive controls through 2 V.

- [SLOT resistor continuation through 2 V, 2026-09-27](validation/slot_ldmos_low_voltage_continuation_2026-09-27.md):
  serial continuation from the recovered 0.8 V state, unchanged device/circuit
  gates, and bounded Windows checkpoint replacement recovery.

- [SLOT candidate failure classification and 0.8 V recovery, 2026-09-27](validation/slot_ldmos_candidate_failure_recovery_2026-09-27.md):
  bound numerical-failure protocol, original-state replay and qualified 0.7–0.8 V
  resistor recovery with unchanged device and circuit gates; full 0-to-BV scan remains open.

- [SLOT full resistor campaign, 2026-09-26](validation/slot_ldmos_full_resistor_campaign_2026-09-26.md):
  in-progress local equilibrium-to-BV campaign, serial execution and checkpoint recovery;
  this work has not yet qualified a complete new-configuration 0-to-BV path.

- [SLOT native transport resistor branch and BVDS comparison, 2026-09-26](validation/slot_ldmos_native_transport_bv_branch_2026-09-26.md):
  qualified 38.2/38.4 V controls, audited resistor branch across the current criterion,
  independent threshold refinement and interpolation comparison with Sentaurus.

- [SLOT native transport 38 V reclosure and continuous path, 2026-09-26](validation/slot_ldmos_native_transport_38v_continuation_2026-09-26.md):
  independent SC/IIC controls, qualified 25-point 35–38 V continuation,
  matching endpoints and external transport resource checks; the subsequent BVDS
  criterion comparison is linked above.

- [SLOT full native transport geometry and qualified 35 V reclosure, 2026-09-26](validation/slot_ldmos_native_transport_35v_2026-09-26.md):
  shared opt-in source guard, assembled/Jacobian checks and independent IIC/SC
  reclosure; the subsequent 38 V and continuous-path evidence is linked above.

- [SLOT near-surface hole source and transport localization, 2026-09-26](validation/slot_ldmos_near_surface_hole_transport_2026-09-26.md):
  same-state source replay, density/field/current factor budget, and frozen
  edge-weight plus mobility-measure control; independent reclosure remains pending.

- [SLOT circumcentric continuation through 38 V, 2026-09-26](validation/slot_ldmos_circumcentric_continuation_38v_2026-09-26.md):
  qualified continuous 35–38 V path, bounded density/QF candidate recovery,
  and remaining local hole-source versus non-Delaunay measure diagnostics.

- [SLOT circumcentric volume localization and 35–38 V controls, 2026-09-26](validation/slot_ldmos_circumcentric_volume_35_38v_2026-09-26.md):
  opt-in charge/recombination measures, independent 35 V gates and high-voltage continuation evidence.

- [SLOT per-cell dielectric operator and independent 35 V reclosure, 2026-09-26](validation/slot_ldmos_per_cell_reclosure_35v_2026-09-26.md):
  opt-in operator, full-device replica acceptance and qualified IIC/SC points;
  source/current improvement with remaining bulk-potential and IIC-source differences.

- [SLOT 35 V IIC nonlinearity and H2 frozen controls, 2026-09-26](validation/slot_ldmos_iic_nonlinearity_h2_35v_2026-09-26.md):
  state-path cancellation checks and charge-volume/permittivity controls;
  promising joint variant, replica-gate limitations, no production policy change.

- [SLOT 35 V adjoint localization and H1 control, 2026-09-26](validation/slot_ldmos_iic_loop_gain_localization_35v_2026-09-26.md):
  verified linear probes, two-segment diagnostic attribution, unsuccessful H1 promotion
  and the remaining IIC linearization discrepancy; no production policy change.

- [SLOT IIC source budget at 35/37/38 V, 2026-09-25](validation/slot_ldmos_iic_source_budget_2026-09-25.md):
  integration-support checks, frozen production-kernel replay, independently
  qualified no-feedback states and a self-consistent source-amplification budget.

- [SLOT 38.2/38.4 V comparison and high-resistance branch, 2026-09-25](validation/slot_ldmos_near_breakdown_branch_2026-09-25.md):
  qualified direct continuation, native field/source comparison and independently
  audited resistor points; the local breakdown-current criterion is not yet crossed.

- [SLOT local-origin continuation through 38 V and resistor start, 2026-09-25](validation/slot_ldmos_local_origin_38v_loadline_2026-09-25.md):
  eight newly qualified voltages, native 37.5/38 V comparisons, and an independently
  audited original-resistor start plus a 0.01 V outer-voltage advance.

- [SLOT local-origin continuation through 37 V, 2026-09-25](validation/slot_ldmos_local_origin_37v_2026-09-25.md):
  57 newly qualified voltages, 63 audited records, native 35/37 V field/source
  comparisons, origin 35 V restart, and tighter-floor checkpoint revalidation.

- [SLOT local-origin continuation through 30 V, 2026-09-25](validation/slot_ldmos_local_origin_30v_2026-09-25.md):
  81 newly qualified voltages, 87 audited records, serial native 25 V reference,
  25/30 V field/source comparisons, origin 15 V restart and residual-floor limits.

- [SLOT local-origin continuation through 20 V, 2026-09-24](validation/slot_ldmos_local_origin_20v_2026-09-24.md):
  48 newly qualified voltages, 52 audited records, guarded Poly restart support,
  native 15/20 V current, field and avalanche-source comparisons, and restart anchors.

- [SLOT local-origin continuation through 10 V, 2026-09-24](validation/slot_ldmos_local_origin_high_voltage_2026-09-24.md):
  36 newly qualified voltages, 40 audited qualification records, serial Release
  continuation in origin 5 V, native 5/10 V comparisons and final restart entry.

- [SLOT reference-frame continuation tools, 2026-09-24](validation/slot_ldmos_reference_frame_tools_2026-09-24.md):
  explicit frame conversion, inactive contact constraints, requalification of both
  predictor anchors, physical-voltage reporting and Release low-bias integration.

- [SLOT 10 V electrical-frame audit, 2026-09-24](validation/slot_ldmos_10v_frame_audit_2026-09-24.md):
  shifted-frame qualification, original-frame revalidation failure, contact-floor
  semantics and inactive contact projection; not continuous cold-origin 10 V or BVDS completion.

- [SLOT local-origin low-voltage completion, 2026-09-24](validation/slot_ldmos_local_origin_low_voltage_2026-09-24.md):
  36 qualified points from local equilibrium to 2.125 V including restart recovery,
  six-point path comparison, serial 92-test tool verification; high-voltage/BVDS remains open.

- [SLOT HDF5 continuation and native low-bias comparison, 2026-09-24](validation/slot_ldmos_hdf5_continuation_2026-09-24.md):
  three-driver state migration, block-probe output fix, Release regression,
  checkpoint confirmation, repeat-chain progress and native 2.025 V evidence.

- [SLOT status and handoff, 2026-09-24](validation/slot_ldmos_status_and_handoff_2026-09-24.md):
  original-deck function and parameter matrix, completed versus unfinished
  validation, current HDF5 anchors and next-stage entry; separates projected-SG
  high-voltage results from the controlled new-source 2.125 V qualification.

- [SLOT local main integration, 2026-09-23](validation/slot_ldmos_main_integration_2026-09-23.md):
  HDF5 restart compatibility, default Jacobian structure reuse, retained SLOT
  physics fixes and bounded Release revalidation; continuation driver migration
  remains separate from high-voltage and BVDS qualification.

- [SLOT mixed-contact guard and escape controls, 2026-09-23](validation/slot_ldmos_contact_guard_and_escape_controls_2026-09-23.md):
  local contact-majority classification fix, T1 update traces, HFS-off and
  single-node residual controls, ULP sensitivity limits and unchanged-gate
  device revalidation. A 0.025 V continuation now qualifies new-source points
  through 2.125 V; the 10 V observed-source checkpoint remains rejected.

- [SLOT comparison analysis for LLM review, 2026-09-23](validation/slot_ldmos_comparison_analysis_for_llm_2026-09-23.md):
  Chinese evidence summary, separate projected-SG and observed-source qualification,
  qualified new-source 2 V, failed 2.125/10 V controls, Jacobian coverage limits,
  units, outstanding hypotheses and concrete questions for independent analysis.

- [SLOT observed element source integration, 2026-09-23](validation/slot_ldmos_observed_source_integration_2026-09-23.md):
  complete fixed-mesh native selector coverage, independent screening-floor
  discrepancy evidence, controlled live source/Jacobian integration and explicit
  qualified 2 V and failed 2.125/10 V coupled controls. Strict original BVDS remains unqualified.

- [SLOT contact avalanche kernel fix, 2026-09-23](validation/slot_ldmos_contact_kernel_fix_2026-09-23.md):
  contact/interface full-gradient correction, analytic regression, full local
  nine-column live-mobility Jacobian checks and native 10/37 V replay.
  Controlled coupled assembly and 2 V qualification are recorded above;
  strict BVDS remains unqualified.

- [SLOT PMI callback follow-up, 2026-09-23](validation/slot_ldmos_pmi_callback_followup_2026-09-23.md):
  complete callback-to-element matching, contact full-gradient current evidence,
  effective current scale, saved-state precision limits and qualified frozen
  AvalDensGradQF-off control. The contact-kernel correction is now verified in
  the contact-kernel report; subsequent controlled assembly is recorded above.

- [SLOT native avalanche identification, 2026-09-23](validation/slot_ldmos_native_avalanche_identification_2026-09-23.md):
  compiler-independent frozen controls, extrema harmonic density, cosine alpha
  smoothing, best-vertex mobility and 10/37 V interior replay with explicit
  geometry-choice uncertainty, native truncated source volumes and a standalone
  nine-column Jacobian kernel. GCC 11 PMI loading and frozen-output parity now
  pass, with direct near-tie vertex-choice evidence. Subsequent fixed-mesh selector
  completion and controlled coupled integration are recorded above.
  Original AvalDensGradQF and complete BVDS remain unqualified.

- [SLOT load-line driver review, 2026-09-23](validation/slot_ldmos_loadline_review_2026-09-23.md):
  guarded independent-device continuation, original all-row/circuit gates,
  reviewed restart identity and budgets, and real 37.075/37.1 V replay evidence.
  Original AvalDensGradQF and complete BVDS remain separate.

- [SLOT completed internal continuation, 2026-09-22](validation/slot_ldmos_continuation_completion_2026-09-22.md):
  both resumed segments complete through 20/37 V, 49 qualified new points,
  and 20/30/35/37 V field/current overlap with independent checkpoints.
  Cold-origin replay and original avalanche/BVDS validation remain separate.

- [SLOT commit handoff and next tasks, 2026-09-22](validation/slot_ldmos_next_tasks_2026-09-22.md):
  committed fixes, a timestamped live-run snapshot, and concrete acceptance
  criteria for continuation, load-line automation, original avalanche physics,
  error budgeting, and final BVDS qualification.

- [SLOT high-voltage continuation, 2026-09-22](validation/slot_ldmos_high_voltage_resume_2026-09-22.md):
  independently qualified projected-SG checkpoints through 37 V, complete
  Jacobian checks, guarded Poly restarts, compensated-QF restart repairs, and
  qualified 37.02/37.05 V external-load-line controls. Original AvalDensGradQF source
  averaging and full BVDS qualification remain unresolved.

- [SLOT energy reference and source volume, 2026-09-21](validation/slot_ldmos_energy_reference_and_source_volume_2026-09-21.md):
  band-derived poly reference, material-local SRH/Auger integration, qualified
  2 V error budget and independent local-equilibrium continuation status.

- [SLOT Auger and poly independent controls, 2026-09-21](validation/slot_ldmos_auger_poly_independent_controls_2026-09-21.md):
  density-enhanced Auger controls at 2/5/10 V, isolated poly transport and
  material-scoped SRH/Auger, qualified 2 V restoration and remaining energy-reference gap.

- [SLOT element IALMob execution, 2026-09-21](validation/slot_ldmos_element_ialmob_execution_2026-09-21.md):
  original isotropic parameter and element replay, strict 2 V checkpoint,
  single-factor error budget and complete mobility/avalanche Jacobian work.
  High-voltage and external-resistor qualification are recorded separately.

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
- [SLOT-LDMOS BVDS validation](validation/slot_ldmos_bvds_validation.md):
  local and global qualification profiles, current low-voltage evidence,
  and remaining physical-model differences.
- [SLOT/main integration and template comparison, 2026-09-21](validation/slot_ldmos_main_integration_2026-09-21.md):
  retained SLOT compatibility, imported template capabilities, and separate
  numerical qualification boundaries.

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
