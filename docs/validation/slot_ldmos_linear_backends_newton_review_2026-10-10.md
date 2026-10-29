# SLOT：线性后端实测与Sentaurus/PETSc非线性算法核查（2026-10-10）

## 结果

已实际测试UMFPACK、Eigen SparseLU、SuperLU_MT、STRUMPACK、GMRES+ILUT、BiCGSTAB+ILUT。六个后端在同一5 mV稳定化系统上均返回方向，物理矩阵、正则化矩阵、闭合替换矩阵及RHS身份全部相同。**返回线性方向不等于器件收敛**；六个方向经过同一有限步处理后均未通过原行比/守恒要求，新增完整合格点0。

本父态不支持“只换一个线性后端即可解决非线性不收敛”。UMFPACK与BiCGSTAB在本次稳定化系统上的分量后向误差更小，但只是单系统结果，不能据此把BiCGSTAB改为生产默认。若继续推进，优先处理实际有限步的全局化、跨块尺度及人口正性协调。

## 共同输入与实验层次

沿用[统一边界与有限步对照](slot_ldmos_boundary_stabilization_finite_steps_2026-10-10.md)的最新兼容配置，Vg=5 V、Vd=0.005 V、300 K、13,881节点、41,643维QF系统。实际AverageBox/Fermi/effective-band SG及完整目标物理不变；原资格门槛不变。

冻结runner为`build/s10four_stage/compat_source/build-release/vela_example_runner.exe`，SHA-256 `f73560bea030e4881ae369bf0347f1344d5a7c58c830687a0664cb50eeab4849`。父态SHA-256 `9db915c19425d3936fab2895751a7569c7381e2cc4a753321e52c7166463a9ff`。使用已有编译后端；SuperLU具体为SuperLU_MT OpenMP版，STRUMPACK为无压缩DIRECT/METIS配置，并非有损压缩或GPU模式。迭代后端为仓库现有Eigen GMRES/BiCGSTAB+IncompleteLUT，容差1e-14、最多4000次，GMRES restart=200。线程请求为1，BLAS请求为1。

本批不启用PETSc库、不连接原生虚拟机、不运行长扫描。Sentaurus依据本地版本绑定手册及既有原生日志；PETSc依据官方文档/源码。线性后端只更换`VELA_LINEAR_SOLVER`及输出路径，物理配置一致。

**层A：生产预限幅方向。** `newton_linear_system_probe`捕获当前正则化及replace_free_row_sum系统，生产实现使用低秩Woodbury对闭合替换处理。每个后端解相同基础矩阵及辅助RHS。accepted_updates=0；工具JSON的converged=true只表示探针成功，其qualified=false始终保留。

**层B：真实有限步。** 从层A导出的方向出发，采用相同0.05/0.1 V载流子限幅、统一UMFPACK物理ψ块线性重修正和接触保护，在同一冻结评估器重算残差/行比/Dirichlet前闭合。固定ψ重修正后端，是为了隔离主方向的后端差异；不把有限试探接受为Newton迭代。

**辅助层P：无正则化物理方向。** 另执行现有`newton_step_probe`。源码核对该接口直接解物理J，再限幅/ψ重修正，**不包含层A的正则化及闭合替换**。故单独报告，不能冒充层A方向的生产首步。未改变该接口源码。

## 层A线性结果

整体相对线性残差为||Ax-b||₂/||b||₂。分量后向误差为max_i |(Ax-b)_i|/(|A||x|+|b|)_i，使用long double重新累加检查；自由行按原生变量身份及配置接触消元确定。方向差以UMFPACK为参考，不将其视为数学真解。

|后端|整次捕获墙钟/s|整体相对线性残差|自由行最大分量后向误差|与UMFPACK最大方向差/V|最大差节点|
|---|---:|---:|---:|---:|---:|
|umfpack|31.25|2.339e-15|1.037e-14|0|0|
|sparselu|36.68|5.688e-15|1.000e+00|6.50097|10057|
|superlu_mt|40.76|5.674e-15|1.000e+00|4.16792|10099|
|strumpack|35.70|3.106e-15|1.000e+00|3.59742|10404|
|gmres_ilut|36.79|9.599e-11|9.455e-01|0.00111409|0|
|bicgstab_ilut|38.42|2.703e-15|1.548e-14|0.00557047|0|

所有后端最大原始QF方向仍约8752.3 V。表中时长含装配、求解、输出，不是纯因子分解性能基准；单轮结果不作速度排名。详细块/节点身份与约束行绝对误差见`comparison.json`。

SparseLU/SuperLU_MT/STRUMPACK的分量相对误差接近1，必须结合绝对尺度解释：其最坏行线性缺陷约1.56e-14、1.44e-14、5.41e-21（内部行单位）。近零RHS/微小行的分母也同样小，不能把该数字解读成整个器件方程100%误差；但它也说明整体残差1e-15不能保证全部弱行方向准确。GMRES最坏行在节点0的ψ行，线性缺陷约1.31e-8，独立计算的整体相对残差约9.60e-11，明显大于UMFPACK。内部迭代器容差和最终Woodbury组合后的整系统误差并不等价。

BiCGSTAB本点得到约1.55e-14的分量后向误差，但原始方向与UMFPACK仍最大相差约5.57mV；病态系统的小后向误差不保证小前向误差。没有用后端自身success标志替代这些外部检查。

## 层B有限步质量

|主方向后端|电子/空穴超标行|冻结组合L2|空穴L2|电子/空穴闭合比|
|---|---:|---:|---:|---:|
|umfpack|1150 / 6949|1.05038713|1.85624e-06|1.000001 / 0.9771879|
|sparselu|2647 / 7132|1.05038713|1.85624e-06|1.000001 / 0.9771879|
|superlu_mt|2584 / 7127|1.05038713|1.85624e-06|1.000001 / 0.9771879|
|strumpack|1804 / 6942|1.05038713|1.85624e-06|1.000001 / 0.9771879|
|gmres_ilut|1187 / 6949|1.05038713|1.85624e-06|1.000001 / 0.9771879|
|bicgstab_ilut|1046 / 6949|1.05038713|1.85624e-06|1.000001 / 0.9771879|

这些试探的整体组合L2仍约1.05，电子闭合比约1，明显未达0.01；故不计合格态。冻结组合范数不是生产normalizedBlockL2停止值。更换后端能够改变部分弱行数量，但不能把这些结果解释为一次完整数值修复。未作多轮非线性收敛速度比较，也未作不同后端的独立合格复求。

辅助层P的原物理方向试探：

|后端|退出码|墙钟/s|取得的试探组合L2|
|---|---:|---:|---:|
|umfpack|0|36.93|1.0503933425621186|
|sparselu|0|38.50|1.0503933425621603|
|superlu_mt|0|41.17|1.0503933425621543|
|strumpack|0|35.66|1.0503933425621395|
|gmres_ilut|0|26.49|1.050393342562196|
|bicgstab_ilut|0|25.50|1.0503933425621448|

退出码124表示单次120 s硬上限，未取得有效输出；若有此情况，仅说明该配置在此上限内未完成，不宣称该算法必不收敛。该层结果不与层B混算。

## Sentaurus手册与当前实际设置

本地手册`build/s09native_newton/sdevice_ug.txt`，T-2022.03，印刷页191–196、1203–1204：Coupled是Newton类方法，每步解JΔz=−F；阻尼决定z_new=z+λΔz。手册讨论Bank–Rose框架，但参数和日志决定某次运行是否实际启用阻尼，不能仅凭章节标题认定每步使用Bank–Rose。

本次原生矩阵来源`build/s10initial05/native/v5/idvd.cmd`及`idvd.log_des.log`明确记录：

|项目|原生实际值|对比较的影响|
|---|---|---|
|线性求解|PARDISO|不是本批测试的四种本地直接后端之一|
|精度|80 bit long double|与本地binary64不同，本批未做原生64位消融|
|LineSearchDamping|1|线搜索默认关闭|
|NotDamped|25|前几次迭代不由Bank–Rose阻尼触发条件主导|
|Digits / ErrRef(e,h)|5 / 1e10 cm⁻³|更新误差按变量及参考人口归一化|
|RhsMin|1e-10|数值单位和停止含义不可直接等同Vela abstol|
|Extrapolate|本矩阵诊断关闭|不把原扫描预测收益混入同态矩阵比较|
|Iterations|0|该版本手册明确仍执行一次迭代；仅第0次系统用于本批|

后续已单独完成 [64/80 位原生成对求解](slot_ldmos_native_precision_2026-10-10.md)：
同一 5 mV 目标两臂均在4次更新后按原生误差条件停止，漏电流相对差2.53e-11；
最终 RHS 仍有明显差异。该结果不更改本表所述原始矩阵来源和本批实验边界。

手册页193说明默认RHS条件或更新条件任一满足可停止；RhsAndUpdateConvergence通常要求两者同时满足，但还存在更小RHS的RhsMinFactor例外。ErrRef平滑连接相对与绝对更新尺度；它不是删除方程或保证弱行残差≤0.001。手册不同章节给出的聚合表达需结合版本实现理解，本报告不把它简化成所有节点逐一通过的最大范数条件。

手册页195：LineSearchDamping是允许的最小阻尼，不是每步固定α；实际阻尼低于此值时，后续处理也不同于Vela当前回溯失败立即返回。页192建议在偏压扫描内若Newton不能快速收敛，倾向缩小外层偏压步而非堆高迭代次数。这支持有界续算与快速诊断，但不是证明所有失败都应盲目缩步。

公开手册不足以重构PARDISO内部选项、全部正性投影和未公开稳定化细节。此前观察到的原生人口处理不能自动当成手册保证的完整算法。

## PETSc官方算法可借鉴之处

1. **区分线性精度与全局化。** [SNESNEWTONLS](https://petsc.org/release/manualpages/SNES/SNESNEWTONLS/)提供可选线搜索；[回溯源码](https://petsc.org/release/src/snes/linesearch/impls/bt/linesearchbt.c.html)计算实际候选F及步方向相关的下降量。对Vela尤其要记录限幅/ψ重修正后的实际方向，而不是继续假定它仍为精确Newton方向。当前`src/numerics/LineSearch.cpp`默认允许任意严格范数下降，调用方也可覆盖下降/合法性条件；因此不能把现有策略简单称为标准Armijo而忽略真实分支。
2. **把幅度约束纳入方向选择。** [SNESNEWTONTR](https://petsc.org/release/manualpages/SNES/SNESNEWTONTR/)按实际/预测改善调整信赖域，并提供Newton、Cauchy、dogleg回退。迁移到本问题需先定义电势/人口尺度，并用完整非线性残差检验模型预测。仅把巨大QF方向统一缩短不能保证越过低密度残差障碍。
3. **投影后仍需检查真正试探。** [SNESLineSearchSetPostCheck](https://petsc.org/release/manualpages/SNES/SNESLineSearchSetPostCheck/)提供方向/候选修改接口；[SNESSetConvergenceTest](https://petsc.org/main/manualpages/SNES/SNESSetConvergenceTest/)允许保留自定义验收。人口有限步、接触保护、ψ重修正之后要重新评估行质量和守恒，不能用投影前的线性预测批准状态。
4. **迭代线性容差应跟随非线性进度。** [Eisenstat–Walker接口](https://petsc.org/release/manualpages/SNES/SNESKSPSetUseEW/)可避免远离解时过度精确求解。本点UMFPACK已达约1e-14分量后向误差，因此继续提高线性精度不是现有证据支持的首要修复；若采用Krylov，要额外检查最终稳定化系统误差。
5. **人口界约束不等于原方程自动成立。** [SNESVINEWTONRSLS](https://petsc.org/release/manualpages/SNES/SNESVINEWTONRSLS/)是变分不等式活动集算法。可借鉴约束方向构造，但不能直接把人口下限处非零残差当成互补解并称满足Vela原方程。最终仍须逐项原资格审计。

PETSc网页本次读取标识为3.26.0；链接为动态release页面，版本身份记录于本报告。未声称已将PETSc接入或运行其SNES求解器。

## 建议的实施优先级

- 保留UMFPACK作基线，BiCGSTAB作为受控备选；先补最终系统的分块/分量线性误差监测，必要时进行同态边界消元、行列均衡或精化对照。不能仅看全局相对残差。
- 优先在当前共同5 mV父态做**尺度明确的信赖域/受约束方向小对照**：载流子人口正性及有限变化进入方向计算，ψ和载流子保持联动；记录真实下降/预测下降、电子闭合、空穴块及新增坏行。最终门槛不变。
- 同时保留限幅后实际方向的斜率检查与有限步重评估。只有得到实际改善，再做2–3次重新装配的连续试验；没有改善即停止该策略，不以多换后端或加时长代替定位。
- 原生更新误差可并行作为影子观测，但不替代Vela逐行及守恒资格。禁止将参考误差/VI活动集作为静默放宽验收。

## 成本与验证

层A六次探针累计219.60 s；辅助层P六次累计204.25 s；层B六态冻结评估及闭合累计237.92 s。均串行、单次上限120 s，接受Newton更新0；有界探针不等于完整求解预算或性能基准。

验证包括全部矩阵文件哈希及RHS逐值相同、扩展精度后向误差重算、输入runner/状态哈希未变、真实有限步行质量/闭合评估与脚本语法。为确认已有工具可用，构建了冻结兼容源码的`linear_solver_replay`目标（成功）；正式比较使用原runner诊断接口，未变更生产实现，未运行全套CTest。前一批边界实验已单独归档，不与本批计数混合。

证据位于`build/s10backend_newton/`：`summary.json`、`comparison.json`、`steps.json`、`stable_finite/summary.json`及各后端矩阵/方向/日志。程序返回成功仅表示诊断完成；所有有限试探均不合格。环境HDF5版本提示保留在日志。
