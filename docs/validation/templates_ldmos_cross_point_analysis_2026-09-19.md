# D5 跨点线性分析复用实验

状态：主 Newton 跨点上下文已实现；三模式×双栅压前 8 点全部通过，共 48 精确点。
相对常驻 worker，候选墙钟减少 8.57%/8.78%，CPU 减少 6.36%/6.49%。
只有一轮短曲线配对，实验开关仍默认关闭；生产默认、物理模型和验收门限不变。

## 实现边界

现有 `--dc-worker` 只缓存网格、材料、掺杂等输入准备。每次 `NewtonSolver::solve`
仍创建局部 `LinearSolver`，因此点内已有分析复用不能跨点传递。

新增 `--dc-worker-linear-reuse`，由顺序运行的 `DCSweep` 持有线性求解上下文，
通过 `NewtonConfig::sequentialLinearSolver` 显式传给主 Newton 求解。每个新请求
调用 `clearNumericCache()`，保留分析但强制重新数值分解。所有非线性状态、
接触和初值仍从当前请求重新准备，不把前一个请求的终态当作隐式初值。
辅助 Poisson/Gummel 线性求解仍走原生命周期；不保证整条曲线只有一次分析。
特别是 `recorrectPoissonStepForClippedQuasiFermi` 在 QF 分量限幅后重新校正
Poisson 步，每次调用仍使用局部求解器；它是 N×N 子块，不能与 3N×3N
主耦合 Jacobian 共用同一个单图缓存，否则会相互驱逐。后续如扩展，应给
每类固定算子独立上下文，而非把所有不同系统都视为同一张图。

缓存约束：

- 准备输入键含网格、单位、材料、掺杂、几何及文件内容；输入缓存失效同时
  释放线性上下文。线性层仍逐次核对矩阵维数、压缩外指针及内索引。
- 图结构变化重建后端；线性异常清空缓存；worker 请求失败或异常释放上下文。
- 跨请求不复用数值因子；同请求内已有的完全相同矩阵因子复用不变。
- 上下文只用于顺序请求，不可在并发求解间共享。
- STRUMPACK 的分析对象还保留值相关匹配/缩放。跨点复用不等同于只缓存
  METIS 排序，必须实测数值轨迹与物理门限；这也是不直接更改默认的原因。

## 对照协议

使用同一冻结 UCRT64 Release (`-O3 -DNDEBUG`、无 `-pg`) 程序，
STRUMPACK 无压缩、METIS、最大对角乘积匹配/缩放、OpenMP=2、OpenBLAS=1、
活跃嵌套层数=1，关闭昂贵因子统计；每个有效点服务核查实际线程设置。

R11 D5 等温 300 K、无 IALMob、Auger 不含生成；原网格 10,241 节点、
19,782 三角形。输入及参考通过
`reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_auger_no_generation_inputs.json`
的哈希合同冻结，Vg=4/8 V 分别从合格零漏压种子检查至前 8 个精确点。
物理单位沿用原配置，电流以 A/μm 输出。计时包含控制器、点服务、状态读写
及原有门限检查，不含生成零压种子的栅压预偏置。

三种模式分别为普通子进程、常驻 worker、常驻 worker＋线性上下文复用。
三者使用相同初值、最大推进步长、外推和失败恢复策略。每档串行执行，
两档轮换顺序；记录系统负载，不宣称独占机器。

逐点重检原块残差、逐行闭合和 KCL 门限；候选对普通子进程的全部精确点
ψ/fn/fp 最大绝对差要求 ≤1e-8 V。保留所有失败和恢复成本。
8 点通过不替代双栅压完整曲线联合评分或完整曲线重复计时。

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build --preset windows-ucrt64-release --target vela_example_runner test_linear_solver test_newton_solver test_dc_sweep --parallel 2
python tests/regression/test_dc_worker.py --runner build-release/vela_example_runner.exe
python scripts/run_templates_ldmos_analysis_reuse.py --output reference_staging/templates_ldmos_cross_point_analysis_20260919
python scripts/analyze_templates_ldmos_analysis_reuse.py reference_staging/templates_ldmos_cross_point_analysis_20260919
```

构建需已实际启用 STRUMPACK 和 `VELA_ENABLE_OPENBLAS_THREAD_CONTROL=ON`；
运行脚本检查 Release 标志，冻结源码、程序、DLL 和构建清单，禁止覆盖证据。

## 结果

所有曲线均完成 Vd=0–9.333333 V 前 8 个精确参考点，各 68 次点服务、
0 次外层减步回退。初始化、内部恢复、未接受候选等成本均保留；
0 次外层回退不表示内部没有恢复。402 个有线性求解的点服务通过运行时线程核验。

| Vg | 模式 | 端到端墙钟 / s | 进程树 CPU / s | Newton 更新 | 全局分析 | 数值分解 |
|---|---|---:|---:|---:|---:|---:|
| 4 V | 普通子进程 | 220.233 | 224.844 | 369 | 125 | 446 |
| 4 V | 常驻 worker | 203.199 | 209.375 | 373 | 129 | 454 |
| 4 V | worker＋主系统分析复用 | **185.789** | **196.062** | 370 | **56** | 448 |
| 8 V | 普通子进程 | 217.830 | 224.797 | 350 | 122 | 427 |
| 8 V | 常驻 worker | 202.955 | 207.641 | 351 | 122 | 428 |
| 8 V | worker＋主系统分析复用 | **185.138** | **194.156** | 349 | **51** | 426 |

相对普通子进程，组合收益为 15.64%/15.01%；相对普通 worker，新增跨点
线性复用收益为 8.57%/8.78%。前者包括进程和输入准备复用，不能全部归于符号分析。
进程树 CPU 为控制器 CPU 与点服务 CPU 的合计，后者沿用 worker 请求前后的
CPU 差分；不含实验冻结、独立后审计及未计入请求的 worker 关闭尾部成本。

| Vg | 模式 | 分析 / s | 数值分解 / s | 回代 / s | 线性总计 / s |
|---|---|---:|---:|---:|---:|
| 4 V | 常驻 worker | 17.749 | 32.317 | 8.182 | 58.662 |
| 4 V | ＋分析复用 | **2.891** | 33.586 | 7.714 | **44.609** |
| 8 V | 常驻 worker | 17.741 | 30.508 | 7.957 | 56.617 |
| 8 V | ＋分析复用 | **2.703** | 32.625 | 7.285 | **43.021** |

数值分解并未变快，主要收益是减少分析准备。Newton 次数只有少量变化，
不应认领非线性算法加速。STRUMPACK 的匹配/缩放复用和有限精度轨迹均可能
影响工作量，本轮不把差异唯一归因于某一个内部机制。

各模式推进目标序列相同，但逐次更新数不完全相同。候选对普通子进程的精确点
ψ/fn/fp 最大差为 4.44e-15/3.55e-15 V；对普通 worker 为 7.99e-15/3.55e-15 V，
均远小于原 1e-8 V 状态一致性门限。所有原块残差、逐行闭合、KCL 检查通过。
这不是完整 D5 双栅压曲线评分，也不继承为 D4/D0 或其他后端的新资格。

## 分析次数如何解释

每档候选均记录一次 `dc.linear_context.misses`、67 次 `hits`，说明同一主上下文
贯穿 68 次点请求。结构核对仍逐次进行；值变化只触发数值分解。

两档各有 54 个有效点服务，其全局 `linear.solve_calls` 恰等于
`newton.linear_l2_row_column.calls`，即没有额外辅助线性求解。这 54 个服务中，
普通子进程和普通 worker 各累计 54 次分析，复用候选均为 **0 次**。
另有一个无需线性求解的零压检查，不计入这 54 个有效服务。

剩余 13 个含辅助系统的服务累计记录候选的 56/51 次全局分析，包含首次主分析。
已确认 QF 限幅后的 Poisson 子块校正、载流子行恢复拥有独立局部求解器。
目前计数器尚未逐算子拆分混合服务，不能断言这 56/51 次全属于某一辅助算子，
也不把全局计数称作“主 Jacobian 分析次数”。

因此，“固定主矩阵只保留一个求解上下文、复用其符号分析”已经落地；
“所有辅助系统也只分析一次”不属于当前已完成的实现。

## 验证和限制

- Release 受影响目标构建成功。线性求解器 **21 测试/541 断言**、
  Newton **123 测试/2,683 断言**、DC sweep **105 测试/3,707 断言**全部通过。
- worker 协议 **2 项**通过，覆盖重复请求的分析复用、强制新数值分解、
  非法请求后的缓存清理、准备输入变化后的失效、与独立进程结果一致。
- linked D5 和后端统计 Python 回归 **23 项**通过，含恢复时禁止改变缓存策略。
  本轮未重跑完整 CTest；没有改动物理公式或验收门限。
- 每配置只有一轮，机器非独占；系统平均忙碌率约 43.9%–48.1%，
  估计后台 CPU 为各曲线 129.9–193.9 CPU·s，负载不完全相同。
  墙钟与 CPU 均有收益，但不能由此认领重复稳定性或同条件优于 UMFPACK。
- 求解器进程峰值工作集：Vg4 普通子进程/worker/候选为 273.5/292.9/264.5 MiB；
  Vg8 为 271.1/287.9/265.4 MiB。这是各求解进程生命周期高水位的最大值，
  不是同时刻整个进程树内存，也不是精确因子内存。

下一步可先给主系统、Poisson 校正和密度恢复分别记录分析次数，
再为固定 Poisson 子块建立独立缓存；保持原门限做单项对照。
当前主系统候选需要完整曲线及交错重复配对后，才决定是否纳入推荐配置。

证据：
[实验汇总](../../reference_staging/templates_ldmos_cross_point_analysis_20260919/summary.json)、
[独立聚合](../../reference_staging/templates_ldmos_cross_point_analysis_20260919/analysis.json)、
[冻结源码/程序清单](../../reference_staging/templates_ldmos_cross_point_analysis_20260919/binary/manifest.json)。
