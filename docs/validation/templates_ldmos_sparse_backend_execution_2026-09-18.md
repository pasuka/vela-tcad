# LDMOS 多直接法后端与 METIS 执行记录

状态：本轮接入及筛选完成（9 月 19 日续接结束）。已完成线性合同测试、两类固定
矩阵筛选、六配置单线程 D5 双栅压前 8 点，以及线程矩阵/STRUMPACK 曲线对照。
新增配置未胜过同批 UMFPACK，因此未触发“优胜候选完整曲线三轮配对”阶段。
遵循[执行计划](templates_ldmos_sparse_backend_plan_2026-09-18.md)，默认仍为 SparseLU。
不改变物理模型、原门限或独立电热求解器配置。

本轮曲线共 20 条、160 个精确点、7,201 次 Newton 更新、0 次回退；
最大精确点电势/QF 差 7.99361e-15 V。计时有背景负载，不能认领独占环境性能。
采用冻结 R11 D5 输入（不启用 IALMob）、300 K，10,241 节点/19,782 三角形，
`unit_scaling` 的微米坐标与 cm^-3 密度口径，端口电流 A/μm。
输入、材料、网格及原始零压种子见
[冻结输入合同](../../reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_auger_no_generation_inputs.json)。

## 已实现的范围

- 可选 `sparselu_metis`、`umfpack_metis`、`mumps`、`mumps_metis`、
  `superlu_mt`、`superlu_mt_metis`、`strumpack`，由共享 `LinearSolver` 管理。
- 保留结构相同的分析复用、值完全相同的数值因子复用；不同 RHS 使用保留因子。
  改值重新分解，改结构重建供应商实例；错误后清除旧缓存，无静默换后端。
- MUMPS 采用非 MPI 的 OpenMP 双精度库、SYM=0；SuperLU_MT 采用低层接口，
  保留符号结构并允许重新选行主元；STRUMPACK 为 METIS、匹配/缩放、无压缩 DIRECT，
  关闭微小主元替换及 GPU。CSC→COO/CSR 和排列恢复均由适配层处理。
- 新后端的 OpenMP 编译选项仅作用于三个适配文件，未传播到物理/Eigen 控制代码。
  正式构建为 UCRT64 Release `-O3 -DNDEBUG`，gprof/LTO 均关闭。
- METIS 使用结构并集图及固定种子；不修改原数值矩阵为对称矩阵。
  UMFPACK 显式检查所请求 METIS 排序确实生效。
- `VELA_LINEAR_THREADS=1/2/4` 仅控制新线程后端；BLAS 线程独立控制。
  因子统计开关扩展至新后端，主计时关闭详细统计。

## 错误处理发现

SuperLU_MT 的 METIS 变体在一个含空行/空列的奇异测试上于库内部发生访问违例；
GDB 定位到库的数值分解调用链。适配层现在在调用供应商代码前拒绝空数值行/列。
另测非空但线性相关的奇异矩阵、非有限输入、失败后重算及改值后的主元变化，均通过。
这不等于证明第三方库对所有病态输入都不会崩溃；保留独立进程超时/失败记录。

## 已完成测试

- 最终 C++ 线性测试：20 项、489 个断言通过（包括七个新增配置及空系统拒绝）。
- 2/4 线程各通过后端合同测试，分别 196 个断言。
- 最终 linked/统计/线性回放 Python 回归：24 项通过（统计模块共 8 项），
  包含新输入坐标、续接归档及串行对照线程固定检查。
  两个已注册 linked/统计 CTest 均通过；未宣称全 CTest。
- PATH 仅保留 Windows System32 时，构建目录中的后端合同测试仍通过，
  证明本次执行文件可使用已复制 DLL 启动，不依赖 UCRT64 bin 进入 PATH。
- 当前原生动态库能力来自本机包及实际调用，不根据仅启用 CMake 开关认领。
- 独立小系统、主元交换、孤立图节点、数值因子复用、结构失效、严重行缩放、
  清缓存及失败恢复纳入测试。尚未宣称全 CTest 或全部物理曲线通过。

## 历史矩阵筛选

原 14 个 30,723 阶、265,211 条目的矩阵，每配置 3 轮、每矩阵 2 个 RHS。
raw/scaled 后向误差 <1e-12、两 RHS 一致性 <1e-10，全通过。
以下为实际 `linear.total` 中位秒数；各小组有各自 UMFPACK 对照，不能跨组当作同时测量。

| 小组 | 配置 | 线性服务中位秒 |
|---|---|---:|
| METIS | SparseLU / UMFPACK | 3.247 / 2.108 |
| METIS | SparseLU+METIS / UMFPACK+METIS | 4.343 / 1.925 |
| MUMPS | UMFPACK / MUMPS / MUMPS+METIS | 2.129 / 3.697 / 3.731 |
| SuperLU_MT | UMFPACK / SuperLU_MT / SuperLU_MT+METIS | 2.711 / 3.025 / 4.490 |
| STRUMPACK | UMFPACK / STRUMPACK | 2.066 / 8.974 |

组内/组间存在后台负载波动，未取得独占空闲环境。
这些数字不等于完整曲线加速，也不用于重新分解此前 UMFPACK 加速的全部原因。

## 当前 R11 代表系统

从已资格的 R11 D5 原尝试输入重跑 Vg=4/8 V，各取 Vd=
0.34621190996608003、4、5.333333333333333、20、40 V。
低压取实际 R11 轨迹点，未重新制造旧序列的 0.375 V。
10 个重跑终态相对对应 R11 终态最大电势/QF 差为 **0 V**。

得到 36 个 `VELALU02` 文件：保存实际进入线性后端的矩阵、RHS 和解，
属于求解器输入坐标，不包含未缩放物理 Jacobian。
故新矩阵仅认领该坐标下的后向误差，原始物理量通过点/曲线原门限检查；
不把同一矩阵重复计算的误差伪装成独立 raw/scaled 两项验证。

9 配置 × 36 系统 × 3 轮，共 972 次矩阵回放、1944 次 RHS 求解全部通过。
最差求解器输入范数型后向误差约 1.25e-17（门限 1e-12）。
单线程、BLAS 单线程、详细统计关闭；表内总服务包含转换、分析、分解和回代。

| 配置 | 线性服务中位秒 |
|---|---:|
| UMFPACK+METIS | 4.617 |
| UMFPACK | 5.978 |
| SuperLU_MT | 7.301 |
| SparseLU | 8.320 |
| MUMPS | 9.303 |
| MUMPS+METIS | 9.602 |
| SparseLU+METIS | 10.888 |
| SuperLU_MT+METIS | 13.625 |
| STRUMPACK | 22.187 |

MUMPS 默认与强制 METIS 的 36 项误差/参考步比较记录相同。
先将 SparseLU、UMFPACK、UMFPACK+METIS、MUMPS、SuperLU_MT、STRUMPACK
列为前 8 点验证配置；其余三个排序控制保留矩阵证据，不重复投入曲线。
是否进入完整曲线阶段，以前 8 点的实际端到端结果及原验收门限决定。

## 9 月 19 日续接

9 月 18 日启动的串行前 8 点任务于次日凌晨中断；检查时相关进程已退出。
Vg=4 V 六配置及 Vg=8 V SparseLU 已完成，中断发生于 Vg=8 V UMFPACK。
续接保留七条完成曲线，校验冻结源码、程序、DLL 和矩阵输入哈希，
将未完成目录及日志另存为 `.interrupted_20260919_160119`，从原零压种子
重跑该条曲线；中断前的部分耗时不计入重跑成绩。原部分日志中的绝对状态路径
仍保留当时文本，查阅归档文件须以该归档目录为根，不用于正式状态资格。

本批次跨两个时间段，均记录背景负载；不能据此认领独占空闲条件下的重复性能。
续接入口新增 `--resume`，须使用原后端列表、轮数、物理配置、线程及统计设置。
后续完整性能配对应在同批次重新测量对照与候选。

## D5 双栅压前 8 点

12 条曲线、96 个精确点全部完成，0 次回退；10 组候选/基准对照通过。
各精确点最大电势/电子 QF/空穴 QF 差为 7.10543e-15 V，门限为 1e-8 V。
范围是 Vd=0–9.333333 V；不转移为新增后端的完整 0–40 V 或 D4 资格。

| 配置 | Vg4 墙钟/s | Vg8 墙钟/s | Vg4/Vg8 Newton 更新 |
|---|---:|---:|---:|
| SparseLU | 264.251 | 276.387 | 374 / 358 |
| UMFPACK | 210.581 | 232.647 | 369 / 349 |
| UMFPACK+METIS | 219.221 | 233.396 | 364 / 349 |
| MUMPS | 255.676 | 268.203 | 364 / 349 |
| SuperLU_MT | 254.702 | 252.340 | 375 / 357 |
| STRUMPACK | 329.760 | 330.077 | 368 / 349 |

本轮是单次曲线筛选，不能据此认领重复稳定性；但没有新增单线程配置在任一
栅压下显示优于 UMFPACK 的端到端耗时。固定矩阵中 UMFPACK+METIS 的收益
没有转化为本轮曲线收益。此后补充的线程对照如下，也未产生优胜候选。

曲线使用 9 月 18 日冻结程序。9 月 19 日最后构建仅补入空系统拒绝、依赖复制及
测试配置清理；非空矩阵算子未改动。线程筛选独立冻结最后构建，不混合二进制计时。

## 1/2/4 线程固定矩阵筛选

每组仍为当前 36 系统、每系统两个 RHS、三轮交错串行；BLAS 固定 1。
表中为 `linear.total` 中位秒数。线程组间背景负载不同，优劣最终仍由同批曲线决定。

| 后端 | 1 线程 | 2 线程 | 4 线程 |
|---|---:|---:|---:|
| UMFPACK 对照 | 5.341 | 6.062 | 6.659 |
| MUMPS | 9.650 | 10.721 | 11.722 |
| SuperLU_MT | 7.767 | 不合格 | 不合格 |
| STRUMPACK | 22.866 | 5.503 | 5.046 |

此处 UMFPACK 是同组对照，适配器未增加 UMFPACK 并行功能；矩阵脚本曾按组统一
设置 OMP 环境，不能把上述三列解释为已证明 UMFPACK 的并行缩放。
后续曲线驱动已显式将 UMFPACK/SparseLU 的 OMP 和后端线程固定为 1，
仅对 MUMPS/SuperLU_MT/STRUMPACK 采用请求的线程数。

SuperLU_MT 2 线程第二轮、4 线程第一轮出现求解器输入后向误差超过 1e-12，
筛选停止该配置的后续轮次。保留失败日志；加入更详细错误文本后的 4 线程
额外四轮未复现，不能覆盖原失败或认领稳定性，原因尚未闭合。
STRUMPACK 各线程三轮均通过；2/4 线程的进程 CPU 中位数约为 8.58/10.30 s，
先选择 CPU 较低的 2 线程进入 D5 双栅压前 8 点，不把矩阵速度直接当作曲线资格。

独立无依赖宏构建还检查了 7 个不可用后端/排列入口，均明确抛出异常。

## STRUMPACK 2/4 线程曲线对照与最终筛选

两组分别重新运行 UMFPACK 单线程基准及 STRUMPACK 候选，均为 D5 双栅压
前 8 点、BLAS=1、无详细因子统计。8 条新增曲线及所有状态差检查通过。
下表“增加”以各自同批 UMFPACK 为基准；未用前一批较有利的时间替代。

| STRUMPACK 线程 | 栅压 | UMFPACK/s | STRUMPACK/s | 墙钟增加 | CPU 增加 | Newton（基准/候选） |
|---|---|---:|---:|---:|---:|---:|
| 2 | 4 V | 228.295 | 240.341 | 5.28% | 16.72% | 369 / 369 |
| 2 | 8 V | 222.541 | 233.618 | 4.98% | 17.62% | 349 / 350 |
| 4 | 4 V | 235.407 | 245.026 | 4.09% | 25.32% | 369 / 371 |
| 4 | 8 V | 226.384 | 243.496 | 7.56% | 25.19% | 349 / 350 |

三批 UMFPACK 对照在每档栅压下的精确状态哈希、目标/更新历史完全一致。
这是控制组复核，不替代新增配置的长曲线重复性资格。
求解子进程峰值工作集最大值：UMFPACK 约 250.45 MiB，STRUMPACK 2/4 线程
约 274.69/276.78 MiB。单线程 SparseLU、UMFPACK+METIS、MUMPS、SuperLU_MT
分别约 278.48、241.09、254.05、286.28 MiB。
这些是含装配和输入的进程峰值，不是 LU 因子独占内存，也不是同时运行的进程树总和。

最终保留 UMFPACK 作为已有性能较优的可选基准，生产默认仍为 SparseLU。
新增后端保留显式实验入口；SuperLU_MT 2/4 线程不得作为已通过精度筛选的配置。
本轮没有新增优胜配置，故按原分阶段筛选计划，不投入 D5/D4 完整曲线三轮配对。
这不是已完成长曲线验收，也不转移 G3/D4/D0 资格。
后续若修复 SuperLU_MT 多线程失败或找到新的端到端收益，须重新通过相应前置门限。

## 证据及续接入口

`reference_staging/` 下：

- `templates_ldmos_metis_screen_20260918`
- `templates_ldmos_mumps_screen_20260918`
- `templates_ldmos_superlu_mt_screen_20260918`
- `templates_ldmos_strumpack_screen_20260918`
- `templates_ldmos_current_linear_captures_20260918_r2`
- `templates_ldmos_current_matrix_screen_20260918`
- `templates_ldmos_multi_backend_first8_20260918`
- `templates_ldmos_threads_{1,2,4}_20260919`
- `templates_ldmos_superlu_mt_failure_detail_20260919`
- `templates_ldmos_superlu_mt_failure_repeat_20260919`
- `templates_ldmos_strumpack_t2_first8_20260919`
- `templates_ldmos_strumpack_t4_first8_20260919`
- `templates_ldmos_sparse_backend_campaign_20260919.json`：汇总及输入摘要文件哈希

最初 `templates_ldmos_current_linear_captures_20260918` 因找不到旧 0.375 V 目标
在运行仿真前终止；保留失败记录，r2 改用真实 R11 低压点。
每次矩阵筛选冻结了工具、源码和 DLL，并核对输入及冻结文件哈希。
当前矩阵的完整分组结果见其 `multi_backend_analysis.json`。

可复用入口：

- `scripts/run_templates_ldmos_matrix_backends.py`：旧/新矩阵、后端、线程数及串行重复。
- `scripts/capture_templates_ldmos_current_linear_inputs.py`：当前 R11 固定尝试诊断采集。
- `scripts/run_templates_ldmos_isothermal_backends.py`：显式后端列表、线程数、8/31 点与轮数。
- `scripts/analyze_templates_ldmos_multi_backend.py`：支持带下划线的多后端名字。

后续不得把新增后端线性测试直接等同于 D5/D4/G3/D0 全曲线资格。
