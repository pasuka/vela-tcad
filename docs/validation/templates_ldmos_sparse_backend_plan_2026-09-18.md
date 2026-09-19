# LDMOS 多稀疏直接法后端与 METIS 排序对比计划

后续接入及筛选已完成；结果与阶段适用范围以[执行记录](templates_ldmos_sparse_backend_execution_2026-09-18.md)为准。
以下保留制定计划时的安装核实和阶段合同。

状态：2026-09-18 已核实本机 UCRT64 安装并完成独立小系统调用检查。
本文是接入及性能验证计划；尚未接入新生产后端，尚未运行新后端 LDMOS 曲线。
保留当前 SparseLU 默认及已验证 UMFPACK 配置，不修改物理模型、Newton 策略和原验收门限。
用户希望先比较实际性能，排序/填充/主元策略的深入归因可后置。

## 1. 已核实的本机环境

工具链根目录 `D:/msys64/ucrt64`；包版本来自本机 `pacman -Q`，
不是根据网上最新版本推测。已检查头文件、导入库、DLL、配置文件及主要 DLL 依赖。

| 组件 | 已安装包版本 | 本机形态与检查结论 |
|---|---|---|
| MUMPS | 5.9.1-1 | `mumps-dso` 顺序、`mumps-dto` OpenMP、`mumps-dmo` MPI+OpenMP；前两者小系统通过 |
| SuperLU_MT | 4.0.2-1 | `libsuperlu_mt_OPENMP.dll`；单线程小系统通过 |
| STRUMPACK | 8.0.0-1 | OpenMP 开启，MPI/CUDA/HIP/SYCL 未启用；无压缩直接解及 METIS 小系统通过 |
| METIS | 5.1.0-7 | `idx_t` 为 32 位；`METIS_NodeND` 与正逆排列检查通过 |
| SuperLU | 7.0.1-1 | 串行包另行安装；本轮仅检查安装/头文件/DLL 依赖，未做求解验证 |
| SuperLU_DIST | 9.2.1-1 | 已安装；分布式路线暂不列入第一批接入 |
| OpenBLAS | 0.3.34-1 | 上述已检查的三个新求解器 DLL 均依赖该库 |
| SCOTCH / ParMETIS | 7.0.12-3 / 4.0.3-5 | 已安装；不等于每个后端均编译启用了这些功能 |
| MS-MPI / ScaLAPACK | 10.1.1-23 / 2.2.3-1 | 已安装包；本轮不认领 MPI 运行资格 |

MUMPS 的三个双精度库使用同一 API 名称，首选 OpenMP 版作为新后端，
线程数设为 1 时建立基线；顺序版保留作独立控制，不在同一适配模块混链两版。
本机 MUMPS、SuperLU_MT 和 METIS 小程序实测索引类型均为 4 字节。
不得将 64 位索引缓冲区强制转换为这些接口的 32 位数组。

需要特别记录的配置差异：

- STRUMPACK `StrumpackConfig.h` 中 `STRUMPACK_USE_METIS` 显示未定义，
  但其 CMake 导出依赖 `METIS::metis`，DLL 实际导入 `libmetis.dll`。
  本轮显式选择 `STRUMPACK_METIS` 后 reorder/factor/solve 均返回成功，
  查询排序枚举为 1。应以编译链接加实际调用检查判定能力，不能仅检查该宏。
- 串行 SuperLU 的配置头显示未启用 `HAVE_METIS`，其直接 DLL 依赖也无 METIS。
  可用外部 METIS 生成排列再传入用户排列接口；不假设内置 METIS 可用。
- SuperLU 和 SuperLU_MT 共享部分全局 C 符号名称，但类型/实现不能混用。
  后续须检查同进程链接的实际符号归属；必要时使用分别链接的适配 DLL。
  独立翻译单元隔离头文件本身不能解决链接符号冲突。

### 本轮实际执行的检查

独立 C++20 小程序以 `-O2 -DNDEBUG` 编译；这只是 API/依赖检查，不是性能测量。
环境：`OMP_NUM_THREADS=1`、`OPENBLAS_NUM_THREADS=1`、`OMP_DYNAMIC=FALSE`。

求解同一非对称系统：

```text
A = [[4,1,0], [2,5,1], [0,3,6]]
b = [6,15,24], x_expected = [1,2,3]
```

| 检查 | 结果 |
|---|---|
| METIS，5 节点路径图 | 排序返回成功，逆排列一致 |
| MUMPS 顺序版，SYM=0、ICNTL(7)=5 | 分析/分解/回代成功，INFOG(7)=5，最大解误差 4.44089e-16 |
| MUMPS OpenMP 版，同上 | 同上 |
| SuperLU_MT，nprocs=1 | info=0，最大解误差 0；该小测试使用库的 AᵀA 最小度排序，未测试 METIS |
| STRUMPACK，NONE 压缩、DIRECT、METIS | 三阶段返回 0，最大解误差 4.44089e-16 |

证据与可重跑脚本位于忽略的构建目录：
`build-release/backend_inventory_20260918/{probe.py,results.json,packages.txt}`。
`results.json` 保存编译参数、返回码、输出及可执行文件 SHA256。
这些检查不包含大型病态矩阵、符号复用、失败恢复、多线程或物理曲线资格。

## 2. 后端接入设计

现有共享入口为 `LinearSolver`，已支持 SparseLU/UMFPACK 的结构和数值缓存。
初期保持调用端 `solve(A,b)` 和 `clearPatternCache()` 不变，后端状态用独立 RAII
对象管理，避免在公开头文件混入多种供应商类型。以下名字均为拟议配置，尚未实现。

| 后端 | 拟议选择名 | 输入与阶段接口 | 第一批范围 |
|---|---|---|---|
| MUMPS | `mumps` | CSC→一基 COO；`dmumps_c`，job=1/2/3，SYM=0，job=-2 释放 | 本机 OpenMP 双精度，不启动 MPI |
| SuperLU_MT | `superlu_mt` | CSC；使用可分离分析、分解、回代的接口 | OpenMP 双精度；不能每次走完整初始化而声称复用 |
| STRUMPACK | `strumpack` | CSC→CSR；set/update matrix、reorder、factor、solve | 无压缩、直接求解；记录 matching、缩放及微小主元处理 |

新增可选 CMake 检测；仅在头文件/链接和运行检查成功后报告可用能力。
MUMPS 优先 `pkg-config mumps-dto`；STRUMPACK 使用 `STRUMPACK::strumpack`；
SuperLU_MT 通过显式查找头文件和 `superlu_mt_OPENMP` 导入库检测。
显式请求未编译的后端应报错，不静默回退。

缓存合同：结构相同且排序/后端参数相同才复用分析；仅 RHS 改变可复用数值因子；
值变化必须重新分解；维度/结构/排序/关键控制变化必须失效；失败后不得复用旧因子。
COO/CSR 的索引转换和元素位置映射按结构缓存，每次更新数值，计入实际线性服务成本。
重用结构不意味着冻结旧的数值主元；MUMPS matching、SuperLU_MT refact 等选项
须按本机版本语义验证，不能强制复用失效的行排列。

## 3. METIS 作为独立排序变量

METIS 优化的是未知量/方程的消去次序以减少 LU 填充，不是按非零数值大小排序。
它可与不同直接求解器组合，不能被当作第四个数值分解后端。

对非对称 Jacobian，先从结构构造无向图：

```text
G = pattern(A) union pattern(transpose(A))
去对角、自环及重复邻接边；保留孤立节点；固定编号、种子与 METIS 选项
(perm, iperm) = METIS_NodeND(G)
检查双射、方向和正逆一致性；按求解器接口转换排列
```

这里是结构并集，不对 A 做数值对称化，也不求解 AᵀA 正规方程。
首先仅以 METIS 输出替换列排列 Q，解 `A Q y=b` 后恢复 `x=Q y`；
若另试双边排列，必须同时变换 RHS 并恢复未知量，作为独立命名配置。
不把仅列排列和双边排列的数据混为一个实验。

安排两类对照：

1. **用户实际使用对照**：每个后端的稳健默认排序，对比其 METIS 选项。
   MUMPS `ICNTL(7)=5`；STRUMPACK `STRUMPACK_METIS`；UMFPACK 检查
   `UMFPACK_ORDERING_METIS` 实际生效；SparseLU 和 SuperLU_MT 使用适配的外部排列。
2. **固定排列对照**：相同外部图/种子生成同一 Q，分别传给可接受用户排列的后端。
   若某后端还做 matching、后排序或内部重排，记录实际排列，不能声称排序完全相同。

先在已验证 SparseLU/UMFPACK 上试 METIS，以较小改动判断排序收益。
排序、图构造成本均计入分析/总服务；固定图重用可单独列出摊销值。
记录填充和内存用于筛选；不要求本阶段解释各算法的所有速度差异。

## 4. 串行执行与验证阶段

| 阶段 | 内容 | 通过条件/交付物 |
|---|---|---|
| P0 已完成 | 安装、依赖、5 个独立调用检查 | 本文及本地检查证据；不代表生产集成完成 |
| P1 | METIS 图/排列适配；SparseLU/UMFPACK 固定矩阵对照 | 排列与还原正确、原系统误差通过；记录排序及总服务成本 |
| P2 | 依次接入 MUMPS → SuperLU_MT → STRUMPACK | 每接入一个就完成小系统、缓存、错误和固定矩阵检查，再做下一个 |
| P3 | 合格后端的 D5 双栅压前 8 个精确点 | 原门限及精确点状态差≤1e-8 V；统计全部尝试，不只成功点 |
| P4 | 有收益候选进行 D5/D4 双栅压完整曲线重复配对 | 各 31 点；每配置至少 3 轮，交错顺序、串行运行，原联合门限 |
| P5 | 优胜后端 1/2/4 线程及适用范围复核 | 单线程基准独立保留；共享路径复核 G3，若推广至电热再验证 D0 |

固定矩阵集先复用已有 14 个 30,723 阶、265,211 条目的历史系统作为烟雾测试，
再补取当前冻结 R11 的低/中/高压与困难点矩阵。历史组通过不代替当前曲线资格。
固定系统测量至少 3 轮、交错后端次序；记录失败/超时，不剔除不利样本。
候选若填充过大或无法满足误差门限，应在固定矩阵阶段止损，不直接投入长曲线。

单元/集成检查包括非对称且需主元交换的已知解、病态缩放、奇异系统失败、
多 RHS、同结构改值、结构变化、清缓存、失败后再次求解及后端不可用报错。
固定矩阵沿用 raw/scaled 范数型后向误差 <1e-12、两 RHS 一致性 <1e-10；
不能要求所有不同后端 Newton 轨迹逐位一致，但须比较轨迹变化、失败恢复及总成本。
每个后端自身重复性单独检查，物理/热学及状态门限不放宽。

## 5. 公平计时与决策

- 正式测量使用同一 UCRT64 Release `-O3 -DNDEBUG`，冻结源码/二进制/依赖哈希，
  单线程先行，保留 SparseLU 和 UMFPACK 两个对照；不复用不同时间段的旧墙钟作唯一基准。
- 单进程串行运行测试；先完成编译/测试再计时。记录负载和电源配置；
  严格空闲条件若不可获得，标记为有背景负载数据，不自动改门限后仍声称空闲。
- 因子详细统计与主要计时分开；`VELA_LINEAR_FACTOR_STATISTICS=0` 的语义扩展至新后端。
  保留必要阶段时间/计数和错误码，禁止在主计时循环遍历所有因子提取 FLOPS。
- 分列准备/格式转换、排序/分析、数值分解、回代/改进、总线性服务及端到端墙钟；
  同时记录 CPU 时间、Newton 更新、装配/候选/回退次数、峰值进程内存和填充。
  库内部内存估计不与进程 RSS 混用；缺失计数标记不可用，不能填 0。
- 多线程先设求解器线程为 1/2/4、BLAS 固定 1，避免嵌套过度订阅；
  BLAS 多线程另做命名实验。线程数由后端 API 与运行环境共同设置并记录实际配置。
- 端到端计时沿用已完成 D5 对照的边界：包含零压种子重闭合、中间推进、输出及
  内部审计，不含历史栅压预偏置和外层联合评分；内存、失败和物理结果同样影响选择。
- 不预设新库一定优于 UMFPACK；依据总耗时、稳定性、资源消耗筛选推荐候选。
  不自动把小矩阵最快、分解最快或多线程最快配置改成生产默认。

串行 SuperLU 可作为第二批低成本参考；SuperLU_DIST、MPI 多进程、压缩 LU、
混合精度和 GPU 暂不加入首轮组合，避免将算法近似或额外通信开销混入直接法基线。

## 6. 依据

- [前轮 D5 低统计开销重复配对](templates_ldmos_umfpack_low_noise_2026-09-18.md)。
- [METIS 官方实现及定位](https://github.com/KarypisLab/METIS)：填充减少排序。
- [MUMPS 官方 FAQ](https://www.mumps-solver.org/index.php?page=faq)：METIS 选择及非对称模式。
- [MUMPS MSYS2 构建配方](https://github.com/msys2/MINGW-packages/blob/master/mingw-w64-mumps/PKGBUILD)：
  具体安装形态仍以本机 `.pc` 与 DLL 为准。
- [SuperLU_MT 官方实现](https://github.com/xiaoyeli/superlu_mt)。
- [STRUMPACK 官方稀疏直接解接口](https://portal.nersc.gov/project/sparse/strumpack/master/sparse.html)：
  网上文档版本可能滞后于本机 8.0.0，函数签名及能力以安装头文件和调用检查为准。
