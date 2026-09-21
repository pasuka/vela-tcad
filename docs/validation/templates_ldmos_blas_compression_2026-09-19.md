# SparseLU BLAS、UMFPACK 线程与 STRUMPACK 压缩验证

## 范围与入口

本轮为用户授权的独立加速实验，生产 `LinearSolver` 和物理模型未改。
两个 `EXCLUDE_FROM_ALL` 诊断程序使用同一源码，分别关闭/开启
`EIGEN_USE_BLAS`，均使用 `EIGEN_DONT_PARALLELIZE`，不链接 `vela_core`，
避免不同 Eigen 宏实例混用及全局修改物理装配。
UCRT64 Release：`-O3 -DNDEBUG`，Eigen 5.0.1、OpenBLAS 0.3.34、
UMFPACK 6.3.8、STRUMPACK 8.0.0。

运行入口：

```powershell
cmake --build --preset windows-ucrt64-release --target linear_solver_acceleration_study linear_solver_acceleration_blas_study --parallel 2
python -m unittest tests.regression.test_linear_acceleration_study
python scripts/run_templates_ldmos_acceleration_study.py --output <new-directory> --captures reference_staging/templates_ldmos_current_linear_captures_20260918_r2 --rounds 3 --configs sparselu_control sparselu_blas1 sparselu_blas2 sparselu_blas4 umfpack_blas1 umfpack_blas2 umfpack_blas4 strumpack_t2 blr6_t2 blr8_t2 hss6_t2 hss8_t2 amalg_t2
python scripts/analyze_templates_ldmos_acceleration_study.py <new-directory>
```

矩阵来自此前 R11 D5 双栅压 10 个代表点的 36 个系统；网格为 10,241 节点、
19,782 三角形，300 K、不含 IALMob。原输入坐标 μm、密度 cm^-3、端口电流
A/μm。`VELALU02` 保存的是实际求解器坐标系统，不是未缩放物理 Jacobian。
每个系统检查两个 RHS；同结构复用分析，改值重新分解。记录准备、分析、分解、
总线性服务时间及进程墙钟/CPU。总线性服务时间不含读文件、误差审计；
它不等于完整曲线端到端时间。保留后台负载，不认领独占机器。

门限保持求解器输入 normwise backward error <1e-12，双 RHS 一致性 <1e-10。
记录分块分量误差和方向差；近零方向的相对差可能很大，另计算绝对分块 L2 差。
本轮不转移 D5/D4/D0 曲线资格。

## 线程口径更正

本机 `libopenblas.dll` 查询为 `USE_OPENMP`，`openblas_get_parallel()==2`。
安装目录与 Release 目录该 DLL 的 SHA-256 相同：
`4b5cbee2aba54a599a031023dcc0d99ebd24988954fe094cf59ed76c004781f6`。
旧脚本的 `OPENBLAS_NUM_THREADS=1` 被此构建忽略；设置 `OMP_NUM_THREADS=4`
时库查询为 4。因此旧多线程报告的“BLAS=1”不是已验证的运行事实。
数值与墙钟记录仍保留，不能将其全部收益归因于求解器自身 OpenMP。

新诊断使用 `openblas_set_num_threads` 并查询验证；STRUMPACK 的 OpenMP
线程和 BLAS 线程分别为 2/1。UMFPACK 本身采用串行稀疏流程，其 OpenMP
线程组大小与 BLAS 的 1/2/4 设置一致。SparseLU 也采用该 BLAS 线程组配置，
但显式关闭 Eigen 自身并行。所有配置禁止嵌套活跃线程组。
这些是运行时配置值，不代表每个小型 BLAS 调用一定启动全部线程。

初版诊断把 UMFPACK 的 OpenMP 组限制为 1，却把 BLAS 设为 2/4，出现 120 s
超时和主动终止；该记录保留在 `templates_ldmos_acceleration_screen_20260919`，
不能作为有效多线程性能结论。修正后的实验独立记录，不覆盖初版证据。

## STRUMPACK 候选

保持 METIS、最大对角乘积匹配/缩放、双精度、无 GPU、无微小主元替换。

- `strumpack`：无压缩 DIRECT 基准。
- `blr6/blr8`：BLR 压缩相对容差 1e-6/1e-8。
- `hss6/hss8`：HSS 压缩相对容差 1e-6/1e-8。
- 压缩候选的绝对压缩容差为 0，最小 separator/front 为 128/256、leaf 为 64；
  使用 PREC_GMRES，外层相对容差 1e-13、绝对容差 0、最多 100 次、restart=30。
  压缩容差不是最终求解容差，也不放宽原误差门限。
- `amalg`：无压缩 DIRECT，启用 MUMPS_SYMQAMD 与 aggressive amalgamation。

这是有限参数筛选，不表示穷尽所有压缩阈值或排序组合。

## 三轮结果与结论

13 个配置均完成 3 轮，共 1,404 个矩阵、2,808 个 RHS，全部通过原门限。
最大输入 normwise backward error 为 1.302e-16；求解器坐标分块 L2 方向差
最大约 1.593e-9，不等同于自洽物理状态误差。所有冻结文件哈希复核通过。

| 配置 | 总线性服务中位数/s | 三轮范围/s | 分解中位数/s | 进程 CPU 中位数/s |
|---|---:|---:|---:|---:|
| sparselu_blas1 | 8.412 | 7.390–8.429 | 7.933 | 8.125 |
| sparselu_blas2 | 7.842 | 7.784–8.158 | 7.407 | 7.797 |
| sparselu_blas4 | 8.478 | 8.446–9.041 | 8.026 | 8.594 |
| umfpack_blas1 | 5.094 | 5.073–5.134 | 3.970 | 5.172 |
| umfpack_blas2 | 5.500 | 5.494–5.508 | 4.360 | 5.766 |
| umfpack_blas4 | 6.601 | 6.395–7.569 | 5.362 | 7.031 |
| sparselu_control | 8.126 | 8.012–9.510 | 7.691 | 8.031 |
| strumpack_t2 | 4.731 | 4.665–4.804 | 2.035 | 7.688 |
| blr6_t2 | 9.144 | 9.114–11.143 | 2.109 | 15.391 |
| blr8_t2 | 8.061 | 8.057–8.096 | 2.150 | 13.031 |
| hss6_t2 | 9.971 | 9.397–10.980 | 2.660 | 16.375 |
| hss8_t2 | 9.205 | 8.716–9.317 | 3.080 | 14.344 |
| amalg_t2 | 4.763 | 4.718–4.896 | 2.051 | 7.797 |

**SparseLU**：GDB 在 `SparseLUImpl::column_bmod → dgemv_` 停下，
证明本机 Eigen 5.0.1 的 `EIGEN_USE_BLAS` 确实覆盖了部分稠密块内核。
此前笼统称该路径不能利用 BLAS 不准确。稀疏符号操作和完整消去流程仍未并行。
BLAS=2 的中位数比控制快约 3.49%，但重复区间重叠；BLAS=1/4 的中位数反而略慢。
初轮小幅收益未成为稳定一致收益，保留为独立候选，不修改生产默认。

**UMFPACK**：DLL 导入包含 OpenBLAS 的 `dgemm_/dgemv_/dtrsm_`；
GDB 在 BLAS=4 时捕获 `libumfpack → libopenblas → libgomp!GOMP_parallel`，
证明实际矩阵分解进入了 BLAS 的 OpenMP 并行调用。BLAS=2/4 比单线程中位数
分别慢约 7.97%/29.58%。开启多线程与取得加速是两件事；本样本继续优先单线程。
旧串行曲线 UMFPACK 的 OMP=1 控制不是已开启多线程加速的证据。

**STRUMPACK**：真正限制 BLAS=1 后，无压缩两线程基准为 4.731 s。
本轮 BLR/HSS 候选约 8.061–9.971 s，增加约 70.4%–110.8%；
每 36 个矩阵、两个 RHS 的累计 GMRES 迭代为 140–212 次。
BLR 因子存储计数总和仅降低约 0.8%，HSS 甚至可能增加；压缩准备与修正成本
超过节省的计算。合并候选为 4.763 s，约慢 0.67%，无明确收益。
各轮 STRUMPACK 的因子条目数与部分 Krylov 次数存在变化，未认领逐比特复现。
库的 factor_nonzeros 是存储计数口径；不同压缩结构不应直接解释成传统精确 LU 填充。

本轮停止于固定矩阵筛选：压缩候选无收益，SparseLU BLAS 的小幅差异受负载波动影响，
未晋级生产配置或进行新的长曲线资格转移。无压缩 STRUMPACK 2/1 线程配置
若继续进入曲线，须使用修正后的线程控制重新做前 8 点；不能继承旧线程口径。

构建及测试：两个 Release 诊断目标通过；新增 2 项 Python 测试包含 16 个
后端/编译变体子用例、改值后的分析复用、线程查询、未知模式拒绝和覆盖保护。
既有 24 项相关 Python 回归通过。生产 C++ 核心未修改。

调试调用栈保存在该证据目录 `debugger/`，不计入性能数据。
生产默认、物理模型及验收门限未变。

## 证据

`reference_staging/` 下：

- `templates_ldmos_acceleration_screen_20260919`：初始诊断及线程配置失败。
- `templates_ldmos_acceleration_repeats_20260919`：修正后冻结程序/源码/DLL/输入，
  13 配置交错串行重复，每个配置 3 轮。

参考：[Eigen BLAS 接口](https://eigen.tuxfamily.org/dox/TopicUsingBlasLapack.html)、
[OpenBLAS 线程说明](https://github.com/OpenMathLib/OpenBLAS)、
[STRUMPACK BLR](https://portal.nersc.gov/project/sparse/strumpack/v4.0.0/BLR_Preconditioning.html)。
实际选项以本机 8.0.0 头文件为准；参考网页版本不代替安装包核实。
