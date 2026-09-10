# LDMOS Vg=8 首轮热点分析：0–1.3333 V（2026-09-09）

本轮在合并后的 `d4285815f648be3f0b20719cbbe4d33e6a1b38b8` 上完成 gprof
采集及同版本 Release 对照。两组均通过原数值门限和完整性审计，24 个接受状态
逐字节一致。首轮确认的主要成本是 **SparseLU 数值分解和 Jacobian 装配**，
其次是残差及连续性检查。尚未修改求解算法，也未验证任何优化的加速收益。

## 实验边界与配置

- Vg=8 V，Vd=0–1.33333333333333 V，沿用已合格 0 V 检查点；不包含栅压预偏置。
- 原外部检查点推进流程：初始步长 0.0025 V，上限 0.1 V，增长因子 1.35，
  直接求解后按需同偏压再求解；预测器关闭、原验收门限不变。
- 网格 10241 节点、19782 三角形、30022 边，其中 5723 个硅节点。
  300 K、Fermi–Dirac、OldSlotboom BGN、SRH/Auger、
  constant_field / transport_cell_vector QF 高场迁移率；雪崩、量子、热关闭。
  几何输入 μm、密度 m⁻³、电势 V、电流 A/μm。
- MSYS2 UCRT64 GCC 16.2.0、C++20、Release `-O3 -DNDEBUG`，实际后端
  Eigen SparseLU/COLAMD、L2 行列均衡。HDF5、SPQR、UMFPACK 被检测到，
  本轮没有选用后两者作为求解后端。
- 独立 `build-gprof-ldmos`：`VELA_ENABLE_GPROF=ON`、LTO 关闭，
  `-pg -fno-omit-frame-pointer`，MinGW 链接禁用动态基址。
  两组均开启内部 `performance_profiling`，不把内部计时的代价视作零。
- gprof 程序 SHA256：`354d01560621b18ae995507da9b16b6d6e16ff2d4e6248d51a508dfa58d94440`。
- Release 对照 SHA256：`8f0b7c0b4e38b032bdac32db5c7782f01efa91ac643cb6361b14fc8f31ff6a1b`。

## 完成与等价性

| 指标 | gprof + 内部计时 | Release + 内部计时 |
| --- | ---: | ---: |
| 子进程墙钟合计 | 299.332 s | 267.118 s |
| 子进程 CPU 合计 | 249.313 s | 224.031 s |
| 内部 dc_sweep.total 合计 | 293.019 s | 263.502 s |
| 接受的物理推进 | 23 | 23 |
| 子进程 | 45 | 45 |
| 实际 Newton 更新 | 230 | 230 |
| 回退 | 0 | 0 |
| 接受状态 | 24 | 24 |

gprof 运行时间为 21:51:55–21:57:01，Release 对照为 21:57:23–22:01:57
（Asia/Shanghai）。控制器间隔另含校验和文件处理，不能与子进程合计混用。
本次 gprof 墙钟比对照多 12.06%；单次顺序运行只能说明观察到的开销，
不用于声称相对历史版本或 Sentaurus 的性能改进。

45 个子步骤的目标、返回码、Newton 更新次数逐项一致；24 个接受状态的 SHA256
全部一致，端点各物理场差异为零，电流同为 `0.00012221004372605 A/μm`。
两个精确区间端点的最大归一化 KCL 为 `4.463896416247188e-12%`。
各接受推进还通过原全局块、局部行和 KCL 检查。此为局部验证，不构成 0–40 V 全曲线验收。

## 以 Release 内部计时为准的热点

下表比例分母为 `dc_sweep.total = 263.502 s`。计时包含子调用，存在嵌套，
**不得将所有行相加**；例如分解包含在线性求解中，边上物理和稀疏模式构造包含在
Jacobian 中，线搜索包含残差及部分连续性检查。

| 范围 | 累计时间 | 占内部总时间 | 解释 |
| --- | ---: | ---: | --- |
| linear.total | 77.945 s | 29.58% | 全部线性求解 |
| linear.factorize | 72.633 s | 27.56% | 占线性总时间 93.19% |
| newton.jacobian | 73.344 s | 27.83% | Jacobian 装配及正则项 |
| jacobian.edge_physics | 49.506 s | 18.79% | 约占 Jacobian 时间 67.50% |
| newton.line_search | 32.039 s | 12.16% | 含试算残差，不能再与残差相加 |
| dd.residual | 31.955 s | 12.13% | 756 次残差计算 |
| dd.continuity_diagnostics | 26.290 s | 9.98% | 630 次连续性项检查 |
| jacobian.pattern_build | 13.605 s | 5.16% | 44 次稀疏模式构造，已包含在 Jacobian 中 |
| dc.record_point | 6.671 s | 2.53% | 点记录及其后处理范围 |

`dc.solve_point = 220.948 s`，其余 `42.554 s` 分布在初始化、记录和未单独
细分的流程中；不能将这部分全部归为文件 I/O。`newton.linear_row_scaling`
和 `newton.linear_l2_row_column` 包含线性求解，不代表缩放本身花费约 77 s。

内部计数确认：251 次 Jacobian、288 次数值分解、711 次线搜索试算，
230 次线搜索接受、21 次线搜索失败。这些失败均由既定同偏压再求解流程处理，
没有物理步回退。251 个迭代阶段含最终未接受更新的阶段，不能写成 251 次实际更新。
符号分析实际执行 81 次、缓存命中 207 次；288 次 `linear.analyze` 阶段还包含缓存检查。
符号分析仅 3.217 s，三角求解 1.878 s，当前主要线性成本在数值分解。

## gprof 交叉证据及限制

45 个子进程各自保留 gmon.out，完整合并为 gmon.sum，并导出 flat profile、
call graph 和 CSV。采样合计 113.06 s，不能当成完整 CPU 或墙钟时间。
SparseLU::factorize 的采样 self 占比为 10.45%，调用数 288，与内部计数吻合；
rebuildFixedJacobianPattern 的 self 占比 5.47%，调用数 44，也与内部计数吻合。
二者的 self 比例与内部含子调用的时间比例不是同一口径。

`_mcount_private`、`__fentry__` 合计 14.62% 的采样 self 时间属于采集开销。
此外出现 `newtonConfigFromJson` 约 239 万次调用、restrictedSparseNorm 占比
9.96%、未启用弧长流程的回调管理符号等异常归属。它们不能证明 JSON 解析、
诊断范数或弧长算法是热点；本轮不据此调整代码。优化后的内联函数、局部符号和
MinGW gprof 归属问题仍需地址/源码级采样验证，call graph 的传播时间不作为优化依据。

## 建议的下一轮定位顺序

1. **拆分 Jacobian 边上通量导数和物理量重算。** 当前边阶段内部计数为
   173,689,980 次 Fermi–Dirac half 调用。源代码的 Fermi 路径仍对边端点变量
   做正负扰动通量评估，值得检查状态不变量和重复计算的缓存机会。这个调用总数
   不等于 Fermi 函数独占时间，也不证明可以直接更换导数公式。
2. **拆分 SparseLU 数值分解。** 这是已确认的最大单项阶段；现有符号分析缓存已生效，
   仅优化 analyzePattern 不会覆盖主要成本。后端/排序的对照需要固定输入和原门限，
   本轮尚未运行替代后端。
3. **检查残差与连续性项的同状态重复计算。** 行权重与局部行验收均调用连续性项计算，
   可以先核对重复状态及可共享结果。必须保留局部行验收；不能为加速关闭检查。
4. **量化跨子进程重建成本。** 45 次启动造成重复初始化，44 次模式构造约 13.6 s。
   若改为保留求解器实例，需验证检查点、严格再求解与缓存失效语义，不能将理论节省
   直接写成已实现收益。外推属于另一组减少迭代次数的 A/B，不与本轮热点结论混合。

源代码入口：
[边上物理及通量导数](../../src/equation/CoupledDDAssembler.cpp)、
[行权重与局部行验收](../../src/solver/NewtonSolver.cpp)、
[SparseLU 阶段](../../src/solver/LinearSolver.cpp)。

## 证据与复现

原始输出保留在本地忽略目录，不提交二进制、gmon 或仿真状态：

- [采集适配器](../../reference_staging/templates_ldmos_gprof_20260909/run_profile.py)
- [汇总脚本](../../reference_staging/templates_ldmos_gprof_20260909/summarize_profile.py)
- [两组等价性审计](../../reference_staging/templates_ldmos_gprof_20260909/comparison.json)
- [Release 阶段汇总](../../reference_staging/templates_ldmos_gprof_20260909/release_control/hotspot_summary.json)
- [gprof 阶段和函数汇总](../../reference_staging/templates_ldmos_gprof_20260909/instrumented/hotspot_summary.json)
- [gprof flat profile](../../reference_staging/templates_ldmos_gprof_20260909/instrumented/gprof_flat.txt)
- [gprof call graph](../../reference_staging/templates_ldmos_gprof_20260909/instrumented/gprof_callgraph.txt)

从当前工作树根目录设置 UCRT64 PATH 后，构建命令为：

```powershell
cmake --preset windows-ucrt64-release -B build-gprof-ldmos -DVELA_ENABLE_GPROF=ON -DVELA_ENABLE_LTO=OFF
cmake --build build-gprof-ldmos --target vela_example_runner --parallel 4
```

运行适配器接受 `--runner`、`--output` 和可选 `--gprof`，输出目录必须不存在。
区间与严格控制流程固定在脚本中。两组运行后分别执行 summarize_profile.py，
再运行 compare_runs.py。此次只构建插桩目标并运行两组真实网格局部验证；
未重复与本轮无关的完整 CTest。此前合并版 741/741 CTest 通过的记录保持独立。
