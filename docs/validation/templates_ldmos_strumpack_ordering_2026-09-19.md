# STRUMPACK METIS / AMD / MMD / AND 排序对照

状态：已完成。承接[2＋1 线程曲线验证](templates_ldmos_strumpack_verified_threads_2026-09-19.md)，
四排序×三分析策略×三轮、1,296 个系统/2,592 个 RHS 均通过原线性门限。
充分复用及按代表点复用时 METIS 最快；逐矩阵重建时 AMD 最快。
没有统一优胜替代项，不更改生产配置，不认领新增曲线资格。

## 固定条件与计时边界

- MSYS2 UCRT64 Release，STRUMPACK 8.0.0；独立诊断程序，不修改生产适配器。
- 无压缩 DIRECT、最大对角乘积匹配/缩放、关闭 tiny-pivot 替换和 GPU；
  OpenMP=2、OpenBLAS=1、活跃嵌套上限=1。通过 API 设置并在每个系统结束核验。
- 四配置依次为 METIS、AMD、MMD、AND，三轮交错串行；不启用额外超节点合并。
- 沿用 `templates_ldmos_current_linear_captures_20260918_r2` 的 36 个 R11 D5
  求解器输入系统，覆盖双栅压既有代表点。输入是 VELALU02 缩放后的矩阵及右端项，
  不能称为原始物理残差矩阵。每个系统检查两组 RHS，沿用后向误差 <1e-12、
  两组 RHS 线性一致性 <1e-10 门限，不放宽任何物理验收标准。
- 三种分析策略：按相邻矩阵结构相同复用；每个矩阵强制重新分析；
  每个代表点目录重建、点内复用。36 个输入结构完全相同，第一种只分析一次。
  逐矩阵重建用于隔离反复准备成本，不代表真实曲线所有点都需要重建。
  重新分析同时刷新值相关的匹配/缩放，可能改变排列与填充，因此三种策略之间
  不能简单用总耗时相减当作纯排序成本；排序优劣以各策略内部的配对为准。
- 通过 STRUMPACK 的虚拟 `perf_counters_start/stop` 回调记录
  `nested dissection`（含排列应用）和 `symbolic factorization`。保留原库回调行为，不重写算法。
  对 AMD/MMD，库仍沿用 `nested dissection` 这个计时标签，它表示所选排序阶段，
  不意味着 AMD/MMD 本身变成嵌套剖分算法。
- `reorder_other_seconds` 是完整 `reorder()` 减去上述两个阶段的余量，包含
  匹配、平衡化、结构对称化及其他未细分开销，不能当作精确的匹配耗时。
- 分析总计还含求解器对象准备、矩阵设置；结构复用时该字段包含矩阵值更新。
  CSC→CSR 转换另记准备时间。线性总计含准备、分析、分解、两次回代及因子条目查询，
  排除捕获文件读取与误差检查；进程墙钟另录。各嵌套阶段不可重复相加。

## 复现

```powershell
cmake --build --preset windows-ucrt64-release --target linear_solver_acceleration_study linear_solver_acceleration_blas_study --parallel 2
python -m unittest tests.regression.test_linear_acceleration_study tests.regression.test_linear_strategy_study tests.regression.test_templates_ldmos_backend_statistics
python scripts/run_templates_ldmos_acceleration_study.py --output reference_staging/templates_ldmos_strumpack_ordering_reuse_20260919 --captures reference_staging/templates_ldmos_current_linear_captures_20260918_r2 --rounds 3 --timeout 180 --configs strumpack_t2 strumpack_amd_t2 strumpack_mmd_t2 strumpack_and_t2
python scripts/run_templates_ldmos_acceleration_study.py --output reference_staging/templates_ldmos_strumpack_ordering_cold_20260919 --captures reference_staging/templates_ldmos_current_linear_captures_20260918_r2 --rounds 3 --timeout 180 --configs strumpack_t2 strumpack_amd_t2 strumpack_mmd_t2 strumpack_and_t2 --cold-analysis
python scripts/run_templates_ldmos_acceleration_study.py --output reference_staging/templates_ldmos_strumpack_ordering_point_20260919 --captures reference_staging/templates_ldmos_current_linear_captures_20260918_r2 --rounds 3 --timeout 180 --configs strumpack_t2 strumpack_amd_t2 strumpack_mmd_t2 strumpack_and_t2 --reset-per-directory
```

两个 Release 诊断目标构建成功；15 项相关 Python 回归通过，包括新增四排序×
两分析策略、目录边界的重建/复用、改值后的复用和已知非对称解核验。
生产核心未改，不宣称重跑完整 CTest。
源文件、输入、可执行文件和 DLL 哈希随实验冻结；后台系统负载单独记录。

每点组是在首两组完成后增加目录重建开关并重新构建的诊断版本；组内四种排序
使用同一冻结程序，首两组原始冻结证据保留。没有修改数值求解步骤或验收门限。
本轮不自动转移曲线资格，也不由单一分析策略的矩阵胜出直接启动长曲线晋级。

## 总成本结果

单位为秒，每格是 36 个输入依次执行的线性服务总计的三轮中位数，包含两次 RHS。
以下不是整条曲线墙钟，也没有运行新的 Newton 轨迹。

| 分析策略 | 每轮实际分析次数 | METIS | AMD | MMD | AND |
|---|---:|---:|---:|---:|---:|
| 相同结构跨输入复用 | 1 | **5.3261** | 7.7458 | 8.1407 | 7.3721 |
| 每个代表点重建、点内复用 | 10 | **6.7416** | 7.6359 | 8.7467 | 7.8083 |
| 每个矩阵重建 | 36 | 11.4057 | **8.6289** | 10.1785 | 9.3650 |

三轮的线性总计范围：

| 分析策略 | METIS | AMD | MMD | AND |
|---|---|---|---|---|
| 跨输入复用 | 4.652–6.850 | 7.227–8.363 | 7.809–9.493 | 7.075–8.318 |
| 按代表点复用 | 6.689–7.282 | 7.219–10.128 | 8.507–10.236 | 7.727–8.440 |
| 逐矩阵重建 | 11.278–12.404 | 8.501–9.355 | 10.010–10.403 | 9.312–10.161 |

METIS 在前两种策略的每一轮都快于其他三种排序，AMD 在逐矩阵重建的每一轮最快。
按代表点复用时 AMD 比 METIS 增加约 13.27% 总成本；逐矩阵重建时减少约 24.35%。
各策略的整机平均忙碌率样本范围依次为 41.18%–79.83%、40.65%–70.98%、
41.35%–57.63%，包含本任务负载；不能把绝对时间波动全部归因于算法。

## 排序、树构造及填充

下表采用按代表点复用组。时间仍为一轮 36 个输入的三轮中位数；填充比为
`factor_nonzeros()/input_nonzeros` 的样本中位数，仅在相同库和定义下比较。

| 项目 | METIS | AMD | MMD | AND |
|---|---:|---:|---:|---:|
| 排序与排列应用/s | 1.8716 | 0.3431 | 0.4001 | 0.4361 |
| 符号树构造/s | 0.0965 | 0.1061 | 0.1386 | 0.1062 |
| reorder 其他开销/s | 0.2779 | 0.2934 | 0.2865 | 0.2732 |
| 分析/改值服务总计/s | 3.0743 | 1.5889 | 1.6656 | 1.6161 |
| 数值分解/s | 2.1319 | 3.1832 | 3.6885 | 5.0350 |
| 两次 RHS 求解/s | 1.4005 | 2.7064 | 3.2392 | 1.0380 |
| 因子填充比 | 7.9749 | 9.4978 | 9.4507 | 15.4112 |

各列取中位数的轮次可能不同，且计时有嵌套，因此不能把表格各行直接相加。
AMD 明显减少排序时间，但其填充比增加，分解与回代成本更高；AND 回代较快，
但填充比接近 METIS 的两倍，分解更慢。符号树构造本身只占本组 METIS 分析
总计的一小部分，排序与排列应用才是已测子阶段中的主要部分。

逐矩阵重建组进一步验证这一点：METIS 排序/树构造分别为 6.2413/0.3001 s，
AMD 为 1.1555/0.3344 s。该组每个矩阵均重新分析，AMD 节省的准备时间
足以抵消分解与回代增加的成本；这不能直接转移到已有分析复用的实际扫压流程。

## 准确性、重复性与限制

- 三组全部线性检查通过，最坏求解器输入后向误差为 **1.9120e-17**，
  两组 RHS 的线性一致性差为 0；所有系统结束时线程核验通过，无失败或超时。
- 36 个输入来自 10 个代表点目录，每点捕获数并不完全相同；按点重建组实际
  执行 10 次分析，与跨输入组的 1 次、逐矩阵组的 36 次分别核对。
- 这批捕获没有原始未缩放矩阵。相对保存参考方向的最大分块 L2 差约
  1.5991e-9，属于求解器坐标，不能替代物理状态差或完整曲线验收。
- 额外重复性检查发现 **METIS 的因子条目数跨轮并非逐值一致**。
  三组每轮因子条目累计量的范围分别为 76.042–76.875、76.070–76.483、
  76.267–76.362 百万（顺序：跨输入、按点、逐矩阵）。逐矩阵组 AMD/MMD
  也有少量条目差异；AND 三组一致。具体差异保存在 `verification.json`。
  尚未定位这些差异来自哪一步，不将其全部归因于 METIS 随机种子或线程调度，
  不宣称三轮因子结构完全确定性。
- 每组 70 项冻结哈希检查通过，共 210 项检查；性能程序已正常退出。
  本轮验证的是线性算法与成本，不是新增物理、Newton 更新数或端到端曲线资格。

## 结论及证据

当前保留 METIS，优先研究跨点上下文与分析复用。AMD 保留为频繁重建场景的
独立候选；若要改变曲线推荐配置，仍需用相同线程和原门限做实际前八点配对，
同时补充排序/匹配及因子结构的重复性诊断。当前不将 MMD/AND 晋级。

三个 `reference_staging/templates_ldmos_strumpack_ordering_{reuse,point,cold}_20260919`
目录分别保存 `summary.json`、逐轮结果、`analysis.json`、`verification.json`、
冻结 `source/` 与 `binary/`。逐轮原始数据保留负载、进程 CPU/墙钟和所有失败记录。
可用 `python scripts/analyze_templates_ldmos_acceleration_study.py <证据目录>` 重建汇总。
