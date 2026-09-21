# LDMOS 边通量导数缓存实施与验证（2026-09-10）

已实现体区迁移率复用、固定几何缓存和边端点状态缓存。Release 全部 742 项测试通过；
Vg=8 V、Vd=0–1.33333333333333 V 的优化前后对照中，24 个接受状态逐字节一致，
两版均为 230 次 Newton 更新。边物理阶段 FermiHalf 调用减少 44.99%，该阶段
累计墙钟减少 30.28%，整个 Jacobian 阶段减少 16.66%。

本次顺序整段对照的总墙钟由 400.783 s 增至 410.621 s（增加 2.45%），CPU 时间
由 275.328 s 降至 266.141 s（减少 3.34%）。因此本报告确认消除了重复计算并保持
该区间数值结果，不将其描述为已确认的整段墙钟加速。

## 实现与数值边界

- **体区迁移率复用**：vector-QF 模式下，全局启用接触电场回退时，对没有实际接触
  回退且没有表面迁移率的体区边复用基态迁移率。接触边及其相邻第三节点的 psi
  扰动仍重新求驱动场与迁移率；原解析高场反馈保留。已给定迁移率时跳过原先无用的
  驱动场求值及表面状态拷贝。
- **固定几何缓存**：装配器初始化时记录接触输运单元掩码、邻接输运面积及 stencil
  节点的面积加权几何敏感度，供电子和空穴共用。每次 Jacobian 装配分别计算一次
  电子/空穴单元 QF 梯度，供基态迁移率和解析反馈复用。实际 QF 梯度依赖当前状态，
  每次装配更新；本次并未消除实际单元梯度求值中的全部几何运算。
- **端点状态缓存**：缓存仅存在于一次 Jacobian 装配的当前边内，电子/空穴及两端
  分开存储。以 psi、局部 QF、边参考系 QF 的 double 位模式作为键，最多保存
  8 组端点状态；换边即清空，容量不足时正常重新求值。保留原参考值/增量精度路径，
  尤其电子的 `electronDensityAt` 与 `Nc*FermiHalf(eta)` 分别缓存，未合并不同
  舍入顺序的密度表达式。有限差分步长和稳定通量公式不变。

代码在 `CoupledDDAssembler.h/.cpp`，新增 `jacobian.endpoint_cache_hits/misses`
聚合计数。新增 Catch2 测试覆盖 Fermi/vector-HFS 在体区和接触回退边上的
Jacobian 与独立残差差分对照，并依次改变、恢复状态，核对缓存生命周期。
未修改 SparseLU 求解器生命周期、外推设置、物理模型或收敛门限。

## 配置与复现证据

- 基础提交 `d4285815f648be3f0b20719cbbe4d33e6a1b38b8`，本轮源码为该提交上的
  工作区变更；不能仅用 HEAD 代表优化版。三份源码 SHA256 和完整差异见下方 manifest/source.diff。
- 优化前程序 SHA256：`8f0b7c0b4e38b032bdac32db5c7782f01efa91ac643cb6361b14fc8f31ff6a1b`。
- 优化后程序 SHA256：`06301e6d44359dffa8889262524401b499c8a619cdb6b2bcb462d8c2f0734b6e`。
- MSYS2 UCRT64 GCC 16.2.0、C++20、Release `-O3 -DNDEBUG`；两版均为内部
  `performance_profiling` 开启、gprof 关闭。实际线性后端 Eigen SparseLU/COLAMD，
  L2 行列均衡；本构建检测到 HDF5、SPQR、UMFPACK，但本算例使用 SparseLU。
- 网格为原 exact topology：10241 节点、19782 三角形、30022 边，5723 个硅节点；
  external average-box 输运系数、barycentric 节点体积。输入几何 μm，密度 m⁻³，
  电势 V，报告电流 A/μm。网格、材料、掺杂、输运系数及种子均冻结并记录哈希。
- 300 K、Fermi–Dirac、OldSlotboom BGN、SRH/Auger、constant_field /
  transport_cell_vector QF 高场迁移率；接触电场回退开启，表面迁移率、雪崩、量子、热关闭。
- 沿用合格 Vd=0 V 种子，不含栅压预偏置。初始步长 0.0025 V、上限 0.1 V、
  增长因子 1.35；原检查点控制器、按需同偏压再求解、预测器关闭。
  原块门限 psi/electron/hole 为 5e-8/1e-11/3e-10，局部行和 KCL 比值门限均 1e-8。

## 完成与等价性

| 指标 | 优化前 | 优化后 |
| --- | ---: | ---: |
| 子进程墙钟合计 | 400.783 s | 410.621 s |
| 子进程 CPU 合计 | 275.328 s | 266.141 s |
| 接受推进 / 接受状态 | 23 / 24 | 23 / 24 |
| 子进程 / Newton 更新 / 回退 | 45 / 230 / 0 | 45 / 230 / 0 |
| Jacobian 次数 / 数值分解次数 | 251 / 288 | 251 / 288 |
| 边物理 FermiHalf 调用 | 173,689,980 | 95,551,170 |
| 端点缓存命中 / 未命中 | 未实施 | 52,092,540 / 63,698,780 |

45 个子步骤的目标、阶段、返回码、Newton 更新数和父状态哈希逐项一致；两个
原完整性审计均通过。24 个接受状态 SHA256 全部一致，端点 psi/phin/phip/n/p
最大差异均为零，电流同为 `0.00012221004372605 A/μm`。
直接求解的部分非零返回码由原严格同偏压再求解流程处理，未引入新的物理步回退。

整段基线时间为 08:58:46–09:05:42，优化版为 09:05:42–09:12:48（Asia/Shanghai）；
此控制器时间还包含校验和文件处理，不能与上表子进程合计混用。

## 分阶段耗时及解释

下表为内部累计墙钟，包含子调用，各行存在嵌套，不能相加。

| 阶段 | 优化前 (s) | 优化后 (s) | 变化 |
| --- | ---: | ---: | ---: |
| dc_sweep.total | 393.296 | 403.689 | +2.64% |
| newton.jacobian | 108.074 | 90.067 | −16.66% |
| jacobian.edge_physics | 71.446 | 49.811 | −30.28% |
| jacobian.pattern_build | 22.184 | 23.783 | +7.21% |
| linear.factorize | 110.746 | 121.641 | +9.84% |
| dd.residual | 43.521 | 47.380 | +8.87% |
| dd.continuity_diagnostics | 36.091 | 41.575 | +15.19% |

FermiHalf 少执行 78,138,810 次，是本轮消除端点重复求值的直接计数证据，
并非该函数独占耗时。缓存命中率为 44.99%。三类优化作为一个组合验证，未做
分别关闭的消融实验，因此不能把全部阶段节省单独归给某一种缓存。

未修改的数值分解、残差和连续性检查阶段也变慢，说明顺序整段计时中存在干扰或
其他运行时变化；这不足以单独证明变化来源。不能用昨晚较快的历史墙钟替代本轮
基线，也不能由本轮 CPU 降幅推算 0–40 V 加速比例。

### 同检查点交错短测

为检查顺序干扰，再固定 `child_00005_direct` 的同一父状态和原配置，求解到
Vd=0.01043125 V。顺序为 ABBA/BAAB（A=基线，B=优化版），每版 4 次，
只替换输出路径与程序。8 次全部收敛，每次 5 次 Newton 更新，状态均与整段
基线该检查点逐字节一致。

| 中位数指标 | 基线 (s) | 优化版 (s) |
| --- | ---: | ---: |
| 整次子进程墙钟 | 10.015 | 11.062 |
| newton.jacobian | 2.963 | 2.598 |
| jacobian.edge_physics | 2.018 | 1.439 |
| linear.factorize | 3.306 | 3.642 |
| dd.residual | 0.520 | 0.684 |

边物理中位耗时减少 28.72%，与整段局部结果方向一致。墙钟范围分别为
8.601–11.658 s 和 8.376–14.190 s；整次墙钟中位数仍增加 10.45%。短测进一步
支持局部优化有效，但没有消除整体耗时的不确定性，且不能代替整段验收或归因实验。

## 检查与证据索引

- `windows-ucrt64-release` 构建通过；Release 全 CTest **742/742 通过**，无跳过，
  总耗时 154.44 s，含新增 Fermi/vector-HFS 测试和既有亚 ULP 参考系及数值回归。
- `git diff --check` 通过。源码变化后的最终重编译和全套测试完成后，才冻结优化版程序。
- 本轮只验证上述低压区间，不构成优化版 Vg=8、0–40 V 全曲线验收。

本地生成物位于忽略目录，未加入版本管理：

- [运行清单与源码哈希](../../reference_staging/templates_ldmos_edge_cache_20260910/manifest.json)
- [冻结源码差异](../../reference_staging/templates_ldmos_edge_cache_20260910/source.diff)
- [整段比较结果](../../reference_staging/templates_ldmos_edge_cache_20260910/comparison.json)
- [顺序运行脚本](../../reference_staging/templates_ldmos_edge_cache_20260910/run_pair.py)
- [比较脚本](../../reference_staging/templates_ldmos_edge_cache_20260910/compare.py)
- [基线内部计时](../../reference_staging/templates_ldmos_edge_cache_20260910/baseline/hotspot_summary.json)
- [优化版内部计时](../../reference_staging/templates_ldmos_edge_cache_20260910/candidate/hotspot_summary.json)
- [交错短测脚本](../../reference_staging/templates_ldmos_edge_cache_20260910/repeat_checkpoint.py)
- [交错短测汇总](../../reference_staging/templates_ldmos_edge_cache_20260910/repeat_summary.json)
- [交错短测逐次记录](../../reference_staging/templates_ldmos_edge_cache_20260910/repeat_records.json)
- [CTest 日志](../../reference_staging/templates_ldmos_edge_cache_20260910/ctest.log)
- [前序重复计算分析](templates_ldmos_stage4_edge_derivative_analysis_2026-09-10.md)
