# R9：UCRT64 Release 前 8 个精确点 gprof 热点复核

日期：2026-09-17。分支 `codex/templates-ldmos-phase-a`。
此前 R9 推荐配置、局部求根及限幅实验已整理提交为 `96bd4b1`；本轮使用
该提交的求解器源码，新增可复现剖析入口，没有修改物理模型、求解策略或门限。

## 冻结范围与数值检查

- 显式 `d0_production_r9_contact.json`，D0 自热，Vg=8 V，原生几何及 SI 物性。
  网格 10241 节点、19782 三角形，其中 5723 个活跃硅节点。
- 原 31 点序列仅截取前 8 项：0、1.3333333333333333、2.6666666666666665、
  4、5.333333333333333、6.666666666666667、8、9.333333333333334 V。
  保留原自适应中间推进、预测器、接触一致性、逐行及块门限。
- 普通 Release 和完整 gprof 均从中性 300 K 独立初始化并预偏置至 Vg=8 V。
  第三次 gprof 从普通 Release 的最终初始化状态开始，只剖析漏压扫描。
  该状态是漏压推进前的合格零压状态，不使用目标点或终点状态。
- 原物性准备、IALMob 屏蔽/局部准备、残差值专用路径和符号分析复用保持开启；
  新中性根缓存、标量 Newton 和局部 QF 裁剪均关闭。

环境：MSYS2 UCRT64 GCC 16.2.0、CMake 4.4.3、GNU gprof 2.47.20260726，
实际启用 UMFPACK 6.3.8；两种构建均 `-O3 -DNDEBUG`、LTO 关闭。
插桩版另加 `-pg -fno-omit-frame-pointer`，链接关闭 dynamicbase，PE 核验通过。
`OMP_NUM_THREADS=1`、`OPENBLAS_NUM_THREADS=1`。三次运行串行，构建和测试在
计时开始前结束；未在 VM 上执行本轮任务。

三组均完成全部 8 个精确点，原终态门限通过。各组漏压尝试均为 13 次，
含零压检查、12 次非零推进，无失败重试。插桩完整运行的初始化及漏压全部
非计时输出与普通 Release **逐值相等**；仅漏压剖析的全部非计时输出亦
逐值相等，四块终态最大差均为零。没有仅凭电流近似相等判定轨迹一致。

## 普通 Release 成本：用于判断整体性能

| 指标 | 初始化/栅压预偏置 | 漏压前 8 点 |
|---|---:|---:|
| 点求解调用 | 27 | 13 |
| Newton 更新 / 数值分解 | 120 / 120 | 137 / 137 |
| 线搜索候选检查 | 120 | 229 |
| 装配 / 其中仅残差 | 147 / 0 | 352 / 201 |
| 符号分析 | 25 | 13 |
| 装配计时 (s) | 4.02 | 72.99 |
| 分解计时，含符号分析 (s) | 3.62 | 27.72 |
| 回代计时 (s) | 0.69 | 4.42 |
| 实际中性根求解 | 8099 | 17794 |
| 中性根密度评估 | 591227 | 1298962 |
| 迭代计数器覆盖的 IALMob 屏蔽极小值求解 | 0 | 4006100 |

IALMob 计数器从点准备后开始取差，表中不包含物理对象构造时的 300 K
屏蔽预计算；中性根计数为该点 assembler 的整个生命周期计数，两者范围不同。

端到端普通 Release 墙钟 **166.93 s**；漏压各次尝试的墙钟合计 **118.97 s**。
以此漏压墙钟为分母：装配 **61.35%**、分解 **23.30%**、回代 **3.72%**，
其余 **13.84 s / 11.64%**。其余包含点准备、门限与结果处理等未被上述
三项包住的工作，不能全部归为 JSON 或 I/O。

端到端减去漏压尝试为 **47.96 s**，包含初始化、栅压预偏置及外层整理；
其中只有 8.34 s 被上述初始化装配/分解/回代计时覆盖。源码显示每次点服务
都会重读网格、构造边映射和物理对象、准备几何，这些发生在迭代阶段计时前。
本轮尚未把剩余时间逐项独立计时，不能把差额全部认领为几何开销。

插桩完整运行 179.79 s；仅漏压插桩运行 135.74 s。它们用于定位热点，
不用于宣称求解加速；本轮仅一组普通 Release 计时，未做性能重复性验收。

## gprof 自身时间采样

下表为可归因采样时间的百分比，**不是端到端墙钟或完整 CPU 百分比**。
完整运行可归因样本 42.42 s，仅漏压 33.46 s，采样单位 0.01 s。
函数简称对应原始 `gprof_flat.txt` 中完整符号；只比较 self time，不把调用图
的父子累计时间相加。

| 函数/路径 | 完整运行样本占比 | 仅漏压样本占比 |
|---|---:|---:|
| `ElectrothermalAssembler::assemble` 主体 | 8.77% | 10.43% |
| IALMob 单元 `evaluateImpl` | 7.36% | 10.34% |
| IALMob `evaluatePrepared<Dual>` | 5.49% | 7.77% |
| IALMob `screeningMinimum` | 5.02% | 6.49% |
| IALMob `evaluate` | 4.27% | 5.26% |
| `evaluateIalHighFieldMobility` | 3.11% | 4.96% |
| IALMob `evaluateWithDerivatives` | 3.32% | 4.27% |
| `updateIalTransportState` | 2.50% | 3.53% |
| `buildIalInterfaceGeometry` | 8.30% | 3.44% |
| `IalScreeningCache::minimum` | 2.90% | 3.32% |
| `thermalSgCurrent` | 2.48% | 2.99% |
| 插桩 `_mcount_private` + `__fentry__` | 11.88% | 8.40% |

按互不重叠的 self 样本合计，IALMob 单元、低/高场、准备、屏蔽缓存及状态
更新约占完整运行 **36.28%**、漏压 **48.98%** 的可归因样本；不含几何准备，
也不含插桩函数。分组明细在 `analysis.json`。这说明需要优先定位迁移率
状态构造和求值链条，而非继续仅盯住中性接触求根。

### 剖析限制与交叉核验

1. `libumfpack.dll`、BLAS 和运行库没有随本轮一起插桩，gprof 不能完整展开
   它们的内部成本；输出中的 Eigen/UMFPACK 适配器条目不等于库内分解时间。
   采用普通 Release 的 27.72 s 内部阶段计时判断线性分解的重要性。
2. MinGW/GCC 本轮调用图存在明显的调用次数归属异常：完整运行的 `assemble`
   显示 8985058 次，而程序内部计数为 499 次；仅漏压显示 8933869 次，内部
   为 352 次。因此不将这些调用次数、递归环或传播的 children 时间作为结论。
   原始调用图保留供排查；符号表和入口反汇编也已保存，未猜测异常的唯一原因。
3. 自身采样受内联、优化和有限采样影响，函数排序只是定位线索；源代码和
   内部计数用于交叉验证，改动收益仍需无插桩的独立配对。GNU 文档也明确
   children 时间依赖调用成本均匀假设，且剖析不能反映未运行期间的等待：
   [实现说明](https://sourceware.org/binutils/docs/gprof/Implementation.html)、
   [flat profile](https://sourceware.org/binutils/docs/gprof/Flat-Profile.html)、
   [children 估计限制](https://sourceware.org/binutils/docs/gprof/Assumptions.html)。

## 对应源路径与下一步优先级

1. **先验证跨点复用不可变几何和输入准备。**
   [点服务](../../src/simulation/ElectrothermalSimulation.cpp)每次重新读网格、
   构造边映射和 assembler；[IALMob 几何](../../src/equation/IalTransport.cpp)
   在每个新对象上重新准备；[界面搜索](../../src/physics/IalInterfaceGeometry.cpp)
   为每个节点遍历界面线段和候选节点。当前复用作用于单次点求解内部，未跨
   40 次点服务共享这些不可变数据。建议先加准备阶段计时，再引入扫压上下文
   保存网格/映射/几何；缓存键须包含网格、单位、晶向、接触和几何来源，
   不跨温度复用旧载流子状态。预期收益尚未测量。
2. **细分 IALMob 单元求值成本，验证残差路径和温度方向的重复工作。**
   [单元核](../../src/equation/IalElementMobility.cpp)的电学/温度方向分别进入
   `evaluateImpl`，虽已复用低场局部偏导，单元场构造和高场求值仍值得核查。
   残差值专用开关已经开启，不应再次提议只开启它；可比较专用标量单元路径
   是否比当前零导数 Dual 表示便宜。屏蔽极小值仍有 4006100 次实际求解，
   但已存在预测保护和精确缓存，须测量缓存失效率和有效温度变化后决定，
   不把它与本轮关闭的中性根 Newton 混为同一算法。
3. **继续定位装配次数与主体 Newton，而非重复做已完成的符号分析复用。**
   漏压每次更新平均 2.57 次装配；151 次完整 Jacobian 装配、201 次仅残差。
   229 次候选检查对应 137 次接受更新，仍有 92 次未接受候选的成本。
   0.375 V 为 13 次更新/43 次候选，4 V 为 18/30，5.333333 V 为 15/24，
   是后续固定目标诊断的优先样本。数值分解已占漏压墙钟 23.30%，减少更新
   可同时减少装配和分解；当前 13 次符号分析对应 137 次分解，不能把二者
   混淆。此前失败的更新映射/恢复候选不因本轮热点结果自动晋级。

上述是本轮分析后的实现建议；本轮没有实施这些新优化。旧 SparseLU 热点
结论也不能直接搬到当前 UMFPACK、电热及 IALMob 配置。

## 复现与证据

使用[生产导出器](../../scripts/export_templates_ldmos_production.py)及本轮新增
[剖析入口](../../scripts/run_templates_ldmos_gprof_first8.py)。后者校验导出清单、
冻结程序和 DLL、逐组保存输入/输出/gmon、校验原门限及跨构建状态，并生成
flat/callgraph 文本和 CSV。需要为 `--output` 选择新目录。

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --preset windows-ucrt64-release -B build-gprof-ldmos-first8 -DVELA_ENABLE_GPROF=ON -DVELA_ENABLE_UMFPACK=ON -DVELA_ENABLE_LTO=OFF
cmake --build build-gprof-ldmos-first8 --parallel 2 --target vela_example_runner
python -X utf8 scripts/run_templates_ldmos_gprof_first8.py --bundle <R9-export-bundle> --release build-release/vela_example_runner.exe --instrumented build-gprof-ldmos-first8/vela_example_runner.exe --gprof D:/msys64/ucrt64/bin/gprof.exe --output <new-run-directory>
```

本轮原始证据：`reference_staging/templates_ldmos_gprof_r9_20260917/`，
包括 `configure.log`、`build.log`、`build_manifest.json`、`pe_imports.txt`、
`tooling_tests.log`、`analysis.json`、`run_matrix.log` 及 `runs/summary.json`。
完整输入和结果留在 ignored 证据目录，不纳入 Git。

| 文件 | SHA256 |
|---|---|
| 普通 Release runner | `4b3fd5816d75a2e63a0d89292bb0d8ce962c07253e410883584414b08dc1c78f` |
| gprof runner | `2a7bf7981f6f6fabedc3714898636e440c21770029fa79e9e7f9ca43547f240d` |
| 完整运行 `gmon.out` | `7894d27efc8a26caeb175b31274918f0508e92298024a5b05b9ee3ef122e9b85` |
| 仅漏压 `gmon.out` | `17a6f88180319606a25d99b8d19386fca68bcf28f02d655d912b58c89dddd329` |

本轮复用了上一轮同源码普通 Release **829/829** 的完整回归，确认 Ninja
无需重建；没有为分析脚本重复全量求解器测试。新运行入口完成上述三组真实
集成运行，复用的导出及 gprof 解析工具 **15 项回归通过**。本轮只认领 Vg8
前 8 点剖析，不认领新的双栅压完整曲线或与 Sentaurus 的重复性能配对。
