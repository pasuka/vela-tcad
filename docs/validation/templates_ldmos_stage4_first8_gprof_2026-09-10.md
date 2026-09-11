# LDMOS Vg=8：前 8 个精确参考点 Release gprof 热点分析

2026-09-10，本轮完成 Vd=0–9.33333333333333 V 的 8 个精确参考点。
**使用 Release `-O3 -DNDEBUG`，主要热点为 SparseLU 数值分解、Jacobian 装配及残差计算。**
108 个接受状态与此前同源码 Release 运行逐字节一致，原适用数值门限全部通过。
本轮增加采集、解析和审计证据，没有修改求解器代码或物理模型。

## 构建与实验范围

- 工作树：`.worktrees/templates-ldmos-phase-a`，HEAD 为
  `d8ff26804a59735ff3f446846e446e0a62add914`，包含上一轮尚未提交的连续性项缓存修改。
  不能仅凭 HEAD 复现；121 个源码/构建输入均已冻结，并与此前生产 Release 快照核对一致。
- 独立目录 `build-gprof-ldmos-first8`，MSYS2 UCRT64 GCC 16.2.0、C++20、
  **Release `-O3 -DNDEBUG -pg -fno-omit-frame-pointer`**，LTO 关闭，
  MinGW 链接参数 `-Wl,--disable-dynamicbase`。实际求解后端为
  Eigen SparseLU/COLAMD、L2 行列均衡。构建检测到 HDF5、SPQR、UMFPACK，
  不表示本次选用了 SPQR/UMFPACK。
- 插桩程序 SHA256：`f7970e9d0aabea552872cded5243191c5764af3ac7e4f24a2b26abce44659368`。
- Vg=8 V，300 K；10241 节点、19782 三角形、30022 边，其中 5723 个硅节点。
  Fermi–Dirac、OldSlotboom BGN、SRH/Auger、constant_field / transport_cell_vector
  QF 高场迁移率，接触电场回退为 contact_node_cell；表面迁移率、雪崩、量子、热关闭。
  外部 average-box 输运系数、barycentric 控制体积。几何 μm、密度 m⁻³、电势 V、电流 A/μm。
- 从已合格的 0 V 检查点重新初始化；不计栅压预偏置。
  精确参考点依次为 **0、1.33333333333333、2.66666666666667、4、
  5.33333333333333、6.66666666666667、8、9.33333333333333 V**。
  内部仍保留初始步长 0.0025 V、最大步长 0.1 V、增长因子 1.35、预测器关闭、
  必要的同偏压严格再求解/密度恢复，故 8 个参考点不等于只启动 8 次求解。
- 原块门限 `5e-8 / 1e-11 / 3e-10`、局部行及 KCL 比率 `1e-8` 不变。
  没有运行 40 V 全曲线或 26.6667 V 换帧验证。

## 完成、耗时与数值验证

运行时间为 14:00:03–14:25:30（Asia/Shanghai，控制器及审计范围）。
214 个子进程均单独保存非空 `gmon.out`，共 453645272 字节，文件逐一计算 SHA256。
0 V 初始化后先验证 gprof 可读，再推进后续参考点。

| 指标 | 本次 Release + gprof | 此前同源码 Release |
| --- | ---: | ---: |
| 子进程墙钟合计 | 1490.740 s（24.85 分钟） | 1235.173 s |
| 子进程 CPU 合计 | 1191.625 s | 1045.969 s |
| 内部 dc_sweep.total | 1467.969 s | 1219.333 s |
| 精确参考点 / 接受状态 | 8 / 108 | 8 / 108 |
| 物理推进 / 子进程 | 107 / 214 | 107 / 214 |
| Newton 更新 / 物理步回退 | 1284 / 0 | 1284 / 0 |

此前 Release 来自上一轮 AB/BA 实验，控制器时间包含等待另一组程序的时间；
这里只比较其子进程累计时间。两次不属于同期配对的插桩开销实验。
墙钟观察值增加 20.69%、CPU 增加 13.93%，不能全部归因于 gprof，也不是优化收益。

独立审计确认 214 个子进程的目标、父状态、返回码及 Newton 更新次数一致，
108 个接受状态 SHA256 一致，所有内部工作计数及阶段调用次数一致。
存在 106 个非零返回的子过程，随后由原流程恢复；其成本全部计入，不把“零回退”写成“每次直接求解成功”。
相对 Sentaurus 的非零电流误差中位数 1.594856%、P95 1.624118%、
最大值 1.627051%，低压微分电阻误差 1.032119%。8 点最大归一化 KCL 为
`9.721109981968602e-12%`，分母为最大绝对端子电流。工程/最终等级的四项适用门限均通过。
终点电流 `0.000312466896071033 A/μm`；0 V 不参与非零电流相对误差统计。

## 内部阶段计时：优先定位的成本

本次比例分母为 `dc_sweep.total = 1467.969 s`。计时**包含子调用且有嵌套**，
不能把所有行相加；下表同时保留此前同源码非插桩 Release 作为性能背景。

| 阶段 | 本次插桩时间 | 本次占比 | 此前 Release 时间 |
| --- | ---: | ---: | ---: |
| SparseLU 数值分解 `linear.factorize` | 439.597 s | 29.95% | 397.748 s |
| Jacobian 装配 `newton.jacobian` | 355.341 s | 24.21% | 290.518 s |
| 残差计算 `dd.residual` | 253.376 s | 17.26% | 184.408 s |
| Jacobian 内的边上物理 `jacobian.edge_physics` | 213.714 s | 14.56% | 167.753 s |
| 连续性诊断 `dd.continuity_diagnostics` | 81.733 s | 5.57% | 60.836 s |
| Jacobian 内的稀疏模式构造 `jacobian.pattern_build` | 73.135 s | 4.98% | 67.162 s |
| 符号分析/缓存检查 `linear.analyze` | 17.143 s | 1.17% | 16.098 s |

数值分解占线性求解总时间 `468.618 s` 的约 93.81%。1612 次数值分解中，
符号分析实际执行 429 次、缓存命中 1183 次，因此只优化符号分析覆盖不到主要成本。
1390 次 Jacobian 中，边上物理约占 60.14%；213 次模式构造反映了跨子进程重建成本。
线搜索范围为 250.875 s，含残差等子调用，不能再与残差相加。

现有连续性缓存正常生效：命中 1763 次、未命中 1525 次；行权重/验收门限保持不变。
共 4528 次残差调用、4312 次线搜索试算。边上物理计数记录
529147968 次 Fermi–Dirac half 调用，这个计数不是该函数的独占耗时。

## gprof 函数热点与符号归属修正

累计 self 采样为 **595.32 s**，每个样本 0.01 s；这是采样覆盖量，不是完整 CPU 或墙钟时间。
默认 PE 符号解析漏掉部分优化克隆函数，将相邻地址的时间和调用归给其他符号。
例如默认输出把 54.60 s 和 13283500 次调用归给 `newtonConfigFromJson`。

本轮使用同一可执行文件的 `nm -n` 文本函数符号，通过 gprof `-S` 完整解析：
排除节标记，对共享地址的 66 组别名保留清单并确定唯一代表，
**从全部原始 gmon.out 重新合并调用弧**，随后重新导出 flat profile 和 call graph。
不能只对默认合并后的 gmon.sum 换符号表，否则已经合并错位的调用弧无法恢复。
两种解析的 self 总量一致；修正后 JSON 配置解析为 214 次、self 样本为 0，
54.60 s 正确落到 `SparseLUImpl::panel_bmod`。原始解析仍作为证据保留。

| 完整符号表下的函数/内核 | self 时间 | self 占比 |
| --- | ---: | ---: |
| Eigen SparseLU `factorize` | 72.32 s | 12.15% |
| Eigen `gebp_kernel` 矩阵运算内核 | 65.48 s | 11.00% |
| Eigen `SparseLUImpl::panel_bmod` | 54.60 s | 9.17% |
| `rebuildFixedJacobianPattern` | 30.72 s | 5.16% |
| Eigen `gemm_pack_lhs` | 25.00 s | 4.20% |
| `fixedJacobianOffset` | 13.96 s | 2.34% |
| `residualImpl` | 13.84 s | 2.32% |

`factorize` 的 1612 次和模式构造的 213 次与内部计数一致。
Eigen 例程合计占 self 样本约 45.60%，其中也包含 SparseLU 以外的 Eigen 调用，
不能全部等同为数值分解。`_mcount_private` / `__fentry__` 合计 81.36 s、13.67%，
属于采集运行时开销；不能把这个采样比例直接当作整程序减速百分比。
call graph 的传播时间依赖调用模型，本轮以 self 样本和内部计时交叉证据作判断。

## 下一步优化建议

1. **优先分解 SparseLU 的实际成本。** 针对 `panel_bmod`、矩阵内核及填充规模，
   在固定输入矩阵上做数值分解专项计时，再选择内核参数或后端对照。
   保留 COLAMD 为现行基准，避免把已生效的符号分析缓存误当成主要优化空间。
   后端或排序变化必须验证线性残差、Newton 轨迹及原 8 点门限，不能直接宣称收益。
2. **继续减少 Jacobian 边上重复物理量计算。** 当前 Fermi 路径仍对端点变量做正负扰动通量评估。
   检查现有端点缓存未覆盖的状态不变量及 `fixedJacobianOffset` 查找，优先保持原计算顺序和导数语义。
   分别量化边阶段与模式阶段，避免用 Fermi 调用总数推断单函数瓶颈。
3. **量化残差重复工作与跨进程重建。** 4528 次残差仍是第三大阶段；检查相同状态下物理辅助量复用。
   模式构造约占 5%，长生命周期求解器可能减少重建，但须先设计检查点、严格再求解和缓存失效语义。
   不通过关闭局部行/KCL 验收或减少参考点来实现加速。

本轮只完成热点定位，不包含上述求解器修改。后续优化统一以同源码普通 Release
固定检查点对照及这 8 个精确参考点验收，不将 gprof 运行时间作为性能收益依据。

## 证据与复现入口

- [Release 构建参数、程序和源码清单](../../reference_staging/templates_ldmos_first8_gprof_20260910/binary/manifest.json)
- [8 点采集脚本](../../reference_staging/templates_ldmos_first8_gprof_20260910/run_first8.py)
- [完整运行计划](../../reference_staging/templates_ldmos_first8_gprof_20260910/instrumented/plan.json)
- [逐文件采样清单](../../reference_staging/templates_ldmos_first8_gprof_20260910/instrumented/profile_inventory.json)
- [独立数值及工作量审计](../../reference_staging/templates_ldmos_first8_gprof_20260910/comparison.json)
- [内部计时与默认符号解析](../../reference_staging/templates_ldmos_first8_gprof_20260910/instrumented/hotspot_summary.json)
- [完整符号表解析与分区间统计](../../reference_staging/templates_ldmos_first8_gprof_20260910/analysis.json)
- [修正后的 flat profile](../../reference_staging/templates_ldmos_first8_gprof_20260910/external_complete_flat.txt)
- [修正后的 call graph](../../reference_staging/templates_ldmos_first8_gprof_20260910/external_complete_callgraph.txt)
- [分析脚本](../../reference_staging/templates_ldmos_first8_gprof_20260910/analyze_hotspots.py)
- [分析文件 SHA256 清单](../../reference_staging/templates_ldmos_first8_gprof_20260910/analysis_inventory.json)

上述原始输出位于本地忽略目录；本报告不提交二进制、状态或采样文件。
从此工作树根目录配置 UCRT64 PATH 后构建：

```powershell
cmake --preset windows-ucrt64-release -B build-gprof-ldmos-first8 -DVELA_ENABLE_GPROF=ON -DVELA_ENABLE_LTO=OFF
cmake --build build-gprof-ldmos-first8 --target vela_example_runner --parallel 4
```

采集脚本接受 `--manifest` 与 `--output`，输出目录必须不存在；程序和源码须先冻结。
本次实际验证为 Release 插桩构建、完整 8 点仿真、原门限、状态/工作量/文件完整性审计及采样交叉核对。
未改动 C++，没有为本轮分析重复全套 CTest；上一轮 743/743 的结果保持独立。
