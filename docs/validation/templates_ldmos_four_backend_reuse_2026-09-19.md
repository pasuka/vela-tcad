# SparseLU / UMFPACK / MUMPS / SuperLU_MT 跨点曲线验证

状态：通过本轮四后端 D5 双栅压前 8 点配对验证。实验默认关闭，生产默认及物理/数值验收门限不变。

## 范围与合同

承接[主 Newton 跨点复用实现](templates_ldmos_cross_point_analysis_2026-09-19.md)，
本轮只扩展实验脚本的后端和对照模式选择，不更改 C++ 求解算法。
四后端分别比较普通常驻 worker 与开启主系统跨点复用的 worker，
各跑 Vg=4/8 V、Vd=0–9.333333 V 的前 8 个精确点，共 16 条短曲线。

- MSYS2 UCRT64 Release，`-O3 -DNDEBUG`、无 `-pg`；实际构建具备
  UMFPACK、MUMPS、SuperLU_MT 和 OpenBLAS 显式线程控制功能。
- 四后端统一请求求解器单线程、OpenBLAS 单线程；MUMPS/SuperLU_MT 适配器
  设置 OpenMP 线程数为 1。每个有效点服务核验 OpenBLAS API 报告值为 1。
  本轮没有新增 MUMPS/SuperLU_MT 内部并行区域线程采样，不能把 BLAS 核验
  称为所有线程池的独立实测。SuperLU_MT 不涉及此前未合格的 2/4 线程配置。
- 排序沿用各后端既有默认，不增加 METIS 变体，不在配对中改变排序或缩放。
- 固定 R11 D5 300 K 等温模型，无 IALMob、Auger 不含生成；10,241 节点、
  19,782 三角形，单位和输入来源沿用冻结 bundle，电流输出 A/μm。
- 输入合同：`reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_auger_no_generation_inputs.json`。
  每组冻结源码、可执行文件、DLL、构建清单及输入哈希；不覆盖既有证据。
- 每后端内部仅改变跨点复用开关。Vg4 先普通 worker 后复用；Vg8 反序。
  全部串行运行，记录后台负载。每配置一轮，不能认领重复稳定性。
- 从已有合格零压种子检查开始计时，含扫压控制器、点服务、状态读写、
  原有门限检查；不含生成种子的栅压预偏置、冻结和独立后审计。
- 所有接受状态检查原块残差、逐行闭合、KCL 门限；两模式每个精确点
  ψ/fn/fp 最大差须 ≤1e-8 V。保留失败、恢复和额外候选的实际成本。
- 全局分析次数含辅助 Poisson 与恢复系统。主系统复用不代表全部辅助系统
  也只分析一次；不把全局计数误称为主 Jacobian 分析次数。

## 复现

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build --preset windows-ucrt64-release --target vela_example_runner --parallel 2
python -m unittest tests.regression.test_analysis_reuse_controls tests.regression.test_templates_ldmos_linked_d5 tests.regression.test_templates_ldmos_backend_statistics
foreach ($backend in @('sparselu','umfpack','mumps','superlu_mt')) {
    $destination = "reference_staging/templates_ldmos_four_backend_reuse_20260919/$backend"
    python scripts/run_templates_ldmos_analysis_reuse.py --output $destination --backend $backend --threads 1 --modes worker reuse
    if ($LASTEXITCODE -eq 0) {
        python scripts/analyze_templates_ldmos_analysis_reuse.py $destination
    }
}
```

Release 目标检查无需重建；26 项相关 Python 测试通过，包括后端/线程约束、
必须具备独立对照、保留原三模式入口和跨点策略恢复合同。未改 C++ 核心，
不重复认领本轮重跑了此前的 C++ 或完整 CTest 验证。

## 实测结果

全部 16 条曲线、128 个精确点通过原有电流、块残差、逐行闭合与 KCL 检查；
1,088 次点请求中，1,072 个有效求解服务的 OpenBLAS API 单线程核验通过。
八组配对推进序列与 Newton 更新序列相同，ψ/fn/fp、电子和空穴密度均逐值一致
（最大绝对差为 0）；外层回退全部为 0。

下表箭头均为普通 worker → 跨点复用。时间单位 s；分析/分解是包含辅助算子的
全局计数。Newton 两模式相同，仅列一个数值。

| 后端 | Vg (V) | 墙钟 (s) | 降幅 | 进程树 CPU (s) | Newton | 分析次数 | 数值分解次数 |
|---|---:|---:|---:|---:|---:|---:|---:|
| SparseLU | 4 | 280.87 → 273.99 | 2.45% | 244.25 → 238.05 | 374 | 126 → 55 | 455 → 455 |
| SparseLU | 8 | 275.04 → 274.41 | 0.23% | 240.22 → 236.44 | 358 | 124 → 53 | 437 → 437 |
| UMFPACK | 4 | 221.56 → 188.29 | 15.02% | 193.36 → 172.22 | 369 | 125 → 55 | 446 → 446 |
| UMFPACK | 8 | 191.23 → 188.36 | 1.50% | 173.97 → 171.19 | 349 | 122 → 51 | 426 → 426 |
| MUMPS | 4 | 228.26 → 215.04 | 5.79% | 208.67 → 197.44 | 364 | 120 → 50 | 436 → 436 |
| MUMPS | 8 | 245.24 → 217.45 | 11.33% | 214.67 → 196.34 | 349 | 122 → 51 | 426 → 426 |
| SuperLU_MT | 4 | 235.72 → 217.92 | 7.55% | 209.39 → 199.61 | 375 | 124 → 54 | 451 → 451 |
| SuperLU_MT | 8 | 226.24 → 216.40 | 4.35% | 205.41 → 197.58 | 357 | 123 → 52 | 435 → 435 |

### 复用范围与计时解释

每条复用曲线均只创建一次 DC 线性上下文，后续 67 次请求命中；数值缓存随新点
失效，数值分解次数没有减少。每个后端/栅压的 54 个纯主系统有效服务，新增分析
从普通 worker 的 54 次降为 0，证明主系统分析对象实际跨点保留。

上下文也由恢复时的 NewtonConfig 副本继承，复用同时覆盖跨点和同点内的 Newton
重试，不能把所有减少的分析次数只理解为跳过每点第一次分析。辅助 Poisson
重闭合及载流子恢复使用独立求解器，未接入本轮主系统上下文；混合服务的聚合
计数不能精确分配给各算子，不能宣称整条曲线的所有系统只分析一次。

| 后端 | Vg (V) | 符号分析时间 (s) | 数值分解时间 (s) | 回代时间 (s) | 求解器进程峰值 (MiB) |
|---|---:|---:|---:|---:|---:|
| SparseLU | 4 | 4.698 → 0.932 | 103.267 → 103.757 | 2.648 → 2.582 | 302.9 → 261.7 |
| SparseLU | 8 | 4.704 → 0.840 | 100.355 → 100.278 | 2.566 → 2.602 | 302.4 → 261.8 |
| UMFPACK | 4 | 3.059 → 0.448 | 49.798 → 42.365 | 7.274 → 6.149 | 273.8 → 248.8 |
| UMFPACK | 8 | 2.631 → 0.428 | 43.024 → 42.960 | 5.849 → 5.883 | 271.9 → 247.2 |
| MUMPS | 4 | 13.976 → 2.179 | 52.989 → 53.438 | 22.165 → 22.496 | 277.9 → 244.6 |
| MUMPS | 8 | 15.453 → 2.258 | 57.861 → 51.993 | 23.498 → 21.972 | 277.4 → 245.5 |
| SuperLU_MT | 4 | 4.596 → 0.811 | 77.291 → 73.145 | 3.019 → 2.898 | 318.1 → 268.3 |
| SuperLU_MT | 8 | 4.772 → 0.758 | 72.371 → 69.860 | 2.831 → 2.743 | 317.3 → 266.7 |

求解器进程峰值是采样到的进程生命周期内存高水位，不是整个进程树同时占用。
CPU 口径为控制器 CPU 加点服务记录 CPU，不包含 worker 关闭尾部及独立后审计。
各曲线平均系统忙碌率约 40.2%–51.6%，保留原始负载采样；不是空载专用机器。

符号分析耗时在全部配对中下降，但单轮总墙钟降幅不能全部归因于符号分析。
例如 UMFPACK Vg4 总墙钟减少 33.28 s，符号分析只减少约 2.61 s，数值分解和
其他阶段也变快。当前证据无法把负载变化、工作空间复用及其他运行波动分开。
因此这些是单轮观测值，不是已经确认的稳定加速比，也不与早先不同线程/负载的
STRUMPACK 批次直接排名。

本轮没有减少主体 Newton 更新，验证的是线性分析对象复用的正确性与成本。
保持实验默认关闭；本轮不新增完整 31 点、D4/D0、跨温度、SuperLU_MT 2/4 线程
或重复稳定性资格。下一阶段可在冻结配置下做完整曲线重复配对，再决定推荐配置。

## 证据与变更边界

- [批次汇总](../../reference_staging/templates_ldmos_four_backend_reuse_20260919/collection_summary.json)：16 条曲线的结果、计数与配对聚合。
- SparseLU：[原始汇总](../../reference_staging/templates_ldmos_four_backend_reuse_20260919/sparselu/summary.json)、[分析](../../reference_staging/templates_ldmos_four_backend_reuse_20260919/sparselu/analysis.json)、[冻结清单](../../reference_staging/templates_ldmos_four_backend_reuse_20260919/sparselu/binary/manifest.json)。
- UMFPACK：[原始汇总](../../reference_staging/templates_ldmos_four_backend_reuse_20260919/umfpack/summary.json)、[分析](../../reference_staging/templates_ldmos_four_backend_reuse_20260919/umfpack/analysis.json)、[冻结清单](../../reference_staging/templates_ldmos_four_backend_reuse_20260919/umfpack/binary/manifest.json)。
- MUMPS：[原始汇总](../../reference_staging/templates_ldmos_four_backend_reuse_20260919/mumps/summary.json)、[分析](../../reference_staging/templates_ldmos_four_backend_reuse_20260919/mumps/analysis.json)、[冻结清单](../../reference_staging/templates_ldmos_four_backend_reuse_20260919/mumps/binary/manifest.json)。
- SuperLU_MT：[原始汇总](../../reference_staging/templates_ldmos_four_backend_reuse_20260919/superlu_mt/summary.json)、[分析](../../reference_staging/templates_ldmos_four_backend_reuse_20260919/superlu_mt/analysis.json)、[冻结清单](../../reference_staging/templates_ldmos_four_backend_reuse_20260919/superlu_mt/binary/manifest.json)。

四组使用同一 Release 可执行文件，SHA-256：
`70054fdaf5ca3b1f2608eec812948e6193acb7df1733fe00085bc6c96be98fe6`。
原始曲线、配置、状态与冻结输入放在上述忽略目录内，不提交生成仿真输出。

运行脚本新增后端、线程及模式选择，保留此前 STRUMPACK 三模式默认；分析脚本
支持普通 worker/复用两模式，防止把缺失 subprocess 当成已执行对照。
新增测试检查有效配对与线程约束。本轮结束后另修正切换曲线时进度显示沿用上条
曲线耗时的问题，只重置显示字段，不改变已冻结证据或实际墙钟计算。

提交整理时，Release 运行器、线性/Newton/DC 测试目标及两个矩阵诊断目标构建
通过；矩阵加速实验 4 项测试、worker 协议 2 项测试通过。此前本轮 26 项 Python
回归及已记录的 C++ 验证结果保持有效；提交整理未重跑完整仿真或完整 CTest。
