# IALMob Halley 默认化集成与复核

日期：2026-09-25。当前状态：默认化及 Codespaces Release 资格复核完成。
完整 CTest 858/858、7 条长曲线/217 个精确点及四组原验收全部通过。

## 范围

在既有四算法完整曲线及 Halley/legacy 三轮重复资格基础上，统一默认入口。
本轮不修改物理公式、屏蔽方程、导数、区间保护、最终残差检查、自动回退、
温度缓存键或原电学/热学验收门限。性能依据仍为
[既有完整曲线与重复计时](templates_ldmos_screening_full_2026-09-25.md)。

- C++ 求根接口、精确缓存、IalMobility 构造使用同一 `defaultIalScreeningMethod`。
- 等温 `solver.mobility.ialmob.screening_method`、电热输入
  `mobility_SI.ialmob.screening_method` 省略时均使用 `halley`。
- 显式 `legacy` 保留；`newton`/`toms748` 仍可用于对照。
- 电热旧诊断字段是显式别名，与正式字段冲突时报错。
- 原生参数及对照测试不能再通过省略算法来暗示 legacy，基线明确指定。

## 验证协议

新目录 `/workspaces/vela-halley-default-20260925`，保留旧目录
`/workspaces/vela-screening-20260924` 的源码、程序及证据不变。
GCC 16 Release、UMFPACK、求解器/BLAS 单线程，HDF5/HighFive 状态接口。
本次独立构建的可选 OpenBLAS 运行时控制 API 未启用；等温启动不设置
`VELA_BLAS_THREADS`，使用 `OPENBLAS_NUM_THREADS=1`，并在新进程加载
实际 `libopenblas.so.0` 核验线程数为 1。没有修改二进制或物理配置来跳过检查。

1. 完整 CTest；补测默认/显式算法、非法值及冲突、精确缓存算法隔离、
   电热点服务回显；保留高精度根、耦合温度和完整装配 Jacobian 检查。
2. 生产 D0 从零初始化，双栅压 0–40 V 各 31 点，配置省略算法；逐点检查
   实际算法、原稳态门限、求根回退，并与上轮显式 Halley 全状态比较（≤1e-8）。
   再比较原生电流、温升场、峰值温升及热平衡。
3. G3 原 31 点 Id–Vg；D5/D4 各双栅压完整 62 点，使用各自物理配置、
   冻结输入散列及原曲线验收门限。合计 7 条长曲线、217 点。

独立复核入口为
[`run_templates_ldmos_halley_default.py`](../../scripts/run_templates_ldmos_halley_default.py)。
Linux 等温 DC worker 的 CPU 计数读取 `/proc/<pid>/stat`，不伪造 Windows
计时。端到端计时只记录本轮成本，不作为新的重复加速结论。

## 结果

首轮完整 CTest 为 855/858：两项历史二分算法测试仍省略算法，导致默认
切换后逐位等价断言比较了不同算法；已显式固定 legacy，保留逐位断言。
新增配置夹具缺少必需的 QF 高场驱动力，已补齐。生产算法及门限未修改。

修正后 GCC16 Release 完整 CTest **858/858 通过**（28.32 s）。
另有屏蔽重复计时协议 Python 测试 4/4 通过。
初次长曲线启动因 `qualification` 目录已存在而安全拒绝写入；未启动求解。
正式批次使用新目录 `evidence/default_20260925`，保留既有目录不变。

D0 双栅压完整 62 点与上轮显式 Halley 的五组保存状态字段最大差均为 0，
原生电流/温升场/热平衡联合验收通过。Vg4/Vg8 墙钟 152.617/140.788 s，
只记录本轮复核成本。G3 原 31 点及六项原门限通过，求解墙钟 35.905 s。

G3 初次准备遇到 Windows 反斜杠路径，修正迁移脚本后继续；其首次求解启动
因请求未启用的可选 OpenBLAS 控制 API 而退出，改用上述已核实环境接口。
这两次均未形成有效 G3 曲线。D0 已完成证据保留，后续脚本版本独立保存，
核对既有冻结源码/程序/摘要散列后复用 D0；不重写旧摘要或移用别的程序资格。

最终批次 `evidence/default_final_20260925/summary.json` 为 `completed`，
D0/G3/D5/D4 均为 `pass`。D0 沿用并重新审核同一程序的已完成证据，
G3/D5/D4 均完整运行，没有继承旧曲线的“通过”状态。

| 复核 | 精确点 | 关键结果 | 本轮墙钟 |
|---|---:|---|---:|
| D0 Vg4 / Vg8 | 62 | 最大电流误差 0.0240871% / 0.00714390%；原热学门限通过 | 152.617 / 140.788 s |
| G3 Id–Vg | 31 | Vth 差 10.7524 mV；峰值 gm 误差 1.61008%；六项门限通过 | 35.905 s |
| D5 Vg4 / Vg8 | 62 | 最大电流误差 1.592417% / 1.625083%；联合验收通过 | 396.345 / 386.845 s |
| D4 Vg4 / Vg8 | 62 | 最大电流误差 1.891761% / 1.801547%；联合验收通过 | 667.358 / 655.747 s |

D0 漏压 Newton 更新 344/303、试算 472/424、装配 784/692；求根回退 0/0，
失败减步 1/0，与已资格化的显式 Halley 对照一致。62 对状态的最大绝对差为 0。
表中等温耗时包含外层逐点控制，与 D0 的生产扫描流程不同，不用于跨物理模型
排名；本轮未追加三轮性能复测，也不替换此前 6.70% 的同程序配对证据。

## 可追溯证据

- 最终 Release 程序 SHA-256：
  `e2105ad756943807becae861cdcbed6bb7a5d25711a6a72a22fec7f8211f7065`。
- 本地摘要：`build/halley_default_20260925/retrieved/`，含 CTest 日志、原生
  联合报告、曲线、扫描台账、源码散列、构建选项及 `retrieval_manifest.json`。
  59 个清单内证据文件的原始字节 SHA-256 全部通过；本地与远端 438 个
  C++/头文件/测试源码、663 个其他脚本经 CRLF/LF 标准化后完全一致，
  最终复核脚本也一致，详情见
  `build/halley_default_20260925/local_source_check.json`。
- 远端完整状态：`/workspaces/vela-halley-default-20260925/evidence/`。
  `default_20260925` 保存 D0；`default_final_20260925` 保存最终等温及汇总。
- 迁移脚本的两个早期版本留在远端，保证已冻结 D0 源码路径的散列不被改写；
  最终执行的 `run_templates_ldmos_halley_default_v3.py` 对应仓库内当前复核脚本。

本轮提交采用 Halley 默认值，显式 legacy 及保护回退保留。
本轮资格限于上述 Codespaces 构建和算例，不额外认领其他平台的性能收益。
确认无遗留仿真/构建任务、证据回收核验完成后，已停止本次启动的 Codespace。
