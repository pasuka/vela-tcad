# 完整曲线后端运行时间对比：执行记录

状态：旧批次 P1 D5 首轮 12 条完整曲线和 6 组双栅压联合验收全部通过；
已在候选名单复核处停止，无仿真进程运行，P2/D4 尚未执行。
用户随后要求统一单线程，现场调度器 T1 已改为 1/1；旧批次 T1 的 2/1 证据
保持原样，禁止直接用修改后的现场脚本恢复旧队列。后续须另建单线程批次。
依据[已批准计划](templates_ldmos_full_curve_timing_plan_2026-09-20.md)执行，默认配置、
物理模型和原验收门限保持不变。

## 运行入口与验证

新增 `scripts/run_templates_ldmos_full_curve_timing.py`，组织 U0/U1/T1/L1/M1/S1
六配置、D5 首轮与重点重复、D4 前八点与完整重复，共 48 条计划曲线。
底层仍使用已验证的 linked D5/D4 扫压程序，不修改 C++ 求解器。

- 新增调度测试 6 项通过：总量/输出点数、栅压反序、轮次轮换、必须通过原 final
  联合门限、后到 U0 对照触发状态核验、状态不一致阻止资格、恢复合同不可变。
- 既有相关 Python 回归 26 项、真实 worker 协议测试 2 项通过。
- UCRT64 Release 运行器目标检查通过，`ninja: no work to do`；程序未重编译。
- 修复新测试最初的局部变量作用域错误后，6 项全部通过；此问题未进入仿真程序。
- 冻结程序/DLL/源码/输入及哈希，另记录参与调度的现场 Python 脚本哈希；恢复时
  检查这些依赖未变。每条曲线独立启动 worker，从合格零压种子开始。
- 每个完整双栅压组执行现有 engineering 和 final 联合评分；每条曲线执行
  独立状态/门限审计和实际 BLAS 线程核验，数值失败停止批次以供定位。
- 每个阶段可明确停止。`PAUSE` 文件只在用户要求暂停时创建；中断会终止本批
  正在运行的子进程树、保存已用时间及日志；下次从原种子重跑该条曲线并追加
  attempt 名称，不覆盖部分轨迹，不将中断耗时混入连续完成计时。

## 当前批次

本轮输出目录：
`reference_staging/templates_ldmos_full_curve_timing_20260920`。

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts/run_templates_ldmos_full_curve_timing.py --output reference_staging/templates_ldmos_full_curve_timing_20260920 --through P1
```

以下是原冻结配置未变时的恢复命令，仅保留作历史说明；由于当前线程合同
已经修订，不应再对本旧目录执行这些命令：

```powershell
python scripts/run_templates_ldmos_full_curve_timing.py --output reference_staging/templates_ldmos_full_curve_timing_20260920 --resume --through P2
python scripts/run_templates_ldmos_full_curve_timing.py --output reference_staging/templates_ldmos_full_curve_timing_20260920 --resume --through P3short
python scripts/run_templates_ldmos_full_curve_timing.py --output reference_staging/templates_ldmos_full_curve_timing_20260920 --resume --through P3full
```

`--resume` 只接受明确 paused 且无活动曲线的批次；failed 或异常遗留锁须先查明，
不自动清锁、掩盖失败或更换配置。每个阶段的结束属于检查边界，不表示用户
取消后续已授权任务。

此前已启用 `supervise_templates_ldmos_full_curve_timing.py` 对本次授权批次进行阶段推进。
它等待当前阶段结束及锁释放，生成独立分析后再恢复下一阶段；数值/联合评分失败、
用户暂停、阶段缺项，或 L1/M1/S1 首轮在任一栅压快于 U1/外层回退更少时停止，
保留 `supervision.json` 的检查原因，不静默排除可能的优胜候选。
新增 3 项阶段推进保护测试和 2 项计时分段/统计测试通过。

```powershell
python scripts/supervise_templates_ldmos_full_curve_timing.py reference_staging/templates_ldmos_full_curve_timing_20260920
```

同一批次只能启动一个阶段推进程序；其独立锁防止重复启动。推进与分析记录分别
保存程序自身及分析器 SHA-256。阶段推进已运行时不要手工同时执行上述 resume 命令。

首条 D5 Vg4 U0 完整曲线已通过：31 点到 40 V，墙钟 647.881 s，Newton 1125 次，
0 次外层回退；全局分析 356 次、数值分解 1271 次，进程树 CPU 569.73 s。
首个 Vg4 完整配对随后通过，31 点 ψ/fn/fp 最大差为 0：

| 配置 | 墙钟 / s | 进程树 CPU / s | Newton | 全局分析 | 数值分解 |
|---|---:|---:|---:|---:|---:|
| U0 普通 worker | 647.881 | 569.734 | 1125 | 356 | 1271 |
| U1 跨点分析复用 | 658.468 | 558.297 | 1125 | 101 | 1271 |

本轮墙钟增加 1.63%，CPU 减少 2.01%，分析次数显著下降；尚不能认领完整曲线
稳定加速。继续按既定顺序运行 T1 及 Vg8，不因本次墙钟变慢而删除或替换记录。
这仅为首档栅压首轮配对，双栅压联合验收和重复性能结论仍待后续完成。

主要墙钟沿用昨日完整漏压扫压范围，排除输入冻结与额外后审计；包含底层扫压
控制器运行内评分。后置双栅压评分另行执行。仅凭短曲线或部分完整点不做排名。

P1 已完成结果如下；旧 T1 为双线程，不能据此排序单线程后端：

| 配置 | Vg4 墙钟 / s | Vg8 墙钟 / s | 求解器/BLAS |
|---|---:|---:|---|
| U0 UMFPACK 普通 worker | 647.881 | 591.959 | 1/1 |
| U1 UMFPACK 复用 | 658.468 | 600.253 | 1/1 |
| T1 STRUMPACK 复用（历史） | 570.753 | 566.076 | 2/1 |
| L1 SparseLU 复用 | 675.972 | 641.169 | 1/1 |
| M1 MUMPS 复用 | 677.144 | 599.498 | 1/1 |
| S1 SuperLU_MT 复用 | 615.640 | 613.787 | 1/1 |

全部 372 个精确点通过原门限，0 次外层回退；配对 ψ/fn/fp 最大差
3.3306690738754696e-14 V。没有三轮重复资格。推进程序因 M1 在 Vg8、S1
在 Vg4 快于 U1 而进入 `review_required`，不是数值失败；定时检查已暂停。

原始证据以批次 `summary.json`、逐曲线
日志/ledger、完整评分及双栅压联合目录为准；不提交生成输出。
