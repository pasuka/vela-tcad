# 五后端统一单线程完整曲线重复对比

## 执行合同

按用户确认顺序，先补跑 STRUMPACK 单线程双栅压完整曲线，通过后自动继续
五后端的统一单线程重复对比。新证据目录为
`reference_staging/templates_ldmos_single_thread_repeats_20260920`。
旧批次及其 STRUMPACK 双线程结果保留，不混入本批次耗时统计。

固定 R11 D5 等温模型、300 K、Vg=4/8 V、Vd=0–40 V 每档 31 个精确点。
五后端全部启用跨点分析复用，求解器/BLAS=1/1，OpenMP 动态线程关闭、
活跃嵌套上限为 1；逐服务核验实际线程。STRUMPACK 保持无压缩、METIS
及既有匹配/缩放设置。使用现有 UCRT64 Release，物理模型、原验收门限及
生产默认均不改变。

| 阶段 | 内容 | 曲线数 |
|---|---|---:|
| S0 | STRUMPACK 单线程 Vg4/Vg8 完整曲线及联合资格检查 | 2 |
| R0 | UMFPACK、SparseLU、MUMPS、SuperLU_MT 第一轮双栅压 | 8 |
| R1 | 五后端第二轮，轮换顺序，Vg8 反序 | 10 |
| R2 | 五后端第三轮，再次轮换顺序，Vg8 反序 | 10 |

合计 30 条完整曲线、930 个精确点、15 组双栅压联合验收。S0 两条曲线
预先定义为 STRUMPACK 第一轮，不额外计算一次。所有曲线串行运行；其他四个
后端也重新采集三轮，避免把不同后台负载下的历史首轮直接并入统计。
原计划 D4 后续资格仍待开展，不由这批 D5 结果自动继承。

## 资格与证据

S0 对照旧批次已通过的 U0 完整状态，冻结其输入、状态、账本及摘要散列，
并要求新旧程序和运行 DLL 散列一致；旧 U0 仅作正确性对照，不作本批计时基准。
S0 两档状态一致性及原双栅压工程/最终门限必须通过才进入 R0。
其后同轮以 U1 为跨后端状态对照，各配置重复轮对照本批自身首轮，
电势/准费米势差异继续采用既有 1e-8 V 门限。失败即停止，不带失败结果晋级。

记录端到端墙钟、CPU、Newton 更新、失败/回退、分析与数值分解计数及后台负载。
计时边界沿用旧批次：曲线运行控制器及其子服务计时，外部汇总验收不计入曲线墙钟。
最终报告三轮原始值、中位数及范围；单轮不宣称稳定加速。

运行入口：[单线程批次](../../scripts/run_templates_ldmos_single_thread_repeats.py)，
复用[原曲线执行/验收逻辑](../../scripts/run_templates_ldmos_full_curve_timing.py)。
该批次由一个串行进程推进，无旧 supervisor；旧 48 条队列不恢复。
`summary.json` 为状态，`analysis.json` 在阶段完成时更新，`batch.lock` 保存 PID；
用户暂停可使用批次 `PAUSE` 文件，保留中断尝试和已完成证据，显式恢复需通过冻结检查。

## 启动前检查

25 项 Python 测试通过：原调度/汇总/监督 12 项、新单线程调度与资格保护 4 项、
后端统计和线程审计 9 项。Release runner 构建检查返回无需重建。
启动前确认无仿真进程、D 盘可用 52.5 GB。
已恢复已有每 20 分钟 heartbeat，检查本批 30 条曲线、真实进程及账本更新时间；
完成或失败/用户暂停后停用检查，不自动启动额外批次。

## 启动状态

2026-09-20 11:14（北京时间）确认串行驱动和仿真进程存活，S0 正在运行
STRUMPACK Vg=4 V，已接受 Vd=0.346212 V，1/31 精确点，0 次外层回退；
全批 0/30 完成。原 heartbeat `ldmos` 已恢复 ACTIVE 并改为监测本新目录。
该记录仅为启动快照，实时状态以证据目录为准。

## 用户要求第二轮结束后暂停（2026-09-20 14:56）

S0/R0/R1 已完成；20 条完整曲线、620 精确点、10 组双栅压联合验收全部通过，
跨后端及重复状态检查通过，0 次外层回退。第三轮 R2 的 10 条曲线未启动。
最后 STRUMPACK 第二轮 Vg8 耗时 810.303 s。

在最后一条 R1 曲线控制器正常退出后设置 PAUSE，驱动完成该曲线验收及
R1 汇总后响应用户暂停，summary.status=paused，active 清空，batch.lock 已移除。
stop_reason=user_interrupt / KeyboardInterrupt('PAUSE requested') 是正常用户暂停，
不是仿真失败。已确认无 python/vela_example_runner 残留进程。

已有每 20 分钟 heartbeat ldmos 已设为 PAUSED。保留原冻结输入、程序、
曲线状态、日志、summary.json 与 R1 analysis.json，不自动恢复。
用户后续明确恢复时，先核验冻结依赖及无进程，删除本批 PAUSE 文件后运行
`python scripts/run_templates_ldmos_single_thread_repeats.py --output reference_staging/templates_ldmos_single_thread_repeats_20260920 --resume`，
仅续跑尚未完成的 R2；是否恢复定时检查遵循后续用户要求。
