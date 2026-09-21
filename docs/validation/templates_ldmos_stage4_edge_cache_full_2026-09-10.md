# 优化版 LDMOS Vg=8、0–40 V 全曲线验证（2026-09-10）

优化版已从合格 0 V 初态完成 Vg=8 V、Vd=0–40 V 全曲线，31/31 个精确参考点齐全，
单栅压工程及最终判据全部通过。累计 5209 次 Newton 更新、0 次物理步回退，
26.6667 V 处的参考系平移等价性检查及最终完整性审计均通过。

运行时间为 2026-09-10 **09:47:32–11:46:04（Asia/Shanghai）**，控制器历时
7111.840 s，约 **1 小时 58 分 32 秒**。本轮补齐了前序缓存优化仅验证低压区间的
覆盖范围；本次未运行优化前的配对全曲线，不能据此宣称全曲线加速比例。

## 版本、配置和单位

- 基础提交 `d4285815f648be3f0b20719cbbe4d33e6a1b38b8` 加体区迁移率复用、几何缓存、
  端点状态缓存的工作区变更。实际程序为前序 742 项 Release 测试通过后冻结的
  `reference_staging/templates_ldmos_edge_cache_20260910/candidate.exe`。
- 程序 SHA256：`06301e6d44359dffa8889262524401b499c8a619cdb6b2bcb462d8c2f0734b6e`。
  三份优化源码的哈希、原差异与程序身份由前序 manifest 及本轮 plan 绑定，不能只用 HEAD 代表本轮版本。
- MSYS2 UCRT64 GCC 16.2.0，C++20，Release；实际后端 Eigen SparseLU/COLAMD，
  L2 行列均衡。内部性能计时开启，gprof 关闭，外推关闭。
- 原 exact-topology 网格：10241 节点、19782 三角形、30022 边，5723 个硅节点；
  external average-box 输运系数、barycentric 节点体积。输入几何 μm、密度 m⁻³、
  电势 V；本报告端电流为 A/μm。
- 300 K、Fermi–Dirac、OldSlotboom BGN、SRH/Auger、constant_field /
  transport_cell_vector QF 高场迁移率，接触电场回退开启；表面迁移率、雪崩、量子、热关闭。
- 使用原合格 Vd=0 V 种子，不包含栅压预偏置。种子 SHA256：
  `d1e48ccd2e089bbe4ef3fa8d3d6a25c904b003d104012c76d848484efc906ccf`。
- 外部检查点控制器初始步长 0.0025 V、上限 0.1 V、增长因子 1.35，逐一到达全部
  31 个参考偏压点；按原流程执行必要的同偏压再求解。psi/electron/hole 原块门限
  分别为 5e-8/1e-11/3e-10，局部行和 KCL 比值门限均为 1e-8，未调整验收门限。

## 完成情况与耗时口径

| 指标 | 结果 |
| --- | ---: |
| 控制器历时（含文件处理及最终审计） | 7111.840 s（1 h 58 min 32 s） |
| 子进程墙钟合计 | 6908.023 s（1 h 55 min 8 s） |
| 子进程 CPU 合计 | 5193.203 s（1 h 26 min 33 s） |
| 内部 dc_sweep.total 合计 | 6803.590 s |
| 精确参考点 / 物理推进 | 31 / 429 |
| 唯一合格状态 | 431（初态 + 429 次推进 + 平移后状态） |
| 子进程 / Newton 更新 / 物理步回退 | 859 / 5209 / 0 |
| Jacobian 装配 / 数值分解 | 5636 / 6479 |

5209 为实际已记账 Newton 更新数，5636 为装配/迭代阶段数，两者不可混用。
部分直接求解子进程返回非零，由既定严格同偏压再求解流程处理；本次没有中断重启，
全部 859 个子步骤的成本均进入账本，0 回退不表示每次直接求解都一次收敛。

## 与 Sentaurus 的单栅压曲线比较

参考为冻结 D5-no-IALMob 的 `IdVd_Vg2_n4_des_drain_curve.csv`，其物理栅压为 8 V。
电流相对误差统计使用 30 个非零且可解析电流点；31 个精确点均参与完整性和 KCL 检查。

| 指标 | 实测 | 工程门限 | 最终门限 |
| --- | ---: | ---: | ---: |
| 电流相对误差中位数 | 1.330933% | 15% | 5% |
| 电流相对误差 P95 | 1.612869% | 25% | 12% |
| 低压电阻误差 | 1.032119% | 20% | 10% |
| 40 V 端点电流误差 | 1.099059% | 20% | 10% |
| 最大归一化 KCL | 1.372946e-10% | 1% | 0.1% |

最大电流相对误差为 1.627051%。40 V 时 Vela 为 `4.2363315296624625e-4 A/μm`，
Sentaurus 为 `4.19027790569749e-4 A/μm`。低压电阻指标沿用评分脚本：首个非零
参考点的 V/I 比较。KCL 分母为最大绝对端电流；零偏压端电流和绝对 KCL 均为零，
没有新增平衡豁免或分母下限。

这里的“工程及最终通过”仅指 Vg=8 的上述五项检查；未合并旧 Vg=4 数据，
未检验本版本的双栅压电流比，不能表述为优化版双栅压 D5 全部通过。

## 26.6667 V 参考系平移

在物理 Vd=26.6666666666667 V 将求解参考系平移 28 V，保持各端物理电压差不变。
平移后的基线子步骤 `child_00578_frame_baseline` 即合格，四项等价性检查全部通过：

| 指标 | 实测 | 原门限 |
| --- | ---: | ---: |
| 最大电势/QF 绝对差 | 7.105427e-15 V | 1e-7 V |
| 最大密度相对差（分母下限 1 m⁻³） | 2.133292e-13 | 1e-6 |
| 电流差 / 最大原端电流 | 3.637536e-14 | 1e-8 |
| 平移前后最大 KCL 比值 | 8.371937e-13 | 1e-8 |

之后继续推进到物理 40 V。本轮确认缓存优化没有阻断该关键平移路径。

## 全曲线内部耗时

计时包含子调用，存在嵌套，各行不能相加；比例分母为内部总时间 6803.590 s。

| 阶段 | 累计时间 | 占内部总时间 |
| --- | ---: | ---: |
| linear.factorize | 2106.529 s | 30.96% |
| newton.jacobian | 1499.598 s | 22.04% |
| dd.residual | 1006.510 s | 14.79% |
| jacobian.edge_physics（包含在 Jacobian 中） | 865.031 s | 12.71% |
| dd.continuity_diagnostics | 692.418 s | 10.18% |

端点缓存命中 1,169,695,440 次、未命中 1,430,304,080 次；边物理阶段 FermiHalf
调用 2,145,523,740 次。符号分析执行 1694 次、缓存命中 4785 次，数值分解 6479 次。
这些为本轮工作量记录，不与低压区间调用总数直接比较。前序短区间的局部缓存收益仍以
配对报告为准；本轮提供全曲线正确性与当前耗时，不提供全曲线前后加速归因。

## 验证与证据

运行器完成原验收及完整性审计后，独立复核又检查了 20 项冻结程序/源码/输入哈希、
431 个唯一合格状态、859 个子步骤的程序/配置/父状态身份与性能文件哈希，重算了
31 点曲线误差和 KCL，结果全部一致。控制器及本次求解进程已退出。
本轮只新增运行适配脚本、审计结果和报告；未修改求解器源码，未重复编译或重跑已通过的 742 项测试。

以下运行证据保留在本地忽略目录，未提交生成的仿真输出：

- [运行配置及冻结输入清单](../../reference_staging/templates_ldmos_edge_cache_full_20260910/vg8_full/plan.json)
- [完整运行账本及平移检查](../../reference_staging/templates_ldmos_edge_cache_full_20260910/vg8_full/fixed/ledger.json)
- [曲线 CSV](../../reference_staging/templates_ldmos_edge_cache_full_20260910/vg8_full/score/curve.csv)
- [曲线判据](../../reference_staging/templates_ldmos_edge_cache_full_20260910/vg8_full/score/summary.json)
- [原完整性审计](../../reference_staging/templates_ldmos_edge_cache_full_20260910/vg8_full/audit_summary.json)
- [独立完成复核](../../reference_staging/templates_ldmos_edge_cache_full_20260910/completion_verification.json)
- [全曲线内部计时](../../reference_staging/templates_ldmos_edge_cache_full_20260910/vg8_full/hotspot_summary.json)
- [全曲线运行脚本](../../reference_staging/templates_ldmos_edge_cache_full_20260910/run_full.py)
- [完成复核脚本](../../reference_staging/templates_ldmos_edge_cache_full_20260910/verify_completed.py)
- [前序缓存优化与低压配对验证](templates_ldmos_stage4_edge_cache_2026-09-10.md)
