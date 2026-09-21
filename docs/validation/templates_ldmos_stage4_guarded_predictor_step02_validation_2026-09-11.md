# LDMOS Stage 4：带保护的预测器与 0.2 V 步长验证（2026-09-11）

本次在原门限下执行“内部 QF 重定心＋带保护的 secant 预测器＋0.2 V 最大步长”验证。先完成最大步长 0.1 V 的保护规则对照，再用同一脚本、二进制、初始种子和门限，仅将上限改为 0.2 V。两个候选串行运行，均从 Vd=0 V 开始，未复用对照的中间解。

## 范围与固定条件

Vg=8 V，精确 Vd 点为 0、1.33333333333333、2.66666666666667、4、5.33333333333333、6.66666666666667、8、9.33333333333333 V；不包含栅压预偏置，也没有运行完整 0–40 V。

使用同一个冻结 Windows UCRT64 Release runner，SHA256 为 `be486a03d52184b2944e435567c42ecfa643f7fee36a27449aa2bb4f810485f2`，Eigen SparseLU/COLAMD，`-O3 -DNDEBUG`，gprof 关闭。本次没有修改 C++、重新编译或改变默认开关。

保留 D5 300 K、Fermi–Dirac、OldSlotboom BGN、SRH/Auger、constant_field vector-QF HFS、contact_node_cell 场回退配置；10,241 节点、19,782 三角形、5,723 个 Si 节点，external-averagebox 输运系数。几何单位 μm、密度 m⁻³、电流 A/μm。

保持原 Poisson/electron/hole 块门限 5e-8/1e-11/3e-10，carrier-row 与 KCL 比率门限 1e-8。内部重定心开启，最多 8 次；Newton 更新预算 160。QF 更新上限随实际推进步长匹配。外层初始步长 0.0025 V、增长因子 1.35，原重闭合、密度恢复和回退流程不变。

## 预测保护

上一成功步若为命中精确参考点而被截短，下一步不做外推；计算出的外推比例超过 2、缺少两个成功历史态、发生回退重试或进入恢复流程时，也不外推。其他 direct 步沿原 secant 公式预测 psi/phin/phip，并保存一致的 QF reference/increment。未外推时沿用原合格父状态，不改变验收门限。

这与上次“比例截到 2 后继续外推”不同。本次选择先跳过风险初值；并没有验证更长历史拟合或新的电场平滑模型。

## 完整前 8 点结果

| 策略 | Newton 更新 | 推进 | 子进程 | 回退 | child wall / s | child CPU / s | 控制器总时间 / s |
|---|---:|---:|---:|---:|---:|---:|---:|
| 上次：重定心＋原外层 secant，0.1 V | 599 | 110 | 117 | 3 | 517.49 | 425.67 | 572.59 |
| 本次：带保护 secant，0.1 V | 436 | 107 | 110 | 0 | 297.18 | 278.86 | 340.00 |
| 本次：带保护 secant，0.2 V | 619 | 67 | 90 | 6 | 498.75 | 433.08 | 538.17 |

当前 0.2 V 相对当前 0.1 V 对照：更新变化 +41.97%，child wall 变化 +67.83%，CPU 变化 +55.30%。

时间为单次顺序运行的实测值，没有交替 AB 重复和置信区间；旧结果来自前一天。child wall 包括所有子进程启动、求解、失败和恢复；控制器总时间还包含外推 CSV、协调和评分。性能审计在计时求解结束后执行。不得将部分曲线时间当成完整曲线比较。

## 精度与逐步审计

- guard01：8 点原工程/最终可用门限全部通过。电流误差中位数 1.594856%，P95 1.624118%，低压差分电阻误差 1.032119%；最大 KCL 比率 1.69275e-13。
  数值分解 479 次，内部重定心 7 次。保护统计：`{"non_direct": 3, "no_history": 1, "predicted": 100, "previous_step_clipped": 6}`。
  成功的完整 0.2 V 步有 0 次，使用密度恢复的成功推进有 1 次。
- guard02：8 点原工程/最终可用门限全部通过。电流误差中位数 1.594856%，P95 1.624118%，低压差分电阻误差 1.032119%；最大 KCL 比率 2.87855e-13。
  数值分解 1003 次，内部重定心 10 次。保护统计：`{"non_direct": 17, "no_history": 1, "predicted": 60, "previous_step_clipped": 6, "retry": 6}`。
  成功的完整 0.2 V 步有 27 次，使用密度恢复的成功推进有 5 次。

| 方案 | 阶段 | 子进程数 | 失败数 | Newton 更新 | child wall / s |
|---|---|---:|---:|---:|---:|
| guard01 | initial | 1 | 0 | 0 | 2.05 |
| guard01 | direct | 107 | 1 | 434 | 289.20 |
| guard01 | reclose | 1 | 1 | 0 | 2.75 |
| guard01 | density | 1 | 0 | 2 | 3.18 |
| guard02 | initial | 1 | 0 | 0 | 0.79 |
| guard02 | direct | 73 | 11 | 568 | 432.34 |
| guard02 | reclose | 11 | 11 | 33 | 47.52 |
| guard02 | density | 5 | 0 | 18 | 18.10 |

独立 audit.py 重新从 transfer/rollback 时间顺序重建每个预测保护判定，核对逐节点外推公式、参考值/增量一致性、输入和输出哈希、原块/逐行/KCL 门限，并核对 profile 的 Newton 计数与 ledger。两个方案的冻结输入、预测保护、初始种子和门限相同，仅 max_step_V 不同。所有检查通过。

## 结论与剩余问题

两组都完整通过原门限，但当前带保护预测器的 0.2 V 方案没有性能收益：推进从 107 次减少到 67 次，更新却从 436 增至 619，child wall 从 297.18 s 增至 498.75 s。因此推荐保留已验证的带保护预测器＋0.1 V 作为当前候选流程。

0.2 V 方案的 6 次物理回退全部发生在 1.333333、2.666667、4、5.333333、6.666667、8 V 六个中间精确参考点之后：上一小步被截短，预测保护触发，直接沿用旧解尝试前进 0.2 V；均失败并缩到 0.1 V。8→8.2 V 的失败 direct 尤其昂贵，耗费 123 次更新，回退后的 8→8.1 V 只需 8 次更新即通过。6 次被回退的尝试（含对应 reclose）合计消耗 241 次更新、196.59 s，占本方案总更新和 child wall 的 38.93% 与 39.42%。

这支持将下一轮保护与步长联动：当预测不可用或主动跳过时，将该次尝试限制为 0.1 V；有可信预测时再允许 0.2 V。另一条途径是为参考点截短保留更合适的历史间距。上述后续改动本次均未执行，不将估计收益作为实测结果。

本次仅新增实验脚本和验证报告，没有修改生产模板或默认门限。相同 C++ 冻结二进制在上次已通过完整 Release 串行 CTest 746/746；本次未重跑该旧测试，实际新增验证为两条完整前 8 点曲线与独立审计。预测器仍位于外层实验适配器，不能声称内置精确点扫压预测器已修复，也不能外推到更大漏压或其他物理配置。

## 证据与复现

原始证据位于工作树忽略目录 `reference_staging/templates_ldmos_guarded_predictor_20260911/`：`run_guarded.py`、`guard01/`、`guard02/`、各自 `plan.json` / `progress.json` / `fixed/ledger.json` / 子进程 profile，以及 `audit.py` / `audit.json`。每次 plan 冻结二进制、脚本和物理输入哈希；实验目录不得覆盖。

原二进制 manifest：`reference_staging/templates_ldmos_qf_recenter_validation_20260910/binary/manifest.json`。上次结果见[重定心及预测器验证报告](templates_ldmos_stage4_qf_recenter_validation_2026-09-10.md)。

复现使用 UCRT64 Python，在工作树根运行下列命令；输出目录必须不存在。

```powershell
D:/msys64/ucrt64/bin/python.exe -X utf8 reference_staging/templates_ldmos_guarded_predictor_20260911/run_guarded.py --manifest reference_staging/templates_ldmos_qf_recenter_validation_20260910/binary/manifest.json --output reference_staging/guard02_repeat_20260911 --recenter --outer-secant --max-step 0.2
```
