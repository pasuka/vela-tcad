# LDMOS Stage 4：预测器与步长联动验证（2026-09-11）

本次验证“无预测时限制为 0.1 V，有可信预测时允许 0.2 V”。在独立输出目录从同一合格零压种子开始，完整运行 Vg=8 V 的前 8 个精确参考点；不复用先前候选的中间解。

## 固定条件与改动范围

精确 Vd 为 0、1.33333333333333、2.66666666666667、4、5.33333333333333、6.66666666666667、8、9.33333333333333 V。没有运行 0–40 V，也不计栅压预偏置。

使用原冻结 UCRT64 Release runner，SHA256 `be486a03d52184b2944e435567c42ecfa643f7fee36a27449aa2bb4f810485f2`；Eigen SparseLU/COLAMD，-O3 -DNDEBUG，gprof 关闭。源和物理输入通过 manifest/plan 的哈希复核。本次仅改变外层实验推进脚本，没有改动 C++、生产配置或求解默认值。

D5 配置保持 300 K、Fermi–Dirac、OldSlotboom BGN、SRH/Auger、constant_field vector-QF HFS 及 contact_node_cell 场回退；没有新增物理平滑或量子、雪崩模型。网格 10,241 节点、19,782 三角形、5,723 个 Si 节点，external-averagebox 输运系数。几何 μm、载流子 m⁻³、电流 A/μm。

原 Poisson/electron/hole 块门限为 5e-8/1e-11/3e-10，carrier-row 和端口 KCL 比率门限为 1e-8；Newton 更新预算 160，内部按需重定心最多 8 次。初始步长 0.0025 V，原增长因子 1.35、重闭合、密度恢复和失败回退流程均保留。

## 联动规则

先对原计划目标检查预测可用性：历史不足、上一步因命中精确点而截短、外推比例超过 2 或处于回退重试时，不使用该预测。预测不可用时，将本次 proposal 限制到 0.1 V，并重新计算目标与匹配的 QF cap；有可信预测时，最大 proposal 为 0.2 V。

成功后从实际采用的 proposal 按原 1.35 倍增长，即主动限制为 0.1 V 后先增长到 0.135、0.18225 V，再到 0.2 V。主动限步单独记录，不计为失败回退；有效 proposal 写入 transfer，避免将主动限步误判成下一次的参考点截短。若降步后重新具备有效预测，可以在较小步长上使用该预测。

## 完整前 8 点结果

| 策略 | 更新 | 推进 | 子进程 | 失败回退 | child wall / s | child CPU / s | 控制器总时间 / s |
|---|---:|---:|---:|---:|---:|---:|---:|
| 前次保护预测器，0.1 V | 436 | 107 | 110 | 0 | 297.18 | 278.86 | 340.00 |
| 前次保护预测器，0.2 V | 619 | 67 | 90 | 6 | 498.75 | 433.08 | 538.17 |
| 本次预测器与步长联动 | 378 | 67 | 78 | 0 | 374.55 | 285.39 | 411.96 |

相对 guard01：Newton 更新变化 -13.30%，child wall 变化 +26.04%，CPU 变化 +2.34%。

相对 guard02：Newton 更新变化 -38.93%，child wall 变化 -24.90%，CPU 变化 -34.10%。

所有时间来自单次串行运行，没有交替 AB 重复。对照为同日稍早的完整实验，环境负载未严格控制，时间变化不具有重复性统计保证。child wall 包含全部子进程启动、失败和恢复；控制器总时间包含外推 CSV、协调和评分。独立性能审计在本次计时求解结束后执行。

## 门限与行为审计

8 个精确点的原块、逐行和 KCL 门限，以及原工程/最终可用曲线门限全部通过。相对 Sentaurus 参考电流误差中位数 1.594856%，P95 1.624118%，低压差分电阻误差 1.032119%；最大 KCL 比率 2.99509e-13。

主动降低 proposal 共 6 次，成功的完整 0.2 V 步共 27 次，密度恢复的成功推进共 5 次。数值分解 491 次，内部重定心 7 次。

| 主动限步父电压 / V | 原计划步长 / V | 有效 proposal / V | 原因 |
|---|---:|---:|---|
| 1.33333333 | 0.2 | 0.1 | previous_step_clipped |
| 2.66666667 | 0.2 | 0.1 | previous_step_clipped |
| 4 | 0.2 | 0.1 | previous_step_clipped |
| 5.33333333 | 0.2 | 0.1 | previous_step_clipped |
| 6.66666667 | 0.2 | 0.1 | previous_step_clipped |
| 8 | 0.2 | 0.1 | previous_step_clipped |

| 阶段 | 子进程 | 失败 | 更新 | child wall / s |
|---|---:|---:|---:|---:|
| initial | 1 | 0 | 0 | 2.36 |
| direct | 67 | 5 | 352 | 326.62 |
| reclose | 5 | 5 | 8 | 23.81 |
| density | 5 | 0 | 18 | 21.77 |

独立 audit.py 核对每个 direct 的计划 guard、主动限步、有效 proposal、目标电压和匹配 QF cap，验证所有无预测的 direct 步均 ≤0.1 V、所有 >0.1 V 的 direct 步均有预测；还核对预测逐节点公式、参考值/增量、原门限、输入/输出哈希和 profile 计数。首个联动变化之前的 28 个子进程与原 0.2 V 对照的更新数、结果及 final_state_hash 一致。所有检查通过。

## 结论与限制

联动验证通过：在 6 个中间精确参考点后主动将 proposal 从 0.2 V 限制为 0.1 V，原 6 次物理回退全部消失，仍保留 27 次成功的完整 0.2 V 推进。相对未联动的 0.2 V 方案，Newton 更新减少 241 次，恰好等于该对照被回退尝试的总更新数。

但本次没有证明它快于带保护的 0.1 V 方案。相比 0.1 V，更新从 436 降到 378，而数值分解从 479 增至 491（+2.51%），CPU 从 278.86 增至 285.39 s（+2.34%），child wall 从 297.18 增至 374.55 s（+26.04%）。墙钟增幅明显大于 CPU 增幅，本次单次顺序测量不能把所有墙钟差异都归因于算法；同时，分解次数表明减少已接受更新并不必然减少线性求解工作。

仍有 5 个 direct 与其 reclose 未通过逐行门限，均在三个全局块已经合格时发生，目标电压为 0.636846、0.836846、1.036846、1.236846 和 1.433333 V；这些点由原密度恢复通过。0.1 V 对照只有 1 次密度恢复。下一项值得验证的是低压启动段继续限制为 0.1 V，或针对这些逐行闭合失败做诊断；本次未执行上述后续修改。

当前仍保留带保护的 0.1 V 作为实测耗时较优的候选，联动策略作为已经完整通过精度和步长规则审计的后续优化基础。

本次没有重跑已有的全量 CTest；该冻结 C++ 二进制在之前已通过 Release 串行 746/746。本次实际完成的是新的完整前 8 点数值验证、独立步长/预测/门限审计和报告检查。实验结果仅适用于当前配置和电压范围，尚未集成到内置扫压流程或生产默认值。

## 证据

忽略目录 `reference_staging/templates_ldmos_predictor_step_link_20260911/` 包含 run_linked.py、linked02/plan.json、progress.json、fixed/ledger.json、各子进程的配置与 profile、audit.py 和 audit.json。二进制 manifest 为 `reference_staging/templates_ldmos_qf_recenter_validation_20260910/binary/manifest.json`。

对照见[带保护预测器及 0.2 V 验证](templates_ldmos_stage4_guarded_predictor_step02_validation_2026-09-11.md)。输出目录必须不存在；从工作树根复现：

```powershell
D:/msys64/ucrt64/bin/python.exe -X utf8 reference_staging/templates_ldmos_predictor_step_link_20260911/run_linked.py --manifest reference_staging/templates_ldmos_qf_recenter_validation_20260910/binary/manifest.json --output reference_staging/linked02_repeat_20260911 --recenter --outer-secant --max-step 0.2
```
