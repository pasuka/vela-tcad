# SLOT 同外压反馈修正（2026-09-30）

## 已完成

已将[反馈检查](slot_ldmos_prediction_feedback_review_2026-09-30.md)中的同外压
评估接入生产驱动，作为显式选项；旧模式仍为默认。没有修改物理模型、求解器、
候选守恒、严格再闭合、有效载流子行、全局闭合、负载线或步长接受守卫。

开启方式是在新批次原有参数上增加：

```text
--port-step-control --port-predictor limited_rational
--port-feedback-mode matched_actual_outer
```

这些是续算驱动参数，不是物理 deck 字段。不要修改正在运行的证据目录。
新模式及回退协议绑定为 `continuation_policy.port_feedback`，版本
`matched_actual_outer_v1`。改变模式必须使用新目录、重新复验导入父态；
普通 resume 和 `--accept-driver-update` 都不能绕过该身份差异。

## 实现

- `scripts/slot_loadline_step_control.py`：提议记录电阻和参考原点；反馈使用
  事先选定的冻结模型在实际接受外压处预测内压和电流。检查历史、系数一致性、
  扩展区间分母/斜率和参考系，不用当前结果重新拟合或重新选模型。
- `scripts/run_slot_ldmos_bv_campaign.py`：增加 CLI 与身份绑定，将模式和
  实际参考原点传入反馈，模型无效时明确阻止增长阶段晋级。
- 原误差保留在 `original_normalized_port_prediction_error`；新模式下
  `normalized_port_prediction_error` 使用同外压比较。额外记录提议外压、
  实际外压、迁移量、冻结模型、检查跨度、同外压预测和晋级资格。
  原有候选/严格电流偏差、候选守恒、行比与成本记录仍保留。
- 历史不足、系数不一致、有效域越界或参考系不匹配时，用原误差回退，反馈
  增长因子最多为 1；电流增长控制器同时判为非 clean，清除晋级连胜并按
  既有规则缩小提议比例。恢复、成本及真实模型误差仍可进一步缩步。
- 电流误差的归一化分母、成本目标、工程质量线和全部最终资格不变；
  没有把候选电流用作严格电流，也没有简单减去电阻压降来制造小误差。

## 验证

执行以下五个轻量 Python 回归模块，**119/119 通过**：

```powershell
D:/msys64/ucrt64/bin/python.exe -B -m unittest tests.regression.test_slot_loadline_step_control tests.regression.test_run_slot_ldmos_bv_campaign tests.regression.test_slot_bounded_recovery tests.regression.test_run_slot_ldmos_loadline_continuation tests.regression.test_slot_candidate_conservation
```

新增 9 个测试覆盖：同一曲线上外压迁移与真实模型误差的区分、成本/恢复限制、
近竖直支路、缺失历史、参考系切换、扩展域越界、冻结系数不一致、其余质量
条件不得绕过，以及 CLI 参数组合、端到端伪后端接入和模式改变拒绝恢复。
这些是控制逻辑回归，不冒充真实器件求解。

生产代码另在 `build/s30feedback_fix/replay.py` 回放已有十八个接受点：

| 数据批次 | 接受点 | 同外压最大归一化误差 |
| --- | ---: | ---: |
| `s30exp6` | 6 | 0.000280527 |
| `s30growth` | 12 | 0.000220040 |

逐点结果与之前独立离线分解一致。默认旧模式返回结果与各批冻结模块逐字典
相等。复用旧提议时，回放程序只从其绑定父态补入新版本所需的参考原点和
固定电阻元数据，没有改写原始账本。回放输出、源代码 SHA-256 和脚本位于
`build/s30feedback_fix/analysis.json`；输入检查点哈希与前次核验一致。

既有 h5py/HDF5 版本不匹配警告仍出现，测试未失败，未安装或升级依赖。
本次没有 C++ 改动，未重建 Release、未重复全量 CTest、未启动 T470p 或本机
重负载仿真，也没有新增合格偏压点。

## 后续边界

建议新目录从最新已合格态进行短批对照，先用约 2.5 mV 初始提议、保留
10 mV 提议上限，并显式记录初始化的增长阶段与倍率；不能把离线回放的
clean 序列当作已执行的新路径来自动晋级。候选误差经验上限、电流增长限制、
最终资格和独立同内压对照继续保留。只有实际成本、拒绝率和复现结果支持后，
才扩大步长。12.5 mV CLI 初始步与动态上限问题不在本次修正范围。

当前仅完成反馈实现和验证，不声称已证明仿真提速或已跨过击穿电流判据。
