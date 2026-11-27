# SLOT IdVd 显式门压初始化集成与全路径复验（2026-10-05）

## 集成内容

承接 [Fermi Poly 初值验证](slot_ldmos_idvd_fermi_poly_initialization_2026-10-05.md)，
在 [IdVd 驱动](../../scripts/run_slot_ldmos_idvd.py) 增加显式
`--gate-initialization fermi_poly_poisson`，默认仍为 `none`。
每个门压点从上一完整合格态构造孤立多晶硅刚性偏压预测，执行
`poisson_only` / `preserve_initial`，再执行完整 Newton 验收。
只有完整验收态进记录；Poisson 状态不继承资格，严格失败后的缩步仍用上一合格父态。
进入漏压段不使用该平衡初始化，保留完整输运 Newton。

没有修改物理模型、C++ 求解器、原行比/闭合门槛或旧证据。
全部子阶段计入同一预算，当前阶段在 summary.active 中可读；
初始化失败与完整资格失败保留各自证据。helper、配置、网格、材料、原生曲线及
runner 的哈希入档，辅助 Python 源码复制到本轮新目录作版本留档。

## 检查

7 项 Python 回归通过：驱动控制流 2 项、Fermi 预测不变量 2 项、IdVd 配置 3 项。
新增控制流检查覆盖：Poisson 输出不可进合格记录，严格失败后不复用被拒态，
漏压阶段不走门压初始化。mock 控制流不当作数值资格。
配置与求解器同前批：300 K、Fermi/Bennett、Si PhuMob/Lombardi/HFS/雪崩、
Poly 常数迁移率，mesh_default 元素盒几何；原生几何等价性及高补偿 BGN 差异仍待验证。

## 全路径运行

新目录：[build/s05idvd_path](../../build/s05idvd_path/)。
本机于北京时间 2026-10-05 22:07:58 启动，控制器初始 PID 16336；
复验预算 10,800 s，普通子求解器上限 600 s，门压 Poisson 上限 240 s。
预算用控制器 monotonic 时钟累计，不自动延长或重启。

路径：新冷平衡 → Vd=0 时 Vg 以 .1 V 步进至 5 V → Vg=5 V 下跟随
归档原生 Vd=0–60 V 样本（必要时插入不超过 1 V 的中间点）。失败最多三层二分。
Fermi helper 的单次门压预测上限仍为 .1 V。
只有 gate 段全部通过后才允许 drain 段开始，未完成前不报告原生 IdVd 定量对比通过。

```powershell
D:/msys64/ucrt64/bin/python.exe -B scripts/run_slot_ldmos_idvd.py --config configs/regression/slot_ldmos_idvd_fermi_control.json --runner build/s05idvd_impl/cold_run3/runner.exe --output build/s05idvd_path --native build/s05idvd/native_drain_path.csv --budget 10800 --child-cap 600 --gate-initialization fermi_poly_poisson
```

此命令是本次执行留档；输出目录已存在，不可原地重复执行。
完整状态及停止原因以 [summary.json](../../build/s05idvd_path/summary.json) 为准。
[status.py](../../build/s05idvd_path/status.py) 使用 Windows FILE_SHARE_DELETE 读取活动文件。
[audit_records.py](../../build/s05idvd_path/audit_records.py) 只核验已保存的门槛记录、父态链和
文件哈希，不代替独立重新求解。

## 结果

北京时间 22:15 的快照：新冷平衡原完整资格通过（278.83 s），
CLI 显式初始化的首个 Vg=.1 V 点也通过；已记录 1 个冷态及 1 个门压态，
正在执行 .2 V 的 Poisson 子阶段。门压阶段 Vd=0，漏电流处于数值零水平，
不据此评价开态输运误差。

**运行中，尚不能给出门压全路径或 IdVd 曲线通过结论。**
后续现场状态以 summary.json 为准，不能把此快照视作后续全部点已经通过。
