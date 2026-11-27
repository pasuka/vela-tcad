# SLOT 消元系统 QF 限幅开关对照

日期：2026-10-11。Vg=5 V、Vd=0.005 V、原点0 V，同一13881节点SLOT网格。

**关闭QF分量限幅未得到新合格点：double与long double均在首个试探残差评估时，
因SG源核人口定义域异常退出，接受更新0次。** 这不是正常完成所有回溯后拒绝，
不能据此断言采用更小阻尼的无分量限幅路径必然不收敛。

## 对照条件

限幅开启基线复用上一批 `build/s11elimination_newton/`，不重复求解；
新增两臂在 `build/s11cap_toggle/`，使用完全相同的冻结runner.exe、
父态 `build/s10initial05/v5_qf.h5`、模型、材料、网格和源项配置。
两种后端分别为 `sparselu_eliminated`、`sparselu_eliminated_long_double`；
后者仅提升共享线性分解/回代精度，装配和非线性状态仍为double。

每种精度的配置差异仅为：

|配置|限幅开|限幅关|
|---|---:|---:|
|quasi_fermi_update_limit_V|0.1 V|0 V|
|quasi_fermi_update_limit_minority_V|0.05 V|0 V|

`max_update` 原本为0，继续关闭全向量上限。原merit回溯线搜索、正则化1e-8、
闭合替换、Fermi、effective_band_sg、AverageBox及体积、接触电离回退均保持。
abstol=1e-9、reltol=1e-8、eps_row=.001、闭合=.01及其它完整资格不变。
每个新增臂最多40次接受更新、软/硬截止1080/1200秒，串行总预算2500秒。

实现上的重要耦合：`componentwise_poisson_recorrect` 仅在QF实际被裁剪时
调用Poisson重修正。将两个限幅值设为0，也就不再触发该重修正。
这是现有“限幅流程开/关”的对照，不能将全部差异独立归因于载流子裁剪，
也不能说两臂都执行了相同的Poisson重修正。

## 结果

|后端|限幅|接受更新数|退出|总超标行|墙钟秒|
|---|---|---:|---|---:|---:|
|double|开启，既有基线|18|载流子行未通过、线搜索拒绝|1512|696.56|
|double|关闭，本次新增|0|SG源核人口定义域异常|未产生新末态|33.20|
|long double|开启，既有基线|24|载流子行未通过、线搜索拒绝|3284|1025.63|
|long double|关闭，本次新增|0|SG源核人口定义域异常|未产生新末态|37.14|

两臂关闭限幅时stderr均为：

```text
vela_example_runner failed: SG source requires finite positive populations
```

新增两次运行合计约70.49秒，均早于预算停止。两个进程退出码均为1，但没有正常返回
结构化Newton结果，status.json为{}；不把其中缺失的iterations/failure_reason填成
求解器实际返回值。接受更新0次由仅保存iteration_0、无后续接受态及退出日志共同核对。
未生成final.h5，父态和首步矩阵/方向保留；原始退出证据未改写。

开启基线虽完成18/24次更新，也均未通过原完整资格。本次四种组合均无新增合格点，
没有将初态或失败试探计入续算历史。

## 首步方向与异常路径

每个开关对照的初态全部数据集、残差、超标行列表一致；physical.bin、
regularized.bin、closure_system.bin、assembled.csv、linear.csv和raw.csv逐文件哈希一致。
关闭限幅时processed.csv与raw.csv也逐文件相同，确认没有隐式保留该分量限幅。

|首步原始方向|double|long double|
|---|---:|---:|
|最大电势位移幅值 V|0.00526864950259|0.00526864950259|
|最大电子QF位移幅值 V|0.00905577646124|0.0101698672841|
|最大空穴QF位移幅值 V|8752.49029221|8752.48962746|

最大空穴方向位于节点10237，符号为负，约−8752.49 V。电压步长5 mV并不意味着
内部载流子QF Newton方向也只有5 mV。两种线性精度都产生了该量级的大方向，
单纯提高线性分解精度没有将其变成可直接施加的有限更新。

实际实现中，LineSearch.cpp先构造candidate，再调用residualFunction(candidate)，
随后才调用hasPositiveFiniteCarriers。PhuMobSgSourceCurrent.inl在对非正或非有限人口
取对数前抛出invalid_argument；本路径未将该异常转为拒绝当前试探并继续减半，
而是传播到runner顶层退出。因此本批没有完成一条正常的无分量限幅回溯路径。
报错不含单元/节点编号，不能把“最大方向位于10237”当作“已定位首个报错人口节点”，
也尚未区分零人口下溢与其它非有限情形。

## 判断和后续

现有分量限幅对保持试探人口处于源核可计算范围具有实际保护作用，但开启限幅
本身仍不足以达到完整资格。简单将限幅值置零，会先暴露试探有效性与异常处理问题。

下一项有辨别力的工作是将人口有效性检查放到源残差求值之前，或只捕获明确的人口
定义域错误并把它转为一次有记录的拒绝试探，然后按原回溯规则减步。
不能吞掉任意异常、修改最终门槛，或把“绕过报错”当作收敛。
完成这类改动并验证后，才能真正比较无分量限幅的回溯轨迹是否有收益。

## 验证与证据

相关quasi-Fermi测试12例103项断言全部通过。未修改源码、编译器、默认配置或旧证据。
配置差异审计、冻结runner及全部输入/旧记录哈希核对、初态字段核对、
初始矩阵和方向绑定均通过。此审计不是独立数值复求。

证据：`build/s11cap_toggle/` 下run.py、configuration_check.json、contract.json、
summary.json、audit.json、cap_tests.log及各臂的stderr、config.log和progress/linear_1。
复用基线见[消元系统5 mV试算](slot_ldmos_eliminated_newton_trial_2026-10-11.md)。
源码位置：[限幅及重修正](../../src/solver/NewtonSolver.cpp)、
[线搜索](../../src/numerics/LineSearch.cpp)、[SG源核](../../include/vela/equation/PhuMobSgSourceCurrent.inl)。
