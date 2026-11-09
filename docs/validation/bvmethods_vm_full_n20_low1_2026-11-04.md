# 关闭bootstrap的N20冷0–1 V验证

本轮计算、终态独立资格与证据核验已完成。关闭bootstrap后，0.066819 V从4轮恢复3轮，与原生一致；但AND全段成本增加，OR未到1V而数值失败，不能推广为全程始终Full更好。

## 为什么会使用bootstrap

C++读取`solver.population_derivative_recovery`的默认值是`full`。此前生产A0/G2实验显式设为`automatic`，不是通用C++默认。自动控制器默认从bootstrap起步，各新目标按原路径重置；默认bootstrapProbeInterval=2，故两次bootstrap实际写回后作Full周期探测。F/load物理门通过或bootstrap线性/模型失败也可提前探测；RMS不直接决定模式选择。Full需要若干可靠观测确认，确认失败可冷却退回bootstrap。

bootstrap保留完整物理RHS及SRH导数，冻结迁移率导数、忽略AV导数，并沿用既定下限相邻边近似处理。它是用于自动恢复的近似Newton路径，不承诺比Full省迭代。来源：`DerivativeRecovery.h`、`CoupledLoadLine.h`、`NewtonSolver.cpp`、`BvMethodsCoupledRuntime.cpp`与当前配置文档。

## 本轮受控范围

Vela VM复用已核验兼容生产核心及二进制，无新构建。2600节点/5020 Tri3、300K/JD/OldSlotboom、总杂质迁移率、SRH+GradQF AV、BTBT关、R1e7、potential_density/QF SG/floor_adjacent_edges、1%/floor1e-100、精确消元/对数外推/full_step保持。生产坐标，不叠加Unit原型。耦合Newton/增长工作量预算N20，Poisson100，原AND/enforce_joint或OR/enforce及独立load/KCL/非有限和report属性保持。最多两solver各PARDISO1线程、每组3600秒自实际启动安全时限。没有新Sentaurus或完整BV。

相对合格automatic N20，只将`population_derivative_recovery`改为`full`。该既有模式停止整个automatic控制器，不仅是替换某一个冻结矩阵；guard选项与报告插桩未另作切换。逐字段差异和所有trace均为Full已断言；冷初始packed/RHS逐SHA与automatic N20相同。不同模式后续预测态不保证相同。

## 0.066819 V逐轮结果

|完成更新|原生error，T-2022.03-SP2参考|automatic N20 RMS|Full-only N20 RMS|
|---|---:|---:|---:|
|1|约1740|1707.8129|1736.0756|
|2|约53.0|67.3551|53.0291|
|3|0.075527|9.33137|0.0755525|
|4|已停止|0.00366691|已停止|

Full-only第三步后F=3.15163e-11、load=-2.25077e-8V，原AND/OR均过；两规则均3轮。其后目标0.100562182394906V也与原生一致。此结果支持bootstrap/Full路径影响该点收敛轮数，但由于从冷态即全程Full，前面接受态/预测态已有微小差异，不能冒称已完成同一预测packed态只改前两次矩阵的严格因果实验。原生版本与Vela U-2022.12兼容常数边界保留。

## 全段终态

|规则/导数模式|接受/尝试/失败|成功轮数|失败轮数|总报告迭代|末接受外压V|
|---|---:|---:|---:|---:|---:|
|AND automatic参考|20/22/2|114|40|154|1|
|AND Full-only|36/47/11|133|166|299|1|
|OR automatic参考|14/14/0|62|0|62|1|
|OR Full-only|16/32/16|92|112|204|0.262404062077|

AND Full-only solver0、非超时、557.704秒；4次预算失败、7次线性资格失败，恢复后到1V。299较154增加94.2%。OR solver1、非超时、386.626秒；1次预算失败、15次线性资格失败，未完成全段，不与完整窗口作成本收益比较。

OR末目标0.262415002612 V、步长10.9405355uV，第7轮分量后向误差1.46445697e-5超过1e-10；F=1.69321e-7、RMS33114.3、load=-5.40183e-10V。下一半步5.47027uV低于Min10uV，数值失败结束。终止alpha0不是此前未更新，不将它归类为20轮耗尽。

36个AND和16个OR接受态均通过适用原门。report载流子/连续性未满足点分别14/15和11/15，保留report属性。OR分支11 residual-only、0 update-only、5 both，含1零更新起点。全程报告迭代不等于分解次数；缺失逐attempt成本不补造，费用attempt+record+other闭合，共享VM墙钟不作为独占加速结论。

## 证据与准备负结果

观察worker/精确监督器退出0，无自有worker/solver；worker0表示观察完成，不代表OR数值通过。6155源文件、原生产core/binary/material与引用物理输入最终SHA不变。46数据/源码成员一次归档逐SHA通过，ZIP1075941字节，SHA256为`9a0294faaffb05c7b4a0debdfa5cce2ae547f9290b93d8e7a38ef5485f7fe36c`。成功日志未读，OR失败日志仅结束后有界一次，marker保留。首次Windows/Linux远端路径格式准备失败发生在Upload，未调用Start、0次求解；修正仅启动路径，新独立目录，保留准备回执。

本任务结束，无待运行队列；未修改生产默认，不重启OR、不放宽任何门或预算。完整配置、逐轮数据和回执保存在ignored任务目录。
