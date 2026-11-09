# BVMethods Quasistationary新模式VM完整扫描（2026-10-30）

用户要求将[新外压步长模式](bvmethods_quasistationary_steps_2026-10-30.md)在vela虚拟机进行一小时完整BV验证。

## 版本与配置

全新Sync源码快照、独立Release构建和独立结果目录。源码包含归一化修复与新Quasistationary外压步长控制，仍有未提交修改；以ignored回执中的源码清单和SHA256为准。不是复用上轮双0.5 V二进制。

核心potential_density、quasi_fermi SG、300 K、JD/OldSlotboom、总杂质迁移率、SRH＋GradQF AV、BTBT关闭、R=1e7 Ω·μm；共同2600节点／5020 Tri3网格。普通Release／实际启用MKL PARDISO单线程，自动导数恢复、局部非正保留1%、浓度下限1e-100 cm⁻³、精确消元、对数人口外推、full_step，原40轮和enforce_joint门。停止电流1.3e-3 A/μm，原各诊断的report/enforce属性保留。

显式coupled_step_control=quasistationary，配置移除旧max_inner_voltage_step_V和coupled_continuation，关闭额外内压预测限制、固定Newton内压修正限制及同目标上一接受态重试。外压初步0.01 V、MinStep=10 μV、MaxStep=2500 V、Increment=1.41、Decrement=2，按实际增广Newton次数调步；不额外设置重试次数截断。器件内压由负载线求解决定。

## 执行和计时

先完成Release构建、八组相关C++回归、Python归档及冷Poisson到0.01 V预检。预检核验实际诊断模式、两项内压限制均关闭、原验收门通过，才开始完整扫描。扫描从冷Poisson重新开始，不从原生或旧Vela接受态续算。

实际完整扫描（含冷Poisson）限3600秒，从status.simulation_started_unix计时；编译、测试、独立冷窗口不计入。远端监督器用单调时钟执行截止，仅终止本轮启动的进程组；SSH断线不取消时限。

初次同步遇SSH连接临时关闭；连接恢复后已完成源码同步和后台Start。源码校验、Release配置、MKL PARDISO配置、构建、八组C++、Python归档及冷窗口全部通过，已进入完整扫描。原两个VM扫描均已结束，本轮为独立新作业。

## 数值失败终态

完整扫描北京时间2026-10-30 22:13:29开始，2026-10-30 22:39:40数值失败结束，耗时1570.36秒（约26分10秒），早于一小时截止23:13:29，非超时。Status已确认监督器退出4，求解器退出1，exit.json同样记录worker退出4。结束后有界读取一次失败日志，明确为耦合Newton更新预算耗尽。

共18个接受点（含零偏），最后接受外压0.8558891781107658 V、器件内压0.8275903757486504 V、电流2.8298845204396072e-9 A/μm，12次更新。未达到停止电流1.3e-3 A/μm，不能给出本轮击穿电压或宣称完整BV通过。

最后失败试探目标0.855906175147395 V，外压增量1.6997036629162743e-5 V（约16.997 μV），retry=0、第40轮、完整导数。器件残差1.0000446060744351e-8仍略高于原1e-8门；负载线残差−3.43260833091108e-8 V、人口RMS=0.6909747255虽满足各自门，联合验收仍不通过。最近记录电子／空穴保护节点均0，不代表整个轨迹未发生保护。按Decrement=2计算下一步约8.499 μV，低于MinStep=10 μV，因此新规则退出，不再额外尝试恰好10 μV的点；此为配置及控制规则的推算，与保存的原始试探严格区分。

最后接受态的实际诊断确认external_step_control=quasistationary、newton_inner_voltage_limit_enabled=false、predicted_inner_step_limit_V=0（表示关闭）。这证明本次使用新规则，但不能仅凭终端残差略超门限就断定剩余非线性失败的唯一成因。没有放宽门限、改变模型或重启。

本地短窗口通过、远端八组C++／Python归档及冷窗口预检通过等结论保持；此前困难路径未闭合的结论也保持。终态回执、配置来源、源码校验清单和有界失败日志已保存，定时监测已停止。私有连接配置、远端源路径和作业标识仅保留在ignored构建目录。
