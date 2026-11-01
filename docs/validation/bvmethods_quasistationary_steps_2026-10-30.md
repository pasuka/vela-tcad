# BVMethods 外接电阻Quasistationary步长接入（2026-10-30）

按用户要求对齐Sentaurus外接电阻法的扫描规则。依据本地T-2022.03用户手册第123–124页及BVMethods外接电阻CMD；不将Continuation的MaxVstep混作Quasistationary的MaxStep。

## 行为

新增`external_circuit.coupled_step_control=quasistationary`。只在外压坐标调步：原始CMD的无量纲步长乘整个外压目标跨度，增长由实际增广Newton更新次数和迭代预算决定，失败乘Decrement的倒数；缩步低于MinStep退出。旧控制额外尝试MinStep的行为留在legacy模式。

该模式关闭固定的Newton内压修正上限，以及0.02 V预测内压默认和secant斜率限步；不在失败后先额外以同一外压重试上一接受态。器件内压由完整负载线求解形成。外压步长保持在外压范围内；不从外压MaxStep推导内压上限。

BVMethods CMD发射器默认选择新模式，JSON及sidecar schema均支持显式选择。普通无显式内压保护的coupled JSON也默认新模式；已有显式内压保护的旧JSON自动解析为legacy_inner_guard，保留复现能力。显式选择新模式并同时要求内压限制时提前报错，要求用户清理冲突设置。nested_scalar保持原算法。

没有改变物理模型、人口坐标、浓度下限／局部保护、RHS或Jacobian、精确消元、线性求解器和原验收门。Newton线搜索仍是单独控制，不能把本次步长调整宣称为完整Sentaurus求解器复现。

## 验证

新增线性电阻器件测试：解析内压2 V的一步Newton方向在旧0.1 V限制下不能在两轮预算内收敛，关闭该限制后直接收敛；保持原器件和负载残差判据。

新增策略冲突／兼容解析、MinStep跨越退出、核心外压迭代反馈及CMD默认配置测试。Windows UCRT64普通Release、实际启用MKL PARDISO：

- Release求解器、负载线、DCSweep和输出契约目标编译通过。
- 注册CTest共5项全部通过：dc_sweep、dc_sweep_output_contract、coupled_load_line、pn2d_config_templates、bvmethods_cmd_profile。
- 另运行CMD/JSON发射与实际runner烟雾测试，共48项全部通过，无跳过。
- git diff --check通过；不为本次改动放宽原数值容差。

## 实际BVMethods冷窗口

保持2600节点／5020 Tri3、300 K、JD/OldSlotboom、总杂质迁移率、SRH＋GradQF AV、BTBT关闭、R=1e7 Ω·μm，核心potential_density、quasi_fermi SG、局部1%正值保护、下限1e-100 cm⁻³、自动导数恢复、精确消元和full_step。使用Release/MKL PARDISO单线程，原40轮预算、器件残差1e-8、负载线残差1e-6 V及人口RMS<1联合门。

从自身冷Poisson开始，新增模式的自主外压步进至0.100562182394906 V，含零偏共6个接受点；进程退出0。

| 接受外压V | 更新次数 | 人口RMS |
|---|---:|---:|
| 0 | 0 | 未定义（初态直接接受） |
| 0.01 | 3 | 0.803045 |
| 0.02382666666666667 | 3 | 0.0988453 |
| 0.042944337777777786 | 3 | 0.896441 |
| 0.06937770436740742 | 4 | 0.00585963 |
| 0.100562182394906 | 4 | 0.103617 |

每个点的实际诊断均为external_step_control=quasistationary、newton_inner_voltage_limit_enabled=false、predicted_inner_step_limit_V=0（表示关闭）。最终内压0.07866268268154383 V、电流2.189948308518218e-9 A/μm；所有接受点满足原联合门。

这些偏压是自主步进生成的，使用当前40轮预算反馈，不是冻结的Sentaurus首五目标，不将两者混作同轨迹对照。数据和新模式完整扫描配置保存在ignored的build-release/bvmethods_quasistationary_20261030目录；完整扫描配置已准备但未启动。

本轮仅验证核心规则与冷短窗口；未新跑完整BV或改变既有困难路径结论，尚不能宣称与Sentaurus完整Newton过程、矩阵或全程I–V完全一致。修改尚未提交。


后续VM完整扫描的实际数值失败终态见[独立VM验证记录](bvmethods_vm_quasistationary_full_2026-10-30.md)，不将本地短窗口通过扩展为全程BV资格。
