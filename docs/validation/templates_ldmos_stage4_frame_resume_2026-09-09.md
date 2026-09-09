# LDMOS Vg=8 精度修复版续跑（2026-09-09）

记录时间：19:14:42，Asia/Shanghai。用户要求继续未完成验证，并明确本次不设半小时
限制，继续至完成或出现需要处理的问题。此前暂停结果和证据保持不变。

## 恢复资格

继续使用精度保持修复版冻结程序，SHA256：
`968ef00ed4cc16447c98853b4044f0895238f77f4972e18399694c7546cad97d`。
此前 739/739 CTest、真实网格局部验证和物理/数值门限保持有效；本次未修改求解器，
未重新构建。恢复检查确认原计划的 30 项冻结输入/代码哈希均未变化。

从此前 Vd=4.9999999999999964 V 的合格状态恢复，父状态 SHA256：
`c8afd25d9b0e89a246c9cf43f5af0074051d8b161b5829c2b5da299c9837068e`。
恢复前再次核验 62 个接受状态、4 个精确参考点、122 个已结束子步骤，以及从原始
0 V 起点延续到 5 V 的父状态链；完整性检查通过。

原入口只支持从精确参考点新建曲线，因此新增独立续跑适配器。它引用旧证据的绝对路径，
在新目录保存后续子步骤与账本；既不修改旧 ledger，也不将 5 V 插入精确参考点列表。
原 31 点列表、步长控制、换参考策略和评分函数沿用冻结版本。

暂停时被终止的旧 `child_00122_direct` 单独标为中断，未接受其输出。该子步骤完整
成本未知，因此后续累计 Newton/CPU/墙钟指标是已完成子步骤的合计，不能写成覆盖
所有中断成本的完整性能总量。

## 启动与首步验证

续跑控制器于 **19:13:37** 启动，初始 PID 3692。
5 V 同偏压复核以 **0 次更新**通过原全局块、局部行和 KCL 门限。
截至 19:14:42，已接受到 **5.1000 V**：累计 680 次已记账更新，其中暂停前 667 次、
本次新增 13 次；精确参考点仍为 4/31，回退数 0。进程路径和启动时间已核对。

当前验收范围为恢复后的同一条 Vg=8、0–40 V 曲线。必须完整取得原 31 个精确参考点，
并通过最终评分和完整性审计后才能宣布单栅压通过；不能与旧二进制 Vg=4 曲线拼接为
最终版本双栅压 D5 通过。既有每分钟聊天监测已恢复，旧暂停截止不再适用于本次运行。

## 证据

- [原暂停确认](../../reference_staging/templates_ldmos_frame_fix_20260909/pause_confirmation.json)
- [续跑适配器](../../reference_staging/templates_ldmos_frame_resume_20260909/run_resume.py)
- [原曲线前段完整性审计](../../reference_staging/templates_ldmos_frame_resume_20260909/vg8_full/prefix_audit.json)
- [冻结计划](../../reference_staging/templates_ldmos_frame_resume_20260909/vg8_full/plan.json)
- [当前账本](../../reference_staging/templates_ldmos_frame_resume_20260909/vg8_full/fixed/ledger.json)
- [当前进度](../../reference_staging/templates_ldmos_frame_resume_20260909/vg8_full/progress.json)
- [运行输出](../../reference_staging/templates_ldmos_frame_resume_20260909/resume_stdout.log)
