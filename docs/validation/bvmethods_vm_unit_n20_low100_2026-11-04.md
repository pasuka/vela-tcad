# N20全变量Unit冷0–100 V及10–100 V终态比较

同一已核验二进制，AND/OR各含关闭适配与全变量Unit配对。耦合Newton及步长工作量参考N20，Poisson100；最多2单线程solver。每组实际启动14400秒安全时限。原物理、保护、独立门保持；不是完整BV，不运行新Sentaurus。

## (10,100] V迭代比较

|组|到100 V|接受/尝试/失败|成功轮数|失败轮数|总报告轮数|实际写回|
|---|---|---:|---:|---:|---:|---:|
|AND_baseline|True|70/91/21|418|279|697|626|
|AND_unit|True|79/103/24|489|230|719|666|
|OR_baseline|True|18/19/1|63|17|80|77|
|OR_unit|True|70/92/22|407|219|626|579|

未覆盖100 V的组不能作为完整区间收益。报告轮数不是分解次数；实际写回与报告迭代口径不同。强制100 V端点贡献见JSON。若路径分叉，累计轮数不等于同态Newton方向单因素比较。

历史参考：T-2022.03-SP2原生(10,100]V为12/12/0、49轮；兼容常数生产A0 N40 AND为72/100/28、954轮，OR为71/98/27、781轮。跨版本/预算/路径不能作唯一因果或独占加速结论。

## 终态资格与限制

AND_baseline: solver=0, timeout=False, 原适用门=True, report载流子/连续性未满足=28/10, 原物理线性复核失败=0, OR分支=None
AND_unit: solver=0, timeout=False, 原适用门=True, report载流子/连续性未满足=50/10, 原物理线性复核失败=10, OR分支=None
OR_baseline: solver=0, timeout=False, 原适用门=True, report载流子/连续性未满足=58/8, 原物理线性复核失败=0, OR分支={'residual_only': 46, 'update_only': 1, 'both': 24, 'zero_update_points': 1, 'failed_points': 0, 'carrier_report_unsatisfied_points': 58, 'continuity_report_unsatisfied_points': 8, 'max_device_residual': 1.0966094808309227e-08, 'max_load_residual_V': 3.3210623229251723e-07, 'max_valid_update_rms': 5622.597534593712}
OR_unit: solver=0, timeout=False, 原适用门=True, report载流子/连续性未满足=40/7, 原物理线性复核失败=3, OR分支={'residual_only': 70, 'update_only': 0, 'both': 44, 'zero_update_points': 1, 'failed_points': 0, 'carrier_report_unsatisfied_points': 40, 'continuity_report_unsatisfied_points': 7, 'max_device_residual': 9.960458647058734e-09, 'max_load_residual_V': 3.694451944369348e-07, 'max_valid_update_rms': 5622.621115041751}

全变量为psi/V0、n/C0、p/C0、Vinner/V0；Unit为隔离诊断原型，既有固定方向门29/30限制保留，production_qualified=false，不切换默认。成功日志未读取，数值失败保留，不自动重跑或调整门/预算/floor。

详细数据：analysis.json、region10_100_iteration_comparison.json、region10_100_accepted_point_iterations.csv、region10_100_attempt_iterations.csv。最终源/core/binary/material身份、精确求解监督器退出及无自有进程由terminal_qualification.json给出；一次终态包SHA与全成员校验由auto_completion.json给出。


## 已核验结果及解释

四组全部solver0、非超时、到100V，全部接受态逐点通过适用原AND或OR及独立门。原求解监督器0、无自有求解进程；自动收尾于北京时间2026-11-04 18:07:15完成，精确收尾监督器0已核对。本机离线不影响VM执行。

(10,100]V的AND报告迭代697→719，增加3.16%；OR80→626，为基线7.825倍。两种规则均未显示本次直接Unit适配减少工作量。扫描路径已经分叉，不能把累计轮数差唯一归为同态Newton方向差或宣称Unit缩放本身无效。

|全0–100V|接受/尝试/失败|报告迭代|
|---|---:|---:|
|AND基线|119/149/30|1194|
|AND Unit|162/206/44|1581|
|OR基线|71/83/12|610|
|OR Unit|114/144/30|1024|

各组强制100V端点均贡献1接受/1成功尝试，轮数依表顺序6/5/3/6；剔除端点后的10–100V轮数为691/714/77/620，结论不变。原物理线性精度复核的拒绝候选数按上述顺序为0/10/0/3，不等于接受态门失败。报告型载流子/连续性未满足数依次28/10、50/10、58/8、40/7，保留report属性。OR基线分支46 residual-only/1 update-only/24 both，OR Unit70/0/44，各含1零更新起点。

费用attempt+record+other顶层均闭合。报告迭代、实际写回和分解次数不同；缺失逐attempt耗时不补造，共享VM墙钟不作为独占提速结论。

6155冻结源文件、生产core/binary/material、复用Unit binary及引用物理输入全部核验；105个数据/源码成员已在VM及本地逐SHA验证，ZIP3132446字节，SHA256为`0200ec5d7f6af705ec9b1f5cfe49ed25afc91f219d2c546f80483eadeb7cd514`。成功日志未读取，无重新求解或放宽门。完整分析、逐点/逐尝试迭代CSV与身份回执保存在ignored任务目录。

本轮任务完成，monitor_finished=true、production_qualified=false；没有遗留等待任务。之前的[Unit固定态验证](bvmethods_vm_unit_basis_2026-11-04.md)29/30方向门限制仍然保留，不改变生产默认。
