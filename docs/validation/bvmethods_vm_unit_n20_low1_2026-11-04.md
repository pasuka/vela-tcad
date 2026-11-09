# 全变量Unit坐标的N20冷0–1 V对照

四组有界计算和终态核验已完成。全变量Unit在本N20窗口未显示明显工作量改善，生产默认不切换。关闭适配AND/OR与既有生产N20的IV、updates、attempts、newton及最终HDF datasets逐位一致。

## 受控范围

同一独立诊断二进制、同冻结兼容常数核心、同物理和配置，只切换baseline/unit。主变量psi/V0、n/C0、p/C0、Vinner/V0；C0=5.089999999999982e20 cm⁻³、V0=0.025851995266484913 V、I0=2.98741982845799e-5 A/µm。器件行沿用原归一化，电压列、端口I/I0、R*I0/V0及负载行同步转换。原物理精度和保护门通过反向适配保持。

2600节点/5020 Tri3、300 K、JD/OldSlotboom、总杂质迁移率、SRH+GradQF AV、BTBT关、R1e7、production A0/G2、potential_density/QF SG/floor_adjacent_edges、1%/floor1e-100、精确消元/对数外推/full_step不变。耦合Newton和增长工作量参考N20，Poisson100；原AND或OR、load/KCL/非有限及报告型属性不变。各MKL PARDISO1线程、最多2solver、每组3600秒实际启动计时。没有新Sentaurus、完整BV或生产核心修改。

使用原因保留修正版Unit诊断源，复用既有30项契约断言通过记录与源码SHA，本轮未冒称重跑测试集。旧固定方向门29/30的缺口仍保留。

## 结果

|模式|坐标|接受/尝试/失败|成功尝试报告迭代|失败尝试报告迭代|报告迭代总数|实际写回|
|---|---|---:|---:|---:|---:|---:|
|AND|baseline|20/22/2|114|40|154|151|
|AND|全变量Unit|20/22/2|114|39|153|148|
|OR|baseline|14/14/0|62|0|62|59|
|OR|全变量Unit|14/14/0|62|0|62|59|

四组均solver0、非超时、到1V，所有接受态适用原门通过。报告迭代不等于实际写回或分解次数；正alpha且非rejected/model_review计实际写回，成功Full重试仍计入。

AND两臂接受外压网格和成功迭代序列相同。153对154并不是成功收敛提速，而是目标0.275176304742091 V在Unit中第19轮bootstrap/Full内部后向精度门失败（1.743288e-8/1.105171e-9，门1e-10）；baseline在20轮耗尽。另一个目标0.238036763592897 V两臂均20轮耗尽，RMS仍未过。两次失败后均缩步恢复到1V。

OR接受网格和逐点迭代序列也相同。20个AND共同外压点最大内压差9.16133e-14 V、电流差7.50200e-15 A/µm；14个OR共同外压点对应2.00092e-11 V、7.32532e-15 A/µm。不把共同外压当作逐位相同的内部Newton轨迹。

原始初始packed和适配前RHS相同；人口经n/C0再乘C0后476分量改变，最大相对误差1.28707e-16。Unit激活已由seed坐标证明：初始电子最大值从约5.09e20变为约1，电势packed保持原有psi/V0，不重复缩放。

原物理线性复核次数：AND baseline154/Unit151、OR各61，所有已执行复核均过；这不能覆盖原固定方向门29/30的失败。AND报告型载流子/连续性两臂均3/10未满足；OR baseline8/8、Unit8/7，保留报告型负结果。OR两臂均9residual-only、0update-only、5both、1零更新点。

包含式费用顶层attempt+record+other均闭合，子项不重复相加，不补造逐attempt成本或分解次数。共享VM及原物理Jacobian复核插桩的墙钟不作为独占性能结论。

## 证据与范围

精确监督器/worker0、无自有进程，6155源文件/core/生产binary/material最终SHA不变。95文件一次归档逐SHA核验；ZIP2660796字节，SHA256为5bea67a9351f4f3f194f890463eab9e116c0832103525f6e6d14b4f11c735f50。原物理输入已打包闭合，成功及数值日志未读。详细结果在ignored任务analysis.json、final_summary.json及accepted_point_iterations.csv。

本结果仅0–1V，不覆盖10–100V。用户随后单独授权N20全变量Unit冷0–100V及10–100V迭代比较，使用独立任务目录继续。原生对照仍为已有T-2022.03-SP2归档，不是新原生运行；本Vela常数profile为sentaurus_u2022_12。
