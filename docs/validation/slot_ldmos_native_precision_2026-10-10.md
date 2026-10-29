# SLOT IdVd 原生 64 位与 80 位精度对照

## 结果

在独立 Sentaurus 配置中把 `ExtendedPrecision` 改为 `-ExtendedPrecision`，
已完成 Vg=5 V、Vd=0.005 V 的双精度求解、导出和成对复核。
64 位和 80 位均在 **4 次 Newton 更新后**按原生更新误差条件收敛，退出码均为 0。
漏电流相对差为 **2.5309955e-11**；64 位最终 RHS 约为 80 位的 **69.36 倍**。

因此，这个低漏压点不需要 80 位才能获得原生认可的收敛结果；但精度影响残差、
近零电流及部分局部源项，不能由端口电流接近推断全部弱行也一致。
两臂都未达到 `RhsMin=1e-10`，原生更新误差停止不能替代 Vela 完整资格。
本批没有新增 Vela 合格点，没有恢复漏压长扫描。

前置：[线性后端和 Newton 方法对照](slot_ldmos_linear_backends_newton_review_2026-10-10.md)。
证据：[build/s10native64](../../build/s10native64/)，
[分析数据](../../build/s10native64/analysis.json)，
[全部场量对照](../../build/s10native64/field_comparison.csv)。

## 配置与单因素控制

使用 Sentaurus T-2022.03-SP2、PARDISO、5 线程；300 K，原点 0 V。
网格为原 SLOT 的 13,881 节点、27,275 三角形，二维单位宽度电流为 A/µm。
从同一 `slot_idvd_20261005/nidvd_VdLin` 零漏压保存态 Load；
它是已有保存态，不是本批重新建立的冷平衡。

沿用最近矩阵比较的物理配置：AverageBox、Fermi、两材料 `DensityIntegral(240)`、
定义积分 Bennett–Wilson 统计修正 PMI、Si PhuMob/Lombardi/HFS、SRH/Auger、
GradQuasiFermi 电离以及现有接触设置。没有关闭源项或更改模型参数。

|设置|64 位臂|80 位臂|
|---|---|---|
|Math 精度选项|`-ExtendedPrecision`|`ExtendedPrecision`|
|日志确认|64 bit double normal precision|80 bit long double extended precision|
|线性求解|PARDISO|PARDISO|
|Iterations / NotDamped|25 / 25|25 / 25|
|Digits / ErrRef(e,h)|5 / 1e10 cm⁻³|相同|
|RhsMin|1e-10|相同|
|Extrapolate|关闭|关闭|
|目标跨度|0→精确 5 mV，一步|相同|
|上限|1200 s，25 次迭代|相同，串行运行|

两份 deck 除精度关键词外逐字一致，par 文件逐字一致。相对于旧 `Iterations=0`
矩阵观测，两臂共同改为 25，并共同加入最终 Plot/Save，才能比较实际收敛。
默认 LineSearchDamping=1，没有额外启用线搜索。每步实际 factor 均为 1。

T-2022.03 用户手册印刷页 223–224 明确：`-ExtendedPrecision` 为默认 double 64 位；
`ExtendedPrecision` 为 Linux long double 80 位。手册也明确扩展精度不保证能解决
一般收敛困难。依据：[本地手册](../../build/s09native_newton/sdevice_ug.txt)。
本次现有标准 BGN PMI 接口本来就是 binary64，两臂使用同一二进制、同一公式，
没有重编插件，也不把 80 位原生模式称为 PMI 全链路 80 位。

## 收敛和端口结果

|指标|64 位|80 位|
|---|---:|---:|
|Newton 接受更新数|4|4|
|最终原生 RHS 日志值|1.63e-2|2.35e-4|
|最终更新误差|5.3150e-2|5.3150e-2|
|停止原因|Error smaller than 1|Error smaller than 1|
|漏极总电流，A/µm|8.04268583993621e-7|8.04268583973265e-7|
|源端空穴电流，A/µm|−2.73567181094911e-10|−2.73567181153941e-10|
|栅极总电流，A/µm|−3.34056885164022e-16|−3.43471525732364e-19|
|程序报告 wallclock，s|91.61|182.64|

两次 wallclock 是本次单样本、含诊断导出的运行成本，不是通用性能基准。
原生内部迭代计时与 wallclock 口径不同，未把它们相加。

|接受更新|64 位 RHS|80 位 RHS|64 位更新误差|80 位更新误差|
|---|---:|---:|---:|---:|
|1|8.73e7|8.73e7|6.10e3|6.10e3|
|2|4.82e4|4.82e4|3.79e2|4.19e2|
|3|1.72e-2|1.69e-2|3.57e2|3.57e2|
|4|1.63e-2|2.35e-4|5.31e-2|5.31e-2|

以上 RHS 与更新误差来自原生日志的有限打印位数，不是 Vela abstol 或逐行 ratio。
没有继续强制迭代探测极限，因此较大的最终 RHS 不能直接命名为已测得的舍入底限。
逐步记录：[64 位](../../build/s10native64/double64_iterations.csv)、
[80 位](../../build/s10native64/extended80_iterations.csv)。

## 第 0 次矩阵和 RHS

两臂均导出 4 套线性系统；维数 39,149×39,149，按原生 ψ/电子人口/空穴人口
块与节点编号对齐。下表比较第 0 次，分母为 80 位对应块范数。
没有行归一化、坐标转换或稳定化再加工，包含边界行；这里不是与 Vela 的新矩阵比较。

|块|Jacobian 相对 Frobenius 差|RHS1 相对 L2 差|SOL 相对 L2 差|矩阵行绝对差最大节点|RHS 绝对差最大节点|
|---|---:|---:|---:|---:|---:|
|ψ|6.9471e-16|2.7064e-10|1.6302e-8|5619|10236|
|电子|1.2065e-9|8.6802e-13|5.4571e-8|10849|5619|
|空穴|2.2410e-8|3.3151e-4|7.1776e-8|10849|13349|

同一个 Load 文件并不意味着两种算术模式重建的 ψ/n/p 逐位相同。
第 0 次原生内部变量相对 L2 差为 ψ 2.50e-15、电子 2.42e-13、空穴 3.10e-13。
因此本对照包含精度切换在装配及状态重建中的作用，不能当成“同一固定矩阵只换线性后端”。
debug 文本本身的打印位数及本地 binary64 解析限制了可见差异，也不是全部 80 位尾数。

## 最终场量

Si 区域按相同节点编号对齐；电场表列的是分量最大绝对差。

|场量|最大绝对差|所在节点|
|---|---:|---:|
|电势，V|6.6058e-14|10344|
|电子 QF，V|1.6507e-11|2854|
|空穴 QF，V|1.5729e-8|9628|
|电场，V/cm|4.8254e-8|10147|

电子/空穴人口的 Si 区域相对 L2 差为 1.5593e-15 / 1.1787e-15；
逐节点最大相对差分别为 6.3922e-10 / 6.0844e-7。
总电离 Plot 场的相对 L2 差为 2.5298e-14，但局部小量仍需单列。
例如节点3598的 SRH Plot 值从 80 位的 −5.5820e7 变为64位的 1.6119e13 cm⁻³s⁻¹；
整体 SRH 场相对 L2 差仍仅 1.4747e-11。不能用整体范数掩盖这种局部相消敏感性，
也不能不经支撑体积和方程分解，就把 Plot 场差直接当成行残差差。
近零 QF、源项或接触电流的相对差尤其容易失真，应同时查看绝对值。

## 导出和审计

- [64 位 deck](../../build/s10native64/native/double64/idvd.cmd) 与
  [原生日志](../../build/s10native64/native/double64/idvd.log_des.log)。
- [64 位最终 TDR](../../build/s10native64/native/double64/final_des.tdr)、
  [最终场量 CSV](../../build/s10native64/native/double64/final_des_export/fields/)、
  [IV](../../build/s10native64/native/double64/idvd.plt)。
- [64 位矩阵和向量 NPZ](../../build/s10native64/native/double64/matrices/)，
  完整原生 debug、每次 Newton TDR 和 final SAV 也已保留。
- [80 位对照](../../build/s10native64/native/extended80/) 使用相同导出结构。
- [完整下载包](../../build/s10native64/results.tgz)，42 个原始结果文件逐一核对远端 SHA256。
  原保存态、circuit SAV 和 PMI 二进制前后哈希相同；两臂输入 deck/par 与上传内容一致。

独立远端目录为 `slot_native64_20261010`；两臂顺序运行，均在上限前完成，
结束后已只读核对无 sdevice 进程。网格 SHA256 为
`ec6cb57f5daa577094c70a5fc3909a2e69a4a6adb69d43c165edb3377fc5258b`。
程序、配置和下载包身份见 [contract](../../build/s10native64/contract.json) 与
[远端收尾核对](../../build/s10native64/remote_final_audit.txt)。

此后比较 Vela 时可以使用本次64位原生参考，减少精度模式这个差异项。
仍须分别审计原生末态的 Vela 资格、同态离散残差及真实更新策略；本批没有证明
长漏压路径、BVDS 或所有弱行对精度不敏感，也没有更改原有验收门槛。
