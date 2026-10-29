# SLOT 5 mV：Sentaurus 与 Vela 收敛条件及实际结果

日期：2026-10-11。仅分析已有记录和保存方向，未新增 Newton 求解、未改验收门槛。

## 结论

原生 double/80 位均在第 4 步由更新误差分支宣布收敛。最新 Vela 消元系统的
double/long-double 限幅开启路径，在第 8 步满足自身相对整体残差条件，第 9 步
满足绝对整体残差条件，第 10 步也满足按原生口径重建的更新误差指标，但仍有
超标载流子行，未通过 Vela 完整资格。

因此两边的“收敛”判定确实不同；不能把本次 Vela 失败概括成整体残差降不下来，
也不能据影子更新指标通过就认定局部误差无害。关闭限幅的两臂在第一次试探时
发生人口定义域异常，尚未进入可比较的接受迭代序列。

## 对象与判据

固定工作点 Vg=5 V、Vd=0.005 V，13881 节点、27275 三角形。原生为
T-2022.03-SP2/PARDISO，分别 normal double 与 extended 80 bit。
Vela 为 `build/s11elimination_newton/runner.exe`，
`sparselu_eliminated` / `sparselu_eliminated_long_double`；后者只扩展线性后端精度，
不是完整装配、状态和非线性流程的 80 位实现。物理采用该批冻结的
AverageBox、Fermi/effective-band SG、PhuMob/HFS、SRH/Auger 和接触电离回退配置。
两端终止结果是各自实际轨迹的比较，不是相同迭代态或相同后端的单因素实验。

| 判据 | Sentaurus 本次有效设置 | Vela 本次有效设置 |
|---|---|---|
| 整体残差 | 原生 RHS 范数 < 1e-10 | 归一化块范数 rho ≤ 1e-9，或 rho/rho0 ≤ 1e-8 |
| 更新误差 | Digits=5，ErrRef(n,p)=1e10 cm^-3；error < 1 可独立终止 | 不是当前配置的独立成功分支；本次仅补算影子指标 |
| 局部载流子行 | 此次未发现与 Vela 同定义、同阈值的逐行残差门槛 | 所有适用载流子行必须满足 eps_row=0.001 |
| 全局连续性闭合 | 未据现有日志证明有相同口径的独立门槛 | 继续要求原相对/绝对闭合条件；相对容差 0.01 |
| 其他约束 | 人口变量、原生边界及求解器保护 | QF/边界等原保护仍生效 |
| 更新与停止策略 | NotDamped=25、Iterations=25；本次 factor 均为 1 | merit 线搜索、QF 0.1/0.05 V 限幅及条件性 Poisson 重修正；最多 40 更新 |

Sentaurus 手册 T-2022.03 第 193–194 页与运行日志一致：默认是 RHS 条件与
更新误差条件的 **OR**；本次未启用 `RhsAndUpdateConvergence`。
`NotDamped` 控制阻尼安排，不是额外的收敛充分条件。

Vela 常规分支为“整体绝对或相对条件”**加上**逐行、闭合及接触 QF 检查。
配置还启用 `carrier_row_qualified_stall_acceptance`，允许满足逐行及其他原条件的
态在 residual floor=1e-8 下退出；该分支仍不能绕过坏行。
源码位置：[NewtonSolver.cpp](../../src/solver/NewtonSolver.cpp) 的
`primaryToleranceConverged` 与 `carrierRowQualifiedStallAcceptance` 分支。

Vela rho 是各变量块 L2 范数按固定块尺度归一化后合成；默认块尺度为
max(初始块范数,1)。行比的尺度是 max(通量绝对值和, |复合净项|, |电离项|,1e-30)。
因此 **原生 RHS 数值与 Vela rho 不能直接相除或套同一个阈值**。

## 实际结束结果

| 路径 | 接受更新数 | 最后整体残差（各自定义） | 更新误差 | 原验收结果 |
|---|---:|---:|---:|---|
| Sentaurus double | 4 | RHS=1.63e-2 | 日志 0.05315 | 更新误差通过，正常收敛；RHS 未通过 |
| Sentaurus 80 bit | 4 | RHS=2.35e-4 | 日志 0.05315 | 更新误差通过，正常收敛；RHS 未通过 |
| Vela 消元 double，限幅开 | 18 | rho=1.65862e-12 | 原始方向影子 E=0.0426143 | 79 电子+1433 空穴坏行，最大行比 1.99298；线搜索拒绝后失败 |
| Vela 消元 long double，限幅开 | 24 | rho=1.10821e-12 | 原始方向影子 E=0.165513 | 58 电子+3226 空穴坏行，最大行比 1.92626；线搜索拒绝后失败 |
| Vela 两种精度，限幅关 | 各 0 | 无新接受末态 | 不适用 | 首次试探抛出 SG 人口必须有限且为正的异常 |

两条限幅开路径的最终电子/空穴闭合比分别为
9.21373e-8 / 1.20736e-7，以及 6.27137e-8 / 5.93949e-8，均已通过原闭合检查。
它们不是预算耗尽退出。新增完整合格点仍为 0。

## 同一 Vela 轨迹上的影子检查

使用原生导出的变量掩码：13881 个 psi、12634 个电子人口、12634 个空穴人口，
总计 39149 项。采用 Digits=5，人口参考 1e10 cm^-3，psi 参考
0.025851995266484913 V。计算

`E = RMS[delta_x / (1e-5*(abs(x_new)+x_ref))]`。

这是依据原生保存方向重建并核验的指标，不是把手册单变量不等式理解成所有节点
的最大值检查。ErrRef 不删除低密度节点的方程。

对于 Vela，使用父态 Fermi 因子 g=F_half'/F_half，把保存的**限幅、线搜索前**
稳定化方向转换为人口线性增量：delta_n=-n*g_n*(d_phin-d_psi)，
delta_p=p*g_p*(d_phip-d_psi)，其中 d 为内部无量纲电势增量。
psi 增量使用该批实际 Vt=0.025851999786435535 V。分母取更新后真实人口。
HDF5 人口由 m^-3 除以 1e6 转为 cm^-3，DOS 也使用 cm^-3。

| 时刻 | double：E / 坏行数 | long double：E / 坏行数 | 含义 |
|---|---:|---:|---|
| 第 8 步 | 203.928 / 1966 | 203.928 / 1361 | 自身相对整体残差首次通过，更新误差不通过 |
| 第 9 步 | 16.9723 / 2004 | 16.9716 / 2154 | 自身绝对整体残差首次通过，更新误差不通过 |
| 第 10 步 | 0.685798 / 2129 | 0.692360 / 1878 | 影子更新误差首次通过，逐行仍失败 |

另同步计算实际有限步人口差指标。末态 double 为 0.0212971；long double 为
2.02300e-5，后者明显受到 alpha=1/1024 缩步影响，不能把它称为未阻尼 Newton
方向很小。原始方向指标用于避免这类混淆；两者均只是诊断，未参与接受。

如果仅采用“Vela 自身整体残差 OR 此影子更新误差”，两臂会在第 8 步停止；
如果仅采用影子更新误差则会在第 10 步停止。这是对既有轨迹的逻辑回放，
**不是已运行的新策略，也不是 Sentaurus 内部 RHS 判据对 Vela 态的直接认证**。

## 核验、证据与后续判断

分析脚本读取两臂共 42 次接受更新及其父态/原始方向，逐步重新计算 Fermi 因子，
同时核对标量逆映射。原生两种精度共 8 次更新的方向误差重建值与日志打印值
相差小于 0.2%，处于打印精度量级。输入文件 SHA-256 前后核对；无原记录修改。
运行环境出现既有 h5py/HDF5 版本警告，所有读取和数值断言正常完成。

- [影子指标 CSV](../../build/s11convergence_comparison/shadow_metrics.csv)
- [分析与输入哈希](../../build/s11convergence_comparison/analysis.json)
- [可重复分析脚本](../../build/s11convergence_comparison/analyze.py)
- [原生 double 日志](../../build/s10native64/native/double64/idvd.log_des.log)
- [原生 80 bit 日志](../../build/s10native64/native/extended80/idvd.log_des.log)
- [本地原生手册文本](../../build/s09native_newton/sdevice_ug.txt)
- [Vela 限幅开启结果](slot_ldmos_eliminated_newton_trial_2026-10-11.md)
- [Vela 限幅关闭结果](slot_ldmos_eliminated_cap_toggle_2026-10-11.md)

本批把“到不了整体小残差”和“整体指标通过但逐行未通过”区分开了。
下一项有意义的验证是对第 10 步等候选态量化端口分载流子电流、源项积分及场量
对剩余坏行的敏感性，再讨论验收是否需要修改；不能用原生已经宣布收敛作为
取消原逐行门槛的依据。限幅关闭的首试探定义域异常应作为另一项更新合法性问题处理。
