# SimpleMOS 高场饱和恢复的只读预检

日期：2026-09-12。Enormal 候选已完成八点双初始化和同扰动校准，十六点控制也已通过，正在计算完整曲线。本页保存等待期间的预检和隔离公式测试，没有启动 HFS 仿真或更改候选生产源代码。

原始 Applications Library 展开脚本使用 `Mobility(PhuMob HighFieldSaturation Enormal)`，没有指定驱动力，也没有 `ComputeGradQuasiFermiAtContacts=UseQuasiFermi`。来源是已冻结的 [M8 合同](../../reference_tcad/simplemos_sentaurus2022/simplemos_m8_original_physics_contract_v1.json)，可核对 [n23 原始展开脚本](../../build-release/reference_tcad/simplemos_sentaurus2022/m8_original_physics/sentaurus_bundle/n23/n23_vd_0p05_des.cmd)。

本地 T-2022.03 用户手册第 439、447–448 页规定：该漂移扩散组合默认是 Canali 高场模型及 GradQuasiFermi 驱动力；但**接触相邻单元默认改用电场**。只有显式指定上述 Math 选项才在这些单元也使用准费米势梯度。因此，全域使用边准费米差或全域使用准费米单元梯度，都不能未经核对就当作原始脚本的等价实现。

当前虚拟机 `sdevice -P:Silicon` 导出的 HighFieldDependence 在 300 K 下为：电子/空穴 vsat0=1.0700e7/8.3700e6 cm/s，beta0=1.109/1.213，alpha=0，ku=kv=1。这些速度和指数与当前 Vela 默认值一致；参数一致不替代离散形成次序校准。

Vela 旧代码已有 `phumob_field_lombardi`、准费米驱动力和向量场导数选项，但当前新验证的 `element_distance_gradient` 与 split 算子没有放行 HFS。下一阶段需要分别核对：非接触单元的三角形准费米梯度、接触单元的电场替代、先合成各顶点低场迁移率再做饱和及 box 平均的具体次序，以及全部邻接电势/准费米列的链式导数。不能用旧模型名的语法支持代替本组合验收。

本地安装 SDK 提供 `PMI_HighFieldMobility` 及 `PMI_HighFieldMobility2` 接口，可用于后续原生输入诊断；须先按接口定义构造零影响或严格等价观察器，并实际验证加入前后的电流和状态一致。Enormal 观察器通过不自动证明 HFS 观察器兼容。实际恢复继续以 Enormal 阶段通过为前置条件。

已准备隔离的 Canali 替换式观察器源码；它不是零加法观察器。300 K 原参数的迁移率及对低场迁移率、驱动力和温度的偏导均有实现。独立 90 位中心差分、两档相对步长检查覆盖 144 组参数，包含零场、零迁移率、弱场和高场，最大有限相对差 9.456e-15。温度范围为测试用的 250/300/350 K，不能据此宣称整个器件模型在这些温度下已通过。

首版对零场一律返回零导数，遗漏电子 250 K 下 beta<1 的奇异右导数；其日志和结果保留。v2 明确返回负无穷并检查三个奇异样本，不把它们算成有限光滑导数。当前器件为 300 K，beta>1。汇总和来源见 [v2 检查](../../reference_tcad/simplemos_sentaurus2022/hfs_observer_preflight_20260912/v2/summary.json) 与 [证据](../../reference_tcad/simplemos_sentaurus2022/hfs_observer_preflight_20260912/v2/evidence.json)。本地 SDK 语法检查因缺少 PDE.h 依赖而未完成，日志保留。随后使用虚拟机完整 SDK 和任务内编译选项适配器，原生 CMI 编译退出 0；[构建证据](../../reference_tcad/simplemos_sentaurus2022/hfs_observer_preflight_20260912/v2/native_build/evidence.json)已取回。尚未执行 SDevice 加载、基准/观察器身份或 HFS 自洽验证；编译成功不等于 ABI 和数值身份通过。

观察输出保存每个顶点最后一次单元上下文调用。后续必须与最终原生态及候选单元的输入逐项匹配；不能将该输出直接解释成唯一节点场，也不能为缺失或歧义样本补零。

后续原生输入准备器为 [prepare_simplemos_hfs_native_20260912.py](../../scripts/prepare_simplemos_hfs_native_20260912.py)，目前仅完成语法检查，尚未执行。它要求 Enormal 完整曲线的完成证据，才会生成八点基准/观察器共 16 个任务；先运行 n19、低 Vd、Vg=0.8 V 的一对身份控制，合格后再运行其余 14 个任务。默认接触驱动力切换保持原脚本设置；全接触准费米梯度另作后续诊断，不替代原脚本。

对应的 [原生身份核验器](../../scripts/check_simplemos_hfs_native_identity_20260912.py)已通过命令行入口检查，尚无 HFS 数据可执行数值验收。核验分为归档与输入身份、原生资格、字段导出、基准/替换观察器对照；使用上述冻结合同门槛，完整字段集合必须匹配，失败记录保留。其结果明确分开“替换观察器身份通过”和“最后调用样本及单元形成方式通过”，前者不自动放行后者。

[样本核验器](../../scripts/inspect_simplemos_hfs_observer_samples_20260912.py)也已通过命令行入口检查。它要求原生身份通过，先核对观察器每个顶点的最后调用是否对应最终电势和密度，再保留所有相邻单元的电场、电子/空穴准费米梯度、接触顶点数、两种接触切换假设及低场迁移率预测。字段形成比较仍为诊断，不从最小误差自动选定唯一单元；遗漏样本及原有 PhuMob 空穴截断资格缺口继续保留。尚未执行 HFS 数据检查。

生产实现前需用这些实际输入区分“顶点低场合成→顶点饱和→box 平均”与“先平均后饱和”，并核对单元仅一个接触顶点时的切换。导数检查须覆盖两种载流子的全部单元邻接列：非接触单元的高场驱动力引入对应准费米势梯度链式项，接触单元的电场替代引入电势梯度项，原有 PhuMob 交叉项和 Enormal 法向场项继续保留。零梯度处须验证组合函数的极限，不直接除以梯度范数。单列或编译通过均不能替代这些检查。

在 300 K、alpha=0 下，若原生数据支持逐顶点饱和，则候选链式导数可以明确写成：令低场合成迁移率为 m，z=mF/vsat，A=1+z^beta，有 μ=m A^(-1/beta)，

```
dμ/dx = A^(-1-1/beta) dm/dx
       - (m²/vsat) z^(beta-1) A^(-1-1/beta) dF/dx
```

非接触单元 F=|Σ φ_i grad(N_i)|，非零梯度时 dF/dφ_j 为该向量与 grad(N_j) 点积除以 F；接触电场分支换为 ψ 的梯度。dm/dx 仍包含 PhuMob 密度交叉项及 Enormal 的全部邻接电势项。当前 300 K 的两个 beta 都大于 1，零驱动力处饱和修正的状态一阶导数极限为零，可直接检查组合函数而避免除零。该公式准备不是原生形成次序验证，也不消除 Enormal 绝对法向场在零点的原有非光滑性。
