# BVMethods人口坐标下物理模型的算术与导数审计，2026-11-02

## 结果和适用范围

除已证实的SRH中间指数溢出外，本轮还复现了AV系数的极小正场AD非有限边界、通用复合限幅区的值/导数不一致，以及派生QF相减丢失微弱SG通量的受控例子。它们的启用路径、触发条件和证据等级不同，不能合并归因于当前BV失败。

本轮是源码审计和本地隔离函数诊断，未修生产源码，未更改/重启VM中的30组窗口，未运行新的Sentaurus任务。当前[窗口消融](bvmethods_bootstrap_trust_windows_2026-11-02.md)仍使用原冻结数值代码。

采用Windows MSYS2 UCRT64 GCC 16.2、C++20、-O2、不启用fast-math，直接编译当前头文件及RecombinationModel、CarrierStatistics、BandToBandTunnelingModel、PerformanceProfiler和Material实现。27项函数输入、3项SG坐标控制已执行并验证预期现象；这不是全模型或全器件回归资格。程序、JSON结果及逐文件SHA256保存在ignored的`build-release/bvmethods_physics_coordinate_audit_20261102/`。

| 项 | 当前BVMethods窗口 | 结果边界 |
|---|---|---|
| SRH | 启用 | 冻结坏态中已确认expm1溢出和inf×0；本轮未修 |
| Auger | 未启用；本profile显式拒绝 | 通用装配路径存在广义excess截断和限幅导数一致性问题，非当前BV根因 |
| 碰撞电离AV | 启用 | 电势依赖真实；人工极小正场复现值0但AD NaN，既有两份失败态AV矩阵仍全部有限 |
| SG输运 | 启用 | 大delta已有避免正指数溢出的分支；另有派生QF相减消减风险 |
| GradQF高场/Lombardi迁移率 | 启用 | 梯度依赖真实；所测Canali样本及既有两坏态输运矩阵有限，不构成全域证明 |
| BTBT | 关闭；本profile拒绝 | 不参与当前窗口；本轮仅查局部E2源代码，未资格化NonLocalPath |

## 1. 物理依赖与坐标依赖必须区分

固定温度/材料/寿命/陷阱参数、固定独立n,p时，局部SRH和通常的Auger源项对ψ的偏导应为0。人口模式中派生`fn=ψ+Cn(n)`、`fp=ψ+Cp(p)`，本地`fp-fn`的共同ψ方向抵消。

碰撞电离和迁移率不同。GradQF驱动包含`grad(fn)=grad(ψ)+grad(Cn(n))`，相邻节点的ψ差是真实驱动。AV还经电流、场依赖迁移率与电离系数依赖ψ；SG电流同样必须保留电势列。不能因SRH抵消问题清零这些模型的ψ偏导。统一增加所有节点的常数电势与单独改变一个节点也是两种不同扰动。

## 2. Auger与通用复合：没有当前人口BV路径，但有共享算术风险

[BvMethodsCoupledRuntime.cpp](../../src/equation/BvMethodsCoupledRuntime.cpp)构造时显式拒绝`augerEnabled()`或BTBT。人口坐标目前又要求该显式算子，因此不能称本轮窗口在计算Auger，也不能把以下通用路径现象当作当前窗口根因。

通用[CoupledDDAssembler.cpp](../../src/equation/CoupledDDAssembler.cpp)的`nodeRecombinationRate`在`generalized_fermi`下复用`generalizedSrhCarrierState`的excess，在`classical_np`下走浓度乘积。前者在[CarrierStatistics.cpp](../../src/physics/CarrierStatistics.cpp)把log(gamma)和准费米分裂分别截断到±500，防止裸指数溢出，但极端态中已不保持浓度恒等式。

用原两个失败节点的冻结参数，把浓度/DOS/ni从cm⁻³统一换成m⁻³，直接调用共用函数：

| 冻结节点 | 共用函数excess，m⁻⁶ | 未独立截断的浓度恒等式，m⁻⁶ | 比值 |
|---|---:|---:|---:|
| 689 | 6.545760695e-45 | 2.028527407e58 | 3.226853467e-103 |
| 331 | 1.677137129e33 | 5.831240869e45 | 2.876123910e-13 |

恒等式为`np-gamma_n gamma_p ni²`，不是把广义模型擅自换成经典`np-ni²`。这些是代入共享辅助函数的反事实诊断，不是原生产失败时启用了Auger。中间量逐个截断可以留下有限值，同时丢失应有的指数补偿。

[RecombinationModel.cpp](../../src/physics/RecombinationModel.cpp)还在n/p的1e30与excess的1e60内部数值限幅之外，继续报告未体现限幅的导数。用SI输入、远离限幅拐点的中心差分：

| 接口及固定量 | 报告偏导 | 差分偏导 |
|---|---:|---:|
| Auger，excess=1e61，n=p=1e22，求excess偏导 | 3.928e-21 | 0 |
| Auger，n=2e30、p=1e20、固定excess=1e40，求n偏导 | 2.9e-3 | 0 |
| SRH，excess=1e61、固定其余量，求excess偏导 | 7.6923e-18 | 0 |

上述数值是接口偏导，不是混合改变n和excess的总导数。它证明报告切线不是饱和值函数的普通导数；如果这是有意的近似/延拓，应明确契约，不能默认为真实完整Jacobian。本轮没有直接决定限幅区必须改成零导数，也没有用这些极端控制证明普通工作范围内Auger已错误。

## 3. 碰撞电离：复现另一种0×inf

当前AV使用[BvMethodsTri3Avalanche.h](../../include/vela/physics/experimental/BvMethodsTri3Avalanche.h)的`avCoefficient`，基本形式为`A exp(-B/E)`。对于有效正A/B和人工有限输入E=1e-160、1e-200 V/cm，系数值为0，而一个AD方向为NaN。

原因是AD先求`-B/E`的导数，`B/E²`超过可表示范围；随后指数已下溢为0，执行`0×inf`得到NaN。它不同于SRH的“大正指数溢出后乘本应抵消的0”，但都属于可避免中间范围损失。数学上该极小正场的系数及其普通导数均极接近0。

测得E=0及1e-140、1e-20、1、1e4、4e5、1e6 V/cm的样本值和AD均有限。精确0场有提前返回，极小正场则进入另一条分支。参数使用有效正Van Overstraeten系数，不声称是当前归档器件的逐字系数配置。

这是人工边界反例，**尚未在真实窗口失败态中复现**。既有两份NaN坏态分别有6/12个Full及bootstrap非有限条目，但独立装配的AV、完整输运与冻结迁移率输运均为0个非有限；不能把它们重新归到AV。

零电流范数的AD另外采用零广义导数约定，当前源码对此有说明；这是非光滑点的线性化契约，不能与上述可避免的NaN混为一谈。

## 4. SG：指数分支已有保护，QF差分仍可能丢失微弱量

`avSgFlux`在abs(delta)>50时改用成对Bernoulli通量，避免直接expm1(delta)的大正指数。对delta=−800、−100、−50、0、50、100、800，等浓度/单位系数控制均返回有限通量及正确的单位delta导数，未复现SRH式正指数NaN。

另外，用非简并极限的受控节点对调用同一QF通量核与`populationSgFlux`：n1=1e10、n2=1e10*(1+1e-14) cm⁻³，ni=1e10、VT=0.0258519997864 V、两节点ψ相同，派生fn=ψ−VT log(n/ni)。去掉公共几何/迁移率系数后：

| ψ共同背景，V | QF核通量 | 直接人口核通量 |
|---|---:|---:|
| 0 | −9.992007222e-5 | −9.918212891e-5 |
| 6.6 | 0 | −9.918212891e-5 |
| 10000 | 0 | −9.918212891e-5 |

在6.6 V处两个派生fn舍入为完全相同的double，微小差被共同背景吞掉。直接人口核保留表示出来的n1−n2。这里是**非简并极限的函数级控制**，不是完整JD器件或实际Jacobian同态对照；通量数也不是终端A/μm。它支持检查QF往返的消减风险，不证明当前1e-8平台由此产生。

已有`band_population`显式实验路径可以绕过SG电流中的QF相减，但GradQF驱动仍有真实物理含义；本轮未切换正在运行的默认quasi_fermi SG，也未把该控制宣称为已资格化的全程修复。

## 5. 迁移率、BTBT及可复用的稳定实现

Canali使用[ElementHighField.h](../../include/vela/physics/ElementHighField.h)中的long double中间运算。低场迁移率1000、vsat=1e7、beta=2的六场样本（0、1e-200、1、1e4、1e100、1e300，使用一致的cm系单位）值及两个报告偏导均有限；范围外、不同beta和复合梯度链未获全域证明。Lombardi/GradQF在非光滑零场的约定与弱人口1/p强切线仍需单独看，不能把大偏导直接判作计算错误。

局部E2 BTBT源有小场指数下溢提前返回；它本身应有电场依赖。当前BV关闭BTBT且拒绝该源，本轮未完整测试E2及NonLocalPath；不能据局部E2代码为非局部路径背书。

仓库[SiliconThermalPhysics.cpp](../../src/physics/SiliconThermalPhysics.cpp)的电热SRH/Auger另有稳定分支：正分裂用`np*(-expm1(-splitting))`，避免单独求巨大正指数。该结构可作为人口SRH修复的参考，但需要同时保持当前统计、单位、陷阱分母、值与导数契约，并做冻结态/平衡附近/浓度方向的验证，不能直接跨模型复制。

## 建议顺序与本轮限制

1. 优先修复已在实际坏态证明的SRH稳定值及人口切线；保留平衡附近精度。
2. 为共享AD/AV原语补极小正场的稳定边界处理，明确源值与导数均有限的契约。
3. 单独审查通用SRH/Auger的指数截断和数值饱和：决定真实导数还是有意近似，不混为Full。
4. 在实际冻结态中测量QF差分损失对SG值/Jv的影响，再决定是否推广直接人口算术。

本轮不修改任何物理模型或验收门，不把人工极端反例当作现有全程失败的实证，也不更改正在运行的窗口实验。
