# BVMethods物理算术修复及T470p复验，2026-11-02

## 当前结果与范围

已实现所选SRH、共享指数/AV边界及通用复合Full导数修复。T470p独立Release/MKL PARDISO单线程复验：八组C++、六冻结态与原AND门0→0.01 V冷窗口全部通过，worker退出0，无自有任务或源目录计算子进程残留。原两个坏态的Full/bootstrap非有限条目由6/12降为0。

本地四组定向C++回归及全量Release构建通过。全量CTest共1217项：初次6项因Windows临时路径过长导致HDF5创建失败，改用短临时目录后原6项全部通过；合并全量及复测为1211项通过、6项环境条件跳过。未运行完整BV，不宣称全程Newton收敛或效率问题已解决；P3仍待用户决定。旧VM的30组窗口继续使用旧冻结源码和二进制，未修改或重启。

前因见[物理模型坐标审计](bvmethods_physics_coordinate_audit_2026-11-02.md)。本轮源文件、测试和可重复诊断入口在工作树中；任务回执、源码清单、原始证据与分析在ignored的`build-release/bvmethods_physics_arithmetic_fix_20261102/`，不把主机连接信息写入仓库。

## 实现及数值契约

1. 人口SRH使用`populationRecombinationExcess`直接构造`np-gamma_n gamma_p ni²`及人口方向，采用补偿宽乘积差保留近乎平衡的差值。人口源项不再经过fn/fp相减及巨大正指数；固定n/p时的ψ方向由计算结构得到0。JD gamma值与切线从稳定log-gamma原语得到。QF模式的正分裂分支避免单独计算exp(+split)。
2. 共享AD的结构零不再执行inf×0；非零指数方向在必要时合并指数后再窄化，保留可表示导数。不可表示的函数值仍不可表示，不以清零掩盖物理/状态失败。AV普通有限路径保留，非有限中间方向用合并指数求值，极小正场不再0×inf。
3. 通用广义复合不再独立截断log(gamma)和分裂到±500，而是保留对数因子并合并指数，近乎平衡仍用expm1。相应QF方向也按同一表达式求导，先合并增长项方向以免大项相消。已有源项饱和值保留，Full偏导与实际饱和值一致，在拐点采用非活跃侧；经典Auger的乘积、Newton链式导数和Gummel实际源率同步。没有给状态新增浓度上限或修改物理参数。
4. SG只做实际冻结态测量。继续使用默认quasi_fermi；band_population仍为显式实验算术。当前证据不足以推广默认路径。

主要源码：[人口复合原语](../../include/vela/physics/experimental/PopulationRecombination.h)、[BVMethods运行时](../../src/equation/BvMethodsCoupledRuntime.cpp)、[AD指数](../../include/vela/equation/Tri3LocalForwardAD.h)、[AV系数](../../include/vela/physics/experimental/BvMethodsTri3Avalanche.h)、[共享统计](../../src/physics/CarrierStatistics.cpp)、[复合值与偏导](../../src/physics/RecombinationModel.cpp)、[通用装配](../../src/equation/CoupledDDAssembler.cpp)。

新增明确排除默认构建的[vela_physics_frozen_audit](../../scripts/diagnostics/physics_frozen_audit.cpp)，只读输入packed状态，不执行Newton更新、不改变接受门；调用为`vela_physics_frozen_audit request.json`，request包含config、packed_binary、inner、out。

## 验证配置与版本身份

2600节点/5020 Tri3，300 K、Joyce–Dixon/OldSlotboom、总杂质迁移率、SRH与GradQF AV，BTBT关闭，R=1e7 Ω·μm。potential_density/quasi_fermi SG、局部非正保留1%、下限1e-100 cm⁻³、精确消元、full_step、40轮，原AND、负载线enforce、F≤1e-8、有效RMS<1、负载线≤1e-6 V及原KCL/其余门的report/enforce属性保留。当前profile仍拒绝Auger/BTBT；共享Auger修复另由单元测试覆盖，没有为本轮BV启用它们。

本地采用MSYS2 UCRT64 Release；T470p采用独立Release、显式MKL PARDISO运行库、每求解进程单线程。T470p最终数值快照清单SHA256为`c4b607b20c207ce777c29efd127e851cbf61f663d500335b67f58c5d8fe631c2`。本地后续仅为全量测试补齐了既有ExactValueBits optional编码支撑及其测试；实际物理求值源码逐文件与T470p清单一致，核对结果保存于`numeric_source_identity.json`。不把旧提交号直接当作修复版身份。

T470p终态证据186文件，压缩包SHA256=`2b9db7c69e65d8cb190fa46ba3d54b8f1d84f7e0be20fa9b3d03c35c82b7dbfd`；已下载并逐文件核验。编译器/运行库与二进制身份保留在私有任务记录中，不以共享主机墙钟推导加速比。

## 原失败态源值与矩阵复验

以下是离散行中的SRH贡献，不是体积复合率。参考用同态物理系数及`np-gamma_n gamma_p ni²`回算；分子与既有110位算术检查一致，不宣称整个材料模型获得110位精度。

| 节点 | 旧有限但错误的行源值 | 修复值 | 浓度恒等式参考 | Full/bootstrap非有限条目 |
|---|---:|---:|---:|---:|
| 689 | 2.79125227e+15 | 56.6213172111 | 56.6213172111 | 6 → 0 / 6 → 0 |
| 331 | 2.64688441e+13 | 1.54346205755e-13 | 1.54346205755e-13 | 12 → 0 / 12 → 0 |

两个坏态现在F、Full及bootstrap均有限，但仍是远离解的坏状态，不能把有限矩阵当作该状态已收敛。另四个高/低压帧同样有限。

## 实际冻结态SG与Jv测量

六帧使用同一packed向量，原始状态SHA逐个核验。比较原quasi_fermi和显式band_population；后者同时改变输运和AV中SG电流的算术，GradQF场继续保留。仅为比较Full，band分支清除与其不兼容的bootstrap专用掩码，不在该分支比较近似矩阵。

固定列尺度为ψ槽1、活跃人口槽n/p，非半导体零人口槽保持正尺度而方向为0。固定方向为小ψ扰动与保正的正人口相对方向；Jv采用宽累加。单侧FD有h=1/.1/.01，行尺度保持算子既有输出，规定行按实际prescribed_rows另行排除。这不是原始Newton方向，也不是所有弱列逐一闭合。

| 帧 | 电子自由行max abs(delta F) | 空穴自由行max abs(delta F) | 全向量max abs(delta Jv) |
|---|---:|---:|---:|
| AND_special_0 | 3.31471e-10 | 2.91038e-11 | 1.25233e-13 |
| AND_special_1 | 6.06207e-09 | 4.18368e-10 | 6.71907e-13 |
| high509 | 9.1812e-11 | 5.52579e-14 | 1.42109e-14 |
| high512 | 6.44399e-11 | 5.62814e-14 | 1.06581e-14 |
| low_near | 2.52363e-11 | 1.01143e-13 | 3.55271e-15 |
| low_stage | 2.62662e-11 | 7.98372e-14 | 7.54952e-15 |

所有帧的两条SG算术路径SRH差为0，Poisson自由行差为0；差主要在输运，AV也有更小差异。低压近解帧电子最差自由行4756：F从−3.92156e-12变为−2.91578e-11；空穴行6435从−1.01481e-13变为−3.38429e-16。不能用全局范数小就断定每个弱行相对差小。

固定方向Jv差很小，FD误差随步长缩小，但未完成多步轨迹和所有弱列验证。因此本轮决定暂不推广直接人口SG默认；不通过切换离散算术掩盖当前剩余非线性问题。

四个普通高/低压帧另与归档旧RHS比较，packed哈希一致，F∞不变、最大分量差约5.3e-13。这是旧VM与新本地的跨主机比较，不能把小差异唯一归因于某个修复或声称逐位轨迹不变。

## 原门冷窗口

T470p原AND冷Poisson 0→0.01 V窗口：接受2点（含零偏），目标点3次更新；内压0.00524425278207732 V，电流4.75577110989e-10 A/μm。设备残差2.25930428e-09、有效RMS=0.803042399、负载线误差2.38919655e-08 V，求解器与worker退出0。逐行原强制门通过。

报告型载流子行未满足1点、全局连续性未满足1点；它们保持report属性，不把它们写成强制门闭合。没有击穿电压资格或完整BV结果。

## 回归与保留记录

- 本地四组定向C++全部通过：recombination、carrier_statistics、bvmethods_av_contract、bvmethods_population_newton。新增性质回归覆盖指数范围、结构零、近平衡补偿差、两个真实强分裂人口输入、饱和Full偏导、经典Auger链式导数与Gummel源值一致性。
- T470p八组全部通过，上述四组加coupled_load_line、mkl_pardiso、joyce_dixon_contract、state_archive；显式PARDISO无静默回退。
- 首轮局部链接因接口补齐时对象未同步失败，最终全量重建后四组通过；未把混合对象当作验证结果。
- 冻结诊断首轮把氧化层人口0用作列尺度，违反正尺度接口；已修正并保存准备失败。T470p R1八组C++通过，随后停在该准备检查，未跑冷窗口；R2只改诊断文件，数值核心逐字不变，并完成全部验证。
- 本地扩大到全量构建时，发现既有ExactValueBits不支持optional字段，阻塞两个旧测试目标。已补齐presence/浮点payload编码与absence、正负零、NaN payload回归。这是测试诊断支撑，不改变求解器。全量构建已通过。全量CTest初次1205通过、6失败、6跳过；6失败均为DCSweep检查点HDF5创建失败（长临时路径），保持同一二进制、门限和测试断言，仅将临时目录缩短到37字符后，6项复测全部通过，耗时7.21秒。合并1211通过、6跳过。跳过项为本地2项PARDISO运行时条件测试及4项未配置归档样本的可选Vela sdevice入口测试；T470p显式PARDISO测试已独立通过。初次与复测JUnit均保留，不把复测写成全量首轮全过。

## 尚未取得的资格

修复说明的是算术和值/导数契约，并不证明旧bootstrap首步失真、低压反复保护、全部困难暖态或完整BV收敛已解决。真实SRH修复与旧VM窗口策略消融保持不同版本；旧P0c若推进需基于旧冻结源隔离，不混入本轮变更。后续完整BV仍须用户依据窗口结果决定，不自动启动P3。
