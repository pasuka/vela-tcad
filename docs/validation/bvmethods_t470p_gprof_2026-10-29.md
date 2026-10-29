# T470p当前核心BVMethods热点与优化（2026-10-29）

当前修改先整理为本地提交854058ae，包含人口对数外推、权威人口归档、直接人口SG实验接口、JD系数及原始量验证。本轮随后在T470p对标准JSON入口开展独立Release／gprof配对，实施两项不改变算术公式的重复计算优化。已有vela VM完整扫描保持原冻结源码，未被替换、停止或重启。

## 显式未知量切换

核心NewtonConfig::stateVariables和JSON solver.state_variables均已存在，默认potential_quasi_fermi；显式potential_density使用电势及电子／空穴浓度作为权威未知量。配置在newtonConfigFromJson及构造器中校验，并进入CoupledDDAssembler、核心人口Newton及外接电阻延续。以下为现有合格BVMethods配置中的关键字段，不是独立完整算例：

```json
{
  "solver": {
    "coupled_operator": "bvmethods_tri3_qf",
    "state_variables": "potential_density",
    "population_sg_current": "quasi_fermi"
  }
}
```

算子名含qf不表示仍求解QF；state_variables独立选择未知量。population_sg_current另行选择SG算术表达式，quasi_fermi表示在SG中使用由权威人口派生的QF。人口模式当前限定300 K、明确BVMethods/JD物理，不代表所有器件和模型均已支持。独立能带原始量的0.221%对照仍为隔离原型，未因本次提交或性能优化变成生产默认。

## 测量条件

T470p为Windows UCRT64 GCC16.2.0、GNU gprof2.47.20260726，8逻辑处理器可用。独立Release与gprof Ninja目录，O3/NDEBUG、LTO关闭；gprof通过VELA_ENABLE_GPROF在编译和链接启用-pg，并沿用MinGW关闭profiling可执行文件ASLR的既有设置。实际使用MKL PARDISO、sequential threading layer和单线性／BLAS线程。

物理输入沿用冻结BVMethods网格：2600节点、2152 Si节点、5020 Tri3，300 K、JD／OldSlotboom及FD-BGN、总杂质迁移率、SRH、自洽GradQF AV、BTBT关闭，R=1e7 Ω·μm。全局材料／物理常数与生产路径不变。核心potential_density、quasi_fermi SG、自动导数恢复、局部非正保留1%、默认精确消元、对数人口外推和原联合验收门。标准延续的预测内压控制开启，无法推导时回退0.02 V。

窗口从自身冷Poisson至外压0.100562182394906 V；共7个接受点（含零偏），更新次数为0、3、3、3、4、4、4。第一项0不包括独立Poisson初始化更新。本轮限制的是扫描终点，未放宽任何残差、更新量、KCL或守恒门；不构成完整BV的加速资格，也不是与Sentaurus同点表的性能比较。

## gprof定位与实现

基线gprof版窗口68.045 s，首次普通Release52.867／52.486 s。gprof采样self合计21.17 s，其中_mcount_private和__fentry__占7.86 s。不能把这些采样占比当作普通Release墙钟占比，也不能用插桩程序的加速倍数代替正式性能结果。GNU说明要求-pg同时参与编译和链接，flat profile的self与调用图children不同，参见[编译说明](https://sourceware.org/binutils/docs/gprof/Compiling.html)和[flat profile说明](https://sourceware.org/binutils/docs/gprof/Flat-Profile.html)。MKL和部分系统库未插桩，其内部不能完整展开。

热点分为两条直接调用链：

1. **同一节点的统计值在相邻三角形反复计算。**BvMethodsTri3Operator的transport／evaluate此前对每单元顶点调用JD人口导数，transport另算统计逆函数。现在在每次装配内部，按节点计算相同的密度单位转换、导数与所需逆函数，再按原局部AD槽散布。没有代换公式、减少导数或跨状态缓存；调用结束即释放临时表。
2. **表面几何构建反复判断相同区域与接口的适用性。**基线updateSurfaceMobilityCellGeometry被调用87次，内部同一区域的每个单元、每条候选接口都重复进行模型和区域字符串判断。现在在一次几何构建内按区域预筛接口，保留原接口次序、最近距离及并列距离选择规则。实时电场继续由当前电势计算，没有跨状态缓存法向场。

| gprof可观测调用 | 基线 | 最终优化 |
|---|---:|---:|
| JD evaluate | 26,757,086 | 12,314,686 |
| JD inverseDensity | 30,622,164 | 7,060,420 |
| surfaceMobilityAppliesToRegionPair | 262,038,980 | 313,194 |
| 主要std::string与字面量比较实例 | 1,679,044,746 | 108,378,786 |
| updateSurfaceMobilityCellGeometry | 87 | 87 |

这些是gprof记录到的实例调用数，不是所有被内联数学操作的完整计数。几何构建次数没有减少，减少的是每次构建内的重复判断。缓存仅限本次装配／几何构建，不增加跨状态或跨owner的失效契约。

单次普通Release阶段计时中，包含式dc.record_point由24.072降至6.610 s；post.contact_current为3.953／3.701 s，linear.total为1.466／1.332 s。父子阶段不相加，后两项的小幅变化不表示本轮修改了端口公式或线性算法。该窗口主要收益来自重复物理／几何准备，PARDISO并非主要优化目标。

## 正式普通Release配对

两版均关闭gprof插桩，使用相同物理、状态、输出、诊断和线程设置。顺序为基线→优化→优化重复→基线，编译和状态比对均在计时区间之外：

| 版本 | 两次墙钟，s | 平均，s |
|---|---|---:|
| 原始基线 | 50.786991、49.635077 | 50.211034 |
| 最终优化 | 32.061930、31.675270 | 31.868600 |

本窗口平均加速约1.576倍、墙钟减少36.53%。两次样本支持该窗口的实际收益，不外推到完整BV、其他网格或其他主机。优化gprof版34.464 s仅用于热点复核，不进入正式加速比。

首次增量尝试曾因补丁保留了早于gprof对象文件的时间戳而未重编译，该组不作为优化剖析证据。最终任务在提取补丁后按远端时间刷新文件时间戳，并核验Release及gprof运行时对象文件SHA256均已改变；两目录分别重新构建，正式结果只取此轮。

## 数值验证与范围

T470p原始Release／重复／gprof及最终基线前后／优化两次／优化gprof共8组，I–V和每次更新CSV一致。用独立C++状态解码器还原双精度值并按IEEE754字节校验，7组对照共63份HDF5、504个场数组与参考逐位一致，包含符号零；输出路径和源配置哈希等归档元数据不要求相同。两项优化没有减少该窗口的Newton更新数。

本机另保留不同延续控制的15点数值窗口，基线与最终优化I–V、更新记录和17份状态／136个场数组一致；它只用于数值回归，不用于性能结论，因为测量期间本机另有构建／CTest负载。零偏及两个有限态检查点还逐字节比较完整7801维矩阵、RHS、6765维精确消元矩阵及其RHS、Newton方向，共15份CSV；最终核验记录见同目录final_numerical_comparison.json。

新增数值用例检查有／无导数请求的输运及AV源值一致、人口变化后本次装配重新计算；新增几何用例检查重复区域单元、区域／接口选择、缺失接口和实时场。BVMethods30例934断言、PhuMob／几何19例219断言通过。完整Release共1169项：1165通过、0失败、4项外部Sentaurus参考e2e测试跳过；Schottky大写名称8项串行执行，其余1161项并行执行。JUnit和各门结果见verification_summary.json及对应日志，未修改容差或参考基线。

PN2D当前模板束、默认掺杂口径、体积策略与物理选择保持原值；共享几何优化保持原候选与浮点距离计算次序，适用回归包含当前已注册PN2D配置／体积／M2门。没有以本轮性能结果重新宣称此前未闭合的原生矩阵或完整BV资格。

## 证据

ignored根build-release/bvmethods_t470p_gprof_20261029保存：源码提交与patch SHA清单、私有部署回执、基线和最终worker、双构建日志、profiles.json、T470p状态位模式摘要、原始和最终状态、算子文件哈希及完整测试结果。远端原始gmon.out、flat／call graph和阶段计时保留在独立任务目录。任务路径、连接资料和大型生成数据不提交。
