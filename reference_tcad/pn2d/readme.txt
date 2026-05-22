本目录包含一个二维 PN 结二极管 Sentaurus TCAD 校准算例，共包含三个脚本：

1. pn2d_sde.cmd
用途：
- 使用 Sentaurus Structure Editor 生成二维 abrupt PN 结二极管结构和网格。
- 器件几何为二维矩形硅区，长度约 2.0 um，高度约 0.5 um。
- PN 结位置位于 x = 1.0 um。
- 左侧为 P 区，右侧为 N 区。
- P 区采用 BoronActiveConcentration 掺杂，浓度约 1e17 cm^-3。
- N 区采用 PhosphorusActiveConcentration 掺杂，浓度约 1e17 cm^-3。
- 左边界设置为 Anode 接触，右边界设置为 Cathode 接触。
- 全局网格较粗，PN 结附近设置局部加密窗口。
- 输出文件前缀为 pn2d，主要生成 pn2d_msh.tdr 等结构网格文件。

关键注意事项：
- 当前 Sentaurus 版本为 T-2022.03。
- 接触设置使用 sdegeo:set-contact，不使用已过时的 sdegeo:define-2d-contact。
- sde:build-mesh 使用空参数：
  (sde:build-mesh "snmesh" "" "pn2d")
- 不在 build-mesh 中直接加入 -quality，网格质量统计如有需要可后处理执行：
  snmesh -quality pn2d_msh.tdr

2. pn2d_sdevice.cmd
用途：
- 对 pn2d_sde.cmd 生成的 PN 结结构进行正向 IV 仿真。
- 输入网格文件为 pn2d_msh.tdr。
- 输出二维场分布文件为 pn2d_des.tdr。
- 输出电流电压曲线文件为 pn2d_iv.plt。
- 输出日志文件为 pn2d_des.log。

电极设置：
- Anode 初始电压为 0 V。
- Cathode 初始电压为 0 V。
- 正向扫描时，将 Anode 从 0 V 扫描到 1.0 V。
- Cathode 保持 0 V。

物理模型：
- Fermi 统计。
- EffectiveIntrinsicDensity(OldSlotboom)。
- Mobility 包含 DopingDep 和 HighFieldSaturation。
- Recombination 包含 SRH 和 Auger。
- 不启用 Avalanche，因此该脚本用于普通正向 IV，不用于 BV。

输出物理量：
- Doping
- DonorConcentration
- AcceptorConcentration
- eDensity
- hDensity
- Potential
- ElectricField
- eCurrentDensity
- hCurrentDensity
- TotalCurrentDensity
- eQuasiFermi
- hQuasiFermi
- SpaceCharge
- SRHRecombination
- AugerRecombination

关键注意事项：
- Plot 段中不要写 x y，因为 Sentaurus Device 不接受 x 和 y 作为 Plot 物理量。
- 坐标和网格信息已经包含在 tdr 文件中。
- Math 段采用兼容 T-2022.03 的保守写法：
  Math {
    Extrapolate
    Derivatives
    RelErrControl
    CNormPrint
    NewtonPlot
  }
- 不要在 Math 段中写 NotDamped、Method=Blocked、SubMethod=ParDiSo、Iterations=30。
- 迭代次数写在 Solve 段的 Coupled(Iterations=...) 中。
- Solve 先求 Poisson 平衡，再求 Poisson Electron Hole 耦合系统，然后进行 Quasistationary 偏压扫描。

3. pn2d_bv_sdevice.cmd
用途：
- 对同一个 PN 结结构进行反向偏压 BV 仿真。
- 输入网格文件为 pn2d_msh.tdr。
- 输出二维场分布文件为 pn2d_bv_des.tdr。
- 输出 BV 曲线文件为 pn2d_bv.plt。
- 输出日志文件为 pn2d_bv_des.log。

电极设置：
- Anode 初始电压为 0 V。
- Cathode 初始电压为 0 V。
- 反向扫描时，将 Cathode 从 0 V 扫描到 50 V。
- Anode 保持 0 V。

物理模型：
- Fermi 统计。
- EffectiveIntrinsicDensity(OldSlotboom)。
- Mobility 包含 DopingDep 和 HighFieldSaturation。
- Recombination 包含 SRH、Auger 和 Avalanche(OkutoCrowell)。
- Avalanche 模型用于观察反向高场下的击穿行为。

输出物理量：
- Doping
- DonorConcentration
- AcceptorConcentration
- eDensity
- hDensity
- Potential
- ElectricField
- eCurrentDensity
- hCurrentDensity
- TotalCurrentDensity
- eQuasiFermi
- hQuasiFermi
- SpaceCharge
- SRHRecombination
- AugerRecombination
- AvalancheGeneration

关键注意事项：
- Plot 段中同样不要写 x y。
- 如某些版本不支持 ImpactIonization 作为 Plot 变量，应删除，仅保留 AvalancheGeneration。
- Math 段采用和 IV 脚本相同的保守写法。
- 不要在 Math 段中放 Iterations=50；迭代次数应放在 Coupled(Iterations=...) 中。
- BV 扫描比 IV 更敏感，因此 MinStep 设置为 1e-8，Increment 设置为 1.2。
- 若仿真发散，可降低 MaxStep，或者将 Drain/Cathode 目标电压分段扫描。

整体执行顺序：
1. 先运行结构和网格生成：
   ../sentaurus/T-2022.03/bin/sde -e -l pn2d_sde.cmd

2. 再运行正向 IV：
   ../sentaurus/T-2022.03/bin/sdevice pn2d_sdevice.cmd

3. 再运行反向 BV：
   ../sentaurus/T-2022.03/bin/sdevice pn2d_bv_sdevice.cmd

4. 可视化结果：
   ../sentaurus/T-2022.03/bin/svisual pn2d_msh.tdr
   ../sentaurus/T-2022.03/bin/svisual pn2d_iv.plt
   ../sentaurus/T-2022.03/bin/svisual pn2d_bv.plt
   ../sentaurus/T-2022.03/bin/svisual pn2d_des.tdr
   ../sentaurus/T-2022.03/bin/svisual pn2d_bv_des.tdr

该算例的目标不是获得真实工艺级 PN 结结果，而是构造一个结构简单、边界清晰、掺杂明确、便于解析和导出数据的基准算例，用于自研二维器件仿真求解器的网格、掺杂、边界条件、泊松方程、漂移扩散方程、IV 和 BV 结果校准。