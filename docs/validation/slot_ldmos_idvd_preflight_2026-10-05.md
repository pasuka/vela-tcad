# SLOT IdVd：原生基线与 Vela 组合能力预检

日期：2026-10-05。前置任务：[IDL 独立材料冷路径复验](slot_ldmos_idl_cold_material_validation_2026-10-05.md)。

## 结论与资格边界

IDL 已完成从新冷平衡到 Vd=0.1 V、Vg=5 V 的全路径复验，89 个合格态，
11 个原生相同门压采样的最大电流差约 0.01111%。该结论只覆盖其已验证的物理组合。

本批按用户授权连接 Sentaurus 虚拟机，在新目录执行**主 `IdVd_des.cmd`**：
冷平衡、Vd=0 下建立 Vg=5 V、固定 Vg 扫描 Vd=0–60 V。原生正常完成，
退出码0，漏压段93条原生接受记录；输出已经下载并核对哈希。

**Vela 同物理 IdVd 对比尚未完成，本批 Vela 合格 IdVd 点为0。** 四次本地有界预检
在配置或装配初始化阶段被拒。阻塞是当前模型组合和材料映射的能力边界，
并非载流子行比资格失败。本批未改源码、物理默认值、IDL 配置或任何验收门槛，
未以关闭 Fermi/雪崩或改用 IdVd1 替代主脚本。

## 原生输入、路径及运行

来源：`D:/gpt-tmp/slot-LDMOS/V3E_NLDMOS_40V_Copy/IdVd_des.cmd`、
`gtree.dat`、`sdevice.par`、`Silicon.par`、`Oxide.par`。
Workbench 实际参数为 Vgate={5}、Vdrain={60}；不能把参数定义中的显示默认值1/20
当成本次运行值。300 K，`model=-`，未启用 Thermodynamic。

虚拟机：SSH别名 `sentaurus`，严格主机密钥核验；目录：
`/home/tcad/sentaurus_runs/vela_oracle_2022/slot_idvd_20261005`。
使用已归档的同算例 `n1_fps.tdr`，SHA-256：
`ec6cb57f5daa577094c70a5fc3909a2e69a4a6adb69d43c165edb3377fc5258b`。
场量导入为13,881个几何节点；原生日志内部顶点计数14,069，不混作不同网格的证据。
冷平衡场量 CSV 保留材料分区和节点映射。

解析 Workbench 占位符、展开参数文件，并增加冷平衡/门压终点和能带场输出；
没有改原物理及扫描步长。程序 T-2022.03-SP2，ParDiSo、ExtendedPrecision、5线程。
原生 Digits=5、RHSmin=1e-10、电子/空穴 ErrRef=1e10、Iterations=25。
原生数值停止判据不等同于 Vela 的完整资格门。

| 路径 | 原脚本步长（归一化） | 本次电压目标 |
|---|---|---|
| 门压建立 | Initial=.001、Increment=1.5、Min=1e-8、Max=.05 | 0→5 V，Vd=0 |
| 漏压扫描 | Initial=1e-4、Increment=1.2、Min=1e-12、Max=.02 | 0→60 V，Vg=5 V |

原生从北京时间13:21:26至13:35:25，包装墙钟839 s，求解器报告835.97 s；
20分钟硬上限未触发。原生内存峰值643 MB。没有遗留的本批求解进程。

## 实际物理模型

以下载的 `native_complete/console.log` 运行回显为准，而非根据脚本缺省项猜测：

| 项目 | 主 IdVd 本次实际启用 | 已验证 IDL |
|---|---|---|
| 统计 | Fermi | Boltzmann |
| 带隙 | Bennett/Wilson，Fermi分支 | Bennett/Wilson，no Fermi |
| Si迁移率 | Philips unified + Lombardi + GradQF高场饱和 | IALMob + GradQF |
| Si源项 | 掺杂依赖SRH、Auger无Generation、van Overstraeten/de Man雪崩 | SRH、Auger；无雪崩 |
| PolySilicon | 常数迁移率，无Si专属源项 | 常数迁移率，无Si专属源项 |
| 偏压 | 固定Vg=5 V，扫描Vd | 固定Vd=.1 V，扫描Vg |

冷平衡场量确认 PolySilicon 电子/空穴迁移率为1417/470.5 cm²/(V·s)，
浮点误差范围内为常数。四个电极为 drain/gate/source/substrate；本脚本没有外电阻。

## 原生曲线归档

电流单位 A/µm。下列全部为 PLT 实际采样，未把近邻点伪称为整数电压精确点。
完整93条记录见 [native_drain_path.csv](../../build/s05idvd/native_drain_path.csv)。

| 实际Vd / V | 原生Id / A·µm⁻¹ |
|---:|---:|
| 0.09068806583 | 1.442064502e-5 |
| 0.9211667051 | 1.283021499e-4 |
| 4.989820418 | 2.879840335e-4 |
| 10.11243379 | 2.909493980e-4 |
| 20.03810264 | 3.200336049e-4 |
| 40.43810264 | 4.772186019e-4 |
| 60.00000000 | 6.496050568e-4 |

这是开态输出曲线，不是关态 BVDS。末点体电流 −3.063106765e-6 A/µm，
不能把漏、源电流直接互为相反数而忽略体端电流。未给出 Vela 相对误差或 IdVd 通过结论。

## Vela 有界预检

程序沿用 IDL 冻结 Release/HDF5/UMFPACK runner，SHA-256：
`18f5544db5d4daaf41552ff50681b25ae96ef30c3a7d9798f63bd29061cc0714`。
每次上限60 s，均在4 s内退出，没有写出合格状态。预检配置是能力探测配置，
其 IDL 材料底稿**不是已验证的 Fermi 材料映射**。

| 预检 | 退出原因 | 可得结论 |
|---|---|---|
| external_full | external AverageBox雪崩仅允许SLOT IALMob专用组合 | 不能直接将既有原生几何路径改成PhuMob+雪崩 |
| mesh_default_full | Intrinsic reference avalanche requires Si-only element IALMob source | 去掉外部几何仍不能处理该非零材料参考势+雪崩组合 |
| element_fermi_noav_isolation | Element box transport requires cell_material geometry | 此次先触及几何先决条件，不能拿它证明Fermi失败 |
| element_fermi_noav_cellmaterial | element_box_phumob Jacobian要求smooth Boltzmann density interval | 完成几何先决条件后仍在初始化被拒；此运行错误不单独证明所有Fermi场景均失败 |

另有直接源码证据：

- [ConfigParsing.cpp](../../src/simulation/ConfigParsing.cpp) 的
  `validateCarrierTransportCouplePhysics` 限制外部AverageBox雪崩组合。
- [CoupledDDAssembler.cpp](../../src/equation/CoupledDDAssembler.cpp) 构造函数
  明确拒绝 `element_box_phumob` 与 Fermi 或雪崩组合；还限制非零材料参考势的雪崩。
- [MobilityModel.cpp](../../src/physics/MobilityModel.cpp) 当前 PhuMob 路径未提供已验证的
  “Si使用PhuMob、Poly使用常数”材料分流接口；IALMob的专用常数材料覆盖不能自动推广。
- [sentaurus_import.py](../../scripts/sentaurus_import.py) 对 BennettWilson+Fermi 导入主动拒绝，
  防止把既有Boltzmann映射当作已验证Fermi映射。

这些是组合限制，不表示 Vela 完全没有 Fermi、PhuMob、IdVd扫描或雪崩功能。
本批没有删除限制条件、把多晶硅参考势置零或更换为简化物理来取得曲线。

## Fermi带隙的补充只读检查

冷平衡原生 `BandgapNarrowing` 场不等于 IDL 使用的纯 Bennett 窄化表达式。
例如硅节点11404：原生0.8459645578 eV，纯 Bennett 为0.3016735689 eV。
在10个按掺杂分位选取的Si/Poly节点，使用独立数值积分计算既有
`fermiStatisticsBandgapCorrection` 所对应公式并与原生场量比较：9点差值绝对值
小于8e-9 eV；高掺杂补偿节点11404仍差6.9216e-4 eV。

证据：[fermi_bgn_spotcheck.json](../../build/s05idvd/fermi_bgn_spotcheck.json)。
这只是公式/场量诊断，尚未定位极端掺杂点差异来源，也未验证密度、接触、广义SG及Jacobian
的一致性，不能据此解除导入器保护或宣称Fermi材料包已经通过。

## 验证与后续边界

- 下载压缩包及18个原生文件SHA-256均通过；原始5份脚本/参数文件哈希未变。
- 冷平衡TDR导入成功；早期单独下载的场量与最终归档TDR哈希一致。
- 原生漏压记录门压均为5 V、漏压单调、末点精确到60 V，解析检查通过。
- 运行现有 `test_solver_physics_maps_explicit_bennett_wilson`，1项通过，
  包含Bennett/Fermi拒绝组合子检查。没有重编译求解器或以测试通过替代数值资格。
- 所有预检、哈希、源文件绑定和退出记录见
  [analysis.json](../../build/s05idvd/analysis.json)、
  [preflight.json](../../build/s05idvd/preflight.json)、
  [evidence_manifest.json](../../build/s05idvd/evidence_manifest.json)。

下一阶段应先完成主IdVd物理组合的实现/验证：Fermi/Bennett材料与接触一致性，
Si PhuMob/Lombardi与Poly常数迁移率分流，以及相应Fermi输运和雪崩组合。
每项先用小算例和本次冻结场量做回归，再做SLOT新冷平衡、5 V门压建立和同漏压采样。
这涉及求解器物理组合扩展，超出本批只用现有能力执行对比验证的范围；不能把解除配置
保护视为完成实现。原资格门槛保留，正式定量对比仍待该阶段完成。
