# SLOT IDL：原生参考与首批 Vela 初始化诊断（2026-10-05）

本批已建立 Sentaurus 完整 IDL 参考：漏压 0.1 V，栅压 0～5 V，77 个点，
退出码 0。**尚未完成同模型 Vela IDL 对比验证**，原因不只是数值收敛：
原生实际启用了默认 Bennett/Wilson 带隙窄化，而首批 Vela 诊断配置使用了
`none`。当前 Vela 接口只支持 `none/slotboom/old_slotboom`，需要先补齐该模型映射。

## 原始脚本与实际模型

原始文件来自 `D:/gpt-tmp/slot-LDMOS/V3E_NLDMOS_40V_Copy/`：
`IDL_des.cmd`、`sdevice.par`、`Silicon.par`、`Oxide.par`、`IDLIN_vis.tcl`、`gtree.dat`。
输入副本和 SHA-256 见 [manifest.json](../../build/s05idl/manifest.json)。

| 项目 | 本次原生 IDL |
| --- | --- |
| 温度/方程 | 300 K；Poisson、Electron、Hole；未启用热方程 |
| 电极 | source、drain、gate、substrate；没有 BVDS 的 SLOT 电极 |
| 扫描 | 零偏压初始化 → 漏压 0.1 V → 栅压 0～5 V |
| Silicon 迁移率 | IALMob(FullPhuMob, PhononCombination=1)，准费米梯度高场饱和 |
| Silicon 复合 | SRH(DopingDep)、Auger without Generation |
| 雪崩 | 未启用；日志中的 avalanche derivatives/Plot 名称不代表启用雪崩源 |
| 统计/带隙窄化 | 非 Fermi；实际默认 Bennett/Wilson BGN |
| 原生求解 | T-2022.03、ParDiSo、ExtendedPrecision；原脚本 Digits=5、RHSmin=1e-10 |

**纠正首轮解释：**原脚本把 `EffectiveIntrinsicDensity(OldSlotboom)` 注释掉，
不能据此推断“无 BGN”。[实际日志](../../build/s05idl/native_resolved/console.log)
明确记录 `Bandgap narrowing model: Bennett/Wilson with bandgap narrowing (no Fermi)`，
Silicon 覆盖段没有改写该项。参数文件 Bennett 的 Ebgn=0.00684 eV、
Nref=3.162e18 cm⁻³、dEg0=0；其浓度口径、适用范围和截断规则还需专项核对。
不得直接用 OldSlotboom 或 none 代替后称为原模型对比。
详见 [映射复核](../../build/s05idl/alignment_review.json)。

本次原生运行复用了此前本器件的工艺网格 `n1_fps.tdr`，SHA-256 为
`ec6cb57f5daa577094c70a5fc3909a2e69a4a6adb69d43c165edb3377fc5258b`。
没有重新运行 SProcess，也没有声称覆盖 Workbench 的全部几何参数分支。
Vela 继续使用对应既有导入网格：13,881 节点、27,275 三角形、7 个区域；
原生内部日志为 14,069 vertices、27,275 elements，与此前原生 BVDS 日志相同。
本批没有重新进行跨格式全网格映射审计。

## Sentaurus 参考结果

独立目录：
`/home/tcad/sentaurus_runs/vela_oracle_2022/slot_idl_20261005_resolved/`。
按原模板展开 300 K、等温方程及文件名；将材料 include 展开，按 gtree 的
`Van {1}` 展开参数宏（IDL 雪崩关闭，仍需满足参数文件语法）。
首轮遗漏该宏引发的语法错误保留在 `build/s05idl/native/`，未覆盖失败证据。

成功参考耗时 **542.59 s**，退出码 0，77 个 IDVG 点，漏压均为 0.1 V。
二维无额外 AreaFactor，按现有参考口径报告电流 A/µm。
下面是实际采样点，不把邻近值写成精确整数栅压：

| 实际栅压 / V | 原生 Id / A·µm⁻¹ |
| --- | --- |
| 0 | 1.22214119e-16 |
| 0.9850623691 | 1.24299072e-16 |
| 1.4850623691 | 3.80130356e-16 |
| 1.9850623691 | 4.09152212e-13 |
| 2.9850623691 | 1.70380361e-6 |
| 3.9850623691 | 1.30427501e-5 |
| 5 | 1.51088979e-5 |

`IDLIN_vis.tcl` 的常电流标准为 `1e-7*L`，当前 L=0.1，对应 1e-8 A/µm。
对本次离散点作**线性电流插值**得到 Vti≈2.60178 V；这不是额外精确偏压求解，
也不是已执行 SVisual 的 ExtractVtgm/全部提取流程。
完整数据：[原生 CSV](../../build/s05idl/native_idl.csv)、
[原始 PLT](../../build/s05idl/native_resolved/IDVG%28LIN%29_idl.plt)。

## Vela 首批诊断及资格边界

使用本工作树合并本地 main 后的 Release 程序，Git HEAD `0b31b408368f69ec8ab4544fc280823bffb8cad1`，
包含合并前保留的本地改动，不能仅凭 Git SHA 重建二进制；实际 runner SHA 及物理资源
哈希单独入档。Windows UCRT64，HDF5，默认直接后端 UMFPACK。
程序冻结为 `build/s05idl/runner.exe`，未改动历史冻结求解器。

Vela 诊断使用四电极、无外电阻、无雪崩、原点 0 V、冷启动；
没有导入 BVDS 高压态。保留元素 IALMob、现有外部 AverageBox 输运资源、
`material_local_circumcentric` 电荷/复合测度、`per_cell` Poisson 介电离散。
旧高压残差归一化未继承，由本点初始化确定。Newton 迭代上限 100，
每子进程墙钟上限 180 s；这些是数值工作量控制，不是验收放宽。

原门槛保持：abstol=1e-9、reltol=1e-8、有效载流子行比 eps_row=0.001、
scale_floor=1e-30、密度资格下限1e-200、全局闭合 tolerance=0.01。
此 IDL 没有外电阻，因此不套用 BVDS 负载线资格。

| 尝试 | 结果 |
| --- | --- |
| 冷启动平衡态 | 18.47 s；配置内完整求解通过，0 条超标行，最大行比1.4284e-12，全局闭合通过 |
| 漏压自适应预偏置 | 180 s 上限停止；0.1 mV 及缩步尝试触发行收敛/全局闭合线搜索拒绝，未接受非零点 |
| 同一平衡态直接求 0.1 V | 124.24 s；37 次接受 Newton 更新后拒绝，3,294 条电子行超标，最大行比1.00011446 |
| 内置 Gummel 密度恢复 | 4.66 s；明确报错：`per_cell Poisson permittivity is not supported by Gummel/DDAssembler`；未切换离散方式 |
| 密度坐标载流子块候选 → 原完整复验 | 候选降低电子块残差，但复验仍有420条超标行；没有把探针成功当作器件资格 |
| 第二轮密度块 → 原完整复验 | 复验27.23 s后拒绝，416条超标行，最大行比1.00065145 |
| 第三轮密度块 → 原完整复验 | 候选完成，完整复验达到180 s上限停止，未取得完整终态资格 |

直接 0.1 V 失败态的3,294条超标行全部为电子：3,123条 Silicon_1 内部，
171条 Silicon_1/Oxide_1 界面；不是少数多晶硅行刚越过门槛。
该次全局归一化残差约9.34e-12，但局部最大行比约1，体现两类指标并不等价。
失败行的 I–V CSV 电流占位值不作为真实器件电流使用。

本地共10个子进程，累计子进程用时约596.34 s（不是并行运行的总墙钟）。
全部已退出，没有遗留本批求解器。各轮有界密度块/完整复验的最终结果与成本见
[analysis.json](../../build/s05idl/analysis.json)；全部候选和被拒状态保留在新目录。
即使某个诊断点满足了配置内数值门槛，也不能消除 Bennett/Wilson 模型缺口。
本批尚无可用于同模型 IDL 比较的 Vela 非零漏压点，未启动 Vela 栅压扫描，
未触碰旧 BVDS 账本或恢复长扫描。

## 下一步

1. 优先核对 Bennett/Wilson 的实际适用材料、浓度定义、截断及能带/本征浓度处理，
   增加显式 Vela 模型和针对性测试；不能把本批 none 配置沿用为正式 IDL 基准。
2. 在模型对齐后复验冷启动与0.1 V预偏置，再判断弱电子行是否仍阻塞。
   如需密度恢复，必须支持现有 per_cell 离散，或沿已兼容的完整 Jacobian 块路径有界修复；
   不以换回旧离散或放宽行比来替代。
3. 预偏置完整合格后，才沿0～5 V栅压扫描，比较同偏压 Id、阈值、跨导、
   场/源项及重复性。原生完整参考已经就绪。

本批只增加独立运行与分析文件、此报告和文档索引，没有改求解器源码或验收规则。
输入/下载/程序哈希、配置门槛检查及全部尝试成本见
[证据清单](../../build/s05idl/evidence_manifest.json)和[分析脚本](../../build/s05idl/analyze.py)。
13项输入/程序/原生输出下载哈希与诊断配置合同检查通过；这不包含模型等效性，
`equivalent_physics=false`、`comparison_passed=false` 均显式保留。
没有将此前 main 合并的全量测试计作本次新验证。
