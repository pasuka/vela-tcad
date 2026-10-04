# SimpleMOS PhuMob 稳定链式导数修复与自洽对照

日期：2026-09-09。工作分支：`codex/simplemos-sdevice-validation`。

**已修复 plain PhuMob legacy 分支的弱载流子链式导数。** 原来 164 个独立非零交叉矩阵项中有 106 项失败、94 项被算成零；本次全部通过，最大相对差 2.71049e-14。四个工作点、两条初始化、修复前后共 16 次 DC 首次尝试均合格，Id 最大相对变化仅 1.66533e-15。这说明已测工作点的剩余电流差并非由这些弱导数的消减造成。

原生温度控制进一步约束了 G(P) 截断行为，但尚未识别其内部极小值搜索算法。反推的 G/P 数值仍只用于诊断。本轮没有启用 PhuMob 单元迁移率平均，没有恢复 Enormal 或高场饱和，没有改变物理常数、残差、电流定义或接受门槛。

## 1. 配置、网格和修复范围

- Vela 对照采用 n19/n23 × Vd=0.05/1 V，Vg=0.8 V，300 K。继承同一 SimpleMOS Si/SiO2 网格、掺杂和匹配 OldSlotboom 基础 ni 的独立材料副本；n23 网格为 1482 节点、2746 三角形。每个输入的文件身份均已冻结。
- 物理模型为 `phumob`、`edge_averaging=legacy`、OldSlotboom、掺杂相关 SRH、Boltzmann 统计。几何仍为已验证的 `transport_edge_geometry=element_box`、`poisson_charge_node_volume=signed_transport`、`poisson_permittivity_policy=cell_material`、`cell_box_policy=delaunay_transfer`。**输运几何使用 element_box 不等于迁移率也使用单元平均。**
- Windows UCRT64 Release 构建，线性后端 UMFPACK；HDF5/TDR 场导出实际运行成功。Vela 输入使用 `unit_scaling`；报告电流单位为 A/μm，电势为 V。PhuMob 标量输入在生产中为 m⁻³、迁移率为 m²/(V·s)，独立本构参考内部换算为 cm⁻³、cm²/(V·s)。
- 修复前保存了旧可执行文件、静态库和相关源码，之后才修改和构建；前后求解使用同一配置、网格、材料和两条种子状态。

生产修改涉及 [MobilityModel.h](../../include/vela/physics/MobilityModel.h)、[MobilityModel.cpp](../../src/physics/MobilityModel.cpp)、[AssemblerUtils.h](../../include/vela/equation/AssemblerUtils.h) 和 [CoupledDDAssembler.cpp](../../src/equation/CoupledDDAssembler.cpp)：

1. 新增 `dμ/dln(n)`、`dμ/dln(p)` 的显式解析链式求值，贯穿 Nsc、P、F、G、有效散射浓度、散射迁移率和 Matthiessen 合成。G 下限分支的导数为零，未截断分支使用解析 dG/dP；F 的求导先做代数消减。
2. 按现有 legacy 边端点载流子平均及相邻有效材料平均形成导数。将同边通量乘以 `(1/μ) dμ/dln(pop)`，再乘端点人口比例及 `±1/Vt`，装配到 ψ、φn、φp 的对应列，并向两端残差行以相反符号散射。
3. 保留固定迁移率 SG 直接导数，合并此前稳定 BGN 电势导数。该解析路径限定为 plain `phumob + legacy + Boltzmann`，且四个载流子指数严格位于 ±500 内；Fermi、表面/高场及指数截断范围仍走原路径。

这部分是手工推导并编码的解析链式导数，**不是自动符号工具生成的全部 Jacobian**。本轮不宣称其他模型、非光滑边界或任意邻接模板已得到完整解析验证。配置说明见 [config_schema.md](../config_schema.md)。

## 2. 独立导数验证

在上轮四个冻结状态的节点 792、1000、1009、1057 上，先关闭 SRH 隔离迁移率交叉项。参考使用独立 100 位 PhuMob 链式导数与实际基态边通量，按守恒符号装配；不以完整残差相减作为弱项的真值。

| 检查 | 修复前 | 修复后 |
|---|---:|---:|
| 164 个独立非零交叉矩阵项未过 1e-5 相对门槛 | 106 | 0 |
| 其中生产矩阵值为零 | 94 | 0 |
| 两个外部步长合计检查数 | 328 | 328 |
| 最大绝对导数差，原始装配残差/V | 3.69701e-16 | 1.80353e-31 |
| 修复后最大相对导数差 | — | 2.71049e-14 |

独立参考非零量级为 2.70534e-29～3.69701e-16。基态通量与残差闭合最大相对差 2.00705e-16，标量 μ 与独立参考最大相对差 5.73756e-16。完整残差双精度差分仍有 177 个零结果，这些没有被用于覆盖独立弱导数资格。

另在自由硅节点每隔三个节点取一组 ψ/φn/φp 方向，使用 1e-4、3e-5、1e-5、3e-6 V 四个幅度。相对范数撤掉探针原有的 1 下限，以 `||Jv−FD||/max(||Jv||,||FD||)` 计算；小三幅度的非弱交叉分块使用原 1e-4 门槛：

| 状态集 | 全部记录 | 正式分块检查 | 失败 | 最大相对差 |
|---|---:|---:|---:|---:|
| 原冻结四态，SRH 关闭 | 144 | 84 | 0 | 7.96607e-7 |
| 修复后四个合格自洽态，恢复 SRH | 144 | 84 | 0 | 7.96607e-7 |

后者不用于宣称弱 SRH 交叉源块已通过独立校准。明细见 [交叉项](../../reference_tcad/simplemos_sentaurus2022/phumob_fix_20260909/cross_analysis/cross_entries.csv)、[固定态分块](../../reference_tcad/simplemos_sentaurus2022/phumob_fix_20260909/blocks.csv)、[最终态分块](../../reference_tcad/simplemos_sentaurus2022/phumob_fix_20260909/post_jvp.csv)。复用的旧交叉分析器输出含 `production_changed=false`，仅代表分析器本身不改源码；本次生产修复的事实与两版身份由 [repair_summary.json](../../reference_tcad/simplemos_sentaurus2022/phumob_fix_20260909/repair_summary.json) 明确记录。

## 3. 同物理模型的自洽对照

修复前后各 8 次首次尝试全部合格，均没有使用合同预留的重载恢复。两端各 4 个双初始化比较全部通过，迭代次数逐工况、逐初始化均未改变。继承的接受条件为全载流子逐行 1e-6、KCL/Id 1e-8、端口重放 1e-8，双初始化电势 1e-6 V、密度相对差 1e-4、Id 相对差 1e-6。

16 个状态的最大全行比为 1.29584e-7，最大 KCL/Id 为 1.49939e-14，最大端口相对差为 1.06581e-14。全局 SRH 源的 32 个分量全部低于继承合同的 1e-10 下限，**全局相对源门槛本轮没有激活**，不能计作独立相对源闭合验证。

下表使用 Vela 初始化的修复后状态，Vg 均为 0.8 V。电流差定义为 `(Id_Vela/Id_Sentaurus−1)×100%`。原生参考使用已核验的原生单元迁移率形成方式，Vela 此处仍为 legacy 平均，因此这些数值是剩余差异的背景，不是单元 PhuMob 恢复的验收。

| NWell | Vd (V) | Id Vela (A/μm) | Id Sentaurus (A/μm) | 电流差 |
|---|---:|---:|---:|---:|
| n19 | 0.05 | 7.82648678718e-7 | 7.83567248147e-7 | −0.11722918% |
| n19 | 1 | 6.58960012622e-6 | 6.62365716117e-6 | −0.51417267% |
| n23 | 0.05 | 3.80835415340e-9 | 3.80999235651e-9 | −0.04299754% |
| n23 | 1 | 2.43630879426e-8 | 2.43716989100e-8 | −0.03533183% |

修复前后 Id 最大相对变化为 1.66533e-15，迭代数也不变。结论仅针对上述四个工作点：稳定弱导数修复提高了 Jacobian 数值正确性，但没有改善剩余 Id 差异。它不是整个 0–1 V 区间的排因结论，也不能与此前 Masetti 简化模型的误差混用。

明细见 [输入合同](../../reference_tcad/simplemos_sentaurus2022/phumob_fix_20260909/contract.json)、[修复前尝试](../../reference_tcad/simplemos_sentaurus2022/phumob_fix_20260909/before/attempts.csv)、[修复后尝试](../../reference_tcad/simplemos_sentaurus2022/phumob_fix_20260909/after/attempts.csv)、[双初始化及电流比较](../../reference_tcad/simplemos_sentaurus2022/phumob_fix_20260909/dc_comparison.csv)、[源门槛激活账本](../../reference_tcad/simplemos_sentaurus2022/phumob_fix_20260909/source_gate_activation.csv)。

## 4. 原生 G 截断：温度控制与开源实现

沿 n19、Vd=0.05 V、Vg=1 V 的原生保存态，新增 299、300、301、350 K 控制。首次仅修改 `Physics Temperature` 后 Load：四次均正常退出，但导出的温度场全部仍为 300 K，因此三个非 300 K 点未通过目标温度预检。原始批次与 [失败账本](../../reference_tcad/simplemos_sentaurus2022/phumob_floor_20260909/temperature_preflight_failure.csv) 保留。

按本地 T-2022.03 用户手册第 126–127 页，另建 `ramped` 批次，使用 `Goal { Model=DeviceTemperature Parameter="Temperature" Value=... }` 从 300 K 做 Quasistationary 渐变。四点的实际硅温度字段均与目标值一致，偏置误差为零，KCL/Id 不超过 1.106e-15。全部任务已经结束并取回；没有遗留后台仿真。

对每个温度、每种载流子，从一个预先选定的全截断单元反推 G 下限，其他单元均作为留出对照。算法比较固定采用同一原生场、掺杂和 box 顶点权重：

- Vela：对 dG/dP 的符号在对数 P 区间二分 80 次，得到精确极小值附近的截断位置。
- Charon：本地 `Charon_Mobility_PhilipsThomas_impl.hpp` 第 599–624 行从 P=0.3246 做 Newton 搜索，更新幅度不超过 1e-5 停止，上限 500 次；本轮四温度的结果与精确根基本一致。
- Genius：本地 `Si_mob_Philips.cc` 第 76–132、166–171 行，以 10 K 温度表和区间宽度 1e-3 的三分搜索产生 P 下限，运行时按温度取表项；该算法也不能复现全部原生空穴单元。

两份开源文件的绝对路径与 SHA 保存在 [source_identity.json](../../reference_tcad/simplemos_sentaurus2022/phumob_floor_20260909/ramped/analysis/source_identity.json)。手册内容取自本地 [sdevice_ug_local_2022.txt](../../build-release/m79_research/sdevice_ug_local_2022.txt)，本轮没有据开源实现推断 Sentaurus 必然使用同一算法。

| T (K) | 电子精确 Gmin | 电子原生反推 G 下限 | 空穴精确 Gmin | 空穴原生反推 G 下限 |
|---|---:|---:|---:|---:|
| 299 | 0.097832686224 | 0.097833586637 | 0.079690921281 | 0.079692786146 |
| 300 | 0.098108696549 | 0.098109650093 | 0.079943109997 | 0.079944873802 |
| 301 | 0.098384133889 | 0.098385141432 | 0.080194802337 | 0.080196468360 |
| 350 | 0.111229007721 | 0.111232831051 | 0.091962604580 | 0.091962695250 |

精确 G 截断对原生空穴单元的留出最大相对差依次为 7.80311e-6、7.35772e-6、6.92865e-6、3.26869e-7，均未过原 1e-7 单元门槛。仅在精确 Pmin 左侧替换 G 下限，299 K 仍有 2.39595e-6、350 K 仍有 1.05811e-7 差异。这证明需要同时考虑截断位置，而不能只替换 G 高度。

进一步对已反推的 G 下限求 `G(P)=Gfloor` 的左右两个根，不另拟合自由参数，比较在左根截断 P、在右根截断 P、单侧对 G 值设下限这三种形式。原 1e-7 单元门槛不变；另外以 1e-12 作为更严格的算法辨别尺度。

| 空穴温度 | 左根 P 截断的留出最大差 | 右根 P 截断的留出最大差 | 可辨别结论 |
|---|---:|---:|---|
| 299 K | 7.77e-16 | 2.39595e-6 | 支持左根；P≈0.2873417141 |
| 300 K | 5.55e-16 | 2.99117e-6 | 支持左根；P≈0.2878542259 |
| 301 K | 7.77e-16 | 7.77e-16 | 网格未采样两根之间，无法区分 |
| 350 K | 1.05811e-7 | 5.55e-16 | 排除左根；右根 P≈0.3126784916 与单侧 G 下限形式均匹配 |

电子四温度均缺少两根之间的样本，三种形式都可匹配到约 6.66e-16，不能独立识别电子截断位置。在当前 G 分支上，右根 P 截断与测试的单侧 G 下限形式等价，因此 350 K 结果也不能识别内部操作顺序。上述仍是基于反推 G 下限的诊断重放，**没有识别原生搜索容差或温度选点算法，更没有把反推数值写入生产模型**。

数据见 [原生工作点](../../reference_tcad/simplemos_sentaurus2022/phumob_floor_20260909/ramped/native_points.csv)、[温度字段核验](../../reference_tcad/simplemos_sentaurus2022/phumob_floor_20260909/ramped/analysis/temperatures.csv)、[算法比较](../../reference_tcad/simplemos_sentaurus2022/phumob_floor_20260909/ramped/analysis/summary.csv)、[截断根](../../reference_tcad/simplemos_sentaurus2022/phumob_floor_20260909/ramped/cutoff/roots.csv)、[留出结果](../../reference_tcad/simplemos_sentaurus2022/phumob_floor_20260909/ramped/cutoff/results.csv)。

## 5. 构建、回归和证据保存

Release 全部目标构建成功。新增三个 Catch2 数值测试覆盖极弱解析导数、单位/温度/As-P 物种控制、实际装配的弱交叉列及守恒符号，18 个 PhuMob 测试全部通过。测试代码见 [test_phumob.cpp](../../tests/test_phumob.cpp)。

完整 `ctest --preset windows-ucrt64-release --parallel 4` 为 **769/785 通过，253.10 秒**。16 项失败名称与前置 BGN 阶段完全一致：13 项历史身份或替代链校验、3 项缺失历史 `tests/test_mos_mixed_material.cpp`。未出现新增失败；完整测试集仍不是全绿，没有修改旧身份基线或删除失败。见 [回归比较](../../reference_tcad/simplemos_sentaurus2022/phumob_fix_20260909/ctest_comparison.json)、[失败分类](../../reference_tcad/simplemos_sentaurus2022/phumob_fix_20260909/ctest_failure_categories.csv)。

本轮脚本和数值结果位于 `phumob_fix_20260909`、`phumob_floor_20260909` 两个证据目录；所有初始温度预检失败、修复前状态与旧二进制均保留。当前源码、可执行文件、报告和状态入口另存于浅层忽略目录 `build-release/phumob_fix_snapshot_20260909/`，避免后续修改覆盖当前可复核版本。历史证据没有被改写成通过。本轮工作保留在当前工作区，未另行提交或推送。

## 6. 下一阶段边界

稳定 legacy 链式导数修复及其四点自洽隔离对照已完成。后续需要实现受保护的“顶点 PhuMob 本构 → box 单元平均 → 有符号边系数”候选，并补齐相邻单元第三顶点的 n/p 列，以及装配器、端口和探针的一致性；不能只移除 element_box 模型保护。

原生 G 搜索算法尚未确定，精确截断的空穴本构仍保留失败资格。应把该差异和单元平均候选分别记录，在定义明确的固定状态/小扰动控制中量化各自影响，不能用反推下限冒充已核实的生产参数。完成这些检查后再进入八点 PhuMob 双初始化对照，随后才恢复 Enormal → 高场饱和。本轮没有重算完整 0–1 V 曲线。
