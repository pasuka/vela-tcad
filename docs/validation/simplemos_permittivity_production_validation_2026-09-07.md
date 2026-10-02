# SimpleMOS 逐材料介电装配生产接入与八点验证

日期：2026-09-07。结论：显式逐单元材料介电装配已接入正式 Poisson、Gummel、耦合 Newton 与耦合接触反应电荷路径；48/48 次 DC、16/16 组双初始化、48/48 组端口一致性检查通过原门槛。新实现复现此前已校准的全部原生 K 结果。

**单独介电修改仍不满足物理放行门槛。** 采用 Vela 自身迁移率时，低 NWell 的四点绝对误差均增大，高 NWell 的 Vg=1 V 两点也变差。加上此前独立校准的原生有效迁移率输入，才能达到高 NWell 0.00539656%–0.00675240% 的绝对电流误差。因此新策略保持显式选择，现有默认和物理接受标准不变；本轮没有将原生迁移率替换写入正式算法。

## 1. 实现范围

新增 `mesh_geometry.poisson_permittivity_policy`：默认 `legacy_average` 保留旧行为；`cell_material` 使用 Gij=Σcell εcell(d/l)cell。逐单元几何先与自身材料 ε 相乘，再累加到公共边。该系数由一个共享函数供给三个装配器，耦合残差与解析 Jacobian 使用同一缓存，接触反应电荷读取同一残差路径。

`cell_material` 无附加数据时使用 BoxGeometryBuilder 实际生成的单元贡献，包含原有负 cotangent 的回退/截断选择。可独立提供 `poisson_cell_edge_coefficients`，每个 Tri3 单元一条记录，明确 cell_id、原序 node_ids 和局部边 (0,1)/(1,2)/(2,0) 的三个无量纲 d/l。输入不包含 ε，材料数据库仍决定 ε。

输入检查拒绝单元数量/节点顺序错误、非有限或负系数、退化单元及冲突策略；不将错误编号或不支持的几何静默平均。n19/n23 分别核对 2742/2746 单元、8226/8238 个系数。已导出的原生系数足够，本轮没有上传或启动新 sdevice 仿真。

原生 box 修改仍由明确输入提供，尚未成为 Vela 自动生成的通用规则。输运几何、三项 Poisson 电荷体积、SRH 体积、常数、SG 公式和接受门槛保持各自原定义。本次未改变 PN2D BV 模板的原子策略或全局默认。

源码：[共享系数](../../include/vela/equation/AssemblerUtils.h)、[几何构建](../../src/mesh/BoxGeometryBuilder.cpp)、[配置解析](../../src/simulation/ConfigParsing.cpp)。接口见[配置说明](../config_schema.md)。

## 2. 隔离对照与电流误差

下表为 (Id/IdSentaurus−1)×100%。前三列均使用 Vela 正式迁移率；最后一列的原生有效迁移率保持为独立诊断输入，其残差、Jacobian 迁移率缓存和端口三处使用同一因子。它用于检查介电修复迁移的一致性，不代表原生迁移率实现已生产化。

| NWell | Vg/V | Vd/V | 原联合方案 | 局部单元几何 | 原生单元几何 | 原生单元几何＋原生有效迁移率 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| n19 | 0.8 | 0.05 | -0.013368% | -0.106682% | -0.110077% | +0.003488% |
| n19 | 0.8 | 1 | -0.287643% | -0.439866% | -0.443506% | +0.002990% |
| n23 | 0.8 | 0.05 | +0.191855% | -0.054132% | -0.058068% | +0.005880% |
| n23 | 0.8 | 1 | +0.351202% | -0.050248% | -0.047084% | +0.006752% |
| n19 | 1 | 0.05 | -0.762563% | -0.772822% | -0.784129% | +0.001006% |
| n19 | 1 | 1 | -1.013748% | -1.074595% | -1.088844% | +0.001242% |
| n23 | 1 | 0.05 | +0.025813% | -0.112740% | -0.115333% | +0.005397% |
| n23 | 1 | 1 | +0.073109% | -0.169243% | -0.171076% | +0.005457% |

完整 DC 分组：旧策略回放 8 次，局部单元几何 8 次，原生单元几何＋Vela 迁移率双初始化 16 次，原生单元几何＋原生有效迁移率双初始化 16 次。完整记录见[电流表](../../reference_tcad/simplemos_sentaurus2022/permittivity_production_20260907/comparison.csv)与[高低配对门槛](../../reference_tcad/simplemos_sentaurus2022/permittivity_production_20260907/physical_gates.csv)。

仅介电改变时，高低 NWell 误差配对与绝对误差没有全部通过，说明原有误差存在抵消；不能依据介电离散公式正确就将单项修改作为整机精度放行结论。原生 K 与原生有效迁移率的组合在四组配对上均通过高 NWell 改善、低 NWell 不恶化和配对 log 误差改善条件。

## 3. 数值资格与迁移精度

- 固定状态 16/16 组通过。新旧原生 K 装配的 Poisson 差异，相对实际介电扰动源最大 1.3453e-13，门槛 1e-8；旧策略回放的残差保持逐位一致。全部载流子残差和输运边通量保持逐位一致。
- 224 个非弱块的小步长 Jv 检查全部通过，最大相对差 5.55854e-09，门槛 1e-4。96 个弱 SRH 交叉块检查保持同态解析值逐位一致，沿用原先独立弱块资格；此次不把弱块差分噪声或单方向通过扩大为任意 Jacobian 完整性结论。
- 48/48 DC 通过全部 1814 个自由载流子行及原全局闭合检查。最差逐行比值 9.10177294e-07，位于 `m65_n19_vd_1p000000_endpoint_vg020/native_cells_native_seed`，仍低于 1e-6，仅约 8.98% 裕量；最大 KCL/Id 为 2.00486e-14，门槛 1e-8。
- 16/16 组双初始化通过。最大势差 1.683e-08 V、密度相对差 6.51174e-07、Id 相对差 1.44345e-11；原门槛依次为 1e-6 V、1e-4、1e-6。
- 新旧原生 K 方案的最大 Id 相对漂移 5.66214e-15；全自由 Si 节点 ψ/φn/φp 最大迁移差 3.4e-16 V。节点 1000/1009 的迁移差单列于[局部状态表](../../reference_tcad/simplemos_sentaurus2022/permittivity_production_20260907/local_state_migration.csv)。这验证生产接入复现旧合格状态，不代表那些状态剩余的原生场误差已经消失。

## 4. 测试与环境

Release 完整构建成功。新增 6 个 Catch2 用例、45 个断言全部通过：不等界面贡献、材料/节点顺序、均匀材料与长度单位、负 cotangent 分支、法向位移和串联电容解析解、三条装配路径与反应电荷、残差/Jacobian 差分、几何重建与错误输入。见[测试源码](../../tests/test_permittivity_assembly.cpp)。

全套 CTest **760/776 通过**，仍为原有同一组 16 项失败：13 项历史源码哈希/替代链检查、3 项引用缺失 `tests/test_mos_mixed_material.cpp`。新增失败为 0；旧哈希及失败记录未改写。详见[逐项审计](../../reference_tcad/simplemos_sentaurus2022/permittivity_production_20260907/ctest_audit.json)。因此不声称全套 CTest 已通过。

Windows MSYS2 UCRT64、C++20 Release；实际 DC 使用 Eigen SparseLU/COLAMD。构建确认 HDF5/TDR、UMFPACK、SPQR 可用，DC 未切换后端。n19/n23：Boron 1e17/2e17 cm⁻³，1480/1482 节点，942 个 Si 节点、907 个自由 Si 节点；网格 μm，状态密度 m⁻³，宽度 1 μm，电流 A/μm。

300 K、Boltzmann/no-BGN、matched ni=1.0750038488844236e10 cm⁻³，Masetti 总杂质与掺杂相关 SRH；HFS、表面迁移率、Auger、雪崩、DG 关闭。max_iter=200、reltol=1e-7、abstol=1e-12、stall floor=1e-9，ψ/准费米更新上限 0.35/0.025 V，contact_basin、原标量线搜索、线性迭代修正 4 次。全局连续性 tolerance=1e-6、source_floor=1e-10；低于净源下限的结果不声明任意相对精度。

## 5. 后续边界与证据

本轮完成逐材料介电装配的显式生产接口与八点迁移验证。下一步应解释并实现原生有效迁移率的单元/边形成方式，验证原生 box 修改的自动生成范围；随后沿已有合格状态对 q/ε0/热电压约定和节点 792/1057 的 SRH 源体积分开展独立同扰动校准。既有毫伏级极少数载流子场差仍需保留检查。整个组合与生产回归接续完成前，维持 16 工况和完整 0–1 V 曲线的原放行门槛。

验证脚本：[执行入口](../../scripts/validate_simplemos_permittivity_production_20260907.py)、[报告入口](../../scripts/report_simplemos_permittivity_production_20260907.py)。阶段顺序为 build、prepare、preflight、run、analyze；prepare 后合同与代码哈希冻结，既有输出拒绝覆盖，重跑需新目录与新合同。生产变更前源码和二进制前像已独立保存，不重写旧证据。完整摘要与关联见[最终证据](../../reference_tcad/simplemos_sentaurus2022/permittivity_production_20260907/final_evidence.json)。本轮未提交或合并代码。
