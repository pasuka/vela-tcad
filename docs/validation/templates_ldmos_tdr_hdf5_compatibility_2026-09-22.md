# TDR 与 Vela HDF5 状态格式兼容性核对

日期：2026-09-22。范围：当前代码及两个实际文件的结构核验，不是新增 TDR 写出或物理状态转换资格。

## 结论

本算例的原生 Sentaurus TDR 确实是 HDF5 文件，可用同一 HDF5 库或 h5py 打开。
HDF5＋HighFive 候选采用独立的 `vela.ddstate/1` 结构；共用容器不等于共用数据结构。
现有代码支持 TDR 到中间网格/场文件的单向导入，不包含 Vela 状态到可由 Sentaurus
读取的 TDR 写出器。更改文件扩展名不能建立兼容性。

## 实际样本

| 内容 | 原生 TDR 样本 | Vela HDF5 样本 |
|---|---|---|
| 来源 | `reference_staging/templates_ldmos_d0_transport_20260912/raw/vg8_low_des.tdr` | T470p 本批 `r0_hdf5_vg4/fixed/child_00252_direct/state.h5` |
| HDF5 检测 | true | true |
| 主结构 | `/collection/geometry_0` | 根级属性 `schema=vela.ddstate/1` |
| 网格 | 10241 个复合记录 `(x,y)`，单位 µm；10 个区域/接触/界面 | 只记录 `node_count=10241`，没有坐标、拓扑、区域或接触 |
| 场组织 | `state_0/dataset_N/values`，按区域、节点/单元/接触支撑组织 | 每个场一个长度10241的未压缩 float64 数组，行号对应外部网格节点 |
| 密度 | `eDensity/hDensity`，硅区域5723项，单位 cm⁻³ | `electrons_m3/holes_m3`，单位 m⁻³ |
| 电势 | `ElectrostaticPotential`、电子/空穴 QF，单位 V | `psi/phin/phip`，另存 QF 参考值与增量，单位 V |
| 温度 | `LatticeTemperature`，单位 K | 当前三方程状态 schema 无温度字段 |
| 其他元数据 | 区域、材料、location type、quantity、单位及换算属性等 | schema、统一单位字符串、节点数、field_names |

两个样本属于不同物理状态，仅比较文件结构，不比较物理场数值。
原生样本 SHA256：`a2879dbc9b91208a6014cf6508c95710bb574b37e706a96f3c8a6fddcd233073`。
Vela 样本 SHA256：`949acb38d7e6e45f621c5f1a7499a54848321b820bdeb2612259e8f3dbf44102`。

## 当前导入器的支持范围

[SentaurusTdrReader.cpp](../../src/io/SentaurusTdrReader.cpp) 使用 `H5Fopen`，
固定读取 `/collection/geometry_0` 及其 `state_0`；坐标只取 x/y，元素仅处理
点、边和三角形。不能由本次验证推断任意三维、多几何、多状态或所有 TDR 版本都支持。

节点场支持全局节点顺序或区域内升序全局节点编号；单元场与接触标量另行映射。
导出中间 CSV 和 field manifest 时保留单位标签。坐标由显式 `coordinateUnit`
选择 µm/cm 后输出 µm；读取单位字符串不代表已经自动执行所有物理单位换算。
当前代码没有对 TDR 的 `conversion factor`、`unit:exponent` 等任意组合做通用转换。
场长度不匹配时部分路径会导出重叠前缀并记录 partial/warning，使用转换结果时
必须检查 manifest，不能仅凭 HDF5 打开成功认定状态完整。

[DDStateFormats.cpp](../../src/io/DDStateFormats.cpp) 和
[dd_state_public.py](../../scripts/dd_state_public.py) 读取 Vela 状态时检查 schema、
单位、节点数、字段、类型及形状，再沿用规范状态解码验证。当前文件内没有网格散列，
节点数相等本身不足以证明节点编号相同；本批依靠外层冻结输入及散列审计保证对应关系。
独立使用恢复文件时仍需保留其配套网格、物性、偏压与来源信息。

## 实测核验

- 使用当前 worktree 的 UCRT64 Release 构建 `sentaurus_import` 和
  `test_sentaurus_tdr_reader`；`ctest -R '^sentaurus_tdr$'` 选中1项并通过，
  Catch2 实际15项测试、174条断言通过。
- C++ 导入器成功读取上述原生样本并生成 inventory JSON。
- h5py 成功读取两个实际文件；Vela Python 读取器成功恢复10241行 Vela状态。
- 将原生 TDR 复制为 `.h5` 后，Vela 读取器因缺少 `schema` 拒绝；
  将 Vela `.h5` 交给 C++ TDR 导入器，因缺少 `/collection/geometry_0` 拒绝。
- 未调用 Sentaurus 对 Vela 文件进行原生加载，故没有认领此方向兼容性。

诊断脚本、样本、inventory 和 JSON 结果保存在
`build/t470p_hdf5_full_repeats_20260922/`，其中
`audit_tdr_schema.py` 可复现本次结构检查，`tdr_schema_audit.json` 保存结果。
本机 h5py 显示编译时 HDF5 2.1.1 与运行时2.2.0版本警告；本次检查通过只证明
这些路径，不扩大为任意版本组合兼容性。诊断在本机进行，没有增加远端仿真计算负载。

## 若后续需要交换文件

1. TDR → Vela：依据已验证的节点/区域映射合并物理场，明确密度
   `n[m^-3]=10^6*n[cm^-3]`；核实电势/QF符号、参考系、绝缘体占位值与缺失字段策略。
   QF参考值/增量不是原生字段的简单改名，必须按Vela表示规则构造并检查恢复残差。
2. Vela → TDR：需补齐网格、区域、接触、场支撑、原生名称和单位及所需元数据，
   再用 Sentaurus 原生工具做加载、显示和逐点数值回读验收。HDF5文件本身合规不够。
3. 当前可继续把 HDF5 用作公开恢复状态容器；HighFive 提供 C++ HDF5接口，
   不负责 TDR 语义转换。本轮不为交换格式扩大生产默认或温度状态资格。
