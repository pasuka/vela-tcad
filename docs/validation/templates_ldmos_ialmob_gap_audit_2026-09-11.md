# Templates/LDMOS IALMob 实现缺口审计

日期：2026-09-11。D5 双栅压原门限通过后开展的第一项 IALMob 工作。
本报告完成输入盘点、版本手册方程映射、实际界面晶向分类与分阶段验证设计；
尚未将 IALMob 接入求解器或完成该物理层验收。
自热、量子势与雪崩不纳入此阶段。

## 已核对的输入与代码

本地 `reference_staging/templates_ldmos_vg8_20260907/ialmob_preparation/Siliconc100.par`
包含默认及 `100/110/111` 四组 IALMob 参数。逐项解析确认默认组与 `100` 相同，
`110/111` 的 B、C、delta 等参数不同。不能从文件名推断网格每个界面均使用 `100`。
D4-classical 的 `IdVd.cmd:21` 显式启用 `Enormal(IALMob(AutoOrientation))`。
该 deck 的第 63–64 行还设置两种载流子的
`RefDens_*GradQuasiFermi_EparallelToInterface=1e12`。
进一步读取原始 bundle 的 `sprocess.cmd:66`，确认 `wafer.orient=100`；
同一 Siliconc100.par 的 `LatticeParameters:2889–2890` 给出
X=(1,0,0)、Y=(0,1,0)，并无生效的 CrystalAxis 反向定义。
因此采用参数文件的显式坐标映射，不凭文件名或工艺默认 flat 方向猜测。

已读取本机
`D:\工作\学习资料\TCAD软件手册\Sentaurus帮助文档\data\sdevice_ug.pdf`，
封面版本为 T-2022.03。核对页 86–87、407–416、426–427、886–889；
其中 408–410 的公式另经页面渲染核对，避免文本抽取丢失括号/指数。
手册与实际晶向信息已经找到，目前无需用户补充这些资料。

源码对照：

| 部分 | 本地输入/实际源码 | 缺口与处理 |
|---|---|---|
| 模型分派 | `MobilityModel.cpp:60` 起支持 Masetti、HFS、增强 Lombardi | 没有 IALMob 分派；须独立配置类型，不改变 D5 `constant_field` |
| 体迁移率 | IALMob 有 mumax/mumin、alpha/theta、n_ref、有效质量参数 | 页 409–410 明确使用完整 PhuMob 库仑部分；不能直接使用现有 Masetti |
| 掺杂输入 | 式 290、300–303 分别使用施主与受主浓度，并含团簇修正 | 当前 mobility 接口只有 netDoping；必须扩展为分别传入 Nd/Na，不能用净掺杂重建补偿区 |
| 声子/粗糙度 | 两者均出现 B/C/delta 等同名参数 | 必须按方程确认量纲、密度指数、场定义和温度因子，不能仅按名称映射 |
| 反型/积累层项 | D1/2_inv、nu*_inv、alpha*_inv 及手册中的 acc 参数 | 现有 `LombardiParameters` 没有对应项；本地 par 未列出的值需按版本默认表补齐并标明来源 |
| 距离衰减 | l_crit 与 l_crit_c | 现有实现仅有一个 criticalLength；不能合并后宣称等价 |
| 指数 | 四组 IALMob 参数的 nu 均为 0 | 现有 `lombardiLimit:376` 要求 nu>0，证明直接套用参数不成立 |
| 晶向 | 页 86–87：按最近界面方向选最近的晶面族，默认不平滑 | 当前模型没有晶向选择；先用显式 X/Y 将几何法向转到晶体系 |
| 高场驱动 | D4 保留 GradQF / EparallelToInterface 密度参考值 | 需要核对开启 Enormal 后的生效范围；当前 D5 的场算子审计不能替代这一项 |

## 已固定的实现合同

- 300 K，FullPhuMob 默认开启，PhononCombination 默认值 1；D4 未指定
  ClusteringEverywhere 或 AsPhPhuMob，不擅自打开这些选项。
- 按式 291–306 实现 2D/3D 库仑散射、反型/积累贡献、屏蔽函数及 G(P) 下限；
  按式 307–311 实现声子组合和粗糙度。l_crit 与 l_crit_c 分开处理。
- 手册页 415 明确 IALMob 自带掺杂依赖并自动关闭常数迁移率贡献；
  不叠加另一个 DopingDependence，以免重复计入体迁移率。
- 页 426–427 的 Enormal 默认取半导体/绝缘体界面的真实几何距离及其梯度法向；
  NormalFieldCorrection 默认 0，本 deck 未覆盖。不以最近网格顶点距离代替真实距离。
- 本 par 的高场 alpha=0、beta0=1.109/1.213、Vsat_Formula=1，保留实际值。
  页 415 中 alpha=1、beta0=2 的 Hänsch 示例是另一种选择，不能照抄来改变本算例。
- 原 D4 已输出 e/hENormal、e/hEparallel、e/hMobility，可用已有 TDR 先做局部对照；
  InterfaceOrientation/NearestInterfaceOrientation 若需直接验证，可另作诊断导出。

按已确认的晶体 X/Y 与手册最近晶面族规则，对实际网格的全部 Si/SiO2
界面边做独立分类：760 条中 618 条属于 {100}、142 条属于 {110}，
对应长度 8.383355 μm / 0.753101 μm；该二维平面下没有 {111} 边。
最小候选点积差 9.77937e-4，没有分类等距情况。这是边法向的几何推断，
尚未与 Sentaurus 节点级 NearestInterfaceOrientation 对照。
该结果说明应保留 100/110 两组参数，而不能把整个结构按 100 组处理。

本地审计结果：
`reference_staging/templates_ldmos_followup_20260911/ialmob_input_audit.json`。
其中保留四组参数、参数文件/deck/工艺脚本/手册/网格/源码 SHA256、
晶向矩阵、界面分类计数，以及无直接 Lombardi 字段的列表。
原始参数文件保持在外部证据目录。历史 Slot-LDMOS 的
`masetti_field_lombardi` 消融不能作为本算例真正 IALMob 的资格证据。

历史 Sentaurus D4→D5 消融的电流差异中位数/P95：Vg4 为 59.84%/60.74%，
Vg8 为 21.96%/46.47%，见
[阶段四消融](templates_ldmos_stage4_idvd_execution_2026-09-03.md)。
该差异说明此模型阶段有实质影响，不宜用经验缩放电流代替实现。

## 下一步实现与验证门槛

同日后续更新：已实现独立 300 K 计算核与 `ialmob_probe`，完成物理极限测试
及 D4 保存状态的初步局部对照。界面和体区长尾差异尚未闭合，尚未接入耦合
求解器；当前进展与剩余门槛以[后续实现报告](templates_ldmos_worker_and_ialmob_execution_2026-09-11.md)为准。

1. 方程与晶向证据已找到。接下来把 par 中显式参数和手册默认表合并为完整
   中性配置，逐字段注明来源与单位；先补齐 Nd/Na 接口和模型项，保持旧模型分派不变。
2. 独立实现体迁移率、表面声子/粗糙度、反型层项与距离因子。先做可检查的
   单点扫描：300 K，实际保存状态覆盖的掺杂/载流子/法向场/界面距离范围。
   测试正值、零场/远界面极限、连续性、单位等价及各晶向独立选择。
3. 在同一保存状态上逐项对照 Sentaurus 的 e/hMobility、Enormal 与高场驱动，
   区分模型差异和离散支持差异；若已有输出缺字段，再补导出相应局部量。
4. 先验 G3 原 31 点六门限，再验 Vg4/Vg8 前八点，最后两条 0–40 V 曲线。
   IALMob 对照采用已具备的 D4-classical 参考，D5 继续作关闭模型的回归控制。
   不以 D5 的参考电流判定开启 IALMob 后的准确性。

资料前置条件已解除。剩余工作是完整模型编码、Jacobian/极限测试、局部状态对照
和 D4 自洽曲线验收；本次未把资料审计或几何分类计作物理模型已实现。
