# BVmethods 原始消费路径与逐边驱动采集（2026-10-15）

T470p 和 Sentaurus VM 均已重新连接。被动追踪确认了实际元素边电流入口，取得
加载首态和更新态的 binary80 Bernoulli 输入/输出及驱动端点项。六次 T470p 审计的
完整首轮打印系统与原始基准逐值一致；跨求解器完整 RHS 和两弱行仍未收口。
接续[原生统计审计](bvmethods_reference_statistics_2026-10-15.md)，
数据见[配套 JSON](bvmethods_passive_consumer_2026-10-15.json)。

后续[端点身份与 RHS1 绑定记录](bvmethods_endpoint_rhs1_binding_2026-10-15.md)已补齐
本报告待办的字段身份、缓存消费和两弱行全支撑绑定；跨求解器全 RHS 仍待办。

## 范围与状态绑定

原 BVmethods 网格 2600 节点、5020 Tri3，Si region 3 为 2152 节点/4156 单元；
300 K、Rext=1e7 ohm*um。沿用源关闭固定态、Fermi/OldSlotboom、DensityIntegral(240)、
原迁移率与接触。Sentaurus T-2022.03-SP2、ExtendedPrecision/binary80、单线程，
既有 blocked decomposition 配置；没有运行新的 Vela Newton 轨迹或更换线性后端。
长度 um、人口/DOS cm^-3、ψ/QF V、带边 eV；本轮底层驱动表按原生内部量记录，
没有给未确认的端点项强加物理单位。

保留旧零源回调，未增加 ReadScalar 统计采集接口。GDB 仅读取调用栈、寄存器和内存，
不调用被调试程序函数、不写求解数据；软件断点正常修改指令。原生退出 2 是单次 Newton
更新预算结束，不是稳态收敛。v1/v4/v5/v6/v7/v8 的变量、Jacobian、RHS1、更新向量和
RHS2 均与原始首轮基准一致。这是原生采集身份核验，不是 Vela/native 残差一致。

实际存在两次 RHS 装配调用。第一次未观察到活动电流条目重算，不能将后一次的
操作数冒充 RHS1 输入。新采集分别标记加载后的 Plot 首态和更新后的第二次 RHS 状态；
首态样本与整个 RHS1 系统共同封存，但缓存消费和目标弱行的逐项绑定仍需补齐。

## 实际消费与字段定义的边界

`e/hFermiDeltaEntry::Compute` 各命中两次：一次在 Plot 的 VertexMobility 依赖链，
一次在 RHS 装配依赖链。因而“原始路径完全不计算 FermiDelta”的解释不成立；
依赖链命中本身不能证明该字段直接作为 SG 扩散因子使用。

只读运算入口检查确认它采用 OneOverTwoIntegrals 的 Inverse/Direct/dDirect 路径；
这与 MultiValley/DensityIntegral 人口路径的实际操作数尚未一一绑定。上轮约 6.37% 的
电子字段差不能直接解释为物理错误。`QuasiFermiPerT` 的原始计算入口没有命中；
其派生运算还涉及 inverse/dInverse 和温度链式项，不能只凭名称解释成 QF/VT。
这些字段的完整物理定义及实际输入身份仍未取得资格。

最初追踪的 `e/hBernoulliEntry` 没有命中。实际路径是
`e/hCurrentWithoutMobilityComby_ElementEdge → CurrentOnEdgeComputations::ComputeOnEdge
→ DES_Bernoulli`。后续采集沿该路径读取，不再依据未执行入口推断本算例的 SG。

## Bernoulli 与端点组合

最后一批包含电子/空穴、首态/更新态各四组，共 16 组有边索引的样本。
T470p 用独立 100 位 Decimal 检查四个输出：B(x)、B′(x)、B(−x)、B′(−x)，
并检查正反方向恒等式。最大绝对差为 1.484552655118454191615161883707662743060046493700774794049293020935017447895840981E-19。
该检查确认参数及返回值读取约定；样本不是完整边集合，也没有覆盖两条目标弱行。

沿实际构造路径取得多个端点项，以精确有理数模拟 binary80 的逐步舍入。
三个 ties-to-even 数学自检通过。完整项重建中，16/16
个驱动逐位一致；先合并端点再相减的对照为 0/16。
后者最大绝对差为 1.376937e-17（原生内部驱动量）。
不得把这种有限样本的运算顺序对照推广为整个 SG 或完整残差通过。

两项解释仅闭合电子样本，空穴调用路径还包含两项额外输入；未通过的两项解释保留。
各端点项暂用 U/W/Z/T 表示，它们的字段身份仍待核对，不能称为已经确认的 QF
reference/increment，也不能据此对 Vela 添加经验偏移或修改材料默认。

对本批有向边 i→j，已验证的逐步 binary80 组合为：电子先求 Uᵢ−Uⱼ、Wᵢ−Wⱼ
再相加；空穴在该和上再加 Zᵢ−Zⱼ、再减 Tᵢ−Tⱼ。每个减法和加法单独舍入。
这里的 U/W/Z/T 只是实际输入槽的中性名称，不是新物理参数或拟合修正。

## 证据与后续

本轮八批原生归档共 220 文件逐项哈希核验通过，包含启动失败证据。
上轮 T470p 补充归档的 8 文件已补取并验证。T470p 六批完整系统身份审计、最后一批
16 组 Bernoulli/导数检查和三项精确舍入自检已完成，分析文件和实际源码哈希已下载核验。
未为私有诊断重建生产求解器，没有声明全 CTest 通过。

v2 在原生可执行文件加载前设置地址断点而启动失败；修正加载时序后继续。
v3 栈追踪成功但寄存器转换阻止操作数读取；后续已修复。私有分析脚本的目录重试和
输出编码失败也有记录，不计入数值资格。任务目录、主机与作业记录只保存在忽略目录
`build-release/bvmethods_reference_increment_20261014/`；公开远端路径使用
`<T470P_DIAGNOSTIC_ROOT>`、`<SENTAURUS_DIAGNOSTIC_ROOT>`。

下一步仍在 Q2/Q3：核对这些实际端点项的字段身份和缓存消费，定向采电子 2156／空穴
2515 的全部支撑边；把保精度驱动、人口、Bernoulli 和逐边电流绑定到原始 RHS1，
再区分表示/舍入与物理本构差异。完整 RHS 对齐后才进入四臂 Jacobian 和 Newton 更新。
生产默认、门限、模板及发行包未改；原始从零完整纯 QS 和 Q4–Q7 仍未通过。
