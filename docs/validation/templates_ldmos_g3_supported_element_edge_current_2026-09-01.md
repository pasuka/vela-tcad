# Templates/LDMOS G3 Sentaurus 支持的 element-edge current 接口审计

## 结论

T-2022.03-SP2 对本问题可用的公开接口已经审计完毕。`Math {
ElementEdgeCurrent }` 在完全相同的高端点 VSV 单模扰动状态上被日志确认激活，
但 minus/plus 两支的 iteration 0、iteration 1 和终态全部归一化导出字段均与
基线逐字节相同，端口电流也在日志精度下相同。它既没有改变本算例可观测的
electron continuity assembly，也没有改变终态顶点 `eCurrentDensity`。

官方 Tcl/PMI 数据接口能够读取 element-edge box coefficient 和
element-vertex Measure，却不能读取 edge 位置的 `eCurrentDensity`。实际能力
探针在完成 Newton iteration 后给出原生错误：

> Tried to read undefined Edge-Vector eCurrentDensity !

因此，约 `33.7806012x` 的 Sentaurus/Vela electron-row 绝对标量无法再通过该
版本的受支持公开输出归属到单边电流或逐行缩放；若必须解释其绝对值，只能
请求 Synopsys 对 T-2022.03-SP2 的内部 assembly/row-scale 约定作版本专属披露。
这不会改变已经确认的局部空间模式，但要求继续禁止把 NewtonPlot
`eContinuityRhs` 的绝对幅值作为跨引擎电流验收门。ledger 保持 `draft`。

## 官方接口合同

T-2022.03-SP2 *Sentaurus Device User Guide* 给出的相关合同如下：

- 第 1254--1256 页：PMI Device Data 的 `ReadCoefficient()` 返回
  element-edge box coefficient，`ReadMeasure()` 返回 element-vertex Measure；
  `ReadScalar()`/`ReadVector()` 只能读取已注册在相应 location 的数据；
  `ReadFlux()` 是顶点场梯度在 box 边界上的面积分再除以 box volume，不是
  element-edge current 导出。
- 第 1483--1493 页：Tcl runtime 暴露相同的 Mesh/Data、Coefficient、Measure、
  Scalar/Vector 和 Flux 访问合同。
- Appendix F 的数据注册表只把 `eCurrentDensity` 注册为 vertex vector，没有
  edge-vector 数据集。
- Math 命令表把 `ElementEdgeCurrent` 定义为电流密度的另一种 element-edge
  approximation；它是装配选项，不是 element-edge current 输出开关。
- 第 208--210 页明确指出 NewtonPlot 的 RHS、error 和 update 是内部、实现相关
  的诊断量；普通 drift-diffusion 没有公开的 general row-scale 字段。

上述合同给出两个相互独立的判据：先用受支持的 Math 选项检查该近似是否改变
当前状态的可观测装配，再用 Tcl/PMI 数据位置请求检查单边量是否可导出。

## 严格同状态 A/B

两组 deck 均从同一 high-endpoint VSV 状态出发，对 node 3721 施加 `+/-1 uV`
electron-QF 对称扰动。唯一的 Sentaurus 输入差异是候选加入
`Math { ElementEdgeCurrent }`；网格、Physics、其余 Math/Solve、初值和输出
字段全部冻结。候选日志均出现 `With ElementEdge Current Density
approximation`，排除了选项未生效。

| 导出阶段 | minus 公共字段 | plus 公共字段 | 结果 |
| --- | ---: | ---: | --- |
| iteration 0 NewtonPlot | 20/20 | 20/20 | 逐字节相同 |
| iteration 1 NewtonPlot | 20/20 | 20/20 | 逐字节相同 |
| 最终状态 | 36/36 | 36/36 | 逐字节相同 |

逐字节相同的关键量包括 NewtonPlot `eContinuityRhs`、第一步 electron density
update、`ElectrostaticPotential`、`eDensity`，以及终态顶点
`eCurrentDensity`。两支 drain electron/hole/conduction current 也分别与基线
在日志精度下完全相同。由于 minus/plus RHS 文件未改变，node-3721 中央差分
标量在 A/B 前后均为 `33.78060115741988x`。

这项零结果的范围很窄但具有决定性：对这个冻结状态，`ElementEdgeCurrent`
不是 `33.78x` 的来源，也不能充当所需的 assembly-level 单边真值出口。它不
证明该选项对所有器件、网格或偏压恒为无效。

## Tcl/PMI 能力探针

能力探针先调用 `ReadCoefficient()` 作为 element-edge location 的正控制，再
调用 `ReadVector $::des_data_edge "eCurrentDensity"`。日志顺序证明 Newton
iteration 与 NewtonPlot 写出已经完成，随后接口才报告 undefined Edge-Vector；
因此失败来自数据注册能力，而不是 deck 未求解、Tcl 未执行或 edge handle
无效。

`ReadFlux()` 没有被伪装成替代答案：它的手册语义是从 vertex field 重构 box
通量散度，不能提供残差装配中每条 element-edge current，也不能辨识内部
row-scale。继续用顶点 `eCurrentDensity` 做几何重构只会回到此前已确认的
后处理闭合地板。

## 冻结解释边界与停止规则

- 受支持公开接口搜索已穷尽：coefficient/Measure 可读，edge current 与普通
  DD row scale 不可读。
- 排除 `ElementEdgeCurrent` 作为本 VSV row 标量的原因或修正项。
- 保留 NewtonPlot 的空间形状、方向、同引擎相对偏压变化和受控扰动响应；禁用
  其跨引擎绝对 RHS 门。
- 不修改 Vela C++ production default，不启用 IALMob/predictor，不调整阈值或
  迁移率参数，也不因该内部标量重跑 31 点曲线。
- 最大 gm 资格仍未通过；本审计只关闭可观测性分支，不把差异批准为固有地板。
  ledger 在器件门槛获得独立处置前继续保持 `draft`。
- 若业务要求解释 `33.78x` 的绝对来源，下一动作是向 Synopsys 请求该版本
  element-edge assembly 或 nonlinear equation row-scale 的内部说明，而不是
  继续盲扫 ErrRef、AreaFactor、Method 或线性预条件参数。

## 可复跑工件

- 审计脚本：`scripts/audit_templates_ldmos_g3_element_edge_current.py`；
- 回归测试：`tests/regression/test_audit_templates_ldmos_g3_element_edge_current.py`；
- ignored oracle：
  `reference_staging/templates_ldmos_g3_element_edge_current_20260901/`；
- 汇总：上述目录的 `analysis/summary.json`；
- 能力探针原始日志：上述目录的
  `remote/probe_high_vsv_node3721_minus_edge_api.log_des.log`。
