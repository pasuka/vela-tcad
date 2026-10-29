# BVmethods 安装库原版：Q0 审计与 Q1 原生参考

日期：2026-10-15。用户要求利用已有原版结果继续，已恢复执行本轮 Q0a、Q0b 审计／工具及 Q1。
**现有原版物理模型在该网格上从零纯 QS 可达到声明停止条件；一次新的独立运行完整复现旧归档。**
Vela 尚缺 NonLocalPath 实现，未取得原版同模型扫描资格；旧 E2／Eparallel 的原生失败不再作为当前原版的失败结论。

## 输入冻结与配置

来源为安装库 `GettingStarted/sdevice/BVmethods` 的 `resistor` 分支，Workbench 节点 5。
旧结果来自 T-2022.03-SP2 的 2026-08-18 归档；新运行复制其实际展开的 cmd、par 和网格，
在独立目录从零初始化，无 Load、瞬态、辅助栅、模式切换或新增精度／PMI／盒文件控制。

| 输入／工具 | SHA-256 |
| --- | --- |
| 实际展开的 cmd | `30f025ed061f0c113cc4e1caa7c00f56003a6ac9ad01a4199e014b08602108a5` |
| 实际参数覆盖 par | `618d0467ee2608a76c7d3bf1dead41d5cc4709c73ad8dcfb59e341636c2630f3` |
| 原网格 | `7cd3968f5dd0e3cd63ef62ed0fcced465dd33bea465167f0be8139c16025902e` |
| 新运行及默认导出使用的 SDevice 可执行文件 | `c786dbd60dd26ca7e3671518256585fd0ae5bd615eaa7c642cbd80743417f71a` |

三份输入与旧归档逐字节一致，新运行前后散列不变。原 `BV_des.cmd`、`sde_dvs.cmd` 与
当前安装库字节相同；实际 SDE 等于原版的节点编号展开加外侧空白变化。旧元数据的轴转置
备注不能直接作为此原版被改写的证据，本轮未改旧清单。

网格有 2,600 个全局点、5,020 个三角形；Si 场数据含 2,152 个节点，不能混用全域与子域数量。
300 K，Fermi＋OldSlotboom，DopingDep／GradQuasiFermi HFS／Enormal，SRH(DopingDep)、
NonLocalPath、Avalanche(GradQuasiFermi)。漏极固定 `R=1e7 Ω·μm`；其他接触外电压为零、
栅 Barrier −0.55。外压目标 100000 V，电流停止声明 `1.3e-3 A/μm`。

新原生日志实际记录 Bank/Rose 非线性求解、blocked decomposition 线性求解、
64-bit double normal precision，装配／求解各一线程；未把 SHEDistribution 的 Super
日志行误作本次 DD 后端。旧扩展精度／DensityIntegral(240) 专项证据不直接成为此默认路径资格。

同一当前可执行文件下导出 Silicon、Oxide、Nitride 默认参数，分别核验并封存，显式 SRH
覆盖的七组数值与安装库参数一致。这些当前默认快照不冒充旧运行当时未保存的完整默认转储。
本轮确认的 Si NonLocalPath 默认包含 Apath `4e14 cm^-3 s^-1`、Bpath `1.9e7 V/cm`、
Ppath `0.037 eV`；原生参数说明以 Ppath>0 区分间接过程，实际路径族和端点算子留待 Q2/Q3。

## 独立轨迹与 BV 提取

| 项目 | 2026-08-18 已有结果 | 2026-10-15 独立复核 |
| --- | ---: | ---: |
| 工具 | T-2022.03-SP2 | T-2022.03-SP2 |
| PLT 点数／数据列 | 66／33 | 66／33 |
| 末外电压 V | 13806.3160273274 | 同左 |
| 末内电压 V | 7.37478070079669 | 同左 |
| 末漏极电流 A/μm | 0.00137989412466267 | 同左 |
| 电流停止日志 | limit=1.3000e-3；正常结束 | 同左；进程退出码 0 |
| 独立线性插值 BV V，abs(Id)=1e-4 | 6.673594818971575 | 同左 |
| 原生 SVisual／原版展开提取器 | 本轮输出 `DOE: BV 6.67` | 本轮输出 `DOE: BV 6.67` |
| SDevice wallclock s | 327.38 | 190.77 |

两份 PLT **逐字节相同**，全部 66×33 数据的最大绝对差为 0。每点从零来源、固定接触、
单调外电压、负载、KCL、位移电流和电流分量闭合均通过；不是仅用正常退出或终点判定。
负载最大差约 `1.2167e-8 V`，四端 KCL 最大绝对差约 `1.67384e-18 A/μm`。
原版提取器按 `%.2f` 输出，原生核验限于该显示精度；不声称已证明其完整精度输出等于
独立线性插值结果。旧 O-2018 参考的 6.37979 V 不用作本套 2022 原版的预设门限。
两次运行的耗时仅作记录，运行环境状态未控制，不作性能提升结论。

新终点 TDR 的七类已检查源项／电离数组与旧终点相同、全部有限；Band2BandGeneration
在 2,152 个 Si 节点中有 840 个非零样本。终点是约 1.38 mA/μm，而非 0.1 mA/μm
击穿提取点，不得拿它代替阈值场量。补充 HDF5 读取存在 h5py 构建／运行库版本提示，
本轮读取成功；仅作为数组身份审计，不认领正式 TDR 导入或跨工具场量资格。

## 工具、测试与能力边界

新增[原版参考校验器](../../scripts/validate_bvmethods_original_reference.py)，显式接收
cmd、par、mesh、PLT、log、冻结参考及网格散列，绑定文件声明、模型／接触／Math／Solve、
实际工具版本及停止记录，输出具名原版资格与输入散列。旧 QS 校验器仅增加显式 stop／resistance
可选参数，默认仍为旧 WP14 契约；没有改全局 STOP，也没有将旧结果重新命名为原版通过。

新增[11 项契约测试](../../tests/regression/test_bvmethods_original_reference.py)，覆盖模型／停止值
替换、重复 Physics、额外精度、重载／辅助栅、中间态失败、参数／网格／文件绑定及工具版本
篡改，并确认 NonLocalPath 仍被前端明确拒绝。新测试注册为 `bvmethods_original_reference`。
本机 70 项相关 unittest 与 Release 配置下五项注册 CTest 均通过；新校验器 CLI 实跑新归档通过。
首稿错误使用 IR 的 code 字段、测试临时目录 ACL 失败均修正，记录保留；最终未放宽预算。

实测解析实际原版 cmd 报 E106、Band2Band 状态 unsupported、无有效 Vela 映射。
[C++ 核心](../../src/physics/BandToBandTunnelingModel.cpp)仅接受 none／e2／sentaurus_e2，
当前数据类型与接口是局部场率和局部导数，尚无此 NonLocalPath 路径算子。
Avalanche(GradQuasiFermi) 有前端映射，实际 SG 电流重建／源支撑／交叉导数仍需原版同态核对。

本轮只配置 CMake 以注册 Python 测试，没有生产 C++ 改动、生产编译、全 CTest、发行包、
提交或推送。配置确认 HDF5、UMFPACK 等可用，但本轮未运行 Vela 线性求解，不能据此认领
某后端求解资格；Windows／Linux Vela 数值重复性留待 Q5。

## 状态与下一步

- Q0a：原件、实际展开、网格、显式参数及当前三材料默认冻结完成。
- Q0b：能力审计、原版专用校验与回归完成；模型真实实现属于 Q3，仍待办。
- Q1：原版已有参考及一次独立从零复核完成，具有原生重复采集证据。
- Q4：仅原生重复轨迹前提已有证据，整体 Vela 同模型资格仍未通过；Q5–Q7 未执行。

下一步先根据手册和冻结参数明确 NonLocalPath 的路径／端点支撑、非局部源积分及导数
契约，形成 Q2/Q3 实现与独立原生控制；获得实际模型能力后再比较完整残差／Jacobian。
旧 E2 弱行两个负向诊断继续暂缓，不把原版成功误解成 Vela 已支持原版。

后续手册与联网研究已形成[NonLocalPath 开发子计划](../superpowers/plans/2026-10-15-nonlocalpath-btbt-development-plan.md)。
手册 p.538 明确 `Band2BandGeneration` 等于空穴源；上述已有数组身份检查不代表电子端验证。
下一步 N0 补双源／版本／导数控制及预算，模型实现与资格仍待办。

所有远端路径／PID、原始商业输入输出、运行和测试日志、失败记录、逐文件哈希及私有脚本
封存在 ignored `build-release/bvmethods_original_q0_q1_20261015/`；旧原版11文件核验见
`build-release/bvmethods_existing_results_20261015/`。新扫描16文件、提取／Si 默认辅助19文件、
其他材料默认7文件均已取回核验，公开文档不保存商业原件或主机身份。
