# BVmethods NonLocalPath N1：接口与单位

日期：2026-10-16。目标是安装库原版 BVmethods 所需的 Si 声子辅助 NonLocalPath。
状态：N1 完成，本机与 T470p Release 验收通过。
本轮只建立私有 C++ 契约；路径搜索、WKB 内核、实际装配和导数计算分别归 N2–N5。

## 接口边界

新增 [NonlocalPathBandToBand.h](../../include/vela/physics/NonlocalPathBandToBand.h) 与
[实现](../../src/physics/NonlocalPathBandToBand.cpp)，均位于 `vela::experimental`。
实验入口是 [Catch2 契约程序](../../tests/test_nonlocalpath_contract.cpp)，CMake 目标
`test_nonlocalpath_contract`／CTest `nonlocalpath_n1_contract`。它直接调用私有 API，
没有新的公开求解命令或发行包能力声明。

- 配置要求显式 Si、区域列表、声子能量、最大长度，并且 A/B 与微观输入恰选一种。
  私有 JSON 严格检查字段及单位；非法整数、非有限数、重复区域、混合输入和未知控制失败。
  参数来源单独保存为 `TunnelingParameterProvenance`，不能混入运行时配置。
- A/B 是 F0=100 V/m 归一化的非局部参数，A 为 m^-3 s^-1，B 为 V/m，
  与局部 E2 的 A 单位不同。参考带隙必须显式给出，不能取有偏置态 BGN 后带隙。
  `mass_ratio` 为 mv/mc，零值选择手册自动质量关系。
- 私有转换可选 CODATA2018 或 N0 已核验的 T-2022.03 常数集；不会修改全局常数。
  微观输入质量为 kg，耦合参数为 J² m/kg；原生单位换算分别是 A×10^6、B×100、
  Cpath×0.01、m/m0×所选 m0、eV×所选 q。声子占据在 A/B 转换中使用参考 300 K，
  实际物理温度来自带状态，后续内核不得将参考温度当作当前温度。
- 带状态包含有效带边、梯度、两类 Fermi 能、人口、DOS、温度、BGN、材料／区域及接触身份，
  并携带这些量对物理 SI 变量的导数。这里只规定实际装配态接口，尚未接入装配器；
  禁止以降精度输出反算后冒充内部态。当前标量类型沿用核心 `Real`。
- 路径记录起点／折点／终点、单元与区域、重心支撑、反射、能量范围、积分区间和终止原因；
  已检查有限值、权重、长度、重复节点、区间重叠和状态身份。
  自主几何生成及 cell/节点拓扑一致性仍由 N2 实现，未宣称当前就是原生离散规则。
- 双源按同一个有符号事件保守分配，空穴在起点、电子在终点，保留负生成（复合）与接触事件。
  节点积分单位 m^-1 s^-1，与节点体积率 m^-3 s^-1 分开；积分除以面积才得到体积率。
  电流换算 `q × sum(integrals) × 10^-6 × AreaFactor` 为 A/μm；此量是源预算，
  后续须结合自由行与接触电流账本验证 KCL，不能直接认领终端电流。
  守恒误差除以绝对事件总量，沿用 N0 的 10^-12 门限，没有 max(1) 单位地板。
  零支撑面积只允许零源，绝缘体不要求虚构半导体 DOS；源校验还会复核逐节点事件账本，
  同时删去两类接触／自由节点源不能仅靠全域电子／空穴总量相等通过。
- 导数明确区分“尚未评估”、原生默认、原生显式及完整数学模式。
  存储导数不等于已经实现计算；未标记计算模式的非空导数被拒绝。
  稀疏支撑独立保存实际行列号，保留跨节点列；排序／去重不改变结构版本，增删或尺寸改变
  均更新版本，缺项抛出重建错误，不静默丢弃。
- 状态身份分别记录网格、材料、参数、接触、人口、温度及试探态版本。
  未来装配适配层必须在 ψ／QF／人口变化及线搜索试探、拒绝、接受时更新 trial／相应版本。
  当前接口检查版本相等，尚未自动追踪外部装配器的修改。

MVP 求解模式只允许 CoupledDC；Gummel、Transient、AC、Noise 明确拒绝。
该检查只是接口适用范围，并不意味着 CoupledDC 已有非局部求解实现。
公开 simulation schema、前端路由及 E2 生产实现均未因本轮修改；NonLocalPath 继续报 E106。

## 验证与证据

本机与 T470p 使用 Windows UCRT64 Release，核验同一最终构建／测试依赖集合。
远程独立作业源路径、PID、源码逐文件 SHA-256、构建／CTest／E2 日志均保存在 ignored
`build-release/bvmethods_nonlocalpath_n1_20261016/`，不把机器身份或仿真生成结果纳入版本控制。
最终快照身份：`31ac515b9ab9118b58c76e71746b67a4f058d1d57ef70cc7d58568969ae39f02`，2122 文件。
基础归档与接口／测试依赖增量共同组成最终快照；远端逐文件核验最终 manifest 后才开始构建。

| 检查 | 本机 Release | T470p Release |
| --- | --- | --- |
| `test_nonlocalpath_contract` | 11 用例／97 断言通过 | 11 用例／97 断言通过 |
| `test_recombination "*E2*"` | 3 用例／13 断言通过 | 3 用例／13 断言通过 |
| 选定 CTest | 6/6 通过 | 6/6 通过 |
| 五组 Python 回归（含在上述 CTest 中） | 81 用例通过 | 80 通过，1 项可选交叉检查跳过 |
| 新 C++ 文件 ASCII／空白与报告链接 | 通过 | 使用相同文件哈希 |

CTest 为 `nonlocalpath_n1_contract`、`bvmethods_nonlocalpath_contract`、
`bvmethods_original_reference`、`sdevice_cmd_parser`、`sdevice_contracts`、`wp14_zero_qs`。
远端未安装可选 `jsonschema`，仅该包的交叉检查跳过，标准库 schema 校验通过。
本轮未跑全量 CTest，未改变生产求解流程；覆盖本接口及 E2、前端能力门和原版契约。

两端配置实际检测到 HDF5/HighFive、SPQR/UMFPACK、STRUMPACK、SuperLU_MT、MUMPS 与 METIS；
这些是可链接能力，本轮没有选择某个后端求解器进行器件仿真。无网格、偏置或新全 BV 结果。
独立参数控制沿用 N0 的 Si 参考带隙 1.1082092307692308 eV、A=4e14 cm^-3 s^-1、
B=1.9e7 V/cm、P=0.037 eV及 T-2022.03 常数，核对微观质量与 Cpath、A 倍增和 B 平方质量响应。

远程初始后台启动未生成状态，后改为保留 SSH 会话执行。裁剪快照首次缺少消息目录文档与
参考清单，补齐后又发现归档时间戳导致 Ninja 复用旧测试对象；最终显式刷新三个变更输入，
日志确认 C++ 实现和测试均重编译，再取得 11 用例／97 断言资格。
E2 启动器的 Windows 相对路径随后修正为绝对路径，仅补跑剩余 E2；所有最终退出码均为零。
这些启动／测试环境失败证据保留在 ignored 目录，未改门限、未跳过必需测试或替换原版物理。

本轮无需网格／求解器后端运行，不认领原版全 RHS、路径、Jacobian 或从零纯 QS 资格。
后续 N2 先验证自主路径几何，再由 N3–N5 提供源项及完整导数。
