# BVmethods NonLocalPath N2：自主路径几何

日期：2026-10-16。状态：私有几何实现及 N2 几何门完成，本机／T470p Release 验证通过。
本轮不包含声子辅助 WKB、成对事件率、实际源装配或 Jacobian；这些仍归 N3–N5。
公开前端继续拒绝 NonLocalPath（E106），没有新增发行包功能声明。

## 实现与物理边界

新增 [NonlocalPathGeometry.h](../../include/vela/physics/NonlocalPathGeometry.h) 和
[实现](../../src/physics/NonlocalPathGeometry.cpp)，保持 `vela::experimental` 命名空间。
依据 T-2022.03 手册 p.530 的动态路径规则：从 Ev−Δv 出发，方向为起点 −∇Ev，
直线穿越三角形，在 Ec+Δc 与起点能量相等处终止；Neumann／Si–绝缘体界面镜面反射。
偏移只改变能量条件，不改变起始方向。

实现使用实际输入网格与二进制64位带边值，不读取原生路径、端点或拟合支撑表。
几何输入为 SI 米、焦耳的共形 Tri3 网格；纯几何拓扑只构造一次，路径由每次传入的状态重算。
校验重复／退化单元、非流形边、共享边两侧局部重叠、外边界分类、材料／选区及状态版本。
不宣称构造器执行了非相邻单元的全局交叉检测；调用方须提供共形网格。

离散约定与手册规定区分如下：

- 内点方向取本单元 P1 价带梯度；自动起点为每个活动 Si 顶点一次，节点方向由邻接
  Si 单元 P1 梯度按面积恢复。默认只使用活动区域，允许跨未选区时使用相邻全部 Si。
  这一恢复及起点采样是 Vela 的候选离散规则，尚未证明与原生逐路径规则相同。
- 通过重心坐标射线交点确定穿越，无空间 epsilon 探针位移。机器精度内的边／顶点归属
  使用固定 128ε 分类；端点与长度／边界重合时使用同尺度舍入判定，优先记录完整端点。
  实际点间弧长用于区间记账，避免理想射线参数与输出坐标舍入不一致。
- 共边掠行优先归属可进入的选中 Si，再按单元 ID 确定唯一归属；编号较小的氧化层／
  未选区不得抢占。默认路径跨未选区明确终止，关闭选区限制后允许穿越同材料 Si；
  不支持的半导体材料仍明确终止。
- 单元穿越保留每段两端重心支撑、材料／区域、累计弧长、反射点及状态版本。
  绝缘体不接收端点支撑。接触／开放外边界、零梯度、无势垒、物理最大长度、反射／遍历
  资源上限、角点歧义及数值失败分别记录，不能都当作零源。
- 节点增加 `additionalRegions`，可表达同材料共享顶点的完整区域成员身份。
  N1 路径校验允许非完整诊断记录仅含一个起点；完整路径仍要求至少两个点及有效区间。
  这是失败审计的数据表达扩展，没有放宽源守恒、残差或物理验收预算。

## 私有审计入口

[vela_nonlocalpath_geometry_probe.cpp](../../src/tools/vela_nonlocalpath_geometry_probe.cpp) 接收
`vela.experimental.nonlocalpath.geometry_input.v1`：SI 坐标、三角形、明确外边界、带状态及
N1 私有配置。Si 材料 ID 为0，其余材料不作为 Si 活动起点。
输出所有起点的 `vela.experimental.nonlocalpath.geometry.v1` 审计，明确
`geometrically_complete`、`resource_failure`、`requires_boundary_resolution`。
程序退出0表示审计完成，不能据此认领所有路径可用于源积分；资源失败退出3，非法输入退出2。
`EXCLUDE_FROM_ALL` 目标未安装到公开命令或发行包中。

```powershell
cmake --preset windows-ucrt64-release
cmake --build --preset windows-ucrt64-release --target test_nonlocalpath_geometry test_nonlocalpath_contract vela_nonlocalpath_geometry_probe
ctest --preset windows-ucrt64-release -R '^nonlocalpath_n[12]_'
./build-release/vela_nonlocalpath_geometry_probe.exe <geometry-input.json> <audit.json>
```

审计只使用带边参与几何求值，QF、占据与其导数不参与；不能将其解释成实际内部态／RHS1绑定。
五态输入的 DOS 在 TDR 的单位标签为 `?`，按手册表195的 cm^-3 转为 SI，并在输入来源中注明；
DOS仅用于状态合法性检查，不影响几何结果。

## 验证

两端均使用 UCRT64 Release，源码快照2126文件 SHA-256：
`3d88f5b7cb20b5925518f8e2ee9a4344aa0a96e31bbbd860af411c54778c06d1`。
远端复用 N1 隔离构建缓存，N2 的作业、源增量、完整文件清单及证据独立保存。
逐文件验证后显式刷新变更输入，构建日志确认实现／测试重新编译，避免归档时间戳复用旧对象。
远端源路径、PID、启动／退出码及原始输入输出保留在 ignored
`build-release/bvmethods_nonlocalpath_n2_20261016/`，不将机器身份或原始仿真输出纳入版本控制。

| 检查 | 本机 Release | T470p Release |
| --- | --- | --- |
| 新增 [几何性质测试](../../tests/test_nonlocalpath_geometry.cpp) | 14用例／8274断言通过 | 14用例／8274断言通过 |
| N1接口回归 | 11用例／97断言通过 | 11用例／97断言通过 |
| E2数值回归 | 3用例／13断言通过 | 3用例／13断言通过 |
| 选定CTest | 7/7通过 | 7/7通过 |
| 五组Python回归（含于上述CTest） | 81用例通过 | 80通过，1项可选jsonschema交叉检查跳过 |
| 实际五态网格审计 | 全部完成 | 全部完成，逐状态文件与本机字节一致 |

性质控制包括线性解析端点、SI尺度／能量规范平移、坐标旋转、网格方向反转、镜像对称、
平面镜面反射、Si／氧化层反射、可分辨掠边、顶点穿越、区域开关与材料限制、无势垒／零梯度、
端点与选区／长度边界重合、无重复起点或间断区间、资源上限与非法网格。
曲带控制 Ev=−x−0.2x²、Ec=Ev+0.7，x=0起点的解析端点为
`(sqrt(1+0.8*0.7)-1)/0.4`；nx=4、8、16、32、64嵌套细化的端点误差严格减小，且小于0.2h²/4。
私有CLI控制确认遍历预算1时退出3（2079条遍历超限）；负预算、未知选项、小数索引、
错误节点数组长度、负区域与负版本均退出2，不静默截断索引或接受非法输入。

CTest为 `nonlocalpath_n2_geometry`、`nonlocalpath_n1_contract`、`bvmethods_nonlocalpath_contract`、
`bvmethods_original_reference`、`sdevice_cmd_parser`、`sdevice_contracts`、`wp14_zero_qs`。
未跑全量CTest；生产装配与求解流程未改。两端实际检测HDF5／HighFive、SPQR／UMFPACK、
STRUMPACK、SuperLU_MT、MUMPS及METIS；本轮未使用这些线性后端求解器件。

## 五态实际网格审计及限制

复用N0冻结的零偏、首非零点、BV前／后、末态TDR：2600节点、5020三角形、2152个Si顶点。
使用导出的带边及明确的 **100nm几何诊断长度上限**；该长度不是原生参数同值资格声明。
没有再次运行Sentaurus，也没有更改原版模型、参考轨迹或N0预算。

| 状态 | 完整路径 | 长度上限 | 接触终止 | 零方向 | 歧义角点 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 零偏 | 76 | 1632 | 404 | 38 | 2 |
| 首非零点 | 79 | 1667 | 404 | 0 | 2 |
| BV前 | 776 | 1011 | 363 | 0 | 2 |
| BV后 | 774 | 1013 | 363 | 0 | 2 |
| 末态 | 790 | 999 | 361 | 0 | 2 |

每态2152个起点恰好覆盖活动Si顶点，无重复或遗漏，数值失败／资源超限均为0。
两端五份审计各自字节一致。独立Python检查从输入网格和每段权重重建位置、逐段弧长与端点能量，
位置相对舍入余差最大1.050e-15，完整端点能量相对舍入余差最大2.719e-15。
这些是几何表达检查，不能替代源项、完整RHS或Jacobian预算。

两条未解决路径始终为节点1287、2567，位置分别为(0,+0.165μm)、(0,−0.165μm)，
都是接触与反射边界交汇处的起点，记录只有1个点、0段。手册未定义该角点冲突的唯一物理规则，
本轮保守标记 `ambiguous_corner`，没有经验法线或接触优先补丁。
N3/N4接入源时须显式处理／拒绝这类未决路径，不能将其默认为零并宣称原版完整物理通过。

原生逐路径内部身份及离散源支撑仍未直接取得；本轮只证明自主几何性质、实际网格覆盖和两机一致性。
下一步为N3声子辅助内核；实际双源／源支撑对照和未决边界对源的影响继续纳入N4能力门。
