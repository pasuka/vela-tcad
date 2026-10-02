# SimpleMOS Enormal 恢复的只读前置核对

日期：2026-09-12。本记录只核对模型定义与当前实现边界；不表示 Enormal 已恢复、校准或自洽通过。当前先完成 [PhuMob 生产迁移](simplemos_split_state_production_validation_2026-09-12.md)后的四条完整曲线。

## 已核对的定义

本地 Sentaurus Device T-2022.03 用户手册（D:/software/ATCNIS01-202203/sdevice_ug.pdf），第 403–407、426–427 页：

- 原脚本 `Enormal` 默认采用 Enhanced Lombardi。其法向场垂直于半导体—绝缘体界面；`ToCurrentEnormal` 是另一种选择，不能静默替换。
- 声学声子与表面粗糙散射按倒迁移率相加，并与 bulk PhuMob 合成。距离衰减为 exp(-distance/lcrit)。保持 300 K 及原生默认参数。
- 默认使用场绝对值；只有显式启用 `MobEnormalSignDependence` 才改成符号相关。
- 默认几何距离为点到界面的真实距离，法向为该距离的梯度；`-GeometricDistances` 才启用最近界面节点的近似。
- 界面顶点的电荷场修正由 `NormalFieldCorrection` 控制，默认系数为零；不要将非零修正当作原脚本隐藏必选项。
- 这些公式不直接规定全部原生节点/单元/边投影次序；仍需用原生有效单元迁移率和固定状态导出来校准。

## 当前源代码的实际边界

`src/physics/MobilityModel.cpp` 已有 `phumob_lombardi` 与 `phumob_field_lombardi`，组合顺序为 PhuMob → Lombardi → 高场饱和；参数默认值也已存在。`include/vela/equation/AssemblerUtils.h` 的旧表面支路使用最近界面几何和实时三角形线性电势梯度，避免把场固定在初始化状态。

但当前合格的 `element_box_phumob` 分支只支持 plain PhuMob，且在 `edgeMobility` 中直接返回 box 结果；新 `SplitDDRuntime` 和宽精度 `SplitDDOperator` 也明确限制 plain PhuMob。因此，仅修改模型字符串既不能得到已校准的表面迁移率形成方式，也不能证明新状态低位和导数保持一致。

旧 M7 合同中的相邻输运单元算术平均是历史实现定义，不能替代本轮已经校准的 box 权重。旧 M7 的语法支持/局部测试通过也不等于本轮生产组合的 Enormal 自洽验收。

## PhuMob 完整曲线通过后的门槛

1. 固定原生与 Vela 同一状态，分别核对界面几何、法向场、距离、Lombardi 分量、有效单元迁移率及 box 边约简，覆盖两个 NWell 和两个 Vd。先查明形成次序，再实施一致的候选；不混入高场饱和。
2. 将场、载流子依赖及边权重纳入统一状态和端口路径，验证低位分拆不变性、亚 ULP 方向、所有场邻接列及载流子交叉块；弱块采用不会被总残差消减吞没的方法。
3. 在相同参数与门槛下进行原生小幅度正负响应校准及双初始化自洽对照，保留全部失败尝试。先控制点再完整曲线，之后才恢复高场饱和。

本记录没有新增 Enormal 输入上传、模型修改或仿真。PhuMob 曲线正在使用冻结源代码和程序，期间不修改求解实现。

## 原生场形成方式的可观测接口

进一步查阅同一手册 Chapter 39 的标准 `PMI_EnormalMobility` 和顶点 Compute Scope：`Compute_muinv(dist,pot,enorm,n,p,t,ct,...)` 接收实际法向场及距离；`ReadCoordinate`、`ReadNearestInterfaceNormal` 和 `ReadDistanceFromSemiconductorInsulatorInterface` 可在 Compute 内提供位置和界面几何。必要时可设计一个返回零附加倒迁移率、所有附加导数也为零的诊断 PMI，并先独立验证其加入前后的状态和 Id 完全一致，再用它观测原生输入。

这是可行接口的文档核对，当前尚未编译、上传或运行该 PMI，也没有将“理论零附加项”当作已验证无扰动。无需仅从有效迁移率反推多个未知的场/距离组合。手册 Math 参数索引也确认 `eMobilityAveraging=Element` 为默认；具体节点、单元、边形成次序仍按实测校准。

## 待验证观察器草稿

官方接口和两份示例已经从 T-2022.03 安装目录下载到忽略目录 `build-release/enormal_restore_20260912/vendor_reference`。自写草稿为 `build-release/enormal_restore_20260912/pmi_vela_enormal_observer.C`，使用标准 PMI 接口：主返回值、六个必需导数及 Na/Nd 可选导数均为零；按坐标和载流子缓存最后一次调用，析构时导出实际场、距离、界面法向和状态。

目前未上传、编译或运行该草稿。后续必须先证明加观察器前后的原生状态/Id 不变，并将其最后一次调用的 potential/n/p 与最终 TDR 逐点匹配；不能直接假设最后一次调用就是合格状态。`ReadCoordinate` 和 `ReadDistanceFromSemiconductorInsulatorInterface` 的文档单位均为微米；Compute 参数中的距离和场暂用 native 名称记录，待独立换算校准。

可优先利用观察到的每个顶点实际法向场及距离，核对节点 Lombardi 合成后的单元平均、再到已验证 box 边权重的形成次序；不预设旧的单元质心场方案与原生相同。仅当这些观测通过资格后，再决定统一宽精度状态与全部邻接导数的实现。
