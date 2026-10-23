# BVmethods 工作区分组归档（2026-10-23）

本记录整理已有修改，不新增物理模型、接受策略或数值资格。当前入口是
[SDevice 交付指南](../sdevice_delivery.md)；历次研究报告中的“下一步”仅为当时记录。
逐文件分组和原始／公开版本哈希见[整理清单](bvmethods_worktree_organization_2026-10-23.json)。

## 三组提交边界

| 组 | 内容 | 边界及依赖 |
|---|---|---|
| 1：交付修复 | cmd 初始化、Fermi/NoFermi 与 OldSlotboom 映射、golden 和相关回归；表面迁移率先于 HFS | 可独立提交的公共修改。保留 D3 已验收行为；不是原版 BVmethods 全轨迹资格 |
| 2：实验实现及共享依赖 | NonLocalPath N1–N4、report/缓存/残差采集、保精度人口与 SG 工具、增广外电阻、共享装配器/Newton/输出修复和 C++ 测试 | 公共修复直接依赖实验采集头文件，作为同一可构建依赖集合保存。不能按整个提交宣布功能已完成；公开 NonLocalPath 仍不支持，report 不参与接受 |
| 3：诊断归档 | 可复现分析脚本及 Python 测试、版本化契约、研究计划、数值摘要和报告索引 | 原始日志片段、机器路径、作业身份不进入公开摘要。大规模运行输出及最新私有 R2 覆盖实现仍留在 ignored 目录 |

第一组单独不足以重建 D3 包；第一、二组构成运行时依赖集合。逐文件检查确认本轮整理前
`src/`、`include/`、`scripts/vela_sdevice/` 及 CMake 的实际文件与冻结快照
`6a56593bf3917cc213a00aeab4f4a614befa52ca` 一致（包含当时未跟踪文件）。
`tests/test_newton_solver.cpp` 另含 D3 已记录的测试分解修订；它不是运行时变更。
CMake 的 Python 诊断注册随第三组落地，第二组只登记对应 C++ 目标。

本轮不改运行时源代码。两处 Fermi/NoFermi 文档说明改为准确描述现有行为：未显式选择时
省略修正键并沿用 C++ 默认 false；显式 NoFermi 写 false。数值门和失败证据保持。

## 成果状态

- D0–D4 的具名交付资格见 [D1–D3](sdevice_delivery_revalidation_2026-10-22.md)及
  [D4](sdevice_delivery_d4_2026-10-23.md)。整理提交不重新授予发行包资格。
- NonLocalPath：接口、几何、声子内核和双源诊断已实现；N4 完整残差、导数、真实接受未收口。
- R2：最新两态、三个细档的指定方向 Jv 通过，实际全步仍拒绝；尚未持续重线性化迭代。
  以[最新 R2 报告](bvmethods_resistor_r2_mobility_contract_2026-10-22.md)为准。
- R3：增广预测／步长历史的严格重启 L2 和从零完整纯 QS 轨迹仍未合格。
- 实验 report、隔离回放及当前方向控制不替代生产行门、全局闭合或物理验收。

## 原始证据保护

整理前 433 个修改／未跟踪文件已逐字节备份到 ignored
`build-release/bvmethods_organization_20261023/pre_organization_sources.zip`。
归档 SHA-256 为 `eecbc9d0e8d27c0bbe612a9c235de9103fc293c77f5452808853bdcc646a8f81`。
公开摘要去除执行元数据和原始日志片段，保留数值结论、门限及原始产物哈希；清单区分
原始文件哈希和整理后哈希，不把脱敏版本宣称为历史原件。

原生商业输入、矩阵、回调、反汇编、机器配置、完整作业记录和大型结果保持在原 ignored
`build-release/bvmethods_*` 目录，不移动、不删除、不纳入本次提交。
原有私有 R2 代码是具名实验覆盖版本，不冒充第二组中的公共实现；其身份仍由各实验清单绑定。

另有六个远程工作流文件在本工作分支未跟踪，但已与本地 main 中的版本逐项比对一致。
它们不属于 BVmethods 三组，本轮原位保留，避免重复提交。

## 本轮验证

- UCRT64 Release 完整构建通过；实际检测到 HDF5 和 UMFPACK。
- 前端、材料、契约及 PN2D 模板：123 项，121 通过、2 跳过。
- 诊断 Python 测试：267 项。首轮 263 通过、4 项因未设置 `VELA_FERMI_AUDIT_EXE`
  报环境错误；由已归档 C++ 源码编译参考程序并设置变量后，4 项全部通过。
- 相关 CTest：475/475 通过，38.68 秒。覆盖 Newton、外电阻、迁移率、源项、report、
  NonLocalPath、外部状态及归档等；未重新运行全量 CTest。
- 未注册到 CTest 的载流子统计可执行程序：20 例／989 断言通过；私有几何：9 例／23 断言通过。
- 887 个原有文档相对链接有效；分组 JSON 可解析；原文件备份 433/433 哈希通过。
- 首次 CTest 的测试清单相对路径解析失败，没有执行测试；改用绝对清单路径后通过，首轮日志保留。
- 本轮没有启动新的 Sentaurus 或远端 BV 仿真，没有重建发行归档。运行时与 D3 快照相符，
  不把这些回归结果解释为原版 BVmethods 全轨迹通过。

最新私有 R2 覆盖实现另有 17 个源码／脚本文件备份至 ignored
`build-release/bvmethods_organization_20261023/r2_private_candidate_sources.zip`，清单记录逐文件哈希。
该版本依赖 `private/bvmethods` 头文件和原始绑定态，因此仅作为私有实验保存，不安装或替换公共装配器。
