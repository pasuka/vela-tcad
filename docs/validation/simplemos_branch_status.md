# SimpleMOS 本地分支整理与验证入口

更新：2026-09-08。工作分支为 `codex/simplemos-sdevice-validation`。

本次整理保存了 M78 以来的生产实现、诊断脚本、测试和验证文档。当前 n23 电流结果见 [Id–Vg 状态表](simplemos_n23_idvg_status_2026-09-08.md)。旧的日期报告保持原文；其中“未提交”等描述只代表当时状态。

## 生产实现范围

- 装配器、接触电流和探针共享有效输运边几何；Poisson 的电子、空穴和掺杂电荷可显式选择独立 signed Si 体积，连续性/SRH 体积保持原策略。
- Poisson、Gummel 的 Poisson 子步、耦合 Newton 和电极反应电荷共享逐单元介电系数装配。
- 显式支持二维 Delaunay box 转移，以及 constant/Masetti 的单元到边迁移率平均；不支持的几何、场依赖模型或冲突配置明确拒绝。Gummel 载流子路径未实现该输运平均。
- 等 ni、无 BGN 的稳定 SG 电势偏导；外部六列状态按物理参考值先相减后缩放；可选高精度线性缺陷迭代修正。
- 补齐从保存态求解到端口提取器的几何配置传递，并增加只读单元 box 和端口响应诊断。

接口与限制见 [配置说明](../config_schema.md)。新几何及迁移率组合均需显式启用；没有改写全局模型默认、物理常数、SRH 默认或收敛门槛。PN2D BV 模板的既有原子组合不变。

已验证实现的主要记录：

- [生产几何、状态精度与数值移植](simplemos_production_migration_and_native_mobility_response_2026-09-07.md)
- [逐材料介电装配](simplemos_permittivity_production_validation_2026-09-07.md)
- [自动 box 与单元/边迁移率](simplemos_generated_box_mobility_validation_2026-09-07.md)

## 仍属于诊断的方案

手册常数、局部 SRH 体积修改、稳定范数差比较及 ψ 小步余量累积在隔离程序中验证，没有合入生产默认。

最近的 58 次 DC 对照中，稳定范数差的控制通过，但高 NWell 原始失败路径仍未通过；ψ 余量累积使原失败路径变差；SRH 校准为 10/20 个 DC 合格、0/8 组响应通过。n23 低 Vd 约 0.104723% 的小源响应偏差几乎由剩余空穴残差解释；这不是 Id–Vg 对 Sentaurus 的电流差。[最近报告](simplemos_stable_merit_step_and_srh_calibration_2026-09-08.md)记录了全部资格边界和后续方向。

当前不能把新版四个 n23 点的约 0.006% 误差外推为完整 0–1 V 曲线，也不能把旧 PhuMob/Lombardi/HFS 全曲线与新 Masetti 子集混成一条曲线。局部响应及初始化资格仍限制进一步物理候选放行。

## 本次提交前检查

- `cmake --build --preset windows-ucrt64-release --parallel 4` 成功；构建已是最新。
- 完整 `ctest --preset windows-ucrt64-release --parallel 4`：764/780 通过，227.54 秒。失败名称与前一轮的 16 项完全相同：13 项历史源码哈希/替代链检查，3 项引用缺失的历史 `tests/test_mos_mixed_material.cpp`；没有删除测试或更新旧哈希。
- 16 个新增诊断 Python 测试模块共 72 项：71 通过，1 项错误。M79b 的 `test_original_failure_record_is_retained` 调用旧 M79 `verify()`，要求当前 `build-release/vela_example_runner.exe` 等于当时二进制哈希；生产程序演进后不满足该历史身份要求。该测试自身也在 M79b 原合同中冻结，本次保留测试、旧身份门槛及失败日志。
- 新增 Python 脚本和测试的 180 个文件语法检查通过；n23 表格重新核对源数据、百分数及偏置覆盖。最近两组独立 C++ 诊断测试已有 5 项、136 条断言通过的冻结记录。

完整测试因此仍非全绿。上述历史身份检查的接续是后续独立工作，不代表已证明所有未测模型均正确。当前构建 HDF5/TDR、UMFPACK/SPQR 可用，SimpleMOS 本次已存 DC 使用 Eigen SparseLU/COLAMD。

本次检查日志位于忽略目录 `build-release/simplemos_commit_20260908/`。这些日志及结果不纳入本次 Git 提交。

## 代码与实验数据的保存边界

本次提交包含源码、CMake 配置、测试、诊断/分析脚本、当前入口和日期报告。新生成的 `reference_tcad/simplemos_sentaurus2022/` 机器证据、CSV、状态和原始仿真输出保留本地；`build-release/` 内二进制、网格、状态、TDR 和日志继续忽略，未添加到 Git。

因此，在只检出本分支的新目录中，部分日期报告的机器证据链接及实验脚本依赖尚不存在。完整重放需要原本地数据包、相同源/二进制前像及对应环境；不能仅凭检出脚本宣称可复现全部实验。脚本中的冻结哈希与拒绝覆盖检查应保留，新的实验应使用新目录和新合同。

本次 n23 文档包含可直接阅读的四点表与两条历史完整曲线的 21 点误差表；其源 CSV 的 SHA256 也随文档保存。未上传代码、推送远端或合并其他分支。
