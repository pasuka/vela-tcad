# LDMOS Vg=8 精度修复版全曲线完成记录（2026-09-09）

Vg=8、Vd=0–40 V 的同版本检查点续跑于 2026-09-09 20:37:50
（Asia/Shanghai）完成。原 31 个精确参考点全部取得，工程和最终单栅压评分均通过，
完整性审计通过。它延续已审计的 0–5 V 前段，不是另一次从零开始的不中断运行。

## 配置与版本

冻结 Release 程序 SHA256：
`968ef00ed4cc16447c98853b4044f0895238f77f4972e18399694c7546cad97d`。
该结果对应合并本地 main 之前的精度修复版本；不能直接作为合并后程序的全曲线证据。
构建使用 MSYS2 UCRT64、C++20，实际线性后端为 Eigen SparseLU/COLAMD，L2 均衡。
网格包含 10241 个节点、19782 个三角形、30022 条边，5723 个硅节点。
温度 300 K，Fermi–Dirac 统计、OldSlotboom BGN、SRH/Auger，
constant_field / transport_cell_vector 准费米梯度高场迁移率；
雪崩、量子、热和预测器关闭。几何输入为 μm，密度为 m⁻³，电势为 V，电流为 A/μm。
原全局、局部行和 KCL 门限保持不变。修复及恢复资格详见
[精度修复记录](templates_ldmos_stage4_frame_fix_2026-09-09.md)和
[续跑记录](templates_ldmos_stage4_frame_resume_2026-09-09.md)。

## 结果

| 指标 | 结果 |
| --- | ---: |
| 精确参考点 | 31/31 |
| 电流相对误差中位数 | 1.330933% |
| 电流相对误差 P95 | 1.612869% |
| 最大电流相对误差 | 1.627051% |
| 低 Vd 微分电阻误差 | 1.032119% |
| 40 V 端点电流误差 | 1.099059% |
| 40 V 电流 | 0.000423633153 A/μm |
| 40 V 参考电流 | 0.000419027791 A/μm |
| 最大归一化 KCL | 1.373226492e-10% |
| 回退数 | 0 |
| 已记账 Newton 更新 | 5065，其中续跑新增 4398 |

审计统计 431 个接受状态、429 次传递和 860 个已完成子步骤。
此前人工暂停中断的 child_00122_direct 成本不完整，单独保留；上述更新总数
不能用于声称覆盖全部中断成本的性能结果。此前该修复版本完整 Release CTest
为 739/739 通过，续跑未改变其程序和冻结输入。

本结果完成 Vg=8 单栅压验收。最终版本双栅压 D5 仍需同版本 Vg=4 全曲线及
两栅压比值验收；旧二进制的 Vg=4 结果不能拼接成最终版本通过。

## 本地原始证据

以下原始输出保留在被 Git 忽略的 reference_staging 中，不随文档提交：

- [评分结果](../../reference_staging/templates_ldmos_frame_resume_20260909/vg8_full/score/summary.json)
- [完整性审计](../../reference_staging/templates_ldmos_frame_resume_20260909/vg8_full/audit_summary.json)
- [完成状态](../../reference_staging/templates_ldmos_frame_resume_20260909/vg8_full/progress.json)
- [最终账本](../../reference_staging/templates_ldmos_frame_resume_20260909/vg8_full/fixed/ledger.json)
