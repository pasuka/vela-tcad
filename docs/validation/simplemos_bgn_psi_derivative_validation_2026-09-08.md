# SimpleMOS BGN 电势导数精度修复与验证（2026-09-08）

工作分支：`codex/simplemos-sdevice-validation`。本轮承接 [OldSlotboom 恢复报告](simplemos_model_restoration_2026-09-08.md) 中的电势方向 Jv 失败，先解决该资格缺口。

已定位到可变本征浓度 SG 电势导数的双精度消减。原展开公式在代数上成立，但近平衡时两个大项相减会产生错误的微小导数。生产修复对状态无关迁移率的 BGN 分支改用稳定边通量求导；原残差、端口电流、材料、几何和接受门槛保持原值。PhuMob、Enormal 和高场饱和仍需要各自的本构及离散校准。

## 独立公式检查

耦合 Newton 从总杂质浓度（施主加受主）建立 BGN 有效 ni，并向本构函数传入零载流子参数。因此 ni 在节点间可变，但该路径中不随 Newton 状态变化。本轮没有缺失的载流子相关 ni 链式项需要补入。

在四个已冻结的 Vg=0.8 V 状态上，对每个载流子选取旧展开式与稳定式双精度差最大的 12 条硅边，检查两个端点，共 192 项。选择按差值和边编号确定，不按检查是否通过筛选。独立参考为 **100 位十进制的未因式分解密度通量中心差分**，步长 1e-25 V。

| 求值方式 | 相对独立高精度差分的最大误差 |
|---|---:|
| 原展开式，100 位求值 | 2.64126302e-53 |
| 原展开式，双精度求值 | 4.36335006e7 |
| 稳定因式形式，双精度求值 | 3.95253359e-8 |

192 项均为非零准费米差；固定门槛分别为高精度公式 1e-30、稳定双精度式 1e-7，全部通过。这里把导出的物理电势、准费米势和 ni 的 binary64 数值精确提升为 Decimal，系数取 1。它用于区分代数遗漏与求值消减，并不是内部拆分准费米参考的逐位重放；生产 Jv 另行验证。

令 `F` 为残差实际使用的稳定边通量，`L=B'/B`，其中电子的 Bernoulli 参数为 `-eta`、空穴为 `+eta`，则固定 ni、固定迁移率的端点电势偏导为：

- 电子：`(F/Vt) * [1+L, -L]`。
- 空穴：`(F/Vt) * [-1-L, L]`。

平坦准费米势下 F=0，两个端点偏导也严格为零。新增装配测试在修复前捕获到应为零的元素约 −2.3579131812e-7；修复后通过。另一项 100 位独立测试覆盖 ni 比值 0.001/3/1000、正负电势差、零及极小准费米差和边反向的符号/端点一致性。

## 生产范围与回归

修改位于 [CoupledDDAssembler.cpp](../../src/equation/CoupledDDAssembler.cpp) 和 [稳定导数接口注释](../../include/vela/discretization/StableSGDerivative.h)。BGN 稳定偏导只用于状态无关迁移率；保留已有 ±500 人口指数保护。未改变 `equal_ni_flux_evaluation` 的无 BGN 语义。

首个较宽候选的 SimpleMOS 数值检查通过，但完整回归发现：它改变了一个高场模型刻意关闭迁移率链式导数的历史近似支路，使冻结的 `max_iterations` 轨迹变成 `line_search_non_decrease`。本轮将修复限定到已验证的状态无关迁移率范围，保持该高场近似协议及冻结期望不变。首个候选的源码、二进制、全部 48 次 DC 和失败回归日志已单独归档，没有覆盖成最终版本。

同时修复了 DCSweep 测试临时目录的并行竞争：先检查目录不存在再使用，可能让两个 Windows 进程选到同一路径；现在使用 `create_directory` 原子占用。它只影响测试隔离，不改变求解器接受条件。

最终 Release 构建成功。完整 `ctest --preset windows-ucrt64-release --parallel 4` 为 **766/782 通过，230.02 秒**。16 个失败名称与此前完全一致：13 项历史源码身份/替代链检查，3 项依赖缺失的历史 `tests/test_mos_mixed_material.cpp`。两个新增 BGN 数值测试及上述高场冻结、并行目录相关测试均通过；完整测试集仍非全绿，没有重写旧哈希或删除失败记录。

## 自洽与 Jv 复核

最终范围为 n19/n23 × Vd=0.05/1 V：OldSlotboom 取 Vg=0/0.2/0.8/1 V，共 16 点；无 BGN 控制取 Vg=0.8/1 V，共 8 点。两条初始化路径各计算一次，**48/48 次首次尝试通过，24/24 点双初始化通过**，没有使用预留的同偏置恢复重载。较宽候选的另 48 次计算独立保留，不混入这组计数。两版累计 96 次 DC，本轮均无 DC 失败；较宽候选的回归失败仍保留。

最差载流子逐行比值 3.76905374e-7，低于冻结的 1e-6 门槛；最大 KCL/Id 为 4.91554685e-12，最大端口提取相对差为 1.13242749e-14。双初始化最大 ψ/φn/φp 差分别为 2.90e-15/3.820699e-11/5.0838023e-10 V，最大载流子相对差为 1.96650265e-8，最大 Id 相对差为 4.92639263e-12。

**全局生成复合源闭合仍有资格边界。** 48 次计算的电子/空穴共 96 个源积分分量全部低于原合同 1e-10 下限，因此都走原有“不激活全局相对源闭合门槛”的分支。不能称这 96 项独立通过了 1e-6 的相对源闭合测试。记录的最大带下限比值为 9.13118815e-5；该项净端口通量为 −3.23432200e-12、源积分为 −3.24345319e-12（诊断残差单位）。本轮保留这些数值和原门槛，依据另外实际执行的逐行及端口检查判定；后续微小 SRH 响应仍需单独资格验证。

Jv 覆盖四个原始冻结状态和 24 个重算状态。对每第三个自由硅节点，分别扰动 ψ/φn/φp，幅度为 1e-4、3e-5、1e-5、3e-6 V，共 336 个方向、1008 个方程分块记录。相对误差按 `||Jv-FD||/max(||Jv||,||FD||)` 计算，撤掉探针原输出中的 1 范数下限。两个小幅度的非载流子交叉块为正式门槛：**392/392 通过，最大 2.91628881e-8，门槛 1e-4**。电子↔空穴弱源交叉块单独标记，不据全残差差分宣称其完整性。

以下是同一原状态、Vg=0.8 V、h=1e-5 V 的直接前后比较（均为相对误差，不是百分数）：

| 工况 | ψ→电子，修复前 | ψ→电子，修复后 | ψ→空穴，修复前 | ψ→空穴，修复后 |
|---|---:|---:|---:|---:|
| n19，Vd=0.05 V | 9.74343e-6 | 2.12539e-8 | 约 1 | 2.53131e-8 |
| n19，Vd=1 V | 6.97116e-7 | 1.92025e-8 | 约 1 | 2.43800e-8 |
| n23，Vd=0.05 V | 1.14180e-3 | 2.23955e-8 | 约 1 | 2.87235e-8 |
| n23，Vd=1 V | 2.15389e-4 | 2.27826e-8 | 约 1 | 2.55120e-8 |

在相同原状态上，h=3e-5/1e-5/3e-6 V 的非交叉分块最大误差依次为 2.56476e-7、2.87235e-8、3.75314e-9，体现中心差分步长收敛；先前较小三幅度的 18/84 项失败现在为 0/84。

最终 Id 相对修复前的最大变化仅 6.70610234e-11。本次修复改善了导数求值精度，没有解释或消除原有深关断漏电差。当前 **OldSlotboom+Masetti、基础 ni 已匹配** 的 n23 结果如下，误差为 `(Id_Vela/Id_Sentaurus−1)×100%`，不是此前无 BGN 完整曲线：

| Vg / V | Vd=0.05 V | Vd=1 V |
|---:|---:|---:|
| 0 | −0.02132498% | −0.88268614% |
| 0.2 | +0.00352636% | −0.02084582% |
| 0.8 | +0.00569486% | +0.00652913% |
| 1.0 | +0.00504174% | +0.00504624% |

高 Vd、Vg=0 V 时，Id 分别为 Vela 3.37084893345e-16、Sentaurus 3.40086792323e-16 A/μm，绝对差约 −3.00190e-18 A/μm。本轮 BGN 只覆盖上述 16 点，未重算 BGN 的完整 51 点曲线。

## 配置和重放边界

复用上一轮原生导出与冻结网格：n19 为 1480 节点/2742 三角形，n23 为 1482 节点/2746 三角形。材料温度 300 K，Boltzmann 统计、Masetti 总杂质依赖、掺杂相关 SRH；BGN 组为 OldSlotboom，硅基础 ni=1.4638805412559193e10 cm⁻³。无 BGN 控制保留其原匹配材料。配置中的历史 `simplemos_m65` 注释标签不是实际材料值，实际值来自冻结的 `materials_file`。

几何为 `element_box` 输运和迁移率平均、`delaunay_transfer`、signed Si Poisson 电荷体积、逐单元介电系数；SRH 体积保持原策略。使用 `unit_scaling`，内部长度/浓度为 cm/cm⁻³；报告电压为 V、电流为 A/μm。网格和诊断文件的显式单位按各自契约解释。构建支持 HDF5/TDR、UMFPACK/SPQR；本次 DC 实际使用 Eigen SparseLU/COLAMD，线性迭代修正 4 次。

Newton 最大 200 次迭代、reltol=1e-7、abstol=1e-12、准费米更新上限 0.025 V；载流子逐行门槛 1e-6，关闭逐行源项/通量人为下限，每次检查 1814 个自由硅载流子行。求解器内的全局闭合开关为 off，后处理依原合同在源积分达到 1e-10 时才要求全局相对闭合不超过 1e-6；本轮该条件全部未激活，详见上文。另检查 KCL/Id 和端口提取一致性均不超过 1e-8。双初始化要求三种势最大差均不超过 1e-6 V、载流子相对差不超过 1e-4、Id 相对差不超过 1e-6。没有降低任何门槛。

本轮采用冻结的原生参考完成局部修复验证，没有新启动 sdevice。后续顺序仍为 PhuMob → Enormal → 高场饱和；其中 PhuMob 首先需要在当前合格 BGN 基线上校准单元载流子状态到迁移率、单元到边的形成方式，以及电子/空穴交叉 Jacobian。当前 `element_box` 对非 constant/Masetti 的保护保留。此次不能据此宣称已恢复原始完整模型或证明任意 Jacobian 列完整。

## 证据入口

- [192 项独立高精度账本](../../reference_tcad/simplemos_sentaurus2022/bgn_psi_derivative_20260908/precision_edges.csv)及[摘要](../../reference_tcad/simplemos_sentaurus2022/bgn_psi_derivative_20260908/precision_summary.json)。
- [最终实验合同](../../reference_tcad/simplemos_sentaurus2022/bgn_psi_derivative_20260908/state_independent/validation_contract.json)、[冻结身份](../../reference_tcad/simplemos_sentaurus2022/bgn_psi_derivative_20260908/state_independent/validation_freeze.json)、[完整回归比较](../../reference_tcad/simplemos_sentaurus2022/bgn_psi_derivative_20260908/state_independent/ctest_comparison.json)。
- [最终摘要](../../reference_tcad/simplemos_sentaurus2022/bgn_psi_derivative_20260908/state_independent/summary.json)、[全部 DC 尝试](../../reference_tcad/simplemos_sentaurus2022/bgn_psi_derivative_20260908/state_independent/attempts.csv)、[双初始化与电流对照](../../reference_tcad/simplemos_sentaurus2022/bgn_psi_derivative_20260908/state_independent/comparison.csv)、[Jv 分块账本](../../reference_tcad/simplemos_sentaurus2022/bgn_psi_derivative_20260908/state_independent/jvp_blocks.csv)。
- [资格边界明细](../../reference_tcad/simplemos_sentaurus2022/bgn_psi_derivative_20260908/state_independent/acceptance_detail.json)、[全局源闭合分量](../../reference_tcad/simplemos_sentaurus2022/bgn_psi_derivative_20260908/state_independent/global_closure_components.csv)、[最终源码与程序归档身份](../../reference_tcad/simplemos_sentaurus2022/bgn_psi_derivative_20260908/state_independent/accepted_snapshot_identity.json)、[最终证据冻结](../../reference_tcad/simplemos_sentaurus2022/bgn_psi_derivative_20260908/state_independent/final_evidence.json)。
- [原版本前像](../../reference_tcad/simplemos_sentaurus2022/bgn_psi_derivative_20260908/preimage.json)和[首个候选归档](../../reference_tcad/simplemos_sentaurus2022/bgn_psi_derivative_20260908/candidate_v1_identity.json)。
- [高精度审计脚本](../../scripts/audit_simplemos_bgn_psi_precision_20260908.py)、[验证脚本](../../scripts/validate_simplemos_bgn_psi_derivative_20260908.py)、[最终限定范围重放入口](../../scripts/validate_simplemos_bgn_psi_scoped_20260908.py)。

原始状态、日志和二进制位于忽略目录 `build-release/simplemos_bgn_psi_derivative_20260908/`；源码前像、首个候选和最终实验分别保存。历史报告和失败记录保持原文。本轮没有提交或推送 Git，生成仿真输出仍保留本地。
