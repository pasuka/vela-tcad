# LDMOS 载流子质量矩阵与投影伪瞬态验证

日期：2026-09-15；分支 `codex/templates-ldmos-phase-a`。承接用户指定的
[质量矩阵验证任务](templates_ldmos_newton_update_execution_2026-09-15.md#5-明确停滞样本后的投影伪瞬态研究)。

结论：质量项的符号、单位、四变量偏导和受控扩散验证通过；三个伪时间尺度均未
恢复两个冻结难点。固定状态隔离证明，1× 质量项消除了初始非正密度目标，但方向
对当前稳态残差范数不再下降，分别冻结迁移率/体复合导数也未改变这一点。
本轮候选不晋级；R7 默认、物理模型和最终电学/热学/逐行门限不变。

## 实现与适用范围

连续性残差为电子 `sum(Jn)-q*A_cont*R`、空穴 `sum(Jp)+q*A_cont*R`。
存储量采用 `rho=(0,-q*A_cont*n,+q*A_cont*p,0)`，单位 C/m，
求解方向满足 `(J + d(rho)/dx / tau) delta = -F`。

- `A_cont` 使用 `recombination_area_m2`，缺省按现有几何合同回退至
  `silicon_area_m2`。不以 Poisson 面积替换显式连续性面积。
- 密度及其 ψ/fn/fp/T 偏导来自同一 `SiliconThermalPhysics`，包括 Fermi、BGN
  和温度链式项。ψ/T 为代数方程，没有新增晶格热容或电势虚假存储项。
- 非硅节点、Dirichlet 和理想接触约束行无质量项；有限空穴交换接触的自由
  空穴连续性行保留存储。当前任务只验证数值伪时间，不认领物理瞬态资格。
- `diagnostic_pseudo_transient` 缺省关闭；与 V1 投影组合时保留 0.01 相对正性下限、
  精确统计反算、ψ/T 原限幅和八次映射失败后 QF 回退。
- 初始 τ 为自由载流子行 `abs(M_ii/J_ii)` 的正有限值中位数乘尺度。
  三档尺度预先固定为 0.1、1、10；两点的 1× τ 分别为
  `7.414123e-14 s`、`7.781222e-14 s`。
- SER 以初始稳态 Jacobian 的行尺度固定衡量后续残差，τ 增长因子限制在
  [0.5,4]，τ 限制在初值的 [1e-6,1e12]。线搜索失败时 τ 除以 10，
  原迭代预算内最多重试三次，全部尝试计入成本。该规则没有放宽最终门限。
- 为隔离原因，增加只读 `diagnostic_pseudo_direction_audit`，只允许单次耦合
  伪瞬态迭代。它比较完整/冻结迁移率导数/冻结体 SRH+Auger 导数/同时冻结的
  方向，并分别加或不加质量项。冻结包含迁移率对边电流及一致热源的链式项；
  接触边界导数不冻结。所有物性值和真实残差保持逐值一致；冻结方向不用于更新。

## 数值验证

1. 独立中心差分覆盖 300/475 K、正负电势、高低掺杂及 40 V 参考平移，检查
   存储量的全部四列偏导；同时验证自由/约束行、面积线性比例和零面积。
2. 用实际 `thermalSgCurrent` 的稀薄等温极限构造两端固定密度的单自由节点扩散。
   对电子和空穴分别验证 `M/tau+J` 的更新等于后向 Euler 解、目标密度为正且
   随 τ 增大趋近两端均值。两组质量项测试共 805 断言通过。
3. 电热装配测试检查导数隔离不改变真实残差/热源，体复合隔离不修改 Poisson、
   温度和接触约束行，完整 Jacobian 仍通过原有限差分测试。
   电热/质量项定向检查共 7 个用例、4962 断言通过。
4. 生产电热回归共 26 项通过，包含全约束质量项无作用、自由载流子质量项、
   实验组合拒绝，以及只读方向审计与无审计更新轨迹逐值一致。
5. 最终完整 Release 构建成功，单线程 CTest **826/826 通过**，总测试时间
   219.37 s。日志保存为 `build-release/ptc_final_build.log` 和 `ptc_final_ctest.log`。

新增测试初次因 `growth_newton=2` 超过单次更新预算及误从 stdout 而非保存日志
读取异常而失败；修正测试调用后通过，没有改变求解器容差。

## 两个冻结难点

Windows MSYS2 UCRT64 GCC 16.2.0、Release `-O3 -DNDEBUG`、UMFPACK，
`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=1`。使用原 R7 同网格、SI 几何、模型、
参考系、输入和最终状态比较门限；Vg=4/8 V，Vd 均为 5.333333333333333 V。
清单中的索引分别为 `vg4_p11`、`vg8_p9`，不是相同推进编号。
点间和变体间串行；开发/构建存在同机负载，墙钟仅为诊断记录。

| 配置 | Vg4 更新/通过 | Vg8 更新/通过 | 总尝试 | 总试算 | 总装配 | 总分解 |
|---|---|---|---:|---:|---:|---:|
| R7 | 15/是 | 15/是 | 30 | 44 | 71 | 30 |
| V1 | 19/否 | 16/否 | 35 | 293 | 330 | 35 |
| PTC+V1，0.1× | 0/否 | 0/否 | 8 | 256 | 258 | 8 |
| PTC+V1，1× | 0/否 | 0/否 | 8 | 256 | 258 | 8 |
| PTC+V1，10× | 60/否 | 2/否 | 66 | 281 | 345 | 66 |
| PTC+QF，1× | 0/否 | 0/否 | 8 | 192 | 194 | 8 |

“0 更新”表示没有接受试步，是失败而非加速。0.1×/1× 均耗尽原预算内的三次
失败重试；10× 在 Vg4 达 60 次预算上限、Vg8 线搜索失败。V1 原失败轨迹仍复现。
Vg4 的 10× 固定尺度残差最终仅降至初始的 0.931834；Vg8 为 0.988609。
前者 Poisson 块已闭合，电子、空穴和逐行门限仍未通过，不能按全局小量提前接受。

R7 墙钟为 16.60/17.70 s；V1 为 47.53/37.05 s；PTC+V1 三档依次为
38.21/34.76、39.69/42.01、70.86/40.09 s。失败状态不能与成功 R7 构成性能资格。
质量装配本身约 0.02 s/次，包含局部物性准备；该嵌套时间不与密度准备重复相加。
本轮主要损失来自被拒试算和未闭合，不能通过只减少质量项准备来解决。

## 固定状态的方向隔离

定义 `s = (W*F)^T (W*J*delta) / ||W*F||²`，W 为初始固定行尺度，
J 是完整真实稳态 Jacobian；s<0 才是该范数的一阶下降方向。
完整 Newton 的 s≈-1，是这一定义的自检。16 个方向的线性方程相对残差最大
`5.87e-11`，不足以解释下面的正斜率。

| 导数处理 | 无质量项的非正目标 n/p：Vg4；Vg8 | 加 1× 质量项后 n/p | 加质量项后 s：Vg4；Vg8 |
|---|---|---|---|
| 完整 | 69/84；71/80 | 两档均 0/0 | +0.06518；+1.10504 |
| 冻结迁移率导数 | 67/80；61/75 | 两档均 0/0 | +0.05928；+1.12644 |
| 冻结体复合导数 | 69/67；72/54 | 两档均 0/0 | +0.08498；+1.14619 |
| 同时冻结 | 67/66；61/49 | 两档均 0/0 | +0.07901；+1.16811 |

非正目标按原始全方向 `n+delta_n<=0` 计数，不等同于 0.01 相对下限触发的
投影数量。完整质量方向使用 ψ/T 上限后，再缩至该步的 1e-4，实际投影候选
残差比仍为 1.000000905 / 1.000005329，与正斜率一致。

因此：质量项已发挥正性正则化作用，但这种方向不保证当前稳态 L2 单调下降。
去掉两类系数导数仍然上升，说明它们不是这些初始失败的充分解释。
这不证明物性导数无关，也不能把受控线性扩散的正性推广为完整耦合系统的收敛保证。

后续应先设计与伪瞬态方向匹配的接受准则或处理代数约束的一致初始化，再做
同样的固定状态斜率和难点验证；不建议仅扩大 τ 扫描或继续放宽终态门限。
本轮没有实现新的非单调接受准则或修改已批准验收口径。
用户随后授权的[实际伪瞬态缺陷接受验证](templates_ldmos_pseudo_acceptance_2026-09-15.md)
已完成：初始拒步被跨过，但新准则与两种时间控制均未闭合稳态，候选不晋级。

## 证据与复现

证据根：`reference_staging/templates_ldmos_extrapolation_20260915/`。

- `ptc_stalls_validated/`：12 组完整控制，含 EXE、16 个 DLL、源码归档、
  输入/输出、所有失败日志及汇总。运行前误填 `vg8_p11` 的目录 `ptc_stalls/`
  标记为 `aborted_before_simulation`，没有仿真结果。入口现先校验清单再建输出目录。
- `ptc_direction_audit_final/`：最终程序的两个固定状态、16 个导数方向；
  单次迭代只检查诊断不扰动更新，不认领完整点资格。
- `ptc_findings_final.json`：独立分析结果；分析入口校验 EXE/DLL/输入/输出哈希
  和矩阵完整性。早期 QF-only 记录中的非正目标计数占位 0 不作正性证据；
  最终程序在未测量时输出 null，分析脚本亦明确排除。

```text
python scripts/run_templates_ldmos_density_projection.py --manifest reference_staging/templates_ldmos_joint_20260914/production/r8_predictor_inputs/manifest.json --probe build-release/electrothermal_probe.exe --output <new-directory> --variants baseline v1 ptc_v1_s01 ptc_v1_s1 ptc_v1_s10 ptc_qf_s1 --cases vg4_p11 vg8_p9
python scripts/run_templates_ldmos_density_projection.py --manifest reference_staging/templates_ldmos_joint_20260914/production/r8_predictor_inputs/manifest.json --probe build-release/electrothermal_probe.exe --output <new-audit-directory> --variants ptc_v1_s1 --cases vg4_p11 vg8_p9 --direction-audit
python scripts/analyze_templates_ldmos_pseudo_transient.py --matrix <matrix-directory> --audit <audit-directory> --output <new-findings.json>
```

独立实验开关与源文件保留用于复查，不修改生产默认，不把旧完整曲线资格转移到候选。
