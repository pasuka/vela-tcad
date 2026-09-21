# LDMOS Vg=8：26.6667 V 换参考失败的分层定位

日期：2026-09-09，Asia/Shanghai。工作树 `templates-ldmos-phase-a`，源码基线
`959c1c9` 加现有未提交修改。本报告是分析证据，不是生产修复或全曲线通过声明。

## 结论

已实测定位到两处 QF 精度损失：外部平移参考值时没有补偿舍入低位，以及
Newton 重启重打包的 `(old_reference + increment) - new_reference` 运算顺序。
CSV 读取本身未改变电势、QF 参考或增量的 binary64 数值。修正这两处的诊断
对照显著减小初始载流子残差，并消除同偏压再次重启引入的局部违规，但仍未
通过电子全局残差门限。因此不能将两项对照宣称为完整修复。

剩余障碍表现为平移后的载流子密度重建／电荷装配误差和小更新尺度下的
舍入敏感性；这支持数值精度导致收敛尾段停滞的解释。尚未用一个通过原门限
的完整反事实修复严格证明其全部因果链，也未逐项隔离缩放、密度计算、SG
通量和 vector HFS 梯度重建各自对最后电子残差的贡献。

## 输入与诊断隔离

- 网格 10,241 节点、19,782 三角形、30,022 边；外部 AverageBox 输运、
  barycentric 体积、material-local Poisson、legacy node-local 接触。
- 300 K，Fermi–Dirac、OldSlotboom BGN、SRH/Auger；`constant_field +
  transport_cell_vector`，完整迁移率 Jacobian。predictor、IALMob、雪崩、量子和热耦合关闭。
- 物理 Vg=8 V、Vd=26.6666666666667 V；frame 0 与整体减 28 V 的 frame 28。
  后者 source/substrate=-28 V、gate=-20 V、drain=-1.3333333333333002 V。
- 输入几何 μm、密度 m⁻³、电势 V、电流 A/μm；unit scaling；UCRT64 GCC
  Release `-O3 -DNDEBUG`，Eigen SparseLU/COLAMD、L2 行列平衡。
- 原全局 psi/e/h 门限 `5e-8 / 1e-11 / 3e-10`；局部 `eps_row=1e-8`、违规零；
  KCL ratio ≤1e-8。诊断单点 solve 使用 40 次上限；所有失败均在线搜索阶段
  提前结束，没有因 40 次预算截断。失败状态没有加入曲线。
- 父状态 SHA256：`f861982783bbdc730d5eb742aa6406430a85488244d15f72574d8c64d6a2dec2`。
  正式 runner SHA256：`52af2e8a2a489f9be299ee45e1b0fdd5bd1ef5f4a61e43693b701ce22e428087`。

所有新增代码、配置、状态和探针输出位于独立的
[诊断目录](../../reference_staging/templates_ldmos_frame_diagnosis_20260909/)。
未编辑正式源码、原控制器、原曲线文件或冻结二进制。诊断编译只替换静态库中的
NewtonSolver 对象，并复用原 runner 对象和其余库；命令、输入哈希及两份局部
diff 分别保存在 `build_manifest.json`、`repack_build_manifest.json`、
`probe_source.diff` 和 `repack_source.diff`。分析脚本核验了 Vg=4 计划内所有冻结文件未变。

## 实际失败与昨天的区别

昨天使用 `edge_projection`，今天使用 `transport_cell_vector` 和补全后的
Jacobian。昨天 Vg=8 的换参考一步通过，完整 31 点数值收敛，但电流精度不通过。
今天 Vg=4 的相同换参考也已通过。不能把昨天的完成记录当作今天新算法的验证。

| 今天 Vg=8 的实际阶段 | psi 残差 | 电子残差 | 空穴残差 | 局部违规 |
| --- | ---: | ---: | ---: | ---: |
| 平移前保存的合格点 | 4.97694e-9 | 9.83988e-12 | 1.53826e-24 | 0 |
| frame baseline，9 次更新后失败 | 4.13662e-8 | 9.64519e-11 | 2.32296e-21 | 0 |
| frame reclose，8 次更新后失败 | 4.13662e-8 | 5.41988e-11 | 4.04722e-14 | 789 |

失败类型为 `line_search_non_decrease`。电子全局块未通过，所以不满足
`density_eligible` 的前提；`change_frame()` 记录 `failed_closure` 并终止。
它尚未进入后续电流／状态 `failed_equivalence` 检查，亦非用尽 160 次正式预算。

## 1. 验证探针对应实际运行路径

现有 `evaluateResidual` / `evaluateDirectionalDerivative` 探针只配置接触 QF
参考，不执行实际 `solve()` 的按初态重定位。未经校准的探针在原合格状态上
给出电子残差 `4.65012e-10`，不能据此否定原合格点。

独立诊断版本增加 saved-reference、initial-reference 和接触投影观测：

- contact 模式在两个参考系下的完整 residual CSV 与正式二进制逐字节相同。
- initial-reference + 接触投影模式的三个初始残差块与实际
  `child_00578_frame_baseline/newton_iterations.csv` 初始行一致。
- 读取后的 psi/phin/phip、两类 QF reference/increment 与输入逐节点浮点相等；
  密度列在两份平移输入中未变，读取时只执行既定单位转换。

因此后续分层比较不是把两个不同的探针坐标策略误当成同一状态。

## 2. 平移参考值丢失低位

原平移脚本对 reference 执行 `float(reference)-28`，increment 原样保留。
以输入 binary64 的精确值作为基准，用 Decimal 80 位计算误差：硅节点中电子
2,936 个、空穴 2,993 个出现差异，最大均约 `1.77636e-15 V`。

诊断补偿仅在新文件中把该舍入余量放回 increment：

```text
new_ref = binary64(old_ref - 28)
new_increment = binary64((exact(old_ref) - 28 - exact(new_ref)) + exact(old_increment))
```

密度和 psi 沿用原平移文件，保持其他条件不变。匹配实际初始化路径时：

| 初始状态 | 电子残差 | 空穴残差 |
| --- | ---: | ---: |
| 原平移 | 2.78871e-8 | 1.19607e-10 |
| 仅补偿参考低位 | 1.46171e-10 | 9.88507e-14 |

电子残差约降低 191 倍，证明这项舍入是实际扰动来源。但原二进制对补偿态
求解 7 次仍失败，电子残差 `4.23824e-11`；再次重启又出现 829 个局部违规。
仅改变外部平移算法不够。

## 3. 重启重打包又吞掉小增量

`packReferencedSolution()` 和 `solve()` 的四处重打包均使用：

```text
(long double(old_ref) + long double(increment)) - long double(new_ref)
```

即使使用 long double，先与约 −28 V 相加也会吞掉足够小的 increment。
例如节点 0 的电子 increment `-1.3114e-29 V` 和空穴 increment
`4.9272e-21 V` 在原平移态被打包成 0。

诊断副本只改为 `(old_ref - new_ref) + increment`，结果：

- 原参考系同一文件的电子残差从重打包后的 `3.02628e-11` 恢复到
  `9.83920e-12`，与保存时 `9.83988e-12` 接近且同在原门限内。
- 与平移低位补偿组合，初始空穴残差进一步从 `9.88507e-14` 降至
  `1.62046e-23`。
- 组合对照的再次重启局部违规为 0，最大 ratio `3.94728e-12`；
  相比仅平移补偿再重启的 829 行、ratio `0.226726`，局部失真消除。
- 但组合对照 5 次更新后电子残差仍为 `7.18131e-11`，再次重启 0 次更新便
  遇到相同线搜索拒绝。局部和 KCL 通过不替代电子全局门限。

## 4. 剩余误差：电荷重建与小尺度装配响应

Poisson 分项探针显示，输入密度不变，但由平移后的电势/QF 重新计算密度时，
最大相对差异从约 `1.4e-14` 增至 `2.7e-13`；补偿 QF 后仍约 `2.13e-13`。
原参考系与补偿平移态的分项差异范数：

| 分项 | L2 差异 |
| --- | ---: |
| 介电通量 | 1.42712e-10 |
| 电子电荷 | 4.91177e-8 |
| 空穴电荷 | 8.20402e-9 |
| 掺杂电荷 | 0 |
| 最终 Poisson 残差 | 4.89040e-8 |

较大电荷项相消后的残差对这些微小变化很敏感。上述实验定位了剩余 Poisson
变化的主要分项，但不能仅凭分项差异宣称已经证明电子全局残差停滞的唯一根因。

Jacobian 共检查 48 个方向，使用 `||Jv-FD||/max(||Jv||,||FD||)` 分块重算，
分母不加 1 的下限：

- 前 30 个为两参考系下节点 2951、3984、4535 的电子 QF，以及节点 4535 的
  psi、4652 的空穴 QF，幅度 `1e-5/1e-6/1e-7 V`。主要误差随步长缩小下降；
  例如原参考系 phin3984→电子块为 `1.2942e-3 → 1.2946e-5 → 1.2930e-7`。
  未复现此前稳定的大幅迁移率导数遗漏，但这不是所有列的穷尽证明。
- 另 18 个以重打包顺序修正的诊断版本比较原态与补偿平移态，幅度
  `1e-9/1e-11/1e-13 V`。很小幅度时误差反而增大，符合舍入影响：
  phin3984→电子块在 `1e-13 V` 为原参考系 `1.4793e-4`、平移系 `7.6043e-3`；
  psi4535→Poisson 块分别为 `1.3323e-4`、`5.1644e-3`。
  这属于特定方向的小幅度差分敏感性，不能外推为整个 Jacobian 有相同比例错误。

## 后续建议与证据

优先把参考值平移和重打包的精度问题变成独立回归测试；生产修复应把低位
保持作为明确契约。随后针对未通过的电子块，继续拆分缩放电势重建、载流子
密度和 SG/vector-field 的残差累加，并核验严格线搜索的实际评价量。不能靠
放宽电子门限或跳过 frame 资格检查继续宣称 D5 通过。

本次仅完成定位和反事实对照，未替换正在运行的求解器，未续跑 Vg=8，未提交 Git。

- [机器可读分析和校验](../../reference_staging/templates_ldmos_frame_diagnosis_20260909/analysis.json)
- [诊断构建及分层探针脚本](../../reference_staging/templates_ldmos_frame_diagnosis_20260909/diagnose.py)
- [外部平移补偿](../../reference_staging/templates_ldmos_frame_diagnosis_20260909/compensated.py)
- [重打包运算顺序对照](../../reference_staging/templates_ldmos_frame_diagnosis_20260909/repack.py)
- [原二进制同偏压验证](../../reference_staging/templates_ldmos_frame_diagnosis_20260909/closures.py)
- [分析脚本](../../reference_staging/templates_ldmos_frame_diagnosis_20260909/analyze.py)
- [昨天的 Vg=8 完整结果](templates_ldmos_stage4_vg8_rerun_2026-09-08.md)
