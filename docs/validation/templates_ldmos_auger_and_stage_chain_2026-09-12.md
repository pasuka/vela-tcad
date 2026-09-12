# LDMOS Auger 密度依赖及后续阶段执行

用户已指定顺序：D4局部差异与Auger H/N0 → D3 → D2 → D1 → D0。
本轮已完成Auger实现、G3/D4/D5完整复核和D3/D2/D1原脚本等效评分；
D0完成独立热传导阶段，完整电热自洽尚未接入，不记作D0曲线通过。

## 修正阶段范围

重新读取原生`prepared/D0-original/IdVd.cmd`和D1--D4脚本发现，
eQuantumPotential(density)/hQuantumPotential(density)仅在Physics声明，
Solve没有对应量子势方程。此前[原生消融报告](templates_ldmos_stage4_idvd_execution_2026-09-03.md)
的D2→D3、D3→D4电流变化仅约1e-12%，且D1→D2的接触分量审计也无可分辨变化。
这些历史观察需核对原始文件及当前参考，不能把该报告早期D5阻塞状态当现状。

因此本轮先重核冻结原生消融证据，再使用最新合格经典解逐级评分；不擅自
往原始Solve增加活跃量子方程。等效曲线验收不代表通用量子或有限速度接触
边界实现通过，也不能自动推定自热条件下接触复合仍不可分辨。

## Auger实现范围

显式新增`solver.auger_density_dependence`，默认关闭。启用时使用
Cn(n)=Cn0[1+Hn exp(-n/N0n)]，Cp(p)同理；Cn0/Cp0仍为原300 K基础系数。
同一接口覆盖复合率、完整密度导数、Newton、Gummel恢复及诊断路径。
原生Siliconc100.par给出Hn/Hp=3.46667/8.25688，N0n/N0p=1e18 cm^-3。
输入的浓度换算遵守既有unit_scaling模式；不改变原残差/逐行/KCL门限。

新独立D4/D5配置命名为`linked_d*_auger_density_config.json`及对应inputs。
证据根为`reference_staging/templates_ldmos_auger_followup_20260912/`；
旧高精度Fermi资格及旧清单保留。Release编译通过，Auger模型3项测试61个断言、
IALMob耦合源Jacobian测试18个断言通过；完整Release CTest 773/773通过（148.92 s）。
随后完成G3及D4/D5器件复核，共155个独立求解的精确参考点通过原门限。

G3新增H/N0后31点、原六门限通过。阈值误差0.0107524 V，最大跨导
相对误差1.61008%，强反型端点误差0.804273%；与高精度Fermi基线一致。

## D4局部差异定位

固定Vg=8 V、Vd=40 V，使用新Release、来源一致300 K材料、原生网格
（10241节点）及原生保存状态。方程项范数为内部单位，端电流换算为A/um；
仅统计5674个有输运通量的节点，排除氧化层准费米势约束占位行。

- H/N0开关新增电子复合源L2=2.23990e-19，约为电子通量绝对量L2的
  1.41964e-21；固定状态端电流完全相同，不支持其为剩余电流差主因。
- 原生状态在Vela中的电子残差L2=0.213344，其中97.7018%的残差平方和
  位于距Si/氧化层界面10 nm以内；99.9967%位于100 nm以内。
- 将原生密度重新表示为Vela统计所需的准费米势后，密度映射最大误差
  降至4.31e-13/3.31e-13，但电子/空穴残差L2增至46.7604/1.19141。
  该操作改变准费米势梯度且可能破坏接触约束，不能作为物理修正或合格初值。
  此前冻结端电流回放因此返回不合格，失败证据保留。

这些结果将后续定位范围缩小到界面输运及统计/通量一致性，尚未宣称所有
局部差异闭合。证据位于`local_r1`、`local_r2_density_mapping`、
`local_r3_density_mapping`。不通过拟合材料常数或放宽原门限消除差异。

## 原生阶段链复核

重新解析原始`.plt`，校验准备脚本、差分文件及实际运行脚本哈希一致。
D3→D4、D2→D3、D1→D2各自最大非零电流变化仍小于4.5e-12%；
源漏分电流最大绝对差1.50162e-17 A/um，未出现可分辨的符号变化。
原始证据重核见`native_chain_r1/provenance.json`及`summary.json`。
新版D4完整通过后，已按D3→D2→D1顺序完成等效曲线评分，三项原工程/最终
联合门限均通过，结果位于`stage_equivalence`。这复用同一组D4的62个合格解，
不额外计为186个新求解点，不代表通用量子或有限速度接触实现。

## 完整等温曲线结果

统一使用UCRT64 Release、Eigen SparseLU/COLAMD，300 K，原精确网格及
AverageBox输运几何；外部坐标um、浓度cm^-3、电流评分A/um。冻结程序SHA256：
`71b79234aa80ec74f647b74e7b7b52a9dfe1f2ef4f9733911d91102de08e87e8`。

| 曲线 | 精确点 | 最大电流误差 | Newton旧→新 | 回退 |
|---|---:|---:|---:|---:|
| D4 Vg4 | 31 | 1.8917615% | 1116→1119 | 0 |
| D4 Vg8 | 31 | 1.8015470% | 1070→1073 | 0 |
| D5 Vg4 | 31 | 1.5924172% | 1133→1126 | 0 |
| D5 Vg8 | 31 | 1.6250834% | 1107→1096 | 0 |

旧值对应已资格高精度Fermi基线，只有新H/N0源项开关不同。
四条曲线相对旧基线最大电流变化≤1.832e-13（相对量），均为数值分辨量级；
Newton更新略有变化，不解释为已获得性能收益。运行期间有测试和本机虚拟机
并行负载，不作独占墙钟性能比较。详细计时、状态清单及哈希见
`completion_summary.json`、`r1/d4_joint/summary.json`、`r1/d5_joint/summary.json`。

保持独立输入合同，例：从worktree根运行以下命令；输出目录必须不存在。

```powershell
D:\msys64\ucrt64\bin\python.exe -X utf8 scripts/run_templates_ldmos_linked_d5.py --physics-profile D4 --bundle reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d4_auger_density_inputs.json --manifest reference_staging/templates_ldmos_auger_followup_20260912/r1/binary/manifest.json --output reference_staging/unique_d4_auger_vg4 --gate 4 --points 31
```

最终Release全部目标构建通过。增加热传导后全量CTest为780项：779项首次通过，
唯一失败是新增热合同进入原三文件黄金合同目录；移至独立`thermal/`目录后，
phase01和热评分两项回归2/2通过（4.21 s），未修改原黄金合同测试或物理门限。
全量运行日志为`final_release_ctest.log`，修复复核为`contract_relocation_ctest.log`。
另有配置/热评分17个Python用例通过，6项热传导Catch2用例通过。
文档38个本地链接检查通过；按仓库正常换行设置执行`git diff --check`通过。

## D0已明确的范围及已确认热门限

本轮用户指定D3→D2→D1→D0顺序，采用原D0作为最终目标。
热边界、材料参数及单位沿用[输入审计](templates_ldmos_self_heating_preflight_2026-09-12.md)：
th_lat=300 K，表面热阻5e-7 m² K/W；Si/Oxide温度相关热参数使用原生来源。
原D0仅在Solve加入Temperature，未启用Thermodynamic或Peltier。
按T-2022.03用户指南p.248方程71逐项映射热源，不用单独J·E代替全部热源。

用户于2026-09-12明确回复“采用建议热门限”，以下标准作为正式D0评分合同：

| 指标 | 已确认门限 |
|---|---|
| 全域稳态热平衡 | 发热与边界净散热的绝对差 / 二者绝对值较大者 ≤0.1%；零源另检绝对残差 |
| 峰值温升误差 | ≤max(1 K, 原生峰值温升的5%) |
| 温升场误差 | 同网格体积加权RMS温升差 ≤max(1 K, 原生温升场加权RMS的5%) |
| 电学指标 | 保留原工程/最终电流门限及数值块、逐行、KCL门限 |

热学零源恒温、解析热传导、Robin边界、界面热流守恒、SI/TCAD换算可以
独立开展，不需要先认领完整D0验收。温度相关统计、复合、迁移率及电热
Jacobian均须一致接入后才可把电流变化归于完整原D0模型。

热门限已保存为受版本管理的
[`d0_thermal_acceptance.json`](../../reference_tcad/templates_ldmos_sentaurus2022/thermal/d0_thermal_acceptance.json)。
独立热传导实现、原生热源对照及其资格边界见
[D0热传导报告](templates_ldmos_d0_thermal_operator_2026-09-12.md)。
