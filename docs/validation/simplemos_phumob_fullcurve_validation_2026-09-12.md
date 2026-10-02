# SimpleMOS PhuMob 四条完整曲线验证

日期：2026-09-12。**204/204 个原生目标、408/408 个 Vela 状态及 204/204 对双初始化全部通过。** 408 次首次求解无失败、无重载。前置 [生产迁移与十六点复验](simplemos_split_state_production_validation_2026-09-12.md)的 32/32 状态、16/16 双初始化及 336 项正式 Jv 检查也已完成。

## 范围与参数

n19/n23 × Vd=0.05/1 V × Vg=0–1 V、步长 0.02 V，共四条 51 点曲线。Vela 分别从合格 0.8 V 状态延续，以及对应原生目标的电势/准费米势重建的相容初值启动。每个失败目标最多允许同偏置重载一次；本轮没有使用。延续路径从 0.8 V 向下到零，并独立向上到 1 V；这不是原始脚本从平衡态启动的完整回归。

物理参数与生产十六点一致：300 K、Boltzmann、OldSlotboom、plain PhuMob 和掺杂 SRH，Si 基础 ni=1.4638805412559193e10 cm^-3。Enormal 和高场饱和仍未启用。独立几何选项保持显式配置：Delaunay box、cell_material 介电系数、硅侧 box 输运/迁移率、三项 Poisson signed_transport 电荷体积；SRH 体积不变。

Vela 使用 UCRT64 Release/Eigen SparseLU/COLAMD、四次线性修正、显式 split 状态路径，默认仍关闭。状态相关残差和端口共享 100 位十进制求值；Jacobian 保留已校准的解析 double 实现，本轮不声称每个曲线点或任意弱交叉源的 Jacobian 都已单独审计。网格坐标 um，状态浓度 m^-3，电流 A/um，宽度 1 um。

Sentaurus 为 T-2022.03-SP2，保留 ExtendedPrecision(128)、Method=Super、Digits=12、ErrRef=1e-2、RhsMin=1e-20、Iterations=40 和 ExitOnFailure。只省去了重复的 Tcl 几何诊断输出，所有目标状态与物理场均已保存并取回。原生完整曲线与前置 16 个控制点的 Id 逐值一致。

## 完整电流对比

误差定义为 `100*(Id_Vela/Id_Sentaurus-1)`。表中使用延续路径；独立初始化路径与其最大相对 Id 差为 7.979680250e-10。

| 工况 | 合格点数 | 最小误差/% | 最大误差/% | 最大绝对误差所在 Vg/V |
| --- | ---: | ---: | ---: | ---: |
| n19，Vd=0.05 V | 51/51 | -0.031456709 | +0.004351757 | 0 |
| n19，Vd=1 V | 51/51 | +0.001221271 | +0.005319955 | 0.16 |
| n23，Vd=0.05 V | 51/51 | -0.018196719 | +0.005725237 | 0 |
| n23，Vd=1 V | 51/51 | -0.808135237 | +0.006529639 | 0 |

因此，这批完整 PhuMob 曲线没有在控制点之间出现更大的高 NWell 电流偏差。n23 高 Vd 深关断仍是最大差异点，约 -0.808135%。本轮没有新增电流百分比门槛，也不把当前 PhuMob 验收扩展为 Enormal/HFS 或原始全部工况已经通过。

![Id-Vg 完整曲线](../../reference_tcad/simplemos_sentaurus2022/phumob_curves_20260912/idvg.png)

![相对电流误差](../../reference_tcad/simplemos_sentaurus2022/phumob_curves_20260912/relative_error.png)

## 原门槛下的数值资格

| 项目 | 最坏值 | 原门槛 |
| --- | ---: | ---: |
| 全部自由载流子行比值 | 3.582585128e-7 | 1e-6 |
| 端口残差积分/提取差 | 4.440892099e-16 | 1e-8 |
| 双初始化 Id 相对差 | 7.979680250e-10 | 1e-6 |
| 双初始化电势最大差 | 1.523753e-11 V | 1e-6 V |
| 双初始化电子准费米势最大差 | 5.6771846e-11 V | 1e-6 V |
| 双初始化空穴准费米势最大差 | 2.96933173e-9 V | 1e-6 V |
| 双初始化密度最大相对差 | 1.148588873e-7 | 1e-4 |

每态均按原规则检查全部 1814 个自由载流子行、独立全局闭合、KCL 和端口；全局源检查仍使用原 1e-10 源下限，不能从带下限的通过推断任意相对 SRH 源精度或局部源扰动已校准。以上场差是两个 Vela 初始化结果之间的差，不是 Vela 对 Sentaurus 的绝对物理场差。

## 执行与失败记录

用户已授权在 `/tmp/vela_simplemos_phumob_curves_20260912` 上传和运行。输入包 SHA-256 为 `e5b79af14ce180f889dde305f2c5f6c2fc43b77c75249484c5b93d64668a2d73`；下载包为 `6345304eec45ac5a2b25e7d4c10a79f9b302c3686d795e721091e486b8040d81`。下载哈希、输入逐字节身份、版本、退出码、目标电压和原生 KCL 均核对通过。全部 204 个物理场已导出并冻结身份。

四条仿真以 0 退出后，控制脚本在旧 Bash 的最终空数组 `pids[@]` 遍历处触发 `set -u` 错误，未执行打包。原 run.sh 和 launcher.log 均保留；另用只归档脚本验证四个退出码、正常结束标记及每工况 51 个状态后补做归档，没有重跑仿真或改动输入。该控制流程缺陷与仿真数值失败分开记录。

本轮完整曲线运行期间没有修改求解器、材料参数、模型或门槛；使用的生产实现回归见前置报告的七项测试、267 个用例及 6875 个断言。没有把本轮称为全库 CTest 或全部物理模型回归。

## 后续与证据

按既定顺序，下一阶段进入 Enormal：先核对原生场、距离及有效迁移率形成方式，并校准零附加观察器，再实现/校准一致的状态及导数路径；其通过后才恢复高场饱和和原始 n17–n24、两个 Vd、0–2.5 V 全部工况。前置资料与未运行的观察器草稿见 [Enormal 前置核对](simplemos_enormal_restore_preflight_2026-09-12.md)。

- [原生合同](../../reference_tcad/simplemos_sentaurus2022/phumob_curves_20260912/native_contract.json)
- [Vela 合同](../../reference_tcad/simplemos_sentaurus2022/phumob_curves_20260912/vela_contract.json)
- [逐点完整比较](../../reference_tcad/simplemos_sentaurus2022/phumob_curves_20260912/comparison.csv)
- [完整比较证据](../../reference_tcad/simplemos_sentaurus2022/phumob_curves_20260912/comparison_evidence.json)
- [数值资格统计](../../reference_tcad/simplemos_sentaurus2022/phumob_curves_20260912/qualification_summary.json)
- [原生控制点重叠检查](../../reference_tcad/simplemos_sentaurus2022/phumob_curves_20260912/native_control_overlap.csv)
- [归档恢复证据](../../reference_tcad/simplemos_sentaurus2022/phumob_curves_20260912/archive_recovery_evidence.json)
- [执行脚本](../../scripts/simplemos_phumob_curves_20260912.py)

本轮尚未提交或推送。先前各阶段失败及本轮中间预览均保留，最终完整比较替代预览结论。
