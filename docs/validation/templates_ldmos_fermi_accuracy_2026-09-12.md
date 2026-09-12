# LDMOS 高精度 Fermi/BGN 联动实现

接续[差异定位报告](templates_ldmos_ialmob_performance_followup_2026-09-11.md)，
本轮按用户指定范围实现高精度Fermi函数、导数、逆函数、BGN联动及材料
常数审计，再复核G3/D5/D4。旧冻结曲线及合同不覆盖。

## 实现及精度

`fermiDiracHalf` 已替换原约0.4%误差的Bednarczyk近似：

- η<-8：六项收敛逸度级数，导数逐项解析计算。
- -8≤η<64：分段Chebyshev展开，Clenshaw计算值及同一展开的解析导数。
- η≥64：六阶Sommerfeld展开及其解析导数。
- 逆函数使用有区间保护的Newton，按相对密度误差停止；极小正密度直接
  使用非简并极限，避免旧绝对停止阈值损害稀薄载流子逆映射。

分段系数由定义积分经能量=t²变换、QUADPACK积分和DCT生成，没有器件
数据拟合。生成网格的最大相对值/导数误差约1.58e-14 / 1.73e-12。
正常C++构建只读系数表，不依赖Python或SciPy运行时。

Catch2使用独立Boost Gauss–Kronrod积分，覆盖[-40,100]、每个分段边界两侧、
从1e-300到1e200的密度逆映射、BGN及SI/TCAD换算。值/导数检查门限为
2e-12 / 2e-11，原器件验收门限不变。12个专项用例、855项断言通过。

BGN原本调用`inverseFermiDiracHalf`，因此直接使用升级后的同一函数，
没有另建低精度BGN分支。旧BGN单测固定的0.13920003585558116 eV来自
Bednarczyk近似；在原Nd/Nc、Na/Nv参数上用定义积分独立求逆得到
0.13961543141335747 eV。更新该参考值，原2e-12相对容差保留；不是放宽
容差来接受失败。修正后完整Release CTest **768/768通过，101.17 s**。
初次767/768的失败记录保留在`r1/release_ctest.log`。
随后将材料解析回归注册为`templates_ldmos_fermi_materials_regression`，
当前CTest共769项；新增项及原linked_d5回归定向2/2通过（1.48 s，内部共13个
Python用例），见`r1/registered_material_tests.log`。此后仅修改测试注册及
文档，没有重编译或替换曲线使用的冻结程序；不把这次定向检查写成769项全跑。

## 材料常数审计

原生来源为`ialmob_difference_location/native/bundle/Siliconc100.par`。
脚本直接读取Bandgap和DOS Formula 1，用当前C++ `PhysicalConstants.h`
中的kB、q、h、m0计算300 K参数。使用完整态密度表达式，不用参数文件
注释里的近似数值前因子2.540e19代替SI常数计算。

| 量 | 历史合同 | 本轮来源一致配置 |
|---|---:|---:|
| Eg / eV | 1.1241592307692307 | 1.12416 |
| 电子亲和势 / eV | 4.05 | 4.0727 |
| Nc / cm^-3 | 2.8566652281655894e19 | 2.856679069182081e19 |
| Nv / cm^-3 | 3.10464190834412e19 | 3.1046570334558958e19 |
| ni / cm^-3 | 10753127231.7754 | 10750016577.953705 |

ni严格由sqrt(Nc Nv) exp(-Eg/(2kT))获得。原生dEg0(OldSlotboom)=0，
活跃BGN路径另行施加窄化及Fermi修正，避免重复计入。新材料文件是独立的
`profiles/materials_fermi_consistent_300k.json`，历史`contracts/materials.json`
不变。材料审计附带输入哈希及每项差异，两个材料脚本回归通过。

## 固定原生状态检查的边界

Vg8/Vd40原生D4状态的10241个节点保留，密度统计只取输入密度>1e12 m^-3：

| 方案 | 电子误差中位数/最大 | 空穴误差中位数/最大 |
|---|---:|---:|
| 旧统计及材料 | 0.04537% / 0.39895% | 0.03595% / 0.21786% |
| 仅升级统计/BGN | 0.02994% / 0.10353% | 0.02997% / 0.07144% |
| 再加入来源一致材料 | 0.001021% / 0.10929% | 0.001037% / 0.07904% |

中位数改善不代表所有节点改善，也不意味着终端电流偏差同比下降。
重掺杂节点687的原生能带η约7.02265；在该原生η上直接计算定义积分，
密度与原生值相差约-0.21784%。该区域BGN逆函数差异会部分补偿密度差，
不能仅调整ni让某个节点恰好吻合。用户指南p.233方程43--50给出定义，
命令参考说明默认Fermi使用Joyce–Dixon；本轮实现定义积分精度，不声称
逐位复制其内部近似。

## 零压种子重新资格

第一次D5直接复用旧种子，在很小的绝对连续性残差下仍有逐行比接近1，
原逐行门限拒绝它。两次重闭合/密度恢复未把这些旧状态变成合格新种子，
失败证据保留在`r1/d5_vg{4,8}_full`，不评分为新曲线。

漏压为零且所有载流子接触均为零偏压时，以常量零QF表达热平衡，分别在
Vg4/Vg8求解Poisson（各3次更新），再用完整Newton原门限独立检查：
两者均0次额外更新通过，电流和连续性残差为零，电势块小于5e-8。
此处没有关闭逐行门限，也没有用Poisson收敛替代完整资格。

版本管理入口`prepare_templates_ldmos_fermi_equilibrium.py`重放得到的
两个合格种子与首次准备SHA256相同。Windows入口显式加入UCRT64 DLL路径；
首次脚本试运行的缺失DLL退出0xC0000135保留为环境失败，修正后成功。
新配置合同指向独立合格种子。种子准备工作与随后漏压曲线的更新数分开。

## G3及完整曲线复核

环境：UCRT64 Release，-O3 -DNDEBUG，Eigen SparseLU/COLAMD，gprof关闭。
原10241节点/19782三角形网格、AverageBox耦合、material_local Poisson，
300 K、Fermi、OldSlotboom、SRH/Auger，D4另开IALMob。密度保存为m^-3，
电流比较为A/um；物理参数变更不等同于数值门限变更。

冻结程序SHA256：`3e69f1a51b739db5c535dea78b4ef4b6260e5a2b1d529e67a4c7d2f78e0b34ee`。
证据根：`reference_staging/templates_ldmos_fermi_accuracy_20260912/`。

G3 Vd=0.1 V、Vg=0--5 V的31点及原六项门限通过：

| 指标 | 旧资格 | 本轮 |
|---|---:|---:|
| 电流log误差中位数 / dex | 0.00396406 | 0.00395654 |
| P95 / dex | 0.09969094 | 0.09967135 |
| 强反型端点相对误差 | 0.00804807 | 0.00804273 |
| Vth绝对误差 / V | 0.01075427 | 0.01075241 |
| 最大gm相对误差 | 0.01591068 | 0.01610083 |
| KCL相对误差 | 0.00123098 | 0.00103888 |

D5/D4从零开始的完整双栅压复核在`r2`执行，使用新种子，不续接旧模型
中间状态。两条同物理曲线并行仅用于资格检查，不构成独占性能基准。

D5完整62点及原工程/最终联合门限通过，零外层回退：

| 指标 | Vg4旧→新 | Vg8旧→新 |
|---|---:|---:|
| 电流相对误差中位数 / % | 1.430678→1.428997 | 1.330933→1.328570 |
| 最大电流相对误差 / % | 1.594271→1.592417 | 1.627051→1.625083 |
| 40 V电流相对误差 / % | 1.441878→1.440470 | 1.099059→1.095416 |
| Newton更新数 | 1139→1133 | 1116→1107 |

40 V栅压电流比误差0.337946%→0.340154%，略有增大，原门限通过。
新版最大归一化KCL为1.40750e-10%。因此局部密度精度提升并未消除约1.5%
的整体电流差异，不能把旧近似误差视为该差异的唯一原因。
来源为`r2/d5_joint/summary.json`和各曲线`fixed/ledger.json`。

D4完整62点及原工程/最终联合门限亦通过，零外层回退：

| 指标 | Vg4旧→新 | Vg8旧→新 |
|---|---:|---:|
| 电流相对误差中位数 / % | 1.857420→1.854387 | 1.691676→1.690462 |
| 最大电流相对误差 / % | 1.894521→1.891761 | 1.801843→1.801547 |
| 40 V电流相对误差 / % | 1.894521→1.891761 | 1.576488→1.574921 |
| 低压差分电阻误差 / % | 1.130280→1.129612 | 1.009554→1.011107 |
| Newton更新数 | 1133→1116 | 1073→1070 |

40 V栅压电流比误差0.312120%→0.310958%；最大归一化KCL为6.47656e-12%。
Vg8低压差分电阻误差略增，仍通过原门限。旧Vg8更新数采用干净完整重跑的
1073次，不混入早期错误参考系试验的19次失败更新。两条新版曲线各253个
子进程，包含零压检查与参考系检查；Vg4/Vg8偏移分别为28 V/4 V，原等价性
检查通过。来源为`r2/d4_joint/summary.json`和各曲线`fixed/ledger.json`。

累计子进程CPU时间在D4为1497.406→1647.063 s、1437.813→1604.547 s；
墙钟为1876.355→2876.875 s、1814.135→2757.407 s。上述为两次运行的原始
观察值，本轮并行及系统竞争条件未控制，不把更新数减少解释为耗时优化，
也不将全部计时差归因于Fermi实现。高精度实现的独占性能比较需单独开展。

本轮G3 31点、D5 62点、D4 62点共155项精确参考点复核全部完成，旧器件
门限未放宽。比较索引为`r2/comparison.json`，包含前后评分及ledger哈希。
高精度统计和材料一致性提升了局部密度精度，但完整曲线仍有约1--2%的
电流差异；剩余局部输运/离散差异及原生统计近似影响仍须分别定位。
收尾核对通过：冻结exe及133个当前C++源码/头文件哈希一致，当前构建为
Release（-O3 -DNDEBUG），四条完整曲线状态均completed、无回退；相关
三份文档24个本地链接和`git diff --check`通过。

## 重现入口

从worktree根运行，仿真及种子输出必须使用新的独立目录；系数和材料
生成器写入指定目标文件：

```powershell
D:/msys64/ucrt64/bin/python.exe scripts/generate_fermi_half_coefficients.py
D:/msys64/ucrt64/bin/python.exe scripts/prepare_templates_ldmos_fermi_materials.py --source-par reference_staging/templates_ldmos_followup_20260911/ialmob_difference_location/native/bundle/Siliconc100.par --output reference_tcad/templates_ldmos_sentaurus2022/profiles/materials_fermi_consistent_300k.json
D:/msys64/ucrt64/bin/python.exe scripts/prepare_templates_ldmos_fermi_equilibrium.py --manifest reference_staging/templates_ldmos_fermi_accuracy_20260912/r1/binary/manifest.json --output reference_staging/new_equilibrium_check
D:/msys64/ucrt64/bin/python.exe scripts/run_templates_ldmos_linked_d5.py --physics-profile D5 --bundle reference_tcad/templates_ldmos_sentaurus2022/profiles/linked_d5_fermi_accurate_inputs.json --manifest reference_staging/templates_ldmos_fermi_accuracy_20260912/r1/binary/manifest.json --output reference_staging/new_d5_vg8_check --gate 8 --points 31
```

D4改用`--physics-profile D4`及`linked_d4_fermi_accurate_inputs.json`，
Vg4改用`--gate 4`。生成系数不是每次运行的必要步骤。重建程序应使用
新的冻结清单，不能把新exe放进旧哈希清单冒充同一资格版本。
独立种子复核不会自动改写输入清单。上述曲线命令仍使用本轮归档种子；
若采用新目录下的种子，应复制生成新的inputs清单，并同步更新`seeds`
路径及`files`哈希，保留原清单和证据。
