# BVmethods NonLocalPath N3：Si 声子辅助内核

日期：2026-10-16。状态：私有 N3 内核及独立数值门完成，本机／T470p 同源 UCRT64 Release 验证通过。
此资格限于均质 Si、等温、分段线性带边、唯一内部 Keldysh 交点和显式横向窗口输入。
实际电子／空穴源装配、原生窗口取样、完整 RHS 和 Jacobian 尚未取得资格；下一步 N4。
公开前端继续拒绝 NonLocalPath（E106），没有修改局部 E2、生产模型默认或发行包能力。

后续[N4非均匀五态对照](bvmethods_nonlocalpath_n4_nonuniform_kernel_2026-10-16.md)
识别出原生带内位置插值与节点系数取样顺序；变带隙时不等于本N3的势垒能量P1数学表示。
本N3资格保持原范围，原生离散匹配候选留在N4诊断，不据此宣告生产模型已对齐。

## 实现及公式核对

新增 [私有接口](../../include/vela/physics/NonlocalPathPhononKernel.h)、
[实现](../../src/physics/NonlocalPathPhononKernel.cpp) 和
[独立数值测试](../../tests/test_nonlocalpath_phonon.cpp)，使用 `vela::experimental`。
依据本地 T-2022.03 手册 pp.531–532 式(508)–(511)，沿用
[N0 参数与预算](bvmethods_nonlocalpath_n0_2026-10-15.md)及
[N1 SI／版本契约](bvmethods_nonlocalpath_n1_2026-10-16.md)。商业手册缓存和原始输出仅保留在 ignored 目录。

- 起点能量为 `E=Ev(0)−Δv`，终点满足 `Ec(l)+Δc=E`。本轮没有自行加入 ±声子能量的端点位移。
  两段虚波数分别由价带／导带有效势垒及对应质量计算，交点满足 `κv=κc`。
- 实现两段 WKB 作用量、两段倒波数积分、沿程前因子积分、Bose 占据和有限横向动量窗口。
  波数使用 `ħ=h/(2π)`，前因子使用完整 Planck 常数 `h`；二者混用会破坏均匀场归一化。
- 沿用 N1 的 A/B 到微观参数转换和 `Rpath=0` 自动质量；新增只读常数访问器，
  不更改全局常数。2022 与 CODATA2018 两种显式常数模式均有归一化控制。
  A/B 的参考带隙与 300K 校准不随运行温度或当前窄化带隙重新拟合。
- `G=−Rnet`：Vela 正值表示生成，负值表示复合。对等温两端的 Fermi 占据差使用稳定的
  对数表达，保留极小非零 Fermi 分离；平衡时精确返回零。
- SI 输入为米、焦耳、千克、开尔文；输出率为 `m^-3 s^-1`，倒波数积分为 `m²`。
  起始价带梯度使用 `J/m`，不再重复乘电荷。

`tunnelingProfileFromGeometry` 从 N2 实际分段重心支撑重建带边与温度，要求路径完整、状态版本匹配、
材料一致及连续支撑。N2 结果新增起始价带梯度幅值和审计字段，几何寻路算法未改。
它拒绝歧义角点与陈旧 trial，不把失败路径转换成零源。

## 稳定积分与窗口边界

P1 带边上的平方根、倒平方根及有效带隙 `Eg,tun^-5/4` 使用稳定的分段解析积分，
在交点处分段，转向端点按其能量契约归零；没有端点 epsilon 截断或奇异点采样。
相邻带边／转向端点的 128 个 binary64 epsilon 判定仅处理序列化舍入，不是放宽 N0 物理预算。
温度只允许 4 个序列化 ULP 的差异，不将非等温模型近似为等温。

运行时用 `long double` 中间量和 `log1p/expm1`；没有新增运行时高精度库依赖。
`estimatedIntegrationRoundoff` 是算术舍入估计，不是严格误差上界，不能替代独立高精度验收。
资源预算及不可达精度、非有限输入、无效窗口、允许区口袋、多重／退化交点均明确拒绝，
不回退到 E2。多重分支或非等温扩展需要另行定义和验证。

调用方必须显式提供两端**带外可达能量**：`εmax−E`、`E−εmin`。
它们不是隧穿势垒段自身的带边极值；直接取转向端点会错误地把窗口置零。
手册公式给出窗口因子，但当前没有足够证据确定原生空间取样和截断身份，本轮没有猜测这套映射。
零窗口返回物理零；非零窗口再小也使用稳定表达。原生窗口和实际路径积分权重继续归 N4。

作用量、完整前因子及源值的对数保留在结果中。真实下溢时输出率为零，
`underflow=true`、符号非零且对数仍有限；平衡／零窗口则为符号零、对数负无穷，二者可区分。
超出可表示范围的结果报错，不饱和截平。

## 独立控制与结果

[测试](../../tests/test_nonlocalpath_phonon.cpp)的 oracle 使用 Boost 十进制 50 位计算、
独立 Gauss–Kronrod 31 点求积（目标精度 `1e-27`）和 180 次二分求交点。
平方坐标变换消除转向端点的积分奇异性；不复用生产解析原函数、求根或参数转换。
oracle 显式使用 N0 冻结的微观质量／耦合和手册常数。
表中误差是在输出 double 层比较，显示零不代表无限精度相等。

| 非均匀 P1 控制：价带窗口 eV（导带为其0.7倍） | 作用量／两段倒波数积分相对误差 | 源值相对误差（对数差反算） |
| --- | --- | ---: |
| 100 | 均在输出 double 层为0 | −7.105427357600977e-15 |
| 0.03 | 均在输出 double 层为0 | 0 |
| 1e-12 | 均在输出 double 层为0 | −1.776356839400249e-15 |

两机上述记录一致，最大源误差约 `7.11e-15`，通过冻结的独立内核 `1e-8` 门限。
另有以下非零控制，均通过：

- 均匀场 `1e7/3e7/1e8/5e8 V/m`、质量比 `0/0.5/2` 恢复
  `A(F/100 V/m)^2.5 exp(−B/F)`；作用量满足 `2S=B/F`。
- 平衡零源、交换两端 Fermi 后生成／复合反号、`1e-30 kT` 分离仍为非零；
  两种有限窗口各自独立改变或归零。
- 300→400K 和 Ppath `0.0185/0.074 eV` 响应，300K A 校准保持；
  400K 比值分别与 N0 原生控制 `1.0459893033702452/0.9007866902069033` 一致于 `1e-8`。
- A 倍增、B 增大、运行带隙改变，能量规范平移及 2/8/32/128 段划分控制。
- 次正规数源与 50 位 oracle 的独立整数 ULP 舍入完全一致；
  极低场下源下溢，保留对数并与均匀场闭式结果比较。
  Boost 十进制直接转 double 会把该次正规数清零，测试改用最小次正规量子为单位舍入，未放宽门限。
- 无效／非等温／多重分支、预算耗尽、N2 未决角点及陈旧状态明确拒绝。

## 构建与回归证据

两机验证源码清单2129文件，SHA-256：
`bb2626b120389cc2f574034b65be2fb3dd771c2f258f6336e64c399c30876726`。
远端复用隔离构建缓存，逐文件校验后刷新变更时间戳，构建日志确认重新编译。
远端独立作业退出0；源路径、PID、各阶段退出码与原始日志保存在 ignored
`build-release/bvmethods_nonlocalpath_n3_20261016/`，机器身份不写入公共文档。

| 检查 | 本机 Release | T470p Release |
| --- | --- | --- |
| N3 内核 | 11用例／108断言通过 | 11用例／108断言通过 |
| N2 几何回归 | 14用例／8274断言通过 | 14用例／8274断言通过 |
| N1 接口回归 | 11用例／97断言通过 | 11用例／97断言通过 |
| E2 数值回归 | 3用例／13断言通过 | 3用例／13断言通过 |
| 选定 CTest | 8/8通过 | 8/8通过 |
| 五组 Python 回归（含于上述CTest） | 81用例通过 | 80通过，1项可选jsonschema检查跳过 |

CTest 为 `nonlocalpath_n3_phonon`、`nonlocalpath_n2_geometry`、`nonlocalpath_n1_contract`、
`bvmethods_nonlocalpath_contract`、`bvmethods_original_reference`、`sdevice_cmd_parser`、
`sdevice_contracts`、`wp14_zero_qs`。未运行全量 CTest：本轮未改变生产装配或器件求解流程。
两机检测 HDF5／HighFive、SPQR／UMFPACK、STRUMPACK、SuperLU_MT、MUMPS 和 METIS；
本轮内核检查没有使用这些线性后端求解器件，不能据此认领后端器件资格。

```powershell
cmake --preset windows-ucrt64-release
cmake --build --preset windows-ucrt64-release --target test_nonlocalpath_phonon test_nonlocalpath_geometry test_nonlocalpath_contract test_recombination --parallel 2
ctest --preset windows-ucrt64-release --output-on-failure -R '^(nonlocalpath_n3_phonon|nonlocalpath_n2_geometry|nonlocalpath_n1_contract|bvmethods_nonlocalpath_contract|bvmethods_original_reference|sdevice_cmd_parser|sdevice_contracts|wp14_zero_qs)$'
./build-release/test_recombination.exe '*E2*'
```

## N4 入口与保留限制

N3 提供自主数学内核及 N2 状态适配，不提供节点实际源或逐路径原生率资格。
N4 应先核对带外窗口、起点积分权重、两端占据身份及电子／空穴远端支撑，再接实际残差和接触账本。
N2 节点1287／2567的接触与反射边界角点、截断／不完整路径和非唯一交点均须显式解析或拒绝，
不得默认为零并宣称原版通过。
本轮没有重跑 Sentaurus、没有更新 BV 参考或 N0 预算；完整 RHS、导数和原版从零纯 QS 仍待 N4–N7。
