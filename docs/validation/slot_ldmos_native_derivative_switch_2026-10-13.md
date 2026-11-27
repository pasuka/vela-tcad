# SLOT 原生导数开关与 Vela 省略项查证（2026-10-13）

## 结论与范围

已在 Sentaurus 同一冻结态上取得显式 `Derivatives`/`-Derivatives` 与
`AvalDerivatives`/`-AvalDerivatives` 的四组矩阵及实际求解 RHS。
**省略迁移率导数有直接原生实验依据；但原生开关比 Vela 当前省略范围更宽。**

- 四组第0次 ψ/n/p、实际求解 RHS1 逐位一致；载入后的171个 Plot 数组也逐位一致。
- 显式全导数有111条负空穴自身对角；关闭 `Derivatives` 后为0。仅关闭
  `AvalDerivatives` 仍为111，不能用雪崩导数解释这组负对角。
- 与 Vela 匹配态相比，111个负对角节点集合完全相同。经独立缩放和 Fermi 链式
  变换，被省略的对角项相对差最小/中位/最大为0.187/0.251/0.413 ppm。
- 空穴自身块的省略矩阵相对 Frobenius 差为8.5846e−7（0.858 ppm）。
  这是局部导数机制的强证据，不是完整 Jacobian 或完整 Newton 路径一致性证明。
- 原生 `-Derivatives` 还改变其他迁移率依赖以及雪崩导数中的一部分；不能用一个
  日志“保留雪崩导数”推断所有源项 Jacobian 元素不变。

本批不改变 Vela 代码、默认配置、物理或停止规则，不启动新的偏压续算；新增合格点0。
只比较各臂第0次装配。原生 `Coupled(Iterations=0)` 实际仍执行一次诊断更新，
该更新及末态单独保留，未被当作冻结态输入或新合格点。

## 输入、版本和实验定义

- 同一 SLOT 网格，13,881节点，原生矩阵39,149维：13,881 ψ + 2×12,634人口。
- 固定 Vg=5 V、Vd=0.005 V；300 K；Sentaurus T-2022.03-SP2；binary64；PARDISO单线程。
- 原生完整 PhuMob/Lombardi/HFS、Fermi定义积分及版本绑定BGN PMI、SRH/Auger/电离。
  日志确认电子和空穴 HFS 均以准费米势梯度为驱动力；未启用低密度驱动力插值。
- 源状态：`build/s10native64/native/double64/newton_1_0.tdr`，即原生首次更新后的状态。
  由保存模板写入该态ψ和QF；按接触名称核对source=0、drain=.005、gate=5、substrate=0。
  SLOT网格标记未被添加为新的外加电压接触。
- 四臂设备保存态SHA256相同：
  `881eba6bd93786d584d1baae5e747c1e8c4d23c305722e2d42bb82226d3c6b41`。
- 原生PMI SHA256：`7cf557d529fee88896100411e0fc71152cbeb97153322ca701e353b9a18cc525`。
  运行时包装器、输入及输出身份见证据清单；包装器哈希不能替代全部运行库身份。

| 臂 | Math开关 | 运行日志 | 墙钟时间 |
|---|---|---|---:|
| full | Derivatives AvalDerivatives | With avalanche derivatives | 78.48 s |
| reduced | -Derivatives AvalDerivatives | Without exact derivatives；With avalanche derivatives | 74.21 s |
| avalanche_off | Derivatives -AvalDerivatives | Without avalanche derivatives | 72.61 s |
| both_off | -Derivatives -AvalDerivatives | Without exact derivatives；Without avalanche derivatives | 69.17 s |

每臂180秒硬上限，四臂串行；13:49:12至13:54:15完成。均退出码0、stderr为空。
这些时间包含初始化、调试矩阵和场量输出，不代表正常偏压步计算成本。

## 同态与 RHS 的审计

四臂原生变量编号、变量值与RHS1完全相同，矩阵差可归因于导数开关。
原生调试输出 `RHS2` **不是本次求解使用的RHS1**：既有连续轨迹的RHS2(k)与
下一次RHS1(k+1)逐位一致，说明它对应更新后的残差。四臂方向不同，RHS2也不同；
不能据此否定同态RHS不变。

重新载入与原历史snapshot1并非位级相同：ψ最大差2.572e−14 V，n/p最大相对差
6.847e−14/6.665e−14；RHS1相对L2差2.825e−10。因此四臂内部的位级因果对照
最强，历史矩阵匹配仍带载入舍入边界。

Vela复用已保存的同态 `native_post_current`/`native_post_frozen` 矩阵，无新本地求解。
这对配置仅切换 `diagnostic_freeze_transport_qf_mobility`，物理残差及solve RHS逐位相同。
Vela导入态与原生新载入人口的差另存 `analysis.json`，近零QF派生量仍有表示敏感性。

## 坐标转换与支撑、符号、量级

比较对象为矩阵差，而非整个矩阵的任意拟合。原生坐标为ψ/Vtn、n/Nref、p/Nref；
Vela为内部缩放ψ、电子QF、空穴QF。使用 Nref=1.4824351e10 cm⁻³，
Vtn=0.025851995266484913 V，Vela电势尺度0.025851999786435535 V。
若 g=F′₁/₂/F₁/₂，则人口列转换为−Nref/(n g_n)、+Nref/(p g_p)；
固定人口改变ψ时，三个电势列同时按Vtn/Vvela变化。

载流子行比例由独立归一化确定：
`C0*mu0*Vvela/(Nref*1e4)=365666721.5415312`，
C0=1.4797799323493304e21 cm⁻³、mu0=1417 cm²/(V·s)。没有用矩阵元素反向拟合比例。
本报告表中矩阵量均为转换后的原生内部坐标，不是安培或伏特本身。

令 Dn=J(Derivatives,AvalDerivatives)−J(-Derivatives,AvalDerivatives)，
Dv为Vela完整输运Jacobian减去窄范围省略Jacobian后的链式变换。

| 行/列块 | ‖Dn−Dv‖/‖Dn‖ | 方向余弦 | 含义 |
|---|---:|---:|---|
| 电子/电子 | 2.0401% | 0.9997926 | 主项接近，原生还删除其他依赖 |
| 空穴/空穴 | 8.5846e−7 | 0.99999999999967 | 低密度空穴自身块高度吻合 |
| 电子/ψ | 12.6764% | 0.991933 | 仍有更宽的电势与源项链 |
| 空穴/ψ | 99.9999% | 0.003699 | 不能称两种开关全块等价 |

Vela窄范围省略项的电子/空穴交叉人口块为零；原生开关在这两个块均非零。
按各块最大元素的1e−10作支撑统计阈值，空穴自身块Dn/Dv分别53,223/48,794项，
共同48,711项，其中48,685项同号。原始未阈值化项、阈值定义与所有块结果完整保留；
不能把这一支撑统计当作精度资格门槛。

| 空穴节点 | 原生完整Jii | 原生省略后Jii | 原生省略项 | Vela省略项（已转换） |
|---|---:|---:|---:|---:|
| 6220 | 27.617265 | 28.221600 | −0.60433418 | −0.60433547 |
| 6132 | −48.977296 | 40.319688 | −89.29698425 | −89.29701873 |
| 7518 | 44.238017 | 44.238498 | −4.8186844e−4 | −4.8185800e−4 |
| 10849 | 1.1856086 | 1.1919010 | −0.0062924631 | −0.0062925332 |

6132原生完整系统的Δp/p=−11.5092；省略后为+118.9981。
7518自身对角两臂都为正，但完整系统方向由−28.0242变为+3757.5780，说明邻接和
跨块耦合仍重要，不能仅用自身对角解释每个节点。
6220在这个**原生首次更新后的5mV态**自身对角为正，不得混同此前Vela6mV失败态的负对角。

## 雪崩导数不是独立于迁移率开关的一整块

四臂差分显示明显交互：

- 保留Derivatives时，开关AvalDerivatives的矩阵差L2=2.851880e6；
- 关闭Derivatives时，同一AvalDerivatives开关差L2仅0.00961676；
- 二阶交互差的电子/空穴对应行几乎相同，成对行差相对L2为6.46e−11，
  支持其属于成对产生导数，而非纯守恒输运散射。

这意味着“`-Derivatives AvalDerivatives`”仍保留某些雪崩导数，但不保留全导数臂的
全部雪崩链。仅据端口电流或日志不能拆出每个内部子项。
关闭AvalDerivatives后再比较迁移率开关，空穴自身块与Vela仍差0.858 ppm；
空穴/ψ块差约99.76%，说明更宽的迁移率电势依赖也存在，差异不全是雪崩。

Vela当前实现仍只去掉输运中的 `flux/mobility * term.quasiFermi`，保留人口、ψ链
及源项AD，见 `src/equation/CoupledDDAssembler.cpp:6030`。这是更窄的非精确Newton方案。

## 默认设置与历史扫描：已确认和未确认的边界

另建 `slot_derivative_context_20261013` 独立目录补充默认设置复核。
不写两个导数关键字的单独 `Coupled` 装配，与显式 `Derivatives AvalDerivatives`
的变量、RHS1和整个矩阵均逐位一致，仍有111条负空穴对角。耗时76.35秒。
这与T-2022.03用户手册p196关于默认全导数的说明一致。
**不能据本实验宣称Sentaurus默认省略迁移率导数，或已证明原生使用Vela同一算法。**

既有原生扫描snapshot1的空穴自身块，反而与本批`both_off`相对Frobenius差约
5.36e−15，与`full`差25.63%。这说明此前“原生没有这些负对角”的观察成立，但
其矩阵对应的迭代上下文和实际有效导数处理，不能仅靠Math日志解释。
旧轨迹与本批Load态带约1e−14的状态舍入差，且QF驱动力/内部历史未位级绑定。
目前只能报告矩阵对应关系，不能据此断言内部存在某项未公开的自动控制。

为检查求解语句上下文，本批也尝试了默认/显式全导数的同目标
`Quasistationary`，但因为载入态已经达到5mV，原生提示“目标已满足、无事可做”，
没有装配矩阵。两组均正常退出，**不计为上下文对照通过**；没有为强制求解而
偷偷改变目标偏压。其日志、输入及哈希保留于 `context_results/`。

合理的后续查证是复核原历史debug在每次装配的完整导数开关/重用状态，并以相同
Quasistationary上下文建立能够实际装配、且逐项验证输入同态的控制。该未决项不
影响本批显式四开关对照对低密度空穴导数贡献的直接证据。

## 证据与验证

本批本地目录 `build/s13native_derivative_switch/`：

- `contract.json`、`prepare.py`、`payload/`：版本与原生输入；
- `results/`：四组debug、矩阵、方向、RHS、载入及更新场量、日志；
- `download_audit.json`：64个原生文件全部通过SHA256核对；
- `analysis.json`、`supplement.json`、`negative_diagonal_audit.json`：定量结论与节点集合；
- `hole_node_deltas.csv`：逐节点空穴对角对照；各差分`.npz`保留完整稀疏支撑；
- `analyze.py`、`supplement.py`、`unpack_parse.py`：复算入口。
- `context_contract.json`、`context_analysis.json`、`context_results/`：默认设置及
  未形成矩阵的两组同目标Quasistationary记录。

远端独立目录：`/home/tcad/sentaurus_runs/vela_oracle_2022/slot_derivative_switch_20261013`。
未覆盖原生历史文件。未修改C++，未为文档单独编译求解器。
