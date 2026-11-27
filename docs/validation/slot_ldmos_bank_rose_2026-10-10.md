# SLOT IdVd Bank–Rose 阻尼实现与验证

本批为耦合 Newton 及其 Poisson 初始化实现了可选的 Bank–Rose Algorithm Global。
默认关闭，原残差、Jacobian、变量限幅、Poisson 重修正及最终器件资格保持不变。
这里的“完整”指公开算法的 K 递推、试探、接受/拒绝及跨迭代状态管理，
不表示复刻 Sentaurus 未公开的实现，也不表示解决了 SLOT 的弱行收敛问题。

同父态 5 mV 有界对照中，两臂各接受 13 次更新后软截止。Bank–Rose 将末态
超标空穴行从基线的 2023 条减至 1259 条，但最大行比反而由 1.62655 增至
2.06986，仍远高于 0.001。新增完整合格点为 0；实现可用，本批未证明器件已修复。

## 算法和接口

依据 Bank 与 Rose 的 *Global Approximate Newton Methods*（1981）第 3 节
[作者公开论文](https://ccom.ucsd.edu/~reb/m270b/newton.pdf) 中 Algorithm Global。
令当前固定尺度残差范数为 rho，当前 Newton 方向为 d，试探为 x + alpha d：

```text
alpha = 1 / (1 + K rho)
接受条件：rho_trial <= (1 - omega alpha) rho，且 rho_trial < rho
拒绝试探：K = seed_K（原 K 为 0），否则 K *= growth
接受试探：下一次迭代 K = decay * K
```

候选须有限且人口有效。K 每次 solve 重置；搜索失败保持原态及搜索前的 K，
若 QF 参考量重整导致整次试探被丢弃，K 同样回滚。限制试探次数、最小阻尼和 K
溢出，失败返回原因，不安装被拒态。K 的单位是残差范数单位的倒数；
改变范数缩放后不能无条件沿用相同 K 参数。

`not_damped=N` 按一次计数：前 N 次更新仅检查有限性和人口有效性，
直接用 alpha=1；第 N+1 次启用上述下降条件。`not_damped=25` 因而从第 26 次
才启用；本批器件对照设为 0，以实际检验阻尼算法。此前原生 64/80 位 5 mV
对照均只用了 4 次更新，不能把那次原生成功归因于 Bank–Rose。
原生版本依据及迭代记录见[精度对照报告](slot_ldmos_native_precision_2026-10-10.md)。

以下片段放在原配置的 `solver` 中；原验收字段保留：

```json
{
  "line_search": false,
  "line_search_mode": "merit",
  "damping_factor": 1.0,
  "stable_merit_comparison": false,
  "diagnostics": true,
  "bank_rose": {
    "enabled": true,
    "not_damped": 0,
    "omega": 0.0001,
    "initial_K": 0.0,
    "seed_K": 1.0,
    "growth": 10.0,
    "decay": 0.1,
    "min_damping": 1e-8,
    "max_trials": 32
  }
}
```

`line_search=false` 关闭旧回溯，由 `bank_rose.enabled=true` 显式选择新算法，
不表示关闭 Bank–Rose 自身的试探。配置冲突、未知字段及非法参数直接报错。
`not_damped` 默认 25；示例的 0 是主动选择。未实现新的停止资格或人口更新规则。

代码位于 `include/vela/numerics/BankRose.h`、`src/numerics/BankRose.cpp`，
接入 `src/solver/NewtonSolver.cpp`，DC 通过耦合 Newton 使用它。
独立的 coupled-load-line/arclength 校正器、Gummel 和线性 Poisson 驱动不在本次范围。
诊断日志记录 K、alpha、当前/目标范数、试探次数、接受标记和拒绝原因；
最终 runner 状态 JSON 也保留失败搜索的试探历史。

公开论文的收敛理论有光滑性、方向质量等前提。现有分量限幅、正则化及源核
可能偏离精确 Newton 假设，不能将该定理作为这个器件必收敛的保证。

## 测试和构建

使用 MSYS2 UCRT64 C++20、`windows-ucrt64-release` 预设。
构建了 `vela_example_runner`、`test_bank_rose`、`test_newton_solver`、
`test_line_search_backtrack_failure`、`test_dc_sweep`；检测到 HDF5 与 UMFPACK。

| 检查 | 结果 |
|---|---|
| Bank–Rose 专项 | 7 项、71 个断言通过 |
| 完整 Newton 测试 | 176 项、4649 个断言通过 |
| 既有回溯失败测试 | 6 项、42 个断言通过 |
| 完整 DC sweep 测试 | 111 项、3953 个断言通过 |
| 筛选 CTest：BankRose、BacktrackingLineSearch、pn2d_config_templates | 17 项全部通过，部分与上述重复 |

覆盖全步过冲后恢复、K 递推及趋近全步、NotDamped 边界、非有限/无效人口拒绝、
失败回滚、范数缩放与 K 的反向补偿、严格下降、零残差与溢出、非法配置、
PN 耦合求解和 Poisson 初始化。默认关闭路径随现有 Newton/DC 回归验证。
首次扩展测试因沙箱临时目录无法原子替换 HDF5 文件失败；改用工作树内新临时目录
后完整重跑通过，未改测试断言或产品门槛。

## 5 mV 同父态对照

证据目录 `build/s10bankrose/`；脚本 `run.py`，输入/求解器 SHA-256 在 `contract.json`。
使用同一提供态 `build/s10initial05/v5_qf.h5`，不是从新冷态建立的续算链。
网格 13881 节点、27275 三角形；Vg=5 V、Vd=0.005 V、300 K、原点 0 V。
电流单位 A/µm，电势 V；内部使用既有 unit_scaling。保留 AverageBox 配套
Poisson/支撑体积、Fermi effective-band SG、PhuMob/Lombardi/HFS、SRH/Auger、
SG 电离及接触电场回退；UMFPACK、binary64、单线性线程。

两臂都由本工作树新构建的相同二进制计算。其常数 q=1.602176634e-19、
kB=1.380649e-23、eps0=8.8541878128e-12；不是先前原生常数兼容 runner，
所以不能把本批与旧报告的迭代差全部归因于阻尼。两臂之间只有阻尼选择及输出路径不同。
原 abstol=1e-9、reltol=1e-8、行比 0.001、闭合 0.01 保留；
carrier regularization=1e-8，`replace_free_row_sum` 和分量限幅后 ψ 重修正保留。

每臂最多 30 次更新、软截止 480 s、硬截止 600 s，两臂串行总上限 1320 s。
预算是本次功能验证保护，不是收敛理论条件；软截止在更新边界检查，允许收尾。
构建/单元测试与部分求解时段重叠，壁钟不可当作严格性能基准。

第 11 次更新开始分叉：原回溯用 2 次试探接受 alpha=0.5，残差
`1.7142875916724228e-12`；Bank–Rose 用 12 次试探，把 K 从 0、1、10 逐级
增至 `1e10`，接受 alpha=`0.979806167342573`，残差
`2.0544804488726007e-12`。两者都满足各自下降条件，但此步 Bank–Rose 的
残差及最大行比并未优于原回溯。

本次 `seed_K=1` 对 `rho≈2e-12` 的尺度较小，最初多个非零 K 仍给出几乎完整的
步长，造成额外装配。这是 K 单位/尺度的实际影响；可以后续单列尺度匹配对照，
不能将本次额外试探成本解释成线性后端变慢，也不能在运行中静默调参。

| 末态指标 | 原回溯 A | Bank–Rose B |
|---|---:|---:|
| 接受更新数 | 13 | 13 |
| 退出原因 | soft_time_limit | soft_time_limit |
| 总墙钟 s | 508.994 | 495.413 |
| normalizedBlockL2 | 1.650941626e-12 | 1.502429851e-12 |
| 超标电子行 | 0 | 0 |
| 超标空穴行 | 2023 | 1259 |
| 最大行比及节点 | 1.626545547，空穴 5068 | 2.069861548，空穴 9874 |
| 电子闭合比 | 8.619698466e-8 | 4.300458609e-8 |
| 空穴闭合比 | 6.443484782e-8 | 4.567785857e-8 |
| SG 漏电流 A/µm，未合格末态 | 8.042624823808215e-7 | 8.042624823753845e-7 |
| 完整资格 | 未通过 | 未通过 |

输入哈希保持不变，第 0 态所有数值 HDF5 数据集完全相同。前 10 次更新的残差、
alpha、最大行比及两载流子闭合逐值相同。Bank–Rose 共记录 24 次试探，其中
13 次接受、11 次拒绝；独立记录审计逐步核对了 K 递推、alpha 公式与下降条件。
第 12、13 次 alpha 分别为 0.9979497317870858、0.9998440284044766。
结果在 `audit.json`、`iterations_comparison.csv`，这是记录审计，不是另一次独立复求。

两臂都因时间保护结束，不能据此断言继续迭代必不收敛。行数减少与最大行比
恶化并存，也不能只用全局残差或端口电流接近宣称成功。Wall time 受同期构建/测试
影响，不据表中 13.6 s 差异主张加速。

本批冻结 runner SHA-256 为
`cd4649d02c1a3961f447f94d6d1d08b8e24347557b4706cc910e01e1fcc10766`。
随后源码补充配置防误用、QF 重整丢弃试探时的 K 回滚及 CLI 试探日志；本批无
QF 重整路径，数值递推未改变。最终 runner 另作单更新冒烟验证，绑定独立哈希，
结果见 `release_smoke/check.json`。未替换本批冻结二进制或既有状态。
最终 runner SHA-256 为
`fb45b5c24313aba1a0f5135b2d8e90d1ce039666bd408d21db6a9b1c128432dc`；
单更新耗时 44.285 s，CLI 试探历史检查通过，首步残差、行比、闭合及块残差
与冻结 B 臂逐值相同。单步上限后返回未收敛是该检查的预期范围，不计合格点。

下一步适合分开验证 K 初始尺度与受限方向的行质量，并在相同更新数预算下做
比较。Bank–Rose 的下降对象是整体范数，不能替代连通弱行的人口/方向协调；
保留原资格。暂不据此恢复长扫描，既有证据和续算账本未修改。
