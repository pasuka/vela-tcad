# SLOT Bank–Rose 人口有效域回退复验

2026-10-11。本批将已识别的 SG 源项人口有效域异常改为有界试探拒绝，随后从同一 5 mV 合格父态尝试 10 mV。

完成 2 个诊断臂，总墙钟 2192.8 秒，新增合格点 0。更新次数指求解器接受的中间更新，不代表合格偏压点。

## 实现范围

SG 源项在取人口对数前继续要求人口有限且严格为正；不截断或替换人口。检查失败现在抛出 `PopulationDomainError`，它继承原先的 `std::invalid_argument`。Bank–Rose 只在试探构造/残差装配范围内捕获此类型，将该试探记为 `bank_rose_population_domain`，不调用该试探的范数或接受回调，也不保存为最佳试探。

阻尼激活时沿原递推增大 K、重试，仍受 max_trials、min_damping、K 溢出限制；未激活的 NotDamped 阶段仍只尝试一次完整步。失败不提交状态或 K。父态装配异常、接受回调异常以及其他配置/程序异常照常传播。原回溯线搜索未改动。

涉及 `include/vela/core/PopulationDomainError.h`、`include/vela/equation/AssemblerUtils.h`、`include/vela/equation/PhuMobSgSourceCurrent.inl`、`src/numerics/BankRose.cpp`；配置语义已补入 `docs/config_schema.md`。

## 固定条件

Vg=5 V，Vd 从 0.005 V 到 0.010 V；300 K、原点 0 V；网格 13881 节点、27275 三角形。AverageBox/Fermi/effective_band_sg、PhuMob/HFS、SRH/Auger/电离与接触回退保持原配置。oneMKL PARDISO binary64、1 线程；正则化 1e-8、replace_free_row_sum/Woodbury 保持。网格、材料、掺杂、父态与运行库路径和哈希见 contract.json。

关闭原线搜索和两种 QF 分量限幅；因此也不触发由截断启动的 Poisson 重修正。Bank–Rose 的 NotDamped=0 从第一步启用，omega=1e-4、seed_K=1、growth=10、decay=0.1、min_damping=1e-14、每步最多32试探。
A_K0 的初始 K=0；B_initial_damped 的初始 K=1762853859.37，对应初始 alpha=4.94635995858e-10。与上一批参数一致；B 初始 K 的估计取原始最大 QF 位移对应 0.05 V，后续没有逐分量限幅或强制 0.05 V 上限。

每臂最多40次接受更新，软截止1080秒、硬截止1200秒，两臂串行总上限2500秒。物理、残差/Jacobian 数值表达及最终资格未修改。

采用用户授权的 SLOT 局部策略：行比只报告；`rho <= 1e-10 OR E_raw < 1` 后仍须通过全局闭合及边界/人口保护。电子和空穴闭合门槛均为0.01；E_raw 根据未限幅、未阻尼的稳定化方向和父态 Fermi 导数计算。实际有限步误差另存，不能以极小阻尼造成的小更新宣布收敛。该 rho/E 口径不是 Sentaurus 内部算法的完全复刻。

## 结果

| 臂 | 接受更新 | 合格 | 墙钟 s | 最后 rho | 最后原始 E | 电子/空穴闭合 |
|---|---:|---|---:|---:|---:|---|
| A_K0 | 31 | False | 1084.2 | 1.1468042056560446 | 54916.13654776624 | 0.999999986951901 / 0.04262216148538171 |
| B_initial_damped | 33 | False | 1108.4 | 1.146824883577348 | 17159.89085690084 | 0.9999999869523085 / 0.042626800303515955 |

### A_K0

退出码 1；求解器原因 `soft_time_limit`。已接受步历史中共记录 32 次人口出域拒绝；不把未接受的最终搜索历史漏计为“没有发生”。
完整 status 历史的人口出域拒绝数为 32；rho 相对初态改善 2.00189878e-05。退出原因由 status 的非空 convergence_reason/failure_reason 确认，保留原 exit.json 中的空标签，不改旧记录。
末个接受步 alpha=8.719884308875231e-11，实际有限更新误差=4.788614102844476e-06，电子/空穴超标行=23/1986，最大行比=1.6383706998119305。
保存态有效载流子人口均为正、字段有限：True。
上一批同臂没有已接受状态可比较。
stderr：`空`。

### B_initial_damped

退出码 1；求解器原因 `soft_time_limit`。已接受步历史中共记录 21 次人口出域拒绝；不把未接受的最终搜索历史漏计为“没有发生”。
完整 status 历史的人口出域拒绝数为 21；rho 相对初态改善 1.9884406e-06。退出原因由 status 的非空 convergence_reason/failure_reason 确认，保留原 exit.json 中的空标签，不改旧记录。
末个接受步 alpha=4.9463697476355535e-09，实际有限更新误差=8.488211455072956e-05，电子/空穴超标行=22/1987，最大行比=1.6319805533469889。
保存态有效载流子人口均为正、字段有限：True。
上一批同臂前 6 个已接受状态逐字段一致：True。
stderr：`空`。

## 判断

两臂都能越过此前直接传播退出的人口出域试探，说明有界拒绝/增K机制按预期工作。但极小统一阻尼未建立完整输运平衡，原始更新误差仍大、电子闭合仍接近1；本次新增合格点为0。不能把“已接受中间更新”或很小的实际更新误差称为合格解。

本结果只支持修复有效域试探的控制流程，不支持把当前无限幅 Bank–Rose 参数作为可靠续算方案，也不能否定 Bank–Rose 在其他变量/方向配置下的适用性。后续分变量人口更新、线性化和停滞保护建议见同期本地 Charon/DEVSIM 核查报告；没有自动启动新的长扫描。

## 验证与证据

UCRT64 Release 预设编译通过。直接运行对应 Catch2 可执行文件：Bank–Rose 10 用例/103断言、元素盒输运25用例/7165断言、Newton177用例/4735断言全部通过。新测试覆盖有效域回退、次数/阻尼下限、NotDamped语义、父态/K保持、其他异常传播，以及实数/AD源项人口检查；既有完整源导数及数值回归同时通过。

停止逻辑10个边界/非有限/闭合用例与物理单位RMS检查通过。两臂首次 physical.bin、regularized.bin、closure_system.bin、raw.csv 均与旧无限幅基线逐字节相同；已接受步 raw 与 processed 向量相同。已冻结的输入、源码、求解器、DLL与脚本哈希未变化。

本机 h5py 报告运行 HDF5 2.2.0 与编译 HDF5 2.1.1 版本不同；已保存状态可读取，字段有限性、人口正性与哈希核对完成。该记录审计不等于独立路径复求。

证据目录 `build/s11bankrose_domain/`：contract.json、summary.json、audit.json、两臂 config/exit/status、progress 状态/矩阵/原始方向及 Bank–Rose 试探历史。本批未生成单独 performance.json；人口出域次数依据完整 status 试探历史统计。`independent_accepted` 是复用控制器的容器字段，本批两臂是 K 初值对照，不是独立电压路径。未将失败末态写入合格续算历史。
