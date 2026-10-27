# BVMethods：gprof 热点与诊断驱动成本定位

日期：2026-10-27。C560 配对运行及热点分析已完成；T470p 独立复核任务正在编译。
本轮新增可选 `cost_profiling` 阶段计时，缺省关闭，只作用于未安装的
`bvmethods_adaptive_bv_probe`。没有优化物理实现、改变验收规则或更新 VM 正在运行的
完整 BV 程序。机器可读证据见[分析摘要](bvmethods_gprof_hotspots_2026-10-27.json)。

## 结论

在当前冷启动前五个非零点窗口，Vela 的主要成本是**反复求值的物理与状态诊断**。
普通 Release 的状态检查占 55.48%，联立线性求解只占 1.15%。完整扫描还存在
大量未接受的更新、人口到边界缩步和外压重试，进一步放大每轮成本。
这两类原因需要分别处理；本轮没有测得同主机、同状态序列的 Sentaurus／Vela 加速比。

当前 VM 采用的是带完整诊断的私有自动恢复扫描驱动。先前把该流程的总墙钟直接与
Sentaurus 扫描比较，会混合求解算法、重复检查、输出和硬件差异。
本轮数据不能用于声称 Vela 的所有算例或生产入口都具有相同比例的瓶颈。

## 冻结条件与校验

- C560：Linux，GCC 15.2.0，GNU gprof 2.46，2 个逻辑处理器可用。
  普通 Release 与 gprof 为独立 Ninja 构建，均为 Release／O3、LTO 关闭；
  gprof 构建通过既有 `VELA_ENABLE_GPROF` 在编译及链接阶段启用 `-pg`。
- 实际线性后端 `pardiso_mkl`，单线性／BLAS 线程，MKL sequential threading layer。
  243 个生产源文件与此前冻结源包逐项哈希一致；改动限定诊断计时和测量脚本。
- 300 K、2600 节点、5020 Tri3、2152 Si 节点、7801 维带端口边框系统；
  R=1e7 Ω·μm、Joyce–Dixon／OldSlotboom、总杂质迁移率、SRH、自洽 GradQF AV，
  BTBT 关闭，与当前自动恢复任务的物理条件相同。
- 冷 Poisson 起步，自动恢复、外推和 report_only 内压预测保持；40 次正确器预算、
  方程门 1e-8、负载线门 1e-6 V、KCL 门及 RMS<1 均不改。
  本轮只把外层尝试预算设为 6，得到零偏加五个非零接受点。
  三组都按预设窗口正常返回程序退出码 3（外层预算耗尽），不属于数值发散，
  也不代表完整 BV 已完成。

| 同主机配对 | 墙钟，s | 联立更新 | 状态场对照 |
|---|---:|---:|---|
| Release，阶段计时关闭 | 128.5534 | 16 | 参考 |
| 同一 Release，阶段计时开启 | 128.7071 | 16 | ψ／fn／fp／n／p 全部逐位相同 |
| gprof 插桩，阶段计时开启 | 148.9720 | 16 | 同上，全部逐位相同 |

三组外压表和更新序列相同，首五点更新次数为 4、3、3、3、3；22 项相关测试通过，
0 失败、0 跳过。开启阶段计时的这次观测开销约 0.12%，gprof 版约增加 15.74%。
各组仅一次运行，不将该开销比例当作重复稳定性资格，也不把插桩墙钟作为正式性能结果。

## 普通 Release 的阶段墙钟

`driver.total` 为 128.6378 s。阶段 timer 为包含式计时，父项与子项不能相加；
尤其 `inspect` 包含自己的 RHS、端口、gates 和 MP ledger，端口总时间又包含其中调用。

| 阶段 | 次数 | 时间，s | 占 driver.total |
|---|---:|---:|---:|
| 状态检查 `driver.inspect` | 36 | 71.3657 | 55.48% |
| └ 物理／守恒门 `inspect.gates` | 36 | 26.4532 | 20.56% |
| └ 高精度账本及其配套重装配 `inspect.mp_ledger` | 36 | 19.3676 | 15.06% |
| 端口电流总求值（跨多个父阶段） | 248 | 22.8468 | 17.76% |
| 每轮完整器件 Jacobian | 16 | 9.5799 | 7.45% |
| 电路参数导数 | 16 | 5.0868 | 3.95% |
| 人口列高精度线性缺陷 | 16 | 3.8011 | 2.95% |
| QF 映射高精度线性缺陷 | 16 | 3.8777 | 3.01% |
| 联立 PARDISO 求解 | 16 | 1.4743 | 1.15% |
| 状态文件写出 | 12 | 0.0476 | 0.04% |

完整 Jacobian 的计时包含其下游工作；生产代码内另有 `dd.jacobian`、
`jacobian.edge_physics` 等计时，不能再次叠加。高精度账本阶段不仅包含 MP 求和，
还包含 `evaluateCarrierTermDiagnostics` 的双精度物理重装配，不能把其 19.37 s
全部归因于高精度算术。状态文件写出很小；JSON 等其他输出留在总时间及 gprof 中，
本轮没有单独量化其全部成本。

## gprof 函数热点

gprof 可归因 self 样本合计约 42.57 s，采样周期 0.01 s。以下比例属于可归因样本，
不是 148.97 s 墙钟的分布。只汇总互不重叠的 self 项，不叠加调用图 children 时间。

| 函数族 | 可归因 self 占比 |
|---|---:|
| Tri3 输运及其局部核 | 24.03% |
| Boost.Multiprecision | 15.39% |
| 迁移率相关函数 | 13.67% |
| Joyce–Dixon 统计 | 10.41% |
| AV／SG 辅助函数 | 9.57% |
| Eigen 稠密向量赋值 | 8.41% |

单函数前列包括 `BvMethodsTri3Operator::transport` 11.79%、100 位十进制 MP 乘法
11.60%、`JoyceDixonStatistics::evaluate` 8.88%、Tri3 输运核 8.43%、Eigen 向量赋值
8.41%。gprof 记录输运入口 753 次、AV 单元入口 483 次、JD evaluate 64,578,390 次。
阶段计时独立记录了 36 次状态检查与 248 次端口求值，支持重复调用链确实很重；
这些统计不能理解为每一条调用都必然需要删除。

MKL 和系统运行库未随程序插桩，其内部成本无法由此完整展开。因此用普通 Release
的阶段计时判断 PARDISO 成本，不用 gprof 外层适配器的 self 时间代替分解成本。
GNU 官方要求 `-pg` 同时参与编译与链接，并说明采样统计和调用图的适用边界：
[编译说明](https://sourceware.org/binutils/docs/gprof/Compiling.html)、
[self profile 定义](https://sourceware.org/binutils/docs/gprof/Flat-Profile.html)、
[采样误差](https://sourceware.org/binutils/docs/gprof/Sampling-Error.html)。

## 对应代码与建议顺序

1. **复用同一候选的验收结果。**
   诊断驱动在完整候选进入控制器前调用 `inspect`，选中后又调用一次。
   建议先缓存完整诊断对象，继续执行相同的全部门，并核验输出逐位一致。
2. **合并同态物理、账本与端口求值。**
   `makeArclengthConvergenceEvaluator` 重新计算 RHS，并分别获取有边界和无边界的
   载流子账本；`BvMethodsCoupledRuntime::terms` 自身又先 `evaluate` 再 `residual`。
   `evaluateCarrierTermDiagnostics` 重新创建 assembler、求 RHS 和账本；
   `makeArclengthContactCurrentFunctional` 每次端口值也重新计算输运残差。
   建议以同一个物理 Evaluation 提供这些数据，保持边界替换、源积分和端口符号契约。
3. **定位多余的矩阵及参数差分。**
   `makeArclengthSystem::parameterDerivative` 为电路参数调用两次完整器件 RHS；
   可在核对接触依赖后只计算必要边界行导数。
   BVMethods Jacobian 当前获取旧 assembler 的完整基矩阵，再替换载流子块；
   应研究只获取所需 Poisson／约束块，保留其电荷与接触导数，不能简单删项。
4. **再讨论审计频率与算术内核。**
   人口线性精度硬门仍需保留。QF 映射缺陷及高精度账本的频率可作为独立受控试验，
   不直接降精度或取消守恒判断。跨状态缓存必须有完整状态／材料／边界版本键。
5. **最后结合失败更新数做端到端配对。**
   本轮低偏五点只需 16 次更新；完整扫描的几十轮失败尝试仍可能支配总时间。
   单轮更快和更多偏压点能收敛是不同改进，不以放宽停止门获得性能结论。

以上是定位与建议，本轮没有实施这些优化。第一步应从消除重复读取同态求值开始，
以相同接受点表、全部状态场、RHS／Jacobian 和守恒报告做普通 Release 配对。
Sentaurus 已有显式 PARDISO 五点控制见[原比较](bvmethods_pardiso_first5_2026-10-25.md)，
其机器、点表和调试输出与本轮不同，不据此给出严格加速倍数。

原始 gmon、flat／call graph、配置与三组完整状态均留在 C560 独立目录。
本地小型证据归档位于 ignored 的 `build-release/bvmethods_gprof_20261027/`。
T470p 的原生 SSH 子进程曾随连接结束而退出；已在确认旧进程结束后使用独立作业组
启动复核，当前已进入 Release 编译。其结果尚未加入本报告，不用于填补 C560 测量。
