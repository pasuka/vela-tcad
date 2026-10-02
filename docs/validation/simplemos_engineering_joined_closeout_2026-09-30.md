# SimpleMOS 完整工程矩阵汇总与回归（2026-09-30）

## 数值结果

本轮将三处已封存结果按 `(case, index)` 汇总：Codespaces 13 条曲线 663 点、T470p 原目录 n17/n18 低 Vd 两条曲线 102 点、T470p 独立重跑 n17 高 Vd 51 点。**16 条曲线、816 个工作点数值资格通过**，键无重复或缺失；与 M60 原生导出对应的六项场量共 4,896 条记录。此结论不等于全 CTest 已通过，也不覆盖新网格、温度或其它模型。

范围为 n17–n24、Vd=0.05/1 V、Vg=0–2.5 V（步长 0.05 V），300 K，原 SimpleMOS 网格。显式工程组合保持 OldSlotboom、PhuMob、Enormal、高场饱和、掺杂依赖 SRH、完整 split 状态、既定几何与源体积策略，SparseLU 后端。电流单位 A/µm；物理场电势为 V，密度对比为 dex。

汇总器重新检查所有逐点通过标志、电流 2% 门槛、双初始化门槛、偏置网格及六项物理场覆盖。载流子行、守恒和端口资格沿用已经由各运行器冻结并验证的记录；不把普通原生字段比较新增为接受门槛。752 个下载证据文件的 SHA256 已核验，Windows 与 Linux 的 204 个 C++ 源/头文件原始字节哈希一致；两个平台的二进制哈希分别保留。

| 指标 | 完整矩阵最大值 | 冻结门槛 |
| --- | ---: | ---: |
| 相对 M60 的绝对电流误差 | 0.186617394672% | 2% |
| 双初始化电势差 | 6.811700358e-11 V | 1e-6 V |
| 双初始化电子准费米势差 | 2.551528712e-8 V | 1e-6 V |
| 双初始化空穴准费米势差 | 5.411943084e-8 V | 1e-6 V |
| 双初始化密度相对差 | 2.093435236e-6 | 1e-4 |
| 双初始化 Id 相对差 | 1.354681034e-9 | 1e-6 |

最大电流误差位于 n23、Vd=0.05 V、Vg=0.05 V，冷启动态为 −0.186617394672%。本次 n17 高 Vd 重跑 51 点全部通过，最大电流误差 0.004708461266%，此前在 1.85 V 的 `std::bad_alloc` 未复现；此次成功不能单独证明旧异常的根因。

## 原生物理场对比

以下统计是当前 Vela 独立原生初态自洽结果对 M60 binary64 导出字段的描述性比较，不是两个 Vela 初态之间的差，也不代表更高精度原生场。

| 场量 | 所有工况的最大指标 | 所在工况 |
| --- | ---: | --- |
| 电势 ψ | 3.524306155e-6 V | n24、高 Vd、Vg=1.00 V |
| 电子准费米势 | 3.080814194e-6 V | n21、高 Vd、Vg=2.50 V |
| 空穴准费米势 | 2.743078200e-6 V | n24、高 Vd、Vg=0.90 V |
| 电子密度绝对 log10 比 | 6.286562226e-5 dex | n23、高 Vd、Vg=0.90 V |
| 空穴密度绝对 log10 比 | 5.918030146e-5 dex | n24、高 Vd、Vg=0.95 V |
| SRH 空间加权相对 L1 | 1.648342076e-5 | n18、低 Vd、Vg=2.50 V |

SRH 指标是空间积分意义下的差，不是过零节点的逐节点相对误差。节点/单元/边迁移率定义的独立覆盖仍需单列，不由此表推断通过。

## 空退出码核查

T470p 重跑的原 `job.exit` 为零字节，保留不变。隔离试验使用相同的 `Start-Process -PassThru`、输出重定向及随后 `WaitForExit()`：实际返回 0 和 7 的短命令均复现读取到空 `ExitCode`。因此可确认包装器存在退出码记录缺陷，不能把空内容解释为 0，也不能据此反推原 Python 进程的真实退出码。

另行执行冻结矩阵驱动的**只读缓存重放**：所有 204 次请求仅准许读取已有零退出状态，任何缓存缺失、启动子进程或改写原结果都会报错。51 点重新计算的比较/场量与原文件一致，162 个 H5 及最终 summary 哈希未变，新校验过程捕获退出码 0。这是新增的独立验证，不冒充恢复原退出码。

新的完整回归驱动使用 Python `subprocess.run().returncode` 记录整数退出码；真实子进程返回 0/7 的两项测试均通过。汇总器另有六项测试，覆盖边界值、重复/缺失点、错误偏置、非有限值、失败资格和物理场覆盖。

## 证据与执行

本地生成文件均在忽略目录 `build/simplemos_engineering_20260929/`：

- `cloud_completion/`：93 份云端证据。
- `closeout/local_cases/`：437 份低 Vd、源码身份与只读重放证据。
- `t470p_retry/completion/`：222 份高 Vd 重跑证据。
- `closeout/joined/`：完整 `summary.json`、`comparison.csv`、`physical_fields.csv`、来源表及哈希封条。

可复核入口：[汇总器](../../scripts/join_simplemos_engineering_20260930.py)、[汇总测试](../../tests/regression/test_simplemos_engineering_join.py)、[回归驱动](../../scripts/run_simplemos_engineering_joined_regression_20260930.py)、[退出码测试](../../tests/regression/test_simplemos_regression_exit_capture.py)。8 项新测试本地通过。

T470p 使用原冻结源码目录 `D:/code-repo/vela-bench/simplemos_engineering_20260929/source`。新回归输出为 `D:/code-repo/vela-bench/simplemos_n17_high_retry_20260930/full_regression_20260930`，不覆盖旧 postmatrix 失败记录。Windows UCRT64 Release 配置、构建及测试清单导出退出码均为 0。实际启用 HDF5/TDR、独立 HDF5 状态接口，以及 SPQR、UMFPACK、STRUMPACK、SuperLU_MT、MUMPS、METIS；工程矩阵沿用显式 SparseLU，不能把可用后端清单当作矩阵实际后端。

北京时间 11:11 完整 CTest 结束：**973/989 通过，16 项失败**，用时 150.56 秒；CTest 退出 8，包装器如实记录 `job.exit=1`。结果下载至 `closeout/regression/`。与 9 月 27 日 `continuation-v2_ctest.log` 比较，失败名称集合完全相同，没有新增或消除的失败项：

- 13 项停在历史冻结源码身份断言：M9、M10、M12、M30、M31、M32、M36、M39、M40、M41、M42、M45、M46。当前源码不等于这些历史证据固定的版本，不能通过直接刷新哈希把历史数值结果宣称为当前实现的复验。
- 3 项停在历史路径读取：M33、M35、M37 引用不存在的 `tests/test_mos_mixed_material.cpp`。这属于历史证据路径缺口，不是本轮求解器的数值失败。

逐项断言、分类及旧日志哈希见 `closeout/regression/failure_audit.json`，完整输出见 `ctest.xml`。本轮未跳过这些测试、修改其历史证据或放宽门槛；**完整回归仍未通过**。工作树与冻结平台的 204 个核心源码文件逐字节哈希一致；回归后独立核验 1,790 个源码、构建配置和测试脚本文件未变，求解器 SHA256 仍为 `5437a84f1eb2b7e7d199386b69e4ca978a644dddd8718d5651adadc21af0dfff`。下载的 7 份回归证据与远端 SHA256 一致，见 `retrieval_verification.json`；控制器、构建、求解器与 CTest 进程均已退出。8 项新增汇总/退出码单测独立执行，不计入远端冻结快照的 989 项。

本轮尚未提交既有大量 WIP。后续先建立历史证据与当前源码之间可追溯的测试合同，处理上述 16 项历史失败，再完成代码审阅与提交；迁移率字段定义覆盖和受控性能工作继续单列。Codespace 保持关闭。
