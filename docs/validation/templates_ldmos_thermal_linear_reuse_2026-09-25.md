# Linked 默认配置与四方程线性对象跨点复用

日期：2026-09-25。基于本地 main `9097c983` 的后续修改。

## 实现范围

- linked D5/D4 驱动默认 `umfpack`、持久 worker、跨点线性分析复用；
  显式 `--linear-solver sparselu --no-worker --no-reuse-linear-analysis`
  保留旧独立子进程运行方式。输入 manifest 仍必须匹配所选后端。
- 四方程 `ElectrothermalPreparationContext` 持有线性对象。
  `reuse_linear_analysis` 和 `reuse_sparselu_symbolic` 均默认开启；后者虽有
  历史名称，实际控制三个支持的后端：UMFPACK、SparseLU COLAMD/AMD。
- 四方程省略后端时优先使用已编译的 UMFPACK；构建未包含 UMFPACK 时采用
  SparseLU COLAMD。显式选择不可用后端仍报错，不静默改换。
- 数值每次重新装配和分解，原有显式 adaptive-Jacobian 实验除外。
  仅在输入身份、后端/排序、矩阵尺寸和完整压缩索引匹配时复用符号分析。
  mesh/材料/几何等准备身份变化、显式关闭、后端变化会重建；普通偏压和温度
  变化更新矩阵数值，不复用旧状态或旧温度物性。
- 点求解异常、未收敛、数值分解或回代失败、外层验收拒绝均使线性缓存失效。
  失败后静态准备仍可保留；没有把非线性失败状态传给下一点。
- 检查点记录有效 `linear_policy`。恢复从新上下文开始，不保存因子；策略改变
  或旧检查点缺少该记录时拒绝恢复，避免混合不同运行配置。
- `performance.symbolic_analyses` 为本点增量，包含单独的切线准备，避免累计
  计数重复相加。`linear_object_reused` 表示对象保留，不代表结构一定命中。

## 验证

环境：本机 MSYS2 UCRT64，`windows-ucrt64-debug`，实际启用 UMFPACK、
HDF5/HighFive；CTest 并行 2，BLAS/OpenMP 环境为 1 线程。
本轮是集成正确性验证，不是 Release 性能计时。

- 构建全部目标成功，最终 CTest **858/858 通过**，总测试时间 222.34 s。
- linked、reference fixture、screening repeat 协议专项 Python **22/22 通过**。
- 三个四方程后端均验证：新状态/偏压/温度下开启和关闭复用结果一致；首点
  分析一次、同结构后续点分析零次，同时仍进行数值分解。
- 覆盖尺寸和索引变化、同非零元数量不同位置、奇异分解后恢复、输入变化、
  显式关闭/重新开启、上下文清空及异常失效。
- 电热入口 50 项 Python 回归（包含在 CTest 中）覆盖初始化及跨点轨迹一致、
  符号分析总量减少、HDF5 恢复、完成检查点保持、修改策略拒绝及外层拒绝失效。
- 早期回归在构建尚未全部完成时运行，出现旧对象缺少新策略校验和 probe
  尚未链接的问题。最终重新编译受影响对象并完整执行 CTest，未修改验收门限。

本地日志：`build/linear-reuse-final-build.log`、`build/linear-reuse-ctest.log`。
本机 h5py 提示其构建 HDF5 版本与已安装运行库版本不同；本轮相关测试仍通过，
未据此改变环境或扩展平台资格。

## 资格边界

本轮未重跑 D0 双栅压 62 点、原生局部场联合验收或三轮 Release 配对计时，
不据单元/入口测试宣称完整曲线加速或继承新默认配置的完整曲线资格。
物理公式、Newton 停止门限、原电学/热学门限和扫压策略没有调整。
