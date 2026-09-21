# LDMOS 性能阶段收尾与默认配置

用户决定结束进一步性能探索，将已验证的复用与后端选择固化，并把当前算例
整理为仓库验证输入。本次不修改物理模型或放宽任何数值/曲线门限。

## 默认行为

共享 `LinearSolver` 按 UMFPACK → SparseLU → STRUMPACK → MUMPS → SuperLU_MT
列出实际可用直接后端。没有显式设置时采用首项；SparseLU 始终可用，因此
未编译 UMFPACK 时自然采用 SparseLU。`VELA_LINEAR_SOLVER` 和构造参数仍可
显式指定后端。显式指定不可用后端或发生数值失败时保持报错，不隐式切换。

默认 `DCSweep` 与普通 `--dc-worker` 为主 Newton 保留同一个线性求解器。
相同结构复用排序/符号分析，每次新点清除数值因子、重新分解新系数；
同次求解中完全相同矩阵的其他右端项仍可复用数值因子。
不改变原网格节点编号，不把旧载流子状态当缓存传给新点。

固定网格通常保持稀疏模式，但方程组合、约束和边界行仍可能改变结构。
因此比较完整压缩稀疏索引，保留结构变化重建、输入文件内容变化失效和失败
请求清理；后端选择变化也会重建。不能将“通常一次分析”实现为无条件复用。
辅助 Poisson/恢复算子的原有生命周期保持不变；不同结构不能共享同一次分析。

对照入口 `--dc-worker-no-linear-reuse` 显式关闭跨点分析复用，旧
`--dc-worker-linear-reuse` 保留为兼容别名；C++ `DCSweep(false, false)` 关闭
准备与线性上下文缓存。历史性能驱动中的关闭选项仍映射到真正关闭的入口，
避免两个对照组在默认变化后实际都复用。

五后端选择作用于共享等温线性路径；独立四方程电热的
`electrothermal_linear_solver` 仍是其自身接口，冻结 R11 已显式指定 UMFPACK。
本轮不把 D5 五后端资格移植为未验证的电热五后端资格。

## 可复现参考算例

[validation_d5](../../reference_tcad/templates_ldmos_sentaurus2022/validation_d5/README.md)
补齐原来散落在忽略目录的九项中性运行依赖，约 6 MB。包括精确网格、掺杂、
AverageBox 权重、材料、配置、双栅压零压种子和完整参考曲线；原始商业工具
二进制文件未入库。只有配置路径重定位，物理设置未改变。
来源与 SHA-256 独立列出，并用 `.gitattributes` 保持已哈希输入的原始字节。

统一入口 `scripts/run_templates_ldmos_reference.py` 默认 UMFPACK、跨点复用、
双栅压各前八点；`--points 31` 才进行完整 62 点及原双栅压联合验收。
旧外部证据配置保留为历史，新的入口不需要 staging 目录或 Sentaurus 虚拟机。
短点验证、完整曲线验收及性能重复资格明确分开。

## 验证记录

MSYS2 UCRT64 Release 最终构建成功，实际启用 UMFPACK、SparseLU、STRUMPACK、
MUMPS 和 SuperLU_MT。全量 CTest **843/843 通过**（136.42 s），覆盖默认后端
选择、显式覆盖、数值失败处理、结构/数值复用和新增参考输入合同。
另行运行的相关 Python 回归 41 项及 worker 协议 2 项均通过；普通 worker 的
默认复用和显式关闭组分别检查，避免以兼容别名代替默认入口验证。

最终二进制 SHA-256：
`9a5942b7452926a19d74d9f3d83c6a52d0377604e88574901ed84533b40c4dd1`。
使用新入库输入、默认 UMFPACK/复用运行双栅压前八个精确参考点，均通过原
状态、逐行闭合与 KCL 审计；每条曲线 68 次点请求只使用一个 worker 进程。

| 栅压 | 精确点数 | Newton 更新 | 回退 | 端到端墙钟 |
| --- | ---: | ---: | ---: | ---: |
| 4 V | 8 | 369 | 0 | 261.52 s |
| 8 V | 8 | 349 | 0 | 194.60 s |

证据：`build/ldmos_closeout/reference_final_first8/summary.json`；构建与全量测试
日志见 `build-release/closeout_rebuild.log`、`build-release/closeout_ctest_final.log`。
本次仅认领前八点数值资格；初段同时运行了回归测试，墙钟不是受控性能配对，
不能据此更新性能排名或宣称新的完整曲线资格。历史完整曲线/重复计时证据
仍按其冻结配置解释，本次未重新执行完整 62 点性能批次。

## 提交与迁移

性能计时记录、驱动和测试与本次默认值/算例整理一起纳入当前分支。
临时图表目录、重复报告压缩包和生成仿真结果不纳入源代码提交。
完成验证后，按顺序将当前分支相对 main 的全部提交 cherry-pick 到本地 main，
保留 main 已有 AC 文档提交及未跟踪汇报文件。不会推送远端。
