# LDMOS worktree 文件清理记录

日期：2026-09-22。工作目录：`D:/code-repo/vela-tcad/.worktrees/templates-ldmos-phase-a`。

当前状态：两个阶段均已完成。第二阶段按用户进一步授权删除历史JSON/CSV输出，
共释放116.09 GB；两阶段累计118.68 GB / 110.53 GiB。下文第一阶段保留情况是
当时快照，历史原始输出的最终保留边界见文末第二阶段记录。

用户授权整理并删除过时、不需要的文件。第一阶段清理8个目标、97个磁盘文件，
共2,588,212,433字节，约2.59 GB / 2.41 GiB。目录统计按文件逻辑长度计，
不是NTFS物理簇或压缩后占用测量。

## 已删除及依据

| 目标 | 字节数 | 删除依据与保留位置 |
|---|---:|---|
| `.chart-data-Bpvi1L/` | 3,795 | 临时图表快照，正式PPTX已归档 |
| `.chart-data-q8irxX/` | 923,233 | 临时PPTX及快照；PPTX解包内容一致，13个嵌入xlsx继续解包后也完全一致，仅打包信息不同 |
| `output/` | 3,327,707 | 仅含重复PDF；与`docs/reports/LDMOS_R11_full_report.pdf`的SHA-256一致 |
| `docs/reports/ldmos_review_2026-09-12_17.zip` | 4,553,369 | 59个成员与保留的同名报告目录逐文件SHA-256一致 |
| `build/public-formats-off/` | 71,926,859 | 已完成的关闭可选格式编译检查产物；配置及编译日志仍在`build/t470p_public_formats_20260922/` |
| `build/t470p_hdf5_full_repeats_20260922/native_tdr_renamed.h5` | 7,554,775 | 原生TDR改名的诊断副本；原文件、核验脚本及结果保留，可重新生成 |
| `reference_staging/templates_ldmos_symbolic_vm_20260917/evidence.tgz` | 1,263,760,520 | 1952个成员与保留的解包文件SHA-256完全一致 |
| `reference_staging/templates_ldmos_preparation_vm_20260917/evidence.tgz` | 1,236,162,175 | 1915个有效成员与保留文件SHA-256一致；仅3个Python 3.6字节码缓存未保留，对应`.py`源码均存在 |

两份大型归档原SHA-256及保留目录：

- symbolic：`3f18dedcf87c62229eae3cf4a11969448f98f32b24e2b0287678892091263dc5`；
  同级`copied_vm/extrapolation_20260915/symbolic_20260917/`。
- preparation：`17cbb305f676badea7888e12d0a38bac37d89bafdf669a5651a3ef9eeaf9d081`；
  同级`copied_vm/extrapolation_20260915/preparation_20260917/`。

历史报告中的归档名称和散列仍是原接收记录。第一阶段删除归档时，解包文件
保持完整；后续用户授权的第二阶段历史JSON/CSV清理将进一步移除部分逐点输出。
因此不能再把残留解包目录当成原归档的完整替代品；保留范围以第二阶段记录为准。

## 保留与目录用途

- `src/`、`include/`、`scripts/`、`schemas/`、`tests/`、`configs/`：当前实现及未提交工作。
- `reference_tcad/`：正式算例及参考输入；本轮没有移除CSV种子或进行HDF5迁移。
- `docs/reports/`：正式会议交付件；`docs/validation/`：数值资格、性能结果及失败调查。
- `build-release/`：当前可用Release程序、依赖库及增量构建状态。
- `build/t470p_*`：近期远端执行脚本、冻结源码包、结果摘要和证据归档。
  最近HDF5完整三轮的归档、104项证据及恢复状态样本保留。
- `reference_staging/`：历史原始状态、局部场、失败轨迹及重现实验依赖；
  日期较早不等于过时，未证实冗余的唯一证据继续保留。

盘点共约190.95 GiB，主要是历史JSON/CSV状态和原生TDR；第一阶段仅删除已确认
重复或可再生成的文件。没有清空历史证据目录，也没有触及其他worktree或远端主机。

## 核验与防止重复堆积

删除前解析并检查所有目标绝对路径位于本worktree内，拒绝已跟踪文件及重解析点，
再次校验待删除文件散列。删除前对2071个保留的已跟踪/未提交文件记录SHA-256；
清理后核对仅本轮文档和忽略规则发生预期变化。

`build/worktree_cleanup_20260922/`保存目录占用清单、删除前Git状态、保护散列、
删除逐文件清单和两份归档的逐成员核验结果。
`.gitignore`增加临时图表目录、PDF暂存目录及该重复报告ZIP的精确忽略规则。
本轮无求解器或物理模型修改，不为清理重新运行仿真。

最终核验：2071个受保护文件无缺失，2068个内容完全不变，另外3个仅为本轮
预期修改的`.gitignore`、文档索引和报告README；最新HDF5归档及104项证据
再次通过散列检查。报告链接检查和`git diff --check`通过。

## 第二阶段：按用户授权删除历史JSON/CSV输出

用户进一步明确建议删除历史JSON和CSV文件。本阶段只处理当前worktree中
`reference_staging/templates_ldmos*/`的未跟踪历史运行输出，没有修改远端主机、
其他worktree、正式算例或历史报告正文。

| 删除类型 | 文件数 |
|---|---:|
| 历史CSV节点状态、导出物理场及逐点诊断 | 274,709 |
| 历史JSON结果、状态及大型诊断输出 | 6,096 |
| 合计 | **280,805** |

释放116,089,149,673字节，即116.09 GB / 108.12 GiB。全部清单目标删除成功，
无变更冲突、跳过项或错误；两阶段累计删除280,902个文件，118,677,362,106字节。

保留范围：

- 全部受版本管理文件与尚未提交的代码、配置、文档；正式`reference_tcad`算例不变。
- 2343个从源码、配置或冻结输入清单识别出的直接依赖，删除前后SHA-256一致。
- 网格、掺杂、材料、几何、种子、输入配置等文件类别，以及汇总、分析、曲线、
  审计、manifest、ledger及小型未分类元数据。保留输入和元数据不等于保留全部旧输出。
- 原生TDR和本阶段未涉及的其他格式；当前Release构建；最新HDF5三轮证据及其归档。

历史报告和ledger可能仍包含已经删除的逐点路径。这些记录只保留历史事实，
不能据此认为原状态仍可读取或直接续跑；对已退役原始场的完整逐点复审需要重新
运行或取得尚存的独立归档。本阶段没有承诺已删除文件都能从其他地方恢复。
特别是第一阶段已删除的两份传输归档，其解包目录如今也只保留上述筛选后的内容。

核验结果：2072个受保护源码/文档文件均存在，仅清理文档有预期修改；2343项
直接依赖及最新HDF5归档、104项证据散列全部通过。清理前后以下22项测试均通过：

```text
python -m unittest tests.regression.test_ldmos_reference_fixture \
  tests.regression.test_templates_ldmos_linked_d5 \
  tests.regression.test_ldmos_production_export
```

详细清单在`build/worktree_cleanup_20260922/`：
`historical_delete_plan.json`、`historical_removed.jsonl`、
`historical_plan_summary.json`、`historical_protected_dependencies.json`及
`historical_verification.json`。删除仅使用PowerShell的LiteralPath文件操作，
逐文件核对目标目录、后缀、Git跟踪状态、文件类型、长度和修改时间。
