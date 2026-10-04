# SimpleMOS 历史证据合同与代码提交审阅

日期：2026-09-30。分支：`codex/simplemos-sdevice-validation`。

## 本轮范围

本轮整理并提交此前已经验证的 SimpleMOS 生产实现、执行器、测试和说明，修复历史证据检查对当前源码的错误依赖。没有修改本轮已通过的 816 点的模型、数值门槛、状态或原生参考；核心源码仅清除一个空白行的 16 个尾部空格，没有数值改动。完整矩阵结论见 [工程矩阵汇总](simplemos_engineering_joined_closeout_2026-09-30.md)。

原先 16 项失败都发生于历史证据的源码身份检查：13 项要求当前文件仍匹配 M44/M45/M46 旧哈希，3 项读取已删除的旧测试路径。历史实验身份与当前实现的数值资格需要分别验证；直接把历史哈希更新成当前源码会错误地给旧结果赋予新身份。

## 修复后的检查语义

- 原 16 份 evidence JSON、产物哈希、数值断言和失败假设均保持不变。
- 增加按 SHA256 寻址的源码归档，恢复 132 组路径/哈希对应的原始内容；记录 Git blob 和历史检出的换行形式。校验时对归档成员直接做字节哈希，不规范化换行、不回退到当前源码。归档约 1.69 MB，无求解状态或仿真结果。
- evidence 原始字节与源码映射也受固定清单约束；Windows 检出用 `.gitattributes` 保留这些 evidence 的字节。
- **仍有一个历史源码溯源缺口**：M33/M35/M37 引用的 `CoupledDDAssembler.cpp`，SHA256 `e525d828d0199dafdea947377dc6faed7cea419b5cb76fc1f51ae05d3dfe64ed`。没有找回其原始内容，仅保留已固定的元数据；测试显式约束这一缺口，未知缺口报错。归档完整性测试通过不代表这个文件已恢复，也不代表所有旧实验均可重跑。
- 新增负面测试覆盖 evidence 改动、未知 evidence、源码映射改动、归档成员缺失/损坏及伪造新缺口。历史报告测试继续检查实际产物及原数值结论。
- 当前实现由 C++ 数值性质、状态接口、执行器门槛测试和独立 816 点工程矩阵约束；这些历史检查不再声称认证当前求解器。

实现入口：[历史证据校验](../../tests/regression/simplemos_evidence_chain.py)、[归档说明](../../tests/fixtures/simplemos_historical_sources/README.md)、[负面测试](../../tests/regression/test_simplemos_historical_provenance.py)。另将工程汇总、退出码、云端任务所有权测试纳入 CTest。

## 提交审阅边界

重点复核了显式数值选项、PhuMob/Enormal/HFS 迁移率与链式导数路径、独立 SRH 源体积、split 三块高低位一致性、HDF5 无损保存、跨偏置参考变换、边界施加、线搜索回退及缓存键。审阅开始时既有生产文件与 816 点资格所冻结的 204 个核心源码文件逐字节一致。提交前仅清除了 `SplitDDOperator.h` 一行空白的尾部空格，逐行去尾空白后完全相等，前后哈希记录于 `provenance_fix/whitespace_only.json`；没有重算或替换数值证据。

已核对的限制保持显式：split 方案是 300 K、Boltzmann、已验证的 PhuMob/Enormal/HFS 组合；未验证的量子、雪崩、弧长或 Gummel 恢复组合被拒绝；独立 SRH 体积不能混入 Auger 等未校准源。仍不将此方案提升为其它器件的全局默认。

源文件、脚本、单元测试和对应文档纳入提交；未跟踪的大量诊断产物及 `build/` 下的 H5、TDR、日志、数值输出继续保留本地，不批量提交。脚本中历史实验入口保留原上下文，不作为自动执行队列。236 个待提交 Python 文件的语法检查通过。

## 验证记录

本地针对 16 个历史模块和新增合同运行 104 项测试通过。第一次 T470p 复验为 991/993：原 16 项全部通过，新发现 ASCII 文档符号和未上传的云端调度脚本依赖两项问题。两者已修复；首轮日志保留于远端 `provenance_regression_20260930`，本地 `build/simplemos_engineering_20260929/closeout/provenance_fix/`。

最终复验使用 T470p Windows UCRT64 Release，输出到新的 `D:/code-repo/vela-bench/simplemos_n17_high_retry_20260930/provenance_regression_v2_20260930`。北京时间 11:50 完成：**993/993 CTest 条目通过**，配置、构建、清单、CTest 和包装器退出码均为 0。求解器二进制前后相同，SHA256 为 `5437a84f1eb2b7e7d199386b69e4ca978a644dddd8718d5651adadc21af0dfff`。11 份结果回传并通过远端字节哈希校验，见本地 `provenance_fix/v2/retrieval_verification.json`。本次完整回归针对冻结核心源码；上述空白清理不改变其程序语义。

另外从暂存区导出独立检出，只包含将被提交的文件，发现并修复 Windows `autocrlf` 对历史产物的换行转换：对 SimpleMOS 冻结数据及 44 个证据引用的报告/图像路径设置 `-text`，原始字节、原哈希不变。复验 **121 项中 119 项通过、2 项按原逻辑跳过**；跳过项是 M9 的可选生成诊断报告，没有把未提交报告冒充检出内容。此检查包含历史模块、证据归档负面测试和工程矩阵合同；日志为 `provenance_fix/index_review_final.log`。Git 暂存差异格式检查通过。

原 816 点、6 项场量每点共 4,896 条记录及最大电流误差 0.186617395% 保持原结论。历史源码归档检查成功，不改变一个原始源码版本仍未恢复的事实。

后续仍独立保留：迁移率字段定义覆盖、受控性能优化，以及上述不可恢复历史源码的归档缺口。Codespaces 保持关闭；本轮不关闭 T470p。
