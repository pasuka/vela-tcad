# SimpleMOS 原始状态与日志归档、main 迁移（2026-10-01）

本批补齐当前 SimpleMOS worktree 内的原始 H5/TDR 和日志保存，承接
[前一批历史证据归档](../local_evidence_20260930/README.md)。不修改模型、门槛或旧仿真结果。

## 原始资产范围

扫描整个当前 worktree（包括 ignored build 目录），使用 Windows 长路径接口，
无扫描错误。共 60,489 份文件、7,772,846,298 字节，包含：

- 48 个 H5、3,033 个 TDR（包含重复副本，按原路径逐份保存）；
- 38,276 份日志/文本记录，包括 stdout、stderr、退出码、成功及失败记录；
- 18,976 份原生曲线和执行上下文，如 PLT、CMD、PAR、状态报告；
- 156 个已有压缩包原样保存，以免遗漏其中尚未展开的状态、输入和日志。

另有结构化日志补充包，收录未进入主包的 465 份 CSV/JSON Newton、尝试、收敛及
失败记录（含相关合同），10,467,811 字节。两包原路径无交集，合计 **60,954 份文件、
7,783,314,109 字节**。

精确筛选规则及分类字节数见 [inventory_summary.json](inventory_summary.json)。
所有成员的原路径、大小、mtime 和 SHA256 保存在 ZIP 内的 `manifest.json`；文件内容位于
`files/<原worktree相对路径>`。压缩包内已有的嵌套归档按整文件哈希保护，不重新解释其内容。

本次是**本地 worktree 的原始资产快照**，不从已关闭的 Codespaces 或 T470p 补拉
仅存在于这些机器的状态，也不是对整个 worktree 所有文件的备份。未单独匹配的旧 CSV
状态或其它二进制不因这份清单自动获得备份资格。原本只在远端的 816 点运行资产继续
由相应运行目录保存。归档生成期间新增的构建/回归日志属于后续验证，不属于盘点快照。

## 远端保存与取回

持久位置、归档 SHA256 和大小见 [archive_location.json](archive_location.json)，
远端完整验证结果见 [remote_verification.json](remote_verification.json)。目录位于
`sentaurus:/root/vela-validation-archive/simplemos/2026-10-01/`，不依赖 worktree，
也不使用 `/tmp`。本地原件未删除，VM 数据仍需用户自己的存储备份策略保护。

从保留这些提交的任意 checkout 根目录，用 Python 3 执行；Windows 建议短根路径。
以下仅下载、校验、恢复，不调用求解器：

```powershell
$meta = Get-Content reference_tcad/simplemos_sentaurus2022/raw_assets_20261001/archive_location.json -Raw | ConvertFrom-Json
New-Item -ItemType Directory -Force build/raw-restore | Out-Null
scp -O -o BatchMode=yes -o ConnectTimeout=10 "sentaurus:$($meta.remote_archive)" build/raw-restore/raw_assets.zip
python scripts/simplemos_raw_archive.py verify --archive build/raw-restore/raw_assets.zip --sha256 $meta.archive_sha256
python scripts/simplemos_raw_archive.py restore --archive build/raw-restore/raw_assets.zip --sha256 $meta.archive_sha256 --root build/raw-restore/restored
scp -O -o BatchMode=yes -o ConnectTimeout=10 "sentaurus:$($meta.structured_logs.remote_archive)" build/raw-restore/structured_logs.zip
python scripts/simplemos_raw_archive.py restore --archive build/raw-restore/structured_logs.zip --sha256 $meta.structured_logs.archive_sha256 --root build/raw-restore/restored
```

恢复前先校验整个压缩包、成员覆盖及逐文件哈希；已有同字节文件跳过，已有不同文件拒绝
覆盖。上述目录仅为恢复示例，历史路径在 `restored` 下重建。全量索引与工具也一并保存
到 VM 归档目录，不依赖本地 build 目录。独立 `verify` 可在 VM 的 Python 3.6 运行。

## main 迁移边界

[cherry_mapping.json](cherry_mapping.json) 对应原分支 45 个普通提交和 2 个历史合并提交，
保留每个来源与目标提交；合并提交以第一父提交重放。用户明确批准冲突文件以已验证
最终版本为准，因此部分后续修复会提前出现在中间重放提交中，**未逐个资格验证中间版本**。

[integration_identity.json](integration_identity.json) 记录原 main、备份分支对应的原提交、
重放后的 main 和整树身份。历史重放最终文件树与 `8680319e` 完全相同；另有一个对齐提交
恢复了 3 行功能等价的重复材料别名判断，以保持证据的字节身份。原 main 中无关的未跟踪
开发计划文档保持不变。本次新增的归档工具、索引和说明另行提交并 cherry-pick 到 main。

当前数值资格仍以 [engineering](../engineering/README.md) 为准。迁移后的构建/回归结果
记录在 [validation.json](validation.json)；同树验证复用此工作树的 Release 构建目录，
不把新编译视为重新运行完整 816 点商业仿真对照。

## 检出字节与验证范围

历史重放期间 `.gitattributes` 晚于证据文件进入 main，795 个文件仍带有早期的
自动换行转换。已从各自 HEAD blob 恢复原字节；main 全部 2,450 个 `-text` 路径
复核无差异。另为 SimpleMOS Python 脚本及测试明确 LF 检出规则，保留早期合同
绑定的源码字节，不改哈希或放宽断言。源码字符检查发现的 README 非 ASCII
连接符已改为 ASCII 连字符。

main 独立运行的 51 项既有、已注册 SimpleMOS CTest 命令全部通过；新归档测试
随本次提交加入。无差别发现所有历史 Python 模块的额外检查仍有失败：525 个
用例中 20 failures、19 errors、4 skips，涉及未注册的后期历史证据脚本、缺少的
外部运行目录及旧 live-source 哈希。这不属于正式 CTest 全通过的范围，不能把
本次迁移表述为所有历史探索脚本在空缓存 checkout 上可重现。原始失败日志与
最终验证日志一并保存到 VM 的 `validation_logs.zip`。
