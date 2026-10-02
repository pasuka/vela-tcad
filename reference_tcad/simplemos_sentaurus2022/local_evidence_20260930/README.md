# SimpleMOS 历史证据归档（2026-09-30）

本目录整理此前未被 Git 跟踪的 1,946 份历史诊断产物。当前工程验收仍以
[engineering](../engineering/README.md) 为准；本次没有运行求解器、修改物理模型、
验收门槛或历史报告，也没有重新认定旧的失败结果。

后续本地 H5/TDR、求解日志和已有压缩包已单独列入
[2026-10-01 原始资产归档](../raw_assets_20261001/README.md)；下面的范围仅描述本批历史诊断产物。

## 存放与范围

| 内容 | 文件数 | 原始字节数 | 存放位置 |
| --- | ---: | ---: | --- |
| 小型合同、身份、摘要及对比表 | 1,424 | 10,450,167 | Git，保留原路径和字节 |
| 详细逐点表、图和历史分析程序备份 | 522 | 356,518,305 | 虚拟机归档；原路径下的本地副本是可选缓存 |
| 完整原件归档 | 1,946 | 366,968,472 | ZIP 共 81,863,221 字节，包含上面两类全部原件 |

持久存储位于 SSH 主机别名 `sentaurus`，不依赖将被删除的 worktree，也不放在 `/tmp`：

```text
/root/vela-validation-archive/simplemos/2026-09-30/88959a459ac9a93e62e6667787ddca92d8c26e6db1b47434a1dbf420e0031d79/simplemos_historical_evidence_20260930.zip
SHA256: 88959a459ac9a93e62e6667787ddca92d8c26e6db1b47434a1dbf420e0031d79
```

[manifest.json](manifest.json) 是逐文件大小、SHA256、存放类别、筛选原因和报告引用的索引；
[catalogue.md](catalogue.md) 按历史阶段导航；
[remote_verification.json](remote_verification.json) 记录虚拟机对归档及全部成员的检查。
远端同时保存 manifest 和独立校验脚本。虚拟机磁盘仍需按用户自己的备份策略保护。

本次范围仅为上述 1,946 份未跟踪产物，**不包括 ignored build 目录中的全量 H5/TDR、
求解日志或可执行文件**，也不是整棵 worktree 的备份。删除 worktree 前，应另外核对
仍需保留的这些原始仿真资产。

## 筛选与历史链接

Git 保留不超过 65,536 字节、被历史报告直接引用或具有合同/身份/摘要/验收意义的
JSON/CSV。具体文件名规则及每一项决定记录在 manifest 中；失败、反例及未通过资格
的摘要同样保留。较大表格、生成图和旧分析备份不作为当前源码提交。

上级 `.gitignore` 只列出 522 个准确路径，未忽略整批历史目录，新文件仍会出现在
`git status`。原件未删除或改写。旧报告保持原文及原链接：新 checkout 中的小型证据
直接可读，指向归档成员的旧链接需要先恢复缓存。不要为了使旧链接可用而重新仿真。

## 校验与恢复

以下命令从未来任意完整 checkout 根目录执行，需 Python 3.9+；下载还需能通过已有
可信 SSH 配置访问 `sentaurus`。不需要 C++ 求解器或 Sentaurus 许可证。

```powershell
python scripts/simplemos_evidence_archive.py verify
New-Item -ItemType Directory -Force build/simplemos_evidence_restore | Out-Null
$archiveDir = '/root/vela-validation-archive/simplemos/2026-09-30/88959a459ac9a93e62e6667787ddca92d8c26e6db1b47434a1dbf420e0031d79'
scp -O -o BatchMode=yes -o ConnectTimeout=10 "sentaurus:${archiveDir}/simplemos_historical_evidence_20260930.zip" build/simplemos_evidence_restore/
python scripts/simplemos_evidence_archive.py verify --archive build/simplemos_evidence_restore/simplemos_historical_evidence_20260930.zip
python scripts/simplemos_evidence_archive.py restore --archive build/simplemos_evidence_restore/simplemos_historical_evidence_20260930.zip
python scripts/simplemos_evidence_archive.py verify
```

不带 `--archive` 的校验只检查 Git 证据及已存在的缓存，并明确输出
`archive_verified: false`，不代表已检查远端。恢复前会核验完整 ZIP 和所有成员的
SHA256，默认只恢复 522 个归档路径；已有相同文件跳过，已有不同文件拒绝覆盖。
`--root <checkout>` 可指定恢复目录，重复 `--path <manifest中的完整相对路径>`
可只恢复指定成员（也可用于取回归档中备份的小型文件）。

注册的 CTest `simplemos_evidence_archive` 检查索引、小型证据、准确忽略范围和
恢复工具的覆盖/路径/冲突保护，无需访问虚拟机或下载完整 ZIP。

## 本轮恢复实测

[restore_verification.json](restore_verification.json) 记录了从 Git 暂存区导出 1,424 份
小型原件、从虚拟机重新下载 ZIP、恢复 522 份归档成员的实测。恢复后 1,946 份文件
哈希全部一致，历史报告引用的 521 个目标均存在。10 项工具单测及三项相关 CTest
（归档、参考交付包、历史来源）通过，无新增求解器运行。

Windows 实测中，过深的临时导出目录触发路径长度限制；改用短目录并对该次 Git 命令
设置 `-c core.longpaths=true` 后通过。建议未来 checkout 使用短根路径，保留文件原名。
