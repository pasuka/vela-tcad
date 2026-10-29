# D4：交付入口、支持矩阵及两机复核

2026-10-23，D4 已完成。用户指定的 T470p（Windows 10 x86_64）和 C560
（Ubuntu 26.04.1 x86_64）均完成新目录中的原包入口验证；结果见
[机器可读账本](sdevice_delivery_d4_2026-10-23.json)。

## 交付内容

- [用户交付指南](../sdevice_delivery.md)：包选择、SHA-256、依赖、解压、材料输入、
  `--check` 发射 JSON、真实运行、退出码和结果检查。
- [支持矩阵](../sdevice_support_matrix.json)：10 个冻结 fixture 的精确哈希、模型和解析结论，
  三个具名 PN2D 数值变体，以及新包原生 check/rejection gates；区分原文与补 Parameter 的变体。
- README、工作流、打包指南和开发计划统一指向当前资格；旧 WP 报告及 fixture 来源记录作为历史证据。
- 配套 `vela-sdevice-delivery-20261023.zip` 包含说明、机器可读矩阵、当前证据和中性示例，
  文件哈希由 `DELIVERY_MANIFEST.json` 记录。运行包仍是 D3 已验收的两个原归档，未重新构建或改包。

运行时快照 `6a56593bf3917cc213a00aeab4f4a614befa52ca`，原包哈希及完整数值资格见
[D1–D3](sdevice_delivery_revalidation_2026-10-22.md)。本轮只有交付文档、支持账本及验证工具工作，
没有修改模型、数值门限、生产代码或模板。

## 实际验证

两台电脑均上传完整原包并核对 SHA-256，在含空格、中文的新路径解压。
检查前后逐项核验包清单文件的大小和哈希；运行环境移除 Vela/Python 覆盖变量，PATH 排除宿主 Python。
验证脚本由随包 Python 执行，真实命令调用随包 launcher、runner 和 importer。

| 检查 | T470p | C560 |
|---|---|---|
| 根入口帮助、版本及 bin 入口帮助 | 通过 | 通过 |
| 不支持的 `device -d` 返回 2 | 通过 | 通过 |
| 新目录 `device -P` 生成默认材料 | 通过 | 通过 |
| `--check` 发射两个 JSON 阶段，不生成求解 CSV | 通过 | 通过 |
| `--threads 1 --exit-on-failure` 真实短算例 | 通过 | 通过 |
| NonLocalPath 输入被 E106/退出 1 拒绝 | 通过 | 通过 |
| 10 个原 fixture 的 LF 哈希及解析退出码 | 10/10 | 10/10 |
| 运行前后包指纹 | 不变 | 不变 |

每平台 8 条 CLI 调用与 10 个 fixture 解析检查符合预期；拒绝用例按预期退出不等于运行失败。
短算例为仓库中性 4 节点/2 三角形、1 µm 方形、施主 `1e15 cm⁻³`、300 K 输入：
初始化 1 个接受点、扫描 6 个接受点、anode 到 1 V，两阶段 `goal_reached=true`。
生成 PLT、VTU/PVD 及 run manifest。长度 µm、人口 cm⁻³、电流 A/µm；使用新包 Release 实际求解器。

这是发行入口演示，未在两台电脑重跑 PN2D 全 BV，也未把原 fixture 的解析通过升级为完整器件数值资格。
完整 PN2D 数值证据由未改动的 D3 原包及其冻结输入保持。
生成 JSON 含本机绝对路径，跨平台原始 JSON 哈希不同是预期现象，不声称 D4 两机 JSON 逐字节相同。

原始输入、作业目录、命令日志和结果保存在 ignored
`build-release/sdevice_d4_delivery_20261023/`；两机结果取回后均与远端 SHA-256 相同。
公开文档不包含连接地址、账户和机器安装路径。初次连接均成功，没有新原生 Sentaurus 运行。

两机随后按交付指南的解包、复制示例、导出材料、检查、运行命令重新演练，均到达 1 V。
Windows 首次演练的 Python 传输脚本解释了 PowerShell 路径中的反斜杠，造成解析错误；
修正为原始字符串后通过。指南命令与发行包未改，失败日志保留。最终说明包仅清理机器路径
并补充此证据，文件清单重新核验；不重跑已通过的数值工作。

## 收口边界

D0–D4 交付主线已完成。当前包的公开数值承诺限于具名 PN2D；ABA、外电阻、电流目标、
Schottky/BJT 等按矩阵各自的解析、发射或 check 层级表述。Math W201 仅记录，不宣称复制
Sentaurus 收敛判据。安装库原版 NonLocalPath、R2/R3、从零纯 QS 继续独立研发。

本轮未移动工作分支、未合并到 main，也未用交付说明覆盖历史发行归档。
