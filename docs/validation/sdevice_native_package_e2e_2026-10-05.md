# SDevice 日志修复、文档补漏与原生包 E2E（2026-10-05）

当前状态：**本次收口完成**。日志修复、文档补漏与两平台部署包原生验收均已完成；
Windows 与 Ubuntu 26.04 各四项验收通过，执行后的包文件指纹保持不变。
原契约与 −20 V 空穴密度已知限制保持不变，WP13/剩余 WP14 按用户决定保留待办。

## 修复与任务状态

- `bdddb838`：保留各段 sidecar 的 `runtime_log.profile`，固定 append 模式；
  CLI `--verbose` 优先并使所有段 debug，否则使用该段 sidecar 或 default。
  真实 C++ runner 回归检查 minimal/debug、未覆盖段和 CLI 优先级。
- 同步计划的 WP14 最小 v2 重跑状态、旧推荐顺序、栅零偏移基准，以及 `-n`/`-P` 的实际 CLI 行为。
- `8dc5498c`：修复 Windows 许可收集器只接受 CMake `FILEPATH` 的问题，
  同时支持 preset 的 `STRING` 编译器项，缺少有效路径仍拒绝打包。
- `run_sdevice_package_e2e.py`：核对解包文件指纹，用包内 CPython 启动验收控制器，
  从平衡态产生全新候选。清除主机 VELA/Python/fixture 覆盖，子进程 PATH 没有系统 Python。
  比对工具来自当前源码树，运行前端、runner、TDR importer 来自部署包；不用历史 candidate。

原 WP0–WP12（含 WP4a/WP6.1）工作包状态不变。WP13 内置材料对齐仍未开工；
WP14 的 PN2D/BJT 最小子集已完成，其余 bvmethods/栅参考不在本次运行范围。
本轮只读连接 Sentaurus VM 确认 v2 归档存在，参考文件已在本地齐备，未新增原生重跑或归档写入。

用户本轮确认先完成当前收口，其他任务保留待办。当前任务表如下；历史工作包的逐项
提交与验收见[计划附录 A](../superpowers/plans/2026-09-26-sdevice-cmd-direct-run-development-plan.md)。

| 任务 | 状态 | 后续动作/边界 |
| --- | --- | --- |
| WP0–WP3、WP4/WP4a、WP5、WP6/WP6.1、WP7–WP9 | 完成 | 接口、前端、发射器和输出已交付 |
| WP10 | 按原契约完成 | 不扩大 IV 比较窗口或 BV 空间资格 |
| WP11 | 完成 | 当前 Linux 为 Ubuntu 26.04；并发按既定平台范围验收 |
| WP12 | 完成 | 当日初始 DoD 缺口由本报告跟踪补齐 |
| 日志覆盖修复与文档补漏 | 完成 | 真实 runner 回归通过；本次同步当前证据入口 |
| 两平台部署包原生 E2E | 完成 | 两平台各三个数值案例 + gates 通过，运行后包指纹通过 |
| WP13 全局内置材料对齐 | 未开工，保留待办 | 会改变无材料文件 deck 的全局默认和基线 |
| WP14 剩余 bvmethods/栅参考 | 未执行，保留待办 | 已完成 PN2D/BJT 子集；其余不在本轮范围 |
| −20 V 空穴密度差异 | 已知限制，保留待办 | 扩大反偏场量资格需另行修复与严格复验 |

## 两平台验收

两包均记录源码提交 `8dc5498cc4562fd7a5327225cd0bbce94ea8f796`，位于忽略目录
`build-release/native_e2e_20261005/packages-final/`。Windows 包内 CPython 3.12.10，
Ubuntu 26.04 包内 CPython 3.12.12；Linux 解包到原生 `/tmp` 文件系统。

| 项目 | Windows | Ubuntu 26.04 |
| --- | --- | --- |
| 原生 PN2D 0V | pass | pass |
| 原生 PN2D IV | pass，206 点到 10 V | pass，206 点到 10 V |
| 原生 PN2D BV | pass，811 点到 −20 V | pass，811 点到 −20 V |
| gates（五项输入检查 + 一项拒绝注入，仅 `--check`） | pass | pass |
| 包执行后指纹 | pass | 修复后 pass；首轮失败保留 |
| 非参考回归 | 全量 CTest 1023/1023，356.14 s | SDevice CTest 13/13，562.61 s |

三个 PN2D 数值案例的原生参考固定为 T-2022.03-SP2 v2；每例均做 manifest 的六项核对、
源 deck/参数/网格校验。`gates` 使用 v1 的具名原生网格核对指纹，验证 Schottky、ABA、
外接电阻及电流目标的输入检查，以及热/IALMob 拒绝；不宣称这些案例已取得原生数值资格，
也不包含尚待 WP14 执行的栅参考重跑。
配置为 300 K、1943 节点/3680 三角形、2×0.5 µm Si PN 结，坐标 µm、密度 cm^-3、
电场 V/cm、电流 A/µm。两平台默认 Eigen SparseLU (COLAMD)、单线程，实际 HDF5/TDR/state 启用。
BV 使用绑定的完整 `element_edge_sg_gss_laux` 原子配置，不修改 C++ 全局默认。

验收范围保持原 [WP10 契约](../../configs/sdevice_acceptance_contract.json)：
IV 曲线比较窗口 0.2…0.3 V、场量点 0.2 V；BV 全窗口曲线/KCL、必过场量点 −13.2 V。
IV 到 10 V 是完整扫描目标，不宣称全正偏范围已取得曲线精度资格。
−13.2 V 不是击穿电压；具名分段 log(I) 膝点为双方 −19.405 V。
−19.5/−19.6/−19.7 V 诊断和 −20 V 空穴密度已知差异按契约原样保留，不放宽门限。

## 首轮完整性失败与复验

Ubuntu 首轮三个数值案例和 gates 报告全部 pass，但控制器以 `-I` 启动 CPython，导致环境变量
`PYTHONDONTWRITEBYTECODE=1` 被忽略。启动时重写包内 encodings 的既有 `.pyc`，
最终指纹检查拒绝该运行；首轮不能计为整体验收通过。

`01fb6687` 修复为控制器显式 `-I -B`，在解释器启动阶段就禁止字节码写入。回归制造无效旧 cache，
验证控制器执行新源码、正确传递参数且旧 cache 字节保持不变；Windows/Ubuntu 均通过。
随后重新解包同一发行档，在新输出目录重新执行 Ubuntu 的四项验收，保留首轮失败目录。
复验进程退出 0，最终摘要为 `status=pass`、`package_unchanged=true`；运行中另一次只读
核对也确认 5055 个清单文件没有变化。BV 求解完成 811 个接受点，runner 的
`elapsed_seconds=5452.182971`，随后曲线、平衡场及必过偏置场量均 pass。
十个偏置场量样本中九个 pass，−20 V 诊断 fail；必过点仍仅 −13.2 V。
Windows 完整四项记录来自修复控制器前的成功运行，其包指纹原本已通过；显式 `-B` 后
另做 Windows 原生 0V 与包完整性复验通过，不重复计入完整四项结果。
这是外部验收控制器修复，发行包内前端/求解器没有变化，不需要改变数值资格或重新标记包版本。

## DoD 与文档核对

WP12 初始账本的第 1 项已由本次两平台原生包验收补齐，按原契约通过。
第 2–5 项的门控、golden、无参考前端、Unicode/中止和并发证据
沿用初始账本，同时加入本轮 Windows 1023 项与 Ubuntu 13 项回归。第 6 项由本轮 Windows
非参考全量回归及文档补漏支撑，三个 PN2D BV 模板门在该全量选集中。
四个原生参考案例单独记录，不以 CTest skip 或合成包冒烟替代。

本轮临时 `audit_followup.py` 核对 67 个登记编号、60 个代码字面编号、32 个解析器模板、
43 个 C++ 字段条目，以及 182 个本地文件链接和 2 个 Markdown 锚点，均通过。
这次文档修改补充证据入口、当前状态与 `-I -B` 控制器说明，没有新增需要重跑的执行示例。

## 证据目录

全部生成产物位于忽略目录 `build-release/native_e2e_20261005/`：

- `windows/package_e2e_summary.json` 与各例 `acceptance_report.json`：Windows 完整四例及包指纹。
- `linux/package_e2e_summary.json`、`linux-e2e.log`：Ubuntu 首轮数值通过、最终完整性失败的原始证据。
- `linux-recheck/`、`linux-recheck.log`：显式 `-B` 后的全新候选和复验。
- `windows-controller-recheck/`：修复后的控制器另做 Windows 原生 0V 与包完整性检查。
- `windows-full-regression.log`、`linux-frontend-final.log`：全量非参考 Windows 与 Linux 前端/并发结果。
- `windows-controller-regression.log`、`linux-controller-test.log`：控制器缓存回归。
- `packages-final/*.sha256`、`*.manifest.json`：发行档哈希、源码身份及包内文件清单。

Windows CTest 的四个参考项（三个数值案例 + gates）不在本轮 1023 个非参考项中重复运行；使用同一
`sdevice_acceptance.py` 从两平台解包入口执行这四个案例，单独记录数值报告。
Ubuntu 没有本轮全仓 CTest 结果；13 项是前端及并发范围。Windows 并发仍按 WP11
约定只执行锁冲突，其余 Linux 矩阵按平台跳过，不将它们伪记为 Windows 已通过。

## 固定指纹

发行档与原始摘要记录的 SHA-256 如下；摘要中的唯一报告路径可定位各原始文件。
生成产物保留于上面的忽略目录，本报告只登记小型指纹证据。

| 平台 | 发行档 SHA-256 | PACKAGE_MANIFEST SHA-256 |
| --- | --- | --- |
| Windows | `ce4d18b7b87d2418937126a91dc32f48b6dca51f84652a429b473bb09fbd1540` | `58c49a048c331eb29963252515ac9a6e9ac3ffaff045e7f7f2781371fdb7a256` |
| Ubuntu 26.04 | `4c245c51e3198c7bc0d44594fc133c5fe1d29219c1384fee2d4a18b29661e944` | `0948ea52ec41f7dff948f661fbaba636b0b3d9ab45ae5cfed284c0ae9e743495` |

| 平台 | 案例 | acceptance_report SHA-256 |
| --- | --- | --- |
| Windows | pn2d_0v | `f6999a2fc3846cd17b4459baa88c49d116c43d4d53fcb9ad5400904e07ba228b` |
| Windows | pn2d_iv | `6f0320bb4ee9f1bf50f19031e4707406bd8e42a6deadf40461b6475fdd8a194e` |
| Windows | pn2d_bv | `c57f36deeb3ef0d7f3bee9af8abcdcb01152374d4dd825f7f3dd68d4d62c0795` |
| Windows | gates | `9715f1da5eb331f248092c78e911410c73f8a5f1287e8ce161c839a3c12381e0` |
| Ubuntu 26.04 | pn2d_0v | `1448902ff3af08628e80c5c2c479c87a3fcd973a94ab3c1b6e2fb0beae80addd` |
| Ubuntu 26.04 | pn2d_iv | `2767428c1be2b7f3db0aef0726b4f1850d74a9dac7eada158da9e8cc7037e97b` |
| Ubuntu 26.04 | pn2d_bv | `00b0a3a11b618591dfb78ca18e581ce35f98a8fe9c07b71a013d34fb0321c026` |
| Ubuntu 26.04 | gates | `4367ec65ec28d2d06c538e3627204a8428cbe9e29c74cc8518fea842bc55afb3` |

验收契约 SHA-256：`88c81240e90182d9d8a07bd07825454eea19439de02dc381357b59db7fd67e3d`。
比对控制器 `sdevice_acceptance.py` SHA-256：`30f8909a770617296c669489a0a20edff66f48dd184d96b1032f89131276eb26`。
显式 `-B` 的入口 `run_sdevice_package_e2e.py` SHA-256：`4a1c864e1b0335dfac2dd6d4718bc6a1582366ca2c6fd3921f9227e51b0f715b`。
