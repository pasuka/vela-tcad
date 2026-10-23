# SDevice 前端发布差分与重新验收

本轮按[范围审查](sdevice_delivery_scope_audit_2026-10-22.md)执行 D1 发布差分、D2 同输入
cmd/JSON 数值对照、D3 两平台打包。资格仍限于具名 PN2D 数值输入及具名 check/rejection
输入；安装库原版 BVmethods 的 NonLocalPath、R2/R3 与从零纯 QS 属于独立研发。

本轮于 2026-10-23 收口：**D1、D2、D3 全部完成，机器可读账本为 pass**。
2026-10-07 包和旧报告保留；本次资格来自新源码、新包和全新数值运行。

## 源码身份与差分决定

审查基线为已验收运行时 `0a89942ef8c487b5dd8487884ec51ade48ba8cdc`；工作分支原 HEAD 为
`a0ed4e073efee0214c39365d3be90a12be51c624`。使用独立临时 Git index 封存实际文件，得到
验证快照 `6a56593bf3917cc213a00aeab4f4a614befa52ca`，tree 为
`73f3a45fff03da7f74fe63b4623641f05868d678`，共 5889 个 Git 可见文件。
分支 HEAD 与用户 index 均未移动；这不是把全部研发修改提交到用户分支。

两平台使用同一原始字节源码归档，其 SHA-256 为
`ff5aac4a1696949349d67943d5dd744d17194360ce9c5ae6da60a9fcbd19e64f`。
逐文件哈希和归档在 ignored `build-release/delivery_revalidation_20261022/`。
Windows 使用 UCRT64 Release；Linux 使用 Ubuntu 26.04.1 原生工具链。可连接的计算 VM 是
另一发行版，因此没有用它构建冒充 Ubuntu 的包。

| 改动组 | 发布审查决定 | 验证边界 |
|---|---|---|
| Poisson/Coupled 初始化 | 保留真实 Poisson 初始化、显式预算和 warm Newton；不能仅凭 golden 接受 | 当前 cmd、当前 JSON、旧合格 JSON 均用新 runner 数值对照 |
| Fermi/NoFermi 与 OldSlotboom | 保留显式统计联动；修复未声明统计时冗余 false 改变冻结策略编码的问题 | 有效默认仍为 false；显式 NoFermi 继续写 false，不修改物理门或指纹门 |
| 表面迁移率与 HFS 顺序 | 保留表面项先于 HFS 的值与导数修复 | Enormal＋HFS 相关测试；PN2D 验收不扩大为所有表面模型资格 |
| 连续性／Poisson 装配 | 保留冻结人口 SRH 身份、场导数与诊断入口 | 常规模式、诊断关闭和相关装配测试；不得把私有装配控制当公共默认 |
| Newton 线性方向与诊断 | 保留自由行正则化边界、同源线性求解封装及只报告采集 | 正则化/精化等显式开关保持原默认；原行门、闭合和线搜索不放宽 |
| 外电阻控制及失败追踪 | 保留 continuation、失败 trace、显式 coupled 验收及电流身份 | 外电阻默认仍 nested_scalar；check 通过不等于完整原版轨迹成功 |
| CSV 与运行输出 | 保留引号、逗号、换行转义及明确接触电流选择 | 输出回归和实际运行链 |
| NonLocalPath 实验实现 | 保留源码及测试，公共 cmd 入口仍拒绝不支持模型 | 不安装私有 probe，不增加公开支持承诺 |
| 打包与依赖 | 沿用 SDeviceRuntime install component、离线 Python、依赖闭包和相对 RUNPATH 审计 | 仅在真实 Ubuntu 26.04 构建 Linux 包；原生包运行另验 |

本轮唯一新增的公共行为修复位于 `deck_emit`：未声明 Fermi/NoFermi 时不再写多余的
`fermi_statistics_correction=false`，保持已合格 PN2D BV 策略编码。对应测试同时检查有效值
和键是否应存在；五个 PN2D stage golden 只去掉此冗余键，初始化变化仍完整保留并接受数值验证。
没有修改验收契约、数值门限、PN2D 原子 profile 或原生参考。

## 同输入数值对照

三条路径使用相同 cmd/par/grid/sidecar 和当前 runner：cmd 转换执行、其发射 JSON 直接执行、
2026-10-07 已合格 JSON 直接执行。前两条比较 CSV/VTU 字节一致性；第三条因初始化流程变化
允许轨迹不同，但仍必须通过原曲线和场量门。每条路径均新建输出，不能复用旧求解结果。

零偏和 IV 三条路径已通过。当前 cmd 与当前 JSON 的 CSV/VTU 逐字节相同。
IV 到达 10 V，但原曲线对照窗口仍仅 0.2–0.3 V；不扩大其数值资格。
BV 的当前 cmd 路径已通过：811 点成功、0 失败、到达 −20 V，求解耗时 6787.485 秒。
409 个曲线比较点的最大相对误差为 0.003336821，最大 KCL 相对量为 0.003376208；
十个固定场量锚点全部通过，包括 −19.5/−19.6/−19.7/−20 V。
契约定义的双段 log 电流拟合 knee 均为 −19.405 V，误差 0 V；该指标不是电流阈值或电离积分
定义的物理击穿电压。两条 JSON 完整对照随后均通过；当前 JSON 与 cmd 的 CSV 和全部
812 个 stage VTU（含零偏初始化）逐字节相同。旧 JSON 与当前 cmd 并非逐字节相同，但同样
满足原曲线、十锚点和 knee 门，支持保留新的初始化发射流程。
PN2D 网格为 1943 节点、3680 三角形；长度 µm、人口 cm⁻³、端口电流 A/µm、
场强 V/cm。BV 的十个场量锚点和曲线 knee 定义保持原契约，−13.2 V 不作为击穿电压。

## 回归与包验收记录

前端 210 项测试完成：202 通过、8 跳过；另一个直接覆盖统计/发射/模板的 83 项集合中
81 通过、2 跳过。这两个集合有重叠，不相加计总数。跳过的外部 runner/reference 项由后续
具名实际运行分别补证，未补齐前不记通过。

IV 首轮 JSON 直跑没有固定 launcher 的线程环境，产生末位差异；补齐四项线程变量及
HDF5 锁设置后，当前 JSON 与 cmd 的结果逐字节一致。失败目录保留，不作为产品失败或通过。
Windows 首轮全量构建曾与 runner 运行争用 exe，链接报文件占用；改用独立冻结运行副本，
后续构建和数值目录隔离。CMake 版本参数的首轮传参错误也保留在日志中。

Ubuntu 首轮完整 CTest：1082 项，1076 通过、4 个外部参考 E2E 跳过、2 项失败。
失败一为 `resistor device and named port derivatives share packed state` 的 SRH 微小跨块响应：
解析方向尺度约 4e−17，完整 double 残差相减得到零；相邻的保精度完整残差测试在原 1e−5 门下通过。
验证修订将前者隔离为无复合源的输运/电荷/端口控制，后者继续保留 SRH 及全部跨块严格检查。
该修订仅修改测试，不改变已冻结生产运行时；独立 patch/测试文件哈希在 JSON 记录。
另一失败为测试环境未将已有 `patchelf` 加入 PATH；已补环境。首轮证据完整保留。
两个导数定向控制均通过，随后 Ubuntu 全套复验完成：**1078 通过、4 跳过、0 失败**，
总计 1082 项，耗时 751.93 秒。四个跳过项为外部参考 PN2D 0V/IV/BV 和 gates，必须由新包
E2E 独立补证；CTest 的“100%”不表示这四项已经运行。

Windows 全套 CTest 同样完成：1078 通过、4 个外部参考 E2E 跳过、0 失败，耗时 348.47 秒。
两平台新包的 PN2D 0V、IV、BV 及 gates 四项 E2E 均通过，补齐上述四个 CTest 跳过项。
两条包内 BV 均完成 811 个扫描接受点并到达 −20 V，十个场量锚点、曲线和 knee 门均通过。
两包曲线最大相对误差同为 0.003336821；最大 KCL 相对量分别为 Windows 0.003376208、
Ubuntu 0.003376207。各原生参考及输入哈希在伴随账本中绑定。
Ubuntu 包首次解压选择了 Windows 挂载盘，Python terminfo 中仅大小写不同的文件与链接发生
路径折叠；归档中的链接无环。这次环境失败完整保留，同一归档改在 Ubuntu 原生文件系统
的新目录重新验收，未修改包内容。
包验收实际使用包内 Python/runner、隔离 PATH、新目录及前后包指纹核验；Linux 每个动态 ELF
的 `$ORIGIN` 相对 RUNPATH 检查通过。[机器可读账本](sdevice_delivery_revalidation_2026-10-22.json)
记录最终通过、平台跳过、初轮失败与产物哈希。重查 1556 个运行时源文件哈希，
与冻结快照无漂移；测试修订独立登记，不改生产二进制。


## 新包身份与验收范围

两包 manifest 均绑定上述 `6a56593b` 运行时快照。Windows 包含 187 个清单文件，Ubuntu 包含
5058 个（含随包 Python）；Linux 的 66 个动态 ELF 均声明基于 `$ORIGIN` 的相对 RUNPATH。
具名新包位于 ignored `build-release/delivery_revalidation_20261022/packages/`：

| 平台 | 文件 | SHA-256 |
|---|---|---|
| Windows x86_64 | `vela-0.1.0-win64.zip` | `1c7835de896c5cb3db01932a9fce7eab8e96de6f1047bc4b255780057dc353d2` |
| Ubuntu 26.04 x86_64 | `vela-0.1.0-ubuntu26.04-x86_64.tar.gz` | `c57ddc0c10de3986b8676c71598948a20d04507447246616fabb87a867381a23` |

验收运行已将包迁移到含空格和中文的目录，真实 PN2D 由随包 Python、launcher、TDR importer 和
runner 执行，运行环境 PATH 排除宿主 Python；随后核验包前后指纹。主机 Python 仅用于编排和
比较工具，不作为发行包运行依赖。两包验收后清单内文件的大小及 SHA-256 均保持一致。

| 检查 | Windows | Ubuntu 26.04 | 边界 |
|---|---|---|---|
| 全量 CTest | 1078 通过、4 跳过 | 1078 通过、4 跳过 | 外部参考四项由新包 E2E 补证 |
| 新包数值/输入 gates | 4/4 通过 | 4/4 通过 | PN2D 0V/IV/BV、具名 check/rejection |
| 包冒烟及依赖控制 | 13 通过、2 跳过 | 15/15 通过 | Windows 跳过 Linux ELF 专项 |
| 并发控制 | 1 通过、3 跳过 | 4/4 通过 | Windows 原范围为锁冲突；完整矩阵在 Linux |
| 空格/中文目录迁移与前后指纹 | 通过 | 通过 | Linux 解压必须保留大小写敏感及符号链接语义 |

Linux 并发覆盖同目录四 cmd、四目录同时运行、线程组合 `[1,2,1,2]`；
PLT/VTU 与串行逐字节一致，观测 solver 线程数不超过对应上限。Windows 的三个平台跳过
不计通过，也不据此新增 Windows 数值并发承诺。

本轮未修改生产接受策略或扩大案例范围，未把用户工作分支合并/提交为发布分支。
D4 随后已完成，见[两机入口验证](sdevice_delivery_d4_2026-10-23.md)及[交付指南](../sdevice_delivery.md)。
扩展案例完整数值承诺、安装库原版 BVmethods
NonLocalPath/R2/R3 继续作为独立任务，不构成本次 PN2D/check-only 包验收的失败。
