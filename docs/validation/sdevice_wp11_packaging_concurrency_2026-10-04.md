# SDevice WP11 打包与并发验收（2026-10-04）

WP11 四项验收已通过：Windows 与 Ubuntu 24.04 独立运行包、Linux 完整并发矩阵、
包内排除参考归档与专有产物、`vela device -P` 导出并复用材料文件。
本次没有修改物理模型、材料默认、离散算法或 WP10 门限。

## 交付与复现

工作树：`codex/sdevice-cmd-frontend`。构建、打包和测试命令见
[当前打包指南](../sdevice_packaging.md)。新增 `SDeviceRuntime` install component、
离线打包器、第三方声明生成器及两组注册到 CTest 的回归测试。

交付目录为当前工作树的 `build-release/wp11/packages/`：

- `vela-0.1.0-win64.zip`
- `vela-0.1.0-ubuntu24.04-x86_64.tar.gz`
- 两个归档各自的 `.sha256` 与 `.manifest.json` 旁文件。

`VERSION` 和 `PACKAGE_MANIFEST.json` 记录完整源码提交号、Python 输入归档哈希及
每个分发文件的大小/哈希。最终归档在本次源码提交后重建，提交号与旁文件是交付身份依据。
构建产物和下载的依赖源码、许可文本保留在忽略的构建目录，未加入 Git。

运行包同时支持根目录与 `bin/` 启动器，携带 zipapp、CPython 3.12、材料、schema、
JSON 模板、运行库和完整许可文本。Windows 的 PE 导入表独立审计覆盖 61 个
EXE/DLL/PYD，导入均落在包内或 Windows 系统库/API set；Linux 在重定位后用 `ldd`
验证所有非 glibc 依赖落在包内 `lib/`，使用可继承的 `$ORIGIN/../lib` RPATH。

Windows Python 输入为官方 3.12.10 embeddable amd64；Linux 为
python-build-standalone 3.12.12、release 20260114 的 x86_64 install-only stripped 包。
两个固定 SHA-256 见打包指南，打包前验证，打包器不联网。
许可声明包含 Python/OpenSSL、实际运行库及编译期依赖的原始许可文本；缺失映射或文本即失败。

## 环境与数值用例

| 项目 | Windows | Linux |
| --- | --- | --- |
| 环境 | 本地 MSYS2 UCRT64 Release | WSL 中隔离的 Ubuntu 24.04.4 x86_64 rootfs，ext4 |
| 编译器 / CMake | GCC 16.2.0 / 4.4.3 | GCC 13.3.0 / 3.28.3 |
| 构建目录 | `build-release/` | `build-linux24/` |
| HDF5 / TDR 导入 | 静态 HDF5，实际启用 | HDF5 1.10.10，实际启用 |
| 状态输出 | HDF5/HighFive 启用 | HDF5/HighFive 启用 |
| 求解后端 | 默认 Eigen SparseLU (COLAMD) | 默认 Eigen SparseLU (COLAMD) |
| 额外可用后端 | SPQR、UMFPACK、METIS、MUMPS、STRUMPACK、SuperLU_MT | SPQR、UMFPACK、METIS |

宿主 WSL 发行版是 Ubuntu 26.04；没有将宿主二进制标成 24.04。
构建与测试均在独立 24.04 用户空间执行。Linux 使用已安装依赖、Catch2 3.7.1，
HighFive 使用 CMake 固定的 `be0ddb3d43ce0f53db2d8b1438e819c5a5cb278a`。

打包冒烟使用 WP5 的合成网格：1 um 方形、4 节点、2 个直角三角形，均匀施主
`1e15 cm^-3`、受主 0；复用原 PN2D IV 命令夹具并扫描至 10 V。
坐标单位 um，电流单位 A/um；各并发命令使用不同输出前缀。
该用例验证打包与并发行为，物理 PN2D 的资格范围仍见
[WP10 收口记录](sdevice_wp10_closeout_2026-10-03.md)。

## 验收结果

| 检查 | Windows | Ubuntu 24.04 |
| --- | --- | --- |
| 非参考 SDevice CTest：`^sdevice_`，排除 `^sdevice_e2e_` | 13/13 通过 | 13/13 通过 |
| 安装包与离线契约 unittest | 9 项，8 通过、1 项 ELF 检查按平台跳过 | 9/9 通过 |
| 并发 unittest | 锁冲突通过，3 项 Linux 矩阵按计划跳过 | 4/4 通过 |
| 实际 runner 的 mesh/deck emitter 回归 | 修复环境注入后追加 2/2 CTest 通过 | 包含在上述 13 项内 |

mesh 模块的无 importer 分支在检测到实际 importer 时按用例设计跳过；真实 mesh 加载与
deck emitter runner 用例均执行通过。

安装包在包含中文和空格的临时路径解包，PATH 只保留系统目录。
根启动器执行 `-P`，导出的 `models.par` 被同一 deck 使用；`bin/` 启动器实际求解至
10 V，PLT/VTU/PVD 完整发布，TDR 导入程序可启动。运行前后包文件清单不变。
归档与内部 zipapp 均检查不含 `reference_staging/`、`reference_tcad/`、Git 元数据、
原生 TDR/TDX/PLT/GRD/DAT/CMD 等专有参考产物。

Linux 完整矩阵结果：

- 同目录 4 个不同 cmd：PLT 逐字节、VTU 全部 ASCII 数值数组与串行一致，零容差；求解子进程最大线程数 1。
- 4 个目录各 1 个 cmd：相同比较全部通过；最大线程数 1。
- 同一 cmd 重复启动：第二个退出 7，保留首个运行的锁所有权；使用受控慢协议夹具保证锁时间窗口。
- 同时运行四组 `(Math, CLI)` 线程配置 `(未指定, 未指定)`、`(2, 未指定)`、`(2, 1)`、`(1, 2)`：
  生效上限依次为 1、2、1、2；用 `ps` 采样实际 C++ 求解子进程，所有进程符合各自上限，最大值 2；输出仍与串行零差异。

比较前移除串行可见输出，再执行并发；数值矩阵使用真实 C++ runner。
慢协议夹具仅用于重复锁测试。Windows 按 WP11 原要求只执行锁冲突，明确跳过
Linux 数值/`ps` 矩阵。未以测试退出码掩盖这些平台跳过项。

可审计日志位于 `build-release/wp11/`：`windows-tests.log`、
`windows-portability-tests.log`、`linux-tests.log`、`linux-frontend-tests.log`、
最终归档的 `windows-final-smoke.log` / `linux-final-smoke.log` 和 `windows-imports.json`。

## 验证中发现并修复的问题

1. 多进程真实 runner 暴露输出发布竞争：后台发布超过 5 秒时，定时 `join` 返回，
   最终发布与后台争用 `.pvd.tmp`，出现 E602。结束阶段改为等待 writer/reader 完成；
   新增 6 秒慢发布回归，确认最终发布不与后台重叠。
2. Windows embeddable 的原计划路径缺少 `.`，无法加载 `_socket` 等原生标准库扩展。
   `python312._pth` 固定为 `python312.zip`、`.`、`..`，不启用 `import site`；路径仍全在包内。
3. Linux CPython 的大小写敏感路径和符号链接不能在 Windows 挂载盘临时目录正确展开。
   Linux staging/冒烟使用原生临时文件系统，tar 保留链接及执行权限；包含空格的 `ldd` 路径正确解析。
4. 随包 Python 禁止生成新的 bytecode cache，避免并发写分发目录；冒烟校验包不被修改。
5. 现有 mesh smoke 的固定 Windows 路径与 deck emitter 的 Nx1 数组断言不适用于 Linux。
   测试改用当前构建 runner、正确读取单分量数据，保留原数值容差。
   Windows 测试环境辅助函数会覆盖 runner 变量，改为在辅助函数之后追加变量，
   原先跳过的真实 runner emitter 用例现已执行通过。

本次运行针对 SDevice 分发与并发，不声称重跑全仓 CTest 或物理 PN2D IV/BV 曲线。
WP12 全文档收尾与总 DoD 举证仍是后续任务。
