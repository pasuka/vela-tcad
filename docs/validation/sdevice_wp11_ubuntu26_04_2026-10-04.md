# WP11 Linux 迁移至 Ubuntu 26.04（2026-10-04）

迁移验收通过。按用户最新决定，Linux 分发目标仅使用 Ubuntu 26.04 x86_64。开发计划、当前工作流、
打包指南和打包器同步更新。没有修改物理模型、材料默认、数值容差或 WP10 资格范围。

## 变更与清理

- `package_vela.py` 要求实际构建环境的 `ID=ubuntu`、`VERSION_ID=26.04` 与 x86_64；
  新回归覆盖接受 26.04、拒绝旧版本/其他发行版/其他架构。
- Linux 独立 Release 构建使用 `build-linux26/`，产出
  `build-release/wp11/packages/vela-0.1.0-ubuntu26.04-x86_64.tar.gz`。
- 删除本任务的旧隔离 rootfs（`noble/`、`noble-ext4/`）、5 GiB `noble.ext4` 镜像、
  `build-linux24/`、旧 Catch2 构建目录、Ubuntu Base 下载及校验文件、旧 Linux 包及
  清单/哈希旁文件，以及指向旧环境的构建/测试辅助脚本。
  删除前验证绝对路径并解除挂载；保留共享源码、参考资产、Windows 构建与共用离线 Python/许可输入。
- [原 WP11 记录](sdevice_wp11_packaging_concurrency_2026-10-04.md)保留当时测试的真实出处，
  已标为历史证据，不再是当前平台验收依据。

## 实际环境

Ubuntu 26.04.1 LTS，现有 WSL `Ubuntu` 发行版；GCC 15.2.0、CMake 4.2.3、
Ninja 1.13.2、系统测试 Python 3.14.4、Catch2 3.7.1、HDF5 1.14.6。
复用已安装依赖，仅补装确认缺失的 `libmetis-dev`，没有系统升级。
TDR 导入和 HDF5/HighFive 状态输出实际启用；SPQR、UMFPACK、METIS 可用，
默认后端为 Eigen SparseLU (COLAMD)。HighFive 沿用固定提交
`be0ddb3d43ce0f53db2d8b1438e819c5a5cb278a`。

分发 Python 仍为离线 CPython 3.12：Linux 的 python-build-standalone 3.12.12 与
Windows 官方 embeddable 3.12.10。输入 SHA-256、许可与命令见
[打包指南](../sdevice_packaging.md)。Linux staging 与解包测试使用原生临时文件系统；
发行包通过 RPATH 在包内解析非 glibc 依赖。

## 验证记录

实际运行日志保留于当前工作树 `build-release/wp11/`：

- `linux26-configure.log` / `linux26-build.log`：新工具链、实际特性与构建。
- `linux26-package.log` / `linux26-frontend-tests.log`：离线包与全套非参考 SDevice 回归。
- `linux26-initial-tests.log` / `linux26-concurrency.log`：保留首轮失败与修复后的并发复验。
- `linux26-final-smoke.log` / `windows26-final-smoke.log`：最终提交对应包的解包复验。
- `windows-imports.json`：最终 Windows 包 61 个 PE 文件的导入闭包审计。

Ubuntu 非参考 SDevice 首轮 CTest 12/13 通过，唯一失败是并发测试未观察到求解子进程。
根因是新构建路径超过 80 列，`ps` 截断 argv 后测试无法识别 runner basename；所有命令本身退出 0。
采样改用 `ps -ww` 后，单独复验该 CTest 通过，4/4 并发 unittest 均通过。
结合首轮已通过的其他 12 项，13 个 SDevice 检查均有通过证据。
没有跳过 Linux 矩阵或放宽任何数值/线程门限。

| 验收项 | Ubuntu 26.04 结果 |
| --- | --- |
| 同目录四个不同 cmd | PLT/VTU 与串行零差异，求解子进程最大线程数 1 |
| 四目录各一个 cmd | 零差异，最大线程数 1 |
| 同 cmd 重复启动 | 第二个退出 7，首个锁所有权保持 |
| Math/CLI 线程四组合 | 生效上限 1、2、1、2，全部符合上限，输出零差异，最大线程数 2 |
| 解包冒烟与离线打包契约 | 10/10 unittest 通过，包含 ELF 依赖闭包与中文/空格路径 |

Windows 最终包执行同一打包契约/冒烟模块（10 项，9 通过、ELF 检查按平台跳过），
并复验锁冲突（通过，3 项 Linux 矩阵按 WP11 约定跳过）。
此前 Windows 13/13 SDevice 基线见原 WP11 记录；本次没有再次运行全部 Windows 前端套件。

数值用例沿用 WP11：1 um 方形、4 节点/2 个直角三角形，施主 `1e15 cm^-3`、受主 0，
电流 A/um，原 PN2D IV cmd 夹具扫描至 10 V。并发结果与串行 PLT 逐字节、VTU
全 ASCII 数值数组以零容差比较；先清除串行可见输出，再执行并发。
重复锁的慢协议夹具只控制时间窗口，其他矩阵均使用真实 C++ runner。

两包的 `VERSION`、内部 zipapp 和 `PACKAGE_MANIFEST.json` 记录一致的完整源码提交号，
外部 `.sha256` 与 `.manifest.json` 是最终交付的校验依据。生成包与数值输出未加入 Git。
此次没有重跑全仓 CTest 或物理 PN2D IV/BV 曲线；WP12 全文档收尾仍待执行。
