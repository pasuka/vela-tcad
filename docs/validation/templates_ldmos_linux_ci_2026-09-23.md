# HDF5 迁移后的 Linux 构建与 CTest

状态：完成，Debug 构建及 850/850 CTest 通过。本次验证当前
`templates-ldmos-phase-a` worktree，包含未提交修改；
不以 Git HEAD 代替源码身份，不改变物理模型及验收门限。

## 环境与证据

- 主机：`qzw-Lenovo-C560.local`，用户 `qzw`；Ubuntu 26.04.1 LTS、x86_64。
- GCC 15.2.0、CMake 4.2.3、Ninja 1.13.2、Catch2 3.7.1、Python 3.14.4。
- HDF5 1.14.6、h5py 3.15.1、numpy 2.3.5；系统 Boost、Eigen、spdlog、json。
- 本地行政证据：`build/c560_ci_20260923/`。
- 远程独立目录：`/home/qzw/vela-ci/hdf5_migration_20260923/`。
- `source.tar.gz` 为 1546 文件的快照，23,210,424 字节，SHA-256：
  `f4631e05207102a13a9a88995493e95c0f7d360f5e378a789e099c9aa9c203c3`。
  上传后逐文件核对 `source_manifest.json`，全部一致。

GitHub workflow 当前为 Ubuntu 24.04 / GCC 16。本轮是 Ubuntu 26.04 / GCC 15
的 Linux 平台资格，不等同于已执行 GitHub Actions。

## 依赖与实际配置

首次检查缺少 Ninja/Catch2，且 sudo 需要交互认证，因而先将发行版 deb 解包到
任务私有 `deps/`。随后用户在 C560 手动安装系统包，复检确认安装成功，正式
配置使用 `/usr/bin/ninja` 和系统 Catch2；私有副本未参与构建。

HighFive 首次 Git 克隆停滞，保留 `configure.log` 和 `status.json`，停止该次
配置进程组后重试。改用同一固定提交的官方归档，不改变依赖版本：

- 提交：`be0ddb3d43ce0f53db2d8b1438e819c5a5cb278a`（项目固定版本）。
- 归档 SHA-256：`2e05a3e085d65689c7ee6c66c239a925dd7c73ef95abf0333a8ca4fbb36e1885`。
- 通过 `FETCHCONTENT_SOURCE_DIR_HIGHFIVE` 指定解包目录；未修改生产 CMake。

第二次配置通过，实际启用 HDF5/HighFive 状态存储及 TDR 导入，线性后端为
Eigen SparseLU/COLAMD。虽然系统已安装 SuiteSparse 7.12.2，当前 CMake 的
pkg-config 检测没有识别 UMFPACK/SPQR；此差异单独记录，不冒称它们已启用。

## 执行方式

隔离 Debug/Ninja 构建，`cmake --build ... --parallel 2`；完整 CTest 串行，
最终使用 `OPENBLAS_NUM_THREADS=1`、`OMP_NUM_THREADS=1`，未设置 Vela 的显式
后端/线程覆盖项。顺序为配置、`ascii_sources`、构建、
注册测试清单、完整 CTest；零测试不认作通过。第二次运行在源码检查时发现
上传包遗漏 `docs/config_schema.md`：439 个源码文件的 ASCII 检查已通过，失败
来自必需的单位文档不存在，未开始编译。补充 487 个文档文本/合同输入，归档
SHA-256 为 `326d792a537911c720bf75c9c3a9d5fdeb1eb8450ef087e40b8602900f02c728`，
新增清单 `extra_manifest.json`，与原源码清单一起逐文件验证。

构建执行者为 `run_ci.py --attempt r3`，状态文件 `r3_status.json`，日志前缀
`r3_`。前两次配置/打包问题及日志保留；失败时停止后续阶段。

第三次配置和源码检查通过，2025 个源文件/文档输入散列核对一致，183 个构建
步骤成功，用时 681.85 s。首次完整 CTest 为 684/850 通过，166 项失败：
验证驱动误设 `VELA_BLAS_THREADS=1`，而默认构建未启用显式 OpenBLAS 线程控制
接口，导致求解前抛出 `OpenBLAS thread control unavailable in this build`。
这属于验证驱动环境问题，未修改生产程序、物理模型或测试容差。

移除 Vela 显式覆盖项后，以相同二进制完整重跑，最终为 **850/850 通过，0 失败，
107.70 s**。最终状态见 `r4_status.json`，日志为 `r4_ctest.log`；构建与测试进程
已退出。关键入口包括 HDF5 C++ 合同、C++/Python 互读、DC worker、四方程
硅电阻自热及电热生产运行器，均通过。

构建、两轮完整测试、测试清单和 CMakeCache 已下载到本地行政证据目录，七份
文件的 SHA-256 逐份与远端一致。最终 CTest 日志 SHA-256：
`0244ffbd1ea5c2e6379f3850ca570ac05bcf4948c92c5298893b0a66533b087a`。

本轮资格为 Ubuntu 26.04 / GCC 15 / Debug / SparseLU 的完整注册回归；不认领
Ubuntu 24.04 / GCC 16 的 GitHub Actions 已执行，也不认领未检测启用的
UMFPACK/SPQR 或 Linux 完整曲线性能已验证。此任务没有新增生产源码修改。
