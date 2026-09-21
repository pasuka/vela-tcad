# RHEL 7.9 与当前 MSYS2 UCRT64 的独立开发环境对齐研究

日期：2026-09-14。范围：只读盘点、公开包元数据核查和依赖解算研究；没有
在虚拟机安装软件、修改系统运行库、启动编译或恢复暂停的 LDMOS 仿真。

## 结论与推荐布局

对齐目标应是实测的 GCC 16.2.0，而非此前举例使用的 GCC 11/13/15。
大多数依赖有相同上游版本的 conda-forge Linux 包。推荐在原 RHEL 7.9 上保留
已资格化的 GCC 11/R7 目录，新增 GCC 16.2.0 独立环境，并将 spdlog、OpenBLAS、
SuiteSparse 构建到专用依赖目录。这样既能对齐主要版本，也能控制稀疏求解与
BLAS 的运行配置。研究结果不是整套环境的运行资格。

拟定目录如下，尚未创建于虚拟机：

| 目录 | 用途 |
|---|---|
| `/opt/vela-align/mamba` | Micromamba 包缓存及管理目录 |
| `/opt/vela-align/envs/core-20260914` | GCC、头文件库、CMake、HDF5、测试工具 |
| `/opt/vela-align/overlays/gcc16-20260914` | 同版本 spdlog/OpenBLAS/SuiteSparse 自建库 |
| `/opt/vela-align/envs/analysis-20260914` | Python 后处理依赖，独立于求解器 BLAS |
| `/home/tcad/vela-build/gcc16-release` | tcad 用户拥有的新 Release 构建目录 |

环境由 root 安装/维护，源码、构建与仿真由 tcad 用户管理。只在 Vela 开发
终端显式启用；不写入全局 LD_LIBRARY_PATH，不替换系统 glibc/Python，
Sentaurus 沿用原启动环境。

## 实测平台与版本

本地通过 `pacman -Q`、`g++ --version`、CMake 缓存、pkg-config 文件、DLL
导入表与 Python 运行时共同检查。虚拟机通过 SSH 只读查询：

- RHEL 7.9 x86_64，glibc 2.17，内核 `3.10.0-1160.el7.x86_64`。
- Windows Release 缓存：GCC 16.2.0，HDF5/UMFPACK=ON，Python API=OFF，
  gprof/LTO=OFF。最终仍须用配置输出及链接结果确认实际功能。
- `libopenblas.dll` 实际导入 `libgomp-1.dll`；UMFPACK 实际导入 OpenBLAS。

下表的“同版可用”仅表示 API 返回主渠道 linux-64/noarch 构建，不表示已经
完成安装或所有传递依赖都经过运行验证。MSYS2 的 `-1/-3` 是其发行包修订，
不能与 conda build number 等同。

| 组件 | 当前 UCRT64 上游版本 | Linux 核查与建议 |
|---|---:|---|
| GCC/G++/GFortran | 16.2.0 | 三个激活包和实现包均同版可用；使用同版工具链 |
| libgcc/libstdc++ | 16.2.0 | 同版可用，包依赖声明 glibc>=2.17 |
| Binutils | 2.47 | 精确版本 API 未找到；已找到 2.46.1，满足所查 GCC 实现包的最低要求 |
| Boost | 1.92.0 | libboost-devel/headers 同版可用 |
| Eigen | 5.0.1 | 同版可用；不要再按 Eigen 3.4 配置 |
| nlohmann/json | 3.12.0 | 同版可用 |
| Catch2 | 3.16.0 | 同版可用 |
| CMake | 4.4.3 | 同版可用 |
| Ninja | 1.13.2 | 同版可用 |
| pkgconf | 3.0.7 | 同版可用，包名与 pkg-config 接口区分 |
| spdlog | 1.17.0 | 同版包存在，但现有构建排除 fmt 12.2；建议源码构建 |
| fmt | 12.2.0 | 同版可用 |
| OpenBLAS | 0.3.34 | 同版存在，但线程运行库有差异；求解器侧建议源码构建 |
| SuiteSparse | 7.12.3 | 精确版本 API 未找到；已找到 7.10.1；建议构建 7.12.3 源码 |
| UMFPACK | 6.3.8 | 精确版本 API 未找到；7.10.1 包组合提供 6.3.5；随 7.12.3 源码构建并复核 |
| HDF5 | 2.2.0 | 同版可用，选择 nompi 构建 |
| Python | 3.14.7 | 同版可用，选择普通 cp314 而非 free-threaded cp314t |
| pybind11 | 3.1.0 | 同版可用 |
| NumPy / SciPy | 2.5.2 / 1.18.1 | 同版可用，包含 cp314 构建 |
| h5py | 3.16.0 | 同版可用，应选择针对 HDF5 2.2 构建的变体 |
| Matplotlib | 3.11.1 | matplotlib-base 同版可用 |
| GDB | 17.2 | 同版可用，包含普通 Python 3.14 变体 |

GCC 实现包依赖的 libgcc、libstdc++、libgomp、libsanitizer 的同版 Linux 包
均声明 glibc>=2.17。`sysroot_linux-64=2.17` 用于编译期头文件/库口径，
不会升级宿主机 glibc，也不能让要求更高 glibc 的程序在旧系统上自动运行。
安装必须继续满足所有传递依赖的真实宿主约束，不可伪造较高 glibc 值绕过。

## 三项需要特别对齐的依赖

### spdlog 与 fmt

本地 `spdlogConfig.cmake` 明确启用外部 fmt。所查 conda-forge spdlog 1.17.0
主渠道构建分别要求 fmt `[12.0,12.1)` 或 `[12.1,12.2)`，因此不能把
`spdlog=1.17.0` 和 `fmt=12.2.0` 同时写入纯预编译包方案。

推荐从上游 v1.17.0 构建 spdlog，链接环境中的 fmt 12.2.0，开启
`SPDLOG_FMT_EXTERNAL=ON`、`SPDLOG_BUILD_SHARED=ON`。MSYS2 对应包包含
`001-spdlog_fmt_external.patch`；实施前核对补丁内容和来源提交，保留跨平台
相关部分，不能仅凭同版本号假定源码完全相同。

### OpenBLAS

本地 `openblas.pc` 实测：

```text
Version: 0.3.34
USE_64BITINT= [未启用]
DYNAMIC_ARCH=ON DYNAMIC_OLDER=OFF NO_AFFINITY=1
USE_OPENMP=ON CORE2 MAX_THREADS=64
```

这是 BLAS 的32位整数接口，不能与64位操作系统或 SuiteSparse 自身的索引
宽度混为一谈。DLL 导入表确认使用 GCC libgomp。

所查 conda-forge 0.3.34 的 OpenMP 构建依赖 LLVM OpenMP；pthreads 变体也
存在。这些包与 Windows 同版本，但线程实现不同。求解器侧推荐用 GCC/
GFortran 16.2.0 构建 OpenBLAS 0.3.34，尽量匹配上述配置，并验证实际链接
GNU libgomp。MSYS2 配方还包含补丁，需按具体来源提交区分通用修复和 Windows
专用修改；不能直接把所有 MinGW 补丁应用到 Linux。

正式计时仍用固定单线程配置并记录实际核函数/CPU能力。即使单线程运行，
不同线程运行库、编译选项和 CPU 分派也可能改变开销。

### SuiteSparse

从上游 v7.12.3 构建 Vela 需要的 UMFPACK、SPQR、CHOLMOD 及它们的完整依赖，
使用同一 OpenBLAS。上游支持：

```text
SUITESPARSE_ENABLE_PROJECTS=umfpack;spqr;cholmod
BLA_VENDOR=OpenBLAS
SUITESPARSE_USE_64BIT_BLAS=OFF
```

保留完整依赖闭包；不要把旧的 AMD/CHOLMOD/SuiteSparse_config 与新 UMFPACK
混用。安装后核查 `pkg-config --modversion umfpack` 为6.3.8、组件实际链接
路径和 OpenMP/排序选项。自建库放入单独 overlay，不直接覆盖包管理器文件。

## 配置、安装与验证顺序

1. 保存当前 Windows 包清单、DLL/工具哈希、宏与构建参数。以本报告版本为
   冻结目标，不随软件源后续更新自动变化。
2. 对 core/analysis 候选做 Linux 平台依赖解算，核对所有 `__glibc` 约束、
   普通 cp314 和 nompi 变体。通过后导出带完整版本/build/hash 的锁定清单。
3. 在虚拟机安装 core，并逐一编译/运行 C++20、线程与 Python/HDF5 冒烟程序；
   这一步检验宿主内核/运行库，而非只看包管理器成功。
4. 按 OpenBLAS → SuiteSparse、以及独立 spdlog 的顺序构建 overlay，固定
   Release、编译器和源码/补丁哈希。设置各库自身的 RPATH，避免依赖全局
   LD_LIBRARY_PATH。若要求所有二进制依赖也由 GCC16 编译，还需扩大源码
   构建范围；同上游版本的 conda 包不保证打包时使用 GCC16。
5. Vela 使用新的 Linux 构建目录，明确 GCC16、C++20、`-O3 -DNDEBUG`，关闭
   gprof/LTO，记录激活脚本注入的 CXXFLAGS/架构优化。CMake 优先搜索 overlay，
   然后 core。必须确认实际 HDF5/UMFPACK 启用，检查链接无缺失且不混入旧R7库。
6. 完成相关数值测试和完整 CTest，再用固定算例比较新旧状态、残差、KCL、
   热平衡与 Newton 轨迹。完整性能对照仍需同VM、同输入、同线程、同输出
   范围的重复配对，不把 Windows/Linux 的时差直接归因于编译器。

当前生产基线的 Python API 为OFF。后处理环境可先使用同版本 NumPy/SciPy，
并明确记录其预编译 BLAS 与求解器自建 BLAS 的区别。若要开启 Python API，
必须额外验证 NumPy 与扩展模块同进程的 BLAS/OpenMP 加载；必要时把自建
OpenBLAS 打成 conda 包或重新构建相关 Python 扩展，不能依靠搜索路径覆盖
同名动态库来宣称一致。

本地还检测到 h5py 3.16.0 的编译 HDF5 为2.1.1、运行 HDF5 为2.2.0的告警。
Linux 方案应选择匹配 HDF5 2.2 构建的 h5py，而不是复刻这一不一致。

## 候选清单与解算边界

- [核心环境候选](environments/vela_ucrt_core_20260914.yml)：不包含待源码构建
  的 spdlog/OpenBLAS/SuiteSparse，不能将它单独当成可编译运行 Vela 的完整环境。
- [数值后处理环境候选](environments/vela_ucrt_analysis_numeric_20260914.yml)：
  暂不包含 Matplotlib，保留普通 cp314、NumPy/SciPy/h5py/HDF5 的目标版本。
- [核心精确包清单](environments/vela_ucrt_core_20260914.explicit.txt)与
  [数值后处理精确包清单](environments/vela_ucrt_analysis_numeric_20260914.explicit.txt)：
  从成功dry-run计划导出，逐包URL/build/MD5固定。
- [解算与SHA256清单](environments/vela_ucrt_20260914_resolution_manifest.json)：
  同时记录计划文件哈希、精确包清单哈希及每个包的SHA256。

本轮在 Windows 研究目录下载了官方便携式 Micromamba，仅用来对 linux-64
做 `--dry-run`。没有执行软件包安装事务，也未在虚拟机安装工具。2.9.0工具
本地启动失败，2.3.2可启动。首个核心解算在240 s超时；固定独立包缓存后，
离线重试已成功返回89个包的核心环境计划，`dry_run=true`、`success=true`。
解算显式模拟Linux/glibc2.17/内核3.10及基础x86_64，未安装这些包。该结果
证明已找到元数据约束满足的核心组合，仍不替代虚拟机实际安装/运行验证。
不含Matplotlib的数值后处理环境也通过相同离线dry-run，共74个包。两套计划
中出现的显式glibc依赖均为`>=2.17,<3.0.a0`。核心包压缩归档大小合计约
348.49 MB，数值后处理约131.85 MB（包含重叠包，不等于安装后占用空间）。
额外离线重试最初因默认缓存目录不存在失败，显式限定`CONDA_PKGS_DIRS`到
已有研究缓存后成功；这属于解算工具缓存问题，不能解释为软件包不兼容。
包含 Matplotlib 的后处理解算明确返回失败：

```text
nothing provides __glibc >=2.28,<3.0.a0 needed by libraqm-0.10.5-h6406941_1
```

进一步读包依赖确认 Matplotlib 3.11.1 的 cp314 构建会引入 libraqm；只核对
Matplotlib 自身的 glibc>=2.17 声明并不足够。因此原“所有分析/绘图包一起装”
清单不能视为可行。较新的 libraqm 变体还需单独核查完整传递依赖。

建议先保留当前 Windows 绘图流程，虚拟机部署核心与数值分析工具。若需要
在虚拟机生成同版本图，单独研究兼容 glibc2.17 的 libraqm/Matplotlib 源码
构建或经确认的较旧绘图组合；不能把 glibc 版本伪造为2.28来强行解算。

实际部署前，在虚拟机将候选清单放入 `/home/tcad/vela-env-plans` 后，用
Linux 原生 Micromamba 执行以下检查（路径为拟定部署路径）：

```bash
/opt/vela-tools/bin/micromamba --no-rc create --dry-run \
  --override-channels -c conda-forge \
  -r /opt/vela-align/mamba -p /opt/vela-align/envs/core-20260914 \
  -f /home/tcad/vela-env-plans/vela_ucrt_core_20260914.yml
```

此命令只解算，core还需要三个源码库补齐。后续真正安装、冒烟测试、源码
依赖构建和 Vela 验收都属于另行执行步骤，本报告没有认领这些步骤已完成。
精确包清单用于复现已选组合，但不能代替原生环境/CPU检查；部署前不设置
较高的`CONDA_OVERRIDE_GLIBC`。本轮只有Windows跨平台解算使用真实VM参数。

## 证据与公开来源

本地研究证据根为
`reference_staging/vela_linux_ucrt_alignment_20260914/`（忽略的研究输出）：
`msys2-packages.txt`、`availability.json`、`followup_availability.json` 及逐版本
原始 API 返回。查询限定 conda-forge 主渠道 linux-64/noarch，保留元数据哈希。
包版本可用性与约束以本轮 API 返回为准；网页检索缓存可能更旧。

- [GCC 16.2.0 激活包元数据](https://api.anaconda.org/release/conda-forge/gcc_linux-64/16.2.0)
- [GCC 实现包元数据](https://api.anaconda.org/release/conda-forge/gcc_impl_linux-64/16.2.0)
- [libgcc 运行约束](https://api.anaconda.org/release/conda-forge/libgcc/16.2.0)
- [spdlog 1.17.0 包依赖](https://api.anaconda.org/release/conda-forge/spdlog/1.17.0)
- [OpenBLAS 0.3.34 运行库依赖](https://api.anaconda.org/release/conda-forge/libopenblas/0.3.34)
- [conda-forge glibc 约束说明](https://conda-forge.org/docs/maintainer/knowledge_base/#requiring-newer-glibc-versions)
- [SuiteSparse v7.12.3 构建说明](https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/v7.12.3/README.md)
- [spdlog v1.17.0 构建选项](https://github.com/gabime/spdlog/blob/v1.17.0/CMakeLists.txt)
- [MSYS2 OpenBLAS 配方](https://github.com/msys2/MINGW-packages/blob/master/mingw-w64-openblas/PKGBUILD)
- [MSYS2 spdlog 配方](https://github.com/msys2/MINGW-packages/blob/master/mingw-w64-spdlog/PKGBUILD)

MSYS2配方的 master URL 是研究定位入口；实际构建前必须固定与已安装包版本
对应的提交和补丁哈希。Windows UCRT 与 Linux glibc、目标ABI、平台补丁不会
因版本号对齐而相同，最终目标是可审计的同版本/同配置和数值资格。
