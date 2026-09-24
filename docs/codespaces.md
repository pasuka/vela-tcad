# 从 Codex 使用 Codespaces 编译 Vela TCAD

本机 Codex 调用 `scripts/Invoke-Codespace.ps1`，脚本通过 GitHub CLI 的
SSH 功能在 Codespace 内执行命令。无需 VS Code 插件，也无需向公网开放 SSH 端口。

## 环境与目录

- `.devcontainer/Dockerfile`：Ubuntu 24.04，GCC 16、CMake、Ninja、Catch2 3，
  Boost、Eigen、JSON、spdlog、HDF5、SuiteSparse、OpenBLAS、METIS，以及 Python
  NumPy、h5py、Pillow。GCC 与 Linux CI 中的 `egor-tensin/setup-gcc@v2`
  使用同一个 `ubuntu-toolchain-r/test` 软件源，安装 `gcc-16` / `g++-16`。
  两者对齐 GCC 主版本，不锁定软件源后续发布的具体补丁版本。
- `.devcontainer/devcontainer.json`：安装 SSH 服务，建议最低 4 CPU / 16 GB。
  创建时只自动配置 CMake，不自动编译或运行仿真。
- 云端源码默认 `/workspaces/vela-tcad`；构建目录 `build-codespaces-gcc16-release/`。
  编译脚本显式指定并检查 GCC 16；原 GCC 13 的 `build-codespaces-release/` 保留，
  切换时不混用 CMake 缓存或旧二进制。
  Windows UCRT64 的二进制和 CMake 缓存不能复制到 Linux 使用。
- Release、Python API 关闭；必须检测到 UMFPACK、SPQR、METIS、HDF5 状态存储。
  MUMPS、SuperLU_MT、STRUMPACK 关闭，和现有 Linux CI 的配置一致。
- 默认编译并行数 4；`OMP_NUM_THREADS`、`OPENBLAS_NUM_THREADS`、
  `VELA_LINEAR_THREADS` 为 1，避免多个算例竞争线程。这不承诺仿真随 vCPU 数线性加速。
- 第一次 CMake 配置会联网获取项目已锁定版本的 HighFive。

## 首次启用

1. 本机已安装 GitHub CLI，并通过 `gh auth login` 登录，授权包含 `codespace`。
   脚本自动复用 `D:\msys64\ucrt64\bin` 中的 Git 和标准安装位置的 GitHub CLI，
   不修改系统 PATH。
2. 将本次配置提交并推送到 GitHub 的目标分支。**Codespaces 使用远端分支，
   不会自动看到本机未提交或未推送的更改。** 可以使用单独的配置分支。
3. 在 GitHub 的账单设置中核对 Codespaces 免费额度、设置每月 30 美元预算，
   如账户支持则启用达到预算后停止使用。脚本不修改账单设置，也不实现费用硬上限。
4. 以下命令会创建并启动实例、开始消耗额度；将分支名替换为包含配置的已推送分支，
   在交互界面选择 4 核机型：

   ```powershell
   $env:Path = "D:\msys64\ucrt64\bin;C:\Program Files\GitHub CLI;$env:Path"
   gh codespace create --repo pasuka/vela-tcad --branch YOUR_PUSHED_BRANCH --idle-timeout 240m
   ```

   如已有 Codespace，先同步配置，再执行 `gh codespace rebuild -c 实例名称`。
   初次创建若进入恢复环境，应先查 `gh codespace logs -c 实例名称`，修复容器构建错误。

## 停止、启动与环境复用

普通停止后再启动同一个 Codespace，不需要重新下载编译器或安装依赖。
源码、构建目录和结果位于持久化的 `/workspaces`，可继续使用已有二进制和 Ninja
增量编译。当前配置只有 `postCreateCommand`，没有每次启动重新安装的命令。
每次 Build 会重新运行一次 CMake 配置检查；这是检查依赖和生成规则，不是重装系统环境，
未变化的目标通常不会重新编译，已经下载的 HighFive 也会复用。

| 操作 | 环境和构建文件如何处理 |
| --- | --- |
| 停止后重新启动同一实例 | 保留已安装工具、源码、编译产物和结果 |
| 修改普通源码并 Sync / Build | 不重装工具链，仅更新代码和增量编译 |
| 修改 Dockerfile / devcontainer 并重建 | 重新构建容器，可复用 Docker 缓存；`/workspaces` 保留 |
| 删除实例后新建，或保留期到期被删除 | 重新创建环境；旧实例未取回的结果不能依赖其继续存在 |

此次从 GCC 13 切换到 GCC 16 需要同步配置并重建一次容器，然后在新的
`build-codespaces-gcc16-release/` 完整编译；以后正常启动无需重复这次迁移。
容器重建不会保留 `/workspaces` 以外手动安装的额外工具，因此长期依赖写入 Dockerfile。
停止实例也不等于暂停并保留仿真进程，长任务恢复仍需要应用自己的检查点。

## 本机 PowerShell / Codex 命令

在本地仓库根目录执行。第一次先查询实例名称，再将名称保存到 `$cs`：

```powershell
.\scripts\Invoke-Codespace.ps1 -Action List
$cs = '替换为列表中的实例名称'

# 重新配置并编译所有目标。第一次可能需要数分钟。
.\scripts\Invoke-Codespace.ps1 -Action Build -Codespace $cs -Jobs 4

# 基础筛选测试；没有匹配测试时返回失败。
.\scripts\Invoke-Codespace.ps1 -Action Test -Codespace $cs -TestRegex 'poisson|pn2d_config_templates' -Jobs 2

# 首次 Linux 环境验收应运行完整 CTest。
.\scripts\Invoke-Codespace.ps1 -Action Test -Codespace $cs -Jobs 2

# 本地代码提交、推送后，快进更新云端当前分支，再执行 Build。
.\scripts\Invoke-Codespace.ps1 -Action Sync -Codespace $cs

# 任意 Bash 命令，从云端仓库根目录执行；这里仅查看构建结果。
.\scripts\Invoke-Codespace.ps1 -Action Run -Codespace $cs -Command 'ls -lh build-codespaces-gcc16-release'

# 下载配置日志，保存到本地被 Git 忽略的构建目录。
New-Item -ItemType Directory -Force build\codespaces-results | Out-Null
.\scripts\Invoke-Codespace.ps1 -Action Fetch -Codespace $cs -RemotePath '/workspaces/vela-tcad/build-codespaces-gcc16-release/configure.log' -Destination '.\build\codespaces-results\'

# 完成并取回结果后停止计算。
.\scripts\Invoke-Codespace.ps1 -Action Stop -Codespace $cs
```

`-DryRun` 只打印命令，不联网。`-RemoteRoot` 可覆盖云端仓库位置。除 List 外要求
明确指定实例，避免误操作其他 Codespace。Sync 遇到云端源码更改或分支分叉时失败，
不会强制覆盖。它只拉取当前分支的上游，不上传本机修改，也不自动切换分支。

本机 GitHub CLI 2.101 与 Windows OpenSSH 的默认文件复制存在路径引号兼容问题。
Fetch 已使用验证通过的兼容参数，并限制远端绝对路径只能包含英文字母、数字、
`/`、`_`、`.`、`-`，拒绝空格、通配符和 shell 表达式；结果目录请使用这些字符命名。

可以直接告诉 Codex：“在指定 Codespace 上编译并运行 Poisson 测试”。执行前应明确
实例名称及代码版本；仿真时还需指定配置、网格、输出目录和需要的求解后端。
使用 Run 提交实际仿真命令，将结果写到 `build-codespaces-gcc16-release/` 下，使用 Fetch
取回结果目录。云端编译失败、测试失败或远程命令失败会传回本机错误，不会作为成功处理。

## 夜间仿真和费用边界

当前脚本提供**前台命令执行**，不是断线续跑队列。SSH 断开时，任务可能中断；
本地终端关闭或电脑睡眠后不能据此保证仿真继续。长任务需要另行配置任务持久化、
退出码记录和模型支持的检查点恢复。

Codespaces 默认空闲超时 30 分钟，可在创建时设置到 240 分钟。后台存在计算进程
不等于用户仍然活跃；`nohup` 或 `tmux` 也不能保证 Codespace 停止后继续运行。
四小时仿真应留意超时和检查点，不能把 240 分钟理解为任务完成保证。

SSH、下载等操作可能启动已停止的实例并产生计算用量；List 不会启动实例。
运行期间即使没有仿真也会消耗计算额度。Stop 后仍占用存储，下载好结果后可以手动
删除不再需要的实例。脚本不会自动删除实例，也不会自动停止可能仍有任务的实例。

## 本次配置验证（2026-09-24）

已检查容器 JSON、PowerShell/Bash 语法，八种动作的命令生成、包含空格和单引号的
路径、测试筛选转义、缺少实例参数、非法并行数以及失败退出码传递。通过新脚本
实际查询 GitHub Codespaces 列表成功。现有 Ubuntu 24.04 WSL 的离线脚本检查通过。

最初 GCC 13 配置在 GitHub Codespaces 的云端验收记录如下；此记录不作为
GCC 16 迁移后的验收证据：

- 配置分支：`codex/codespaces-compute`。
- 实例：`vela-tcad-compute-69r6rj7pvvjc54g9`，区域 SouthEastAsia，
  `standardLinux32gb`（4 vCPU、16 GB 内存、32 GB 工作区磁盘）。
- 被编译的提交：`a90d7903a1e3b940f5ea2f3dcdbfbb3a26baf182`。
- GCC 13.3.0、CMake 3.28.3、Python 3.12.3、h5py 3.10.0、NumPy 1.26.4、
  HDF5 1.10.10；UMFPACK、SPQR、METIS、HDF5 状态存储实际启用。
- Release 全目标编译：`Build -Jobs 4`，183/183 个 Ninja 步骤完成。
- 完整回归：`Test -Jobs 2`，**852/852 通过，0 失败**，CTest 实际用时 27.83 秒。
- 本地证据目录：`build/codespaces-validation/`，包含编译、CTest 和原始远端日志，
  按项目约定不提交生成日志。

这是构建及完整 CTest 验收，不是与 T470p 的性能对比，也不代表已经运行用户的
15 分钟至 4 小时生产仿真。此实例的使用示例：

```powershell
$cs = 'vela-tcad-compute-69r6rj7pvvjc54g9'
.\scripts\Invoke-Codespace.ps1 -Action Build -Codespace $cs -Jobs 4
```

实例停止后，上述 SSH 操作可重新启动它并消耗额度。实例仍受 GitHub 保留期策略管理，
如果以后已被自动删除，需要从已推送的配置分支重新创建。

## 参考

- [GitHub 开发容器配置](https://docs.github.com/en/codespaces/setting-up-your-project-for-codespaces/adding-a-dev-container-configuration/introduction-to-dev-containers)
- [GitHub CLI SSH 与 SSH 服务要求](https://cli.github.com/manual/gh_codespace_ssh)
- [创建 Codespace](https://cli.github.com/manual/gh_codespace_create)
- [复制文件](https://cli.github.com/manual/gh_codespace_cp)
- [空闲超时与计费期间](https://docs.github.com/en/codespaces/setting-your-user-preferences/setting-your-timeout-period-for-github-codespaces)
- [重建与文件持久化](https://docs.github.com/en/codespaces/developing-in-a-codespace/rebuilding-the-container-in-a-codespace)
- [CI 使用的 GCC 安装动作源码](https://github.com/egor-tensin/setup-gcc/blob/v2/action.yml)
