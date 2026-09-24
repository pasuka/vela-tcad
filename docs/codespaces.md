# 从 Codex 使用 Codespaces 编译 Vela TCAD

本机 Codex 调用 `scripts/Invoke-Codespace.ps1`，脚本通过 GitHub CLI 的
SSH 功能在 Codespace 内执行命令。无需 VS Code 插件，也无需向公网开放 SSH 端口。

## 环境与目录

- `.devcontainer/Dockerfile`：Ubuntu 24.04，系统 GCC 13、CMake、Ninja、Catch2 3，
  Boost、Eigen、JSON、spdlog、HDF5、SuiteSparse、OpenBLAS、METIS，以及 Python
  NumPy、h5py、Pillow。与 Linux CI 使用相同的依赖类别，但 CI 单独安装 GCC 16。
- `.devcontainer/devcontainer.json`：安装 SSH 服务，建议最低 4 CPU / 16 GB。
  创建时只自动配置 CMake，不自动编译或运行仿真。
- 云端源码默认 `/workspaces/vela-tcad`；构建目录 `build-codespaces-release/`。
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
.\scripts\Invoke-Codespace.ps1 -Action Run -Codespace $cs -Command 'ls -lh build-codespaces-release'

# 下载配置日志，保存到本地被 Git 忽略的构建目录。
New-Item -ItemType Directory -Force build\codespaces-results | Out-Null
.\scripts\Invoke-Codespace.ps1 -Action Fetch -Codespace $cs -RemotePath '/workspaces/vela-tcad/build-codespaces-release/configure.log' -Destination '.\build\codespaces-results\'

# 完成并取回结果后停止计算。
.\scripts\Invoke-Codespace.ps1 -Action Stop -Codespace $cs
```

`-DryRun` 只打印命令，不联网。`-RemoteRoot` 可覆盖云端仓库位置。除 List 外要求
明确指定实例，避免误操作其他 Codespace。Sync 遇到云端源码更改或分支分叉时失败，
不会强制覆盖。它只拉取当前分支的上游，不上传本机修改，也不自动切换分支。

可以直接告诉 Codex：“在指定 Codespace 上编译并运行 Poisson 测试”。执行前应明确
实例名称及代码版本；仿真时还需指定配置、网格、输出目录和需要的求解后端。
使用 Run 提交实际仿真命令，将结果写到 `build-codespaces-release/` 下，使用 Fetch
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

本机没有 Docker，WSL 缺少完整 Python 回归依赖，本次没有构建容器镜像或执行
云端 C++ 编译/CTest；首次实例启动后仍须完成上面的 Build 和完整 Test 验收。
容器配置尚不能作为 Linux 数值回归通过的证据。

## 参考

- [GitHub 开发容器配置](https://docs.github.com/en/codespaces/setting-up-your-project-for-codespaces/adding-a-dev-container-configuration/introduction-to-dev-containers)
- [GitHub CLI SSH 与 SSH 服务要求](https://cli.github.com/manual/gh_codespace_ssh)
- [创建 Codespace](https://cli.github.com/manual/gh_codespace_create)
- [复制文件](https://cli.github.com/manual/gh_codespace_cp)
- [空闲超时与计费期间](https://docs.github.com/en/codespaces/setting-your-user-preferences/setting-your-timeout-period-for-github-codespaces)
