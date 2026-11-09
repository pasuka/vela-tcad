# Sentaurus SDE → Gmsh 离线工作流交接（2026-09-05）

## 1. 交接目标

本文件用于在新对话中继续维护 `Sentaurus SDE/SDevice → Device IR → Gmsh →
Vela IV/BV` 工作流。当前功能开发 M0–M5 已完成并通过已有验收；后续工作的重点不是
重新实现功能，而是同步最新 `main`、清理仓库卫生问题、重新执行最终回归并完成集成。

## 2. 仓库与分支状态

- 主仓库：`D:\code-repo\vela-tcad`
- 开发 worktree：`D:\code-repo\copilot-worktrees\vela-tcad-dev`
- 开发分支：`codex/sentaurus2022-sde-workflow`
- 本交接编写前开发分支 HEAD：`61bac71`（`ci: validate Gmsh integration on Ubuntu`）
- 本地 `main` HEAD：`8e99e57`（`Promote validated TCAD fixtures and remove examples`）
- 2026-09-05 实测 `git rev-list --left-right --count main...HEAD`：`22 19`，即
  `main` 独有 22 个提交、开发分支独有 19 个提交。
- 两个 worktree 在本交接编写前均无未提交修改。

注意：计划及 M5 验收文档中的“分支落后 0 个提交”只代表 2026-08-27 的同步快照，
已经不是当前状态。不要直接根据旧结论执行 fast-forward，也不要在未审查冲突和测试前
把开发分支合入 `main`。

## 3. 已完成能力

### M0：Gmsh 运行时与合同

- 固定并探测 Gmsh Python API 4.15.2。
- 支持受控 runtime staging、严格缺依赖失败报告。
- Device IR v2、v1→v2 迁移及 `vela.mesh_field_bundle.v1` schema 已落地。

### M1：最小跨格式闭环

- 公开默认路径为 `offline-static + gmsh`。
- 从同一规范化 bundle 输出 `device.msh`、`device.vtu`、可选 `device.vtk`、
  `device.grd/device.dat`、`mesh.json`、`doping.csv`。
- MSH Physical Groups 和五个 NodeData 掺杂字段均支持独立回读验证。

### M2：SDE 核心语义

- 支持逐 body 的 ABA/BAB 几何语义。
- 支持二维直线 `insert-vertex`、完整边段 contact、创建态 ownership。
- 支持 `NoReplace`、`LocalReplace`、`Replace` 及 line-baseline Gaussian。
- Sentaurus 2022 oracle 已冻结共享节点、placement 顺序和二维 Positive 方向语义。

### M3：refinement

- 支持 Line、Rectangle、Polygon Ref/Eval window。
- 支持 region、material、window placement，以及 MaxLenInt、DoubleSide、
  UseRegionNames。
- 各向异性和重叠约束以可审计 Gmsh field DAG 表达。

### M4：IV/BV 和 DF-ISE/TDX 闭环

- 单命令能够把 SDE/SDevice 输入传递到预处理、Vela deck、求解曲线和验收报告。
- 内部 TDX T-2022.03-SP2 已成功将 GRD/DAT 转换为 TDR，拓扑、接触及节点场严格回读。
- 完整 PN2D IV 已达到 `accepted`。

### M5：Transfinite Rectangle 与最终 BV 验收

- 支持 `auto | transfinite-rect | unstructured`。
- 支持 `Left | Right | AlternateLeft | AlternateRight` 三角形对角布置。
- PN2D `auto` 选择 Transfinite Rectangle，得到 1482 个节点、2800 个 Tri3；每个
  三角形包含直角且满足非钝角门槛。
- PN2D BV 使用原子化
  `element_edge_sg_gss_laux + element_vertex_box_measure + bernoulli + mixed_voronoi`
  组合，406 个收敛点覆盖至 -20 V。
- 端子电流混合容差最大比值为 0.966，小于合同门槛 1.0。
- 与冻结 Sentaurus BV 曲线比较 213 个点，最大差异 0.0561 decade、最大相对误差
  13.78%，通过 0.08 decade / 20% 门槛。
- M5 TDX/TDR live gate：1482 节点、2800 个体 Tri3、7 个节点数据集均零差异。

权威完成记录：

- `docs/superpowers/plans/2026-08-24-sde-gmsh-offline-mesh-and-field-export.md`
- `docs/validation/sentaurus_gmsh_m5_acceptance_2026-08-26.md`

## 4. 已有验证证据

- Windows MSYS2 UCRT64 Release：最终 M5/主线同步树全量 CTest 712/712 通过。
- 提交 `69b259b` 后再次执行 Windows Release 全量 CTest：712/712 通过。
- 提交 `61bac71` 后 Windows 定向 Gmsh 测试：6/6 通过。
- 提交 `61bac71` 的 WSL Ubuntu 24.04 Debug 构建成功，全量 CTest 697/697 通过。
- GitHub-hosted Ubuntu workflow 已配置 Gmsh 4.15.2 安装和 probe，但尚无本交接可确认的
  GitHub Actions 实际运行记录。

上述证据属于当时提交树。同步当前 `main` 后必须重新构建和测试，不能直接把旧结果视为
新合并树的验收结果。

## 5. 当前已知问题

### P0：开发分支与最新 main 已分叉

`main...HEAD` 当前为 `22 19`。`main` 的新提交包括 AC 计划、Genius/Sentaurus BJT
参考、BJT 电流恢复与 TCAD fixture 整理，可能与 reference、测试注册、文档和 CI 路径
发生交叉。应先审阅 `git log --left-right --cherry-pick main...HEAD` 和三点差异，再选择
merge 或 rebase；保留清晰审计历史时优先普通 merge。

### P1：oracle 文本被 `.gitattributes` 隐藏为二进制差异

`reference_tcad/pn2d_sentaurus2022_oracle/.gitattributes` 当前对 `artifacts/**` 和
`input/**` 同时设置 `-diff -text`，导致 `.cmd`、`.scm`、`.txt` 和 `.log.eval` 也不能
进行文本审查。建议保留 `-text` 以避免 Git 自动换行改变受哈希保护的原始字节，同时只
对真正的 `.tdr` 保持 `-diff`，并对已确认的文本后缀显式设置 `diff`。不要执行
`git add --renormalize`，否则可能破坏 manifest 中的原始 SHA-256。

### P1：oracle 原始日志的公开发布策略未决定

提交 `2a6e9a7` 中的 evaluated log、stdout、stderr、版本及退出码属于可追溯 oracle，
不是普通构建日志；但原始 stdout/version 可能含商业软件声明、主机名、构建主机或进程
信息。若仓库将公开，应由维护者决定：内部原始、哈希封印的 oracle 与公开脱敏结构化
fixture 分离。不要直接改写原始文件后继续宣称其为原始 oracle。

### P1：CI YAML 行尾导致默认 `git diff --check` 失败

当前 `git diff --check main...HEAD` 仅报告 `.github/workflows/ci.yml` 新增行带 CRLF。
`git -c core.whitespace=cr-at-eol diff --check main...HEAD` 可通过，但建议在根
`.gitattributes` 为 YAML 明确 LF，并只规范化该 workflow。必须避免全仓机械换行，尤其
不能改动受哈希保护的 oracle。

### P1：旧验收文档存在阶段性结论歧义

`docs/sentaurus2022_e2e_validation_report.md`、
`docs/validation/sentaurus2022_remaining_tasks_20260822.json` 和 M4 验收记录含 BV
尚未通过的历史状态。M5 文档是当前结论。建议给旧文件增加“历史快照、已被 M5 验收
取代”的显著标识，而不是删除失败证据。

### P1：最新合并树尚未执行最终双平台回归

当前缺少“合入 22 个最新 main 提交后”的 Windows UCRT64 Release 全量 CTest、WSL
Ubuntu 24.04 构建/CTest、Gmsh probe/M1–M5 定向测试以及 GitHub Actions 结果。

## 6. 明确的 fail-closed 边界

以下能力仍未实现或未完成等价性资格认证，遇到时必须明确失败，不能静默近似：

- AB、ABiA、ABiB、XX 布尔模式。
- Polygon body、3D、ACIS、process emulation。
- body contact、`remove`、三维 `imprint`。
- MaxTransDiff、MaxGradient、MaxInterval、multibox、offsetting。
- 非零 constant/analytical restriction decay 的 Gmsh 等价实现。
- 完整 Elk Scheme 及任意控制流求值。
- 通用 Newton carrier/global continuity `enforce` 策略；PN2D 当前使用 report-only
  求解并由独立强制 acceptance 合同验收。

## 7. 建议的新对话执行顺序

1. 在开发 worktree 确认分支和干净状态，不触碰其他 worktree 的用户修改。
2. 审阅最新 `main` 的 22 个独有提交与开发分支 19 个独有提交，识别 reference、CTest、
   workflow、schema 和 solver profile 的交叉点。
3. 在开发分支普通 merge 最新 `main`；解决冲突时保持 PN2D BV 原子 profile 整体移动，
   不拆分 SG/GSS-Laux、element-vertex box source、Bernoulli midpoint、mixed-Voronoi 和
   non-obtuse qualification。
4. 修复 oracle `.gitattributes` 文本可审查性、CI YAML LF 规则及旧文档状态标识；对
   oracle 发布策略先形成明确决定，再改动原始 fixture。
5. 使用 MSYS2 UCRT64 Release 增量构建，先运行 Gmsh probe、M1–M5、PN2D IV/BV 和
   Sentaurus E2E 定向测试，再运行全量 CTest。
6. 使用 WSL Ubuntu 24.04 重新配置、构建并运行全量 CTest。
7. 运行严格 JSON/schema、`git diff --check`、`git status`，确认无生成产物、临时日志或
   商业工具输出被意外提交。
8. 若 GitHub Actions 可用，推送候选分支并取得 Ubuntu workflow 的真实通过证据。
9. 更新计划和最终验收报告中的 HEAD、分叉计数、测试数量与耗时，再由维护者决定合入
   `main` 的方式。

## 8. 建议验证命令

PowerShell 中先初始化 UCRT64：

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
Set-Location "D:\code-repo\copilot-worktrees\vela-tcad-dev"

git status --short --branch
git rev-list --left-right --count main...HEAD
git log --oneline --left-right --cherry-pick main...HEAD
git diff --stat main...HEAD

cmake --preset windows-ucrt64-release
cmake --build --preset windows-ucrt64-release --parallel
ctest --test-dir build-release --output-on-failure `
  -R "sentaurus_(device_e2e_contract|gmsh_m[1-5]|sample_integration)"
ctest --test-dir build-release --output-on-failure

git diff --check main...HEAD
git status --short --branch
```

WSL Ubuntu 24.04 应使用独立构建目录，避免混用 Windows 产物：

```bash
cmake -S . -B build-wsl-ubuntu2404 -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-wsl-ubuntu2404 --parallel
ctest --test-dir build-wsl-ubuntu2404 --output-on-failure
```

## 9. 可直接交给新对话的任务提示词

请在 `D:\code-repo\copilot-worktrees\vela-tcad-dev` 的
`codex/sentaurus2022-sde-workflow` 分支继续收尾。先完整阅读
`docs/validation/sentaurus_gmsh_handoff_2026-09-05.md`、M5 验收记录和开发计划。当前
开发分支与本地 `main` 的实测分叉为 `22 19`，不得沿用旧文档中的“落后 0”结论。
请先审阅双方提交和三点差异，再把最新 `main` 以可审计方式同步到开发分支；保留
PN2D BV 的原子 `element_edge_sg_gss_laux` profile。随后处理 oracle 文本 diff 属性、
CI YAML LF、历史验收文档标识和 oracle 公开发布策略，且不得对哈希封印 oracle 执行
批量 renormalize。完成后依次运行 Windows UCRT64 Release 定向及全量 CTest、WSL
Ubuntu 24.04 全量测试、Gmsh probe/M1–M5/PN2D IV-BV 合同和 `git diff --check`。
更新最终验收文档，提交修改，并报告提交哈希、测试结果、剩余风险和安全合入 main 的
命令。不要触碰其他 worktree 的未提交内容，也不要提交 build、临时日志或未审查的商业
工具输出。
