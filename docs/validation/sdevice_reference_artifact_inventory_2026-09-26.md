# SDevice cmd 直跑：参考产物盘点与归档提案（2026-09-26）

状态：**只读盘点与首批归档已完成（A-1/A-2 已确认）；重跑未执行（A-3 延后）**。
本文是 [开发计划](../superpowers/plans/2026-09-26-sdevice-cmd-direct-run-development-plan.md)
第 12 节（产物管理）与 WP0 卡片“VM 目录 + 首份 manifest”的执行前输入。

## 1. 盘点方法

- 本地：`reference_tcad/**`、`reference_staging/`、`build*/reference_tcad/` 搜索
  `*.tdr *.plt *.grd *.dat *.sav`。
- 虚拟机：SSH 别名 `sentaurus`（`~/.ssh/config`：192.168.119.131，登录用户 **root**，
  `HOME=/root`；[工作流文档](../sentaurus_vm_ssh_workflow.md) 写的 `tcad` 是数据目录的
  拥有者名，不是登录用户）。只读命令：`find`/`stat`/`sha256sum`/`head`。
  注意 MSYS2 的 `ssh` 不读 Windows 用户目录的 `~/.ssh/config`，需用
  `C:\Windows\System32\OpenSSH\ssh.exe`。
- 判定“deck 与仓库一致”：VM 上运行时使用的 cmd 文件与 `reference_tcad/<case>/source/`
  中同名文件 SHA-256 相同（前 16 位见下表）。
- Sentaurus 版本：`/atctools/Synopsys/tcad/T-2022.03`，日志头 `T-2022.03-SP2 (0.7745337)`。
  磁盘剩余 111 GB。sdevice 每次启动的 license 签出约 3 分钟（`pn2d_0v`：wallclock 192 s，
  CPU 1.1 s）。

## 2. 本地结论

仓库与本机构建目录中**没有任何** Sentaurus 原始产物；`reference_staging/`、
`build-release/`、`build-python/` 在当前检出中不存在。仓库只含 A 类源文件（cmd/par）、
B 类曲线 CSV 与派生 JSON/CSV/PNG（`templates_ldmos_sentaurus2022` 另含 6 个 Vela `.h5` 种子）。

## 3. 虚拟机逐案例结果

VM 根目录：`/root/sentaurus_runs/vela_oracle_2022`（记 R）、
`/home/tcad/sentaurus_runs/vela_oracle_2022`（记 H）、`/root/sentaurus_runs/vela_sde_fixtures`。

| 案例 / deck | VM 目录 | 网格 | `_des.tdr` | `.plt` | `_des.log` | deck 一致 | 结论 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| pn2d `pn2d_0v_sdevice.cmd` | `R/sentaurus2022_license_smoke/source` | `pn2d_msh.tdr`（61 021 B） | `pn2d_0v_des.tdr` | `pn2d_0v.plt` | `pn2d_0v.log_des.log` | ✓ cmd `0e699077`、`pn2d_sde.cmd` `407d1326`、`pn2d_msh.cmd` `beb6cb22`、`models.par` `b4b3ebfd` | 可直接归档 |
| pn2d `pn2d_iv_sdevice.cmd` | — | 同上网格 | ✗ | ✗ | ✗ | — | 需重跑（deck 内 `Plot(Time=Intervals=200)` → 201 个 TDR，约 60 MB） |
| pn2d `pn2d_bv_sdevice.cmd` | — | 同上网格 | ✗ | ✗ | ✗ | — | 需重跑（`Intervals=400` → 401 个 TDR，>100 MB） |
| pn2d variants P0–P3/D1/D2/G1/G2（`forward/reverse_sdevice.cmd`） | `H/pn2d_variants_8ff9b76aac2e/<V>/<M>` | ✓ | ✓ 多偏置 | ✓ | ✓ | ✓（P0/original：forward `6acee155`、reverse `98f018f9`） | 备用，非主验收 deck |
| schottky `schottky_n_des.cmd` | `R/schottky_t2022_refresh_20260818` | `schottky_n_msh.tdr` | `schottky_n_des.tdr` | `forward_schottky_n_iv_des.plt`（+ 初始 `schottky_n_iv_des.plt`） | `schottky_n_des.log` | ✓ des `2d27fafd`、sde `c381e21d` | 可直接归档 |
| genius_bjt `bjt_m0_des.cmd` | `R/genius_bjt_wp0_wp2_20260902/run_v2` | `bjt_msh.tdr`（coarse） | `bjt_m0_des.tdr` + `bjt_m0_{equilibrium,vbe070}_des.tdr` + `bjt_m0_vce_000N_des.tdr` | `bjt_m0_{base,collector}_bjt_m0_des.plt` | `bjt_m0_des.log` | ✓ `417f2914`；sde `2cd887d2` | 可归档 |
| genius_bjt `bjt_m1_des.cmd` | `R/genius_bjt_refined_multibias_20260904/source`（同 `H/genius_bjt_mesh_sensitivity_20260903/local_refined/source`） | `bjt_msh.tdr`（local_refined） | 同上模式 | 同上模式 | `bjt_m1_des.log` | ✓ `1cfb8cbe`；`bjt_sde_local_refined.cmd` `1c3b0f0d`（run_v2 内的 m1 是旧版 `40d4321a`，勿用） | 可归档 |
| bvmethods 6 个 deck | `R/bvmethods_t2022_refresh_20260818`（Workbench 节点 n3–n8） | `n1_msh.tdr` | `n3..n8_des.tdr` | `n3..n8_des.plt` | `n3..n8_des.log` | ✗ **物理/设置不同**：预处理后的 `pp3–pp8_des.cmd` 与仓库整理版比，aba_poisson `BreakAtIonIntegral(1 1.1)` vs 仓库 `(3 1.)`（击穿判据不同 → 扫描终点不同）、Plot 字段名不同；external_resistor `Band2Band(Model=NonLocalPath)`/`Avalanche(GradQuasiFermi)`/`AbsVal=1.3e-3` vs 仓库 `E2`/`Eparallel`/`1.443e-3`（隧穿模型、雪崩驱动力、停止电流三处物理/判据不同）。参数文件 `pp*_des.par` 与 `models.par` 只有格式差异（数值相同）。网格：VM `n1_msh.tdr` 由 `pp1_dvs.cmd` 生成，仓库 `bvmethods_nmos_sde.cmd` 是同一器件的**坐标轴转置重写**（VM x=深度、仓库 y=深度，`PolySilicon`→`PolySi`，`Tox 20e-4`→`0.002`），几何等价但节点坐标不同 | VM 结果是“另一组设置”的结果，不能作为仓库 deck 的验收基准；若要 deck 一致的参考，需按仓库 deck 重跑（含用 `bvmethods_nmos_sde.cmd` 重建网格） |
| bvmethods `full_physics_ablation_sdevice.cmd` | — | — | ✗ | ✗ | ✗ | 仓库文件是模板（`SRH_MODEL`/`ENORMAL_MODEL`/`SRH_PARAMETER_FILE` 占位符） | 需先物化 4 个 deck；用户决定暂缓 |
| singledevice | `R/singledevice_t2022_refresh_20260818` | `n2_msh.tdr` | `singledevice_des.tdr`、`n2_{Lin,Sat}_des.tdr` | `IdVgs{Lin,Sat}_singledevice.plt` | `singledevice.log_des.log` | ✗（形式）：VM 版是手工替换占位符后的 CRLF 文件（`@tdr@`→`n2_msh.tdr`、`n@node@`→`n2` 等），另有 Plot 名 `ImpactIonization`→仓库 `AvalancheGeneration`；物理与 Solve 序列相同 | 仓库版本本身会被 E101/E102 拒绝，只作“预处理指令应被拒绝”夹具；VM 结果仅供人工参考 |
| transportmodels pp6/pp7/pp12/pp13 | `R/transportmodels_dd_dg_t2022_20260818_run02` | `n1_msh.tdr` | `n6/n7/n12/n13_des.tdr` | `IdVgs_n6/n7_des.plt`、`IdVds_n12/n13_des.plt` | `n6..n13_des.log` | ✓（归一化后）：VM 版 2 301 B vs 仓库 2 265 B 的差异**只是行尾空白与末尾空行**（`pp6` 逐行 diff 验证），剥除空白后相同 | 可作参考；当前因 `Thermode{}`、`System(...)` 为 fail-closed 夹具，无需重跑 |
| templates_ldmos | `R/templates_ldmos/phase01_original_20260826_02/work` | `n1_fps.tdr`（Sentaurus Process 输出，9.5 MB） | `n2/n4/n6_des.tdr` | `IdVg_n2_des.plt`、`IdVd_Vg{1,2}_n4_des.plt`、`n6_des.plt` | `n2/n4/n6_des.log` | 仓库**没有**任何 sdevice cmd（只有 JSON 契约与 `.h5` 种子） | fail-closed 夹具需先把 `IdVg.cmd/IdVd.cmd/BVdss.cmd` 从 VM 拷回 |

附带事实：

- TDR 文件内含时间戳：同一 `pn2d_sde.cmd`/`pn2d_msh.cmd` 三次运行得到三个不同哈希的
  `pn2d_msh.tdr`（`920b55f7`/`23af0444`/`1012c7bc`，字节数相同）。因此 manifest 的哈希只对
  “归档的那份拷贝”成立，不能作为可复现性证据；可复现性由“deck 哈希 + 版本 + 结构摘要”表达。
- **`File{Parameter}` 缺失是常态**：`pn2d_{0v,iv,bv}`、`schottky_n_des.cmd`、`bjt_m0/m1_des.cmd`
  都没有 `Parameter=` 行，Sentaurus 使用内置默认参数运行（手册第 2 章“Default Parameters”）；
  运行目录里的 `models.par`（107 158 B，与仓库 `pn2d_sentaurus2018/source/models.par` 同哈希）是
  `sdevice -P` 导出的默认参数副本。按计划决定 4（缺 `Parameter` → 致命 `E301`），这 5 个主验收
  deck 的**原文**都会在校验阶段被拒绝；验收时需使用追加了
  `Parameter="models.par"` 的夹具副本。见第 6 节待决项。

## 4. 用户决定（2026-09-26）

1. 重跑范围：**暂不重跑**（`pn2d_iv`/`pn2d_bv`/bvmethods 留待后续批次）；只归档现有产物。
2. `full_physics_ablation`：暂缓，不物化、不运行。
3. VM 归档写入：先只读盘点，布局确认后再复制（**已于同日确认并执行**，见第 5 节）。

## 5. 归档布局（A-1 已确认并执行，2026-09-26）

按计划第 12.2 节（决定 18，方案 S3）。根目录 `/home/tcad/vela_reference_artifacts/`，
每案例一个版本子目录，`cp -p` 复制（不移动、不删除原目录），登记后 `chmod -R a-w`。
类别目录在前、变体在后（与 C7 schema 的 `path` 模式 `^(source|mesh|plt|tdr|log)/` 一致）：

```text
/home/tcad/vela_reference_artifacts/<case>/v1/
  source/[<variant>/]   A 类：产生结果的 cmd/.par（VM 上实际使用的副本）
  mesh/[<variant>/]     C 类：*_msh.tdr、*_bnd.tdr、sde/snmesh 日志
  plt/[<variant>/]      B 类
  tdr/[<variant>/]      D 类：*_des.tdr、多偏置序列、*.sav 状态
  log/[<variant>/]      E 类：*_des.log、Workbench/脚本 *.out
  artifact_manifest.json   ← 与仓库 reference_tcad/<case>/artifact_manifest.json 相同
```

首批归档结果（206 个文件 + 7 份 manifest，113 MB，均为 T-2022.03-SP2 (0.7745337)；
所有 sdevice 日志以 `Good Bye` 结束）：

| 案例 | 文件数 | `deck_matches_repository` | 来源目录 |
| --- | --- | --- | --- |
| `pn2d_sentaurus2018`（仅 0v） | 12 | true | `R/sentaurus2022_license_smoke/source` |
| `schottky_charon_sentaurus2018` | 13 | true | `R/schottky_t2022_refresh_20260818` |
| `genius_bjt_sentaurus2022`（变体 `m0`/`m1`） | 42 | true | `R/genius_bjt_wp0_wp2_20260902/run_v2`（m0）、`R/genius_bjt_refined_multibias_20260904/source`（m1） |
| `transportmodels_sentaurus2022` | 36 | true（归一化后一致） | `R/transportmodels_dd_dg_t2022_20260818_run02` |
| `bvmethods_sentaurus2018` | 47 | **false** | `R/bvmethods_t2022_refresh_20260818` |
| `singledevice_sentaurus2018` | 27 | **false** | `R/singledevice_t2022_refresh_20260818` |
| `templates_ldmos_sentaurus2022` | 29 | **false**（仓库无 deck） | `R/templates_ldmos/phase01_original_20260826_02/work` |

manifest 每条记录：相对名、类别 A–E、VM 副本的 SHA-256 与字节数、生产者（`sde`/`snmesh`/`sdevice`/
`sprocess`/`repository`/`other`，含 cmd 名、cmd 的 LF 归一化哈希、参数文件、运行目录、完成时间）、
`verify.status=registered`；sdevice 日志附 `summary.exit_status`。结构摘要（网格节点/接触、`.plt`
数据组与首末偏置）由 WP9 的 `register` 补充。生成脚本为一次性工具（`build/make_artifact_manifests_v1.py`，
不入库）；WP9 交付 `scripts/reference_artifacts.py` 后以其 `verify` 重新核对这 7 份 manifest。
`bvmethods`、`singledevice`、`templates_ldmos` 标注 `deck_matches_repository: false`，只作夹具或
“部分参考”，不进入第 8 节一致性验收。

“deck 一致”的判定规则与计划第 12.4 节第 2 条相同：剥除注释、折叠空白后逐 token 比较，而不是
字节相等（TDR/日志的哈希才按字节）。按此规则：pn2d_0v、schottky、genius_bjt m0/m1、
transportmodels 一致；bvmethods（物理与判据不同、网格坐标转置）与 singledevice（占位符替换、
Plot 名）不一致；templates_ldmos 仓库无 deck 可比。

## 6. 待决项（登记到计划附录 A）

| 编号 | 事项 | 影响 | 建议 |
| --- | --- | --- | --- |
| A-1 | 归档布局（第 5 节）与是否复制 `*.out`/`*.log`（含主机名等标识） | WP0 VM 部分、WP9 manifest | **已确认并执行（2026-09-26）**：按第 5 节布局复制 7 个案例，deck 不一致者标 `deck_matches_repository: false`；日志只存 VM，仓库 manifest 只存哈希 |
| A-2 | `File{Parameter}` 缺失的处理：现决定 4 为致命 `E301`，但 5/6 个主验收 deck 原文无 `Parameter` | WP1 golden、WP5 门 4、WP10 验收 | **已确认（2026-09-26）**：保持 `E301`；验收夹具用追加了 `Parameter="models.par"` 的副本（物理不变：该文件就是 Sentaurus 默认参数导出） |
| A-3 | `pn2d_iv`/`pn2d_bv`/bvmethods 重跑时机 | WP10 曲线与场量验收、决定 20 的符号因子（预期取自 `pn2d_iv.plt`） | 在 WP5 可运行 `--check` 之后、WP10 之前的一个批次内完成（约 4 次 sdevice + 1 次 sde/snmesh） |
| A-4 | `.plt` 多段扫描的 `time` 语义（每段重置或累加）需含两个 `Quasistationary` 的参考 `.plt` | C6 输出契约 | 现成证据：`bjt_m0/m1_{base,collector}_*.plt`（两段各自文件，因 `NewCurrentPrefix` 切换）与 transportmodels `IdVgs_n6_des.plt`（`init_` 段被 `System("rm ...")` 删除）；在 WP8 读取这些文件确定 |

## 7. 证据文件

本次盘点的原始命令输出保存在本机 `build/vm_inventory{2,3,4,5}.txt`（已忽略目录，不入库）。
关键哈希已写入第 3 节表格；deck 的完整 SHA-256 以 `tests/fixtures/sdevice/<case>/*.ir.json`
的 `source.sha256`（LF 归一化后）为准。
