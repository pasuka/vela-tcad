# SDevice 前端交付入口与支持矩阵

更新：2026-10-23。本文是本轮具名发行包的用户入口；详细字段见
[工作流](sentaurus_sdevice_workflow.md)，构建发行包见[打包指南](sdevice_packaging.md)。
运行资格以[本轮 D1–D3 账本](validation/sdevice_delivery_revalidation_2026-10-22.json)
和 [D4 两机验证](validation/sdevice_delivery_d4_2026-10-23.md)为准。

## 选择发行包

运行时快照：`6a56593bf3917cc213a00aeab4f4a614befa52ca`。版本号同为 `0.1.0` 的旧包不能仅凭
文件名区分，交付时应核对以下 SHA-256 及包内 `PACKAGE_MANIFEST.json` 的 `commit`。

| 平台 | 包文件 | SHA-256 |
|---|---|---|
| Windows x86_64 | `vela-0.1.0-win64.zip` | `1c7835de896c5cb3db01932a9fce7eab8e96de6f1047bc4b255780057dc353d2` |
| Ubuntu 26.04 x86_64 | `vela-0.1.0-ubuntu26.04-x86_64.tar.gz` | `c57ddc0c10de3986b8676c71598948a20d04507447246616fabb87a867381a23` |

包内带 Python、前端 pyz、runner、网格 importer、材料、schema、模板、动态依赖和许可。
使用包不需要编译器、CMake、宿主 Python、pybind11 或 Sentaurus 安装；原生 TDR/材料输入需由用户提供。
Linux 仍依赖目标系统的 glibc/加载器；不能把 Ubuntu 包视作任意 Linux 通用包。
全部动态 ELF 使用各自 `$ORIGIN` 相对 RUNPATH，不需配置 `LD_LIBRARY_PATH`。
Linux 应在原生、大小写敏感且支持符号链接的文件系统解压，不能先在 Windows 文件系统展开再复制。

本轮文档与中性示例另以 `vela-sdevice-delivery-20261023.zip` 提供，内含文件哈希清单。
它配套上述两个原包，不改动其二进制或包哈希。原包内早期文档的完成状态以本指南及绑定账本更新；
仓库中示例源为 [examples/sdevice](../examples/sdevice)。商用参考数据不随交付说明包分发。

## 从解包到运行

将运行包和配套说明包放在同一空目录。以下 PowerShell 命令从该目录执行：

```powershell
Get-FileHash .\vela-0.1.0-win64.zip -Algorithm SHA256
Expand-Archive .\vela-0.1.0-win64.zip -DestinationPath .\runtime
Expand-Archive .\vela-sdevice-delivery-20261023.zip -DestinationPath .\delivery
$vela = (Resolve-Path .\runtime\vela\vela.exe).Path
& $vela device -h
& $vela device -v
Copy-Item .\delivery\examples\sdevice .\example-run -Recurse
Push-Location .\example-run
& $vela device -P
& $vela device --check .\device_des.cmd
& $vela device --threads 1 --exit-on-failure .\device_des.cmd
if ($LASTEXITCODE -ne 0) { throw 'Vela run failed; inspect run_des.log' }
Get-Content .\device_des_vela\run_manifest.json
Pop-Location
```

Ubuntu 26.04 对应命令（使用随包 Python 解开说明包）：

```bash
sha256sum vela-0.1.0-ubuntu26.04-x86_64.tar.gz
mkdir runtime delivery
tar -xzf vela-0.1.0-ubuntu26.04-x86_64.tar.gz -C runtime
./runtime/vela/bin/python/bin/python3 -I -B -m zipfile -e vela-sdevice-delivery-20261023.zip delivery
vela="$(pwd)/runtime/vela/vela"
"$vela" device -h
"$vela" device -v
cp -R delivery/examples/sdevice example-run
cd example-run
"$vela" device -P
"$vela" device --check device_des.cmd
"$vela" device --threads 1 --exit-on-failure device_des.cmd
cat device_des_vela/run_manifest.json
```

先将显示的归档哈希与表格核对，再解压。每次验证用新目录；`-P` 在当前目录写 `models.par`，
会覆盖同名文件，只用于此新建中性示例。实际器件应在 `File{Parameter="..."}` 中声明自己的 `.par`。
所有 `File{}` 相对路径以 cmd 所在目录解析；sidecar 为 `device_des.vela.json`。
`--allow-lossy` 只在逐项理解材料报告后使用，不能解除不支持模型的拒绝。

`--check` 既校验也发射 JSON，但不求解。检查 `device_des_vela/` 中的 `stage_plan.json`、
`simulation_stage00.json`、`simulation_stage01.json`、`materials_report.json` 和网格资格文件。
不需要另加不存在的 `--emit-json` 选项。需要调试 JSON 时保留整个生成目录及依赖；阶段链由
`vela device` 编排，不能把任意一个中间 JSON 当作独立初始态。

短算例是 1 µm 方形、4 节点/2 三角形、均匀施主 `1e15 cm⁻³`、300 K，anode 扫描至 1 V。
它用于入口演示，不代表 PN 结物理验证。输出含 `current_des.plt`、`final_des.vtu`、
`final_des.pvd`、分步 VTU 和 `run_des.log`；电流为 A/µm。必须同时检查退出码和
`run_manifest.json` 各阶段的 `goal_reached`。使用 `--exit-on-failure` 可使未达目标返回失败；
不启用该选项时退出 0 不足以证明所有目标达到。错误定位见工作流第 2、4 节。

## 具名输入的资格上限

精确原文、LF 哈希、模型及当前资格见[机器可读支持矩阵](sdevice_support_matrix.json)。
`reference_tcad` 的 43 个 cmd 包含变体和消融，不能根据文件存在推断全部支持。

| 输入/版本 | 资格 | 边界 |
|---|---|---|
| PN2D `pn2d_0v_sdevice.cmd`，补 Parameter 的 v2 | 两平台完整 cmd 数值 E2E | 具名材料/网格/输入；原文缺声明会报 E301 |
| PN2D `pn2d_iv_sdevice.cmd`，v2 | 两平台完整 cmd 数值 E2E | 到 10 V；原曲线验证窗口仅 0.2–0.3 V |
| PN2D `pn2d_bv_sdevice.cmd`，v2＋冻结 sidecar/profile | 两平台完整 cmd 数值 E2E | 到 −20 V、十个场量锚点；资格不推广到任意网格 |
| BVmethods `aba_poisson_sdevice.cmd` | 解析/发射/真实输入 `--check` | 未承诺完整 ABA 数值对齐 |
| BVmethods `aba_coupled_sdevice.cmd` | 真实输入 `--check` | 未承诺完整轨迹 |
| BVmethods `external_resistor_sdevice.cmd`（原 fixture E2/Eparallel） | 解析/发射/真实输入 `--check` | 从零纯 QS 未获数值资格 |
| BVmethods `voltage_to_current_sdevice.cmd` | 真实输入 `--check` | 未承诺完整电流控制轨迹 |
| Schottky `schottky_n_des.cmd` | 补 Parameter 后真实输入 `--check` | 独立 JSON IV/重启证据不等于本包 cmd 数值 E2E |
| BJT `bjt_m0_des.cmd` | 补 Parameter 后解析/发射、阶段/时间契约 | 无本轮新包 BJT 数值 E2E |
| `singledevice_sdevice.cmd` | 正确拒绝 | 未展开预处理/占位符 |
| `full_physics_ablation_sdevice.cmd` | 正确拒绝 | 白名单外 Physics |
| transportmodels `pp6_des.cmd` | 正确拒绝 | 超范围结构 |
| T-2022.03-SP2 安装库原版 BVmethods | 独立研发；公开入口拒绝 NonLocalPath | 与原 fixture 的 E2/Eparallel 目标不同 |

PN2D 数值门保持不变，网格 1943 节点/3680 非钝角三角形，单位 µm、cm⁻³、A/µm、V/cm。
−13.2 V 是场量锚点。当前拟合 knee −19.405 V 是契约定义的曲线指标，不等同于任意定义的击穿电压。
Math 的 `Digits`、`RelErrControl`、`ErrRef` 解析并记录 W201，不自动替代 Vela 收敛规则；
接受过程不声称与 Sentaurus 等价。外电阻默认仍为 `nested_scalar`；实验 `coupled_newton` 不是完整资格。

## 平台及证据边界

本轮 D1–D3 已在 Windows/UCRT64 与 Ubuntu 26.04 完成新包 0V/IV/BV、输入 gates、迁移/指纹、
相对 RUNPATH 及原平台范围的并发验收。Linux 完整并发矩阵通过；Windows 原范围只验证锁冲突。
D4 在 T470p 和 C560 上补充原包实际用户入口、短算例和 10 个 fixture 的解析身份检查，
不将此短算例冒充新的全 BV 或商业参考重复计算。

历史 WP 报告及 fixture 的 `reference_artifacts` 文本是当时的来源记录；本指南、支持矩阵和
绑定的新包验收账本共同给出本次状态。扩展 Schottky/BJT 等数值资格以及 NonLocalPath/R2/R3
需按各自计划验收，不能静默升级为已支持。
