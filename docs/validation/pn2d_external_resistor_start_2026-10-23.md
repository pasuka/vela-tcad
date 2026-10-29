# PN2D 外接电阻法起步验证（2026-10-23）

## 结果与范围

本地 `main`（`179aa4e0`）已合并至 `codex/sdevice-cmd-frontend`，合并提交为
`5ece7f33`，父提交为 `4e777034` 和 `179aa4e0`。合并后的受跟踪文件树与
main 一致。六个重复未跟踪文件先核对内容并备份，再由 main 纳入跟踪；原有
ignored 私有实验和商业输入原位保留。私有候选 zip SHA-256 与交接记录一致。

按用户要求改从 PN2D 开始，电阻经用户确认固定为 **1e7 Ω·μm**。
直接定压和 nested_scalar 完成零偏到 −10 mV 的六点控制；coupled_newton
未达到首个 −2 mV 目标。**这不是完整 BV、原生外电阻对照或 BVmethods R2/R3 收口。**

## 输入与执行

- Windows UCRT64 Release；配置检测到 Eigen 5.0.1、Boost 1.92、HDF5/HighFive、
  UMFPACK。运行环境指定 `VELA_LINEAR_SOLVER=umfpack`，OMP/BLAS 各一线程。
- PN2D：1943 节点、3680 三角形；Anode 扫描、Cathode 接地，300 K。
  网格、掺杂、材料来自本 worktree 原有 `build/iv_bv_atomic_probe/`，复制到新实验目录并绑定哈希。
- 通过当前 `pn2d_bv` v3 模板生成，保留 Masetti/HFS、SRH、OldSlotboom、
  VanOverstraeten 及完整 `element_edge_sg_gss_laux` 原子组合；模板原本无 BTBT。
  不把 BVmethods NMOS 的 Fermi/Lombardi 配置套用到 PN2D。
- 三路分别从自身零偏初始化，无原生映射态或非零重启种子。
  请求外压/直接定压节点为 `0, -0.002, -0.004, -0.006, -0.008, -0.01 V`。
- 外电阻使用 Anode 具名电流，`current_direction=+1`，
  `L = Vinner + 1e7 * IAnode - Vouter`；电流为 A/μm。
  增广采用 `direct_bordered`，设备标量门 `1e-8`，负载门 `1e-6 V`，
  载流子行门 `0.001`、全局闭合门 `0.01` 保持模板/配置原值。

最初外电阻入口禁止递减外压，因而在装配前拒绝 PN2D 负向目标。
本轮将该检查改为允许任一方向的单调序列，仍拒绝途中换向；现有有符号步长
控制和负载方程保持。新增 C++ 回归同时覆盖 nested_scalar 和 coupled_newton。

## 数值账本

| 路径 | 退出码 | 接受点数 | 最后请求坐标 | 最大 KCL 误差 A/μm | 最大重算负载误差 V |
|---|---:|---:|---:|---:|---:|
| 直接定压 | 0 | 6 | −0.01 V | 1.947e-23 | 不适用 |
| nested_scalar | 0 | 6 | −0.01 V | 1.947e-23 | 2.331e-12 |
| coupled_newton | 1 | 2 | −1.953125e-6 V | 7.099e-23 | 1.517e-17 |

直接定压和 nested_scalar 六点均通过生产接受状态和独立 KCL/负载检查。
前两路的最大电子/空穴全局闭合比约为 `8.391e-5 / 1.804e-6`。
这只是低压控制，电阻压降很小，不能建立击穿区负载反馈资格。

增广首次 −2 mV 候选的负载与 KCL 已闭合，但载流子行最大比
`0.00692535 > 0.001`，空穴全局闭合比 `0.0135964 > 0.01`，均不通过。
缩步后只接受约 −1.953 µV；该点的空穴比约 `0.5864`，但 source-floor
资格为 false，故生产规则接受。**不能把它称为严格空穴闭合通过。**

最后拒绝发生在外压约 `−2.953125e-6 V`，内压
`−2.9531249974405186e-6 V`：设备残差 `1.3180182e-11`、负载残差
`−3.19543e-18 V`、载流子行最大比 `7.47914e-5`；空穴全局闭合
已具资格，比例 `0.587142`，由此拒绝。标量残差收敛不能代替物理接受。
后续应在同态诊断为何标量停止先于物理闭合，并验证继续校正策略；不降低原门限。

## 证据和复现

所有路径相对本 worktree。原始证据位于 ignored
`build-release/bvmethods_resume_20261023/`：

- `pre_merge_untracked_*.zip`、`pre_merge_inventory.json`：合并前备份。
- `configure.log`、`build*.log`、`ctest*.log`：环境和构建/回归账本。
- `pn2d_reverse_enabled/`：三路配置、输入哈希、执行摘要、CSV 和运行日志。
- `pn2d_failure_capture/`：仅增广的重复捕获，保存 13 个拒绝态 HDF5；
  两轮接受 CSV 逐字节一致。仅打开拒绝态输出，不改变接受策略。
- `evidence_manifest.json`：实验 JSON/CSV/HDF5/备份的 SHA-256 清单。

入口为 [控制脚本](../../scripts/run_pn2d_external_resistor_control.py)。例如在 worktree 根目录：

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts/run_pn2d_external_resistor_control.py `
  --runner build-release/vela_example_runner.exe `
  --mesh build/iv_bv_atomic_probe/mesh.json `
  --doping build/iv_bv_atomic_probe/doping.csv `
  --materials build/iv_bv_atomic_probe/materials.json `
  --output build-release/pn2d_resistor_fresh --resistance 1e7 --timeout 120
```

输出目录必须不存在；`--arms coupled_newton` 可单独捕获增广失败。
脚本退出 2 表示有一路执行失败或未完成六点物理检查；各子进程真实退出码另存。
此复现依赖本地 ignored 网格/材料，不能仅凭 Git 克隆宣称完整可复现。

## 验证与运行边界

Release 全量初次构建及最终受影响目标构建通过；没有运行全量 CTest。
新增反向扫描 C++ 测试通过，1 例/24 断言；Python 独立物理检查 4 例通过。
最终相关 CTest **7/7 通过**（14.57 秒），包括完整 dc_sweep、输出契约、
coupled_load_line、三个 PN2D 既有配置/原子组合回归及新增物理检查。
Newton 的 arclength/terminal-current/resistor_contract 控制 **4 例/75 断言通过**。
报告相对链接有效；ignored 证据清单 **123 项哈希全部通过**。
首次临时目录权限失败属 Windows 沙箱环境问题；将 TEMP/TMP 指向本 worktree
ignored 目录后重跑。最初 runner 参数遗漏和递减序列拒绝均保留原始失败日志。

T470p、C560 只读连通性检查通过；Vela VM 两次握手关闭，未派发远端作业。
本轮数值控制在本地完成，机器身份及原始连接记录仅保留在 ignored 目录。
没有完整 BV、原生新仿真、严格重启资格或生产默认迁移。
