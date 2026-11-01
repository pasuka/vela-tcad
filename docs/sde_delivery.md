# SDE 当前交付与验收

核对日期：2026-11-01。建模与交付主范围已完成；验收源码为
`72a5b1294ab29bdc55731f4920329cbcee0093ac`。本次文档更新不改变该版本身份。

## 功能与 CI

八份原始 SDE 脚本（PN2D、Schottky、gate_reference、两份 BJT、SingleDevice、
TransportModels、BVmethods）均支持离线建模，输出同源的 `device.grd`、`device.dat`、
`device.msh`。DAT 的必交付为掺杂场。TDR 与原始日志仅用于内部验证。

[CI 运行 41193813343](https://github.com/pasuka/vela-tcad/actions/runs/41193813343)
于 2026-10-31 成功完成，单一 [ci.yml](../.github/workflows/ci.yml) 中两个 job 均运行于 Ubuntu 26.04：

| Job | 实际结果 |
| --- | --- |
| build-and-test | Debug 构建成功；CTest 选择 1190 项，1186 通过、0 失败、4 条件跳过 |
| SDE modeling qualification | Release 运行程序构建、离线包组装和迁移成功；八例各两次通过；专项 CTest 18/18 通过 |

全量 job 的四个跳过项是 `sdevice_e2e_pn2d_0v`、`sdevice_e2e_pn2d_iv`、
`sdevice_e2e_pn2d_bv`、`sdevice_e2e_gates`，不能记为已执行通过。
SDE 专项另外执行了自产网格的原生导入及八份命令检查；这些不等于完整 IV/BV 数值资格。
完整统计与逐次输出哈希见 [main 结构化验收](validation/reference_sde_main_acceptance_2026-11-01.json)。

## 当前 Ubuntu 离线包

| 项目 | 核对值 |
| --- | --- |
| 文件 | `vela-sde-ubuntu26.04.tar.gz` |
| 平台 | Ubuntu 26.04 x86-64 |
| 字节数 | 34185184 |
| SHA-256 | `33a0b2540c64890c2a61a1247e55afcf521196eaab671aa3a6a498cc642b16e6` |
| Python | 隔离 CPython 3.12.12 |
| Gmsh | 外部显式提供 4.15.2，包内不捆绑 |
| 源码状态 | 上述 main 提交，打包时工作树干净 |
| 本次归档核对 | GitHub ZIP 摘要、包摘要及 sidecar、内嵌 manifest、4996 项载荷哈希全部一致 |

GitHub artifact 名为 `sde-qualification-ubuntu26.04`，ID 为 `13040705264`。
ZIP SHA-256 为 `0e6b0a0ee5533885a3b104da9c7dec7a53892f10913fbe73555fcaa81e248ce7`。
CI 保留期为 14 天，该产物到期日为 2026-11-14；本次已将原 ZIP 及展开内容保存至 ignored
`build-release/sde-delivery/2026-11-01-main-72a5b129/`，独立于 GitHub 保留期保存。

本地目录包含原 ZIP、`packages/` 下的包与 manifest/SHA-256 sidecar、
`package-corpus/qualification.json`、`baseline-corpus/case_summary.json`，以及
本轮生成的 `ARCHIVE_INVENTORY.json` 文件清单。二进制、完整 manifest、原始路径及运行输出
不入库，公开仓库仅保存结构化摘要与哈希。本轮不创建 GitHub Release 或远程推送。

在 GitHub 保留期内，可从运行页下载，或下载到新的本地目录后校验：

```bash
gh run download 41193813343 --repo pasuka/vela-tcad \
  --name sde-qualification-ubuntu26.04 --dir build-release/sde-delivery/from-ci-72a5b129
cd build-release/sde-delivery/from-ci-72a5b129/packages
sha256sum --check vela-sde-ubuntu26.04.tar.gz.sha256
```

包的使用和显式 Gmsh 路径设置见 [P6 使用说明](sde_p6_qualification.md#build-and-supply-the-package)。

## 证据边界与后续工作

历史 Windows/Ubuntu 包、32 次商业网格对照及 8 次 TDX 回转仍以 `e1b83fd3` 和
[2026-10-30 验收](validation/reference_sde_p6_acceptance_2026-10-30.json) 为准，原记录保持不变。
当前 CI 没有调用商业 SDE/TDX，也没有生成当前 main 的 Windows 包；不将历史资格换标为当前版本。

P1–P6 建模与交付主范围已收尾。以下项目仍需另外确定范围和验收合同：

- SimpleMOS/LDMOS 原始 SDE、Polygon body、其他布尔模式、非零 decay。
- 自产网格的完整 IV/BV、NonLocalPath、外电阻 R2/R3 数值资格。
- 原 MOS 命令的 PolySilicon 参数与 Thermode/System 前端支持，保留现有明确拒绝边界。

3D、ACIS、完整 Elk Scheme 仍在本轮承诺之外。需要当前 main 的 Windows 正式交付时，
应另外组装并验证该版本，不影响目前 Ubuntu CI 和八例建模的完成状态。
