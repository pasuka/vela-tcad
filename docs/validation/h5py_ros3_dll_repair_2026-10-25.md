# Windows h5py DLL 导入修复（2026-10-25）

原先由 h5py DLL 导入导致失败的 **8 项 CTest 全部通过**；修复后完整 Release CTest
为 **1113 项：1109 通过、0 失败、4 跳过**，退出码 0。机器可读证据见
[配套 JSON](h5py_ros3_dll_repair_2026-10-25.json)。本轮仅修复 Python 包环境并补充文档，
没有修改 Vela 求解器、物理模型、验收门限或测试跳过规则。

## 4 项跳过测试

用户引用的原完整测试记录为 `build-release/branch_organization_20261025/ctest.xml`，
当时共 1105 项。以下四项的 system-out 均为：

```text
SKIP: VELA_REFERENCE_ARTIFACTS_DIR is not set; native references unavailable.
```

| 测试名称 | 范围 |
|---|---|
| `sdevice_e2e_pn2d_0v` | PN2D 零偏原生参考端到端对照 |
| `sdevice_e2e_pn2d_iv` | PN2D I–V 原生参考端到端对照 |
| `sdevice_e2e_pn2d_bv` | PN2D BV 原生参考端到端对照 |
| `sdevice_e2e_gates` | Schottky / BVMethods 的 check 验收及不支持模型的拒绝检查 |

`scripts/sdevice_acceptance.py` 在未提供参考目录时返回 77，CMake 将 77 注册为
`SKIP_RETURN_CODE`。修复后的完整测试仍保留相同四项跳过；它们没有因 h5py 修复而获得
原生参考数据，也没有被计为已通过的器件数值验收。

## 根因

CTest 实际使用 UCRT64 Python 3.14.9，已安装 NumPy 2.5.4、h5py 3.16.0-2、
HDF5 2.2.0-4。对原 h5py 与 HDF5 执行包完整性检查，均为 0 altered files。
失败是 DLL 已找到但缺少被引用入口，不是单纯缺少 DLL 文件。

逐个核对所有 h5py `.pyd` 的 PE 导入表与实际 DLL 导出符号，发现仅以下 ROS3 符号缺失：

| 扩展模块 | 从 `libhdf5-320.dll` 引用但未导出的符号 |
|---|---|
| `defs` | `H5Pget_fapl_ros3`、`H5Pget_fapl_ros3_token`、`H5Pset_fapl_ros3`、`H5Pset_fapl_ros3_token` |
| `h5fd` | `H5FD_ROS3_id_g` |

当前 HDF5 的 `H5_HAVE_ROS3_VFD` 未定义，而旧 h5py 二进制包含这些 ROS3 引用，
因此导入 `h5py.defs` 时即失败。三个电热相关测试栈中的相对导入异常位于导入回退过程，
最终根因仍是同一个 h5py DLL 错误；本轮没有改写这些导入回退逻辑。

## 修复

使用 [MSYS2 官方 h5py 配方](https://github.com/msys2/MINGW-packages/tree/master/mingw-w64-python-h5py)
及相同版本的 PyPI 3.16.0 源码，按官方配方校验源码和两个 MinGW 补丁的 SHA256。
按 [h5py 的源码构建方式](https://docs.h5py.org/en/stable/build.html#downstream-packagers)
复用当前 Python / NumPy / HDF5，仅补齐缺少的构建依赖。

重建时由 h5py 自动探测现有 HDF5 的功能，得到 `ros3=False`；未重编译或替换 HDF5。
本地包版本标为 `mingw-w64-ucrt-x86_64-python-h5py 3.16.0-2.1`，通过 pacman 安装，
与官方配方的差别仅为本地 pkgrel。原 3.16.0-2 包及签名已备份在 ignored 证据目录。
构建中的 Windows 长路径问题通过临时短路径映射解决，映射已移除。

安装后验证：

- 默认 CTest 解释器正常 `import h5py`，编译与运行 HDF5 版本均为 2.2.0。
- 再次检查所有 `.pyd` 的 DLL 导入符号，缺失列表为空。
- gzip 和 LZF 压缩数组的实际 HDF5 写入 / 回读精确一致。
- h5py 包完整性为 242 文件、0 altered files。
- HDF5 DLL 的 SHA256 修复前后相同，见 JSON；没有更换 C++ HDF5 库。

## 测试结果

先按原失败清单精确选择八项，8/8 通过、0 跳过，用时约 80.41 s：

| 测试 | 结果 |
|---|---|
| `state_archive_python_interop` | 通过 |
| `state_archive_dc_worker` | 通过 |
| `reference_tcad_regression` | 通过 |
| `simplemos_engineering_codespaces` | 通过 |
| `simplemos_engineering_matrix_contract` | 通过 |
| `templates_ldmos_electrothermal_input_regression` | 通过 |
| `templates_ldmos_joint_local_fields_regression` | 通过 |
| `electrothermal_dc_runner_regression` | 通过 |

随后运行 `ctest --preset windows-ucrt64-release --output-on-failure --parallel 4`，
1113 项中 1109 通过、0 失败、上述 4 项跳过，用时约 216.63 s。
较原 1105 项多出的 8 项来自此前新增的 BV 步长控制、sweep 和 CMD profile 测试，
具体名称记录在 JSON 中。

原故障、包备份、缺失符号、源码校验、构建 / 安装日志、读写检查与两轮 CTest XML / 回执
均保存在 `build-release/h5py_dll_repair_20261025/`。历史失败记录没有改写。
