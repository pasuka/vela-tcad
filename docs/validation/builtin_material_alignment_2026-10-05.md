# WP13 内置 Si/SiO2 材料对齐与基线审计（2026-10-05）

状态：**完成 WP13 原卡片的 Windows 验收**。在 `codex/sdevice-cmd-frontend` 当前工作树上，
以 `9f9d2f5d` 为变更前状态，对齐无 BGN 的 300 K 内置标量，并将历史校准路径依赖的旧默认显式冻结。
Windows Release 最终全量非外部参考 CTest **1026/1026**；三个 PN2D 门全部通过；没有放宽数值容差。
额外发现的 Schottky 历史重启问题经对齐前发行包复现，作为独立待修复项保留。Ubuntu 额外验证未完成，原因见下文。

## 对齐值与边界

参数依据为现有 WP3 默认 `share/vela/materials/models.par`，由 WP3 发射器以 `active_models=[]`、
`temperature_K=300`、`mesh_materials=[Si, SiO2]` 生成
[`default_materials.json`](../../share/vela/materials/default_materials.json)。
`DefaultLibraryTest.test_silicon_matches_the_sentaurus_default_dump` 继续逐字段核对已归档的 Sentaurus 默认转储；
本次不需要新的 VM 文件传输或 sdevice 重跑。原 `models.par` 内容未变。
其 SHA-256 为 `7399cbd4e3df62f1fe3ef4fe7842f88cdb83c6d97c5ff7e399bcd3bac7765983`。

| 材料/标量 | 变更前 | 变更后，300 K 无 BGN | 相对变化 |
| --- | ---: | ---: | ---: |
| Si εᵣ | 11.7 | 11.7 | 0 |
| Si E_g，eV | 1.12 | 1.1241592307692307 | +0.37136% |
| Si χ，eV | 4.05 | 4.05 | 0 |
| Si N_c，cm⁻³ | 2.8e19 | 2.8566791451911786e19 | +2.02426% |
| Si N_v，cm⁻³ | 1.04e19 | 3.104657033455896e19 | +198.52471% |
| Si nᵢ，cm⁻³ | 1e10 | 10750176656.471613 | +7.50177% |
| Si μ_n，cm²/(V·s) | 1350 | 1417 | +4.96296% |
| Si μ_p，cm²/(V·s) | 480 | 470.5 | −1.97917% |
| SiO2 εᵣ / E_g，eV | 3.9 / 9 | 3.9 / 9 | 0 |
| SiO2 χ，eV | 0.95 | 0.9 | −0.05 eV；−5.26316% |

C++ 内置值使用 SI：密度乘 `1e6`、迁移率乘 `1e-4`；`ni` 在构造器中由
`sqrt(Nc*Nv)*exp(-Eg*q/(2*kB*T))` 推导，而非另设独立旧常数。
Catch2 在 SI 与 `unit_scaling` 两种内部单位下，以 **1e-12 相对容差**逐字段核对快照，并验证快照实际经加载器读入后的等价性。

本次仅改内置标量：不激活 BGN、量子、复合或电离模型，不修改 WP3 发射器物理实现和任何求解算法。
`Material::atTemperature` 保持固定带隙、`Nc/Nv ∝ T^(3/2)`、迁移率 `∝ T^(-2.2)`；新增质量作用关系测试覆盖
250/300/350 K，但这不等于其他温度也匹配 `.par` 的 Varshni/有效质量律。无材料文件 deck 与部分覆盖中省略的字段
会使用新默认；显式材料字段继续优先。

`generate_default_models_par.py --check` 同时检查 `.par` 与生成的 JSON，能拒绝被修改的快照。
JSON 随 zipapp 资源及 `SDeviceRuntime` 安装组件提供。当前工作流与配置文档已同步；早期报告保留历史数值并增加后续状态链接。
`SDeviceRuntime` 安装冒烟通过，安装后的 JSON 与源文件字节一致，SHA-256 为
`c52a40007c339d83451366c1bcecbd06e9314b97e1847309df97584f236b9251`。
`.gitattributes` 为默认参数库及快照固定 LF，防止 Windows `core.autocrlf` 改变 `.par` 原始哈希而使生成器门失效。

## 依赖盘点

变更前共有 **25 个直接使用 MaterialDatabase 的 C++ 测试文件**，不是旧计划的 20 个。
新增 `test_builtin_material_alignment.cpp` 后为 26 个。间接通过 DCSweep/Poisson/runner 的依赖也由全量 CTest 覆盖。
文件名下表均相对于 `tests/`；性质断言与明确指定的独立常数没有批量改写。

| 文件 | 分类与处理 |
| --- | --- |
| test_bv_process_probe.cpp | 物理/算子性质；原断言通过 |
| test_cell_reconstructed_avalanche.cpp | 电离源守恒、离散性质；原断言通过 |
| test_dd_gummel.cpp | DD 平衡与输运性质；原断言通过 |
| test_device_stability.cpp | 稳定性与温度性质；原断言通过 |
| test_element_box_transport.cpp | 输运离散性质；原断言通过 |
| test_element_edge_gss_laux_avalanche.cpp | SG/GSS-Laux 离散与守恒；原断言通过 |
| test_external_state.cpp | 状态输入与单位性质；原断言通过 |
| test_fixed_state_operator_audit.cpp | 固定状态算子性质；原断言通过 |
| test_gummel_high_doping.cpp | 高掺杂正性与稳定性；原断言通过 |
| test_ialmob_transport.cpp | 界面迁移率/输运性质；原断言通过 |
| test_impact_ionization.cpp | 物理/VTK 数值一致性；一处隐含旧 ni 改为读取数据库 |
| test_material_local_poisson_charge.cpp | 局部材料 Poisson 电荷性质；原断言通过 |
| test_mesh.cpp | 材料与单位性质；内置精确数值门迁到生成快照对照，显式旧覆盖仍保留 |
| test_mobility.cpp | 迁移率模型性质；原断言通过 |
| test_mos_interface_charge.cpp | MOS 电荷、界面性质；原断言通过 |
| test_newton_solver.cpp | Newton 方程/Jacobian/接触性质；原断言通过 |
| test_permittivity_assembly.cpp | 介电装配性质；原断言通过 |
| test_phumob.cpp | PhuMob 模型性质；原断言通过 |
| test_poisson.cpp | Poisson/部分材料覆盖性质；保留未覆盖的原迁移率，改为捕获真实原值 |
| test_poisson_charge_volume.cpp | 电荷体积性质；原断言通过 |
| test_poisson_neumann.cpp | Neumann 边界性质；原断言通过 |
| test_recombination.cpp | 复合/统计/温度性质；原断言通过 |
| test_scaling.cpp | SI/TCAD 单位性质；显式给定的旧比例常数保留 |
| test_schottky_contact.cpp | 势垒与接触解析性质；显式 Nc/Nv 常数保留，原断言通过 |
| test_sg_flux.cpp | SG 通量、KCL 与单位性质；原断言通过 |

对 Git 跟踪的 `configs/`、`examples/`、`tests/fixtures/`、`reference_tcad/` JSON 进行结构盘点：
识别 `mesh_file` 且含 `solver`/`sweep` 的 deck，包括模板的 `config` 容器，共 **77 个**：

| 目录/用途 | deck 数 | 依赖与处理 |
| --- | ---: | --- |
| configs/templates | 4 | PN2D IV/BV 与 bvmethods 两模板均显式声明材料文件；模板、原子 profile 与门限未变 |
| reference_tcad/genius_bjt_sentaurus2022/vela/configs | 5 | 完整显式材料；未改 |
| reference_tcad/pn2d_sentaurus2018/variants/inputs/{D1,D2,G1,G2,P0,P1,P2,P3} | 16 | 各自完整材料；未改 |
| reference_tcad/schottky_charon_sentaurus2018/vela | 8 | 两个无材料文件基准与四个部分覆盖变体改为完整旧快照；两个 full_materials 变体未改 |
| reference_tcad/simplemos_sentaurus2022/engineering/inputs/n17…n24 | 8 | 完整显式材料；未改 |
| reference_tcad/templates_ldmos_sentaurus2022/profiles 与 validation_d5 | 10 | 显式材料与既有绑定；未改 |
| reference_tcad/transportmodels_sentaurus2022/vela/configs | 12 | 完整显式材料；未改 |
| tests/fixtures/sdevice/**/expected_stage*.json | 14 | WP4 生成 deck golden 始终带材料文件；未改 |

其余 JSON 是 schema、合同、网格、报告或元数据，不按仿真 deck 计数；通用导出器生成的无材料 deck 按新全局默认运行。
`examples/sdevice` 是命令文件示例，经前端生成完整材料。全部 77 路径和原始/部分材料键盘点保留于
`build-release/wp13/dependency-inventory.{json,txt}`（变更前快照）。

## 历史校准路径与重新基线化

只有显式文件中已列出的字段会覆盖内置值，因此仅有 `materials_file` 还不足以保证旧校准路径不变。
[`generate_legacy_material_baselines.py`](../../scripts/generate_legacy_material_baselines.py) 生成以下完整旧快照，
保留原文件字段，只补省略的 pre-WP13 标量；**原始部分文件、封存 manifest、网格、状态和原 CSV 不变**。
冻结文件的 `ni/mun/mup/Nc_m3/Nv_m3` 使用历史 `unit_scaling` 原生单位，即 cm⁻³、cm²/(V·s)，
不能当成 SI JSON 直接加载。`_comment` 标明其旧基线用途。它们不是新默认库，也不要求旧 ni 与旧 Eg/DOS 自洽。

| 原输入 | 冻结输出/实际消费者 | 保持的有效参数 |
| --- | --- | --- |
| Schottky 两个无材料文件基准 | vela/materials_legacy_builtin_300k.json；simulation_iv 与 arclength | 完整旧 Si/SiO2 |
| bvmethods/vela/materials_sentaurus2018.json | 同目录 materials_sentaurus2018_legacy_300k.json；boundary_methods 与 iic_current_support 脚本 | 原 ni=1.4638914958767616e10，其余旧 Si；Nitride 原字段 |
| pn2d/source/pn2d_sentaurus2018_iv_materials.json | 同目录 *_legacy_300k.json；五个 reference 导入仿真及 sensitivity/localization 脚本 | 原 ni=1.4638914958767616e10，其余旧 Si |
| Schottky/vela/materials_mu_n_sentaurus2018.json | 同目录 *_legacy_300k.json；两个 mu_n deck | 原 mun=1417，其余旧 Si |
| Schottky/vela/materials_mu_n_nc_eg_sentaurus2018.json | 同目录 *_legacy_300k.json；两个 mu_n_nc_eg deck | 原 mun/Nc/Eg，其余旧 Si |
| singledevice/vela/materials_sentaurus2018.json | 生成 source 与 vela 两份 *_legacy_300k.json；两条导入分支 | 原量子/亲和势与其他材料字段，补旧 Si 基本标量 |
| tests/fixtures/pn2d_bv_rejected_transition/materials.json | 同目录 materials_legacy_300k.json；sealed transition Catch2 | 原 ni + 旧其余标量，原封存文件与哈希不动 |

新的 Catch2 将“完整旧内置材料 + 原部分文件”与“新的完整冻结文件”通过实际加载器逐字段比较，覆盖上述六个部分文件。
`--check` 验证八个冻结输出没有漂移。Node7 历史诊断生成器同样在复制后的输出里补旧标量；
SingleDevice 固定状态脚本改为复制 deck 实际声明的材料文件，避免继续复制旧固定文件名。
Schottky 两段工作流将声明的材料路径绝对化，避免隔离输出目录中的相对路径失效。

以下每个测试预期变化均有语义依据，容差未改；没有重新生成任何参考曲线来掩盖差异。

| 修改文件 | 旧预期 → 新预期及原因 |
| --- | --- |
| tests/test_mesh.cpp | 旧内置 ni/μ 数字 → 基本正性/载流子排序；新增独立快照门在 1e-12 下覆盖所有精确数值，显式材料旧值断言不变 |
| tests/test_poisson.cpp | 未覆盖 mun=0.135 → 与覆盖前读取的 mun=0.1417 相同；保持部分覆盖语义 |
| tests/test_impact_ionization.cpp | VTK 性质测试写死 ni=1e16 m⁻³ → 读取当前数据库 ni；构造状态与输出器使用同一材料，原源积分容差不变 |
| tests/test_dc_sweep.cpp | sealed parent 继续用原材料标量，由完整旧快照显式表达；40 次迭代及 max_iterations 失败证据断言不变 |
| tests/regression/test_reference_tcad_tools.py | PN2D 路径断言使用完整旧快照名称，并新增旧 ni/μ/Eg/Nc/Nv 断言；物理和门限不变 |
| scripts/check_sdevice_gate_contact_mapping.py 与对应测试 | 默认材料从新快照读取；Φ_S,i=4.6228018935… → 4.611003610676079 eV；n-type 1e17 cm⁻³ 平带预期下降约1.869972 mV；公式、1 mV 门限及 refinement 门不变，显式旧参数公式仍测试 |
| tests/regression/test_sdevice_materials.py | “与 C++ 遗留值不同” → 生成快照与 emitter 完全相同；增加 JSON 漂移拒绝和旧快照生成门 |
| tests/regression/test_schottky_reference_workflow.py | 增加冻结材料的存在/绝对路径检查；额外执行时发现既有重启后缀断言仍为 CSV，按当前工作流改为 HDF5 |

SiO2 χ 降低 0.05 eV 会影响启用绝缘体量子/带边作用的路径；本次没有更改相关方程。
`test_mos_interface_charge`、`test_newton_solver` 中 MOS/量子测试、已有 SingleDevice/TransportModels 材料合同与回归均在全量测试中通过。
它们分别检验性质或使用已显式给出的材料；SingleDevice 的显式 oxide χ=0.9 已在原文件中，不发生二次变化。
Poisson 介电装配依赖 ε，ε 未变；实际 runner 的栅平带三种 Barrier 复验通过，严格门限未放宽。

## 验证与证据

Windows 环境：MSYS2 UCRT64 GCC **16.2.0**，C++20，`windows-ucrt64-release`，Python API 关闭。
配置实际检测到 HDF5 静态库与独立 HDF5/HighFive 状态接口、SPQR、UMFPACK、METIS，以及
STRUMPACK/SuperLU_MT/MUMPS；默认 Eigen SparseLU 使用 COLAMD。可选后端不因此成为所有案例的默认。
初次配置后的 `-isystem D:/msys64/ucrt64/bin/../include` 触发 GCC `include_next <math.h>` 错误；
本地 cache 使用 `-DCMAKE_NO_SYSTEM_FROM_IMPORTED=ON` 后完整构建成功，没有修改全局工具链或升级依赖。

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --preset windows-ucrt64-release -DCMAKE_NO_SYSTEM_FROM_IMPORTED=ON
cmake --build --preset windows-ucrt64-release --parallel 4
ctest --preset windows-ucrt64-release -j 4 -E "^sdevice_e2e_"
ctest --preset windows-ucrt64-release -R "pn2d_config_templates|pn2d_node_volume_policy_default_acceptance|pn2d_bv_m2_mixed_voronoi_self_consistent_control"
python scripts/generate_default_models_par.py --check
python scripts/generate_legacy_material_baselines.py --check
python -m unittest tests.regression.test_schottky_reference_workflow
```

临时目录在工作树的 `build-release/wp13/test-tmp`，运行时设置 `TEMP/TMP` 指向该目录。
首轮 1026 项中三项失败：旧 PN2D 路径、VTK 隐含 ni、栅检查器旧参数；修复后相关六项 **6/6**，
最终全量 **1026/1026**，459.21 s。三个 PN2D 门包含在两轮全量及具名复验中。
全量包括新增三项 Catch2：双单位快照、温度质量作用关系、历史部分覆盖保持。WP3 材料 Python 回归为 **32 项**；
额外 Schottky 工作流 Python 测试 **3/3**。源码 ASCII 门通过，源默认库内容不变，物理模板/验收合同/容差 diff 已核对。
四个 `sdevice_e2e_*` 外部参考项未在此次本地全量中重跑；此前 `8dc5498c` 两平台原生包验收保留为历史包证据，
不标成 WP13 新二进制验收。本次没有重新发布部署包，也没有扩大 PN2D 场量资格。

额外 Ubuntu 26.04.1 验证尝试：配置成功（GCC15.2，HDF5::HDF5、SPQR、UMFPACK、METIS，SparseLU/COLAMD），
编译了新材料实现与对照测试对象。首次全目标构建包含许多此前未构建的测试目标，改为受影响目标时停止了本任务的 ninja；
随后 WSL 启动返回 `HCS_E_CONNECTION_TIMEOUT`，未完成链接与测试。**不声称 Linux WP13 已通过**；后续可在恢复的 Ubuntu
环境复验受影响目标。环境失败日志和被停止的构建均保留，没有将其当作产品测试失败或通过。
WSL 失败的工具调用结果摘要保存在 `build-linux26/wp13/environment-failures.txt`，不冒充原始标准输出。

生成证据位于忽略目录 `build-release/wp13/` 与 `build-linux26/wp13/`：

- configure/build*.log：实际功能检测、初次环境失败和成功构建。
- full-ctest.log：首轮失败；final-ctest.log 与 final-ctest.xml：最终 1026 项通过。
- fixes-ctest.log：三个修复项与三个 PN2D 门的 6/6。
- materials-python.log、materials-workflow-python.log、schottky-workflow-test.log：生成器与额外工作流测试；前者保留修改过程中的历史结果。
- dependency-inventory.json/txt：77 个 deck、25 个测试与部分材料输入清单。
- docs-and-probe.log、install-smoke.log：文档链接、JUnit 汇总、历史 Node7 材料完整化与安装资源检查。
- schottky-rerun/、schottky-control/、schottky-control-comparison.json：以下独立问题的完整对照。

最终全量日志 SHA-256：`97e997e28c265a41ed546fc5a1a696cffce8417661640f51ebf477eb11639e98`。

## 独立待验证项：Schottky 0.82 V 重启

额外重放历史 Schottky 完整 1 V 曲线，阶段 A 的 **42 点、终点0.8200000000000004 V**成功，
阶段 B 在0.82 V 初始点失败：`line_search_non_decrease`，51 次迭代，电子连续性残差
`2.3756890637732237e-9`，高于固定 `abstol=2e-9`。没有调大阈值或更改物理。

为区分旧问题与 WP13 回归，使用之前通过包完整性门的 Windows 发行包作为对照：
`PACKAGE_MANIFEST.commit=8dc5498cc4562fd7a5327225cd0bbce94ea8f796`，发行 ZIP SHA-256
`ce4d18b7b87d2418937126a91dc32f48b6dca51f84652a429b473bb09fbd1540`。
同一工作流和冻结材料分别执行新旧 runner：

| 对照 | 结果 |
| --- | --- |
| 阶段 A 两份 CSV | 字节完全一致；SHA-256 `270c31891437e86a712484e340e90af26900b07f02fea49a0f756ac9eb7e88c6` |
| 阶段 B | 同一初始点、失败类别、51 次迭代、全部数值列和残差；只有失败证据文件路径不同 |

这证明该重启失败在对齐前也存在，冻结材料保持有效参数和可重放的当前行为。
历史 2026-08-17 的完整曲线成功记录不能据此标为当前已复验；需要独立修复/重新验证当前 HDF5 重启路径。
本次材料对齐未顺带放宽其2e-9门限。WP14 剩余参考案例与 PN2D −20 V 空穴密度差异也继续保留待办。
