# `vela device`：SDevice 命令文件直跑工作流

状态：**WP5 编排器与 WP6.1 FieldData 已交付**（2026-10-02）。本文是
[开发计划](superpowers/plans/2026-09-26-sdevice-cmd-direct-run-development-plan.md)
第 4–6 节契约与第 10.0 节接口契约的可执行版本：schema 文件、夹具、消息目录以本文为准；
计划文档保留决策依据。WP10 已按原契约收口，资格范围与已知限制见第 3.3 节；
WP11 分发与并发验收已完成（2026-10-04），见[打包指南](sdevice_packaging.md)与
[验收记录](validation/sdevice_wp11_packaging_concurrency_2026-10-04.md)。

## 1. 契约文件索引

| 契约 | 文件 | 生产者 → 消费者 | 状态 |
| --- | --- | --- | --- |
| C1 执行 IR v2 | `schemas/vela.sentaurus_execution_ir.v2.schema.json` | WP1 → WP3/WP4/WP5 | 冻结 |
| C2 网格包 | `schemas/vela.mesh_qualification.v1.schema.json` + 现有 `mesh.json`/`doping.csv` 契约 | WP2 → WP4/WP5 | 冻结 |
| C3 材料文件 | `schemas/vela.templates_ldmos.materials.v1.schema.json`（WP3 不另设别名，见第 2.4 节） | WP3 → WP4 | 沿用现有 |
| C4 stage 计划 | `schemas/vela.stage_plan.v1.schema.json` | WP4 → WP5/WP8 | 发射器已交付 |
| C5 求解器协议 | `docs/config_schema.md` + `configs/schema/vela-simulation.schema.json` | WP6 → WP4/WP5/WP8 | WP6 已交付 |
| C6 输出文件契约 | 本文第 3 节 | WP8/WP6 → WP10 | 冻结 |
| C7 产物清单 | `schemas/vela.reference_artifact_manifest.v1.schema.json` | WP9 → WP10 | 冻结 |
| C8 CLI/退出码/消息目录 | 本文第 2、4 节 | WP0 → WP5/WP7/WP12 | 冻结 |
| C9 启动器协议 | 本文第 2.3 节 | WP7 → WP5 | 启动器与 zipapp 已交付 |
| C10 sidecar | `schemas/vela.sdevice_sidecar.v1.schema.json` | WP0 → WP4 | 校验与合并已交付 |

校验工具：`python scripts/vela_sdevice/schema_validator.py <schema> <json>...`
（纯标准库的 draft-07 子集校验器；安装了 `jsonschema` 时单测会交叉比对两者结论）。
解析器（WP1）：`python -m vela_sdevice.cmd_parser <deck.cmd>`（`PYTHONPATH=scripts`）打印 IR v2 并以
`verdict.exit_code` 退出；库接口 `vela_sdevice.cmd_parser.parse_cmd_file/parse_cmd_text`。
夹具：`tests/fixtures/sdevice/`（索引 `fixtures.json`，比较规则见第 5 节）。
契约单测：`tests/regression/test_sdevice_contracts.py`（CTest 名 `sdevice_contracts`）；解析器单测
`tests/regression/test_sdevice_cmd_parser.py`（CTest 名 `sdevice_cmd_parser`）。

## 2. CLI 与退出码（C8）

### 2.1 调用形式

```text
vela device [选项] <name_des.cmd | name_des | name>
```

命令文件解析顺序：字面路径存在 → `<name>.cmd` → `<name>_des.cmd`。未给命令文件或要求从标准输入
读取 → 用法错误（退出码 2）。工作目录 = 命令文件所在目录；`File{}` 内相对路径相对该目录。

镜像的 `sdevice` 选项：`-h`、`-v`、`-q`、`-n`、`-i`、`-P`、`--exit-on-failure`、
`--threads <int>`、`--max_threads <int>`、`--verbose`。拒绝（退出码 2，E202）：`-d -M -S -r -L*
-ver -rel -versions --tcl --xml --field-names --parameter-names --compiler-version`。
Vela 自有选项：`--check`、`--allow-lossy`、`--keep-intermediate/--no-keep-intermediate`、
`--vtu-binary`、`--plot-final-only`。

### 2.2 退出码

| 码 | 含义 | 消息编号段 |
| --- | --- | --- |
| 0 | 流程完成（含无 `ExitOnFailure` 时某 stage 未达 Goal） | — |
| 1 | 输入致命错误（解析、不支持的命令/模型、预处理指令） | E1xx |
| 2 | 用法错误 | E2xx |
| 3 | 输入文件缺失或未声明 | E3xx |
| 4 | 校验失败 | E4xx |
| 5 | 求解失败 | E5xx |
| 6 | 输出写入失败 | E6xx |
| 7 | 并发冲突 | E7xx |
| 8 | 被外部中止 | E8xx |

取得工作目录锁后的每个非 0 退出同时写 log 尾部 `[FATAL]` 摘要与
`<stem>_vela/fail_report.json{stage,message_id,message,details,exit_code}`。
未取得锁时保留持有者文件；命令文件定位和 argparse 用法错误发生在创建本次工作目录之前，只写 stderr。
CLI 输出和用户日志均为 UTF-8，包括 Windows 本地化系统错误详情。

### 2.3 启动器协议（C9，WP7）

`vela device <args>` → 子进程 `<python> bin/vela_sdevice.pyz device <args>`，cwd 不变，退出码透传。
`<python>` 查找顺序：`VELA_PYTHON` → `bin/python/` → PATH 中 `python3`/`python`（≥ 3.10，log 警告）。
Windows 用宽字符 API 透传参数；`bin/` 前置到 `PATH`。

CMake 默认构建生成 `<build>/bin/vela_sdevice.pyz`；`vela` 目标依赖该打包目标。每次构建重新记录当前 Git 提交。
zipapp 包含标准库之外的四个仓库模块（`sentaurus_parameter_ir`、`sentaurus_parameter_map`、
`compare_sentaurus_tdr_tdx`、`convert_tcad_export`）、默认材料库、schema 和版本记录；运行不依赖源码目录或第三方 Python 包。
`VELA_BIN` 可覆盖 runner/importer 目录；默认查找 zipapp 目录及其父目录，再查 PATH。
`VELA_EXAMPLE_RUNNER` / `VELA_SENTAURUS_IMPORT` 可指定对应工具。

在 UCRT64 环境的仓库根运行：

```powershell
cmake --build --preset windows-ucrt64-release
.\build-release\vela.exe device -h
.\build-release\vela.exe device --check .\case\device_des.cmd
.\build-release\vela.exe device --threads 4 .\case\device_des.cmd
```

### 2.4 材料发射器（C3，WP3）

库接口 `vela_sdevice.materials.emit_materials(par, active_models, temperature_K, mesh_materials,
electrode_materials=(), allow_lossy=False, band_edge_sensitive=False)`；命令行
`python -m vela_sdevice.materials --par <file> --ir execution_ir.json --materials Si SiO2 --out <dir>`。
输出 `materials.json`（C3；`schema` 恒为 `vela.templates_ldmos.materials.v1`，C++ 装载器只接受该 id，
因此不发布通用别名）与 `materials_report.json`（`vela.sdevice_materials_report.v1`：每个参数的
矩阵状态、来源文件/行、`origin ∈ {user_material, user_global, vela_default}`；`solver_parameters[]`
把 Masetti/Scharfetter/OldSlotboom/van Overstraeten/HighField 参数按矩阵 `target` 换成 SI 值交给 WP4；
`electrodes[]` 给出 `Electrode{Material=}` 的 `Bandgap.WorkFunction`；`inert[]`/`lossy[]` 列出被判定
无效果或有损的参数及理由）。

C3 的单位由自身 schema 与 `unit_system` 定义；装载到 `format_version: 2` deck 时不做 deck 键名转换。
它仍接受原生 `_cm3`、`_cm2_V_s` 键并拒绝混用 `_m3` 等未声明的键。

规则：

- 激活模型来自 C1 `physics.effective[].key`（`vela_sdevice.materials.active_models_from_ir`），映射到
  `scripts/sentaurus_parameter_map.py` 的模型名；参数状态由该矩阵给出，发射器只评估冻结定律。
- `share/vela/materials/models.par`（`scripts/generate_default_models_par.py` 生成，即 `vela device -P` 的输出）是默认库：
  Silicon 采用与 Sentaurus 默认参数集相同的文献值（Green 1990 带隙/态密度质量、Masetti 1983、Canali 1975、
  Scharfetter、Dziewior–Schmid Auger、Slotboom–de Graaff 1976、van Overstraeten–de Man 1970、Wettstein 2001），
  Oxide/Nitride 采用仓库已对齐的 Sentaurus 介电参数（ε 3.9/7.5，E_g 9/5 eV，χ 0.9/1.9 eV）；发射器支持的每个模型在库中都有节。
  Vela 不能表示的 Sentaurus 默认值改为 Vela 自身定律（`ConstantMobility.Exponent=2.2,2.2`，无 `mutunnel`，`Auger.H=0`，
  `HighFieldDependence.betaexp=vsatexp=0`）。用户 `.par` 未给出的节/参数逐参数回退到该库并记 `origin=vela_default`；
  网格材料在两者中都没有 → E403；激活模型所需的节在两者中都缺失 → E404 `(missing_section)`。
  C++ `MaterialDatabase` 内置 Si（nᵢ 1e10、N_c 2.8e19、N_v 1.04e19 cm⁻³、E_g 1.12 eV）仍是无材料文件 JSON deck 的
  遗留默认，与该库不同；`vela device` 始终经由材料文件，不使用内置值。
- 半导体标量：`Eg(T)` 为 Varshni；`Nc/Nv` 用 Formula 1 质量律与 `2(2π m kB T/h²)^{3/2}`（CODATA 2018，与
  `include/vela/core/PhysicalConstants.h` 一致）；`ni = sqrt(Nc Nv) exp(−(Eg(T)+dEg0)/2kT)`，`dEg0(OldSlotboom)`
  只在 OldSlotboom 激活时折入 `ni`，不进入 `bandgap_eV`；`mun/mup` 取 `ConstantMobility.mumax`。
- 绝缘体只需 `Epsilon`（可带 `Bandgap`）；`Region=`/`*Interface=` 作用域含激活参数 → E404 `(unsupported_scope)`。
- `approximated` 需 `--allow-lossy`；`unsupported_*`/`unmapped` 致命。发射器对以下矩阵行做有据可查的“无效果”
  降级并写入 `inert[]`：`Bandgap.alpha2=0`、`beta2`（alpha2=0）、`EgMin=0`、`dEgMin`（EgMin=0）、
  `Bgn2Chi`（无带边敏感接触/异质界面时；`--band-edge-sensitive` 保持有损）、`ConstantMobility.mutunnel`
  （无隧穿模型）、`ConstantMobility.Exponent`（两载流子均为 2.2，或 T=300 K）、
  `HighFieldDependence.betaexp/vsatexp`（T=300 K）；`Auger.A/B/C` 在冻结温度上精确合成为常数。
- 已知偏差：历史参考 `pn2d_sentaurus2018_iv_materials.json` 的 `ni=1.4638914958767616e10` 无法由任何记录在案的
  公式复现（最接近 −2.7e-4），发射器以闭式公式为准并在测试中以 5e-4 核对；C++ 内置 Si 的 `ni=1e10` 与其
  `Nc/Nv/Eg` 不自洽（推导值 6.68e9），默认库往返因此只对 `ni` 以外的字段闭合。

单测 `tests/regression/test_sdevice_materials.py`（CTest 名 `sdevice_materials`）。

### 2.5 栅接触参考能级（C4 发射规则，WP4a）

手册 Eq. 104 `ψ = φ_F − Φ_MS` 与 Vela `metal_gate`（`ψ = bias − work_function_eV`）共用本征半导体参考能级，
映射精确（平带电压夏具偏差 2e-16 V，见
[栅接触核对](validation/sdevice_gate_contact_reference_2026-09-29.md)）：

- `Electrode{ Barrier=Φ_MS }`（缺省 0）→ `contacts[]: {type:"metal_gate", work_function_eV: Φ_MS}`；
- `Electrode{ WorkFunction=Φ_M }` → `work_function_eV = Φ_M − (χ + E_g/2 + (kT/2)·ln(N_c/N_v))`，参数取 C3 中参考半导体（Si）的值；
- `Electrode{ Material="X" }` → `Φ_M` 取 WP3 报告 `electrodes[].work_function_eV`（`.par` `Bandgap{WorkFunction}`），再按上条换算；
- `Material="Silicon"(N=...)` 多晶硅栅与氧化层固定电荷不在本规则内（WP4 处理）。

当前 WP3 默认库在 300 K 给出 `Φ_S,i=4.611003610676079 eV`；2026-10-02 用户同意将 WP4
零偏移验收输入改为该值。遗留内置 Si 的 `4.6228018936 eV` 用作输入时应得到
`work_function_eV=0.011798282923921377 eV`，测试同时覆盖这两个结果；公式和材料参数未改。

### 2.6 求解器新增键（C5，WP6）

`format_version: 2` deck 的 `sweep` 块新增（定义与默认值见 `docs/config_schema.md`）：`vtk_format` (`legacy`|`vtu`)、
`vtk_fields`（TDR 数据名，`Name` / `Name/Vector`）、`vtk_binary`、`vtk_bias_points`（必经点，1e-12 V 落点，不改变
其它自适应点）、`csv_all_contacts`（每接触 `contact_<name>_{electron,hole,total}_A_per_um` 列），以及
`breakdown.current_{abs,min,max}_A_per_um` + `breakdown.current_contact`（电流终止判据，所有 sweep 模式；触发点的
`criterion` 为 `current_abs|current_min|current_max`，CSV 增加 `stop_criterion`/`stop_current_A_per_um` 列）。
`.vtu` 字段目录与单位由 `tdrFieldCatalogue()`（`include/vela/io/VtuWriter.h`）定义；不可提供的请求字段记入
`DCSweepPoint::skippedVtkFields` 并写 runtime log `vtk_fields_skipped`，由 WP5 转为 W202。

WP6.1 增加可选 `sweep.vtk_field_data` 对象：ASCII 标识符键、字符串或有限数值；拒绝布尔、嵌套值和内置键重名。
例如 `{"CmdSHA256":"abc","Stage":1,"t":0.25}` 原样写入 FieldData。提供该对象时 runner 另写
`GoalReached`/`goal_reached`（本段目标到达为 1，否则 0）以及被扫电极的内外电压；未给 `t` 时按本段电压区间
计算，Current 模式按实测初始电流到最终电流目标计算。二进制 VTU 的 FieldData 同样可读；不提供该对象保持原有输出字节。

### 2.7 stage/deck 发射器（C4/C10，WP4）

库模块 `vela_sdevice.deck_emit` 只消费 accepted C1、C2 网格资格报告与 WP3 材料报告；不读 cmd 文本、不运行求解器：

```python
emit_stage_plan(ir, *, mesh_qualification, materials_report, work_dir,
                cmd_basename, cmd_sha256, options=None, sidecar=None)
# -> StagePlanResult(stage_plan, decks, decisions)
rebase_stage(result, stage_index, reached_voltages_V)
write_stage_plan(out_dir, result)
```

`EmitOptions` 包含 `exit_on_failure`、`threads`、`plot_final_only`、`vtu_binary`、
`keep_intermediate`、`allow_lossy`。`write_stage_plan` 先校验再写 `stage_plan.json`、
`simulation_stageNN.json` 与 `decisions.json`。路径为交给 WP5 暂存的 ASCII 基名。
独立 CLI（从仓库根设置 `PYTHONPATH=scripts`）为：

```text
python -m vela_sdevice.deck_emit --ir execution_ir.json --mesh-qualification mesh_qualification.json --materials-report materials_report.json --out device_vela
```

- 连续初始 Coupled 合并为平衡 stage；后段通过 `sweep.initial_state_file`/`write_state_file` 串接 HDF5。
  Save 改写前段状态文件名并保存电压表；Load 用 `frozen_state` 恢复快照。同一段多个 Save 别名暂不支持，报 E407。
- 电压步长乘以本段 `|V1−V0|`；反向扫描保留 `step` 的负号。缺省采用手册 T-2022.03 Table 348：
  InitialStep=0.1、MinStep=0.001、MaxStep=1、Increment=2、Decrement=2。
  C4 `shrink_factor` 保留 Decrement 除数；runner deck `shrink_factor=1/Decrement`。
  零电压跨度变为单点平衡并记录 W211。
- Plot/CurrentPlot 的 Time 范围展开后取排序去重并集；`vtk_bias_points` 保证电压必经点，C4 记录写出选择和全局 plot 序号。
  Current 目标使用 `voltage_to_current`；达到后，WP5 必须将实测接触电压传给 `rebase_stage`。
  纯 WP4 发射器在没有实测初始电流时仍对显式 Time 中间点报 E407；WP5 用 `defer_current_times=True`
  先生成可校验计划，在运行 Current 段前用 `bind_current_start` 绑定前段 CSV 的实测电流与电压。
  `--check` 不求解，Current Time 映射会留待实际执行时绑定。
  未达 Goal 或 Current 后续接时，WP5 调用 `rebase_stage(..., update_saved_state=True)`，
  使前段 Save 对应的电压表也记录实测值；默认库调用继续保留原快照语义。
- 数值参数取 WP3 报告，按公开 v2 单位键写出；`band_to_band` 保留其 SI 约定，Auger 合成值由 SI 转 cm⁶/s。
  Scharfetter `taumax` 同时写入常数 SRH 的 `tau_n_s/tau_p_s`，避免回退到 C++ 内置寿命。
  材料参数声明不激活模型：普通 `SRH` 显式写 `srh_doping_dependence.enabled=false`，
  `SRH(DopingDependence)` 才启用掺杂寿命。普通欧姆接触固定两种载流子的准费米电势，
  前端显式关闭 Newton 的 minority-electron relaxation。
  缺省高场迁移率和雪崩驱动力按 SDevice DD 语义取准费米梯度；显式驱动力优先。
  雪崩采用 `generation=current_density` 与通用 SG edge `current_approximation=density_gradient`。
  此映射不启用 PN2D BV 模板专用的 SG/GSS-Laux 与 mixed-Voronoi 原子配置。
  不同材料要求互异的全局求解参数时报 E407。材料界面固定电荷展开为网格 region 对。
- 初始平衡保留 Gummel/Newton 初始化；已有状态的 Coupled 段用 warm-start Newton，
  `quasi_fermi_reference=contact_majority`，避免逐点 Gummel 破坏延续初值。
  同时启用现有局部连续性行检查（0.001）、行缩放和全局源/接触守恒检查（1%），
  单次准费米更新上限 0.1 V；不放宽已有验收门限，`min_newton_max_iter` 保持源 Iterations。
- metal_gate 使用第 2.5 节公式；Schottky 使用绝对功函数或 `barrier_eV`，采用 `thermionic_robin` 与 Newton。
  串联电阻仅允许配置在本段扫描接触；与电流目标互斥。
- sidecar 先按 C10 白名单校验，再递归合并对象、替换数组与标量；每个覆盖键记录 W207。
  不允许覆盖材料物理参数。`mixed_voronoi` 要求 C2 的非钝角资格、无钝角，
  最大角不超过 90°（仅容许 1e−8° 浮点误差），与当前运行时 non-obtuse 检查一致。
- sidecar 可显式选 `avalanche_current_support_profile=element_edge_sg_gss_laux` 或
  `legacy_cell_reconstructed`，应用完整的电流支撑、源映射、中点密度、节点体积和网格资格配置；
  拒绝拆分几何覆盖，必须已有 GradQF avalanche。它不激活模型或覆盖系数，记录 W207。
  前者依赖非钝角资格；后者完整恢复 barycentric/no-guard 配置。源导数按 Math 使用 local AD 或 frozen，
  接触元素保留原生 electric-field fallback。
  PN2D BV 参考案例随附显式 sidecar；验收报告记录其哈希，复用旧输出也核验它。
- C4 和公开 deck schema 均校验；额外遍历拒绝非有限值与 v1 键。错误使用 `DeckEmitError(message_id, **fields)`。

CTest `sdevice_deck_emit` 对七个 accepted 夹具比较完整 golden，并验证真实 runner 的平衡、电压扫描、电流目标与 HDF5 串接。
夹具测试补 `Parameter="models.par"` 后重新解析，不改变原始 IR golden。

### 2.8 编排器执行与输出（WP5）

`python -m vela_sdevice device ...`（源码运行时 `PYTHONPATH=scripts`）与原生入口使用同一 CLI。
`--check/-n` 完成 C1/C2/C3/C4/C10 gates 后写 manifest，不调用 runner；`-i` 仅执行首段平衡。
每个接受点回显一行，`-q` 关闭进度；线程数同时传入计划和四个 OMP/BLAS/Eigen 环境变量。
HDF5 子进程固定 `HDF5_USE_FILE_LOCKING=FALSE`。第一段失败、异常进程结束或 ExitOnFailure 返回 5；
runner 退出 2 映射为生成 deck 被拒的退出 4，并保留其输出详情。

工作目录通常为 `<stem>_vela`；非 ASCII 命令文件基名使用 `device_<sha256前8位>_vela`，保持 runner 暂存基名为 ASCII。

锁采用排他创建，记录 PID、主机、时间及所有权 token。存活进程、异地主机或无法读出所有权的锁返回 7；
只有同主机已死进程的锁自动回收，记录 W701。SIGINT/SIGTERM/Windows CTRL_BREAK 中止子进程、
保留最后接受点的 PLT/VTU/PVD，释放本进程锁，返回 8。

PLT 只发布完整且收敛的 CSV 行，通过临时文件替换增量更新；读取时使用稳定行快照。
VTU 在文件闭合后按 C4 全局编号及 FilePrefix 复制，PVD 原子更新，最终文件始终来自最后接受点。
`--plot-final-only` 只发布最终 VTU；成功后的 `--no-keep-intermediate` 移除工作目录阶段 CSV/VTU，
保留状态文件、输入、计划和诊断记录。失败与中止保留中间产物。FieldData 的 goal 状态针对该段；
整个流程的完成状态见 `run_manifest.json` 的退出码和各段记录。

runner 以 ASCII 暂存名读写输入与运行日志，编排器在段结束后将 runtime 日志合并至唯一用户日志，
以支持 Windows 中文和空格路径。用户日志含版本/提交/哈希、cmd 回显、映射决策、每段结果与结束摘要。
`.plt` 电流符号因子全局固定为 +1，多段 time 每段重置；WP10 的原生参考证据见第 3.3 节。

## 3. 输出契约（C6）

### 3.1 文件命名

手册附录 C（Table 194）规则：基名含 `.` 的字符串按字面使用，否则追加 `_des.<ext>`
（命令 `_des.cmd`/`.cmd`，日志 `_des.log`，参数 `.par`，网格 `.tdr`，Save/Load `_des.sav`，
器件 Plot `_des.tdr`，电流 `_des.plt`）；Quasistationary 内的 plot/save 按全局序号编号
`<prefix>_NNNN_des.<ext>`；`Output` 缺省为 `output`；`NewCurrentPrefix` 前置到电流文件名。

**观测到的偏差（T-2022.03-SP2，按实测行为实现）**：`File{Output}` 不遵守“含 `.` 即字面”规则——
`Output="pn2d_0v.log"` 生成 `pn2d_0v.log_des.log`，`Output="singledevice.log"` 生成
`singledevice.log_des.log`，而 `Output="schottky_n_des.log"` 生成 `schottky_n_des.log`。
因此 Vela 的规则为：`Output` 值已以 `_des.log` 结尾则字面使用，否则追加 `_des.log`。
`Current`/`Plot`/`Grid`/`Parameter` 遵守手册规则（证据：`Current="pn2d_0v.plt"` → `pn2d_0v.plt`；
`Current="bjt_m0"` + `NewCurrentPrefix="bjt_m0_base_"` → `bjt_m0_base_bjt_m0_des.plt`；
`Current="schottky_n_iv"` + `NewCurrentPrefix="forward_"` → `forward_schottky_n_iv_des.plt`）。
证据目录见 [产物盘点](validation/sdevice_reference_artifact_inventory_2026-09-26.md)。

Vela 替换：器件 plot `.tdr` → `.vtu`；`.sav` → `.h5`；另生成 `<plot stem>.pvd` 序列索引
（timestep = 全局偏置点序号）。

### 3.2 目录布局

```text
<cmd 目录>/
  name_des.cmd
  <Output 解析名>                      # 单一 log，顺序追加
  <NewCurrentPrefix><Current 解析名>   # DF-ISE xyplot，每接受步追加（写 .tmp 再 rename）
  <Plot 解析名 .tdr→.vtu>             # 末态；失败/中止时为最后接受的状态，FieldData goal_reached=false
  <FilePrefix>_NNNN_des.vtu ...        # Solve 内 Plot(Time=...)；全局序号
  <Plot stem>_NNNN_des.vtu ...         # 默认每接受偏置点一个（--plot-final-only 关闭）
  <Plot stem>.pvd
  name_des_vela/                       # 中间产物：execution_ir.json mesh.json doping.csv materials.json
                                       # mesh_qualification.json stage_plan.json simulation_stageNN.json
                                       # stageNN_state.h5 stageNN.csv run_manifest.json fail_report.json .lock
```

### 3.3 `.plt` 与 `.vtu`

- `.plt` 数据组：`time` + 每个电极的 `OuterVoltage`、`InnerVoltage`、`eCurrent`、`hCurrent`、
  `TotalCurrent`（Table 179）；2-D 电流单位 A/µm；`time` = 当前 Quasistationary 的 `t`；
  WP10 已用通过六项参考核对的 T-2022.03-SP2 v2 BJT m0 base/collector 文件确认：
  **每段重置为 0～1，不累加**（36/31 行；SHA-256 分别为
  `dd05b1888149055604614241312f51e8ba3da32028f5cf7b40d131d5aa391c26`、
  `0aff50fa989d04249b606e9b811a1fd655a7ea545ab4dcefe2e8293bfae0ed9b`）。
  电流符号全局固化为 `CURRENT_SIGN_FACTOR=+1.0`；依据同样通过核对的 v2 `pn2d_iv.plt`
  首三个正向 Anode 点及 Vela 原始 Anode 电流方向；PLT SHA-256
  `a5f712e3ba0bfeea221e0de681ce380b3ad4ab14eb42c63b1324eac431709b25`。
  符号核对采用已接受原始电流的 log-current 插值，只确定方向；数值幅值另行验收。
  当前结果见 [WP10 收口报告](validation/sdevice_wp10_closeout_2026-10-03.md)：原契约验收完成。
  IV 到达 10 V 并通过既定窗口，BV 到达 −20 V 且全窗口曲线/KCL 通过；BV 独立场量资格
  仅限绑定输入/配置的 −13.2 V。该点不是击穿电压；现有分段 log(I) 膝点拟合给出原生与
  Vela 均 −19.405 V，−19.5/−19.6/−19.7 V 场量诊断通过，−20 V 空穴密度诊断失败。
  附加诊断失败随 JSON 保留，不推导全反偏范围的空间资格。后续核对可复用已完成候选，
  cmd/par/grid/sidecar 必须通过指纹检查；配置漂移不能自动沿用 BV 场量资格。
- `.vtu`：XML UnstructuredGrid，Float64，ASCII 默认；字段名/单位用手册附录 F 的 TDR 数据名
  （计划附录 B），向量为同名 3 分量数组；`region_id` CellData；`FieldData` 含 `VelaVersion`、
  `CmdSHA256`、`Stage`、`t`、各电极 `OuterVoltage/InnerVoltage/TotalCurrent`、`goal_reached`。
  只写 `Plot{}` 请求的字段，隐式追加 `ElectrostaticPotential`、`eDensity`、`hDensity`、
  `DopingConcentration`（`PlotExplicit` 时不追加）；无法提供的字段 W202 跳过。

## 4. 消息目录（C8）

编号稳定、文案英文、每条带源文件行号（无行号时 `line=null`）。平台按编号匹配。
`{...}` 为替换位；`text` 字段保存原文 token。

| 编号 | 级别 | 退出码 | 模板 | 触发 |
| --- | --- | --- | --- | --- |
| W101 | warning | — | `Math keyword '{kw}' is unknown and ignored` | `Math{}` 中不在 Table 211 的关键字 |
| W102 | warning | — | `Math keyword '{kw}' is recorded and ignored (no Vela equivalent)` | `NotDamped`、`Transient`、`Method`、`SubMethod`、`CNormPrint`、`WallClock`、`CheckUndefinedModels`、`ParameterInheritance` 等 |
| W103 | warning | — | `Coupled option '{kw}' is recorded and ignored` | `Coupled(LineSearchDamping= Method= ...)` 中除 `Iterations` 外的选项 |
| W104 | warning | — | `Plot field '{name}' is not a known data name or alias; skipped` | `Plot{}` 未识别名字 |
| W105 | warning | — | `CurrentPlot coordinate sample '{text}' is ignored` | `CurrentPlot{ <Field>((x y)) }` |
| W106 | warning | — | `CurrentPlot integral or extremum '{text}' is ignored` | `CurrentPlot{ X(Integrate(...)) / Maximum(...) }` |
| W107 | warning | — | `File keyword '{kw}' is recorded and ignored (extra output only)` | `NewtonPlot=`、`ACExtract=`、`Extraction=`、`NonLocalPlot=` 等纯附加输出 |
| W108 | warning | — | `Plot option '{kw}' is ignored` | `-Loadable`、`DatasetsFromGrid`、`-PlotLoadable` 等 |
| W109 | warning | — | `Keyword '{kw}' appears more than once; the last occurrence is used` | 同段重复关键字 |
| W110 | warning | — | `Quasistationary step parameter '{kw}' is absent; sdevice default assumed` | 缺 `InitialStep/MinStep/MaxStep/Increment/Decrement` |
| W111 | warning | — | `Extrapolate option '{kw}' is ignored` | `Extrapolate(Order=...)` |
| W112 | warning | — | `Electrode attribute 'Resistor' accepted as an alias of 'Resist' (Table 206)` | `Resistor=` |
| W113 | warning | — | `Electrode attribute '{kw}={value}' equals the default and is dropped` | `Poisson=Dirichlet`、`AreaFactor=1` |
| W114 | warning | — | `Physics model '{key}' is recorded as metadata only` | 有效集中 `status=metadata_only` 的模型 |
| W201 | warning | — | `Math '{kw}={value}' mapped heuristically to '{vela}'` | `Digits/ErrRef/Iterations` → `solver.*` |
| W202 | warning | — | `Plot field '{data_name}' is not available in Vela; skipped` | 如 `eTemperature` |
| W203 | warning | — | `Implicit plot fields appended: {list}` | 无 `PlotExplicit` 时 |
| W204 | warning | — | `Device plot '{tdr}' is written as '{vtu}'` | `.tdr`→`.vtu` |
| W205 | warning | — | `State file '{sav}' is written as '{h5}'` | `.sav`→`.h5` |
| W206 | warning | — | `Stage {n} stopped at {contact}={V} V before reaching its goal {goal} V; continuing from the last accepted bias` | 未达 Goal 且无 `ExitOnFailure` |
| W207 | warning | — | `Sidecar override applied to stage {n}: {key} = {value}` | `<stem>.vela.json` |
| W208 | warning | — | `Load restores electrode voltages to {table}` | `Load(FilePrefix=)` |
| W209 | warning | — | `Schottky contact '{name}': default recombination velocities {e}/{h} cm/s injected` | 缺 `e/hRecVelocity` |
| W211 | warning | — | `Stage {n} has zero voltage span; emitted as a single equilibrium point` | WP4 退化扫描 |
| W210 | warning | — | `Thread count {n} taken from Math{NumberOfThreads}` | 无 `--threads` |
| W701 | warning | — | `Stale lock from PID {pid} on {host} replaced` | 锁文件进程不存在 |
| E101 | fatal | 1 | `Workbench preprocessor directive '{text}' is not supported; preprocess the deck first` | `#if #elif #else #endif #setdep #define #include #noprint #header #endheader #split #endsplit #exit` |
| E102 | fatal | 1 | `Workbench placeholder '{text}' is not expanded` | `@…@` |
| E103 | fatal | 1 | `Tcl block '!( ... )!' is not supported` | `!(` … `)!` |
| E104 | fatal | 1 | `Syntax error: {detail}` | 括号不匹配、未知段名、非法赋值 |
| E105 | fatal | 1 | `Section '{name}' is not supported` | `System{}`、`Thermode{}`、`Device{}`、`Optics{}`、`Transient`（段级） |
| E106 | fatal | 1 | `Physics model '{key}' is not supported` | `IALMob`、`Thermodynamic`、`Hydrodynamic`、`eTemperature`、`hQuantumPotential`、非均匀 `FixedCharge` 以外的 `Trap(s)`、`Okuto/Lackner/UniBo` |
| E107 | fatal | 1 | `Physics token '{text}' is not recognised` | 未知模型名/占位符 |
| E108 | fatal | 1 | `Electrode attribute '{kw}' is not supported` | Table 206 中未映射属性（`Charge= Current= DistResist EqOhmic Neumann AreaFactor≠1 Material=…(N=|P=)`） |
| E109 | fatal | 1 | `Solve statement '{text}' is not supported` | `System(...)`、`Transient`、`ACCoupled`、`Plugin`、`Continuation`、`Unset`、`Set` 非 `mode current` |
| E110 | fatal | 1 | `Goal form '{text}' is not supported` | 多 `Goal`、`Regexp`、`Model=/Parameter=`、`DopingWell` |
| E111 | fatal | 1 | `Electrode attribute '{kw}' is not defined in Table 206` | 未知属性 |
| E112 | fatal | 1 | `File keyword '{kw}' changes inputs or physics and is not supported` | `Load= Save= MobilityDoping= LifeTime= OpticalGenerationInput=` 等 |
| E113 | fatal | 1 | `Effective Physics model sets differ between semiconductor regions: {detail}` | 作用域展开后不一致 |
| E114 | fatal | 1 | `Grid file '{name}' has an unsupported extension (expected .tdr or .grd)` | |
| E115 | fatal | 1 | `Command file cannot be decoded as UTF-8/ASCII text or is empty` | |
| E116 | fatal | 1 | `Coupled equation '{name}' is not supported` | `Contact Circuit Temperature eTemperature hTemperature hQuantumPotential` |
| E117 | fatal | 1 | `Electrode Name=Regexp(...) is not supported` | |
| E201 | fatal | 2 | `No command file given (reading from stdin is not supported)` | |
| E202 | fatal | 2 | `Option '{opt}' is not supported by vela device` | |
| E301 | fatal | 3 | `File{Parameter} is not declared: vela device requires a parameter file (add Parameter="<name>.par" or export the defaults with 'vela device -P')` | 缺 `Parameter`/`Parameters` |
| E302 | fatal | 3 | `Grid file '{name}' not found` | |
| E303 | fatal | 3 | `Parameter file '{name}' not found` | |
| E304 | fatal | 3 | `DF-ISE companion file '{name}.dat' not found` | |
| E305 | fatal | 3 | `Command file '{name}' not found` | |
| E401 | fatal | 4 | `Electrode '{name}' has no matching contact in the mesh (contacts: {list})` | |
| E402 | fatal | 4 | `Mesh is not a 2-D triangle mesh ({detail})` | |
| E403 | fatal | 4 | `Material '{name}' is neither built in nor defined in the parameter file` | |
| E404 | fatal | 4 | `Parameter '{path}' of active model '{model}' cannot be imported ({status})` | `approximated` 无 `--allow-lossy`、`unsupported_*` |
| E405 | fatal | 4 | `Mesh coordinates are implausible (bounding box {box} um)` | |
| E406 | fatal | 4 | `Generated deck '{deck}' fails schema validation: {detail}` | |
| E407 | fatal | 4 | `Stage chain is inconsistent: {detail}` | 状态文件依赖、`t` 越界、起始电压不一致 |
| E408 | fatal | 4 | `Node-volume policy '{policy}' requires a non-obtuse mesh ({n} obtuse elements)` | |
| E409 | fatal | 4 | `Sidecar '{file}' is invalid: {detail}` | 白名单外键、schema 失败、合并后 deck 失败 |
| E410 | fatal | 4 | `Parameter file '{name}' cannot be parsed: {detail}` | `.par` 语法无法分类（`ParParseError`） |
| E501 | fatal | 5 | `Equilibrium stage failed: {detail}` | 首个 stage 不收敛 |
| E502 | fatal | 5 | `Stage {n} failed with ExitOnFailure: {detail}` | |
| E503 | fatal | 5 | `Solver terminated abnormally (exit code {code})` | runner 崩溃 |
| E601 | fatal | 6 | `Cannot write output file '{name}': {os_error}` | |
| E602 | fatal | 6 | `Output writer failed for '{name}': {detail}` | |
| E701 | fatal | 7 | `Another vela device run holds '{lock}' (PID {pid} on {host})` | |
| E801 | fatal | 8 | `Interrupted by signal {sig} at stage {n}, {contact}={V} V; partial outputs kept` | |

## 5. 执行 IR v2 约定（C1 补充）

schema 之外、golden 夹具所依赖的确定性规则：

1. **注释与预处理扫描**：先在原文上扫描 `#<directive>`（仅第 4 节 E101 列出的指令名，
   `#Thermal`/`#Plot(` 等视为注释）、`@…@`、`!( … )!`；再剥离 `*`/`#` 到行尾的注释（`*`/`#`
   必须是 token 起始；引号内文本不受影响）。
2. **哈希**：`source.sha256` = 去 UTF-8 BOM、CRLF→LF 后的文本 SHA-256；BOM/CRLF 变体产生相同 IR。
3. **文件名**：按第 3.1 节；`files.current.raw` 为 null 时 `resolved` = `<cmd stem>_des.plt`
   （stem 去掉 `_des.cmd`/`.cmd`）；`files.plot.vela` = `resolved` 的 `.tdr`→`.vtu`。
4. **电极属性**：只记录 deck 写出的属性；`Resistor=` 折叠为 `resist_ohm_um` 并发 W112；
   `Poisson=Dirichlet`/`AreaFactor=1` 丢弃并发 W113；其它 Table 206 属性 → E108，未知 → E111。
5. **Physics 规范键**（`effective[].key`，按 ASCII 升序排列）：
   `EffectiveIntrinsicDensity.OldSlotboom`、`Fermi`、`Mobility.DopingDependence`（别名 `DopingDep`）、
   `Mobility.HighFieldSaturation`（别名 `HighFieldsaturation`；`e/hHighFieldSaturation` 保留前缀为
   `Mobility.eHighFieldSaturation`/`Mobility.hHighFieldSaturation`）、`Mobility.Enormal`、
   `Recombination.SRH`（`parameters.doping_dependence`/`temperature_dependence`）、`Recombination.Auger`、
   `Recombination.Avalanche`（`parameters.model`、`parameters.driving_force`，均只在写出时存在）、
   `Recombination.Band2Band`（`parameters.model`）、`eQuantumPotential`。`Temperature=` 不进入
   `global[]`，写入 `temperature_K`（缺省 300）。`vela` 字段为说明性映射字符串，WP4 可修正其内容
   并同步更新 golden。
6. **Plot 别名 → 数据名**：`Potential→ElectrostaticPotential`、`eQuasiFermi→eQuasiFermiPotential`、
   `Doping→DopingConcentration`、`eCurrent→eCurrentDensity`（`TotalCurrent→TotalCurrentDensity`）、
   `SRH/SRHRecombination→srhRecombination`、`Auger→AugerRecombination`、`Band2Band→Band2BandGeneration`、
   `ConductionBand→ConductionBandEnergy`、`ValenceBand→ValenceBandEnergy`、`Affinity→ElectronAffinity`、
   `BandGapNarrowing→BandgapNarrowing`、`eENormal→eEnormal`；其余与附录 F 同名。`/Vector`、`/Element`
   后缀决定 `component`，`name` 不含后缀。
7. **Math 默认值**：`iterations=50`、`digits=5`、`rel_err_control=true`、`extrapolate=false`、
   `derivatives=true`、`aval_derivatives=true`、`compute_ionization_integrals=false`、
   `break_at_ion_integral=null`、`exit_on_failure=false`、`number_of_threads=1`；`ErrRef` 缺省 null。
   记录并忽略的关键字进入 `math.ignored`（关键字按原拼写）并各发一条 W102。
8. **Solve 语句与 stage**：`Poisson` 裸语句 → `statements[].kind="poisson"`，stage 视同
   `Coupled{Poisson}`；`Coupled` 无 `Iterations` → `iterations=null`；顶层 `Plot(FilePrefix=)`
   折叠到其前一个 stage 的 `outputs.plot_files`（`time=null`，文件名 `<prefix>_des.tdr`）；
   Quasistationary 内 `Plot(... Time=...)` 文件名 `<prefix>_NNNN_des.tdr`；`NewCurrentPrefix` 更新
   后续 stage 的 `current_prefix`；`Set("c" mode current)` 记入后续 stage 的 `current_mode_contacts`。
   可选键（`goal`、`t_control`、`options`、`current_mode_contacts`）只在有内容时出现。
9. **analysis.kind 规则**：R1 无 Quasistationary → `equilibrium`；R2 单个电压扫描且存在
   `Recombination.Avalanche`、`ComputeIonizationIntegrals`/`BreakAtIonIntegral` 或电流
   `BreakCriteria` → `bv`；R3 单个电压扫描（其它） → `iv`（单个电流目标扫描也记 `iv`，
   `reason` 为 `"R3: one Quasistationary with a current goal"`）；R4 多个 Quasistationary → `multi_stage`；
   R5 在门 1 被拒（预处理指令/占位符/Tcl 块）或文件不可解码（E115）→ `rejected`，各节不再求值，
   `electrodes`/`solve.statements` 为空；其它致命发现仍按 R1–R4 推导（`verdict` 单独表达拒绝）。
   `reason` 字符串固定为 `"R<n>: <本节描述>"`（见 golden 与 `scripts/vela_sdevice/cmd_parser/solve.py`）。
10. **verdict**：`accepted` = `errors` 与 `preprocessor_findings` 皆空；`exit_code` 取致命发现中
    最小的退出码段（E1xx→1，E301→3）；`fatal_ids` 按首次出现顺序去重。

## 6. 夹具与比较规则

`tests/fixtures/sdevice/fixtures.json` 列出每个 deck 的来源、LF 归一化哈希、golden 与比较模式：

- `compare="full"`：解析器输出与 `*.ir.json` 作为解析后的 JSON 相等（键序无关）。
- `compare="findings"`：只比较 `verdict`（完全相等）以及 golden 中出现的
  `preprocessor_findings`/`errors` 列表（逐项 `{id,line,text}` 相等，顺序一致）。

主验收 deck（pn2d、schottky、genius_bjt）原文都没有 `Parameter=`，golden 因而带 `E301`；
验收运行使用追加了 `Parameter="models.par"` 的副本（该文件即 `sdevice -P` 的默认参数导出），
见盘点文档待决项 A-2。

## 7. 待决与后续

实施状态、剩余工作包票据与交接说明见 [2026-09-30 交接文档](validation/sdevice_cmd_frontend_handoff_2026-09-30.md)。
该交接是当日快照；WP4 的后续实现、基准修正与完成验证见
[2026-10-02 WP4 验证记录](validation/sdevice_deck_emit_validation_2026-10-02.md)。

- A-1 VM 归档布局、A-2 `Parameter` 缺失处理、A-3 重跑时机、A-4 多段 `.plt` 的 `time` 语义：
  旧盘点见 [产物盘点](validation/sdevice_reference_artifact_inventory_2026-09-26.md) 第 6 节；
  用户已授权并完成 PN2D 0V/IV/BV 和 BJT m0 的最小 v2 重跑。WP10 已固化 A-4 与 C6，
  当前状态见 [WP10 收口报告](validation/sdevice_wp10_closeout_2026-10-03.md)；其余 WP14 案例未执行。
- A-5 C++ `MaterialDatabase` 内置 Si/SiO2 与本文第 2.4 节默认库的对齐及重新基线化：开发方案 WP13（用户 2026-09-29 要求列入）；
  完成前第 2.4 节的“遗留默认”说明保持有效。
- 惰性 `Thermode{}`（未求解温度方程）是否从 E105 降级为警告：transportmodels DD deck 的唯一阻塞项之一，
  阶段 1 维持拒绝。
- C5 新增键已由 WP6 回填第 2.6 节；C6 的 `time` 语义与符号因子由 WP10 回填第 3.3 节（C3 见第 2.4 节）。
