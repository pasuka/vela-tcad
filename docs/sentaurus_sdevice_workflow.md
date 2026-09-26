# `vela device`：SDevice 命令文件直跑工作流（契约初稿）

状态：**WP0 契约冻结稿**（2026-09-26）。本文是
[开发计划](superpowers/plans/2026-09-26-sdevice-cmd-direct-run-development-plan.md)
第 4–6 节契约与第 10.0 节接口契约的可执行版本：schema 文件、夹具、消息目录以本文为准；
计划文档保留决策依据。尚未实现的部分（解析器、编排器、C++ 增量）以“待 WPn”标注。

## 1. 契约文件索引

| 契约 | 文件 | 生产者 → 消费者 | 状态 |
| --- | --- | --- | --- |
| C1 执行 IR v2 | `schemas/vela.sentaurus_execution_ir.v2.schema.json` | WP1 → WP3/WP4/WP5 | 冻结 |
| C2 网格包 | `schemas/vela.mesh_qualification.v1.schema.json` + 现有 `mesh.json`/`doping.csv` 契约 | WP2 → WP4/WP5 | 冻结 |
| C3 材料文件 | `schemas/vela.templates_ldmos.materials.v1.schema.json`（通用别名由 WP3 发布） | WP3 → WP4 | 沿用现有 |
| C4 stage 计划 | `schemas/vela.stage_plan.v1.schema.json` | WP4 → WP5/WP8 | 冻结 |
| C5 求解器协议 | `docs/config_schema.md` + `configs/schema/vela-simulation.schema.json` | WP6 → WP4/WP5/WP8 | 新增键待 WP6 |
| C6 输出文件契约 | 本文第 3 节 | WP8/WP6 → WP10 | 初稿 |
| C7 产物清单 | `schemas/vela.reference_artifact_manifest.v1.schema.json` | WP9 → WP10 | 冻结 |
| C8 CLI/退出码/消息目录 | 本文第 2、4 节 | WP0 → WP5/WP7/WP12 | 冻结 |
| C9 启动器协议 | 本文第 2.3 节 | WP7 → WP5 | 初稿 |
| C10 sidecar | `schemas/vela.sdevice_sidecar.v1.schema.json` | WP0 → WP4 | 冻结 |

校验工具：`python scripts/vela_sdevice/schema_validator.py <schema> <json>...`
（纯标准库的 draft-07 子集校验器；安装了 `jsonschema` 时单测会交叉比对两者结论）。
夹具：`tests/fixtures/sdevice/`（索引 `fixtures.json`，比较规则见第 5 节）。
契约单测：`tests/regression/test_sdevice_contracts.py`（CTest 名 `sdevice_contracts`）。

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

每个非 0 退出同时写 log 尾部 `[FATAL]` 摘要与 `<stem>_vela/fail_report.json{stage,message,details,exit_code}`。

### 2.3 启动器协议（C9，待 WP7）

`vela device <args>` → 子进程 `<python> bin/vela_sdevice.pyz device <args>`，cwd 不变，退出码透传。
`<python>` 查找顺序：`VELA_PYTHON` → `bin/python/` → PATH 中 `python3`/`python`（≥ 3.10，log 警告）。
Windows 用宽字符 API 透传参数；`bin/` 前置到 `PATH`。

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
  多段扫描的 `time` 语义（重置或累加）待 WP8 用 `bjt_m0_{base,collector}_*.plt` 与
  transportmodels `IdVgs_n6_des.plt` 确定后回填本节。电流符号因子：全局常量，取自首个通过
  第 12.4 节核对的 `pn2d_iv.plt`（待重跑），确定后连同文件哈希写入本节。
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
   `BreakCriteria` → `bv`；R3 单个电压扫描（其它） → `iv`；R4 多个 Quasistationary → `multi_stage`；
   R5 存在致命发现 → 仍按 R1–R4 推导（`verdict` 单独表达拒绝）。`reason` 字符串固定为
   `"R<n>: <本节描述>"`（见 golden）。
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

- A-1 VM 归档布局、A-2 `Parameter` 缺失处理、A-3 重跑时机、A-4 多段 `.plt` 的 `time` 语义：
  见 [产物盘点](validation/sdevice_reference_artifact_inventory_2026-09-26.md) 第 6 节。
- 惰性 `Thermode{}`（未求解温度方程）是否从 E105 降级为警告：transportmodels DD deck 的唯一阻塞项之一，
  阶段 1 维持拒绝。
- C3 通用别名 `vela.materials.v1`、C5 新增键、C6 的 `time` 语义与符号因子由对应 WP 回填本文。
