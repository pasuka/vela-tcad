# LDMOS 单位修正完整验收、常驻进程与 IALMob 计算核

2026-09-11，工作分支 `codex/templates-ldmos-phase-a`。
本报告接续[前轮执行](templates_ldmos_followup_execution_2026-09-11.md)。
单位修正后的 D5 双栅压完整验收通过；常驻进程数值等价但本次未提速，保持显式实验选项；
IALMob 独立计算核已实现，尚未完成局部对齐或接入自洽求解器。

## 1. 正确 Auger 单位的完整曲线

使用 `linked_d5_auger_units_inputs.json`，显式系数按配置的 cm⁶/s 语义为
Cn=2.9e-31、Cp=1.028e-31。两条曲线均从各自合格 Vd=0 种子开始，排除栅压预偏置。
保持可信预测最大 0.2 V、无预测最大 0.1 V、联动恢复和参考系平移原策略。

| 指标 | Vg=4 V | Vg=8 V |
| --- | ---: | ---: |
| 精确参考点 / 电压区间 | 31 / 0–40 V | 31 / 0–40 V |
| 非零推进 / 求解请求 | 251 / 253 | 251 / 253 |
| Newton 更新（包含参考系闭合） | 1139 | 1116 |
| 回退 | 0 | 0 |
| 请求累计墙钟时间 | 1088.165 s | 1154.914 s |
| 子进程累计 CPU 时间 | 844.547 s | 905.672 s |
| Id 相对误差中位数 / P95 | 1.43068% / 1.58331% | 1.33093% / 1.61287% |
| Ron / 40 V 端点误差 | 1.10245% / 1.44188% | 1.03212% / 1.09906% |

双栅压 40 V 电流比误差 **0.337945884635%**，原工程及最终联合评分均通过。
独立审计重计全部 trace/attempt/profiler 更新、原块/逐行/KCL 门限、状态及输入哈希。
这次是从零压完整推进，补齐前轮“62 点同偏压重闭合”的资格边界。

原历史系数全曲线的 Vg4/Vg8 更新数为 1144/1089，本次为 -5/+27 次。
这些完整曲线运行期间存在其他构建/验证负载，表中墙钟时间不可解释为单位修正导致的性能变化。

据此将受版本管理运行入口的默认 bundle 晋级为正确单位配置。
旧 `linked_d5_inputs.json` 和配置保持原值，通过显式 `--bundle` 用于历史回放。
尚未实现 Auger H/N0 密度增强；本次通过不代表所有 Sentaurus 模型差异已消除。

## 2. 常驻 DC 进程及输入准备缓存

实现位置：`DCSweep`、`vela_example_runner --dc-worker`、`scripts/dc_worker_client.py`、
`scripts/run_templates_ldmos_linked_d5.py --worker`。

- 一行 JSON 对应一个顺序求解请求，响应包含请求 ID、状态和捕获日志；单次失败不结束进程。
- 复用网格读取、box 几何、外部 AverageBox 耦合系数、材料和掺杂的准备结果。
- 每次比较相关配置及外部文件的完整字节和解析后的绝对路径；同名、同大小、同时间戳的文件变化也会失效。
- 每次复制准备模板后再解析接触、求解设置和显式种子；接受/拒绝的非线性状态不隐式跨请求继承。
- `DCSweep()` 默认不缓存；仅 `DCSweep(true)` 顺序调用启用，不用于并发共享。

### 同一冻结 Release 的前八点串行对照

固定 Vg=8 V、正确 Auger 单位、Vd=0–9.33333333333333 V 的前八个精确点。
先运行子进程入口，再运行常驻进程入口。二进制相同，参数和数值控制相同。
完整曲线作业已结束；这是一组顺序测量，未做反序重复或环境频率控制。

| 指标 | 逐请求子进程 | 常驻进程 |
| --- | ---: | ---: |
| 求解请求 / 进程 | 68 / 68 | 68 / 1 |
| 非零推进 / Newton 更新 / 回退 | 67 / 384 / 0 | 67 / 384 / 0 |
| 请求累计墙钟时间 | 290.4145 s | 316.8178 s |
| CPU 时间 | 250.0000 s | 258.0938 s |
| 控制器完整墙钟时间 | 319.9783 s | 355.6665 s |
| 输入准备 build | 16.0483 s | 0.2244 s |
| 输入字节检查 | 0 | 3.7062 s |
| SparseLU 数值分解 | 109.1059 s | 126.0627 s |
| Jacobian | 68.3436 s | 79.1154 s |

每次请求的目标电压、父状态哈希、更新数一致，八个精确点状态文件 SHA256 一致；
两组均通过独立原门限审计。常驻进程首次准备、后续 67 次缓存命中。
历史 `children` 字段在 worker 模式表示求解请求数；新增 `solver_requests` 与
`process_count` 明确区分请求和进程。阶段计时可能嵌套，不能相加当作总耗时。

本次请求墙钟时间 **增加 9.09%**，CPU **增加 3.24%**。
准备缓存节约未转化为整体提速：分解和 Jacobian 等求解阶段同时变慢。
仅凭这一组计时不能区分机器负载/频率变化与长期进程内存分配行为，不将变慢归因于已证实的根因。
结论：保留 `--worker` 为可复现实验入口，**不晋级默认，不宣称性能收益**。
后续若继续该方向，应反序/交错复测并记录 RSS、峰值内存与分配器行为，再决定是否保留复杂度。

### 高压和失败隔离

同一 worker 独立重放四个保存请求：26.666667 V 平移前、同偏压平移后、
26.766667 V、40 V；更新数分别 3/1/8/3，状态哈希均与完整曲线对应请求一致。
中间插入缺失种子的失败请求，后续三次仍通过，缓存命中且没有状态污染。
这验证高压参考系与异常隔离，**不是常驻入口完整 31 点曲线的资格**。

## 3. IALMob 独立计算核及局部探针

新增 `include/vela/physics/IalMobility.h`、`src/physics/IalMobility.cpp`、
`src/tools/ialmob_probe.cpp`。依据本机 T-2022.03 User Guide 第 408–414 页
方程 290–311，独立实现 3D/2D Coulomb、屏蔽/团簇、声子、粗糙度和界面距离衰减。
Nd、Na、n、p 分别输入；状态和输出边界采用 SI，系数采用手册 cm 单位。
G(P) 的最小值位置只在模型构造时计算一次。

范围明确限定为 **300 K、FullPhuMob、PhononCombination=1**：
没有应力/薄层修正，ClusteringEverywhere 关闭，声子/粗糙度中杂质加权系数及其子指数固定为默认 1。
该范围覆盖实际 par 中已审计的相应设置；温度指数保留但在 300 K 无作用。
独立晶向函数识别立方晶体 {100}/{110}/{111} 族，包含符号和轴置换对称性。
尚未把晶向选择、Nd/Na 和载流子依赖完整接入主求解器/Jacobian。

`ialmob_probe INPUT.json` 在 stdout 输出局部项，不运行曲线。
必须明确 `temperature_K`、`carrier`、`states_SI`；可用 `parameters_cm` 覆盖手册默认参数。
未知参数、错误单位键和非 300 K 输入拒绝运行。无限大的单项迁移率以 JSON null 表示无对应散射。
`crystal_normal` 只输出晶面族标签，**不自动选择参数组**；调用方必须传入已确认的组。
未把这个工具当作可用的 `solver.mobility=ialmob` 配置。

已完成迁移率测试 36 项、6375 条断言，其中新增 7 项 IALMob 测试，覆盖：
零场本征 SI 极限、远界面多数载流子掺杂曲线、界面场/距离变化、补偿掺杂、
异类载流子散射、物理范围正值/有限性、晶向对称性及错误输入。
单独探针协议检查覆盖参数单位与拒绝未支持模式。

### D4 原生保存状态的初步对照

从已有 D4-classical `n4_des.tdr` 导出 Vg=8 V、Vd=40 V 的原生字段，
核对节点坐标、300 K 温度及字段单位。采用实际 `Siliconc100.par` 的 {100} 参数。

- 严格平直、法向唯一的 {100} 界面只筛得 **2 个节点**，样本很小。
  按输出 Eparallel 估算高场修正后，电子相对偏差约 -10.61%/-11.94%，空穴约 -12.64%/-13.50%。
  这组估算尚不是实际高场驱动的精确重构，不用于计算核或曲线验收。
- 另筛得距 Si/SiO2 界面至少 0.2 μm（20 倍 l_crit）的 **2150 个体区节点**。
  第一次用 Eparallel 筛“弱驱动”出现长尾，随后查手册确认：deck 未显式指定时，
  HighFieldSaturation 默认驱动是 **GradQuasiFermi**；输出 Eparallel 不能直接代替。
- 再使用原生 QF 梯度，并筛载流子密度 ≥1e16 cm⁻³、无量纲 QF 驱动 <0.01，
  电子/空穴样本为 253/311 个。加入 QF 高场估算后的相对误差中位数约
  -1.90e-9 / -4.11e-4，但绝对相对误差 P95 仍为 **24.93% / 5.23%**。
  因此不能仅凭接近零的中位数宣布体区通过。

下一步需逐节点核对长尾区域、输出字段与迁移率求值的离散位置，以及
RefDens=1e12 的界面平行场插值。界面 signed Enormal 的处理也需明确后再扩大探针范围。
现有证据尚不能确定差异全部来自离散支持，也不能排除模型实现差异。
保持原 D5 关闭 IALMob 回归，通过局部对照及 Jacobian 测试后才开展 D4 自洽曲线。
现有本地手册和参数可继续使用，目前无需用户重复提供资料。

## 4. 环境、复现与证据

本次曲线均为 Windows UCRT64 **Release，-O3 -DNDEBUG，gprof 关闭**，
Eigen SparseLU/COLAMD（环境显式 `VELA_LINEAR_SOLVER=sparselu`）。
网格 10241 节点、19782 三角形，5723 个 Si 节点；16237 条外部 AverageBox 耦合，
barycentric 节点体积及 material_local Poisson 电荷体积。
300 K、Fermi、OldSlotboom、SRH/Auger；几何 μm，状态密度 m⁻³，端电流 A/μm。

完整曲线冻结 runner SHA256：
`1e7464bb5df420cf02b37527e19098ede0b0300bc64100b0a9137837f79b1331`。
常驻/子进程成对对照共同 runner SHA256：
`6a4a275cb8d565345b04800810978a14149eef6aa01c09313fd3ce76ac9dc14d`。
前者不含 worker 入口；后者包含 worker/cache，尚未包含后续新增独立 IALMob 工具。

从当前 worktree 根复现（输出目录必须不存在）：

```powershell
D:\msys64\ucrt64\bin\python.exe -X utf8 scripts/run_templates_ldmos_linked_d5.py --manifest reference_staging/templates_ldmos_followup_20260911/worker_binary/manifest.json --output reference_staging/unique_worker_vg8_first8 --gate 8 --points 8 --worker
D:\msys64\ucrt64\bin\python.exe -X utf8 scripts/audit_templates_ldmos_linked_d5.py reference_staging/unique_worker_vg8_first8
```

删除 `--worker` 即同二进制子进程控制；`--preflight` 仅检查输入。
默认单位修正 bundle 与早期冻结历史 bundle 必须明确区分。

本机证据根目录 `reference_staging/templates_ldmos_followup_20260911/`（不纳入 Git）：

- `vg4_corrected_auger_full/`、`vg8_corrected_auger_full/`、`corrected_full_joint_score/`。
- `cold_first8/`、`worker_first8/`、`worker_pair_comparison.json`、`worker_binary/`。
- `worker_frame_probe/summary.json` 和有意失败请求的原始响应。
- `ialmob_input_audit.json`、`d4_inventory.json`、`d4_export/`。
- `ialmob_local_probe_r2/`、`ialmob_bulk_probe/`，含逐节点 CSV、输入参数和探针输出。
- `ialmob_bulk_probe/driving_field_audit.json` 明确修正首轮 Eparallel 驱动假设。
- `ialmob_final_probe/manifest.json` 冻结最终探针及源码；四组已有局部输入重新运行，输出与前述对照逐项一致。
- `ialmob_tests.log`、`worker_ctest_full.log` 及最终 `final_build.log`、`final_ctest.log`。

最终按 `windows-ucrt64-release` 构建成功；完整 Release CTest **758/758 通过**，
124.63 s，包含 IALMob 计算核/探针、DC worker 协议、模板回归及 ASCII 检查。
`git diff --check` 无内容错误（仅 Windows LF/CRLF 提示）。
此前 worker/cache 阶段的 750/750 结果保留在独立日志中。
