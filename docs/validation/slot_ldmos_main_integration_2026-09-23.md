# SLOT 合并本地 main：检查点与结构复用影响（2026-09-23）

## 范围与提交

本次从本地 `main` 合并，不 fetch 远端。合并前 SLOT 分支为 `a6930b39`，
main 为 `f6ff6cc517bbfb508036448266e72218bcd14601`。
目标是集成与兼容性分析，没有修改物理模型、收敛门槛或原生 Sentaurus deck。

| 提交 | 内容 | 对当前 SLOT 的作用 |
| --- | --- | --- |
| `b8e92343` | 清理旧 LDMOS 会议报告及图片 | 文档清理；保留本分支独立 SLOT 日报及验证索引，不删除本地 build 证据 |
| `f6ca2e3e` | CI 合并为单一 Linux 工作流 | 改善持续集成维护；不能替代本机 Windows Release 验证 |
| `45292d2e` | 生产检查点统一 HDF5、状态校验、性能实验基础设施 | 保留 binary64 与分离 QF 增量，增加网格/原点/单位/来源校验；同时引入旧续算脚本不兼容问题 |
| `cb51ce53` | 经典 DC 默认复用 Jacobian 稀疏结构 | 同一 DCSweep/worker 可复用装配位置；数值及所有偏导仍重新计算 |
| `f6ff6cc5` | 兼容较旧 HDF5 的 link 迭代 API | 降低旧 Linux/HDF5 环境构建限制，不改变离散物理 |

两处文档冲突已按内容整合：`docs/README.md` 保留双方新条目；
`docs/reports/README.md` 保留 SLOT 日报入口，同时接受 main 对旧 LDMOS 材料的删除。
源代码自动合并；Newton 接触局部掺杂分类及新源装配保留。

## HDF5：有助于可靠续算，但需要迁移驱动

当前契约见 [state_archive.md](../state_archive.md)。生产状态必须为 `.h5`、
`vela.state/2`；标准构建要求 HDF5/HighFive，TDR 导入仍是独立可选项。
网格身份、节点数、字段单位和势参考原点必须匹配。源配置/物理输入哈希用于
追溯，不等于跨模型的合格证明；显式初值允许来自不同偏压。

特别适合 SLOT 的两点：四个 QF 参考/增量数组独立保留，避免把小于参考值 ULP
的增量折回绝对 QF；状态先写同目录临时文件再替换，失败路径具有状态保护。
HDF5 本身不提高求解算术精度，不能消除当前 Poisson 舍入敏感度。

生产求解器不再回退读取 CSV/VDS/NPY。主线已迁移部分旧 SLOT 准备脚本，但
本分支较新的脚本仍依赖 CSV 状态，需要后续单独适配：

- `run_slot_ldmos_noavalanche_continuation.py`：状态输出、割线预测、carrier-only 修正。
- `run_slot_ldmos_loadline_continuation.py`：状态读取、确认态比较及输出路径。
- `slot_poly_equilibrium.py`：多晶硅重启状态的读取和改写。

因此不能直接以合并后的程序重启旧自动续算命令。应统一调用 `state_archive.py`，
保留源状态、校验 mesh/origin/单位、记录预测或修正来源，并用解码后的字段比较
确认态；不同 HDF5 文件的容器字节不必相同。普通 I-V、项分解、网格及掺杂 CSV
不属于状态迁移范围。本次没有增加生产 CSV 后备读取，也没有覆盖历史检查点。

## Jacobian 结构复用：作用及边界

`solver.reuse_jacobian_structure` 默认开启；显式设置优先于旧
`diagnostic_reuse_jacobian_structure`。共享对象只含零数值稀疏模式、索引和
edge/cell scatter；残差、迁移率及雪崩源偏导、数值 Jacobian 和数值分解仍更新。
结构键检查网格/节点顺序、边/单元支持、约束行、温度相关量和 stencil 等。

同一进程的连续 DCSweep/worker 才能保留跨请求缓存；现有 SLOT 外部驱动频繁
启动独立进程，跨点收益有限。应测命中次数和实际耗时，不能将 D5 的历史加速比
直接写成 SLOT 收益。Fermi 节点缓存仍默认关闭，且当前 Boltzmann 配置不受益。

这不是复用旧数值 Jacobian，不能代替完整偏导；也不会直接修复少数载流子逃逸、
HFS 负对角或接触/Poisson 残差问题。新源 10 V 或 BVDS 资格不能从该优化推导。

## 本地验证

证据目录：`build/slot_main_integration_20260923/`（生成物，不提交）。
合并前 runner 已独立保存，SHA256 为
`c5c7c6847bdc2cfee36d15a57955cdd691cfa3d63628c420f7188640ddaf5878`。
合并后 runner SHA256 为
`bd12915dfeef016aa0149a8cdedb325c3c3d86e5e8b02af06ed509d3487d2707`。
采用 Windows UCRT64 Release，实际配置启用 HDF5 静态库、HighFive、TDR 导入
及 UMFPACK；器件复验通过 `VELA_LINEAR_SOLVER=umfpack` 选择后端。

器件网格 13,881 节点，unit_scaling（长度内部单位 µm），状态密度 m⁻³，
本报告端电流 A/µm（原始 I-V CSV 为 A/m，乘 1e−6）。受控配置为 300 K、Boltzmann、元素 IALMob 完整偏导、
`ialmob_element_density_grad_qf` 雪崩、Si 内 SRH/Auger 及恢复掺杂的多晶硅。
不是旧投影 SG，也不代表原始完整 BVDS 已通过。

三个旧状态（2 V、2.125 V、10 V 未合格候选）已通过迁移工具转为独立 HDF5：
每份 10 个字段、13,881 节点，全部与 CSV 解析后的 binary64 值逐位一致。
输出各 1,120,720 bytes；原 CSV 约 2.05–2.10 MB。该比较针对数值，不代表
恢复 CSV 写入之前已经丢失的精度。见 `migration.json`。

本机 Python h5py 报告构建 HDF5 2.1.1、运行 HDF5 2.2.0 的版本警告；
不将此环境表述为完全匹配。迁移字段一致性已检查，C++ 读取及回归结果另列。

HDF5 输入的完整器件复验及独立全有效载流子行审计结果：

| 状态 | 电流 A/µm | ψ 块范数 | 最大载流子行比 | 结果 |
| --- | ---: | ---: | ---: | --- |
| 2 V | 9.264787215720656e−16 | 3.593676240276348e−10 | 5.393425429606717e−5 | 合格，电流保持 |
| 2.125 V | 9.589273548958036e−16 | 3.532364657664011e−10 | 2.082955501614698e−8 | 合格，电流保持 |
| 10 V 候选 | 不报告失败行零电流为物理结果 | 1.276409970746016e−9 | 2.398930310365803e−12 | 未合格，仍为 line_search_non_decrease |

三者全局闭合均通过，独立审计各覆盖 12,634 个有效电子行及 12,634 个有效空穴行。
10 V 接触多数载流子 QF drop 仍为 5.35167645226901e−11 V，高于 5e−11 V
保护阈值；没有跳过保护或调宽阈值。证据：`compatibility.json`、各 `verify_*`。

另外从同一 2.1 V 状态向 2.125 V 执行**恰好一次 Newton 更新**：
合并前 CSV 程序、合并后 HDF5 默认结构缓存、合并后显式关闭缓存三者的
10 个保存状态字段解码后全部逐位一致。性能计数确认每次实际装配一次 Jacobian、
接受一次 Newton 更新。因为人为限制一次迭代，三次均为未收敛候选，不能称为新合格点。
开启缓存的一次独立调用为 1 次 miss，无跨请求 hit；此实验支持装配及开关兼容性，
不证明 SLOT 的缓存命中收益。证据：`one_step_comparison.json`、`one_step_*`。

Release 全量构建成功；`ctest --preset windows-ucrt64-release -j 2`
**912/912 通过，146.95 s**。包括 HDF5 状态契约、Python/C++ 互操作、worker、
旧格式拒绝、缓存默认开关/失效、IALMob/新源 Jacobian 及混合接触保护回归。
见 `configure.log`、`build.log`、`ctest.log`。相关脚本回归通过不表示旧 CSV
续算流程已经适配新程序；上述器件验证使用独立 HDF5 配置，并未调用仍待迁移
的 CSV 状态改写链。文件链接与 `git diff --check` 检查通过。

## 对后续验证顺序的影响

1. 优先完成三类 SLOT 状态驱动的 HDF5 适配和对应回归，再恢复自动续算。
2. 沿用已验证 0.025 V 路径，整理从本地检查点开始的可重复续算；不能把分段
   carrier-only 修正后复验称为未经修正的冷启动连续扫描。
3. 在相同状态和门槛下测 worker 结构缓存收益，再决定是否改造驱动以减少进程启动。
4. 独立解决 10 V Poisson/接触保护阻塞，继续高压同模型场量及电流对照。
   旧投影 SG 的内部 37 V / 外部 37.1 V 资格不能沿用到新源。
5. 最后完成严格原模型、外电阻与击穿后支路验证；本次合并不改变其未完成状态。

合并前结果及物理问题详见
[接触保护与逃逸对照](slot_ldmos_contact_guard_and_escape_controls_2026-09-23.md)。
