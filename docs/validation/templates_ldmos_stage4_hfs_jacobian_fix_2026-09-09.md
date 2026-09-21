# Stage-4 D5：二维 HFS 连续性 Jacobian 修复与验证（2026-09-09）

## 当前结论

二维 `transport_cell_vector` HFS 连续性 Jacobian 的迁移率反馈遗漏已修复。
最终候选第二版已通过 Release CTest **736/736**，同初态局部 HFS 开启控制的
更新数从 **79→4**、**130→11**，HFS 关闭控制的状态和迭代记录与旧版逐字节一致。
198 个 JVP 方向和两栅压低偏压检查已完成。两条独立完整曲线的状态见下方更新。
**目前不据局部结果宣布 D5 通过，不进入 IALMob。**

工作树为 `.worktrees/templates-ldmos-phase-a`，分支
`codex/templates-ldmos-phase-a`，基底 HEAD `959c1c91513c62caeb700c2f71d8cef194879419`。
保留交接时已有步长控制修改和所有历史证据，未提交 Git。

## 修复内容与判别性测试

- 在 `CoupledDDAssembler::assembleJacobian` 中对面积加权后的二维梯度矢量模长
  应用迁移率链式导数，覆盖端点及相邻输运单元的第三顶点；只计入输运材料。
  电子和空穴分别反馈到本载流子的连续性列，边两端贡献符号相反。
- 准费米参考值保持为常量坐标变换，梯度使用物理准费米势；梯度和列导数
  保留原 `fieldFactor`、势缩放路径。接触电场回退生效时仍使用 psi 模板，
  不重复加入 QF 反馈。关闭场导数时保留匹配冻结的 Jacobian 行为。
- 稀疏模板显式覆盖 vector HFS 的同载流子邻域列。第一版采用更宽的单元模板，
  引入的跨载流子结构零改变了 HFS-off 的排序／浮点求解路径；第二版收窄模板，
  已恢复下述 HFS-off 对照的逐字节一致性。第一版证据单独保留，不混入最终曲线。
- 新增两个 Catch2 用例，覆盖非沿边且非均匀的二维 QF 梯度、电子／空穴分块列、
  单独移除一条边以隔离第三顶点贡献、接触回退开／关，以及 Fermi–Dirac
  live／匹配 frozen 方向差分。分块归一化不使用单位下限，不由 Poisson 掩盖。
  修复前第一项测试有 12 个断言失败，其中第三顶点两载流子导数均遗漏 100%；
  最终两个测试的 64 个断言全部通过。

证据：[修复前测试](../../reference_staging/templates_ldmos_hfs_jacobian_20260909/test_before.log)、
[最终专项测试](../../reference_staging/templates_ldmos_hfs_jacobian_20260909/test_after_v2.log)、
[最终全量 CTest](../../reference_staging/templates_ldmos_hfs_jacobian_v2_20260909/ctest.log)。
全量测试包括 `pn2d_config_templates`、`pn2d_node_volume_policy_default_acceptance`、
`pn2d_bv_m2_mixed_voronoi_self_consistent_control`。

## 环境、物理与血缘

UCRT64 GCC 16.2.0，CMake `windows-ucrt64-release`，`-O3 -DNDEBUG`。
实际数值后端为 Eigen SparseLU/COLAMD、L2 行列平衡；配置同时检测到
HDF5-static、SuiteSparse SPQR 和 UMFPACK，可用不等于本算例实际启用。

- 旧 runner SHA-256：`9638411ec0d38434f05afc765d573060adde0d406c6818f3428c76777a066391`。
- 最终 runner SHA-256：`52af2e8a2a489f9be299ee45e1b0fdd5bd1ef5f4a61e43693b701ce22e428087`。
- 两个二进制均固定保存在独立 ignored 实验目录；使用 manifest 核验 SHA，
  不改写旧脚本的硬编码断言。[最终目录](../../reference_staging/templates_ldmos_hfs_jacobian_v2_20260909/)
  包含 `binaries.json`、`environment.json`、`source.diff`、控制脚本及状态血缘。

保持交接中的 10,241 节点／19,782 三角形二维网格、外部 AverageBox 输运耦合、
barycentric 体积、material-local Poisson 电荷及 legacy node-local 接触。
300 K，Fermi–Dirac、OldSlotboom、SRH/Auger；predictor、IALMob、雪崩、量子与热耦合关闭。
HFS 开为 `constant_field + transport_cell_vector`，关为 `constant`。
电压 V、几何 μm、密度 m⁻³、电流 A/μm。

不放宽原始门限：psi/e/h 块分别 `5e-8`／`1e-11`／`3e-10`，
局部 carrier-row `eps_row=1e-8`、违规数零，端口 KCL 比值 `<=1e-8`；
零偏单独保留既有评分器的规则。

## 同初态局部对照

Vg=8 V、Vd=4 V 的原合格 HFS-on／off checkpoint 在新旧 runner 上均无需更新
即可通过同偏压验收。每个电压控制都从完全相同的原 checkpoint 独立开始。
原始状态 SHA：

- HFS-on：`27d90ee9be95a323d431501d6a1e6dbc88d41479aba05a68e2612e44b9b3e8c1`。
- HFS-off：`4e83f6cd233929505ef65a1ca441ce0bc7386c94ebf70a31991a6336cfd33dda`。

| 控制 | 旧版更新数 | 最终修复版更新数 | 修复版 Id（A/μm） |
| --- | ---: | ---: | ---: |
| HFS 开，4→4.0025 V，QF 上限 0.0025 V | 79 | 4 | 2.3068573526169577e-4 |
| HFS 开，4→4.1 V，QF 上限 0.1 V | 129+1=130 | 10+1=11 | 2.3303736801698474e-4 |
| HFS 关，4→4.0025 V | 5 | 5 | 4.417251955549226e-4 |
| HFS 关，4→4.1 V | 8+1=9 | 8+1=9 | 4.520297630292511e-4 |

更新数包含失败的直接步及同偏压重收敛成本。HFS-off 两个最终状态和全部迭代
CSV 与旧版逐字节一致。100 mV 步的共有尾段重收敛需求仍存在，不能宣布其根因已解决。
控制期间与 CTest、其他工作树仿真存在并发，墙钟及 CPU 均保留为成本记录，
不作为隔离环境下的端到端加速结论。

## 全量 JVP 与低偏压资格

最终 runner 完成 HFS-off 90、HFS-on 90、匹配 frozen 18 个方向。
原来稳定的大电子反馈偏差在 `h=1e-7 V` 下为：节点 2951 `1.51858e-8`、
2949 `6.26824e-9`、3432 `3.27825e-10`，使用
`||Jv-FD||/max(||Jv||,||FD||)` 分块归一化。2951 的 `h=1e-5/1e-6/1e-7 V`
误差为约 `1.77119e-4/1.76955e-6/1.51858e-8`，与中心差分截断误差的二阶下降一致。
节点 5569 仍有差分幅度敏感性（电子 `h=1e-6/1e-7 V` 为
`1.50977e-4/7.49914e-6`），完整步长序列保留，不把它混同于原先稳定的大反馈遗漏。

两栅压零漏压旧 checkpoint 均在最终模型下零次更新通过；0→0.0025 V 的
Vg=8／Vg=4 控制分别需要 6+1／6 次更新，均通过原全局块、局部行与 KCL 门限。
Id 分别为 `2.883623841621566e-7`／`2.315620250270672e-7 A/μm`。
所有四个 4 V 起点控制的新旧 source/drain/gate/substrate 电流均完成核对；
差值按漏电流尺度归一化小于 `1e-8`，KCL 均不超过 `1e-8`。
详见 [局部审计](../../reference_staging/templates_ldmos_hfs_jacobian_v2_20260909/local_analysis.json)
及 [控制结果与初态血缘](../../reference_staging/templates_ldmos_hfs_jacobian_v2_20260909/local_progress.json)。

## 独立完整曲线状态

2026-09-09 09:08（Asia/Shanghai）已启动 `run_full_queue.py`，queue PID 20400、
Vg=8 controller PID 8372；实际进程与输出同时核验。顺序运行 Vg=8、Vg=4
独立 0→40 V 曲线，保留原 31 个精确点及 26.666667 V 处的坐标参考变换。
输出位于最终实验目录的 `vg8_full/`、`vg4_full/`，新旧二进制不混用。
当前完整曲线尚未完成，D5 最终资格仍待两条曲线和两栅压比评分。
实时状态：[队列](../../reference_staging/templates_ldmos_hfs_jacobian_v2_20260909/queue_progress.json)。
