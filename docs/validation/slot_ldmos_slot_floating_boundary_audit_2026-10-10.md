# SLOT 金属结构与浮置边界核查（2026-10-10）

## 结论

工艺脚本支持 SLOT 是金属结构：填充 Tungsten 后以接触替换金属。
但所检查的原生 IdVd 配置与原生保存态没有体现浮置金属等势条件；
BVDS 模板反而显式指定 SLOT 电压为 0 V。两者均不能当作已验证的浮置电极配置。
Vela 当前支持解析 floating 类型名称，但 Poisson、DC sweep 和诊断入口拒绝执行该类型。
因此物理设计意图与已有仿真边界存在需要澄清的差距，尚不能归因于 Vela 单方面遗漏，
也不能据此认定它是低密度行不收敛的原因。

## 证据范围

只读现有源文件、安装手册的本地提取副本和保存态；没有连接远程、启动求解、
修改物理配置或验收门槛。新增内容仅为本报告及 `build/s10slot_boundary_audit/`
中的逐节点统计、输入 SHA-256 绑定。HDF5 读取工具报告编译/运行版本不一致警告；
原生电势统计直接从已有 CSV 导出读取，不依赖本次 HDF5 读取。

### 工艺与原生电极声明

- 原目录：`D:/gpt-tmp/slot-LDMOS/V3E_NLDMOS_40V_Copy/`。
- `sprocess_fps.cmd:481`：填充 Tungsten。
- `sprocess_fps.cmd:489`：`contact name= SLOT point x= -0.4 y= 3.0 Tungsten replace`。
- `IdVd_des.cmd:14`：只声明 drain/gate/source/substrate；没有 SLOT 或其 Charge 条件。
- `BVDS_des.cmd:10`：`{ Name="SLOT" Voltage=0.0 }`。这是所读模板的固定电压声明，
  本次没有将结论扩展到所有历史 BVDS 运行的展开配置。
- 已实际运行的 `build/s09pmi_source/native/stdout.log:241` 电极清单也只有四个电极。
  所读原生运行版本为 T-2022.03，末态为 Vg=5 V、Vd=0.006 V。

### 手册与保存态

本地手册 `build/s09native_newton/sdevice_ug.txt`，T-2022.03，
印刷页 115–116、280：电压、电流和电荷是不同边界类型；浮置金属通过 Charge
指定总电荷，并保持表面电势恒定，以电通量积分与指定电荷相等确定浮置电位。
不能把省略电极声明等同于启用 Charge=0；零净电荷也不能仅由“浮置”二字推定。

从 `build/s10qualification_direction/native_final_export/contacts.csv` 的 SLOT 边集
核对：62 个节点构成同一个连通分量。按全局节点编号读取各材料区域的
ElectrostaticPotential CSV，62/62 节点均覆盖，单位 V。

| 状态 | SLOT 最小 ψ / V | 最大 ψ / V | 电势跨度 / V |
|---|---:|---:|---:|
| 原生已保存收敛末态 | 0.427280799097 | 4.542077053100 | 4.114796254003 |
| Vela 本批共同初态 P | 0.427339594167 | 4.542077704884 | 4.114738110716 |
| Vela A 第一步 | 0.427281044598 | 4.542076839109 | 4.114795794511 |

单个连通 SLOT 标记上的原生电势相差约 4.115 V，不是等势金属浮置边界的结果。
该证据来自实际场量，而非仅从没有 SLOT 的电极列表推测。
网格不包含 Tungsten 体区域；接触替换后 SLOT 节点邻接介质，
54 个仅邻接 SiO2、6 个仅邻接 Nitride、2 个邻接二者。介质邻接本身不否定金属设计意图。

### Vela 实现

- `docs/config_schema.md:1058` 说明 floating 在 Poisson 路径被拒绝。
- `src/simulation/PoissonSimulation.cpp:274`、`src/simulation/DCSweep.cpp:2807`、
  `src/tools/vela_example_runner.cpp:187` 均显式拒绝 Floating 类型。
- `src/solver/NewtonSolver.cpp` 的接触遍历对不在 contactBiases 中的名称直接跳过。
- 本批配置 `build/s10qualification_direction/prepared_initial_terms/config.json`
  仅含四个电极；已有 parent_poisson 行探针中 SLOT 62 行的 constrained 均为 0。

因此当前 SLOT 不是通过浮置电位未知量、等势约束和总电荷方程求解。
A 第一步的最大 93.043134 µV 变化来自这些自由电势节点，不能称为浮置电位整体移动，
也不是已配置四个电极的固定边界漂移。

## 后续边界

应先确认目标器件确实要求孤立金属，以及初始净电荷、是否存在额外电容耦合或泄漏路径。
随后在独立原生配置中显式建立浮置金属条件；Vela 需实现对应等势和电荷约束，
不能使用固定 0 V、仅提高介电常数或省略边界代替。
修正后属于新的物理配置，需从新冷态验证电荷、等势性及门/漏压路径，
旧链保留为原配置证据，不能直接合并资格。
