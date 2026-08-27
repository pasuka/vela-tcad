# Sentaurus Gmsh M5 验收记录（2026-08-26）

更新日期：2026-08-27。

## 范围

M5 在 M3/M4 语义和导出闭环之上增加可选择的 Gmsh Transfinite Rectangle 网格。公开
默认仍为 `offline-static + gmsh`，网格策略默认为 `auto`；商业 Sentaurus/TDX 仅用于
内部门禁。

## 已通过证据

- PN2D `auto` 选择 `transfinite-rect`，对角布置为 `AlternateRight`。
- refinement 分区：x 坐标 `[0, 0.9, 1.1, 2]`，区间数 `[18, 20, 18]`；y 坐标
  `[0, 0.5]`，区间数 `[25]`。
- 网格：1482 节点、2800 个 Tri3、单连通、全部 CCW；Anode/Cathode 各 26 个节点。
- 每个 Tri3 均含直角，最大直角误差 `7.96e-13` 度；最大角
  `90.00000000000077` 度，非钝角门槛通过。
- `device.msh/.vtu/.vtk/.dat` 的 Donor、Acceptor、Net、Abs 和 Log10Abs 五个节点场相对
  bundle 的最大绝对差均为 0；GRD/DAT 和 Physical Groups 回读通过。
- 完整 PN2D IV 使用 Release 求解器覆盖 0 至 0.3 V，状态 `accepted`。
- 完整 PN2D BV 最终使用原子化 `element_edge_sg_gss_laux +
  element_vertex_box_measure + bernoulli + mixed_voronoi` profile 生成 406 个收敛点并覆盖至
  `-19.999999999999996 V`；最大端子电流混合容差比为 `0.9662744481`，小于门槛 1.0，
  状态为 `accepted`。
- 冻结 Sentaurus BV 曲线已纳入 E2E contract：比较 213 个参考点，最大差异
  `0.0560836635 decade`，最大相对误差 `0.1378464621`，均通过 `0.08 decade / 0.20`
  门槛；曲线 LF 归一化 SHA-256 为
  `9ede4233957834fbb42b37c0465f3d581c7e1e0c321f20335dc775de30cd225b`。
- PN2D BV 模板升级至 v4：移除默认 `contact_majority`，连续性行缩放禁用，carrier row
  和 global continuity closure 使用 `report`；最终 acceptance 仍强制执行端子守恒与
  Sentaurus 曲线合同。对应提交为 `69b259b`。
- `69b259b` 后的 UCRT64 Release 全量 CTest 为 712/712 通过，总耗时 798.37 秒；其中
  `sentaurus_device_e2e_contract`、`sentaurus_gmsh_m1` 至 `m5` 和
  `sentaurus_sample_integration` 均通过。
- M1-M5 定向 Release CTest 全部通过；增量 Release 构建成功，全量 CTest 712/712 通过。
- 四种对角布置 `Left/Right/AlternateLeft/AlternateRight` 均通过逐单元直角和非钝角验证。
- 局部 Rectangle refinement window 的轴向投影测试通过；人工分区切点不会截断原始完整
  contact 边。

## 内部 TDX/TDR live gate

- 授权日期：2026-08-27；目标：SSH `sentaurus`；远端隔离目录：
  `~/sentaurus_runs/vela_gmsh_tdx/pn2d_gmsh_m5_20260827`。
- Sentaurus Data Explorer：`T-2022.03-SP2`；`tdx -d device.grd device.dat device.tdr`
  退出码为 0，stderr 为空。
- 返回 TDR SHA-256：
  `6ae5439f0562c9c6180470a762cc9e80029361d76300cb5be43292ac367d0ad1`。
- 拓扑严格匹配：1482 节点、2850 总单元、2800 个体 Tri3；区域集合为
  `Anode/Cathode/R.Si`，Anode/Cathode 各 25 条接触边。
- 7 个数据集均比较 1482 个节点：Donor、Acceptor、Doping、AbsDoping、Log10AbsDoping、
  BoronActive 和 PhosphorusActive；每个数据集的最大绝对差和最大相对差均为 0。
- 机器可读清单位于忽略的 Release 产物目录
  `build-release/m5-tdx-live-20260827/tdx_gate_manifest.json`，状态为 `passed`。

## 已知边界

- 通用 Newton `carrier/global continuity enforce` 策略尚未作为 PN2D 默认求解策略完成
  数值资格认证；模板 v4 通过 report-only 求解和独立强制 acceptance 合同隔离该风险。
- `MaxTransDiff`、`MaxGradient` 等 M3 4B 能力仍保持 source-located fail-closed，不属于
  本次 M5 收尾范围。
