# Sentaurus Gmsh M5 验收记录（2026-08-26）

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
- 完整 PN2D BV 生成 406 个收敛点并覆盖至 `-19.999999999999996 V`；网格合同通过，
  但最大端子电流混合容差比为 `1.2157959514`，超过门槛 1.0，因此状态如实为 `failed`。
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

## 已知未通过项

- 完整 PN2D BV 已到达 -20 V，但端子电流混合守恒合同仍失败；不得把 TDX 文件兼容性
  门禁通过解释为 BV 数值合同通过。
