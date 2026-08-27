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

## 待完成门禁

- 将 M5 GRD/DAT 上传到已授权的 Sentaurus 2022 主机，执行 TDX 转换并回读 TDR 对比。
  该步骤需要对 M5 新文件单独授予出站授权；M4 授权不会被扩大解释。
