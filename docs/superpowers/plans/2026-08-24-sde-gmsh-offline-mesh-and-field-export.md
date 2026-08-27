# Sentaurus SDE 到 Gmsh 离线网格与场数据导出开发计划

## 计划状态

- 状态：已完成；M0-M5 公开 Gmsh→Vela IV/BV 工作流、本地 Release 回归和内部 live
  TDX/TDR 最终门禁均已执行并形成可审计结果
- 更新日期：2026-08-27
- 计划源文件：`docs/superpowers/plans/2026-08-24-sde-gmsh-offline-mesh-and-field-export.md`
- 实施 worktree：`D:\code-repo\copilot-worktrees\vela-tcad-dev`
- 实施分支：`codex/sentaurus2022-sde-workflow`
- 审计起点：`efb590f`，当时分支相对 `main` 落后 9 个提交、独有 3 个提交。
- 同步结果：计划提交 `9d0c77a` 后，通过 merge commit `e1987ef` 合入本地 `main`
  `f66dbd1`；`main` 现为实施分支祖先，既有分支提交均保留，并已完成合入后的 Release
  基线回归。
- 最终收尾同步：M5 文档收尾提交为 `605cd6c`；经确认 234 项工作区差异全部仅为
  CRLF/LF 行尾差异后恢复至索引版本。随后通过标准双父 merge commit `af01ea2` 合入
  本地 `main` `270f525` 的 4 个后续提交，合并无内容冲突，分支相对 `main` 为
  0 个落后、17 个独有提交。合并后 UCRT64 Release 增量构建 104/104 成功，全量
  CTest 712/712 通过，耗时 795.58 秒。
- 阶段 -1 验证：目标 worktree 独立 Release 配置和 137 个构建步骤成功；Sentaurus 定向
  回归 10/10 通过；全量 CTest 705/705 通过，耗时 246.90 秒。该结果是 M0 之前的冻结
  基线，后续新增失败必须相对此基线归因。
- M0 结果：已新增 Gmsh 4.15.2 Python API probe、UCRT64 构建目录 runtime staging、
  严格 JSON fail report、Device IR v1-to-v2 失真感知迁移、Device IR v2 schema 和
  `vela.mesh_field_bundle.v1` schema。真实 probe 验证 wrapper/API/CLI 均为 4.15.2，
  实际加载 staged `gmsh-4.15.dll`，并生成 3 节点、1 个 Tri3；拓扑相关选项包含单线程、
  `Mesh.Algorithm=6`、`Mesh.RandomFactor=0`、禁用重组/高阶/优化及 MSH 4.1 ASCII。
  UCRT64 Release 构建成功，最终全量 CTest 707/707 通过，耗时 234.09 秒。
- M1 结果：已实现独立的公开离线预处理 CLI、严格分离 Rectangle 的 OCC/Tri3 建模、
  region/material/contact Physical Groups、constant `NoReplace` 节点掺杂以及单一 bundle
  驱动的 MSH/VTU/VTK/GRD/DAT/mesh.json/doping.csv 导出。所有格式均由独立读回路径与
  bundle 比较；MSH 严格包含一份网格段和五份 NodeData。重叠、接触界面、ABA/BAB、
  insert-vertex、Gaussian、replacement 与 refinement 继续 fail-closed。验收说明见
  `docs/validation/sentaurus_gmsh_m1_acceptance_2026-08-25.md`。最终 UCRT64 Release
  构建成功，全量 CTest 708/708 通过，耗时 209.44 秒。
- M2 结果：已实现 per-body ABA/BAB OCC fragment ownership、source-map 与几何包含双重
  审计、创建顺序界面命名、straight-edge insert-vertex、完整 split-edge contact、曲线
  邻接 owner、有序 NoReplace/LocalReplace/Replace 和 line-baseline Gaussian。M2 节点场
  按 Sentaurus region-local 值在共享节点求和，不使用 tag、centroid、last-owner 或单元
  majority。T-2022.03-SP2 官方 oracle 退出码为 0，得到 7 点/5 单元、`R.P=1 um2`、
  `R.N=2 um2`、`R.P+R.N` 界面，并将 Positive 左法向、共享节点与 placement 参数顺序
  冻结到机器可读 golden。验收说明见
  `docs/validation/sentaurus_gmsh_m2_acceptance_2026-08-25.md`。UCRT64 Release 构建
  成功；全量 CTest 首轮 708/709 通过，唯一失败为新增 provenance 引起的旧 PN2D IR
  golden 漂移，更新受哈希封印的 golden 后该 CTest 独立复跑通过，因此 709 项均有
  Release 通过证据。
- M3 结果：阶段 4A 已完成。Line/Rectangle/闭合逆时针 Polygon Ref/Eval window、
  region/material/window placement、各向异性 max_x/max_y、重叠 MinAniso 最细方向约束及
  MaxLenInt（含 DoubleSide、UseRegionNames 和 factor）均生成可审计 Gmsh field DAG。
  各向异性网格遵循 Gmsh 4.15.2 官方 t17，使用 BAMG/Algorithm 7；min spacing 只作为
  Sentaurus refinement 停止目标记录，不作为 Gmsh 硬下限。完整 PN2D 代表性运行得到
  1044 点/1954 个 Tri3，全局和结区 P95 边长分量通过 1.35 倍门槛，全部格式节点场差异
  为 0。MaxTransDiff/MaxGradient 作为 4B 继续 source-located fail-closed。验收说明见
  `docs/validation/sentaurus_gmsh_m3_acceptance_2026-08-25.md`。UCRT64 Release 构建
  成功，全量 CTest 710/710 通过，耗时 303.47 秒。
- M4 阶段结果：`sentaurus_import.py device` 已以 `offline-static + gmsh` 为公开默认，
  完成 2×3 前端/网格后端组合矩阵、Gmsh 专属 barycentric/钝角质量合同、完整跨格式产物
  向 Vela deck 的单命令传递以及独立的内部 TDX gate。完整 PN2D IV 使用真实 Gmsh 网格
  得到 16 个收敛点并覆盖 0–0.3 V，状态 accepted；完整 BV 得到 406 个收敛点并到达
  −20 V，但当时端子电流混合守恒比 1.043 超过门槛 1，状态如实为 failed；该阶段性
  缺口已由 M5 后续的原子化 SG/GSS-Laux BV profile 修复。TDX gate 的
  hash-sealed dry-run 已通过；获明确授权后，live TDX T-2022.03-SP2 转换退出码为 0，
  返回 TDR 的 1062 节点、2010 总单元、1990 个体 Tri3、两个各含 10 条边的 contact 及
  7 个节点掺杂数据集均通过严格对比，场值最大绝对/相对差为 0。M4 影响范围 Python
  回归 114/114、UCRT64 Release 全量 CTest 711/711 通过（157.26 秒）。验收记录见
  `docs/validation/sentaurus_gmsh_m4_acceptance_2026-08-25.md`。
- M5 最终结果：公开 Gmsh 路径新增 `auto | transfinite-rect | unstructured` 网格策略及
  `Left | Right | AlternateLeft | AlternateRight` 对角布置。`auto` 只对单个轴对齐矩形、
  Rectangle refinement window、无 refinement function/insert-vertex 的 Device IR 选择
  Transfinite；其他合法 M3 输入回退非结构 BAMG，显式请求不满足条件时 fail-closed。
  PN2D 的 0.9/1.1 um 结区窗口被分为 18/20/18 个 x 向区间，y 向为 25 个区间，实测
  1482 节点、2800 个 Tri3、两端各 26 个接触节点；每个单元均有一个直角，最大角
  90.00000000000077 度，MSH/VTU/VTK/DAT 五个节点场差值为 0。完整 IV 已 accepted。
  完整 BV 随后使用原子化 `element_edge_sg_gss_laux + element_vertex_box_measure +
  bernoulli + mixed_voronoi` 组合重新验收：406 个点到达 -20 V，端子电流混合容差比
  0.966，小于门槛 1.0；与冻结 Sentaurus BV 曲线比较 213 个点，最大差异 0.0561 decade、
  最大相对误差 13.78%，最终状态 accepted。参考曲线、LF 归一化 SHA-256 和误差门槛已
  正式写入 E2E contract，提交为 `69b259b`。M1-M5 定向 CTest 全部通过；原 M5 基线的
  UCRT64 Release 全量 CTest 712/712 通过。提交 `69b259b` 后再次执行全量 Release
  CTest，712/712 通过，总耗时 798.37 秒。
  经单独授权后，内部 TDX T-2022.03-SP2 live gate 退出码为 0；返回 TDR 的 1482 节点、
  2850 总单元、2800 个体 Tri3、Anode/Cathode 各 25 条接触边以及 7 个节点数据集全部与
  GRD/DAT 严格匹配，所有场的最大绝对/相对差均为 0。
- 产品默认：`offline-static + gmsh`
- 内部验证：`sentaurus2022-evaluated + tdx`，仅限具有许可证的开发环境
- 取代范围：本计划取代“离线生成器使用正交张量网格并固定对角线切分三角形”的网格部分；
  已有 SDE/SDevice 解析、Device IR、求解和 IV/BV 验收能力在核实后复用。

## 背景与架构纠正

目标不是让 Vela 自行维护一个简化版二维网格器，而是读取 Sentaurus SDE `.cmd`
脚本，解析几何、区域、材料、接触、网格控制和掺杂定义，再通过 Gmsh Python API
完成建模与三角形网格生成。生成后的节点应被赋予掺杂浓度，并同时导出可供 Gmsh、
ParaView/VisIt、Vela 和 Sentaurus/TDX 工具链使用的结果。

当前 `scripts/sentaurus_mesh_builder.py` 构造正交张量网格，并把每个矩形沿固定对角线
切成两个直角三角形。该实现不满足本目标，后续只保留为显式兼容后端
`tensor-legacy`，不得继续作为默认离线实现或 Gmsh 验收替代品。

## 本次代码与环境审计结论

以下结论来自实施 worktree 的当前代码和本机环境，不得在实施时重新降级为假设：

1. 当前公开入口已经使用 `scripts/sentaurus_import.py device`，前端参数名称是
   `--sde-frontend`，默认值已经是 `offline-static`。后续必须保持该名称兼容，不能另造
   `--frontend` 破坏现有命令。
2. 当前离线路径直接调用 `sentaurus_mesh_builder.build_mesh_and_doping()`，输出
   `mesh.json/doping.csv`，并将来源写成 `offline-static-approximation`。Gmsh 后端必须替换
   这条默认分支，而不是在旁边增加一个无人调用的脚本。
3. 当前 evaluated 前端只要发现 `official_mesh_tdr` 就隐式优先导入官方网格。新架构必须
   将前端和网格后端正交化；选择 evaluated 前端不能隐式改变 `--mesh-backend`。
4. 当前 Device IR schema 是 `vela.sentaurus_device_ir.v1`。本计划需要增加 body 创建态、
   contact entity、refinement 和 doping provenance，因此首个不兼容版本固定升级为 v2，
   不允许继续向 v1 偷加必选字段。
5. 当前 `mesh.json` 只保存 contact node set，不能表达 contact line element。必须引入一份
   exporter 共用的规范化 mesh/field bundle；`mesh.json` 只是从该 bundle 派生的 Vela 视图。
6. 本机 UCRT64 环境为 Python 3.14.7、Gmsh 4.15.2。`gmsh.exe -version` 可用；默认
   `import gmsh` 失败。系统存在 `D:\msys64\ucrt64\lib\gmsh.py` 和
   `D:\msys64\ucrt64\bin\libgmsh.dll`，但 wrapper 寻找 `gmsh-4.15.dll`，即使设置
   `PYTHONPATH` 仍无法初始化。阶段 0 必须先解决并自动检测这项包装布局问题。
7. `scripts/pn2d_minimal6_topology.py` 已有固定拓扑 DF-ISE writer/parser 素材，但它不能
   直接扩展成通用实现；通用 writer 必须以规范化 bundle 为输入并单独测试。

## 已冻结的架构决策

- **前端与网格后端正交**：`--sde-frontend` 只决定如何得到 Device IR；
  `--mesh-backend` 只决定如何得到网格。支持值为公开默认 `gmsh`、显式兼容
  `tensor-legacy` 和内部专用 `official-tdr`。`official-tdr` 仅在 evaluated manifest 提供
  官方网格时可用。
- **IR 明确升级**：新增 `vela.sentaurus_device_ir.v2`；提供独立的 v1-to-v2 迁移函数，
  迁移无法恢复的 provenance 必须标为 unknown 并在需要该信息的后端 fail-closed。
- **单一规范化数据源**：新增 versioned `vela.mesh_field_bundle.v1`，保存原始 Gmsh tags、
  连续 Vela/DF-ISE IDs、Tri3、边界 Line2、region/material/contact/interface membership、
  节点 species/donor/acceptor/net fields 和单位。所有 exporter 只消费此 bundle。
- **稳定命名与编号**：Physical Name 使用命名空间 `region:<name>`、`material:<name>`、
  `contact:<name>`、`interface:<left>|<right>`；Physical tag 由排序后的 canonical name 和
  固定保留区间分配，并写入 manifest，不使用 Python hash 或 OCC 返回顺序。
- **确定性范围**：第一里程碑固定 Gmsh 4.15.2、MSH 4.1 ASCII、单线程、固定随机种子、
  一阶 Tri3。跨平台不要求字节相同，但同一环境重复运行必须满足 ID 映射和语义不变量。
- **质量与求解解耦**：通用预处理始终计算并报告钝角；是否拒绝由显式
  `--mesh-quality-profile` 或 E2E contract 决定。不得因 Gmsh 产生钝角而暗中修改
  PN2D BV 原子 solver profile。
- **商业工具隔离**：任何公开命令、测试和默认路径都不能探测或调用 `sde/tdx`；
  `official-tdr` 与 TDX round-trip 只存在于内部验证 job。

## 总体目标

提供一条无需商业 Sentaurus 运行时的公开命令：

```text
Sentaurus SDE .cmd
        -> offline-static 解析
        -> versioned Device IR
        -> Gmsh OCC 几何与 Physical Groups
        -> Gmsh 二维一阶三角形网格
        -> 节点掺杂求值
        -> mesh.json / doping.csv
        -> device.msh（含 Physical Groups 和 NodeData）
        -> device.vtu / device.vtk
        -> device.grd / device.dat（DF-ISE 文本）
        -> manifest、编号映射和质量报告
```

公开离线模式直接写出兼容 DF-ISE 的 `.grd/.dat`，不调用商业 `tdx`。内部许可证环境
使用 `tdx -d device.grd device.dat device.tdr` 验证兼容性，不能把内部验证命令变成
公开运行时依赖。

## 产品与安全边界

1. 默认前端为 `offline-static`；`sentaurus2022-evaluated` 只能显式启用，并在帮助、
   manifest 和文档中标记为内部开发验证模式。
2. 默认网格后端为 `gmsh`；旧实现只能通过 `--mesh-backend tensor-legacy` 显式选择。
3. 无法精确解释的 SDE 命令、布尔模式、接触选择、掺杂语义或网格控制必须
   fail-closed，并报告源文件、行号、命令和未支持能力。
4. 不允许把未知 SDE 布尔模式近似成 ABA，也不允许按接触邻接单元数量猜测 owner。
5. 不要求重新实现完整 Elk Scheme、ACIS、Sentaurus Mesh、三维几何、process
   emulation 或任意外部 TDR profile。

## 目标输出合同

### `device.msh`

- 使用 Gmsh MSH 4.1。
- 保存几何实体、二维三角形、边界线、Physical Names 和 Physical Groups。
- 二维 Physical Groups 表示 SDE region/material；一维 Physical Groups 表示 contact
  和需保留的材料/区域界面。
- 使用模型关联的 `$NodeData` 保存以下标量场，单位均记录在 manifest：
  - `DonorConcentration`，`cm^-3`；
  - `AcceptorConcentration`，`cm^-3`；
  - `DopingConcentration`，有符号净掺杂，`cm^-3`；
  - `AbsDopingConcentration`，派生可视化量；
  - `Log10AbsDopingConcentration`，派生可视化量，并在 manifest 记录零值截断下限。
- MSH 中使用 Gmsh 原始 node/element tag，不假定 tag 从 0 开始或连续。

### `device.vtu` 与 `device.vtk`

- `device.vtu` 为推荐的 VTK XML 输出；`device.vtk` 为 Legacy ASCII 兼容输出。
- `PointData` 至少包含与 MSH 一致的五个掺杂字段。
- 三角形 `CellData` 至少包含 `RegionId`、`MaterialId`、
  `PhysicalSurfaceTag` 和 `ElementId`。
- 接触边界作为 line cells 输出，并包含 `ContactId` 和 `PhysicalCurveTag`；非适用单元
  使用明确的哨兵值，不得丢失接触信息。
- 字段名、单位和符号约定与 MSH、DF-ISE 完全一致。

### `device.grd` 与 `device.dat`

- `.grd` 保存顶点、唯一边、带方向边引用的三角形、region/material ownership、
  contact 边界元素和 Locations。
- `.dat` 保存 `DopingConcentration`、donor/acceptor 总量，以及 IR 中能保留物种名称的
  active species 数据集。
- DF-ISE 使用从 0 开始的连续顶点、边和元素编号。
- 现有 `scripts/pn2d_minimal6_topology.py` 的 writer 只作为格式验证素材；实现必须抽取
  成与 PN2D 节点数、坐标和区域名无关的通用 writer。

### Vela 中间文件

- `mesh.json`：从 Gmsh 网格规范化得到的 Vela 网格。
- `doping.csv`：每个连续 Vela node ID 的 donor/acceptor，单位 `cm^-3`。
- `id_mapping.json`：Gmsh tag、Vela ID、DF-ISE ID 的双向映射及哈希。
- `mesh_qualification.json`：拓扑、方向、角度、退化单元、区域、接触和网格策略检查。
- `run_manifest.json`：输入和输出哈希、Gmsh/Python/Vela 版本、全部 Gmsh options、前端、
  后端、单位、能力清单和生成状态。

## 实施阶段

## 阶段 -1：同步实施基线并落地计划

### 任务

- 将本计划纳入实施分支，确保后续提交和验收报告引用同一份版本化方案。
- 在不丢弃实施分支 3 个既有提交的前提下合入最新 `main`；记录 merge base、合入方式和
  冲突解决结果。禁止用 hard reset 覆盖任一侧工作。
- 合入后重新审计 `sentaurus_import.py`、SDE frontend、mesh builder、TDR importer 和
  PN2D IV/BV contract，更新基线测试数及已知失败清单。
- 使用 Release 配置建立干净基线；不得复用另一个 worktree 的 build 目录或未提交产物。

### 验收

- 实施 worktree 除计划落地提交外保持干净，`git diff --check` 通过。
- 最新 `main` 与实施分支的既有 3 个提交均可从提交图追溯。
- Release 构建、现有 Sentaurus Python regression 和全量 CTest 通过；若存在基线失败，
  必须在首个实现提交前形成机器可读豁免清单，不能归因于后续 Gmsh 修改。

## 阶段 0：冻结基线与配置 Gmsh Python API

### 任务

- 在实施目标分支上重新审计现有 SDE IR、geometry/contact/refinement/doping 支持，
  不以旧报告代替代码核查。
- 新增 `scripts/probe_gmsh_python.py`；成功条件不是简单 `import gmsh`，而是能够
  `initialize()`、创建模型、生成一个最小三角形、查询 nodes/elements 并 `finalize()`。
- 修复已确认的 UCRT64 wrapper/DLL 布局问题，但不得修改系统安装或把 26 MB DLL 提交
  到仓库。允许的方案是构建目录内的受控 runtime staging：复制 `gmsh.py` 与共享库到
  ignored 目录，并使用 wrapper 期望的 DLL 名；记录源路径、版本和 SHA-256。也可使用
  已验证兼容 Python 3.14 的官方 4.15.2 wheel，但两种布局必须通过同一 probe。
- 增加统一的 Gmsh 能力探测和清晰错误信息。缺少模块时不得静默回退
  `tensor-legacy`。
- 首轮只接受 Gmsh API 4.15.2；CLI 与 Python API 版本必须一致。固定 MSH 4.1 ASCII、
  二维算法、随机种子、线程数、重组/高阶和优化选项。
- 在 manifest 保存 Gmsh 版本和全部影响拓扑的选项；跨 Gmsh 版本的网格连接变化
  不视为字节级回归，但必须满足几何和质量不变量。

### 验收

- 有 Gmsh 环境能生成单矩形 Tri3 网格。
- 无 Gmsh 环境以非零退出并生成机器可读 fail report，包含安装/配置建议。
- 不产生未声明的 CLI 后端回退。
- probe 同时验证 Python wrapper 版本、实际加载的共享库和 `gmsh.exe` 版本，三者不一致
  时 fail-closed。

## 阶段 1：扩展并冻结 Device IR 合同

### 任务

- 为每个 body 保存创建顺序、源位置、创建时有效的 boolean mode、region 和 material。
- 明确支持的二维形状：第一里程碑至少为 Rectangle、Polygon 和作为 refinement/
  Gaussian baseline 的 Line；其他形状 fail-closed。
- contact IR 保存 contact set、显示参数、current contact、pick point/entity list、
  `insert-vertex` 和源位置。
- refinement IR 保存 window 类型、placement target、max/min spacing、函数和接口参数。
- doping placement 保存 target 类型、replacement、decay、方向、顺序、物种和源位置。
- 输出 `vela.sentaurus_device_ir.v2` JSON schema、v1-to-v2 migration 和 capability
  manifest；v2 中每个 provenance 字段都必须有 schema 约束。

### 验收

- IR golden fixture 覆盖 ABA、BAB、局部 contact、region/material/window placement、
  replacement 和 line-baseline Gaussian。
- 未支持命令报告精确到文件和行号。
- v1 fixture 经迁移后可以继续由 `tensor-legacy` 消费；需要 v2 provenance 的 Gmsh
  操作遇到 unknown 时精确 fail-closed。

## 阶段 2：Gmsh OCC 几何与区域/材料 Physical Groups

### 任务

- 新增 `scripts/sentaurus_gmsh_backend.py`，所有 Gmsh API 调用集中在该模块。
- 将支持的 SDE body 转成 OCC 实体；使用 OCC fragment 建立共形、无重叠二维 faces。
- 必须消费 `occ.fragment()` 返回的 source-to-fragment 映射建立候选 ownership，再用面积
  和几何包含检查复核；不得用 surface tag 大小、返回数组位置或单一质心猜测归属。
- 在 fragment 后按创建顺序和每个 body 的 boolean mode 对碎片做所有权分类：
  - ABA：新 body 保持完整，从所有旧 body 删除重叠部分；
  - BAB：旧 body 保持完整，从新 body 删除与旧 body 重叠部分。
- 对无法唯一分类、零面积或非流形碎片 fail-closed；不得依赖 OCC 返回 tag 顺序。
- `occ.synchronize()` 后，为每个 region 建立二维 Physical Group；另保存 material 映射。
- 同一 region 的多个碎片属于同一 Physical Group，不能错误合并跨 region 表面。
- 对同一 material 的多个 region 另建 material Physical Group；region 与 material 双重
  membership 必须能从 MSH 4.1 Entities/PhysicalNames round-trip。

### 验收

- ABA/BAB 双向部分重叠、完全包含、多旧区域覆盖和零面积接触 fixture 通过。
- 面积、区域数、材料归属和内部共形界面均有数值 golden 检查。
- Gmsh 重新打开 `.msh` 后 Physical Names、稳定 tag 和实体成员保持一致。

## 阶段 3：接触、边界和接口 Physical Groups

### 任务

- 在 OCC 层执行 `insert-vertex` 对直线边的真实分割，再进行 contact 选择。
- `set-contact` 的 pick point 必须唯一落在一条完整 curve 上；端点、多边重合或内部界面
  歧义必须 fail-closed。
- 支持显式 contact 名称和 current contact；名称处理遵循已冻结的 SDE 语义。
- 每个 contact 建立一维 Physical Group，并保存完整 curve，而不是只保存采样节点。
- region/material 界面建立稳定的内部 curve group，供 `MaxLenInt` 和验证使用。
- Gmsh 网格提取后，从 Physical Curves 导出接触 line elements 和有序节点集合。
- contact owner 由 curve 邻接的器件 surface membership 决定；内部界面、多个邻接器件面或
  无邻接面必须 fail-closed，禁止恢复现有的“接触节点邻接单元多数”推断。

### 验收

- `insert-vertex + set-contact` 只绑定分割后的目标完整边段。
- 多接触不共享边界元素；共享端点是否允许由明确合同控制。
- contact curve、line elements、Vela contact nodes、VTK line cells 和 DF-ISE contact
  elements 相互一致。

## 阶段 4：SDE refinement 到 Gmsh mesh fields

### 任务

- 支持 Line、Rectangle、Polygon Ref/Eval window。
- 支持 region、material 和 window refinement placement。
- 将区域尺寸、窗口尺寸、距离/阈值和重叠最小尺寸组合成可追踪的 Gmsh field graph。
- 阶段 4A 先实现可直接对应几何的 region/material/window size 与 `MaxLenInt`，包括
  `DoubleSide` 和 `UseRegionNames`，并输出可审计的 field DAG。
- `MaxTransDiff` 和 `MaxGradient` 列入阶段 4B：必须先从 Sentaurus 2022 手册冻结变量、
  变换和阈值公式。若 Gmsh 原生 field 无精确等价，则实现“试网格 -> 节点掺杂/梯度 ->
  PostView background field -> 重网格”的有界迭代，并记录每轮收敛指标；不得用全区域
  最细尺寸冒充函数支持。
- max/min spacing 的定义遵循 SDE 语义；min spacing 不能误当成绝对不可突破下限。
- 对各向异性 max_x/max_y 只有在已验证 metric 映射后才宣称支持；在此之前遇到不等值
  必须 fail-closed，不能静默取两者最小值近似。
- 强制二维一阶三角形，不允许 recombination 或高阶单元。

### 验收

- refinement fixture 比较目标区域、窗口和界面附近的实际边长分布。
- 重叠 placement 的有效最大尺寸使用最细约束。
- 未实现的 MaxInterval、multibox、offsetting 明确 fail-closed 并进入 capability report。
- 4B 的迭代达到最大轮数仍不满足边长/场变化门槛时返回 failed，不输出成功 manifest。

## 阶段 5：网格提取、连续编号与质量门槛

### M5 补充目标：可选择 Transfinite Rectangle 直角三角形网格

- `sentaurus_gmsh_pipeline.py` 的 M5 默认策略为 `auto`；设备入口保持公开默认
  `offline-static + gmsh + auto`。
- `transfinite-rect` 将 Rectangle refinement window 边界作为结构化分区线，对每段按
  所有适用 region/material/window placement 的最小有效 max spacing 计算区间数，并以
  `ceil(length/spacing)` 冻结节点数。
- 每个矩形 tile 调用 Gmsh `setTransfiniteCurve` 和 `setTransfiniteSurface`；固定一阶
  Tri3、`Mesh.RecombineAll=0`，禁止生成四边形或高阶单元。
- 对角布置提供 `Left`、`Right`、`AlternateLeft`、`AlternateRight`；默认
  `AlternateRight`。所有布置均须通过“每个 Tri3 有且仅需至少一个 90 度角”和非钝角
  门槛，不能只凭命令名宣称为直角网格。
- 输出 `transfinite_mesh_contract.json`，记录分区坐标、每段区间数、预期/实际节点与
  单元数、对角布置、重组和单元阶次；同样信息进入 bundle provenance、资格报告和
  run manifest；格式由 `schemas/vela.gmsh_transfinite_mesh_contract.v1.schema.json` 冻结。
- 第一版资格范围刻意保守：多 body、Polygon body、insert-vertex、Line/Polygon
  refinement window 或 refinement function 不得进入 Transfinite 路径。`auto` 对其中
  仍属 M3 能力范围的输入回退 `unstructured`；显式 `transfinite-rect` 返回含原因列表的
  fail report。

### M5 补充验收

- PN2D `auto` 必须选择 `transfinite-rect`，节点/单元数严格为 1482/2800，接触节点不少于
  11，最大角不超过 90.000001 度且全部 Tri3 含直角。
- PN2D refinement 分段必须为 x=[18,20,18]、y=[25]；不得以全域最细网格替代窗口约束。
- MSH、VTU、VTK、GRD/DAT 和 CSV 的节点掺杂字段继续逐节点一致。
- 完整 PN2D IV/BV 必须读取策略感知 E2E contract；Transfinite 使用
  `mixed_voronoi + require_non_obtuse`，非结构 Gmsh 保持 `barycentric`。
- 在已授权的 Sentaurus 2022 内部环境执行 TDX 转换，回读 TDR 后核对拓扑、接触、区域
  和全部节点场；公开运行时不得依赖 TDX。

### 任务

- 使用 `gmsh.model.mesh.getNodes/getElements` 提取原始 tag、坐标和 Tri3 connectivity。
- 保留 `.msh` 原始正整数 tag；另生成从 0 开始且连续的 Vela/DF-ISE ID。
- 输出 `id_mapping.json` 并验证映射双射、无缺口、无悬空引用。
- 先构造 `vela.mesh_field_bundle.v1`，再由它派生 `mesh.json`、`doping.csv` 和各格式输出；
  禁止 exporter 各自调用 Gmsh 查询并形成不同排序。
- 将三角形统一成逆时针顺序，拒绝零面积、重复节点、非 Tri3 和非流形网格。
- 计算最小/最大角、面积、连通分量、区域覆盖和接触覆盖。
- 非钝角门槛由下游 solver profile 显式要求；Gmsh 不能保证时必须报告失败或由调用方
  明确选择兼容的 volume policy，禁止自动拆散 PN2D BV 的
  `element_edge_sg_gss_laux` 原子策略包。

### 验收

- Vela node/cell ID 均从 0 开始、连续且引用有效。
- `.msh` 原始 tag 经映射后与 `mesh.json` 坐标、连接和区域完全一致。
- 所有三角形 CCW，所有 Physical Surface 的三角形只属于一个器件区域。
- Line2 contact/interface elements、Tri3 和节点的原始 tag 与连续 ID 均存在双向映射；
  相同原始 tag 不得在不同 exporter 中映射成不同连续 ID。

## 阶段 6：在 Gmsh 节点上执行掺杂赋值

### 任务

- 将掺杂求值从张量网格 builder 中拆成与网格后端无关的模块。
- 按 Gmsh 原始节点坐标和 placement 源顺序计算 donor/acceptor。
- 实现并验证 `NoReplace`、`LocalReplace`、`Replace`；donor 与 acceptor 分别累计，
  最后生成 signed net doping。
- 实现已规划的 constant 与 Gaussian 子集，包括 PeakPos、PeakVal/Dose、结深参数、
  lateral Gauss/Erf、line baseline、Positive/Negative/Both、Eval/NoEval 和 target 限制。
- 二维 Positive 方向按 T-2022.03-SP2 oracle 定义为 baseline 起点到终点向量的左法向；
  这取代了实施前未经 oracle 验证的“右侧”假设。
- 保留 species 数据；可视化派生量不能回写或替代权威 donor/acceptor/net 数据。
- 对共享 region/material 界面节点先记录 incident membership set。区域定向 placement 在
  界面产生冲突值时，第一版必须 fail-closed 并输出冲突 placement/source location；只有
  Sentaurus oracle fixture 冻结了边界规则后才能实现确定性归并，禁止随意选择某个邻接
  region owner。

### 验收

- replacement、Gaussian 峰值、结深、方向、横向衰减和 NoEval 具有解析值测试。
- `doping.csv` 每行坐标和 node ID 可追溯到 Gmsh tag。
- 节点位于区域/窗口边界时遵循冻结的边界包含规则，不受浮点偶然性影响。
- semiconductor/insulator 界面上的掺杂字段策略有独立 golden；insulator-only 节点明确
  写 donor=acceptor=net=0，而不是缺失值或 NaN。

## 阶段 7：MSH NodeData、VTK/VTU 与通用 DF-ISE 导出

### 任务

- 先用 `gmsh.write()` 写基础 MSH，再用 `view.addModelData(..., "NodeData", ...)` 和
  `view.write(path, append=True)` 追加五个掺杂场；具体 API 调用必须由 4.15.2 round-trip
  测试验证，不能只检查文件中出现 `$NodeData` 字符串。
- 重新打开 MSH 并通过 Gmsh view API 读取字段，验证 tag 和数值。
- 新增与预处理流水线解耦的自有 VTK exporter：默认写 `.vtu`，同时支持 Legacy `.vtk`；
  可复用现有字段命名和 `VTKWriter` 测试经验，但不得遗漏 contact line cells。
- 将 minimal6 DF-ISE writer 抽取为通用 `scripts/sentaurus_dfise_writer.py`，支持任意
  合法二维 Tri3 网格、区域、材料、接触和节点数据。
- 每个 exporter 都消费同一份规范化 mesh/field bundle，不得各自重新计算几何或掺杂。
- 输出文件采用临时文件加原子替换；任一导出失败时状态保持 failed，不留下旧成功
  manifest。

### 验收

- `.msh` 在 Gmsh 中显示五个独立掺杂 view，并可 probe 节点值。
- `.vtu/.vtk` 中 PointData、CellData 和 contact line cells 完整且能 round-trip 解析。
- `.grd/.dat` 通过仓库 parser round-trip，区域、接触、拓扑、方向和数据集全部一致。
- 对每个节点，MSH、VTU、VTK、DAT 和 CSV 的 donor/acceptor/net 在序列化容差内一致。
- 每个文件先写同目录临时文件并完成回读，再原子替换正式文件；`run_manifest.json` 最后
  写入。任务开始时删除旧 accepted marker，任一 exporter 失败后不得残留旧成功状态。

## 阶段 8：CLI、manifest 和默认路径迁移

### 任务

- 扩展 `scripts/sentaurus_import.py device`，保持现有参数名：

```text
--sde-frontend offline-static             # 默认，现有兼容参数
--mesh-backend gmsh                       # 默认
--mesh-backend tensor-legacy              # 显式兼容
--mesh-backend official-tdr               # 仅内部 evaluated 模式
--mesh-quality-profile report-only        # 仅报告通用质量
--mesh-quality-profile vela-non-obtuse    # 求解合同可显式选择
--msh-format 4.1
--write-msh / --no-write-msh
--write-vtu / --no-write-vtu
--write-vtk                               # 可选兼容
--write-dfise / --no-write-dfise
```

- 默认生成 `device.msh`、`device.vtu`、`device.grd`、`device.dat`、Vela 中间文件和
  manifest；Legacy `.vtk` 可由配置或 CLI 开启。
- 帮助文本必须说明 Gmsh 为开源离线依赖，Sentaurus/tdx 不属于公开运行依赖。
- `simulation.json` 继续引用规范化的 `mesh.json` 和 `doping.csv`，避免求解器直接依赖
  Gmsh Python 运行时。
- manifest 状态区分 `generated`、`solved`、`accepted`、`failed`；只完成预处理时不得
  写成 solved/accepted。

### 验收

- 单条公开命令从 SDE/SDevice 脚本生成全部默认文件，并可选择继续运行 Vela。
- 默认执行路径不会探测或调用 `sde`、`tdx` 或其他商业二进制。
- `tensor-legacy` 产物明确带兼容/非权威标识。
- `sentaurus2022-evaluated + gmsh` 仍使用 Gmsh；只有显式
  `sentaurus2022-evaluated + official-tdr` 才导入官方 TDR。完整的 2×3 前端/后端组合矩阵
  都有 CLI 测试，其中非法组合必须 fail-closed。

## 阶段 9：公开 CI 与内部 Sentaurus oracle

### 公开 CI

- 测试分三层：不需要 Gmsh 的 IR/doping/exporter bundle 单元测试；固定 4.15.2 API 的
  OCC/mesh 集成测试；完整 `sentaurus_import.py device` E2E。缺少 API 时只有明确标记的
  API 测试可 skip，默认 gmsh CLI smoke test 必须失败而不是回退。
- CI 安装固定版本的开源 Gmsh Python API，运行 parser、OCC、mesh、doping 和 exporter
  回归测试，并保存 probe manifest。
- fixture 至少包含：
  1. 单区域矩形与两个接触；
  2. ABA/BAB 多区域几何；
  3. `insert-vertex + set-contact` 局部接触；
  4. 三种 replacement；
  5. line-baseline Gaussian 正/负/双向及 NoEval；
  6. region/window/interface refinement；
  7. 完整 PN2D SDE 预处理。
- 不对不同 Gmsh 版本要求 `.msh` 字节一致或完全相同的内部对角线；要求面积、区域、
  Physical Groups、接触、字段、质量和输出 round-trip 不变量一致。

### 内部许可证环境

- 使用 Sentaurus 2022 evaluated log 核对解析后的几何、区域、接触和掺杂定义。
- 执行 `tdx -d device.grd device.dat device.tdr`，再由现有 TDR importer 读回。
- 比较节点坐标、三角形拓扑、region/material、contact 和掺杂字段；差异生成机器可读
  report，不允许只保存截图或人工结论。
- 官方 SDE mesh 只作为语义/结果 oracle，不要求 Gmsh 生成完全相同的节点连接。
- TDX 命令参数先由目标 VM 的 `tdx -h` 冻结并写入内部 fixture 生成说明；在完成 smoke
  test 前，不把未经验证的命令拼写当作公开合同。

## 阶段 10：端到端 IV/BV 回归与文档

### 任务

- 使用 Gmsh 生成的完整 PN2D 网格重新运行 IV/BV 合同。
- 将网格差异与求解器差异分开报告，不得用已有张量网格或官方 TDR 结果替代 Gmsh
  E2E 验收。
- 更新 `docs/sentaurus_import.md`、安装说明、capability matrix、VM fixture 生成说明和
  输出格式文档。
- 保留旧后端迁移说明和弃用周期，不能无提示删除历史命令。

### 验收

- Gmsh PN2D 预处理产物满足所有跨格式一致性门槛。
- Vela 能读取生成的 `mesh.json/doping.csv` 并完成规定的 IV/BV 验收；若某个 BV
  目标仍未完成，报告最后接受偏压和明确原因，不伪装为完成。
- Release 构建、相关 Python regression、全量 `ctest --output-on-failure` 和
  `git diff --check` 通过。

## 总体验收标准

只有同时满足以下条件，才可宣布本计划完成：

1. 默认公开路径确实调用 Gmsh Python API，而非内部张量网格器或 Gmsh CLI 替代物。
2. SDE 几何、区域、材料和接触被映射为正确的 OCC 实体与 Physical Groups。
3. 输出仅包含一阶二维三角形，并通过方向、退化、区域和接触质量检查。
4. Gmsh tag 被保留，Vela/DF-ISE 编号从 0 开始且连续，映射可逆。
5. donor、acceptor、net doping 在所有输出格式中逐节点一致。
6. `device.msh` 含可由 Gmsh 直接显示的 NodeData 掺杂场。
7. `device.vtu` 和可选 `device.vtk` 含节点场、区域/材料 cell data 和接触 line cells。
8. 通用 `.grd/.dat` writer 不依赖 PN2D minimal6 固定拓扑，并通过离线 round-trip。
9. 内部 `tdx -d` 能读取 `.grd/.dat` 并生成可回读 TDR。
10. 缺少 Gmsh、遇到不支持 SDE 语义或导出失败时均 fail-closed。
11. 默认路径不依赖 Sentaurus 许可证、开发者绝对路径或未提交 fixture。
12. Gmsh 生成的 PN2D 网格实际参与 IV/BV 端到端验收。
13. `vela.mesh_field_bundle.v1` 是所有输出的唯一语义来源，并通过 schema 校验。
14. 目标分支已同步实施起点的 `main`，不存在用旧基线掩盖的测试结果。
15. Windows UCRT64、公开 CI 至少两个环境均通过 Gmsh initialize/mesh/finalize probe。

## 里程碑与优先级

### M0 / P0：依赖与合同可执行

- 完成阶段 -1、阶段 0 和 Device IR v2/schema/bundle schema 的最小骨架。
- 交付可重复的 Gmsh 4.15.2 probe、缺依赖 fail report 和 CI 安装路径。
- 此里程碑未通过前，不进入 OCC、refinement 或 exporter 实现。

### M1 / P0：最小真实跨格式闭环

- 状态：已完成（2026-08-25）。
- 两个不重叠 Rectangle region、两个 contact、constant doping、无 refinement function。
- Gmsh OCC + Tri3 + Physical Groups + NodeData。
- 从同一 bundle 输出 MSH/VTU/VTK/GRD/DAT/mesh.json/doping.csv 并全部回读比较。
- M1 是首个可演示版本，但不宣称支持 ABA/BAB、Gaussian 或完整 PN2D。

### M2 / P0：SDE 核心语义

- 状态：已完成（2026-08-25）。
- ABA/BAB fragment ownership、insert-vertex/contact、ordered replacement 和 Gaussian。
- region/material/window 与界面边界规则经过 Sentaurus 2022 oracle 冻结。

### M3 / P1：refinement 与完整 PN2D 预处理

- 状态：已完成（2026-08-25）。
- 完成阶段 4A；4B 未证明等价的功能继续 fail-closed。
- 完整 PN2D 生成、跨格式一致性和 mesh quality contract 通过。

### M4 / P1：IV/BV 与内部 TDX 收口

- 状态：已完成（2026-08-26）；公开 IV/BV 已实跑，内部 live TDX/TDR 严格往返通过。
- Gmsh PN2D 网格参与实际 IV/BV；内部环境完成 TDX/TDR round-trip。
- M0-M4 已全部验收，总体状态改为完成；PN2D BV 的端子守恒超限继续按 failed 报告，
  不通过放宽合同伪装为 accepted。

## 建议提交拆分

1. `build(gmsh): add Python API capability probe and pinned options`
2. `refactor(sde): freeze geometry contact refinement and doping IR`
3. `feat(mesh): add Gmsh OCC geometry and physical groups`
4. `feat(mesh): map SDE contacts and refinement fields to Gmsh`
5. `feat(mesh): normalize Gmsh Tri3 tags and qualify topology`
6. `feat(doping): evaluate ordered profiles on Gmsh nodes`
7. `feat(export): write MSH NodeData and VTK VTU fields`
8. `feat(export): add generic DF-ISE grid and dataset writer`
9. `feat(cli): make offline Gmsh the default SDE device path`
10. `test(e2e): add Gmsh cross-format and internal TDX gates`
11. `docs(sde): document offline Gmsh workflow and capability matrix`

每个提交必须有对应测试；不得把 IR 语义、网格算法、导出格式和求解器物理修改混在
同一个提交中。

## 主要风险与应对

- **Gmsh Python 模块与本机 `gmsh.exe` 版本不一致**：启动时比较版本并 fail-closed。
- **已确认的 MSYS2 wrapper/DLL 命名不一致**：使用 build 目录 runtime staging 或经验证
  wheel；不得复制 DLL 进源码树、修改系统安装或只靠设置 `PYTHONPATH` 宣称解决。
- **OCC fragment tag 不稳定**：按几何分类和源 IR 重建 ownership，不依赖返回顺序。
- **各向异性 refinement 无直接等价映射**：未验证前 fail-closed，不静默近似。
- **非结构网格出现钝角**：提供明确质量 gate；不得自动改变 PN2D BV 原子 solver
  profile 规避失败。
- **不同 Gmsh 版本连接不同**：固定版本并基于不变量测试，不维护脆弱的全文件 golden。
- **MSH、VTK、DF-ISE 字段错位**：所有 exporter 消费同一规范化 bundle，并实施逐节点
  跨格式比较。
- **TDX 属于商业环境**：公开 CI 只做 DF-ISE parser round-trip，TDX 只在内部 gate 使用。

## 首轮实施终点

首轮不以完整 IV/BV 为开始条件。第一里程碑必须先交付一个最小但真实的闭环：读取含
两个区域和两个接触的 SDE fixture，通过 Gmsh Python API 生成 Physical Groups 和
Tri3 网格，在节点上写入 donor/acceptor/net doping，并产出能重新读取且相互一致的
`.msh`、`.vtu`、`.vtk`、`.grd/.dat`、`mesh.json` 和 `doping.csv`。该闭环通过后，
再扩展完整 refinement/Gaussian 和 PN2D IV/BV。首轮 fixture 不使用重叠 body 或界面
冲突 doping，从而验证真实 API 与 exporter 闭环，但不能被误报为 M2 或总体完成。
