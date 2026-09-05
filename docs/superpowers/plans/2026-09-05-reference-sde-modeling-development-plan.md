# reference_tcad SDE 读取、建模与 GRD/DAT/MSH 导出开发计划

日期：2026-09-05。状态：基于合并树的实测计划，新增建模能力尚未实施。

## 目标与本次基线

让 `reference_tcad` 中的二维 SDE 器件源脚本，经过离线读取、语义求值、几何建模、
网格生成和节点掺杂计算，直接输出 `device.grd`、`device.dat`、`device.msh`。
正常建模流程不需要 SDevice 输入、商业 SDE、TDX 或预先存在的 TDR。
沿用现有 `offline-static → Device IR v2 → Gmsh → mesh-field bundle → exporters`。
C++20 求解核心保持原架构；本计划主要扩展已有 Python 预处理层。

- 工作目录：`D:/code-repo/copilot-worktrees/vela-tcad-dev`。
- 分支：`codex/sentaurus2022-sde-workflow`。
- 交接提交：`abc2d4a`；合入本地 main：`8e99e57`。
- 本次普通双父合并：`505e1ea`，无文本冲突；main 的 22 个新提交已纳入。
- 原 M0–M5 是已实现基础，历史 PN2D IV/BV 验收不能代替本次合并树回归。
- 本次用户要求已经确定发布策略：TDR 和原始日志仅用于本地验证，不纳入当前提交树。

依据：`docs/validation/sentaurus_gmsh_handoff_2026-09-05.md`、原 M0–M5 计划及
M5 验收记录；当前源码和本次逐脚本离线建模结果优先于旧能力矩阵。

## 实测覆盖范围

使用 Git 跟踪清单查找含几何创建的 `.cmd/.scm`，而非只匹配 `*sde*.cmd`，因此包括
TransportModels 的 `nMOS_dvs.cmd`。使用未删减、未改写的原脚本逐份执行默认 M5 管线。

| 器件源脚本（相对 reference_tcad） | 本次实际结果 | 后续工作 |
| --- | --- | --- |
| `pn2d_sentaurus2018/source/pn2d_sde.cmd` | 导出及跨格式回读通过，transfinite | 纳入批量稳定回归 |
| `pn2d_sentaurus2018_coarse7x3/source/pn2d_sde.cmd` | 导出及跨格式回读通过，transfinite | 验证脚本语义；目录名称不代表生成网格必须为 7×3 |
| `nmos2d_sentaurus2018/source/nmos2d_sde.cmd` | 导出及跨格式回读通过，非结构网格 | 多材料、共享界面和四接触基线 |
| `schottky_charon_sentaurus2018/source/schottky_n_sde.cmd` | 导出及跨格式回读通过，transfinite | Schottky 结构基线 |
| `genius_bjt_sentaurus2022/source/bjt_sde.cmd` | 第 6 行：define 内创建几何不支持 | body 返回值绑定、General 解析掺杂，然后 MaxTransDiff |
| `genius_bjt_sentaurus2022/source/bjt_sde_local_refined.cmd` | 第 7 行：同上 | 同上，另验证局部加密效应 |
| `pn2d_sentaurus2018_skewed_tri3/source/pn2d_sde.cmd` | 第 27 行：create-polygon 不支持 | Polygon body 与共形界面 |
| `singledevice_sentaurus2018/source/singledevice_sde.cmd` | 第 57 行：函数式 define 不支持 | 有界函数、body 绑定、圆角、MaxTransDiff、镜像 |
| `transportmodels_sentaurus2022/source/nMOS_dvs.cmd` | 第 57 行：函数式 define 不支持 | 辅助函数、圆角、refine-box/doping/interface 快捷命令和镜像 |
| `bvmethods_sentaurus2018/source/bvmethods_nmos_sde.cmd` | 第 32 行：PolySi 材料名不支持 | 别名、圆角、body contact/remove、MaxTransDiff、multibox、镜像 |

器件源脚本合计 10 份：4 份导出成功、6 份失败。失败行只代表首个实际阻塞；后续缺口
根据源脚本及实现审阅列出，尚未冒充逐项动态验证。

额外检查的 `pn2d_sentaurus2022_oracle/input/pn2d_sde.cmd` 也通过，但它是已有 PN2D
源的验证副本，不计入新增器件覆盖数。`vela_sde_evaluated_wrapper.scm` 的 `set!`
在离线模式明确失败；它属于内部日志采集包装器，不应计入公共器件失败率。
`pn2d_sentaurus2018_minimal6` 只有固定网格/SDevice 用例，不凭空宣称有待导入 SDE 源。

机器可读审计摘要：`docs/validation/reference_sde_audit_2026-09-05.json`。
完整本地生成结果：`build-release/reference-sde-audit-20260905/`。

## 输出合同

三个主文件从同一规范化 bundle 导出，不能分别重建几何或重新采样掺杂。

| 输出 | 必须保留的信息 |
| --- | --- |
| GRD（DF-ISE ASCII） | 坐标、唯一边、带方向的单元引用、区域/材料归属、接触边界 |
| DAT（与 GRD 成对） | 节点 donor/acceptor、带符号净掺杂、可追踪的 active species、有效区域 |
| MSH（现有 4.1 ASCII） | 三角形及边界线、Physical Groups/Names、节点掺杂 NodeData |
| 辅助 JSON/CSV | IR、bundle、ID 映射、单位、输入/输出哈希、网格质量与失败原因 |

保持长度 um、掺杂 cm^-3 的既有约定；统一保留 Gmsh tag、Vela ID、DF-ISE ID 映射。
网格必须有正面积、正确方向、无裂缝和重复体单元；接触不得丢失或跨越错误边段。
保留共享节点的 region-local 掺杂语义及 placement 顺序，不用简单平均替代。
对每个新几何操作，region-local 数据、边界 ownership 和字段传播必须一起实现。
所有成功结果都必须经过现有独立回读路径；失败返回非零退出码和带源位置的报告。
失败重跑不得留下可误认成新成功结果的旧主文件，此项加入批量工具的失败恢复验收。

这里的 DAT 是建模阶段掺杂场；电势、电流、载流子浓度等求解场属于后续仿真输出，
不作为本计划的隐含交付要求。三种格式成功也不自动等于 Sentaurus 物理曲线一致。

## 实施顺序与里程碑

### P0：合并基线和验证文件隔离（本次收尾）

1. 合入本地 main，检查自动合并的 solver、CMake、CTest 及模板交叉点。
2. TDR、`.log`、`.log.eval`、stdout/stderr、版本和退出码原始文件保留本地并退出跟踪；
   保留小型源文件、IR golden、原始哈希及结构化验证摘要。
3. 公共 oracle 测试不依赖私有原件；有原件则校验哈希，内部严格模式缺文件必须失败。
4. 恢复 oracle 源文件文本 diff；仅为 YAML 规范 LF，不批量改写哈希封印原件。
5. Windows Release 构建和回归；Linux/GitHub 的新合并树证据单独列出状态。

验收：本地 main 是祖先；当前 Git 索引没有 TDR/原始日志；公共无原件环境能够验证
源脚本与 frozen IR；内部模式不会把缺证据当成成功；`git diff --check` 通过。

### P1：把已通过的四份源脚本变成公共建模基线（约 2–3 工作日）

- 复用 `scripts/sentaurus_gmsh_pipeline.py` 作为仅需 SDE 的公开入口，并补充使用文档。
- 新增版本化 case manifest：源脚本、源哈希、类别、期望区域/材料/接触、验收门槛。
- 新增批量 runner 和参数化 regression，输出 parse/model/export/reference 四级状态。
- 区分器件源、验证副本、内部包装器、只有固定网格的用例，防止覆盖统计失真。
- 新增仓库卫生门禁：检查 tracked 文件，禁止 TDR、原始日志和采集目录混入提交。

验收：四个真实脚本从干净输出目录一条命令生成三个文件；跨格式拓扑、标签和节点场
比较通过；同输入重复运行满足确定性合同；失败报告能定位脚本、行、命令和缺失能力。
现有入口示例：

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python scripts/sentaurus_gmsh_pipeline.py `
  --sde reference_tcad/nmos2d_sentaurus2018/source/nmos2d_sde.cmd `
  --output-dir build-release/reference-modeling/nmos2d `
  --runtime-dir build-release/gmsh-python-runtime
```

### P2：有界 Scheme 语义与 BJT 完整建模（约 9–13 工作日）

- 扩展 `sentaurus_sexp.py` / `sentaurus_device_ir.py`：有类型的 body handle、
  `(define body (sdegeo:create-rectangle ...))`、源脚本实际使用的纯辅助函数和表达式。
- 维护语句执行顺序、作用域与 source location；限制递归/执行步数；任意 shell 调用
  继续拒绝。只实现实际所需子集，不承诺完整 Elk Scheme。
- 在 `sentaurus_device_ir.py` / `sentaurus_doping.py` 新增 BJT 实际使用的
  `sdedr:define-analytical-profile ... "General"`：解析参数赋值串及表达式中的
  x/y、abs、exp、乘除和幂运算。采用受限表达式 AST，不使用 Python eval；保存设备坐标
  求值语义、默认值及源位置，并用源合同的解析点和截面对比。不能用现有 line-baseline
  Gaussian 替代 General 表达式。同步扩展 IR/schema 与能力清单；旧 v1/v2 输入显式迁移，
  不把未知表达式类型降级为 constant。
- 在 `sentaurus_gmsh_backend.py` 实现 MaxTransDiff：先用小型 oracle 冻结数据变换、
  阈值、掺杂近零/符号变化及 min-size 停止语义，再选择场采样和迭代重网格策略。
- 把每轮指标、停止原因、最小尺寸、最大迭代数和节点数上限写入报告；未达门槛不能
  静默以固定尺寸网格当作等价结果。MaxGradient 作为独立扩展，非这两份 BJT 的前置。
- 扩展并复用 placement 求值，按 General 表达式评价连续掺杂，不从参考 TDR 拷贝场完成建模。

验收：BJT 两份原脚本不改写即可输出三种格式；局部加密确实作用于指定窗口；几何、
三个接触、掺杂截面/积分、结位置和网格尺寸合同通过。累计 6/10 器件源脚本通过。
新增用例覆盖 body 返回值、函数参数、错误作用域、超限以及自适应停止失败。

### P3：Polygon body 与偏斜网格用例（约 3–5 工作日）

- 扩展 Device IR/schema 和 Gmsh OCC 建模，支持源脚本实际使用的二维简单多边形。
- 验证轮廓闭合、自交、退化、方向和容差；布尔分割后保留 region/material ownership。
- 共形共享边、边段接触选择、节点掺杂与三种 exporter 全链路覆盖；明确区分已有
  Polygon refinement window 和新增 Polygon body。

验收：skewed_tri3 原脚本通过，区域面积和邻接关系正确，完整接触边界可回读。
不要求 Gmsh 复刻固定参考网格的节点编号或三角形连接；按几何和场比较。
一般建模允许合格的钝角网格；PN2D BV 的非钝角资格门槛仍单独强制执行。
累计 7/10 器件源脚本通过。

### P4：SingleDevice / TransportModels 的完整 MOS（约 9–13 工作日）

- 实现 `fillet-2d` 的二维曲边几何、圆角半径校验和边界 chord-error 网格合同。
  不将圆角无声忽略为直角，接触选择与圆角后的拓扑保持一致。
- 复用 P2 的辅助函数/body 绑定/MaxTransDiff；验证细氧化层和材料界面加密。
- TransportModels 还实际使用 `sdedr:refine-box`、`sdedr:refine-doping`、
  `sdedr:refine-interface`，当前 HANDLERS 未注册。先冻结其参数、缺省值、作用区域和
  叠加次序，再规范化到 IR 的 window/placement/field 合同，不能因为它与 SingleDevice
  外形相似就漏掉这三种快捷命令。
- 为 `@node@` 等实际使用占位符提供显式参数和缺参数诊断。
- 对脚本末尾已知 `tdx -mtt ... -ren drain=source ...` 形式做严格语义识别，转换为
  bundle 中的镜像操作，在本地执行坐标反射、轴上节点合并、三角形方向修正、
  region/contact 重命名和所有场传播。不得执行原始 `system:command`。
- 先按半器件验收几何，再验收完整对称器件，防止输出只有 drain 的半器件却报完成。

验收：SingleDevice 和 TransportModels 两份原脚本的完整器件均生成 GRD/DAT/MSH；
源漏接触、对称轴无裂缝、区域命名、圆角和掺杂镜像正确。累计 9/10 通过。
未知 TDX 参数或任意其他外部命令明确失败，不转去依赖商业 TDX 完成正常建模。

### P5：BVmethods 与最终参考覆盖（约 5–8 工作日）

- 在材料映射中显式处理 `PolySi`/`PolySilicon` 别名并保留原名称；不泛化为未知材料
  一律映射到 Si。
- 实现 body contact + `remove`：明确创建态接触边界、移除的 body、保留邻接区域、
  后续引用已移除区域的掺杂/网格命令语义；以小型 oracle 固定行为。
- 实现 multibox 尺寸、正负 grading ratio、方向、窗口裁剪及与已有 refinement 的组合。
  检查 BVmethods 的 Erf 横向 Gaussian 与 Symm 方向现有实现，不重复宣称它们尚缺。
- 使用 P4 本地镜像操作生成完整 BVmethods 器件；保证 gate 删除后接触仍完整有效。

验收：10/10 器件源脚本离线生成三个主文件，逐用例满足几何、接触、掺杂和网格合同。
如果某项官方语义尚未确认，应保留明确阻塞，不删命令或改写简化版原脚本来报覆盖率。

### P6：资格验证、文档与集成（约 3–5 工作日）

- 对所有 case 运行独立 GRD/DAT 和 MSH 回读；在公共 CI 中运行完整离线 corpus。
- 内部 TDX 将自产 GRD/DAT 转成临时 TDR，验证格式拓扑和节点场；同时对官方 SDE
  几何、接触、掺杂及尺寸约束做参考对比。两种网格不要求节点数/节点编号完全相同。
- 几何面积、边界长度、结位置、掺杂截面/积分、尺寸分位数分别设 case 合同。
  门槛在实施前基于解析几何、小型 oracle 和网格收敛确定，不能为单次结果放宽。
- Windows UCRT64 Release 和 Ubuntu 使用独立 build 目录完成全量回归，获取 GitHub
  Actions 实际记录后再标记 CI 通过；更新能力矩阵与用户使用文档。
- 保留 PN2D IV/BV、BJT 现有 solver 回归，防止预处理扩展破坏既有仿真；不把所有 MOS
  的 IV/BV 物理精度扩展成当前建模项目的必交付项。

单人串行估算约 31–47 工作日，不含商业参考环境排队；P2 自适应等价性、P4 圆角与镜像、
P5 删除/渐变语义是主要不确定性。建议每个里程碑独立提交，先交付现有四例的稳定入口。

## 仓库与验证证据政策

- 公共仓库：源码、SDE 原脚本、schema、配置、手工或解析构造的小型 fixture、
  脱离原始输出的结构化验收摘要及哈希。
- 忽略目录：生成 GRD/DAT/MSH、TDR、日志、可执行文件、运行时 DLL、大规模仿真结果。
  GRD/DAT/MSH 是产品输出，并不意味着所有生成结果都应提交；回归优先现场生成回读。
- 内部证据：TDR、`.log*` 和 stdout/stderr/版本输出仅本地验证；不得以改扩展名或
  把完整日志嵌入 JSON 的方式重新提交。小型原始 mesh command 可作为可审阅文本保留。
- 已从索引移除的历史 TDR/日志仍存在于旧 Git 提交。此次不重写分支历史；若需要
  对外发布且要求历史对象也不包含这些文件，应另行安排历史清理及协作迁移。
- 本次没有推送分支，也没有把开发分支合回 main。

## 新合并树验证记录

逐脚本实测已完成，见上表及 JSON。SDE/Gmsh 定向 CTest 9/9 通过（206.22 秒）；
本地 oracle 3/3 通过。无原始证据的公共 checkout 模拟：2 项通过、1 项明确跳过；
内部严格模式确认缺证据失败。

- Windows UCRT64 Release 配置、136/136 编译与链接步骤全部完成。
- 另 4 项定向 CTest 全部通过：ASCII、Gmsh runtime、Sentaurus device E2E contract、
  sample integration（363.97 秒）。两组共 13/13 项通过。
- 修改 JSON 的严格解析、暂存差异检查及相对 main 的完整差异检查均通过。
- 当前索引无 TDR/原始日志；原始文件哈希校验通过；YAML 仅改变行尾。
- 本次未运行全量 Windows CTest、WSL/Ubuntu、GitHub Actions 或 live TDX/SDE；
  也未重新生成完整 PN2D IV/BV 物理曲线。这些新合并树证据仍按各自合同待验收。
历史文档的 712/712 等数字属于旧提交，不能作为新合并树的测试结果。

## 下一项具体开发任务

实施 P1：为四份已通过的源脚本加入版本化 case manifest、批量回归和单命令使用文档，
同时将本次其余六份脚本的失败类别固定成可定位的基线。随后进入 P2，先修复 BJT 的
body 变量绑定和 General 解析掺杂，再实现并验证 MaxTransDiff。无需重写已有 GRD/DAT/MSH exporter。
