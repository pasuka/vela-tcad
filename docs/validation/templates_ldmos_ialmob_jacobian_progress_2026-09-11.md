# IALMob 单元 Jacobian 接口与暂停交接 — 2026-09-11

## 当前完成边界

本轮实现独立距离/晶向构造、六变量低场偏导、单元顶点迁移率的九电势导数，
并将这些 C++ 接口用于固定原生状态全场对照。**尚未把新接口接入生产
CoupledDDAssembler 的全局输运残差/Jacobian、ContactCurrent 或恢复求解路径；
尚未启动 D4 自洽曲线，不能认领 D4 通过。** 用户要求在北京时间 16:35
保存并暂停，当前阶段的交接位置是生产全局接入之前。

原 D5/G3 资格保持其已记录的程序和配置范围；本轮未重跑它们，未改原门限。
没有把 IALMob 别名映射到 Masetti 或 Lombardi，也没有开启新的生产模型选项。

## 已实现及核对

- `IalMobility::evaluateWithDerivatives`：Nd、Na、n、p、法向场绝对值、距离
  六个 SI 输入的正向自动微分。标量值与导数共用公式，避免两套实现漂移。
  无散射的无穷分量不参与无穷乘零的导数运算。零场/零密度分数幂拐点采用
  明确的零分支斜率，不将其声称为经典导数；正输入使用两档步长的中心差分。
- `evaluateIalElementMobility`：每个顶点分别保留 Nd/Na、n/p、距离和晶向模型；
  使用未归一化 P1 距离梯度。密度响应由调用方按实际 Fermi/Boltzmann 统计提供。
  HFS 在顶点施加后进行 Measure 加权，包含 RefDens 标量插值、界面切向
  QF 投影和有效电极电场回退。返回各节点 psi/phin/phip 共九列的导数。
- `buildIalInterfaceGeometry`：距离仍取到真实线段的最短距离；晶向取最近界面
  顶点的等权入射单位法向，经显式晶体坐标变换选择晶面族。该二维规则在当前
  网格全部 5723 个 Si 节点上与原生标签一致，不输入原生标签。相同网格上，
  最近线段法向有 15 处误配；插值顶点法向有 7 处。两种最近顶点候选都吻合，
  本次接口明确采用全局最近界面顶点；不能据此认领任意三维网格通用规则。
- 增加 `ialmob_geometry_probe`、`ialmob_element_probe`。二者为固定状态诊断，
  不是完整曲线运行入口。原 `ialmob_probe` 增加可选 `derivatives: true`，偏导
  顺序为 Nd/Na/n/p/Enormal/distance，保持原默认输出字段兼容。

## 实测证据

环境为 Windows UCRT64、C++20、Release，`-O3 -DNDEBUG`，gprof 关闭。
本次配置检测到 HDF5、SuiteSparse SPQR/UMFPACK；默认排序为 Eigen SparseLU
COLAMD。本轮局部接口和固定状态对照不调用稀疏线性求解器，不能据此比较 Newton 耗时。
几何、状态及导数接口均为 SI；原生浓度、迁移率、坐标转换分别使用
cm^-3 → m^-3、cm²/(V s) → m²/(V s)、μm → m。

网格为 10241 节点、19782 三角形，其中 10515 个 Si 单元；对照固定于
300 K、Vg=8 V/Vd=40 V 的 D4 保存状态。没有新增 Sentaurus Newton 作业。

| 检查 | 结果 |
| --- | --- |
| 四组原局部输入重放 | 63090 个载流子样本；迁移率相对漂移最大 1.23e-15 |
| 六变量偏导有限性 | 378540 个偏导全部有限 |
| 原生样本抽样中心差分 | 594 项；最大归一化误差 9.77e-9 |
| 独立 C++ 晶向 | 5723 节点，0 处不匹配；{100}/{110} = 4968/755 |
| 全部 Si 单元低场最大相对误差 | 电子 3.93e-6；空穴 5.07e-6 |
| 全部 Si 单元高场最大相对误差 | 电子 1.023e-5；空穴 3.97e-6 |
| 聚焦 Release 测试 | 12 个 IALMob Catch2 用例及 2 个工具协议回归，共 14/14 通过 |
| 源码及差异检查 | ascii_sources 1/1 通过，git diff --check 无内容错误 |

误差表为无量纲相对值，不是百分数。九列 Jacobian 单元测试覆盖两种统计、
HFS 开关、PartialLayer 开关、接触回退开关、混合晶向模型、两档扰动步长和
规范平移不变性。原生全场重放用于验证函数值；其中密度响应显式置零，
不能把该重放解释为原生全局 Jacobian 比较，导数资格来自独立活密度差分测试。

全场重放已经改用由电势构造的 C++ P1 场和独立 C++ 晶向/距离。
**Measure 仍来自已有原生几何证据**，沿用已验证的局部槽映射；尚未解决
生产外部 AverageBox 边耦合与单元顶点迁移率之间的完整装配合同。
有效电极为本例实际 source/drain/substrate，接触边集合另外包含热接触边；
生产接入时须按接触类型构造，不能依据 `th_lat` 名称在库代码中硬编码。

## 保存位置和复现

证据目录：`reference_staging/templates_ldmos_followup_20260911/ialmob_jacobian_r3/`。

- `replay/summary.json`：原值漂移、偏导有限性、抽样差分及程序哈希。
- `orientation_summary.json`：几何候选对照；原生标签仅用于评分。
- `element/summary.json`：独立 C++ 几何加单元接口的全场对照。
- `checkpoint/`：暂停保存的本轮涉及源码/文档副本及清单；包含工作树当前内容，
  不替代 Git 历史，也不表示此前其它未提交修改已提交。
- 父目录 `ialmob_r3_final_tests.log`、`ialmob_element_build.log`、
  `ialmob_r3_configure.log` 保存本轮构建/验证日志。
- 父目录 `audit_ialmob_orientation_r3.py`、`replay_ialmob_jacobian_r3.py`、
  `replay_ialmob_element_r3.py` 为实验复现脚本，依赖上轮外部证据；重放脚本拒绝
  覆盖已有输出目录。再次运行前使用新输出目录，保留本次证据。

从当前 worktree 根复查代码：

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build --preset windows-ucrt64-release --target test_mobility ialmob_probe ialmob_geometry_probe ialmob_element_probe -j 4
ctest --preset windows-ucrt64-release -R '^(IALMob|ialmob_.*protocol)' --output-on-failure
```

## 恢复后的执行顺序

1. 明确单元 Measure 与外部 AverageBox 输运边系数的装配合同。现有
   `AssemblerUtils::edgeMobility` 向模型传入的 n/p 固定为零；
   `CoupledDDAssembler::cachedEdgeMobility` 缓存低场值，均不适用于活 IALMob。
   新增显式配置、参数组、几何缓存和随状态刷新的单元迁移率/导数缓存；
   在同一残差状态中复用，不从原生参考载流子、标签或迁移率取值。
2. 同步接入全局残差、固定稀疏结构及九列链式导数、端电流、行闭合和密度恢复。
   以完整装配的定向差分、守恒及单位缩放等价性检查接入正确性。
   局部 Jacobian 通过不能替代这一步。
3. 明确生产有效电极/边界分类、单元 Measure 输入依赖及零场分支处理，验证
   接口的非生产选项不会静默生效；对 D5/G3 做必要回归。
4. 使用新冻结 Release 获取合格 D4 零压种子，先 Vg=8 V 从零压的前8个精确点，
   再两栅压0–40 V完整曲线及电流比，逐项沿用 D4 原门限。不得把 D5 种子
   的原资格直接当作 D4 零压已闭合，也不得在全局接口缺失时启动替代模型曲线。

本次没有 Git 提交或推送；暂停后需要用户明确要求继续才恢复工作。
保存时已检查本轮编译、诊断和测试会话全部退出，相关活动进程数为零。
