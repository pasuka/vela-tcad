# Si/SiO₂ 介电系数：Sentaurus、Charon 与 Genius 源码核对

2026-09-07。结论：三者在本次核对的普通、各向同性、无界面偶极层的材料贴合网格上，都保留各材料自己的介电系数，再装配其几何/通量贡献。它们的网格几何和装配组织不同，不能据此宣称所有边系数数值相同。尤其 Genius 中的端点算术平均不是跨 Si/SiO₂ 的材料平均。

本轮是手册及源码阅读，没有修改求解器或运行 Charon/Genius 仿真。读取版本：Charon `7cc38745625a6011ae3584ed111ec7ee74fb890e`，Genius `543da8452d5dfd33e6f8c457f962f6f670f0fce7`；两库 tracked 工作区均无改动。Sentaurus 依据本地 T-2022.03 用户手册及既有原生几何导出。

## 1. 界面条件与要区分的几何

普通界面满足电势连续；以法向从 Si 指向氧化层，`n·(D_ox−D_Si)=σ_interface`，其中 `D=εE`。无界面自由电荷时，法向 D 连续，法向 E 可因 ε 不同而跳变。界面电荷作为源处理，不通过修改一个平均 ε 模拟。

对本算例沿界面的共边 i–j，Sentaurus box 系数可写为：

```
G_ij = Σ_e ε_e κ_ij^e
κ_ij^e = d_ij^e / l_ij                 （2D）
F_ij = G_ij (ψ_i − ψ_j)               （取此通量方向约定）
```

若只有一个 Si 和一个氧化层相邻单元，则 `G=ε_Si g_Si+ε_ox g_ox`。在总几何系数非零、权重适合平均的情况下，也可写成几何加权有效值 `ε_eff=(ε_Si g_Si+ε_ox g_ox)/(g_Si+g_ox)`。对于可能带符号或被原生算法修正的 box 几何，直接累加 ε×几何贡献比笼统称为“平均”更准确。

这与法向穿过两层材料的一维串联路径不同：后者等效 `ε_eff=(l_Si+l_ox)/(l_Si/ε_Si+l_ox/ε_ox)`。不能把串联路径的调和平均直接替换为本例共边的装配规则。

## 2. Sentaurus

T-2022.03 手册 Chapter 38、pp.1177–1178、式 1264–1266：先给出 box 离散，再明确方程实际按单元累加；Table 188 中 Poisson 的边物理项为 `ε(u_i−u_j)`。因此结合单元 box 系数，得到 `Σ_e ε_e κ_ij^e`。

来源：[本地手册：elementwise assembly](D:/code-repo/vela-tcad/.worktrees/simplemos-sdevice-validation/build-release/m79_research/sdevice_ug_local_2022.txt:62247)。介电参数定义在 Chapter 7、pp.229–230 的 `Epsilon/epsilon`；特殊偶极层模型可允许势跳跃，不能套用普通界面连续条件到这些扩展模型。

本轮前一阶段的原生系数导出与该规则相符：两个网格各 20 条 Si/SiO₂ 共界面边全部通过公式核对，Vela 旧算法的系数偏高约 42.7%。参见[逐边数据](../../reference_tcad/simplemos_sentaurus2022/remaining_poisson_20260907/interface_formula.csv)与[定位验证](simplemos_remaining_poisson_dielectric_localization_2026-09-07.md)。这是本算例的实测证据，不是获取或声称审阅了 Sentaurus 私有源码。

## 3. Charon

Charon 官方手册 §1.6 说明其支持多种离散方法；Poisson/Laplace 包括 FEM 与 CVFEM-SG。[官方用户手册](https://www.sandia.gov/app/uploads/sites/106/2022/08/Charon_UserManual-v2_2.pdf#page=21)

材料参数绑定在 material model/physics block 上。`ClosureModelFactory` 对当前材料创建常数场，并分别提供 IP 和 BASIS 布局；Si 单元和 SiO₂ 单元因此可以在同一几何节点位置保留不同的 ε 数据。

- [材料常数创建](D:/code-repo/tcad-charon/src/evaluators/Charon_ClosureModel_Factory_impl.hpp:667)：显式参数创建当前材料常数；1164 行附近处理材料默认值；2006 行附近创建 Cell/IP 与 Cell/BASIS 布局。
- [FEM 通量](D:/code-repo/tcad-charon/src/evaluators/Charon_PotentialFlux_impl.hpp:63)：`phi_flux(cell,ip,dim)=Lambda2*rel_perm(cell,ip)*gradphi`。随后通过 [FEM Poisson 装配](D:/code-repo/tcad-charon/src/equation_sets/Charon_EquationSet_NLPoisson_impl.hpp:115) 积分 `ε∇ψ·∇N`。物理单位下其刚度为 `A_ij=Σ_e ∫_e ε_e ∇N_i·∇N_j dΩ`。
- [CVFEM 通量](D:/code-repo/tcad-charon/src/evaluators/Charon_SGCVFEM_PotentialFlux_impl.hpp:105)：每个单元内插值 ε、求形函数梯度对应的 ∇ψ，再形成通量；[子控制面通量积分](D:/code-repo/tcad-charon/src/evaluators/Charon_Integrator_SubCVFluxDotNorm_impl.hpp:100) 用加权法向积分，并向单元边两端装配正负残差。
- [旧 SGCharon1 通量](D:/code-repo/tcad-charon/src/evaluators/Charon_SGCharon1_PotentialFlux_impl.hpp:151) 虽有 `(rel_perm(cell,node0)+rel_perm(cell,node1))/2`，但两值仍来自同一 cell，不能读成将两种材料的 ε 混合。

默认连续电势字段在相接单元块间汇合，两侧通量进入共享电势方程。代码允许额外的 discontinuous-field 配置，本结论限定普通 Si/SiO₂ 配置。现成 [nMOS 输入](D:/code-repo/tcad-charon/test/nightlyTests/trapSRH/nmosfet.nlp.inp:14) 明确分为 silicon 的 NLP-CVFEM 与 sio2 的 Laplace-CVFEM 两个 block，两侧使用 `ELECTRIC_POTENTIAL`。

Charon 的 CVFEM 使用单元梯度、子控制面加权法向；其一般形式不能直接视为 Sentaurus 的逐边二点 box 系数。两者一致的是材料贡献在局部计算后累加的原则，而非已经验证其全部几何离散等价。

## 4. Genius-TCAD-Open

Genius 为同一几何界面位置保留不同区域的 FVM_Node。每个区域的节点 ε 由该区域自己的材料模型初始化：[Si 区域](D:/code-repo/Genius-TCAD-Open/src/solution/semiconductor_region.cc:259)、[绝缘区域](D:/code-repo/Genius-TCAD-Open/src/solution/insulator_region.cc:169)。FVM 邻接表限定当前区域，参见 [FVM_Node](D:/code-repo/Genius-TCAD-Open/include/solution/fvm_node_info.h:386)。

区域内 Poisson 边通量为：

```
eps_edge = 0.5*(eps1+eps2)
flux = eps_edge * cv_surface_area / distance * (V2−V1)
```

源码：[半导体 Poisson](D:/code-repo/Genius-TCAD-Open/src/solver/poisson/poisson_semiconductor.cc:149)、[绝缘体 Poisson](D:/code-repo/Genius-TCAD-Open/src/solver/poisson/poisson_insulator.cc:118)。`eps1/eps2` 是当前区域中边的两个端点值；均匀 Si 区域两者都是 ε_Si，均匀 SiO₂ 区域两者都是 ε_ox。因此这处平均没有产生 `(ε_Si+ε_ox)/2`。

界面处理分两步：

1. [预处理](D:/code-repo/Genius-TCAD-Open/src/solver/poisson/poisson_boundary_is_interface.cc:72) 把绝缘侧 Poisson 行加入半导体对应行，并清空原绝缘行。全局求解器实际执行 [VecAddClearRow](D:/code-repo/Genius-TCAD-Open/src/solver/poisson/poisson.cc:342)。
2. 绝缘侧被替换为 [V_ins−V_semi=0](D:/code-repo/Genius-TCAD-Open/src/solver/poisson/poisson_boundary_is_interface.cc:196)，使电势连续；半导体行另加界面固定电荷源。Jacobian 同样合并两侧行后添加连续性约束。

消去重复电势后，均匀材料区域对同一共边的贡献为 `ε_Si A_Si/l + ε_ox A_ox/l`，A 为各自区域内控制面贡献。Genius 的实际 A 仍受自身网格/控制体处理影响，未证明与 Sentaurus 原生 box 系数逐项相同。

还核对了漂移扩散生产路径：[DDM1 半导体通量](D:/code-repo/Genius-TCAD-Open/src/solver/ddm1/ddm1_semiconductor.cc:228)、[DDM1 绝缘侧通量](D:/code-repo/Genius-TCAD-Open/src/solver/ddm1/ddm1_insulator.cc:117)、[DDM1 界面行合并](D:/code-repo/Genius-TCAD-Open/src/solver/ddm1/ddm1_boundary_is_interface.cc:76) 采用相同区域划分原则，因此结论不只来自初始化 Poisson 求解器。

## 5. 对 Vela 的含义与参数差异

Vela 旧算法相当于 `((ε_Si+ε_ox)/2)*(g_Si+g_ox)`，与逐材料贡献相加一般不等。误差正是 `0.5*(ε_Si−ε_ox)*(g_ox−g_Si)`；只有介电常数相同或几何权重相同等特殊情况才消失。

后续实现应保留材料对应的局部介电几何贡献，并确保残差和 Jacobian 使用同一系数；仍需分别验证网格几何、非 Delaunay/钝角处理及其他 Poisson 诊断分支。不能只复制 Genius 中那一行端点平均，也不能只将平均方式改成调和平均。

材料数值还需单独统一：本地 Charon 的 [Si 默认 εr=11.8](D:/code-repo/tcad-charon/src/Charon_Material_Properties.cpp:1061)、[SiO₂ 默认 εr=3.9](D:/code-repo/tcad-charon/src/Charon_Material_Properties.cpp:2211)；Genius 的 [Si 默认 εr=11.7](D:/code-repo/Genius-TCAD-Open/src/material/Si/Si_basic.cc:41)、[SiO₂ 默认 εr=3.9](D:/code-repo/Genius-TCAD-Open/src/material/SiO2/SiO2_basic.cc:45)。Charon 上述 nMOS 输入又显式将 Si 覆盖为 11.9，因此不能把库默认值当作某个算例最终生效的参数。
