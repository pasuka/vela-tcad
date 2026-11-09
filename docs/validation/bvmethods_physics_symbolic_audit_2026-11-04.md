# 当前BV算例物理模型与电势/浓度偏导审计（2026-11-04）

## 结论和资格范围

本轮完成模型清单、手工符号推导、34项精确有理多项式恒等式检查，以及92项当前C++原语值/偏导对照。后者最大缩放误差1.09912e-14（误差定义为abs(actual-reference)/max(1,abs(reference))，并非所有弱分量的相对误差）。16个被审计物理/装配/求解源码与最近VM Full-only N20任务冻结清单逐SHA相同。这只证明所列16文件与该冻结版本一致，不证明传递源码/编译身份闭包，也不自动覆盖其他运行。

**本轮人工符号核对未发现所列公式的符号/链式错误；数值检查只覆盖下列92项局部原语，不授予装配后7801维Jacobian的光滑分支一致性资格。** 这不等于所有状态/分支/完整稀疏Jacobian已被证明正确，也不证明与Sentaurus内部模型一致。没有修改生产源码、材料、保护、验收门或默认导数策略，没有启动新的Newton/Sentaurus计算。

确认一项归档元数据讹误：兼容常数材料JSON顶层physical_constants_profile为sentaurus_u2022_12，Si的provenance文字仍称CODATA2018。它不能证明数值用错；已核验的Nc/Nv等数值和二进制配置是兼容profile。当前转换器源码已经按profile分支生成说明，故不能把归档旧说明直接认定为当前转换器仍有同样缺陷。本轮只记录，不修改冻结归档。

## 1. 当前算例、材料和实际分派

使用最近Full-only N20冷0–1V任务的已核验配置/网格/材料，物理与此前automatic A0相同。2600节点、5020 Tri3：

|材料|区域|单元数|相对介电常数|载流子方程|
|---|---|---:|---:|---|
|Si|R.Substrate|4156|11.7|电子、空穴|
|SiO2|R.Gateox及两个R.PolyReox区域|548|3.9|无，非活性槽位约束为0|
|Nitride|两个R.Spacer区域|316|7.5|无，非活性槽位约束为0|

300K、经典等温DD；启用Joyce–Dixon FD近似、OldSlotboom及静态统计校正、Masetti(total_impurity)+Enhanced Lombardi+Canali、掺杂相关SRH、自洽Van Overstraeten GradQF AV。BTBT关闭，Auger未启用，无量子势、热方程或动态陷阱占据方程。源/漏/衬底为欧姆接触，gate为metal_gate。串联电阻R=1e7 Ω·µm。几何配置含node_volume_policy=barycentric，但当前Poisson电荷被poisson_charge_node_volume=signed_transport覆盖：电子、空穴、掺杂三项均使用同一套带符号AverageBox体积。连续性、SRH、AV及迁移率平均沿用cellBoxNodeMeasures/相应单元box支持，结合delaunay_transfer后的非负边系数；它们不能统称为Poisson的barycentric体积。不套用PN2D模板另一个atomic profile的mixed-Voronoi默认。

兼容常数q=1.602192e-19 C，kB=1.380662e-23 J/K，eps0=8.8542e-12 F/m，VT=kBT/q=0.025851995266484913V。300K材料Nc=2.8582520073866568e19、Nv=3.1046289324748386e19 cm^-3；材料基础ni=1.4638805412559162e10 cm^-3，节点实际ni_eff还包含掺杂BGN。

当前生产packed变量为(psi/V0,n_internal,p_internal)，外接未知量lambda=Vd为物理V；浓度internal为本算例cm^-3。全变量Unit原型是独立路径，不把其默认资格混入本审计。以下先用物理变量(psi,n,p,lambda)推导，最后作单位/坐标变换。

## 2. 连续方程与符号约定

令Nd、Na、介电常数、材料参数和温度为固定数据，C=Nd-Na，rho=q(p-n+C)+rho_fixed，电场E=-grad(psi)。常规电流为Jn、Jp：

\[
-\nabla\cdot(\epsilon\nabla\psi)-\rho=0,\quad
-\nabla\cdot J_n/q+U-G=0,\quad
+\nabla\cdot J_p/q+U-G=0.
\]

在均匀材料的光滑连续极限，Jn=-q mu_n n grad(phi_n)，Jp=-q mu_p p grad(phi_p)。离散代码以SG粒子通量而非简单中心差分实现；电子残差对应负电子常规电流散度，空穴残差对应正空穴常规电流散度。因此两行均加SRH、减AV是相容的，不能只看电子/空穴同号源项就认定符号错。

## 3. 统计、BGN与准费米势反变换

令r_n=n/Nc，r_p=p/Nv，H为当前JD反函数，eta_n=H(r_n)、eta_p=H(r_p)，k_n=H'(r_n)/Nc、k_p=H'(r_p)/Nv。再令F为与H配对的正向函数，则k_n=1/[Nc F'(eta_n)]。

\[
\phi_n=\psi-V_T[H(n/N_c)-\ln(n_{i,e}/N_c)],\quad
\phi_p=\psi+V_T[H(p/N_v)-\ln(n_{i,e}/N_v)].
\]

在材料、掺杂和温度固定时：

\[
\partial_\psi\phi_n=\partial_\psi\phi_p=1,\quad
\partial_n\phi_n=-V_T k_n,\quad \partial_p\phi_p=+V_T k_p,
\quad\partial_p\phi_n=\partial_n\phi_p=0.
\]

这与CarrierStatistics.cpp及Tri3人口jet的正负号一致。浓度是主变量并不意味着不再需要QF：QF仍用于GradQF迁移率、AV、SG及欧姆边界的物理构成。

JD稀分支实际为H(r)=ln r+P(r)，P=a1 r+a2 r²+a3 r³+a4 r⁴，系数分别0.3535533905932738、-0.004950089729875255、0.00014838577128872333、-4.425630118996707e-6。令K=1+rP'(r)，则

\[
F'(\eta)=r/K,\quad F''(\eta)=r[K-rK']/K^3.
\]

高密度反分支H=sqrt(u²-c)-shift，u=(r/A)^(2/3)，A=4/(3sqrt(pi))，c=pi²/6，H'=(2/3)u²/[r sqrt(u²-c)]。连接为C1，F二阶导在拼接点可跳变。跨拼接点不能按同一光滑二阶公式验收，应使用分支内或单侧检查；本轮未对拼接附近完成方向导数收敛检验。

BGN实际节点路径：

\[
\Delta E_g(N)=\max\{0,E_0+E_{ref}[\ell+\sqrt{\ell^2+s}]\},\quad \ell=\ln(N/N_{ref}),\ N=N_D+N_A.
\]

Eref=.009eV、Nref=1e17cm^-3、s=.5；静态FD校正为VT乘以供受体各自max(0,H(N_species/N_DOS)-ln(N_species/N_DOS))。最终ni_eff=ni_material exp(DeltaEg_total/(2VT))。材料dEg0(OldSlotboom)已在转换阶段处理，与节点掺杂BGN区分。虽然通用deltaEg函数签名包含n,p，**本算例buildEffectiveNodeNi调用传入n=p=0**，故ni_eff对当前coupled Newton psi,n,p偏导均为0；仅在这一调用路径中不加BGN人口偏导不是漏项。Gummel路径可能传入n、p，其模型语义另计，不推广此结论。

## 4. Poisson离散行：耦合状态与初始化不可混同

物理自由行（忽略统一正行尺度）：

\[
R_{\psi i}=\sum_j k_{ij}(\psi_i-\psi_j)+q(W_{n i}n_i-W_{p i}p_i-W_{D i}C_i)-Q_{fixed,i}.
\]

k_ij从各材料单元介电系数构造。通式保留Wn/Wp/WD记号，但**当前signed_transport配置中Wn=Wp=WD=W_signed**；实际需区分的是Poisson带符号电荷体积与连续性/源项体积。W_signed不假定为正：若节点总带符号体积为负，人口列符号随之翻转；若只是单元局部份额为负，不自动意味着节点总和为负。本轮未给出节点总负体积数量或最小值，相关统计仍待做。网格26个钝角单元已按点积与舍入界独立确认，但这不能替代节点体积总和统计。偏导：

\[
J_{\psi_i\psi_i}=\sum_jk_{ij},\quad J_{\psi_i\psi_j}=-k_{ij},\quad
J_{\psi_i n_i}=+qW_{n i},\quad J_{\psi_i p_i}=-qW_{p i}.
\]

全耦合人口变量中，对psi求导固定n,p，因此没有额外Debye项。BvMethodsCoupledRuntime::populationBase直接这样装配，而非从QF矩阵相减重建。

冷Poisson-only阶段不同：固定phi_n/phi_p，n,p由psi重建，n_psi=Nc F'(eta_n)/VT、p_psi=-Nv F'(eta_p)/VT，因此对角上有额外统计链式项q[Wn Nc F'+Wp Nv F']/VT；当前带符号W不保证它为正。初始化确实转入QF统计的Poisson-only路径；不能用这一块替代后续人口Jacobian的固定浓度偏导。

## 5. SG输运的完整偏导

对三角形有向边i->j，载流子c=n或p，s_n=+1、s_p=-1。设l=ln(c_j/c_i)，eta=H(c/N_DOS)，广义Einstein系数

\[
g=(\eta_j-\eta_i)/l,\quad
\delta=s(\phi_j-\phi_i)/(V_Tg),\quad a=l+\delta,\quad
K_e=C_{edge}\bar\mu V_Tg.
\]

当前输运的mu_bar为固定单元box面积加权的顶点迁移率，不是AV的代表顶点迁移率。物理粒子通量为

\[
\Phi=K_e[c_iB(-a)-c_jB(a)],\quad B(a)=a/(e^a-1).
\]

等价稳定值路径为K_e c_j B(a) expm1(delta)，由B(-a)=exp(a)B(a)得到。近平衡不能直接用两个大Bernoulli项的普通双精度相减作为独立真值。

\[
B'(a)=\frac{e^a-1-ae^a}{(e^a-1)^2},\quad B'(0)=-1/2.
\]

令T=c_iB(-a)-c_jB(a)，对任意物理未知量u：

\[
\Phi_u=T(K_e)_u+K_e\{B(-a)(c_i)_u-B(a)(c_j)_u-[c_iB'(-a)+c_jB'(a)]a_u\},
\]
\[
(K_e)_u=C_{edge}V_T(g\bar\mu_u+\bar\mu g_u),\quad
 a_u=l_u+\frac{s\Delta\phi_u}{V_Tg}-\frac{s\Delta\phi}{V_Tg^2}g_u.
\]

电势列l_psi=g_psi=0，a_psi_i=-s/(VTg)、a_psi_j=+s/(VTg)，但mu_psi通常不为0。浓度列必须同时保留显式c项、QF人口反变换、g和迁移率链式项。

\[
g_{c_i}=\frac{-k_i l+(\eta_j-\eta_i)/c_i}{l^2},\quad
g_{c_j}=\frac{k_j l-(\eta_j-\eta_i)/c_j}{l^2}.
\]

当abs(l)<=1e-10，代码使用mid_eta处max(1,F/F')极限分支，其对mid_eta偏导为1-F F''/(F')²（值大于1时），再乘eta人口链式项。该小分支不是在所有密度区都强行令g=1。

每条边将同一个Phi以+/-加入两节点，因此自由输运内部边对残差和Jacobian的总和均消去。单元组装后还需用各节点是否为接触行作边界替换，不能把边界行拿来做自由连续性守恒求和。

## 6. 迁移率及对psi/n/p的链式项

Masetti使用固定N=Nd+Na：mu_bulk=mu_min1 exp(-Pc/N)+(mu0-mu_min2)/(1+(N/Cr)^alpha)-mu1/(1+(Cs/N)^beta)，并取非负值。固定掺杂/T下其psi,n,p偏导为0；本算例未用PhuMob的载流子散射依赖。

Enhanced Lombardi实际alpha=0。令Fperp=abs(grad psi dot grad distance)*坐标换算，d=exp(-distance/Lcrit)，Z=C_T,N Fperp^(2/3)，则逆散射项

\[
I(F)=d\left[a_{ac}\frac{F}{B+Z}+a_{sr}\left(\frac{(F/100)^A}{\Delta}+\frac{F^3}{\eta}\right)\right],\quad
\mu_s=(1/\mu_{bulk}+I)^{-1}.
\]

\[
I'=d\left[a_{ac}\frac{B+Z/3}{(B+Z)^2}+a_{sr}\left(\frac{A(F/100)^A}{\Delta F}+\frac{3F^2}{\eta}\right)\right],\quad
(\mu_s)_F=-\mu_s^2 I',\quad (\mu_s)_{\mu_b}=\mu_s^2/\mu_b^2.
\]

距离、掺杂、几何投影是固定数据，不能误加其Newton人口导数。法向场使用距离场梯度，代码未假设该梯度已单位归一；psi的偏导按实际几何梯度计算。

高场Canali：F为实际单元GradQF范数，接触相邻规定单元使用电势梯度作为驱动。F<1V/cm时返回mu_s；其他分支令z=mu_s F/vsat、T=1+z^beta：

\[
\mu=\mu_sT^{-1/\beta},\quad
\mu_{\mu_s}=T^{-1-1/\beta},\quad
\mu_F=-\mu_s^2 z^{\beta-1}T^{-1-1/\beta}/v_{sat}.
\]

完整状态偏导为mu_u=mu_mu_s (mu_s)_Fperp (Fperp)_u+mu_F F_u。若F=norm(A phi)，非零时F_u=unit(A phi) dot A phi_u。psi/n/p链式项由第3节给出；alpha=0且Masetti静态时电子mu对p偏导为0、空穴mu对n偏导为0。实际C++同时保留法向场和GradQF驱动这两项。

风险：1V/cm处为刻意采用的分段阈值，存在极小值跳变；Fperp=0的abs及零梯度范数不是普通可微点，代码取零广义导数。不能在跨分支/尖点的有限差分失败后直接认定光滑公式错。

## 7. 广义FD-SRH：不能只用Boltzmann商式

掺杂寿命tau_c=tau_min+(tau_max-tau_min)/(1+(N/Nref)^gamma)，本算例tau_max=1e-7s、tau_min=0、Nref=1e16cm^-3、gamma=1。N和T固定，因此tau的psi,n,p偏导为0。

令gamma_n=(n/Nc)exp(-H(n/Nc))、gamma_p同理；稳定原语使用log gamma，不在低密度用两个巨大倒数相减。设

\[
\gamma_n'=\gamma_n(1/n-k_n),\quad \gamma_p'=\gamma_p(1/p-k_p),
\]
\[
X=np-\gamma_n\gamma_p n_{i,e}^2,\quad
D=\tau_p(n+\gamma_n n_{i,e})+\tau_n(p+\gamma_p n_{i,e}),\quad U=X/D.
\]

\[
U_n=\frac{p-n_{i,e}^2\gamma_p\gamma_n'}{D}-\frac{X\tau_p(1+n_{i,e}\gamma_n')}{D^2},
\]
\[
U_p=\frac{n-n_{i,e}^2\gamma_n\gamma_p'}{D}-\frac{X\tau_n(1+n_{i,e}\gamma_p')}{D^2},\quad U_\psi=0\ \text{(固定n,p)}.
\]

源码RecombinationModel提供固定X和固定gamma的偏导，Runtime再组装X_n/X_p与gamma链式项；将其dRateDn单独当作总导数会误报缺项。populationRecombinationExcess使用补偿乘积差，避免旧QF splitting经expm1溢出。极值饱和和abs(D)<1e-100分支导数须按返回值分段解释，本轮普通分支符号推导不证明所有饱和边界。

## 8. AV产生率及全部状态偏导

300K分支alpha_c(F)=A_c exp(-B_c/F)，F>0时alpha'=alpha B/F²，F<=0返回0。参数electron low/high A=703000cm^-1、B=1231000V/cm；hole low A=1582000、B=2036000，high A=671000、B=1693000。临界Fsw=400000V/cm附近[.95Fsw,1.05Fsw]采用余弦混合：

\[
\alpha=(1-w)\alpha_l+w\alpha_h,\quad
\alpha'=(1-w)\alpha_l'+w\alpha_h'+w'(\alpha_h-\alpha_l).
\]

代码确实包含最后的混合权重导数，不能只检查两端指数项。300K温度guard使phonon温度因子为参考温度值，不把T变体视为已验证。

AV的alpha驱动场与迁移率驱动共享既定stencil选择：接触相邻规定单元同样从GradQF切换到grad(psi)，不只是迁移率如此。当前AV电流通过代表顶点的两条QF-SG边恢复二维向量，再乘固定几何投影。其场F及代表顶点迁移率mu_rep均保持实时；输运行使用顶点迁移率平均，因此不能以同一mu_bar替代这里mu_rep。令j_c为不含mu/q的恢复向量：

\[
G=\sum_{c=n,p}\alpha_c(F_c)\mu_{c,rep}\|j_c\|,
\]
\[
G_u=\sum_c\left[(\alpha_c'F_{c,u}\mu_c+\alpha_c\mu_{c,u})\|j_c\|+
\alpha_c\mu_c\frac{j_c\cdot j_{c,u}}{\|j_c\|}\right].
\]

j的每条SG边导数仍按第5节、使用该AV原语自己的系数和g=1来求，再作固定矩阵线性变换。AV行按box面积积分，在电子和空穴残差均为-G，其Jacobian均为-G_u，q在终端电流换算阶段出现，不能在粒子源中重复乘q。

**零电流处norm无Frechet导数。** 代码在j=0选择零广义切线；真实方向导数为norm(dj)，不是零线性映射。近零时j方向对QF最后几位舍入敏感，是此前六个空穴行偏差的未闭合模型/数值约定之一。92项原语检查包含零切线实现选择，但绝不把它计为零点普通可微性证明。与Sentaurus内部AV切线是否完全一致仍未确认，不能据此全局翻转AV符号。

## 9. 接触/界面、电路与Schur符号

欧姆接触同时规定psi=Vcontact+psi_built_in、phi_n=phi_p=Vcontact。换成人口主变量后，约束行仍通过QF反变换表达，并非没有浓度边界：

\[
b_n=(\phi_n-V_c)/V0,\ b_p=(\phi_p-V_c)/V0,
\]

其psi列均1/V0，人口列分别-VT*k_n/V0、+VT*k_p/V0，Vcontact列均-1/V0。电势边界行同样为(psi-psi_bi-Vc)/V0。金属gate仅施加电势，其work-function偏置先进入有效电势；不在氧化物中凭空施加输运n/p。界面通过共享电势、各单元eps和积分电荷装配；本案未启用动态界面陷阱。

终端电流使用未替换边界行的**transport-only**节点通量，符号I=scale*q*sum_contact(-r_n_transport+r_p_transport)。其人口/电势梯度来自同一输运Jacobian，不把AV/SRH积分再次作为接触传导电流相加。在固定完整packed态时I_lambda=0；沿解支dI/dlambda仍非零，来自x(lambda)的响应。

\[
L(x,\lambda)=c_V\lambda+R I(x)-V_{outer},\quad
L_x=R I_x,\quad L_\lambda=c_V+R I_\lambda.

\]

本算例c_V=voltageCoefficient=1、固定packed态I_lambda=0。一般完整Newton块为[[J,F_lambda],[R I_x,c_V+R I_lambda]]，RHS为[-F,-L]。a=-J^-1 F，z=J^-1 F_lambda时：

\[
s=c_V+R I_\lambda-RI_xz,\quad
\delta\lambda=(-L-RI_xa)/s,\quad\delta x=a-z\delta\lambda.
\]

符号与CoupledLoadLine一致，方向不能只看s的符号。当前人口F_lambda列实际由中心差分边界残差取得；相关makeArclengthSystem默认及本次既有探针调用h=1e-4 V，差分分母2h。应逐具体运行核对调用参数，不把默认视为所有路径固定值；解析预期只在受激接触行有-1/V0。本次核查确认其数学结构，但**没有新增对全部真实状态的该列近根绝对精度资格**，不能把所有端口/电路弱列问题宣布解决。

## 10. 坐标、行尺度和近似策略的界限

物理到当前packed的电势列需乘V0；浓度局部jet为适应末尾V0散射，先种入dc/V0，最后乘V0，两个因子必须成对。直接列方向和极低人口相对方向在局部抵消后再撤销，不改变残差状态定义。

若从当前packed再作全变量Unit，D=diag(1,C0,C0,V0)，残差行变换为S，则Jhat=S J D、bhat=S b；从完全物理变量出发的D的首块才是V0，不能重复除电势尺度。负载行、端口I/I0、电阻R I0/V0须同步。固定预条件/求解坐标尺度不当成新的材料物理偏导。

bootstrap保留完整物理RHS、完整SRH，但冻结迁移率导数、去掉AV导数；floor_adjacent_edges还冻结相关边的统计人口反馈。这是**有意的近似Jacobian**，不等于完整残差的精确导数，不能拿Full符号检查通过去替bootstrap资格背书。另一个细节是电路端口I_x使用完整transport梯度，设备块在bootstrap时却可近似，Schur组合是近似Newton系统。偏导审计不能证明它在高压首步给出可靠方向。

Full模式数学光滑分支检查通过，并不意味着全程稳健；最近N20 Full-only OR在0.262404V失败，AND成本增加。此前Unit固定物理方向29/30门缺口、原生首帧六行差异及报告型守恒未满足均继续保留。

## 11. 验证、发现和限制清单

|对象|本轮结论|
|---|---|
|固定人口Poisson电荷/介电项、Poisson-only链式项|符号与装配一致|
|JD低/高分支、QF人口反变换、gamma|推导一致，光滑分支原语抽样通过|
|SG Bernoulli正负号、显式人口与a/mu/g链式项|代数推导核对；数值只测avSgFlux的g=1、常系数原语，未测完整transport装配|
|Masetti静态掺杂、Lombardi、Canali链式|当前300K/alpha0分支一致；阈值/abs非光滑需单列|
|广义SRH|人工核对X/D完整链式；C++仅数值测X，未测完整商式|
|AV指数及混合导数、乘积和范数链式|非零光滑分支一致；零点约定与原生内部等价未闭合|
|bootstrap|确定性省略导数，不能声称等于Full残差Jacobian|
|接触、端口、电路、Schur|结构与符号一致；F_lambda数值列完整精度资格未新增|
|冻结材料Si provenance|确认文字与顶层profile矛盾，未据此判数值错|
|完整7801维所有状态Jacobian/原生等价|本轮不授予这一总资格|

34项符号检查使用标准库有理多项式表示、商式/链式微分及交叉相乘零多项式验证，未用浮点采样拟合；不是依赖缺失的外部CAS。92项C++检查直接包含当前生产头文件，但覆盖必须按第12节的实际分类解释。参考式同样来自本轮推导，存在共模风险。未运行完整输运装配、SRH完整商式、Poisson/边界行、全局Jv或列和不变量测试；未冒称全Ctest/完整网格MP检查。

探针首版用双精度两个大SG项相减作参考，在近平衡处出现参考消去错误，失败源码保留；改成长精度expm1稳定参考后92项通过。工具链首个调用缺少UCRT64 PATH而未编译，随后采用仓库规定PATH；两者均不计为生产数值失败。

证据位于ignored `build-release/bvmethods_physics_symbolic_audit_20261104`：case_inventory.json、source_identity.json、symbolic_checks.py及结果、primitive_derivative_audit.cpp及结果、输入/源码SHA与审计包回执。原始算例归档SHA为9a0294faaffb05c7b4a0debdfa5cce2ae547f9290b93d8e7a38ef5485f7fe36c。未修改求解代码或冻结VM材料。


## 12. 范围更正：92项、误差口径和装配资格

本节收窄此前较强措辞。92项实际为：Bernoulli16、JD16、avSgFlux(g=1且k常数)24、AV系数10、Canali12、Lombardi9、SRH乘积X2（gamma切线手动给定）、norm3。SG小delta测试固定l=ln(1.2)，**没有测试小abs(l)**。34项符号恒等式的exp/log等链式因子由人工输入；它们验证代数推导，不独立认证源码和装配。

未数值覆盖：evaluateBvMethodsTri3Transport、广义g及其近等密度分支、K_e完整链、迁移率平均、场stencil/dualAbs、AV向量恢复、SRH完整商式、Poisson/接触行装配、行列散射/V0缩放、F_lambda以及输运列和为零。此前将这些推导统称为“原语对照通过”的表述撤回。

原误差abs(a-ref)/max(1,abs(ref))与5e-11门会弱化小量验证。Canali场切线可能仅有5e-4至5e-2相对约束；Lombardi约1e-11切线低于该绝对门，B(60)约5e-25也不能据此认为有效相对精度通过。AV低场下溢为0属于平凡通过，须单列。**1.09912e-14是该混合口径的最大值，不能代表所有小量的相对精度。** 尚未生成逐项自然尺度误差与非平凡覆盖统计，不将其写成已完成。

新增数值风险：g=(eta_j-eta_i)/l在abs(l)>1e-10时使用两个独立eta的差。近等密度时，eta差的舍入误差可被1/l、偏导中的1/l²放大；其风险与极低人口floor尾部强切线不同。误差估计依赖eta、驱动、边几何与实际状态，不把示例0.2~0.3作为当前网格已测结果。真实冻结态边abs(l)直方图、小l扫描及高精度输运Jv仍未完成。

模型语义待核实：Lombardi使用未归一化grad(d_h)；接触相邻AV的alpha也使用grad(psi)；AV临界带余弦混合是否与原生一致未核实；AV恢复电流g=1而输运用广义g；原生Fermi-SG离散、Enormal、BGN Fermi修正仍须同版本手册核实。这些与“代码对其自身函数求导一致”不同。

源码身份原16项中已包含AssemblerUtils.h，但未建立传递闭包；缺失PopulationSg/PrimitiveState/FermiSgDifferential/ElementLombardi/MobilityModel/DerivativeRecovery及几何实现等依赖。Full-only冻结清单不能在没有另行证据时直接替代automatic/bootstrap运行身份。独立探针曾使用C++20/-O2、兼容常数宏和UCRT64 include路径，但尚未归档所有传递依赖与完整工具链/编译记录，原16项SHA不授予这一闭包资格。

建议的装配层多态高精度Jv、输运列和、缺失原语及自然尺度误差、小l直方图、signed AverageBox/连续性体积差与grad(d_h)分布、同版本手册核实和传递源码闭包，均列为后续待办，本报告不把它们当已授权执行或已完成。生产门限和代码未改。
