"""Seal the bounded dielectric localization report and its reproducible evidence."""
import ast, math, re
from pathlib import Path
import validate_simplemos_remaining_dielectric_20260907 as t

a=t.a;d=t.d;l=t.l;OUT=t.OUT
REPORT=t.p.REPO/'docs/validation/simplemos_remaining_poisson_dielectric_localization_2026-09-07.md'

def main():
    for path in (l.OUT/'ledger_evidence.json',l.OUT/'interface_formula_evidence.json',OUT/'validation_evidence.json',OUT/'fields_evidence.json',OUT/'post_K_ledger/evidence.json',t.q.OUT/'final_evidence.json',t.v.OUT/'validation_evidence.json'):a.verify(path)
    rows=a.rows(OUT/'comparison.csv');small=a.rows(OUT/'small_dc.csv');finite=a.rows(OUT/'finite_dc.csv');response=a.rows(OUT/'response.csv');dual=a.rows(OUT/'dual_init.csv');ports=a.rows(OUT/'ports.csv');summary=a.read(OUT/'summary.json')
    passed=lambda rows:sum(r['qualified']=='True' for r in rows)
    maximum=lambda rows,key:max(float(r[key]) for r in rows)
    minimum=lambda rows,key:min(float(r[key]) for r in rows)
    pct=lambda x:f'{float(x)*100:+.6f}%'
    def row_for(c,axis):return next(r for r in rows if r['key']==c['key'] and r['axis']==axis)
    cases=sorted(a.read(OUT/'contract.json')['cases'],key=lambda c:(c['vg'],c['device'],c['vd']))
    table=['| NWell | Vg/V | Vd/V | 联合几何生产基线 | 再加原生 μ | 再改全部介电 K | 仅改 Si/SiO₂ 边 |','| --- | --- | --- | --- | --- | --- | --- |']
    for c in cases:
        r=row_for(c,'all_native_K');other=row_for(c,'Si_SiO2_only')
        table.append(f"| {c['device']} | {c['vg']} | {c['vd']} | {pct(r['production_relative_error'])} | {pct(r['mobility_relative_error'])} | {pct(r['candidate_relative_error'])} | {pct(other['candidate_relative_error'])} |")
    formula=a.rows(l.OUT/'interface_formula.csv');share=[r for r in a.rows(l.OUT/'interface_share.csv') if '_n23_' in r['key']]
    preflight=a.rows(OUT/'preflight.csv');jvp=[r for r in a.rows(OUT/'jvp.csv') if r['unchanged_SRH_cross_block']=='False' and float(r['step_V'])<1e-5]
    field=a.rows(OUT/'local_fields.csv');localtable=['| Vg/V | 节点 | Δψ 前 → 后 / μV | Δφn 前 → 后 / μV | Δφp 前 → 后 / μV | SRH 相对差前 → 后 |','| --- | --- | --- | --- | --- | --- |']
    for r in sorted((r for r in field if r['device']=='n23' and float(r['vd'])==1 and r['axis']=='all_native_K'),key=lambda r:(float(r['vg']),int(r['node_id']))):
        vals=[f"{float(r[f+'_old_delta_V'])*1e6:+.6f} → {float(r[f+'_new_delta_V'])*1e6:+.6f}" for f in ('psi','phin','phip')]
        localtable.append(f"| {r['vg']} | {r['node_id']} | {' | '.join(vals)} | {pct(r['old_SRH_relative'])} → {pct(r['new_SRH_relative'])} |")
    metrics=a.rows(OUT/'full_Si_field_metrics.csv');globaltable=[r'| NWell | Vg/V | Vd/V | Si 区 max\|Δψ\| / μV | max\|Δφn\| / μV | max\|Δφp\| / μV | SRH 相对 L2 差 |','| --- | --- | --- | --- | --- | --- | --- |']
    for c in cases:
        data={r['field']:r for r in metrics if r['key']==c['key'] and r['state']=='all_native_K'}
        vals=[f"{float(data[f]['max_absolute'])*1e6:.6g}" for f in ('psi','phin','phip')]
        globaltable.append(f"| {c['device']} | {c['vg']} | {c['vd']} | {' | '.join(vals)} | {float(data['SRH']['relative_l2'])*100:.6g}% |")
    worst=[]
    for c in cases:
        geo,mask=t.v.m.previous.prior.support(c);state=d.ordered(t.LOCAL/c['key']/'all_native_K/replacement/state.csv',geo.count)
        for field in ('phin','phip'):
            r=next(r for r in metrics if r['key']==c['key'] and r['state']=='all_native_K' and r['field']==field);i=int(r['worst_node'])
            worst.append(dict(key=c['key'],field=field,node_id=i,x_um=geo.coords[i][0],y_um=geo.coords[i][1],max_absolute_difference_V=float(r['max_absolute']),
                electrons_m3=float(state[i]['electrons_m3']),holes_m3=float(state[i]['holes_m3']),signed_Si_over_current_SRH_volume=float(geo.volumes['signed_si'][i]/geo.volumes['all_cell'][i]),native_SRH_quadrature_calibrated=False))
    a.write_csv(OUT/'worst_quasi_fermi_nodes.csv',worst)
    post=a.rows(OUT/'post_K_ledger/ledger.csv');postchecks=a.rows(OUT/'post_K_ledger/closure.csv');posttable=['| Vg/V | Vd/V | Poisson 电荷常数项 | ε₀ 项 | 连续性状态差项 | 同状态端口差 | 实际剩余 Id 误差 |','| --- | --- | --- | --- | --- | --- | --- |']
    for c in cases:
        if c['device']!='n23':continue
        data={r['component']:float(r['relative_to_native_Id']) for r in post if r['key']==c['key']}
        continuity=math.fsum(data[k] for k in ('electron_flux_state_difference','electron_SRH_state_difference','hole_flux_state_difference','hole_SRH_state_difference','continuity_other_and_roundoff'))
        vals=[data['Poisson_charge_constant_convention'],data['Poisson_epsilon0_convention'],continuity,data['mapped_port_minus_native_port'],data['actual_remaining_Id_error']]
        posttable.append(f"| {c['vg']} | {c['vd']} | {' | '.join(pct(x) for x in vals)} |")
    replay=a.rows(l.OUT/'ledger_checks.csv');screen=a.rows(l.OUT/'ledger.csv');constants=[]
    for c in cases:
        if c['device']!='n23':continue
        data={r['component']:r for r in screen if r['key']==c['key'] and r['state']=='native_mu'}
        constants.append(f"| {c['vg']} | {c['vd']} | {pct(data['candidate_epsilon0']['relative_to_Id'])} | {pct(data['candidate_charge_constant']['relative_to_Id'])} |")
    allrows=[r for r in rows if r['axis']=='all_native_K'];high=[r for r in allrows if r['device']=='n23']
    promotion='通过' if summary['all_native_K_promoted'] else '未通过'
    text=f'''# SimpleMOS 剩余 Poisson / 介电耦合误差定位

2026-09-07。**本轮把主要剩余误差定位到 Si/SiO₂ 共界面边的介电系数形成方式，并完成同扰动校准和双初始化自洽对照。** 在已合格的联合几何加原生有效迁移率基线上，全部原生介电几何替换后的高 NWell 四个工作点 Id 绝对相对误差范围为 {min(abs(float(r['candidate_relative_error'])) for r in high)*100:.6g}%–{max(abs(float(r['candidate_relative_error'])) for r in high)*100:.6g}%。本轮八点候选比较门槛{promotion}；没有修改正式默认值，没有扩展 16 工况或完整曲线。

承接[生产移植与原生有效迁移率响应验证](simplemos_production_migration_and_native_mobility_response_2026-09-07.md)。本轮使用既有 Sentaurus 2022 网格、原生 box 系数、物理场和电流；无新增上传下载或 sdevice 仿真。全部新增 DC 在 Vela 隔离程序中执行。

## 1. 定位到的系数差异

生产 `edgeEpsilon` 调用 `edgeAvgMaterialProp`，对相邻单元的 ε 作算术平均；`buildEdgeAssemblyKernels` 再乘全边 couple/length。对这里两侧各一个单元的界面边，记无量纲几何系数为 g：

```
K_old    = (ε_Si + ε_ox)/2 × (g_Si + g_ox)
K_native = ε_Si × g_Si + ε_ox × g_ox
K_old − K_native = (ε_Si − ε_ox) × (g_ox − g_Si)/2
```

n19、n23 各 20 条 Si/SiO₂ 共界面边的三个恒等式均通过：{passed(formula)}/40。旧边系数比原生几何系数高 {minimum(formula,'old_over_native_minus_one')*100:.6f}%–{maximum(formula,'old_over_native_minus_one')*100:.6f}%。本轮 K_native 使用 **Vela 的 ε₀**，因此该差异来自介电几何加权方式，不混入 ε₀ 的版本约定，也没有拟合材料参数。

这些界面边的 `region_local` 系数也与原生对应边一致。不过全器件还有原生对部分单元 box 几何的修改，不能据此把现有全局 `poisson_edge_coupling` 开关等同于本轮全部原生 K 方案。源代码见 [AssemblerUtils.h](../../include/vela/equation/AssemblerUtils.h) 和 [CoupledDDAssembler.cpp](../../src/equation/CoupledDDAssembler.cpp)；逐边证据见[界面公式核对](../../reference_tcad/simplemos_sentaurus2022/remaining_poisson_20260907/interface_formula.csv)。

## 2. 从空间归因到独立因果验证

先在八点的联合基线和原生 μ 替换合格态上重新求完整 DD 漏端伴随。高 NWell 的 20 条界面边贡献占全部介电修正带符号响应的 {minimum(share,'interface_over_total')*100:.4f}%–{maximum(share,'interface_over_total')*100:.4f}%。比例可超过 100%，因为 Si、Nitride 等其他区域存在抵消；这不是互斥误差百分比的相加。

两条候选方向分别为全部原生介电几何、仅 Si/SiO₂ 共界面边。均使用 `K(α)=K_old+α(K_target−K_old)`，每条边成对加入相反源项。迁移率固定为上一轮原生有效 μ，输运几何、三项 signed Si Poisson 电荷体积、SRH 原体积、Vela q/kb/ε₀ 和所有边界保持不变。

两个方向分别作 α=±0.001、±0.0005 的自洽计算，另有每点一次 α=0 重放；之后才作 α=1 双初始化。Poisson 缓存系数同时进入残差与解析 Jacobian；固定态连续性行、SG 边通量、端口电流逐值不变，直接端口项为零，预测为 `−λ_PᵀδR_P`。

- 原生 Poisson 重放和 Vela live-state 重放：{passed(replay)}/8；最大重放误差分别为 {maximum(replay,'native_replay_over_charge'):.6g}（相对电荷范数）和 {maximum(replay,'live_replay_relative'):.6g}（相对 Vela 映射残差）。使用算子导出的实际 ψ/n/p，避免六列输入与打包状态的混淆。
- 两方向固定态预检：{passed(preflight)}/16，最大 Poisson 源差分误差 {maximum(preflight,'Poisson_source_relative'):.6g}，门槛 1e-8。
- 受影响及非弱块 Jv 在两个较小步长的检查：{passed(jvp)}/{len(jvp)}，最大相对误差 {maximum(jvp,'relative_error'):.6g}，门槛 1e-4。弱 SRH 交叉块保持相同解析值；没有把消减噪声当作新的弱列认证。
- 小扰动 DC：{passed(small)}/{len(small)}；幅度校准：{passed(response)}/{len(response)}。最大预测误差 {maximum(response,'prediction_relative_error'):.6g}、双幅度差 {maximum(response,'two_amplitude_relative'):.6g}、偶/奇分量比 {maximum(response,'even_over_odd'):.6g}，最小信号/零漂移 {minimum(response,'signal_over_zero_drift'):.6g}。门槛依次为 0.001、0.001、0.01、100。
- 有限替换 DC：{passed(finite)}/{len(finite)}；双初始化：{passed(dual)}/{len(dual)}；独立端口提取一致性：{passed(ports)}/{len(ports)}。全部 DC 最大逐行比值 {maximum(small+finite,'max_row_ratio'):.6g}，门槛 1e-6；最大 KCL/Id {maximum(small+finite,'kcl_over_Id'):.6g}，门槛 1e-8。

见[冻结合同](../../reference_tcad/simplemos_sentaurus2022/remaining_poisson_20260907/dielectric/contract.json)、[幅度校准](../../reference_tcad/simplemos_sentaurus2022/remaining_poisson_20260907/dielectric/response.csv)、[双初始化](../../reference_tcad/simplemos_sentaurus2022/remaining_poisson_20260907/dielectric/dual_init.csv)。

## 3. 自洽 Id 对比

误差统一为 `Id_Vela / Id_Sentaurus − 1`；电流 A/μm，宽度 1 μm。后两列都以“联合几何加原生 μ”为起点。

{chr(10).join(table)}

全部 K 方案相对“仅加原生 μ”和“联合几何生产基线”的高 NWell 绝对误差、低 NWell 不恶化与高低 log 电流配对条件，逐项记录在[配对表](../../reference_tcad/simplemos_sentaurus2022/remaining_poisson_20260907/dielectric/pairing.csv)。界面子集只用于定位，不能推广为任意器件的生产策略。有限替换响应与初始切线预测的差异也保留在[电流表](../../reference_tcad/simplemos_sentaurus2022/remaining_poisson_20260907/dielectric/comparison.csv)，不以小扰动线性校准替代有限幅度验证。

## 4. 剩余势与 SRH 场差

下表是 n23、Vd=1 V，原生 μ 基线到全部 K 替换的带符号差；势单位 **μV**，不是 mV。

{chr(10).join(localtable)}

此外对全部 907 个自由 Si 节点作直接原生物理场比较。下表是全部 K 替换后的范围统计，场误差是诊断结果，没有新增或放宽接受门槛：

{chr(10).join(globaltable)}

SRH 由实际装配源、保持不变的 SRH 体积和边通量物理单位转换还原为 cm⁻³ s⁻¹；原生侧为节点 Plot 速率。相对 L2 是节点速率误差范数除以原生速率范数，不能当作原生求解器导出的积分残差。32 项局部势差链与直接原生比较、SRH 单位转换的独立核对见[场核对表](../../reference_tcad/simplemos_sentaurus2022/remaining_poisson_20260907/dielectric/field_reconciliation.csv)。全部字段见[Si 场统计](../../reference_tcad/simplemos_sentaurus2022/remaining_poisson_20260907/dielectric/full_Si_field_metrics.csv)和[选定节点](../../reference_tcad/simplemos_sentaurus2022/remaining_poisson_20260907/dielectric/selected_node_fields.csv)。

**不能把节点 1000/1009 的微伏级改善推广到全部准费米场。** n23、Vg=0.8 V、Vd=1 V 的最差 φn/φp 分别位于节点 792/1057，差约 3.804/3.901 mV。节点 792 的 n/p 仅约 5.85e8/7.55e9 m⁻³，是强耗尽区域；节点 1057 的 n/p 约 3.71e19/3.92e5 m⁻³，空穴为少数载流子。其 signed Si 体积与当前 SRH 体积之比分别约 1.1796/0.2173。这给出剩余连续性/SRH 几何的明确对照位置，尚不证明直接换 SRH 体积正确；[最差节点表](../../reference_tcad/simplemos_sentaurus2022/remaining_poisson_20260907/dielectric/worst_quasi_fermi_nodes.csv)保留所有八点。

## 5. 仍未归零的部分及后续边界

在全部 K 替换后的合格态重新求伴随，并对原生映射态导出完整残差，建立 `Id_strict−Id_native ≈ −λᵀ(R_mapped−R_strict)+(Id_mapped−Id_native)` 的一阶账本。这里保留接触边界行；映射原生状态可能不满足 Vela 的边界约定，不能沿用固定边界参数扰动时将接触权重置零的做法。

下表各贡献均除以原生 Id。它描述实际状态差的一阶分配，**不是新常数候选的正负自洽校准**；表中未单列的热密度映射、边界、重放及有限状态余项保留在完整账本中。

{chr(10).join(posttable)}

八点一阶账本与实际余差的最大差为 {max(abs(float(r['finite_state_remainder_over_Id'])) for r in postchecks):.6g}×Id，最大占当前余差 {max(abs(float(r['finite_state_remainder_over_gap'])) for r in postchecks)*100:.6g}%。这将最后的端口余量约束到 Poisson 常数约定和同状态输运/端口差异的组合；SRH 状态差对 Id 的投影很小，但局部少数载流子场仍值得独立验证。见[剩余误差账本](../../reference_tcad/simplemos_sentaurus2022/remaining_poisson_20260907/dielectric/post_K_ledger/ledger.csv)和[闭合检验](../../reference_tcad/simplemos_sentaurus2022/remaining_poisson_20260907/dielectric/post_K_ledger/closure.csv)。

本轮分开保留物理常数约定的切线筛查：Vela q=1.602176634e-19 C、ε₀=8.8541878128e-12 F/m；原生重放采用既有 2022 手册示例常数 q=1.602192e-19 C、ε₀=8.8542e-12 F/m。八点原生 Poisson 重放支持该约定，无拟合。以下数值是在替换介电 K **之前**的原生 μ 合格态上投影，仅用于量级与方向判断：

| Vg/V | Vd/V | 单独 ε₀ 约定的预测 ΔId/Id | 单独 Poisson 电荷 q 约定的预测 ΔId/Id |
| --- | --- | --- | --- |
{chr(10).join(constants)}

这两项没有做本轮同方向自洽校准，不能据此声称常数修改已经作为独立物理候选通过，或已经消除了最后的 Id 差。尤其 q/kb/热电压约定同时涉及密度映射、SG 通量和端口单位，不能只全局替换一个常数后按 Id 拟合。SRH 节点场仍保留差异，不能因 Poisson/Id 改善就宣布 SRH 离散完全一致。

后续应利用本轮合格全部 K 状态的伴随，对一致的常数约定方向作正负自洽校准，并在节点 792/1057 对剩余 SRH 源及其体积另做独立校准；生产化则需把逐材料介电几何来源与其他 Poisson/电荷诊断分支一致地接入，验证原生 box 修改的适用范围。当前任务完成的是八点误差定位与隔离验证，没有把原生导出的逐边表设为全局默认，也没有启动 16 工况或逐 Vg 扫描。

## 6. 环境、复现与证据

Windows MSYS2 UCRT64，C++20 Release，实际 DC 后端 Eigen SparseLU/COLAMD；依赖构建中 HDF5/TDR、UMFPACK、SPQR 可用，不代表本次 DC 使用 UMFPACK。n19/n23：Boron 1e17/2e17 cm⁻³，节点 1480/1482，三角形 2742/2746；300 K，Boltzmann/no-BGN，matched ni=1.0750038488844236e10 cm⁻³，Masetti 总杂质与掺杂相关 SRH；HFS、表面迁移率、Auger、雪崩、DG 关闭。网格 μm，状态密度 m⁻³。

max_iter=200、reltol=1e-7、abstol=1e-12、stall floor=1e-9，ψ/准费米步长上限 0.35/0.025 V，contact_basin、原标量线搜索，线性迭代修正 4 次；1814 载流子行 eps=1e-6，无行排除或新尺度下限。全局连续性闭合 tolerance=1e-6、source_floor=1e-10；低于 source_floor 的净源不声称任意相对精度。双初始化势差 ≤1e-6 V、密度相对差 ≤1e-4、Id 相对差 ≤1e-6，保持原门槛。

新增隔离编译成功，几何守恒、固定态残差/Jacobian、同源正负响应、有限 DC 和字段核对如上。Python 脚本语法检查通过。本轮未修改正式源码、未重跑全套 CTest；上一轮 754/770 与 16 项历史证据失败继续保留，没有更改旧哈希或重新标记通过。

复现入口：[残差账本脚本](../../scripts/localize_simplemos_remaining_poisson_20260907.py)、[界面公式审计](../../scripts/audit_simplemos_dielectric_interface_formula_20260907.py)、[冻结自洽验证](../../scripts/validate_simplemos_remaining_dielectric_20260907.py)、[完整 Si 场比较](../../scripts/analyze_simplemos_remaining_poisson_fields_20260907.py)。这些脚本和输出按阶段冻结，已有目录拒绝覆盖；新重跑应使用新实验目录并记录新合同。

本报告与全部阶段证据关联在[最终证据](../../reference_tcad/simplemos_sentaurus2022/remaining_poisson_20260907/final_evidence.json)。本轮未提交或合并代码。
'''
    assert not REPORT.exists();REPORT.write_text(text,encoding='utf-8',newline='\n')
    scripts=[Path(__file__).resolve()]+[t.p.REPO/'scripts'/name for name in ('localize_simplemos_remaining_poisson_20260907.py','validate_simplemos_remaining_dielectric_20260907.py','audit_simplemos_dielectric_interface_formula_20260907.py','analyze_simplemos_remaining_poisson_fields_20260907.py','audit_simplemos_post_dielectric_residual_20260907.py')]
    for script in scripts:ast.parse(script.read_text(encoding='utf-8'))
    for link in re.findall(r'\]\(([^)]+)\)',text):
        target=(REPORT.parent/link.split('#')[0]).resolve()
        if target.name!='final_evidence.json':assert target.is_file(),target
    files=scripts+[REPORT,l.OUT/'ledger_evidence.json',l.OUT/'interface_formula_evidence.json',OUT/'validation_evidence.json',OUT/'fields_evidence.json',OUT/'post_K_ledger/evidence.json',OUT/'worst_quasi_fermi_nodes.csv',t.q.OUT/'final_evidence.json',t.v.OUT/'validation_evidence.json']
    a.write(l.OUT/'final_evidence.json',dict(status='completed_bounded_localization',date='2026-09-07',new_DC=len(small)+len(finite),qualified_DC=passed(small)+passed(finite),
        production_changes=False,acceptance_changes=False,full_curve_started=False,summary=summary,input_hashes={a.rel(f):a.sha(f) for f in files}))
    a.verify(l.OUT/'final_evidence.json');print('Report and evidence sealed',REPORT,flush=True)

if __name__=='__main__':main()
