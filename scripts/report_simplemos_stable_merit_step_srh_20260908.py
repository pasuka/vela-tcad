"""Seal the three isolated diagnostic batches without changing their contracts."""
import ast
import re
from pathlib import Path

import validate_simplemos_stable_merit_20260908 as m
import validate_simplemos_compensated_step_20260908 as u
import calibrate_simplemos_srh_hotspots_20260908 as s
import check_simplemos_hotspot_contract_20260908 as q
import audit_simplemos_srh_response_remainder_20260908 as rem

a=m.a; d=m.d; p=m.p
OUT=m.OUT/'combined'
REPORT=p.REPO/'docs/validation/simplemos_stable_merit_step_and_srh_calibration_2026-09-08.md'


def yes(value):
    return value is True or value=='True'


def fmt(value):
    return f'{float(value):.9g}'


def key(row):
    return ('n19' if '_n19_' in row['key'] else 'n23')+' / '+('0.05' if '0p050000' in row['key'] else '1')


def table(headers, rows):
    return '\n'.join(['| '+' | '.join(headers)+' |', '| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,r))+' |' for r in rows])


def link(path,label):
    return f'[{label}](../../{a.rel(path)})'


def result(row):
    return ('通过' if yes(row['qualified']) else '失败')+'；'+fmt(row['max_row_ratio'])+'；'+str(row['iterations'])


def main():
    manifests=[m.OUT/'build_evidence.json',m.OUT/'freeze.json',m.OUT/'evidence.json',
               u.OUT/'build_evidence.json',u.OUT/'freeze.json',u.OUT/'evidence.json',
               s.OUT/'build_evidence.json',s.OUT/'freeze.json',s.OUT/'preflight_evidence.json',s.OUT/'evidence.json',
               q.OUT/'freeze.json',q.OUT/'guards_evidence.json',q.OUT/'evidence.json',rem.OUT/'freeze.json',rem.OUT/'evidence.json']
    for path in manifests:a.verify(path)
    mr=a.rows(m.OUT/'dc.csv');ur=a.rows(u.OUT/'dc.csv');sr=a.rows(s.OUT/'dc.csv')
    mc=a.rows(m.OUT/'comparison.csv');uc=a.rows(u.OUT/'comparison.csv');sp=a.rows(s.OUT/'preflight.csv');resp=a.rows(s.OUT/'response.csv')
    ms=a.read(m.OUT/'summary.json');us=a.read(u.OUT/'summary.json');ss=a.read(s.OUT/'summary.json');qs=a.read(q.OUT/'summary.json')
    remainder=a.rows(rem.OUT/'metrics.csv');remainder_blocks=a.rows(rem.OUT/'blocks.csv')
    assert len(mr)==24 and len(ur)==14 and len(sr)==20
    failure=[]
    for name,label in [('failure_original_restart','原始重启输入'),('failure_saved_state','重新加载失败输出')]:
        for dev in ('n19','n23'):
            rows=[r for r in mr if r['name']==name and '_'+dev+'_' in r['key']]
            failure.append([dev,label]+[result(next(r for r in rows if r['mode']==mode)) for mode in ('legacy','stable')])
    carry=[]
    for row in ur:
        if not yes(row['carry']):continue
        base=next(r for r in mr if all(r[x]==row[x] for x in ('key','name','mode')))
        carry.append([key(row),'原失败路径' if row['name'].startswith('failure') else '联合控制',row['mode'],result(base),result(row)])
    srhtable=[]
    for case in a.read(s.OUT/'contract.json')['cases']:
        rows=[r for r in sr if r['key']==case['key']]
        srhtable.append([key(rows[0])]+[result(next(r for r in rows if r['label']==label)) for label in ('zero','plus_full','minus_full','plus_half','minus_half')])
    response_table=table(['NWell / Vd','幅度','切向预测相对差','双幅度相对差','偶/奇响应','信号/零漂移','资格'],
        [[key(r),r['amplitude'],fmt(r['prediction_relative']),fmt(r['two_amplitude_relative']),fmt(r['even_over_odd']),fmt(r['signal_over_zero_drift']),'通过' if yes(r['qualified']) else '失败'] for r in resp])
    controls=[r for r in mc if r['name'].startswith('control')]
    carrycontrols=[r for r in uc if yes(r['carry']) and r['name'].startswith('control')]
    statuses=[]
    for module in (m,u):
        for job in a.read(module.OUT/'contract.json')['jobs']:
            dest=Path(job['dest']);status=a.read(dest/'config.status.json')
            rowgate=status.get('carrier_row_convergence',{})
            statuses.append(dict(branch='merit' if module is m else 'carry',key=job['case']['key'],name=job['name'],mode=job['mode'],carry=job.get('carry',False),
                iterations=status['iterations'],converged=status['converged'],convergence_reason=status.get('convergence_reason',''),failure_reason=status.get('failure_reason',''),elapsed_seconds=status['elapsed_seconds'],
                worst_node=rowgate.get('max_ratio_node'),worst_carrier=rowgate.get('max_ratio_carrier'),max_ratio=rowgate.get('max_ratio')))
    for case in a.read(s.OUT/'contract.json')['cases']:
        for job in case['jobs']:
            status=a.read(Path(job['dest'])/'config.status.json');rowgate=status.get('carrier_row_convergence',{})
            statuses.append(dict(branch='srh',key=case['key'],name=job['label'],mode='legacy',carry=False,
                iterations=status['iterations'],converged=status['converged'],convergence_reason=status.get('convergence_reason',''),failure_reason=status.get('failure_reason',''),elapsed_seconds=status['elapsed_seconds'],
                worst_node=rowgate.get('max_ratio_node'),worst_carrier=rowgate.get('max_ratio_carrier'),max_ratio=rowgate.get('max_ratio')))
    a.write_csv(OUT/'solver_status.csv',statuses)
    expected_errors=dict(wrong_mesh='Invalid finite SRH node manifest',wrong_count='Invalid finite SRH node manifest',
        missing_record='Invalid finite SRH volume record',duplicate='Invalid finite SRH volume record',zero_ratio='Invalid finite SRH volume record',
        out_of_range_node='Invalid finite SRH volume record',extra='Extra finite SRH data',bad_alpha='SRH amplitude outside frozen interpolation range',missing_alpha='Missing finite SRH amplitude')
    guard_reasons=[]
    for name,expected in expected_errors.items():
        message=(q.LOCAL/name/'config.stderr.txt').read_text(encoding='utf-8').strip()
        assert expected in message,(name,message)
        guard_reasons.append(dict(guard=name,expected=expected,actual=message,matched=True))
    a.write_csv(OUT/'guard_reasons.csv',guard_reasons)
    srhfail=[r for r in statuses if r['branch']=='srh' and not r['converged']]
    srh_linefail=sum(r['failure_reason']=='carrier_row_convergence_line_search_rejected' for r in srhfail)
    srh_iterfail=sum(r['failure_reason']=='carrier_row_convergence' and r['iterations']==200 for r in srhfail)
    residual_rejections=[]
    target=next(j for j in a.read(m.OUT/'contract.json')['jobs'] if j['case']['device']=='n23' and j['name']=='failure_original_restart' and j['mode']=='stable')
    trace=a.rows(Path(target['dest'])/'merit_trace.csv')
    if trace:
        # Preserve the last attempt using the actual trace schema, not a guessed iteration count.
        iter_key=next(k for k in trace[0] if k in ('iteration','iter'))
        last=max(int(r[iter_key]) for r in trace)
        residual_rejections=[r for r in trace if int(r[iter_key])==last]
        a.write_csv(OUT/'remaining_rejection.csv',residual_rejections)
    summary=dict(DC_total=58,stable=ms,compensated_step=us,srh=ss,srh_supplementary=qs,
        stable_control_strict_pairs=sum(yes(r['strict_pair_qualified']) for r in controls),stable_control_pairs=len(controls),
        carry_control_state_close=sum(yes(r['state_close']) for r in carrycontrols),carry_control_comparisons=len(carrycontrols),
        qualified_srh_remainder_max_relative=max(float(r['remaining_relative']) for r in remainder),
        production_changes=False,native_new_simulations=0,finite_srh_expansion=False)
    a.write(OUT/'summary.json',summary)
    source_error=max(float(r['source_relative']) for r in sp)
    weak_error=max(float(r['weak_J_relative']) for r in sp)
    linear_error=max(float(r['linear_relative']) for r in sp)
    port_max=max(float(r['relative']) for module in (m,u) for r in a.rows(module.OUT/'ports.csv'))
    body=f'''# SimpleMOS 稳定范数差、电势小更新及邻域 SRH 校准

日期：2026-09-08。完成 **58 次有界 DC 对照**、两组诊断单元测试（5 项、136 条断言）、四工况固定源/Jacobian 检查和九项输入拒绝检查。三组候选分别冻结输入及原接受门槛，未进行失败后的自动重启。

稳定范数差比较能修复部分原失败路径，但未解决高 NWell 的原始失败路径；简单累积被舍入吸收的 ψ 小更新使两条原失败路径变差。邻域 SRH 的源项与 Jacobian 插入检查通过，自洽响应仍有资格缺口。**本轮不放行生产算法替换、邻域 SRH 有限替换或完整曲线扩展。**

## 1. 配置、隔离及门槛

接续[上一轮精度与邻域审计](simplemos_linesearch_precision_and_srh_neighborhood_2026-09-08.md)。n19/n23 的 NWell Boron 为 1e17/2e17 cm⁻³，Vg=0.8 V，Vd=0.05/1 V。300 K，Boltzmann、无 BGN、总杂质 Masetti、掺杂 SRH；HFS、表面迁移率、Auger、雪崩、DG 关闭。沿用已验证的输运 element-box、独立 signed Si Poisson 电荷体积、单元介电系数及两节点联合 SRH 锚点。原失败路径只改变此前的第二 SRH 节点，不能与联合控制混作同一物理状态。

n19/n23 网格分别 1480/1482 节点、2742/2746 单元；各有 907 个自由 Si 节点，逐行检查全部 1814 个自由载流子行。坐标 μm，物理密度 m⁻³，Id A/μm，宽度 1 μm。Windows UCRT64/C++20 Release，实际 DC 后端 Eigen SparseLU/COLAMD，四次线性修正；HDF5/TDR、UMFPACK/SPQR 编译可用，本轮未切换后端。独立切向求解使用 SciPy SuperLU 和 long-double 残差修正。

沿用诊断程序的 q=1.602192e-19、kb=1.380662e-23、ε₀=8.8542e-12、Vt(300 K)=0.025851995266484913 V，未把手册末位数差当作本轮待修原因，也未修改生产常数。

所有 DC 保持 max_iter=200、rel_tol=1e-7、abs_tol=1e-12、stall residual floor=1e-9，ψ/准费米步长上限 0.35/0.025 V、contact_basin。逐行比值门槛 1e-6、KCL/Id 1e-8，独立全局闭合 eps=1e-6/source_floor=1e-10。求解时 global_continuity_closure.mode=off，线搜索为三块权重/尺度均为 1 的 L2 残差；全局闭合在输出态独立复核。合格还要求求解器返回 converged、退出码成功。状态对照门槛为势 1e-6 V、密度相对差 1e-4、Id 相对差 1e-6。没有排除新的自由载流子行。

## 2. 只替换范数差比较：24 次 DC

诊断比较器计算 Σ(Rtrial−Rbase)(Rtrial+Rbase)，先用补偿 long double；若符号不能由误差界确定，再用精确二进制整数累加平方差。只判断与原数学规则一致的严格 L2 下降，保留零对零例外，拒绝非有限输入。残差公式、Jacobian、状态形成、物理体积及资格门槛均不变；仅支持已经冻结的无加权 L2/无全局惩罚配置。

独立测试通过 3 项、96 条断言，覆盖双精度总范数看不见的小项下降/上升、极端指数、次正规数、精确消减、非有限输入，并复核上轮 26 个实际候选的高精度符号。两条关闭候选的原轨迹逐位复现检查：{ms['legacy_identity']}。

以下单元格依次为「资格；最差载流子行比值；迭代数」，不能只看总残差或 Id：

{table(['NWell','输入','原比较器','稳定比较器'],failure)}

从原始重启输入看，n19 获得改善，n23 仍失败。重新加载已导出的失败输出时，原比较器本身也能通过，因此这两次恢复不能归因于稳定比较器。n23 稳定分支最后一次尝试的 13 个实际候选平方差均为正；该次拒绝不是比较器丢失下降符号，详见 {link(OUT/'remaining_rejection.csv','最后拒绝记录')}。

稳定分支总资格 {ms['stable_qualified']}/{ms['stable_total']}，其中联合控制 {ms['stable_controls_qualified']}/{ms['stable_controls_total']}；原/稳定联合控制的严格状态配对 {summary['stable_control_strict_pairs']}/{len(controls)}。控制最大 Id 相对变化 {fmt(max(float(r['Id_relative']) for r in controls))}，最大势差 {fmt(max(float(r[k]) for r in controls for k in ('psi_max_V','phin_max_V','phip_max_V')))} V，最大密度相对变化 {fmt(max(float(r['density_max_relative']) for r in controls))}。

部分稳定分支走满 200 次迭代，按原有 max_iter_stall_residual_floor 接受并独立通过逐行门槛；这不是本轮新增接受条件，也不是迭代效率改善。全部退出原因保留于 {link(OUT/'solver_status.csv','求解器状态表')}。{link(m.OUT/'dc.csv','逐点 DC')}、{link(m.OUT/'comparison.csv','状态配对')}、{link(m.OUT/'trace_summary.csv','比较轨迹汇总')}。

## 3. 单独累积 ψ 小更新：14 次 DC

在同一初始状态上，分别以原/稳定比较器启用 ψ 更新余量：使用 TwoSum 与 FMA 保存被浮点加法吸收的部分，只在候选被接受后提交余量；拒绝候选不改变它。仅 ψ 有余量，载流子更新保持原样，每次 Newton 求解开始清零。

**这只是一次求解内的步长补偿，不是持久的电势参考值+增量表示，也没有提高 Poisson 残差求值精度。** 残差、Jacobian 和输出仍使用实际可表示的候选状态，余量不作为另一个物理状态导出。独立测试通过 2 项、40 条断言，对照 100 位累加并检查拒绝候选不污染已接受状态。

{table(['NWell / Vd','输入','比较器','不累积','累积 ψ 更新'],carry)}

累积候选共 {us['carry_qualified']}/{us['carry_total']} 合格：原失败路径 {us['carry_failures_qualified']}/{us['failure_cases']}、联合控制 {us['carry_controls_qualified']}/{us['control_cases']}。关闭累积的两条身份检查 {us['identity_passed']}/2。控制状态接近门槛通过 {summary['carry_control_state_close']}/{len(carrycontrols)}；该数字与 DC 资格分别记录，不能用接近但未闭合的状态替代合格状态。两条原失败路径均变差，不能将该简单累积方案迁入生产。{link(u.OUT/'dc.csv','逐点 DC')}、{link(u.OUT/'comparison.csv','同输入对照')}。

两种数值候选共 38 次端口对照中，独立 functional、contact extractor 和求解结果的最大相对差 {fmt(port_max)}，通过 {ms['ports_qualified']+us['ports_passed']}/38。端口一致只能排除本批次端口提取分支不一致，不能代替内部行闭合。

## 4. 四节点 SRH 体积方向：20 次 DC

保持原比较器和原 ψ 更新，沿用合格两节点联合基态；四个新增节点同时沿一个预定义几何方向移动：Vsrh(α)=Vold[1+α(VSi/Vold−1)]。n23 节点 794/795/1056/1087 对应 n19 793/794/1055/1086，空间映射最大偏差低于 1e-8 μm，与原两节点锚点不重叠。原锚点固定，Poisson 体积和输运几何保持不变。{link(s.OUT/'scope.csv','节点、坐标及几何体积比')}。

四工况固定状态检查全部通过：零扰动残差逐位相同，源公式相对差最大 {fmt(source_error)}（门槛 1e-5），弱交叉 Jacobian 相对差 {fmt(weak_error)}（1e-12），独立线性切向缺陷 {fmt(linear_error)}（1e-8），强 Jv 检查通过，边通量/迁移率/几何逐位不变。α=1 仅用于固定状态源/Jacobian 插入检查，没有进行有限自洽替换。九个错误清单/幅度输入均被拒绝，stderr 原因逐项核实为预期输入校验，非任意非零退出。{link(s.OUT/'preflight.csv','固定状态检查')}、{link(q.OUT/'guards.csv','输入拒绝检查')}、{link(OUT/'guard_reasons.csv','拒绝原因核实')}。

自洽批次每工况运行 α=0、±0.001、±0.0005，初始状态相同，不追加重启。单元格仍为「资格；最差行比值；迭代数」：

{table(['NWell / Vd','零扰动','+0.001','−0.001','+0.0005','−0.0005'],srhtable)}

DC 合格 {ss['qualified_DC']}/{ss['DC']}。响应使用全部 907 个自由 Si 节点的物理 ψ/φn/φp 向量；要求正负状态及零态都合格、切向预测差和双幅度差均 ≤1e-3、偶/奇 ≤0.01、信号/零漂移 ≥100。另独立要求零扰动 Id 漂移 ≤1e-5 dex。结果：

{response_table}

失败中 {srh_linefail} 点在载流子行尚未闭合时线搜索拒绝，{srh_iterfail} 点达到 200 次仍未通过载流子行门槛。最差空穴行集中于 n19 低 Vd 节点 794，以及 n19/n23 高 Vd 节点 1055/1056，详见求解器状态表；不能将这些退出标为软件崩溃或已校准的物理响应。

原响应资格 {ss['qualified_responses']}/{ss['responses']}；零扰动电流漂移通过 {qs['zero_drift_passed']}/4，九项输入保护通过，组合总资格为 {qs['qualified']}。{link(q.OUT/'summary.json','最终 SRH 资格')}。失败状态仍保留，仅作诊断，不用于宣称响应吻合；不能把固定状态公式正确当作自洽校准通过。本轮也未执行 Sentaurus 同源正负扰动，原生响应资格仍未取得。{link(s.OUT/'response.csv','响应数据')}、{link(s.OUT/'dc.csv','全部 DC 与失败原因')}。

## 5. 已合格 n23/低 Vd 的响应偏差进一步定位

该工况五个 DC 状态均合格，但全幅/半幅的切向预测相对差为 0.001047230271/0.001047230367，即约 **0.104723%**，略高于冻结的 0.1% 门槛。双幅度差仅 1.5728e-7，偶/奇响应和零漂移也通过；不能把该失败与其他未收敛点混为一谈。

对四个已有扰动态追加只读 functional 残差导出，没有追加 DC。令 Δxodd=(x+−x−)/2、Rodd=(R+−R−)/2，使用冻结基态的完整 J、相同源方向 s 和原切向 t=−J⁻¹s，检查误差 e=Δxodd−αt 与 J⁻¹Rodd。全部比较仍使用自由 Si 上的物理 ψ/φn/φp；这是残差误差的线性归因，不把修正后向量写回状态，也不重新标记原校准结果。

{table(['幅度','原响应相对差','剩余残差投影相对量','扣除投影后的相对差'],[[r['amplitude'],fmt(r['original_prediction_relative']),fmt(r['residual_projected_relative']),fmt(r['remaining_relative'])] for r in remainder])}

响应差几乎由空穴残差投影解释：全幅/半幅的空穴投影约 0.00104723119/0.00104723391，Poisson 投影约 2.50e-13/3.40e-13，电子投影约 4.47e-7/5.09e-7，均相对于预测响应范数。各块范数不能直接相加；沿实际误差方向的有符号贡献另列于 {link(rem.OUT/'blocks.csv','分块投影')}。原 J 线性解缺陷约 1.90e-14，残差投影求解缺陷低于 1.7e-16；扣除投影后的误差降至 9.90e-9/2.33e-8，支持当前局部切向与残差差分一致。

四个扰动均在 200 次以原 max_iter_stall_residual_floor 接受，最大行是节点 795 的电子行，约 1.078e-7/5.389e-8；但**响应误差由空穴残差主导**。这说明最大逐行比值与响应敏感方向不是同一个量。通过统一逐行 1e-6 仍可能不足以满足特定小源响应的 0.1% 精度，不能靠端口或最大行单一指标推断校准精度。

{link(rem.OUT/'metrics.csv','投影指标')}、{link(rem.OUT/'fields.csv','逐节点场误差与投影')}、{link(rem.OUT/'evidence.json','只读归因证据')}。本结果支持先处理剩余闭合精度，尚不能证明任意源方向、有限替换或原生同扰动响应已通过。

## 6. 结论与下一步边界

稳定范数差解决了一个可独立复现的数值障碍，但高 NWell 原始路径仍未闭合；简单 ψ 余量累积的反例说明，需要在残差与状态表示之间保持更完整的一致性。本轮不能据此断言完整的电势参考坐标方案失败，也不能声称已实现该方案。

下一候选应先独立验证 Poisson 残差的补偿求值/电势参考坐标，并用同一组原失败输入和联合控制检查实际候选残差、完整 Jacobian 及端口一致性；保持本轮失败记录和门槛。邻域 SRH 一方面应恢复未通过的正负小扰动态，另一方面需使 n23/低 Vd 已合格状态的残余空穴误差进一步下降，直到原始响应本身通过冻结门槛。不能用本轮只读残差投影代替自洽结果，也不应在失败态上拟合体积。高 NWell/高 Vd 反例未消除前，不进入全局默认或完整 0–1 V 曲线放行。

本轮只新增诊断头文件、独立 Catch2 测试、冻结构建/执行/分析脚本及证据，未修改生产 src/include、常数默认或模型默认。保留工作树已有生产修改。没有新增远程上传或 sdevice 仿真；未运行全量 CTest，不能将此前已有失败描述为本轮全绿。{link(OUT/'summary.json','汇总')}、{link(OUT/'final_evidence.json','最终证据')}。
'''
    REPORT.parent.mkdir(parents=True,exist_ok=True)
    assert not REPORT.exists()
    REPORT.write_text(body,encoding='utf-8')
    scripts=[p.REPO/'scripts'/name for name in (
        'build_simplemos_stable_merit_20260908.py','validate_simplemos_stable_merit_20260908.py',
        'build_simplemos_compensated_step_20260908.py','validate_simplemos_compensated_step_20260908.py',
        'calibrate_simplemos_srh_hotspots_20260908.py','check_simplemos_hotspot_contract_20260908.py','audit_simplemos_srh_response_remainder_20260908.py',Path(__file__).name)]
    for path in scripts:ast.parse(path.read_text(encoding='utf-8'),filename=str(path))
    links=[]
    for target in re.findall(r'\]\(([^)]+)\)',body):
        path=(REPORT.parent/target).resolve()
        if path==OUT/'final_evidence.json':continue
        assert path.exists(),str(path)
        links.append(path)
    a.write(OUT/'verification.json',dict(manifests_verified=len(manifests),python_syntax_files=len(scripts),local_links_checked=len(links),unit_tests=5,unit_assertions=136,DC=58))
    d.matrix.freeze(OUT/'final_evidence.json',manifests+scripts+[REPORT,OUT/'summary.json',OUT/'solver_status.csv',OUT/'remaining_rejection.csv',OUT/'guard_reasons.csv',OUT/'verification.json'])
    a.verify(OUT/'final_evidence.json')
    print(summary,flush=True)
    print(REPORT,flush=True)


if __name__=='__main__':main()
