"""Report completed precision and source-support audits with unchanged gates."""
import ast
import re
from decimal import Decimal as D
from pathlib import Path
import run_simplemos_linesearch_precision_20260908 as run
import audit_simplemos_srh_neighborhood_20260908 as neighbor

core=run.core; audit=core.audit; a=core.a; d=core.d; p=audit.p
OUT=audit.OUT
REPORT=p.REPO/'docs/validation/simplemos_linesearch_precision_and_srh_neighborhood_2026-09-08.md'


def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+['| '+' | '.join(map(str,row))+' |' for row in rows])


def main():
    manifests=[OUT/'build_evidence.json',OUT/'freeze.json',OUT/'replay_evidence.json',core.OUT/'runtime_evidence.json',core.OUT/'evidence.json',neighbor.OUT/'evidence.json']
    for path in manifests:a.verify(path)
    jobs=a.read(OUT/'contract.json')['jobs'];identity=a.rows(OUT/'identity.csv');merits=a.rows(core.OUT/'merit.csv');checks=a.rows(core.OUT/'checks.csv');linear=a.rows(core.OUT/'linear_defect.csv');coords=a.rows(core.OUT/'coordinates.csv')
    nn=a.rows(neighbor.OUT/'nodes.csv');nr=a.rows(neighbor.OUT/'rows.csv');nsummary=a.rows(neighbor.OUT/'summary.csv')
    assert len(checks)==28 and all(x['identity']=='True' and x['repeat_residual_identity']=='True' for x in identity)
    assert all(x['raw_equals_capped']=='True' for x in linear)
    effective=[];modecounts=[];normtable=[];steptable=[];holetable=[];roots=[]
    for job in jobs:
        key=job['case']['key'];dev=job['case']['device'];root=Path(job['dest']);base=a.rows(root/'base_merit.csv')[0]
        for label in ['base']+['trial_'+str(i) for i in range(13)]:
            m=a.rows(root/(label+'_merit.csv'))[0];alpha=float(m['alpha']);double=float(m['merit']);b=float(base['merit']);passed=label!='base' and (double<b or double<=(1-1e-4*alpha)*b)
            effective.append(dict(key=key,label=label,alpha=alpha,actual_double_merit=double,actual_double_relative_change=double/b-1,actual_line_search_would_decrease=passed))
        assert not any(x['actual_line_search_would_decrease'] for x in effective if x['key']==key)
        for mode in ('production','rounded','kernel','split','ideal'):
            rows=[x for x in merits if x['key']==key and x['mode']==mode]
            modecounts.append(dict(key=key,mode=mode,decreasing_trials=sum(x['would_decrease']=='True' for x in rows),total_trials=13,
                note='production means high-precision norm of exported production residual, not actual double norm'))
        for label in ('trial_0','trial_1'):
            values=[]
            for mode in ('kernel','split','ideal'):
                x=next(x for x in merits if x['key']==key and x['mode']==mode and x['label']==label);values.append(f"{100*float(x['relative_change']):+.8g}%")
            actual=next(x for x in effective if x['key']==key and x['label']==label)
            normtable.append([dev,actual['alpha'],f"{100*actual['actual_double_relative_change']:+.8g}%"]+values)
        for label in ('trial_0','trial_1','trial_12'):
            x=next(x for x in coords if x['key']==key and x['label']==label and x['field']=='psi')
            steptable.append([dev,2**(-int(label[6:])),x['nonzero_intended'],x['nonzero_actual'],f"{float(x['max_rounding_error_V']):.5e}"])
        for mode in ('rounded','kernel','split'):
            x=next(x for x in merits if x['key']==key and x['label']=='base' and x['mode']==mode)
            holetable.append([dev,job['node'],mode,f"{float(x['hole_row_ratio']):.10g}"])
    for n in nn:
        if n['role']!='root':continue
        # Use the carrier corresponding to the frozen phin/phip hotspot, never
        # an unresolved majority row to infer the source volume.
        carrier='electron' if n['vd']=='1.0' and ((n['device']=='n19' and n['node']=='793') or (n['device']=='n23' and n['node']=='794')) else 'hole'
        r=next(x for x in nr if x['key']==n['key'] and x['node']==n['node'] and x['carrier']==carrier and x['mode']=='qf')
        assert r['source_support_resolved']=='True'
        roots.append([n['device'],n['vd'],n['node'],carrier,f"{float(n['Si_over_current_volume']):.9f}",f"{float(r['current_volume_balance']):.6g}",f"{float(r['signed_Si_balance']):.6g}"])
    a.write_csv(OUT/'actual_line_search_decisions.csv',effective);a.write_csv(OUT/'precision_mode_counts.csv',modecounts)
    summary=dict(identity_replays=2,identity_replays_passed=2,snapshots=28,actual_trials=26,actual_decreasing_trials=0,
        max_precision_agreement=max(float(x['precision_agreement']) for x in checks),max_reconstruction=max(float(x['rounded_reconstruction']) for x in checks),
        max_norm_replay_relative=max(float(x['norm_replay_relative']) for x in checks),
        neighborhood_case_nodes=len(nn),native_row_representations=len(nr),source_resolved_representations=sum(x['source_support_resolved']=='True' for x in nr),
        qf_rows=len(nr)//2,qf_source_resolved=sum(x['source_support_resolved']=='True' and x['mode']=='qf' for x in nr),
        hotspot_carrier_rows_resolved=len(roots),original_failed_states_still_unqualified=2,
        independent_global_acceptance_unchanged=True,global_closure_in_line_search='off in both actual solver configurations',
        production_changes=False,new_native_DC=False,new_finite_volume_substitution=False,new_DC_qualification=False,
        conclusion='Residual evaluation, scalar norm comparison and quantization of the actual state update all affect last-step acceptance; full linear solve and clipping are not the limiting causes at these two states. All six selected hotspot source rows support signed Si integration in fixed native exports. Native same-source DC response and expanded finite substitution remain unqualified.')
    a.write(OUT/'summary.json',summary)
    report=f"""# SimpleMOS 线搜索残差精度与剩余 SRH 热点邻域审计

日期：2026-09-08。本轮完成两项只读验证：两个失败过程逐位重现，并对基态及实际 13 个回溯候选做 60/100 位精度扫描；四工况共 38 个热点/邻域节点、152 个原生行表达完成体积对照。**剩余失败同时受残差求值、标量范数比较和状态更新舍入影响；六个预选热点的相应载流子行均支持 signed Si SRH 体积。**

本轮未修改生产求解器、常数、接受门槛或 SRH 策略，没有新增有限替换或 Sentaurus 仿真。两个原失败态在高精度下仍未满足逐行 1e-6 门槛，不能把只读诊断的下降候选当作新的合格 DC 状态。

## 1. 范围与复现资格

接续[有限自洽验证](simplemos_constants_and_srh_finite_validation_2026-09-08.md)。材料、网格和配置保持该报告原样：n19/n23（NWell Boron 1e17/2e17 cm⁻³）、Vg=0.8 V、Vd=0.05/1 V；300 K、Boltzmann/no-BGN、Masetti、掺杂 SRH，HFS/表面/Auger/雪崩/DG 关闭。网格 n19/n23 分别 1480/1482 节点，每例 907 个自由 Si 节点、1814 个自由载流子行，坐标 μm，物理密度 m⁻³，Id A/μm、宽度 1 μm。沿用手册常数的隔离程序、已验证自动 box/介电/输运配置。

线搜索对象仅为 Vd=1、只改第二节点、原生来源初始化的一次重启失败态；分别是 n19 节点 1056、n23 节点 1057。隔离程序只增加输出，链接同一常数约定的 core 对象。Windows UCRT64/C++20 Release、Eigen SparseLU/COLAMD、四次线性修正均不变。HDF5/TDR、UMFPACK、SPQR 编译可用，本轮未切换后端。

从原重启输入重放后，输出状态 SHA256、接受迭代数（81/18）、退出码和失败原因全部相同；最终分别在第 82/19 次尝试拒绝。28 个快照的重复残差逐位一致。导出完整 x、原始/截断步、Jacobian、边通量、SRH 核函数输入、Poisson 系数/电荷/边界及实际范数。[复现表](../../reference_tcad/simplemos_sentaurus2022/linesearch_srh_audit_20260908/identity.csv)。

严格逐行 1e-6、端口 KCL/Id 1e-8 和独立全局闭合资格不变。特别核实：**这两个实际配置中的 global_continuity_closure.mode=off，线搜索没有加入全局源和惩罚。** 先前的全局闭合是独立资格检查，不应把其 source_floor 误认为当前 merit 的放大因子。实际 residual_norm 为 block，三块 scale/weight 都为 1，因此当前标量为三块残差平方和开根号。公式及实际开关见[原判定回放](../../reference_tcad/simplemos_sentaurus2022/linesearch_srh_audit_20260908/actual_line_search_decisions.csv)。

## 2. 精度层次及校验

所有高精度计算以导出值对应的精确二进制浮点数为输入，固定双精度几何、材料参数、派生系数和热电压。分别检查：

1. CSV `production`：只用高精度重算已导出的生产残差范数；**不是实际双精度范数**。实际值另存 actual_line_search_decisions.csv。
2. `rounded`：使用已舍入的边通量、Poisson 电荷乘积及 SRH 速率，提高累加、SRH 积分乘积和缩放的精度。初始合同对 source products 的简写不表示 SRH 整个乘积已经预先舍入；具体操作以此说明和冻结分析代码为准。
3. `kernel`：按原核函数输入高精度重算 Poisson 乘积/差、SG、SRH，保持已经转换到核函数坐标的输入值。
4. `split`：从**实际已舍入候选 x** 高精度形成物理电势和参考值+增量，再重算密度、SG、SRH、Poisson。
5. `ideal`：同 split，但用高精度形成 x₀+α·step，不发生双精度更新吸收；这是数学方向对照，未写回生产状态。

各表示用自身基态比较，避免跨表示的绝对残差混比。60/100 位结果的最大全行尺度归一化差为 {summary['max_precision_agreement']:.6g}（门槛 1e-20）；已舍入项对原残差重建最大差 {summary['max_reconstruction']:.6g}（1e-10）；原范数公式回放相对差最大 {summary['max_norm_replay_relative']:.6g}。不同表示间可有显著残差差异，高精度位数一致只说明参考计算稳定，不能证明其已成为新的生产离散方程。

分析脚本第一次读取时因共享 CSV 接口要求 Path 而停止，未产生数值结论。随后用单独冻结的路径类型适配器运行，原脚本与合同保留，未改公式或阈值。[路径适配记录](../../reference_tcad/simplemos_sentaurus2022/linesearch_srh_audit_20260908/precision/path_compatibility/adapter.json)、[精度检查](../../reference_tcad/simplemos_sentaurus2022/linesearch_srh_audit_20260908/precision/path_compatibility/checks.csv)。

## 3. 失败原因的分层定位

原线搜索接受条件为范数严格下降，或满足 (1−1e-4·α) 的充分下降；26 个实际候选均未通过，全部载流子有效性检查为真。下面比较同一方向，百分比为相对各自基态的范数变化：

{table(['NWell','α','实际双精度','kernel 高精度','split 实际候选高精度','ideal 未舍入方向'],normtable)}

两点的完整原始步与截断步在所有列完全相同，不仅最差行相同。100 位 J·dx+R 的 L2 分别约 4.10e-27、3.36e-27，最差空穴行线性缺陷/该行原残差分别约 2.43e-15、6.50e-16。因此，当前更新上限与线性方程求解精度不能解释这两个最后拒绝事件。

残差求值精度确实会改变判定：例如 n19 半步在实际范数中增加 0.09067%，split 重算则下降 0.38416%；n23 完整步由增加 4.76193% 变为下降 6.37632%。但只提高求和精度不足以恢复所有数学方向，核函数坐标形成与实际状态舍入也有影响。

电势更新吸收很明显：

{table(['NWell','α','拟更新 ψ 节点数','实际改变 ψ 节点数','更新舍入最大差/V'],steptable)}

n19 完整步在 split 仍增大，却在 ideal 下降 16.2801%；n23 ideal 完整步下降 26.6167%。两点的 13 个 ideal 步均下降，实际 ψ 经双精度加法后很多更新消失。进一步缩小 α 会丢失更多 ψ 更新，因此单纯加深回溯无法保证恢复该方向。

另一个独立因素是范数分辨率：n19 的 α=1/128 与 n23 的 α=1/4096，生产双精度总范数与基态相等；仅用 100 位重算同一组生产残差的范数，却分别能看到约 −8.86e-34、−7.02e-33 的相对下降。它来自极少数载流子残差改善，被约 1e-10 的 Poisson 主导范数淹没。该数值是微小数学下降，不表示一次足以越过逐行资格门槛。

最差空穴行的基态重算仍未闭合：

{table(['NWell','节点','表示','行比值（门槛 1e-6）'],holetable)}

高精度求值不能把这两个失败态直接重新标为合格。[全部范数扫描](../../reference_tcad/simplemos_sentaurus2022/linesearch_srh_audit_20260908/precision/path_compatibility/merit.csv)、[线性缺陷](../../reference_tcad/simplemos_sentaurus2022/linesearch_srh_audit_20260908/precision/path_compatibility/linear_defect.csv)、[状态更新](../../reference_tcad/simplemos_sentaurus2022/linesearch_srh_audit_20260908/precision/path_compatibility/coordinates.csv)、[最差 Poisson 行](../../reference_tcad/simplemos_sentaurus2022/linesearch_srh_audit_20260908/precision/path_compatibility/worst_poisson.csv)。

## 4. 剩余热点邻域的 SRH 体积

使用此前合格的两节点联合有限状态，在四工况分别冻结 7、11、7、13 个节点：高 Vd 的 n19 793/794、n23 794/1056，低 Vd 的 n19 1086、n23 1087 及各自一圈相邻自由 Si 节点，同时保留原先两体积节点作为锚点。38 为按工况计的节点数，包含跨偏置重复位置。

原生回放从原生场和独立 box 单元加权 Masetti 出发，保留准费米势式与密度式 SG，比较原生体积率乘当前积分体积或 signed Si 体积；未以 Vela 场替代原生场。预先要求源信号/(64ε×绝对边通量和)≥100、两种 SG 净流出差/源≤1e-4、signed Si 闭合≤1e-5，才将该行标为源体积可分辨。强流出消减的多数载流子行保留为未分辨，未据其异常推断物理体积。

152 个行表达中 105 个满足这一诊断条件，准费米势式为 53/76；所有六个预选热点的相应载流子行均满足。表内体积比来自几何；闭合指标分母统一为 max(|净流出|, |原生 SRH×signed Si 体积|)。

{table(['NWell','Vd/V','热点节点','载流子','signed Si/当前体积','当前体积不闭合','signed Si 不闭合'],roots)}

这些位置的未替换体积差与剩余准费米势热点一致，支持把邻域 SRH 支撑作为下一候选。两个已替换锚点在每个工况的体积比为 1，提供 8 个按工况计的体积不变控制。原生多数载流子净通量很小时，双表达差与源尺度比仍可能很大；不能以端口电流吻合否定这些未分辨记录，也不能把它们用于反推体积。

本表是**固定原生导出场的源定义审计**，尚未得到原生同源正负自洽响应资格；没有执行邻域扩大后的有限替换。Vela 联合基态仍保留原严格资格，当前体积策略的自由行最大比值低于 1e-6。[节点/场/体积](../../reference_tcad/simplemos_sentaurus2022/linesearch_srh_audit_20260908/neighborhood/nodes.csv)、[全部行与可分辨性](../../reference_tcad/simplemos_sentaurus2022/linesearch_srh_audit_20260908/neighborhood/rows.csv)。

## 5. 下一步应隔离验证的实现

先验证数值上等价的**范数差比较**：从同一原始方向比较残差平方和的增量，避免将微小载流子下降加到较大 Poisson 范数后丢失；保留原逐行及端口资格。这是最小候选，但不能单独解决 ψ 更新吸收，应先用本轮两失败态及合格控制做独立对照。

另行验证 Poisson 残差的补偿求值及电势参考值+增量表示，分别检查残差/Jacobian/端口一致性和完整步是否保留。不能把高精度理想步直接写回双精度后宣称已实现，也不能同时改范数、状态和物理体积以至失去归因。

物理方向应先形成几何定义的 SRH 支撑候选，在上述热点邻域做独立小幅正负源响应校准，覆盖两个 NWell、两个 Vd 后再有限自洽对照。保留现有 Poisson/输运独立体积；不得用全局 mesh 体积开关替代 SRH 独立策略。原生同源响应、Vg=1、16 工况和完整曲线仍待后续验证。

本轮验证包括两次原轨迹逐位复现、28 个双精度/60 位/100 位快照、全行重建与范数公式校验、完整 Jacobian 线性缺陷、152 个原生源行表达和证据哈希检查。无生产代码变更，未重跑全量 CTest。脚本：[只读重放](../../scripts/audit_simplemos_linesearch_precision_20260908.py)、[高精度分析](../../scripts/analyze_simplemos_linesearch_precision_20260908.py)、[运行适配](../../scripts/run_simplemos_linesearch_precision_20260908.py)、[邻域审计](../../scripts/audit_simplemos_srh_neighborhood_20260908.py)、[报告生成](../../scripts/report_simplemos_linesearch_srh_audit_20260908.py)。[汇总](../../reference_tcad/simplemos_sentaurus2022/linesearch_srh_audit_20260908/summary.json)、[最终证据](../../reference_tcad/simplemos_sentaurus2022/linesearch_srh_audit_20260908/final_evidence.json)。
"""
    assert not REPORT.exists();REPORT.write_text(report,encoding='utf-8')
    scripts=[Path(x.__file__).resolve() for x in (audit,core,run,neighbor)]+[Path(__file__).resolve()]
    for script in scripts:ast.parse(script.read_text(encoding='utf-8'),filename=str(script))
    links=re.findall(r'\]\(([^)]+)\)',report)
    for link in links:
        if link.endswith('/final_evidence.json'):continue
        assert (REPORT.parent/link).resolve().exists(),link
    a.write(OUT/'report_checks.json',dict(syntax_checked=5,local_links_checked=len(links),stage_hashes_verified=len(manifests),production_defaults_changed=False))
    d.matrix.freeze(OUT/'final_evidence.json',manifests+scripts+[REPORT,OUT/'summary.json',OUT/'report_checks.json',OUT/'actual_line_search_decisions.csv',OUT/'precision_mode_counts.csv'])
    a.verify(OUT/'final_evidence.json');print(summary,flush=True);print(a.rel(REPORT),flush=True)


if __name__=='__main__':main()
