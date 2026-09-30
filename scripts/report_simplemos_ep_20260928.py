"""Summarize measured EP/projection/initialization evidence, without new gates."""
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'build/simplemos_ep_20260928'
def rows(p):
    with p.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def main():
    comparisons=rows(OUT/'comparison.csv')
    summary=json.loads((OUT/'summary.json').read_text())
    checks=rows(OUT/'run_checks.csv')
    results=OUT/'probe_results'
    overlap=json.loads((results/'complete_cold_overlap.json').read_text())
    field=rows(OUT/'field_comparison.csv')
    logs=[]
    for c in checks:
        if c['exit_code']!='0':continue
        text=(OUT/'raw/bundle'/c['name'].split('_')[2]/(c['name']+'.console.log')).read_text()
        verified=all(s in text for s in ('Use 128 bit (double-double)', 'Relative error : 15 digits',
            'Mininum |rhs| : 1.0000e-15', 'Linear solver : blocked decomposition', 'Good Bye !'))
        assert verified,c['name']
        logs.append(dict(name=c['name'],precision_and_backend_verified=True,
            rhs_stop_count=text.count('|RHS| less than'),update_stop_count=text.count('Error smaller than')))
    (OUT/'precision_backend_checks.json').write_text(json.dumps(logs,indent=2)+'\n')
    summary['actual_precision_and_backend_require_log_review']=len(logs)!=8
    summary['logs_precision_and_reported_backend_verified']=len(logs)
    summary['max_abs_native_KCL_relative']=max(abs(float(r['default_KCL_relative'])) for r in comparisons)
    summary['max_abs_direct_default_Id_relative']=max(abs(float(r['direct_default_Id_relative'])) for r in comparisons)
    summary['max_abs_substrate_gap_relative']=max(abs(float(r['substrate_gap_relative'])) for r in comparisons)
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    residuals=[]
    for dev in ('n23','n24'):
        free=set(json.loads((ROOT/'build/outlier_analysis_20260928/inputs'/dev/'geometry.json').read_text())['free_si'])
        for arm in ('a','b','c'):
            rr=[r for r in rows(results/f'{dev}_{arm}_residual.csv') if int(r['node_id']) in free and float(r['net_doping_m3'])<0]
            residuals.append(dict(device=dev,arm=arm,nodes=len(rr),mask='free_Si_and_net_doping_negative',
                unit='raw_solver_scaled_residual_not_A_per_um',**{f'{k}_sum':math.fsum(float(r[k+'_residual']) for r in rr) for k in ('psi','phin','phip')}))
    (OUT/'auxiliary_p_region_residual_sums.json').write_text(json.dumps(residuals,indent=2)+'\n')
    lines=['# SimpleMOS 八次 EP 对照、三态探针和双初始化校核', '', '2026-09-28。用户批准后执行；原参考、2% 电流标准和原数值/双初始化门槛均未修改。', '',
        '## 原生高精度收敛包', '',
        f"已完成 {summary['completed_runs']}/8 次，形成 {summary['paired_points']} 个 default/Direct 配对偏置点。",
        '范围：n23/n24/n17/n19，Vd=0.05 V，各 Vg=0–2.5 V、间隔 0.05 V。',
        '保持 M60 网格、接触、300 K、Boltzmann、OldSlotboom、PhuMob、Enormal、HFS、SRH(DopingDependence) 和求解路径；只加入 ExtendedPrecision(128)、Digits=15、RhsMin=1e-15。n17 增加与其它器件一致的五份输出快照。',
        '8/8 Physics、Electrode 和除输出 Plot 外的 Solve 块一致。14/14 上传文件哈希匹配。',
        f'{len(logs)}/8 日志确认实际使用 double-double 128 bit、Digits=15、RhsMin=1e-15；后端回显为 blocked decomposition。日志中 SHEDistribution 的 Super 不能充当本 DD 方程内层后端的证明。未显式改动 Method。',
        '本次是三项设置共同改变的收敛包，不能单独量化 EP、Digits、RhsMin 各自贡献。停止条件仍可能是 RHS 或更新误差，不声称每个接受步都满足 RHS 阈值。', '',
        '| 器件 | Vg/V | Vela/M60−1，% | Vela/EP−1，% | EP KCL/abs(Id) | 衬底电子 default−Direct /abs(Id) |',
        '|---|---:|---:|---:|---:|---:|']
    for r in comparisons:
        if (r['device'] in ('n23','n24') and int(r['index'])<=2) or (r['device']=='n17' and int(r['index'])==3) or (r['device']=='n19' and int(r['index'])==1):
            err=100*(float(r['vela_Id'])/float(r['m60_default_Id'])-1)
            lines.append(f"| {r['device']} | {float(r['vg']):g} | {err:+.9f} | {float(r['vela_ep_error_percent']):+.9f} | {float(r['default_KCL_relative']):.4e} | {float(r['substrate_gap_relative']):.4e} |")
    lines+=['',f"当前配对曲线最大绝对 Vela/EP 电流误差为 {summary['max_abs_vela_ep_error_percent']:.9f}%。零观测差仅指 PLT 输出精度内相同，不宣称无限精度下严格为零。",
        f"全部配对点最大 abs(KCL)/abs(Id)={summary['max_abs_native_KCL_relative']:.6e}；最大 Direct/default 漏电流相对差={summary['max_abs_direct_default_Id_relative']:.6e}；最大衬底电子观测差/abs(Id)={summary['max_abs_substrate_gap_relative']:.6e}。",
        '该实验支持原生有限精度求值/收敛对此前差异的解释，未分离原生内部端口与状态项；不能把 KCL 当电流误差条，也不把四器件结论外推为全部 816 点 EP 验证。', '',
        '## 场导出检查', '',
        '每运行在 Vg=0/.05/.10/.15/.8 V 导出五份 TDR。同节点对比区域为全部 Si 节点（含接触）；速率、密度使用原生 cm 单位；势为 V。L1 指共同节点未加体积权的相对范数，不是源项积分。',
        '| 器件，Vg=.05 | EP−M60 max Δψ/V | max Δφn/V | max Δφp/V | n 的相对 L1 | SRH 的相对 L1 |',
        '|---|---:|---:|---:|---:|---:|']
    for dev in ('n23','n24'):
        ff={r['field']:r for r in field if r['device']==dev and float(r['vg'])==.05 and r['comparison']=='EP_default_minus_M60_default'}
        if ff:
            vals=[float(ff[k]['max_abs']) for k in ('ElectrostaticPotential','eQuasiFermiPotential','hQuasiFermiPotential')]
            vals += [float(ff[k]['unweighted_L1_relative']) for k in ('eDensity','srhRecombination')]
            lines.append('| '+dev+' | '+' | '.join(f'{x:.5e}' for x in vals)+' |')
    pairs=[r for r in field if r['comparison']=='EP_default_minus_EP_direct']
    lines += ['', f"已检验 {len(pairs)} 个字段/状态对的 EP default–Direct 最大绝对差为 {max(float(r['max_abs']) for r in pairs):.5e}（按 CSV 输出）。这些显示场近似不变，不等于原生内部状态位级相同。",'',
        '## 三态生产固定求值', '',
        '使用冻结 Linux 二进制 SHA256 `87b934c79451110467cee69b62f6ec2e312a8aea51ab85c07bf2c0577ad05ea3`，SparseLU/COLAMD 配置、统一 split 状态契约及冻结 HFS/几何/SRH 配置。',
        'n23/n24、Vd=Vg=.05 V，各 a/b/c × drain/source/substrate 共 18 次只读求值；全部成功 restore。a 为原始完整 split H5；b 从其物理 binary64 字段重新打包；c 为 M60 导出字段重新打包。并未以 EP 导出替换 c。',
        '打包用 C++ cpp_bin_float_100；物理字段和 QF 增量从最终 hi+lo 重建。两器件均逐节点核对坐标，兼容元数据保留，派生来源单独记录。a 的三端电流均通过 1e-12 相对复现检查，另外记录绝对差。17,784 个打包坐标另经 Python Fraction 精确有理数独立核对，hi/lo、物理场及增量的舍入全部一致；本地两个 split 状态头文件哈希与冻结云端源码一致。',
        '提取器内部约定 total=electron−hole；报告同时输出 hole_conventional=−hole_internal，与原生 electron+hole=total 对齐。', '',
        '| 器件 | (b−a) 漏电流 /abs(Id_a) | (c−b) 漏电流 /abs(Id_a) |',
        '|---|---:|---:|']
    for dev in ('n23','n24'):
        s=json.loads((results/f'{dev}_three_state_summary.json').read_text())
        d={r['meaning']:r for r in s['differences'] if r['contact']=='drain'}
        lines.append(f"| {dev} | {d['projection_repacking']['total_delta_over_original_drain']:.12g} | {d['export_precision_state_response']['total_delta_over_original_drain']:.12g} |")
    lines += ['', '**结论：投影态不具备重放原深关断接触电流的精度资格。** b 的漏电流丢失约99.75%，源/衬底基本不变；c 的接触电流偏移远大于真实 Id。它们是已测得的导出/投影状态求值结果，不能用于真实内部状态项归因，不能解释成模型误差，也不能用来替换自洽电流。',
        'p 型自由 Si 节点的辅助残差和已保存，但只是明确掺杂掩码下的 Vela 原始缩放残差，未转换为 A/µm，也未证明等同于 Sentaurus 原生阱域；不拿它与 default−Direct 电流作恒等式。', '',
        '## 双初始化与覆盖修正', '',
        '旧比较器构造低 Vd 目录时使用 0.05，但冻结目录使用 0p05，导致静默跳过低 Vd。保留第一次仅高 Vd 的16组输出；本轮独立比较器改为严格检查全部32组目录存在。生产求解器和历史报告未改。',
        f"最终 {overlap['qualified']}/{overlap['expected']} 组通过：n19/n23 × Vd=.05/1 × Vg=0/.2/.8/1 × baseline/native 两初始化。比较冷启动全曲线状态与既有两初始化控制态，不是新做816点双初始化。",
        '原门槛：三个势最大差≤1e-6 V、密度最大相对差≤1e-4、Id 相对差≤1e-6，同时两侧均已有数值资格。', '',
        '## 证据与未覆盖范围', '',
        '本地输出：`build/simplemos_ep_20260928/`。主要文件：`manifest.json`、`deck_changes.diff`、`uploaded.sha256`、`comparison.csv`、`terminal_components.csv`、`precision_backend_checks.json`、`field_comparison.csv`、`probe_results/*three_state_summary.json`、`probe_results/complete_cold_overlap.csv`。',
        '远端原生目录：`sentaurus:/tmp/vela_simplemos_ep_20260928`；云端独立目录：`/workspaces/simplemos-hfs-merged-20260926/ep_review_20260928`。Codespace 校核结束且无相关进程后已停止，状态确认 Shutdown。原始大状态保留，仅回传校核输出。',
        '未验证范围：n20/n21/n22 的 EP 原生场、其余 Vd 的 EP、全816点高精度参考、原生内部全精度状态跨程序重放。没有据此替换参考、放宽门槛或修改生产电流算法。','']
    (ROOT/'docs/validation/simplemos_ep_and_three_state_validation_2026-09-28.md').write_text('\n'.join(lines),encoding='utf-8')
    print(json.dumps(dict(completed_runs=summary['completed_runs'],report_written=True)))

if __name__=='__main__':main()
